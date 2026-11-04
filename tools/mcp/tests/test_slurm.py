"""Jobs on a Slurm host, against a stand-in Slurm: sbatch runs the batch script at
once (or once a hold file goes), squeue and sacct answer from what it left, and a
test can make sacct say how a job ended. Here and through the stand-in ssh."""

import os
import time

import pytest
from test_hosts import ATOMS, FAKE_BIN, FAKE_SSH, INPUT, executable

from turbogap_mcp import jobs
from turbogap_mcp.config import HostConfig, load_hosts
from turbogap_mcp.calculator import CalcError, Registry
from turbogap_mcp.hosts import HostError, make_host

FAKES = {
    "sbatch": """#!/bin/bash
for last; do :; done
n=$(cat {S}/next 2>/dev/null || echo 100); echo $((n + 1)) > {S}/next
out=$(sed -n 's/^#SBATCH --output=//p' "$last")
setsid bash -c "while [ -f {S}/hold ]; do sleep 0.1; done; touch {S}/$n.started
bash '$last' >> '$out' 2>&1; echo \\$? > {S}/$n.code" > /dev/null 2>&1 < /dev/null &
echo $! > {S}/$n.pid
cp "$last" {S}/$n.script
echo "$n;cluster"
""",
    "squeue": """#!/bin/bash
while [ $# -gt 0 ]; do case $1 in -j) id=$2; shift ;; -o) fmt=$2; shift ;; esac; shift; done
if ! kill -0 "$(cat {S}/$id.pid 2>/dev/null)" 2>/dev/null; then
    echo "slurm_load_jobs error: Invalid job id specified" >&2; exit 1; fi
state=RUNNING; reason=None
[ -f {S}/$id.started ] || {{ state=PENDING; reason=Priority; }}
if [ "$fmt" = %T ] || [ -z "$fmt" ]; then echo $state; else echo "$state|$reason|node1|0:05"; fi
""",
    "sacct": """#!/bin/bash
while [ $# -gt 0 ]; do case $1 in -j) id=$2; shift ;; esac; shift; done
if [ -f {S}/$id.sacct ]; then line=$(cat {S}/$id.sacct)
elif [ -f {S}/$id.cancelled ]; then line="CANCELLED by 1|0:15"
elif [ -f {S}/$id.code ]; then c=$(cat {S}/$id.code)
    if [ "$c" = 0 ]; then line="COMPLETED|0:0"; else line="FAILED|$c:0"; fi
else exit 0; fi
echo "$line|node1|00:00:05"
""",
    "scancel": """#!/bin/bash
touch {S}/$1.cancelled
kill -TERM -- -$(cat {S}/$1.pid) 2>/dev/null || true
""",
    "srun": """#!/bin/bash
exec "$@"
""",
}


@pytest.fixture
def slurm_dir(tmp_path):
    state = tmp_path / "slurm"
    state.mkdir()
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    for name, text in FAKES.items():
        executable(bindir / name, text.replace("{S}", str(state)).replace("{{", "{").replace("}}", "}"))
    return state, bindir


@pytest.fixture(params=["local", "ssh"])
def host(request, tmp_path, monkeypatch, slurm_dir):
    monkeypatch.setenv("TURBOGAP_MCP_STATE", str(tmp_path / "state"))
    _, bindir = slurm_dir
    config = HostConfig(name=request.param, bin=executable(tmp_path / "fake_turbogap", FAKE_BIN),
                        runs_root=str(tmp_path / "runs"), max_ranks=4, omp_threads=2, timeout=20,
                        setup=f"export PATH={bindir}:$PATH; echo noise from setup",
                        slurm={"account": "acct", "partition": "small", "time": "00:05:00",
                               "options": ["--mem-per-cpu=1G"]})
    if request.param == "ssh":
        config.ssh = "nowhere"
        config.ssh_command = [executable(tmp_path / "fake_ssh", FAKE_SSH)]
    return make_host(config)


def start(host, table, tree, first_line="", ranks=2):
    out = jobs.start(host, table, tree, mode="predict", input_text=first_line + INPUT, atoms_xyz=ATOMS,
                     files={"pot.gap": "gap_beg distance_2b\ngap_end\n"}, ranks=ranks)
    return jobs.load(out["job_id"])


def wait(host, job, until=("done", "failed", "cancelled", "lost"), seconds=20):
    for _ in range(int(seconds * 10)):
        st = jobs.status(host, job)
        if st["state"] in until:
            return st
        time.sleep(0.1)
    raise AssertionError(f"still {st['state']}")


def slurm_id(job):
    return open(f"{job.dir}/slurm_id").read().strip()


def test_queued_then_done_with_the_batch_script_asked_for(host, table, tree, slurm_dir):
    state, _ = slurm_dir
    (state / "hold").touch()
    job = start(host, table, tree)
    st = wait(host, job, until=("queued",))
    assert (st["slurm_id"], st["slurm_state"], st["reason"]) == ("100", "PENDING", "Priority")
    (state / "hold").unlink()
    st = wait(host, job)
    assert (st["state"], st["exit_code"], st["slurm_state"]) == ("done", 0, "COMPLETED")
    assert "mode predict threads 2" in st["out_log_tail"]
    script = (state / "100.script").read_text()
    for line in ["#SBATCH --account=acct", "#SBATCH --partition=small", "#SBATCH --time=00:05:00",
                 "#SBATCH --ntasks=2", "#SBATCH --cpus-per-task=2", f"#SBATCH --output={job.dir}/out.log",
                 "#SBATCH --mem-per-cpu=1G", "export OMP_NUM_THREADS=2", "echo $? > rc"]:
        assert line in script.splitlines(), line
    assert jobs.results(host, job)["frame"]["energy"] == -1.5


def test_a_failure_keeps_its_code(host, table, tree):
    st = wait(host, start(host, table, tree, "! fake fail\n"))
    assert (st["state"], st["exit_code"]) == ("failed", 3)


def test_scancel_while_queued_and_running(host, table, tree, slurm_dir):
    state, _ = slurm_dir
    (state / "hold").touch()
    queued = start(host, table, tree, "! fake sleep\n")
    wait(host, queued, until=("queued",))
    assert jobs.cancel(host, queued)["state"] == "cancelled"
    (state / "hold").unlink()
    running = start(host, table, tree, "! fake sleep\n")
    wait(host, running, until=("running",))
    assert jobs.cancel(host, running)["state"] == "cancelled"
    assert not os.path.exists(f"{running.dir}/rc")


def test_a_job_slurm_ended_is_failed_with_slurms_reason(host, table, tree, slurm_dir):
    state, _ = slurm_dir
    job = start(host, table, tree, "! fake sleep\n")
    wait(host, job, until=("running",))
    n = slurm_id(job)
    (state / f"{n}.sacct").write_text("TIMEOUT|0:0")
    os.killpg(int((state / f"{n}.pid").read_text()), 9)
    st = wait(host, job)
    assert (st["state"], st["slurm_state"]) == ("failed", "TIMEOUT")


def test_no_calculator_on_a_batch_host(host, table, tree):
    with pytest.raises(CalcError, match="login node"):
        Registry().start(host, table, tree, INPUT, ATOMS)


def test_an_expired_certificate_is_named(tmp_path):
    key = tmp_path / "id"
    os.system(f"ssh-keygen -q -t ed25519 -N '' -f {key}")
    os.system(f"ssh-keygen -q -s {key} -I test -n me -V 20200101:20200102 {key}.pub")
    config = HostConfig(name="roihu", bin="x", runs_root="/x", ssh="roihu",
                        certificate=f"{key}-cert.pub", renew="roihucer")
    with pytest.raises(HostError, match=r"expired at 2020-01-02.*renew it with `roihucer`"):
        make_host(config).sh("true")


def test_slurm_config_needs_account_partition_and_time(tmp_path):
    path = tmp_path / "hosts.toml"
    path.write_text('[r]\nbin = "x"\nruns_root = "/r"\n[r.slurm]\naccount = "a"\npartition = "p"\n')
    with pytest.raises(ValueError, match="needs time"):
        load_hosts(path)
    path.write_text('[r]\nbin = "x"\nruns_root = "/r"\n[r.slurm]\naccount = "a"\npartition = "p"\n'
                    'time = "1:00"\nqos = "x"\n')
    with pytest.raises(ValueError, match="qos"):
        load_hosts(path)
