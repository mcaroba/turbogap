"""Job mechanics on a stand-in binary, here and through a stand-in ssh that hands
the remote command string to a shell the way sshd does. Fast; no TurboGAP needed."""

import os
import subprocess
import time

import pytest

from turbogap_mcp import jobs
from turbogap_mcp.config import HostConfig, load_hosts
from turbogap_mcp.hosts import HostError, make_host

FAKE_BIN = """#!/bin/bash
case "$(head -1 input)" in
  "! fake sleep") sleep 300 ;;
  "! fake fail") echo boom; exit 3 ;;
  *) printf '2\\nLattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3:forces:R:3 energy=-1.5\\n' > trajectory_out.xyz
     printf 'C 0 0 0 1 0 0\\nO 1.2 0 0 -1 0 0\\n' >> trajectory_out.xyz
     printf '#     Step             Time      Temperature                E_kin                E_pot             Pressure\\n' > thermo.log
     printf '0 0.0 300.0 0.1 -1.5 2.0\\n1 1.0 310.0 0.2 -1.4 2.1\\n' >> thermo.log
     echo "mode $1 threads $OMP_NUM_THREADS" ;;
esac
"""

# sshd runs the remote command through the user's shell; so does this.
FAKE_SSH = """#!/bin/bash
for last; do :; done
exec bash -c "$last"
"""

INPUT = """atoms_file = 'atoms.xyz'
pot_file = 'pot.gap'
n_species = 2
species = C O
"""
ATOMS = "2\nLattice=\"5 0 0 0 5 0 0 0 5\"\nC 0 0 0\nO 1.2 0 0\n"


def executable(path, text):
    path.write_text(text)
    path.chmod(0o755)
    return str(path)


@pytest.fixture(params=["local", "ssh"])
def host(request, tmp_path, monkeypatch):
    monkeypatch.setenv("TURBOGAP_MCP_STATE", str(tmp_path / "state"))
    config = HostConfig(name=request.param, bin=executable(tmp_path / "fake_turbogap", FAKE_BIN),
                        runs_root=str(tmp_path / "runs with space"), max_ranks=2, omp_threads=3,
                        read_roots=[str(tmp_path / "shared")], timeout=20)
    if request.param == "ssh":
        config.ssh = "nowhere"
        config.ssh_command = [executable(tmp_path / "fake_ssh", FAKE_SSH)]
    return make_host(config)


def start(host, table, tree, first_line="", **kw):
    return jobs.start(host, table, tree, mode="predict", input_text=first_line + INPUT,
                      atoms_xyz=ATOMS, files={"pot.gap": "gap_beg distance_2b\ngap_end\n"}, **kw)


def wait(host, job, until=("done", "failed", "cancelled", "lost"), seconds=20):
    for _ in range(int(seconds * 10)):
        st = jobs.status(host, job)
        if st["state"] in until:
            return st
        time.sleep(0.1)
    raise AssertionError(f"still {st['state']}")


def test_run_to_completion(host, table, tree):
    started = start(host, table, tree)
    job = jobs.load(started["job_id"])
    st = wait(host, job)
    assert (st["state"], st["exit_code"]) == ("done", 0)
    assert "mode predict threads 3" in st["out_log_tail"]
    res = jobs.results(host, job)
    assert res["frame"]["energy"] == -1.5
    assert res["frame"]["species"] == {"C": 1, "O": 1}
    assert res["frame"]["force_max"] == 1.0
    assert res["thermo"]["last"]["Temperature"] == 310.0


def test_failure_is_reported_with_its_code(host, table, tree):
    job = jobs.load(start(host, table, tree, "! fake fail\n")["job_id"])
    st = wait(host, job)
    assert (st["state"], st["exit_code"]) == ("failed", 3)
    assert "boom" in st["out_log_tail"]


def test_cancel_leaves_no_process(host, table, tree):
    job = jobs.load(start(host, table, tree, "! fake sleep\n")["job_id"])
    assert wait(host, job, until=("running",))["state"] == "running"
    group = (open(f"{job.dir}/pid").read().strip())
    assert jobs.cancel(host, job)["state"] == "cancelled"
    left = subprocess.run(["pgrep", "-g", group], capture_output=True, text=True).stdout
    assert left == ""


def test_a_killed_job_is_lost_not_running(host, table, tree):
    job = jobs.load(start(host, table, tree, "! fake sleep\n")["job_id"])
    wait(host, job, until=("running",))
    group = open(f"{job.dir}/pid").read().strip()
    os.killpg(int(group), 9)
    assert wait(host, job)["state"] == "lost"


def test_invalid_input_starts_nothing(host, table, tree):
    out = jobs.start(host, table, tree, mode="predict", input_text=INPUT + "md_stpe = 1\n")
    assert out["started"] is False
    assert not os.path.exists(host.config.runs_root)


def test_missing_named_file_refuses_to_launch(host, table, tree):
    with pytest.raises(jobs.JobError, match="not there"):
        jobs.start(host, table, tree, mode="predict", input_text=INPUT, atoms_xyz=ATOMS)


def test_links_outside_the_roots_are_refused(host, table, tree, tmp_path):
    for source in ["/etc/passwd", str(tmp_path / "shared" / ".." / "x"), "relative"]:
        with pytest.raises(jobs.JobError, match="outside"):
            start(host, table, tree, links=[{"source": source, "target": "x"}])


def test_links_inside_a_root_are_made(host, table, tree, tmp_path):
    (tmp_path / "shared").mkdir()
    (tmp_path / "shared" / "data.dat").write_text("1 2\n")
    job = jobs.load(start(host, table, tree, links=[{"source": str(tmp_path / "shared" / "data.dat"),
                                                     "target": "exp.dat"}])["job_id"])
    assert open(f"{job.dir}/exp.dat").read() == "1 2\n"


def test_ranks_are_capped(host, table, tree):
    with pytest.raises(jobs.JobError, match="max_ranks"):
        start(host, table, tree, ranks=3)


def test_put_is_verbatim(host, tmp_path):
    text = "a 'quoted' \"line\" with $HOME and `cmd` and TURBOGAP_MCP_EOF\n\\n end\n"
    target = str(tmp_path / "put.txt")
    host.put(text, target)
    assert open(target).read() == text


def test_a_command_reading_stdin_does_not_eat_the_script(host):
    assert host.sh("cat\necho after\n").strip() == "after"


def test_unreachable_host_times_out(tmp_path):
    config = HostConfig(name="dead", bin="x", runs_root="/x", ssh="dead", timeout=1,
                        ssh_command=[executable(tmp_path / "hang", "#!/bin/bash\nsleep 30\n")])
    with pytest.raises(HostError, match="no answer within 1 s"):
        make_host(config).sh("true")


def test_config_rejects_unknown_keys(tmp_path):
    path = tmp_path / "hosts.toml"
    path.write_text('[a]\nbin = "x"\nruns_root = "/r"\nmax_rank = 3\n')
    with pytest.raises(ValueError, match="max_rank"):
        load_hosts(path)
