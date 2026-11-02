"""Real TurboGAP runs through the job layer, against the regression suite's own run
of the same case with the same binary: the outputs must be byte-identical.

TURBOGAP_MCP_TEST_BIN names the binary (default: the tree's bin/turbogap) and
TURBOGAP_DATA_ROOT the turbogap_tests checkout (default: beside the tree). Either
missing skips these tests, loudly.
"""

import filecmp
import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest

from turbogap_mcp import examples, jobs
from turbogap_mcp.config import HostConfig
from turbogap_mcp.hosts import make_host

CASES = ["co_predict", "co_predict_mpi4", "co_md"]


def binary(tree):
    return Path(os.environ.get("TURBOGAP_MCP_TEST_BIN", tree / "bin" / "turbogap"))


def data_root(tree):
    return Path(os.environ.get("TURBOGAP_DATA_ROOT", tree.parent / "turbogap_tests"))


@pytest.fixture(scope="module")
def setup(tree):
    if not os.access(binary(tree), os.X_OK):
        pytest.skip(f"no TurboGAP binary at {binary(tree)}: set TURBOGAP_MCP_TEST_BIN")
    if not (data_root(tree) / "CO").is_dir():
        pytest.skip(f"no test data at {data_root(tree)}: set TURBOGAP_DATA_ROOT")
    if shutil.which("mpirun") is None:
        pytest.skip("no mpirun")


def suite_run(tree, case, work):
    """The suite refuses a reference that is the binary under test, so it gets a copy.
    Threads as the job has them (omp_threads = 1): the sums depend on the count."""
    work.mkdir(parents=True)
    reference = work / "reference_turbogap"
    shutil.copy2(binary(tree), reference)
    env = dict(os.environ, OMP_NUM_THREADS="1", TURBOGAP_KEEP="1", TMPDIR=str(work), TURBOGAP_BIN=str(binary(tree)),
               TURBOGAP_REF_BIN=str(reference), TURBOGAP_DATA_ROOT=str(data_root(tree)))
    done = subprocess.run([str(tree / "tests" / "regression" / "run.sh"), "--cpu", case],
                          cwd=tree, env=env, capture_output=True, text=True, timeout=900)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    kept = sorted(work.glob(f"turbogap_regression.*/{case}.test"))
    assert kept, f"the suite kept no staging directory under {work}"
    return kept[0]


@pytest.mark.parametrize("case", CASES)
def test_case_matches_the_suite(setup, tree, table, tmp_path, monkeypatch, case):
    monkeypatch.setenv("TURBOGAP_MCP_STATE", str(tmp_path / "state"))
    host = make_host(HostConfig(name="local", bin=str(binary(tree)), runs_root=str(tmp_path / "runs"),
                                test_data=str(data_root(tree)), max_ranks=4, timeout=60))
    started = jobs.start(host, table, tree, case=case)
    assert started["started"], started
    job = jobs.load(started["job_id"])
    for _ in range(900):
        st = jobs.status(host, job)
        if st["state"] != "running":
            break
        time.sleep(1)
    assert (st["state"], st["exit_code"]) == ("done", 0), st

    reference = suite_run(tree, case, tmp_path / "suite")
    outputs = examples_outputs(tree, case)
    for name in outputs:
        assert filecmp.cmp(reference / name, Path(job.dir) / name, shallow=False), name

    frames, energy = last_frame_energy(reference / "trajectory_out.xyz")
    res = jobs.results(host, job)
    assert res["frames_walked"] == frames
    assert res["frame"]["energy"] == energy


def examples_outputs(tree, case):
    return examples.read_case(tree, case)["outputs"]


def last_frame_energy(path):
    """Frame count and the last frame's energy, read without the server's code."""
    lines = path.read_text().splitlines()
    i = 0
    frames = 0
    header = ""
    while i < len(lines):
        header = lines[i + 1]
        i += int(lines[i]) + 2
        frames += 1
    return frames, float(header.split("energy=")[1].split()[0])
