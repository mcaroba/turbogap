"""The calculator on real TurboGAP, checked the way tests/ipi_socket/run.sh checks
`turbogap ipi`: two geometries in a triclinic cell, without and with a supercell,
against `turbogap predict` on the geometries actually sent, to the precision of
trajectory_out.xyz; and the transposed cell must disagree.

TURBOGAP_MCP_TEST_BIN names the binary and TURBOGAP_DATA_ROOT the turbogap_tests
checkout; either missing skips these tests, loudly.
"""

import importlib.util
import os
import subprocess
import sys

import numpy as np
import pytest

from turbogap_mcp.calculator import Registry, load_refserver
from turbogap_mcp.config import HostConfig
from turbogap_mcp.hosts import make_host

COMMON = """pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345
"""


def binary():
    return os.environ.get("TURBOGAP_MCP_TEST_BIN", "")


def co(tree):
    return os.path.join(os.environ.get("TURBOGAP_DATA_ROOT", str(tree.parent / "turbogap_tests")), "CO")


def module(tree, name):
    path = tree / "tests" / "ipi_socket" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"ipi_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def configs(tree, tmp_path_factory):
    if not os.access(binary(), os.X_OK):
        pytest.skip(f"no TurboGAP binary at '{binary()}': set TURBOGAP_MCP_TEST_BIN")
    if not os.path.isdir(os.path.join(co(tree), "gap_files")):
        pytest.skip(f"no CO potential under {co(tree)}: set TURBOGAP_DATA_ROOT")
    work = tmp_path_factory.mktemp("configs")
    subprocess.run([sys.executable, str(tree / "tests" / "ipi_socket" / "make_configs.py"),
                    os.path.join(co(tree), "atoms_897.xyz")], cwd=work, check=True, capture_output=True)
    return work


def predict(tree, work, xyz_text, tag):
    """`turbogap predict` on the given frames, run directly, not through the server."""
    d = work / f"predict_{tag}"
    d.mkdir()
    os.symlink(os.path.join(co(tree), "gap_files"), d / "gap_files")
    (d / "atoms.xyz").write_text(xyz_text)
    (d / "input").write_text(COMMON + 'atoms_file = "atoms.xyz"\nwrite_xyz = 1\n')
    done = subprocess.run([binary(), "predict"], cwd=d, capture_output=True, text=True,
                          env=dict(os.environ, OMP_NUM_THREADS="1"))
    assert done.returncode == 0, done.stdout[-2000:]
    return module(tree, "compare").read_traj(d / "trajectory_out.xyz")


def effective(ref, path, out):
    species, pos, cell = ref.read_extxyz(path)
    ref.write_effective(out, species, pos, cell)
    return out.read_text()


def transposed(xyz_text):
    """The same text with the lattice matrix transposed, as run.sh makes its control."""
    lines = xyz_text.splitlines()
    head, rest = lines[1].split('Lattice="', 1)
    lattice, tail = rest.split('"', 1)
    cell = np.array(lattice.split()).reshape(3, 3)
    lines[1] = head + 'Lattice="' + " ".join(cell.T.ravel()) + '"' + tail
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("tag", ["big", "small"])
def test_socket_path_agrees_with_predict(tree, table, configs, tmp_path, monkeypatch, tag):
    monkeypatch.setenv("TURBOGAP_MCP_STATE", str(tmp_path / "state"))
    host = make_host(HostConfig(name="local", bin=binary(), runs_root=str(tmp_path / "runs"),
                                read_roots=[co(tree)], timeout=60))
    ref = load_refserver(tree)
    tol = module(tree, "compare")
    frames = [configs / f"{tag}_{k}.xyz" for k in range(2)]

    registry = Registry()
    try:
        started = registry.start(host, table, tree, COMMON + 'atoms_file = "atoms.xyz"\n',
                                 frames[0].read_text(),
                                 links=[{"source": os.path.join(co(tree), "gap_files"), "target": "gap_files"}])
        assert started["started"], started
        calc = registry.get(started["calc_id"])
        got = [calc.evaluate(f.read_text()) for f in frames]
        log = (tmp_path / "runs" / started["calc_id"] / "out.log").read_text()
        stopped = registry.stop(started["calc_id"])
    finally:
        registry.stop_all()
    assert (stopped["state"], stopped["exit_code"]) == ("done", 0), log[-2000:]
    assert ("smaller" in log) == (tag == "small"), "the supercell is built in the small cell only"

    sent = "".join(effective(ref, f, tmp_path / f"eff_{k}.xyz") for k, f in enumerate(frames))
    expected = predict(tree, tmp_path, sent, tag)
    assert len(expected) == len(frames)
    for k, (E, f, v) in enumerate(expected):
        assert abs(got[k]["energy"] - E) <= tol.TOL_E, k
        assert np.abs(np.array(got[k]["forces"]) - f).max() <= tol.TOL_F, k
        assert np.abs(np.array(got[k]["virial"]).reshape(3, 3) - v).max() <= tol.TOL_V, k

    if tag == "big":
        (_, fT, _), = predict(tree, tmp_path, transposed((tmp_path / "eff_0.xyz").read_text()), "T")
        assert np.abs(np.array(got[0]["forces"]) - fT).max() > tol.DISCRIMINATE
