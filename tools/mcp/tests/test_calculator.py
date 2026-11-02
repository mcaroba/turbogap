"""The calculator's lifecycle on a stand-in `turbogap ipi`: an i-PI client written
from the protocol whose potential is E = |r|^2 / 2 in atomic units, so the units
coming back can be checked against a formula. Here and through a stand-in ssh."""

import os
import sys

import pytest

from turbogap_mcp import jobs
from turbogap_mcp.calculator import CalcError, Registry
from turbogap_mcp.config import HostConfig
from turbogap_mcp.hosts import make_host

BOHR = 0.5291772109
HARTREE = 27.2113862460

FAKE_IPI = f"""#!{sys.executable}
import re, socket, struct, sys, time
text = open("input").read()
behaviour = text.splitlines()[0]
if "never connect" in behaviour:
    sys.exit(2)
name = re.search(r"ipi_address = 'UNIX:(\\w+)'", text).group(1)
s = socket.socket(socket.AF_UNIX)
s.connect("/tmp/ipi_" + name)
def recvn(n):
    b = b""
    while len(b) < n:
        c = s.recv(n - len(b))
        if not c:
            sys.exit(5)
        b += c
    return b
init = have = False
while True:
    head = recvn(12).decode().strip()
    if head == "STATUS":
        s.sendall((b"HAVEDATA" if have else b"READY" if init else b"NEEDINIT").ljust(12))
    elif head == "INIT":
        recvn(4); recvn(struct.unpack("<i", recvn(4))[0]); init = True
    elif head == "POSDATA":
        recvn(144); nat = struct.unpack("<i", recvn(4))[0]
        pos = struct.unpack(f"<{{3 * nat}}d", recvn(24 * nat))
        if "crash" in behaviour:
            sys.exit(4)
        have = True
    elif head == "GETFORCE":
        s.sendall(b"FORCEREADY".ljust(12))
        s.sendall(struct.pack("<d", 0.5 * sum(x * x for x in pos)))
        s.sendall(struct.pack("<i", nat))
        s.sendall(struct.pack(f"<{{3 * nat}}d", *(-x for x in pos)))
        s.sendall(struct.pack("<9d", *range(9)))
        s.sendall(struct.pack("<i", 0))
        have = False
    elif head == "EXIT":
        sys.exit(0)
"""

# sshd runs the remote command through the user's shell. A forward (-N) is the
# same machine here, so the socket is already where the host looks for it.
FAKE_SSH = """#!/bin/bash
for a; do [ "$a" = -N ] && exec sleep 3600; done
for last; do :; done
exec bash -c "$last"
"""

INPUT = """atoms_file = 'atoms.xyz'
pot_file = 'pot.gap'
n_species = 2
species = C O
masses = 12.01 15.99
"""


def frame(positions, species=("C", "O"), cell="6 0 0 0.5 7 0 0.3 0.2 8"):
    lines = [str(len(species)), f'Lattice="{cell}" Properties=species:S:1:pos:R:3']
    lines += [f"{s} {x} {y} {z}" for s, (x, y, z) in zip(species, positions)]
    return "\n".join(lines) + "\n"


START = frame([(0.1, 0.2, 0.3), (1.2, 0.0, 0.1)])


def executable(path, text):
    path.write_text(text)
    path.chmod(0o755)
    return str(path)


@pytest.fixture(params=["local", "ssh"])
def host(request, tmp_path, monkeypatch):
    monkeypatch.setenv("TURBOGAP_MCP_STATE", str(tmp_path / "state"))
    config = HostConfig(name=request.param, bin=executable(tmp_path / "fake_ipi", FAKE_IPI),
                        runs_root=str(tmp_path / "runs"), timeout=20)
    if request.param == "ssh":
        config.ssh = "nowhere"
        config.ssh_command = [executable(tmp_path / "fake_ssh", FAKE_SSH)]
    return make_host(config)


@pytest.fixture
def registry():
    reg = Registry()
    yield reg
    reg.stop_all()


def start(registry, host, table, tree, first_line=""):
    return registry.start(host, table, tree, first_line + INPUT, START,
                          files={"pot.gap": "gap_beg distance_2b\ngap_end\n"}, start_timeout=20, eval_timeout=20)


def test_evaluations_in_ev_and_angstrom(registry, host, table, tree):
    started = start(registry, host, table, tree)
    calc = registry.get(started["calc_id"])
    for positions, cell in [([(0.5, -1.0, 2.0), (1.5, 0.25, 0.0)], "6 0 0 0.5 7 0 0.3 0.2 8"),
                            ([(0.4, -1.1, 2.2), (1.4, 0.35, 0.1)], "6.2 0 0 0.4 7.1 0 0.3 0.1 8.3")]:
        out = calc.evaluate(frame(positions, cell=cell))
        bohr = [[x / BOHR for x in r] for r in positions]
        assert out["energy"] == pytest.approx(HARTREE * 0.5 * sum(x * x for r in bohr for x in r), rel=1e-14)
        assert sum(out["forces"], []) == pytest.approx([-x * HARTREE / BOHR for r in bohr for x in r], rel=1e-14)
        assert out["virial"] == pytest.approx([HARTREE * k for k in range(9)], rel=1e-14)
    stopped = registry.stop(started["calc_id"])
    assert (stopped["state"], stopped["exit_code"], stopped["evaluations"]) == ("done", 0, 2)
    assert not os.path.exists(calc.server.path)


def test_other_atoms_are_refused_and_the_calculator_survives(registry, host, table, tree):
    calc = registry.get(start(registry, host, table, tree)["calc_id"])
    with pytest.raises(CalcError, match="species"):
        calc.evaluate(frame([(0, 0, 0), (1, 0, 0)], species=("O", "C")))
    with pytest.raises(CalcError, match="Lattice"):
        calc.evaluate("2\nno cell\nC 0 0 0\nO 1 0 0\n")
    assert calc.evaluate(START)["energy"] > 0


def test_a_crash_mid_evaluation_ends_the_calculator(registry, host, table, tree):
    calc_id = start(registry, host, table, tree, "! crash\n")["calc_id"]
    calc = registry.get(calc_id)
    with pytest.raises(CalcError, match="exit 4"):
        calc.evaluate(START)
    with pytest.raises(CalcError, match="stopped"):
        calc.evaluate(START)
    assert host.poll(calc.job.dir)[0] in ("done", "cancelled")


def test_a_client_that_never_connects_fails_the_start_fast(registry, host, table, tree):
    with pytest.raises(CalcError, match="exit 2") as caught:
        start(registry, host, table, tree, "! never connect\n")
    assert "before connecting" in str(caught.value)
    assert not [p for p in os.listdir("/tmp") if p.startswith("ipi_mcp_") and
                os.stat(f"/tmp/{p}").st_uid == os.getuid() and not registry.calcs]


def test_ipi_address_is_the_calculators_to_set(registry, host, table, tree):
    with pytest.raises(CalcError, match="ipi_address"):
        registry.start(host, table, tree, INPUT + "ipi_address = 'UNIX:mine'\n", START)


def test_invalid_input_starts_nothing(registry, host, table, tree):
    out = registry.start(host, table, tree, INPUT + "md_stpe = 1\n", START)
    assert out["started"] is False
    assert not jobs.all_jobs()
