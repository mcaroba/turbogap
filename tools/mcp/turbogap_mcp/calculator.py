"""The live calculator: `turbogap ipi` kept running behind an i-PI socket, so each
geometry costs one force call and not a potential load.

The server side is tests/ipi_socket/refserver.py, imported, because that is the
tested reading of the protocol. TurboGAP connects to /tmp/ipi_<name>; on a remote
host `ssh -R` puts that path there and forwards it to the listener here.
"""

import importlib.util
import math
import secrets
import shlex
import socket
import subprocess
import threading
import time

import numpy as np

from . import jobs
from .hosts import Host, HostError
from .results import MAX_ARRAY_ATOMS, columns, parse_header
from .validate import parse_entries

POLL_S = 0.5
FORWARD_WAIT_S = 15.0
# How long a client whose socket broke gets to exit, so its code and log can be reported.
DYING_WAIT_S = 10.0


class CalcError(ValueError):
    pass


def load_refserver(root):
    path = root / "tests" / "ipi_socket" / "refserver.py"
    spec = importlib.util.spec_from_file_location("turbogap_refserver", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_frame(text):
    """Species, positions (A) and the cell with its vectors as rows, from one extxyz frame."""
    lines = text.strip("\n").splitlines()
    try:
        n = int(lines[0].split()[0])
        header = parse_header(lines[1])
        cell = np.array([float(x) for x in header["Lattice"].split()]).reshape(3, 3)
        cols = columns(header.get("Properties", "species:S:1:pos:R:3"))
        s, _ = cols["species"]
        p, _ = cols["pos"]
        atoms = [line.split() for line in lines[2:2 + n]]
        species = [a[s] for a in atoms]
        positions = np.array([[float(x) for x in a[p:p + 3]] for a in atoms])
    except (IndexError, KeyError, ValueError) as e:
        raise CalcError(f"not an extxyz frame with a Lattice, species and pos: {e}") from None
    if len(species) != n:
        raise CalcError(f"the frame says {n} atoms and has {len(species)}")
    return species, positions, cell


def with_address(input_text, name):
    entries, _ = parse_entries(input_text)
    if any(e.keyword == "ipi_address" for e in entries):
        raise CalcError("leave ipi_address out: the calculator sets it")
    return input_text.rstrip("\n") + f"\nipi_address = 'UNIX:{name}'\n"


class Forward:
    """ssh -R from /tmp/ipi_<name> on the host to the same path here."""

    def __init__(self, host: Host, path):
        c = host.config
        self.host = host
        self.path = path
        self.proc = subprocess.Popen(
            c.ssh_command + ["-N", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8",
                             "-o", "ControlMaster=no", "-o", "ControlPath=none",
                             "-o", "ExitOnForwardFailure=yes", "-o", "StreamLocalBindUnlink=yes",
                             "-R", f"{path}:{path}", c.ssh],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        q = shlex.quote(path)
        tries = int(FORWARD_WAIT_S / 0.2)
        try:
            host.sh(f"for i in $(seq {tries}); do [ -S {q} ] && exit 0; sleep 0.2; done; exit 1\n",
                    timeout=FORWARD_WAIT_S + 10)
        except HostError:
            self.close()
            raise CalcError(f"ssh -R to {host.name} did not come up: {self.stderr()}") from None

    def stderr(self):
        if self.proc.poll() is None:
            return "still running, but the socket never appeared"
        return self.proc.stderr.read().decode(errors="replace").strip()[-300:]

    def close(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(5)
            except subprocess.TimeoutExpired:
                self.proc.kill()
        try:
            self.host.sh(f"rm -f {shlex.quote(self.path)}\n")
        except HostError:
            pass


class Calculator:
    def __init__(self, host: Host, job: jobs.Job, server, forward, species, eval_timeout):
        self.host = host
        self.job = job
        self.server = server
        self.forward = forward
        self.species = species
        self.eval_timeout = eval_timeout
        self.lock = threading.Lock()
        self.evaluations = 0
        self.closed = False

    def handshake(self, start_timeout):
        """Wait for TurboGAP to connect and take INIT. It answers its first STATUS only
        after the evaluation it discards, so this includes loading the potential."""
        deadline = time.monotonic() + start_timeout
        self.server.listener.settimeout(POLL_S)
        while self.server.conn is None:
            try:
                self.server.conn, _ = self.server.listener.accept()
            except socket.timeout:
                self.check_alive("before connecting")
                if time.monotonic() > deadline:
                    raise CalcError(f"TurboGAP did not connect within {start_timeout:.0f} s")
        self.server.conn.settimeout(start_timeout)
        self.protocol(lambda: (self.server.expect_status("NEEDINIT", "before INIT"),
                               self.server.init(payload=b"turbogap-mcp")))

    def check_alive(self, when):
        state, code = self.host.poll(self.job.dir)
        if state != "running":
            tail = self.host.tail(f"{self.job.dir}/out.log", 15)
            raise CalcError(f"TurboGAP is {state} (exit {code}) {when}:\n{tail}")

    def protocol(self, step):
        """A protocol step; any failure leaves the stream out of step, so it ends the calculator."""
        try:
            return step()
        except (SystemExit, RuntimeError, OSError) as e:
            detail = str(e) or type(e).__name__
            for _ in range(int(DYING_WAIT_S / POLL_S)):
                try:
                    self.check_alive("after a failed exchange")
                except CalcError as dead:
                    detail = f"{detail}\n{dead}"
                    break
                time.sleep(POLL_S)
            self.close(kill=True)
            raise CalcError(f"calculator {self.job.id} ended: {detail}") from None

    def evaluate(self, xyz, arrays=True):
        species, positions, cell = parse_frame(xyz)
        if species != self.species:
            raise CalcError(f"the species list must be the start frame's, atom for atom "
                            f"({len(self.species)} atoms); i-PI sends positions only")
        with self.lock:
            if self.closed:
                raise CalcError(f"calculator {self.job.id} is stopped")
            self.server.conn.settimeout(self.eval_timeout)
            t0 = time.perf_counter()

            def exchange():
                self.server.expect_status("READY", "before a geometry")
                self.server.posdata(cell, positions)
                self.server.expect_status("HAVEDATA", "after a geometry")
                return self.server.getforce()

            energy, forces, virial = self.protocol(exchange)
            self.evaluations += 1
        return summarise(energy, forces, virial, time.perf_counter() - t0, arrays)

    def close(self, kill=False):
        if self.closed:
            return
        self.closed = True
        if not kill and self.server.conn is not None:
            try:
                self.server.conn.settimeout(10)
                self.server.expect_status("READY", "before EXIT")
                self.server.exit()
            except (SystemExit, RuntimeError, OSError):
                kill = True
        self.server.close()
        if self.forward is not None:
            self.forward.close()
        if kill:
            self.host.kill(self.job.dir)
            return
        for _ in range(int(30 / POLL_S)):
            if self.host.poll(self.job.dir)[0] != "running":
                return
            time.sleep(POLL_S)
        self.host.kill(self.job.dir)


def summarise(energy, forces, virial, seconds, arrays):
    norms = np.linalg.norm(forces, axis=1)
    out = {"energy": float(energy), "virial": [float(v) for v in virial.ravel()],
           "force_max": float(norms.max()) if len(norms) else 0.0,
           "force_rms": float(math.sqrt((norms ** 2).mean())) if len(norms) else 0.0,
           "wall_s": round(seconds, 4), "units": "eV, eV/A, eV"}
    if arrays and len(forces) <= MAX_ARRAY_ATOMS:
        out["forces"] = forces.tolist()
    elif arrays:
        out["arrays_omitted"] = f"{len(forces)} atoms is over {MAX_ARRAY_ATOMS}"
    return out


class Registry:
    """The calculators of this server process. They end with it: TurboGAP aborts
    when its socket closes."""

    def __init__(self):
        self.calcs = {}
        self.lock = threading.Lock()

    def start(self, host: Host, table, root, input_text, atoms_xyz, links=(), files=None, ranks=1,
              start_timeout=600.0, eval_timeout=600.0):
        species, _, _ = parse_frame(atoms_xyz)
        name = f"mcp_{secrets.token_hex(4)}"
        text = with_address(input_text, name)
        server = load_refserver(root).Server(name)
        forward = None
        try:
            if host.config.ssh:
                forward = Forward(host, server.path)
            started = jobs.start(host, table, root, mode="ipi", input_text=text, atoms_xyz=atoms_xyz,
                                 links=links, files=files, ranks=ranks)
        except BaseException:
            server.close()
            if forward is not None:
                forward.close()
            raise
        if not started["started"]:
            server.close()
            if forward is not None:
                forward.close()
            return started
        calc = Calculator(host, jobs.load(started["job_id"]), server, forward, species, eval_timeout)
        try:
            calc.handshake(start_timeout)
        except CalcError:
            calc.close(kill=True)
            raise
        with self.lock:
            self.calcs[calc.job.id] = calc
        return {"started": True, "calc_id": calc.job.id, "host": host.name, "dir": calc.job.dir,
                "n_atoms": len(species), "warnings": started["warnings"]}

    def get(self, calc_id):
        with self.lock:
            calc = self.calcs.get(calc_id)
        if calc is None:
            raise CalcError(f"no running calculator '{calc_id}' in this server")
        return calc

    def stop(self, calc_id):
        calc = self.get(calc_id)
        calc.close()
        with self.lock:
            self.calcs.pop(calc_id, None)
        state, code = calc.host.poll(calc.job.dir)
        return {"calc_id": calc_id, "evaluations": calc.evaluations, "state": state, "exit_code": code}

    def stop_all(self):
        with self.lock:
            calcs = list(self.calcs.values())
            self.calcs.clear()
        for calc in calcs:
            try:
                calc.close()
            except Exception:
                calc.close(kill=True)
