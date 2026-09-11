#!/usr/bin/env python3
"""An i-PI SERVER, written from the protocol, to drive `turbogap ipi`.

This is the reference the test compares against, and the point of it is that
it is not TurboGAP. It was written from the wire protocol -- 12-byte
space-padded ASCII headers, little-endian doubles and 32-bit integers, atomic
units -- as the i-PI documentation and the reference Fortran driver state it,
not from src/ipi_socket.f90. If the two agree on a triclinic cell, in atomic
units, over two consecutive geometries, then the client's reading of the
protocol is not merely self-consistent.

What it does:

  1. Listens on a UNIX socket and waits for TurboGAP to connect.
  2. Runs the handshake, ASSERTING the answers: NEEDINIT before INIT, READY
     after INIT and before coordinates, HAVEDATA once forces exist. A client
     that answered HAVEDATA too early would hand back forces for the geometry
     in atoms_file rather than the one it was sent, which is a wrong
     trajectory and not a crash, so the order is checked rather than assumed.
  3. Sends each geometry in bohr with the cell as i-PI holds it -- lattice
     vectors in the COLUMNS of h, sent row-major -- and reads back energy,
     forces and virial in Hartree and Hartree/bohr.
  4. Writes what came back, converted to eV and Angstrom, for run.sh to
     compare against `turbogap predict` on the same geometries.
  5. Writes the geometry it ACTUALLY SENT, back-converted to Angstrom at full
     precision. Angstrom -> bohr -> Angstrom is not the identity in floating
     point, and on this system the last-bit difference moves a force by up to
     1e-7 eV/A -- ten times the file format's own granularity, which would
     otherwise have to be absorbed by a tolerance nobody could justify.
     Feeding `turbogap predict` the effective geometry instead makes the two
     paths see bit-identical input, so what is left is the print format and
     nothing else.

Usage:  refserver.py SOCKETNAME CONFIG.xyz [CONFIG2.xyz ...] --out PREFIX
"""
from __future__ import annotations

import argparse
import os
import socket
import struct
import sys

import numpy as np

MSGLEN = 12
BOHR = 0.5291772109          # Angstrom per bohr
HARTREE = 27.2113862460      # eV per Hartree


def read_extxyz(path):
    """Species, positions in Angstrom, and the 3x3 cell with vectors as ROWS."""
    with open(path) as fh:
        nat = int(fh.readline().split()[0])
        comment = fh.readline()
        key = 'Lattice="'
        i = comment.index(key) + len(key)
        lat = np.array([float(x) for x in comment[i:comment.index('"', i)].split()])
        cell_rows = lat.reshape(3, 3)
        sp, pos = [], []
        for _ in range(nat):
            f = fh.readline().split()
            sp.append(f[0])
            pos.append([float(f[1]), float(f[2]), float(f[3])])
    return sp, np.array(pos), cell_rows


def write_effective(path, species, pos_ang, cell_rows_ang):
    """The geometry the client will reconstruct: (x/BOHR)*BOHR, not x.

    17 significant digits so the text round trip is itself exact -- a
    shorter format would put back the very error this is here to remove.
    """
    pos = (pos_ang / BOHR) * BOHR
    cell = (cell_rows_ang / BOHR) * BOHR
    with open(path, "w") as fh:
        fh.write(f"{len(species)}\n")
        fh.write('Lattice="' + " ".join(f"{v:.17g}" for v in cell.ravel())
                 + '" Properties=species:S:1:pos:R:3\n')
        for s, r in zip(species, pos):
            fh.write(f"{s} " + " ".join(f"{v:.17g}" for v in r) + "\n")


class Server:
    def __init__(self, name):
        self.path = f"/tmp/ipi_{name}"
        if os.path.exists(self.path):
            os.unlink(self.path)
        self.listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.listener.bind(self.path)
        self.listener.listen(1)
        self.conn = None

    def accept(self, timeout):
        self.listener.settimeout(timeout)
        self.conn, _ = self.listener.accept()
        self.conn.settimeout(timeout)

    def close(self):
        if self.conn is not None:
            self.conn.close()
        self.listener.close()
        if os.path.exists(self.path):
            os.unlink(self.path)

    # -- wire ------------------------------------------------------------
    def put_header(self, text):
        assert len(text) <= MSGLEN
        self.conn.sendall(text.ljust(MSGLEN).encode("ascii"))

    def get_header(self):
        return self.recvn(MSGLEN).decode("ascii").strip()

    def recvn(self, n):
        buf = b""
        while len(buf) < n:
            chunk = self.conn.recv(n - len(buf))
            if not chunk:
                raise RuntimeError("the client closed the connection")
            buf += chunk
        return buf

    def status(self):
        self.put_header("STATUS")
        return self.get_header()

    def expect_status(self, want, when):
        got = self.status()
        if got != want:
            raise SystemExit(f"FAIL: STATUS {when}: expected {want}, got '{got}'")
        return got

    def init(self, bead=0, payload=b"reference server"):
        self.put_header("INIT")
        self.conn.sendall(struct.pack("<i", bead))
        self.conn.sendall(struct.pack("<i", len(payload)))
        self.conn.sendall(payload)

    def posdata(self, cell_rows_ang, pos_ang):
        # i-PI's h holds the lattice vectors as COLUMNS and is sent row-major.
        # cell_rows_ang holds them as rows, as extxyz writes them, so h is its
        # transpose. Getting this backwards is invisible for a cubic cell,
        # which is why the test cell is triclinic.
        h = (cell_rows_ang.T / BOHR).astype(np.float64)
        self.put_header("POSDATA")
        self.conn.sendall(np.ascontiguousarray(h).tobytes())
        self.conn.sendall(np.ascontiguousarray(np.linalg.inv(h)).tobytes())
        self.conn.sendall(struct.pack("<i", len(pos_ang)))
        self.conn.sendall(np.ascontiguousarray(pos_ang / BOHR, dtype=np.float64).tobytes())

    def getforce(self):
        self.put_header("GETFORCE")
        head = self.get_header()
        if head != "FORCEREADY":
            raise SystemExit(f"FAIL: expected FORCEREADY, got '{head}'")
        pot = struct.unpack("<d", self.recvn(8))[0]
        nat = struct.unpack("<i", self.recvn(4))[0]
        forces = np.frombuffer(self.recvn(8 * 3 * nat), dtype=np.float64).reshape(nat, 3)
        vir = np.frombuffer(self.recvn(72), dtype=np.float64).reshape(3, 3)
        nextra = struct.unpack("<i", self.recvn(4))[0]
        if nextra:
            self.recvn(nextra)
        # to eV and Angstrom
        return pot * HARTREE, forces * HARTREE / BOHR, vir * HARTREE

    def exit(self):
        self.put_header("EXIT")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("configs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--timeout", type=float, default=600.0)
    ap.add_argument("--emit-effective", metavar="PREFIX",
                    help="write the geometry actually sent, in Angstrom, at 17 digits")
    args = ap.parse_args()

    srv = Server(args.name)
    print(f"reference server listening on {srv.path}", flush=True)
    try:
        srv.accept(args.timeout)
        print("client connected", flush=True)

        srv.expect_status("NEEDINIT", "before INIT")
        srv.init()

        for k, cfg in enumerate(args.configs):
            sp, pos, cell = read_extxyz(cfg)
            if args.emit_effective:
                write_effective(f"{args.emit_effective}_{k}.xyz", sp, pos, cell)
            srv.expect_status("READY", f"before geometry {k}")
            srv.posdata(cell, pos)
            srv.expect_status("HAVEDATA", f"after geometry {k}")
            pot, forces, vir = srv.getforce()
            with open(f"{args.out}_{k}.dat", "w") as fh:
                fh.write(f"# energy_eV {pot:.16e}\n")
                fh.write("# virial_eV " + " ".join(f"{v:.16e}" for v in vir.ravel()) + "\n")
                for row in forces:
                    fh.write(" ".join(f"{v:.16e}" for v in row) + "\n")
            print(f"geometry {k}: E = {pot:.9f} eV, max|f| = {np.abs(forces).max():.9f} eV/A",
                  flush=True)

        srv.expect_status("READY", "before EXIT")
        srv.exit()
        print("sent EXIT", flush=True)
    finally:
        srv.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
