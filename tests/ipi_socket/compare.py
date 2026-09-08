#!/usr/bin/env python3
"""Socket path against file path, for both phases, plus the cell control.

Run from the working directory run.sh built. Exits non-zero on any failure and
says which check failed and by how much.
"""
from __future__ import annotations

import re
import sys

import numpy as np

# Two effects set the floor, and both are accounted for rather than absorbed.
#
#   The PRINT FORMAT. trajectory_out.xyz carries forces and the virial to 8
#   decimals and the energy to 6, so values that agree exactly can still differ
#   by half a unit in the last place printed.
#
#   THE WRAP. `turbogap ipi` wraps what i-PI sends, because every other mode in
#   the tree keeps its positions wrapped and a long run would otherwise let the
#   coordinates grow without bound. `turbogap predict` does not wrap at all. So
#   the two paths evaluate coordinate sets that differ by exact lattice
#   translations: the same configuration and the same forces, but not the same
#   last bits, because minimum-image arithmetic takes a different branch on
#   each. Feeding predict the pre-wrapped geometry does not close it, since
#   TurboGAP's wrap follows get_distance's shift convention and that is not
#   reproducible from outside the code without copying it.
#
# 5e-7 eV/A is under 1e-9 of the largest force in this system, and transposing
# the cell moves those forces by hundreds of eV/A, so the margin between
# "agrees" and "a convention is wrong" is nine orders of magnitude wide.
TOL_F = 5.0e-7
TOL_V = 5.0e-7
TOL_E = 1.0e-5
# What "different" means for the control.
DISCRIMINATE = 1.0


def read_traj(path):
    out = []
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                return out
            nat = int(line.split()[0])
            c = fh.readline()
            E = float(re.search(r"energy=(\S+)", c).group(1))
            v = np.array([float(x) for x in
                          re.search(r'virial="([^"]*)"', c).group(1).split()]).reshape(3, 3)
            f = np.array([[float(x) for x in fh.readline().split()[4:7]] for _ in range(nat)])
            out.append((E, f, v))


def read_socket(path):
    lines = open(path).read().splitlines()
    E = float(lines[0].split()[2])
    v = np.array([float(x) for x in lines[1].split()[2:]]).reshape(3, 3)
    f = np.array([[float(x) for x in ln.split()] for ln in lines[2:]])
    return E, f, v


def main() -> int:
    fail = False

    ref_big = read_traj("ref_big.xyz")
    refT = read_traj("refT.xyz")
    d = np.abs(ref_big[0][1] - refT[0][1]).max()
    print(f"  control: transposing the cell moves the forces by {d:.3f} eV/A "
          f"(needs > {DISCRIMINATE})")
    if d <= DISCRIMINATE:
        print("  FAIL: the triclinic cell does not discriminate; the agreement "
              "checks would pass with the index order reversed")
        fail = True

    for tag, label in (("big", "no supercell, indices = 1"),
                       ("small", "supercell, indices > 1")):
        ref = read_traj(f"ref_{tag}.xyz")
        print(f"  --- {tag} ({label}) ---")
        for k in range(len(ref)):
            E, f, v = ref[k]
            Es, fs, vs = read_socket(f"sock_{tag}_{k}.dat")
            dE, df, dv = abs(Es - E), np.abs(fs - f).max(), np.abs(vs - v).max()
            print(f"  geometry {k}: dE = {dE:.2e} eV (tol {TOL_E:.0e}), "
                  f"max|df| = {df:.2e} eV/A (tol {TOL_F:.0e}), "
                  f"max|dvirial| = {dv:.2e} eV (tol {TOL_V:.0e})")
            print(f"              scale |E| = {abs(E):.2f}, "
                  f"max|f| = {np.abs(f).max():.2f}, max|virial| = {np.abs(v).max():.1f}")
            if dE > TOL_E or df > TOL_F or dv > TOL_V:
                print(f"  FAIL: {tag} geometry {k} disagrees beyond the print granularity")
                fail = True

    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
