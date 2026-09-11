#!/usr/bin/env python3
"""Build the structure for the topology regression case.

An oxidised graphene sheet with free CO molecules above it: the case
mc_mol_identify = topology exists for, in the elements the C/O potential in the
test data covers. The sheet's own oxygens are bonded into it, so a match on
composition alone would find CO molecules that are not there; the free ones are
separate components of exactly two atoms.

Written out rather than generated at run time so the case is reproducible, and
kept in the case directory rather than the test-data repository because it is
eighty atoms.
"""

import math

A = 2.46          # graphene lattice constant, Angstrom
CO_BOND = 1.128   # the CO bond length, as in co_molecule.xyz
CELLS = 4         # 4x4 unit cells -> 32 carbons
Z_HEIGHT = 20.0   # vacuum, so the sheet does not see its own image
N_CO = 3          # free molecules above the sheet


def sheet():
    """A 4x4 graphene sheet, two atoms per cell."""
    a1 = (A, 0.0)
    a2 = (A / 2.0, A * math.sqrt(3.0) / 2.0)
    basis = [(0.0, 0.0), (A / 2.0, A / (2.0 * math.sqrt(3.0)))]
    atoms = []
    for i in range(CELLS):
        for j in range(CELLS):
            for bx, by in basis:
                x = i * a1[0] + j * a2[0] + bx
                y = i * a1[1] + j * a2[1] + by
                atoms.append(("C", x, y, 0.0))
    return atoms


def epoxides(carbons):
    """Two epoxy oxygens bridging C-C bonds, 1.22 A above the sheet.

    Bonded to the sheet, so they belong to its component and must not be
    matched -- which is the point of putting them here.
    """
    out = []
    for index in (0, 10):
        _, x1, y1, _ = carbons[index]
        _, x2, y2, _ = carbons[index + 1]
        out.append(("O", (x1 + x2) / 2.0, (y1 + y2) / 2.0, 1.22))
    return out


def free_molecules():
    """CO molecules well above the sheet and well apart from each other."""
    out = []
    for k in range(N_CO):
        x = 1.5 + k * 3.0
        y = 1.5 + k * 1.7
        out.append(("C", x, y, 6.0))
        out.append(("O", x, y, 6.0 + CO_BOND))
    return out


def main():
    carbons = sheet()
    atoms = carbons + epoxides(carbons) + free_molecules()
    ax = CELLS * A
    ay = CELLS * A * math.sqrt(3.0) / 2.0

    with open("intercalated.xyz", "w", encoding="utf-8") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            f'Lattice="{ax:.6f} 0.0 0.0 {A * CELLS / 2.0:.6f} {ay:.6f} 0.0 0.0 0.0 {Z_HEIGHT:.6f}" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for name, x, y, z in atoms:
            handle.write(f"{name:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

    print(f"intercalated.xyz: {len(atoms)} atoms "
          f"({len(carbons)} sheet C, {len(atoms) - len(carbons) - 2 * N_CO} epoxy O, "
          f"{N_CO} free CO)")


main()
