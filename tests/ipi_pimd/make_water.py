#!/usr/bin/env python3
"""Write a small box of water molecules for the PIMD test.

Generated rather than taken from the test-data repository so the test says what
it is running on: four molecules at the experimental gas-phase geometry, far
enough apart not to start on top of each other, in a box big enough that the
ring polymers do not see their own images at the cutoff.
"""

import math
import sys

BOX = 12.0
O_H = 0.9572          # Angstrom
ANGLE = 104.52        # degrees
SPACING = 5.0


def molecule(origin, turn):
    """One water, rotated about z so the four are not all aligned."""
    half = math.radians(ANGLE) / 2.0
    local = [
        ("O", 0.0, 0.0, 0.0),
        ("H", O_H * math.sin(half), O_H * math.cos(half), 0.0),
        ("H", -O_H * math.sin(half), O_H * math.cos(half), 0.0),
    ]
    out = []
    for name, x, y, z in local:
        rx = x * math.cos(turn) - y * math.sin(turn)
        ry = x * math.sin(turn) + y * math.cos(turn)
        out.append((name, origin[0] + rx, origin[1] + ry, origin[2] + z))
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "water.xyz"
    origins = [(2.5, 2.5, 2.5), (2.5 + SPACING, 2.5, 2.5),
               (2.5, 2.5 + SPACING, 2.5), (2.5 + SPACING, 2.5 + SPACING, 2.5)]
    atoms = []
    for index, origin in enumerate(origins):
        atoms.extend(molecule(origin, index * math.pi / 5.0))

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(f'Lattice="{BOX} 0.0 0.0 0.0 {BOX} 0.0 0.0 0.0 {BOX}" '
                     'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for name, x, y, z in atoms:
            handle.write(f"{name:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    print(f"{path}: {len(atoms)} atoms, {len(origins)} water molecules, {BOX} A box")


main()
