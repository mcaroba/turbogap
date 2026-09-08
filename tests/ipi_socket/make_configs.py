#!/usr/bin/env python3
"""The four test configurations, written at 17 significant digits.

TWO CELLS, and the difference between them is the point.

  big    the CO system's own 897 atoms in a cell sheared out of its original
         orthorhombic one. Every perpendicular width stays well above the 4.5 A
         cutoff, so TurboGAP does not replicate: indices is 1 and a_box is the
         primitive cell.

  small  a sub-box of the same system, small enough that TurboGAP builds a
         supercell. indices is then 2 and a_box is the SUPERCELL vector, while
         i-PI still sends the primitive one -- the case where the two
         conventions part company.

BOTH ARE TRICLINIC, and their transposes are different cells. That is what
makes the index order of i-PI's cell matrix testable at all: with a cubic box
a client that reads h transposed gets the right answer and the test proves
nothing.

The sub-box is carved by fractional coordinate, so atoms that were far apart in
the original cell can end up as periodic neighbours across the new faces. Any
pair closer than MIN_SEP is thinned out -- one atom of the pair removed, worst
first. The result is not a physical structure and does not need to be: both
paths evaluate the identical configuration and the test is about the plumbing.
"""
from __future__ import annotations

import sys

import numpy as np

# Rows are the lattice vectors, so the transpose -- which is what goes on the
# wire as i-PI's h -- is upper triangular, as i-PI requires.
BIG = np.array([[17.5589570312, 0.0, 0.0],
                [2.0, 17.7712514289, 0.0],
                [1.5, 2.5, 17.2175376634]])
# Every perpendicular width here is below 2*rcut = 9.0 A, which is what forces
# the replication. run.sh asserts that it actually happened.
SMALL = np.array([[8.5, 0.0, 0.0],
                  [1.2, 8.6, 0.0],
                  [0.9, 1.4, 8.3]])
MIN_SEP = 1.0
NOISE = 0.05          # every force moves, no bond breaks
SEED = 7


def read_xyz(path):
    with open(path) as fh:
        nat = int(fh.readline().split()[0])
        fh.readline()
        sp, pos = [], []
        for _ in range(nat):
            f = fh.readline().split()
            sp.append(f[0])
            pos.append([float(x) for x in f[1:4]])
    return sp, np.array(pos)


def write(path, species, pos, cell):
    with open(path, "w") as fh:
        fh.write(f"{len(species)}\n")
        fh.write('Lattice="' + " ".join(f"{v:.17g}" for v in cell.ravel())
                 + '" Properties=species:S:1:pos:R:3\n')
        for s, r in zip(species, pos):
            fh.write(f"{s} " + " ".join(f"{v:.17g}" for v in r) + "\n")


def min_image_distances(pos, cell):
    """All pair distances under the minimum image of `cell` (rows = vectors)."""
    inv = np.linalg.inv(cell)
    d = pos[:, None, :] - pos[None, :, :]
    s = d @ inv
    s -= np.round(s)
    d = s @ cell
    return np.sqrt((d ** 2).sum(-1))


def carve(species, pos, cell):
    """Atoms inside `cell`, thinned until nothing is closer than MIN_SEP."""
    s = pos @ np.linalg.inv(cell)
    keep = np.all((s >= 0.0) & (s < 1.0), axis=1)
    idx = np.flatnonzero(keep)
    sp = [species[i] for i in idx]
    p = pos[idx]
    while True:
        d = min_image_distances(p, cell)
        np.fill_diagonal(d, np.inf)
        if d.min() >= MIN_SEP:
            return sp, p
        i, j = np.unravel_index(d.argmin(), d.shape)
        # Drop whichever of the offending pair is more crowded.
        drop = i if (d[i] < MIN_SEP).sum() >= (d[j] < MIN_SEP).sum() else j
        sp = [x for k, x in enumerate(sp) if k != drop]
        p = np.delete(p, drop, axis=0)


def main() -> int:
    species, pos = read_xyz(sys.argv[1])
    rng = np.random.default_rng(SEED)

    write("big_0.xyz", species, pos, BIG)
    write("big_1.xyz", species, pos + rng.normal(0.0, NOISE, pos.shape), BIG)
    print(f"  big:   {len(species)} atoms, triclinic, no supercell expected")

    sp, p = carve(species, pos, SMALL)
    if len(sp) < 20:
        print(f"  ERROR: the sub-box kept only {len(sp)} atoms; widen SMALL")
        return 1
    write("small_0.xyz", sp, p, SMALL)
    write("small_1.xyz", sp, p + rng.normal(0.0, NOISE, p.shape), SMALL)
    d = min_image_distances(p, SMALL)
    np.fill_diagonal(d, np.inf)
    print(f"  small: {len(sp)} atoms, triclinic, closest approach "
          f"{d.min():.3f} A, supercell expected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
