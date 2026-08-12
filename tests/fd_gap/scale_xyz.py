#!/usr/bin/env python3
"""Scale an extxyz frame's cell and positions by a common factor.

The core potential is the short-range repulsion a GAP is not fitted to
describe, so on an equilibrated structure it is identically zero -- every pair
sits outside its table. A finite-difference check of a contribution that is
zero passes without checking anything, so the core_pot leg compresses the cell
until the contribution switches on.

Fractional coordinates are unchanged, so this is a uniform strain and nothing
else.

Usage:  scale_xyz.py <in.xyz> <factor> <out.xyz>
"""
import re
import sys


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    src, factor, dst = sys.argv[1], float(sys.argv[2]), sys.argv[3]

    lines = open(src).read().splitlines()
    nat = int(lines[0].split()[0])
    comment = lines[1]

    m = re.search(r'Lattice="([^"]+)"', comment)
    if m is None:
        raise SystemExit("no Lattice= in the comment line")
    lattice = [float(v) * factor for v in m.group(1).split()]
    comment = (comment[:m.start(1)]
               + " ".join("%.10f" % v for v in lattice)
               + comment[m.end(1):])

    out = [lines[0], comment]
    for line in lines[2:2 + nat]:
        w = line.split()
        w[1:4] = ["%.10f" % (float(x) * factor) for x in w[1:4]]
        out.append(" ".join(w))

    open(dst, "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
