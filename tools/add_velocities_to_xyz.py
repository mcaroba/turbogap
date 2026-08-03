#!/usr/bin/env python3
"""Build a TurboGAP input xyz that carries explicit initial velocities.

Reads frame 0 of a TurboGAP trajectory_out.xyz (which contains
velocities:R:3) and writes an atoms file with species, positions and
velocities only. Supplying velocities explicitly means TurboGAP never has to
randomize them, so a CPU and a GPU run start from exactly the same state and
their trajectories can be compared directly.

Usage: add_velocities_to_xyz.py <trajectory_out.xyz> <out.xyz>
"""
import re, sys


def main(src, dst):
    with open(src) as f:
        nat = int(f.readline().split()[0])
        comment = f.readline()
        rows = [f.readline().split() for _ in range(nat)]

    props = re.search(r'Properties=(\S+)', comment)
    if not props:
        sys.exit("no Properties field in frame 0")
    fields, col = [], 0
    for i in range(0, len(props.group(1).split(':')), 3):
        name, _kind, n = props.group(1).split(':')[i:i + 3]
        fields.append((name, col, int(n)))
        col += int(n)
    idx = {name: (start, n) for name, start, n in fields}
    for need in ("species", "pos", "velocities"):
        if need not in idx:
            sys.exit(f"frame 0 has no '{need}' column (found {list(idx)})")

    lat = re.search(r'Lattice="([^"]*)"', comment)
    pbc = re.search(r'pbc="([^"]*)"', comment)

    head = []
    if lat:
        head.append('Lattice="%s"' % lat.group(1))
    head.append("Properties=species:S:1:pos:R:3:velocities:R:3")
    head.append('pbc="%s"' % (pbc.group(1) if pbc else "T T T"))

    sp_i = idx["species"][0]
    p_i = idx["pos"][0]
    v_i = idx["velocities"][0]
    with open(dst, "w") as o:
        o.write("%d\n" % nat)
        o.write(" ".join(head) + "\n")
        for r in rows:
            pos = r[p_i:p_i + 3]
            vel = r[v_i:v_i + 3]
            o.write("%-2s %20.12f %20.12f %20.12f %20.12e %20.12e %20.12e\n"
                    % (r[sp_i], *[float(x) for x in pos], *[float(x) for x in vel]))
    print(f"wrote {dst}: {nat} atoms with velocities")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
