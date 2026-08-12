#!/usr/bin/env python3
"""Compare two extxyz trajectories frame by frame, atom by atom.

The left file is what a Monte Carlo run recorded as it went; the right is a
single-shot prediction on those same coordinates. If the neighbour lists the MC
path carried were right, the two agree to round-off. If a rank was left holding
a stale list across an insertion or a removal, the atoms whose neighbourhood
changed disagree by far more than that, and this says which ones.

Frames are matched in order and must have the same species and coordinates;
what is compared is the per-atom forces and local energies.
"""
import sys
import argparse


def frames(path):
    """Yield (natoms, comment, [ (symbol, [fields...]) ]) per frame."""
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if not line.strip():
                continue
            nat = int(line.split()[0])
            comment = fh.readline()
            rows = [fh.readline().split() for _ in range(nat)]
            yield nat, comment, rows


def columns(comment):
    """Map an extxyz Properties= spec to {name: (first_column, width)}."""
    spec = None
    for tok in comment.split():
        if tok.startswith("Properties="):
            spec = tok[len("Properties="):]
            break
    if spec is None:
        raise SystemExit("no Properties= in comment line")
    out, col = {}, 0
    parts = spec.split(":")
    for i in range(0, len(parts), 3):
        name, _kind, width = parts[i], parts[i + 1], int(parts[i + 2])
        out[name] = (col, width)
        col += width
    return out


def get(rows, cols, name):
    first, width = cols[name]
    return [[float(r[first + k]) for k in range(width)] for r in rows]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mc")
    p.add_argument("predict")
    p.add_argument("--force-tol", type=float, default=1e-6,
                   help="max absolute force component difference, eV/A")
    p.add_argument("--energy-tol", type=float, default=1e-7,
                   help="max absolute local energy difference, eV")
    p.add_argument("--quiet", action="store_true")
    a = p.parse_args()

    worst_f = worst_e = 0.0
    worst_f_at = worst_e_at = None
    n = 0

    for i, (left, right) in enumerate(zip(frames(a.mc), frames(a.predict))):
        nat_l, com_l, rows_l = left
        nat_r, com_r, rows_r = right
        if nat_l != nat_r:
            raise SystemExit(
                "frame %d: %d atoms in %s but %d in %s -- the prediction did "
                "not follow the same frames" % (i, nat_l, a.mc, nat_r, a.predict))
        cl, cr = columns(com_l), columns(com_r)

        if [r[0] for r in rows_l] != [r[0] for r in rows_r]:
            raise SystemExit("frame %d: species differ" % i)

        # Positions are the input to both sides; a difference here would mean
        # the prediction was not run on these coordinates at all.
        pl, pr = get(rows_l, cl, "pos"), get(rows_r, cr, "pos")
        dp = max(abs(x - y) for u, v in zip(pl, pr) for x, y in zip(u, v))
        if dp > 1e-6:
            raise SystemExit("frame %d: positions differ by %.3g A" % (i, dp))

        fl, fr = get(rows_l, cl, "forces"), get(rows_r, cr, "forces")
        for j, (u, v) in enumerate(zip(fl, fr)):
            for x, y in zip(u, v):
                if abs(x - y) > worst_f:
                    worst_f, worst_f_at = abs(x - y), (i, j)

        if "local_energy" in cl and "local_energy" in cr:
            el, er = get(rows_l, cl, "local_energy"), get(rows_r, cr, "local_energy")
            for j, (u, v) in enumerate(zip(el, er)):
                if abs(u[0] - v[0]) > worst_e:
                    worst_e, worst_e_at = abs(u[0] - v[0]), (i, j)
        n += 1

    if n == 0:
        raise SystemExit("no frames compared")

    ok = worst_f <= a.force_tol and worst_e <= a.energy_tol
    if not a.quiet or not ok:
        print("    %d frames" % n)
        print("    max |dF| = %.3e eV/A (tol %.0e)%s"
              % (worst_f, a.force_tol,
                 "" if worst_f_at is None else "  at frame %d atom %d" % worst_f_at))
        print("    max |dE_i| = %.3e eV   (tol %.0e)%s"
              % (worst_e, a.energy_tol,
                 "" if worst_e_at is None else "  at frame %d atom %d" % worst_e_at))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
