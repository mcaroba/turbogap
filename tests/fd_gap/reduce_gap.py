#!/usr/bin/env python3
"""Keep only one kind of block from a .gap file.

A potential file is a flat sequence of `gap_beg <type> ... gap_end` blocks, and
dropping all but one type leaves a potential whose total energy, forces and
virial belong to that type alone. That is what lets a finite-difference check
of, say, the 3b forces be a plain check of the total forces, with no
per-contribution diagnostic to add to the Fortran and nothing to difference.

The `desc_sparse` / `alphas_sparse` paths inside the blocks are copied
unchanged, so the reduced file has to be used from a directory where they still
resolve -- in practice, next to the same `gap_files` symlink as the original.

Usage:  reduce_gap.py <in.gap> <type> <out.gap>
"""
import sys


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    src, want, dst = sys.argv[1:]

    out, keep = [], False
    for line in open(src):
        stripped = line.strip()
        if stripped.startswith("gap_beg"):
            parts = stripped.split()
            keep = len(parts) > 1 and parts[1] == want
            if keep:
                out.append(line)
            continue
        if stripped.startswith("gap_end"):
            if keep:
                out.append(line)
            keep = False
            continue
        if keep:
            out.append(line)

    if not out:
        raise SystemExit(f"reduce_gap.py: {src} has no '{want}' block")

    open(dst, "w").writelines(out)
    print(f"    {dst}: {sum(1 for l in out if l.strip().startswith('gap_beg'))} "
          f"{want} block(s)")


if __name__ == "__main__":
    main()
