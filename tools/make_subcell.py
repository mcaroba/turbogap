#!/usr/bin/env python3
"""Cut a smaller periodic cell out of a large orthorhombic extended-XYZ file.

    tools/make_subcell.py atoms_1000000.xyz --cells 25 --out atoms_125000.xyz
    tools/make_subcell.py atoms_1000000.xyz --fraction 0.5 --out half.xyz
    tools/make_subcell.py atoms_1000000.xyz --info

Why this exists rather than "take the first N atoms": the first N lines of a
supercell are not a cell. They are a slab with two open faces, and running it
periodic puts atoms at bonding distance across a boundary that is not there,
which changes the energy, the neighbour counts, and therefore every performance
number you were trying to measure. A scaling ladder built that way measures the
wrong systems.

This takes a sub-BOX instead: keep every atom with 0 <= r_i < L_i', and set the
lattice to L'. For that to be a valid periodic cell, L' has to be a whole number
of the underlying crystal repeat, which is what --cells expresses. The
1,000,000-atom diamond in this project is 50x50x50 conventional cells of 3.5 A,
so --cells 25 gives 125,000 atoms in a real 87.5 A cell.

Only orthorhombic lattices are handled, and the script refuses anything else
rather than silently producing a cell whose volume is wrong.
"""

import argparse
import sys


def parse_lattice(comment):
    """Pull Lattice="..." out of an extended-XYZ comment line."""
    key = 'Lattice="'
    i = comment.find(key)
    if i < 0:
        raise SystemExit('no Lattice="..." in the comment line; not extended XYZ')
    j = comment.index('"', i + len(key))
    v = [float(x) for x in comment[i + len(key):j].split()]
    if len(v) != 9:
        raise SystemExit(f'Lattice has {len(v)} numbers, expected 9')
    return v


def check_orthorhombic(v):
    off = [v[1], v[2], v[3], v[5], v[6], v[7]]
    if any(abs(x) > 1e-8 for x in off):
        raise SystemExit('lattice is not orthorhombic; this script would produce '
                         'a cell with the wrong volume. Refusing.')
    return v[0], v[4], v[8]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('xyz')
    p.add_argument('--out', help='output file (required unless --info)')
    p.add_argument('--cells', type=int,
                   help='keep this many crystal repeats per side (needs --repeat)')
    p.add_argument('--repeat', type=float,
                   help='length of one crystal repeat in Angstrom. Default: infer '
                        'by assuming the input is a cubic n^3 supercell of a '
                        '4-atom-per-cell fcc lattice (diamond, zincblende).')
    p.add_argument('--fraction', type=float,
                   help='keep this fraction of each side instead (0 < f <= 1). '
                        'Does not check that the result is a whole crystal repeat.')
    p.add_argument('--shift', type=float,
                   help='move the cut planes by this many Angstrom before '
                        'cutting. Default: one eighth of a repeat, which puts '
                        'the boundary between lattice planes instead of '
                        'through one.')
    p.add_argument('--info', action='store_true',
                   help='report the cell and the sizes --cells can produce')
    a = p.parse_args()

    with open(a.xyz) as f:
        n_atoms = int(f.readline().split()[0])
        comment = f.readline()
        lat = parse_lattice(comment)
        lx, ly, lz = check_orthorhombic(lat)

        # A cubic supercell of a conventional cubic cell with 8 atoms (diamond)
        # has n^3 * 8 atoms. Solve for n, then repeat = L/n.
        n_cells = round((n_atoms / 8.0) ** (1.0 / 3.0))
        inferred = lx / n_cells if n_cells else 0.0

        if a.info:
            print(f'{a.xyz}')
            print(f'  atoms   {n_atoms}')
            print(f'  cell    {lx:.4f} x {ly:.4f} x {lz:.4f} A')
            print(f'  density {n_atoms / (lx * ly * lz):.5f} atoms/A^3')
            print(f'  looks like {n_cells}^3 conventional cells of '
                  f'{inferred:.4f} A ({n_cells ** 3 * 8} atoms)')
            print('  sizes --cells can give:')
            for k in (8, 12, 16, 20, 25, 32, 40, n_cells):
                if k <= n_cells:
                    print(f'    --cells {k:3d}  ->{k ** 3 * 8:9d} atoms, '
                          f'{k * inferred:8.3f} A box')
            return 0

        if not a.out:
            p.error('--out is required unless --info')

        repeat = a.repeat if a.repeat is not None else inferred
        if a.fraction is not None:
            fx = fy = fz = a.fraction
            if not 0 < a.fraction <= 1:
                p.error('--fraction must be in (0, 1]')
        elif a.cells is not None:
            side = a.cells * repeat
            fx, fy, fz = side / lx, side / ly, side / lz
            if side > lx + 1e-8:
                p.error(f'--cells {a.cells} needs {side:.3f} A but the cell is '
                        f'{lx:.3f} A')
        else:
            p.error('give --cells or --fraction (or --info)')

        # Half-open interval on each axis: an atom exactly on the far face
        # belongs to the next image and keeping it would double-count a plane.
        bx, by, bz = lx * fx, ly * fy, lz * fz

        # Two corrections, both of which showed up as a density that came out
        # BELOW the parent cell's -- i.e. atoms lost, and lost one-sidedly.
        #
        # 1. Wrap into [0, L) first. This is a thermal snapshot, so an atom that
        #    started at x = 0 can be at x = -0.05, and the parent file keeps it
        #    there. A raw `0 <= x` test drops it, and its periodic image at
        #    x ~ L is outside the sub-box too, so it is lost from both ends.
        #
        # 2. Offset the cut plane off the lattice plane. Cutting at a whole
        #    number of repeats puts the boundary exactly through a plane of
        #    atoms, where thermal displacement decides each one's fate. Shifting
        #    by an eighth of a repeat moves the boundary into the gap between
        #    planes, where there is nothing to be ambiguous about. It does not
        #    change the box, only which periodic image of each atom is taken.
        shift = a.shift if a.shift is not None else repeat / 8.0

        def wrap(v, L):
            v = (v + shift) % L
            return v

        kept = []
        for _ in range(n_atoms):
            line = f.readline()
            if not line:
                break
            parts = line.split()
            x = wrap(float(parts[1]), lx)
            y = wrap(float(parts[2]), ly)
            z = wrap(float(parts[3]), lz)
            if x < bx and y < by and z < bz:
                # Rewrite the position: the wrap and shift moved it, and writing
                # the original would put atoms outside the cell we declare.
                parts[1], parts[2], parts[3] = f'{x:.8f}', f'{y:.8f}', f'{z:.8f}'
                kept.append(' '.join(parts))

    props = 'Properties=species:S:1:pos:R:3'
    if 'Properties=' in comment:
        i = comment.index('Properties=')
        props = comment[i:].split()[0]

    with open(a.out, 'w') as g:
        g.write(f'{len(kept)}\n')
        g.write(f'Lattice="{bx:.10f} 0.0 0.0 0.0 {by:.10f} 0.0 0.0 0.0 {bz:.10f}" '
                f'{props} pbc="T T T"\n')
        g.write('\n'.join(kept))
        g.write('\n')

    print(f'{a.out}: {len(kept)} atoms in {bx:.3f} x {by:.3f} x {bz:.3f} A '
          f'(density {len(kept) / (bx * by * bz):.5f} atoms/A^3, '
          f'was {n_atoms / (lx * ly * lz):.5f})')
    if abs(len(kept) / (bx * by * bz) - n_atoms / (lx * ly * lz)) > 1e-3:
        print('  WARNING: the density moved. The cut did not land on a whole '
              'crystal repeat, so the sub-cell has a seam at its boundary.',
              file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
