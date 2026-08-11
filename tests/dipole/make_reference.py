#!/usr/bin/env python3
"""Extract QUIP's dipole predictions from a `quip` run log into a reference file.

QUIP dumps each frame as an extended-xyz block whose every line is prefixed with
"AT ". The header line of a block carries

    mu="..."       the DFT target the model was fitted against
    dipole="..."   the model's prediction          <- what we compare to
    energy=...     the model's fictitious scalar

and the per-atom lines end with the three local_dipole columns, whose position
is given by the Properties string.

Only `dipole` and `local_dipole` are extracted; `mu` is the training target, not
a prediction, and comparing against it would measure the quality of the fit
rather than whether TurboGAP agrees with QUIP.

Usage: make_reference.py out_quip_predict reference.dat
"""

import re
import sys


def parse_properties(header):
    """Return the column index at which local_dipole starts, or None.

    Properties is a colon-separated list of (name, type, count) triples; the
    column offset is the sum of the counts that precede local_dipole.
    """
    m = re.search(r'Properties=(\S+)', header)
    if not m:
        return None
    fields = m.group(1).split(':')
    col = 0
    for i in range(0, len(fields) - 2, 3):
        name, _, count = fields[i], fields[i + 1], int(fields[i + 2])
        if name == 'local_dipole':
            return col
        col += count
    return None


def parse_vector(header, key):
    m = re.search(key + r'="([^"]*)"', header)
    if not m:
        raise ValueError(f'no {key}= in header: {header[:120]}')
    parts = m.group(1).split()
    if len(parts) != 3:
        raise ValueError(f'{key}= is not a 3-vector: {m.group(1)!r}')
    return [float(x) for x in parts]


def read_quip_log(path):
    """Yield (dipole, [local_dipole per atom]) per frame.

    Each block is a complete extended-xyz frame with every line prefixed "AT ",
    the atom-count line included:

        AT 6
        AT Lattice="..." Properties=...:local_dipole:R:3 dipole="..." ...
        AT O  0.0 0.0 0.0  8  0 0 0  5  0.61 1.16 -1.26
        ... (6 atom lines)
    """
    with open(path) as f:
        lines = [ln[3:] for ln in f if ln.startswith('AT ')]

    i = 0
    while i < len(lines):
        tok = lines[i].split()
        if len(tok) != 1 or not tok[0].isdigit():
            i += 1
            continue
        n_atoms = int(tok[0])

        header = lines[i + 1]
        if 'dipole=' not in header:
            i += 1
            continue
        col = parse_properties(header)
        if col is None:
            raise ValueError('header has dipole= but no local_dipole in Properties')
        dipole = parse_vector(header, 'dipole')

        local = []
        for row in lines[i + 2:i + 2 + n_atoms]:
            cols = row.split()[col:col + 3]
            if len(cols) != 3:
                raise ValueError('atom line has no local_dipole at column %d: %r'
                                 % (col, row[:120]))
            local.append([float(x) for x in cols])
        if len(local) != n_atoms:
            raise ValueError('frame truncated: expected %d atoms, got %d'
                             % (n_atoms, len(local)))

        yield dipole, local
        i += 2 + n_atoms


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]

    frames = list(read_quip_log(src))
    if not frames:
        sys.exit(f'ERROR: no frames with dipole= found in {src}. '
                 'Was quip run with calc_args="dipole=dipole local_dipole=local_dipole"?')

    with open(dst, 'w') as out:
        out.write('# QUIP dipole reference, extracted from %s\n' % src)
        out.write('# frame  n_atoms\n#   dipole_x dipole_y dipole_z\n')
        out.write('#   local_dipole_x local_dipole_y local_dipole_z  (one line per atom)\n')
        for n, (dipole, local) in enumerate(frames):
            out.write('%d %d\n' % (n, len(local)))
            out.write('  %.10f %.10f %.10f\n' % tuple(dipole))
            for mu in local:
                out.write('  %.10f %.10f %.10f\n' % tuple(mu))

    n_atoms = {len(local) for _, local in frames}
    print('wrote %s: %d frames, %s atoms/frame'
          % (dst, len(frames), '/'.join(str(x) for x in sorted(n_atoms))))


if __name__ == '__main__':
    main()
