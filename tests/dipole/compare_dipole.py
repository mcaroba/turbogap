#!/usr/bin/env python3
"""Compare TurboGAP's dipole prediction against a QUIP reference.

Unlike tests/regression, this is not a bit-identity check: it asks whether two
independent codes agree on the same model, so it needs a tolerance. Both call
the same soap_turbo library, so the residual should sit at print precision
(TurboGAP writes F16.8, i.e. ~5e-9 of rounding on these magnitudes), not at the
level of a real numerical difference.

Checks three things:
  1. per-atom local_dipole against QUIP
  2. per-frame dipole against QUIP
  3. that TurboGAP's own dipole equals the sum of its local_dipole columns,
     which QUIP's output also satisfies

Usage: compare_dipole.py reference.dat trajectory_out.xyz [tolerance]
"""

import re
import sys


def read_reference(path):
    frames = []
    with open(path) as f:
        rows = [ln for ln in f if not ln.lstrip().startswith('#') and ln.strip()]
    i = 0
    while i < len(rows):
        _, n_atoms = (int(x) for x in rows[i].split())
        dipole = [float(x) for x in rows[i + 1].split()]
        local = [[float(x) for x in rows[i + 2 + k].split()] for k in range(n_atoms)]
        frames.append((dipole, local))
        i += 2 + n_atoms
    return frames


def properties_column(header, name):
    m = re.search(r'Properties=(\S+)', header)
    if not m:
        return None
    fields = m.group(1).split(':')
    col = 0
    for i in range(0, len(fields) - 2, 3):
        if fields[i] == name:
            return col
        col += int(fields[i + 2])
    return None


def read_turbogap(path):
    frames = []
    with open(path) as f:
        while True:
            count = f.readline()
            if not count.strip():
                break
            n_atoms = int(count)
            header = f.readline()
            m = re.search(r'dipole="([^"]*)"', header)
            if not m:
                raise SystemExit(
                    'ERROR: %s has no dipole= on the comment line. Is dipole_model '
                    '= .true. set in the .gap file?' % path)
            dipole = [float(x) for x in m.group(1).split()]
            col = properties_column(header, 'local_dipole')
            if col is None:
                raise SystemExit('ERROR: %s has no local_dipole in Properties' % path)
            local = []
            for _ in range(n_atoms):
                local.append([float(x) for x in f.readline().split()[col:col + 3]])
            frames.append((dipole, local))
    return frames


def stats(pairs):
    """Return (max_abs_error, rms_error) over a list of (a, b) scalars."""
    if not pairs:
        return 0.0, 0.0
    diffs = [abs(a - b) for a, b in pairs]
    rms = (sum(d * d for d in diffs) / len(diffs)) ** 0.5
    return max(diffs), rms


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    ref_path, tg_path = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) == 4 else 1e-6

    ref = read_reference(ref_path)
    tg = read_turbogap(tg_path)

    if len(ref) != len(tg):
        sys.exit('ERROR: frame count differs: reference %d, turbogap %d'
                 % (len(ref), len(tg)))

    local_pairs, total_pairs, sum_pairs = [], [], []
    worst_frame = (0.0, -1)
    for n, ((r_dip, r_loc), (t_dip, t_loc)) in enumerate(zip(ref, tg)):
        if len(r_loc) != len(t_loc):
            sys.exit('ERROR: frame %d atom count differs: reference %d, turbogap %d'
                     % (n, len(r_loc), len(t_loc)))
        frame_worst = 0.0
        for rm, tm in zip(r_loc, t_loc):
            for a, b in zip(rm, tm):
                local_pairs.append((a, b))
                frame_worst = max(frame_worst, abs(a - b))
        for a, b in zip(r_dip, t_dip):
            total_pairs.append((a, b))
            frame_worst = max(frame_worst, abs(a - b))
        # TurboGAP internal consistency: mu = sum_i mu_i
        for k in range(3):
            sum_pairs.append((t_dip[k], sum(m[k] for m in t_loc)))
        if frame_worst > worst_frame[0]:
            worst_frame = (frame_worst, n)

    local_max, local_rms = stats(local_pairs)
    total_max, total_rms = stats(total_pairs)
    sum_max, _ = stats(sum_pairs)

    print('frames compared            %d  (%d atoms/frame)' % (len(ref), len(ref[0][1])))
    print('local_dipole vs QUIP       max %.3e   rms %.3e' % (local_max, local_rms))
    print('dipole       vs QUIP       max %.3e   rms %.3e' % (total_max, total_rms))
    print('dipole vs sum(local)       max %.3e   (turbogap internal)' % sum_max)
    print('worst frame                %d' % worst_frame[1])
    print('tolerance                  %.3e' % tol)

    bad = []
    if local_max > tol:
        bad.append('local_dipole')
    if total_max > tol:
        bad.append('dipole')
    # The internal sum is over n_atoms terms of printed values, so it accumulates
    # a little more rounding than a single comparison does.
    if sum_max > tol * max(1, len(ref[0][1])):
        bad.append('dipole/sum(local) consistency')

    if bad:
        print('FAIL: %s outside tolerance' % ', '.join(bad))
        return 1
    print('PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
