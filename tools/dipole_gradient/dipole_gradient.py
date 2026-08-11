#!/usr/bin/env python3
"""Cartesian gradient of the dipole, d(mu_a)/d(r_jb), by central differences.

WHY FINITE DIFFERENCE. A GAP dipole model predicts mu_i = dE_i/dr_i, so its
Cartesian gradient is a second derivative of the fitted scalar:

    d(mu_ia)/d(r_jb) = sum_de (d2E_i/dq_d dq_e)(dq_d/dr_ia)(dq_e/dr_jb)   [A]
                     + sum_d  (dE_i/dq_d)(d2q_d/dr_ia dr_jb)             [B]

Term A is computable from quantities TurboGAP already forms. Term B needs the
SECOND derivative of the SOAP descriptor, which soap_turbo does not provide --
it returns soap_cart_der and nothing beyond it. So an exact analytic gradient
is not available without re-differentiating soap_turbo_radial/angular/compress.
Central differences need none of that, are exact to O(h^2), and are the
reference any future analytic implementation has to reproduce.

COST. 6N+1 configurations, but TurboGAP evaluates every frame of an xyz file in
one invocation, so this is a single run, not 6N runs. That is affordable for
validation and for small systems; it is not a production route for MD on large
cells. See README.md for the charge-based reformulation that makes the gradient
first-order.

CHECK. mu is invariant under a rigid translation of the whole system, so

    sum_j d(mu_a)/d(r_jb) = 0   for every a, b

That sum rule is computed and reported: it is a property of the model, not of
the differencing, so a large residual means something is wrong beyond step size.

Usage:
    dipole_gradient.py --gap gap_files/water_dipole.gap --atoms frame.xyz \\
                       --turbogap /path/to/bin/turbogap [-h-step 1e-4] [--out J.dat]

`--atoms` must hold ONE frame. Species and masses are taken from `--species`
(default "H O" with water masses) or from an existing `input` via --template.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile


def read_frame(path):
    with open(path) as f:
        n = int(f.readline())
        comment = f.readline().rstrip('\n')
        species, pos = [], []
        for _ in range(n):
            t = f.readline().split()
            species.append(t[0])
            pos.append([float(x) for x in t[1:4]])
    return species, pos, comment


def lattice_of(comment):
    m = re.search(r'Lattice="([^"]*)"', comment)
    if not m:
        raise SystemExit('ERROR: the frame has no Lattice="..." in its comment line')
    return m.group(1)


def write_frames(path, species, frames, lattice):
    with open(path, 'w') as out:
        for pos in frames:
            out.write('%d\n' % len(species))
            out.write('Lattice="%s" Properties=species:S:1:pos:R:3 pbc="T T T"\n' % lattice)
            for s, p in zip(species, pos):
                out.write('%s %.12f %.12f %.12f\n' % (s, p[0], p[1], p[2]))


def read_dipoles(path):
    """Return one dipole 3-vector per frame of a trajectory_out.xyz."""
    out = []
    with open(path) as f:
        while True:
            count = f.readline()
            if not count.strip():
                break
            n = int(count)
            header = f.readline()
            m = re.search(r'\bdipole="([^"]*)"', header)
            if not m:
                raise SystemExit(
                    'ERROR: no dipole= in %s. Does the .gap set dipole_model = .true.?' % path)
            out.append([float(x) for x in m.group(1).split()])
            for _ in range(n):
                f.readline()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gap', required=True, help='path to the dipole .gap file')
    ap.add_argument('--atoms', required=True, help='xyz holding exactly one frame')
    ap.add_argument('--turbogap', required=True, help='turbogap binary')
    # 5e-4 is the measured optimum for this output format: below it the F16.8
    # rounding of dipole= dominates (error ~ 1/h), above it the O(h^2)
    # truncation does. See README.md for the scan.
    ap.add_argument('--step', type=float, default=5e-4, help='displacement h in Angstrom')
    ap.add_argument('--species', default='H O')
    ap.add_argument('--masses', default='1.008 15.999')
    ap.add_argument('--out', default='dipole_gradient.dat')
    ap.add_argument('--keep', action='store_true', help='keep the staging directory')
    args = ap.parse_args()

    species, pos, comment = read_frame(args.atoms)
    lattice = lattice_of(comment)
    n = len(species)
    h = args.step

    work = tempfile.mkdtemp(prefix='dipole_grad.')
    try:
        gap_dir = os.path.dirname(os.path.abspath(args.gap)) or '.'
        shutil.copytree(gap_dir, os.path.join(work, 'gap_files'))
        gap_name = os.path.basename(args.gap)

        # Frame 0 is the undisplaced reference; then +h/-h for every (atom, cart).
        frames, order = [list(map(list, pos))], []
        for j in range(n):
            for b in range(3):
                for sgn in (+1, -1):
                    p = [list(q) for q in pos]
                    p[j][b] += sgn * h
                    frames.append(p)
                    order.append((j, b, sgn))
        write_frames(os.path.join(work, 'atoms.xyz'), species, frames, lattice)

        with open(os.path.join(work, 'input'), 'w') as f:
            f.write('atoms_file = "atoms.xyz"\n')
            f.write('pot_file = "gap_files/%s"\n' % gap_name)
            f.write('n_species = %d\n' % len(args.species.split()))
            f.write('species = %s\n' % args.species)
            f.write('masses = %s\n' % args.masses)
            f.write('e0 = %s\n' % ' '.join('0.' for _ in args.species.split()))
            f.write('random_seed = 12345\n\nwrite_xyz = 1\n')

        # One invocation for all 6N+1 frames.
        r = subprocess.run([os.path.abspath(args.turbogap), 'predict'],
                           cwd=work, capture_output=True, text=True)
        if r.returncode != 0:
            sys.stderr.write(r.stdout[-3000:])
            raise SystemExit('ERROR: turbogap exited %d' % r.returncode)

        mu = read_dipoles(os.path.join(work, 'trajectory_out.xyz'))
        if len(mu) != len(frames):
            raise SystemExit('ERROR: expected %d frames back, got %d' % (len(frames), len(mu)))

        # J[a][j][b] = d(mu_a)/d(r_jb)
        J = [[[0.0] * 3 for _ in range(n)] for _ in range(3)]
        plus, minus = {}, {}
        for idx, (j, b, sgn) in enumerate(order, start=1):
            (plus if sgn > 0 else minus)[(j, b)] = mu[idx]
        for j in range(n):
            for b in range(3):
                for a in range(3):
                    J[a][j][b] = (plus[(j, b)][a] - minus[(j, b)][a]) / (2 * h)

        with open(args.out, 'w') as f:
            f.write('# d(mu_a)/d(r_jb) by central differences, h = %g\n' % h)
            f.write('# dipole = %s\n' % ' '.join('%.10f' % x for x in mu[0]))
            f.write('# atom  species   a   b        d(mu_a)/d(r_jb)\n')
            for j in range(n):
                for a in range(3):
                    for b in range(3):
                        f.write('%6d %8s %3s %3s %22.10f\n'
                                % (j, species[j], 'xyz'[a], 'xyz'[b], J[a][j][b]))

        print('dipole            %s' % '  '.join('%.8f' % x for x in mu[0]))
        print('wrote             %s  (%d atoms, h = %g, %d frames in one run)'
              % (args.out, n, h, len(frames)))

        # Translational sum rule: a rigid shift cannot change mu.
        worst = 0.0
        scale = max(abs(J[a][j][b]) for a in range(3) for j in range(n) for b in range(3))
        print('\ntranslational sum rule  sum_j d(mu_a)/d(r_jb) = 0')
        for a in range(3):
            row = []
            for b in range(3):
                s = sum(J[a][j][b] for j in range(n))
                worst = max(worst, abs(s))
                row.append('%12.3e' % s)
            print('  a=%s   %s' % ('xyz'[a], ' '.join(row)))
        print('  max |residual| %.3e   (largest |J| element %.3e, ratio %.1e)'
              % (worst, scale, worst / scale if scale else 0.0))
    finally:
        if args.keep:
            print('\nstaging kept in %s' % work)
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
