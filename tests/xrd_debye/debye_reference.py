#!/usr/bin/env python3
"""An independent Debye scattering sum, for checking TurboGAP's xrd_debye path.

Reads an extended-xyz frame and an xrd_prediction.dat / nd_prediction.dat and
compares them. Deliberately naive: an explicit double loop over pairs, no
reuse of anything from the Fortran beyond the tabulated constants, so that a
sign, a factor, a missing self term or a wrong window derivative shows up
rather than cancelling on both sides.

Usage:
    debye_reference.py <atoms.xyz> <prediction.dat> [options]

Exit status 0 if they agree to the tolerance, 1 otherwise.
"""
import argparse
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Scattering power
# ---------------------------------------------------------------------------
# X-ray form factor, Cromer-Mann form:  f(s) = c + sum_j a_j exp(-b_j s^2),
# with s = sin(theta)/lambda. The numbers are copied from
# get_scattering_factor_params in src/exp_utils.f90.
#
# They are then rounded through float32, because that is what the Fortran
# effectively does: the table is declared real(dp), but its constants are
# written 2.31000 with no d0 suffix, so each is a default-real literal that is
# rounded before it is widened. Without this the reference disagrees with
# TurboGAP in the 8th significant digit for that reason alone, and the
# comparison would have to be loosened to a tolerance that hides real
# mistakes.
def _f32(x):
    return float(np.float32(x))


_XRAY = {
    "C": ([2.31000, 1.02000, 1.58860, 0.865000],
          [20.8439, 10.2075, 0.568700, 51.6512], 0.215600),
    "O": ([3.04850, 2.28680, 1.54630, 0.867000],
          [13.2771, 5.70110, 0.323900, 32.9089], 0.250800),
}
XRAY = {e: ([_f32(v) for v in a], [_f32(v) for v in b], _f32(c))
        for e, (a, b, c) in _XRAY.items()}

# Neutron scattering lengths, from get_neutron_scattering_length. That table
# does write its literals with a d0, so there is no rounding to mimic here.
NEUTRON = {"C": 6.6460, "O": 5.803}


def form_factor(element, s):
    a, b, c = XRAY[element]
    f = np.full_like(s, c)
    for aj, bj in zip(a, b):
        f = f + aj * np.exp(-bj * s * s)
    return f


def scattering_table(elements, x, neutron):
    """f_i(q) for every atom, on the grid x = 2 sin(theta)/lambda."""
    if neutron:
        return np.array([np.full_like(x, NEUTRON[e]) for e in elements])
    return np.array([form_factor(e, x / 2.0) for e in elements])


# ---------------------------------------------------------------------------
# The Debye sum
# ---------------------------------------------------------------------------
def debye(elements, positions, x, neutron=False, window=False, r_cut=0.0):
    """I(q) = 1/N sum_i sum_j f_i f_j sinc(Q r_ij) w(r_ij),  Q = 2 pi x.

    The i == j terms contribute f_i^2 -- the self term. No minimum image and
    no periodic images, matching get_xrd_debye.
    """
    n = len(elements)
    f = scattering_table(elements, x, neutron)
    Q = 2.0 * np.pi * x

    d = positions[:, None, :] - positions[None, :, :]
    r = np.sqrt((d * d).sum(-1))

    y = np.zeros_like(x)
    for i in range(n):
        for j in range(n):
            ff = f[i] * f[j]
            if i == j:
                y += ff
                continue
            if window:
                if r[i, j] > r_cut:
                    continue
                aw = np.pi * r[i, j] / r_cut
                w = np.sin(aw) / aw
            else:
                w = 1.0
            qr = Q * r[i, j]
            y += ff * np.sin(qr) / qr * w
    return y / n


def apply_output(y, elements, x, output, neutron):
    """The xrd_output conventions, as get_xrd_from_partial_structure_factors
    defines them: its y is the interference part I - sum_i c_i f_i^2."""
    if output == "xrd":
        return y
    kinds = sorted(set(elements))
    conc = {e: elements.count(e) / len(elements) for e in kinds}
    ftab = {e: scattering_table([e], x, neutron)[0] for e in kinds}
    self_term = sum(conc[e] * ftab[e] ** 2 for e in kinds)
    sth = sum(conc[e] * ftab[e] for e in kinds)
    if output in ("q*i(q)", "q*F(q)"):
        return 2.0 * np.pi * x * (y - self_term) / sth ** 2
    if output in ("F(q)", "i(q)"):
        return (y - self_term) / sth ** 2 + 1.0
    raise SystemExit(f"unknown output convention {output!r}")


def lorentz_polarization(x, wavelength, polarization, sin_theta_min=1e-3):
    """(1 + P cos^2 2theta) / (sin^2 theta cos theta), zero where there is no
    usable scattering angle."""
    sth = 0.5 * wavelength * x
    lp = np.zeros_like(x)
    ok = (sth >= sin_theta_min) & (sth < 1.0)
    s = sth[ok]
    cth = np.sqrt(1.0 - s * s)
    c2 = 1.0 - 2.0 * s * s
    lp[ok] = (1.0 + polarization * c2 * c2) / (s * s * cth)
    return lp


# ---------------------------------------------------------------------------
def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    elements, positions = [], []
    for line in lines[2:2 + n]:
        w = line.split()
        elements.append(w[0])
        positions.append([float(v) for v in w[1:4]])
    return elements, np.array(positions)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xyz")
    p.add_argument("prediction")
    p.add_argument("--qmin", type=float, default=1.0)
    p.add_argument("--qmax", type=float, default=10.0)
    p.add_argument("--n", type=int, default=51)
    p.add_argument("--units", default="q", choices=["q", "twotheta"])
    p.add_argument("--output", default="xrd")
    p.add_argument("--neutron", action="store_true")
    p.add_argument("--window", action="store_true")
    p.add_argument("--rcut", type=float, default=4.0)
    p.add_argument("--lp", action="store_true")
    p.add_argument("--polarization", type=float, default=1.0)
    p.add_argument("--wavelength", type=float, default=1.5405981)
    # xrd_prediction.dat is written F20.8, so 5e-9 per sample is the floor and
    # nothing here can be checked tighter than that.
    p.add_argument("--atol", type=float, default=1e-7)
    p.add_argument("--rtol", type=float, default=1e-10)
    a = p.parse_args()

    grid = np.linspace(a.qmin, a.qmax, a.n)
    if a.units == "twotheta":
        x = 2.0 * np.sin(np.pi * grid / 180.0 / 2.0) / a.wavelength
    else:
        x = grid / 2.0 / np.pi

    elements, positions = read_xyz(a.xyz)
    y = debye(elements, positions, x, neutron=a.neutron,
              window=a.window, r_cut=a.rcut)
    y = apply_output(y, elements, x, a.output, a.neutron)
    if a.lp:
        y = lorentz_polarization(x, a.wavelength, a.polarization) * y

    got = np.loadtxt(a.prediction, comments="#")
    if got.ndim != 2 or got.shape[0] != a.n:
        print(f"FAIL: {a.prediction} has "
              f"{0 if got.ndim != 2 else got.shape[0]} samples, expected {a.n}")
        return 1
    gx, gy = got[:, 0], got[:, 1]

    dx = np.abs(gx - grid).max()
    scale = max(float(np.abs(y).max()), 1e-30)
    dy = float(np.abs(gy - y).max())
    tol = a.atol + a.rtol * scale
    bad = dy > tol or dx > 1e-8

    print(f"    grid  n = {a.n}  max|dq| = {dx:.2e}")
    print(f"    max|y| = {scale:.6e}  max|dy| = {dy:.2e}  tol = {tol:.2e}")
    if bad:
        k = int(np.argmax(np.abs(gy - y)))
        print(f"    worst sample {k}: q = {gx[k]:.6f}  "
              f"reference = {y[k]:.10e}  turbogap = {gy[k]:.10e}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
