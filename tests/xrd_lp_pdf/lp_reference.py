#!/usr/bin/env python3
"""Check that a predicted pattern carries the Lorentz-polarization factor.

Given the same pattern predicted twice, once with xrd_lorentz_polarization off
and once with it on, this asserts

    y_lp(q) == y_raw(q) * LP(theta),
    LP(theta) = ( 1 + P cos^2(2 theta) ) / ( sin^2(theta) cos(theta) ),

with LP set to zero below sin_theta_min and at or beyond the Ewald limit
sin(theta) >= 1. LP is recomputed here from the q grid in the file and from
nothing else, so a wrong wavelength, a dropped polarization or a factor
applied on the wrong grid all show up.

The point of differencing two runs rather than predicting y_raw from scratch
is that it isolates the factor: whatever the pdf -> structure factor -> I(q)
chain produces cancels, so this says nothing about that chain and everything
about the weight applied on top of it.

Usage: lp_reference.py <raw.dat> <lp.dat> --wavelength L --polarization P
"""
import argparse
import sys

import numpy as np


def read_pattern(path):
    d = np.loadtxt(path, comments="#")
    if d.ndim != 2 or d.shape[1] < 2:
        raise SystemExit(f"{path}: expected two columns of numbers")
    return d[:, 0], d[:, 1]


def working_grid(x, units, wavelength):
    """The file holds the user's q units; LP needs x = 2 sin(theta)/lambda."""
    if units in ("q", "saxs"):
        return x / (2.0 * np.pi)
    if units in ("twotheta", "xrd"):
        return 2.0 * np.sin(np.radians(x) / 2.0) / wavelength
    raise SystemExit(f"unknown q_units {units!r}")


def lorentz_polarization(x, wavelength, polarization, sin_theta_min):
    sth = 0.5 * wavelength * x
    ok = (sth >= sin_theta_min) & (sth < 1.0)
    lp = np.zeros_like(sth)
    s = sth[ok]
    cth = np.sqrt(1.0 - s * s)
    c2th = 1.0 - 2.0 * s * s
    lp[ok] = (1.0 + polarization * c2th * c2th) / (s * s * cth)
    return lp


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("raw", help="prediction with the factor off")
    p.add_argument("weighted", help="prediction with the factor on")
    p.add_argument("--wavelength", type=float, default=1.5405981)
    p.add_argument("--polarization", type=float, default=1.0)
    p.add_argument("--sin-theta-min", type=float, default=1e-3)
    p.add_argument("--q-units", default="q")
    # The floor is the q column, not the pattern: both are written F20.8, and
    # LP falls off as 1/sin^2(theta), so a 5e-9 error in q moves LP by
    # 2 LP / q times that -- an order of magnitude more than the 5e-9 in y
    # once LP itself is around 100. Measured worst case is 6e-8 relative;
    # 1e-6 is loose enough not to be flaky and tight enough that a wrong
    # wavelength or polarization cannot hide.
    p.add_argument("--tol", type=float, default=1e-6)
    a = p.parse_args()

    x_raw, y_raw = read_pattern(a.raw)
    x_lp, y_lp = read_pattern(a.weighted)

    if x_raw.shape != x_lp.shape or not np.allclose(x_raw, x_lp, rtol=0, atol=1e-9):
        print("    FAIL: the two runs did not use the same q grid")
        return 1

    lp = lorentz_polarization(working_grid(x_raw, a.q_units, a.wavelength),
                              a.wavelength, a.polarization, a.sin_theta_min)
    expected = y_raw * lp

    scale = max(float(np.abs(y_lp).max()), 1e-12)
    err = float(np.abs(y_lp - expected).max())
    rel = err / scale

    n_zero = int((lp == 0.0).sum())
    print(f"    LP over the grid: min {lp.min():.4f}  max {lp.max():.4f}  "
          f"zeroed samples {n_zero}/{lp.size}")
    print(f"    max |y_lp - y_raw * LP| = {err:.3e}   "
          f"relative to max|y_lp| = {scale:.4e}: {rel:.2e}   "
          f"(tolerance {a.tol:.1e})")

    if lp.max() <= 0.0:
        print("    FAIL: LP is zero over the whole grid, the check is vacuous")
        return 1
    if np.allclose(y_lp, y_raw, rtol=0, atol=1e-12):
        print("    FAIL: the two runs are identical, the factor was not applied")
        return 1
    if rel > a.tol:
        worst = int(np.argmax(np.abs(y_lp - expected)))
        print(f"    FAIL: worst at q = {x_raw[worst]:.6f}: "
              f"got {y_lp[worst]:.8f}, expected {expected[worst]:.8f} "
              f"(y_raw {y_raw[worst]:.8f} * LP {lp[worst]:.8f})")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
