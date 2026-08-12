#!/usr/bin/env python3
"""Check the shape of the velocity distribution in an extxyz frame.

randomize_velocities offers two draws. "uniform" makes each component uniform
on [0,1), removes the centre-of-mass velocity and then scales everything so the
kinetic energy is exactly (3/2)(N-1)kT. "maxwell" draws each component from
N(0, kT/m) and removes the centre-of-mass velocity, with no rescaling.

Both give the right mean kinetic energy, so a temperature check cannot tell
them apart -- which is the whole reason the wrong one went unnoticed. What
separates them is the shape, and the cheapest statistic that sees it is the
excess kurtosis of the mass-scaled components sqrt(m_i) v_ij:

    Gaussian          0
    uniform on a box  -1.2

Scaling, centre-of-mass removal and the choice of temperature all leave those
values alone, so the test needs no reference data and no tolerance tuning.
"""
import sys
import argparse
import math

KB = 8.6173303e-5
AMU = 103.6426965268  # amu -> eV fs^2 / A^2, the same factor xyz.f90 divides by


def read_frame(path):
    with open(path) as fh:
        nat = int(fh.readline().split()[0])
        comment = fh.readline()
        rows = [fh.readline().split() for _ in range(nat)]
    spec = None
    for tok in comment.split():
        if tok.startswith("Properties="):
            spec = tok[len("Properties="):]
    if spec is None:
        raise SystemExit("no Properties= in comment line")
    cols, col = {}, 0
    parts = spec.split(":")
    for i in range(0, len(parts), 3):
        cols[parts[i]] = (col, int(parts[i + 2]))
        col += int(parts[i + 2])
    return cols, rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("xyz")
    p.add_argument("--shape", choices=["gaussian", "box"], required=True)
    p.add_argument("--temperature", type=float, required=True)
    p.add_argument("--temperature-tol", type=float, default=0.30,
                   help="fractional tolerance on the kinetic temperature")
    p.add_argument("--kurtosis-tol", type=float, default=0.35)
    a = p.parse_args()

    cols, rows = read_frame(a.xyz)
    for need in ("velocities", "masses"):
        if need not in cols:
            raise SystemExit("frame has no %s column" % need)

    vc, _ = cols["velocities"]
    mc, _ = cols["masses"]

    v, m = [], []
    for r in rows:
        m.append(float(r[mc]) * AMU)
        v.append([float(r[vc + k]) for k in range(3)])

    n = len(v)
    e_kin = sum(0.5 * m[i] * sum(x * x for x in v[i]) for i in range(n))
    temperature = 2.0 / 3.0 / (n - 1) / KB * e_kin

    # sqrt(m) v has the same distribution for every species, so the components
    # can be pooled across a multi-species cell.
    x = [math.sqrt(m[i]) * v[i][k] for i in range(n) for k in range(3)]
    mean = sum(x) / len(x)
    m2 = sum((u - mean) ** 2 for u in x) / len(x)
    m4 = sum((u - mean) ** 4 for u in x) / len(x)
    kurtosis = m4 / (m2 * m2) - 3.0

    want = 0.0 if a.shape == "gaussian" else -1.2

    ok = True
    dt = abs(temperature - a.temperature) / a.temperature
    if dt > a.temperature_tol:
        ok = False
    dk = abs(kurtosis - want)
    if dk > a.kurtosis_tol:
        ok = False

    print("    %d atoms" % n)
    print("    T = %.1f K (want %.0f +- %.0f%%)" % (temperature, a.temperature,
                                                    100 * a.temperature_tol))
    print("    excess kurtosis = %+.3f (want %+.1f +- %.2f)" % (kurtosis, want,
                                                                a.kurtosis_tol))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
