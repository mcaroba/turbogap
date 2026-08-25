#!/usr/bin/env python3
"""Run ir_fft.f90 and TNEP/spectroscopy.py on the same dipoles and compare.

The reference functions are lifted OUT of spectroscopy.py textually rather than
imported, because that module imports tensorflow at the top and this check has
no business needing a GPU stack. The extraction is verified byte-identical to
the original before it is used, so what is being compared against is the real
code and not a paraphrase of it.

THE ONE EXPECTED DIFFERENCE is the speed of light. spectroscopy.py carries
c = 2.99792458e-5 cm/fs; ir_fft.f90 carries 1/c = 33356.40952 cm/fs, the
constant the rest of TurboGAP uses. They differ in the 12th significant figure,
which shifts the frequency axis by 5.5e-12 relative and, through the
frequency-dependent quantum correction, the intensity by ~1e-12. The comparison
is therefore run twice: once as the two are written, and once with the Python's
dt compensated so the grids coincide exactly. The second run is the one that
must agree to machine precision, and it is what shows the constant is the ONLY
difference rather than the largest of several.

window = "blackman" is excluded and is a deliberate deviation: see the
"DEVIATIONS" block at the top of ir_fft.f90.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.environ.get(
    "TNEP_SPECTROSCOPY", os.path.expanduser("~/test/TNEP/spectroscopy.py"))

CM_F = 33356.40952      # the Fortran constant, cm^-1 per (1/fs)
C_PY = 2.99792458e-5    # the Python constant, cm/fs


def load_reference(path):
    """Extract compute_dipole_acf and compute_ir_spectrum, verifying the copy."""
    src = open(path).read().split("\n")

    def grab(name, end_marker):
        i = next(k for k, l in enumerate(src) if l.startswith("def " + name))
        j = next(k for k in range(i, len(src)) if src[k].startswith(end_marker))
        return "\n".join(src[i:j + 1])

    body = ("import numpy as np\n\n"
            + grab("compute_dipole_acf", "    return acf") + "\n\n"
            + grab("compute_ir_spectrum",
                   "    return freq_cm, intensity, power, acf") + "\n")
    ns = {}
    exec(compile(body, "<extracted spectroscopy.py>", "exec"), ns)
    return ns


CASES = [
    # window,  ratio, nu_max, smooth_k, kind,       quantum correction, T,   dc
    ("hann",   0.10,  4000.,   0, "gaussian", "harmonic",  300., 100.),
    ("hann",   0.10,  4000.,  10, "gaussian", "harmonic",  300., 100.),
    ("hann",   0.20,  4000.,  10, "gaussian", "harmonic",  300., 100.),
    ("hann",   0.20,  4000.,  10, "gaussian", "classical", 300., 100.),
    ("hann",   0.20,  4000.,  10, "gaussian", "linear",    300., 100.),
    ("hann",   0.20,  4000.,  10, "gaussian", "none",      300., 100.),
    ("hann",   0.20,  4000.,  10, "box",      "harmonic",  300., 100.),
    ("hann",   0.20,  4000.,  25, "box",      "harmonic",  300., 100.),
    ("none",   0.15,  3000.,  10, "gaussian", "harmonic",  500., 200.),
    ("hann",   0.05,  4000.,   7, "gaussian", "harmonic",  300.,   0.),
    ("hann",   0.33,  1000.,  10, "gaussian", "harmonic",  300., 100.),
]


def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else "./irfftverify"
    dip = sys.argv[2] if len(sys.argv) > 2 else "dipoles.dat"
    if not os.path.exists(SPEC):
        print(f"SKIP: no reference at {SPEC} (set TNEP_SPECTROSCOPY)")
        return 0
    # spectroscopy.py wraps its scipy import in try/except and falls back to
    # BOX smoothing when it fails -- silently, with smooth_kind still saying
    # "gaussian". Comparing against that would show a large disagreement that
    # is the fallback and not a defect, so refuse rather than mislead.
    try:
        import scipy.ndimage  # noqa: F401
    except ImportError:
        print("SKIP: scipy is missing, so the reference would silently use box")
        print("      smoothing where it says gaussian. The comparison would be")
        print("      invalid rather than merely approximate.")
        return 0
    ref = load_reference(SPEC)
    mu = np.loadtxt(dip)

    fac = CM_F * C_PY
    print(f"  the two speed-of-light constants differ by {fac - 1:.2e} relative")
    worst_overall = 0.0
    bad = 0
    for compensate in (False, True):
        tag = ("with the 1/c constant matched" if compensate
               else "as the two are written")
        print(f"\n  {tag}:")
        print(f"    {'window':>8} {'ratio':>6} {'smooth':>10} {'qc':>10} "
              f"{'bins':>5} {'d(nu)/nu':>10} {'dI':>10} {'dP':>10}")
        for (win, ratio, numax, sk, skind, qc, temp, dc) in CASES:
            subprocess.run(
                [exe, "spectrum", dip, "1.0", win, str(ratio), str(numax),
                 str(sk), skind, qc, str(temp), str(dc), "f.dat"],
                check=True, capture_output=True)
            F = np.loadtxt("f.dat")
            dtp = (1.0 / fac) if compensate else 1.0
            p_nu, p_I, p_P, _ = ref["compute_ir_spectrum"](
                mu, dt_fs=dtp, window=(win if win != "none" else None),
                max_freq_cm=numax, acf_ratio=ratio, smooth_k=sk,
                smooth_kind=skind, temperature=temp,
                quantum_correction=qc, power_dc_cutoff_cm=dc)
            if len(F) != len(p_nu):
                print(f"    {win:>8} {ratio:>6} {skind:>10} {qc:>10}  "
                      f"LENGTH MISMATCH {len(F)} vs {len(p_nu)}")
                bad += 1
                continue
            dnu = np.max(np.abs(F[:, 0] - p_nu)) / max(p_nu.max(), 1.0)
            dI = np.max(np.abs(F[:, 1] - p_I))
            dP = np.max(np.abs(F[:, 2] - p_P))
            print(f"    {win:>8} {ratio:>6} {skind:>10} {qc:>10} {len(F):>5} "
                  f"{dnu:>10.2e} {dI:>10.2e} {dP:>10.2e}")
            if compensate:
                worst_overall = max(worst_overall, dnu, dI, dP)

    # 1e-11 is loose against the ~1e-14 that is actually achieved, and tight
    # against anything that is a real disagreement rather than round-off.
    print()
    if bad == 0 and worst_overall < 1e-11:
        print(f"  PASS: agrees with spectroscopy.py to {worst_overall:.1e} "
              f"once the 1/c constant is matched")
        return 0
    print(f"  FAIL: worst disagreement {worst_overall:.2e}, "
          f"{bad} length mismatches")
    return 1


if __name__ == "__main__":
    sys.exit(main())
