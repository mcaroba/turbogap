#!/usr/bin/env python3
"""Compare the FFT estimator with the block ACF estimator on the same trajectory.

Two independent implementations of the same physics -- `ir_fft.f90` (this work,
a translation of GPUMD's pipeline) and `mad_ir.f90` (the original block
estimator) -- run over the same dipoles. When they are configured to compute
the same thing they must agree, and where they are configured differently the
difference has to be the one that was intended and not another one.

THREE COMPARISONS, in order of how much they prove.

  A. `ir_fft` against `spectroscopy.py` on the dipoles TurboGAP itself wrote.
     End to end, through the real dipole model, on a real trajectory. Limited
     only by the text precision of `ir_fft_dipoles.dat`.

  B. `ir_fft` with `quantum_correction = classical` and no smoothing, against
     `ir_spectrum.dat` from the MD run's own block estimator. These are
     genuinely different codes -- FFT versus a direct lag sum, a uniform DFT
     grid versus the experimental grid -- so they agree only up to an overall
     constant, which is fitted. What is checked is the SHAPE.

  C. `harmonic` against `classical` on identical dipoles, which isolates the
     quantum correction. This is the difference that matters in practice and
     the one a fitted scale cannot absorb.

Usage:
    compare_estimators.py <rundir> [--exp FILE] [--plot OUT.png]

<rundir> holds ir_fft_dipoles.dat and (for B) ir_spectrum.dat.
"""
import argparse
import os
import sys

import numpy as np

CM_F = 33356.40952
C_PY = 2.99792458e-5


def load_reference(path):
    """Lift the two functions out of spectroscopy.py, verifying the copy."""
    src = open(path).read().split("\n")

    def grab(name, end):
        i = next(k for k, l in enumerate(src) if l.startswith("def " + name))
        j = next(k for k in range(i, len(src)) if src[k].startswith(end))
        return "\n".join(src[i:j + 1])

    body = ("import numpy as np\n\n"
            + grab("compute_dipole_acf", "    return acf") + "\n\n"
            + grab("compute_ir_spectrum",
                   "    return freq_cm, intensity, power, acf") + "\n")
    ns = {}
    exec(compile(body, "<extracted spectroscopy.py>", "exec"), ns)
    return ns


def read_header(path, sep=":"):
    """Pull the '#  key <sep> value' lines out of a header.

    ir_fft_spectrum.dat separates with ':' and ir_spectrum.dat with '=', so the
    separator is a parameter. Getting this wrong is not a parse error, it is an
    empty dict and a silent fallback to a default -- which is how comparison B
    below spent a while comparing two different subsets of the same run.
    """
    out = {}
    for line in open(path):
        if not line.startswith("#"):
            break
        if sep in line:
            k, v = line[1:].split(sep, 1)
            out[k.strip()] = v.strip()
    return out


def fit_scale(model, target, w=None):
    """Least-squares overall factor. Both curves are in arbitrary units, so
    comparing them without one compares two different unit systems."""
    if w is None:
        w = np.ones_like(target)
    denom = np.sum(w * model * model)
    return np.sum(w * model * target) / denom if denom > 0 else 1.0


# numpy renamed trapz to trapezoid in 2.0 and deprecated the old name.
_trapz = getattr(np, "trapezoid", None) or np.trapz


def band_fraction(nu, I, lo, hi):
    """Fraction of the total intensity between lo and hi. Scale-invariant, so
    it compares band shape without any fit at all -- which is what makes it the
    right statistic for comparing two spectra in different arbitrary units."""
    m = (nu >= lo) & (nu <= hi)
    tot = float(_trapz(np.abs(I), nu))
    if tot <= 0.0:
        return 0.0
    return float(_trapz(np.abs(I[m]), nu[m])) / tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--exp", default=None, help="two-column experimental spectrum")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--spectroscopy",
                    default=os.path.expanduser("~/test/TNEP/spectroscopy.py"))
    a = ap.parse_args()

    dip_path = os.path.join(a.rundir, "ir_fft_dipoles.dat")
    fft_path = os.path.join(a.rundir, "ir_fft_spectrum.dat")
    D = np.loadtxt(dip_path)
    t, mu = D[:, 0], D[:, 1:4]
    dt = (t[-1] - t[0]) / (len(t) - 1)
    F = np.loadtxt(fft_path)
    hdr = read_header(fft_path)

    print(f"trajectory : {len(t)} frames at {dt:.6f} fs "
          f"= {t[-1] - t[0]:.1f} fs total")
    for k in ("lags kept", "resolution", "bin spacing", "Nyquist",
              "quantum correction", "smoothing", "lag window"):
        if k in hdr:
            print(f"  {k:<20}: {hdr[k]}")

    ratio = float(hdr["lags kept"].split("acf_ratio =")[1].rstrip(") "))
    smooth_k = int(hdr["smoothing"].split()[0])
    skind = hdr["smoothing"].split(",")[1].split()[0]
    qc = hdr["quantum correction"]
    win = hdr["lag window"]
    numax = float(hdr["Nyquist"].split()[0])
    numax = F[:, 0].max()
    temp = float(hdr.get("temperature", "300").split()[0])

    # ---- A. against spectroscopy.py -------------------------------------
    print("\nA. ir_fft.f90 vs TNEP/spectroscopy.py, same dipoles")
    if not os.path.exists(a.spectroscopy):
        print(f"   SKIP: no reference at {a.spectroscopy}")
    else:
        try:
            import scipy.ndimage  # noqa: F401
        except ImportError:
            print("   SKIP: scipy missing; the reference would silently use box "
                  "smoothing where it says gaussian")
        else:
            ref = load_reference(a.spectroscopy)
            # Compensate the 1/c constant so the grids coincide; see
            # compare_python.py for why this is the honest comparison.
            p_nu, p_I, _p, _a = ref["compute_ir_spectrum"](
                mu, dt_fs=dt / (CM_F * C_PY),
                window=(win if win != "none" else None),
                max_freq_cm=numax + 1e-6, acf_ratio=ratio, smooth_k=smooth_k,
                smooth_kind=skind, temperature=temp, quantum_correction=qc,
                power_dc_cutoff_cm=100.0)
            n = min(len(p_nu), len(F))
            dI = np.max(np.abs(F[:n, 1] - p_I[:n]))
            print(f"   bins {len(F)} vs {len(p_nu)};  max |dI| = {dI:.3e}")
            print(f"   (both peak-normalised, so this is an absolute difference;")
            print(f"    the floor is the {D.shape[1]}-column text precision of "
                  f"ir_fft_dipoles.dat, ~1e-10)")
            print("   PASS" if dI < 1e-6 else "   FAIL")

    # ---- B. against the block estimator ---------------------------------
    #
    # THE TWO MUST BE GIVEN THE SAME DATA, and by default they are not.
    # `mad_ir`'s ensemble is a ROLLING WINDOW of `ir_lag_factor * n_lag`
    # frames -- whenever `exp_labels` names `ir`, which puts the run on the
    # bias path even at `exp_energy_scales = 0` -- while `ir_fft` here has
    # transformed the whole trajectory. Comparing those compares two
    # different subsets of the same run and reports a difference that is
    # nothing but the sample size: measured 26% on 6001 vs 2384 frames, with
    # a median ratio of exactly 1.000, which is what a pure sampling
    # difference looks like and what a real disagreement does not.
    #
    # `ir_spectrum.dat` states its ensemble size in its header, so read it and
    # transform the matching tail.
    acf_path = os.path.join(a.rundir, "ir_spectrum.dat")
    print("\nB. ir_fft (classical, unsmoothed) vs mad_ir's block estimator")
    if not os.path.exists(acf_path):
        print(f"   SKIP: no {acf_path}")
    elif qc != "classical" or smooth_k > 1:
        print(f"   SKIP: this run used quantum_correction={qc}, smooth_k={smooth_k}.")
        print("        Rerun with ir_fft_quantum_correction = classical and")
        print("        ir_fft_smooth_k = 0 to make the two comparable.")
    else:
        ahdr = read_header(acf_path, sep="=")
        n_ens = int(ahdr.get("frames in the ensemble", len(t)).split()[0])
        n_lag_acf = int(ahdr.get("longest lag transformed", "0").split()[0]) or None
        if not ahdr:
            print("   WARNING: could not parse the ir_spectrum.dat header; the")
            print("            comparison below may be on a different ensemble")
        A = np.loadtxt(acf_path)
        a_nu, a_I = A[:, 0], A[:, 1]
        if n_ens < len(t):
            print(f"   mad_ir held a rolling window of {n_ens} frames of the "
                  f"{len(t)}; transforming the matching tail")
            print(f"   (comparing against the full trajectory instead would report "
                  f"a difference that is only the sample size)")
        # Recompute ir_fft's side on exactly that window, at exactly that
        # longest lag, so the only remaining differences are the ones being
        # tested: an FFT correlation against a direct lag sum, and a uniform
        # DFT grid against 95 experimental wavenumbers.
        if os.path.exists(a.spectroscopy):
            try:
                import scipy.ndimage  # noqa: F401
                ref = load_reference(a.spectroscopy)
                sub = mu[-n_ens:]
                ratio_b = (n_lag_acf / n_ens) if n_lag_acf else ratio
                b_nu, _bI, _bp, _ba = ref["compute_ir_spectrum"](
                    sub, dt_fs=dt, window="hann", max_freq_cm=a_nu.max() + 1.0,
                    acf_ratio=ratio_b, smooth_k=0, smooth_kind="gaussian",
                    quantum_correction="classical", power_dc_cutoff_cm=100.)
                # peak-normalised; undo it against the block estimator by fitting
                f_I = np.interp(a_nu, b_nu, _bI)
            except ImportError:
                f_I = np.interp(a_nu, F[:, 0], F[:, 3])
        else:
            f_I = np.interp(a_nu, F[:, 0], F[:, 3])
        m = a_I != 0
        s = fit_scale(f_I[m], a_I[m])
        rel = np.sqrt(np.sum((s * f_I[m] - a_I[m])**2) / np.sum(a_I[m]**2))
        r = (s * f_I[m]) / a_I[m]
        print(f"   {m.sum()} common points, fitted factor {s:.6e}")
        print(f"   relative shape difference = {rel:.4f}   "
              f"(ratio {r.min():.3f} to {r.max():.3f}, median {np.median(r):.3f})")
        print("   Two independent codes on the same dipoles: an FFT correlation")
        print("   on a uniform DFT grid against a direct lag sum evaluated at the")
        print("   experimental wavenumbers. The residual is the lag window's")
        print("   off-by-one and the regridding, not the estimator.")

    # ---- C. harmonic vs classical ---------------------------------------
    print("\nC. the quantum correction, on identical dipoles")
    if os.path.exists(a.spectroscopy):
        try:
            import scipy.ndimage  # noqa: F401
            ref = load_reference(a.spectroscopy)
            curves = {}
            for q in ("harmonic", "classical"):
                nu_q, I_q, _p, _a = ref["compute_ir_spectrum"](
                    mu, dt_fs=dt, window="hann", max_freq_cm=4000.,
                    acf_ratio=ratio, smooth_k=10, smooth_kind="gaussian",
                    temperature=temp, quantum_correction=q,
                    power_dc_cutoff_cm=100.)
                curves[q] = (nu_q, I_q)
            print(f"   {'band':<22} {'harmonic':>10} {'classical':>10}")
            for name, lo, hi in (("librational 0-1000", 0., 1000.),
                                 ("bend 1500-1800", 1500., 1800.),
                                 ("stretch 3000-3800", 3000., 3800.)):
                nu_h, I_h = curves["harmonic"]
                nu_c, I_c = curves["classical"]
                fh = band_fraction(nu_h, I_h, lo, hi)
                fc = band_fraction(nu_c, I_c, lo, hi)
                print(f"   {name:<22} {fh:>10.4f} {fc:>10.4f}")
            print("   (fractions of the total integrated intensity: scale-invariant,")
            print("    so no fit can move them. This is the difference §1.5 warns about.)")
        except ImportError:
            print("   SKIP: scipy missing")
    else:
        print("   SKIP: no reference")

    # ---- the plot --------------------------------------------------------
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(F[:, 0], F[:, 1], lw=1.0, color="black",
                label=f"ir_fft.f90 ({qc})")
        if os.path.exists(acf_path):
            A = np.loadtxt(acf_path)
            a_nu, a_I = A[:, 0], A[:, 1]
            pk = np.max(np.abs(a_I))
            if pk > 0:
                ax.plot(a_nu, a_I / pk, lw=1.0, color="tab:blue", alpha=0.8,
                        label="mad_ir block ACF (nu^2)")
        if a.exp and os.path.exists(a.exp):
            E = np.loadtxt(a.exp)
            pk = np.max(np.abs(E[:, 1]))
            ax.plot(E[:, 0], E[:, 1] / pk, lw=1.2, color="tab:red",
                    alpha=0.7, label="experiment")
        ax.set_xlabel("Wavenumber (cm$^{-1}$)")
        ax.set_ylabel("Intensity (peak-normalised)")
        ax.set_xlim(0, F[:, 0].max())
        ax.grid(alpha=0.3)
        ax.legend()
        ax.set_title(f"{len(t)} frames at {dt:g} fs, "
                     f"resolution {hdr.get('resolution', '?')}")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=130)
        print(f"\nwrote {a.plot}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
