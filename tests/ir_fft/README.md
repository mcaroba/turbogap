# `ir_fft` — the FFT IR estimator

Tests for `src/ir_fft.f90` and `src/ir_fft_io.f90`: the Wiener–Khinchin IR
pipeline translated from `TNEP/spectroscopy.py`, selected with
`ir_bias_mode = fft`, and the `turbogap predict` route that turns a trajectory
on disk into a spectrum.

Operational documentation is `docs/mad_ir_matching.md` §1.5, §1.6, §2.4, §2.5
and §4.4.

## What makes these tests different from the other IR ones

`acf` and `xl` can only be checked against their own derivations. `fft` has a
**reference implementation**, so agreement is checkable to the last bit rather
than to a tolerance somebody chose. That is what `run.sh` step 4 does, and it
is the strongest statement available about this module: **machine precision,
3.5e-14 worst case over 11 configurations**, once the two speed-of-light
constants are matched.

## The scripts

| script | needs | what it does |
|---|---|---|
| `run.sh` | gfortran, numpy, scipy | the mathematics, with no potential and no data |
| `traj_run.sh` | a built `turbogap`, a trajectory, a dipole model | the same code through an actual run |
| `compare_python.py` | numpy, scipy | called by `run.sh`; the head-to-head with the Python |
| `compare_estimators.py` | numpy (matplotlib to plot) | called by `traj_run.sh`; also usable on any run directory |
| `gen_dipoles.py` | numpy | the synthetic trajectory `run.sh` uses |
| `irfftverify.f90` | — | the driver; `spectrum`, `acf`, `grad`, `fft` modes |

```
./run.sh                                  # from anywhere; defaults to ../..
./traj_run.sh traj.xyz model.gap /tmp/w  ../../bin/turbogap  exp_ir.dat
```

`run.sh` and `traj_run.sh` honour `TURBOGAP_PYTHON`, because on some hosts the
system interpreter has no numpy.

**`compare_python.py` needs scipy, and skips loudly without it.** This is not
fussiness: `spectroscopy.py` wraps its scipy import in `try/except` and falls
back to **box** smoothing when it fails, silently, while still reporting
`smooth_kind = "gaussian"`. Comparing against that would show a large
disagreement that is the fallback and not a defect.

## What each check is actually for

**1. The transform.** Round trip and forward-against-the-definition on the
radix-2 FFT. Nothing else can be right if this is not.

**2. `C(tau)` against the definition, not only against the Python.** The Python
computes the correlation by FFT too, so an error in the zero padding would
agree with itself perfectly and be wrong in both. The direct `O(T·L)` sum is
the only independent statement of what the correlation is.

**3. The pipeline against `spectroscopy.py`,** over every window, every quantum
correction, both smoothers and several `acf_ratio`s. The two functions are
lifted out of the module textually — verifying the copy is byte-identical
first — so the comparison is against the real code, and so it does not import
TensorFlow.

The comparison runs **twice**: as the two are written, and with the Python's
`dt` compensated so the frequency grids coincide exactly. `spectroscopy.py`
carries `c = 2.99792458e-5` cm/fs and TurboGAP carries `1/c = 33356.40952`,
which differ in the twelfth significant figure (5.5e-12 relative). Running it
both ways is what shows that constant is the **only** difference rather than
the largest of several.

**4. The gradient, by h-scan.** Not a single `h`. A gradient wrong by a
constant factor passes any single-`h` check with a loose enough tolerance and
fails a scan immediately, so what is asserted is the **ratio**: the error must
fall by a factor of four per halving of `h`, which is the central-difference
truncation error and nothing else. Measured 4.00, 4.00, 4.00.

## The one deliberate deviation from the Python

`window = "blackman"`, and it is excluded from the comparison.

`np.blackman(L)` is a **symmetric** window. Applied to lags `0..L-1` it is
≈0 at lag 0, peaks at `L/2`, and returns to ≈0 at `L-1`. As a lag window that
deletes `C(0)` — the largest and best determined term, the one carrying the
total intensity — and weights the middle of the correlation most. `ir_fft` uses
the half window instead,

    w(tau) = 0.42 + 0.5 cos(pi tau/L) + 0.08 cos(2 pi tau/L)

which is 1 at lag 0 and 0 at lag `L`, and is what "Blackman lag window" means
everywhere else. `hann` and `none` agree with the Python exactly.

## `traj_run.sh`, and why the dipole check comes first

A trajectory written by TurboGAP carries `dipole="..."` on each comment line.
Predicting the dipoles again from the same configurations with the same model
must reproduce those numbers, and the check is first because if it fails
nothing downstream is worth reading — the spectrum would be a perfectly correct
transform of the wrong input. Measured agreement on 64 waters: exact to the
8 decimals the tag is written with.

Then the spectrum is checked against `spectroscopy.py` on the dipoles TurboGAP
itself wrote (`ir_fft_dipoles.dat` exists for exactly this), and against
`mad_ir.f90`'s block estimator with the settings that make the two comparable
— `ir_fft_quantum_correction = classical`, `ir_fft_smooth_k = 0`. Those are two
genuinely different codes, one an FFT and one a direct lag sum, on different
frequency grids, so agreement there is a cross-validation rather than a
tautology. They agree only up to an overall factor, which is fitted; what is
compared is the shape.
