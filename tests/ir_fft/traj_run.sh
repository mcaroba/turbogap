#!/usr/bin/env bash
#
# The FFT estimator end to end: predict an IR spectrum from a trajectory that
# is already on disk, and check it three ways.
#
# `run.sh` tests the mathematics with no potential and no data. This tests the
# other half — the keywords, the `time=` parsing, the dipole model, the frame
# loop, the files — and it fails in different ways, so it is a separate script
# that needs a case to run.
#
# WHAT IT CHECKS
#
#   1. The dipoles TurboGAP predicts from the trajectory agree with the ones
#      the MD run that produced it recorded on each frame's comment line, to
#      within the round-off the trajectory file's F16.8 positions force -- and,
#      more to the point, DIFFER ONLY LIKE ROUND-OFF: zero-mean, white, and
#      uncorrelated with the dipole. If they do not, the predict-mode
#      descriptor path is wrong and nothing downstream is worth looking at:
#      the spectrum would be a perfectly correct transform of the wrong input.
#   2. `ir_fft_spectrum.dat` matches what `TNEP/spectroscopy.py` gives for the
#      dipoles in `ir_fft_dipoles.dat`. End to end, through the real model.
#   3. With `classical` weighting and no smoothing, it matches the block ACF
#      estimator's `ir_spectrum.dat` in SHAPE — two independent codes, so
#      agreement here is a genuine cross-validation rather than a tautology.
#
# and then reports the harmonic-vs-classical band weights, which is the
# difference that matters when reading either spectrum as physics.
#
# Usage:
#   traj_run.sh <trajectory.xyz> <potential.gap> <workdir> [turbogap] [exp.dat]
#
# The trajectory must carry `time=` on each comment line — TurboGAP's own
# `trajectory_out.xyz` does — and the potential must contain a dipole model.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

TRAJ="$(cd "$(dirname "${1:?trajectory.xyz}")" && pwd)/$(basename "$1")"
POT="$(cd "$(dirname "${2:?potential.gap}")" && pwd)/$(basename "$2")"
WORK="${3:?workdir}"
TG="${4:-$HERE/../../bin/turbogap}"
EXPFILE="${5:-}"
PY=${PY:-${TURBOGAP_PYTHON:-python3}}
NP=${NP:-1}

mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
cd "$WORK"

nframes=$(grep -c 'Lattice=' "$TRAJ")
echo "==> $nframes frames in $(basename "$TRAJ")"

ln -sf "$TRAJ" traj.xyz

# A .gap file names its companion tables (alphas_*, cutoffs, ...) by paths that
# are relative to the WORKING DIRECTORY, not to itself -- typically
# "gap_files/energy/...". So the potential's own directory has to appear under
# its own name here, and pot_file has to go through that name. Pointing
# pot_file straight at an absolute path opens the .gap perfectly well and then
# dies on the first companion table.
POTDIR="$(dirname "$POT")"
POTBASE="$(basename "$POTDIR")"
rm -f "$POTBASE"
ln -s "$POTDIR" "$POTBASE"

# `classical` and no smoothing, so that check 3 compares like with like: those
# are the two settings that make the FFT estimator compute the same quantity
# the block estimator does. Everyday use would leave them at their defaults.
cat > input <<EOF
atoms_file = "traj.xyz"
pot_file   = "$POTBASE/$(basename "$POT")"
n_species  = 3
species    = H C O
masses     = 1.008 12.011 15.999
e0         = 0. 0. 0.
soap_radial_legacy_filter = .false.

do_prediction = .true.
do_forces     = .false.

do_ir        = .true.
ir_bias_mode = fft
ir_nu_max    = 4000.0

ir_fft_acf_ratio          = 0.2
ir_fft_quantum_correction = classical
ir_fft_smooth_k           = 0
ir_fft_temperature        = 300.0
EOF

echo "==> turbogap predict  (this evaluates the dipole model on every frame)"
if [ "$NP" -gt 1 ]; then
    mpirun -np "$NP" "$TG" predict > predict.log 2>&1
else
    "$TG" predict > predict.log 2>&1
fi
grep -A10 'IR spectrum from the trajectory' predict.log | tr '\r' '\n' | head -11

echo
echo "==> 1. predicted dipoles vs the trajectory's own dipole= tags"
$PY - "$TRAJ" ir_fft_dipoles.dat <<'EOF'
import re, sys
import numpy as np

# WHAT THE TOLERANCE HAS TO BE, AND WHY IT IS NOT A GUESS.
#
# write_extxyz writes positions with F16.8, so the configurations read back are
# not the ones the MD integrated: each coordinate is quantised to 1e-8 Ang,
# uniform, sigma = 1e-8/sqrt(12) = 2.9e-9. The dipole model is a smooth
# function of those coordinates, so it inherits
#
#     sigma(dmu) ~ sqrt(3 N) * sigma_q * |dmu/dr|
#
# with |dmu/dr| a Born charge, 0.5-1.7 e for water. On 192 atoms that is
# ~7e-8 e.Ang, and the MAXIMUM over 3T samples is another factor sqrt(2 ln 3T).
# A fixed tolerance would either pass a real defect on a small system or fail
# this one on a large one.
#
# So the check is not only on the size. Position round-off is ZERO-MEAN, WHITE,
# and UNCORRELATED WITH THE DIPOLE; every way the predict path could actually
# be wrong -- a different descriptor, a botched MPI reduce, a frame offset, a
# scale error -- shows up as a bias, a drift, a correlation, or all three.
# Those are what separate noise from a defect, and they are what is asserted.
tags = []
for line in open(sys.argv[1]):
    m = re.search(r'\sdipole="([^"]+)"', line)
    if m:
        tags.append([float(x) for x in m.group(1).split()])
if not tags:
    print("   SKIP: the trajectory carries no dipole= tags to compare against")
    sys.exit(0)
ref = np.array(tags)
got = np.loadtxt(sys.argv[2])[:, 1:4]
n = min(len(ref), len(got))
ref, got = ref[:n], got[:n]
natoms = int(open(sys.argv[1]).readline().split()[0])

d = (got - ref).ravel()
sig_q = 1e-8 / np.sqrt(12.0)                 # one ulp of F16.8, uniform
expected = np.sqrt(3 * natoms) * sig_q       # per component, at |dmu/dr| = 1 e
born = d.std() / expected                    # the implied Born charge
bound = 5.0 * expected * np.sqrt(2 * np.log(max(d.size, 2)))

print(f"   {n} frames, {natoms} atoms; |mu| ~ {np.sqrt((ref**2).sum(1).mean()):.3f} e.Ang")
print(f"   max |difference|  = {np.abs(d).max():.3e}   (bound from F16.8 "
      f"positions: {bound:.1e})")
print(f"   rms difference    = {d.std():.3e}  ->  implied |dmu/dr| = {born:.2f} e")
print(f"                                          (Born charges in water: 0.5-1.7 e)")

checks = [
    ("zero-mean",              abs(d.mean()) / d.std(),                    0.05),
    ("white (lag-1 autocorr)", abs(np.corrcoef(d[:-1], d[1:])[0, 1]),      0.05),
    ("not a scale error",      abs(np.corrcoef(d, ref.ravel())[0, 1]),     0.05),
    ("no drift",               abs(d[:d.size // 2].std() / d[d.size // 2:].std() - 1), 0.10),
]
ok = np.abs(d).max() < bound and 0.1 < born < 5.0
for name, val, tol in checks:
    good = val < tol
    ok = ok and good
    print(f"   {name:<24} {val:8.4f} < {tol:<5}  {'ok' if good else 'FAIL'}")
print("   PASS: the difference is position round-off, not a disagreement"
      if ok else "   FAIL: the predict-mode dipole path disagrees")
sys.exit(0 if ok else 1)
EOF

echo
echo "==> 2/3. against spectroscopy.py and against the block estimator"
# ir_spectrum.dat, if the run that made the trajectory wrote one, is the block
# estimator's answer for the same dipoles.
for cand in "$(dirname "$TRAJ")/ir_spectrum.dat"; do
    [ -f "$cand" ] && [ ! -f ir_spectrum.dat ] && ln -sf "$cand" ir_spectrum.dat
done
$PY "$HERE/compare_estimators.py" . ${EXPFILE:+--exp "$EXPFILE"} \
    --plot spectrum.png

echo
echo "==> outputs in $WORK"
ls -la ir_fft_spectrum.dat ir_fft_dipoles.dat spectrum.png 2>/dev/null || true
