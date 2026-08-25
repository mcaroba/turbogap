#!/usr/bin/env bash
#
# Validation of the FFT IR estimator, src/ir_fft.f90.
#
# ir_fft.f90 is a translation of TNEP/spectroscopy.py, so the question it has
# to answer is not "does this look like a spectrum" but "is it the same
# spectrum". That makes the tests unusually sharp: there is a reference
# implementation, and agreement with it is checkable to the last bit.
#
# Four checks, in order of how little they assume:
#
#   1. THE TRANSFORM. A round trip through the radix-2 FFT must be the
#      identity, and the forward transform must equal the DFT written out as
#      a sum. Nothing else in the module can be right if this is not.
#
#   2. THE CORRELATION. C(tau) from the FFT against C(tau) from the
#      definition. This is checked against the DEFINITION and not only
#      against the Python, because the Python computes it by FFT too -- an
#      error in the zero padding would agree with itself perfectly and be
#      wrong in both.
#
#   3. THE PIPELINE, against spectroscopy.py, over a grid of configurations:
#      every window, every quantum correction, both smoothers, several
#      acf_ratios. Needs numpy and scipy; skipped with a message if they are
#      not there, because a skipped check that says so is better than a check
#      that quietly does not run.
#
#   4. THE GRADIENT, by an h-scan against central differences. Not a single h:
#      a gradient wrong by a constant factor passes any single-h check with a
#      loose enough tolerance, and fails a scan immediately. What is being
#      looked for is the error falling as h^2 until round-off takes over.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: every path below is used after a cd into the build directory, so a
# relative root stops resolving the moment we move.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g}
# Honours TURBOGAP_PYTHON, as the other test scripts in this tree do: on some
# hosts the system interpreter has no numpy and the one that does is in a venv.
PY=${PY:-${TURBOGAP_PYTHON:-python3}}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> building irfftverify"
# ir_fft.f90 depends on nothing but kinds.f90 -- no MPI, no LAPACK, no types --
# which is what makes it drivable from a bare program like this one.
$FC $FFLAGS -o irfftverify \
    "$SRC/kinds.f90" "$SRC/ir_fft.f90" "$HERE/irfftverify.f90"

fail=0
note() { echo; echo "---- $* ----"; }

note "1. the radix-2 transform"
./irfftverify fft 4096 || fail=1

note "2. the dipole trajectory"
$PY "$HERE/gen_dipoles.py" dipoles.dat || {
    echo "SKIP: $PY has no numpy; set TURBOGAP_PYTHON to one that does"
    exit 0
}

note "3. C(tau): FFT against the definition"
./irfftverify acf dipoles.dat acf_check.dat || fail=1

note "4. the pipeline against TNEP/spectroscopy.py"
# scipy is not optional here even though spectroscopy.py treats it as optional:
# without it the reference SILENTLY falls back from gaussian smoothing to box,
# and the comparison would report a large disagreement that is the fallback
# rather than a defect. Better to skip loudly. See compare_python.py.
if $PY -c "import numpy, scipy.ndimage" 2>/dev/null; then
    $PY "$HERE/compare_python.py" ./irfftverify dipoles.dat || fail=1
else
    echo "SKIP: numpy and scipy are needed to run the reference implementation"
    echo "      (scipy especially: without it spectroscopy.py quietly switches"
    echo "       gaussian smoothing to box, and the comparison is then invalid)"
fi

note "5. the MAD gradient, by h-scan"
$PY - <<'EOF'
import numpy as np
nu = np.linspace(400., 4000., 95)
I  = (np.exp(-0.5*((nu-3400.)/180.)**2) + 0.5*np.exp(-0.5*((nu-1650.)/90.)**2)
      + 0.3*np.exp(-0.5*((nu-650.)/220.)**2))
np.savetxt('exp_target.dat', np.column_stack([nu, I]), fmt='%.10e')
EOF
for cfg in "harmonic gaussian mean" "classical box mean" "harmonic gaussian nomean"; do
    echo "  configuration: $cfg"
    ./irfftverify grad dipoles.dat 1.0 exp_target.dat 1e-1 $cfg > grad.out || fail=1
    grep -E '^  -> ' grad.out
    # The scan must fall like h^2 over the first few halvings. Checking the
    # RATIO rather than an absolute tolerance is what makes this a test of the
    # gradient rather than of the tolerance: a wrong gradient gives a ratio of
    # 1 (the error is the constant discrepancy, not the truncation).
    $PY - grad.out <<'EOF'
import sys, re
errs = [float(m.group(1)) for m in
        (re.search(r'worst rel\.err =\s+(\S+)', l) for l in open(sys.argv[1]))
        if m]
ok = True
for i in range(3):
    r = errs[i]/errs[i+1]
    tag = "ok" if 3.0 < r < 5.0 else "NOT h^2"
    if tag != "ok": ok = False
    print(f"    err({i}) / err({i+1}) = {r:6.2f}   (h^2 predicts 4)   {tag}")
print("    PASS" if ok else "    FAIL: the error does not fall as h^2")
sys.exit(0 if ok else 1)
EOF
done

echo
if [ "$fail" -eq 0 ]; then
    echo "ir_fft: all checks passed"
else
    echo "ir_fft: FAILURES above"
    exit 1
fi
