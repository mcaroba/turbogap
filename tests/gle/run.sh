#!/usr/bin/env bash
#
# Validation of the generalized Langevin thermostat.
#
#   gleverify   the gle module alone: the matrix exponential, the
#               semi-definite Cholesky, the closed-form ns = 0 propagator, the
#               stationary distribution across three masses, the velocity
#               autocorrelation against exp(-A t) C, the refusals, the restart,
#               and a harmonic-oscillator h-scan of the splitting bias.
#               Needs only kinds.f90 and gle.f90.
#
#   md_run.sh   the same thermostat through an actual turbogap run: keywords,
#               the report, the free-particle temperature, restart resume and
#               refusal, and the input guards. A different question from the
#               mathematics, and it fails in different ways.
#
# The two are separate because they can be run at different times: gleverify
# needs no potential and no test data and takes a couple of minutes, md_run.sh
# needs a case to run and takes longer.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: every path below is used after a cd into the build directory,
# so a relative root stops resolving the moment we move.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g}
# gle.f90 uses no other module in this tree, but it does call LAPACK's dgeev to
# report the kernel's relaxation times. Nothing the propagator does needs it.
LIBS=${LIBS:--llapack -lblas}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> gleverify (the gle module alone)"
$FC $FFLAGS -c "$SRC/kinds.f90" -o kinds.o
$FC $FFLAGS -c "$SRC/gle.f90"   -o gle.o
$FC $FFLAGS -o gleverify "$HERE/gleverify.f90" kinds.o gle.o $LIBS
./gleverify

echo
echo "==> done. The statistical checks are sampling estimates and carry a"
echo "    tolerance; everything else is an identity and is checked to"
echo "    round-off. For the wiring, run ./md_run.sh."
