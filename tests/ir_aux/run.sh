#!/usr/bin/env bash
#
# Validation of the envelope-targeted resonator bank, ir_bias_mode = aux.
#
#   iraux_verify   the module alone, against things that are not it: the
#                  analytic period of a 1000 cm^-1 resonator, the closed-form
#                  steady state of a damped driven oscillator, the linearised
#                  control loop's own 2*pi*tau, an analytic dipole Jacobian by
#                  finite differences, and a restart round trip.
#
# The three checks that carry the design are 3, 4 and 5. Check 3 is that the
# controller moves R to R_target at all -- a bias POTENTIAL in R returns the
# same amplitude for every target, which is why the scheme this replaced could
# not work. Check 4 asserts both that the loop rings at 2*pi*tau without the
# proportional term and that it does not ring with it, which pins the loop's
# transfer function rather than only its fixed point. Check 5 is that active
# control does not move the resonant frequency, which is what would corrupt
# the observable.
#
# Needs only kinds.f90, mad_ir.f90 and ir_auxiliary_dynamics.f90 -- no library,
# no MPI, no GAP.
#
# Usage:  ./run.sh [path/to/turbogap]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: every path below is used after a cd into the build directory.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> building iraux_verify (kinds + mad_ir + ir_auxiliary_dynamics)"
$FC $FFLAGS -c "$SRC/kinds.f90" -o kinds.o
$FC $FFLAGS -c "$SRC/mad_ir.f90" -o mad_ir.o
$FC $FFLAGS -c "$SRC/ir_auxiliary_dynamics.f90" -o ir_auxiliary_dynamics.o
$FC $FFLAGS -o iraux_verify "$HERE/iraux_verify.f90" \
    kinds.o mad_ir.o ir_auxiliary_dynamics.o

echo
./iraux_verify
