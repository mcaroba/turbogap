#!/usr/bin/env bash
#
# Validation of the MAD IR observable and the force it produces.
#
#   irverify   the mad_ir module alone: window sizing, Nyquist aliasing, peak
#              placement, and lambda = dL/dmu(newest) against finite
#              differences of the loss.
#
#   madverify  the whole chain, linked against lib/libturbogap.a so the gap.f90
#              routines under test are the ones that ship:
#              get_soap -> get_soap_dipole_weights -> get_soap_central_hessian
#              -> accumulate_dmu_dr -> mad_ir_evaluate -> mad_ir_forces,
#              against finite differences of the loss in an atom's position.
#
# madverify needs the library, so build the tree first (make) if it is not
# there. irverify needs only kinds.f90 and mad_ir.f90.
#
# Usage:  ./run.sh [path/to/turbogap]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$(cd "$HERE/../.." && pwd)}"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
MPIFC=${MPIFC:-mpif90}
FFLAGS=${FFLAGS:--O2 -g}
LIBS=${LIBS:--llapack -lblas}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> irverify (mad_ir alone)"
$FC $FFLAGS -c "$SRC/kinds.f90"  -o kinds.o
$FC $FFLAGS -c "$SRC/mad_ir.f90" -o mad_ir.o
$FC $FFLAGS -o irverify "$HERE/irverify.f90" kinds.o mad_ir.o
./irverify

echo
echo "==> madverify (full chain, against lib/libturbogap.a)"
if [ ! -f "$ROOT/lib/libturbogap.a" ]; then
  echo "    lib/libturbogap.a is missing; run make in $ROOT first." >&2
  exit 1
fi
cp -f "$HERE/../soap_derivatives/gharness.f90" .
$MPIFC $FFLAGS -I"$ROOT/include" -c gharness.f90 -o gharness.o
$MPIFC $FFLAGS -I"$ROOT/include" -o madverify "$HERE/madverify.f90" gharness.o \
       "$ROOT/lib/libturbogap.a" $LIBS
./madverify

echo
echo "==> done. Every h-scan should fall as h^2 and then turn round on"
echo "    round-off; anything flat in h is a systematic error."
