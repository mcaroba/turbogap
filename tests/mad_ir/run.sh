#!/usr/bin/env bash
#
# Validation of the MAD IR observable and the force it produces.
#
#   irverify   the mad_ir module alone: window sizing, Nyquist aliasing, peak
#              placement, and lambda = dL/dmu(newest) against finite
#              differences of the loss.
#
#   auxverify  the same module in AUXILIARY-VARIABLE form (ir_acf_mode =
#              exponential), where the running correlation is integrated
#              alongside the atoms rather than recomputed from the buffer. The
#              gradient is entirely different code from irverify's and gets its
#              own h-scan, over both the mean subtraction and the fitted scale.
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
# Absolute: every path below is used after a cd into the build directory,
# so a relative root stops resolving the moment we move.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
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
echo "==> auxverify (mad_ir with auxiliary variables)"
$FC $FFLAGS -o auxverify "$HERE/auxverify.f90" kinds.o mad_ir.o
./auxverify

echo
echo "==> madverify (full chain, against lib/libturbogap.a)"
# LIBDIR/INCDIR let a variant object tree (lib-dbg, lib-gle, ...) be tested
# without rebuilding the default one that everything else names.
LIBDIR=${LIBDIR:-$ROOT/lib}
INCDIR=${INCDIR:-$ROOT/include}
if [ ! -f "$LIBDIR/libturbogap.a" ]; then
  echo "    $LIBDIR/libturbogap.a is missing; run make in $ROOT first." >&2
  exit 1
fi
cp -f "$HERE/../soap_derivatives/gharness.f90" .
$MPIFC $FFLAGS -I"$INCDIR" -c gharness.f90 -o gharness.o
$MPIFC $FFLAGS -I"$INCDIR" -o madverify "$HERE/madverify.f90" gharness.o \
       "$LIBDIR/libturbogap.a" $LIBS
./madverify

echo
echo "==> done. Every h-scan should fall as h^2 and then turn round on"
echo "    round-off; anything flat in h is a systematic error."
