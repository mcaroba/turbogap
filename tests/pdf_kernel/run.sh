#!/usr/bin/env bash
#
# The pdf kernel's recurrence (exp_utils:pair_kde_row with step) and the exp
# it replaces against the Gaussian in quad precision, on the grids of a 12 A and
# a 24.4 A pdf at three kernel widths: see pdf_kernelverify.f90. Linked against
# lib/libturbogap.a, which make builds.
#
# Usage:  ./run.sh [path/to/tree]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
BUILD="$HERE/build"

MPIFC=${MPIFC:-mpif90}
FFLAGS=${FFLAGS:--O2 -g}
LIBS=${LIBS:--llapack -lblas}
LIBDIR=${LIBDIR:-$ROOT/lib}
INCDIR=${INCDIR:-$ROOT/include}

if [ ! -f "$LIBDIR/libturbogap.a" ]; then
  echo "    $LIBDIR/libturbogap.a is missing; run make in $ROOT first." >&2
  exit 1
fi
mkdir -p "$BUILD"
cd "$BUILD"

echo "==> pdf_kernelverify (the recurrence and exp against quad precision)"
# -fopenmp: the library is threaded unless built with OPENMP=0.
$MPIFC $FFLAGS -I"$INCDIR" -o pdf_kernelverify "$HERE/pdf_kernelverify.f90" "$LIBDIR/libturbogap.a" $LIBS -fopenmp
./pdf_kernelverify
