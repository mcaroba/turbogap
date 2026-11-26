#!/usr/bin/env bash
#
# The pdf's sample grid and Gaussian kernel (exp_utils:make_pdf_kernel,
# pair_kde_row) against quad precision, on the grids of a 12 A and a 24.4 A pdf
# at three kernel widths: see pdf_kernelverify.f90. Linked against
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

echo "==> pdf_kernelverify (the grid and the kernel against quad precision)"
# -fopenmp: the library is threaded unless built with OPENMP=0.
$MPIFC $FFLAGS -I"$INCDIR" -o pdf_kernelverify "$HERE/pdf_kernelverify.f90" "$LIBDIR/libturbogap.a" $LIBS -fopenmp
./pdf_kernelverify
