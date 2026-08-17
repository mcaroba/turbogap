#!/usr/bin/env bash
#
# Finite-difference validation of the SOAP first and second derivatives.
#
# These drive get_soap and get_soap_central_hessian directly, on a synthetic
# non-periodic system, rather than going through a fitted potential. That is
# deliberate: it reaches the uncompressed branch (no GAP in turbogap_tests was
# fitted without compression, so an end-to-end run cannot exercise it), it puts
# atoms at awkward geometries, and it isolates the descriptor from everything
# downstream of it.
#
# Every check is an h-scan. A single step size proves nothing -- a wrong
# derivative can hide inside truncation error at one h. What identifies a
# correct derivative is the shape of the scan: h^2 decay until round-off takes
# over, and a clear minimum.
#
# Usage:   ./run.sh [path/to/turbogap]
# Default is two directories up from here.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$(cd "$HERE/../.." && pwd)}"
SRC="$ROOT/src/soap_turbo/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g}
LIBS=${LIBS:--llapack -lblas}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> building against $SRC"
for f in soap_turbo_functions soap_turbo_radial soap_turbo_angular soap_turbo; do
  $FC $FFLAGS -c "$SRC/$f.f90" -o "$f.o"
done
$FC $FFLAGS -c "$HERE/gharness.f90" -o gharness.o
OBJ="gharness.o soap_turbo.o soap_turbo_radial.o soap_turbo_angular.o soap_turbo_functions.o"
for p in rdcheck dverify hverify dipverify; do
  $FC $FFLAGS -o "$p" "$HERE/$p.f90" $OBJ $LIBS
done

run () {
  echo
  echo "################################################################"
  echo "# $1"
  echo "################################################################"
  "./$1"
}

run rdcheck    # radial expansion d/dr, d2/dr2, and F_l(r)
run dverify    # soap_cart_der
run hverify    # d2 soap / d(central) d(neighbour)
run dipverify  # d(dipole)/d(r), both kernel-curvature and descriptor-curvature terms

echo
echo "==> done. What the numbers should look like:"
echo "    every h-scan falls roughly as h^2 and then turns round; the minimum"
echo "    is the answer. Anything flat in h is a systematic error, not noise."
