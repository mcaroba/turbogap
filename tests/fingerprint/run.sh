#!/usr/bin/env bash
#
# Checks for src/fingerprint.f90, the array hash domain_sync_state uses to
# decide what to broadcast.
#
#   1. The arithmetic: the Fortran result for fixed inputs against the same
#      algorithm in Python's unbounded integers. A product that overflowed, or
#      a sign extension that leaked into a lane, disagrees here.
#   2. Properties: every single-bit flip, an atom swap, a reshape, a changed
#      species name or fix flag must change the fingerprint; unallocated must
#      differ from zero-size.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
BUILD="$HERE/build"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g}
PY=${PY:-${TURBOGAP_PYTHON:-python3}}

mkdir -p "$BUILD"
cd "$BUILD"
echo "==> building fingerprintverify"
$FC $FFLAGS -o fingerprintverify "$SRC/kinds.f90" "$SRC/fingerprint.f90" "$HERE/fingerprintverify.f90"
./fingerprintverify > out.txt

fail=0
echo "---- arithmetic against Python ----"
grep -E '^(REF|HASH)' out.txt | "$PY" "$HERE/reference.py" || fail=1
echo "---- properties ----"
grep -E '^(FAIL|PROPERT)' out.txt
grep -q '^PROPERTIES OK$' out.txt || fail=1

if [ "$fail" -ne 0 ]; then
    echo "FAILED"
    exit 1
fi
echo "PASSED"
