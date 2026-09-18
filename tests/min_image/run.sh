#!/usr/bin/env bash
#
# Validation of the nearest-image shortcut in get_distance, src/neighbors.f90.
#
# The routine used to brute-force the 27 images around the rounded fractional
# offset. It now returns the rounded one directly whenever that is closer than
# half the smallest perpendicular width of the cell, which makes it provably
# the nearest. The claim is that this is bit-exact, so the reference is the
# same routine with the shortcut switched off: two binaries from one source,
# one dump each, compared byte for byte.
#
# Negative control, run on 2026-09-18: taking the shortcut for every pair
# instead (d_near = 1.d5) makes the two dumps differ on all six cells, from
# 118 pairs in 5000 on the barely sheared one to 3106 on the skewed one. So
# the threshold is what the agreement rests on, not the rounding being right
# by accident.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: the paths below are used after a cd into the build directory.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
ST="$SRC/soap_turbo/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g -ffree-line-length-none}

SHORTCUT='res_near = d_near\*d_near'

# One source, two copies: the reference is this code with the shortcut off.
# If the anchor ever stops matching, the reference would silently become the
# same binary twice, so the counts are checked rather than assumed.
prepare_sources() {
    rm -rf "$BUILD"
    mkdir -p "$BUILD/fast" "$BUILD/full"
    cp "$SRC/neighbors.f90" "$BUILD/fast/neighbors.f90"
    sed "s/$SHORTCUT/res_near = 0.d0/" "$SRC/neighbors.f90" > "$BUILD/full/neighbors.f90"
    local n_fast n_full
    n_fast=$(grep -c "$SHORTCUT" "$BUILD/fast/neighbors.f90" || true)
    n_full=$(grep -c "$SHORTCUT" "$BUILD/full/neighbors.f90" || true)
    if [ "$n_fast" != "1" ] || [ "$n_full" != "0" ]; then
        echo "FAIL: could not switch the shortcut off (found $n_fast, left $n_full);" >&2
        echo "      get_distance has changed shape and this test needs updating" >&2
        exit 1
    fi
}

# neighbors.f90 guards its `use mpi` with _GPU, so without that define the
# chain is kinds, soap_turbo_functions and timing -- and nothing needs MPI.
build_one() {
    local variant=$1
    ( cd "$BUILD/$variant" && $FC $FFLAGS -cpp -o minimageverify \
        "$SRC/kinds.f90" "$ST/soap_turbo_functions.f90" "$SRC/nvtx.f90" "$SRC/timing.f90" \
        neighbors.f90 "$HERE/minimageverify.f90" )
}

run_one() {
    local variant=$1
    ( cd "$BUILD/$variant" && ./minimageverify )
}

compare_dumps() {
    if cmp -s "$BUILD/fast/dump.txt" "$BUILD/full/dump.txt"; then
        echo
        echo "PASS: the shortcut returns what the 27-image search returns, on every pair"
    else
        echo
        echo "FAIL: the shortcut and the full search disagree"
        diff "$BUILD/full/dump.txt" "$BUILD/fast/dump.txt" | head -10
        exit 1
    fi
}

main() {
    prepare_sources
    build_one fast
    build_one full
    echo "  pairs per cell that the shortcut can take:"
    run_one fast
    run_one full > /dev/null
    compare_dumps
}

main "$@"
