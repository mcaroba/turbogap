#!/usr/bin/env bash
#
# Validation of the cell list in build_neighbors_list, src/neighbors.f90.
#
# The routine has two branches. Orthorhombic boxes larger than the cutoff
# sphere take one; every other cell -- sheared, hexagonal, or smaller than its
# own cutoff -- takes the other, which until recently was a loop over every
# site against every position.
#
# What this checks is that replacing that loop with a cell list found the same
# neighbours. The reference is the loop itself, written out in the driver and
# calling the module's own get_distance, so the two differ in one thing only:
# which pairs are tested. Any disagreement is a binning error.
#
# The regression suite covers two cell shapes. This covers seven, chosen for
# the bin counts they produce: one bin per axis, two bins per axis, and enough
# bins to wrap, on a hexagonal cell, two degrees of shear and a slab.
#
# Usage:  ./run.sh [path/to/turbogap_root]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: the paths below are used after a cd into the build directory.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2 -g -ffree-line-length-none}

mkdir -p "$BUILD"
cd "$BUILD"

# neighbors.f90 guards its `use mpi` with _GPU, so without that define the
# chain is kinds, soap_turbo_functions and timing -- and nothing needs MPI.
compile() {
    $FC $FFLAGS -cpp -o neighbors_cellverify \
        "$SRC/kinds.f90" "$SRC/soap_turbo/src/soap_turbo_functions.f90" "$SRC/nvtx.f90" "$SRC/timing.f90" \
        "$SRC/neighbors.f90" "$HERE/neighbors_cellverify.f90"
}

run() {
    ./neighbors_cellverify
}

main() {
    echo "=== the cell list finds what the all-pairs loop found ==="
    compile
    run
}

main "$@"
