#!/usr/bin/env bash
#
# Validation of the Verlet skin, src/neighbors_skin.f90.
#
# The module decides when a neighbours list built with a padded cutoff has
# stopped holding every pair inside the real cutoff. That decision is checkable
# without any of the list machinery: build both pair sets by brute force over
# all N^2 pairs and ask whether one contains the other. The reference is the
# definition of a neighbour, not a second copy of this code and not a stored
# baseline.
#
# Three checks:
#
#   1. THE SELECTION. skin_two_largest picks the two largest displacements
#      without sorting. Compared against a full descending sort, on fields
#      built to contain a tie and a lone large mover so the branch that keeps
#      the runner-up is actually taken.
#
#   2. THE ACCUMULATOR. An atom walking steadily across a periodic boundary
#      must accumulate its path length, not a box length. The check also
#      reports what the accumulator gives WITHOUT the minimum image, and fails
#      if that agrees -- if it does, no wrap was crossed and the check has
#      stopped testing anything.
#
#   3. SAFETY. Build the padded set once, walk the atoms, and for as long as
#      skin_needs_rebuild says the list still stands, require every pair inside
#      rcut to be in that set. Run for orthorhombic and triclinic cells.
#
#      The same trajectory is also run past the criterion this module replaces
#      -- a squared displacement compared against half the buffer, unsquared --
#      and THAT one has to be caught losing pairs. A safety check no criterion
#      can fail is not a check. The run also fails if the skin survives fewer
#      than two steps on average: always rebuilding is safe and useless.
#
# Usage:  ./run.sh [path/to/turbogap_root] [seed]

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Absolute: the paths below are used after a cd into the build directory.
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
SEED="${2:-20260904}"
SRC="$ROOT/src"
BUILD="$HERE/build"

FC=${FC:-gfortran}
# -cpp because error.f90 guards its MPI use with #ifdef _MPIF90, which is not
# defined here: without the preprocessor the guard is a comment and "use mpi"
# is compiled, which fails on a build that has no MPI module.
FFLAGS=${FFLAGS:--O2 -g -cpp}

mkdir -p "$BUILD"
cd "$BUILD"

echo "==> building neighbors_skinverify"
# neighbors_skin.f90 depends on kinds.f90 and error.f90 and nothing else --
# no types, no neighbours builder, and no MPI outside _MPIF90 -- which is what
# lets a bare program drive it.
$FC $FFLAGS -o neighbors_skinverify \
    "$SRC/kinds.f90" "$SRC/error.f90" "$SRC/neighbors_skin.f90" "$HERE/neighbors_skinverify.f90"

./neighbors_skinverify "$SEED"
