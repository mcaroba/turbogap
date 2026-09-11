#!/usr/bin/env bash
#
# Topology: does the bond graph find the molecules that are there, and only
# those?
#
# Compiles kinds.f90 + elements.f90 + topology.f90 and the driver, and nothing
# else. No GAP file and no potential is involved, which is the point: topology
# is geometry and graph theory, so a failure here names the topology code
# rather than something it was standing next to.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

F90=${F90:-gfortran}
work=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_topology.XXXXXX")
trap 'rm -rf "$work"' EXIT

printf 'building against %s/src\n' "$repo"
"$F90" -ffree-line-length-none -J "$work" -c "$repo/src/kinds.f90"    -o "$work/kinds.o"
"$F90" -ffree-line-length-none -J "$work" -c "$repo/src/elements.f90" -o "$work/elements.o"
"$F90" -ffree-line-length-none -J "$work" -c "$repo/src/topology.f90" -o "$work/topology.o"
"$F90" -ffree-line-length-none -I "$work" -J "$work" \
       "$here/topologyverify.f90" "$work/kinds.o" "$work/elements.o" "$work/topology.o" \
       -o "$work/topologyverify"

printf '\n'
"$work/topologyverify"
