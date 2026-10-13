#!/usr/bin/env bash
#
# The OpenMP regions in the SOAP descriptor, its derivatives and the GAP
# prediction (src/soap_turbo, src/gap.f90, src/local_properties.f90,
# src/gap_interface.f90) are written so the result does not depend on the
# thread count: every thread writes only its own sites and pairs, and the sums
# into forces and virial stay serial and in pair order.
#
# So the check is exact: a single point with forces, Hirshfeld volumes and
# SOAP derivatives, at one thread and at several, must produce byte-identical
# output. BLAS is pinned to one thread, because OpenBLAS splits dgemm by its
# own thread count and that moves the last digit independently of this code.
# A race shows up here as a difference that changes from run to run, so the
# threaded run is repeated.
#
# Usage:  ./run.sh [path/to/turbogap_root]   (needs bin/turbogap built)

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${1:-$HERE/../..}" && pwd)"
BIN=${TURBOGAP_BIN:-$ROOT/bin/turbogap}
DATA=${TURBOGAP_DATA:-$ROOT/../turbogap_tests/CO}
THREADS=${THREADS:-"2 5 8"}
REPEATS=${REPEATS:-3}
WORK="$HERE/work"

[ -x "$BIN" ] || { echo "no binary at $BIN"; exit 1; }
[ -f "$DATA/atoms_7176.xyz" ] || { echo "SKIPPED: no CO data at $DATA"; exit 0; }

run_at() {
    local threads=$1 dir=$2
    rm -rf "$dir"
    mkdir -p "$dir"
    ln -s "$DATA/gap_files" "$dir/gap_files"
    ln -s "$DATA/atoms_7176.xyz" "$dir/atoms.xyz"
    cat > "$dir/input" <<INPUT
atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
e0 = -.16138053 0.
do_forces = .true.
INPUT
    (cd "$dir" && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=$threads "$BIN" predict > out.log 2>&1)
}

echo "==> one thread"
run_at 1 "$WORK/t1"
fail=0
for t in $THREADS; do
    for r in $(seq 1 "$REPEATS"); do
        if ! run_at "$t" "$WORK/t${t}_$r"; then
            echo "    $t threads, run $r: CRASHED (see $WORK/t${t}_$r/out.log)"
            fail=1
        elif cmp -s "$WORK/t1/trajectory_out.xyz" "$WORK/t${t}_$r/trajectory_out.xyz"; then
            echo "    $t threads, run $r: identical"
        else
            echo "    $t threads, run $r: DIFFERS"
            fail=1
        fi
    done
done
[ "$fail" -eq 0 ] && echo "PASSED" || { echo "FAILED"; exit 1; }
