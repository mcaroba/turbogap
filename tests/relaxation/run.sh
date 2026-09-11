#!/usr/bin/env bash
#
# Relaxation: does it converge, to the tolerance it was given?
#
# The regression cases run a relaxation for fifteen steps and compare the
# output byte for byte. That detects change and says nothing about whether the
# thing converges, so this drives it to convergence at several force tolerances
# and checks what should be true of the result.
#
# Needs the test data, because a relaxation needs a potential. It runs about
# seven relaxations of a 512-atom cell and takes a few minutes.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

PYTHON=${TURBOGAP_PYTHON:-python3}
BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}

# shellcheck source=../data_root.sh
. "$repo/tests/data_root.sh"
DATA=$DATA_ROOT/xps_opt

check_inputs() {
    [ -x "$BIN" ] || { printf 'ERROR: no binary at %s\n' "$BIN" >&2; exit 2; }
    [ -d "$DATA" ] || { printf 'ERROR: no test data at %s\n' "$DATA" >&2; exit 2; }
    [ -f "$DATA/gap_files/CO.gap" ] || {
        printf 'ERROR: %s/gap_files/CO.gap is missing\n' "$DATA" >&2; exit 2; }
}

# Set by run_checks, and read by the EXIT trap after it has returned -- so it
# cannot be local to it.
WORK=""
cleanup() {
    [ -n "$WORK" ] && rm -rf "$WORK"
}
trap cleanup EXIT

run_checks() {
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_relaxation.XXXXXX")
    cd "$WORK"
    TURBOGAP_BIN=$BIN \
    TURBOGAP_TOOLS=$repo/tools \
    RELAX_DATA=$DATA \
        "$PYTHON" "$here/relaxcheck.py"
}

check_inputs
printf 'binary %s\ndata   %s\n' "$BIN" "$DATA"
run_checks
