#!/usr/bin/env bash
#
# Does the Kokkos backend change any result?
#
# A plain Kokkos-against-CUDA diff cannot answer that on its own. The device
# binary is not reproducible run to run -- run estat_gsf twice with ONE
# unchanged binary and energy_soap moves in the tenth digit -- so several cases
# differ from themselves, and a Kokkos build would be blamed for all of them.
#
# So this runs the suite twice:
#
#   control   the plain device binary against a copy of itself. Whatever fails
#             here is device non-determinism and nothing else.
#   test      the Kokkos binary against the plain device binary.
#
# and reports the difference of the two failure sets. A case failing in the
# test pass but not in the control is the port's doing; anything in both is
# noise that was there before.
#
# Usage:
#     tools/verify_kokkos.sh                 # every case
#     tools/verify_kokkos.sh estat_gsf xrd_predict
#
# Environment:
#     TURBOGAP_GPU_BIN      plain device binary  (default bin-gpu/turbogap)
#     TURBOGAP_KOKKOS_BIN   Kokkos binary        (default bin-kokkos-gpu/turbogap)
set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

# ---------------------------------------------------------------- parameters

CUDA_BIN=${TURBOGAP_GPU_BIN:-$repo/bin-gpu/turbogap}
KOKKOS_BIN=${TURBOGAP_KOKKOS_BIN:-$repo/bin-kokkos-gpu/turbogap}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_kokkos_verify.XXXXXX")

# ----------------------------------------------------------------- functions

say() { printf '%s\n' "$*"; }

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

cleanup() {
    rm -rf "$WORK"
}

check_binaries() {
    [ -x "$CUDA_BIN" ] ||
        die "no device binary at $CUDA_BIN -- build it with ./compile_gpu.sh"
    [ -x "$KOKKOS_BIN" ] ||
        die "no Kokkos binary at $KOKKOS_BIN -- build it with KOKKOS=1 ./compile_gpu.sh"
    return 0
}

# Both passes must run the SAME binary from two paths. run.sh refuses to
# compare a binary with itself, and copying is also what stops one pass from
# reading a binary the other is rebuilding.
stage_binaries() {
    cp "$CUDA_BIN" "$WORK/tg_cuda"
    cp "$CUDA_BIN" "$WORK/tg_cuda_copy"
    cp "$KOKKOS_BIN" "$WORK/tg_kokkos"
    return 0
}

# --cpu is the bit-exact mode. That is what is wanted here even though these
# are device binaries: the question is whether ANY digit moved, and the
# control pass is what separates a moved digit from a meaningful one.
run_pass() {
    local label=$1
    local test_bin=$2
    local ref_bin=$3
    local log=$4
    shift 4

    say "=== $label ==="
    TURBOGAP_BIN="$test_bin" TURBOGAP_REF_BIN="$ref_bin" \
        "$repo/tests/regression/run.sh" --cpu "$@" > "$log" 2>&1
    say "  rc=$?  $(grep -E 'passed:' "$log" | tail -1)"
    return 0
}

failing_cases() {
    local log=$1
    local out=$2
    grep -B1 'FAIL' "$log" | grep '^== ' | tr -d '= ' | sort -u > "$out"
}

report() {
    local n_control n_test n_attributable

    failing_cases "$WORK/control.log" "$WORK/set_control.txt"
    failing_cases "$WORK/test.log" "$WORK/set_test.txt"
    comm -13 "$WORK/set_control.txt" "$WORK/set_test.txt" > "$WORK/set_attributable.txt"

    n_control=$(grep -c . < "$WORK/set_control.txt")
    n_test=$(grep -c . < "$WORK/set_test.txt")
    n_attributable=$(grep -c . < "$WORK/set_attributable.txt")

    say ""
    say "cases differing device-to-device with no code change at all: $n_control"
    sed 's/^/    /' "$WORK/set_control.txt"
    say ""
    say "cases differing Kokkos against CUDA: $n_test"
    sed 's/^/    /' "$WORK/set_test.txt"
    say ""

    if [ "$n_attributable" -eq 0 ]; then
        say "ATTRIBUTABLE TO THE PORT: none."
        say "Every case that differs under Kokkos also differs without it."
        return 0
    fi

    say "ATTRIBUTABLE TO THE PORT: $n_attributable"
    sed 's/^/    /' "$WORK/set_attributable.txt"
    say ""
    say "These differ under Kokkos and NOT in the control. Logs: $WORK"
    return 1
}

# --------------------------------------------------------------------- wiring

trap cleanup EXIT

check_binaries
stage_binaries

run_pass "control: CUDA against itself" \
    "$WORK/tg_cuda" "$WORK/tg_cuda_copy" "$WORK/control.log" "$@"
run_pass "test: Kokkos against CUDA" \
    "$WORK/tg_kokkos" "$WORK/tg_cuda" "$WORK/test.log" "$@"

# report reads the logs, so keep them if it has something to say.
if report; then
    exit 0
fi
trap - EXIT
exit 1
