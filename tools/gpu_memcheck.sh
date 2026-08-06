#!/usr/bin/env bash
#
# Run a regression case under compute-sanitizer.
#
# The bugs this branch has actually shipped are all in the sanitizer's scope
# and none of them announced themselves:
#
#   * a device buffer allocated and never filled, then used to index another
#     -- negative index, out-of-bounds atomics, wrong answers (the missing
#     cpy_htod for neighbors_list_d, found by a structural bisect after an
#     argument-by-argument comparison could not find it, because every
#     argument was byte-identical and the bug was a MISSING STATEMENT);
#   * a host-side overrun that counted pairs with < and filled with <=, which
#     on the 7176-atom CO system lands in memory the allocator already handed
#     out and corrupts the heap silently;
#   * a debugging exit(0) that terminated the process with status 0.
#
# initcheck is the tool for the first, memcheck for the second. Neither is
# reachable from the regression suite, which only compares outputs.
#
# Usage:
#   tools/gpu_memcheck.sh                        # memcheck, vdw_ts (small cell)
#   tools/gpu_memcheck.sh --tool initcheck
#   tools/gpu_memcheck.sh --case CO_predict --tool racecheck
#   tools/gpu_memcheck.sh --leak
#   tools/gpu_memcheck.sh --list
#
# Tools:
#   memcheck    out-of-bounds and misaligned device accesses    (start here)
#   initcheck   device reads of UNINITIALISED device memory     (the bisect bug)
#   racecheck   shared-memory data races within a block
#   synccheck   invalid __syncthreads / barrier usage
#
# PREFER A SMALL CELL. The default case is the P4 dimer, deliberately: a suite
# of one system size cannot see a whole class of bug, and the overrun above
# aborted instantly on 8 atoms while corrupting the heap invisibly on 7176.
# compute-sanitizer is also 10-100x slower, so CO_predict under memcheck is a
# coffee break.

set -eu

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
DATA_ROOT=${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}
TOOL=memcheck
CASE=vdw_ts
LEAK=0
KEEP=${TURBOGAP_KEEP:-0}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --tool) TOOL=$2; shift 2 ;;
    --case) CASE=$2; shift 2 ;;
    --leak) LEAK=1; shift ;;
    --keep) KEEP=1; shift ;;
    --list) printf 'cases: vdw_ts (P4 dimer, 8 atoms)  CO_predict (7176 atoms)  XRD_mad\n'
            printf 'tools: memcheck initcheck racecheck synccheck\n'; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

# Prefer the copy inside the CUDA toolkit. The /usr/bin/compute-sanitizer on
# alt is a stub that cannot find libsanitizer-collection.so and fails with
# "Unable to find injection library" -- and it fails while still producing
# plausible-looking output, so a PATH lookup that takes the first hit reports
# a clean run that never ran.
SANITIZER=""
for c in /usr/local/cuda/bin /usr/local/cuda-*/bin; do
  [ -x "$c/compute-sanitizer" ] && { SANITIZER="$c/compute-sanitizer"; break; }
done
[ -n "$SANITIZER" ] || SANITIZER=$(command -v compute-sanitizer 2>/dev/null || true)
[ -n "$SANITIZER" ] || die "compute-sanitizer not found (looked in /usr/local/cuda*/bin and on PATH)"
[ -x "$BIN" ] || die "binary not found: $BIN (build with DEBUG=0)"

WORK=${TMPDIR:-/tmp}/turbogap_memcheck.$$
mkdir -p "$WORK"
[ "$KEEP" = 1 ] || trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# ---- stage the case, mirroring tests/gpu/run_regression.sh -----------------
case "$CASE" in
  vdw_ts)
    d=$DATA_ROOT/vdw_P
    [ -d "$d" ] || die "test data not found: $d (run tests/fetch_test_data.sh)"
    cp "$d/p4_dimer.xyz" atoms.xyz
    ln -sf "$d/gap_files" gap_files
    cat > input <<'EOF'
atoms_file = "atoms.xyz"
pot_file = "gap_files/phosphorus.gap"
n_species = 1
species = P
e0 = -0.52375977
masses = 30.97
random_seed = 12345
vdw_rcut = 25.
vdw_r0_ref = 2.12
vdw_alpha0_ref = 3.7046
vdw_c6_ref = 110.54
vdw_buffer = 0.5
vdw_type = ts
write_xyz = 1
EOF
    MODE=predict ;;
  CO_predict)
    d=$DATA_ROOT/CO
    [ -d "$d" ] || die "test data not found: $d"
    cp "$d/atoms.xyz" atoms.xyz 2>/dev/null || cp "$d"/atoms_*.xyz atoms.xyz
    ln -sf "$d/gap_files" gap_files
    cat > input <<'EOF'
atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
e0 = 0. 0.
masses = 12.01 16.00
random_seed = 12345
write_xyz = 1
EOF
    MODE=predict ;;
  *) die "unknown case: $CASE (try --list)" ;;
esac

SAN_ARGS="--tool $TOOL --print-limit 40"
[ "$TOOL" = memcheck ] && [ "$LEAK" = 1 ] && SAN_ARGS="$SAN_ARGS --leak-check full"
# --padding pads every allocation so a small overrun lands in guard bytes
# rather than in the next live buffer. This is what makes an overrun VISIBLE
# on a large system instead of silently corrupting a neighbouring array.
[ "$TOOL" = memcheck ] && SAN_ARGS="$SAN_ARGS --padding 64"

printf '=== compute-sanitizer --tool %s : case %s (%s) ===\n' "$TOOL" "$CASE" "$MODE"
printf '    sanitizer %s\n    binary %s\n    workdir %s\n\n' "$SANITIZER" "$BIN" "$WORK"

set +e
"$SANITIZER" $SAN_ARGS "$BIN" "$MODE" > sanitizer.log 2>&1
rc=$?
set -e

# The error COUNT is the verdict, not the exit status: compute-sanitizer
# reports errors and still exits 0 in some configurations, which is the same
# "reports failure, exits successfully" trap as the bare `stop` and the
# debugging exit(0) that this branch has hit twice.
if grep -q 'Unable to find injection library' sanitizer.log; then
  cat sanitizer.log
  die "compute-sanitizer could not inject; use the CUDA toolkit copy, not a distro stub"
fi

# Take the verdict from compute-sanitizer's own ERROR SUMMARY, not from a
# hand-written grep over the record headers. The first version of this script
# matched '^====== ' where the records actually begin with nine '=', so it
# counted zero errors on a run that reported one and printed "clean" -- the
# exact "reports success while being wrong" shape catalogued in
# session_handoff_2026-08.md section 1. If the summary line is absent at all,
# that is a failure too: it means the run did not complete.
summary=$(grep -E '^=+ ERROR SUMMARY' sanitizer.log | tail -1 || true)
if [ -z "$summary" ]; then
  tail -30 sanitizer.log
  die "no ERROR SUMMARY line -- the sanitizer run did not complete"
fi
total=$(printf '%s' "$summary" | sed -E 's/.*ERROR SUMMARY: ([0-9]+).*/\1/')

# ERROR SUMMARY counts two different things and only one of them is ours.
#
#   real     Invalid / Uninitialized / Race / Barrier -- a fault in code we
#            wrote, or in a kernel we launched with bad arguments.
#   benign   "CUDA API Error: Kernel (...cublasGemv...) cannot be loaded",
#            raised from cuLibraryLoadData inside libcublasLt. That is the
#            sanitizer failing to instrument a closed-source cuBLAS kernel,
#            not a memory error, and it fires on every run that touches
#            cublasDgemv -- which is every SOAP prediction.
#
# Counting them together means a genuine fault is one line among the noise and
# the run "has always had one error anyway". Separate them, fail only on the
# first, and still PRINT the second rather than filtering it away.
real=$(grep -cE '^=+ (Invalid|Uninitialized|Race|Barrier) ' sanitizer.log || true)
benign=$(grep -cE '^=+ CUDA API Error: Kernel .* cannot b' sanitizer.log || true)

grep -E '^=+ (Invalid|Uninitialized|Race|Barrier|Program hit|CUDA API Error)' sanitizer.log \
  | cut -c1-160 | head -40 || true
printf -- '\n--- %s\n' "$summary"
printf -- '--- device memory faults: %s   cuBLAS instrumentation notices: %s   process exit: %s\n' \
  "$real" "$benign" "$rc"
[ "$benign" -gt 0 ] && [ "$real" = 0 ] && printf -- '--- the cuBLAS notices are expected; see the comment in this script\n'
[ "$KEEP" = 1 ] && printf -- '--- full log kept at %s/sanitizer.log\n' "$WORK"

if [ "$real" -gt 0 ]; then
  printf 'FAIL: %s device memory fault(s)\n' "$real"
  exit 1
fi
if [ "$total" -gt "$((real + benign))" ]; then
  printf 'FAIL: ERROR SUMMARY reports %s but only %s were classified -- read the log\n' \
    "$total" "$((real + benign))"
  exit 1
fi
if [ "$rc" -ne 0 ]; then
  printf 'NOTE: no sanitizer errors, but the process exited %s -- a host-side\n' "$rc"
  printf '      fault. Rebuild with DEBUG=1 (-fcheck=bounds) and rerun; gdb\n'
  printf '      cannot attach on alt (ptrace_scope = 2).\n'
  exit 1
fi
printf 'clean\n'
