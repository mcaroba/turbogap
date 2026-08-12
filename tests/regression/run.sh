#!/usr/bin/env bash
#
# Refactor regression suite.
#
# The contract for the turbogap.f90 modularization is that output is
# BIT-IDENTICAL to the pre-refactor baseline. Extracting a block into a
# procedure must not reorder floating-point operations, so there is no reason
# for a single digit to move. That makes plain diff the right comparator: it is
# exact, needs no tolerance tuning, and cannot be quietly loosened. A change
# that legitimately alters results is not a refactor and does not belong on
# this branch.
#
# Each case runs twice: once with the reference binary, once with the binary
# under test, in separate staging directories. Every file listed in the case's
# OUTPUTS is diffed. Wall-clock is recorded for both so that a performance
# regression shows up on the commit that caused it rather than at the end.
#
# Environment:
#   TURBOGAP_BIN        binary under test  (default: <repo>/bin/turbogap)
#   TURBOGAP_REF_BIN    reference binary   (default: <repo>/tests/regression/baseline/turbogap.e6eb1aa)
#   TURBOGAP_DATA_ROOT  directory holding the test systems. Default: the
#                       turbogap_tests clone beside this repository, which
#                       tests/fetch_test_data.sh creates on the first run
#                       (override where it goes with TURBOGAP_TESTS_DIR).
#   TURBOGAP_TIME_TOL   fail if test/ref wall-clock exceeds this ratio
#                       (default: unset -> timing is reported, not enforced)
#   TURBOGAP_CASE_TIMEOUT  seconds a single run may take before it is killed
#                       and the case failed (default: 1800)
#   TURBOGAP_KEEP       keep staging directories for inspection
#   TURBOGAP_BLESS      regenerate the expected/ output of golden cases
#
# Most cases compare the binary under test against the frozen baseline. A case
# whose case.conf sets REFERENCE=golden instead compares against a checked-in
# expected/ directory, for the situation where the baseline cannot produce a
# reference at all because it crashes on that input. Those are characterization
# tests: they pin behaviour against future drift and assert nothing about
# whether the physics is right.
#
# Usage: run.sh [case ...]      (no arguments runs every case)
#        run.sh --list

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
REF_BIN=${TURBOGAP_REF_BIN:-$here/baseline/turbogap.e6eb1aa}
TIME_TOL=${TURBOGAP_TIME_TOL:-}
CASE_TIMEOUT=${TURBOGAP_CASE_TIMEOUT:-1800}
WORK=${TMPDIR:-/tmp}/turbogap_regression.$$

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

if [ "${1:-}" = --list ]; then
  for d in "$here"/cases/*/; do basename "$d"; done
  exit 0
fi

[ -x "$BIN" ] || die "binary under test not found or not executable: $BIN"
[ -x "$REF_BIN" ] || die "reference binary not found or not executable: $REF_BIN"
# The test systems live in their own repository, fetched on first use rather
# than making every new checkout a two-step setup.
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
[ -d "$DATA_ROOT" ] || die "data root not found: $DATA_ROOT (set TURBOGAP_DATA_ROOT)"

if [ "$BIN" -ef "$REF_BIN" ]; then
  die "binary under test and reference are the same file; nothing would be tested"
fi

mkdir -p "$WORK"
if [ -z "${TURBOGAP_KEEP:-}" ]; then
  trap 'rm -rf "$WORK"' EXIT
else
  trap 'printf "staging kept in %s\n" "$WORK"' EXIT
fi

if [ $# -gt 0 ]; then
  cases=("$@")
else
  cases=()
  for d in "$here"/cases/*/; do cases+=("$(basename "$d")"); done
fi

# stage <case> <label> -> prints the staging directory
stage() {
  local name=$1
  local label=$2
  local dir=$WORK/$name.$label
  mkdir -p "$dir"
  local spec src dst
  for spec in $LINKS; do
    src=${spec%%:*}
    dst=${spec#*:}
    [ "$dst" = "$spec" ] && dst=$src
    [ -e "$DATA_ROOT/$DATA/$src" ] || die "$name: missing input $DATA_ROOT/$DATA/$src"
    ln -sf "$DATA_ROOT/$DATA/$src" "$dir/$dst"
  done
  cp "$here/cases/$name/input" "$dir/input"
  # Pin the SOAP batch split. gpu_memory_budget_init sizes
  # max_Gbytes_per_process from the node's memory, and the batch count it
  # produces changes the order the per-batch force contributions are summed --
  # so without this the stored references would only reproduce on a machine
  # with the same amount of RAM. 1.0 is the historical default the baselines
  # were recorded with. Pinned here rather than in each cases/*/input because
  # it is a property of comparing against stored text, not of any one case.
  printf '\nmax_Gbytes_per_process = 1.0\n' >>"$dir/input"
  printf '%s' "$dir"
}

# run <dir> <binary> <label> -> sets REPLY to elapsed seconds
#
# Wrapped in a timeout because a defect that hangs the run is a real failure
# mode here, not a hypothetical: KNOWN_ISSUES #6 made the input-file loop
# re-read the same record for ever, and without this the whole suite waits
# behind it with nothing on the terminal to say which case is stuck. A timed-out
# case reports exit 124 and the suite moves on.
run() {
  local dir=$1 bin=$2 label=$3 rc t0 t1
  t0=$(date +%s.%N)
  if [ "$RANKS" -gt 1 ]; then
    (cd "$dir" && timeout "$CASE_TIMEOUT" mpirun -np "$RANKS" "$bin" "$MODE") >"$dir/run.log" 2>&1
  else
    (cd "$dir" && timeout "$CASE_TIMEOUT" "$bin" "$MODE") >"$dir/run.log" 2>&1
  fi
  rc=$?
  t1=$(date +%s.%N)
  REPLY=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
  if [ $rc -eq 124 ]; then
    printf '    %s run timed out after %ss (TURBOGAP_CASE_TIMEOUT)\n' "$label" "$CASE_TIMEOUT"
    tail -15 "$dir/run.log" | sed 's/^/      /'
    return 1
  fi
  if [ $rc -ne 0 ]; then
    printf '    %s run failed (exit %d):\n' "$label" "$rc"
    tail -15 "$dir/run.log" | sed 's/^/      /'
    return 1
  fi
  return 0
}

pass=0
fail=0
skip=0
printf 'binary   %s\n' "$BIN"
printf 'baseline %s\n' "$REF_BIN"
printf 'data     %s\n\n' "$DATA_ROOT"

for name in "${cases[@]}"; do
  conf=$here/cases/$name/case.conf
  [ -f "$conf" ] || die "no such case: $name"

  DATA=
  MODE=
  RANKS=1
  LINKS=
  OUTPUTS=
  NEEDS=
  REFERENCE=baseline
  # shellcheck source=/dev/null
  . "$conf"

  printf '== %s ==\n' "$name"

  missing=
  for n in $NEEDS; do
    [ -e "$DATA_ROOT/$DATA/$n" ] || missing="$missing $n"
  done
  if [ -n "$missing" ]; then
    printf '    SKIP: missing test data:%s\n' "$missing"
    skip=$((skip + 1))
    continue
  fi
  if [ "$RANKS" -gt 1 ] && ! command -v mpirun >/dev/null; then
    printf '    SKIP: mpirun not available (case needs %d ranks)\n' "$RANKS"
    skip=$((skip + 1))
    continue
  fi

  tdir=$(stage "$name" test)

  if [ "$REFERENCE" = golden ]; then
    # The baseline binary cannot produce a reference for this case -- it
    # crashes on it. Compare against a checked-in expected output instead.
    # This is a characterization test: it pins behaviour so later refactors
    # cannot drift it, and says nothing about whether the physics is right.
    rdir=$here/cases/$name/expected
    if [ ! -d "$rdir" ] && [ -z "${TURBOGAP_BLESS:-}" ]; then
      printf '    SKIP: no expected/ directory (regenerate with TURBOGAP_BLESS=1)\n'
      skip=$((skip + 1))
      continue
    fi
    if ! run "$tdir" "$BIN" test; then
      fail=$((fail + 1))
      continue
    fi
    ttime=$REPLY
    rtime=$ttime
    if [ -n "${TURBOGAP_BLESS:-}" ]; then
      mkdir -p "$rdir"
      for out in $OUTPUTS; do cp "$tdir/$out" "$rdir/$out"; done
      printf '    BLESSED expected output (%ss)\n' "$ttime"
      pass=$((pass + 1))
      continue
    fi
  else
    rdir=$(stage "$name" ref)
    if ! run "$rdir" "$REF_BIN" reference; then
      fail=$((fail + 1))
      continue
    fi
    rtime=$REPLY
    if ! run "$tdir" "$BIN" test; then
      fail=$((fail + 1))
      continue
    fi
    ttime=$REPLY
  fi

  bad=
  for out in $OUTPUTS; do
    if [ ! -f "$rdir/$out" ] && [ ! -f "$tdir/$out" ]; then
      bad="$bad $out(absent-from-both)"
      continue
    fi
    if [ ! -f "$rdir/$out" ] || [ ! -f "$tdir/$out" ]; then
      bad="$bad $out(produced-by-only-one)"
      continue
    fi
    if ! diff -q "$rdir/$out" "$tdir/$out" >/dev/null; then
      bad="$bad $out"
    fi
  done

  ratio=$(awk -v t="$ttime" -v r="$rtime" 'BEGIN{ if (r>0) printf "%.2f", t/r; else print "n/a" }')
  if [ -z "$bad" ]; then
    printf '    PASS   ref %ss  test %ss  (x%s)\n' "$rtime" "$ttime" "$ratio"
    pass=$((pass + 1))
  else
    printf '    FAIL   differing outputs:%s\n' "$bad"
    for out in $bad; do
      case $out in *\(*\)) continue ;; esac
      printf '      --- first differences in %s ---\n' "$out"
      diff "$rdir/$out" "$tdir/$out" | head -12 | sed 's/^/      /'
    done
    fail=$((fail + 1))
    continue
  fi

  if [ -n "$TIME_TOL" ]; then
    if awk -v x="$ratio" -v tol="$TIME_TOL" 'BEGIN{exit !(x>tol)}'; then
      printf '    SLOW   %sx slower than baseline (tolerance %sx)\n' "$ratio" "$TIME_TOL"
      fail=$((fail + 1))
      pass=$((pass - 1))
    fi
  fi
done

printf '\npassed: %d   failed: %d   skipped: %d\n' "$pass" "$fail" "$skip"
[ "$fail" -eq 0 ]
