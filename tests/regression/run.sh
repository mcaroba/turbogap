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
#   TURBOGAP_DATA_ROOT  directory holding the test systems
#                       (default: $HOME/work/cpu_vs_gpu_tests/input)
#   TURBOGAP_TIME_TOL   fail if test/ref wall-clock exceeds this ratio
#                       (default: unset -> timing is reported, not enforced)
#   TURBOGAP_KEEP       keep staging directories for inspection
#
# Usage: run.sh [case ...]      (no arguments runs every case)
#        run.sh --list

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
REF_BIN=${TURBOGAP_REF_BIN:-$here/baseline/turbogap.e6eb1aa}
DATA_ROOT=${TURBOGAP_DATA_ROOT:-$HOME/work/cpu_vs_gpu_tests/input}
TIME_TOL=${TURBOGAP_TIME_TOL:-}
WORK=${TMPDIR:-/tmp}/turbogap_regression.$$

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

if [ "${1:-}" = --list ]; then
  for d in "$here"/cases/*/; do basename "$d"; done
  exit 0
fi

[ -x "$BIN" ]     || die "binary under test not found or not executable: $BIN"
[ -x "$REF_BIN" ] || die "reference binary not found or not executable: $REF_BIN"
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
  printf '%s' "$dir"
}

# run <dir> <binary> <label> -> sets REPLY to elapsed seconds
run() {
  local dir=$1 bin=$2 label=$3 rc t0 t1
  t0=$(date +%s.%N)
  if [ "$RANKS" -gt 1 ]; then
    ( cd "$dir" && mpirun -np "$RANKS" "$bin" "$MODE" ) > "$dir/run.log" 2>&1
  else
    ( cd "$dir" && "$bin" "$MODE" ) > "$dir/run.log" 2>&1
  fi
  rc=$?
  t1=$(date +%s.%N)
  REPLY=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
  if [ $rc -ne 0 ]; then
    printf '    %s run failed (exit %d):\n' "$label" "$rc"
    tail -15 "$dir/run.log" | sed 's/^/      /'
    return 1
  fi
  return 0
}

pass=0; fail=0; skip=0
printf 'binary   %s\n' "$BIN"
printf 'baseline %s\n' "$REF_BIN"
printf 'data     %s\n\n' "$DATA_ROOT"

for name in "${cases[@]}"; do
  conf=$here/cases/$name/case.conf
  [ -f "$conf" ] || die "no such case: $name"

  DATA=; MODE=; RANKS=1; LINKS=; OUTPUTS=; NEEDS=
  # shellcheck source=/dev/null
  . "$conf"

  printf '== %s ==\n' "$name"

  missing=
  for n in $NEEDS; do
    [ -e "$DATA_ROOT/$DATA/$n" ] || missing="$missing $n"
  done
  if [ -n "$missing" ]; then
    printf '    SKIP: missing test data:%s\n' "$missing"
    skip=$((skip+1)); continue
  fi
  if [ "$RANKS" -gt 1 ] && ! command -v mpirun >/dev/null; then
    printf '    SKIP: mpirun not available (case needs %d ranks)\n' "$RANKS"
    skip=$((skip+1)); continue
  fi

  rdir=$(stage "$name" ref)
  tdir=$(stage "$name" test)

  if ! run "$rdir" "$REF_BIN" reference; then fail=$((fail+1)); continue; fi
  rtime=$REPLY
  if ! run "$tdir" "$BIN" test; then fail=$((fail+1)); continue; fi
  ttime=$REPLY

  bad=
  for out in $OUTPUTS; do
    if [ ! -f "$rdir/$out" ] && [ ! -f "$tdir/$out" ]; then
      bad="$bad $out(absent-from-both)"; continue
    fi
    if [ ! -f "$rdir/$out" ] || [ ! -f "$tdir/$out" ]; then
      bad="$bad $out(produced-by-only-one)"; continue
    fi
    if ! diff -q "$rdir/$out" "$tdir/$out" >/dev/null; then
      bad="$bad $out"
    fi
  done

  ratio=$(awk -v t="$ttime" -v r="$rtime" 'BEGIN{ if (r>0) printf "%.2f", t/r; else print "n/a" }')
  if [ -z "$bad" ]; then
    printf '    PASS   ref %ss  test %ss  (x%s)\n' "$rtime" "$ttime" "$ratio"
    pass=$((pass+1))
  else
    printf '    FAIL   differing outputs:%s\n' "$bad"
    for out in $bad; do
      case $out in *\(*\)) continue;; esac
      printf '      --- first differences in %s ---\n' "$out"
      diff "$rdir/$out" "$tdir/$out" | head -12 | sed 's/^/      /'
    done
    fail=$((fail+1))
    continue
  fi

  if [ -n "$TIME_TOL" ]; then
    if awk -v x="$ratio" -v tol="$TIME_TOL" 'BEGIN{exit !(x>tol)}'; then
      printf '    SLOW   %sx slower than baseline (tolerance %sx)\n' "$ratio" "$TIME_TOL"
      fail=$((fail+1)); pass=$((pass-1))
    fi
  fi
done

printf '\npassed: %d   failed: %d   skipped: %d\n' "$pass" "$fail" "$skip"
[ "$fail" -eq 0 ]
