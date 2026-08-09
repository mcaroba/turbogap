#!/usr/bin/env bash
#
# GPU regression tests: run TurboGAP on a small set of cases and check the
# results against a reference produced by the CPU code.
#
# The reference can come from either
#   * a CPU TurboGAP binary (TURBOGAP_REF_BIN), which is run here, or
#   * a stored reference trajectory (tests/gpu/reference/<case>.xyz),
# so the suite works both next to a CPU build and on a machine that only has
# the GPU build.
#
# Environment:
#   TURBOGAP_BIN       GPU binary under test      (default: <repo>/bin/turbogap)
#   TURBOGAP_REF_BIN   CPU binary for reference   (default: unset -> stored reference)
#   TURBOGAP_TEST_DATA directory holding the CO test inputs
#                      (default: $HOME/work/cpu_vs_gpu_tests/input/CO)
#   TURBOGAP_MPI_RANKS if set, run the GPU binary under mpirun with this many ranks
#
# With no arguments every case in tests/gpu/cases.sh runs; name cases on the
# command line to run a subset. That file is the single definition of what the
# cases ARE -- tools/profile_gpu.sh sources it too, so a profile and a test run
# are of the same input.
#
# Exit status is non-zero if any case fails.

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
REF_BIN=${TURBOGAP_REF_BIN:-}
DATA=${TURBOGAP_TEST_DATA:-$repo/../turbogap_tests/CO}
COMPARE=$repo/tools/compare_xyz.py
REFDIR=$here/reference
WORK=${TMPDIR:-/tmp}/turbogap_gpu_regression.$$

fail=0
pass=0
xfail=0
xpass=0

log() { printf '%s\n' "$*"; }
die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[ -x "$BIN" ] || die "GPU binary not found or not executable: $BIN"
[ -f "$COMPARE" ] || die "comparison script not found: $COMPARE"
# The test systems live in their own repository. $DATA points at one system
# inside it, so fetch the whole clone and then re-check.
if [ ! -d "$DATA" ]; then
  "$repo/tests/fetch_test_data.sh" "${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}" \
    || die "could not fetch test data"
fi
[ -d "$DATA" ] || die "test data directory not found: $DATA (set TURBOGAP_TEST_DATA)"
command -v python3 >/dev/null || die "python3 is required"

mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# The case table. tools/profile_gpu.sh sources the same file, so a profile is
# always of the input this suite checks.
TG_TESTS_DIR=${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}
TG_DATA=$DATA
. "$here/cases.sh"

# run_case <dir> <binary> <mode> <label>
run_case() {
  local dir=$1 bin=$2 mode=$3 label=$4
  # </dev/null is load-bearing. mpirun reads its own stdin, and the case loop
  # at the bottom of this file feeds names to `while read` on the shell's
  # stdin -- so under TURBOGAP_MPI_RANKS the first mpirun swallowed the rest of
  # the case list and the suite reported "passed: 1" as though it had run
  # everything and been happy.
  (cd "$dir" &&
    if [ -n "${TURBOGAP_MPI_RANKS:-}" ] && [ "$label" = gpu ]; then
      mpirun --oversubscribe -np "$TURBOGAP_MPI_RANKS" "$bin" "$mode"
    else
      "$bin" "$mode"
    fi) >"$dir/$label.log" 2>&1 </dev/null
  local rc=$?
  if [ $rc -ne 0 ]; then
    log "    $label run failed (exit $rc); tail of log:"
    tail -15 "$dir/$label.log" | sed 's/^/      /'
    return 1
  fi
  if grep -q GPUassert "$dir/$label.log"; then
    log "    $label run reported a GPU error:"
    grep -m3 GPUassert "$dir/$label.log" | sed 's/^/      /'
    return 1
  fi
  return 0
}

check_case() {
  local name=$1
  log "== $name =="
  tg_case_def "$name" || die "no such case: $name"
  local data=$TG_DATA_DIR
  if [ ! -d "$data" ]; then
    log "    SKIP: data directory not present: $data"
    return
  fi
  if [ ! -e "$data/$TG_ATOMS" ]; then
    log "    SKIP: $data/$TG_ATOMS not present"
    [ "$name" = CO_md ] &&
      log "          (create it with tools/add_velocities_to_xyz.py)"
    return
  fi
  local mode=$TG_MODE
  local extras=(${TG_EXTRAS+"${TG_EXTRAS[@]}"})
  local dir
  dir=$(tg_stage_case "$name" "$WORK")

  local xreason
  xreason=$(tg_case_xfail "$name")
  if ! run_case "$dir" "$BIN" "$mode" gpu; then
    if [ -n "$xreason" ]; then
      log "    XFAIL: $xreason"
      xfail=$((xfail + 1))
    else
      fail=$((fail + 1))
    fi
    return
  fi
  if [ -n "$xreason" ] && [ ! -f "$dir/trajectory_out.xyz" ]; then
    log "    XFAIL: $xreason"
    log "           (ran to exit 0 but produced no trajectory)"
    xfail=$((xfail + 1))
    return
  fi

  local ref=$REFDIR/$name.xyz
  if [ -n "$REF_BIN" ]; then
    local rdir=$WORK/$name.ref
    mkdir -p "$rdir"
    ln -sf "$data/$TG_ATOMS" "$rdir/atoms.xyz"
    ln -sf "$data/gap_files" "$rdir/gap_files"
    local extra
    for extra in ${extras+"${extras[@]}"}; do ln -sf "$data/$extra" "$rdir/$extra"; done
    cp "$dir/input" "$rdir/input"
    if ! run_case "$rdir" "$REF_BIN" "$mode" cpu; then
      log "    reference run failed"
      fail=$((fail + 1))
      return
    fi
    ref=$rdir/trajectory_out.xyz
  elif [ ! -f "$ref" ]; then
    log "    SKIP: no reference (set TURBOGAP_REF_BIN or add $ref)"
    return
  fi

  if python3 "$COMPARE" "$ref" "$dir/trajectory_out.xyz" >"$dir/compare.txt" 2>&1; then
    if [ -n "$xreason" ]; then
      log "    XPASS: expected to fail but did not -- remove it from xfail_reason"
      xpass=$((xpass + 1))
    else
      log "    PASS"
      pass=$((pass + 1))
    fi
  elif [ -n "$xreason" ]; then
    log "    XFAIL: $xreason"
    xfail=$((xfail + 1))
  else
    log "    FAIL"
    sed 's/^/      /' "$dir/compare.txt"
    fail=$((fail + 1))
  fi
}

# Run the whole table, or the cases named on the command line.
if [ $# -gt 0 ]; then
  for name in "$@"; do check_case "$name"; done
else
  while read -r name; do check_case "$name"; done < <(tg_case_list)
fi

log ""
log "passed: $pass   failed: $fail   xfail: $xfail   xpass: $xpass"
if [ "$xpass" -gt 0 ]; then
  log ""
  log "An expected failure passed. Remove it from tg_case_xfail() in cases.sh"
  log "and update"
  log "docs/gpu_fixes_handoff.md before this masks a real regression."
fi
[ "$fail" -eq 0 ]
