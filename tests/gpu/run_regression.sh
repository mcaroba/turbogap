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
# Exit status is non-zero if any case fails.

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
REF_BIN=${TURBOGAP_REF_BIN:-}
DATA=${TURBOGAP_TEST_DATA:-$HOME/work/cpu_vs_gpu_tests/input/CO}
COMPARE=$repo/tools/compare_xyz.py
REFDIR=$here/reference
WORK=${TMPDIR:-/tmp}/turbogap_gpu_regression.$$

fail=0
pass=0

log()  { printf '%s\n' "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

[ -x "$BIN" ]      || die "GPU binary not found or not executable: $BIN"
[ -f "$COMPARE" ]  || die "comparison script not found: $COMPARE"
[ -d "$DATA" ]     || die "test data directory not found: $DATA (set TURBOGAP_TEST_DATA)"
command -v python3 >/dev/null || die "python3 is required"

mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# stage_case <case> <atoms-file> <input-body>
stage_case() {
  local name=$1 atoms=$2 body=$3 dir=$WORK/$name
  mkdir -p "$dir"
  ln -sf "$DATA/$atoms"   "$dir/atoms.xyz"
  ln -sf "$DATA/gap_files" "$dir/gap_files"
  printf '%s\n' "$body" > "$dir/input"
  printf '%s' "$dir"
}

# run_case <dir> <binary> <mode> <label>
run_case() {
  local dir=$1 bin=$2 mode=$3 label=$4
  ( cd "$dir" && \
    if [ -n "${TURBOGAP_MPI_RANKS:-}" ] && [ "$label" = gpu ]; then
      mpirun -np "$TURBOGAP_MPI_RANKS" "$bin" "$mode"
    else
      "$bin" "$mode"
    fi ) > "$dir/$label.log" 2>&1
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
  local name=$1 atoms=$2 mode=$3 body=$4
  log "== $name =="
  local dir; dir=$(stage_case "$name" "$atoms" "$body")

  if ! run_case "$dir" "$BIN" "$mode" gpu; then
    fail=$((fail+1)); return
  fi

  local ref=$REFDIR/$name.xyz
  if [ -n "$REF_BIN" ]; then
    local rdir=$WORK/$name.ref
    mkdir -p "$rdir"
    ln -sf "$DATA/$atoms" "$rdir/atoms.xyz"
    ln -sf "$DATA/gap_files" "$rdir/gap_files"
    cp "$dir/input" "$rdir/input"
    if ! run_case "$rdir" "$REF_BIN" "$mode" cpu; then
      log "    reference run failed"; fail=$((fail+1)); return
    fi
    ref=$rdir/trajectory_out.xyz
  elif [ ! -f "$ref" ]; then
    log "    SKIP: no reference (set TURBOGAP_REF_BIN or add $ref)"
    return
  fi

  if python3 "$COMPARE" "$ref" "$dir/trajectory_out.xyz" > "$dir/compare.txt" 2>&1; then
    log "    PASS"
    pass=$((pass+1))
  else
    log "    FAIL"
    sed 's/^/      /' "$dir/compare.txt"
    fail=$((fail+1))
  fi
}

common='atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345'

# Single point: exercises the SOAP, 2b and 3b energy, force and virial paths.
check_case CO_predict atoms_7176.xyz predict "$common
write_xyz = 1"

# Short MD from explicit velocities, so no randomization is involved and the
# trajectory is reproducible between the CPU and GPU builds.
if [ -f "$DATA/atoms_7176_vel.xyz" ]; then
  check_case CO_md atoms_7176_vel.xyz md "$common
md_nsteps = 5
thermostat = berendsen
write_xyz = 1"
else
  log "== CO_md =="
  log "    SKIP: $DATA/atoms_7176_vel.xyz not present"
  log "          (create it with tools/add_velocities_to_xyz.py)"
fi

log ""
log "passed: $pass   failed: $fail"
[ "$fail" -eq 0 ]
