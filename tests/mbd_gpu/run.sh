#!/usr/bin/env bash
#
# Does MBD on the device (vdw_mbd_device, src/gpu/gap_mbd.cu) give what MBD on
# the host threads gives, in the same GPU-build binary?
#
# Two systems with the vdw_mbd_cell_mpi2 settings: the dilute 108-atom P4
# cell, a few atoms in each MBD cluster, and a dense black-phosphorus cell,
# hundreds. Each runs with vdw_mbd_device off and on, on several threads (each
# thread hands its own sites to the device on its own stream), and the outputs
# must agree to the precision they are written with. The device sums in a
# different order, so this is a tolerance, not a byte comparison.
#
# Needs a GPU build and a device; skips loudly without either.
#
# Environment:
#   TURBOGAP_GPU_BIN    device binary (default: <repo>/bin-gpu/turbogap)
#   TURBOGAP_DATA_ROOT  directory holding the test systems (see ../data_root.sh)
#   TURBOGAP_PYTHON     interpreter for the cell and the comparison
#   TURBOGAP_KEEP       keep the staging directory for inspection

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

BIN=${TURBOGAP_GPU_BIN:-$repo/bin-gpu/turbogap}
PYTHON=${TURBOGAP_PYTHON:-python3}
# shellcheck source=../data_root.sh
. "$here/../data_root.sh"
DATA=$DATA_ROOT/vdw_P
WORK=${TMPDIR:-/tmp}/turbogap_mbd_gpu.$$
CASE_INPUT=$repo/tests/regression/cases/vdw_mbd_cell_mpi2/input
DENSE_CELL="2 1 1"
THREADS=4
# Energies and forces are written to 8 decimals.
RTOL=1e-7
ATOL=1e-7

skip() {
  printf '    SKIP: %s\n' "$1"
  exit 0
}

check_preconditions() {
  [ -x "$BIN" ] || skip "no GPU build at $BIN"
  command -v nvidia-smi > /dev/null || skip "no nvidia-smi on PATH"
  nvidia-smi -L > /dev/null 2>&1 || skip "no usable device"
  [ -e "$DATA/p4_cell.xyz" ] || { printf 'ERROR: missing %s\n' "$DATA/p4_cell.xyz" >&2; exit 2; }
}

run_one() {
  local dir=$1
  local atoms=$2
  local device=$3
  mkdir -p "$dir"
  ln -sfn "$DATA/gap_files" "$dir/gap_files"
  cp "$atoms" "$dir/cell.xyz"
  sed 's/p4_cell.xyz/cell.xyz/' "$CASE_INPUT" > "$dir/input"
  echo "vdw_mbd_device = $device" >> "$dir/input"
  (cd "$dir" && OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=1 "$BIN" predict > out.log 2>&1) || {
    printf '    FAIL: %s did not run; see %s/out.log\n' "$(basename "$dir")" "$dir" >&2
    return 1
  }
}

# The device run must say MBD ran on the device, and the host run must not:
# otherwise the comparison below could pass with both on the host.
check_where() {
  local name=$1
  local said='on the device'
  grep -q "$said" "$WORK/${name}_device/out.log" || {
    printf '    FAIL  %s: the device run does not report MBD on the device\n' "$name"
    return 1
  }
  if grep -q "$said" "$WORK/${name}_host/out.log"; then
    printf '    FAIL  %s: the host run reports MBD on the device\n' "$name"
    return 1
  fi
}

compare() {
  local name=$1
  local host=$WORK/${name}_host/trajectory_out.xyz
  local device=$WORK/${name}_device/trajectory_out.xyz
  if "$PYTHON" "$repo/tests/regression/compare_tol.py" "$host" "$device" --rtol $RTOL --atol $ATOL; then
    printf '    PASS  %s: device agrees with the host\n' "$name"
  else
    printf '    FAIL  %s: device and host differ beyond rtol %s atol %s\n' "$name" $RTOL $ATOL
    return 1
  fi
}

main() {
  printf '\n== MBD on the device against MBD on the host threads, same GPU build\n'
  check_preconditions
  mkdir -p "$WORK"
  "$PYTHON" "$here/black_p.py" $DENSE_CELL "$WORK/black_p.xyz"
  local fails=0
  for name in p4_cell black_p; do
    local atoms=$DATA/p4_cell.xyz
    [ "$name" = black_p ] && atoms=$WORK/black_p.xyz
    run_one "$WORK/${name}_host" "$atoms" .false. || fails=$((fails + 1))
    run_one "$WORK/${name}_device" "$atoms" .true. || fails=$((fails + 1))
    check_where "$name" || fails=$((fails + 1))
    compare "$name" || fails=$((fails + 1))
  done
  [ -n "${TURBOGAP_KEEP:-}" ] || rm -rf "$WORK"
  [ "$fails" -eq 0 ]
}

main
