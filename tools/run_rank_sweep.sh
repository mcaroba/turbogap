#!/usr/bin/env bash
#
# Run one system at several MPI rank counts on the SAME GPU, and tabulate what
# it buys.
#
#     tools/run_rank_sweep.sh --atoms atoms_124959.xyz
#     tools/run_rank_sweep.sh --ranks '1 2 4 8' --input input.md --mode md
#     tools/run_rank_sweep.sh --mps --ranks '1 4'
#
# The question this exists for
# ----------------------------
# "Is each GPU tied to one rank?" No. src/gpu/gpu_memory.cu does
#
#     hipSetDevice(my_rank % num_gpus)
#
# so on a one-GPU node EVERY rank binds to device 0. One rank per GPU is a
# convention the submission scripts follow, not something the code enforces.
# Oversubscribing is therefore already supported, and this measures what it
# actually gives you.
#
# What to expect, and why the --mps row matters
# ---------------------------------------------
# TurboGAP's MPI decomposition splits SITES (turbogap.f90: i_beg/i_end per
# rank), and build_neighbors_list is called with a do_list mask, so each rank
# builds only its own atoms' neighbours. All the host-side work therefore
# parallelises across ranks for free.
#
# The device side does not, by default. Each rank gets its own CUDA context, and
# the driver TIME-SLICES between contexts: kernels from different ranks do not
# run at the same time, they take turns, with a context switch between. So
# without MPS the gain comes from the host half of the run -- which the
# profiling in docs/PROFILING.md measured at 38-50% of the SOAP phase and 100%
# of the neighbour build, so it is not a small half.
#
# With MPS (--mps) the ranks share one context and their kernels really do run
# concurrently on the SMs, the same way streams within one process would. That
# is the row that separates "the host work parallelised" from "the GPU is now
# being shared properly".
#
# Memory is the constraint on how far this goes: N contexts on one card, each
# with its own buffers. Each rank holds ~1/N of the atoms so the buffers shrink,
# but the ~300 MB per-context overhead does not. gpu_mem_fraction divides the
# budget by ranks-per-device for exactly this reason.

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

DIR=${TURBOGAP_LARGE_DIR:-$repo/../large_systems/diamond_1M}
ATOMS=atoms_124959.xyz
INPUT=input.gap
MODE=predict
BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
WORK=${TMPDIR:-/tmp}/turbogap_rank_sweep
RANKS="1 2 4 6"
USE_MPS=0

log() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case $1 in
  --dir) DIR=$2; shift 2 ;;
  --atoms) ATOMS=$2; shift 2 ;;
  --input) INPUT=$2; shift 2 ;;
  --mode) MODE=$2; shift 2 ;;
  --bin) BIN=$2; shift 2 ;;
  --ranks) RANKS=$2; shift 2 ;;
  --work) WORK=$2; shift 2 ;;
  --mps) USE_MPS=1; shift ;;
  -h | --help) sed -n '2,45p' "$0" | sed 's/^#\{0,1\} \{0,1\}//'; exit 0 ;;
  *) die "unknown option: $1" ;;
  esac
done

[ -x "$BIN" ] || die "no binary at $BIN"
[ -f "$DIR/$ATOMS" ] || die "no $ATOMS in $DIR"
[ -f "$DIR/$INPUT" ] || die "no $INPUT in $DIR"
command -v mpirun >/dev/null || die "mpirun not found"

mkdir -p "$WORK"

mps_started=0
if [ "$USE_MPS" = 1 ]; then
  command -v nvidia-cuda-mps-control >/dev/null ||
    die "--mps asked for but nvidia-cuda-mps-control is not installed"
  export CUDA_MPS_PIPE_DIRECTORY=$WORK/mps-pipe
  export CUDA_MPS_LOG_DIRECTORY=$WORK/mps-log
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  if echo get_server_list | nvidia-cuda-mps-control 2>/dev/null | grep -q '[0-9]'; then
    log "MPS: already running, leaving it alone"
  else
    log "MPS: starting a per-user daemon (pipes in $CUDA_MPS_PIPE_DIRECTORY)"
    nvidia-cuda-mps-control -d || die "could not start the MPS daemon"
    mps_started=1
  fi
  # Stopped on the way out, including on an interrupt: a stray MPS daemon
  # silently changes how every later CUDA process on this machine behaves,
  # which is a horrible thing to leave behind.
  trap 'if [ "$mps_started" = 1 ]; then
          echo "MPS: stopping the daemon"
          echo quit | nvidia-cuda-mps-control >/dev/null 2>&1
        fi' EXIT INT TERM
fi

log ""
log "$ATOMS, $INPUT, mode $MODE, MPS $([ "$USE_MPS" = 1 ] && echo on || echo off)"
log ""
printf '%7s  %9s  %9s  %10s  %12s  %s\n' ranks wall_s speedup host_GB dev_peak_GB result
printf '%7s  %9s  %9s  %10s  %12s  %s\n' ------- --------- --------- ---------- ------------ ------

base=

for n in $RANKS; do
  d=$WORK/r$n
  rm -rf "$d"
  mkdir -p "$d"
  ln -sf "$DIR/$ATOMS" "$d/atoms.xyz"
  ln -sf "$DIR/gap_files" "$d/gap_files"
  cp "$DIR/$INPUT" "$d/input"

  # --oversubscribe: more ranks than the machine claims cores for is normal
  # here; the point is to share one GPU, not to fill the node.
  ( cd "$d" && /usr/bin/time -f '%e %M' -o "$d/time.txt" \
      mpirun --oversubscribe -np "$n" "$BIN" "$MODE" >"$d/run.log" 2>&1 )
  rc=$?

  # /usr/bin/time prepends "Command exited with non-zero status N" on a failure,
  # so the numbers are on the LAST line, not the first. Reading $1 of line 1 on a
  # failed run picked up the word "Command" and fed it to the speedup division.
  wall=$(tail -1 "$d/time.txt" 2>/dev/null | awk '{print $1+0}')
  rss_kb=$(tail -1 "$d/time.txt" 2>/dev/null | awk '{print $2+0}')
  host_gb=$(awk -v k="${rss_kb:-0}" 'BEGIN{printf "%.2f", k/1048576}')
  # Rank 0's ledger only. Every rank has its own; the sum is what the card
  # sees, and rank 0's is the per-rank figure the budget is set from.
  dev_peak=$(grep -m1 'this rank' "$d/run.log" 2>/dev/null |
    sed -n 's/.*peak *\([0-9.]*\) GB.*/\1/p')

  if [ $rc -eq 0 ]; then result=ok
  elif grep -q 'GPUmem: device allocation' "$d/run.log" 2>/dev/null; then result="DEVICE OOM"
  elif grep -q 'illegal memory access' "$d/run.log" 2>/dev/null; then
    # Known and pre-existing at e70d896: the GPU branch takes an illegal device
    # address with more than one rank, on inputs that are fine with one. Named
    # here rather than left as "exit 188" so a sweep that hits it says so.
    result="ILLEGAL ADDRESS (multi-rank bug)"
  else result="failed (exit $rc)"; fi

  # Only a SUCCESSFUL run sets the baseline. A crash that happens to be quick
  # would otherwise become the reference every later row is measured against.
  [ -z "$base" ] && [ "$rc" -eq 0 ] && [ -n "$wall" ] && base=$wall
  speedup=$(awk -v b="${base:-0}" -v w="${wall:-0}" -v ok="$rc" \
    'BEGIN{ if (ok==0 && w>0 && b>0) printf "%.2fx", b/w; else printf "-" }')

  printf '%7s  %9s  %9s  %10s  %12s  %s\n' \
    "$n" "${wall:-?}" "$speedup" "$host_gb" "${dev_peak:-?}" "$result"
done

log ""
log "Working directories under $WORK"
