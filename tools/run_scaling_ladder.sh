#!/usr/bin/env bash
#
# Walk a system-size ladder and record, for each size, whether it ran, how long
# it took, and what it cost in host and device memory.
#
#     tools/run_scaling_ladder.sh
#     tools/run_scaling_ladder.sh --max 262202
#     tools/run_scaling_ladder.sh --input input.gap --dir /path/to/systems
#
# This is not a test. There is no reference to compare a million-atom result
# against, and producing one is not the point -- the questions are "does it fit"
# and "how does it grow". It answers both in a table.
#
# Why it stops at the first failure by default: past the point where the device
# runs out, every larger size fails the same way for the same reason, and each
# one costs longer than the last to find that out. --keep-going overrides it.
#
# Peak host RSS comes from /usr/bin/time, peak device memory from the ledger
# line TurboGAP prints at the end of a run. Both matter and they fail
# differently: the host runs out during the neighbour build, the device during
# the descriptor loop.

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

DIR=${TURBOGAP_LARGE_DIR:-$repo/../large_systems/diamond_1M}
INPUT=input.gap
BIN=${TURBOGAP_BIN:-$repo/bin/turbogap}
MODE=predict
WORK=${TURBOGAP_LADDER_WORK:-${TMPDIR:-/tmp}/turbogap_ladder}
MAX=0
KEEP_GOING=0
RANKS=

log() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case $1 in
  --dir) DIR=$2; shift 2 ;;
  --input) INPUT=$2; shift 2 ;;
  --bin) BIN=$2; shift 2 ;;
  --mode) MODE=$2; shift 2 ;;
  --work) WORK=$2; shift 2 ;;
  --max) MAX=$2; shift 2 ;;
  --ranks) RANKS=$2; shift 2 ;;
  --keep-going) KEEP_GOING=1; shift ;;
  -h | --help) sed -n '2,25p' "$0" | sed 's/^#\{0,1\} \{0,1\}//'; exit 0 ;;
  *) die "unknown option: $1" ;;
  esac
done

[ -x "$BIN" ] || die "no binary at $BIN"
[ -d "$DIR" ] || die "no system directory at $DIR"
[ -f "$DIR/$INPUT" ] || die "no $INPUT in $DIR"

mkdir -p "$WORK"

# The ladder is whatever atoms_<N>.xyz files exist, smallest first. Deriving it
# from the directory rather than hardcoding sizes means tools/make_subcell.py
# can add a rung without this script knowing.
sizes=$(ls "$DIR"/atoms_*.xyz 2>/dev/null |
  sed 's/.*atoms_\([0-9]*\)\.xyz/\1/' | sort -n)
[ -n "$sizes" ] || die "no atoms_<N>.xyz files in $DIR"

printf '%12s  %9s  %10s  %12s  %10s  %s\n' \
  atoms wall_s host_GB dev_peak_GB budget_GB result
printf '%12s  %9s  %10s  %12s  %10s  %s\n' \
  ------------ --------- ---------- ------------ ---------- ------

for n in $sizes; do
  [ "$MAX" -gt 0 ] && [ "$n" -gt "$MAX" ] && continue

  d=$WORK/n$n
  rm -rf "$d"
  mkdir -p "$d"
  ln -sf "$DIR/atoms_$n.xyz" "$d/atoms.xyz"
  ln -sf "$DIR/gap_files" "$d/gap_files"
  cp "$DIR/$INPUT" "$d/input"

  # TG_LADDER_EXTRAS: an input.mad names an experimental data file, which
  # lives beside the systems and is not otherwise staged.
  for _x in "$DIR"/*; do
    case ${_x##*/} in atoms_*.xyz|input*|gap_files) continue ;; esac
    [ -f "$_x" ] && ln -sf "$_x" "$d/${_x##*/}"
  done

  # %M is peak RSS in kilobytes. GNU time only; the shell builtin has no
  # equivalent, which is why this calls /usr/bin/time explicitly.
  pre=""
  [ -n "$RANKS" ] && pre="mpirun -np $RANKS"
  ( cd "$d" && /usr/bin/time -f '%e %M' -o "$d/time.txt" \
      $pre "$BIN" "$MODE" >"$d/run.log" 2>&1 )
  rc=$?

  wall=$(awk '{print $1}' "$d/time.txt" 2>/dev/null)
  rss_kb=$(awk '{print $2}' "$d/time.txt" 2>/dev/null)
  host_gb=$(awk -v k="${rss_kb:-0}" 'BEGIN{printf "%.2f", k/1048576}')

  # The ledger line TurboGAP prints at the end of a run:
  #   GPUmem   this rank      X.XXX GB live in N buffers, peak Y.YYY GB
  dev_peak=$(grep -m1 'this rank' "$d/run.log" 2>/dev/null |
    sed -n 's/.*peak *\([0-9.]*\) GB.*/\1/p')
  # The budget, not the batch count: nothing prints how many SOAP batches
  # get_number_of_atom_pairs settled on, and a column headed "batches" showing a
  # number of gigabytes is worse than no column.
  budget=$(grep -m1 'max_Gbytes_per_process set from the device' "$d/run.log" 2>/dev/null |
    awk '{printf "%.2f", $NF}')
  [ -n "$budget" ] || budget=$(grep -m1 'budgeting is OFF' "$d/run.log" 2>/dev/null |
    awk '{printf "%.2f*", $(NF-1)}')

  if [ $rc -eq 0 ]; then
    result=ok
  elif grep -q 'GPUmem: device allocation' "$d/run.log" 2>/dev/null; then
    result="DEVICE OOM"
  elif grep -qiE 'Allocation would exceed|cannot allocate|Killed|std::bad_alloc' "$d/run.log" 2>/dev/null; then
    result="HOST OOM"
  else
    result="failed (exit $rc)"
  fi

  printf '%12s  %9s  %10s  %12s  %10s  %s\n' \
    "$n" "${wall:-?}" "$host_gb" "${dev_peak:-?}" "${budget:-?}" "$result"

  if [ $rc -ne 0 ] && [ "$KEEP_GOING" = 0 ]; then
    log ""
    log "Stopped at $n atoms. Tail of $d/run.log:"
    tail -20 "$d/run.log" | sed 's/^/  /'
    log ""
    log "Re-run with --keep-going to try the larger sizes anyway."
    exit 1
  fi
done

log ""
log "Working directories kept under $WORK"
