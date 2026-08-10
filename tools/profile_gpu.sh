#!/usr/bin/env bash
#
# Build TurboGAP for profiling, run a test-suite case under the NVIDIA
# profilers, and print the analysis.
#
#     tools/profile_gpu.sh --list
#     tools/profile_gpu.sh CO_predict
#     tools/profile_gpu.sh --tool ncu --kernel 'soap|gap' CO_predict
#     tools/profile_gpu.sh --no-build --out /scratch/prof XRD_mad
#
# The cases come from tests/gpu/cases.sh, which is the same file
# tests/gpu/run_regression.sh reads. A profile of an input that is not the one
# the suite checks is worth very little, and a pasted copy of an input starts
# ageing the moment it is pasted, so there is exactly one definition.
#
# Three things this arranges that are easy to get wrong by hand:
#
#   1. It builds with PROFILE=1 DEBUG=0. The default for this tree is DEBUG=1,
#      which compiles device code with -G: about 2.1x slower and, worse,
#      differently ranked, so a DEBUG profile does not describe the build you
#      run. The Makefile refuses the combination outright.
#
#   2. The profiling build has its own object tree, tagged by the flags it was
#      built with (bin-profile/, bin-omp-profile/). make rebuilds on timestamps,
#      not on flags, so a PROFILE=1 build over an existing build/ would relink
#      the old objects, produce a binary with neither line info nor NVTX ranges,
#      and report success.
#
#   3. It runs in a staged directory and leaves it there. TurboGAP writes its
#      output next to its input, and a profile you cannot go back and re-query
#      is a profile you will capture again.
#
# Output lands in <out>/<case>/ :
#     input, atoms.xyz, gap_files    the staged case
#     run.log                        TurboGAP's own output, including its timers
#     <case>.nsys-rep                the timeline; open it in the Nsight GUI
#     <case>.sqlite                  the same data, queryable
#     stats/*.csv                    nsys's own summaries
#     report.txt                     the digest tools/nsys_report.py prints
#     <case>.ncu-rep                 (--tool ncu) per-kernel counters

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

TOOL=nsys
BUILD=1
OUT=${TURBOGAP_PROF_OUT:-$repo/profiling}
KERNEL_FILTER=
LAUNCH_COUNT=3
# basic, not speedoflight. `ncu --list-sets` offers basic / detailed / full /
# roofline; SpeedOfLight is a SECTION inside them, not a set. Given a set name it
# does not recognise ncu does NOT fail: it warns once, "No metrics to collect
# found in sections", and then writes a normal ~1 MB report containing no metric
# values at all. Every kernel appears with empty columns and the CSV export is
# empty, so the run looks like it worked. This default cost a whole GH200 job
# before it was noticed.
NCU_SET=${TURBOGAP_NCU_SET:-basic}
JOBS=${TURBOGAP_MAKE_JOBS:-12}
OPENMP=0
OMP_THREADS=
SUFFIX=
overrides=()
MPI_RANKS=${TURBOGAP_MPI_RANKS:-}
cases=()

log() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

usage() {
  sed -n '2,45p' "$0" | sed 's/^#\{0,1\} \{0,1\}//'
  cat <<'EOF'

Options:
  --list                 print the case names and exit
  --tool nsys|ncu|both   which profiler to run          (default: nsys)
  --no-build             use the existing bin[-omp]-profile/turbogap
  --out DIR              where results go                (default: <repo>/profiling)
  --kernel REGEX         ncu: profile only matching kernels
  --launch-count N       ncu: GLOBAL cap on profiled launches, not per kernel
                         (default: 3). Pair it with --kernel, or you profile
                         only whatever launches first.
  --ncu-set SET          ncu: metric set (basic|detailed|full|roofline)  (default: basic)
  --jobs N               make -j                         (default: 12)
  --ranks N              run under mpirun with N ranks
  --openmp [N]           build OPENMP=1 and run with N host threads. This is
                         the switch that makes the per-thread GPU streams real;
                         without it n_omp is 1 and the batched loops are serial
                         however many batches the input asks for.
  --input 'k = v'        override one keyword in the staged input; repeatable.
                         gpu_n_batches is the one to reach for: it defaults to
                         1, gpu_context_init clamps n_omp to it, and with one
                         batch there is nothing for a second stream to do. So
                         --openmp without --input 'gpu_n_batches = N' measures
                         a serial loop with extra threads sitting idle.
  --suffix S             write results to <out>/<case>S/, so two variants of
                         one case can be compared side by side.
EOF
}

while [ $# -gt 0 ]; do
  case $1 in
  --list) LIST=1; shift ;;
  --tool) TOOL=$2; shift 2 ;;
  --no-build) BUILD=0; shift ;;
  --build) BUILD=1; shift ;;
  --out) OUT=$2; shift 2 ;;
  --kernel) KERNEL_FILTER=$2; shift 2 ;;
  --launch-count) LAUNCH_COUNT=$2; shift 2 ;;
  --ncu-set) NCU_SET=$2; shift 2 ;;
  --jobs) JOBS=$2; shift 2 ;;
  --ranks) MPI_RANKS=$2; shift 2 ;;
  --input) overrides+=("$2"); shift 2 ;;
  --suffix) SUFFIX=$2; shift 2 ;;
  --openmp)
    OPENMP=1
    case ${2:-} in [0-9]*) OMP_THREADS=$2; shift 2 ;; *) shift ;; esac
    ;;
  -h | --help) usage; exit 0 ;;
  -*) die "unknown option: $1 (try --help)" ;;
  *) cases+=("$1"); shift ;;
  esac
done

case $TOOL in nsys | ncu | both) ;; *) die "--tool must be nsys, ncu or both" ;; esac

# ---------------------------------------------------------------- the cases
DATA=${TURBOGAP_TEST_DATA:-$repo/../turbogap_tests/CO}
TG_TESTS_DIR=${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}
TG_DATA=$DATA
. "$repo/tests/gpu/cases.sh"

if [ -n "${LIST:-}" ]; then
  show() {
    while read -r name; do
      tg_case_def "$name"
      mark=' '
      [ -e "$TG_DATA_DIR/$TG_ATOMS" ] || mark='!'
      printf ' %s%-14s %-8s %s\n' "$mark" "$name" "$TG_MODE" "$TG_DATA_DIR/$TG_ATOMS"
    done
  }
  echo "regression cases (also run by tests/gpu/run_regression.sh)"
  tg_case_list | show
  echo
  echo "large systems (profiling only -- no reference exists at these sizes)"
  tg_case_list_large | show
  echo
  echo " ! = data file not present"
  exit 0
fi

[ ${#cases[@]} -gt 0 ] || { usage; exit 2; }

if [ ! -d "$DATA" ]; then
  "$repo/tests/fetch_test_data.sh" "$TG_TESTS_DIR" || die "could not fetch test data"
fi

# ---------------------------------------------------------------- the build
# The Makefile tags the object and bin trees with the flag combination, so this
# has to agree with it or --no-build picks up the other variant's binary.
TAG=
[ "$OPENMP" = 1 ] && TAG=$TAG-omp
TAG=$TAG-profile
BIN=$repo/bin$TAG/turbogap

if [ "$BUILD" = 1 ]; then
  [ -n "${HOP_ROOT:-}" ] || die "HOP_ROOT is not set"
  [ -n "${TURBOGAP_ARCH:-}" ] || die "TURBOGAP_ARCH is not set"
  log "== building PROFILE=1 DEBUG=0 OPENMP=$OPENMP (-> bin$TAG/turbogap) =="
  make -C "$repo" -j"$JOBS" PROFILE=1 DEBUG=0 OPENMP="$OPENMP" \
    >"${TMPDIR:-/tmp}/tg_profile_build.log" 2>&1 ||
    { tail -30 "${TMPDIR:-/tmp}/tg_profile_build.log"; die "build failed"; }
  log "   ok"
fi
[ -x "$BIN" ] || die "no profiling binary at $BIN (drop --no-build, or build it)"

# The check that the build is the one we asked for. A binary with no NVTX
# symbol was linked from stale objects, and every NVTX-keyed number in the
# report would then be silently empty rather than wrong -- the failure mode that
# is hardest to notice.
if ! nm -D "$BIN" 2>/dev/null | grep -q nvtxRangePush; then
  die "$BIN has no NVTX symbol: it was not built with PROFILE=1.
    Rebuild:  make -C $repo PROFILE=1 DEBUG=0"
fi

# Require only the profiler actually being run. Demanding nsys unconditionally
# made --tool ncu impossible on Roihu, where nsys is not installed on the GPU
# nodes at all -- and ncu is the only profiler that works there.
if [ "$TOOL" != ncu ]; then
  command -v nsys >/dev/null || die "nsys not found; run tools/setup_profiling_env.sh --check"
fi
if [ "$TOOL" != nsys ]; then
  command -v ncu >/dev/null || die "ncu not found; run tools/setup_profiling_env.sh --check"
fi

# ---------------------------------------------------------------- run a case
run_under_nsys() {
  local dir=$1 name=$2
  log "-- nsys"
  # -t cuda,nvtx,osrt: the device work, our phase markers, and the OS calls the
  #    host blocks in. osrt is what shows a synchronize as a wait rather than as
  #    a gap with no explanation.
  # --cuda-memory-usage: device allocation over time. TurboGAP mallocs and frees
  #    per batch and per descriptor, and whether that is a cost is a question a
  #    timeline alone cannot answer.
  # --capture-range=none: profile the whole run. These cases are seconds long,
  #    and the setup phase is exactly where one-off costs hide.
  local pre=()
  [ -n "$MPI_RANKS" ] && pre=(mpirun -np "$MPI_RANKS")
  [ -n "$OMP_THREADS" ] && export OMP_NUM_THREADS=$OMP_THREADS
  (cd "$dir" && nsys profile \
    -t cuda,nvtx,osrt \
    --cuda-memory-usage=true \
    --force-overwrite=true \
    -o "$dir/$name" \
    "${pre[@]}" "$BIN" "$TG_MODE") >"$dir/run.log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    log "   run failed (exit $rc); tail of run.log:"
    tail -20 "$dir/run.log" | sed 's/^/     /'
    return 1
  fi
  grep -q GPUassert "$dir/run.log" && {
    log "   NOTE: the run reported a GPU error; the profile is of a broken run"
    grep -m3 GPUassert "$dir/run.log" | sed 's/^/     /'
  }

  mkdir -p "$dir/stats"
  # Export once and report from the export. `nsys stats` re-exports the .nsys-rep
  # on every invocation otherwise, which on the bigger cases costs more than the
  # profiled run did.
  nsys export --type sqlite --force-overwrite=true \
    -o "$dir/$name.sqlite" "$dir/$name.nsys-rep" >>"$dir/run.log" 2>&1 ||
    { log "   nsys export failed"; return 1; }

  #   cuda_gpu_kern_sum   which kernels, and for how long
  #   cuda_gpu_mem_time_sum / _size_sum   the transfers, in time and in bytes
  #   cuda_api_sum        the host side of the same story: where it blocks
  #   nvtx_pushpop_sum    our phases, on the host
  #   nvtx_gpu_proj_sum   our phases PROJECTED ONTO THE DEVICE -- the one that
  #                       answers "how busy is the GPU during SOAP", which is
  #                       the question the whole stream discussion turns on
  #   cuda_kern_exec_sum  launch overhead vs execution, per kernel
  local reports=(cuda_gpu_kern_sum cuda_gpu_mem_time_sum cuda_gpu_mem_size_sum
                 cuda_api_sum nvtx_pushpop_sum nvtx_gpu_proj_sum cuda_kern_exec_sum)
  local r
  for r in "${reports[@]}"; do
    nsys stats --report "$r" --format csv --force-export=false \
      -o "$dir/stats/$name" "$dir/$name.sqlite" >>"$dir/run.log" 2>&1 ||
      log "   (report $r unavailable)"
  done

  python3 "$here/nsys_report.py" "$dir/$name.sqlite" | tee "$dir/report.txt"
}

run_under_ncu() {
  local dir=$1 name=$2
  log "-- ncu (set=$NCU_SET, launch-count=$LAUNCH_COUNT${KERNEL_FILTER:+, kernel~/$KERNEL_FILTER/})"
  # --launch-count is not a nicety. ncu serialises the context and replays each
  # kernel once per metric pass; on CO_predict, which launches thousands of
  # kernels, an unbounded run does not finish in a useful time. Profiling a few
  # launches of each distinct kernel gives the same per-kernel numbers.
  local args=(--set "$NCU_SET" --launch-count "$LAUNCH_COUNT"
              --target-processes all --force-overwrite -o "$dir/$name")
  [ -n "$KERNEL_FILTER" ] && args+=(--kernel-name-base function --kernel-name "regex:$KERNEL_FILTER")
  (cd "$dir" && ncu "${args[@]}" "$BIN" "$TG_MODE") >"$dir/ncu.log" 2>&1
  local rc=$?
  if grep -q ERR_NVGPUCTRPERM "$dir/ncu.log" 2>/dev/null; then
    log "   ncu was refused access to the GPU performance counters."
    log "   Fix, then reboot or reload the nvidia module:"
    log "     echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' \\"
    log "       | sudo tee /etc/modprobe.d/nvidia-profiling.conf"
    return 1
  fi
  [ $rc -eq 0 ] || { log "   ncu failed (exit $rc):"; tail -20 "$dir/ncu.log" | sed 's/^/     /'; return 1; }

  ncu --import "$dir/$name.ncu-rep" --page details --csv >"$dir/ncu_details.csv" 2>/dev/null
  log "   wrote $dir/$name.ncu-rep"
  # Speed-of-light in one line per kernel: the first question about any kernel
  # is whether it is compute bound, memory bound, or neither (i.e. latency
  # bound), and that is what decides whether streams would help it.
  ncu --import "$dir/$name.ncu-rep" --page details 2>/dev/null |
    grep -E '^\s+(Compute \(SM\) Throughput|Memory Throughput|Duration|void |[a-zA-Z_].*\(.*\))' |
    head -60 | sed 's/^/     /'
}

status=0
for name in "${cases[@]}"; do
  log ""
  log "======================================================================"
  log "== $name"
  log "======================================================================"
  tg_case_def "$name" || { status=1; continue; }
  if [ ! -e "$TG_DATA_DIR/$TG_ATOMS" ]; then
    log "   SKIP: no data at $TG_DATA_DIR/$TG_ATOMS"
    continue
  fi
  dir=$(tg_stage_case "$name" "$OUT")
  if [ -n "$SUFFIX" ]; then
    rm -rf "$dir$SUFFIX" && mv "$dir" "$dir$SUFFIX" && dir=$dir$SUFFIX
  fi
  # Keyword overrides. Replaced in place when the keyword is already there and
  # appended when it is not, so `--input 'gpu_n_batches = 4'` works whether or
  # not the case sets it -- and the file keeps the case's own value visible in
  # the diff rather than carrying two conflicting lines, which TurboGAP resolves
  # by last-one-wins and is a trap to read back.
  for ov in ${overrides+"${overrides[@]}"}; do
    key=${ov%%=*}; key=${key%% }
    if grep -qE "^[[:space:]]*$key[[:space:]]*=" "$dir/input"; then
      sed -i -E "s|^[[:space:]]*$key[[:space:]]*=.*|$ov|" "$dir/input"
    else
      printf '%s\n' "$ov" >>"$dir/input"
    fi
    log "   input override: $ov"
  done
  log "   staged in $dir"

  case $TOOL in
  nsys) run_under_nsys "$dir" "$name" || status=1 ;;
  ncu) run_under_ncu "$dir" "$name" || status=1 ;;
  both)
    run_under_nsys "$dir" "$name" || status=1
    run_under_ncu "$dir" "$name" || status=1
    ;;
  esac
done

log ""
log "Results under $OUT"
exit $status
