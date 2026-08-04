#!/usr/bin/env bash
#
# Scheduled GPU regression run, intended to be driven by cron.
#
# cron gives us almost no environment, so everything is set explicitly here.
# By default this tests the binary that is already in bin/ rather than
# rebuilding: this is a live working tree, and a scheduled job should not
# quietly replace the binary you are using. Set REBUILD=1 to build first.
#
# Overlapping runs are skipped, and the run is skipped entirely if another
# turbogap process is already using the GPU, so the job never competes with
# interactive work.
#
# Usage:  cron_regression.sh          test the current bin/turbogap
#         REBUILD=1 cron_regression.sh  rebuild (DEBUG=0) and then test

set -u

REPO=${TURBOGAP_REPO:-/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_gpu_commit_mahti}
LOGDIR=${TURBOGAP_LOGDIR:-$REPO/tests/gpu/cron-logs}
KEEP=${TURBOGAP_LOG_KEEP:-30}
LOCK=${TMPDIR:-/tmp}/turbogap_gpu_regression.lock

export PATH=/usr/local/bin:/usr/bin:/bin:/usr/local/cuda/bin
export HOP_ROOT=${HOP_ROOT:-/u/74/zarrout1/unix/work/hop}
export TURBOGAP_ARCH=${TURBOGAP_ARCH:-Aalto_gfortran_openblas_hip_cuda}
export TURBOGAP_TEST_DATA=${TURBOGAP_TEST_DATA:-/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/input/CO}
export TURBOGAP_REF_BIN=${TURBOGAP_REF_BIN:-/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_master_2026/bin/turbogap}

mkdir -p "$LOGDIR"
stamp=$(date +%Y%m%d-%H%M%S)
log=$LOGDIR/$stamp.log

# Only one scheduled run at a time.
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "$(date -Is) another regression run holds the lock; skipping" >> "$LOGDIR/skipped.log"
  exit 0
fi

# Don't fight an interactive job for the GPU.
if pgrep -x turbogap >/dev/null 2>&1; then
  echo "$(date -Is) turbogap already running; skipping" >> "$LOGDIR/skipped.log"
  exit 0
fi

# Run in a subshell so that an early exit inside it still leaves the
# bookkeeping below (latest.log, history.log, log pruning) to run.
(
  echo "=== turbogap GPU regression: $(date -Is) on $(hostname) ==="
  cd "$REPO" || { echo "cannot cd to $REPO"; exit 2; }
  echo "repo   : $REPO"
  echo "branch : $(git rev-parse --abbrev-ref HEAD 2>/dev/null)"
  echo "commit : $(git rev-parse --short HEAD 2>/dev/null) $(git log -1 --pretty=%s 2>/dev/null)"
  echo "submod : $(git -C src/soap_turbo rev-parse --short HEAD 2>/dev/null)"

  if [ "${REBUILD:-0}" = 1 ]; then
    echo "--- rebuilding (DEBUG=0) ---"
    make deepclean >/dev/null 2>&1
    if ! make DEBUG=0 -j4; then
      echo "RESULT: BUILD FAILED"
      exit 1
    fi
  fi
  echo "binary : bin/turbogap ($(date -r bin/turbogap -Is 2>/dev/null))"

  echo "--- comparator self-test ---"
  python3 tests/gpu/test_compare_xyz.py || { echo "RESULT: COMPARATOR SELF-TEST FAILED"; exit 1; }

  echo "--- regression suite ---"
  if tests/gpu/run_regression.sh; then
    echo "RESULT: PASS"
  else
    echo "RESULT: FAIL"
    exit 1
  fi
) > "$log" 2>&1
rc=$?

ln -sfn "$log" "$LOGDIR/latest.log"
printf '%s %s %s\n' "$(date -Is)" "$([ $rc -eq 0 ] && echo PASS || echo FAIL)" "$(basename "$log")" \
  >> "$LOGDIR/history.log"

# Keep the log directory bounded.
ls -1t "$LOGDIR"/*.log 2>/dev/null | grep -v 'latest.log$' | tail -n +$((KEEP+1)) | xargs -r rm -f

exit $rc
