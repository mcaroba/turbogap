#!/usr/bin/env bash
#
# Make sure the regression test data is present.
#
# The test systems live in their own repository rather than in the source tree:
# they are ~35 MB of trajectories and potentials, they are shared by the CPU and
# GPU branches, and they change on a different cadence from the code.
#
#   https://github.com/TiganyZ/turbogap_tests
#
# Usage:
#   tests/fetch_test_data.sh [<target-dir>]
#
# With no argument the target is $TURBOGAP_TESTS_DIR, else <repo>/../turbogap_tests
# -- the same default the two run.sh scripts use, so a plain `run.sh` finds what
# this leaves behind.
#
# Environment:
#   TURBOGAP_TESTS_URL   clone from somewhere else (a fork, or a local mirror)
#   TURBOGAP_TESTS_REF   branch or tag to check out (default: the remote HEAD).
#                        Implies a full clone, since a shallow one has no other
#                        ref to check out.
#   TURBOGAP_TESTS_DEPTH clone depth (default: 1). The history of a data
#                        repository is of no use to a test run and is most of
#                        what there is to transfer, so the default is shallow.
#                        Set to 0 for a full clone.
#   TURBOGAP_TESTS_UPDATE=1
#                        also `git pull` when the clone already exists. Off by
#                        default: a test run should not silently change the data
#                        it is testing against.
#
# Exit status is non-zero if the data could not be made available.

set -u

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

URL=${TURBOGAP_TESTS_URL:-https://github.com/TiganyZ/turbogap_tests.git}
REF=${TURBOGAP_TESTS_REF:-}
DEPTH=${TURBOGAP_TESTS_DEPTH:-1}
TARGET=${1:-${TURBOGAP_TESTS_DIR:-$repo/../turbogap_tests}}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

command -v git >/dev/null || die "git is required to fetch the test data"

# A ref only exists to be checked out if the history came with it.
if [ -n "$REF" ] || [ "$DEPTH" = 0 ]; then
  depth_opt=""
else
  depth_opt="--depth $DEPTH"
fi

if [ -d "$TARGET/.git" ]; then
  if [ "${TURBOGAP_TESTS_UPDATE:-0}" = 1 ]; then
    printf 'updating test data in %s\n' "$TARGET"
    git -C "$TARGET" pull --ff-only || die "could not update $TARGET"
  else
    printf 'test data already present: %s\n' "$TARGET"
    printf '  (set TURBOGAP_TESTS_UPDATE=1 to pull)\n'
  fi
elif [ -d "$TARGET" ]; then
  # A directory that is not a clone. Could be a manual copy or a symlink into a
  # shared location; either way it is not ours to overwrite.
  if [ -n "$(ls -A "$TARGET" 2>/dev/null)" ]; then
    printf 'test data directory exists and is not a git clone: %s\n' "$TARGET"
    printf '  leaving it alone; remove it or set TURBOGAP_TESTS_DIR to fetch afresh\n'
  else
    rmdir "$TARGET" && git clone $depth_opt "$URL" "$TARGET" || die "clone into $TARGET failed"
  fi
else
  printf 'cloning %s\n  into %s\n' "$URL" "$TARGET"
  git clone $depth_opt "$URL" "$TARGET" || die "clone into $TARGET failed"
fi

if [ -n "$REF" ]; then
  git -C "$TARGET" checkout --quiet "$REF" || die "no such ref in the test data: $REF"
  printf 'checked out %s\n' "$REF"
fi

[ -d "$TARGET" ] || die "test data still missing after fetch: $TARGET"
printf 'test data ready: %s\n' "$TARGET"
