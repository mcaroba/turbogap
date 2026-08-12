#!/usr/bin/env bash
#
# Rebuild the reference binary the regression suite compares against.
#
# The binary itself is deliberately not committed: it is 1.6 MB, it is
# architecture-specific, and a stale one committed by accident would make the
# whole suite pass vacuously. baseline/COMMIT records what it must be built
# from, and this script rebuilds exactly that in a detached worktree so the
# branch you are working on is never touched.
#
# Usage: tests/regression/make_baseline.sh

set -eu

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)
commit_file=$here/baseline/COMMIT

[ -f "$commit_file" ] || { echo "ERROR: $commit_file missing" >&2; exit 2; }

parent=$(sed -n 1p "$commit_file")
submod=$(sed -n 2p "$commit_file")
out=$here/baseline/turbogap.$(printf '%.7s' "$parent")

if [ -x "$out" ]; then
  echo "baseline already present: $out"
  echo "(delete it to force a rebuild)"
  exit 0
fi

# Which arch makefile to build with. The baseline and the binary under test
# have to differ only in source, so this follows whatever the working tree is
# configured with -- but the two Makefiles say it differently, and reading the
# wrong one produces a build that fails several minutes later with
# `makefiles/Makefile.: No such file or directory`:
#
#   current   TURBOGAP_ARCH ?= Ubuntu_gfortran_mpi
#             include makefiles/Makefile.$(TURBOGAP_ARCH)
#   baseline  include makefiles/Makefile.Ubuntu_gfortran_mpi
#
# So: the environment first, then the variable, then the literal include, and
# reject anything still holding an unexpanded make variable.
arch=${TURBOGAP_ARCH:-}
if [ -z "$arch" ]; then
  arch=$(sed -n 's|^ *TURBOGAP_ARCH *[?:]\{0,1\}= *\([^ ]*\) *$|\1|p' "$repo/Makefile" | head -1)
fi
if [ -z "$arch" ]; then
  arch=$(sed -n 's|^include makefiles/Makefile\.\(.*\)$|\1|p' "$repo/Makefile" | head -1)
fi
case ${arch:-} in
'' | *'$('*)
  echo "ERROR: could not determine the arch from $repo/Makefile" >&2
  echo "       set TURBOGAP_ARCH to one of:" >&2
  ls "$repo/makefiles" | sed 's|^Makefile\.|         |' >&2
  exit 2
  ;;
esac
[ -f "$repo/makefiles/Makefile.$arch" ] || {
  echo "ERROR: no such arch makefile: makefiles/Makefile.$arch" >&2
  exit 2
}

# Both commits have to be in the object store before anything is built. Without
# this the failure is git's bare `fatal: invalid reference: <sha>`, which says
# nothing about which of the two it means or what to do about it. A shallow or
# freshly cloned checkout is the usual reason.
cd "$repo"
git rev-parse --verify --quiet "$parent^{commit}" >/dev/null || {
  echo "ERROR: commit $parent is not in this clone" >&2
  echo "       It is the pre-refactor master HEAD the suite compares against," >&2
  echo "       recorded in $commit_file. Fetch it:" >&2
  echo "           git fetch --unshallow    # if this is a shallow clone" >&2
  echo "           git fetch origin" >&2
  exit 2
}
if [ -e "$repo/src/soap_turbo/.git" ]; then
  git -C "$repo/src/soap_turbo" rev-parse --verify --quiet "$submod^{commit}" >/dev/null || {
    echo "ERROR: soap_turbo commit $submod is not in the submodule clone" >&2
    echo "       Fetch it:" >&2
    echo "           git -C src/soap_turbo fetch origin" >&2
    exit 2
  }
fi

wt=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_baseline.XXXXXX")
cleanup() { cd "$repo"; git worktree remove --force "$wt" >/dev/null 2>&1 || true; rm -rf "$wt"; }
trap cleanup EXIT

echo "building baseline from $parent (soap_turbo $submod), arch $arch"
cd "$repo"
# A worktree left registered by a run that was killed rather than trapped makes
# the next `worktree add` fail for a reason that has nothing to do with this one.
git worktree prune
git worktree add --detach "$wt" "$parent" >/dev/null

cd "$wt"
git submodule update --init --recursive >/dev/null 2>&1 || true
( cd src/soap_turbo && git checkout --quiet "$submod" )

sed -i "s|^include makefiles/Makefile\..*|include makefiles/Makefile.$arch|" Makefile
grep -q "^include makefiles/Makefile\.$arch\$" Makefile || {
  echo "ERROR: could not point the baseline Makefile at makefiles/Makefile.$arch" >&2
  exit 2
}

make clean >/dev/null
make all >/dev/null
cp bin/turbogap "$out"
echo "wrote $out"
