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

wt=$(mktemp -d "${TMPDIR:-/tmp}/turbogap_baseline.XXXXXX")
cleanup() { cd "$repo"; git worktree remove --force "$wt" >/dev/null 2>&1 || true; rm -rf "$wt"; }
trap cleanup EXIT

echo "building baseline from $parent (soap_turbo $submod)"
cd "$repo"
git worktree add --detach "$wt" "$parent" >/dev/null

cd "$wt"
git submodule update --init --recursive >/dev/null 2>&1 || true
( cd src/soap_turbo && git checkout --quiet "$submod" )

# Use the same arch makefile the working tree is configured with, so the
# baseline and the binary under test differ only in source.
arch=$(sed -n 's|^include makefiles/Makefile\.\(.*\)$|\1|p' "$repo/Makefile" | head -1)
[ -n "$arch" ] || { echo "ERROR: could not determine arch from $repo/Makefile" >&2; exit 2; }
sed -i "s|^include makefiles/Makefile\..*|include makefiles/Makefile.$arch|" Makefile

make clean >/dev/null
make all >/dev/null
cp bin/turbogap "$out"
echo "wrote $out"
