#!/usr/bin/env bash
#
# Create the tooling environment this repo's hooks and scripts expect, and
# install the pre-commit hook.
#
# Why this exists as a script rather than a line in a README: the environment
# is not reconstructible from `pip install pre-commit`. Three things about it
# are load-bearing and all three are easy to get wrong.
#
#   1. fprettify is PINNED to 0.3.7. The pre-commit hook pins it too, but a
#      manual `fprettify` run must agree with the hook or the two fight over
#      every file. Formatting drift between the CPU and GPU branches is not
#      cosmetic here -- it inflated the measured merge surface of
#      exp_interface.f90 by an order of magnitude (757 diverged lines reported
#      against a real 67), which is the whole reason the trees are normalised.
#
#   2. The host may be PEP-668 externally managed (alt is). `pip install
#      --user` fails there, and `--break-system-packages` is the wrong answer
#      -- it puts an unpinned fprettify on PATH ahead of nothing in particular
#      and makes the version depend on when you installed it. A venv is the
#      only arrangement that is both permitted and reproducible.
#
#   3. `pre-commit install` must run inside the clone. It is per-clone state in
#      .git/hooks, not something the committed .pre-commit-config.yaml provides,
#      so a fresh checkout has no hook at all and silently commits unformatted
#      code.
#
# Usage:
#     tools/setup_dev_env.sh              # default venv location
#     TURBOGAP_VENV=/path tools/setup_dev_env.sh
#     tools/setup_dev_env.sh --check      # verify only, change nothing
#
# Then add the venv's bin to PATH (the script prints the exact line).

set -eu

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

VENV=${TURBOGAP_VENV:-$HOME/.venvs/turbogap-tools}
FPRETTIFY_VERSION=0.3.7

check_only=0
[ "${1:-}" = "--check" ] && check_only=1

say() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- check mode
if [ "$check_only" = 1 ]; then
  rc=0
  if [ -x "$VENV/bin/fprettify" ]; then
    v=$("$VENV/bin/fprettify" --version 2>&1 | awk '{print $NF}')
    if [ "$v" = "$FPRETTIFY_VERSION" ]; then
      say "  ok    fprettify $v"
    else
      say "  WRONG fprettify $v (expected $FPRETTIFY_VERSION)"; rc=1
    fi
  else
    say "  MISS  fprettify -- run tools/setup_dev_env.sh"; rc=1
  fi

  if [ -x "$VENV/bin/pre-commit" ]; then
    say "  ok    pre-commit $("$VENV/bin/pre-commit" --version | awk '{print $NF}')"
  else
    say "  MISS  pre-commit -- run tools/setup_dev_env.sh"; rc=1
  fi

  if [ -f "$repo/.git/hooks/pre-commit" ]; then
    say "  ok    hook installed in $repo/.git/hooks/pre-commit"
  else
    say "  MISS  hook NOT installed -- commits will not be formatted"; rc=1
  fi

  case ":$PATH:" in
    *":$VENV/bin:"*) say "  ok    $VENV/bin is on PATH" ;;
    *)               say "  MISS  $VENV/bin is NOT on PATH"; rc=1 ;;
  esac

  exit $rc
fi

# ---------------------------------------------------------------- create venv
command -v python3 >/dev/null || die "python3 not found"

if [ ! -x "$VENV/bin/python" ] && [ ! -x "$VENV/bin/python3" ]; then
  say "creating venv at $VENV"
  python3 -m venv "$VENV" || die "python3 -m venv failed (is python3-venv installed?)"
else
  say "venv already at $VENV"
fi

PY="$VENV/bin/python"
[ -x "$PY" ] || PY="$VENV/bin/python3"

say "installing pinned tooling"
"$PY" -m pip install --quiet --upgrade pip
"$PY" -m pip install --quiet "fprettify==$FPRETTIFY_VERSION" pre-commit

# --------------------------------------------------------------- install hook
# pre-commit needs fprettify's own pinned copy, which it builds into its cache
# from additional_dependencies in .pre-commit-config.yaml -- that is separate
# from the venv copy above, and both are wanted: the hook uses its cache, a
# manual run uses the venv, and they must be the same version.
if [ -d "$repo/.git" ]; then
  say "installing the git hook"
  (cd "$repo" && PATH="$VENV/bin:$PATH" "$VENV/bin/pre-commit" install >/dev/null)
else
  say "NOTE: $repo is not a git checkout, skipping hook install"
fi

cat <<EOF

Done. Add this to your shell profile, or run it in each session:

    export PATH="$VENV/bin:\$PATH"

Verify at any time with:

    tools/setup_dev_env.sh --check

Notes:

  * fprettify rewrites files DURING 'git commit', and when it does the commit
    ABORTS. That is not a failure -- rebuild, re-run the suite, 'git add' and
    commit again. Confirm the rewrite was whitespace-only first:

        git diff -w --ignore-blank-lines    # must report nothing

  * To reformat everything after changing the settings (rare, and never on one
    branch alone -- see .pre-commit-config.yaml):

        pre-commit run --all-files

  * The hooks are ALLOW-LISTED to src/ tools/ docs/ makefiles/. They are not
    exclude-listed, because an exclude list fails open: the first version of
    the config excluded tests/*/expected and the fixers then rewrote
    silicon.gap, silicon.xml, the .beta tables and a case's input deck. Those
    are parsed inputs, where a stripped trailing space changes what the code
    reads. Only ever widen 'files:', never narrow 'exclude:'.
EOF
