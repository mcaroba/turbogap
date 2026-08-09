#!/usr/bin/env bash
#
# Set up everything tools/profile_gpu.sh needs, and say clearly what is missing
# when something is not there.
#
# This is the profiling counterpart to tools/setup_dev_env.sh, and it is a
# separate script on purpose: that one builds the *formatting* environment
# (fprettify pinned to 0.3.7, the pre-commit hook) which every commit needs and
# which has nothing to do with the GPU. Profiling needs a different set of
# things, most of which are not Python at all, and merging the two would make
# the common case -- a contributor who only wants their commit formatted --
# depend on a CUDA installation.
#
# It manages the Python side with uv rather than venv+pip because uv resolves
# and installs into a fresh venv in a couple of seconds, and `uv pip compile`
# gives a lock file with hashes, so a profile analysed a year from now uses the
# same pandas that produced the numbers in the handoff docs. It is also the only
# way to get a *specific* Python here: alt's system interpreter is 3.12 and
# PEP-668 externally managed, so `pip install --user` fails outright (the same
# trap tools/setup_dev_env.sh documents).
#
# Usage:
#     tools/setup_profiling_env.sh              # create/refresh the venv
#     tools/setup_profiling_env.sh --check      # verify only, change nothing
#     TURBOGAP_PROF_VENV=/path tools/setup_profiling_env.sh
#
# Then put the venv on PATH; the script prints the exact line.

set -eu

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

VENV=${TURBOGAP_PROF_VENV:-$HOME/.venvs/turbogap-profiling}
PYTHON_VERSION=${TURBOGAP_PROF_PYTHON:-3.12}

check_only=0
[ "${1:-}" = "--check" ] && check_only=1

rc=0
ok()   { printf '  ok    %s\n' "$*"; }
warn() { printf '  WARN  %s\n' "$*"; }
miss() { printf '  MISS  %s\n' "$*"; rc=1; }
say()  { printf '%s\n' "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

# The analysis packages. Deliberately short:
#
#   tools/nsys_report.py reads the nsys SQLite export with the standard
#   library's sqlite3 and needs nothing from this list -- that is intentional,
#   so the core report cannot be blocked by a failed install on a machine where
#   the profile has already been captured.
#
#   These are for the parts where a library genuinely pays: --plot draws the
#   phase and kernel breakdown with matplotlib, and pandas is what you want in
#   a REPL when a report raises a question the canned queries do not answer.
PACKAGES="pandas matplotlib"

# ------------------------------------------------------------------ tool check
#
# Reported whether or not the venv is being built, because these are the things
# that actually stop a profile, and every one of them has a specific fix.
check_tools() {
  say ""
  say "Profiling tools"

  if command -v nsys >/dev/null 2>&1; then
    ok "nsys $(nsys --version 2>&1 | sed -n 's/.*version //p' | head -1)"
  else
    miss "nsys not on PATH -- install Nsight Systems, or load a cuda module"
  fi

  if command -v ncu >/dev/null 2>&1; then
    ok "ncu $(ncu --version 2>&1 | sed -n 's/^Version //p' | head -1)"
  else
    warn "ncu not on PATH -- per-kernel counters unavailable (nsys still works)"
  fi

  if command -v nvcc >/dev/null 2>&1; then
    ok "nvcc $(nvcc --version 2>&1 | sed -n 's/.*release \([0-9.]*\).*/\1/p' | head -1)"
  else
    miss "nvcc not on PATH -- the GPU build needs it"
  fi

  say ""
  say "Device"

  # nvidia-smi is NOT the test. It talks to the driver through NVML, and NVML
  # has its own version handshake: on a host whose userspace libraries have been
  # updated under a running kernel module, nvidia-smi fails with
  # "Driver/library version mismatch" while the CUDA runtime keeps working
  # perfectly. That is the state alt is in, and taking nvidia-smi's word for it
  # would mean refusing to profile a machine that profiles fine. Ask the CUDA
  # runtime instead, which is the thing TurboGAP actually uses.
  local probe=${TMPDIR:-/tmp}/tg_cuda_probe.$$
  if command -v nvcc >/dev/null 2>&1; then
    cat >"$probe.cu" <<'EOF'
#include <cstdio>
int main() {
  int n = 0;
  cudaError_t e = cudaGetDeviceCount(&n);
  if (e != cudaSuccess) { printf("ERR %s\n", cudaGetErrorString(e)); return 1; }
  if (n == 0) { printf("ERR no CUDA device visible\n"); return 1; }
  cudaDeviceProp p;
  cudaGetDeviceProperties(&p, 0);
  printf("OK %s sm_%d%d %d SMs %.1f GB\n", p.name, p.major, p.minor,
         p.multiProcessorCount, p.totalGlobalMem / 1073741824.0);
  return 0;
}
EOF
    if nvcc -o "$probe" "$probe.cu" >/dev/null 2>&1 && out=$("$probe" 2>&1); then
      ok "CUDA runtime: ${out#OK }"
    else
      miss "CUDA runtime cannot open a device: ${out:-compile failed}"
    fi
    rm -f "$probe" "$probe.cu"
  fi

  if command -v nvidia-smi >/dev/null 2>&1 && ! nvidia-smi >/dev/null 2>&1; then
    warn "nvidia-smi fails but the CUDA runtime works -- an NVML userspace/kernel"
    warn "      version mismatch. Harmless for profiling; nsys and ncu do not use NVML."
  fi

  # ncu reads hardware performance counters, and on most distributions that is
  # admin-only by default (NVIDIA's ERR_NVGPUCTRPERM). nsys does not need it, so
  # this is a warning: it costs you --tool ncu, not the whole workflow.
  if command -v ncu >/dev/null 2>&1; then
    local rp=/proc/driver/nvidia/params
    if [ -r $rp ] && grep -q '^RmProfilingAdminOnly: 0' $rp 2>/dev/null; then
      ok "GPU performance counters readable by non-root (ncu will work)"
    elif [ -r $rp ]; then
      warn "GPU performance counters are admin-only. ncu will fail with"
      warn "      ERR_NVGPUCTRPERM. Fix (needs root, then reboot or reload the module):"
      warn "        echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' \\"
      warn "          | sudo tee /etc/modprobe.d/nvidia-profiling.conf"
      warn "      nsys is unaffected."
    fi
  fi

  say ""
  say "Build environment"
  if [ -n "${HOP_ROOT:-}" ] && [ -d "${HOP_ROOT:-}" ]; then
    ok "HOP_ROOT=$HOP_ROOT"
  else
    miss "HOP_ROOT unset or not a directory -- the GPU build refuses to start"
  fi
  if [ -n "${TURBOGAP_ARCH:-}" ] && [ -f "$repo/makefiles/Makefile.${TURBOGAP_ARCH:-}" ]; then
    ok "TURBOGAP_ARCH=$TURBOGAP_ARCH"
  else
    miss "TURBOGAP_ARCH unset or no makefiles/Makefile.\$TURBOGAP_ARCH"
  fi
}

# ------------------------------------------------------------------ check mode
if [ "$check_only" = 1 ]; then
  say "Python analysis environment"
  if [ -x "$VENV/bin/python" ]; then
    ok "venv $VENV ($("$VENV/bin/python" -V 2>&1))"
    for p in $PACKAGES; do
      if "$VENV/bin/python" -c "import ${p%%[*}" >/dev/null 2>&1; then
        ok "$p $("$VENV/bin/python" -c "import ${p%%[*} as m; print(m.__version__)" 2>/dev/null)"
      else
        miss "$p not importable"
      fi
    done
  else
    miss "no venv at $VENV -- run tools/setup_profiling_env.sh"
  fi
  case ":$PATH:" in
    *":$VENV/bin:"*) ok "$VENV/bin is on PATH" ;;
    *)               warn "$VENV/bin is NOT on PATH (only needed for --plot and ad-hoc analysis)" ;;
  esac
  check_tools
  say ""
  [ $rc -eq 0 ] && say "Ready." || say "Some checks failed (see MISS above)."
  exit $rc
fi

# ------------------------------------------------------------------- build mode
if ! command -v uv >/dev/null 2>&1; then
  # Not installed silently. uv's installer writes to $HOME and edits shell
  # profiles; that is the user's decision, not this script's.
  die "uv is not on PATH. Install it with one of

    curl -LsSf https://astral.sh/uv/install.sh | sh     # then restart your shell
    pipx install uv
    python3 -m pip install --user uv

then re-run this script. (On alt it is already at ~/.local/bin/uv -- that
directory may simply not be on PATH.)"
fi

say "uv $(uv --version | awk '{print $2}')"
say "Creating $VENV with Python $PYTHON_VERSION"

# --python <version> makes uv fetch and manage that interpreter if the system
# does not have it, which is what keeps this reproducible across the login node,
# a compute node and a laptop.
uv venv --python "$PYTHON_VERSION" "$VENV"

# A lock file, refreshed only when it does not exist, so re-running this script
# never silently upgrades the packages under a set of numbers already recorded
# in the docs. Delete it deliberately to move versions.
lock=$repo/tools/profiling-requirements.txt
if [ ! -f "$lock" ]; then
  say "Resolving $PACKAGES -> $lock"
  # shellcheck disable=SC2086
  printf '%s\n' $PACKAGES |
    uv pip compile --quiet --generate-hashes --python "$VENV/bin/python" - -o "$lock"
  say "Wrote $lock. Commit it: it is what makes two profiles comparable."
else
  say "Using existing $lock (delete it to re-resolve)"
fi

VIRTUAL_ENV=$VENV uv pip install --quiet --require-hashes -r "$lock"

say ""
say "Add the venv to PATH:"
say ""
say "    export PATH=\"$VENV/bin:\$PATH\""
say ""

check_tools
say ""
if [ $rc -eq 0 ]; then
  say "Ready. Next:  tools/profile_gpu.sh --list"
else
  say "The Python environment is built, but some tools are missing (MISS above)."
fi
exit $rc
