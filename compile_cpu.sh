#!/usr/bin/env bash
#
# Build the host binary, from a fresh clone or an existing one.
#
# Produces bin/turbogap. Nothing here is required -- `make -j` does the same --
# but it fetches the submodules first and checks the few things whose absence
# produces a failure much later and further away than its cause.
set -euo pipefail

# ---------------------------------------------------------------- parameters

ARCH=${TURBOGAP_ARCH:-Ubuntu_gfortran_mpi}
JOBS=${JOBS:-12}
# 1 gives bin-dbg/turbogap, bounds-checked, in its own object tree.
DEBUG=${DEBUG:-0}

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ----------------------------------------------------------------- functions

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

check_arch() {
    [ -f "$repo/makefiles/Makefile.$ARCH" ] || {
        printf 'no such arch: %s\navailable:\n' "$ARCH" >&2
        ls "$repo/makefiles" | sed 's|^Makefile\.|    |' >&2
        exit 1
    }
    grep -q '^GPU = 1' "$repo/makefiles/Makefile.$ARCH" &&
        die "$ARCH is a device arch; use compile_gpu.sh"
    return 0
}

check_compiler() {
    command -v mpif90 >/dev/null || die "mpif90 not found; load an MPI module or install openmpi"
}

fetch_submodules() {
    # src/soap_turbo carries the SOAP routines; without it the build fails on a
    # missing source file rather than on a missing checkout.
    git -C "$repo" submodule update --init --recursive
}

build() {
    make -C "$repo" -j "$JOBS" TURBOGAP_ARCH="$ARCH" DEBUG="$DEBUG"
}

report() {
    local bin=$repo/bin/turbogap
    [ "$DEBUG" = 0 ] || bin=$repo/bin-dbg/turbogap
    [ -x "$bin" ] || die "build reported success but $bin is not there"
    printf '\nbuilt %s\n' "$bin"
    printf 'test with: tests/regression/run.sh\n'
}

# ------------------------------------------------------------------- the run

check_arch
check_compiler
fetch_submodules
build
report
