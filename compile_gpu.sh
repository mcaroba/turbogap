#!/usr/bin/env bash
#
# Build the device binary, from a fresh clone or an existing one.
#
# Produces bin-gpu/turbogap, beside the host binary rather than instead of it:
# the regression suite's --both mode runs the two against each other and needs
# them at the same time.
set -euo pipefail

# ---------------------------------------------------------------- parameters

ARCH=${TURBOGAP_ARCH:-Aalto_gfortran_openblas_hip_cuda}
JOBS=${JOBS:-8}
# The device arch makefiles default DEBUG to 1, which is -G and about 2.1x
# slower. Never profile or time a DEBUG=1 build; see docs/PROFILING.md.
DEBUG=${DEBUG:-0}
# Its own object tree, so the host build is not overwritten.
TAG=${BUILD_TAG_EXTRA:--gpu}

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ----------------------------------------------------------------- functions

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

check_arch() {
    [ -f "$repo/makefiles/Makefile.$ARCH" ] || {
        printf 'no such arch: %s\navailable device arches:\n' "$ARCH" >&2
        grep -l '^GPU = 1' "$repo"/makefiles/Makefile.* | sed 's|.*/Makefile\.|    |' >&2
        exit 1
    }
    grep -q '^GPU = 1' "$repo/makefiles/Makefile.$ARCH" ||
        die "$ARCH is not a device arch; use compile_cpu.sh"
    return 0
}

check_toolchain() {
    command -v mpif90 >/dev/null || die "mpif90 not found; load an MPI module"
    command -v nvcc >/dev/null || command -v hipcc >/dev/null ||
        die "neither nvcc nor hipcc found; load a CUDA or ROCm module"
}

fetch_submodules() {
    # A device build takes its SOAP routines from src/soap_turbo_gpu, whose
    # get_soap signature differs from the host one -- it is a second
    # implementation, not a later version, so both submodules are needed.
    git -C "$repo" submodule update --init --recursive
}

build() {
    make -C "$repo" -j "$JOBS" TURBOGAP_ARCH="$ARCH" DEBUG="$DEBUG" BUILD_TAG_EXTRA="$TAG"
}

report() {
    local bin=$repo/bin$TAG/turbogap
    [ -x "$bin" ] || die "build reported success but $bin is not there"
    printf '\nbuilt %s\n' "$bin"
    printf 'test with: tests/regression/run.sh --gpu     (or --both)\n'
}

# ------------------------------------------------------------------- the run

check_arch
check_toolchain
fetch_submodules
build
report
