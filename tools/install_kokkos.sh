#!/usr/bin/env bash
#
# Build and install Kokkos, for the KOKKOS=1 device build.
#
# Why this is a script and not a README line: Kokkos is a header-and-library
# dependency that has to be compiled against the same CUDA and the same device
# architecture as TurboGAP itself, and three of those settings are easy to get
# silently wrong.
#
#   1. THE ARCHITECTURE IS BAKED IN. Kokkos compiles its CUDA backend for one
#      -arch, chosen at configure time by Kokkos_ARCH_<NAME>=ON. Pick a
#      different one from the CUDA_ARCH the arch makefile passes and the link
#      succeeds, the binary runs, and every Kokkos kernel launch fails at
#      runtime with "no kernel image is available". So the name is derived here
#      from the same sm_XX the makefile uses, rather than typed twice.
#
#   2. KOKKOS MUST BE BUILT WITH nvcc_wrapper, not with g++ or nvcc directly.
#      Kokkos ships it in bin/; it is a g++ shim that routes .cpp through nvcc
#      when the CUDA backend is on. Configuring with plain g++ produces a
#      host-only library that links and then has no device code in it.
#
#   3. cmake IS NOT ON THIS BOX. It is not a module and not in the base image,
#      but PyPI ships official binary wheels of it, so it installs into the
#      same venv as fprettify and i-PI with no root and no system change. That
#      is what setup_dev_env.sh does; this script uses it if it is there and
#      installs it if it is not.
#
# Usage:
#     tools/install_kokkos.sh                  # build and install
#     tools/install_kokkos.sh --check          # verify only, change nothing
#     KOKKOS_ROOT=/path tools/install_kokkos.sh
#     CUDA_ARCH=sm_80 tools/install_kokkos.sh
#
# Then export the line it prints. The KOKKOS=1 build reads KOKKOS_ROOT the same
# way the device build reads HOP_ROOT.

set -eu

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/.." && pwd)

# ---------------------------------------------------------------- parameters

VENV=${TURBOGAP_VENV:-$HOME/.venvs/turbogap-tools}
# Pinned for the same reason the formatters are: Kokkos changes which CUDA and
# host-compiler combinations it accepts between minor releases.
KOKKOS_VERSION=${KOKKOS_VERSION:-4.7.04}
KOKKOS_ROOT=${KOKKOS_ROOT:-$HOME/.local/kokkos-$KOKKOS_VERSION}
# Must match CUDA_ARCH in the device arch makefile; see note 1 above.
CUDA_ARCH=${CUDA_ARCH:-sm_86}
KOKKOS_SRC=${KOKKOS_SRC:-$HOME/.cache/turbogap/kokkos-$KOKKOS_VERSION}
JOBS=${JOBS:-8}

check_only=0
[ "${1:-}" = "--check" ] && check_only=1

# ----------------------------------------------------------------- functions

say() { printf '%s\n' "$*"; }

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

# Kokkos names the device generation, nvcc names the compute capability. One
# setting, two spellings; keep them derived from each other rather than typed.
kokkos_arch_name() {
    case "$1" in
        sm_60) echo PASCAL60 ;;
        sm_61) echo PASCAL61 ;;
        sm_70) echo VOLTA70 ;;
        sm_72) echo VOLTA72 ;;
        sm_75) echo TURING75 ;;
        sm_80) echo AMPERE80 ;;
        sm_86) echo AMPERE86 ;;
        sm_89) echo ADA89 ;;
        sm_90) echo HOPPER90 ;;
        *) return 1 ;;
    esac
}

find_cmake() {
    if [ -x "$VENV/bin/cmake" ]; then
        echo "$VENV/bin/cmake"
    elif command -v cmake >/dev/null 2>&1; then
        command -v cmake
    else
        return 1
    fi
}

install_cmake() {
    [ -d "$VENV" ] || die "no venv at $VENV; run tools/setup_dev_env.sh first"
    say "installing cmake into $VENV"
    "$VENV/bin/python" -m pip install --quiet cmake
}

fetch_kokkos() {
    if [ -d "$KOKKOS_SRC/.git" ]; then
        say "kokkos $KOKKOS_VERSION source already at $KOKKOS_SRC"
        return 0
    fi
    say "fetching kokkos $KOKKOS_VERSION"
    mkdir -p "$(dirname "$KOKKOS_SRC")"
    rm -rf "$KOKKOS_SRC"
    git clone --quiet --depth 1 --branch "$KOKKOS_VERSION" \
        https://github.com/kokkos/kokkos.git "$KOKKOS_SRC" ||
        die "could not clone kokkos $KOKKOS_VERSION"
}

build_kokkos() {
    local cmake_bin=$1
    local arch_name
    arch_name=$(kokkos_arch_name "$CUDA_ARCH") ||
        die "no Kokkos architecture known for CUDA_ARCH=$CUDA_ARCH; add it to kokkos_arch_name"

    command -v nvcc >/dev/null 2>&1 || die "nvcc is not on PATH; the CUDA backend needs it"

    say "building kokkos for $CUDA_ARCH (Kokkos_ARCH_$arch_name) into $KOKKOS_ROOT"
    rm -rf "$KOKKOS_SRC/build"
    "$cmake_bin" -S "$KOKKOS_SRC" -B "$KOKKOS_SRC/build" \
        -DCMAKE_INSTALL_PREFIX="$KOKKOS_ROOT" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER="$KOKKOS_SRC/bin/nvcc_wrapper" \
        -DCMAKE_CXX_STANDARD=17 \
        -DKokkos_ENABLE_SERIAL=ON \
        -DKokkos_ENABLE_CUDA=ON \
        -DKokkos_ENABLE_CUDA_LAMBDA=ON \
        -DKokkos_ENABLE_CUDA_CONSTEXPR=ON \
        -DKokkos_ARCH_"$arch_name"=ON \
        >"$KOKKOS_SRC/configure.log" 2>&1 ||
        { tail -30 "$KOKKOS_SRC/configure.log" >&2; die "cmake configure failed; full log in $KOKKOS_SRC/configure.log"; }

    "$cmake_bin" --build "$KOKKOS_SRC/build" -j "$JOBS" \
        >"$KOKKOS_SRC/build.log" 2>&1 ||
        { tail -30 "$KOKKOS_SRC/build.log" >&2; die "kokkos build failed; full log in $KOKKOS_SRC/build.log"; }

    "$cmake_bin" --install "$KOKKOS_SRC/build" >>"$KOKKOS_SRC/build.log" 2>&1 ||
        die "kokkos install failed; full log in $KOKKOS_SRC/build.log"
}

# An install is usable when the CMake package, the header and the static
# library are all present. Checking only the directory would pass on a build
# that failed halfway through installing.
check_install() {
    [ -f "$KOKKOS_ROOT/include/Kokkos_Core.hpp" ] || return 1
    [ -f "$KOKKOS_ROOT/lib/libkokkoscore.a" ] || [ -f "$KOKKOS_ROOT/lib64/libkokkoscore.a" ] || return 1
    return 0
}

report_check() {
    local rc=0
    if check_install; then
        say "  ok    kokkos $KOKKOS_VERSION at $KOKKOS_ROOT"
    else
        say "  MISS  no kokkos install at $KOKKOS_ROOT"
        rc=1
    fi
    if find_cmake >/dev/null 2>&1; then
        say "  ok    cmake $("$(find_cmake)" --version | head -1 | awk '{print $3}')"
    else
        say "  MISS  cmake is not in $VENV and not on PATH"
        rc=1
    fi
    return $rc
}

print_env() {
    say ""
    say "kokkos is installed. Add this to your environment:"
    say ""
    say "    export KOKKOS_ROOT=$KOKKOS_ROOT"
    say ""
    say "then build the device binary with Kokkos:"
    say ""
    say "    KOKKOS=1 ./compile_gpu.sh"
    say ""
}

# --------------------------------------------------------------------- wiring

if [ "$check_only" = 1 ]; then
    say "kokkos:"
    report_check
    exit $?
fi

if check_install; then
    say "kokkos $KOKKOS_VERSION is already installed at $KOKKOS_ROOT"
    print_env
    exit 0
fi

cmake_bin=$(find_cmake) || {
    install_cmake
    cmake_bin=$(find_cmake) || die "cmake still not found after installing it"
}
say "using cmake at $cmake_bin"

fetch_kokkos
build_kokkos "$cmake_bin"

check_install || die "kokkos build reported success but $KOKKOS_ROOT is incomplete"
print_env
