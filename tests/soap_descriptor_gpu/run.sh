#!/bin/bash
# cuda_get_soap_p and cuda_get_soap_der_one, spread over (site or pair, n, l)
# instead of one thread per site or pair, against the original loop orders kept
# in descverify.cu: the same sums in the same order, so the outputs must agree
# bit for bit. gpu_soap_pair_contract against the derivative chain it replaces,
# to a tolerance. Needs nvcc and a device; skips loudly without either.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)

CUDA_ARCH=${CUDA_ARCH:-sm_90}
HOP_ROOT=${HOP_ROOT:-$repo/src/hop}
GL_STUB=${GL_STUB:-$repo/makefiles/roihu-gl-stub}
WORK=${TURBOGAP_WORK:-/tmp/soap_descriptor_gpu.$$}
N_SITES=${N_SITES:-1000}
NEIGH=${NEIGH:-40}
REPS=${REPS:-5}

skip() { printf '    SKIP: %s\n' "$1"; exit 0; }

check_preconditions() {
    command -v nvcc > /dev/null || skip "no nvcc on PATH"
    command -v nvidia-smi > /dev/null || skip "no nvidia-smi on PATH"
    nvidia-smi -L > /dev/null 2>&1 || skip "no usable device"
}

build() {
    mkdir -p "$WORK"
    nvcc -x cu -O3 -arch="$CUDA_ARCH" --expt-extended-lambda --expt-relaxed-constexpr \
        -I"$HOP_ROOT" -I"$HOP_ROOT/source/hip" -DHOP_TARGET_CUDA -DCUDA -I"$GL_STUB" -I"$repo/src/gpu" \
        "$here/descverify.cu" -x cu "$repo/src/gpu/gpu_memory.cu" -lcurand -lcublas -o "$WORK/descverify"
}

main() {
    printf '\n== cuda_get_soap_p and _der_one against the originals, bit for bit; the contraction against the chain\n'
    check_preconditions
    build
    "$WORK/descverify" "$N_SITES" "$NEIGH" "$REPS"
    rm -rf "$WORK"
}

main
