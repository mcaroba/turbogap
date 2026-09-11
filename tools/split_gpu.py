#!/usr/bin/env python3
"""Split the monolithic GPU sources into src/gpu/.

The split is done by exact line ranges out of the three original files, so the
body of every kernel and wrapper is carried over byte for byte. Only the
preambles (includes, macros, the two copies of gpuAssert) are rewritten, and
they are rewritten into shared headers.

Run from the directory that holds src/.
"""

import os
import re
import sys

SRC = 'src'
OUT = os.path.join(SRC, 'gpu')

CW = os.path.join(SRC, 'cuda_wrappers.cu')
GE = os.path.join(SRC, 'gpu_exp.cu')
TB = os.path.join(SRC, '3b_final.cc')

def read(p):
    with open(p) as f:
        return f.read().split('\n')

src = {CW: read(CW), GE: read(GE), TB: read(TB)}

# ---------------------------------------------------------------- the split
#
# (file, first_line, last_line, destination). Inclusive, 1-based. Every line of
# every source file appears exactly once below, either with a destination or
# with None meaning "preamble, re-expressed in a header". check_coverage()
# enforces that, so nothing can be dropped by accident.

PLAN = [
    # ---- cuda_wrappers.cu -------------------------------------------------
    (CW,    1,   29, None),                  # includes + tuning defines
    (CW,   30,  174, 'gpu_memory.cu'),       # device memory ledger
    (CW,  175,  214, None),                  # dead `counter`, gpuErrchk -> gpu_common.h
    (CW,  215,  221, 'gap_soap_radial.cu'),  # N_a
    (CW,  222,  271, 'gap_2b.cu'),           # gpu_spline, gpu_spline_der
    (CW,  272,  306, 'gpu_memory.cu'),       # vect_dble
    (CW,  307,  533, 'gpu_memory.cu'),       # malloc / free / memcpy wrappers
    (CW,  534,  553, 'gpu_blas.cu'),         # cublas handle lifetime
    (CW,  554,  570, 'gap_predict.cu'),      # gpu_pow
    (CW,  571,  604, 'gpu_blas.cu'),         # mmul_t_n, mvmul_n
    (CW,  605,  620, 'gap_predict.cu'),      # axpc
    (CW,  621,  693, 'gpu_memory.cu'),       # device enumeration and selection
    (CW,  694,  742, 'gap_predict.cu'),      # matvect_kernels, matvect_qs
    (CW,  743,  772, 'gpu_blas.cu'),         # mmul_n_t, dgemm_n_n
    (CW,  773,  964, 'gap_soap_forces.cu'),  # soap forces/virial, local props
    (CW,  965, 1452, 'gap_soap_descriptor.cu'),  # soap_p, soap_der, normalise
    (CW, 1453, 1960, 'gap_soap_angular.cu'), # cnk, plm, exp coefficients
    (CW, 1961, 1970, 'gap_soap_radial.cu'),  # check_nan
    (CW, 1971, 2498, 'gap_soap_radial.cu'),  # global scaling, poly3, poly3gauss
    (CW, 2499, 2649, 'gap_2b.cu'),           # 2b and the core potential
    (CW, 2650, 2684, 'gpu_memory.cu'),       # sync, pinned host memory
    # ---- gpu_exp.cu -------------------------------------------------------
    (GE,    1,   65, None),                  # includes, macros, gpuAssert
    (GE,   66,  389, 'gpu_scan.cu'),         # warp/block reduce, inclusive scan
    (GE,  390,  480, 'mad_electrostatics.cu'),
    (GE,  481,  486, 'gpu_scan.cu'),         # gpu_inclusive_scan_int
    (GE,  487,  860, 'mad_electrostatics.cu'),
    (GE,  861, 1474, 'mad_pdf.cu'),
    (GE, 1475, 1491, 'gpu_memory.cu'),       # gpu_meminfo
    (GE, 1492, 1702, 'mad_xrd.cu'),
    (GE, 1703, 1709, 'gpu_memory.cu'),       # gpu_print_pointer_*
    # ---- 3b_final.cc ------------------------------------------------------
    (TB,    1,   14, None),                  # includes
    (TB,   15,  886, 'gap_3b.cc'),
]

def check_coverage():
    ok = True
    for f, lines in src.items():
        covered = [0] * (len(lines) + 1)
        for pf, a, b, _ in PLAN:
            if pf != f:
                continue
            for i in range(a, b + 1):
                covered[i] += 1
        n = len(lines)
        # the final element of split('\n') is the empty string after the last \n
        for i in range(1, n):
            if covered[i] != 1:
                print(f'COVERAGE {f}:{i} counted {covered[i]}: {lines[i-1][:70]}')
                ok = False
    return ok

if not check_coverage():
    sys.exit('plan does not tile the sources exactly')

# ------------------------------------------------------------------ headers
#
# Per-file preamble. tpb differs between the two originals (64 in
# cuda_wrappers, 512 in gpu_exp) and every kernel that reads it was tuned
# against its own file's value, so it is set per destination from the file the
# code came from -- never unified.

COMMON = '''// Shared by every translation unit under src/gpu.
//
// This is the preamble the two monolithic sources used to carry a copy of
// each: the HIP/CUDA includes, and the error check. The two copies of
// gpuAssert had drifted -- cuda_wrappers.cu dumped a host backtrace, gpu_exp.cu
// did not -- and the backtrace is what lets addr2line name the Fortran call
// site behind a thin C shim, so that is the version kept here.
#ifndef TURBOGAP_GPU_COMMON_H
#define TURBOGAP_GPU_COMMON_H

#include "hip/hip_runtime.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <assert.h>
#include <hip/hip_complex.h>
#include <execinfo.h>

// Launch geometry shared across the MAD kernels and the scan primitives. Note
// that `tpb` is deliberately NOT here: the two original files disagreed on it
// (64 against 512) and every kernel was tuned against its own file's value, so
// each source under src/gpu sets it locally. These two never disagreed.
#define WARP_SIZE 32
#define BLOCK_SIZE 512 // Number of threads per block

#define gpuErrchk(ans)                                                                                                             \\
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char* file, int line, bool abort = true) {
  if (code != hipSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\\n", hipGetErrorString(code), file, line);
    // The wrapper file/line above only names the thin C shim; dump the native
    // stack so the offending Fortran call site can be identified with addr2line.
    void* bt_[64];
    int nbt_ = backtrace(bt_, 64);
    fprintf(stderr, "GPUassert: host backtrace (%d frames):\\n", nbt_);
    fflush(stderr);
    backtrace_symbols_fd(bt_, nbt_, fileno(stderr));
    fflush(stderr);
    if (abort)
      exit(code);
  }
}

// Report a launch geometry the device will refuse, and say which launch it was.
//
// CUDA rejects any launch with a zero extent as "invalid configuration
// argument" -- an error that names neither the kernel nor the dimension, and is
// raised at the next error check rather than at the launch. Both properties
// cost real time here, because the batched experimental paths compute their
// extents from per-batch, per-species-pair counts that go to zero as soon as
// gpu_n_batches > 1 splits the cell finely enough. This turns that into one
// line saying which kernel and what it asked for, immediately.
//
// Costs nothing when the launch is valid: three integer comparisons on the
// host, off the critical path of a kernel that is about to run.
#define TG_CHECK_LAUNCH(name_, grid_, block_, extent_name_, extent_)                                                               \\
  do {                                                                                                                             \\
    if ((grid_).x == 0 || (grid_).y == 0 || (grid_).z == 0 || (block_).x == 0) {                                                   \\
      fprintf(stderr,                                                                                                              \\
              "GPUlaunch: %s would launch grid(%u,%u,%u) block(%u,%u,%u) -- a zero extent.\\n"                                      \\
              "GPUlaunch:   %s = %ld. An empty batch reaches this kernel; guard the call.\\n",                                      \\
              (name_), (grid_).x, (grid_).y, (grid_).z, (block_).x, (block_).y, (block_).z, (extent_name_), (long) (extent_));     \\
      fflush(stderr);                                                                                                              \\
    }                                                                                                                              \\
  } while (0)

#endif // TURBOGAP_GPU_COMMON_H
'''

# What each generated .cu needs on top of gpu_common.h: its own tpb, the
# tuning macros its kernels read, and the headers declaring anything it calls
# that now lives in another translation unit.
PREAMBLE = {
'gpu_memory.cu': '''// Device memory: the allocation ledger, the malloc/free/memcpy wrappers, the
// pinned host allocator, device enumeration and selection, and the debug
// pointer printers. Used by both the GAP and the MAD paths.
#include "gpu_common.h"
#include "gpu_memory.h"

#include <hiprand/hiprand.h>
''',
'gpu_blas.cu': '''// hipBLAS/cuBLAS handle lifetime and the dense products the GAP prediction
// needs. Nothing here is specific to a descriptor.
#include "gpu_common.h"
#include "gpu_blas.h"
''',
'gpu_scan.cu': '''// Parallel primitives shared by the MAD and electrostatics neighbour-counting
// paths: a block reduction, a recursive reduction and an inclusive scan over
// per-pair flags.
//
// These used to sit at the top of gpu_exp.cu, where the pdf and the
// electrostatics entry points could launch kernel_multiply_flags directly.
// They no longer share a translation unit, and a <<<>>> launch cannot cross
// one without -rdc=true, so that kernel is now reached through the host
// launcher gpu_multiply_flags -- which computes the same geometry both call
// sites computed inline.
#include "gpu_common.h"
#include "gpu_scan.h"

#define tpb 512

//#define LOG_BLOCK_SIZE 10 // Log base 2 of BLOCK_SIZE
#define NUM_BANKS 32    // Define the number of shared memory banks
#define LOG_NUM_BANKS 5 // Logarithm base 2 of NUM_BANKS
//#define CONFLICT_FREE_OFFSET(index) ((index) >> LOG_NUM_BANKS)
#ifdef ZERO_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(n) ((n) >> (LOG_NUM_BANKS) + (n) >> (2 * LOG_NUM_BANKS))
#else
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)
#endif
''',
'gap_predict.cu': '''// GAP prediction arithmetic that is not tied to one descriptor: the kernel
// power, the energy offset, and the two mat-vec products over the sparse set.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
''',
'gap_soap_radial.cu': '''// The soap_turbo radial expansion: poly3 and poly3gauss basis coefficients and
// their derivatives, plus the global scaling applied to them.
//
// This is the single most expensive kernel family in a GAP run.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
#define mode_polynomial 1
''',
'gap_soap_angular.cu': '''// The soap_turbo angular expansion: the associated Legendre arrays, the
// exp(i m phi) coefficients and their derivatives, and the cnk contraction of
// the radial and angular parts.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
#define tpbcnk 64 // this is because k_max is 45???
''',
'gap_soap_descriptor.cu': '''// The SOAP descriptor itself, formed from cnk: the power spectrum, its
// normalisation, the radial/azimuthal/polar derivatives and their conversion
// to Cartesian derivatives, plus the transposes those need.
#include "gpu_common.h"
#include "gap_gpu.h"

#include <cstddef>

#define tpb 64
#define tpb_get_soap_der_one 128
#define tpbcnk 64 // this is because k_max is 45???

//For tiled transpositions
constexpr std::size_t TRANSPOSE_TILE_DIM = 32;
constexpr std::size_t TRANSPOSE_BLOCK_ROWS = 8;
''',
'gap_soap_forces.cu': '''// Contraction of the SOAP Cartesian derivatives into forces and the virial,
// and the same contraction for local properties.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
''',
'gap_2b.cu': '''// The two-body descriptor and the core (repulsive) potential, both of which
// evaluate a cubic spline over pair distances.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
''',
'mad_pdf.cu': '''// MAD: the pair distribution function.
//
// Counting the pairs inside the pdf cutoff, compacting them into the k-indexed
// buffers, the kernel density estimate of g(r) and its derivative, and the
// reduction over pairs.
#include "gpu_common.h"
#include "gpu_scan.h"
#include "mad_gpu.h"

#define tpb 512
''',
'mad_xrd.cu': '''// MAD: the X-ray/neutron structure factor.
//
// Gk and its derivative with respect to atomic positions, the Hadamard product
// with the scattering factors, the dgemv that turns that into per-pair forces,
// and the collection of those into forces and the virial.
#include "gpu_common.h"
#include "mad_gpu.h"

#define tpb 512
''',
'mad_electrostatics.cu': '''// Electrostatics: the shifted-force Coulomb sum (GSF).
//
// Grouped with the MAD sources rather than the GAP ones because it is reached
// from the same experimental-prediction driver and shares their pair
// compaction, but it is a separate physical model from pdf and xrd.
#include "gpu_common.h"
#include "gpu_scan.h"
#include "mad_gpu.h"

#define tpb 512
''',
'gap_3b.cc': '''// The three-body descriptor.
//
// Kept as .cc, and so compiled with $(CC) rather than $(CU), because it needs
// -std=c++20 for <numbers> and <bit>.
#include "gpu_common.h"
#include "gap_gpu.h"

#include <numbers>
#include <iostream>
#include <iomanip>
#include <bit>
#include <chrono>
''',
}

# --------------------------------------------------------------- assemble
os.makedirs(OUT, exist_ok=True)

bodies = {}
for f, a, b, dest in PLAN:
    if dest is None:
        continue
    bodies.setdefault(dest, []).extend(src[f][a - 1:b])

for dest, lines in bodies.items():
    with open(os.path.join(OUT, dest), 'w') as fh:
        fh.write(PREAMBLE[dest] + '\n' + '\n'.join(lines).rstrip('\n') + '\n')

with open(os.path.join(OUT, 'gpu_common.h'), 'w') as fh:
    fh.write(COMMON)

# ------------------------------------------------- declaration headers
#
# Generated from the definitions rather than typed, so a header can only ever
# disagree with its source by a change made after this ran -- and then it
# disagrees at compile time, which is the point of having it.

def signatures(dest, want):
    """Pull `extern "C"` signatures out of a generated file, in order.

    Two spellings appear: `extern "C" void f(...) {` per definition, and a
    single `extern "C" { ... }` block wrapping several plain definitions (3b
    uses the latter).
    """
    path = os.path.join(OUT, dest)
    lines = open(path).read().split('\n')
    out = []
    i = 0
    in_block = False
    while i < len(lines):
        line = lines[i]
        if line.rstrip() == 'extern "C" {':
            in_block = True
            i += 1
            continue
        if in_block and line.startswith('}'):
            in_block = False
            i += 1
            continue
        prefix = 'extern "C" ' if line.startswith('extern "C" ') else ('' if in_block else None)
        if prefix is None or not re.match(r'(extern "C" )?(void|int|double|bool|size_t|int\*|double\*) ', line):
            i += 1
            continue
        buf = [line]
        while not buf[-1].rstrip().endswith('{'):
            i += 1
            if i >= len(lines):
                break
            buf.append(lines[i])
        # A couple of definitions carry whole-line comments between the closing
        # paren and the brace. Carried into a declaration they would swallow
        # the semicolon, so drop them.
        buf = [b for b in buf if not b.lstrip().startswith('//')]
        sig = '\n'.join(buf).rstrip()
        if not sig.endswith('{'):
            i += 1
            continue
        sig = sig[:-1].rstrip() + ';'
        if prefix:
            sig = sig[len(prefix):]
        name = re.search(r'([A-Za-z_][A-Za-z_0-9]*)\s*\(', sig)
        if name and (want is None or name.group(1) in want):
            out.append(sig)
        i += 1
    return out

HEADERS = {
 'gpu_memory.h': ('''// Device and pinned-host memory, and device selection.
//
// Every one of these is called from Fortran through bind(C) in
// src/fortran_cuda_interfaces.f90; the declarations exist so the C++ side
// cannot drift from the definitions in gpu_memory.cu without a compile error.
''', ['gpu_memory.cu'], None),
 'gpu_blas.h': ('''// hipBLAS handle lifetime and the dense products used by GAP prediction.
''', ['gpu_blas.cu'], None),
 'gap_gpu.h': ('''// The GAP side of the GPU backend: 2b, 3b and soap_turbo.
//
// Grouped by descriptor, in the order a prediction uses them --
//   radial + angular -> cnk -> descriptor -> forces
// with 2b, 3b and the prediction arithmetic alongside.
''', ['gap_predict.cu', 'gap_soap_radial.cu', 'gap_soap_angular.cu',
      'gap_soap_descriptor.cu', 'gap_soap_forces.cu', 'gap_2b.cu',
      'gap_3b.cc'], None),
 'mad_gpu.h': ('''// The MAD side of the GPU backend: the pair distribution function, the
// structure factor, and the shifted-force electrostatics that shares their
// pair-compaction machinery.
''', ['mad_pdf.cu', 'mad_xrd.cu', 'mad_electrostatics.cu'], None),
}

for hdr, (blurb, files, want) in HEADERS.items():
    guard = 'TURBOGAP_' + hdr.upper().replace('.', '_')
    body = [blurb, f'#ifndef {guard}', f'#define {guard}', '',
            '#include "gpu_common.h"', '', 'extern "C" {', '']
    for f in files:
        body.append(f'// ---- {f}')
        for s in signatures(f, want):
            body.append(s)
        body.append('')
    body += ['} // extern "C"', '', f'#endif // {guard}', '']
    with open(os.path.join(OUT, hdr), 'w') as fh:
        fh.write('\n'.join(body))


# ------------------------------------------------ gpu_scan.h, written by hand
#
# Not generated: these are the internal C++ helpers the two MAD entry points
# share, not bind(C) entry points, so there is no definition spelling to copy.

SCAN_H = '''// Parallel primitives over per-pair flag arrays, shared by the pdf and the
// electrostatics neighbour counting. Implemented in gpu_scan.cu.
#ifndef TURBOGAP_GPU_SCAN_H
#define TURBOGAP_GPU_SCAN_H

#include "gpu_common.h"

// Sum n ints into d_out, reducing recursively until one block remains.
void recursiveReduce(int* d_in, int* d_out, int n, hipStream_t* stream);

// In-place inclusive prefix sum over n ints.
void inclusiveScan(int* d_data_out, int n, hipStream_t* stream);

// Zero every prefix-sum entry whose flag is zero, so the result is a compaction
// index for the flagged pairs and nothing else.
//
// This is a host launcher rather than a kernel the callers launch themselves.
// The kernel lives in gpu_scan.cu, and a <<<>>> launch cannot cross a
// translation unit without -rdc=true; both call sites computed the same
// geometry inline, so it moved here with them.
void gpu_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream);

// Print the pending error on a stream without clearing the sticky one.
void gpu_peek_stream_error(hipStream_t* stream);

extern "C" {
void gpu_inclusive_scan_int(int size, int* n_neigh_index_d, hipStream_t* stream);
}

#endif // TURBOGAP_GPU_SCAN_H
'''
with open(os.path.join(OUT, 'gpu_scan.h'), 'w') as fh:
    fh.write(SCAN_H)

# The launcher itself, appended to gpu_scan.cu next to the kernel it launches.
LAUNCHER = '''
// See gpu_scan.h: the two callers used to launch kernel_multiply_flags
// directly, with exactly this geometry, when they shared a translation unit
// with it.
void gpu_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_multiply_flags<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, nk_flags_d, nk_sum_flags_d);
}
'''
with open(os.path.join(OUT, 'gpu_scan.cu'), 'a') as fh:
    fh.write(LAUNCHER)

# and the two call sites now go through it
OLD_CALL = '''  nblocks = dim3((n_pairs + tpb) / tpb, 1, 1);
  nthreads = dim3(tpb, 1, 1);
'''
NEW_CALL = ''
LAUNCH = '  kernel_multiply_flags<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, nk_flags_d, nk_flags_sum_d);\n'
REPL = '  gpu_multiply_flags(n_pairs, nk_flags_d, nk_flags_sum_d, stream);\n'
for dest in ('mad_pdf.cu', 'mad_electrostatics.cu'):
    p = os.path.join(OUT, dest)
    t = open(p).read()
    assert t.count(LAUNCH) == 1, (dest, t.count(LAUNCH))
    assert t.count(OLD_CALL) == 1, (dest, t.count(OLD_CALL))
    t = t.replace(OLD_CALL, NEW_CALL).replace(LAUNCH, REPL)
    open(p, 'w').write(t)

print('wrote', sorted(os.listdir(OUT)))
