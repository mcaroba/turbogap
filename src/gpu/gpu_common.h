// Shared by every translation unit under src/gpu.
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

#define gpuErrchk(ans)                                                                                                             \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char* file, int line, bool abort = true) {
  if (code != hipSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
    // The wrapper file/line above only names the thin C shim; dump the native
    // stack so the offending Fortran call site can be identified with addr2line.
    void* bt_[64];
    int nbt_ = backtrace(bt_, 64);
    fprintf(stderr, "GPUassert: host backtrace (%d frames):\n", nbt_);
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
#define TG_CHECK_LAUNCH(name_, grid_, block_, extent_name_, extent_)                                                               \
  do {                                                                                                                             \
    if ((grid_).x == 0 || (grid_).y == 0 || (grid_).z == 0 || (block_).x == 0) {                                                   \
      fprintf(stderr,                                                                                                              \
              "GPUlaunch: %s would launch grid(%u,%u,%u) block(%u,%u,%u) -- a zero extent.\n"                                      \
              "GPUlaunch:   %s = %ld. An empty batch reaches this kernel; guard the call.\n",                                      \
              (name_), (grid_).x, (grid_).y, (grid_).z, (block_).x, (block_).y, (block_).z, (extent_name_), (long) (extent_));     \
      fflush(stderr);                                                                                                              \
    }                                                                                                                              \
  } while (0)

#endif // TURBOGAP_GPU_COMMON_H
