// One kernel, two backends, chosen at compile time.
//
// Every kernel under src/gpu used to be a raw <<<>>> launch of a __global__
// function. That is one implementation of a parallel-for, welded to one
// vendor. This header adds a second -- Kokkos -- and lets the SAME kernel body
// serve both, so a kernel is written once and neither backend can drift from
// the other.
//
//     -D_KOKKOS absent   the launch is the <<<>>> it always was
//     -D_KOKKOS present  the body runs under Kokkos, on the same stream
//
// Why Kokkos rather than a second hand-written path: LAMMPS's KOKKOS package is
// Kokkos, and a pair_style that calls TurboGAP has to hand it LAMMPS's
// neighbour data in LAMMPS's memory. Expressing the kernels in Kokkos is what
// makes that possible later without writing them a third time. See
// docs/LAMMPS_KOKKOS.md.
//
// THE STREAM IS NOT REPLACED. Kokkos adopts the stream gpu_context already
// owns, so a Kokkos kernel and a <<<>>> kernel issued next to each other are
// still ordered with respect to one another. That is what lets the port
// proceed file by file instead of all at once: an unconverted kernel keeps
// working beside a converted one.
//
// WHAT IS PRESERVED, AND WHY IT MATTERS. A backend swap must not move a
// single digit, or the regression suite cannot tell a port from a bug. Two
// properties give that:
//
//   - the converted primitives are INTEGER (pair counts, prefix sums, bucket
//     indices). Integer addition is associative, so a different summation
//     order is the same answer, exactly.
//   - the converted floating-point kernels sum SEQUENTIALLY INSIDE ONE
//     THREAD, in an order fixed by the neighbour list. gpu_scatter.cu was
//     written that way on purpose (see its header); this header only changes
//     how threads are handed their indices, never the order a thread adds in.
//
// So the converted primitives return the same values under either backend,
// exactly. Do not convert a kernel whose reduction Kokkos would reassociate --
// a parallel_reduce over doubles -- without saying so at the call site.
//
// That does NOT make whole-output bit-equality a usable test, because the
// device binary is not reproducible against ITSELF: run estat_gsf twice with
// one unchanged CUDA binary and energy_soap moves in the tenth digit. The
// check that does work is the difference of two failure sets -- Kokkos vs
// CUDA, against CUDA vs CUDA on the same cases. Anything failing in the first
// and not the second is the port's doing. tools/verify_kokkos.sh runs both.
#ifndef TURBOGAP_GPU_BACKEND_H
#define TURBOGAP_GPU_BACKEND_H

#include "gpu_common.h"

#ifdef _KOKKOS
#include <Kokkos_Core.hpp>
#endif

// The annotation a kernel body carries. Both spellings need nvcc's
// --expt-extended-lambda, which the Makefile adds for a GPU build.
#ifdef _KOKKOS
#define TG_LAMBDA KOKKOS_LAMBDA
#else
#define TG_LAMBDA [=] __device__
#endif

// The annotation a device HELPER function carries -- a spline evaluation, an
// integer power -- as opposed to a kernel body. Same reasoning as TG_LAMBDA:
// __device__ is CUDA's spelling and does not exist for a Kokkos host backend.
#ifdef _KOKKOS
#define TG_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define TG_INLINE_FUNCTION __device__ __forceinline__
#endif

// An atomic add usable from a body that compiles under either backend.
//
// atomicAdd would in fact compile both ways today, because the Kokkos backend
// here is CUDA and the compiler is still nvcc. It is spelled through a macro
// anyway so that a body does not silently stop being portable the day Kokkos
// is built against a host backend, where atomicAdd does not exist.
#ifdef _KOKKOS
#define TG_ATOMIC_ADD(addr, val) Kokkos::atomic_add((addr), (val))
#else
#define TG_ATOMIC_ADD(addr, val) atomicAdd((addr), (val))
#endif

// The same, where the caller needs the value the address held BEFORE the add
// -- handing out slots in a counting sort, say. These are two different
// functions in Kokkos: atomic_add returns void and only atomic_fetch_add
// returns the old value, whereas CUDA's atomicAdd always returns it. Using
// TG_ATOMIC_ADD for a slot would compile under CUDA and fail under Kokkos.
#ifdef _KOKKOS
#define TG_ATOMIC_FETCH_ADD(addr, val) Kokkos::atomic_fetch_add((addr), (val))
#else
#define TG_ATOMIC_FETCH_ADD(addr, val) atomicAdd((addr), (val))
#endif

#ifdef _KOKKOS

// Named to match LAMMPS's LMPDeviceType, which a pair_style templates on.
#if defined(HOP_TARGET_CUDA) || defined(CUDA)
using TGDeviceType = Kokkos::Cuda;
#else
using TGDeviceType = Kokkos::HIP;
#endif

// Kokkos must come up on the card THIS RANK selected, not on device 0.
// cuda_set_device hands out cards round-robin (my_rank % num_gpus), so a
// default-initialised Kokkos would put every rank's kernels on one device
// while the rest of the code used another -- which is not an error, just
// silently wrong results and a serialised run. Called from cuda_set_device,
// once, after the card is chosen.
inline void tg_backend_initialize(int device_id) {
  if (Kokkos::is_initialized() || Kokkos::is_finalized())
    return;
  Kokkos::InitializationSettings settings;
  settings.set_device_id(device_id);
  settings.set_disable_warnings(true);
  Kokkos::initialize(settings);
  // A backstop only, for a path that never reaches gpu_context_finalize. The
  // ordered teardown is tg_backend_finalize below; by the time this runs the
  // CUDA context is usually already gone.
  std::atexit([]() {
    if (Kokkos::is_initialized() && !Kokkos::is_finalized())
      Kokkos::finalize();
  });
}

// Kokkos must be torn down BEFORE the CUDA context it holds handles into.
//
// Leaving this to atexit does not work: gpu_context_finalize calls
// cuda_device_reset, which destroys the context, and Kokkos::finalize then
// fails freeing its own pinned staging buffer --
//
//   cudaFreeHost(Kokkos::Impl::CudaInternal::constantMemHostStagingPerDevice)
//   error( cudaErrorInvalidValue )
//
// thrown as a std::runtime_error out of an exit handler, which is an abort and
// a SIGABRT backtrace on a run that had already finished its work correctly.
// Called from cuda_device_reset, ahead of the reset.
inline void tg_backend_finalize() {
  if (Kokkos::is_initialized() && !Kokkos::is_finalized())
    Kokkos::finalize();
}

// An execution space instance wrapping a stream the Fortran context owns.
//
// Constructing one is not free -- it queries the stream's context through the
// driver API -- and the pdf and electrostatics paths launch dozens of kernels
// per step, so they are cached by stream handle. The cache LEAKS its entries
// deliberately: a Kokkos::Cuda destructor running at static-destruction time,
// after Kokkos::finalize, is undefined. Eight handles is the whole cost, and
// TurboGAP uses two.
inline TGDeviceType& tg_exec(hipStream_t* stream) {
  constexpr int NCACHE = 8;
  static hipStream_t keys[NCACHE];
  static TGDeviceType* vals[NCACHE];
  static int n_cached = 0;

  hipStream_t s = (stream == nullptr) ? (hipStream_t) 0 : stream[0];
  for (int i = 0; i < n_cached; i++)
    if (keys[i] == s)
      return *vals[i];

  if (n_cached == NCACHE) {
    // Not fatal: fall back to the first entry's space, which is still a valid
    // place to run, and say so once rather than growing without bound.
    static bool warned = false;
    if (!warned) {
      fprintf(stderr, "GPUbackend: more than %d streams; reusing the first. Raise NCACHE.\n", NCACHE);
      warned = true;
    }
    return *vals[0];
  }

  keys[n_cached] = s;
  vals[n_cached] = new TGDeviceType(s);
  return *vals[n_cached++];
}

#else // plain CUDA/HIP

inline void tg_backend_initialize(int) {
}
inline void tg_backend_finalize() {
}

// The generic launcher the <<<>>> path needs: one instantiation per kernel
// body, which is what a hand-written __global__ per kernel was anyway.
template <class Functor> __global__ void tg_range_kernel(int n, Functor f) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    f(i);
  }
}

template <class Functor> __global__ void tg_range_kernel_long(long long n, Functor f) {
  long long i = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    f(i);
  }
}

// The 3-D form keeps the geometry the hand-written launches used: the fast
// index tiled across a block, the other two carried on the grid's y and z.
template <class Functor> __global__ void tg_range_kernel_3d(int n0, int n1, int n2, Functor f) {
  int i0 = blockIdx.x * blockDim.x + threadIdx.x;
  int i1 = blockIdx.y;
  int i2 = blockIdx.z;
  if (i0 < n0 && i1 < n1 && i2 < n2) {
    f(i0, i1, i2);
  }
}

#endif // _KOKKOS

// ---------------------------------------------------------------- parallel_for
//
// block is the <<<>>> block size, ignored by Kokkos, which picks its own.
// Pass the value the file used before conversion so an unconverted and a
// converted build launch identical geometry.
template <class Functor>
inline void tg_parallel_for(const char* name, int n, hipStream_t* stream, const Functor& f, int block = 256) {
  if (n <= 0)
    return;
#ifdef _KOKKOS
  (void) block;
  Kokkos::parallel_for(name, Kokkos::RangePolicy<TGDeviceType>(tg_exec(stream), 0, n), f);
#else
  (void) name;
  dim3 nblocks((n + block - 1) / block, 1, 1);
  dim3 nthreads(block, 1, 1);
  tg_range_kernel<<<nblocks, nthreads, 0, stream[0]>>>(n, f);
#endif
}

// The same, over a range that does not fit in an int.
//
// The structure-factor kernels index n_k * n_samples, a product of two counts
// each of which is comfortably an int and whose product need not be. They were
// written with a `long long` thread index for that reason, and converting them
// to the int form above would reintroduce exactly the overflow they avoid.
template <class Functor>
inline void tg_parallel_for_long(const char* name, long long n, hipStream_t* stream, const Functor& f, int block = 256) {
  if (n <= 0)
    return;
#ifdef _KOKKOS
  (void) block;
  Kokkos::parallel_for(name, Kokkos::RangePolicy<TGDeviceType, Kokkos::IndexType<long long>>(tg_exec(stream), 0, n), f);
#else
  (void) name;
  long long n_blocks = (n + block - 1) / block;
  dim3 nblocks((unsigned int) n_blocks, 1, 1);
  dim3 nthreads(block, 1, 1);
  tg_range_kernel_long<<<nblocks, nthreads, 0, stream[0]>>>(n, f);
#endif
}

// ------------------------------------------------------------------ 3-D range
//
// For a kernel whose work is naturally three indices -- site, radial n, angular
// k -- rather than one. Flattening it into a single index by hand would work
// and would read badly, and Kokkos has MDRangePolicy for exactly this.
//
// ONLY for a pure map: Kokkos tiles an MDRange differently from the y/z grid
// the <<<>>> path uses, so the order the triples are visited is not the same.
// That is immaterial when each triple writes its own output and reads nothing
// another triple writes, and wrong the moment it does not.
template <class Functor>
inline void tg_parallel_for_3d(const char* name, int n0, int n1, int n2, hipStream_t* stream, const Functor& f, int block = 256) {
  if (n0 <= 0 || n1 <= 0 || n2 <= 0)
    return;
#ifdef _KOKKOS
  (void) block;
  Kokkos::parallel_for(name, Kokkos::MDRangePolicy<TGDeviceType, Kokkos::Rank<3>>(tg_exec(stream), {0, 0, 0}, {n0, n1, n2}), f);
#else
  (void) name;
  dim3 nblocks((n0 + block - 1) / block, n1, n2);
  dim3 nthreads(block, 1, 1);
  tg_range_kernel_3d<<<nblocks, nthreads, 0, stream[0]>>>(n0, n1, n2, f);
#endif
}

#ifdef _KOKKOS

// ------------------------------------------------------------- inclusive scan
//
// In place over n ints, the result landing back in data. Integers, so this is
// exact whichever backend runs it.
//
// Only the Kokkos backend defines these. The <<<>>> path keeps the hand-written
// scan and recursive reduction that gpu_scan.cu has always carried, reached
// through the same public names -- so no caller learns which backend it got.
inline void tg_inclusive_scan_int(const char* name, int* data, int n, hipStream_t* stream) {
  if (n <= 0)
    return;
  Kokkos::parallel_scan(
      name, Kokkos::RangePolicy<TGDeviceType>(tg_exec(stream), 0, n), KOKKOS_LAMBDA(const int i, int& running, const bool final) {
        // Read before the write: Kokkos runs the body twice
        // for the same i, and only the final pass stores, so
        // data[i] is still the input value both times.
        const int v = data[i];
        running += v;
        if (final)
          data[i] = running;
      });
}

// Sum n ints into out_d, which is a DEVICE address: the callers keep the total
// on the card and copy it back themselves when they need it.
inline void tg_reduce_sum_int(const char* name, const int* in, int n, int* out_d, hipStream_t* stream) {
  if (n <= 0)
    return;
  Kokkos::View<int, TGDeviceType::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>> result(out_d);
  Kokkos::parallel_reduce(
      name, Kokkos::RangePolicy<TGDeviceType>(tg_exec(stream), 0, n),
      KOKKOS_LAMBDA(const int i, int& partial) { partial += in[i]; }, result);
}

#endif // _KOKKOS

#endif // TURBOGAP_GPU_BACKEND_H
