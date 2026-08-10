// Device memory: the allocation ledger, the malloc/free/memcpy wrappers, the
// pinned host allocator, device enumeration and selection, and the debug
// pointer printers. Used by both the GAP and the MAD paths.
#include "gpu_common.h"
#include "gpu_memory.h"

#include <hiprand/hiprand.h>

//  Device memory ledger
// ==========================================================================
//
// Why a ledger and not just hipMemGetInfo
// ---------------------------------------
// hipMemGetInfo answers "how much does the DEVICE have left", which includes
// every other process on the card and the driver's own reservations. It cannot
// answer "how much has TurboGAP asked for" -- the number you need when deciding
// whether to raise gpu_n_batches -- and it cannot answer it after the fact,
// because by the time an allocation fails the interesting quantity is the
// high-water mark, not the current one.
//
// The Fortran side has total_gpu_memory() already, but it is hand maintained:
// each call site adds its own byte count, only the experimental paths call it,
// and nothing ever checks it against the device. This sits underneath every
// wrapper instead, so it cannot disagree with what was actually allocated.
//
// Cost is one hash insert or erase per allocation, around 100 ns against the
// 2-3 us hipMallocAsync itself takes. Sampling hipMemGetInfo per allocation
// would have been far more expensive and less informative.
//
// Thread safety: OPENMP=1 runs the batched loops on several threads, each with
// its own stream, all allocating. One mutex, never held across a CUDA call.
#include <mutex>
#include <unordered_map>

namespace {

std::unordered_map<void*, size_t> tg_mem_live_map;
std::mutex tg_mem_mutex;
size_t tg_mem_live = 0;      // bytes currently held by wrapper allocations
size_t tg_mem_peak = 0;      // high-water mark of the above
size_t tg_mem_requested = 0; // cumulative, including memory since freed
size_t tg_mem_n_alloc = 0;
size_t tg_mem_n_free = 0;
size_t tg_mem_n_untracked_free = 0;

// What the input asked for, so a failure can name the keyword that fixes it.
// Registered from Fortran; negative means "never told".
double tg_mem_hint_max_gb = -1.0;
int tg_mem_hint_n_batches = -1;

const double TG_GB = 1024.0 * 1024.0 * 1024.0;

void tg_mem_note_alloc(void* p, size_t n) {
  if (p == nullptr)
    return;
  std::lock_guard<std::mutex> lock(tg_mem_mutex);
  tg_mem_live_map[p] = n;
  tg_mem_live += n;
  tg_mem_requested += n;
  tg_mem_n_alloc += 1;
  if (tg_mem_live > tg_mem_peak)
    tg_mem_peak = tg_mem_live;
}

void tg_mem_note_free(void* p) {
  if (p == nullptr)
    return;
  std::lock_guard<std::mutex> lock(tg_mem_mutex);
  auto it = tg_mem_live_map.find(p);
  if (it == tg_mem_live_map.end()) {
    // A pointer this ledger never saw. Counted rather than warned about: the
    // scan/reduce helpers in gpu_exp.cu and the vendored soap_turbo code
    // allocate directly, so this is expected and only interesting in bulk.
    tg_mem_n_untracked_free += 1;
    return;
  }
  tg_mem_live -= it->second;
  tg_mem_n_free += 1;
  tg_mem_live_map.erase(it);
}

} // namespace

// The whole picture. Called on an allocation failure, and from Fortran via
// gpu_memory_report at phase boundaries.
extern "C" void gpu_mem_report(const char* label) {
  size_t dev_free = 0, dev_total = 0;
  hipMemGetInfo(&dev_free, &dev_total);
  size_t live, peak, requested, n_alloc, n_free, n_untracked, n_map;
  {
    std::lock_guard<std::mutex> lock(tg_mem_mutex);
    live = tg_mem_live;
    peak = tg_mem_peak;
    requested = tg_mem_requested;
    n_alloc = tg_mem_n_alloc;
    n_free = tg_mem_n_free;
    n_untracked = tg_mem_n_untracked_free;
    n_map = tg_mem_live_map.size();
  }
  fprintf(stderr, "GPUmem [%s]\n", label ? label : "");
  fprintf(stderr, "GPUmem   device      %8.3f GB free of %8.3f GB total (%.1f%% in use)\n", dev_free / TG_GB, dev_total / TG_GB,
          dev_total ? 100.0 * (double) (dev_total - dev_free) / (double) dev_total : 0.0);
  fprintf(stderr, "GPUmem   this rank   %8.3f GB live in %lu buffers, peak %8.3f GB\n", live / TG_GB, (unsigned long) n_map,
          peak / TG_GB);
  fprintf(stderr, "GPUmem   cumulative  %8.3f GB requested over %lu allocations, %lu frees (%lu untracked)\n", requested / TG_GB,
          (unsigned long) n_alloc, (unsigned long) n_free, (unsigned long) n_untracked);
  fflush(stderr);
}

// Numbers for the Fortran side, in bytes. Any pointer may be null.
extern "C" void gpu_mem_stats(size_t* dev_free, size_t* dev_total, size_t* live, size_t* peak) {
  size_t f = 0, t = 0;
  hipMemGetInfo(&f, &t);
  if (dev_free)
    *dev_free = f;
  if (dev_total)
    *dev_total = t;
  std::lock_guard<std::mutex> lock(tg_mem_mutex);
  if (live)
    *live = tg_mem_live;
  if (peak)
    *peak = tg_mem_peak;
}

// Let the input's batching keywords appear in a failure message.
extern "C" void gpu_mem_set_hint(double max_gb, int n_batches) {
  std::lock_guard<std::mutex> lock(tg_mem_mutex);
  tg_mem_hint_max_gb = max_gb;
  tg_mem_hint_n_batches = n_batches;
}

// What to say when there is no memory left.
//
// The generic gpuErrchk path prints "out of memory", the wrapper's own file and
// line, and a host backtrace -- none of which says how much was wanted, how
// much there was, or which input keyword controls it. On a million-atom cell
// that is the difference between a five-minute fix and an afternoon.
static void tg_mem_report_oom(size_t wanted, hipError_t code) {
  fprintf(stderr, "\nGPUmem: device allocation of %.3f GB FAILED: %s\n", wanted / TG_GB, hipGetErrorString(code));
  gpu_mem_report("at the failure");
  fprintf(stderr, "GPUmem   The batching keywords are what fix this:\n");
  if (tg_mem_hint_max_gb > 0.0)
    fprintf(stderr, "GPUmem     max_Gbytes_per_process = %.3f  <- splits the SOAP loop; HALVE it\n", tg_mem_hint_max_gb);
  else
    fprintf(stderr, "GPUmem     max_Gbytes_per_process       <- splits the SOAP loop (default 1.0); lower it\n");
  if (tg_mem_hint_n_batches > 0)
    fprintf(stderr, "GPUmem     gpu_n_batches          = %d     <- splits the pdf/xrd/estat loops; DOUBLE it\n",
            tg_mem_hint_n_batches);
  else
    fprintf(stderr, "GPUmem     gpu_n_batches                <- splits the pdf/xrd/estat loops (default 1); raise it\n");
  fprintf(stderr, "GPUmem     gpu_mem_fraction             <- sets both from the device automatically\n");
  fflush(stderr);
}
__global__ void vect_dble(double* a, int N) {
  int idx = threadIdx.x + blockIdx.x * gridDim.x;
  if (idx < N)
    printf(" %lf \n", a[idx]);
}


// Stream-ordered device allocation.
//
// This used to end with
//
//     hipError_t err;
//     hipDeviceSynchronize();
//     err = hipGetLastError();      // and nothing read err
//
// which was the single most expensive statement in the program. Every device
// buffer in TurboGAP is allocated through here, so a whole-device barrier ran
// on every allocation: nsys measured 3982 cudaDeviceSynchronize calls costing
// 60.5% of the CO_predict span, 1838 at 31.0% on XRD_mad, 1025 at 49.7% on
// estat_gsf. It bought nothing -- the error it collected was assigned to a
// local and discarded, with the test commented out -- and gpuErrchk on the line
// above already checks the only status hipMallocAsync reports.
//
// It was also a correctness ceiling, not just a cost. hipDeviceSynchronize
// waits for every stream, so no two streams could ever overlap while it was
// here, whatever the batched loops did with them: measured peak concurrency
// was 1 on all four test cases, including the two that create more than one
// stream. Removing it is the precondition for the OpenMP/stream work, not an
// independent optimisation.
//
// What still makes this safe: hipMallocAsync is stream-ORDERED. The returned
// pointer is valid for work enqueued after it on `stream`, which is the only
// way any caller here uses it. Handing such a pointer to a DIFFERENT stream
// does need an explicit event, and that is a constraint on new stream code, not
// something this barrier was correctly providing.
extern "C" void gpu_malloc_async(void** a_d, size_t Np, hipStream_t* stream) {
  hipError_t code = hipMallocAsync((void**) a_d, Np, stream[0]);

  if (code == hipErrorOutOfMemory) {
    // Retry once, and this is not optimism.
    //
    // hipFreeAsync is STREAM-ORDERED: it returns a buffer to the pool when the
    // stream reaches the free, not when the call is made. A loop that allocates
    // and frees per batch therefore has, at any instant, a pile of memory that
    // is logically free and physically still held. Draining the stream converts
    // it, and trimming the pool hands anything the pool is hoarding back to the
    // driver so a large contiguous request can be met.
    //
    // Both are expensive, which is exactly why they are on the failure path
    // only: the normal path is one hipMallocAsync and nothing else. This is the
    // one place a synchronise buys something -- unlike the whole-device barrier
    // that used to run here on every successful allocation.
    hipStreamSynchronize(stream[0]);
    // The CURRENT device, not device 0. cuda_set_device assigns cards
    // round-robin, so on a multi-GPU node a hardcoded 0 would drain some other
    // rank's pool and leave this one exactly as full as it was.
    int dev = 0;
    hipMemPool_t pool;
    if (hipGetDevice(&dev) == hipSuccess && hipDeviceGetDefaultMemPool(&pool, dev) == hipSuccess)
      hipMemPoolTrimTo(pool, 0);
    code = hipMallocAsync((void**) a_d, Np, stream[0]);
    if (code == hipSuccess)
      fprintf(stderr,
              "GPUmem: note -- a %.3f GB allocation needed a pool drain to succeed. The run is\n"
              "GPUmem:         close to the device limit; raise gpu_n_batches or lower\n"
              "GPUmem:         max_Gbytes_per_process before it stops fitting at all.\n",
              Np / (1024.0 * 1024.0 * 1024.0));
  }

  if (code != hipSuccess) {
    tg_mem_report_oom(Np, code);
    gpuAssert(code, __FILE__, __LINE__);
  }

  tg_mem_note_alloc(*a_d, Np);
  return;
}

//  extern "C" void gpu_malloc_async(void **a_d, size_t Np, hipStream_t *stream )
// {


//   //gpuErrchk(hipMallocAsync((void **) a_d,  Np ,stream[0]));
//   gpuErrchk(hipMalloc((void **) a_d,  Np ));
//   hipError_t err;
//   hipDeviceSynchronize();
//   err = hipGetLastError();
// //  if (err != hipSuccess) {
// //}
//    return;
// }

extern "C" void cuda_memset_async(void* a_d, int value, size_t Np, hipStream_t* stream) {
  gpuErrchk(hipMemsetAsync(a_d, value, Np, stream[0]));
}
extern "C" void cuda_malloc_all_blocking(void** a_d, size_t Np) {
  hipError_t code = hipMalloc((void**) a_d, Np);
  if (code != hipSuccess) {
    tg_mem_report_oom(Np, code);
    gpuAssert(code, __FILE__, __LINE__);
  }
  tg_mem_note_alloc(*a_d, Np);
  return;
}
extern "C" void cuda_device_reset() {
  hipDeviceReset();
}
extern "C" void cuda_free(void** a_d) {
  tg_mem_note_free(*a_d);
  gpuErrchk(hipFree(*a_d));
  //printf("GPU memory freed \n");
  //hipDeviceReset();
  return;
}


// extern "C" void cuda_free_async(void **a_d, hipStream_t *stream )
// {
//   //gpuErrchk(hipFreeAsync(*a_d, stream[0]));
//   gpuErrchk(hipFree(*a_d));
//    return;
// }

extern "C" void cuda_free_async(void** a_d, hipStream_t* stream) {
  // Noted before the call, not after: hipFreeAsync is stream-ordered, so the
  // memory is not back in the pool yet, but the pointer is dead to this code
  // the moment the call is made and must leave the ledger with it. The ledger
  // therefore tracks what TurboGAP believes it holds, which is the quantity the
  // batching keywords control -- hipMemGetInfo, reported next to it, is what
  // the device actually has.
  tg_mem_note_free(*a_d);
  gpuErrchk(hipFreeAsync(*a_d, stream[0]));
  return;
}

/*extern "C" void GPU_fill_rand(double *A, int N, int ccc) {
// Create a pseudo-random number generator
hiprandGenerator_t prng;
hiprandCreateGenerator(&prng, HIPRAND_RNG_PSEUDO_DEFAULT);

// Set the seed for the random number generator using the system clock
hiprandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock() + (unsigned long long)  ccc * N);

// Fill the array with random numbers on the device
hiprandGenerateUniformDouble(prng, A,N);
//vect_dble<<<(N+128-1)/128,128>>>(A,N);
//hipDeviceSynchronize();
printf("\n Filled \n");
}
 */

extern "C" void cuda_cpy_htod(void* a, void* a_d, size_t N, hipStream_t* stream) {
  gpuErrchk(hipMemcpyAsync(a_d, a, N, hipMemcpyHostToDevice, stream[0]));
  //gpuErrchk(hipMemcpy(a_d, a, N, hipMemcpyHostToDevice));
  return;
}
extern "C" void cuda_cpy_htod_blocking(void* a, void* a_d, size_t N) {
  ;
  gpuErrchk(hipMemcpy(a_d, a, N, hipMemcpyHostToDevice));
  return;
}

extern "C" void cuda_cpy_dtod(void* b_d, void* a_d, size_t N, hipStream_t* stream) {
  gpuErrchk(hipMemcpyAsync(a_d, b_d, N, hipMemcpyDeviceToDevice, stream[0]));
  //gpuErrchk(hipMemcpy( a_d, b_d, N, hipMemcpyDeviceToDevice));
  return;
}

extern "C" void cuda_cpy_dtoh(void* a_d, void* a, size_t N, hipStream_t* stream) {
  gpuErrchk(hipMemcpyAsync(a, a_d, N, hipMemcpyDeviceToHost, stream[0]));
  //gpuErrchk(hipMemcpy(a, a_d,  N, hipMemcpyDeviceToHost));
  return;
}

extern "C" void cuda_cpy_dtoh_event(void* a_d, void* a, size_t N, hipStream_t* stream) {
  // Create a CUDA event
  hipEvent_t copyComplete;
  hipEventCreate(&copyComplete);

  gpuErrchk(hipMemcpyAsync(a, a_d, N, hipMemcpyDeviceToHost, stream[0]));

  // Record the event after the asynchronous copy
  hipEventRecord(copyComplete);

  // Wait for the event to complete
  hipEventSynchronize(copyComplete);

  // Clean up
  hipEventDestroy(copyComplete);

  //gpuErrchk(hipMemcpy(a, a_d,  N, hipMemcpyDeviceToHost));
  return;
}


extern "C" void cuda_cpy_dtoh_blocking(void* a_d, void* a, size_t N) {
  gpuErrchk(hipMemcpy(a, a_d, N, hipMemcpyDeviceToHost));
  return;
}


/* extern "C" void cuda_cpy_double_htod(double *a, double *a_d, int N)
   {
   gpuErrchk(hipMemcpy(a_d, a, sizeof(double) * N, hipMemcpyHostToDevice));
   return;
   } */
/*
   extern "C" void cuda_cpy_bool_htod(bool *a, double *a_d, int N)
   {
   gpuErrchk(hipMemcpy(a_d, a, sizeof(bool) * N, hipMemcpyHostToDevice));
//gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(bool) * N, hipMemcpyHostToDevice));
return;
} */
/*
   extern "C" void cuda_cpy_bool_dtoh(bool *a_d, bool *a, int N)
   {
   gpuErrchk(hipMemcpy(a, a_d, sizeof(bool) * N, hipMemcpyDeviceToHost));
   return;
   } */

/* extern "C" void cuda_cpy_double_complex_htod(hipDoubleComplex *a, hipDoubleComplex *a_d, int N)
   {

   gpuErrchk(hipMemcpy(a_d, a, sizeof(hipDoubleComplex) * N, hipMemcpyHostToDevice));
//gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(hipDoubleComplex) * N, hipMemcpyHostToDevice));
return;
} */


/*
   extern "C" void cuda_cpy_double_complex_dtoh(hipDoubleComplex *a_d, hipDoubleComplex *a ,int N)
   {
//hipMemcpyAsync( a, a_d, sizeof(double) * N, hipMemcpyDeviceToHost );
gpuErrchk(hipMemcpy( a, a_d, sizeof(hipDoubleComplex) * N, hipMemcpyDeviceToHost ));
//printf("\nTest cpy D to H \n");

return;
} */

/* extern "C" void cuda_cpy_int_htod(int *a, int *a_d, int N)
   {

   gpuErrchk(hipMemcpy(a_d, a, sizeof(int) * N, hipMemcpyHostToDevice ));
//gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(int) * N, hipMemcpyHostToDevice ));
return;
}
 */
/*
   extern "C" void cuda_cpy_double_dtoh(double *a_d, double *a ,int N)
   {
   gpuErrchk(hipMemcpy( a, a_d, sizeof(double) * N, hipMemcpyDeviceToHost ));

   return;
   } */
/*
   extern "C" void cuda_cpy_double_dtod(double *b_d, double *a_d,int N)
   {
   gpuErrchk(hipMemcpy( a_d, b_d, sizeof(double) * N, hipMemcpyDeviceToDevice ));
   return;
   }
 */


// extern "C" void wrappers_all(double *soap, double *kernels, double *kernels_copy, double *Qs, double *energies, double delta, double zeta, double e0, int n_sites, int n_soap, int n_sparse, int size_kernels, int size_soap, int size_Qs, int size_alphas, int  size_energies)
// {
//   int ntpb=256;
//   int nblocks=(size_kernels+ntpb-1)/ntpb;
//   // Create a handle for CUBLAS
// 	hipblasHandle_t handle;
// 	hipblasCreate(&handle);
//   double *kernels_d, *kernels_copy_d, *soap_d, *Qs_d, *energies_d;
//   hipMalloc( &kernels_d, sizeof(double) * size_kernels );
//   hipMalloc( &kernels_copy_d, sizeof(double) * size_kernels );
//   hipMalloc( &soap_d, sizeof(double) * size_soap );
//   hipMalloc( &Qs_d, sizeof(double) * size_Qs );
//   hipMalloc( &energies_d, sizeof(double)*size_energies);


//   const double alf = 1;
//   const double bet = 0;

//   hipMemcpy(kernels_d, kernels, sizeof(double) * size_kernels, hipMemcpyHostToDevice );
//   hipMemcpy(soap_d, soap, sizeof(double) * size_soap, hipMemcpyHostToDevice );
//   hipMemcpy(Qs_d, Qs, sizeof(double) * size_Qs, hipMemcpyHostToDevice );
//   // Do the actual multiplication

//   hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d, n_sites);
// //hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N,  nBy, nAx, nAy, alpha, B, nAy, A, nAy, beta, C, nBy);
//     //printf("\n hipblasDgemm \n");
//   // gpu_blas_mmul_t_n(cubhandle,     A,     B,      C,         nAx,      nAy,       nBy,             bb, zeta, N)
//   // gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites, kernels_copy_d, zeta, size_kernels)

//   hipMemcpy( kernels, kernels_d, sizeof(double) * size_kernels, hipMemcpyDeviceToHost );
//   hipLaunchKernelGGL(gpu_pow, dim3(nblocks,1,1), dim3(ntpb,1,1), 0, 0, kernels_d,kernels_copy_d, zeta, size_kernels);
//   hipMemcpy( kernels_copy, kernels_copy_d, sizeof(double) * size_kernels, hipMemcpyDeviceToHost );
// 	// Destroy the handle
// 	hipblasDestroy(handle);
//   hipFree(kernels_d);
//   hipFree(kernels_copy_d);
//   hipFree(soap_d);
//   hipFree(Qs_d);
//   hipFree(energies_d);
//   //printf("\n %d %d %d %d %d %d %d %d  \n", n_sites, n_soap, n_sparse, size_kernels,  size_soap,  size_Qs,  size_alphas,  size_energies);
//   //printf("\n %d %d %d\n", nblocks,ntpb, size_kernels);
//   //exit(0);
//  return;
// }

// How many devices this process can see.
//
// Needed to divide the memory budget correctly. cuda_set_device below hands out
// devices round-robin (my_rank % num_gpus), so the number of ranks sharing one
// card is ceil(ntasks / num_gpus) -- not ntasks, and not one. Getting this
// wrong in either direction is silent: too small a divisor and every rank sizes
// its batches to memory another rank is already using; too large and the run is
// needlessly slow.
extern "C" int gpu_device_count() {
  int n = 0;
  if (hipGetDeviceCount(&n) != hipSuccess)
    return 1;
  return n > 0 ? n : 1;
}

extern "C" void cuda_set_device(int my_rank) {
  int num_gpus = 0;
  int mygpuid;
  gpuErrchk(hipGetDeviceCount(&num_gpus));
  gpuErrchk(hipSetDevice(my_rank % num_gpus));
  gpuErrchk(hipGetDevice(&mygpuid));
  /*gpuErrchk(hipSetDevice(0));*/
  //  printf("\n Seta Aset at %d %d %d %d\n", num_gpus, my_rank%num_gpus,my_rank, mygpuid);
  //exit(0);
  return;
}

extern "C" void gpu_device_sync() {
  gpuErrchk(hipDeviceSynchronize());
}

extern "C" void gpu_stream_sync(hipStream_t* stream) {
  gpuErrchk(hipStreamSynchronize(stream[0]));
}


extern "C" int* host_alloc(size_t size) {
  int* a;
  size_t fm, gm;
  hipMemGetInfo(&fm, &gm);
  printf("Host GPU alloc:  memory usage: %lu/%lu MB\n", fm / 1024 / 1024, gm / 1024 / 1024);
  hipHostMalloc((void**) &(a), size);
  printf("%s\n", hipGetErrorString(hipGetLastError()));
  return a;
}

extern "C" void host_free(void* ptr) {
  hipHostFree(ptr);
  printf("%s\n", hipGetErrorString(hipGetLastError()));
}

extern "C" void cpy_htoh_pinned(void* src, void* dest, size_t size) {
  gpuErrchk(hipMemcpy(dest, src, size, hipMemcpyHostToHost));
  printf("%s\n", hipGetErrorString(hipGetLastError()));
}


extern "C" void gpu_check_error() {
  hipError_t code = hipDeviceSynchronize();
  printf("\n %s \n", hipGetErrorString(code));
  gpuErrchk(code);
}
extern "C" void gpu_meminfo() {
  // Initialize variables to store memory info
  size_t freeMem, totalMem;

  // Get memory information
  hipMemGetInfo(&freeMem, &totalMem);

  // Calculate used memory
  size_t usedMem = totalMem - freeMem;

  // Print out the memory information
  printf("--- Total Memory: %lu bytes\n", totalMem);
  printf("--- Free Memory: %lu bytes\n", freeMem);
  printf("--- Used Memory: %lu bytes\n", usedMem);
}


extern "C" void gpu_print_pointer_int(int* p) {
  printf(" address:  %p \n", p);
}
extern "C" void gpu_print_pointer_double(double* p) {
  printf(" address:  %p \n", p);
}
