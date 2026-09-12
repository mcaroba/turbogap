// Parallel primitives shared by the MAD and electrostatics neighbour-counting
// paths: a reduction and an inclusive scan over per-pair flags.
//
// These used to sit at the top of gpu_exp.cu, where the pdf and the
// electrostatics entry points could launch kernel_multiply_flags directly.
// They no longer share a translation unit, and a <<<>>> launch cannot cross
// one without -rdc=true, so that kernel is now reached through the host
// launcher gpu_multiply_flags -- which computes the same geometry both call
// sites computed inline.
//
// TWO BACKENDS. Under -D_KOKKOS the reduction and the scan are Kokkos
// primitives; otherwise they are the hand-written block reduction and Blelloch
// scan below. The public names are the same either way, so mad_pdf.cu and
// mad_electrostatics.cu do not know which they got.
//
// The two agree EXACTLY, not approximately: both sum ints, and integer
// addition is associative, so no reassociation Kokkos performs can move a
// result. That is the property that makes the bit-exact regression suite a
// real check on this port rather than a tolerance to be tuned. See
// src/gpu/gpu_backend.h.
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gpu_scan.h"

#define tpb 512

#define NUM_BANKS 32    // Define the number of shared memory banks
#define LOG_NUM_BANKS 5 // Logarithm base 2 of NUM_BANKS
#ifdef ZERO_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(n) ((n) >> (LOG_NUM_BANKS) + (n) >> (2 * LOG_NUM_BANKS))
#else
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)
#endif

#ifdef _KOKKOS

// Kokkos owns the recursion and the block geometry; these are the same two
// primitives, asked for by name.
void recursiveReduce(int* d_in, int* d_out, int n, hipStream_t* stream) {
  tg_reduce_sum_int("turbogap_reduce_flags", d_in, n, d_out, stream);
}

void inclusiveScan(int* d_data_out, int n, hipStream_t* stream) {
  tg_inclusive_scan_int("turbogap_scan_flags", d_data_out, n, stream);
}

#else // the hand-written CUDA/HIP primitives

//------------------------------------------------------------//
//-------------------   Reduction Kernel   -------------------//
//------------------------------------------------------------//
__global__ void blockReduceKernel(int* d_in, int* d_out, int n) {
  // Shared memory for partial results
  __shared__ int sharedData[BLOCK_SIZE];

  // Calculate global thread index
  int tid = threadIdx.x;
  int globalIndex = blockIdx.x * blockDim.x + tid;

  // Load data into shared memory
  sharedData[tid] = (globalIndex < n) ? d_in[globalIndex] : 0;

  // Perform reduction in shared memory
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    __syncthreads();
    if (tid < s) {
      sharedData[tid] += sharedData[tid + s];
    }
  }

  // Write the result of the block to global memory
  if (tid == 0) {
    d_out[blockIdx.x] = sharedData[0];
  }
}

// Recursive function to perform reduction
void recursiveReduce(int* d_in, int* d_out, int n, hipStream_t* stream) {
  int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

  // If only one block remains, no further recursion is needed
  if (numBlocks == 1) {
    blockReduceKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_in, d_out, n);
    //        hipDeviceSynchronize();
  } else {
    // Allocate memory for intermediate results
    int* d_intermediate;
    //        gpuErrchk(hipMalloc(&d_intermediate, numBlocks * sizeof(int)));
    gpuErrchk(hipMallocAsync((void**) &d_intermediate, numBlocks * sizeof(int), stream[0]));

    // Perform the reduction on the blocks
    blockReduceKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_in, d_intermediate, n);
    //hipDeviceSynchronize();

    // Recursively reduce the intermediate results
    recursiveReduce(d_intermediate, d_out, numBlocks, stream);

    // Free intermediate memory
    gpuErrchk(hipFreeAsync(d_intermediate, stream[0]));
  }
}


//------------------------------------------------------------//
//---------------------   Scan Kernel   ----------------------//
//------------------------------------------------------------//
// Kernel for the inclusive scan using shared memory with dynamic padding
__global__ void inclusiveScanKernel(int* d_data, int* d_blockSums, int n) {
  // Shared memory with padding to avoid bank conflicts
  __shared__ int sharedData[BLOCK_SIZE * 2 + BLOCK_SIZE / NUM_BANKS];

  // Calculate thread and global indices
  int tid = threadIdx.x;
  int globalIndexA = 2 * blockIdx.x * blockDim.x + tid;
  int globalIndexB = globalIndexA + blockDim.x;

  // Calculate bank offsets to avoid conflicts
  int bankOffsetA = tid + (tid >> LOG_NUM_BANKS);                             // Bank offset for the first element
  int bankOffsetB = bankOffsetA + blockDim.x + (blockDim.x >> LOG_NUM_BANKS); // Bank offset for the second element

  // Load data into shared memory with padding
  sharedData[bankOffsetA] = (globalIndexA < n) ? d_data[globalIndexA] : 0;
  sharedData[bankOffsetB] = (globalIndexB < n) ? d_data[globalIndexB] : 0;

  int ai;
  int bi;

  // Up-sweep (reduce) phase
  int offset = 1;
  for (int d = blockDim.x; d > 0; d >>= 1) {
    __syncthreads();
    if (tid < d) {
      ai = (offset * (2 * tid + 1)) - 1;
      bi = (offset * (2 * tid + 2)) - 1;

      // Adjust indices for bank conflicts
      ai += ai >> LOG_NUM_BANKS;
      bi += bi >> LOG_NUM_BANKS;

      sharedData[bi] += sharedData[ai];
    }
    offset *= 2;
  }

  // Clear the last element for the down-sweep phase
  if (tid == 0) {
    int lastIndex = (2 * blockDim.x - 1) + ((2 * blockDim.x - 1) >> LOG_NUM_BANKS);
    d_blockSums[blockIdx.x] = sharedData[lastIndex];
    sharedData[lastIndex] = 0;
  }

  // Down-sweep phase
  for (int d = 1; d < 2 * blockDim.x; d *= 2) {
    offset >>= 1;
    __syncthreads();
    if (tid < d) {
      ai = (offset * (2 * tid + 1)) - 1;
      bi = (offset * (2 * tid + 2)) - 1;

      // Adjust indices for bank conflicts
      ai += ai >> LOG_NUM_BANKS;
      bi += bi >> LOG_NUM_BANKS;

      int temp = sharedData[ai];
      sharedData[ai] = sharedData[bi];
      sharedData[bi] += temp;
    }
  }

  __syncthreads();

  // Write the results back to global memory
  if (globalIndexA < n)
    d_data[globalIndexA] = sharedData[bankOffsetA];
  if (globalIndexB < n)
    d_data[globalIndexB] = sharedData[bankOffsetB];
}

// Kernel to add block sums to each element
// out[n-1] = out[n-2] + input[n-1]. See the note in inclusiveScan.
__global__ void fixLastInclusiveKernel(int* d_data_out, int n) {
  if (threadIdx.x == 0 && blockIdx.x == 0 && n >= 2) {
    d_data_out[n - 1] += d_data_out[n - 2];
  }
}

__global__ void addBlockSumsKernel(int* d_data, int* d_blockSums, int n) {
  int globalIndex = threadIdx.x + blockIdx.x * blockDim.x * 2;
  if (blockIdx.x > 0) {
    int blockSum = d_blockSums[blockIdx.x - 1];
    if (globalIndex < n)
      d_data[globalIndex] += blockSum;
    if (globalIndex + blockDim.x < n)
      d_data[globalIndex + blockDim.x] += blockSum;
  }
}

// Function to perform an inclusive scan on an array
void inclusiveScan(int* d_data_out, int n, hipStream_t* stream) {
  // Calculate the size needed for padding
  int paddedN = ((n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2)) * (BLOCK_SIZE * 2); //+ 1;
  size_t size = paddedN * sizeof(int);

  // Allocate memory on the device
  int* d_data;
  int* d_blockSums;
  gpuErrchk(hipMallocAsync(&d_data, size, stream[0]));
  gpuErrchk(hipMallocAsync(&d_blockSums, ((n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2)) * sizeof(int), stream[0]));

  // Copy data into buffer
  gpuErrchk(hipMemcpyAsync(d_data, d_data_out, n * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));
  gpuErrchk(hipMemsetAsync(d_data + n, 0, (paddedN - n) * sizeof(int), stream[0])); // Zero out padding

  // Calculate the number of blocks needed
  int numBlocks = (n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2);

  //  printf("\n In recursive scan, numBlocks = %d", numBlocks);


  // Launch kernel for the main scan operation


  inclusiveScanKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_data, d_blockSums, n);

  // If there are multiple blocks, perform a scan on the block sums
  if (numBlocks > 1) {
    inclusiveScan(d_blockSums, numBlocks, stream);
    //    hipStre();
    // Add block sums to each element in the array
    addBlockSumsKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_data, d_blockSums, n);
  }

  // Shifted by one, which turns the kernel's exclusive scan into an inclusive
  // one -- for every element but the last, whose value the kernel never wrote.
  gpuErrchk(hipMemcpyAsync(d_data_out, d_data + 1, (n - 1) * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));

  // d_data_out[n-1] still holds the caller's own last input at this point, so
  // adding the element before it -- now the last exclusive value -- completes
  // the scan. n == 1 needs nothing: one element is its own inclusive scan.
  fixLastInclusiveKernel<<<1, 1, 0, stream[0]>>>(d_data_out, n);

  // Free device memory
  gpuErrchk(hipFreeAsync(d_data, stream[0]));
  gpuErrchk(hipFreeAsync(d_blockSums, stream[0]));
}


#endif // _KOKKOS

void gpu_peek_stream_error(hipStream_t* stream) {
  hipError_t code = hipDeviceSynchronize();
  printf("\n %s \n", hipGetErrorString(code));
  gpuErrchk(code);

  gpuErrchk(hipStreamSynchronize(stream[0]));
  gpuErrchk(hipPeekAtLastError());
}

extern "C" void gpu_inclusive_scan_int(int size, int* n_neigh_index_d, hipStream_t* stream) {
  //  printf("doing gpu inclusive scan, size = %d ", size  ) ;
  inclusiveScan(n_neigh_index_d, size, stream);
  //  gpu_peek_stream_error( stream );
}

// See gpu_scan.h: the two callers used to launch kernel_multiply_flags
// directly, with exactly this geometry, when they shared a translation unit
// with it. One body now, dispatched to whichever backend is compiled in.
void gpu_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_multiply_flags", n_pairs, stream, TG_LAMBDA(const int tid) { nk_sum_flags_d[tid] *= nk_flags_d[tid]; }, tpb);
}
