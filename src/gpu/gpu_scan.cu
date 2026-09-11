// Parallel primitives shared by the MAD and electrostatics neighbour-counting
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


// Warp reduction function to sum values within a warp
__inline__ __device__ int warpReduceSum(int val) {
  // Use shuffle down to reduce across the warp
  for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
#ifdef CUDA
    val += __shfl_down_sync(0xffffffff, val, offset, warpSize);
#else
    val += __shfl_down(val, offset, warpSize);
#endif
  }
  return val;
}


__inline__ __device__ double warpReduceSumDouble(double val) {
  // Use shuffle down to reduce across the warp
  for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
#ifdef CUDA
    val += __shfl_down_sync(0xffffffff, val, offset, warpSize);
#else
    val += __shfl_down(val, offset, warpSize);
#endif
  }
  return val;
}


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

// // Kernel for inclusive scan with shared memory padding to avoid bank conflicts
// __global__ void inclusive_scan_kernel_diff(int *d_in, int *d_out, int n, int block_size) {
//     // Dynamically allocate shared memory with padding to avoid bank conflicts
//     extern __shared__ int temp[];

//     int tid = threadIdx.x;
//     int gid = blockIdx.x * block_size + tid;

//     // Calculate the effective index in padded shared memory
//     int warp_size = 32; // Warp size for CUDA architectures
//     int effective_index = tid + tid / warp_size; // Add padding every warp_size elements

//     // Load input into shared memory
//     if (gid < n) {
//         temp[effective_index] = d_in[gid];
//     } else {
//         temp[effective_index] = 0; // Pad with zero if outside array bounds
//     }

//     __syncthreads();

//     // Perform the scan
//     for (int offset = 1; offset < block_size; offset *= 2) {
//         int val = 0;
//         if (tid >= offset) {
//             val = temp[effective_index - offset];
//         }
//         __syncthreads(); // Synchronize before updating
//         temp[effective_index] += val;
//         __syncthreads(); // Synchronize after updating
//     }

//     // Write result to output array
//     if (gid < n) {
//         d_out[gid] = temp[effective_index]; // Inclusive: No adjustment needed
//     }
// }

// void inclusive_scan_diff(int *h_in, int *h_out, int n, int block_size) {
//     int *d_in, *d_out;
//     size_t size = n * sizeof(int);

//     // Allocate device memory
//     cudaMalloc(&d_in, size);
//     cudaMalloc(&d_out, size);

//     // Copy input to device
//     cudaMemcpy(d_in, h_in, size, cudaMemcpyHostToDevice);

//     // Launch kernel
//     int threads_per_block = block_size;
//     int blocks_per_grid = (n + threads_per_block - 1) / threads_per_block;

//     // Calculate shared memory size with padding
//     int warp_size = 32;
//     int shared_memory_size = (block_size + block_size / warp_size) * sizeof(int);

//     inclusive_scan_kernel_diff<<<blocks_per_grid, threads_per_block, shared_memory_size>>>(d_in, d_out, n, block_size);

//     // Copy result back to host
//     cudaMemcpy(h_out, d_out, size, cudaMemcpyDeviceToHost);

//     // Free device memory
//     cudaFree(d_in);
//     cudaFree(d_out);
// }

// int main() {
//     // Example input
//     int h_in[] = {1, 2, 3, 4, 5};
//     int n = sizeof(h_in) / sizeof(h_in[0]);
//     int h_out[n];

//     // Set block size (can be any power of 2, up to maximum threads per block)
//     int block_size = 128;

//     // Perform inclusive scan
//     inclusive_scan(h_in, h_out, n, block_size);

//     // Print result
//     std::cout << "Input: ";
//     for (int i = 0; i < n; i++) {
//         std::cout << h_in[i] << " ";
//     }
//     std::cout << "\nOutput: ";
//     for (int i = 0; i < n; i++) {
//         std::cout << h_out[i] << " ";
//     }
//     std::cout << std::endl;

//     return 0;
// }


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

  // Copy result back to host
  gpuErrchk(hipMemcpyAsync(d_data_out, d_data + 1, (n - 1) * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));

  // Free device memory
  gpuErrchk(hipFreeAsync(d_data, stream[0]));
  gpuErrchk(hipFreeAsync(d_blockSums, stream[0]));
}

__global__ void kernel_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < n_pairs) {
    nk_sum_flags_d[tid] *= nk_flags_d[tid];
  }
}

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
// with it.
void gpu_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_multiply_flags<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, nk_flags_d, nk_sum_flags_d);
}
