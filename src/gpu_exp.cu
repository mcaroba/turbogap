#include "hip/hip_runtime.h"
// wrappers file
// compile with:
// rm cuda_wrappers.o; nvcc -lcublas -lcurand -arch=sm_80 src/cuda_wrappers.cu -c;
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <cmath>

#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
//#include <hipsolver.h>
#include <hiprand/hiprand.h>
#include <assert.h>
#include <hip/hip_complex.h>

#define tpb 256 // optimize for best performance & check the effect on each kernel which used tpb for the shared memory
#define tpbcnk 64 // this is because k_max is 45???

// #ifdef CUDA
// #define WARP_SIZE 32
// #else
// #define WARP_SIZE 64
// #endif

#define WARP_SIZE 32

#define BLOCK_SIZE 512 // Number of threads per block
//#define LOG_BLOCK_SIZE 10 // Log base 2 of BLOCK_SIZE
#define NUM_BANKS 32      // Define the number of shared memory banks
#define LOG_NUM_BANKS 5   // Logarithm base 2 of NUM_BANKS
//#define CONFLICT_FREE_OFFSET(index) ((index) >> LOG_NUM_BANKS)
#ifdef ZERO_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(n) ((n) >> (LOG_NUM_BANKS) + (n) >> (2 * LOG_NUM_BANKS))
#else
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)
#endif

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Now define the kernels for the pair distribution function and xrd

// __device__ double warp_red(double data) {

//    double res = data;
//    for (int i = 32; i!=0; i=i>>1) {
//       res += __shfl_down(res, i, 64);
//    }
//    return res;
// }


// __device__ double warp_red(double data) {

//    double res = data;
//    for (int i =warpSize/2; i!=0; i=i>>1) {
//       res += __shfl_down_sync(0xffffffff,res, i,warpSize);
//    }
//    return res;
// }

//  __device__ void warpReduce(volatile int *sdata, unsigned int tid) 
//   {

// #ifdef CUDA
// #define WARP_SIZE 32
// #else
// #define WARP_SIZE 64
// #endif
    

//     sdata[tid] = sdata[tid] + sdata[tid+16];
//     sdata[tid] = sdata[tid] + sdata[tid+8];
//     sdata[tid] = sdata[tid] + sdata[tid+4];
//     sdata[tid] = sdata[tid] + sdata[tid+2];
//     sdata[tid] = sdata[tid] + sdata[tid+1];
      
//   }  

// // Function to perform reduction within a warp
// __device__ int warpReduceSum(int val) {
//     // Reduce across a warp using shuffle instructions
//     for (int offset = WARP_SIZE / 2; offset > 0; offset /= 2) {
//         val += __shfl_down(val, offset);
//     }
//     return val;
// }


// Warp reduction function to sum values within a warp
__inline__ __device__ int warpReduceSum(int val) {
    // Use shuffle down to reduce across the warp
    for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
      val += __shfl_down_sync(0xffffffff,val, offset,warpSize);
    }
    return val;
}


__inline__ __device__ double warpReduceSumDouble(double val) {
    // Use shuffle down to reduce across the warp
    for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
      val += __shfl_down_sync(0xffffffff,val, offset,warpSize);
    }
    return val;
}



// Kernel to perform a block-wide reduction using warpReduceSum
__global__ void blockReduceSum(int* d_in, int* d_out, int n) {
    // Shared memory for block-level reduction
    __shared__ int sharedData[WARP_SIZE];

    // Calculate global thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int lane = threadIdx.x % WARP_SIZE; // Lane index within the warp
    int warpId = threadIdx.x / WARP_SIZE; // Warp index within the block

    // Load data from global memory and perform warp reduction
    int val = (tid < n) ? d_in[tid] : 0; // Check bounds
    val = warpReduceSum(val);

    // Store reduced value of each warp in shared memory
    if (lane == 0) {
        sharedData[warpId] = val;
    }

    // Synchronize to ensure all warp reductions are done
    __syncthreads();

    // Use first warp to reduce values in shared memory
    if (warpId == 0) {
        val = (threadIdx.x < blockDim.x / WARP_SIZE) ? sharedData[lane] : 0;
        val = warpReduceSum(val);
        // The final sum is in sharedData[0] and saved to global memory by the first thread
        if (lane == 0) {
            d_out[blockIdx.x] = val;
        }
    }
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
void recursiveReduce(int* d_in, int* d_out, int n) {
    int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // If only one block remains, no further recursion is needed
    if (numBlocks == 1) {
        blockReduceKernel<<<numBlocks, BLOCK_SIZE>>>(d_in, d_out, n);
        hipDeviceSynchronize();
    } else {
        // Allocate memory for intermediate results
        int* d_intermediate;
        hipMalloc(&d_intermediate, numBlocks * sizeof(int));

        // Perform the reduction on the blocks
        blockReduceKernel<<<numBlocks, BLOCK_SIZE>>>(d_in, d_intermediate, n);
        hipDeviceSynchronize();

        // Recursively reduce the intermediate results
        recursiveReduce(d_intermediate, d_out, numBlocks);

        // Free intermediate memory
        hipFree(d_intermediate);
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
    int bankOffsetA = tid + (tid >> LOG_NUM_BANKS);  // Bank offset for the first element
    int bankOffsetB = bankOffsetA + blockDim.x + (blockDim.x >> LOG_NUM_BANKS); // Bank offset for the second element

    // Load data into shared memory with padding
    sharedData[bankOffsetA] = (globalIndexA < n) ? d_data[globalIndexA] : 0;
    sharedData[bankOffsetB] = (globalIndexB < n) ? d_data[globalIndexB] : 0;

    int ai;
    int bi;
    // [  0 ,2,  3 ,4 ,5 ,6  ]
    // [  0  0  2   5  9 ]
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
      // BLOCK_SIZE * 2 + BLOCK_SIZE / NUM_BANKS
      //        int lastIndex = (2 * blockDim.x - 1) + ((2 * blockDim.x - 1) >> LOG_NUM_BANKS);
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
    if (globalIndexA < n) d_data[globalIndexA] = sharedData[bankOffsetA];
    if (globalIndexB < n) d_data[globalIndexB] = sharedData[bankOffsetB];
}

// Kernel to add block sums to each element
__global__ void addBlockSumsKernel(int* d_data, int* d_blockSums, int n) {
    int globalIndex = threadIdx.x + blockIdx.x * blockDim.x * 2;
    if (blockIdx.x > 0) {
        int blockSum = d_blockSums[blockIdx.x - 1];
	//	printf("\n globalIndex = %d, globalIndex+dim = %d, blockSum = %d\n", globalIndex, globalIndex+blockDim.x, blockSum);
        if (globalIndex < n) d_data[globalIndex] += blockSum;
        if (globalIndex + blockDim.x < n) d_data[globalIndex + blockDim.x] += blockSum;
    }
}

// Function to perform an inclusive scan on an array
void inclusiveScan(int* d_data_out, int n) {
    // Calculate the size needed for padding
    int paddedN = ((n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2)) * (BLOCK_SIZE * 2) + 1;
    size_t size = paddedN * sizeof(int);

    // Allocate memory on the device
    int* d_data;
    int* d_blockSums;
    hipMalloc(&d_data, size);
    hipMalloc(&d_blockSums, ((n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2)) * sizeof(int));

    // Copy data from host to device
    hipMemcpy(d_data, d_data_out, n * sizeof(int), hipMemcpyDeviceToDevice);
    hipMemset(d_data + n, 0, (paddedN - n) * sizeof(int));  // Zero out padding

    // Calculate the number of blocks needed
    int numBlocks = (n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2);

    printf("\n In recursive scan, numBlocks = %d", numBlocks);
    
    // Launch kernel for the main scan operation
    inclusiveScanKernel<<<numBlocks, BLOCK_SIZE>>>(d_data, d_blockSums, n);
    hipDeviceSynchronize();

    // If there are multiple blocks, perform a scan on the block sums
    if (numBlocks > 1) {
      inclusiveScan(d_blockSums, numBlocks);
      //      inclusiveScanKernel<<<1, BLOCK_SIZE>>>(d_blockSums, d_blockSums, numBlocks);
      hipDeviceSynchronize();

      //--- Sum block sums on host for ease of implementation, see above comment
      // int* h_block_sums;
      // gpuErrchk(hipHostMalloc((void**) &h_block_sums, numBlocks*sizeof(int)));
      // gpuErrchk(hipMemcpy( h_block_sums, d_blockSums, numBlocks*sizeof(int), hipMemcpyDeviceToHost));

      // for( int i = 1; i < numBlocks; ++i ){
      // 	h_block_sums[i] += h_block_sums[i-1];
      // }
      // gpuErrchk(hipMemcpy( d_blockSums, h_block_sums, numBlocks*sizeof(int), hipMemcpyHostToDevice));
    
      // gpuErrchk(hipHostFree(h_block_sums));

      // Add block sums to each element in the array
      addBlockSumsKernel<<<numBlocks, BLOCK_SIZE>>>(d_data, d_blockSums, n);
      hipDeviceSynchronize();
    }

    // Copy result back to host
    //    int* h_data;
    //    hipHostMalloc(&h_data, n * sizeof(int));    
    hipMemcpy(d_data_out, d_data+1, n * sizeof(int), hipMemcpyDeviceToDevice);    
    //    h_data[n-1] = d_data[n-1] + h_data[n-2];
    //    hipMemcpy(d_data_out, h_data, n * sizeof(int), hipMemcpyHostToDevice);

    // Free device memory
    hipFree(d_data);
    hipFree(d_blockSums);
    //    hipHostFree( h_data );
}




  
// __global__ void prescan(float *g_idata, float *g_odata, int n)
// {
//   extern __shared__ float temp[2*BLOCK_SIZE + (2*BLOCK_SIZE) / NUM_BANKS];// allocated on invocation
//   int thid = threadIdx.x;
//   int offset = 1;

//   int ai = thid;
//   int bi = thid + (n/2);
//   int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
//   int bankOffsetB = CONFLICT_FREE_OFFSET(bi);
//   temp[ai + bankOffsetA] = g_idata[ai];
//   temp[bi + bankOffsetB] = g_idata[bi];
//   // temp[2*thid] = g_idata[2*thid]; // load input into shared memory
//   // temp[2*thid+1] = g_idata[2*thid+1];
//   for (int d = n>>1; d > 0; d >>= 1) // build sum in place up the tree
//     {
//       __syncthreads();
//       if (thid < d)
// 	{
// 	  int ai = offset*(2*thid+1)-1;
// 	  int bi = offset*(2*thid+2)-1;
// 	  ai += CONFLICT_FREE_OFFSET(ai);
// 	  bi += CONFLICT_FREE_OFFSET(bi);
// 	  temp[bi] += temp[ai];
// 	}
//       offset *= 2;
//     }
//   //  if (thid == 0) { temp[n - 1] = 0; } // clear the last element
//   if (thid==0) { temp[n – 1 + CONFLICT_FREE_OFFSET(n - 1)] = 0; }
//   for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
//     {
//       offset >>= 1;
//       __syncthreads();
//       if (thid < d)
// 	{
// 	  int ai = offset*(2*thid+1)-1;
// 	  int bi = offset*(2*thid+2)-1;
// 	  ai += CONFLICT_FREE_OFFSET(ai);
// 	  bi += CONFLICT_FREE_OFFSET(bi);
// 	  float t = temp[ai];
// 	  temp[ai] = temp[bi];
// 	  temp[bi] += t;
// 	}
//     }
//   __syncthreads();

//   g_odata[ai] = temp[ai + bankOffsetA];
//   g_odata[bi] = temp[bi + bankOffsetB];
// }





















//------------------------------------------------------------//
//----------------------   Setup PDF   -----------------------//
//------------------------------------------------------------//

__global__
void kernel_get_pair_distribution_nk(int i_beg, int i_end, int n_sites0, int* neighbors_list_d,
				     int* n_neigh_d, int* neighbor_species_d,
				     int* species_d, double* rjs_d, double* xyz,
				     double r_min, double r_max, double r_cut, double buffer,
				     int* nk_flags_d, int sp1, int sp2){

  int i_site=i_beg-1+threadIdx.x+blockIdx.x*blockDim.x;
  int k_val=threadIdx.x+blockIdx.x*blockDim.x;  
  int i,j,k,s,k1,k2;
  double r;
  int tid = threadIdx.x;
  int lane = tid % WARP_SIZE;
  int warpId = tid / WARP_SIZE;
  
  // nk_array_d[blockIdx.x] = 0;
  
  // __shared__ int nk_shared[tpb];
  
  //  printf(" In exp count kernel \n");
  int nk_loc = 0;
  if(i_site<i_end){
    if( species_d[i_site] == sp1 ||  species_d[i_site] == sp2 ) {
      k=0;
      for (i=i_beg-1; i<i_site; i++) k += n_neigh_d[i];
      for (j=1; j<n_neigh_d[i_site]; j++) {
        k +=1;
        if( !(
	      (species_d[i_site]==sp1 && neighbor_species_d[k]==sp2) ||
	      (species_d[i_site]==sp2 && neighbor_species_d[k]==sp1))   ) continue;
        r = rjs_d[k];
	if (r < 0.001 || r > r_cut || r < r_min || r > r_max + buffer ) continue;
	nk_loc += 1;
	nk_flags_d[k] = 1;
      }
    }
  }
    
}


extern "C" void  gpu_get_pair_distribution_nk(int i_beg, int i_end, int n_pairs,  int n_sites0,  int* neighbors_list,
						  int* n_neigh, int* neighbor_species,
					      int* species, double* rjs, double* xyz, double r_min, double r_max, double r_cut,
						  double buffer, 
					      int* nk_out_d, int* nk_flags_d, int* nk_flags_sum_d,  int species_1, int species_2,
						  hipStream_t *stream ){

  // This function is to set the k_index array for the partial pair distributions 

  
  dim3 nblocks=dim3((i_end - i_beg + tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);


  kernel_get_pair_distribution_nk<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, n_sites0, neighbors_list,
								      n_neigh, neighbor_species,
								      species, rjs, xyz, r_min, r_max, r_cut, buffer,
								       nk_flags_d, species_1, species_2);

  hipStreamSynchronize(stream[0]);

  // Perform recursive reduction
  recursiveReduce(nk_flags_d, nk_out_d, n_pairs);
  gpuErrchk(hipStreamSynchronize(stream[0]));
  

  printf("\n Starting recursive scan  \n");  
  gpuErrchk(hipMemcpy( nk_flags_sum_d, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToDevice));

  // int* nk_flags_sum_temp_d;
  // hipMalloc( &nk_flags_sum_temp_d, n_pairs * sizeof( int ) );
  // gpuErrchk(hipMemcpy( nk_flags_sum_temp_d, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToDevice));    
  
  // // Perform recursive inclusive scan
  // inclusiveScan(nk_flags_sum_temp_d, nk_flags_sum_d, n_pairs);

  // hipFree(nk_flags_sum_temp_d);


  inclusiveScan(nk_flags_sum_d, n_pairs);
  
  printf("\n Finished recursive scan  \n");  

  // Print results for debugging of the scan
  //---------------------------------------------//
  // --- Check that the host can do the same --- //

  // int* nk_sum_check;
  // int* nk_sum_arr_check;  

  // gpuErrchk(hipHostMalloc((void**) &nk_sum_check, n_pairs*sizeof(int)));
  // gpuErrchk(hipMemcpy( nk_sum_check, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToHost));

  // gpuErrchk(hipHostMalloc((void**) &nk_sum_arr_check, n_pairs*sizeof(int)));  
  // gpuErrchk(hipMemcpy( nk_sum_arr_check, nk_flags_sum_d, n_pairs*sizeof(int), hipMemcpyDeviceToHost));  
  // hipDeviceSynchronize();

  // // printf("CHECK: i = %d, nk_sum_check %d, nk_flag_sum_check %d \n", 0, nk_sum_check[0], nk_sum_arr_check[0]);  
  // // for(int i = 1; i < n_pairs; ++i){
  // //   nk_sum_check[i] += nk_sum_check[i-1];
  // //   printf("CHECK: i = %d, nk_sum_check %d, nk_flag_sum_check %d \n", i, nk_sum_check[i], nk_sum_arr_check[i]);      
  // //   //    printf("CHECK: i = %d, nk_sum_check %d \n", i, nk_sum_check[i]);
  // // }
  // gpuErrchk(hipHostFree(nk_sum_check));
  // gpuErrchk(hipHostFree(nk_sum_arr_check));  

  //---------------------------------------------//

  
  
  
}



// Function to find the next power of two
int nextPowerOfTwo(int x) {
    if (x < 1) return 1;
    int power = 1;
    while (power < x) {
        power <<= 1;
    }
    return power;
}



__global__
void kernel_set_pair_distribution_k_index(int i_beg, int i_end, int n_pairs,  int n_sites0, int* neighbors_list_d,
					  double* rjs, double* xyz,
					  int* k_index_d, int* j2_index_d,
					  double* rjs_index_d, double* xyz_index_d, int* nk_flags_d, int* nk_sum_flags_d){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int i,j2,nk;

  if( tid < n_pairs){
    //    printf("tid = %d  nk_flags = %d  nk_sum_flags %d \n", tid, nk_flags_d[tid], nk_sum_flags_d[tid]);
    
    if( nk_flags_d[tid] == 1 ){

      nk = nk_sum_flags_d[tid]-1;      
      if( tid == n_pairs-1 ){
	nk = nk_sum_flags_d[tid-1]; 
      }
      
      //      printf(" tid = %d, nk = %d\n", tid, nk);
      k_index_d[nk] = tid;

      j2 = ( (neighbors_list_d[tid] - 1) % n_sites0 );
      j2_index_d[nk] = j2;

      rjs_index_d[nk] = rjs[tid];      
      
      xyz_index_d[3*nk    ] = xyz[3*tid    ];
      xyz_index_d[3*nk + 1] = xyz[3*tid + 1];
      xyz_index_d[3*nk + 2] = xyz[3*tid + 2];            
    }
  }
}



extern "C" void  gpu_set_pair_distribution_k_index(int i_beg, int i_end, int n_pairs, int n_sites0,  int* neighbors_list,
						   double* rjs, double* xyz,
						   int* k_index_d, int* j2_index_d, double* rjs_index_d,
						   double* xyz_k_d, int* nk_flags_d, int* nk_sum_flags_d,
						   hipStream_t *stream ){

  // First do an inclusive scan to a temporary array which is of the
  // size of the flags so we get the nk index efficiently
  
    
  
  // Allocate temporary array and copy device data
  // int n = n_pairs; 
  // int* nk_sum_flags_d;

  // int paddedN = ((n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2)) * (BLOCK_SIZE * 2);
  // size_t size = paddedN * sizeof(int);


  // gpuErrchk(hipMalloc((void**)&nk_sum_flags_d, n*sizeof(int)));
  // //  printf("-- Copying nk_flags_d to nk_sum_flags_d \n");
  // gpuErrchk(hipMemcpy( nk_sum_flags_d, nk_flags_d, n*sizeof(int), hipMemcpyDeviceToDevice));
  //  hipMemset(nk_sum_flags_d + n , 0, (paddedN - n ) * sizeof(int));  // Set padding elements to 0

  // Perform recursive inclusive scan
  //  recursiveInclusiveScan(nk_sum_flags_d, n);

  
  
  // // Allocate the array 
  // gpuErrchk(hipMalloc((void**)&nk_sum_flags_d, paddedN*sizeof(int)));
  // //  printf("-- Copying nk_flags_d to nk_sum_flags_d \n");
  // gpuErrchk(hipMemcpy( nk_sum_flags_d, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToDevice));
  // hipMemset(nk_sum_flags_d + n , 0, (paddedN - n ) * sizeof(int));  // Set padding elements to 0

  //---------------------------------------------//
  // --- Check that the host can do the same --- //

  // int* nk_sum_check;
  // gpuErrchk(hipHostMalloc((void**) &nk_sum_check, n_pairs*sizeof(int)));
  // gpuErrchk(hipMemcpy( nk_sum_check, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToHost));
  // hipDeviceSynchronize();

  // printf("CHECK: i = %d, nk_sum_check %d \n", 1, nk_sum_check[0]);  
  // for(int i = 1; i < n_pairs; ++i){
  //   nk_sum_check[i] += nk_sum_check[i-1];
  //   printf("CHECK: i = %d, nk_sum_check %d \n", i, nk_sum_check[i]);
  // }
  // gpuErrchk(hipFree(nk_sum_check));

  //---------------------------------------------//

  
  printf("> Starting set pair distribution arrays \n");

  // int numBlocks = (n + BLOCK_SIZE * 2 - 1) / (BLOCK_SIZE * 2);

  // // Arrays for the block sums 
  // int* d_block_sums;

  // if (numBlocks > 1) {
  //   gpuErrchk(hipMalloc(&d_block_sums, numBlocks * sizeof(int)));
  // }
  
  // hipDeviceSynchronize();
  // // Perform the scan on each block
  // printf(">-- Starting inclusive scan \n");
  // //  inclusiveScanout(nk_flags_d, nk_sum_flags_d, n  );
  
  // inclusiveScanKernel<<<numBlocks, BLOCK_SIZE>>>(nk_sum_flags_d, d_block_sums, n);
  // gpuErrchk( hipPeekAtLastError() );
  // hipDeviceSynchronize();//  hipStreamSynchronize(stream[0]);  

  // if (numBlocks > 1) {
  //   // // Inclusive scan the block sums to accumulate the total sum 
  //   // inclusiveScanKernel<<<1, BLOCK_SIZE>>>(d_block_sums, d_block_sums, numBlocks);
  //   // hipDeviceSynchronize();

  //   // We would have to recursively call the kernel as the number of
  //   // blocks may exceed the block size, I trie dthis but couldnt get it fully working so for now I'm just going to
  //   // sum them on the host

    
  //   //--- Sum block sums on host for ease of implementation, see above comment
  //   int* h_block_sums;
  //   gpuErrchk(hipHostMalloc((void**) &h_block_sums, numBlocks*sizeof(int)));
  //   gpuErrchk(hipMemcpy( h_block_sums, d_block_sums, numBlocks*sizeof(int), hipMemcpyDeviceToHost));

  //   for( int i = 1; i < numBlocks; ++i ){
  //     h_block_sums[i] += h_block_sums[i-1];
  //   }
  //   gpuErrchk(hipMemcpy( d_block_sums, h_block_sums, numBlocks*sizeof(int), hipMemcpyHostToDevice));
    
  //   gpuErrchk(hipHostFree(h_block_sums));

  //   addBlockSumsKernel<<<numBlocks, BLOCK_SIZE>>>(nk_sum_flags_d, d_block_sums, n);
  //   gpuErrchk( hipPeekAtLastError() );
  //   hipDeviceSynchronize();    
  //   //    hipStreamSynchronize(stream[0]);      
  // }

  // //  gpuErrchk(hipMemcpy(h_data, nk_sum_flags_d, n * sizeof(int), hipMemcpyDeviceToHost));

  // //gpuErrchk(hipFree(d_data));
  // //  if (numBlocks > 1) {

  // gpuErrchk(hipFree(d_block_sums));


  // //    gpuErrchk(hipHostFree(h_block_sums));    
  //   //}

  // // // fflush(stdin);

  
  dim3 nblocks=dim3((n_pairs + tpb-1)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  // Can do a cumulative reduction of nk to get the nk indices for the flags 
  gpuErrchk( hipPeekAtLastError() );
  hipDeviceSynchronize();

  kernel_set_pair_distribution_k_index<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, n_pairs, n_sites0, neighbors_list,
									   rjs, xyz, k_index_d, j2_index_d, rjs_index_d,
									   xyz_k_d, nk_flags_d,
									   nk_sum_flags_d);
  printf(">> Set k index arrays  \n");
  gpuErrchk( hipPeekAtLastError() );
  //gpuErrchk(hipFreeAsync(nk_sum_flags_d, stream[0]) );  
  //  hipStreamSynchronize(stream[0]);  
  hipDeviceSynchronize();
}






__global__
void kernel_setup_matrix_forces(int i_beg, int i_end, int n_pairs,  int n_sites0, int* neighbors_list_d,
					   double* xyz_all_d,
					  int* k_index_d, int* j2_index_d,
					   double* xyz_index_d, int* nk_flags_d, int* nk_sum_flags_d){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int i,j2,nk;
  
  if( tid < n_pairs){
    //    printf("tid = %d  nk_flags = %d  nk_sum_flags %d \n", tid, nk_flags_d[tid], nk_sum_flags_d[tid]);

    if( nk_flags_d[tid] == 1 ){

      nk = nk_sum_flags_d[tid]-1;      
      if( tid == n_pairs-1 ){
	nk = nk_sum_flags_d[tid-1]; 
      }

      printf(" tid = %d, nk = %d\n", tid, nk);
      k_index_d[nk] = tid;

      j2 = ( (neighbors_list_d[tid] - 1) % n_sites0 )+1;
      j2_index_d[nk] = j2;

      //      rjs_index_d[nk] = rjs[tid];      
      //      (1:3, Nk_max), (i-1) + 3 * nk
      xyz_index_d[3*nk    ] = xyz_all_d[3*tid    ];
      xyz_index_d[3*nk + 1] = xyz_all_d[3*tid + 1];
      xyz_index_d[3*nk + 2] = xyz_all_d[3*tid + 2];            
    }
  }
}



extern "C" void  gpu_setup_matrix_forces(int i_beg, int i_end, int n_pairs, int n_sites0,  int* neighbors_list,
					 double* xyz_all_d, int* k_index_d, int* j2_index_d, 
					 double* xyz_k_d, int* nk_flags_d, int* nk_flags_sum_d,
						   hipStream_t *stream ){

  // First do an inclusive scan to a temporary array which is of the
  // size of the flags so we get the nk index efficiently
  
    
  
  // Allocate temporary array and copy device data
  //  int n = n_pairs; 
  //---------------------------------------------//
  // --- Check that the host can do the same --- //

  // int* nk_sum_check;
  // gpuErrchk(hipHostMalloc((void**) &nk_sum_check, n_pairs*sizeof(int)));
  // gpuErrchk(hipMemcpy( nk_sum_check, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToHost));
  // hipDeviceSynchronize();

  // printf("CHECK: i = %d, nk_sum_check %d \n", 1, nk_sum_check[0]);  
  // for(int i = 1; i < n_pairs; ++i){
  //   nk_sum_check[i] += nk_sum_check[i-1];
  //   printf("CHECK: i = %d, nk_sum_check %d \n", i, nk_sum_check[i]);
  // }
  // gpuErrchk(hipFree(nk_sum_check));

  //---------------------------------------------//

  
  dim3 nblocks=dim3((n_pairs + tpb-1)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  // Can do a cumulative reduction of nk to get the nk indices for the flags 
  hipDeviceSynchronize();
  printf(">> Set matrix force arrays  \n");
  kernel_setup_matrix_forces<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, n_pairs, n_sites0, neighbors_list,
									   xyz_all_d, k_index_d, j2_index_d, 
									   xyz_k_d, nk_flags_d,
									   nk_flags_sum_d);

  //  gpuErrchk(hipFreeAsync(nk_sum_flags_d, stream[0]) );  
  //  hipStreamSynchronize(stream[0]);  
  hipDeviceSynchronize();
  printf(">> Finished set matrix force arrays  \n");  
}



// ----------------------------------------------------- //
// --- Pair distribution calculation and derivatives --- //
// ----------------------------------------------------- //


  /*
--- Thinking behind this kernel ---



- We normally iterate over each r and the for each r compute the gaussian 
- Each thread can get a different r
- Each thread 
- n_samples will generally be less than the number of threads, but also what we could do is just work on on


- I could actually try and set it up that each block works on a single r value and we then warp reduce to get the sum 
- The question then becomes how many blocks do I go over
- So we have 
- warp_size = 32 lets say 
- So 32 k values are sent to each thread in the warp 
- nblocks = ( nk + 32 - 1) / 32 
- This is for one n_samples l index... 

- So then maybe I actually just send all of them 
- nblocks = n_samples * ( nk + 32 - 1) / 32 
- if blockIdx.x / ( ( nk + 32 - 1) / 32  ) = 1 then we are on the next l value 
- So make an array which stores the sum which is the size of the number of blocks 
- Then do a parallel reduce on sections of the sum to get the final result
- These sections will be of size ((nk + 32 -1)/ 32  +  1) 


  */  

  // do l = 1,n_samples
  //    kde(l) = kde(l) +  exp( -( (x(l) - r) / kde_sigma )**2 / 2.d0 )
  // end do
  // pair_distribution(1:n_samples) = pair_distribution(1:n_samples) + &
  //      & kde(1:n_samples)

  // if (do_derivatives)then
  //    ! reuse kde to get the derivatives
  //    do l = 1,n_samples
  //       kde(l) = kde(l)  *  ( (x(l) - r) / kde_sigma**2 ) / r / dV(l)
  //    end do

  //    pair_distribution_der(1:n_samples, n_dim_idx,  k) = &
  //         & pair_distribution_der(1:n_samples, n_dim_idx, &
  //         & k) + kde(1:n_samples)

  // end if


  // __shared__ double      x_pdf[ n_samples ];
  // __shared__ double     dV_pdf[ n_samples ];    
  // __shared__ double shared_pdf[ n_samples ];



  // if( tid < n_samples ){
  //   x_pdf[ tid ] = x[ tid ];
  //   dV_pdf[ tid ] = dV[ tid ];    
  //   shared_pdf[ tid ] = 0.0;

    

  // }

  // int k = k_index_d[tid];
  // double r = rjs_d[tid];
      
  // __syncthreads();

    
  // if(tid < n_samples){
  //   for( int i = 0; i < n_samples; ++i ){
      

  //   }
  // }

  

__global__
void kernel_get_pair_distribution_kde( double* pdf_out, double* pdf_der_out,
				       int n_k, int n_samples,
				       double kde_sigma, double* x_d, double* dV_d, double* rjs_d, double pdf_factor){

  int utid=threadIdx.x+blockIdx.x*blockDim.x;
  // int k_index = blockIdx.x * blockDim.x + threadIdx.x;  // Index in the j dimension (columns)
  // int l_index = blockIdx.y * blockDim.y + threadIdx.y;  // Index in the l dimension (rows)
  
  int i;
  double x, dV, r, x_gauss, pdf_temp, pdf_der_temp;

  int k_index = utid / n_samples;    
  int l_index = utid % n_samples;

  if ( k_index < n_k && l_index < n_samples ){
    
    x  =  x_d[l_index];
    dV = dV_d[l_index];

    r = rjs_d[k_index];        

    x_gauss = ( (x - r) / kde_sigma );
    
    pdf_temp =  pdf_factor * exp( - 0.5 * x_gauss * x_gauss  ) / dV ;

    
    if( utid % 10000 == 0 ){
      double e = exp( - 0.5 * (x_gauss * x_gauss)  );
      printf(" pdf_temp %lf, x %lf, dV %lf, x_gauss %lf, exp %lf,  k %d/ %d, l %d / %d\n", pdf_temp, x, dV, x_gauss, e,  k_index, n_k, l_index, n_samples);
    }

    
    pdf_der_temp = (pdf_temp * x_gauss) / (kde_sigma * r);// / dV;

    //  l_index + n_samples * k_index
    pdf_out[  k_index + l_index * n_k    ] = pdf_temp;    

    pdf_der_out[ l_index + n_samples * k_index ] = pdf_der_temp;
  }
    
}

__global__
void kernel_reduce_pair_distribution( double* pdf_in, double* pdf_out, int n_k, int n_samples){

  int tid = threadIdx.x;
  int i, stride ;
  double pdf_temp;
  // Now we want to reduce over the n_samples portions of continuous memory 

  int n_strides = ( n_k + BLOCK_SIZE - 1) / BLOCK_SIZE;

  __shared__ double sharedData[BLOCK_SIZE];


  pdf_temp = 0.0;
  
  for( i = 0; i < n_strides; ++i ){
    // int k_index = (utid + i) / n_samples;    
    // int l_index = (utid + i) % n_samples;

    int k_index = tid + i*BLOCK_SIZE  ;//(utid + i) / n_samples;    
    int l_index = blockIdx.x;

					  
    if ( k_index < n_k && l_index < n_samples ){
      // l_index + n_samples * k_index
      pdf_temp += pdf_in[ k_index + l_index * n_k ];
    }
  }

  sharedData[tid] = pdf_temp;

  // Perform reduction in shared memory
  for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
    __syncthreads();
    if (tid < s) {
      sharedData[tid] += sharedData[tid + s];
    }
  }

  // Write the result of the block to global memory
  if (tid == 0) {
    pdf_out[blockIdx.x] =  sharedData[0] ;//      d_out[blockIdx.x] = val;
  }
}
  


__global__ void reduce_k_index(double* pdf_in, double* pdf_out, int n_samples, int k_size) {
    int l_index = blockIdx.x * blockDim.x + threadIdx.x;

    // Ensure we're within bounds
    if (l_index >= n_samples) return;

    float sum = 0.0f;

    // Reduce over k_index
    for (int k_index = 0; k_index < k_size; ++k_index) {
        sum += pdf_in[l_index + n_samples * k_index];
    }

    pdf_out[l_index] = sum ;
}




extern "C" void  gpu_get_pair_distribution_and_ders( double* pair_distribution_d, double* pair_distribution_der_d,
						     int n_k, int n_samples, double kde_sigma, double* x_d, double* dV_d, 
						     double* rjs_d, double pdf_factor,  hipStream_t *stream ){

  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  int threads = BLOCK_SIZE;
  dim3 nblocks=dim3((n_k * n_samples + threads-1)/threads,1,1);
  dim3 nthreads=dim3(threads,1,1);

  double* pdf_to_reduce;
  hipMalloc(&pdf_to_reduce, n_k * n_samples * sizeof(double) );


  // dim3 blockDim(16, 16); // 16x16 threads per block
  // dim3 gridDim((n_k + blockDim.x - 1) / blockDim.x, (n_samples + blockDim.y - 1) / blockDim.y);

  
  gpuErrchk( hipPeekAtLastError() );
  printf("> pdf evaluation kernel starting\n");
  kernel_get_pair_distribution_kde<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce,  pair_distribution_der_d,
									n_k, n_samples, kde_sigma, x_d, dV_d, rjs_d, pdf_factor );
  hipDeviceSynchronize();
  gpuErrchk( hipPeekAtLastError() );
  printf("> pdf evaluation kernel finished\n");

  // Then we need to reduce over the size of the blocks for the pair distribution 

  //  // More complex reduction
  //  //----------------------------------------
  nblocks=dim3(n_samples,1,1);
  nthreads=dim3(threads,1,1);
  printf("> pdf reduction kernel starting\n");
  kernel_reduce_pair_distribution<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_d, n_k, n_samples );
  printf("> pdf reduction kernel finished\n");
  //  //----------------------------------------

  // Simple reduction
  // nblocks=dim3((n_samples + threads-1)/threads,1,1);
  // nthreads=dim3(threads,1,1);
  // printf("> pdf reduction kernel starting\n");
  // reduce_k_index<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_d, n_samples, n_k);  
  // printf("> pdf reduction kernel finished\n");
  
  
  hipFree(pdf_to_reduce);
  hipDeviceSynchronize();  
  gpuErrchk( hipPeekAtLastError() );
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


__global__
void kernel_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int j, k, j2;
  k = tid / n_samples;
  j = tid % n_samples;  

  // Gka(j,k) = xyz_k(i, k) * Gk(j, k)
  // Gka = ( n_samples * n_k )
  // (2,3) = 1,1,,2,1,,3,1,,4,1,, n_samples * (k-1) + (j-1),    (4,1), (4-1) + (1-1)*n_k = 3 
  // (j,k) == ( j - 1 ) + n_samples * (k - 1)
  //xyz = (1,1), (2,1), (3,1)
  if( k < n_k && j < n_samples ){
    Gka_d[ j + k * n_samples ] = Gk_d[ j + k * n_samples ] * xyz_k_d[ (i-1) + k * 3 ]; 
  }
  
  
}


extern "C" void  gpu_get_Gka(int i, int n_k, int n_samples,
			     double* Gka_d, double* Gk_d, double* xyz_k_d, hipStream_t* stream ){

  int threads = BLOCK_SIZE;
  dim3 nblocks=dim3((n_k * n_samples + threads-1)/threads,1,1);
  dim3 nthreads=dim3(threads,1,1);

  kernel_get_Gka<<<nblocks, nthreads, 0, stream[0]>>>(i, n_k, n_samples, Gka_d, Gk_d, xyz_k_d );

  
}



__global__
void kernel_hadamard_vec_mat_product(int n_samples_sf, int n_k,
				     double* all_scattering_factors_d, double* dermat_d){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int j, l, j2;
  j = tid / n_samples_sf;
  l = tid % n_samples_sf;  

  if( j < n_k && l < n_samples_sf ){
    dermat_d[ l + j*n_samples_sf  ] *= all_scattering_factors_d[ l ]; 
  }
  
  
}

__global__ void dermat_kernel(double *dermat, double *all_scattering_factors, int n_samples_sf, int n_k) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;  // Index in the j dimension (columns)
    int l = blockIdx.y * blockDim.y + threadIdx.y;  // Index in the l dimension (rows)

    if (j < n_k && l < n_samples_sf) {

        dermat[l + j * n_samples_sf] *= all_scattering_factors[l];
    }
}

extern "C" void  gpu_hadamard_vec_mat_product(int n_samples_sf, int n_k,
					      double* all_scattering_factors_d, double* dermat_d, hipStream_t* stream ){

  // int threads = BLOCK_SIZE;
  // dim3 nblocks=dim3((n_k * n_samples_sf + threads-1)/threads,1,1);
  // dim3 nthreads=dim3(threads,1,1);

  // kernel_hadamard_vec_mat_product<<<nblocks, nthreads, 0, stream[0]>>>(n_samples_sf, n_k, all_scattering_factors_d, dermat_d );


  
  // Define grid and block dimensions
  dim3 blockDim(16, 16); // 16x16 threads per block
  dim3 gridDim((n_k + blockDim.x - 1) / blockDim.x, (n_samples_sf + blockDim.y - 1) / blockDim.y);

  
    // Launch the kernel
  dermat_kernel<<<gridDim, blockDim, 0, stream[0]>>>(dermat_d, all_scattering_factors_d, n_samples_sf, n_k);

  
  
}



extern "C" void gpu_get_fi_dgemv(const int i, const int n_samples_sf, const int n_k, double* dermat_d,
				 double* prefactor_d, double* fi_d, hipblasHandle_t handle, hipStream_t* stream ){

  	const double alf = 1;
	const double bet = 0;
	const double *alpha = &alf;
	const double *beta = &bet;
	double *ptr = fi_d + (i - 1) * n_k;

	// Do the actual multiplication
       // ! Now we take the dot products by matrix vector
       // call dgemv("T", n_samples_sf,  n_k, 1.d0 , dermat, n_samples_sf,&
       //      &  prefactor, 1, 0.d0, fi(:,i), 1)
	
	hipblasDgemv(handle, HIPBLAS_OP_T, n_samples_sf, n_k, alpha,
		     dermat_d, n_samples_sf, prefactor_d, 1, beta,
		     ptr, 1);


}


__global__
void kernel_exp_force_virial_collection(int n_k, double3* forces0, double energy_scale, double* fi,
						int* j2_list, double* virial, double3* xyz){
  // tid == some n_k value 
  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int j2;
  
  double3 this_force;

  if( tid < n_k ){
    j2 = j2_list[ tid ]-1;    
    //    this_force = forces0[ j2 ];
    double fi_tmp[ 3 ];

    fi_tmp[0] = fi[ tid ];
    fi_tmp[1] = fi[ tid +     n_k ];
    fi_tmp[2] = fi[ tid + 2 * n_k ];        
    
    this_force.x = energy_scale * fi_tmp[0];
    this_force.y = energy_scale * fi_tmp[1];
    this_force.z = energy_scale * fi_tmp[2];

    atomicAdd(&forces0[j2].x, this_force.x);
    atomicAdd(&forces0[j2].y, this_force.y);
    atomicAdd(&forces0[j2].z, this_force.z);

    //    forces0(1:3, j2_list(j)) = forces0(1:3, j2_list(j)) + this_force(1:3)

    double tmp_this_force[3];
    tmp_this_force[0] = this_force.x;
    tmp_this_force[1] = this_force.y;
    tmp_this_force[2] = this_force.z;
    
    
    double3 tmp_xyz;
    tmp_xyz=xyz[tid];
    double this_xyz[3];
    this_xyz[0]=tmp_xyz.x;
    this_xyz[1]=tmp_xyz.y;
    this_xyz[2]=tmp_xyz.z;

    for(int k1=0;k1<3;k1++){
      for(int k2=0;k2<3;k2++){
        double loc_viri=0.5*(tmp_this_force[k1]*this_xyz[k2] + tmp_this_force[k2]*this_xyz[k1]); 
        atomicAdd(&virial[k2+3*k1], loc_viri);
      }
    }

    
  }
  
}

extern "C" void gpu_exp_force_virial_collection(int n_k, double3* forces0, double energy_scale, double* fi,
						int* j2_list, double* virial, double3* xyz, hipStream_t *stream ) {
  
  // We can have a kernel go over the values of nk, which furnish us
  // with the j index for fi, and we can pass it ot the forces
  // j2_list.

  dim3 nblocks=dim3((n_k + tpb-1)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_exp_force_virial_collection<<<nblocks, nthreads, 0, stream[0]>>>( n_k, forces0, energy_scale, fi,
									   j2_list, virial, xyz);
  
  
}


extern "C" void gpu_print_pointer_int(int* p) {  printf(" address:  %p \n", p);  }
extern "C" void gpu_print_pointer_double(double* p) {  printf(" address:  %p \n", p);  }


// do j = 1, n_k
//        this_force(1:3) = energy_scale * fi(j,1:3)
//        forces0(1:3, j2_list(j)) = forces0(1:3, j2_list(j)) + this_force(1:3)

//        do k1 = 1, 3
//           do k2 =1, 3
//              virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
//           end do
//        end do
//     end do




// extern "C" void gpu_stream_synchronize(hipStream_t *stream){
//   hipStreamSynchronize(stream[0]);
// }