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
void recursiveReduce(int* d_in, int* d_out, int n, hipStream_t *stream) {
    int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // If only one block remains, no further recursion is needed
    if (numBlocks == 1) {
      blockReduceKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_in, d_out, n);
      //        hipDeviceSynchronize();
    } else {
        // Allocate memory for intermediate results
        int* d_intermediate;
        hipMalloc(&d_intermediate, numBlocks * sizeof(int));

        // Perform the reduction on the blocks
        blockReduceKernel<<<numBlocks, BLOCK_SIZE, 0, stream[0]>>>(d_in, d_intermediate, n);
        //hipDeviceSynchronize();

        // Recursively reduce the intermediate results
        recursiveReduce(d_intermediate, d_out, numBlocks, stream);

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
    if (globalIndexA < n) d_data[globalIndexA] = sharedData[bankOffsetA];
    if (globalIndexB < n) d_data[globalIndexB] = sharedData[bankOffsetB];
}

// Kernel to add block sums to each element
__global__ void addBlockSumsKernel(int* d_data, int* d_blockSums, int n) {
    int globalIndex = threadIdx.x + blockIdx.x * blockDim.x * 2;
    if (blockIdx.x > 0) {
        int blockSum = d_blockSums[blockIdx.x - 1];
        if (globalIndex < n) d_data[globalIndex] += blockSum;
        if (globalIndex + blockDim.x < n) d_data[globalIndex + blockDim.x] += blockSum;
    }
}



// Function to perform an inclusive scan on an array
void inclusiveScan(int* d_data_out, int n, hipStream_t *stream) {
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
    inclusiveScanKernel<<<numBlocks, BLOCK_SIZE,0, stream[0]>>>(d_data, d_blockSums, n);
    hipDeviceSynchronize();

    // If there are multiple blocks, perform a scan on the block sums
    if (numBlocks > 1) {
      inclusiveScan(d_blockSums, numBlocks, stream);

      // Add block sums to each element in the array
      addBlockSumsKernel<<<numBlocks, BLOCK_SIZE,0,stream[0]>>>(d_data, d_blockSums, n);
    }

    // Copy result back to host
    hipMemcpy(d_data_out, d_data+1, n * sizeof(int), hipMemcpyDeviceToDevice);    

    // Free device memory
    hipFree(d_data);
    hipFree(d_blockSums);
}



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

  // Perform recursive reduction
  recursiveReduce(nk_flags_d, nk_out_d, n_pairs, stream);

  gpuErrchk(hipMemcpy( nk_flags_sum_d, nk_flags_d, n_pairs*sizeof(int), hipMemcpyDeviceToDevice));

  // Perform inclusive scan to get the nk indexes 
  inclusiveScan(nk_flags_sum_d, n_pairs, stream);
  
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
      j2_index_d[nk] = j2+1;

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


  
  dim3 nblocks=dim3((n_pairs + tpb-1)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_set_pair_distribution_k_index<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, n_pairs, n_sites0, neighbors_list,
									   rjs, xyz, k_index_d, j2_index_d, rjs_index_d,
									   xyz_k_d, nk_flags_d,
									   nk_sum_flags_d);

}



// ----------------------------------------------------- //
// --- Pair distribution calculation and derivatives --- //
// ----------------------------------------------------- //

__global__
void kernel_get_pair_distribution_kde( double* pdf_out, double* pdf_der_out,
				       int n_k, int n_samples,
				       double kde_sigma, double* x_d, double* dV_d, double* rjs_d,
				       double pdf_factor, double der_factor){

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
    
    pdf_der_temp = (der_factor * pdf_temp * x_gauss) / (kde_sigma * r);// / dV;

    // Reversed the indexing of this temporary array for memory coalescence 
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

    int k_index = tid + i*BLOCK_SIZE;
    int l_index = blockIdx.x;

					  
    if ( k_index < n_k && l_index < n_samples ){
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
    pdf_out[blockIdx.x] =  sharedData[0];
  }
}
  



extern "C" void  gpu_get_pair_distribution_and_ders( double* pair_distribution_d, double* pair_distribution_der_d,
						     int n_k, int n_samples, double kde_sigma, double* x_d, double* dV_d, 
						     double* rjs_d, double pdf_factor, double der_factor,  hipStream_t *stream ){

  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  int threads = BLOCK_SIZE;
  dim3 nblocks=dim3((n_k * n_samples + threads-1)/threads,1,1);
  dim3 nthreads=dim3(threads,1,1);

  double* pdf_to_reduce;
  hipMalloc(&pdf_to_reduce, n_k * n_samples * sizeof(double) );

  
  printf("> pdf evaluation kernel starting\n");
  kernel_get_pair_distribution_kde<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce,  pair_distribution_der_d,
									n_k, n_samples, kde_sigma, x_d, dV_d,
									rjs_d, pdf_factor, der_factor );
  //  hipDeviceSynchronize();
  // gpuErrchk( hipPeekAtLastError() );
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
  
  
  hipFree(pdf_to_reduce);
  //  hipDeviceSynchronize();  
  //gpuErrchk( hipPeekAtLastError() );
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
void kernel_set_Gka( int nk, int n_samples, int* k_index_d, double* Gk_d,
		     double*  pair_distribution_partial_der_d, double c_factor){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int l, k;
  k = tid / n_samples;
  l = tid % n_samples;  

  // Gk(1:n_samples, n_k) =  -2.d0 *  c_factor * pair_distribution_der(1:n_samples,  k )

  if( k < nk && l < n_samples ){
    Gk_d[ l + k * n_samples ] = - 2.0 * c_factor * pair_distribution_partial_der_d[ l + k * n_samples ]; 
  }
}


extern "C" void gpu_set_Gk(int nk, int n_samples, int* k_index_d, double* Gk_d,
			   double*  pair_distribution_partial_der_d, double c_factor,  hipStream_t *stream ){

  
  int threads = BLOCK_SIZE;
  dim3 nblocks=dim3((nk * n_samples + threads-1)/threads,1,1);
  dim3 nthreads=dim3(threads,1,1);

  kernel_set_Gka<<<nblocks, nthreads, 0, stream[0]>>>( nk, n_samples, k_index_d,
						       Gk_d, pair_distribution_partial_der_d, c_factor );

}

  

__global__
void kernel_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d){

  int tid=threadIdx.x+blockIdx.x*blockDim.x;  
  int j, k, j2;
  k = tid / n_samples;
  j = tid % n_samples;  

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
  int j, l;
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

  int threads = BLOCK_SIZE;
  dim3 nblocks=dim3((n_k * n_samples_sf + threads-1)/threads,1,1);
  dim3 nthreads=dim3(threads,1,1);

  kernel_hadamard_vec_mat_product<<<nblocks, nthreads, 0, stream[0]>>>(n_samples_sf, n_k, all_scattering_factors_d, dermat_d );


  
  // // Define grid and block dimensions
  // dim3 blockDim(16, 16); // 16x16 threads per block
  // dim3 gridDim((n_k + blockDim.x - 1) / blockDim.x, (n_samples_sf + blockDim.y - 1) / blockDim.y);

  
  //   // Launch the kernel
  // dermat_kernel<<<gridDim, blockDim, 0, stream[0]>>>(dermat_d, all_scattering_factors_d, n_samples_sf, n_k);
  
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

