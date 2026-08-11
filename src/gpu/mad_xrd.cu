// MAD: the X-ray/neutron structure factor.
//
// Gk and its derivative with respect to atomic positions, the Hadamard product
// with the scattering factors, the dgemv that turns that into per-pair forces,
// and the collection of those into forces and the virial.
#include "gpu_common.h"
#include "mad_gpu.h"
#include "gpu_scatter.h"

#define tpb 512

__global__ void kernel_set_Gka(int nk, int n_samples, int* k_index_d, double* Gk_d, double* pair_distribution_partial_der_d,
                               double c_factor) {
  long long tid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  long long l, k;
  k = tid / n_samples;
  l = tid % n_samples;

  // Gk(1:n_samples, n_k) =  -2.d0 *  c_factor * pair_distribution_der(1:n_samples,  k )

  if (k < nk && l < n_samples) {
    // if( isnan( pair_distribution_partial_der_d[l + k*n_samples  ] ) ){
    //   printf("pair_distribution_partial_der_d[ l + k*n_samples  ] is NaN! l = %d, k = %d, l + k*n_samples = %d \n", l, k, l + k*n_samples);
    // }


    Gk_d[l + k * n_samples] = -2.0 * c_factor * pair_distribution_partial_der_d[l + k * n_samples];
  }
}


extern "C" void gpu_set_Gk(int nk, int n_samples, int* k_index_d, double* Gk_d, double* pair_distribution_partial_der_d,
                           double c_factor, hipStream_t* stream) {
  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) nk * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  kernel_set_Gka<<<nblocks, nthreads, 0, stream[0]>>>(nk, n_samples, k_index_d, Gk_d, pair_distribution_partial_der_d, c_factor);
}


__global__ void kernel_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d) {
  long long tid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  long long j, k, j2;
  k = tid / n_samples;
  j = tid % n_samples;

  if (k < n_k && j < n_samples) {
    Gka_d[j + k * n_samples] = Gk_d[j + k * n_samples] * xyz_k_d[(i - 1) + k * 3];


    // if( isnan( Gk_d[ j + k*n_samples  ] ) ){
    //   printf("Gk_d  setup [ l + j*n_samples  ] is NaN!    j = %d, k = %d, j + k*n_samples = %d, n_k %d\n", j, k, j + k*n_samples, n_k);
    // }

    // if( isnan( Gka_d[ j + k*n_samples  ] ) ){
    //   printf("xyz_k_d is probably nan!                    j = %d, k = %d, j + k*n_samples = %d, n_k %d\n", j, k, j + k*n_samples, n_k);
    // }
  }
}


extern "C" void gpu_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d, hipStream_t* stream) {
  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  kernel_get_Gka<<<nblocks, nthreads, 0, stream[0]>>>(i, n_k, n_samples, Gka_d, Gk_d, xyz_k_d);
}


__global__ void kernel_get_Gka_inplace(int i, int n_k, int n_samples, double* Gk_d, double* xyz_k_d) {
  long long tid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  long long j, k, j2;
  k = tid / n_samples;
  j = tid % n_samples;

  if (k < n_k && j < n_samples) {
    Gk_d[j + k * n_samples] *= xyz_k_d[(i - 1) + k * 3];
  }
}


extern "C" void gpu_get_Gka_inplace(int i, int n_k, int n_samples, double* Gk_d, double* xyz_k_d, hipStream_t* stream) {
  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  kernel_get_Gka_inplace<<<nblocks, nthreads, 0, stream[0]>>>(i, n_k, n_samples, Gk_d, xyz_k_d);
}


__global__ void kernel_hadamard_vec_mat_product(int n_samples_sf, int n_k, double* all_scattering_factors_d, double* dermat_d) {
  long long tid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  long long j, l;
  j = tid / n_samples_sf;
  l = tid % n_samples_sf;

  if (j < n_k && l < n_samples_sf) {
    // if( isnan( dermat_d[ l + j*n_samples_sf  ] ) ){
    //   printf("dermat_d     [ l + j*n_samples_sf  ] is NaN! l = %d, j = %d, l + j*n_samples_sf = d \n", l, j, l + j*n_samples_sf);
    // }

    dermat_d[l + j * n_samples_sf] *= all_scattering_factors_d[l];
  }
}

__global__ void dermat_kernel(double* dermat, double* all_scattering_factors, int n_samples_sf, int n_k) {
  long long j = blockIdx.x * blockDim.x + threadIdx.x; // Index in the j dimension (columns)
  long long l = blockIdx.y * blockDim.y + threadIdx.y; // Index in the l dimension (rows)

  if (j < n_k && l < n_samples_sf) {
    dermat[l + j * n_samples_sf] *= all_scattering_factors[l];
  }
}

extern "C" void gpu_hadamard_vec_mat_product(int n_samples_sf, int n_k, double* all_scattering_factors_d, double* dermat_d,
                                             hipStream_t* stream) {
  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples_sf + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  kernel_hadamard_vec_mat_product<<<nblocks, nthreads, 0, stream[0]>>>(n_samples_sf, n_k, all_scattering_factors_d, dermat_d);


  // // Define grid and block dimensions
  // dim3 blockDim(16, 16); // 16x16 threads per block
  // dim3 gridDim((n_k + blockDim.x - 1) / blockDim.x, (n_samples_sf + blockDim.y - 1) / blockDim.y);


  //   // Launch the kernel
  // dermat_kernel<<<gridDim, blockDim, 0, stream[0]>>>(dermat_d, all_scattering_factors_d, n_samples_sf, n_k);
}


extern "C" void gpu_get_fi_dgemv(const int i, const int n_samples_sf, const int n_k, double* dermat_d, double* prefactor_d,
                                 double* fi_d, hipblasHandle_t handle, hipStream_t* stream) {
  const double alf = 1;
  const double bet = 0;
  const double* alpha = &alf;
  const double* beta = &bet;
  double* ptr = fi_d + (i - 1) * n_k;

  // Do the actual multiplication
  // ! Now we take the dot products by matrix vector
  // call dgemv("T", n_samples_sf,  n_k, 1.d0 , dermat, n_samples_sf,&
  //      &  prefactor, 1, 0.d0, fi(:,i), 1)

  hipblasDgemv(handle, HIPBLAS_OP_T, n_samples_sf, n_k, alpha, dermat_d, n_samples_sf, prefactor_d, 1, beta, ptr, 1);
}


// The pair's own force, staged at its own index. It used to be added straight
// into forces0[j2] and the virial into nine addresses, both with atomicAdd,
// which made the result depend on the order the threads arrived; the sum is
// gpu_pair_scatter_reduce's now. See gpu_scatter.h.
__global__ void kernel_exp_force_virial_collection(int n_k, double* pair_force, double energy_scale, double* fi,
                                                   int* j2_list) {
  // tid == some n_k value
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  double3 this_force;

  if (tid < n_k) {
    double fi_tmp[3];

    fi_tmp[0] = fi[tid];
    fi_tmp[1] = fi[tid + n_k];
    fi_tmp[2] = fi[tid + 2 * n_k];

    // if(isnan( fi_tmp[0] )){ printf("> tid = %d, n_k = %d, fi_tmp[0] = %lf\n", tid, n_k,  fi_tmp[0]); }
    // if(isnan( fi_tmp[1] )){ printf("> tid = %d, n_k = %d, fi_tmp[1] = %lf\n", tid, n_k,  fi_tmp[1]); }
    // if(isnan( fi_tmp[2] )){ printf("> tid = %d, n_k = %d, fi_tmp[2] = %lf\n", tid, n_k,  fi_tmp[2]); }
    // if(isnan( energy_scale )){ printf("> tid = %d, n_k = %d,  energy_scale = %lf\n",  tid, n_k,  energy_scale); }

    this_force.x = energy_scale * fi_tmp[0];
    this_force.y = energy_scale * fi_tmp[1];
    this_force.z = energy_scale * fi_tmp[2];

    //    forces0(1:3, j2_list(j)) = forces0(1:3, j2_list(j)) + this_force(1:3)
    pair_force[3 * tid] = this_force.x;
    pair_force[3 * tid + 1] = this_force.y;
    pair_force[3 * tid + 2] = this_force.z;
  }
}

extern "C" void gpu_exp_force_virial_collection(int n_k, int n_sites, double3* forces0, double energy_scale, double* fi,
                                                int* j2_list, double* virial, double3* xyz, hipStream_t* stream) {
  // We can have a kernel go over the values of nk, which furnish us
  // with the j index for fi, and we can pass it ot the forces
  // j2_list.

  if (n_k <= 0)
    return;

  dim3 nblocks = dim3((n_k + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  double* pair_force_d;
  gpuErrchk(hipMallocAsync(&pair_force_d, (size_t) 3 * n_k * sizeof(double), stream[0]));

  kernel_exp_force_virial_collection<<<nblocks, nthreads, 0, stream[0]>>>(n_k, pair_force_d, energy_scale, fi, j2_list);

  // 0.25 = 0.5 to symmetrise, times 0.5 because k_list holds ordered pairs and
  // the reversed one reproduces the term rather than adding a new one.
  gpu_pair_scatter_reduce(n_k, n_sites, j2_list, pair_force_d, (const double*) xyz, (double*) forces0, virial, 0.25, stream);

  gpuErrchk(hipFreeAsync(pair_force_d, stream[0]));
}
