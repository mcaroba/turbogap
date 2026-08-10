// MAD: the pair distribution function.
//
// Counting the pairs inside the pdf cutoff, compacting them into the k-indexed
// buffers, the kernel density estimate of g(r) and its derivative, and the
// reduction over pairs.
#include "gpu_common.h"
#include "gpu_scan.h"
#include "mad_gpu.h"

#define tpb 512


__global__ void kernel_get_pair_distribution_nk(int i_beg, int i_end, int n_sites0, int* neighbors_list_d, int* n_neigh_d,
                                                int* neighbor_species_d, int* species_d, double* rjs_d, double* xyz, double r_min,
                                                double r_max, double r_cut, double buffer, int* nk_flags_d, int sp1, int sp2) {
  int i_site = i_beg - 1 + threadIdx.x + blockIdx.x * blockDim.x;
  int k_val = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j, k, s, k1, k2;
  double r;
  int tid = threadIdx.x;
  int lane = tid % WARP_SIZE;
  int warpId = tid / WARP_SIZE;


  int nk_loc = 0;
  if (i_site < i_end) {
    if (species_d[i_site] == sp1 || species_d[i_site] == sp2) {
      k = 0;
      for (i = i_beg - 1; i < i_site; i++)
        k += n_neigh_d[i];
      for (j = 1; j < n_neigh_d[i_site]; j++) {
        k += 1;
        if (!((species_d[i_site] == sp1 && neighbor_species_d[k] == sp2) ||
              (species_d[i_site] == sp2 && neighbor_species_d[k] == sp1)))
          continue;
        r = rjs_d[k];
        if (r < 0.001 || r > r_cut || r < r_min || r > r_max + buffer)
          continue;
        nk_loc += 1;
        nk_flags_d[k] = 1;
      }
    }
  }
}


extern "C" void gpu_get_pair_distribution_nk(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, int* n_neigh,
                                             int* neighbor_species, int* species, double* rjs, double* xyz, double r_min,
                                             double r_max, double r_cut, double buffer, int* nk_out_d, int* nk_flags_d,
                                             int* nk_flags_sum_d, int species_1, int species_2, hipStream_t* stream) {
  // This function is to set the k_index array for the partial pair distributions


  dim3 nblocks = dim3((i_end - i_beg + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  // gpuErrchk( hipPeekAtLastError() );

  // "invalid configuration argument" names no dimension and no kernel, which
  // makes a zero-extent launch one of the more expensive errors here to track
  // down -- and splitting a cell into batches produces empty extents readily.
  // Say which launch and what it asked for.
  TG_CHECK_LAUNCH("kernel_get_pair_distribution_nk", nblocks, nthreads, "n_pairs", n_pairs);
  kernel_get_pair_distribution_nk<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, n_sites0, neighbors_list, n_neigh,
                                                                       neighbor_species, species, rjs, xyz, r_min, r_max, r_cut,
                                                                       buffer, nk_flags_d, species_1, species_2);

  hipStreamSynchronize(stream[0]);
  gpuErrchk(hipPeekAtLastError());
  fflush(stdout);

  // Perform recursive reduction
  recursiveReduce(nk_flags_d, nk_out_d, n_pairs, stream);
  hipStreamSynchronize(stream[0]);
  gpuErrchk(hipPeekAtLastError());
  fflush(stdout);

  gpuErrchk(hipMemcpyAsync(nk_flags_sum_d, nk_flags_d, n_pairs * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));
  hipStreamSynchronize(stream[0]);
  gpuErrchk(hipPeekAtLastError());
  fflush(stdout);

  // Perform inclusive scan to get the nk indexes
  inclusiveScan(nk_flags_sum_d, n_pairs, stream);

  hipStreamSynchronize(stream[0]);
  gpuErrchk(hipPeekAtLastError());
  fflush(stdout);
  // Now multiply the flags_sum_d and flags_d to get final nk_array

  gpu_multiply_flags(n_pairs, nk_flags_d, nk_flags_sum_d, stream);
  hipStreamSynchronize(stream[0]);
  gpuErrchk(hipPeekAtLastError());
  fflush(stdout);
}


__global__ void kernel_set_pair_distribution_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list_d,
                                                     double* rjs, double* xyz, int* k_index_d, int* j2_index_d, double* rjs_index_d,
                                                     double* xyz_index_d, // int* nk_flags_d,
                                                     int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j2, nk, nk_temp;

  if (tid < n_pairs) {
    if (nk_sum_flags_d[tid] > 0) {
      nk = nk_sum_flags_d[tid] - 1;

      if (tid == n_pairs - 1) {
        // Search for the last non-zero value
        i = 0;
        while (nk == 0) {
          i += 1;
          nk = nk_sum_flags_d[tid - i];
        }
      }

      k_index_d[nk] = tid;

      j2 = ((neighbors_list_d[tid] - 1) % n_sites0);
      j2_index_d[nk] = j2 + 1;

      rjs_index_d[nk] = rjs[tid];

      xyz_index_d[3 * nk] = xyz[3 * tid];
      xyz_index_d[3 * nk + 1] = xyz[3 * tid + 1];
      xyz_index_d[3 * nk + 2] = xyz[3 * tid + 2];
    }
  }
}


extern "C" void gpu_set_pair_distribution_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, double* rjs,
                                                  double* xyz, int* k_index_d, int* j2_index_d, double* rjs_index_d,
                                                  double* xyz_k_d, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_set_pair_distribution_k_index<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, n_pairs, n_sites0, neighbors_list, rjs,
                                                                            xyz, k_index_d, j2_index_d, rjs_index_d,
                                                                            xyz_k_d, //nk_flags_d,
                                                                            nk_sum_flags_d);
}


__global__ void kernel_set_pair_distribution_k_index_only(int n_pairs, int* k_index_d, int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int nk;
  if (tid < n_pairs) {
    if (nk_sum_flags_d[tid] > 0) {
      nk = nk_sum_flags_d[tid] - 1;
      if (tid == n_pairs - 1) {
        nk = nk_sum_flags_d[tid - 1];
      }
      k_index_d[nk] = tid;
    }
  }
}

extern "C" void gpu_set_pair_distribution_k_index_only(int n_pairs, int* k_index_d, int* nk_sum_flags_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_set_pair_distribution_k_index_only<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, k_index_d, nk_sum_flags_d);
}


__global__ void kernel_set_pair_distribution_j2_only(int n_pairs, int n_sites0, int* neighbors_list_d, int* j2_index_d,
                                                     int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int j2, nk;
  if (tid < n_pairs) {
    if (nk_sum_flags_d[tid] > 0) {
      nk = nk_sum_flags_d[tid] - 1;
      if (tid == n_pairs - 1) {
        nk = nk_sum_flags_d[tid - 1];
      }
      j2 = ((neighbors_list_d[tid] - 1) % n_sites0);
      j2_index_d[nk] = j2 + 1;
    }
  }
}

extern "C" void gpu_set_pair_distribution_j2_only(int n_pairs, int n_sites0, int* neighbors_list_d, int* j2_index_d,
                                                  int* nk_sum_flags_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_set_pair_distribution_j2_only<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, n_sites0, neighbors_list_d, j2_index_d,
                                                                            nk_sum_flags_d);
}


__global__ void kernel_set_pair_distribution_rjs_only(int n_pairs, double* rjs, double* rjs_index_d, int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int nk;
  if (tid < n_pairs) {
    if (nk_sum_flags_d[tid] > 0) {
      nk = nk_sum_flags_d[tid] - 1;
      if (tid == n_pairs - 1) {
        nk = nk_sum_flags_d[tid - 1];
      }
      rjs_index_d[nk] = rjs[tid];
    }
  }
}

extern "C" void gpu_set_pair_distribution_rjs_only(int n_pairs, double* rjs, double* rjs_index_d, int* nk_sum_flags_d,
                                                   hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_set_pair_distribution_rjs_only<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, rjs, rjs_index_d, nk_sum_flags_d);
}


__global__ void kernel_set_pair_distribution_xyz_only(int n_pairs, double* xyz, double* xyz_index_d, int* nk_sum_flags_d) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int nk;
  if (tid < n_pairs) {
    if (nk_sum_flags_d[tid] > 0) {
      nk = nk_sum_flags_d[tid] - 1;
      if (tid == n_pairs - 1) {
        nk = nk_sum_flags_d[tid - 1];
      }
      xyz_index_d[3 * nk] = xyz[3 * tid];
      xyz_index_d[3 * nk + 1] = xyz[3 * tid + 1];
      xyz_index_d[3 * nk + 2] = xyz[3 * tid + 2];
    }
  }
}

extern "C" void gpu_set_pair_distribution_xyz_only(int n_pairs, double* xyz, double* xyz_index_d, int* nk_sum_flags_d,
                                                   hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  kernel_set_pair_distribution_xyz_only<<<nblocks, nthreads, 0, stream[0]>>>(n_pairs, xyz, xyz_index_d, nk_sum_flags_d);
}


// ----------------------------------------------------- //
// --- Pair distribution calculation and derivatives --- //
// ----------------------------------------------------- //

__global__ void kernel_get_pair_distribution_kde(double* pdf_out, double* pdf_der_out, int n_k, int n_samples, double kde_sigma,
                                                 double* x_d, double* dV_d, double* rjs_d, double pdf_factor, double der_factor) {
  int utid = threadIdx.x + blockIdx.x * blockDim.x;
  // int k_index = blockIdx.x * blockDim.x + threadIdx.x;  // Index in the j dimension (columns)
  // int l_index = blockIdx.y * blockDim.y + threadIdx.y;  // Index in the l dimension (rows)

  int i;
  double x, dV, r, x_gauss, pdf_temp, pdf_der_temp;

  int k_index = utid / n_samples;
  int l_index = utid % n_samples;

  if (k_index < n_k && l_index < n_samples) {
    x = x_d[l_index];
    dV = dV_d[l_index];

    r = rjs_d[k_index];

    x_gauss = ((x - r) / kde_sigma);

    pdf_temp = pdf_factor * exp(-0.5 * x_gauss * x_gauss) / dV;

    pdf_der_temp = (der_factor * pdf_temp * x_gauss) / (kde_sigma * r); // / dV;

    // Reversed the indexing of this temporary array for memory coalescence
    pdf_out[k_index + l_index * n_k] = pdf_temp;

    pdf_der_out[l_index + n_samples * k_index] = pdf_der_temp;
  }
}

__global__ void kernel_reduce_pair_distribution(double* pdf_in, double* pdf_out, int n_k, int n_samples) {
  int tid = threadIdx.x;
  int i, stride;
  double pdf_temp;

  // Now we want to reduce over the n_samples portions of continuous memory
  int n_strides = (n_k + BLOCK_SIZE - 1) / BLOCK_SIZE;

  __shared__ double sharedData[BLOCK_SIZE];


  pdf_temp = 0.0;

  for (i = 0; i < n_strides; ++i) {
    long long k_index = tid + i * BLOCK_SIZE;
    long long l_index = blockIdx.x;


    if (k_index < n_k && l_index < n_samples) {
      // Remember the reversed indexing!
      pdf_temp += pdf_in[k_index + l_index * n_k];
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
    pdf_out[blockIdx.x] = sharedData[0];

    //    i = blockIdx.x;
    //    printf("> after pdf reduce, blockidx = %d, pdf_out = %lf\n", i, pdf_out[i] );
  }
}


extern "C" void gpu_get_pair_distribution_and_ders(double* pair_distribution_d, double* pair_distribution_der_d, int n_k,
                                                   int n_samples, double kde_sigma, double* x_d, double* dV_d, double* rjs_d,
                                                   double pdf_factor, double der_factor, hipStream_t* stream) {
  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  size_t freeMem, totalMem;

  // Get memory information
  hipMemGetInfo(&freeMem, &totalMem);

  // Calculate used memory
  size_t usedMem = totalMem - freeMem;

  // Print out the memory information
  // printf("\n--- b4 pdf Total Memory: %lu bytes\n", totalMem);
  // printf("--- b4 pdf Free Memory: %lu bytes\n", freeMem);
  // printf("--- b4 pdf Used Memory: %lu bytes\n", usedMem);


  double* pdf_to_reduce;
  hipMallocAsync(&pdf_to_reduce, n_k * n_samples * sizeof(double), stream[0]);

  // Get memory information
  hipMemGetInfo(&freeMem, &totalMem);
  usedMem = totalMem - freeMem;

  // Print out the memory information
  // printf("\n--- after alloc pdf Total Memory: %lu bytes\n", totalMem);
  // printf("--- after alloc pdf Free Memory: %lu bytes\n", freeMem);
  // printf("--- after alloc pdf Used Memory: %lu bytes\n", usedMem);


  // printf("> pdf evaluation kernel starting\n");
  kernel_get_pair_distribution_kde<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_der_d, n_k, n_samples,
                                                                        kde_sigma, x_d, dV_d, rjs_d, pdf_factor, der_factor);
  //  hipDeviceSynchronize();
  // gpuErrchk( hipPeekAtLastError() );
  // printf("> pdf evaluation kernel finished\n");


  // Get memory information
  hipMemGetInfo(&freeMem, &totalMem);
  usedMem = totalMem - freeMem;

  // Print out the memory information
  // printf("\n--- after kernel pdf Total Memory: %lu bytes\n", totalMem);
  // printf("--- after kernel pdf Free Memory: %lu bytes\n", freeMem);
  // printf("--- after kernel pdf Used Memory: %lu bytes\n", usedMem);


  // Then we need to reduce over the size of the blocks for the pair distribution

  //  // More complex reduction
  //  //----------------------------------------
  nblocks = dim3(n_samples, 1, 1);
  nthreads = dim3(threads, 1, 1);
  // printf("> pdf reduction kernel starting\n");
  kernel_reduce_pair_distribution<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_d, n_k, n_samples);
  // printf("> pdf reduction kernel finished\n");
  //  //----------------------------------------


  hipFreeAsync(pdf_to_reduce, stream[0]);
  //  hipDeviceSynchronize();
  //gpuErrchk( hipPeekAtLastError() );
}


__global__ void kernel_get_pair_distribution_kde_der_only(double* pdf_der_out, int n_k, int n_samples, double kde_sigma,
                                                          double* x_d, double* dV_d, double* rjs_d, double pdf_factor,
                                                          double der_factor) {
  long long utid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  // int k_index = blockIdx.x * blockDim.x + threadIdx.x;  // Index in the j dimension (columns)
  // int l_index = blockIdx.y * blockDim.y + threadIdx.y;  // Index in the l dimension (rows)

  int i;
  double x, dV, r, x_gauss, pdf_temp, pdf_der_temp;

  long long k_index = utid / n_samples;
  long long l_index = utid % n_samples;

  if (k_index < n_k && l_index < n_samples) {
    x = x_d[l_index];
    dV = dV_d[l_index];

    r = rjs_d[k_index];

    x_gauss = ((x - r) / kde_sigma);

    pdf_temp = pdf_factor * exp(-0.5 * x_gauss * x_gauss) / dV;

    pdf_der_temp = (der_factor * pdf_temp * x_gauss) / (kde_sigma * r); // / dV;


    // if( isnan( pdf_der_temp ) ){
    //   // printf("pdf_der_temp is NaN! l = %d, k = %d, n_k = %d, r %lf, kde_sigma %lf\n", l_index, k_index, n_k,  r, kde_sigma);
    // }

    // Reversed the indexing of this temporary array for memory coalescence
    //    pdf_out[  k_index + l_index * n_k    ] = pdf_temp;

    pdf_der_out[l_index + n_samples * k_index] = pdf_der_temp;
  }
}


extern "C" void gpu_get_pair_distribution_der_only(double* pair_distribution_der_d, int n_k, int n_samples, double kde_sigma,
                                                   double* x_d, double* dV_d, double* rjs_d, double pdf_factor, double der_factor,
                                                   hipStream_t* stream) {
  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  //  // printf("> pdf der evaluation kernel starting\n");
  kernel_get_pair_distribution_kde_der_only<<<nblocks, nthreads, 0, stream[0]>>>(pair_distribution_der_d, n_k, n_samples, kde_sigma,
                                                                                 x_d, dV_d, rjs_d, pdf_factor, der_factor);
}


__global__ void kernel_get_pair_distribution_kde_only(double* pdf_out, int n_k, int n_samples, double kde_sigma, double* x_d,
                                                      double* dV_d, double* rjs_d, double pdf_factor, double der_factor) {
  long long utid = (long long) blockIdx.x * blockDim.x + threadIdx.x;
  // int k_index = blockIdx.x * blockDim.x + threadIdx.x;  // Index in the j dimension (columns)
  // int l_index = blockIdx.y * blockDim.y + threadIdx.y;  // Index in the l dimension (rows)

  int i;
  double x, dV, r, x_gauss, pdf_temp, pdf_der_temp;

  long long k_index = utid / n_samples;
  long long l_index = utid % n_samples;

  if (k_index < n_k && l_index < n_samples) {
    x = x_d[l_index];
    dV = dV_d[l_index];

    r = rjs_d[k_index];

    x_gauss = ((x - r) / kde_sigma);

    //    // printf(" pdf only kernel r = %lf, x = %lf, x_gauss = %lf\n", r, x, x_gauss);

    pdf_temp = pdf_factor * exp(-0.5 * x_gauss * x_gauss) / dV;

    // x =  (- 0.5 * x_gauss * x_gauss) ;

    //    if( pdf_temp > 0.000 ){
    // // printf("> checking calc has data, k_index %d, l_index %d, pdf_factor = %lf, x = %lf, pdf_temp = %lf \n", k_index, l_index, pdf_factor, x, pdf_temp);
    //      }
    // //    pdf_der_temp = (der_factor * pdf_temp * x_gauss) / (kde_sigma * r);// / dV;

    // Reversed the indexing of this temporary array for memory coalescence
    pdf_out[k_index + l_index * n_k] = pdf_temp;

    //pdf_der_out[ l_index + n_samples * k_index ] = pdf_der_temp;
  }
}


// An empty species pair contributes nothing, and asking the device for it is an
// error rather than a no-op.
//
// n_k is the number of pairs of one species combination inside one batch, so it
// is zero whenever a batch happens to contain none -- which cannot happen with
// the default gpu_n_batches = 1 over a whole cell, and happens readily as soon
// as the cell is split. The grid below is then dim3(0,1,1), and CUDA rejects a
// zero-extent launch with "invalid configuration argument". The allocation is
// zero bytes too, so rjs_d is NULL.
//
// This is why XRD_mad aborts inside gpu_get_pair_distribution_only below with
// `gpu_n_batches = 4` (it was gpu_exp.cu:1397 before the src/gpu split): it is
// not a stream problem and not a memory problem, it is the batched experimental
// path never having been run with more than one batch. estat_gsf, the only case
// in the suite that sets gpu_n_batches > 1, exhausts device memory before it
// reaches here, so nothing exercised this.
//
// Returning early is the whole fix: the caller memsets pair_distribution_d
// before calling, and the answer for no pairs is zero.
static inline bool pdf_batch_is_empty(int n_k, int n_samples) {
  return n_k <= 0 || n_samples <= 0;
}

extern "C" void gpu_get_pair_distribution_only(double* pair_distribution_d, int n_k, int n_samples, double kde_sigma, double* x_d,
                                               double* dV_d, double* rjs_d, double pdf_factor, double der_factor,
                                               hipStream_t* stream) {
  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  if (pdf_batch_is_empty(n_k, n_samples)) {
    if (n_samples > 0)
      gpuErrchk(hipMemsetAsync(pair_distribution_d, 0, n_samples * sizeof(double), stream[0]));
    return;
  }

  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  // gpuErrchk( hipPeekAtLastError() );

  // size_t freeMem, totalMem;

  // // Get memory information
  // hipMemGetInfo(&freeMem, &totalMem);

  // // Calculate used memory
  // size_t usedMem = totalMem - freeMem;

  // // Print out the memory information
  // // // printf("\n--- b4 pdf Total Memory: %lu bytes\n", totalMem);
  // // // printf("--- b4 pdf Free Memory: %lu bytes\n", freeMem);
  // // // printf("--- b4 pdf Used Memory: %lu bytes\n", usedMem);


  double* pdf_to_reduce;
  hipMallocAsync(&pdf_to_reduce, n_k * n_samples * sizeof(double), stream[0]);
  // gpuErrchk( hipPeekAtLastError() );

  // Get memory information
  // hipMemGetInfo(&freeMem, &totalMem);
  // usedMem = totalMem - freeMem;

  // // Print out the memory information
  // // // printf("\n--- after alloc pdf Total Memory: %lu bytes\n", totalMem);
  // // // printf("--- after alloc pdf Free Memory: %lu bytes\n", freeMem);
  // // // printf("--- after alloc pdf Used Memory: %lu bytes\n", usedMem);


  //  // printf("> pdf evaluation kernel starting\n");
  kernel_get_pair_distribution_kde_only<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, n_k, n_samples, kde_sigma, x_d, dV_d,
                                                                             rjs_d, pdf_factor, der_factor);
  //  hipDeviceSynchronize();
  // gpuErrchk( hipPeekAtLastError() );
  // printf("> pdf evaluation kernel finished\n");
  //   // Get memory information
  //   hipMemGetInfo(&freeMem, &totalMem);
  //   usedMem = totalMem - freeMem;

  //   // Print out the memory information
  //   printf("\n--- after kernel pdf Total Memory: %lu bytes\n", totalMem);
  //   printf("--- after kernel pdf Free Memory: %lu bytes\n", freeMem);
  //   printf("--- after kernel pdf Used Memory: %lu bytes\n", usedMem);

  // Then we need to reduce over the size of the blocks for the pair distribution

  //  // More complex reduction
  //  //----------------------------------------
  nblocks = dim3(n_samples, 1, 1);
  nthreads = dim3(threads, 1, 1);
  //  printf("> pdf reduction kernel starting\n");
  kernel_reduce_pair_distribution<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_d, n_k, n_samples);
  //printf("> pdf reduction kernel finished\n");
  //  //----------------------------------------


  hipFreeAsync(pdf_to_reduce, stream[0]);
  //  hipDeviceSynchronize();
  //gpuErrchk( hipPeekAtLastError() );
}


extern "C" void gpu_get_pair_distribution_only_falloc(double* pair_distribution_d, double* pdf_to_reduce, int n_k, int n_samples,
                                                      double kde_sigma, double* x_d, double* dV_d, double* rjs_d, double pdf_factor,
                                                      double der_factor, hipStream_t* stream) {
  // We want to evaluate the smoothed pair distribution in a quick manner using these kernels
  // 1. Evaluate the exponential over the threads

  // See pdf_batch_is_empty above: with the cell split into batches a species
  // pair can be absent from a batch, and dim3(0,1,1) is a launch error.
  if (pdf_batch_is_empty(n_k, n_samples)) {
    if (n_samples > 0)
      gpuErrchk(hipMemsetAsync(pair_distribution_d, 0, n_samples * sizeof(double), stream[0]));
    return;
  }

  int threads = BLOCK_SIZE;
  dim3 nblocks = dim3((unsigned int) (((long long) n_k * n_samples + threads - 1) / threads), 1, 1);
  dim3 nthreads = dim3(threads, 1, 1);

  // gpuErrchk( hipPeekAtLastError() );
  //  printf("> pdf evaluation kernel starting\n");
  gpuErrchk(hipPeekAtLastError());
  kernel_get_pair_distribution_kde_only<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, n_k, n_samples, kde_sigma, x_d, dV_d,
                                                                             rjs_d, pdf_factor, der_factor);
  //printf("> pdf evaluation kernel finished\n");
  //  gpuErrchk( hipPeekAtLastError() );
  hipDeviceSynchronize();
  gpuErrchk(hipPeekAtLastError());
  nblocks = dim3(n_samples, 1, 1);
  nthreads = dim3(threads, 1, 1);

  //printf("> pdf reduction kernel starting\n");
  kernel_reduce_pair_distribution<<<nblocks, nthreads, 0, stream[0]>>>(pdf_to_reduce, pair_distribution_d, n_k, n_samples);
  hipDeviceSynchronize();
  gpuErrchk(hipPeekAtLastError());
  //printf("> pdf reduction kernel finished\n");
  //----------------------------------------
  //  gpuErrchk( hipPeekAtLastError() );
}
