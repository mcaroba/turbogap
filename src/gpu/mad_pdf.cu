// MAD: the pair distribution function.
//
// Counting the pairs inside the pdf cutoff, compacting them into the k-indexed
// buffers, the kernel density estimate of g(r) and its derivative, and the
// reduction over pairs.
#include "gpu_backend.h"
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
  // Report and clear anything left over from earlier in the stream.
  // hipPeekAtLastError does not reset the flag, so without this the first
  // failure anywhere is re-reported at every later check and every one of them
  // blames its own line -- which is how a launch that was fine came to be
  // accused of "invalid configuration argument" on every rank but zero.
  // Clearing without saying so would instead throw away the evidence that
  // something earlier failed, so it is printed here and named as inherited.
  {
    hipError_t stale = hipGetLastError();
    if (stale != hipSuccess) {
      fprintf(stderr,
              "GPUlaunch: gpu_get_pair_distribution_nk entered with an error already set: %s.\n"
              "GPUlaunch:   It came from something earlier in this stream that did not check.\n",
              hipGetErrorString(stale));
      fflush(stderr);
    }
  }
  kernel_get_pair_distribution_nk<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, n_sites0, neighbors_list, n_neigh,
                                                                       neighbor_species, species, rjs, xyz, r_min, r_max, r_cut,
                                                                       buffer, nk_flags_d, species_1, species_2);
  {
    hipError_t launch_err = hipGetLastError();
    if (launch_err != hipSuccess) {
      fprintf(stderr,
              "GPUlaunch: kernel_get_pair_distribution_nk grid(%u,%u,%u) block(%u,%u,%u) i_beg=%d i_end=%d n_pairs=%d: %s\n",
              nblocks.x, nblocks.y, nblocks.z, nthreads.x, nthreads.y, nthreads.z, i_beg, i_end, n_pairs,
              hipGetErrorString(launch_err));
      fflush(stderr);
    }
  }

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


extern "C" void gpu_set_pair_distribution_k_index_only(int n_pairs, int* k_index_d, int* nk_sum_flags_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_pdf_k_index", n_pairs, stream,
      TG_LAMBDA(const int tid) {
        if (nk_sum_flags_d[tid] <= 0)
          return;
        int nk = nk_sum_flags_d[tid] - 1;
        if (tid == n_pairs - 1) {
          nk = nk_sum_flags_d[tid - 1];
        }
        k_index_d[nk] = tid;
      },
      tpb);
}


extern "C" void gpu_set_pair_distribution_j2_only(int n_pairs, int n_sites0, int* neighbors_list_d, int* j2_index_d,
                                                  int* nk_sum_flags_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_pdf_j2", n_pairs, stream,
      TG_LAMBDA(const int tid) {
        if (nk_sum_flags_d[tid] <= 0)
          return;
        int nk = nk_sum_flags_d[tid] - 1;
        if (tid == n_pairs - 1) {
          nk = nk_sum_flags_d[tid - 1];
        }
        j2_index_d[nk] = ((neighbors_list_d[tid] - 1) % n_sites0) + 1;
      },
      tpb);
}


extern "C" void gpu_set_pair_distribution_rjs_only(int n_pairs, double* rjs, double* rjs_index_d, int* nk_sum_flags_d,
                                                   hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_pdf_rjs", n_pairs, stream,
      TG_LAMBDA(const int tid) {
        if (nk_sum_flags_d[tid] <= 0)
          return;
        int nk = nk_sum_flags_d[tid] - 1;
        if (tid == n_pairs - 1) {
          nk = nk_sum_flags_d[tid - 1];
        }
        rjs_index_d[nk] = rjs[tid];
      },
      tpb);
}


extern "C" void gpu_set_pair_distribution_xyz_only(int n_pairs, double* xyz, double* xyz_index_d, int* nk_sum_flags_d,
                                                   hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_pdf_xyz", n_pairs, stream,
      TG_LAMBDA(const int tid) {
        if (nk_sum_flags_d[tid] <= 0)
          return;
        int nk = nk_sum_flags_d[tid] - 1;
        if (tid == n_pairs - 1) {
          nk = nk_sum_flags_d[tid - 1];
        }
        xyz_index_d[3 * nk] = xyz[3 * tid];
        xyz_index_d[3 * nk + 1] = xyz[3 * tid + 1];
        xyz_index_d[3 * nk + 2] = xyz[3 * tid + 2];
      },
      tpb);
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

// NOT converted to the backend layer, deliberately. This is a __shared__ block
// reduction over DOUBLES, and Kokkos would sum the block in a different order.
// Reassociating a floating-point sum moves the result, so converting it would
// be a physics change wearing a refactor's clothes. It needs a TeamPolicy that
// reproduces this tree order, worth doing separately and checking on its own.
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


// ------------------------------------------------------------------------- //
// --- Windowed pdf: every species pair at once, no n_k x n_samples array --- //
// ------------------------------------------------------------------------- //
//
// A pair's Gaussian is evaluated only where the host evaluates it, at the
// samples lo..hi within kde_cut sigma (pdf_window_lo/hi copy the host's index
// arithmetic in pair_kde_row), so the device sums the same terms as the host.
// The pairs of all species pairs are tagged with their channel and compacted
// once into one list, ordered by channel and, within a channel, by pair index;
// the list stays on the device for the force pass.

// chan_d(k) for every entry of the batch's neighbour list: the pair's channel,
// or -1 for a site's own entry and for a pair outside [r_lo, r_hi].
extern "C" void gpu_pdf_tag_pairs(int n_sites_b, int n_pairs_b, int* n_neigh_d, int* species_d, int* neighbor_species_d,
                                  double* rjs_d, int n_species, double r_lo, double r_hi, int* chan_d, hipStream_t* stream) {
  if (n_sites_b <= 0 || n_pairs_b <= 0)
    return;
  int* site_end;
  gpuErrchk(hipMallocAsync(&site_end, (size_t) n_sites_b * sizeof(int), stream[0]));
  gpuErrchk(hipMemcpyAsync(site_end, n_neigh_d, (size_t) n_sites_b * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));
  inclusiveScan(site_end, n_sites_b, stream);

  tg_parallel_for(
      "turbogap_pdf_tag_pairs", n_sites_b, stream,
      TG_LAMBDA(const int i) {
        const int end = site_end[i];
        const int beg = end - n_neigh_d[i];
        const int si = species_d[i];
        chan_d[beg] = -1;
        for (int k = beg + 1; k < end; k++) {
          const double r = rjs_d[k];
          chan_d[k] = (r >= r_lo && r <= r_hi) ? pdf_channel(si, neighbor_species_d[k], n_species) : -1;
        }
      },
      tpb);

  gpuErrchk(hipFreeAsync(site_end, stream[0]));
}

// Compacts the tagged pairs, channel by channel, into r_out, xyz_out (3 per
// pair), j2_out (1-based atom) and chan_out, each of capacity n_pairs_b.
// off_h(0:n_channels) receives each channel's first position on the host.
// seg_out(c * (n_sites_b + 1) + i) is where site i's pairs of channel c begin,
// and (c, n_sites_b) is the channel's end: within a channel the pairs keep
// neighbour-list order, so each site's are one run.
extern "C" void gpu_pdf_compact_pairs(int n_pairs_b, int n_sites0, int n_channels, int* chan_d, int* neighbors_list_d,
                                      double* rjs_d, double* xyz_d, double* r_out, double* xyz_out, int* j2_out, int* chan_out,
                                      int* off_h, int n_sites_b, int* n_neigh_d, int* seg_out, hipStream_t* stream) {
  for (int c = 0; c <= n_channels; c++)
    off_h[c] = 0;
  if (n_pairs_b <= 0 || n_channels <= 0)
    return;

  const size_t n = (size_t) n_pairs_b;
  int* rank;
  int* off_d;
  gpuErrchk(hipMallocAsync(&rank, n * n_channels * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&off_d, (size_t) (n_channels + 1) * sizeof(int), stream[0]));

  for (int c = 0; c < n_channels; c++) {
    int* rank_c = rank + c * n;
    tg_parallel_for(
        "turbogap_pdf_channel_flags", n_pairs_b, stream, TG_LAMBDA(const int k) { rank_c[k] = (chan_d[k] == c) ? 1 : 0; }, tpb);
    inclusiveScan(rank_c, n_pairs_b, stream);
  }

  tg_parallel_for(
      "turbogap_pdf_channel_offsets", 1, stream,
      TG_LAMBDA(const int) {
        off_d[0] = 0;
        for (int c = 0; c < n_channels; c++)
          off_d[c + 1] = off_d[c] + rank[c * n + n - 1];
      },
      1);

  tg_parallel_for(
      "turbogap_pdf_compact", n_pairs_b, stream,
      TG_LAMBDA(const int k) {
        const int c = chan_d[k];
        if (c < 0)
          return;
        const int p = off_d[c] + rank[c * n + k] - 1;
        r_out[p] = rjs_d[k];
        xyz_out[3 * (size_t) p] = xyz_d[3 * (size_t) k];
        xyz_out[3 * (size_t) p + 1] = xyz_d[3 * (size_t) k + 1];
        xyz_out[3 * (size_t) p + 2] = xyz_d[3 * (size_t) k + 2];
        j2_out[p] = ((neighbors_list_d[k] - 1) % n_sites0) + 1;
        chan_out[p] = c;
      },
      tpb);

  // A site's run in channel c starts after the channel's pairs that come
  // before the site's own entry in the list: rank, an inclusive count, at the
  // entry before it.
  int* site_end;
  gpuErrchk(hipMallocAsync(&site_end, (size_t) n_sites_b * sizeof(int), stream[0]));
  gpuErrchk(hipMemcpyAsync(site_end, n_neigh_d, (size_t) n_sites_b * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));
  inclusiveScan(site_end, n_sites_b, stream);
  const int n_seg = n_sites_b + 1;
  tg_parallel_for(
      "turbogap_pdf_site_segments", n_channels * n_seg, stream,
      TG_LAMBDA(const int idx) {
        const int c = idx / n_seg;
        const int i = idx - c * n_seg;
        if (i == n_sites_b) {
          seg_out[idx] = off_d[c + 1];
          return;
        }
        const int beg = site_end[i] - n_neigh_d[i];
        seg_out[idx] = off_d[c] + ((beg > 0) ? rank[c * n + beg - 1] : 0);
      },
      tpb);

  gpuErrchk(hipMemcpyAsync(off_h, off_d, (size_t) (n_channels + 1) * sizeof(int), hipMemcpyDeviceToHost, stream[0]));
  gpuErrchk(hipStreamSynchronize(stream[0]));
  gpuErrchk(hipFreeAsync(site_end, stream[0]));
  gpuErrchk(hipFreeAsync(rank, stream[0]));
  gpuErrchk(hipFreeAsync(off_d, stream[0]));
}

// Pairs per block of the gather. Their r and windows sit in shared memory.
#define PDF_GATHER_CHUNK 1024

// NOT converted to the backend layer, deliberately: it stages a chunk of pairs
// in __shared__ memory for the whole block, which needs a TeamPolicy the layer
// does not carry.
//
// The chunk's pairs are bucketed by the first sample of their window: a count
// per bucket, the buckets' starts, each pair dropped into its bucket (the slot
// from an atomic, so in any order) and each bucket then sorted by pair index,
// which takes the scheduler back out -- a bucket holds a pair or two, as the
// chunk has fewer pairs than the pdf samples. A sample's pairs are then the
// buckets of the w_max samples up to it, one run of the sorted list, and one
// thread per sample sums them in that order. Scanning the whole chunk for
// each sample instead checked ~13x more pairs than it summed at 25 A.
// Dynamic shared memory: 2 n_samples + 1 ints.
__global__ void kernel_pdf_window_gather(const double* r_d, int n_k, int n_samples, double r_min, double dr, double kde_sigma,
                                         double kde_cut, int w_max, const double* x_d, double* partial_d) {
  __shared__ double r_s[PDF_GATHER_CHUNK];
  __shared__ int lo_s[PDF_GATHER_CHUNK];
  __shared__ int hi_s[PDF_GATHER_CHUNK];
  __shared__ int order_s[PDF_GATHER_CHUNK];
  extern __shared__ int bucket_s[];
  int* start_s = bucket_s;
  int* cursor_s = bucket_s + n_samples + 1;

  const int beg = blockIdx.x * PDF_GATHER_CHUNK;
  const int n = (n_k - beg < PDF_GATHER_CHUNK) ? n_k - beg : PDF_GATHER_CHUNK;

  for (int l = threadIdx.x; l <= n_samples; l += blockDim.x)
    start_s[l] = 0;
  for (int l = threadIdx.x; l < n_samples; l += blockDim.x)
    cursor_s[l] = 0;
  __syncthreads();

  for (int p = threadIdx.x; p < n; p += blockDim.x) {
    const double r = r_d[beg + p];
    const int lo = pdf_window_lo(r, r_min, dr, kde_sigma, kde_cut);
    r_s[p] = r;
    lo_s[p] = lo < n_samples - 1 ? lo : n_samples - 1;
    hi_s[p] = pdf_window_hi(r, r_min, dr, kde_sigma, kde_cut, n_samples);
    atomicAdd(&start_s[lo_s[p] + 1], 1);
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    for (int l = 1; l <= n_samples; l++)
      start_s[l] += start_s[l - 1];
  }
  __syncthreads();

  for (int p = threadIdx.x; p < n; p += blockDim.x) {
    const int b = lo_s[p];
    order_s[start_s[b] + atomicAdd(&cursor_s[b], 1)] = p;
  }
  __syncthreads();

  for (int b = threadIdx.x; b < n_samples; b += blockDim.x) {
    for (int s = start_s[b] + 1; s < start_s[b + 1]; s++) {
      const int v = order_s[s];
      int t = s - 1;
      while (t >= start_s[b] && order_s[t] > v) {
        order_s[t + 1] = order_s[t];
        t--;
      }
      order_s[t + 1] = v;
    }
  }
  __syncthreads();

  for (int l = threadIdx.x; l < n_samples; l += blockDim.x) {
    const double x = x_d[l];
    const int b0 = (l - w_max + 1 > 0) ? l - w_max + 1 : 0;
    double acc = 0.0;
    for (int s = start_s[b0]; s < start_s[l + 1]; s++) {
      const int p = order_s[s];
      if (hi_s[p] < l)
        continue;
      const double g = (x - r_s[p]) / kde_sigma;
      acc += exp(-0.5 * g * g);
    }
    partial_d[(size_t) blockIdx.x * n_samples + l] = acc;
  }
}

// pdf_out(l, c) = pdf_factor / dV(l) * sum over channel c's pairs of
// exp(-(x_l - r)^2 / 2 sigma^2), for the compacted list r_d with channel
// offsets off_h(0:n_channels). Blocks of pairs are summed in block order.
extern "C" void gpu_pdf_window_kde(int n_channels, int* off_h, int n_samples, double r_min, double dr, double kde_sigma,
                                   double kde_cut, double* x_d, double* dV_d, double* r_d, double pdf_factor, double* pdf_out_d,
                                   hipStream_t* stream) {
  if (n_channels <= 0 || n_samples <= 0)
    return;

  int* blk_h = (int*) malloc((size_t) (n_channels + 1) * sizeof(int));
  blk_h[0] = 0;
  for (int c = 0; c < n_channels; c++) {
    const int n_c = off_h[c + 1] - off_h[c];
    blk_h[c + 1] = blk_h[c] + (n_c + PDF_GATHER_CHUNK - 1) / PDF_GATHER_CHUNK;
  }
  const int n_blk = blk_h[n_channels];

  int* blk_d;
  double* partial = nullptr;
  gpuErrchk(hipMallocAsync(&blk_d, (size_t) (n_channels + 1) * sizeof(int), stream[0]));
  gpuErrchk(hipMemcpyAsync(blk_d, blk_h, (size_t) (n_channels + 1) * sizeof(int), hipMemcpyHostToDevice, stream[0]));
  if (n_blk > 0)
    gpuErrchk(hipMallocAsync(&partial, (size_t) n_blk * n_samples * sizeof(double), stream[0]));

  int threads = ((n_samples < 1024 ? n_samples : 1024) + 31) / 32 * 32;
  // The widest window, as the host's accumulate_pair_kde sizes it.
  const int w_max = (int) (2.0 * kde_cut * kde_sigma / dr) + 5;
  const size_t shared = (size_t) (2 * n_samples + 1) * sizeof(int);
  for (int c = 0; c < n_channels; c++) {
    const int n_c = off_h[c + 1] - off_h[c];
    const int nb = blk_h[c + 1] - blk_h[c];
    if (nb == 0)
      continue;
    kernel_pdf_window_gather<<<nb, threads, shared, stream[0]>>>(r_d + off_h[c], n_c, n_samples, r_min, dr, kde_sigma, kde_cut,
                                                                 w_max, x_d, partial + (size_t) blk_h[c] * n_samples);
  }

  tg_parallel_for(
      "turbogap_pdf_window_sum_blocks", n_channels * n_samples, stream,
      TG_LAMBDA(const int idx) {
        const int c = idx / n_samples;
        const int l = idx - c * n_samples;
        double s = 0.0;
        for (int b = blk_d[c]; b < blk_d[c + 1]; b++)
          s += partial[(size_t) b * n_samples + l];
        pdf_out_d[idx] = pdf_factor * s / dV_d[l];
      },
      tpb);

  gpuErrchk(hipStreamSynchronize(stream[0]));
  if (partial != nullptr)
    gpuErrchk(hipFreeAsync(partial, stream[0]));
  gpuErrchk(hipFreeAsync(blk_d, stream[0]));
  free(blk_h);
}

// ---------------------------------------------------------------------------
// The pdf's pairs built on the device from the positions, rather than uploaded
// from the host's neighbour list every step. For an orthorhombic cell at least
// three cutoffs wide: a cell list whose cells are at least rcut / 2 wide, the
// atoms copied into cell order and, within a cell, index order; then each
// site's neighbours within rcut among the 5 x 5 x 5 cells around its own,
// nearest image. One warp per site walks those cells' atoms 32 at a time, so
// the loads of a cell are one sweep; a ballot counts the hits and places them
// in the order walked, which is fixed. The output is a batch's list in the
// host's layout -- each site's own entry, then its neighbours -- so everything
// that read the uploaded list reads it unchanged.
//
// NOT converted to the backend layer, deliberately: the warp-wide ballot is a
// team pattern the layer does not carry. The cell list itself is flat maps.

#define PDF_CELL_REACH 2
#define PDF_WARPS 8

struct pdf_cells_t {
  int n;
  int nc[3];
  double L[3];
  double rcut;
  double* x_d;
  double* y_d;
  double* z_d;
  int* atom_d;
  int* species_d;
  int* species_idx_d;
  int* cell_start_d;
  int* cell_of_site_d;
};

TG_INLINE_FUNCTION int pdf_cell_of(double x, double L, int nc) {
  double s = x / L;
  s -= floor(s);
  const int c = (int) (s * nc);
  return c < nc ? c : nc - 1;
}

// The cell list of n atoms (pos_h 3 x n, species_h 1-based) in an lx x ly x lz
// cell. The caller has checked it is at least three rcut a side, so there are
// at least six cells a side and the 125 around a cell are distinct.
extern "C" void* gpu_pdf_cells_build(int n, double* pos_h, int* species_h, double lx, double ly, double lz, double rcut,
                                     hipStream_t* stream) {
  pdf_cells_t* cl = new pdf_cells_t;
  cl->n = n;
  cl->L[0] = lx;
  cl->L[1] = ly;
  cl->L[2] = lz;
  cl->rcut = rcut;
  for (int d = 0; d < 3; d++)
    cl->nc[d] = (int) (cl->L[d] / (rcut / PDF_CELL_REACH));
  const int nc0 = cl->nc[0];
  const int nc1 = cl->nc[1];
  const int nc2 = cl->nc[2];
  const int ncell = nc0 * nc1 * nc2;
  const size_t nn = (size_t) (n > 0 ? n : 1);

  double* pos;
  int* cell_of;
  int* cursor;
  int* order;
  gpuErrchk(hipMallocAsync(&pos, 3 * nn * sizeof(double), stream[0]));
  gpuErrchk(hipMallocAsync(&cell_of, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cursor, (size_t) ncell * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&order, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->x_d, nn * sizeof(double), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->y_d, nn * sizeof(double), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->z_d, nn * sizeof(double), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->atom_d, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->species_d, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->species_idx_d, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->cell_start_d, (size_t) (ncell + 1) * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cl->cell_of_site_d, nn * sizeof(int), stream[0]));
  gpuErrchk(hipMemcpyAsync(pos, pos_h, 3 * (size_t) n * sizeof(double), hipMemcpyHostToDevice, stream[0]));
  gpuErrchk(hipMemcpyAsync(cl->species_idx_d, species_h, (size_t) n * sizeof(int), hipMemcpyHostToDevice, stream[0]));
  gpuErrchk(hipMemsetAsync(cl->cell_start_d, 0, (size_t) (ncell + 1) * sizeof(int), stream[0]));
  gpuErrchk(hipMemsetAsync(cursor, 0, (size_t) ncell * sizeof(int), stream[0]));

  int* cell_start = cl->cell_start_d;
  int* cell_of_site = cl->cell_of_site_d;
  double* x = cl->x_d;
  double* y = cl->y_d;
  double* z = cl->z_d;
  int* atom = cl->atom_d;
  int* spec = cl->species_d;
  const int* spec_in = cl->species_idx_d;
  const double L0 = lx;
  const double L1 = ly;
  const double L2 = lz;

  tg_parallel_for(
      "turbogap_pdf_cells_bin", n, stream,
      TG_LAMBDA(const int a) {
        const int c = pdf_cell_of(pos[3 * (size_t) a], L0, nc0) +
                      nc0 * (pdf_cell_of(pos[3 * (size_t) a + 1], L1, nc1) + nc1 * pdf_cell_of(pos[3 * (size_t) a + 2], L2, nc2));
        cell_of[a] = c;
        cell_of_site[a] = c;
        TG_ATOMIC_ADD(&cell_start[c + 1], 1);
      },
      tpb);
  inclusiveScan(cell_start + 1, ncell, stream);
  tg_parallel_for(
      "turbogap_pdf_cells_place", n, stream,
      TG_LAMBDA(const int a) {
        const int c = cell_of[a];
        order[cell_start[c] + TG_ATOMIC_FETCH_ADD(&cursor[c], 1)] = a;
      },
      tpb);
  // The atomics left each cell in any order; index order fixes it.
  tg_parallel_for(
      "turbogap_pdf_cells_order", ncell, stream,
      TG_LAMBDA(const int c) {
        for (int s = cell_start[c] + 1; s < cell_start[c + 1]; s++) {
          const int v = order[s];
          int t = s - 1;
          while (t >= cell_start[c] && order[t] > v) {
            order[t + 1] = order[t];
            t--;
          }
          order[t + 1] = v;
        }
      },
      tpb);
  tg_parallel_for(
      "turbogap_pdf_cells_gather", n, stream,
      TG_LAMBDA(const int t) {
        const int a = order[t];
        x[t] = pos[3 * (size_t) a];
        y[t] = pos[3 * (size_t) a + 1];
        z[t] = pos[3 * (size_t) a + 2];
        atom[t] = a;
        spec[t] = spec_in[a];
      },
      tpb);

  gpuErrchk(hipFreeAsync(order, stream[0]));
  gpuErrchk(hipFreeAsync(cursor, stream[0]));
  gpuErrchk(hipFreeAsync(cell_of, stream[0]));
  gpuErrchk(hipFreeAsync(pos, stream[0]));
  gpuErrchk(hipStreamSynchronize(stream[0]));
  return cl;
}

// Site i's neighbours within sqrt(rc2), one warp: counted, and written from
// entry k on when rjs is not null (1-based atom, species, distance,
// separation). The site's position comes from the cell-ordered arrays through
// its slot, found by scanning its own cell.
__device__ int pdf_warp_site(const int i, const pdf_cells_t cl, const double rc2, int k, int* nl, int* nspec, double* rjs,
                             double* xyz) {
  const int lane = threadIdx.x & 31;
  const unsigned lanes_below = (1u << lane) - 1u;
  const int ci = cl.cell_of_site_d[i];
  int slot = -1;
  for (int t0 = cl.cell_start_d[ci]; t0 < cl.cell_start_d[ci + 1]; t0 += 32) {
    const int t = t0 + lane;
    const unsigned found = __ballot_sync(0xffffffffu, t < cl.cell_start_d[ci + 1] && cl.atom_d[t] == i);
    if (found != 0u) {
      slot = t0 + __ffs(found) - 1;
      break;
    }
  }
  const double xi = cl.x_d[slot];
  const double yi = cl.y_d[slot];
  const double zi = cl.z_d[slot];
  const int nc0 = cl.nc[0];
  const int nc1 = cl.nc[1];
  const int nc2 = cl.nc[2];
  const int cx = ci % nc0;
  const int cy = (ci / nc0) % nc1;
  const int cz = ci / (nc0 * nc1);
  int n = 0;
  for (int dz = -PDF_CELL_REACH; dz <= PDF_CELL_REACH; dz++) {
    const int zc = (cz + dz + nc2) % nc2;
    for (int dy = -PDF_CELL_REACH; dy <= PDF_CELL_REACH; dy++) {
      const int yc = (cy + dy + nc1) % nc1;
      for (int dx = -PDF_CELL_REACH; dx <= PDF_CELL_REACH; dx++) {
        const int xc = (cx + dx + nc0) % nc0;
        const int c = xc + nc0 * (yc + nc1 * zc);
        const int beg = cl.cell_start_d[c];
        const int end = cl.cell_start_d[c + 1];
        for (int t0 = beg; t0 < end; t0 += 32) {
          const int t = t0 + lane;
          bool hit = false;
          double ax = 0.0;
          double ay = 0.0;
          double az = 0.0;
          double r2 = 0.0;
          if (t < end && t != slot) {
            ax = cl.x_d[t] - xi;
            ay = cl.y_d[t] - yi;
            az = cl.z_d[t] - zi;
            ax -= cl.L[0] * rint(ax / cl.L[0]);
            ay -= cl.L[1] * rint(ay / cl.L[1]);
            az -= cl.L[2] * rint(az / cl.L[2]);
            r2 = ax * ax + ay * ay + az * az;
            hit = r2 < rc2;
          }
          const unsigned hits = __ballot_sync(0xffffffffu, hit);
          if (rjs != nullptr && hit) {
            const int p = k + n + __popc(hits & lanes_below);
            nl[p] = cl.atom_d[t] + 1;
            nspec[p] = cl.species_d[t];
            rjs[p] = sqrt(r2);
            xyz[3 * (size_t) p] = ax;
            xyz[3 * (size_t) p + 1] = ay;
            xyz[3 * (size_t) p + 2] = az;
          }
          n += __popc(hits);
        }
      }
    }
  }
  return n;
}

__global__ void kernel_pdf_cells_count(const pdf_cells_t cl, int site_beg, int n_site, double rc2, int* counts) {
  const int s = blockIdx.x * PDF_WARPS + threadIdx.x / 32;
  if (s >= n_site)
    return;
  const int n = pdf_warp_site(site_beg + s, cl, rc2, 0, nullptr, nullptr, nullptr, nullptr);
  if ((threadIdx.x & 31) == 0)
    counts[s] = n;
}

__global__ void kernel_pdf_cells_fill(const pdf_cells_t cl, int site_beg, int n_site, double rc2, const int* n_neigh,
                                      const int* site_end, int* species_b, int* nl, int* nspec, double* rjs, double* xyz) {
  const int s = blockIdx.x * PDF_WARPS + threadIdx.x / 32;
  if (s >= n_site)
    return;
  const int i = site_beg + s;
  const int k = site_end[s] - n_neigh[s];
  if ((threadIdx.x & 31) == 0) {
    species_b[s] = cl.species_idx_d[i];
    nl[k] = i + 1;
    nspec[k] = cl.species_idx_d[i];
    rjs[k] = 0.0;
    xyz[3 * (size_t) k] = 0.0;
    xyz[3 * (size_t) k + 1] = 0.0;
    xyz[3 * (size_t) k + 2] = 0.0;
  }
  pdf_warp_site(i, cl, rc2, k + 1, nl, nspec, rjs, xyz);
}

// counts_h(s): the neighbours of site site_beg + s (0-based) within rcut.
extern "C" void gpu_pdf_cells_count(void* h, int site_beg, int n_site, int* counts_h, hipStream_t* stream) {
  if (n_site <= 0)
    return;
  const pdf_cells_t* cl = (const pdf_cells_t*) h;
  int* counts_d;
  gpuErrchk(hipMallocAsync(&counts_d, (size_t) n_site * sizeof(int), stream[0]));
  const int blocks = (n_site + PDF_WARPS - 1) / PDF_WARPS;
  kernel_pdf_cells_count<<<blocks, 32 * PDF_WARPS, 0, stream[0]>>>(*cl, site_beg, n_site, cl->rcut * cl->rcut, counts_d);
  gpuErrchk(hipMemcpyAsync(counts_h, counts_d, (size_t) n_site * sizeof(int), hipMemcpyDeviceToHost, stream[0]));
  gpuErrchk(hipStreamSynchronize(stream[0]));
  gpuErrchk(hipFreeAsync(counts_d, stream[0]));
}

// A batch's list, sites site_beg .. site_beg + n_site - 1 (0-based), in the
// host's layout. n_neigh_h(s) is the site's entries, its own included, as
// gpu_pdf_cells_count gave them plus one; n_pairs their sum. The arrays are
// allocated here on the stream, for gpu_free_neighbors to release.
extern "C" void gpu_pdf_cells_fill(void* h, int site_beg, int n_site, int* n_neigh_h, int n_pairs, void** n_neigh_out,
                                   void** species_out, void** nl_out, void** nspec_out, void** rjs_out, void** xyz_out,
                                   hipStream_t* stream) {
  const pdf_cells_t* cl = (const pdf_cells_t*) h;
  const size_t ns = (size_t) (n_site > 0 ? n_site : 1);
  const size_t np = (size_t) (n_pairs > 0 ? n_pairs : 1);
  int* n_neigh;
  int* species_b;
  int* nl;
  int* nspec;
  double* rjs;
  double* xyz;
  gpuErrchk(hipMallocAsync(&n_neigh, ns * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&species_b, ns * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&nl, np * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&nspec, np * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&rjs, np * sizeof(double), stream[0]));
  gpuErrchk(hipMallocAsync(&xyz, 3 * np * sizeof(double), stream[0]));
  *n_neigh_out = n_neigh;
  *species_out = species_b;
  *nl_out = nl;
  *nspec_out = nspec;
  *rjs_out = rjs;
  *xyz_out = xyz;
  if (n_site <= 0)
    return;

  int* site_end;
  gpuErrchk(hipMallocAsync(&site_end, ns * sizeof(int), stream[0]));
  gpuErrchk(hipMemcpyAsync(n_neigh, n_neigh_h, ns * sizeof(int), hipMemcpyHostToDevice, stream[0]));
  gpuErrchk(hipMemcpyAsync(site_end, n_neigh, ns * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));
  inclusiveScan(site_end, n_site, stream);
  const int blocks = (n_site + PDF_WARPS - 1) / PDF_WARPS;
  kernel_pdf_cells_fill<<<blocks, 32 * PDF_WARPS, 0, stream[0]>>>(*cl, site_beg, n_site, cl->rcut * cl->rcut, n_neigh, site_end,
                                                                  species_b, nl, nspec, rjs, xyz);
  gpuErrchk(hipFreeAsync(site_end, stream[0]));
}

extern "C" void gpu_pdf_cells_free(void* h, hipStream_t* stream) {
  pdf_cells_t* cl = (pdf_cells_t*) h;
  gpuErrchk(hipFreeAsync(cl->x_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->y_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->z_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->atom_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->species_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->species_idx_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->cell_start_d, stream[0]));
  gpuErrchk(hipFreeAsync(cl->cell_of_site_d, stream[0]));
  gpuErrchk(hipStreamSynchronize(stream[0]));
  delete cl;
}
