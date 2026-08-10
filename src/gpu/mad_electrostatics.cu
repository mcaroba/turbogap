// Electrostatics: the shifted-force Coulomb sum (GSF).
//
// Grouped with the MAD sources rather than the GAP ones because it is reached
// from the same experimental-prediction driver and shares their pair
// compaction, but it is a separate physical model from pdf and xrd.
#include "gpu_common.h"
#include "gpu_scan.h"
#include "mad_gpu.h"

#define tpb 512

__global__ void kernel_get_electrostatics_nk(int i_beg, int i_end, int* n_neigh_d, int* n_neigh_index_d, double* rjs_d, double* xyz,
                                             double r_cut, int* nk_flags_d) {
  int i_site = i_beg - 1 + threadIdx.x + blockIdx.x * blockDim.x;
  int k_val = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j, k, s, k1, k2;
  double r;
  int tid = threadIdx.x;
  int lane = tid % WARP_SIZE;
  int warpId = tid / WARP_SIZE;


  int nk_loc = 0;
  if (i_site < i_end) {
    k = 0;
    for (i = i_beg - 1; i < i_site; i++)
      k += n_neigh_d[i];

    for (j = 0; j < n_neigh_d[i_site]; j++) {
      r = rjs_d[k];
      if (r > r_cut)
        continue;
      nk_loc += 1;
      nk_flags_d[k] = 1;

      k += 1;
    }
    n_neigh_index_d[k_val] = nk_loc;
    //      printf(" - site %d  nk_local %d\n", i_site, nk_loc);
  }
}


extern "C" void gpu_get_electrostatics_nk(int i_beg, int i_end, int n_pairs, int* n_neigh, int* n_neigh_index_d, double* rjs,
                                          double* xyz, double r_cut, int* nk_out_d, int* nk_flags_d, int* nk_flags_sum_d,
                                          hipStream_t* stream) {
  // This function is to set the k_index array for the partial pair distributions


  dim3 nblocks = dim3((i_end - i_beg + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  // gpuErrchk( hipPeekAtLastError() );

  kernel_get_electrostatics_nk<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, n_neigh, n_neigh_index_d, rjs, xyz, r_cut,
                                                                    nk_flags_d);

  recursiveReduce(nk_flags_d, nk_out_d, n_pairs, stream);

  gpuErrchk(hipMemcpyAsync(nk_flags_sum_d, nk_flags_d, n_pairs * sizeof(int), hipMemcpyDeviceToDevice, stream[0]));

  printf("After recursive reduce");
  gpu_peek_stream_error(stream);
  hipError_t err;
  hipDeviceSynchronize();
  err = hipGetLastError();

  inclusiveScan(nk_flags_sum_d, n_pairs, stream);

  printf("After incluseive scan reduce");
  gpu_peek_stream_error(stream);
  hipDeviceSynchronize();
  err = hipGetLastError();


  //  gpu_peek_stream_error( stream );
  // hipError_t err;
  // hipDeviceSynchronize();
  // err = hipGetLastError();

  /* hipError_t code=hipDeviceSynchronize() ; */
  /* printf("\n %s \n", hipGetErrorString(code)); */
  /* gpuErrchk( code ); */

  /* gpuErrchk( hipStreamSynchronize( stream[0] ) ); */
  /* gpuErrchk( hipPeekAtLastError() ); */

  //  inclusiveScan(n_neigh_index_d, i_end - i_beg + 1, stream);

  // Now multiply the flags_sum_d and flags_d to get final nk_array


  gpu_multiply_flags(n_pairs, nk_flags_d, nk_flags_sum_d, stream);

  printf("After multiply flage reduce");
  gpu_peek_stream_error(stream);
  hipDeviceSynchronize();
  err = hipGetLastError();
}

__global__ void kernel_set_electrostatics_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list_d,
                                                  double* rjs, double* xyz, double* charges_d, double* neighbor_charges_index_d,
                                                  double* charge_gradients_d, double* charge_gradients_index_d, int* k_index_d,
                                                  int* j2_index_d, double* rjs_index_d, double* xyz_index_d, int* nk_sum_flags_d) {
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
      j2_index_d[nk] = j2;

      rjs_index_d[nk] = rjs[tid];

      xyz_index_d[3 * nk] = xyz[3 * tid];
      xyz_index_d[3 * nk + 1] = xyz[3 * tid + 1];
      xyz_index_d[3 * nk + 2] = xyz[3 * tid + 2];

      neighbor_charges_index_d[nk] = charges_d[j2];

      charge_gradients_index_d[3 * nk] = charge_gradients_d[3 * tid];
      charge_gradients_index_d[3 * nk + 1] = charge_gradients_d[3 * tid + 1];
      charge_gradients_index_d[3 * nk + 2] = charge_gradients_d[3 * tid + 2];
    }
  }
}


extern "C" void gpu_set_electrostatics_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, double* rjs,
                                               double* xyz, double* charges_d, double* neighbor_charges_index_d,
                                               double* charge_gradients_d, double* charge_gradients_index_d, int* k_index_d,
                                               int* j2_index_d, double* rjs_index_d, double* xyz_k_d, int* nk_sum_flags_d,
                                               hipStream_t* stream) {
  dim3 nblocks = dim3((n_pairs + tpb - 1) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_set_electrostatics_k_index<<<nblocks, nthreads, 0, stream[0]>>>(
      i_beg, i_end, n_pairs, n_sites0, neighbors_list, rjs, xyz, charges_d, neighbor_charges_index_d, charge_gradients_d,
      charge_gradients_index_d, k_index_d, j2_index_d, rjs_index_d, xyz_k_d, nk_sum_flags_d);

  // hipError_t err;
  // hipDeviceSynchronize();
  // err = hipGetLastError();

  /* hipError_t code=hipDeviceSynchronize() ; */
  /* printf("\n %s \n", hipGetErrorString(code)); */
  /* gpuErrchk( code ); */

  /* gpuErrchk( hipStreamSynchronize( stream[0] ) ); */
  /* gpuErrchk( hipPeekAtLastError() ); */
}


__device__ double estat_B0(double rij, double alpha) {
  return erfc(alpha * rij) / rij;
}


__device__ double estat_B0_der_pre(double rij, double alpha, double B0) {
  const double TWO_OVER_SQRT_PI = 1.1283791670955126;
  return -1.0 * B0 / rij - TWO_OVER_SQRT_PI * alpha * exp(-1.0 * alpha * alpha * rij * rij) / rij;
}


__device__ double v_01_gsf(double rij, double alpha, double rcut, double B0, double B0_rcut, double B0_rcut_der) {
  return B0 - B0_rcut - (rij - rcut) * B0_rcut_der;
}
//v_01 = B0 - B0_rcut - (rij - rcut) * B0_rcut_der

__device__ double damping_function_cosine(double distance, double r_inner, double r_outer) {
  double damping = 0.0;
  double arg;
  const double PI = 3.14159265358979323846;

  if (distance < r_inner) {
    damping = 0.0;
  } else if (distance > r_outer) {
    damping = 1.0;
  } else {
    arg = (distance - r_inner) * PI / (r_outer - r_inner);
    damping = 0.5 - 0.5 * cos(arg);
  }
  return damping;
}


__device__ double damping_function_cosine_der(double distance, double r_inner, double r_outer) {
  double damping = 0.0;
  double arg;
  const double PI = 3.14159265358979323846;

  if (distance < r_inner) {
    damping = 0.0;
  } else if (distance > r_outer) {
    damping = 0.0;
  } else {
    arg = (distance - r_inner) * PI / (r_outer - r_inner);
    damping = 0.5 / (r_outer - r_inner) * sin(arg);
  }
  return damping;
}


__global__ void kernel_electrostatics_gsf(const int i_beg, const int nk_max, double* energies_d, double* forces_d, double* virial_d,
                                          int* j2_index_d, const int n_sites, const int this_n_sites, const int this_n_pairs,
                                          int* n_neigh_index_d, double* charges_d, double* charge_gradients_d,
                                          double* neighbor_charges_index_d, double* rjs_index_d, double* xyz_index_d,
                                          const double alpha, const double rcut, const double rcut_in, const double rcut_width,
                                          const double B0_rcut, const double B0_rcut_der, const bool do_cosine_damping,
                                          const bool do_forces) {
  int i_site = i_beg - 1 + blockIdx.x;
  int temp_i_site = blockIdx.x;
  int i, n_strides;
  double rij, neigh_charge, loc_viri;
  int tid = threadIdx.x;
  int pair_idx;
  double this_force[3];
  double this_xyz[3];
  const double K = 14.399645478380902;
  double inner_damp = 1.0;
  double pair_energy, rcut_soft;
  double w_a;
  // In this kernel, each atom is controlled by one block of
  // threads. The threads in the block evaluate each pair and sum them
  // to give the sum of the pairwise energies and the force prefactor.


  // Create the pair offset
  // - The pair offset is given by the inclusive scan of n_neigh_index_d, which is already done beforehand.
  // - Therefore, the offset is given summation before, and the number
  // - of neighbours for the site is given by the subtraction of this
  // - from the current index,
  int n_site_pair_offset = 0;
  if (temp_i_site > 0) {
    n_site_pair_offset = n_neigh_index_d[temp_i_site - 1];
  }

  int n_site_pairs = n_neigh_index_d[temp_i_site] - n_site_pair_offset;

  if (temp_i_site + 1 == this_n_sites) {
    n_site_pairs = n_neigh_index_d[temp_i_site];
  }

  int neigh_idx;

  // Each thread works on a pair, however the number of pairs may
  // exceed the number of threads, hence we set a loop over the number
  // of strides, which is the number of pairs we need to work on
  // divided by the number of threads, with a minimum of 1 stride.
  n_strides = (n_site_pairs - 1 + BLOCK_SIZE) / BLOCK_SIZE;

  // if( tid == 0 ){
  //   printf("- centeridx = %d, cslocal %d,  n_site_pairs = %d, n_strides = %d\n\n", i_site, temp_i_site,  n_site_pairs, n_strides);
  // }

  // The charge of the given site is given by the atom index offset, i_beg, plus
  // the blockidx.x which is the index of the atom
  double charge_site = charges_d[i_site];
  double center_term = K * charge_site;

  __shared__ double shared_energies[BLOCK_SIZE];
  __shared__ double shared_prefactor[BLOCK_SIZE];


  double temp_energy = 0.0;
  double temp_prefactor = 0.0;


  rcut_soft = rcut_in - rcut_width;


  for (int stride_idx = 0; stride_idx < n_strides; ++stride_idx) {
    // The pair_index is offset by the total number of pairs
    // beforehand. Over each stride, it is incremented by the block
    // size.
    pair_idx = n_site_pair_offset + tid + BLOCK_SIZE * stride_idx;

    if (pair_idx - n_site_pair_offset < n_site_pairs && pair_idx - n_site_pair_offset > 0) {
      neigh_charge = neighbor_charges_index_d[pair_idx];

      rij = rjs_index_d[pair_idx];

      double B0 = estat_B0(rij, alpha);

      double f = v_01_gsf(rij, alpha, rcut, B0, B0_rcut, B0_rcut_der);

      // if ( i_site == 3390 - 1 ){
      // 	printf(">eeval site %d, pairid %d / %d tid  %d  rij %lf alpha %lf, B0 %lf f %lf, BOrc, %lf BOrcd %lf\n",
      // 	       i_site,
      // 	       pair_idx, nk_max,
      // 	       tid, rij, alpha,  B0, f, B0_rcut, B0_rcut_der );
      // }


      if (do_cosine_damping) {
        if (rij < rcut_soft) {
          inner_damp = damping_function_cosine(rij, rcut_soft, rcut_in);
        }
      }


      pair_energy = f * center_term * neigh_charge * inner_damp;


      // Adding the energy to the shared data
      temp_energy += 0.5 * pair_energy;

      // Atomic add to the forces_prefactor_d
      temp_prefactor += pair_energy;


      if (do_forces) {
        neigh_idx = j2_index_d[pair_idx];
        double inner_damp_der = 1.0;
        if (do_cosine_damping) {
          if (rij < rcut_soft) {
            inner_damp_der = damping_function_cosine_der(rij, rcut_soft, rcut_in);
          }
        }


        w_a = ((estat_B0_der_pre(rij, alpha, B0) - B0_rcut_der) / rij) * center_term * neigh_charge * inner_damp_der;

        this_xyz[0] = xyz_index_d[3 * pair_idx];
        this_xyz[1] = xyz_index_d[3 * pair_idx + 1];
        this_xyz[2] = xyz_index_d[3 * pair_idx + 2];

        this_force[0] = -this_xyz[0] * w_a;
        this_force[1] = -this_xyz[1] * w_a;
        this_force[2] = -this_xyz[2] * w_a;

        atomicAdd(&forces_d[3 * neigh_idx], this_force[0]);
        atomicAdd(&forces_d[3 * neigh_idx + 1], this_force[1]);
        atomicAdd(&forces_d[3 * neigh_idx + 2], this_force[2]);

        // Virial --------------------------------
#pragma unroll
        for (int k1 = 0; k1 < 3; k1++) {
#pragma unroll
          for (int k2 = 0; k2 < 3; k2++) {
            loc_viri = 0.5 * (this_force[k1] * this_xyz[k2] + this_force[k2] * this_xyz[k1]);
            atomicAdd(&virial_d[k2 + 3 * k1], loc_viri);
          }
        }
        //----------------------------------------
      }
    }
  }


  shared_energies[tid] = temp_energy;
  shared_prefactor[tid] = temp_prefactor;

  // Now we have gone over the strides, we can now reduce over the
  // block in shared memory to get the energy


  // Perform reduction in shared memory
  for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
    __syncthreads();
    if (tid < s) {
      shared_energies[tid] += shared_energies[tid + s];
      shared_prefactor[tid] += shared_prefactor[tid + s];
    }
  }

  // Write the result of the block to global memory
  __syncthreads();
  if (tid == 0) {
    energies_d[blockIdx.x] = shared_energies[0];
  }

  // Now we have the prefactor for each site, we can add the force
  // from the charge gradients

  // Syncing threads so every thread can get the prefactor value
  __syncthreads();
  temp_prefactor = shared_prefactor[0];

  if (do_forces) {
    for (int stride_idx = 0; stride_idx < n_strides; ++stride_idx) {
      // Now we evaluate the variable-charge term
      pair_idx = n_site_pair_offset + tid + BLOCK_SIZE * stride_idx;

      if (pair_idx < n_site_pair_offset + n_site_pairs) {
        // Note we are going over all the gradients now

        neigh_idx = j2_index_d[pair_idx];

        this_force[0] = -temp_prefactor / charge_site * charge_gradients_d[3 * pair_idx];
        this_force[1] = -temp_prefactor / charge_site * charge_gradients_d[3 * pair_idx + 1];
        this_force[2] = -temp_prefactor / charge_site * charge_gradients_d[3 * pair_idx + 2];

        atomicAdd(&forces_d[3 * neigh_idx], this_force[0]);
        atomicAdd(&forces_d[3 * neigh_idx + 1], this_force[1]);
        atomicAdd(&forces_d[3 * neigh_idx + 2], this_force[2]);

        // atomicAdd(&forces_d[neigh_idx].x, this_force[0]);
        // atomicAdd(&forces_d[neigh_idx].y, this_force[1]);
        // atomicAdd(&forces_d[neigh_idx].z, this_force[2]);

        // Virial --------------------------------
#pragma unroll
        for (int k1 = 0; k1 < 3; k1++) {
#pragma unroll
          for (int k2 = 0; k2 < 3; k2++) {
            loc_viri = 0.5 * (this_force[k1] * this_xyz[k2] + this_force[k2] * this_xyz[k1]);
            atomicAdd(&virial_d[k2 + 3 * k1], loc_viri);
          }
        }
        //----------------------------------------
      }
    }
  }
}


extern "C" void gpu_get_electrostatics_energies(const int i_beg, const int nk_max, double* energies_d, double* forces_d,
                                                double* virial_d, int* j2_index_d, const int n_sites, const int this_n_sites,
                                                const int this_n_pairs, int* n_neigh_index_d, double* charges_d,
                                                double* charge_gradients_d, double* neighbor_charges_index_d, double* rjs_index_d,
                                                double* xyz_index_d, const double alpha, const double rcut, const double rcut_in,
                                                const double rcut_width, const double B0_rcut, const double B0_rcut_der,
                                                const bool do_damping, const bool do_forces, hipStream_t* stream) {
  // We want to reduce over the total number of pairs up to the maximum number of pairs
  // Each charge can be associated with a thread
  //
  // So for each site, we need to go up to the number of available pairs
  // This is given by
  // This means we can have the number of blocks as the number of /atoms/
  // Then we should have enough threads to deal with the number of neighbours for that atom.

  // dim3 nblocks=dim3((this_n_sites - 1 + BLOCK_SIZE)/BLOCK_SIZE,1,1);
  // dim3 nthreads=dim3(BLOCK_SIZE,1,1);


  // Can do proper inclusive scan later, for now just sum on cpu
  // // Perform inclusive scan to get the nk indexes
  // inclusiveScan(n_neigh_index_d, this_n_sites, stream);

  // int* n_neigh_index_sum_d;

  // hipHostMalloc((void **) &(n_neigh_index_sum_d), size(int) * this_n_sites);

  // hipMemcpy(  )


  dim3 nblocks = dim3(this_n_sites, 1, 1);
  dim3 nthreads = dim3(BLOCK_SIZE, 1, 1);


  kernel_electrostatics_gsf<<<nblocks, nthreads, 0, stream[0]>>>(
      i_beg, nk_max, energies_d, forces_d, virial_d, j2_index_d, n_sites, this_n_sites, this_n_pairs, n_neigh_index_d, charges_d,
      charge_gradients_d, neighbor_charges_index_d, rjs_index_d, xyz_index_d, alpha, rcut, rcut_in, rcut_width, B0_rcut,
      B0_rcut_der, do_damping, do_forces);
}


//------------------------------------------------------------//
//----------------------   Setup PDF   -----------------------//
//------------------------------------------------------------//
