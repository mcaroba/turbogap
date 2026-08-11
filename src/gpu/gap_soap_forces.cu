// Contraction of the SOAP Cartesian derivatives into forces and the virial,
// and the same contraction for local properties.
#include "gpu_common.h"
#include "gap_gpu.h"
#include "gpu_scatter.h"

#define tpb 64

// The block's own pair, and nothing about where it lands. Both the force and
// the virial used to be added straight from here with atomicAdd, and the order
// the blocks got there is what made the build disagree with itself run to run;
// the sum now happens in gpu_pair_scatter_reduce, in pair order. See
// gpu_scatter.h.
__global__ void cuda_soap_forces_virial_two(int n_sites, double* Qss_d, int n_soap, int* l_index_d, double3* soap_der_d,
                                            double* pair_force_d) {
  int l_nn = blockIdx.x;
  int tid = threadIdx.x;
  int i_site = l_index_d[l_nn] - 1;

  __shared__ double shxthis_block_force[tpb];
  __shared__ double shythis_block_force[tpb];
  __shared__ double shzthis_block_force[tpb];

  shxthis_block_force[tid] = 0;
  shythis_block_force[tid] = 0;
  shzthis_block_force[tid] = 0;

  double locx_this_force = 0;
  double locy_this_force = 0;
  double locz_this_force = 0;

  for (int ii = tid; ii < n_soap; ii = ii + tpb) {
    int i_Qss = i_site + ii * n_sites;      // --> (i, 1:n_soap)
    double loc_this_Qss = Qss_d[i_Qss];     // this read  seems OK
    int in_soap_der = (l_nn * n_soap + ii); // (k,:,l) l is pair index, soap_der(3,n_soap,n_pairs)
    double3 loc_soap_der = soap_der_d[in_soap_der];
    /*     if(isnan( loc_soap_der.x)|| isnan( loc_soap_der.y)||isnan( loc_soap_der.z)){
	   printf("\n loc_soap_der is nan\n");
	   } */
    // if(isnan( loc_this_Qss)){
    //   printf("\n loc_this_Qss is nan  %lf %lf %lf %lf\n", loc_this_Qss, loc_soap_der.x, loc_soap_der.y, loc_soap_der.z);
    // }
    locx_this_force += loc_this_Qss * loc_soap_der.x;
    locy_this_force += loc_this_Qss * loc_soap_der.y;
    locz_this_force += loc_this_Qss * loc_soap_der.z;
  }

  shxthis_block_force[tid] = locx_this_force;
  shythis_block_force[tid] = locy_this_force;
  shzthis_block_force[tid] = locz_this_force;
  /*   if(isnan(locx_this_force)||isnan(locy_this_force)||isnan(locz_this_force)){
       printf("\n loc_this_force is nan\n");
       }   */
  __syncthreads();
  //reduction
  for (int s = tpb / 2; s > 0; s >>= 1) // s=s/2'
  {
    if (tid < s) {
      shxthis_block_force[tid] += shxthis_block_force[tid + s];
      shythis_block_force[tid] += shythis_block_force[tid + s];
      shzthis_block_force[tid] += shzthis_block_force[tid + s];
    }
    __syncthreads();
  }

  //  at this point this_force is computed
  if (tid == 0) {
    /*     if(isnan(shxthis_block_force[0])||isnan(shythis_block_force[0])||isnan(shzthis_block_force[0])){
	   printf("\n this_force is nan\n");
	   } */
    pair_force_d[3 * l_nn] = shxthis_block_force[0];
    pair_force_d[3 * l_nn + 1] = shythis_block_force[0];
    pair_force_d[3 * l_nn + 2] = shzthis_block_force[0];
  }
}

extern "C" void gpu_final_soap_forces_virial(int n_sites, double* Qss_d, int n_soap, int* l_index_d, int* j2_index_d,
                                             double3* soap_der_d, double3* xyz_d, double* virial_d, int n_sites0, double* forces_d,
                                             int n_pairs, hipStream_t* stream) {
  dim3 nblocks(n_pairs, 1);

  hipMemsetAsync(forces_d, 0, 3 * n_sites0 * sizeof(double), stream[0]);
  hipMemsetAsync(virial_d, 0, 9 * sizeof(double), stream[0]);

  // The kernel now stages its result per pair and gpu_pair_scatter_reduce sums
  // it, so that the answer does not depend on the order the blocks finish in.
  double* pair_force_d;
  gpuErrchk(hipMallocAsync(&pair_force_d, (size_t) 3 * n_pairs * sizeof(double), stream[0]));

  cuda_soap_forces_virial_two<<<nblocks, tpb, 0, stream[0]>>>(n_sites, Qss_d, n_soap, l_index_d, soap_der_d, pair_force_d);

  // 0.5 symmetrises; this route does not halve again, its l_nn running over the
  // neighbour list rather than over ordered pairs of a k_list.
  gpu_pair_scatter_reduce(n_pairs, n_sites0, j2_index_d, pair_force_d, (const double*) xyz_d, forces_d, virial_d, 0.5,
                          stream);

  gpuErrchk(hipFreeAsync(pair_force_d, stream[0]));

  return;
}


__global__ void cuda_local_property_derivatives(int n_sites, double* Qss_d, int n_soap, int* l_index_d, double3* soap_der_d,
                                                double* local_property_cart_der_d) {
  int l_nn = blockIdx.x;
  int tid = threadIdx.x;
  int i_site = l_index_d[l_nn] - 1;

  __shared__ double shxthis_block_local_property_cart_der[tpb];
  __shared__ double shythis_block_local_property_cart_der[tpb];
  __shared__ double shzthis_block_local_property_cart_der[tpb];

  shxthis_block_local_property_cart_der[tid] = 0;
  shythis_block_local_property_cart_der[tid] = 0;
  shzthis_block_local_property_cart_der[tid] = 0;

  double locx_this_local_property_cart_der = 0;
  double locy_this_local_property_cart_der = 0;
  double locz_this_local_property_cart_der = 0;

  // Here, every thread acts on a different pair
  // > We sum over the thread id as long as its less than the n_soap
  // > Then we increment by the number of threads in the block
  // > So, i_site is the pair index which comes from blockidx, this is over n_pairs
  // > Each thread goes over n_soap here, such that we can compute
  //   the dot_product( this_Qss(i,1:n_soap), soap_der( cart, 1:n_soap, n_pairs ))
  for (int ii = tid; ii < n_soap; ii = ii + tpb) {
    int i_Qss = i_site + ii * n_sites; // --> (i, 1:n_soap)
    double loc_this_Qss = Qss_d[i_Qss];
    int in_soap_der = (l_nn * n_soap + ii); // (k,:,l) l is pair index, soap_der(3,n_soap,n_pairs)
    double3 loc_soap_der = soap_der_d[in_soap_der];

    locx_this_local_property_cart_der += loc_this_Qss * loc_soap_der.x;
    locy_this_local_property_cart_der += loc_this_Qss * loc_soap_der.y;
    locz_this_local_property_cart_der += loc_this_Qss * loc_soap_der.z;
  }

  shxthis_block_local_property_cart_der[tid] = locx_this_local_property_cart_der;
  shythis_block_local_property_cart_der[tid] = locy_this_local_property_cart_der;
  shzthis_block_local_property_cart_der[tid] = locz_this_local_property_cart_der;

  __syncthreads();

  // Reduction for the dot product
  for (int s = tpb / 2; s > 0; s >>= 1) // s=s/2'
  {
    if (tid < s) {
      shxthis_block_local_property_cart_der[tid] += shxthis_block_local_property_cart_der[tid + s];
      shythis_block_local_property_cart_der[tid] += shythis_block_local_property_cart_der[tid + s];
      shzthis_block_local_property_cart_der[tid] += shzthis_block_local_property_cart_der[tid + s];
    }
    __syncthreads();
  }

  if (tid == 0) {
    // No need for this to be atomic as they should be independent
    local_property_cart_der_d[3 * l_nn] = -shxthis_block_local_property_cart_der[0];
    local_property_cart_der_d[3 * l_nn + 1] = -shythis_block_local_property_cart_der[0];
    local_property_cart_der_d[3 * l_nn + 2] = -shzthis_block_local_property_cart_der[0];
  }
}


extern "C" void gpu_local_property_derivatives(int n_sites, double* Qss_d, int n_soap, int* l_index_d, double3* soap_der_d,
                                               double* local_property_cart_der_d, int n_pairs, hipStream_t* stream) {
  dim3 nblocks(n_pairs, 1);

  /*double *this_force_d;
    hipMalloc((void**)&this_force_d,sizeof(double)*n_pairs*3);*/
  hipMemsetAsync(local_property_cart_der_d, 0, 3 * n_pairs * sizeof(double), stream[0]);

  cuda_local_property_derivatives<<<nblocks, tpb, 0, stream[0]>>>(n_sites, Qss_d, n_soap, l_index_d, soap_der_d,
                                                                  local_property_cart_der_d);

  /*gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );*/

  return;
}


// Local dipoles from a dipole GAP.
//
// mu_i = dE_i/dr_i, the fitted scalar differentiated with respect to the
// central atom's OWN position. soap_turbo builds that self derivative as the
// j == 1 pair of each site, so this is the same contraction
// cuda_local_property_derivatives performs, restricted to one pair per site
// instead of every pair. The grid is therefore n_sites blocks rather than
// n_pairs, which is where nearly all of the saving is: a site has of order 50
// neighbours and only its own entry is wanted.
//
// beg_index_d[i] is the 1-based index of site i's self pair, i.e. the Fortran
// neighbors_beg. The caller supplies Qss_d already carrying +zeta*delta^2, so
// no sign is applied here -- mu is a gradient, not a force.
__global__ void cuda_soap_dipole(int n_sites, double* Qss_d, int n_soap, int* beg_index_d, double3* soap_der_d,
                                 double* dipoles_d) {
  int i_site = blockIdx.x;
  int tid = threadIdx.x;
  int l_self = beg_index_d[i_site] - 1;

  __shared__ double shx[tpb];
  __shared__ double shy[tpb];
  __shared__ double shz[tpb];

  double locx = 0;
  double locy = 0;
  double locz = 0;

  // Each thread strides over the SOAP dimension; the block reduces the three
  // dot products dot(Qss(i,1:n_soap), soap_der(cart,1:n_soap,l_self)).
  for (int ii = tid; ii < n_soap; ii = ii + tpb) {
    int i_Qss = i_site + ii * n_sites;  // Qss is (n_sites, n_soap), column major
    double loc_this_Qss = Qss_d[i_Qss];
    int in_soap_der = (l_self * n_soap + ii);  // soap_der is (3, n_soap, n_pairs)
    double3 loc_soap_der = soap_der_d[in_soap_der];

    locx += loc_this_Qss * loc_soap_der.x;
    locy += loc_this_Qss * loc_soap_der.y;
    locz += loc_this_Qss * loc_soap_der.z;
  }

  shx[tid] = locx;
  shy[tid] = locy;
  shz[tid] = locz;

  __syncthreads();

  for (int s = tpb / 2; s > 0; s >>= 1) {
    if (tid < s) {
      shx[tid] += shx[tid + s];
      shy[tid] += shy[tid + s];
      shz[tid] += shz[tid + s];
    }
    __syncthreads();
  }

  if (tid == 0) {
    // One block per site, so these writes are independent and need no atomics.
    dipoles_d[3 * i_site] = shx[0];
    dipoles_d[3 * i_site + 1] = shy[0];
    dipoles_d[3 * i_site + 2] = shz[0];
  }
}

extern "C" void gpu_soap_dipole(int n_sites, double* Qss_d, int n_soap, int* beg_index_d, double3* soap_der_d,
                                double* dipoles_d, hipStream_t* stream) {
  dim3 nblocks(n_sites, 1);

  hipMemsetAsync(dipoles_d, 0, 3 * n_sites * sizeof(double), stream[0]);

  cuda_soap_dipole<<<nblocks, tpb, 0, stream[0]>>>(n_sites, Qss_d, n_soap, beg_index_d, soap_der_d, dipoles_d);

  return;
}
