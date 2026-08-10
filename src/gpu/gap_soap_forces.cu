// Contraction of the SOAP Cartesian derivatives into forces and the virial,
// and the same contraction for local properties.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64

__global__ void cuda_soap_forces_virial_two(int n_sites, double* Qss_d, int n_soap, int* l_index_d, int* j2_index_d,
                                            double3* soap_der_d, double3* xyz_d, double* virial_d, int n_sites0, double* forces_d) {
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
    int j2 = j2_index_d[l_nn] - 1;
    atomicAdd(&forces_d[j2 * 3], shxthis_block_force[0]);
    atomicAdd(&forces_d[j2 * 3 + 1], shythis_block_force[0]);
    atomicAdd(&forces_d[j2 * 3 + 2], shzthis_block_force[0]);

    // now the virial
    double this_force[3];
    this_force[0] = shxthis_block_force[0];
    this_force[1] = shythis_block_force[0];
    this_force[2] = shzthis_block_force[0];
    /*     if(isnan(shxthis_block_force[0])||isnan(shythis_block_force[0])||isnan(shzthis_block_force[0])){
	   printf("\n this_force is nan\n");
	   } */
    /* if(isnan(this_force[0])||isnan(this_force[1])||isnan(this_force[2])){
       printf("\n this_force is nan\n");
       } */

    double3 tmp_xyz;
    tmp_xyz = xyz_d[l_nn];
    double this_xyz[3];
    this_xyz[0] = tmp_xyz.x;
    this_xyz[1] = tmp_xyz.y;
    this_xyz[2] = tmp_xyz.z;

    for (int k1 = 0; k1 < 3; k1++) {
      for (int k2 = 0; k2 < 3; k2++) {
        double loc_viri = 0.5 * (this_force[k1] * this_xyz[k2] + this_force[k2] * this_xyz[k1]);
        /*         if(isnan(loc_viri)){
		   printf("\n locviri is nan\n");
		   } */
        atomicAdd(&virial_d[k2 + 3 * k1], loc_viri);
      }
    }
    /*     if(isnan(tmp_xyz.x)||isnan(tmp_xyz.y)||isnan(tmp_xyz.z)){
	   printf("\n tmp is nan\n");
	   } */
  }
}

extern "C" void gpu_final_soap_forces_virial(int n_sites, double* Qss_d, int n_soap, int* l_index_d, int* j2_index_d,
                                             double3* soap_der_d, double3* xyz_d, double* virial_d, int n_sites0, double* forces_d,
                                             int n_pairs, hipStream_t* stream) {
  dim3 nblocks(n_pairs, 1);

  /*double *this_force_d;
    hipMalloc((void**)&this_force_d,sizeof(double)*n_pairs*3);*/
  hipMemsetAsync(forces_d, 0, 3 * n_sites0 * sizeof(double), stream[0]);
  hipMemsetAsync(virial_d, 0, 9 * sizeof(double), stream[0]);

  cuda_soap_forces_virial_two<<<nblocks, tpb, 0, stream[0]>>>(n_sites, Qss_d, n_soap, l_index_d, j2_index_d, soap_der_d, xyz_d,
                                                              virial_d, n_sites0, forces_d);

  /*gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );*/

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
