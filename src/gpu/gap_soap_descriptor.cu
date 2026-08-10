// The SOAP descriptor itself, formed from cnk: the power spectrum, its
// normalisation, the radial/azimuthal/polar derivatives and their conversion
// to Cartesian derivatives, plus the transposes those need.
#include "gpu_common.h"
#include "gap_gpu.h"

#include <cstddef>

#define tpb 64
#define tpb_get_soap_der_one 128
#define tpbcnk 64 // this is because k_max is 45???

//For tiled transpositions
constexpr std::size_t TRANSPOSE_TILE_DIM = 32;
constexpr std::size_t TRANSPOSE_BLOCK_ROWS = 8;

__global__ void cuda_get_soap_p(double* soap_d, double* sqrt_dot_p_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                                bool* skip_soap_component_d, int n_sites, int n_soap, int n_max, int l_max) {
  int i_site = threadIdx.x + blockIdx.x * blockDim.x;
  int k_max = 1 + l_max * (l_max + 1) / 2 + l_max;
  double my_sqrt_dot_p = 0.0;
  if (i_site < n_sites) {
    int counter = 0;
    int counter2 = 0;
    //int ssc_counter=0;
    for (int n = 1; n <= n_max; n++) {
      for (int np = n; np <= n_max; np++) {
        for (int l = 0; l <= l_max; l++) {
          //if(!skip_soap_component_d[ssc_counter]){ //if( skip_soap_component(l, np, n) )cycle
          bool my_skip = skip_soap_component_d[l + (l_max + 1) * (np - 1 + (n - 1) * n_max)];
          if (!(my_skip)) { //if( skip_soap_component(l, np, n) )cycle

            counter++;
            double my_soap = 0.0; //soap_d[counter-1+i_site*n_soap];
            for (int m = 0; m <= l; m++) {
              int k = 1 + l * (l + 1) / 2 + m; //k = 1 + l*(l+1)/2 + m
              counter2++;
              hipDoubleComplex tmp_1_cnk_d =
                  cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)]; //cnk_d[k-1+k_max*(n-1 +i_site*n_max)];
              hipDoubleComplex tmp_2_cnk_d =
                  cnk_d[i_site + n_sites * ((k - 1) + (np - 1) * k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
              my_soap += multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_d.x * tmp_2_cnk_d.x + tmp_1_cnk_d.y * tmp_2_cnk_d.y);
              /*               if(isnan(my_soap)){
			       printf("\n my_soap is nan %lf %lf %lf %lf %lf!!\n", multiplicity_array_d[counter2-1], tmp_1_cnk_d.x, tmp_1_cnk_d.y,tmp_2_cnk_d.x,tmp_2_cnk_d.y);
			       } */
              //soap(counter, i) = soap(counter, i) + multiplicity * real(cnk(k, n, i) * conjg(cnk(k, np, i)))
            }
            soap_d[counter - 1 + i_site * n_soap] = my_soap;
            my_sqrt_dot_p += my_soap * my_soap;
          }
        }
      }
    }
    my_sqrt_dot_p = sqrt(my_sqrt_dot_p);
    if (my_sqrt_dot_p < 1.0e-5) {
      my_sqrt_dot_p = 1.0;
    }
    sqrt_dot_p_d[i_site] = my_sqrt_dot_p;
  }
}

extern "C" void gpu_get_sqrt_dot_p(double* sqrt_dot_d, double* soap_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                                   bool* skip_soap_component_d, int n_sites, int n_soap, int n_max, int l_max,
                                   hipStream_t* stream) {
  dim3 nblocks = dim3((n_sites - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  cuda_get_soap_p<<<nblocks, nthreads, 0, stream[0]>>>(soap_d, sqrt_dot_d, multiplicity_array_d, cnk_d, skip_soap_component_d,
                                                       n_sites, n_soap, n_max, l_max);
  return;
}


__global__ void cuda_get_soap_der_one(double* soap_rad_der_d, double* soap_azi_der_d, double* soap_pol_der_d,
                                      double* multiplicity_array_d, double* trans_soap_rad_der_d, double* trans_soap_azi_der_d,
                                      double* trans_soap_pol_der_d, hipDoubleComplex* cnk_d, hipDoubleComplex* cnk_rad_der_d,
                                      hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, int* k2_i_site_d,
                                      bool* skip_soap_component_d, int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max,
                                      int l_max) {
  int k2 = threadIdx.x + blockIdx.x * blockDim.x;
  if (k2 < n_atom_pairs) {
    int i_site = k2_i_site_d[k2] - 1;
    int counter = 0;
    int counter2 = 0;
    for (int n = 1; n <= n_max; n++) {
      for (int np = n; np <= n_max; np++) {
        for (int l = 0; l <= l_max; l++) {
          if (!skip_soap_component_d
                  [l +
                   (l_max + 1) *
                       (np - 1 +
                        (n - 1) *
                            n_max)]) { //if( skip_soap_component(l, np, n) )cycle // if it happens lots of time, do it in reverse
            counter++;
            double my_soap_rad_der = 0; //trans_soap_rad_der_d[k2+(counter-1)*n_atom_pairs]; //soap_rad_der_d[counter-1+k2*n_soap];
            double my_soap_azi_der = 0; //trans_soap_azi_der_d[k2+(counter-1)*n_atom_pairs]; //soap_azi_der_d[counter-1+k2*n_soap];
            double my_soap_pol_der = 0; //trans_soap_pol_der_d[k2+(counter-1)*n_atom_pairs]; //soap_pol_der_d[counter-1+k2*n_soap];
            for (int m = 0; m <= l; m++) {
              int k = 1 + l * (l + 1) / 2 + m;
              counter2++;
              /* if(threadIdx.x==121 && blockIdx.x==154){
		   printf("\n Pair  %d \n" , k2, i_site);
		   } */
              hipDoubleComplex tmp_1_cnk_d =
                  cnk_d[i_site +
                        n_sites *
                            (k - 1 +
                             (n - 1) *
                                 k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(n-1)*k_max)];  //cnk_d[k-1+ k_max*(n-1 +i_site*n_max)];
              hipDoubleComplex tmp_2_cnk_d =
                  cnk_d[i_site +
                        n_sites *
                            (k - 1 +
                             (np - 1) *
                                 k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(np-1)*k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
              hipDoubleComplex tmp_1_cnk_rad_d = cnk_rad_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (n - 1) *
                            k_max)]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; // cnk_rad_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_rad_d = cnk_rad_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (np - 1) *
                            k_max)]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; // cnk_rad_der_d[k-1+k_max*(np-1+k2*n_max)];
              hipDoubleComplex tmp_1_cnk_azi_d = cnk_azi_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (n - 1) *
                            k_max)]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_azi_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_azi_d = cnk_azi_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (np - 1) *
                            k_max)]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_azi_der_d[k-1+k_max*(np-1+k2*n_max)];
              hipDoubleComplex tmp_1_cnk_pol_d = cnk_pol_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (n - 1) *
                            k_max)]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_pol_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_pol_d = cnk_pol_der_d
                  [k2 +
                   n_atom_pairs *
                       (k - 1 +
                        (np - 1) *
                            k_max)]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_pol_der_d[k-1+k_max*(np-1+k2*n_max)];
              my_soap_rad_der +=
                  multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_rad_d.x * tmp_2_cnk_d.x + tmp_1_cnk_rad_d.y * tmp_2_cnk_d.y +
                                                        tmp_1_cnk_d.x * tmp_2_cnk_rad_d.x + tmp_1_cnk_d.y * tmp_2_cnk_rad_d.y);
              my_soap_azi_der +=
                  multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_azi_d.x * tmp_2_cnk_d.x + tmp_1_cnk_azi_d.y * tmp_2_cnk_d.y +
                                                        tmp_1_cnk_d.x * tmp_2_cnk_azi_d.x + tmp_1_cnk_d.y * tmp_2_cnk_azi_d.y);
              my_soap_pol_der +=
                  multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_pol_d.x * tmp_2_cnk_d.x + tmp_1_cnk_pol_d.y * tmp_2_cnk_d.y +
                                                        tmp_1_cnk_d.x * tmp_2_cnk_pol_d.x + tmp_1_cnk_d.y * tmp_2_cnk_pol_d.y);
            }
            trans_soap_rad_der_d[k2 + (counter - 1) * n_atom_pairs] =
                my_soap_rad_der; //soap_rad_der_d[counter-1+k2*n_soap]=my_soap_rad_der;
            trans_soap_azi_der_d[k2 + (counter - 1) * n_atom_pairs] =
                my_soap_azi_der; //soap_azi_der_d[counter-1+k2*n_soap]=my_soap_azi_der;
            trans_soap_pol_der_d[k2 + (counter - 1) * n_atom_pairs] =
                my_soap_pol_der; //soap_pol_der_d[counter-1+k2*n_soap]=my_soap_pol_der;
          }
        }
      }
    }
  }
}


__global__ void cuda_get_soap_der_two_one(double* soap_d, double* sqrt_dot_p_d, double* soap_rad_der_d, double* soap_azi_der_d,
                                          double* soap_pol_der_d, double* trans_soap_rad_der_d, double* trans_soap_azi_der_d,
                                          double* trans_soap_pol_der_d, double* tdotoprod_der_rad, double* tdotoprod_der_azi,
                                          double* tdotoprod_der_pol, int* k2_i_site_d, int n_sites, int n_atom_pairs, int n_soap,
                                          int k_max, int n_max, int l_max) {
  int k2 = blockIdx.x;
  int tid = threadIdx.x;
  int i_site = k2_i_site_d[k2] - 1;
  __shared__ double sh_soap_rad_der_dot[tpb];
  __shared__ double sh_soap_azi_der_dot[tpb];
  __shared__ double sh_soap_pol_der_dot[tpb];
  double this_dotprod_rad = 0.0;
  double this_dotprod_azi = 0.0;
  double this_dotprod_pol = 0.0;

  for (int s = tid; s < n_soap; s = s + tpb) {
    this_dotprod_rad += soap_d[s + i_site * n_soap] * soap_rad_der_d[s + k2 * n_soap];
    this_dotprod_azi += soap_d[s + i_site * n_soap] * soap_azi_der_d[s + k2 * n_soap];
    this_dotprod_pol += soap_d[s + i_site * n_soap] * soap_pol_der_d[s + k2 * n_soap];
  }
  sh_soap_rad_der_dot[tid] = this_dotprod_rad;
  sh_soap_azi_der_dot[tid] = this_dotprod_azi;
  sh_soap_pol_der_dot[tid] = this_dotprod_pol;
  __syncthreads();

  //reduction
  for (int s = tpb / 2; s > 0; s >>= 1) // s=s/2
  {
    if (tid < s) {
      sh_soap_rad_der_dot[tid] += sh_soap_rad_der_dot[tid + s];
      sh_soap_azi_der_dot[tid] += sh_soap_azi_der_dot[tid + s];
      sh_soap_pol_der_dot[tid] += sh_soap_pol_der_dot[tid + s];
    }
    __syncthreads();
  }
  for (int s = tid; s < n_soap; s = s + tpb) {
    tdotoprod_der_rad[s + k2 * n_soap] = sh_soap_rad_der_dot[0];
    tdotoprod_der_azi[s + k2 * n_soap] = sh_soap_azi_der_dot[0];
    tdotoprod_der_pol[s + k2 * n_soap] = sh_soap_pol_der_dot[0];
  }
}

__global__ void cuda_get_soap_der_two_two(double* soap_d, double* sqrt_dot_p_d, double* soap_rad_der_d, double* soap_azi_der_d,
                                          double* soap_pol_der_d, double* tdotoprod_der_rad, double* tdotoprod_der_azi,
                                          double* tdotoprod_der_pol, int* k2_i_site_d, int n_sites, int n_atom_pairs, int n_soap,
                                          int k_max, int n_max, int l_max) {
  int k2 = blockIdx.x;
  int tid = threadIdx.x;
  int i_site = k2_i_site_d[k2] - 1;
  double loc_sqrt_dot_p = sqrt_dot_p_d[i_site];
  for (int s = tid; s < n_soap; s = s + tpb) {
    double my_soap = soap_d[s + i_site * n_soap];

    double my_soap_rad_der = soap_rad_der_d[s + k2 * n_soap];
    double my_soap_azi_der = soap_azi_der_d[s + k2 * n_soap];
    double my_soap_pol_der = soap_pol_der_d[s + k2 * n_soap];

    double myprod_der_rad = tdotoprod_der_rad[s + k2 * n_soap];
    double myprod_der_azi = tdotoprod_der_azi[s + k2 * n_soap];
    double myprod_der_pol = tdotoprod_der_pol[s + k2 * n_soap];


    soap_rad_der_d[s + k2 * n_soap] =
        my_soap_rad_der / loc_sqrt_dot_p - my_soap / (loc_sqrt_dot_p * loc_sqrt_dot_p * loc_sqrt_dot_p) * myprod_der_rad;
    soap_azi_der_d[s + k2 * n_soap] =
        my_soap_azi_der / loc_sqrt_dot_p - my_soap / (loc_sqrt_dot_p * loc_sqrt_dot_p * loc_sqrt_dot_p) * myprod_der_azi;
    soap_pol_der_d[s + k2 * n_soap] =
        my_soap_pol_der / loc_sqrt_dot_p - my_soap / (loc_sqrt_dot_p * loc_sqrt_dot_p * loc_sqrt_dot_p) * myprod_der_pol;
  }
}


__global__ void cuda_get_soap_der_thr_one(double3* soap_cart_der_d, double* soap_rad_der_d, double* soap_azi_der_d,
                                          double* soap_pol_der_d, double* thetas, double* phis, double* rjs, int* k3_index,
                                          int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max) {
  int k2 = blockIdx.x;
  int tid = threadIdx.x;
  int k3 = k3_index[k2] - 1;

  double my_theta = thetas[k2];
  double my_phi = phis[k2];
  double my_rj = rjs[k2];
  for (int s = tid; s < n_soap; s = s + tpb) {
    if (k3 != k2) {
      double my_soap_rad_der = soap_rad_der_d[s + k2 * n_soap];
      double my_soap_azi_der = soap_azi_der_d[s + k2 * n_soap];
      double my_soap_pol_der = soap_pol_der_d[s + k2 * n_soap];
      double3 my_soap_cart_der;
      my_soap_cart_der.x = sin(my_theta) * cos(my_phi) * my_soap_rad_der - cos(my_theta) * cos(my_phi) / my_rj * my_soap_pol_der -
                           sin(my_phi) / my_rj * my_soap_azi_der;
      my_soap_cart_der.y = sin(my_theta) * sin(my_phi) * my_soap_rad_der - cos(my_theta) * sin(my_phi) / my_rj * my_soap_pol_der +
                           cos(my_phi) / my_rj * my_soap_azi_der;
      my_soap_cart_der.z = cos(my_theta) * my_soap_rad_der + sin(my_theta) / my_rj * my_soap_pol_der;
      soap_cart_der_d[s + k2 * n_soap] = my_soap_cart_der;
    }
  }
}

__global__ void cuda_get_soap_der_thr_two(double3* soap_cart_der_d, double* soap_rad_der_d, double* soap_azi_der_d,
                                          double* soap_pol_der_d, double* thetas, double* phis, double* rjs, int* n_neigh_d,
                                          int* i_k2_start_d, int* k2_i_site_d, int* k3_index_d, int n_sites, int n_atom_pairs,
                                          int n_soap, int k_max, int n_max, int l_max, int maxneigh) {
  int i_site = blockIdx.x;
  int tid = threadIdx.x;
  int my_start = i_k2_start_d[i_site] - 1;
  int k3 = my_start;
  int my_n_neigh = n_neigh_d[i_site];

  for (int s = tid; s < n_soap; s = s + tpb) {
    double3 loc_sum;
    loc_sum.x = 0, loc_sum.y = 0, loc_sum.z = 0;
    int k2 = my_start + 1;
    for (int j = 1; j < my_n_neigh; j++) {
      double3 my_soap_cart_der = soap_cart_der_d[s + k2 * n_soap];
      loc_sum.x -= my_soap_cart_der.x;
      loc_sum.y -= my_soap_cart_der.y;
      loc_sum.z -= my_soap_cart_der.z;
      k2++;
    }
    soap_cart_der_d[s + k3 * n_soap] = loc_sum;
  }
}


__global__ void naive_transpose_soap_rad_azi_pol(double* soap_rad_der_d, double* tran_soap_rad_der_d, int n_soap,
                                                 int n_atom_pairs) {
  // Handles bank conflicts too
  __shared__ double tile[TRANSPOSE_TILE_DIM][TRANSPOSE_TILE_DIM + 1];
  int x = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.y;

  for (int j = 0; j < TRANSPOSE_TILE_DIM; j += TRANSPOSE_BLOCK_ROWS) {
    if (x < n_soap && (y + j) < n_atom_pairs) {
      tile[threadIdx.y + j][threadIdx.x] = soap_rad_der_d[(y + j) * n_soap + x];
    }
  }

  __syncthreads();

  x = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.x;
  y = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.y;
  for (int j = 0; j < TRANSPOSE_TILE_DIM; j += TRANSPOSE_BLOCK_ROWS) {
    if (x < n_atom_pairs && (y + j) < n_soap) {
      tran_soap_rad_der_d[(y + j) * n_atom_pairs + x] = tile[threadIdx.x][threadIdx.y + j];
    }
  }
}

__global__ void naive_transpose_cnk_arrays(hipDoubleComplex* C, hipDoubleComplex* tran_C, int k_max, int n_max, int n_sites) {
  // in Fortran is cnk( 1:k_max, 1:n_max, 1:n_sites) --> (1:n_sites,1:k_max, 1:n_max)
  //       cnk_rad_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
  int i_g = threadIdx.x + blockIdx.x * blockDim.x;
  if (i_g < k_max * n_max * n_sites) {
    hipDoubleComplex loc_C = C[i_g]; // i_g=i_k+k_max*(i_n+i_site*n_max)
    int i_k = i_g % k_max;
    int i_z = i_g / k_max;
    int i_n = i_z % n_max;
    int i_site = i_z / n_max;
    int new_i_g = i_site + n_sites * (i_k + i_n * k_max);
    tran_C[new_i_g] = loc_C;
  }
}

extern "C" void gpu_get_soap_der(double* soap_d, double* sqrt_dot_d, double3* soap_cart_der_d, double* soap_rad_der_d,
                                 double* soap_azi_der_d, double* soap_pol_der_d, double* thetas_d, double* phis_d, double* rjs_d,
                                 double* multiplicity_array_d, hipDoubleComplex* cnk_d, hipDoubleComplex* cnk_rad_der_d,
                                 hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, int* n_neigh_d,
                                 int* i_k2_start_d, int* k2_i_site_d, int* k3_index_d, bool* skip_soap_component_d, int n_sites,
                                 int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max, int maxneigh, hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  dim3 nblocks_get_soap_der_one = dim3((n_atom_pairs - 1 + tpb_get_soap_der_one) / tpb_get_soap_der_one, 1, 1);
  dim3 nthreads_get_soap_der_one = dim3(tpb_get_soap_der_one, 1, 1);
  //size_t mf, ma;
  //hipMemGetInfo(&mf, &ma);
  //printf("\n free: %zu total: %zu", mf, ma);
  double *tdotoprod_der_rad, *tdotoprod_der_azi, *tdotoprod_der_pol;
  hipMallocAsync((void**) &tdotoprod_der_rad, sizeof(double) * n_atom_pairs * n_soap, stream[0]);
  hipMallocAsync((void**) &tdotoprod_der_azi, sizeof(double) * n_atom_pairs * n_soap, stream[0]);
  hipMallocAsync((void**) &tdotoprod_der_pol, sizeof(double) * n_atom_pairs * n_soap, stream[0]);


  double *trans_soap_rad_der_d, *trans_soap_azi_der_d, *trans_soap_pol_der_d;
  hipMallocAsync((void**) &trans_soap_rad_der_d, sizeof(double) * n_atom_pairs * n_soap, stream[0]);
  hipMallocAsync((void**) &trans_soap_azi_der_d, sizeof(double) * n_atom_pairs * n_soap, stream[0]);
  hipMallocAsync((void**) &trans_soap_pol_der_d, sizeof(double) * n_atom_pairs * n_soap, stream[0]);


  cuda_get_soap_der_one<<<nblocks_get_soap_der_one, nthreads_get_soap_der_one, 0, stream[0]>>>(
      soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, multiplicity_array_d, trans_soap_rad_der_d, trans_soap_azi_der_d,
      trans_soap_pol_der_d, cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, k2_i_site_d, skip_soap_component_d, n_sites,
      n_atom_pairs, n_soap, k_max, n_max, l_max);

  dim3 transpose_block(TRANSPOSE_TILE_DIM, TRANSPOSE_BLOCK_ROWS, 1);
  dim3 transpose_grid((n_atom_pairs + TRANSPOSE_TILE_DIM - 1) / TRANSPOSE_TILE_DIM,
                      (n_soap + TRANSPOSE_TILE_DIM - 1) / TRANSPOSE_TILE_DIM, 1);

  naive_transpose_soap_rad_azi_pol<<<transpose_grid, transpose_block, 0, stream[0]>>>(trans_soap_rad_der_d, soap_rad_der_d,
                                                                                      n_atom_pairs, n_soap);

  naive_transpose_soap_rad_azi_pol<<<transpose_grid, transpose_block, 0, stream[0]>>>(trans_soap_azi_der_d, soap_azi_der_d,
                                                                                      n_atom_pairs, n_soap);

  naive_transpose_soap_rad_azi_pol<<<transpose_grid, transpose_block, 0, stream[0]>>>(trans_soap_pol_der_d, soap_pol_der_d,
                                                                                      n_atom_pairs, n_soap);

  cuda_get_soap_der_two_one<<<n_atom_pairs, nthreads, 0, stream[0]>>>(
      soap_d, sqrt_dot_d, soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, trans_soap_rad_der_d, trans_soap_azi_der_d,
      trans_soap_pol_der_d, tdotoprod_der_rad, tdotoprod_der_azi, tdotoprod_der_pol, k2_i_site_d, n_sites, n_atom_pairs, n_soap,
      k_max, n_max, l_max);


  cuda_get_soap_der_two_two<<<n_atom_pairs, nthreads, 0, stream[0]>>>(
      soap_d, sqrt_dot_d, soap_rad_der_d, soap_azi_der_d, soap_pol_der_d,
      //trans_soap_rad_der_d, trans_soap_azi_der_d, trans_soap_pol_der_d,
      tdotoprod_der_rad, tdotoprod_der_azi, tdotoprod_der_pol, k2_i_site_d, n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max);

  cuda_get_soap_der_thr_one<<<n_atom_pairs, nthreads, 0, stream[0]>>>(soap_cart_der_d, soap_rad_der_d, soap_azi_der_d,
                                                                      soap_pol_der_d, thetas_d, phis_d, rjs_d, k3_index_d, n_sites,
                                                                      n_atom_pairs, n_soap, k_max, n_max, l_max);
  cuda_get_soap_der_thr_two<<<n_sites, nthreads, 0, stream[0]>>>(
      soap_cart_der_d, soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, thetas_d, phis_d, rjs_d, n_neigh_d, i_k2_start_d,
      k2_i_site_d, k3_index_d, n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh);
  hipFreeAsync(tdotoprod_der_rad, stream[0]);
  hipFreeAsync(tdotoprod_der_azi, stream[0]);
  hipFreeAsync(tdotoprod_der_pol, stream[0]);
  hipFreeAsync(trans_soap_rad_der_d, stream[0]);
  hipFreeAsync(trans_soap_azi_der_d, stream[0]);
  hipFreeAsync(trans_soap_pol_der_d, stream[0]);
  /* hipFreeAsync(trans_cnk_d,0); */
  /* hipFreeAsync(trans_cnk_rad_der_d,0);hipFreeAsync(trans_cnk_azi_der_d,0);hipFreeAsync(trans_cnk_pol_der_d,0); */

  // hipError_t code=hipDeviceSynchronize() ;
  // printf("\n %s \n", hipGetErrorString(code));
  // gpuErrchk( code );
  return;
}

__global__ void cuda_soap_normalize(double* soap_d, double* sqrt_dot_d, int n_soap, int n_sites) {
  int i_site = blockIdx.x;
  int tid = threadIdx.x;
  double my_sqrt_dot_p = sqrt_dot_d[i_site];
  for (int s = tid; s < n_soap; s = s + tpb) {
    double my_soap_final = soap_d[s + i_site * n_soap] / my_sqrt_dot_p;
    soap_d[s + i_site * n_soap] = my_soap_final;
  }
}
extern "C" void gpu_soap_normalize(double* soap_d, double* sqrt_dot_d, int n_soap, int n_sites, hipStream_t* stream) {
  cuda_soap_normalize<<<n_sites, tpb, 0, stream[0]>>>(soap_d, sqrt_dot_d, n_soap, n_sites);
}
// Two earlier versions of the derivative kernel (cuda_get_derivatives and
// cuda_get_derivatives_new) used to live here, inside a block comment and
// with their only call site commented out as well. They also wrote
// cnk_*_der_d with a different memory layout than the live kernel below and
// its consumer use, so they were a trap for anyone editing this path.
// Removed.
__global__ void cuda_get_derivatives_new_new(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d,
                                             double* radial_exp_coeff_der_d, hipDoubleComplex* angular_exp_coeff_rad_der_d,
                                             hipDoubleComplex* angular_exp_coeff_azi_der_d,
                                             hipDoubleComplex* angular_exp_coeff_pol_der_d, hipDoubleComplex* cnk_rad_der_d,
                                             hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, double* rjs_d,
                                             double rcut_max, int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max,
                                             int l_max) {
  int k2 = threadIdx.x + blockDim.x * blockIdx.x;
  int n = blockIdx.y + 1;
  int k = blockIdx.z + 1;
  double Pi = 4.0 * acos(-1.0);
  if (k2 < n_atom_pairs) {
    double my_rjs = rjs_d[k2];
    if (my_rjs < rcut_max) {
      double my_radial_exp_c;
      my_radial_exp_c = radial_exp_coeff_d[n - 1 + k2 * n_max];
      double my_radial_exp_c_der;
      my_radial_exp_c_der = radial_exp_coeff_der_d[n - 1 + k2 * n_max];
      //int k=1+l*(l+1)/2+m;
      hipDoubleComplex my_cnk_rad_der;
      hipDoubleComplex my_cnk_azi_der;
      hipDoubleComplex my_cnk_pol_der;
      hipDoubleComplex my_ang_exp_c;
      my_ang_exp_c = angular_exp_coeff_d[k2 + n_atom_pairs * (k - 1)]; //angular_exp_coeff_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_rad_der;
      my_ang_exp_c_rad_der = angular_exp_coeff_rad_der_d[k2 + (k - 1) * n_atom_pairs]; //angular_exp_coeff_rad_der_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_azi_der;
      my_ang_exp_c_azi_der = angular_exp_coeff_azi_der_d[k2 + (k - 1) * n_atom_pairs]; //angular_exp_coeff_azi_der_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_pol_der;
      my_ang_exp_c_pol_der = angular_exp_coeff_pol_der_d[k2 + (k - 1) * n_atom_pairs]; //angular_exp_coeff_pol_der_d[k-1+k2*k_max];

      my_cnk_rad_der.x = Pi * (my_ang_exp_c.x * my_radial_exp_c_der + my_ang_exp_c_rad_der.x * my_radial_exp_c);
      my_cnk_rad_der.y = Pi * (my_ang_exp_c.y * my_radial_exp_c_der + my_ang_exp_c_rad_der.y * my_radial_exp_c);

      my_cnk_azi_der.x = Pi * (my_ang_exp_c_azi_der.x * my_radial_exp_c);
      my_cnk_azi_der.y = Pi * (my_ang_exp_c_azi_der.y * my_radial_exp_c);

      my_cnk_pol_der.x = Pi * (my_ang_exp_c_pol_der.x * my_radial_exp_c);
      my_cnk_pol_der.y = Pi * (my_ang_exp_c_pol_der.y * my_radial_exp_c);

      //i_site+n_sites*(i_k+i_n*k_max);
      cnk_rad_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] =
          my_cnk_rad_der; //cnk_rad_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_rad_der;
      cnk_azi_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] =
          my_cnk_azi_der; //cnk_azi_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_azi_der;
      cnk_pol_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] =
          my_cnk_pol_der; //cnk_pol_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_pol_der;
    }
  }
}


extern "C" void gpu_get_derivatives(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d,
                                    double* radial_exp_coeff_der_d, hipDoubleComplex* angular_exp_coeff_rad_der_d,
                                    hipDoubleComplex* angular_exp_coeff_azi_der_d, hipDoubleComplex* angular_exp_coeff_pol_der_d,
                                    hipDoubleComplex* cnk_rad_der_d, hipDoubleComplex* cnk_azi_der_d,
                                    hipDoubleComplex* cnk_pol_der_d, double* rjs_d, double rcut_max, int n_atom_pairs, int n_sites,
                                    int n_soap, int k_max, int n_max, int l_max, hipStream_t* stream) {
  cuda_get_derivatives_new_new<<<dim3((n_atom_pairs + tpbcnk - 1) / tpbcnk, n_max, k_max), tpbcnk, 0, stream[0]>>>(
      radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d, angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d,
      angular_exp_coeff_pol_der_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, rjs_d, rcut_max, n_atom_pairs, n_sites, n_soap,
      k_max, n_max, l_max);

  /*}
      hipEventRecord(stop);
      hipEventSynchronize(stop);
      milliseconds = 0.0;
      hipEventElapsedTime(&milliseconds, start, stop);
      printf("\n Time of the second kernel in s %f\n", milliseconds/1000.0);

      exit(0);*/
}
