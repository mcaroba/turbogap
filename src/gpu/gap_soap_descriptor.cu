// The SOAP descriptor itself, formed from cnk: the power spectrum, its
// normalisation, the radial/azimuthal/polar derivatives and their conversion
// to Cartesian derivatives, plus the transposes those need.
// Partly converted. cuda_get_soap_p, _der_one and _index are flat maps and
// go through tg_parallel_for; the rest do not, and are left as <<<>>> launches
// on purpose:
//
//   cuda_get_soap_der_two_one
//       closes with a __shared__ TREE REDUCTION OVER DOUBLES. Kokkos would sum
//       the block in a different order and the result would move, so this one
//       is not merely unconverted -- it must not be converted without a
//       TeamPolicy reproducing this exact tree order. See gap_soap_forces.cu.
//   naive_transpose_soap_rad_azi_pol
//       __shared__ used as a transpose tile: shared data, one barrier, no
//       accumulation. TeamPolicy with team_scratch, no numerical hazard.
//   cuda_get_soap_der_two_two, _thr_one, _thr_two, cuda_soap_normalize,
//   cuda_get_derivatives_new_new
//       the block index is part of the addressing rather than a tiling of one
//       range. MDRangePolicy, or an explicit flattening. The last now runs
//       only for the TURBOGAP_DUMP_CNK dump: _der_one forms the same values.
//
// See docs/LAMMPS_KOKKOS.md.
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gap_gpu.h"
#include "gpu_memory.h"

#include <cstddef>

#define tpb 64
#define tpb_get_soap_der_one 128
#define tpbcnk 64 // this is because k_max is 45???

//For tiled transpositions
constexpr std::size_t TRANSPOSE_TILE_DIM = 32;
constexpr std::size_t TRANSPOSE_BLOCK_ROWS = 8;

// Where (n, np, l) lands: .x the output column (-1 if skipped), .y the first
// multiplicity index of its m run. The order is the one the loops used to walk;
// each entry adds the kept ones of the rows n before it, counted once per row,
// to those before it in its own row.
static void cuda_get_soap_index(int2* index_d, bool* skip_soap_component_d, int n_max, int l_max, hipStream_t* stream) {
  int2* row_d;
  hipMallocAsync((void**) &row_d, sizeof(int2) * n_max, stream[0]);
  tg_parallel_for(
      "turbogap_soap_soap_index_rows", n_max, stream,
      TG_LAMBDA(const int r) {
        int n = r + 1;
        int counter = 0;
        int counter2 = 0;
        for (int np = n; np <= n_max; np++) {
          for (int l = 0; l <= l_max; l++) {
            if (!skip_soap_component_d[l + (l_max + 1) * (np - 1 + (n - 1) * n_max)]) {
              counter++;
              counter2 += l + 1;
            }
          }
        }
        row_d[r] = make_int2(counter, counter2);
      },
      tpb);
  tg_parallel_for(
      "turbogap_soap_soap_index", n_max * n_max * (l_max + 1), stream,
      TG_LAMBDA(const int s) {
        int l = s % (l_max + 1);
        int np = (s / (l_max + 1)) % n_max + 1;
        int n = s / (l_max + 1) / n_max + 1;
        if (np < n)
          return;
        if (skip_soap_component_d[s]) {
          index_d[s] = make_int2(-1, 0);
          return;
        }
        int counter = 0;
        int counter2 = 0;
        for (int r = 0; r < n - 1; r++) {
          counter += row_d[r].x;
          counter2 += row_d[r].y;
        }
        for (int np1 = n; np1 <= np; np1++) {
          for (int l1 = 0; l1 <= l_max; l1++) {
            if (np1 == np && l1 == l) {
              index_d[s] = make_int2(counter, counter2);
              return;
            }
            if (!skip_soap_component_d[l1 + (l_max + 1) * (np1 - 1 + (n - 1) * n_max)]) {
              counter++;
              counter2 += l1 + 1;
            }
          }
        }
      },
      tpb);
  hipFreeAsync(row_d, stream[0]);
}

// One thread per (site, n, l), then one per site for the norm, summed over the
// components in output order as before; one thread per site left most of the
// device idle.
static void cuda_get_soap_p(double* soap_d, double* sqrt_dot_p_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                            int2* index_d, int n_sites, int n_soap, int n_max, int l_max, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_soap_soap_p", n_sites * n_max * (l_max + 1), stream,
      TG_LAMBDA(const int t) {
        int i_site = t % n_sites;
        int l = (t / n_sites) % (l_max + 1);
        int n = t / n_sites / (l_max + 1) + 1;
        int k_max = 1 + l_max * (l_max + 1) / 2 + l_max;
        for (int np = n; np <= n_max; np++) {
          int2 idx = index_d[l + (l_max + 1) * (np - 1 + (n - 1) * n_max)];
          if (idx.x < 0)
            continue;
          double my_soap = 0.0;
          for (int m = 0; m <= l; m++) {
            int k = 1 + l * (l + 1) / 2 + m;
            hipDoubleComplex tmp_1_cnk_d = cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)];
            hipDoubleComplex tmp_2_cnk_d = cnk_d[i_site + n_sites * ((k - 1) + (np - 1) * k_max)];
            my_soap += multiplicity_array_d[idx.y + m] * (tmp_1_cnk_d.x * tmp_2_cnk_d.x + tmp_1_cnk_d.y * tmp_2_cnk_d.y);
          }
          soap_d[idx.x + i_site * n_soap] = my_soap;
        }
      },
      tpb);
  tg_parallel_for(
      "turbogap_soap_sqrt_dot_p", n_sites, stream,
      TG_LAMBDA(const int i_site) {
        double my_sqrt_dot_p = 0.0;
        for (int c = 0; c < n_soap; c++) {
          double my_soap = soap_d[c + i_site * n_soap];
          my_sqrt_dot_p += my_soap * my_soap;
        }
        my_sqrt_dot_p = sqrt(my_sqrt_dot_p);
        if (my_sqrt_dot_p < 1.0e-5) {
          my_sqrt_dot_p = 1.0;
        }
        sqrt_dot_p_d[i_site] = my_sqrt_dot_p;
      },
      tpb);
}

extern "C" void gpu_get_sqrt_dot_p(double* sqrt_dot_d, double* soap_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                                   bool* skip_soap_component_d, int n_sites, int n_soap, int n_max, int l_max,
                                   hipStream_t* stream) {
  int2* index_d;
  hipMallocAsync((void**) &index_d, sizeof(int2) * n_max * n_max * (l_max + 1), stream[0]);
  cuda_get_soap_index(index_d, skip_soap_component_d, n_max, l_max, stream);
  cuda_get_soap_p(soap_d, sqrt_dot_d, multiplicity_array_d, cnk_d, index_d, n_sites, n_soap, n_max, l_max, stream);
  hipFreeAsync(index_d, stream[0]);
}


// cnk_rad_der and cnk_azi_der/cnk_pol_der of one pair, n and k, from the
// radial and angular expansion coefficients they are the product of. Spelled
// as the operations the stand-alone kernel compiled to, so that inlined into
// another kernel they cannot be contracted differently. The rounding
// intrinsics are device-only; Kokkos also compiles these for the host.
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define CNK_MUL(a, b) __dmul_rn(a, b)
#define CNK_FMA(a, b, c) __fma_rn(a, b, c)
#else
#define CNK_MUL(a, b) ((a) * (b))
#define CNK_FMA(a, b, c) fma(a, b, c)
#endif

TG_INLINE_FUNCTION hipDoubleComplex cnk_rad_der_of(hipDoubleComplex ang, hipDoubleComplex ang_rad_der, double rad, double rad_der,
                                                   double Pi) {
  hipDoubleComplex c;
  c.x = CNK_MUL(Pi, CNK_FMA(ang.x, rad_der, CNK_MUL(ang_rad_der.x, rad)));
  c.y = CNK_MUL(Pi, CNK_FMA(ang.y, rad_der, CNK_MUL(ang_rad_der.y, rad)));
  return c;
}

TG_INLINE_FUNCTION hipDoubleComplex cnk_ang_der_of(hipDoubleComplex ang_der, double rad, double Pi) {
  hipDoubleComplex c;
  c.x = CNK_MUL(Pi, CNK_MUL(ang_der.x, rad));
  c.y = CNK_MUL(Pi, CNK_MUL(ang_der.y, rad));
  return c;
}

// One thread per (pair, n, l) rather than per pair: one per pair left about ten
// warps on an SM and the kernel waited on its loads. Each output sums over m in
// the same order and with the same expression as before, so it is
// bit-identical; the n-side coefficients are loaded once per block of
// DER_ONE_NP np values instead of once per np. cnk_*_der are formed here from
// the expansion coefficients rather than read from arrays three times the
// size of the descriptor derivatives.
#define DER_ONE_NP 4

static void cuda_get_soap_der_one(double* multiplicity_array_d, double* trans_soap_rad_der_d, double* trans_soap_azi_der_d,
                                  double* trans_soap_pol_der_d, hipDoubleComplex* cnk_d, double* radial_exp_coeff_d,
                                  double* radial_exp_coeff_der_d, hipDoubleComplex* angular_exp_coeff_d,
                                  hipDoubleComplex* angular_exp_coeff_rad_der_d, hipDoubleComplex* angular_exp_coeff_azi_der_d,
                                  hipDoubleComplex* angular_exp_coeff_pol_der_d, double* rjs_d, double rcut_max, int* k2_i_site_d,
                                  int2* index_d, int n_sites, int n_atom_pairs, int k_max, int n_max, int l_max,
                                  hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_soap_soap_der_one", n_atom_pairs * n_max * (l_max + 1), stream,
      TG_LAMBDA(const int t) {
        int k2 = t % n_atom_pairs;
        int l = (t / n_atom_pairs) % (l_max + 1);
        int n = t / n_atom_pairs / (l_max + 1) + 1;
        int i_site = k2_i_site_d[k2] - 1;
        int k0 = l * (l + 1) / 2;
        double Pi = 4.0 * acos(-1.0);
        bool in_cut = rjs_d[k2] < rcut_max;
        double rad_n = radial_exp_coeff_d[n - 1 + k2 * n_max];
        double rad_der_n = radial_exp_coeff_der_d[n - 1 + k2 * n_max];
        for (int np0 = n; np0 <= n_max; np0 += DER_ONE_NP) {
          int2 idx[DER_ONE_NP];
          double rad[DER_ONE_NP];
          double azi[DER_ONE_NP];
          double pol[DER_ONE_NP];
          double rad_np[DER_ONE_NP];
          double rad_der_np[DER_ONE_NP];
#pragma unroll
          for (int b = 0; b < DER_ONE_NP; b++) {
            int np = np0 + b;
            idx[b] = make_int2(-1, 0);
            rad_np[b] = 0;
            rad_der_np[b] = 0;
            if (np <= n_max) {
              idx[b] = index_d[l + (l_max + 1) * (np - 1 + (n - 1) * n_max)];
              rad_np[b] = radial_exp_coeff_d[np - 1 + k2 * n_max];
              rad_der_np[b] = radial_exp_coeff_der_d[np - 1 + k2 * n_max];
            }
            rad[b] = 0;
            azi[b] = 0;
            pol[b] = 0;
          }
          for (int m = 0; m <= l; m++) {
            int k = 1 + k0 + m;
            hipDoubleComplex zero = make_hipDoubleComplex(0.0, 0.0);
            hipDoubleComplex ang = angular_exp_coeff_d[k2 + n_atom_pairs * (k - 1)];
            hipDoubleComplex ang_rad = angular_exp_coeff_rad_der_d[k2 + n_atom_pairs * (k - 1)];
            hipDoubleComplex ang_azi = angular_exp_coeff_azi_der_d[k2 + n_atom_pairs * (k - 1)];
            hipDoubleComplex ang_pol = angular_exp_coeff_pol_der_d[k2 + n_atom_pairs * (k - 1)];
            hipDoubleComplex tmp_1_cnk_d = cnk_d[i_site + n_sites * (k - 1 + (n - 1) * k_max)];
            hipDoubleComplex tmp_1_cnk_rad_d = in_cut ? cnk_rad_der_of(ang, ang_rad, rad_n, rad_der_n, Pi) : zero;
            hipDoubleComplex tmp_1_cnk_azi_d = in_cut ? cnk_ang_der_of(ang_azi, rad_n, Pi) : zero;
            hipDoubleComplex tmp_1_cnk_pol_d = in_cut ? cnk_ang_der_of(ang_pol, rad_n, Pi) : zero;
#pragma unroll
            for (int b = 0; b < DER_ONE_NP; b++) {
              if (idx[b].x < 0)
                continue;
              int np = np0 + b;
              double mult = multiplicity_array_d[idx[b].y + m];
              hipDoubleComplex tmp_2_cnk_d = cnk_d[i_site + n_sites * (k - 1 + (np - 1) * k_max)];
              hipDoubleComplex tmp_2_cnk_rad_d = in_cut ? cnk_rad_der_of(ang, ang_rad, rad_np[b], rad_der_np[b], Pi) : zero;
              hipDoubleComplex tmp_2_cnk_azi_d = in_cut ? cnk_ang_der_of(ang_azi, rad_np[b], Pi) : zero;
              hipDoubleComplex tmp_2_cnk_pol_d = in_cut ? cnk_ang_der_of(ang_pol, rad_np[b], Pi) : zero;
              rad[b] += mult * (tmp_1_cnk_rad_d.x * tmp_2_cnk_d.x + tmp_1_cnk_rad_d.y * tmp_2_cnk_d.y +
                                tmp_1_cnk_d.x * tmp_2_cnk_rad_d.x + tmp_1_cnk_d.y * tmp_2_cnk_rad_d.y);
              azi[b] += mult * (tmp_1_cnk_azi_d.x * tmp_2_cnk_d.x + tmp_1_cnk_azi_d.y * tmp_2_cnk_d.y +
                                tmp_1_cnk_d.x * tmp_2_cnk_azi_d.x + tmp_1_cnk_d.y * tmp_2_cnk_azi_d.y);
              pol[b] += mult * (tmp_1_cnk_pol_d.x * tmp_2_cnk_d.x + tmp_1_cnk_pol_d.y * tmp_2_cnk_d.y +
                                tmp_1_cnk_d.x * tmp_2_cnk_pol_d.x + tmp_1_cnk_d.y * tmp_2_cnk_pol_d.y);
            }
          }
#pragma unroll
          for (int b = 0; b < DER_ONE_NP; b++) {
            if (idx[b].x < 0)
              continue;
            trans_soap_rad_der_d[k2 + idx[b].x * n_atom_pairs] = rad[b];
            trans_soap_azi_der_d[k2 + idx[b].x * n_atom_pairs] = azi[b];
            trans_soap_pol_der_d[k2 + idx[b].x * n_atom_pairs] = pol[b];
          }
        }
      },
      tpb);
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


extern "C" void gpu_get_soap_der(double* soap_d, double* sqrt_dot_d, double3* soap_cart_der_d, double* soap_rad_der_d,
                                 double* soap_azi_der_d, double* soap_pol_der_d, double* thetas_d, double* phis_d, double* rjs_d,
                                 double* multiplicity_array_d, hipDoubleComplex* cnk_d, double* radial_exp_coeff_d,
                                 double* radial_exp_coeff_der_d, hipDoubleComplex* angular_exp_coeff_d,
                                 hipDoubleComplex* angular_exp_coeff_rad_der_d, hipDoubleComplex* angular_exp_coeff_azi_der_d,
                                 hipDoubleComplex* angular_exp_coeff_pol_der_d, double rcut_max, int* n_neigh_d, int* i_k2_start_d,
                                 int* k2_i_site_d, int* k3_index_d, bool* skip_soap_component_d, int n_sites, int n_atom_pairs,
                                 int n_soap, int k_max, int n_max, int l_max, int maxneigh, hipStream_t* stream) {
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


  int2* der_one_index_d;
  hipMallocAsync((void**) &der_one_index_d, sizeof(int2) * n_max * n_max * (l_max + 1), stream[0]);
  cuda_get_soap_index(der_one_index_d, skip_soap_component_d, n_max, l_max, stream);
  cuda_get_soap_der_one(multiplicity_array_d, trans_soap_rad_der_d, trans_soap_azi_der_d, trans_soap_pol_der_d, cnk_d,
                        radial_exp_coeff_d, radial_exp_coeff_der_d, angular_exp_coeff_d, angular_exp_coeff_rad_der_d,
                        angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d, rjs_d, rcut_max, k2_i_site_d, der_one_index_d,
                        n_sites, n_atom_pairs, k_max, n_max, l_max, stream);
  hipFreeAsync(der_one_index_d, stream[0]);

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

// The 1-based site of each pair from each site's 0-based first pair, and the
// species mask from Fortran logicals: per-pair loops get_soap used to run.
extern "C" void gpu_soap_pair_site(int n_sites, int* k2_start_d, int* n_neigh_d, int* k2_i_site_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_soap_pair_site", n_sites, stream,
      TG_LAMBDA(const int i) {
        for (int k = k2_start_d[i]; k < k2_start_d[i] + n_neigh_d[i]; k++)
          k2_i_site_d[k] = i + 1;
      },
      tpb);
}

extern "C" void gpu_logical_to_bool(int* logical_d, bool* bool_d, int n, hipStream_t* stream) {
  tg_parallel_for("turbogap_soap_logical_to_bool", n, stream, TG_LAMBDA(const int k) { bool_d[k] = logical_d[k] != 0; }, tpb);
}

// The expansion of one batch, kept by get_soap in place of the descriptor
// derivatives, for gpu_soap_pair_contract. Everything but soap_d and the
// multiplicities was get_soap's and is freed by gpu_soap_expansion_release.
struct soap_expansion {
  double* soap_d;
  double* sqrt_dot_d;
  hipDoubleComplex* cnk_d;
  double* rad_d;
  double* rad_der_d;
  hipDoubleComplex* ang_d;
  hipDoubleComplex* ang_rad_d;
  hipDoubleComplex* ang_azi_d;
  hipDoubleComplex* ang_pol_d;
  double* rjs_d;
  double* thetas_d;
  double* phis_d;
  int* k2_i_site_d;
  int* k2_start_d;
  int* n_neigh_d;
  double* multiplicity_d;
  int2* index_d;
  double rcut_max;
  int n_sites;
  int n_pairs;
  int n_soap;
  int k_max;
  int n_max;
  int l_max;
};

static soap_expansion kept = {};

extern "C" void gpu_soap_keep_expansion(double* soap_d, double* sqrt_dot_d, hipDoubleComplex* cnk_d, double* rad_d,
                                        double* rad_der_d, hipDoubleComplex* ang_d, hipDoubleComplex* ang_rad_d,
                                        hipDoubleComplex* ang_azi_d, hipDoubleComplex* ang_pol_d, double* rjs_d, double* thetas_d,
                                        double* phis_d, int* k2_i_site_d, int* k2_start_d, int* n_neigh_d, double* multiplicity_d,
                                        bool* skip_soap_component_d, double rcut_max, int n_sites, int n_pairs, int n_soap,
                                        int k_max, int n_max, int l_max, hipStream_t* stream) {
  kept = {soap_d,    sqrt_dot_d, cnk_d,    rad_d,   rad_der_d,   ang_d,      ang_rad_d, ang_azi_d,
          ang_pol_d, rjs_d,      thetas_d, phis_d,  k2_i_site_d, k2_start_d, n_neigh_d, multiplicity_d,
          nullptr,   rcut_max,   n_sites,  n_pairs, n_soap,      k_max,      n_max,     l_max};
  hipMallocAsync((void**) &kept.index_d, sizeof(int2) * n_max * n_max * (l_max + 1), stream[0]);
  cuda_get_soap_index(kept.index_d, skip_soap_component_d, n_max, l_max, stream);
}

extern "C" void gpu_soap_expansion_release(hipStream_t* stream) {
  // All but index_d came from gpu_malloc_async and leave its ledger here.
  void* owned[] = {kept.sqrt_dot_d, kept.cnk_d, kept.rad_d,    kept.rad_der_d, kept.ang_d,       kept.ang_rad_d,  kept.ang_azi_d,
                   kept.ang_pol_d,  kept.rjs_d, kept.thetas_d, kept.phis_d,    kept.k2_i_site_d, kept.k2_start_d, kept.n_neigh_d};
  for (void* p : owned)
    if (p)
      cuda_free_async(&p, stream);
  if (kept.index_d)
    hipFreeAsync(kept.index_d, stream[0]);
  kept = {};
}

// pair_d(1:3, k) = sign * sum_s G(i, s) dq_s(i)/dr_k, i the site of pair k, for
// G(n_sites, n_soap) a gradient with respect to the normalised descriptor q.
// The same contraction as G against soap_cart_der, but dq/dr is never formed:
// G is carried back through the normalisation and the power spectrum to one
// complex vector per site and expansion coefficient, H(n, lm), and each pair
// then needs only its own coefficient derivatives against H.
extern "C" void gpu_soap_pair_contract(double* G_d, double sign, double* pair_d, hipStream_t* stream) {
  soap_expansion e = kept;
  int n_sites = e.n_sites;
  int n_pairs = e.n_pairs;
  int n_soap = e.n_soap;
  int k_max = e.k_max;
  int n_max = e.n_max;
  int l_max = e.l_max;

  // Through q = p/|p|: dE/dp = (G - q (q.G)) / |p|.
  double* Gp_d;
  hipMallocAsync((void**) &Gp_d, sizeof(double) * n_sites * n_soap, stream[0]);
  tg_parallel_for(
      "turbogap_soap_contract_norm", n_sites, stream,
      TG_LAMBDA(const int i) {
        double qg = 0.0;
        for (int s = 0; s < n_soap; s++)
          qg += e.soap_d[s + i * n_soap] * G_d[i + s * n_sites];
        double inv = 1.0 / e.sqrt_dot_d[i];
        for (int s = 0; s < n_soap; s++)
          Gp_d[i + s * n_sites] = (G_d[i + s * n_sites] - e.soap_d[s + i * n_soap] * qg) * inv;
      },
      tpb);

  // Through p(n, np, l) = sum_m mult cnk(n).cnk(np): the diagonal is the
  // product of one coefficient with itself, so it carries twice.
  hipDoubleComplex* H_d;
  hipMallocAsync((void**) &H_d, sizeof(hipDoubleComplex) * n_sites * k_max * n_max, stream[0]);
  tg_parallel_for(
      "turbogap_soap_contract_h", n_sites * n_max * (l_max + 1), stream,
      TG_LAMBDA(const int t) {
        int i = t % n_sites;
        int l = (t / n_sites) % (l_max + 1);
        int n = t / n_sites / (l_max + 1) + 1;
        for (int m = 0; m <= l; m++) {
          int k = 1 + l * (l + 1) / 2 + m;
          hipDoubleComplex h = make_hipDoubleComplex(0.0, 0.0);
          for (int np = 1; np <= n_max; np++) {
            int lo = n < np ? n : np;
            int hi = n < np ? np : n;
            int2 idx = e.index_d[l + (l_max + 1) * (hi - 1 + (lo - 1) * n_max)];
            if (idx.x < 0)
              continue;
            double w = Gp_d[i + idx.x * n_sites] * e.multiplicity_d[idx.y + m];
            if (n == np)
              w *= 2.0;
            hipDoubleComplex c = e.cnk_d[i + n_sites * (k - 1 + (np - 1) * k_max)];
            h.x += w * c.x;
            h.y += w * c.y;
          }
          H_d[n - 1 + n_max * (k - 1 + k_max * i)] = h;
        }
      },
      tpb);

  // Each pair's coefficient derivatives against H, one thread per (pair, n)
  // with n fastest so that rad and H are read contiguously across a warp; then
  // per pair, the sum over n in order and the Cartesian conversion of
  // cuda_get_soap_der_thr_one. The central pair is redone below.
  double* part_d;
  hipMallocAsync((void**) &part_d, sizeof(double) * 3 * n_pairs * n_max, stream[0]);
  tg_parallel_for_long(
      "turbogap_soap_contract_pair_n", (long long) n_pairs * n_max, stream,
      TG_LAMBDA(const long long t) {
        int n = t % n_max + 1;
        int k2 = t / n_max;
        int i = e.k2_i_site_d[k2] - 1;
        double rad = e.rad_d[n - 1 + k2 * n_max];
        double rad_der = e.rad_der_d[n - 1 + k2 * n_max];
        double f_rad = 0.0;
        double f_azi = 0.0;
        double f_pol = 0.0;
        for (int k = 1; k <= k_max; k++) {
          hipDoubleComplex h = H_d[n - 1 + n_max * (k - 1 + k_max * i)];
          hipDoubleComplex ang = e.ang_d[k2 + n_pairs * (k - 1)];
          hipDoubleComplex ang_rad = e.ang_rad_d[k2 + n_pairs * (k - 1)];
          hipDoubleComplex ang_azi = e.ang_azi_d[k2 + n_pairs * (k - 1)];
          hipDoubleComplex ang_pol = e.ang_pol_d[k2 + n_pairs * (k - 1)];
          f_rad += rad_der * (ang.x * h.x + ang.y * h.y) + rad * (ang_rad.x * h.x + ang_rad.y * h.y);
          f_azi += rad * (ang_azi.x * h.x + ang_azi.y * h.y);
          f_pol += rad * (ang_pol.x * h.x + ang_pol.y * h.y);
        }
        part_d[3 * t] = f_rad;
        part_d[3 * t + 1] = f_azi;
        part_d[3 * t + 2] = f_pol;
      },
      tpb);
  tg_parallel_for(
      "turbogap_soap_contract_pair", n_pairs, stream,
      TG_LAMBDA(const int k2) {
        double f_rad = 0.0;
        double f_azi = 0.0;
        double f_pol = 0.0;
        if (e.rjs_d[k2] < e.rcut_max) {
          double Pi = 4.0 * acos(-1.0);
          for (int n = 0; n < n_max; n++) {
            f_rad += part_d[3 * (n + (long long) n_max * k2)];
            f_azi += part_d[3 * (n + (long long) n_max * k2) + 1];
            f_pol += part_d[3 * (n + (long long) n_max * k2) + 2];
          }
          f_rad *= Pi;
          f_azi *= Pi;
          f_pol *= Pi;
        }
        double theta = e.thetas_d[k2];
        double phi = e.phis_d[k2];
        double rj = e.rjs_d[k2];
        pair_d[3 * k2] = sign * (sin(theta) * cos(phi) * f_rad - cos(theta) * cos(phi) / rj * f_pol - sin(phi) / rj * f_azi);
        pair_d[3 * k2 + 1] = sign * (sin(theta) * sin(phi) * f_rad - cos(theta) * sin(phi) / rj * f_pol + cos(phi) / rj * f_azi);
        pair_d[3 * k2 + 2] = sign * (cos(theta) * f_rad + sin(theta) / rj * f_pol);
      },
      tpb);
  hipFreeAsync(part_d, stream[0]);

  // The central pair is minus the sum of the site's others, in their order.
  tg_parallel_for(
      "turbogap_soap_contract_self", n_sites, stream,
      TG_LAMBDA(const int i) {
        int k3 = e.k2_start_d[i];
        for (int c = 0; c < 3; c++) {
          double s = 0.0;
          for (int j = 1; j < e.n_neigh_d[i]; j++)
            s -= pair_d[3 * (k3 + j) + c];
          pair_d[3 * k3 + c] = s;
        }
      },
      tpb);

  hipFreeAsync(Gp_d, stream[0]);
  hipFreeAsync(H_d, stream[0]);
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

      my_cnk_rad_der = cnk_rad_der_of(my_ang_exp_c, my_ang_exp_c_rad_der, my_radial_exp_c, my_radial_exp_c_der, Pi);
      my_cnk_azi_der = cnk_ang_der_of(my_ang_exp_c_azi_der, my_radial_exp_c, Pi);
      my_cnk_pol_der = cnk_ang_der_of(my_ang_exp_c_pol_der, my_radial_exp_c, Pi);

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
