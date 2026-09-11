// The soap_turbo radial expansion: poly3 and poly3gauss basis coefficients and
// their derivatives, plus the global scaling applied to them.
//
// This is the single most expensive kernel family in a GAP run.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
#define mode_polynomial 1

// A note on why there is almost no pow() left in this file.
//
// nvcc does NOT fold pow(x, 2) into x * x, even with a literal exponent: the
// PTX for this translation unit carried 27 real calls to
// __internal_accurate_pow -- a log2, an exp2 and their correction terms, a few
// hundred FP64 instructions, to square a number. The 51 literal-exponent calls
// are now multiplies and pow(pi, 0.25) is sqrt(sqrt(pi)).
//
// Every replacement is PARENTHESISED. `1.0 / pow(x, 2)` written as
// `1.0 / x * x` is 1.0, because / and * have equal precedence and associate
// left. The unparenthesised version of this change built cleanly, ran, and put
// an 80% error in energy_soap; the regression suite is what found it.
//
// The four that remain cannot be done this way: pow(tmp1, amplitude_scaling)
// takes its exponent from a per-species input.


__device__ double N_a(double rcut, int a) {
  const int b = 2 * a + 5;
  return sqrt(rcut / static_cast<double>(b));
}


__global__ void check_nan(double* G, int Nt) {
  int id = threadIdx.x + threadIdx.x + blockIdx.x * blockDim.x;
  {
    if (isnan(G[id]) && id < Nt) {
      printf("Is nan %lf at %d", G[id], id);
    }
  }
}


__global__ void cuda_global_scaling(double* radial_exp_coeff_d, int* i_beg_d, int* i_end_d, double* global_scaling_d, int n_max,
                                    int n_atom_pairs, int n_species, double* rcut_hard_d, int* k2_i_site_d, int* k2_start_d,
                                    int divide) {
  int i_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (i_ij < n_atom_pairs) {
    //int i_site=k2_i_site_d[i_ij];
    //int k2start=k2_start_d[i_site-1];
    //if(i_ij!=k2start)
    {
      int i_one = 0;
      for (int i = 0; i < n_species; i++) {
        for (int ii = i_beg_d[i]; ii <= i_end_d[i]; ii++) {
          double loc_rad_exp_coeff =
              radial_exp_coeff_d[i_one + i_ij * n_max] *
              global_scaling_d[i]; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]*global_scaling_d[i];
          // The value and its radial derivative come out of the radial
          // kernels in different units (the derivative is with respect to
          // the reduced coordinate r/rcut_hard), hence the opposite powers
          // of sqrt(rcut_hard) here. Verified: forcing both to the same
          // factor makes cnk_rad_der far worse.
          if (divide == 0) {
            loc_rad_exp_coeff *= sqrt(rcut_hard_d[i]);
          }
          if (divide == 1) {
            loc_rad_exp_coeff *= 1.0 / sqrt(rcut_hard_d[i]);
          }
          radial_exp_coeff_d[i_one + i_ij * n_max] =
              loc_rad_exp_coeff; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]=loc_rad_exp_coeff;

          i_one++;
        }
      }
    }
  }
}


// cuda_poly3gauss_one used to be launched here. Its body had been commented
// out, so it read nothing, wrote nothing and only ran an empty loop nest.
// Removed along with the kernel; the real work is done by cuda_global_scaling
// below.

extern "C" void gpu_get_radial_exp_coeff_poly3gauss(double* radial_exp_coeff_d, double* radial_exp_coeff_der_d, int* i_beg_d,
                                                    int* i_end_d, double* global_scaling_d, int size_radial_exp_coeff_one,
                                                    int size_radial_exp_coeff_two, int n_species, bool c_do_derivatives,
                                                    int bintybint, double* rcut_hard_d, int* k2_i_site_d, int* k_2start_d,
                                                    hipStream_t* stream) {
  dim3 nblocks = dim3((size_radial_exp_coeff_two - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  int divide;
  divide = 0;
  cuda_global_scaling<<<nblocks, nthreads, 0, stream[0]>>>(radial_exp_coeff_d, i_beg_d, i_end_d, global_scaling_d,
                                                           size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species,
                                                           rcut_hard_d, k2_i_site_d, k_2start_d, divide);
  /* gpuErrchk( hipPeekAtLastError() );
       gpuErrchk( hipDeviceSynchronize() ); */
  if (c_do_derivatives) {
    divide = 1;
    cuda_global_scaling<<<nblocks, nthreads, 0, stream[0]>>>(radial_exp_coeff_der_d, i_beg_d, i_end_d, global_scaling_d,
                                                             size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species,
                                                             rcut_hard_d, k2_i_site_d, k_2start_d, divide);
  }
  /* gpuErrchk( hipPeekAtLastError() );
       gpuErrchk( hipDeviceSynchronize() );   */
}

__global__ void kernel_get_radial_poly3gauss(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d,
                                             int n_sites, int* n_neigh_d, int n_max, int n_temp, bool do_derivatives,
                                             double* exp_coeff_d, double* exp_coeff_der_d, double* rcut_soft_d,
                                             double* atom_sigma_d, double* exp_coeff_temp1_d, double* exp_coeff_temp2_d,
                                             double* exp_coeff_der_temp_d, int* i_beg_d, int* i_end_d, double* atom_sigma_scaling_d,
                                             int mode, int radial_enhancement, double* amplitude_scaling_d, int* alpha_max_d,
                                             double* nf_d, int n_temp_der, double* W_d) {
  int n, d, k;
  double ampli_tude, ampli_tude_der, atom_sigma_scaled, amplitude_scaling, C1, C2, W_exp, nf, atom_sigma_f, rj_f, sf2;
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (k_ij < n_atom_pairs) {
    double rjs = rjs_d[k_ij];
    for (int i_sp = 0; i_sp < n_species; i_sp++) {
      if (mask_d[k_ij + i_sp * n_atom_pairs]) {
        double rcut_hard_in = rcut_hard_d[i_sp];
        if (rjs <= rcut_hard_in) {
          int alpha_max = alpha_max_d[i_sp];
          int alpha_max_der = alpha_max;
          if (do_derivatives)
            alpha_max_der += 2;
          int j = 0;
          int j1 = 0;
          for (n = 0; n < n_sites; n++) {
            j1 += n_neigh_d[n];
            if (k_ij < j1)
              break;
            j = j1;
          }
          if (k_ij > j) {
            double pi = acos(-1.0);
            double sq2 = sqrt(2.0);
            double rcut_soft_in = rcut_soft_d[i_sp];
            double rcut_soft = rcut_soft_in / rcut_hard_in;
            double rcut_hard = 1.0;
            double atom_sigma_in = atom_sigma_d[i_sp];
            double atom_sigma = atom_sigma_in / rcut_hard_in;
            double dr = 1.0 - rcut_soft_in / rcut_hard_in;
            double N_gauss = sqrt(2.0 / atom_sigma) / sqrt(sqrt(pi));
            double pref_f = 0.0;
            for (n = 0; n < alpha_max_der; n++) {
              exp_coeff_temp1_d[k_ij * n_temp + n] = 0.0;
              exp_coeff_temp2_d[k_ij * n_temp + n] = 0.0;
            }
            if (do_derivatives)
              for (n = 0; n < alpha_max; n++)
                exp_coeff_der_temp_d[k_ij * n_temp_der + n] = 0.0;
            double rj = rjs / rcut_hard_in;
            double atom_sigma_scaling = atom_sigma_scaling_d[i_sp];
            double atom_sigma_scaled = atom_sigma + atom_sigma_scaling * rj;
            double s2 = (atom_sigma_scaled * atom_sigma_scaled);
            double amplitude_scaling = amplitude_scaling_d[i_sp];
            tmp1 = 1.0 + rj * rj * (2.0 * rj - 3.0);
            tmp2 = atom_sigma_scaling / atom_sigma_scaled;
            tmp3 = 6.0 / atom_sigma_scaled * rj * (rj - 1.0);
            if (mode == mode_polynomial) {
              if (amplitude_scaling == 0.0) {
                ampli_tude = 1.0 / atom_sigma_scaled;
                ampli_tude_der = -atom_sigma_scaling / s2;
              } else if (tmp1 <= 1.e-10) {
                ampli_tude = 0.0;
                ampli_tude_der = 0.0;
              } else {
                if (amplitude_scaling == 1.0) {
                  ampli_tude = 1.0 / atom_sigma_scaled * tmp1;
                  ampli_tude_der = tmp3 - tmp2 * ampli_tude;
                } else {
                  ampli_tude = 1.0 / atom_sigma_scaled * pow(tmp1, amplitude_scaling);
                  // The trailing term is a product, matching the CPU:
                  //   - atom_sigma_scaling/atom_sigma_scaled * amplitude
                  // It was written as "- tmp2 - ampli_tude", subtracting the two
                  // factors instead of multiplying them. Only this branch was
                  // affected; the amplitude_scaling == 1 case above is correct.
                  ampli_tude_der = amplitude_scaling * tmp3 * pow(tmp1, amplitude_scaling - 1.0) - tmp2 * ampli_tude;
                }
              }
            }
            tmp3 = rj + sqrt(2.0 / pi) * atom_sigma_scaled;
            tmp4 = sqrt(8.0 / pi) * atom_sigma_scaled;
            tmp5 = rj * rj + s2 + tmp4 * rj;
            if (radial_enhancement == 1) {
              ampli_tude_der = ampli_tude * (1.0 + sqrt(2.0 / pi) * atom_sigma_scaling) + ampli_tude_der * tmp3;
              ampli_tude = ampli_tude * tmp3;
            } else if (radial_enhancement == 2) {
              ampli_tude_der = ampli_tude * (2.0 * rj + 2.0 * atom_sigma_scaled * atom_sigma_scaling + tmp4 +
                                             sqrt(8.0 / pi) * rj * atom_sigma_scaling) +
                               ampli_tude_der * tmp5;
              ampli_tude = ampli_tude * tmp5;
            }
            double I_n = 0.0;
            double N_n = 1.0;
            double N_np1 = N_a(rcut_hard, -2);
            double I_np1 = sqrt(pi / 2.0) * atom_sigma_scaled *
                           (erf((rcut_soft - rj) / sq2 / atom_sigma_scaled) - erf((-rj) / sq2 / atom_sigma_scaled)) / N_np1;
            double I_np2, N_np2;
            C1 = (rcut_hard_in == rcut_soft_in) ? 0.0 : s2 / dr * exp(-0.5 * ((rcut_soft - rj) * (rcut_soft - rj)) / s2);
            C2 = s2 / rcut_hard * exp(-0.5 * (rj * rj) / s2);
            for (n = -1; n <= alpha_max_der - 1; n++) {
              C1 = C1 * dr;
              C2 = C2 * rcut_hard;
              N_np2 = N_a(rcut_hard, n);
              I_np2 = s2 * double(n + 1) * N_n / N_np2 * I_n - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 + C1 / N_np2 - C2 / N_np2;
              if (n > 0)
                exp_coeff_temp1_d[k_ij * n_temp + n - 1] = I_np2;
              N_n = N_np1;
              N_np1 = N_np2;
              I_n = I_np1;
              I_np1 = I_np2;
            }
            if (do_derivatives) {
              tmp1 = atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled;
              tmp2 = (rj - rcut_hard) / s2 * (tmp1 - 1.0);
              tmp3 = rcut_hard * (2.0 * tmp1 - 1.0) / s2;
              tmp4 = atom_sigma_scaling * rcut_hard * rcut_hard / (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled);
              tmp5 = exp_coeff_temp1_d[k_ij * n_temp];
              tmp6 = exp_coeff_temp1_d[k_ij * n_temp + 1];
              for (n = 1; n <= alpha_max - 1; n++) {
                tmp7 = exp_coeff_temp1_d[k_ij * n_temp + n + 1];
                exp_coeff_der_temp_d[k_ij * n_temp_der + n - 1] = tmp2 * tmp5 +
                                                                  tmp3 * N_a(rcut_hard, n + 1) / N_a(rcut_hard, n) * tmp6 +
                                                                  tmp4 * N_a(rcut_hard, n + 2) / N_a(rcut_hard, n) * tmp7;
                tmp5 = tmp6;
                tmp6 = tmp7;
              }
            }
            if (false || (rcut_soft - rj) < 4.0 * atom_sigma_scaled) {
              nf = nf_d[i_sp];
              tmp1 = dr * dr / nf / nf;
              atom_sigma_f = atom_sigma_scaled * dr / nf / sqrt(s2 + tmp1);
              rj_f = (s2 * rcut_soft + tmp1 * rj) / (s2 + tmp1);
              sf2 = (atom_sigma_f * atom_sigma_f);
              pref_f = exp(-0.5 * ((rcut_soft - rj) * (rcut_soft - rj)) / (s2 + tmp1));
              I_n = 0.0;
              N_n = 1.0;
              N_np1 = N_a(rcut_hard, -2);
              I_np1 = sqrt(pi / 2.0) * atom_sigma_f *
                      (erf((rcut_hard - rj_f) / sq2 / atom_sigma_f) - erf((rcut_soft - rj_f) / sq2 / atom_sigma_f)) / N_np1;
              C2 = sf2 / dr * exp(-0.5 * ((rcut_soft - rj_f) * (rcut_soft - rj_f)) / sf2);
              for (n = -1; n <= alpha_max_der - 1; n++) {
                C2 *= dr;
                double N_np2 = N_a(rcut_hard, n);
                double I_np2 = sf2 * double(n + 1) * N_n / N_np2 * I_n - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 - C2 / N_np2;
                if (n > 0)
                  exp_coeff_temp2_d[k_ij * n_temp + n - 1] += I_np2;
                N_n = N_np1;
                N_np1 = N_np2;
                I_n = I_np1;
                I_np1 = I_np2;
              }
              if (do_derivatives) {
                double denom = s2 + tmp1;
                double der_pref_f = pref_f * ((rcut_soft - rj) / denom + ((rcut_soft - rj) * (rcut_soft - rj)) / (denom * denom) *
                                                                             atom_sigma_scaled * atom_sigma_scaling);
                double der_rjf_rj = (2.0 * atom_sigma_scaled * rcut_soft * atom_sigma_scaling + tmp1) / denom -
                                    (s2 * rcut_soft + tmp1 * rj) * 2.0 * atom_sigma_scaled * atom_sigma_scaling / (denom * denom);
                double der_sjf_rj =
                    atom_sigma_scaling * dr / nf / sqrt(denom) * (1.0 - (atom_sigma_scaled * atom_sigma_scaled) / denom);
                tmp2 = (rj_f - rcut_hard) / sf2 * (der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj);
                tmp3 = rcut_hard / sf2 * (2.0 * der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj);
                tmp4 = der_sjf_rj * rcut_hard * rcut_hard / (atom_sigma_f * atom_sigma_f * atom_sigma_f);
                tmp5 = exp_coeff_temp2_d[k_ij * n_temp];
                tmp6 = exp_coeff_temp2_d[k_ij * n_temp + 1];
                for (n = 1; n <= alpha_max - 1; n++) {
                  tmp7 = exp_coeff_temp2_d[k_ij * n_temp + n + 1];
                  exp_coeff_der_temp_d[k_ij * n_temp_der + n - 1] +=
                      pref_f * (tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n + 1) / N_a(rcut_hard, n) * tmp6 +
                                tmp4 * N_a(rcut_hard, n + 2) / N_a(rcut_hard, n) * tmp7) +
                      der_pref_f * tmp5;
                  tmp5 = tmp6;
                  tmp6 = tmp7;
                }
              }
            }
            exp_coeff_temp1_d[k_ij * n_temp + alpha_max - 1] = 0.0;
            exp_coeff_temp2_d[k_ij * n_temp + alpha_max - 1] = 0.0;
            if (false || rj < 4.0 * (atom_sigma + atom_sigma_scaled)) {
              double sigma_star = sqrt((atom_sigma * atom_sigma) + s2);
              exp_coeff_temp1_d[k_ij * n_temp + alpha_max - 1] =
                  exp(-0.5 * (rj * rj) / (sigma_star * sigma_star)) * sqrt(pi / 2.0) * atom_sigma_scaled * atom_sigma / sigma_star *
                  (1.0 + erf(atom_sigma / atom_sigma_scaled * rj / sq2 / sigma_star)) * N_gauss;
              if (do_derivatives)
                exp_coeff_der_temp_d[k_ij * n_temp_der + alpha_max - 1] =
                    ((rj * rj) * atom_sigma_scaling / (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled) -
                     rj / (sigma_star * sigma_star) +
                     atom_sigma_scaling * (rj * rj) * (atom_sigma * atom_sigma * atom_sigma * atom_sigma) /
                         (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled) /
                         (sigma_star * sigma_star * sigma_star * sigma_star) +
                     atom_sigma_scaling * (atom_sigma * atom_sigma) / atom_sigma_scaled / (sigma_star * sigma_star) -
                     2.0 * (rj * rj) * atom_sigma_scaling * (atom_sigma * atom_sigma) /
                         (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled) / (sigma_star * sigma_star)) *
                        exp_coeff_temp1_d[k_ij * n_temp + alpha_max - 1] +
                    (1. / s2 - 2.0 * rj * atom_sigma_scaling / (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled)) * s2 *
                        (atom_sigma * atom_sigma) / (sigma_star * sigma_star) * sqrt(2.0 / atom_sigma) / sqrt(sqrt(pi)) *
                        exp(-0.5 * (rj * rj) / (sigma_star * sigma_star) * (1.0 + (atom_sigma * atom_sigma) / s2)) +
                    sqrt(2.0 / atom_sigma) / sqrt(sqrt(pi)) *
                        exp(-0.5 * (rj * rj) / (sigma_star * sigma_star) * (1.0 + (atom_sigma * atom_sigma) / s2)) *
                        atom_sigma_scaling / atom_sigma_scaled * rj * (atom_sigma * atom_sigma * atom_sigma * atom_sigma) /
                        (sigma_star * sigma_star * sigma_star * sigma_star);
            }
            if (do_derivatives) {
              for (n = 0; n < alpha_max; n++)
                exp_coeff_der_temp_d[k_ij * n_temp_der + n] =
                    ampli_tude * exp_coeff_der_temp_d[k_ij * n_temp_der + n] +
                    ampli_tude_der * (exp_coeff_temp1_d[k_ij * n_temp + n] + pref_f * exp_coeff_temp2_d[k_ij * n_temp + n]);
              for (d = i_beg_d[i_sp]; d <= i_end_d[i_sp]; d++) {
                W_exp = 0.0;
                k = 0;
                for (n = i_beg_d[i_sp]; n <= i_end_d[i_sp]; n++) {
                  W_exp += W_d[(n - 1) * n_max + d - 1] * exp_coeff_der_temp_d[k_ij * n_temp_der + k];
                  k += 1;
                }
                exp_coeff_der_d[k_ij * n_max + d - 1] = W_exp;
              }
            }
            for (d = i_beg_d[i_sp]; d <= i_end_d[i_sp]; d++) {
              W_exp = 0.0;
              k = 0;
              for (n = i_beg_d[i_sp]; n <= i_end_d[i_sp]; n++) {
                W_exp += W_d[(n - 1) * n_max + d - 1] *
                         (exp_coeff_temp1_d[k_ij * n_temp + k] + pref_f * exp_coeff_temp2_d[k_ij * n_temp + k]);
                k += 1;
              }
              exp_coeff_d[k_ij * n_max + d - 1] = ampli_tude * W_exp;
            }
          }
        }
      }
    }
  }
}

extern "C" void gpu_radial_poly3gauss(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d,
                                      int n_sites, int* n_neigh_d, int n_max, int n_temp, bool do_derivatives, double* exp_coeff_d,
                                      double* exp_coeff_der_d, double* rcut_soft_d, double* atom_sigma_d, double* exp_coeff_temp1_d,
                                      double* exp_coeff_temp2_d, double* exp_coeff_der_temp_d, int* i_beg, int* i_end,
                                      double* atom_sigma_scaling_d, int mode, int radial_enhancement, double* amplitude_scaling_d,
                                      int* alpha_max_d, double* nf_d, int n_temp_der, double* W_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_get_radial_poly3gauss<<<nblocks, nthreads, 0, stream[0]>>>(
      n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, n_temp, do_derivatives, exp_coeff_d,
      exp_coeff_der_d, rcut_soft_d, atom_sigma_d, exp_coeff_temp1_d, exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end,
      atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, nf_d, n_temp_der, W_d);
}

__global__ void kernel_get_radial_poly3(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d,
                                        int n_sites, int* n_neigh_d, int n_max, int n_temp, bool do_derivatives,
                                        double* exp_coeff_d, double* exp_coeff_der_d, double* rcut_soft_d, double* atom_sigma_d,
                                        double* exp_coeff_temp1_d, double* exp_coeff_temp2_d, double* exp_coeff_der_temp_d,
                                        int* i_beg_d, int* i_end_d, double* atom_sigma_scaling_d, int mode, int radial_enhancement,
                                        double* amplitude_scaling_d, int* alpha_max_d, double* nf_d, int n_temp_der, double* W_d,
                                        bool* do_central_d, double* central_weight_d) {
  int n, d, k;
  double ampli_tude, ampli_tude_der, atom_sigma_scaled, amplitude_scaling, C1, C2, W_exp, nf, atom_sigma_f, rj_f, sf2;
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (k_ij < n_atom_pairs) {
    double rjs = rjs_d[k_ij];
    for (int i_sp = 0; i_sp < n_species; i_sp++) {
      if (mask_d[k_ij + i_sp * n_atom_pairs]) {
        double rcut_hard_in = rcut_hard_d[i_sp];
        if (rjs <= rcut_hard_in) {
          int alpha_max = alpha_max_d[i_sp];
          int alpha_max_der = alpha_max;
          if (do_derivatives)
            alpha_max_der += 2;
          int j = 0;
          int j1 = 0;
          for (n = 0; n < n_sites; n++) {
            j1 += n_neigh_d[n];
            if (k_ij < j1)
              break;
            j = j1;
          }
          if (k_ij > j || do_central_d[i_sp]) {
            double pi = acos(-1.0);
            double sq2 = sqrt(2.0);
            double rcut_soft_in = rcut_soft_d[i_sp];
            double rcut_soft = rcut_soft_in / rcut_hard_in;
            double rcut_hard = 1.0;
            double atom_sigma_in = atom_sigma_d[i_sp];
            double atom_sigma = atom_sigma_in / rcut_hard_in;
            double dr = 1.0 - rcut_soft_in / rcut_hard_in;
            double pref_f = 0.0;
            for (n = 0; n < alpha_max_der; n++) {
              exp_coeff_temp1_d[k_ij * n_temp + n] = 0.0;
              exp_coeff_temp2_d[k_ij * n_temp + n] = 0.0;
            }
            if (do_derivatives)
              for (n = 0; n < alpha_max; n++)
                exp_coeff_der_temp_d[k_ij * n_temp_der + n] = 0.0;
            double rj = rjs / rcut_hard_in;
            double atom_sigma_scaling = atom_sigma_scaling_d[i_sp];
            double atom_sigma_scaled = atom_sigma + atom_sigma_scaling * rj;
            double s2 = (atom_sigma_scaled * atom_sigma_scaled);
            double amplitude_scaling = amplitude_scaling_d[i_sp];
            tmp1 = 1.0 + rj * rj * (2.0 * rj - 3.0);
            tmp2 = atom_sigma_scaling / atom_sigma_scaled;
            tmp3 = 6.0 / atom_sigma_scaled * rj * (rj - 1.0);
            if (mode == mode_polynomial) {
              if (amplitude_scaling == 0.0) {
                ampli_tude = 1.0 / atom_sigma_scaled;
                ampli_tude_der = -atom_sigma_scaling / s2;
              } else if (tmp1 <= 1.e-10) {
                ampli_tude = 0.0;
                ampli_tude_der = 0.0;
              } else {
                if (amplitude_scaling == 1.0) {
                  ampli_tude = 1.0 / atom_sigma_scaled * tmp1;
                  ampli_tude_der = tmp3 - tmp2 * ampli_tude;
                } else {
                  ampli_tude = 1.0 / atom_sigma_scaled * pow(tmp1, amplitude_scaling);
                  ampli_tude_der = tmp3 * amplitude_scaling * pow(tmp1, amplitude_scaling - 1.0) - tmp2 * ampli_tude;
                }
              }
            }
            if (k_ij == j) {
              ampli_tude = central_weight_d[i_sp] * ampli_tude;
              ampli_tude_der = central_weight_d[i_sp] * ampli_tude_der;
            }
            tmp3 = rj + sqrt(2.0 / pi) * atom_sigma_scaled;
            tmp4 = sqrt(8.0 / pi) * atom_sigma_scaled;
            tmp5 = rj * rj + s2 + tmp4 * rj;
            if (radial_enhancement == 1) {
              ampli_tude_der = ampli_tude * (1.0 + sqrt(2.0 / pi) * atom_sigma_scaling) + ampli_tude_der * tmp3;
              ampli_tude = ampli_tude * tmp3;
            } else if (radial_enhancement == 2) {
              ampli_tude_der = ampli_tude * (2.0 * rj + 2.0 * atom_sigma_scaled * atom_sigma_scaling + tmp4 +
                                             sqrt(8.0 / pi) * rj * atom_sigma_scaling) +
                               ampli_tude_der * tmp5;
              ampli_tude = ampli_tude * tmp5;
            }
            double I_n = 0.0;
            double N_n = 1.0;
            double N_np1 = N_a(rcut_hard, -2);
            double I_np1 = sqrt(pi / 2.0) * atom_sigma_scaled *
                           (erf((rcut_soft - rj) / sq2 / atom_sigma_scaled) - erf((-rj) / sq2 / atom_sigma_scaled)) / N_np1;
            double I_np2, N_np2;
            C1 = (rcut_hard_in == rcut_soft_in) ? 0.0 : s2 / dr * exp(-0.5 * ((rcut_soft - rj) * (rcut_soft - rj)) / s2);
            C2 = s2 / rcut_hard * exp(-0.5 * (rj * rj) / s2);
            for (n = -1; n <= alpha_max_der; n++) {
              C1 *= dr;
              C2 *= rcut_hard;
              N_np2 = N_a(rcut_hard, n);
              I_np2 = s2 * double(n + 1) * N_n / N_np2 * I_n - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 + C1 / N_np2 - C2 / N_np2;
              if (n > 0)
                exp_coeff_temp1_d[k_ij * n_temp + n - 1] = I_np2;
              N_n = N_np1;
              N_np1 = N_np2;
              I_n = I_np1;
              I_np1 = I_np2;
            }
            if (do_derivatives) {
              tmp1 = atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled;
              tmp2 = (rj - rcut_hard) / s2 * (tmp1 - 1.0);
              tmp3 = rcut_hard * (2.0 * tmp1 - 1.0) / s2;
              tmp4 = atom_sigma_scaling * rcut_hard * rcut_hard / (atom_sigma_scaled * atom_sigma_scaled * atom_sigma_scaled);
              tmp5 = exp_coeff_temp1_d[k_ij * n_temp];
              tmp6 = exp_coeff_temp1_d[k_ij * n_temp + 1];
              for (n = 1; n <= alpha_max - 1; n++) {
                tmp7 = exp_coeff_temp1_d[k_ij * n_temp + n + 1];
                exp_coeff_der_temp_d[k_ij * n_temp_der + n - 1] = tmp2 * tmp5 +
                                                                  tmp3 * N_a(rcut_hard, n + 1) / N_a(rcut_hard, n) * tmp6 +
                                                                  tmp4 * N_a(rcut_hard, n + 2) / N_a(rcut_hard, n) * tmp7;
                tmp5 = tmp6;
                tmp6 = tmp7;
              }
            }
            if (false || (rcut_soft - rj) < 4.0 * atom_sigma_scaled) {
              nf = nf_d[i_sp];
              tmp1 = dr * dr / nf / nf;
              atom_sigma_f = atom_sigma_scaled * dr / nf / sqrt(s2 + tmp1);
              rj_f = (s2 * rcut_soft + tmp1 * rj) / (s2 + tmp1);
              sf2 = (atom_sigma_f * atom_sigma_f);
              pref_f = exp(-0.5 * ((rcut_soft - rj) * (rcut_soft - rj)) / (s2 + tmp1));
              I_n = 0.0;
              N_n = 1.0;
              N_np1 = N_a(rcut_hard, -2);
              I_np1 = sqrt(pi / 2.0) * atom_sigma_f *
                      (erf((rcut_hard - rj_f) / sq2 / atom_sigma_f) - erf((rcut_soft - rj_f) / sq2 / atom_sigma_f)) / N_np1;
              C2 = sf2 / dr * exp(-0.5 * ((rcut_soft - rj_f) * (rcut_soft - rj_f)) / sf2);
              for (n = -1; n <= alpha_max_der - 1; n++) {
                C2 *= dr;
                double N_np2 = N_a(rcut_hard, n);
                double I_np2 = sf2 * double(n + 1) * N_n / N_np2 * I_n - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 - C2 / N_np2;
                if (n > 0)
                  exp_coeff_temp2_d[k_ij * n_temp + n - 1] += I_np2;
                N_n = N_np1;
                N_np1 = N_np2;
                I_n = I_np1;
                I_np1 = I_np2;
              }
              if (do_derivatives) {
                double denom = s2 + tmp1;
                double der_pref_f = pref_f * ((rcut_soft - rj) / denom + ((rcut_soft - rj) * (rcut_soft - rj)) / (denom * denom) *
                                                                             atom_sigma_scaled * atom_sigma_scaling);
                double der_rjf_rj = (2.0 * atom_sigma_scaled * rcut_soft * atom_sigma_scaling + tmp1) / denom -
                                    (s2 * rcut_soft + tmp1 * rj) * 2.0 * atom_sigma_scaled * atom_sigma_scaling / (denom * denom);
                double der_sjf_rj =
                    atom_sigma_scaling * dr / nf / sqrt(denom) * (1.0 - (atom_sigma_scaled * atom_sigma_scaled) / denom);

                tmp2 = (rj_f - rcut_hard) / sf2 * (der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj);
                tmp3 = rcut_hard / sf2 * (2.0 * der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj);
                tmp4 = der_sjf_rj * rcut_hard * rcut_hard / (atom_sigma_f * atom_sigma_f * atom_sigma_f);
                tmp5 = exp_coeff_temp2_d[k_ij * n_temp];
                tmp6 = exp_coeff_temp2_d[k_ij * n_temp + 1];
                for (n = 1; n <= alpha_max - 1; n++) {
                  tmp7 = exp_coeff_temp2_d[k_ij * n_temp + n + 1];
                  exp_coeff_der_temp_d[k_ij * n_temp_der + n - 1] +=
                      pref_f * (tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n + 1) / N_a(rcut_hard, n) * tmp6 +
                                tmp4 * N_a(rcut_hard, n + 2) / N_a(rcut_hard, n) * tmp7) +
                      der_pref_f * tmp5;
                  tmp5 = tmp6;
                  tmp6 = tmp7;
                }
              }
            }
            if (do_derivatives) {
              for (n = 0; n < alpha_max; n++)
                exp_coeff_der_temp_d[k_ij * n_temp_der + n] =
                    ampli_tude * exp_coeff_der_temp_d[k_ij * n_temp_der + n] +
                    ampli_tude_der * (exp_coeff_temp1_d[k_ij * n_temp + n] + pref_f * exp_coeff_temp2_d[k_ij * n_temp + n]);
              for (d = i_beg_d[i_sp]; d <= i_end_d[i_sp]; d++) {
                W_exp = 0.0;
                k = 0;
                for (n = i_beg_d[i_sp]; n <= i_end_d[i_sp]; n++) {
                  W_exp += W_d[(n - 1) * n_max + d - 1] * exp_coeff_der_temp_d[k_ij * n_temp_der + k];
                  k += 1;
                }
                exp_coeff_der_d[k_ij * n_max + d - 1] = W_exp;
              }
            }
            for (d = i_beg_d[i_sp]; d <= i_end_d[i_sp]; d++) {
              W_exp = 0.0;
              k = 0;
              for (n = i_beg_d[i_sp]; n <= i_end_d[i_sp]; n++) {
                W_exp += W_d[(n - 1) * n_max + d - 1] *
                         (exp_coeff_temp1_d[k_ij * n_temp + k] + pref_f * exp_coeff_temp2_d[k_ij * n_temp + k]);
                k += 1;
              }
              exp_coeff_d[k_ij * n_max + d - 1] = ampli_tude * W_exp;
            }
          }
        }
      }
    }
  }
}

extern "C" void gpu_radial_poly3(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d, int n_sites,
                                 int* n_neigh_d, int n_max, int n_temp, bool do_derivatives, double* exp_coeff_d,
                                 double* exp_coeff_der_d, double* rcut_soft_d, double* atom_sigma_d, double* exp_coeff_temp1_d,
                                 double* exp_coeff_temp2_d, double* exp_coeff_der_temp_d, int* i_beg, int* i_end,
                                 double* atom_sigma_scaling_d, int mode, int radial_enhancement, double* amplitude_scaling_d,
                                 int* alpha_max_d, double* nf_d, int n_temp_der, double* W_d, bool* do_central_d,
                                 double* central_weight_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_get_radial_poly3<<<nblocks, nthreads, 0, stream[0]>>>(
      n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, n_temp, do_derivatives, exp_coeff_d,
      exp_coeff_der_d, rcut_soft_d, atom_sigma_d, exp_coeff_temp1_d, exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end,
      atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, nf_d, n_temp_der, W_d, do_central_d,
      central_weight_d);
}
