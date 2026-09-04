// The soap_turbo "poly3operator" radial expansion on the device.
//
// The atomic density is a piecewise-cubic bump of half-width
//   width_j = 2*sqrt(2*ln2) * sigma_j
// centred on the neighbour, cut against the soft cutoff by a second cubic
// filter, and projected onto the poly3 radial basis analytically. Every
// integral is therefore a polynomial one, and closes as
//
//   c_alpha = sum_q  s_q  (1 - l_q)^(alpha+3)  sum_k A(alpha,k) M_q(k)
//
// over the three integration limits l_1 <= l_2 <= l_3 of a region, with A the
// atom-independent basis coefficients (get_constant_poly_coeff on the host)
// and M_q the density's polynomial coefficients evaluated at l_q. The soft
// region needs k = 1..4, the buffer region k = 1..7 because the filter raises
// the degree from 3 to 6.
//
// The CPU reference (get_radial_expansion_coefficients_poly3operator in
// soap_turbo_radial.f90) writes this as a per-site sweep over dynamically
// allocated (n_neigh, alpha_max, 3) arrays, which is how it vectorises. None
// of that structure is inherent: every pair is independent, so here it is one
// thread per (i,j) pair, with the per-alpha accumulators in shared memory and
// A staged there once per block. No global scratch, no allocation.
//
// Everything is in units of rcut_hard (rcut_hard = 1), as on the CPU. The
// global sqrt(rcut_hard) and the per-species global_scaling are applied
// afterwards by cuda_global_scaling, exactly as for poly3 and poly3gauss.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb_op 64
#define mode_polynomial 1

// (p-1)! / (p-c)!, the entries of the monomial-derivative matrix R* of degree
// 6 stripped of their powers of r. Row p, column c, both 1-based; the matrix
// is upper triangular. This is M_radial_monomial() in the Fortran.
__device__ const double c_fact_ratio[49] = {
    1.0, 0.0, 0.0,  0.0,   0.0,   0.0,   0.0,  //
    1.0, 1.0, 0.0,  0.0,   0.0,   0.0,   0.0,  //
    1.0, 2.0, 2.0,  0.0,   0.0,   0.0,   0.0,  //
    1.0, 3.0, 6.0,  6.0,   0.0,   0.0,   0.0,  //
    1.0, 4.0, 12.0, 24.0,  24.0,  0.0,   0.0,  //
    1.0, 5.0, 20.0, 60.0,  120.0, 120.0, 0.0,  //
    1.0, 6.0, 30.0, 120.0, 360.0, 720.0, 720.0 //
};

// g_aux(), left and right pieces: the cubic bump and its first three
// derivatives with respect to r, at x = (r - r0)/width. The Fortran array
// form multiplies by 1/width, so do the same -- the scalar form divides, and
// the two differ in the last bit.
__device__ __forceinline__ void op_g_left(double x, double iw, double* g) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double iw2 = iw * iw;
  g[0] = 1.0 - 3.0 * x2 - 2.0 * x3;
  g[1] = -6.0 * (x2 + x) * iw;
  g[2] = -3.0 * (2.0 * x + 1.0) * iw2;
  g[3] = -2.0 * iw2 * iw;
}

__device__ __forceinline__ void op_g_right(double x, double iw, double* g) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double iw2 = iw * iw;
  g[0] = 1.0 - 3.0 * x2 + 2.0 * x3;
  g[1] = 6.0 * (x2 - x) * iw;
  g[2] = 3.0 * (2.0 * x - 1.0) * iw2;
  g[3] = 2.0 * iw2 * iw;
}

// d/dr0 of the above, with width = width0 + width_scaling * r0.
__device__ __forceinline__ void op_g_der_left(double x, double iw, double ws, double* g) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double iw2 = iw * iw;
  const double iw3 = iw2 * iw;
  g[0] = 6.0 * (x + (ws + 1.0) * x2 + ws * x3) * iw;
  g[1] = 6.0 * (1.0 + 2.0 * (ws + 1.0) * x + 3.0 * ws * x2) * iw2;
  g[2] = 6.0 * (ws + 1.0 + 3.0 * ws * x) * iw3;
  g[3] = 6.0 * ws * iw3 * iw;
}

__device__ __forceinline__ void op_g_der_right(double x, double iw, double ws, double* g) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double iw2 = iw * iw;
  const double iw3 = iw2 * iw;
  g[0] = 6.0 * (x + (ws - 1.0) * x2 - ws * x3) * iw;
  g[1] = 6.0 * (1.0 + 2.0 * (ws - 1.0) * x - 3.0 * ws * x2) * iw2;
  g[2] = 6.0 * (ws - 1.0 - 3.0 * ws * x) * iw3;
  g[3] = -6.0 * ws * iw3 * iw;
}

// One soft-region limit: accumulate sign * (1-l)^(alpha+3) * sum_k A(alpha,k)
// * (1-l)^(k-1) * g(k) * vect(k) into the alpha accumulators. vect comes from
// integrating the cubic against the basis by parts.
__device__ __forceinline__ void op_accumulate_soft(const double* g, double sign, double p, const double* sA, int ib, int alpha_max,
                                                   double* acc, int stride) {
  const double v0 = -g[0];
  const double v1 = -p * g[1];
  const double v2 = -2.0 * p * p * g[2];
  const double v3 = -6.0 * p * p * p * g[3];
  double pa = p * p;
  pa = pa * pa; // (1-l)^4, the alpha = 1 power
  for (int a = 0; a < alpha_max; ++a) {
    const double* Ar = sA + (ib - 1 + a) * 7;
    acc[a * stride] += sign * (Ar[0] * v0 + Ar[1] * v1 + Ar[2] * v2 + Ar[3] * v3) * pa;
    pa *= p;
  }
}

// One buffer-region limit. B holds the degree-6 coefficients of the filtered
// density (get_constant_poly_filter_coeff_array); contracting it with R*(l)
// gives the density and its derivatives there.
__device__ __forceinline__ void op_accumulate_buffer(const double* B, double r, double sign, double p, const double* sA, int ib,
                                                     int alpha_max, double* acc, int stride) {
  double M[7];
  double pc = 1.0;
  for (int c = 0; c < 7; ++c) {
    double s = 0.0;
    double rp = 1.0;
    for (int q = c; q < 7; ++q) {
      s += B[q] * rp * c_fact_ratio[q * 7 + c];
      rp *= r;
    }
    M[c] = -s * pc;
    pc *= p;
  }
  double pa = p * p;
  pa = pa * pa;
  for (int a = 0; a < alpha_max; ++a) {
    const double* Ar = sA + (ib - 1 + a) * 7;
    double t = 0.0;
    for (int c = 0; c < 7; ++c) {
      t += Ar[c] * M[c];
    }
    acc[a * stride] += sign * t * pa;
    pa *= p;
  }
}

// Polynomial product of the two cubics: the Toeplitz matmul the Fortran builds
// explicitly, which is a length-7 convolution and nothing more.
__device__ __forceinline__ void op_conv4(const double* a, const double* b, double* out) {
  for (int r = 0; r < 7; ++r) {
    double s = 0.0;
    for (int i = 0; i < 4; ++i) {
      const int j = r - i;
      if (j >= 0 && j < 4) {
        s += a[j] * b[i];
      }
    }
    out[r] = s;
  }
}

__global__ void kernel_get_radial_poly3operator(int n_atom_pairs, int n_species, const bool* mask_d, const double* rjs_d,
                                                const double* rcut_hard_d, const double* rcut_soft_d, const double* atom_sigma_d,
                                                const double* atom_sigma_scaling_d, const double* amplitude_scaling_d,
                                                const double* central_weight_d, const int* alpha_max_d, const int* i_beg_d,
                                                const int* i_end_d, const int* k2_start_d, const int* k2_i_site_d,
                                                const double* A_d, const double* W_d, int n_max, int max_alpha, int mode,
                                                int radial_enhancement, bool do_derivatives, double* exp_coeff_d,
                                                double* exp_coeff_der_d) {
  extern __shared__ double smem[];
  double* sA = smem;
  double* sC = smem + n_max * 7;
  double* sD = sC + blockDim.x * max_alpha;

  for (int t = threadIdx.x; t < n_max * 7; t += blockDim.x) {
    sA[t] = A_d[t];
  }
  __syncthreads();

  const int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (k_ij >= n_atom_pairs) {
    return;
  }

  const double pi = acos(-1.0);
  const double fwhm = 2.0 * sqrt(2.0 * log(2.0));
  const double rjs_in = rjs_d[k_ij];
  const bool is_central = (k_ij == k2_start_d[k2_i_site_d[k_ij] - 1]);
  const int stride = blockDim.x;
  double* C = sC + threadIdx.x;
  double* D = sD + threadIdx.x;

  for (int i_sp = 0; i_sp < n_species; ++i_sp) {
    if (!mask_d[k_ij + i_sp * n_atom_pairs]) {
      continue;
    }
    const double rcut_hard_in = rcut_hard_d[i_sp];
    const double rcut_soft_in = rcut_soft_d[i_sp];
    const double atom_sigma_in = atom_sigma_d[i_sp];
    const double sig_scaling = atom_sigma_scaling_d[i_sp];

    // Region membership, decided in the same unnormalised units as the CPU so
    // that a pair sitting exactly on a boundary lands on the same side. Both
    // boundaries are harmless anyway: there the three limits coincide and the
    // region contributes zero.
    const double width_in = fwhm * (atom_sigma_in + sig_scaling * rjs_in);
    const bool in_soft = (rjs_in - width_in < rcut_soft_in);
    const bool in_buffer = (rcut_soft_in < rcut_hard_in) && (rjs_in + width_in > rcut_soft_in);
    if (!in_soft && !in_buffer) {
      continue;
    }

    const int alpha_max = alpha_max_d[i_sp];
    const int ib = i_beg_d[i_sp];
    const double rcut_soft = rcut_soft_in / rcut_hard_in;
    const double rj = rjs_in / rcut_hard_in;
    const double ass = atom_sigma_in / rcut_hard_in + sig_scaling * rj;
    const double s2 = ass * ass;
    const double width = fwhm * ass;
    const double iwidth = 1.0 / width;
    const double width_scaling = fwhm * sig_scaling;

    double amp = 0.0;
    double amp_der = 0.0;
    if (mode == mode_polynomial) {
      const double as = amplitude_scaling_d[i_sp];
      if (as == 0.0) {
        amp = 1.0 / ass;
        amp_der = -sig_scaling / s2;
      } else {
        const double t1 = 1.0 + 2.0 * rj * rj * rj - 3.0 * rj * rj;
        if (as == 1.0) {
          amp = 1.0 / ass * t1;
          amp_der = 6.0 / ass * (rj * rj - rj) - sig_scaling / ass * amp;
        } else {
          amp = 1.0 / ass * pow(t1, as);
          amp_der = 6.0 * as / ass * (rj * rj - rj) * pow(t1, as - 1.0) - sig_scaling / ass * amp;
        }
      }
    }
    // The central atom is scaled, never skipped: the CPU call site passes
    // do_central = .true. unconditionally for this basis and lets
    // central_weight do the work.
    if (is_central) {
      amp *= central_weight_d[i_sp];
      amp_der *= central_weight_d[i_sp];
    }
    if (radial_enhancement == 1) {
      const double t3 = rj + sqrt(2.0 / pi) * ass;
      amp_der = amp * (1.0 + sqrt(2.0 / pi) * sig_scaling) + amp_der * t3;
      amp = amp * t3;
    } else if (radial_enhancement == 2) {
      const double t4 = sqrt(8.0 / pi) * ass;
      const double t5 = rj * rj + s2 + t4 * rj;
      amp_der = amp * (2.0 * rj + 2.0 * ass * sig_scaling + t4 + sqrt(8.0 / pi) * rj * sig_scaling) + amp_der * t5;
      amp = amp * t5;
    }

    for (int a = 0; a < alpha_max; ++a) {
      C[a * stride] = 0.0;
      D[a * stride] = 0.0;
    }

    if (in_soft) {
      const double l1 = fmax(0.0, rj - width);
      const double l2 = fmin(rj, rcut_soft);
      const double l3 = fmin(rcut_soft, rj + width);
      const double p1 = 1.0 - l1;
      const double p2 = 1.0 - l2;
      const double p3 = 1.0 - l3;
      double g[4];
      op_g_left((l1 - rj) * iwidth, iwidth, g);
      op_accumulate_soft(g, -1.0, p1, sA, ib, alpha_max, C, stride);
      op_g_left((l2 - rj) * iwidth, iwidth, g);
      op_accumulate_soft(g, 1.0, p2, sA, ib, alpha_max, C, stride);
      op_g_right((l2 - rj) * iwidth, iwidth, g);
      op_accumulate_soft(g, -1.0, p2, sA, ib, alpha_max, C, stride);
      op_g_right((l3 - rj) * iwidth, iwidth, g);
      op_accumulate_soft(g, 1.0, p3, sA, ib, alpha_max, C, stride);
      if (do_derivatives) {
        op_g_der_left((l1 - rj) * iwidth, iwidth, width_scaling, g);
        op_accumulate_soft(g, -1.0, p1, sA, ib, alpha_max, D, stride);
        op_g_der_left((l2 - rj) * iwidth, iwidth, width_scaling, g);
        op_accumulate_soft(g, 1.0, p2, sA, ib, alpha_max, D, stride);
        op_g_der_right((l2 - rj) * iwidth, iwidth, width_scaling, g);
        op_accumulate_soft(g, -1.0, p2, sA, ib, alpha_max, D, stride);
        op_g_der_right((l3 - rj) * iwidth, iwidth, width_scaling, g);
        op_accumulate_soft(g, 1.0, p3, sA, ib, alpha_max, D, stride);
      }
    }

    if (in_buffer) {
      const double l1 = fmax(rcut_soft, rj - width);
      const double l2 = fmax(rj, rcut_soft);
      const double l3 = fmin(1.0, rj + width);
      const double p1 = 1.0 - l1;
      const double p2 = 1.0 - l2;
      const double p3 = 1.0 - l3;

      // The filter is the right piece of a cubic centred on rcut_soft with
      // width 2*sqrt(2*ln2)*(rcut_hard - rcut_soft), evaluated at r = 0.
      // g_aux() divides by the width here, so match that.
      const double filter_width = fwhm * (1.0 - rcut_soft);
      const double xf = -rcut_soft / filter_width;
      double cf[4];
      cf[0] = 1.0 - 3.0 * xf * xf + 2.0 * xf * xf * xf;
      cf[1] = 6.0 * (xf * xf - xf) / filter_width;
      cf[2] = 3.0 * (2.0 * xf - 1.0) / (filter_width * filter_width);
      cf[3] = 2.0 / (filter_width * filter_width * filter_width);

      const double xj = -rj * iwidth;
      double cp[4];
      double B[7];
      op_g_left(xj, iwidth, cp);
      op_conv4(cp, cf, B);
      op_accumulate_buffer(B, l1, -1.0, p1, sA, ib, alpha_max, C, stride);
      op_accumulate_buffer(B, l2, 1.0, p2, sA, ib, alpha_max, C, stride);
      op_g_right(xj, iwidth, cp);
      op_conv4(cp, cf, B);
      op_accumulate_buffer(B, l2, -1.0, p2, sA, ib, alpha_max, C, stride);
      op_accumulate_buffer(B, l3, 1.0, p3, sA, ib, alpha_max, C, stride);
      if (do_derivatives) {
        op_g_der_left(xj, iwidth, width_scaling, cp);
        op_conv4(cp, cf, B);
        op_accumulate_buffer(B, l1, -1.0, p1, sA, ib, alpha_max, D, stride);
        op_accumulate_buffer(B, l2, 1.0, p2, sA, ib, alpha_max, D, stride);
        op_g_der_right(xj, iwidth, width_scaling, cp);
        op_conv4(cp, cf, B);
        op_accumulate_buffer(B, l2, -1.0, p2, sA, ib, alpha_max, D, stride);
        op_accumulate_buffer(B, l3, 1.0, p3, sA, ib, alpha_max, D, stride);
      }
    }

    // Amplitude first, then the change of basis, in that order: it is what the
    // CPU does, and W mixes the alphas so the two orders do not round alike.
    for (int a = 0; a < alpha_max; ++a) {
      const double c = C[a * stride];
      if (do_derivatives) {
        D[a * stride] = amp * D[a * stride] + amp_der * c;
      }
      C[a * stride] = amp * c;
    }

    const int ie = i_end_d[i_sp];
    for (int d = ib; d <= ie; ++d) {
      double wc = 0.0;
      double wd = 0.0;
      for (int a = 0; a < alpha_max; ++a) {
        const double wv = W_d[(ib - 1 + a) * n_max + d - 1];
        wc += wv * C[a * stride];
        wd += wv * D[a * stride];
      }
      exp_coeff_d[k_ij * n_max + d - 1] = wc;
      if (do_derivatives) {
        exp_coeff_der_d[k_ij * n_max + d - 1] = wd;
      }
    }
  }
}

extern "C" void gpu_radial_poly3operator(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d,
                                         double* rcut_soft_d, double* atom_sigma_d, double* atom_sigma_scaling_d,
                                         double* amplitude_scaling_d, double* central_weight_d, int* alpha_max_d, int* i_beg_d,
                                         int* i_end_d, int* k2_start_d, int* k2_i_site_d, double* A_d, double* W_d, int n_max,
                                         int max_alpha, int mode, int radial_enhancement, bool do_derivatives, double* exp_coeff_d,
                                         double* exp_coeff_der_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb_op) / tpb_op, 1, 1);
  dim3 nthreads = dim3(tpb_op, 1, 1);
  size_t shmem = (static_cast<size_t>(n_max) * 7 + 2 * static_cast<size_t>(tpb_op) * max_alpha) * sizeof(double);

  kernel_get_radial_poly3operator<<<nblocks, nthreads, shmem, stream[0]>>>(
      n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, rcut_soft_d, atom_sigma_d, atom_sigma_scaling_d, amplitude_scaling_d,
      central_weight_d, alpha_max_d, i_beg_d, i_end_d, k2_start_d, k2_i_site_d, A_d, W_d, n_max, max_alpha, mode,
      radial_enhancement, do_derivatives, exp_coeff_d, exp_coeff_der_d);
}
