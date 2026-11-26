// MAD: the X-ray/neutron structure factor.
//
// The per-pair forces from the pdf derivative, and their collection into
// forces and the virial.
//
// Every kernel here is a flat map, written once through tg_parallel_for and
// dispatched to CUDA or Kokkos by src/gpu/gpu_backend.h.
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gpu_scatter.h"
#include "mad_gpu.h"

#define tpb 512

// fi(k, c) = -2 c_factor * sum_s w_s dg_s/dr_k * xyz_k(c): pair k's force before
// the energy scale, with w = S^T (residual x scattering factor) formed on the
// host. pdf_der is (n_samples, nk), fi is (nk, 3).
extern "C" void gpu_xrd_pair_forces(int nk, int n_samples, double* w_d, double* pdf_der_d, double c_factor, double* xyz_k_d,
                                    double* fi_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_xrd_pair_forces", nk, stream,
      TG_LAMBDA(const int k) {
        double a = 0.0;
        for (int s = 0; s < n_samples; s++)
          a += w_d[s] * pdf_der_d[s + (long long) n_samples * k];
        a *= -2.0 * c_factor;
        for (int c = 0; c < 3; c++)
          fi_d[k + (long long) nk * c] = a * xyz_k_d[c + 3 * (long long) k];
      },
      tpb);
}

// w(s) = sum_q S(q, s) residual(q) f(q), f = 1 when f_d is null: the weights
// gpu_xrd_pair_forces contracts the pdf derivative with, where S, the residual
// and the scattering factors are already on the device. S is (n_samples_sf,
// n_samples).
extern "C" void gpu_xrd_weights(int n_samples_sf, int n_samples, double* S_d, double* residual_d, double* f_d, double* w_d,
                                hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_xrd_weights", n_samples, stream,
      TG_LAMBDA(const int s) {
        double a = 0.0;
        for (int q = 0; q < n_samples_sf; q++)
          a += S_d[q + (long long) n_samples_sf * s] * residual_d[q] * (f_d ? f_d[q] : 1.0);
        w_d[s] = a;
      },
      tpb);
}

extern "C" void gpu_exp_force_virial_collection(int n_k, int n_sites, double3* forces0, double energy_scale, double* fi,
                                                int* j2_list, double* virial, double3* xyz, hipStream_t* stream) {
  if (n_k <= 0)
    return;

  double* pair_force_d;
  gpuErrchk(hipMallocAsync(&pair_force_d, (size_t) 3 * n_k * sizeof(double), stream[0]));

  // The pair's own force, staged at its own index. It used to be added straight
  // into forces0[j2] and the virial into nine addresses, both with atomicAdd,
  // which made the result depend on the order the threads arrived; the sum is
  // gpu_pair_scatter_reduce's now. See gpu_scatter.h.
  tg_parallel_for(
      "turbogap_exp_force_collection", n_k, stream,
      TG_LAMBDA(const int tid) {
        pair_force_d[3 * tid] = energy_scale * fi[tid];
        pair_force_d[3 * tid + 1] = energy_scale * fi[tid + n_k];
        pair_force_d[3 * tid + 2] = energy_scale * fi[tid + 2 * n_k];
      },
      tpb);

  // 0.25 = 0.5 to symmetrise, times 0.5 because k_list holds ordered pairs and
  // the reversed one reproduces the term rather than adding a new one.
  gpu_pair_scatter_reduce(n_k, n_sites, j2_list, pair_force_d, (const double*) xyz, (double*) forces0, virial, 0.25, stream);

  gpuErrchk(hipFreeAsync(pair_force_d, stream[0]));
}

// ------------------------------------------------------------------------ //
// --- Windowed forces: every species pair at once, no stored derivative --- //
// ------------------------------------------------------------------------ //

// W(l, c) = fac / dV(l) * sum_q S(q, l) residual(q) f(q) for channel c, f = 1
// when f_d is null. fac carries everything constant over the channel: the
// energy scale, -2 c_factor, the kde normalisation and the partial's scale.
extern "C" void gpu_xrd_window_weights(int n_samples_sf, int n_samples, int channel, double* S_d, double* residual_d, double* f_d,
                                       double* dV_d, double fac, double* w_d, hipStream_t* stream) {
  double* w_c = w_d + (size_t) channel * n_samples;
  tg_parallel_for(
      "turbogap_xrd_window_weights", n_samples, stream,
      TG_LAMBDA(const int s) {
        double a = 0.0;
        for (int q = 0; q < n_samples_sf; q++)
          a += S_d[q + (size_t) n_samples_sf * s] * residual_d[q] * (f_d ? f_d[q] : 1.0);
        w_c[s] = fac * a / dV_d[s];
      },
      tpb);
}

// The forces and the pair half of the virial from the compacted pairs of
// gpu_pdf_compact_pairs. Each pair rebuilds its kernel over the samples the
// pdf used and contracts it with its channel's W:
//
//   F_p = xyz_p / (sigma r_p) sum_l W(l, c_p) exp(-g^2 / 2) g,  g = (x_l - r_p) / sigma
//
// The kernel is built from the sample c nearest r outwards: on the pdf's evenly
// spaced grid the ratio of two neighbouring values changes by exp(-h^2) from
// one sample to the next (h = dr / sigma), so one exp at c and two for the
// ratios start it and each further sample is two multiplications. Above
// h = 2 (a kernel narrower than half a sample) each sample takes its exp.
//
// F_p acts on the pair's neighbour. The list is full and the reversed pair has
// the same r and channel and the opposite xyz, so its force is -F_p: a site's
// force is minus the sum of its own pairs' F, which are one run per channel
// (seg_d, from gpu_pdf_compact_pairs) summed in order -- no scatter. Site i of
// the batch is atom i_off + i. The same loop sums its pairs' virial, and the
// sites' are then summed in site order. forces_d (3 x n_sites) and virial_d (9)
// are accumulated into.
extern "C" void gpu_xrd_window_forces(int n_k, int n_sites, int n_samples, double r_min, double dr, double kde_sigma,
                                      double kde_cut, double* x_d, double* w_d, double* r_d, double* xyz_d, int* chan_d,
                                      int n_channels, int n_sites_b, int i_off, int* seg_d, double* forces_d, double* virial_d,
                                      hipStream_t* stream) {
  if (n_k <= 0)
    return;

  double* pair_force;
  gpuErrchk(hipMallocAsync(&pair_force, (size_t) 3 * n_k * sizeof(double), stream[0]));

  tg_parallel_for(
      "turbogap_xrd_window_pair_forces", n_k, stream,
      TG_LAMBDA(const int p) {
        const double r = r_d[p];
        const double* w_c = w_d + (size_t) chan_d[p] * n_samples;
        const int lo = pdf_window_lo(r, r_min, dr, kde_sigma, kde_cut);
        const int hi = pdf_window_hi(r, r_min, dr, kde_sigma, kde_cut, n_samples);
        const double h = dr / kde_sigma;
        double a = 0.0;
        if (lo <= hi && h <= 2.0) {
          int c = lo + (int) rint((r - x_d[lo]) / dr);
          c = (c < lo) ? lo : ((c > hi) ? hi : c);
          const double u = (x_d[c] - r) / kde_sigma;
          const double q = exp(-h * h);
          const double g0 = exp(-0.5 * u * u);
          a = w_c[c] * g0 * u;
          double g = g0;
          double t = exp(-u * h - 0.5 * h * h);
          for (int l = c + 1; l <= hi; l++) {
            g *= t;
            t *= q;
            a += w_c[l] * g * (u + (l - c) * h);
          }
          g = g0;
          t = exp(u * h - 0.5 * h * h);
          for (int l = c - 1; l >= lo; l--) {
            g *= t;
            t *= q;
            a += w_c[l] * g * (u - (c - l) * h);
          }
        } else {
          for (int l = lo; l <= hi; l++) {
            const double g = (x_d[l] - r) / kde_sigma;
            a += w_c[l] * exp(-0.5 * g * g) * g;
          }
        }
        a /= kde_sigma * r;
        pair_force[3 * (size_t) p] = a * xyz_d[3 * (size_t) p];
        pair_force[3 * (size_t) p + 1] = a * xyz_d[3 * (size_t) p + 1];
        pair_force[3 * (size_t) p + 2] = a * xyz_d[3 * (size_t) p + 2];
      },
      tpb);

  double* site_virial;
  gpuErrchk(hipMallocAsync(&site_virial, (size_t) 9 * (n_sites_b > 0 ? n_sites_b : 1) * sizeof(double), stream[0]));

  // 0.25: see gpu_exp_force_virial_collection.
  const int n_seg = n_sites_b + 1;
  tg_parallel_for(
      "turbogap_xrd_window_site_forces", n_sites_b, stream,
      TG_LAMBDA(const int i) {
        double f[3] = {0.0, 0.0, 0.0};
        double v[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int c = 0; c < n_channels; c++) {
          for (int p = seg_d[c * n_seg + i]; p < seg_d[c * n_seg + i + 1]; p++) {
            const double* fp = pair_force + 3 * (size_t) p;
            const double* xp = xyz_d + 3 * (size_t) p;
            for (int k1 = 0; k1 < 3; k1++) {
              f[k1] += fp[k1];
              for (int k2 = 0; k2 < 3; k2++)
                v[3 * k1 + k2] += 0.25 * (fp[k1] * xp[k2] + fp[k2] * xp[k1]);
            }
          }
        }
        const size_t a = 3 * (size_t) (i_off + i);
        forces_d[a] -= f[0];
        forces_d[a + 1] -= f[1];
        forces_d[a + 2] -= f[2];
        for (int k = 0; k < 9; k++)
          site_virial[9 * (size_t) i + k] = v[k];
      },
      tpb);

  tg_parallel_for(
      "turbogap_xrd_window_virial", 9, stream,
      TG_LAMBDA(const int k) {
        double s = 0.0;
        for (int i = 0; i < n_sites_b; i++)
          s += site_virial[9 * (size_t) i + k];
        virial_d[k] += s;
      },
      9);

  gpuErrchk(hipFreeAsync(site_virial, stream[0]));

  gpuErrchk(hipFreeAsync(pair_force, stream[0]));
}
