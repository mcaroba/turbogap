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
// then gpu_pair_scatter_reduce sums them onto atoms in neighbour-list order.
// forces_d (3 x n_sites) and virial_d (9) are accumulated into.
extern "C" void gpu_xrd_window_forces(int n_k, int n_sites, int n_samples, double r_min, double dr, double kde_sigma,
                                      double kde_cut, double* x_d, double* w_d, double* r_d, double* xyz_d, int* j2_d, int* chan_d,
                                      double* forces_d, double* virial_d, hipStream_t* stream) {
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
        double a = 0.0;
        for (int l = lo; l <= hi; l++) {
          const double g = (x_d[l] - r) / kde_sigma;
          a += w_c[l] * exp(-0.5 * g * g) * g;
        }
        a /= kde_sigma * r;
        pair_force[3 * (size_t) p] = a * xyz_d[3 * (size_t) p];
        pair_force[3 * (size_t) p + 1] = a * xyz_d[3 * (size_t) p + 1];
        pair_force[3 * (size_t) p + 2] = a * xyz_d[3 * (size_t) p + 2];
      },
      tpb);

  // 0.25: see gpu_exp_force_virial_collection.
  gpu_pair_scatter_reduce(n_k, n_sites, j2_d, pair_force, xyz_d, forces_d, virial_d, 0.25, stream);

  gpuErrchk(hipFreeAsync(pair_force, stream[0]));
}
