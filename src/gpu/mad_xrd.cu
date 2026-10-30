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
