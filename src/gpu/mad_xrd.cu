// MAD: the X-ray/neutron structure factor.
//
// Gk and its derivative with respect to atomic positions, the Hadamard product
// with the scattering factors, the dgemv that turns that into per-pair forces,
// and the collection of those into forces and the virial.
//
// Every kernel here is a flat map, written once through tg_parallel_for and
// dispatched to CUDA or Kokkos by src/gpu/gpu_backend.h.
//
// Most index n_k * n_samples, a product of two ints that need not be an int,
// so they take the long-indexed launcher. The in-kernel `if (k < n_k && ...)`
// guards the hand-written launches carried are gone with them: those existed
// because a <<<>>> grid is rounded up to whole blocks and over-runs the extent.
// tg_parallel_for_long passes the exact extent, so the guard was dead weight.
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gpu_scatter.h"
#include "mad_gpu.h"

#define tpb 512

// Gk(1:n_samples, n_k) = -2 * c_factor * pair_distribution_der(1:n_samples, k)
extern "C" void gpu_set_Gk(int nk, int n_samples, int* k_index_d, double* Gk_d, double* pair_distribution_partial_der_d,
                           double c_factor, hipStream_t* stream) {
  tg_parallel_for_long(
      "turbogap_set_Gk", (long long) nk * n_samples, stream,
      TG_LAMBDA(const long long tid) { Gk_d[tid] = -2.0 * c_factor * pair_distribution_partial_der_d[tid]; }, BLOCK_SIZE);
}

extern "C" void gpu_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d, hipStream_t* stream) {
  tg_parallel_for_long(
      "turbogap_get_Gka", (long long) n_k * n_samples, stream,
      TG_LAMBDA(const long long tid) {
        const long long k = tid / n_samples;
        Gka_d[tid] = Gk_d[tid] * xyz_k_d[(i - 1) + k * 3];
      },
      BLOCK_SIZE);
}

extern "C" void gpu_get_Gka_inplace(int i, int n_k, int n_samples, double* Gk_d, double* xyz_k_d, hipStream_t* stream) {
  tg_parallel_for_long(
      "turbogap_get_Gka_inplace", (long long) n_k * n_samples, stream,
      TG_LAMBDA(const long long tid) {
        const long long k = tid / n_samples;
        Gk_d[tid] *= xyz_k_d[(i - 1) + k * 3];
      },
      BLOCK_SIZE);
}

extern "C" void gpu_hadamard_vec_mat_product(int n_samples_sf, int n_k, double* all_scattering_factors_d, double* dermat_d,
                                             hipStream_t* stream) {
  tg_parallel_for_long(
      "turbogap_hadamard_vec_mat", (long long) n_k * n_samples_sf, stream,
      TG_LAMBDA(const long long tid) {
        const long long l = tid % n_samples_sf;
        dermat_d[tid] *= all_scattering_factors_d[l];
      },
      BLOCK_SIZE);
}

extern "C" void gpu_get_fi_dgemv(const int i, const int n_samples_sf, const int n_k, double* dermat_d, double* prefactor_d,
                                 double* fi_d, hipblasHandle_t handle, hipStream_t* stream) {
  const double alf = 1;
  const double bet = 0;
  const double* alpha = &alf;
  const double* beta = &bet;
  double* ptr = fi_d + (i - 1) * n_k;

  // call dgemv("T", n_samples_sf, n_k, 1.d0, dermat, n_samples_sf,
  //            prefactor, 1, 0.d0, fi(:,i), 1)
  hipblasDgemv(handle, HIPBLAS_OP_T, n_samples_sf, n_k, alpha, dermat_d, n_samples_sf, prefactor_d, 1, beta, ptr, 1);
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
