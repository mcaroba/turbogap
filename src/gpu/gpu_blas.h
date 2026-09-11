// hipBLAS handle lifetime and the dense products used by GAP prediction.

#ifndef TURBOGAP_GPU_BLAS_H
#define TURBOGAP_GPU_BLAS_H

#include "gpu_common.h"

extern "C" {

// ---- gpu_blas.cu
void create_cublas_handle(hipblasHandle_t* handle, hipStream_t* stream);
void destroy_cublas_handle(hipblasHandle_t* handle, hipStream_t* stream);
void gpu_blas_mmul_t_n(hipblasHandle_t handle, const double* Qs_d, const double* soap_d, double* kernels_d, const int n_sparse,
                       const int n_soap, const int n_sites);
void gpu_blas_mvmul_n(hipblasHandle_t handle, double* kernels_copy_d, const double* alphas_d, double* energies_d, const int n_sites,
                      const int n_sparse);
void gpu_blas_mmul_n_t(hipblasHandle_t handle, const double* kernels_der_d, const double* Qs_copy_d, double* Qss_d,
                       const int n_sparse, const int n_soap, const int n_sites, double cdelta);
void gpu_dgemm_n_n(int m, int n, int k, double alpha, double* A, int lda, double* B, int ldb, double beta, double* C, int ldc,
                   hipblasHandle_t handle);

} // extern "C"

#endif // TURBOGAP_GPU_BLAS_H
