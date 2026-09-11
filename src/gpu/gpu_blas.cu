// hipBLAS/cuBLAS handle lifetime and the dense products the GAP prediction
// needs. Nothing here is specific to a descriptor.
#include "gpu_common.h"
#include "gpu_blas.h"

extern "C" void create_cublas_handle(hipblasHandle_t* handle, hipStream_t* stream) {
  hipblasCreate(handle);
  hipStreamCreate(stream);
  hipblasSetStream(*handle, *stream);
  //    hipsolverCreate(handle);
  //    hipsolverSetStream(*handle,*stream);
  /*printf("\n cublas handle created \n");
    exit(0);*/

  return;
}

extern "C" void destroy_cublas_handle(hipblasHandle_t* handle, hipStream_t* stream) {
  // Destroy the handle
  hipblasDestroy(*handle);
  hipStreamDestroy(*stream);
  //printf("\n cublas handle destroyed. \n The End? \n");
  return;
}

extern "C" void gpu_blas_mmul_t_n(hipblasHandle_t handle, const double* Qs_d, const double* soap_d, double* kernels_d,
                                  const int n_sparse, const int n_soap, const int n_sites)
//                                                           const double *A,     const double *B,         double *C,       const int nAx,
// const int nAy,      const int nBy,double *b, double zeta, int N)
{
  // (hipblasHandle_t handle, const double *Qs_d, const double *soap_d, double *kernels_d, const int n_sparse, const int n_soap, const int n_sites,double *b, double zeta, int N)
  const double alf = 1;
  const double bet = 0;

  // soap(n_soap,n_sites)
  // Qs(1:n_soap, 1:n_sparse)
  // kernels(1:n_sites, 1:n_sparse)
  // call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)

  // Do the actual multiplication
  hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d,
               n_sites);

  return;
}

extern "C" void gpu_blas_mvmul_n(hipblasHandle_t handle, double* kernels_copy_d, const double* alphas_d, double* energies_d,
                                 const int n_sites, const int n_sparse) {
  const double alf = 1;
  const double bet = 0;
  const double* alpha = &alf;
  const double* beta = &bet;

  // Do the actual multiplication
  hipblasDgemv(handle, HIPBLAS_OP_N, n_sites, n_sparse, alpha, kernels_copy_d, n_sites, alphas_d, 1, beta, energies_d, 1);
  return;
}


extern "C" void gpu_blas_mmul_n_t(hipblasHandle_t handle, const double* kernels_der_d, const double* Qs_copy_d, double* Qss_d,
                                  const int n_sparse, const int n_soap, const int n_sites, double cdelta) {
  const double alf = cdelta;
  const double bet = 0;
  const double* alpha = &alf;
  const double* beta = &bet;
  // soap(n_soap,n_sites)
  // Qs(1:n_soap, 1:n_sparse)
  // kernels(1:n_sites, 1:n_sparse)
  // call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)
  // hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, alpha, soap_d, n_soap, Qs_d, n_soap, beta, kernels_d, n_sites);

  // allocate( kernels_der(1:n_sites, 1:n_sparse)
  // allocate( Qs_copy(1:n_soap, 1:n_sparse) ))
  // allocate( Qss(1:n_sites, 1:n_soap) )
  // call dgemm("n", "t", n_sites, n_soap, n_sparse, cdelta, kernels_der, n_sites, Qs_copy, n_soap, 0.d0, Qss, n_sites)
  hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, n_sites, n_soap, n_sparse, alpha, kernels_der_d, n_sites, Qs_copy_d, n_soap,
               beta, Qss_d, n_sites);
}


extern "C" void gpu_dgemm_n_n(int m, int n, int k, double alpha, double* A, int lda, double* B, int ldb, double beta, double* C,
                              int ldc, hipblasHandle_t handle) {
  const double* alpha_a = &alpha;
  const double* beta_a = &beta;

  hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, m, n, k, alpha_a, A, lda, B, ldb, beta_a, C, ldc);
}
