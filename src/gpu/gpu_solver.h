/* ---------------------------------------------------------------------------
   hipSOLVER spellings on cuSOLVER, for the CUDA target.

   hop maps the HIP runtime, hipBLAS, hipFFT, hipRAND and hipRTC onto their
   CUDA equivalents, but ships no solver header -- so a file written against
   hipSOLVER does not compile under HOP_TARGET_CUDA. orthonormalization_kernels.cc
   is the only such file here, and its own #include is commented out, which is
   part of why it sits disabled in the Makefile.

   Most of the surface maps by name. Three entry points genuinely differ in
   shape and need a wrapper rather than a #define:

     hipsolverDgesvd_bufferSize takes jobu and jobvt; cuSOLVER's does not.
     hipsolverDgetrf takes a workspace SIZE; cuSOLVER's takes only the buffer.
     hipsolverDgetrs takes a workspace; cuSOLVER's does not.

   The wrappers accept the hipSOLVER argument list and drop what cuSOLVER has
   no parameter for, so call sites stay in one spelling.

   THIS IS NECESSARY BUT NOT SUFFICIENT, and the reason matters.

   orthonormalization_kernels.cc declares its helpers as taking a
   hipblasHandle_t and passes that straight to the hipsolver* calls:

       void sqrt_mat(double* W, const int alpha_max, hipStream_t& s, hipblasHandle_t& handle)
           ... hipsolverDgesvd_bufferSize(handle, ...)

   On ROCm that compiles because hipSOLVER's handle is a void*, so a BLAS
   handle passes the type check. On CUDA cublasHandle_t and cusolverDnHandle_t
   are distinct types and it does not.

   Making them agree here -- typedefing hipsolverHandle_t to void* and casting
   inside the wrappers -- would compile and then fail at runtime, because
   cuSOLVER needs a handle from cusolverDnCreate and would be handed a cuBLAS
   one. A compile error is the better outcome of the two, so this header does
   not paper over it.

   Re-enabling the file on CUDA therefore needs a real change: a cuSOLVER
   handle created and threaded through those five helpers and their Fortran
   call sites, beside the cuBLAS one. This header is the API mapping that work
   will need; it is not the work.

   Compile-tested against CUDA 13.4 headers. Nothing here has been run.
--------------------------------------------------------------------------- */

#ifndef TG_GPU_SOLVER_H
#define TG_GPU_SOLVER_H

#if defined(HOP_TARGET_CUDA) || defined(CUDA)

#include <cusolverDn.h>

typedef cusolverDnHandle_t hipsolverHandle_t;
typedef cusolverStatus_t hipsolverStatus_t;
typedef cublasOperation_t hipsolverOperation_t;
typedef cublasFillMode_t hipsolverFillMode_t;

#define HIPSOLVER_STATUS_SUCCESS CUSOLVER_STATUS_SUCCESS
#define HIPSOLVER_FILL_MODE_UPPER CUBLAS_FILL_MODE_UPPER
#define HIPSOLVER_FILL_MODE_LOWER CUBLAS_FILL_MODE_LOWER
#define HIPSOLVER_OP_N CUBLAS_OP_N
#define HIPSOLVER_OP_T CUBLAS_OP_T

/* Same argument list on both sides. */
#define hipsolverDnDgesvd cusolverDnDgesvd
#define hipsolverDpotrf_bufferSize cusolverDnDpotrf_bufferSize
#define hipsolverDpotrf cusolverDnDpotrf
#define hipsolverDpotri_bufferSize cusolverDnDpotri_bufferSize
#define hipsolverDpotri cusolverDnDpotri
#define hipsolverDgetrf_bufferSize cusolverDnDgetrf_bufferSize
#define hipsolverCreate cusolverDnCreate
#define hipsolverDestroy cusolverDnDestroy
#define hipsolverSetStream cusolverDnSetStream

/* cuSOLVER sizes the gesvd workspace from the extents alone. */
inline cusolverStatus_t hipsolverDgesvd_bufferSize(cusolverDnHandle_t handle, signed char /*jobu*/, signed char /*jobvt*/, int m,
                                                   int n, int* lwork) {
  return cusolverDnDgesvd_bufferSize(handle, m, n, lwork);
}

/* cuSOLVER's getrf takes the workspace but not its size. */
inline cusolverStatus_t hipsolverDgetrf(cusolverDnHandle_t handle, int m, int n, double* A, int lda, double* work, int /*lwork*/,
                                        int* devIpiv, int* devInfo) {
  return cusolverDnDgetrf(handle, m, n, A, lda, work, devIpiv, devInfo);
}

/* cuSOLVER's getrs needs no workspace at all. */
inline cusolverStatus_t hipsolverDgetrs(cusolverDnHandle_t handle, cublasOperation_t trans, int n, int nrhs, double* A, int lda,
                                        int* devIpiv, double* B, int ldb, double* /*work*/, int /*lwork*/, int* devInfo) {
  return cusolverDnDgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb, devInfo);
}

/* No cuSOLVER counterpart: getrs allocates nothing, so the answer is zero. */
inline cusolverStatus_t hipsolverDgetrs_bufferSize(cusolverDnHandle_t /*handle*/, cublasOperation_t /*trans*/, int /*n*/,
                                                   int /*nrhs*/, double* /*A*/, int /*lda*/, int* /*devIpiv*/, double* /*B*/,
                                                   int /*ldb*/, int* lwork) {
  *lwork = 0;
  return CUSOLVER_STATUS_SUCCESS;
}

#else

#include <hipsolver/hipsolver.h>

#endif

#endif
