// GAP prediction arithmetic that is not tied to one descriptor: the kernel
// power, the energy offset, and the two mat-vec products over the sparse set.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64

__global__ void gpu_pow(double* a, double* b, double zeta, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < N) {
    double loca = a[idx];
    b[idx] = pow(loca, zeta);
  }
}

extern "C" void gpu_kernels_pow(double* a, double* b, double zeta, int size, hipStream_t* stream) {
  int ntpb = 256;
  int nblocks = (size + ntpb - 1) / ntpb;
  gpu_pow<<<nblocks, ntpb, 0, stream[0]>>>(a, b, zeta, size);
  // gpuErrchk( hipPeekAtLastError() );
  // gpuErrchk( hipDeviceSynchronize() );
  return;
}

__global__ void gpu_simpleaxpc(double* a, double dccc, double e0, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < N) {
    double loca = a[idx];
    a[idx] = dccc * loca + e0;
  }
}

extern "C" void gpu_axpc(double* a, double dccc, double e0, int size, hipStream_t* stream) {
  int ntpb = 256;
  int nblocks = (size + ntpb - 1) / ntpb;
  gpu_simpleaxpc<<<nblocks, ntpb, 0, stream[0]>>>(a, dccc, e0, size);
  /*gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );*/
  return;
}

__global__ void matvect_kernels(double* kernels_d, double* alphas_d, int n_sites, int n_sparse) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int ispa = idx / n_sites;
  int isite = idx % n_sites;
  if (ispa < n_sparse && isite < n_sites) {
    double lock = kernels_d[idx] * alphas_d[ispa];
    kernels_d[idx] = lock;
  }
}

extern "C" void cuda_matvect_kernels(double* kernels_d, double* alphas_d, int n_sites, int n_sparse, hipStream_t* stream) {
  int ntpb = 256;
  int nblocks = (n_sites * n_sparse + ntpb - 1) / ntpb;
  matvect_kernels<<<nblocks, ntpb, 0, stream[0]>>>(kernels_d, alphas_d, n_sites, n_sparse);
  /*gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );*/
  return;
}


__global__ void matvect_qs(double* qs_d, double* qs_copy_d, double* alphas_d, int n_soap, int n_sparse) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int ispa = idx / n_soap;
  int isoap = idx % n_soap;
  if (ispa < n_sparse && isoap < n_soap) {
    double lock = qs_d[idx] * alphas_d[ispa];
    qs_copy_d[idx] = lock;
  }
}

extern "C" void cuda_matvect_qs(double* qs_d, double* qs_copy_d, double* alphas_d, int n_soap, int n_sparse, hipStream_t* stream) {
  /*
     alphas(n_sparse)
     allocate( Qs_copy(1:n_soap, 1:n_sparse) )
     do i = 1, n_soap
     Qs_copy(i,:) = Qs(i,:)*alphas(:)
     end do
   */
  int ntpb = 256;
  int nblocks = (n_soap * n_sparse + ntpb - 1) / ntpb;
  matvect_qs<<<nblocks, ntpb, 0, stream[0]>>>(qs_d, qs_copy_d, alphas_d, n_soap, n_sparse);
  /*gpuErrchk( hipPeekAtLastError() );
    gpuErrchk( hipDeviceSynchronize() );*/
  return;
}


// gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, n_soap, n_sites, cdelta)
