// GAP prediction arithmetic that is not tied to one descriptor: the kernel
// power, the energy offset, and the two mat-vec products over the sparse set.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64

// Exponentiation by squaring. Handles a negative base correctly, which matters:
// the kernel matrix is a dot product and is routinely negative, and pow() with
// an exactly-integral exponent is defined there too.
__device__ __forceinline__ double tg_pow_int(double x, int n) {
  double r = 1.0;
  double b = x;
  while (n > 0) {
    if (n & 1)
      r *= b;
    b *= b;
    n >>= 1;
  }
  return r;
}

// Raise every element of the kernel matrix to zeta.
//
// zeta is a hyperparameter read from the potential file, and every published
// TurboGAP potential uses a small whole number -- 6 in carbon.gap. The device
// pow(double, double) cannot know that: it is the general exp2(zeta*log2(x))
// routine with its argument reduction and correction terms, a few hundred FP64
// instructions per element.
//
// Measured on 124,959 atoms before this change: gpu_pow raised 8.26 M elements
// per launch in 13.4 ms. That is 132 MB of traffic in 13.4 ms -- 9.9 GB/s
// against the A2000's ~288, i.e. 3% of bandwidth. The kernel was not moving
// memory, it was evaluating pow(). Over the 260 launches of a single point it
// came to 3.47 s of a 35.8 s run, 9.7%, making it the third most expensive
// kernel in the whole profile behind the three-body kernel and the SOAP
// derivative.
//
// The Fortran already knew this and did nothing with it: gap.f90 and
// local_properties.f90 both derive `zeta_int` and `is_zeta_int` from the input
// and then pass the double anyway, while the CPU paths in vdw.f90 and
// local_properties.f90 use `K**zeta_int`. So the integer path also brings the
// two backends closer together rather than further apart.
//
// Squaring rounds a different number of times than pow does, so results move in
// the last one or two ulp -- far inside the ~1e-10 at which the GPU and CPU
// already differ. The regression suite is the check.
__global__ void gpu_pow(double* a, double* b, double zeta, int zeta_int, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < N) {
    double loca = a[idx];
    // zeta_int is the same for every thread in the grid, so this branch is
    // uniform: no warp diverges on it.
    b[idx] = (zeta_int >= 0) ? tg_pow_int(loca, zeta_int) : pow(loca, zeta);
  }
}

extern "C" void gpu_kernels_pow(double* a, double* b, double zeta, int size, hipStream_t* stream) {
  int ntpb = 256;
  int nblocks = (size + ntpb - 1) / ntpb;
  // -1 means "not a small non-negative whole number, use pow()". The 1e-5
  // tolerance is the one gap.f90 already uses to set is_zeta_int, so the two
  // cannot disagree about whether a given potential's zeta is an integer. Note
  // that the force path calls this with zeta-1, which is why 0 must be allowed.
  int zeta_int = -1;
  double r = round(zeta);
  if (r >= 0.0 && r <= 64.0 && fabs(zeta - r) < 1.0e-5)
    zeta_int = (int) r;
  gpu_pow<<<nblocks, ntpb, 0, stream[0]>>>(a, b, zeta, zeta_int, size);
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
