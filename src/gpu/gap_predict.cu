// GAP prediction arithmetic that is not tied to one descriptor: the kernel
// power, the energy offset, and the two mat-vec products over the sparse set.
//
// Four flat maps, written once through tg_parallel_for and dispatched to CUDA
// or Kokkos by src/gpu/gpu_backend.h. The in-kernel bounds guards the
// hand-written launches carried are gone: a <<<>>> grid is rounded up to whole
// blocks and over-runs the extent, tg_parallel_for passes the exact extent.
#include "gap_gpu.h"
#include "gpu_backend.h"
#include "gpu_common.h"

#define tpb 64

// Exponentiation by squaring. Handles a negative base correctly, which matters:
// the kernel matrix is a dot product and is routinely negative, and pow() with
// an exactly-integral exponent is defined there too.
TG_INLINE_FUNCTION double tg_pow_int(double x, int n) {
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
extern "C" void gpu_kernels_pow(double* a, double* b, double zeta, int size, hipStream_t* stream) {
  // -1 means "not a small non-negative whole number, use pow()". The 1e-5
  // tolerance is the one gap.f90 already uses to set is_zeta_int, so the two
  // cannot disagree about whether a given potential's zeta is an integer. Note
  // that the force path calls this with zeta-1, which is why 0 must be allowed.
  int zeta_int = -1;
  double r = round(zeta);
  if (r >= 0.0 && r <= 64.0 && fabs(zeta - r) < 1.0e-5)
    zeta_int = (int) r;

  tg_parallel_for(
      "turbogap_kernels_pow", size, stream,
      TG_LAMBDA(const int idx) {
        const double loca = a[idx];
        // zeta_int is the same for every thread in the grid, so this branch is
        // uniform: no warp diverges on it.
        b[idx] = (zeta_int >= 0) ? tg_pow_int(loca, zeta_int) : pow(loca, zeta);
      },
      256);
}

extern "C" void gpu_axpc(double* a, double dccc, double e0, int size, hipStream_t* stream) {
  tg_parallel_for("turbogap_axpc", size, stream, TG_LAMBDA(const int idx) { a[idx] = dccc * a[idx] + e0; }, 256);
}

// Clamp a local property at zero: the device half of local_property_predict's
// zero_trunc. A Hirshfeld volume or a core-electron binding energy cannot be
// negative, so a negative prediction is the model being wrong rather than a
// value to keep. n_floored_d counts them for the caller's warning and must be
// zeroed before the call.
extern "C" void gpu_zero_trunc(double* v, int* n_floored_d, int size, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_zero_trunc", size, stream,
      TG_LAMBDA(const int idx) {
        if (v[idx] < 0.0) {
          v[idx] = 0.0;
          TG_ATOMIC_ADD(n_floored_d, 1);
        }
      },
      256);
}

// Zero the gradient of every clamped site.
//
// Qss is (n_sites, n_soap) column-major and each of a site's pairs contracts
// the same row, so zeroing the row zeroes all of that site's derivatives at
// once -- what local_properties.f90 does pair by pair on the host. A sum of
// zero products is exactly 0.0, so the two routes agree bit for bit.
//
// The test is v == 0.0 rather than a separate mask because the host's is too:
// it zeroes the gradient of any site whose value is exactly zero, floored or
// not. Matching that quirk is the point -- the two backends must not disagree.
extern "C" void gpu_zero_trunc_der(double* Qss_d, const double* v, int n_sites, int n_soap, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_zero_trunc_der", n_sites * n_soap, stream,
      TG_LAMBDA(const int idx) {
        if (v[idx % n_sites] == 0.0)
          Qss_d[idx] = 0.0;
      },
      256);
}

extern "C" void cuda_matvect_kernels(double* kernels_d, double* alphas_d, int n_sites, int n_sparse, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_matvect_kernels", n_sites * n_sparse, stream,
      TG_LAMBDA(const int idx) { kernels_d[idx] *= alphas_d[idx / n_sites]; }, 256);
}

// alphas(n_sparse); Qs_copy(1:n_soap, 1:n_sparse), Qs_copy(i,:) = Qs(i,:)*alphas(:)
extern "C" void cuda_matvect_qs(double* qs_d, double* qs_copy_d, double* alphas_d, int n_soap, int n_sparse, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_matvect_qs", n_soap * n_sparse, stream,
      TG_LAMBDA(const int idx) { qs_copy_d[idx] = qs_d[idx] * alphas_d[idx / n_soap]; }, 256);
}

// gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, n_soap, n_sites, cdelta)
