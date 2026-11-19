// Many-body dispersion on the device: the energy series and the 2b force
// pass of one site's cluster, for get_mbd_energies_and_forces (src/vdw.f90)
// under the central approximation with the default coefficients. The host
// builds the cluster, its nonzeros and the SCS polarisabilities; each host
// thread hands its site here on its own stream, so the threads keep the device
// busy while they build the next sites.
//
// Every frequency's matrix is the first one scaled by per-atom factors,
// A_w = D_w A_1 D_w with D_w = diag(sqrt(1/(1 + (w/o_p)^2))), so only A_1 is
// stored and one power step of the series is a single product of A_1 with
// 3 x n_freq columns. Sums are in a fixed order (sequential per lane, then a
// fixed shuffle or shared-memory tree), so a run reproduces itself; they do
// not follow the host's order, so device and host agree to round-off.
#include "gap_gpu.h"
#include "gpu_backend.h"
#include "gpu_common.h"

#include <mutex>
#include <vector>

#define MBD_LANES 32
#define MBD_ROWS_PER_BLOCK 4
#define MBD_MAX_COLS 39
#define MBD_FORCE_TPB 256
#define MBD_DOT_TPB 256

#if defined(__HIP_PLATFORM_AMD__)
#define MBD_SHFL_DOWN(v, off) __shfl_down(v, off, MBD_LANES)
#else
#define MBD_SHFL_DOWN(v, off) __shfl_down_sync(0xffffffffu, v, off, MBD_LANES)
#endif

namespace mbd_gpu {

// One site's cluster on the device. Indices into host arrays keep their
// 1-based values where the host stores them so.
struct mbd_site_d {
  int n_atoms; // n_mbd_sites
  int n_pairs; // entries of the cluster pair list, own entries and duplicates included
  int nnz;
  int n_freq;
  int n_order;
  int i0;
  const int* atom_of_pair; // 0-based row atom of each entry
  const int* own_entry;    // 0-based own entry of each atom
  const int* neigh;        // mbd_neighbors_list
  const int* p_mbd;        // 1-based column atom, -1 for a duplicate
  const int* k4_start;     // 0-based first nonzero of each pair in COO order, -1 if none
  const double* omegas;
  const double* a_mbd;
  const double* o_mbd;
  const double* r0_scs;
  const double* f_damp;
  const double* t_lr_mult;
  const double* rjs;
  const double* xyz;   // 3 x n_pairs, pair vectors
  const double* xyz0;  // 3 x n_pairs, positions from the central atom
  const double* t_mbd; // 9 x n_pairs
  const double* pcut;
  const double* v1; // val_sym(:, 1), COO order
  const double* v1_csr;
  const int* i_csr;       // 1-based
  const int* j_csr;       // 1-based
  const double* t_lr_sym; // 3 x 3N
};

struct mbd_consts {
  double d, sR, Bohr, rcut_mbd, rcut_mbd2, r_buf_ij;
};

TG_INLINE_FUNCTION double mbd_pow5(double x) {
  double x2 = x * x;
  return x2 * (x * x2);
}

TG_INLINE_FUNCTION double mbd_pow7(double x) {
  double x2 = x * x;
  return (x * x2) * (x2 * x2);
}

// s(p, w) = sqrt(1/(1 + (w/o)^2)) and q(p, w) = sqrt(a/(1 + (w/o)^2)), atom-major.
void mbd_freq_factors(mbd_site_d s, double* s_d, double* q_d, hipStream_t* stream) {
  int n = s.n_atoms;
  tg_parallel_for(
      "mbd_freq_factors", n * s.n_freq, stream, TG_LAMBDA(const int idx) {
        int p = idx % n;
        int w = idx / n;
        int k = s.own_entry[p];
        double t = s.omegas[w] / s.o_mbd[k];
        double den = 1.0 + t * t;
        s_d[idx] = sqrt(1.0 / den);
        q_d[idx] = sqrt(s.a_mbd[k] / den);
      });
}

// AT_sym(r, c, w) = q(p(r), w) * (q(1, w) * T_LR_sym(c, r)); layout r + 3N*(c + 3*w).
void mbd_at_sym(mbd_site_d s, const double* q_d, double* at_d, hipStream_t* stream) {
  int n3 = 3 * s.n_atoms;
  int n = s.n_atoms;
  tg_parallel_for(
      "mbd_at_sym", n3 * 3 * s.n_freq, stream, TG_LAMBDA(const int idx) {
        int r = idx % n3;
        int c = (idx / n3) % 3;
        int w = idx / (3 * n3);
        at_d[idx] = q_d[r / 3 + n * w] * (q_d[0 + n * w] * s.t_lr_sym[c + 3 * r]);
      });
}

void mbd_csr_values(mbd_site_d s, const int* perm_d, double* v1_csr_d, hipStream_t* stream) {
  tg_parallel_for("mbd_csr_values", s.nnz, stream, TG_LAMBDA(const int k) { v1_csr_d[k] = s.v1[perm_d[k] - 1]; });
}

// Z(r, c, w) = s(p(r), w) X(r, c, w).
void mbd_scale_rows(mbd_site_d s, const double* s_d, const double* x_d, double* z_d, hipStream_t* stream) {
  int n3 = 3 * s.n_atoms;
  int n = s.n_atoms;
  tg_parallel_for(
      "mbd_scale_rows", n3 * 3 * s.n_freq, stream, TG_LAMBDA(const int idx) {
        int r = idx % n3;
        int w = idx / (3 * n3);
        z_d[idx] = s_d[r / 3 + n * w] * x_d[idx];
      });
}

// Columns c0.. of Y = D A_1 Z: 32 lanes per row, each summing every 32nd entry
// of the row, then a shuffle tree. Then pol -= k2*res(k2+1)*Y and, past the
// first power, the diagonal kept for the integrand.
__global__ void mbd_spmm_cols_kernel(mbd_site_d s, int k2, int c0, double pol_coeff, const double* s_d, const double* z_d,
                                     double* y_d, double* pol_d, double* diag_d) {
  int n3 = 3 * s.n_atoms;
  int n = s.n_atoms;
  int ncol = min(MBD_MAX_COLS, 3 * s.n_freq - c0);
  int lane = threadIdx.x % MBD_LANES;
  int r = blockIdx.x * MBD_ROWS_PER_BLOCK + threadIdx.x / MBD_LANES;
  double acc[MBD_MAX_COLS];
#pragma unroll
  for (int c = 0; c < MBD_MAX_COLS; c++)
    acc[c] = 0.0;
  if (r < n3) {
    for (int k = s.i_csr[r] - 1 + lane; k < s.i_csr[r + 1] - 1; k += MBD_LANES) {
      int j = s.j_csr[k] - 1;
      double v = s.v1_csr[k];
#pragma unroll
      for (int c = 0; c < MBD_MAX_COLS; c++)
        if (c < ncol)
          acc[c] = acc[c] + v * z_d[j + (size_t) n3 * (c0 + c)];
    }
  }
  for (int off = MBD_LANES / 2; off > 0; off /= 2) {
#pragma unroll
    for (int c = 0; c < MBD_MAX_COLS; c++)
      acc[c] = acc[c] + MBD_SHFL_DOWN(acc[c], off);
  }
  if (lane == 0 && r < n3) {
#pragma unroll
    for (int c = 0; c < MBD_MAX_COLS; c++) {
      if (c >= ncol)
        continue;
      int gc = c0 + c;
      int w = gc / 3;
      int cc = gc % 3;
      size_t at = r + (size_t) n3 * gc;
      double y = s_d[r / 3 + n * w] * acc[c];
      y_d[at] = y;
      pol_d[at] = pol_d[at] - pol_coeff * y;
      if (k2 > 3 && r == cc)
        diag_d[(k2 - 4) * 3 + r + 3 * (s.n_order - 3) * w] = y;
    }
  }
}

void mbd_series_step(mbd_site_d s, int k2, const double* s_d, const double* x_d, double* z_d, double* y_d, double* pol_d,
                     double* diag_d, hipStream_t* stream) {
  double pol_coeff = (double) k2 * (-1.0 / (double) k2);
  mbd_scale_rows(s, s_d, x_d, z_d, stream);
  int blocks = (3 * s.n_atoms + MBD_ROWS_PER_BLOCK - 1) / MBD_ROWS_PER_BLOCK;
  for (int c0 = 0; c0 < 3 * s.n_freq; c0 += MBD_MAX_COLS)
    mbd_spmm_cols_kernel<<<blocks, MBD_ROWS_PER_BLOCK * MBD_LANES, 0, stream[0]>>>(s, k2, c0, pol_coeff, s_d, z_d, y_d, pol_d,
                                                                                   diag_d);
}

// The series' last terms, one block per (w, c, kind): kind 0 AT.X_last, kind 1 AT.AT.
__global__ void mbd_series_dots_kernel(mbd_site_d s, const double* at_d, const double* x_d, double* dots_d) {
  __shared__ double sh[MBD_DOT_TPB];
  int n3 = 3 * s.n_atoms;
  int c = blockIdx.x % 3;
  int w = (blockIdx.x / 3) % s.n_freq;
  int kind = blockIdx.x / (3 * s.n_freq);
  const double* a = at_d + (size_t) n3 * (c + 3 * w);
  const double* b = (kind == 0 ? x_d : at_d) + (size_t) n3 * (c + 3 * w);
  double t = 0.0;
  for (int r = threadIdx.x; r < n3; r += MBD_DOT_TPB)
    t = t + a[r] * b[r];
  sh[threadIdx.x] = t;
  __syncthreads();
  for (int off = MBD_DOT_TPB / 2; off > 0; off /= 2) {
    if (threadIdx.x < off)
      sh[threadIdx.x] = sh[threadIdx.x] + sh[threadIdx.x + off];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    dots_d[blockIdx.x] = sh[0];
}

// integrand_sym(w) from the kept diagonals and the dots, in the host's order.
void mbd_series_finish(mbd_site_d s, const double* at_d, const double* x_d, const double* diag_d, double* dots_d,
                       double* integrand_d, hipStream_t* stream) {
  mbd_series_dots_kernel<<<6 * s.n_freq, MBD_DOT_TPB, 0, stream[0]>>>(s, at_d, x_d, dots_d);
  tg_parallel_for(
      "mbd_series_integrand", s.n_freq, stream, TG_LAMBDA(const int w) {
        double integ = 0.0;
        for (int k2 = 4; k2 <= s.n_order; k2++) {
          double res = -1.0 / (double) (k2 - 1);
          for (int c = 0; c < 3; c++)
            integ = integ + res * diag_d[(k2 - 4) * 3 + c + 3 * (s.n_order - 3) * w];
        }
        double res_last = -1.0 / (double) (s.n_order + 1 - 1);
        for (int c = 0; c < 3; c++)
          integ = integ + res_last * dots_d[c + 3 * w];
        double res3 = -1.0 / (double) (3 - 1);
        for (int c = 0; c < 3; c++)
          integ = integ + res3 * dots_d[3 * s.n_freq + c + 3 * w];
        integrand_d[w] = integ;
      });
}

// fq_t(p, w) for this c3: da/(1 + (w/o)^2) + a*(2 w^2 o)*do/(o^2 + w^2)^2; da and dO per atom.
void mbd_freq_derivs(mbd_site_d s, const double* da, const double* dO, double* t_d, hipStream_t* stream) {
  int n = s.n_atoms;
  tg_parallel_for(
      "mbd_freq_derivs", n * s.n_freq, stream, TG_LAMBDA(const int idx) {
        int p = idx % n;
        int w = idx / n;
        int k = s.own_entry[p];
        double om = s.omegas[w];
        double o = s.o_mbd[k];
        double t = om / o;
        double den = o * o + om * om;
        t_d[idx] = da[p] / (1.0 + t * t) + ((s.a_mbd[k] * ((2.0 * (om * om)) * o)) * dO[p]) / (den * den);
      });
}

// One thread per pair entry: the frequency-independent part once (dT_LR_val,
// and dT_LR_sym for the central atom's pairs), then for each frequency the
// pair's dot, virial and damped virial, reduced over each warp by a shuffle
// tree and over the block's warps in order into part_d (8 per block and w;
// the 8th is the dot's magnitude, kept for tests).
__global__ void mbd_force_kernel(mbd_site_d s, mbd_consts k, int c3, const double* dr0, const double* s_d, const double* q_d,
                                 const double* t_d, double* part_d, double* dtlr_sym_d) {
  __shared__ double sh[MBD_FORCE_TPB / MBD_LANES][8];
  int k2 = blockIdx.x * MBD_FORCE_TPB + threadIdx.x;
  int lane = threadIdx.x % MBD_LANES;
  int warp = threadIdx.x / MBD_LANES;
  int n = s.n_atoms;
  int c3m = c3 - 1;
  bool live = k2 < s.n_pairs;
  int p = 0, q = 0, kb = 0, k4 = 0;
  if (live) {
    p = s.atom_of_pair[k2];
    kb = s.own_entry[p];
    live = !(k2 == kb || s.p_mbd[k2] == -1);
  }
  double dtv[9];
  double f = 0.0, tlr = 0.0, pcut_i = 0.0, pcut_next = 0.0;
  double xi[3] = {0.0, 0.0, 0.0}, xn[3] = {0.0, 0.0, 0.0};
  bool dup = false;
  if (live) {
    q = s.p_mbd[k2] - 1;
    k4 = s.k4_start[k2];
    int i2 = s.neigh[kb];
    int j = s.neigh[k2];
    double rjs = s.rjs[k2];
    f = s.f_damp[k2];
    tlr = s.t_lr_mult[k2];
    const double* x = s.xyz + 3 * k2;
    double S = k.sR * (s.r0_scs[kb] + s.r0_scs[k2]);
    double dS = k.sR * (dr0[p] + dr0[q]);
    double e = exp(-k.d * (rjs / S - 1.0));
    double fdds = -((((k.d * rjs) / (S * S)) * (f * f)) * e) * dS;
    for (int c = 0; c < 9; c++)
      dtv[c] = 0.0 + (s.t_mbd[9 * k2 + c] * fdds) * tlr;
    if (s.i0 == i2 || s.i0 == j) {
      double fdm = ((((k.d / S) * (f * f)) * e) * x[c3m]) / rjs;
      double dtm = 0.0;
      if (p == 0) {
        double rb;
        if (rjs > (k.rcut_mbd - k.r_buf_ij) / k.Bohr && rjs <= k.rcut_mbd / k.Bohr) {
          rb = (rjs * k.Bohr - k.rcut_mbd + k.r_buf_ij) / k.r_buf_ij;
          if (s.i0 == i2)
            dtm = (-6.0 * rb + 6.0 * (rb * rb)) * (-x[c3m] / rjs / (k.r_buf_ij / k.Bohr));
          else if (s.i0 == j)
            dtm = (-6.0 * rb + 6.0 * (rb * rb)) * (x[c3m] / rjs / (k.r_buf_ij / k.Bohr));
        } else if (rjs > (k.rcut_mbd2 - k.r_buf_ij) / k.Bohr && rjs <= k.rcut_mbd2 / k.Bohr) {
          rb = (rjs * k.Bohr - k.rcut_mbd2 + k.r_buf_ij) / k.r_buf_ij;
          if (s.i0 == i2)
            dtm = (-6.0 * rb + 6.0 * (rb * rb)) * (-x[c3m] / rjs / (k.r_buf_ij / k.Bohr));
          else if (s.i0 == j)
            dtm = (-6.0 * rb + 6.0 * (rb * rb)) * (x[c3m] / rjs / (k.r_buf_ij / k.Bohr));
        }
      }
      double r5 = mbd_pow5(rjs);
      double r7 = mbd_pow7(rjs);
      for (int c1 = 0; c1 < 3; c1++)
        for (int c2 = 0; c2 < 3; c2++) {
          int c = 3 * c1 + c2;
          double dt = (((-15.0 * x[c1]) * x[c2]) * x[c3m]) / r7;
          if (c1 == c2)
            dt = dt + (3.0 / r5) * x[c3m];
          if (c2 == c3m)
            dt = dt + (3.0 / r5) * x[c1];
          if (c1 == c3m)
            dt = dt + (3.0 / r5) * x[c2];
          double tm = s.t_mbd[9 * k2 + c];
          if (s.i0 == i2)
            dtv[c] = dtv[c] - (tm * fdm + dt * f) * tlr;
          else if (s.i0 == j)
            dtv[c] = dtv[c] + (tm * fdm + dt * f) * tlr;
          dtv[c] = dtv[c] + (f * tm) * dtm;
        }
    }
    if (p == 0)
      for (int c1 = 0; c1 < 3; c1++)
        for (int c2 = 0; c2 < 3; c2++)
          dtlr_sym_d[c1 + 3 * (3 * q + c2)] = dtv[3 * c1 + c2];
    for (int c = 0; c < 3; c++)
      xi[c] = s.xyz0[3 * kb + c];
    pcut_i = s.pcut[kb];
    dup = (k2 + 2 < s.n_pairs) && (s.p_mbd[k2 + 1] == -1);
    if (dup) {
      for (int c = 0; c < 3; c++)
        xn[c] = s.xyz0[3 * (k2 + 1) + c];
      pcut_next = s.pcut[k2 + 1];
    }
  }
  for (int w = 0; w < s.n_freq; w++) {
    double acc[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (live) {
      double fq_si = q_d[p + n * w], fq_sj = q_d[q + n * w];
      double fq_ss = fq_si * fq_sj;
      double fq_gi = ((0.5 / fq_si) * t_d[p + n * w]) * fq_sj;
      double fq_gj = ((0.5 / fq_sj) * t_d[q + n * w]) * fq_si;
      double s_pq = s_d[p + n * w] * s_d[q + n * w];
      int kk = k4;
      for (int c1 = 0; c1 < 3; c1++)
        for (int c2 = 0; c2 < 3; c2++) {
          double tm = s.t_mbd[9 * k2 + 3 * c1 + c2];
          if (!(fabs(f * tm * tlr) > 1.0e-20))
            continue;
          double dval = fq_ss * dtv[3 * c1 + c2] + ((fq_gi * f) * tm) * tlr + ((fq_gj * f) * tm) * tlr;
          double val = s_pq * s.v1[kk];
          for (int c = 0; c < 3; c++)
            acc[1 + c] = acc[1 + c] - ((0.5 * xi[c]) * dval) * val;
          for (int c = 0; c < 3; c++)
            acc[4 + c] = acc[4 + c] - (((0.5 * xi[c]) * dval) * val) * pcut_i;
          acc[0] = acc[0] + val * dval;
          acc[7] = acc[7] + fabs(val * dval);
          kk++;
          if (dup) {
            double val2 = s_pq * s.v1[kk];
            for (int c = 0; c < 3; c++)
              acc[1 + c] = acc[1 + c] - ((0.5 * xn[c]) * dval) * val2;
            for (int c = 0; c < 3; c++)
              acc[4 + c] = acc[4 + c] - (((0.5 * xn[c]) * dval) * val2) * pcut_next;
            acc[0] = acc[0] + val2 * dval;
            acc[7] = acc[7] + fabs(val2 * dval);
            kk++;
          }
        }
    }
    for (int off = MBD_LANES / 2; off > 0; off /= 2) {
#pragma unroll
      for (int c = 0; c < 8; c++)
        acc[c] = acc[c] + MBD_SHFL_DOWN(acc[c], off);
    }
    if (lane == 0)
#pragma unroll
      for (int c = 0; c < 8; c++)
        sh[warp][c] = acc[c];
    __syncthreads();
    if (threadIdx.x < 8) {
      double t = 0.0;
      for (int wp = 0; wp < MBD_FORCE_TPB / MBD_LANES; wp++)
        t = t + sh[wp][threadIdx.x];
      part_d[8 * (blockIdx.x + gridDim.x * w) + threadIdx.x] = t;
    }
    __syncthreads();
  }
}

// Block sums in block order: out[8*w + c].
void mbd_force_reduce(mbd_site_d s, int n_blocks, const double* part_d, double* out_d, hipStream_t* stream) {
  tg_parallel_for(
      "mbd_force_reduce", 8 * s.n_freq, stream, TG_LAMBDA(const int idx) {
        int c = idx % 8;
        int w = idx / 8;
        double t = 0.0;
        for (int b = 0; b < n_blocks; b++)
          t = t + part_d[8 * (b + n_blocks * w) + c];
        out_d[idx] = t;
      });
}

template <class T> void grow(T*& p, size_t& cap, size_t n, hipStream_t* stream) {
  if (n <= cap)
    return;
  if (p)
    gpuErrchk(hipFreeAsync(p, stream[0]));
  gpuErrchk(hipMallocAsync((void**) &p, sizeof(T) * n, stream[0]));
  cap = n;
}

template <class T> void put(T*& p, size_t& cap, const T* host, size_t n, hipStream_t* stream) {
  grow(p, cap, n > 0 ? n : 1, stream);
  if (n)
    gpuErrchk(hipMemcpyAsync(p, host, sizeof(T) * n, hipMemcpyHostToDevice, stream[0]));
}

// A host thread's stream and buffers, kept between sites and grown when a
// site needs more.
struct mbd_ctx {
  hipStream_t stream;
  mbd_site_d s;
  mbd_consts k;
  int n_blocks;
  std::vector<int> own, atom_of, k4s;
  std::vector<double> da_atom, do_atom, dr0_atom, out;
  int *own_d = nullptr, *atom_of_d = nullptr, *neigh_d = nullptr, *p_mbd_d = nullptr, *k4s_d = nullptr;
  int *i_csr_d = nullptr, *j_csr_d = nullptr, *perm_d = nullptr;
  double *omegas_d = nullptr, *a_d = nullptr, *o_d = nullptr, *r0_d = nullptr, *f_d = nullptr, *tlr_d = nullptr;
  double *rjs_d = nullptr, *xyz_d = nullptr, *xyz0_d = nullptr, *t_mbd_d = nullptr, *pcut_d = nullptr;
  double *v1_d = nullptr, *v1_csr_d = nullptr, *t_lr_sym_d = nullptr;
  double *s_d = nullptr, *q_d = nullptr, *t_d = nullptr, *at_d = nullptr, *x_d = nullptr, *y_d = nullptr, *z_d = nullptr;
  double *pol_d = nullptr, *diag_d = nullptr, *dots_d = nullptr, *integrand_d = nullptr;
  double *da_d = nullptr, *do_d = nullptr, *dr0_d = nullptr, *part_d = nullptr, *out_d = nullptr, *dtlr_d = nullptr;
  size_t c_own = 0, c_atom_of = 0, c_neigh = 0, c_p_mbd = 0, c_k4s = 0, c_i_csr = 0, c_j_csr = 0, c_perm = 0;
  size_t c_omegas = 0, c_a = 0, c_o = 0, c_r0 = 0, c_f = 0, c_tlr = 0, c_rjs = 0, c_xyz = 0, c_xyz0 = 0, c_t_mbd = 0;
  size_t c_pcut = 0, c_v1 = 0, c_v1_csr = 0, c_t_lr_sym = 0, c_s = 0, c_q = 0, c_t = 0, c_at = 0, c_x = 0, c_y = 0;
  size_t c_z = 0, c_pol = 0, c_diag = 0, c_dots = 0, c_integrand = 0, c_da = 0, c_do = 0, c_dr0 = 0, c_part = 0;
  size_t c_out = 0, c_dtlr = 0;
};

int g_device = 0;
std::mutex g_mutex;
std::vector<mbd_ctx*> g_ctx;

mbd_ctx* context(int tid) {
  std::lock_guard<std::mutex> lock(g_mutex);
  if ((int) g_ctx.size() <= tid)
    g_ctx.resize(tid + 1, nullptr);
  if (!g_ctx[tid]) {
    g_ctx[tid] = new mbd_ctx();
    gpuErrchk(hipSetDevice(g_device));
    gpuErrchk(hipStreamCreateWithFlags(&g_ctx[tid]->stream, hipStreamNonBlocking));
  }
  return g_ctx[tid];
}

// Own entries, row atoms, and the first nonzero of each pair in the order the
// host walks them (the filter the host applied when it kept nonzeros).
void index_pairs(mbd_ctx& c, int n_atoms, int n_pairs, int nnz, const int* n_neigh, const int* p_mbd, const double* f_damp,
                 const double* t_mbd, const double* tlr) {
  c.own.assign(n_atoms, 0);
  c.atom_of.assign(n_pairs, 0);
  c.k4s.assign(n_pairs, -1);
  int k2 = 0;
  for (int p = 0; p < n_atoms; p++) {
    c.own[p] = k2;
    for (int j = 0; j < n_neigh[p]; j++)
      c.atom_of[k2 + j] = p;
    k2 += n_neigh[p];
  }
  int k4 = 0;
  for (int p = 0; p < n_atoms; p++) {
    for (int k = c.own[p] + 1; k < c.own[p] + n_neigh[p]; k++) {
      if (p_mbd[k] == -1)
        continue;
      bool dup = (k + 2 < n_pairs) && p_mbd[k + 1] == -1;
      for (int cc = 0; cc < 9; cc++) {
        if (!(fabs(f_damp[k] * t_mbd[9 * (size_t) k + cc] * tlr[k]) > 1.0e-20))
          continue;
        if (c.k4s[k] < 0)
          c.k4s[k] = k4;
        k4 += dup ? 2 : 1;
      }
    }
  }
  if (k4 != nnz) {
    fprintf(stderr, "ERROR: MBD on the device counts %d nonzeros where the host kept %d\n", k4, nnz);
    exit(1);
  }
}

} // namespace mbd_gpu

using namespace mbd_gpu;

extern "C" {

void mbd_gpu_init() {
  gpuErrchk(hipGetDevice(&g_device));
}

void mbd_gpu_site(int tid, int n_atoms, int n_pairs, int nnz, int n_freq, int n_order, int i0, double d, double sR, double Bohr,
                  double rcut_mbd, double rcut_mbd2, double r_buf_ij, const int* n_neigh, const int* neigh, const int* p_mbd,
                  const double* omegas, const double* a_mbd, const double* o_mbd, const double* r0_scs, const double* f_damp,
                  const double* t_lr_mult, const double* rjs, const double* xyz, const double* xyz0, const double* t_mbd,
                  const double* pcut, const double* v1, const int* i_csr, const int* j_csr, const int* perm, const double* t_lr_sym,
                  double* integrand, double* pol) {
  mbd_ctx& c = *context(tid);
  gpuErrchk(hipSetDevice(g_device));
  hipStream_t* st = &c.stream;
  index_pairs(c, n_atoms, n_pairs, nnz, n_neigh, p_mbd, f_damp, t_mbd, t_lr_mult);
  size_t n3f = (size_t) 9 * n_atoms * n_freq;
  size_t naf = (size_t) n_atoms * n_freq;
  put(c.own_d, c.c_own, c.own.data(), n_atoms, st);
  put(c.atom_of_d, c.c_atom_of, c.atom_of.data(), n_pairs, st);
  put(c.k4s_d, c.c_k4s, c.k4s.data(), n_pairs, st);
  put(c.neigh_d, c.c_neigh, neigh, n_pairs, st);
  put(c.p_mbd_d, c.c_p_mbd, p_mbd, n_pairs, st);
  put(c.i_csr_d, c.c_i_csr, i_csr, 3 * (size_t) n_atoms + 1, st);
  put(c.j_csr_d, c.c_j_csr, j_csr, nnz, st);
  put(c.perm_d, c.c_perm, perm, nnz, st);
  put(c.omegas_d, c.c_omegas, omegas, n_freq, st);
  put(c.a_d, c.c_a, a_mbd, n_pairs, st);
  put(c.o_d, c.c_o, o_mbd, n_pairs, st);
  put(c.r0_d, c.c_r0, r0_scs, n_pairs, st);
  put(c.f_d, c.c_f, f_damp, n_pairs, st);
  put(c.tlr_d, c.c_tlr, t_lr_mult, n_pairs, st);
  put(c.rjs_d, c.c_rjs, rjs, n_pairs, st);
  put(c.xyz_d, c.c_xyz, xyz, 3 * (size_t) n_pairs, st);
  put(c.xyz0_d, c.c_xyz0, xyz0, 3 * (size_t) n_pairs, st);
  put(c.t_mbd_d, c.c_t_mbd, t_mbd, 9 * (size_t) n_pairs, st);
  put(c.pcut_d, c.c_pcut, pcut, n_pairs, st);
  put(c.v1_d, c.c_v1, v1, nnz, st);
  put(c.t_lr_sym_d, c.c_t_lr_sym, t_lr_sym, 9 * (size_t) n_atoms, st);
  grow(c.v1_csr_d, c.c_v1_csr, nnz > 0 ? nnz : 1, st);
  grow(c.s_d, c.c_s, naf, st);
  grow(c.q_d, c.c_q, naf, st);
  grow(c.t_d, c.c_t, naf, st);
  grow(c.at_d, c.c_at, n3f, st);
  grow(c.x_d, c.c_x, n3f, st);
  grow(c.y_d, c.c_y, n3f, st);
  grow(c.z_d, c.c_z, n3f, st);
  grow(c.pol_d, c.c_pol, n3f, st);
  grow(c.diag_d, c.c_diag, (size_t) 3 * (n_order - 3) * n_freq, st);
  grow(c.dots_d, c.c_dots, (size_t) 6 * n_freq, st);
  grow(c.integrand_d, c.c_integrand, n_freq, st);

  mbd_site_d& s = c.s;
  s = {n_atoms,   n_pairs,  nnz,        n_freq,   n_order, i0,         c.atom_of_d, c.own_d,   c.neigh_d,
       c.p_mbd_d, c.k4s_d,  c.omegas_d, c.a_d,    c.o_d,   c.r0_d,     c.f_d,       c.tlr_d,   c.rjs_d,
       c.xyz_d,   c.xyz0_d, c.t_mbd_d,  c.pcut_d, c.v1_d,  c.v1_csr_d, c.i_csr_d,   c.j_csr_d, c.t_lr_sym_d};
  c.k = {d, sR, Bohr, rcut_mbd, rcut_mbd2, r_buf_ij};
  c.n_blocks = (n_pairs + MBD_FORCE_TPB - 1) / MBD_FORCE_TPB;

  mbd_freq_factors(s, c.s_d, c.q_d, st);
  mbd_at_sym(s, c.q_d, c.at_d, st);
  mbd_csr_values(s, c.perm_d, c.v1_csr_d, st);
  gpuErrchk(hipMemcpyAsync(c.x_d, c.at_d, sizeof(double) * n3f, hipMemcpyDeviceToDevice, st[0]));
  gpuErrchk(hipMemsetAsync(c.pol_d, 0, sizeof(double) * n3f, st[0]));
  for (int k2 = 3; k2 <= n_order; k2++) {
    mbd_series_step(s, k2, c.s_d, c.x_d, c.z_d, c.y_d, c.pol_d, c.diag_d, st);
    double* t = c.x_d;
    c.x_d = c.y_d;
    c.y_d = t;
    size_t tc = c.c_x;
    c.c_x = c.c_y;
    c.c_y = tc;
  }
  mbd_series_finish(s, c.at_d, c.x_d, c.diag_d, c.dots_d, c.integrand_d, st);
  gpuErrchk(hipMemcpyAsync(integrand, c.integrand_d, sizeof(double) * n_freq, hipMemcpyDeviceToHost, st[0]));
  gpuErrchk(hipMemcpyAsync(pol, c.pol_d, sizeof(double) * n3f, hipMemcpyDeviceToHost, st[0]));
  gpuErrchk(hipStreamSynchronize(st[0]));
}

void mbd_gpu_force(int tid, int c3, const double* da, const double* dO, const double* dr0, double* dot_2b, double* virial_2b,
                   double* damped_2b, double* dt_lr_sym) {
  mbd_ctx& c = *context(tid);
  gpuErrchk(hipSetDevice(g_device));
  hipStream_t* st = &c.stream;
  mbd_site_d& s = c.s;
  c.da_atom.resize(s.n_atoms);
  c.do_atom.resize(s.n_atoms);
  c.dr0_atom.resize(s.n_atoms);
  for (int p = 0; p < s.n_atoms; p++) {
    c.da_atom[p] = da[c.own[p]];
    c.do_atom[p] = dO[c.own[p]];
    c.dr0_atom[p] = dr0[c.own[p]];
  }
  put(c.da_d, c.c_da, c.da_atom.data(), s.n_atoms, st);
  put(c.do_d, c.c_do, c.do_atom.data(), s.n_atoms, st);
  put(c.dr0_d, c.c_dr0, c.dr0_atom.data(), s.n_atoms, st);
  grow(c.part_d, c.c_part, (size_t) 8 * c.n_blocks * s.n_freq, st);
  grow(c.out_d, c.c_out, (size_t) 8 * s.n_freq, st);
  grow(c.dtlr_d, c.c_dtlr, 9 * (size_t) s.n_atoms, st);
  gpuErrchk(hipMemsetAsync(c.dtlr_d, 0, sizeof(double) * 9 * s.n_atoms, st[0]));
  mbd_freq_derivs(s, c.da_d, c.do_d, c.t_d, st);
  mbd_force_kernel<<<c.n_blocks, MBD_FORCE_TPB, 0, st[0]>>>(s, c.k, c3, c.dr0_d, c.s_d, c.q_d, c.t_d, c.part_d, c.dtlr_d);
  mbd_force_reduce(s, c.n_blocks, c.part_d, c.out_d, st);
  c.out.resize(8 * (size_t) s.n_freq);
  gpuErrchk(hipMemcpyAsync(c.out.data(), c.out_d, sizeof(double) * 8 * s.n_freq, hipMemcpyDeviceToHost, st[0]));
  gpuErrchk(hipMemcpyAsync(dt_lr_sym, c.dtlr_d, sizeof(double) * 9 * s.n_atoms, hipMemcpyDeviceToHost, st[0]));
  gpuErrchk(hipStreamSynchronize(st[0]));
  for (int w = 0; w < s.n_freq; w++) {
    dot_2b[w] = c.out[8 * w];
    for (int cc = 0; cc < 3; cc++) {
      virial_2b[cc + 3 * w] = c.out[8 * w + 1 + cc];
      damped_2b[cc + 3 * w] = c.out[8 * w + 4 + cc];
    }
  }
}

} // extern "C"
