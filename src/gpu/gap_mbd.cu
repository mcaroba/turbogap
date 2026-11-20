// Many-body dispersion on the device: the energy series and the 2b force
// pass of one site's cluster, for get_mbd_energies_and_forces (src/vdw.f90)
// under the central approximation with the default coefficients. The host
// builds the cluster's pair list, its atoms' polarisabilities and the central
// atom's row of the dipole matrix; the device forms every pair's values. Each
// host thread hands its site here on its own stream, so the threads keep the
// device busy while they build the next sites.
//
// Every frequency's matrix is the first one scaled by per-atom factors,
// A_w = D_w A_1 D_w with D_w = diag(sqrt(1/(1 + (w/o_p)^2))), so only A_1 is
// stored, as one 3 x 3 block per pair, and one power step of the series is a
// single product of A_1 with 3 x n_freq columns. Sums are in a fixed order (sequential per lane, then a
// fixed shuffle or shared-memory tree), so a run reproduces itself; they do
// not follow the host's order, so device and host agree to round-off.
#include "gap_gpu.h"
#include "gpu_backend.h"
#include "gpu_common.h"

#include <mutex>
#include <vector>

#define MBD_LANES 32
#define MBD_ROWS_PER_BLOCK 4
#define MBD_COL_SLOTS 2
#define MBD_REDUCE_TPB 256
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
  int n_freq;
  int n_order;
  int i0;
  const int* atom_of_pair; // 0-based row atom of each entry
  const int* own_entry;    // 0-based own entry of each atom
  const int* gidx;         // mbd_neighbors_list at each atom's own entry
  const int* p_mbd;        // 1-based column atom, -1 for a duplicate
  const int* bptr;         // block row p is blocks bptr[p] .. bptr[p+1]-1
  const int* bcol;         // a block's column atom
  const int* bpos;         // a pair entry's block, -1 for own entries and duplicates
  const int* bdup;         // the block a central pair's duplicate fills, else -1
  const double* omegas;
  const double* a_own; // a_mbd, o_mbd, r0_ii_SCS at each atom's own entry
  const double* o_own;
  const double* r0_own;
  const double* a_nb; // a_mbd and r0_ii_SCS of the atom as a pair's column atom
  const double* r0_nb;
  const double* pos;  // 3 x n_atoms, positions from the central atom
  const double* pcut; // per atom
  double* f_damp;     // per pair entry, formed here
  double* t_lr_mult;
  double* blk;            // 9 per block, in block order: frequency-1 values, 0 where the host drops one
  int* mask;              // the nonzeros, bit 3*c1 + c2
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
        double t = s.omegas[w] / s.o_own[p];
        double den = 1.0 + t * t;
        s_d[idx] = sqrt(1.0 / den);
        q_d[idx] = sqrt(s.a_own[p] / den);
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

// The pair vector, length and dipole tensor (the host's order, no fused
// multiply-adds).
TG_INLINE_FUNCTION double mbd_pair_tensor(const mbd_site_d& s, int p, int q, double* x, double* tm9) {
  for (int c = 0; c < 3; c++)
    x[c] = __dsub_rn(s.pos[3 * q + c], s.pos[3 * p + c]);
  double rjs = sqrt(__dadd_rn(__dadd_rn(__dmul_rn(x[0], x[0]), __dmul_rn(x[1], x[1])), __dmul_rn(x[2], x[2])));
  double r5 = mbd_pow5(rjs);
  for (int c1 = 0; c1 < 3; c1++)
    for (int c2 = 0; c2 < 3; c2++)
      tm9[3 * c1 + c2] = c1 == c2 ? __dsub_rn(__dmul_rn(__dmul_rn(3.0, x[c1]), x[c1]), __dmul_rn(rjs, rjs)) / r5
                                  : __dmul_rn(__dmul_rn(3.0, x[c1]), x[c2]) / r5;
  return rjs;
}

// Each pair's damping, cutoff weight and frequency-1 values, as the host's
// build forms them: val = sqrt(a_p) sqrt(a_q) f T tlr where |f T tlr| > 1e-20.
// The central atom's pairs are cut at rcut_mbd, the others at rcut_mbd2.
void mbd_pair_values(mbd_site_d s, mbd_consts k, hipStream_t* stream) {
  tg_parallel_for(
      "mbd_pair_values", s.n_pairs, stream, TG_LAMBDA(const int k2) {
        int p = s.atom_of_pair[k2];
        int msk = 0;
        if (k2 != s.own_entry[p] && s.p_mbd[k2] != -1) {
          int q = s.p_mbd[k2] - 1;
          double x[3];
          double tm9[9];
          double rjs = mbd_pair_tensor(s, p, q, x, tm9);
          double rcut = p == 0 ? k.rcut_mbd : k.rcut_mbd2;
          double r_b = __dmul_rn(rjs, k.Bohr);
          double tlr = 1.0;
          if (r_b > rcut - k.r_buf_ij && r_b <= rcut) {
            double rb = __dadd_rn(__dsub_rn(r_b, rcut), k.r_buf_ij) / k.r_buf_ij;
            double rb2 = __dmul_rn(rb, rb);
            tlr = __dadd_rn(__dsub_rn(1.0, __dmul_rn(3.0, rb2)), __dmul_rn(2.0, __dmul_rn(rb2, rb)));
          } else if (r_b > rcut) {
            tlr = 0.0;
          }
          double f = 1.0 / (1.0 + exp(-k.d * (rjs / (k.sR * (s.r0_own[p] + s.r0_nb[q])) - 1.0)));
          double sa = sqrt(s.a_own[p]) * sqrt(s.a_nb[q]);
          double* vb = s.blk + 9 * (size_t) s.bpos[k2];
          double* vd = s.bdup[k2] >= 0 ? s.blk + 9 * (size_t) s.bdup[k2] : nullptr;
          for (int c = 0; c < 9; c++) {
            bool keep = fabs((f * tm9[c]) * tlr) > 1.0e-20;
            vb[c] = keep ? ((sa * f) * tm9[c]) * tlr : 0.0;
            if (vd)
              vd[c] = vb[c];
            if (keep)
              msk |= 1 << c;
          }
          s.f_damp[k2] = f;
          s.t_lr_mult[k2] = tlr;
        }
        s.mask[k2] = msk;
      });
}

// Z(r, c, w) = s(p(r), w) X(r, c, w), stored row by row: the 3 n_freq
// columns of a row are contiguous.
void mbd_scale_rows(mbd_site_d s, const double* s_d, const double* x_d, double* z_d, hipStream_t* stream) {
  int n3 = 3 * s.n_atoms;
  int n = s.n_atoms;
  int ncols = 3 * s.n_freq;
  tg_parallel_for(
      "mbd_scale_rows", n3 * ncols, stream, TG_LAMBDA(const int idx) {
        int gc = idx % ncols;
        int r = idx / ncols;
        int w = gc / 3;
        z_d[idx] = s_d[r / 3 + n * w] * x_d[r + (size_t) n3 * gc];
      });
}

// Columns c0.. of Y = D A_1 Z, one warp per atom p (rows 3p .. 3p+2): lane
// l sums columns c0 + l and c0 + l + 32 over p's blocks in order, each block
// one contiguous read of Z's rows 3q .. 3q+2. Then pol -= k2*res(k2+1)*Y and,
// past the first power, the diagonal kept for the integrand.
__global__ void mbd_spmm_blocks_kernel(mbd_site_d s, int k2, int c0, double pol_coeff, const double* s_d, const double* z_d,
                                       double* y_d, double* pol_d, double* diag_d) {
  int n = s.n_atoms;
  int n3 = 3 * n;
  int ncols = 3 * s.n_freq;
  int lane = threadIdx.x % MBD_LANES;
  int p = blockIdx.x * MBD_ROWS_PER_BLOCK + threadIdx.x / MBD_LANES;
  if (p >= n)
    return;
  double acc[3][MBD_COL_SLOTS];
#pragma unroll
  for (int c1 = 0; c1 < 3; c1++)
#pragma unroll
    for (int m = 0; m < MBD_COL_SLOTS; m++)
      acc[c1][m] = 0.0;
#pragma unroll 2
  for (int b = s.bptr[p]; b < s.bptr[p + 1]; b++) {
    const double* vb = s.blk + 9 * (size_t) b;
    const double* zr = z_d + (size_t) ncols * 3 * s.bcol[b];
    double v9[9];
#pragma unroll
    for (int c = 0; c < 9; c++)
      v9[c] = vb[c];
#pragma unroll
    for (int m = 0; m < MBD_COL_SLOTS; m++) {
      int gc = c0 + lane + MBD_LANES * m;
      if (gc < ncols) {
        double z0 = zr[gc];
        double z1 = zr[ncols + gc];
        double z2 = zr[2 * ncols + gc];
#pragma unroll
        for (int c1 = 0; c1 < 3; c1++) {
          acc[c1][m] = acc[c1][m] + v9[3 * c1] * z0;
          acc[c1][m] = acc[c1][m] + v9[3 * c1 + 1] * z1;
          acc[c1][m] = acc[c1][m] + v9[3 * c1 + 2] * z2;
        }
      }
    }
  }
#pragma unroll
  for (int m = 0; m < MBD_COL_SLOTS; m++) {
    int gc = c0 + lane + MBD_LANES * m;
    if (gc >= ncols)
      continue;
    int w = gc / 3;
    int cc = gc % 3;
#pragma unroll
    for (int c1 = 0; c1 < 3; c1++) {
      int r = 3 * p + c1;
      size_t at = r + (size_t) n3 * gc;
      double y = s_d[p + n * w] * acc[c1][m];
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
  int blocks = (s.n_atoms + MBD_ROWS_PER_BLOCK - 1) / MBD_ROWS_PER_BLOCK;
  for (int c0 = 0; c0 < 3 * s.n_freq; c0 += MBD_LANES * MBD_COL_SLOTS)
    mbd_spmm_blocks_kernel<<<blocks, MBD_ROWS_PER_BLOCK * MBD_LANES, 0, stream[0]>>>(s, k2, c0, pol_coeff, s_d, z_d, y_d, pol_d,
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
        double om = s.omegas[w];
        double o = s.o_own[p];
        double t = om / o;
        double den = o * o + om * om;
        t_d[idx] = da[p] / (1.0 + t * t) + ((s.a_own[p] * ((2.0 * (om * om)) * o)) * dO[p]) / (den * den);
      });
}

// One thread per pair entry: the frequency-independent part once (dT_LR_val,
// and dT_LR_sym for the central atom's pairs), then for each frequency the
// pair's dot, virial and damped virial, reduced over each warp by a shuffle
// tree and over the block's warps in order into part_d (8 per block and w;
// the 8th is the dot's magnitude, kept for tests).
__global__ void mbd_force_kernel(mbd_site_d s, mbd_consts k, int c3, const double* da, const double* dO, const double* dr0,
                                 const double* s_d, const double* q_d, const double* t_d, double* part_d, double* dtlr_sym_d) {
  __shared__ double sh[MBD_FORCE_TPB / MBD_LANES][8];
  int k2 = blockIdx.x * MBD_FORCE_TPB + threadIdx.x;
  int lane = threadIdx.x % MBD_LANES;
  int warp = threadIdx.x / MBD_LANES;
  int n = s.n_atoms;
  int c3m = c3 - 1;
  bool live = k2 < s.n_pairs;
  int p = 0, q = 0, kb = 0;
  if (live) {
    p = s.atom_of_pair[k2];
    kb = s.own_entry[p];
    live = !(k2 == kb || s.p_mbd[k2] == -1);
  }
  // Two atoms whose polarisabilities do not move with atom i0, neither of
  // them i0: the pair adds exactly nothing.
  if (live) {
    int qq = s.p_mbd[k2] - 1;
    bool still_p = da[p] == 0.0 && dO[p] == 0.0 && dr0[p] == 0.0;
    bool still_q = da[qq] == 0.0 && dO[qq] == 0.0 && dr0[qq] == 0.0;
    live = !(still_p && still_q && s.gidx[p] != s.i0 && s.gidx[qq] != s.i0);
  }
  double dtv[9];
  double tm9[9];
  double f = 0.0, tlr = 0.0, pcut_i = 0.0, pcut_next = 0.0;
  double xi[3] = {0.0, 0.0, 0.0}, xn[3] = {0.0, 0.0, 0.0};
  bool dup = false;
  int msk = 0;
  if (live) {
    q = s.p_mbd[k2] - 1;
    msk = s.mask[k2];
    int i2 = s.gidx[p];
    int j = s.gidx[q];
    f = s.f_damp[k2];
    tlr = s.t_lr_mult[k2];
    double x[3];
    double rjs = mbd_pair_tensor(s, p, q, x, tm9);
    double S = k.sR * (s.r0_own[p] + s.r0_nb[q]);
    double dS = k.sR * (dr0[p] + dr0[q]);
    double e = exp(-k.d * (rjs / S - 1.0));
    double fdds = -((((k.d * rjs) / (S * S)) * (f * f)) * e) * dS;
    for (int c = 0; c < 9; c++)
      dtv[c] = 0.0 + (tm9[c] * fdds) * tlr;
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
          double tm = tm9[c];
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
      xi[c] = s.pos[3 * p + c];
    pcut_i = s.pcut[p];
    dup = (k2 + 2 < s.n_pairs) && (s.p_mbd[k2 + 1] == -1);
    if (dup) {
      for (int c = 0; c < 3; c++)
        xn[c] = s.pos[3 * q + c];
      pcut_next = s.pcut[q];
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
      const double* v = s.blk + 9 * (size_t) s.bpos[k2];
      for (int c1 = 0; c1 < 3; c1++)
        for (int c2 = 0; c2 < 3; c2++) {
          double tm = tm9[3 * c1 + c2];
          if (!((msk >> (3 * c1 + c2)) & 1))
            continue;
          double dval = fq_ss * dtv[3 * c1 + c2] + ((fq_gi * f) * tm) * tlr + ((fq_gj * f) * tm) * tlr;
          double val = s_pq * v[3 * c1 + c2];
          for (int c = 0; c < 3; c++)
            acc[1 + c] = acc[1 + c] - ((0.5 * xi[c]) * dval) * val;
          for (int c = 0; c < 3; c++)
            acc[4 + c] = acc[4 + c] - (((0.5 * xi[c]) * dval) * val) * pcut_i;
          acc[0] = acc[0] + val * dval;
          acc[7] = acc[7] + fabs(val * dval);
          if (dup) {
            double val2 = val;
            for (int c = 0; c < 3; c++)
              acc[1 + c] = acc[1 + c] - ((0.5 * xn[c]) * dval) * val2;
            for (int c = 0; c < 3; c++)
              acc[4 + c] = acc[4 + c] - (((0.5 * xn[c]) * dval) * val2) * pcut_next;
            acc[0] = acc[0] + val2 * dval;
            acc[7] = acc[7] + fabs(val2 * dval);
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

// Block sums, one block per out[8*w + c]: thread t sums blocks t, t + TPB, ...
// in that order, then a fixed shared-memory tree.
__global__ void mbd_force_reduce_kernel(int n_blocks, const double* part_d, double* out_d) {
  __shared__ double sh[MBD_REDUCE_TPB];
  int idx = blockIdx.x;
  int c = idx % 8;
  int w = idx / 8;
  double t = 0.0;
  for (int b = threadIdx.x; b < n_blocks; b += MBD_REDUCE_TPB)
    t = t + part_d[8 * (b + n_blocks * w) + c];
  sh[threadIdx.x] = t;
  __syncthreads();
  for (int off = MBD_REDUCE_TPB / 2; off > 0; off /= 2) {
    if (threadIdx.x < off)
      sh[threadIdx.x] = sh[threadIdx.x] + sh[threadIdx.x + off];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    out_d[idx] = sh[0];
}

void mbd_force_reduce(mbd_site_d s, int n_blocks, const double* part_d, double* out_d, hipStream_t* stream) {
  mbd_force_reduce_kernel<<<8 * s.n_freq, MBD_REDUCE_TPB, 0, stream[0]>>>(n_blocks, part_d, out_d);
}

// Byte offsets of a site's arrays in one allocation, each 256-byte aligned.
struct mbd_layout {
  size_t bytes = 0;
  template <class T> size_t add(size_t n) {
    size_t at = bytes;
    bytes += (sizeof(T) * (n > 0 ? n : 1) + 255) / 256 * 256;
    return at;
  }
};

// A host thread's stream and its one device allocation, kept between sites
// and grown, with a quarter to spare, when a site needs more.
struct mbd_ctx {
  hipStream_t stream;
  mbd_site_d s;
  mbd_consts k;
  int n_blocks;
  char* arena_d = nullptr;
  size_t c_arena = 0;
  std::vector<int> own, atom_of, gidx, bptr, bcol, bpos, bdup;
  std::vector<double> a_own, o_own, r0_own, pos, pcut, da_atom, out;
  int *own_d, *atom_of_d, *gidx_d, *p_mbd_d, *bptr_d, *bcol_d, *bpos_d, *bdup_d, *mask_d;
  double *omegas_d, *a_own_d, *o_own_d, *r0_own_d, *a_nb_d, *r0_nb_d, *pos_d, *pcut_d, *f_d, *tlr_d, *blk_d;
  double *t_lr_sym_d, *s_d, *q_d, *t_d, *at_d, *x_d, *y_d, *z_d, *pol_d, *diag_d, *dots_d, *integrand_d;
  double *da_d, *do_d, *dr0_d, *part_d, *out_d, *dtlr_d;
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

template <class T> void up(T* dev, const T* host, size_t n, hipStream_t* stream) {
  if (n)
    gpuErrchk(hipMemcpyAsync(dev, host, sizeof(T) * n, hipMemcpyHostToDevice, stream[0]));
}

// Own entries, row atoms, each atom's own-entry values, and the block rows:
// atom p's pairs in the host's order, and for each of the central atom's
// pairs (0, q) followed by its duplicate, a block (q, 0) holding the same
// values, as the host's duplicate entries do.
void index_pairs(mbd_ctx& c, int n_atoms, int n_pairs, const int* n_neigh, const int* neigh, const int* p_mbd, const double* a_mbd,
                 const double* o_mbd, const double* r0_scs, const double* xyz0, const double* pcut) {
  c.own.assign(n_atoms, 0);
  c.atom_of.assign(n_pairs, 0);
  int k2 = 0;
  for (int p = 0; p < n_atoms; p++) {
    c.own[p] = k2;
    for (int j = 0; j < n_neigh[p]; j++)
      c.atom_of[k2 + j] = p;
    k2 += n_neigh[p];
  }
  c.gidx.resize(n_atoms);
  c.a_own.resize(n_atoms);
  c.o_own.resize(n_atoms);
  c.r0_own.resize(n_atoms);
  c.pos.resize(3 * (size_t) n_atoms);
  c.pcut.resize(n_atoms);
  for (int p = 0; p < n_atoms; p++) {
    int kb = c.own[p];
    c.gidx[p] = neigh[kb];
    c.a_own[p] = a_mbd[kb];
    c.o_own[p] = o_mbd[kb];
    c.r0_own[p] = r0_scs[kb];
    for (int cc = 0; cc < 3; cc++)
      c.pos[3 * p + cc] = xyz0[3 * (size_t) kb + cc];
    c.pcut[p] = pcut[kb];
  }
  c.bptr.assign(n_atoms + 1, 0);
  for (int p = 0; p < n_atoms; p++)
    for (int k = c.own[p] + 1; k < c.own[p] + n_neigh[p]; k++) {
      if (p_mbd[k] == -1)
        continue;
      c.bptr[p + 1]++;
      if (k + 1 < n_pairs && p_mbd[k + 1] == -1)
        c.bptr[p_mbd[k]]++;
    }
  for (int p = 0; p < n_atoms; p++)
    c.bptr[p + 1] += c.bptr[p];
  c.bcol.resize(c.bptr[n_atoms]);
  c.bpos.assign(n_pairs, -1);
  c.bdup.assign(n_pairs, -1);
  std::vector<int> fill(c.bptr.begin(), c.bptr.end() - 1);
  for (int p = 0; p < n_atoms; p++)
    for (int k = c.own[p] + 1; k < c.own[p] + n_neigh[p]; k++) {
      if (p_mbd[k] == -1)
        continue;
      int q = p_mbd[k] - 1;
      c.bcol[fill[p]] = q;
      c.bpos[k] = fill[p]++;
      if (k + 1 < n_pairs && p_mbd[k + 1] == -1) {
        c.bcol[fill[q]] = p;
        c.bdup[k] = fill[q]++;
      }
    }
}

// Every array of a site and of its force passes, in c's one allocation.
void place_arrays(mbd_ctx& c, int n_atoms, int n_pairs, int n_freq, int n_order, hipStream_t* stream) {
  size_t n = n_atoms;
  size_t n3f = 9 * n * n_freq;
  size_t naf = n * n_freq;
  size_t n_blk = c.bcol.size();
  c.n_blocks = (n_pairs + MBD_FORCE_TPB - 1) / MBD_FORCE_TPB;
  mbd_layout L;
  size_t o_own = L.add<int>(n), o_atom_of = L.add<int>(n_pairs), o_gidx = L.add<int>(n);
  size_t o_p_mbd = L.add<int>(n_pairs), o_bptr = L.add<int>(n + 1), o_bcol = L.add<int>(n_blk);
  size_t o_bpos = L.add<int>(n_pairs), o_bdup = L.add<int>(n_pairs), o_mask = L.add<int>(n_pairs);
  size_t o_omegas = L.add<double>(n_freq), o_a_own = L.add<double>(n), o_o_own = L.add<double>(n);
  size_t o_r0_own = L.add<double>(n), o_a_nb = L.add<double>(n), o_r0_nb = L.add<double>(n);
  size_t o_pos = L.add<double>(3 * n), o_pcut = L.add<double>(n), o_f = L.add<double>(n_pairs);
  size_t o_tlr = L.add<double>(n_pairs), o_blk = L.add<double>(9 * n_blk), o_t_lr_sym = L.add<double>(9 * n);
  size_t o_s = L.add<double>(naf), o_q = L.add<double>(naf), o_t = L.add<double>(naf);
  size_t o_at = L.add<double>(n3f), o_x = L.add<double>(n3f), o_y = L.add<double>(n3f), o_z = L.add<double>(n3f);
  size_t o_pol = L.add<double>(n3f), o_diag = L.add<double>(3 * (size_t) (n_order - 3) * n_freq);
  size_t o_dots = L.add<double>(6 * (size_t) n_freq), o_integrand = L.add<double>(n_freq);
  size_t o_da = L.add<double>(3 * n), o_part = L.add<double>(8 * (size_t) c.n_blocks * n_freq);
  size_t o_out = L.add<double>(8 * (size_t) n_freq), o_dtlr = L.add<double>(9 * n);
  if (L.bytes > c.c_arena) {
    if (c.arena_d)
      gpuErrchk(hipFreeAsync(c.arena_d, stream[0]));
    c.c_arena = L.bytes + L.bytes / 4;
    gpuErrchk(hipMallocAsync((void**) &c.arena_d, c.c_arena, stream[0]));
  }
  char* a = c.arena_d;
  c.own_d = (int*) (a + o_own);
  c.atom_of_d = (int*) (a + o_atom_of);
  c.gidx_d = (int*) (a + o_gidx);
  c.p_mbd_d = (int*) (a + o_p_mbd);
  c.bptr_d = (int*) (a + o_bptr);
  c.bcol_d = (int*) (a + o_bcol);
  c.bpos_d = (int*) (a + o_bpos);
  c.bdup_d = (int*) (a + o_bdup);
  c.mask_d = (int*) (a + o_mask);
  c.omegas_d = (double*) (a + o_omegas);
  c.a_own_d = (double*) (a + o_a_own);
  c.o_own_d = (double*) (a + o_o_own);
  c.r0_own_d = (double*) (a + o_r0_own);
  c.a_nb_d = (double*) (a + o_a_nb);
  c.r0_nb_d = (double*) (a + o_r0_nb);
  c.pos_d = (double*) (a + o_pos);
  c.pcut_d = (double*) (a + o_pcut);
  c.f_d = (double*) (a + o_f);
  c.tlr_d = (double*) (a + o_tlr);
  c.blk_d = (double*) (a + o_blk);
  c.t_lr_sym_d = (double*) (a + o_t_lr_sym);
  c.s_d = (double*) (a + o_s);
  c.q_d = (double*) (a + o_q);
  c.t_d = (double*) (a + o_t);
  c.at_d = (double*) (a + o_at);
  c.x_d = (double*) (a + o_x);
  c.y_d = (double*) (a + o_y);
  c.z_d = (double*) (a + o_z);
  c.pol_d = (double*) (a + o_pol);
  c.diag_d = (double*) (a + o_diag);
  c.dots_d = (double*) (a + o_dots);
  c.integrand_d = (double*) (a + o_integrand);
  c.da_d = (double*) (a + o_da);
  c.do_d = c.da_d + n;
  c.dr0_d = c.da_d + 2 * n;
  c.part_d = (double*) (a + o_part);
  c.out_d = (double*) (a + o_out);
  c.dtlr_d = (double*) (a + o_dtlr);
}

} // namespace mbd_gpu

using namespace mbd_gpu;

extern "C" {

void mbd_gpu_init() {
  gpuErrchk(hipGetDevice(&g_device));
}

void mbd_gpu_site(int tid, int n_atoms, int n_pairs, int n_freq, int n_order, int i0, double d, double sR, double Bohr,
                  double rcut_mbd, double rcut_mbd2, double r_buf_ij, const int* n_neigh, const int* neigh, const int* p_mbd,
                  const double* omegas, const double* a_mbd, const double* o_mbd, const double* r0_scs, const double* xyz0,
                  const double* pcut, const double* a_in, const double* r0_in, const double* t_lr_sym, double* integrand,
                  double* pol) {
  mbd_ctx& c = *context(tid);
  gpuErrchk(hipSetDevice(g_device));
  hipStream_t* st = &c.stream;
  index_pairs(c, n_atoms, n_pairs, n_neigh, neigh, p_mbd, a_mbd, o_mbd, r0_scs, xyz0, pcut);
  place_arrays(c, n_atoms, n_pairs, n_freq, n_order, st);
  size_t n3f = (size_t) 9 * n_atoms * n_freq;
  up(c.own_d, c.own.data(), n_atoms, st);
  up(c.atom_of_d, c.atom_of.data(), n_pairs, st);
  up(c.gidx_d, c.gidx.data(), n_atoms, st);
  up(c.p_mbd_d, p_mbd, n_pairs, st);
  up(c.bptr_d, c.bptr.data(), n_atoms + 1, st);
  up(c.bcol_d, c.bcol.data(), c.bcol.size(), st);
  up(c.bpos_d, c.bpos.data(), n_pairs, st);
  up(c.bdup_d, c.bdup.data(), n_pairs, st);
  up(c.omegas_d, omegas, n_freq, st);
  up(c.a_own_d, c.a_own.data(), n_atoms, st);
  up(c.o_own_d, c.o_own.data(), n_atoms, st);
  up(c.r0_own_d, c.r0_own.data(), n_atoms, st);
  up(c.a_nb_d, a_in, n_atoms, st);
  up(c.r0_nb_d, r0_in, n_atoms, st);
  up(c.pos_d, c.pos.data(), 3 * (size_t) n_atoms, st);
  up(c.pcut_d, c.pcut.data(), n_atoms, st);
  up(c.t_lr_sym_d, t_lr_sym, 9 * (size_t) n_atoms, st);

  mbd_site_d& s = c.s;
  s = {n_atoms,   n_pairs,  n_freq,   n_order,  i0,         c.atom_of_d, c.own_d,   c.gidx_d,    c.p_mbd_d,
       c.bptr_d,  c.bcol_d, c.bpos_d, c.bdup_d, c.omegas_d, c.a_own_d,   c.o_own_d, c.r0_own_d,  c.a_nb_d,
       c.r0_nb_d, c.pos_d,  c.pcut_d, c.f_d,    c.tlr_d,    c.blk_d,     c.mask_d,  c.t_lr_sym_d};
  c.k = {d, sR, Bohr, rcut_mbd, rcut_mbd2, r_buf_ij};

  mbd_pair_values(s, c.k, st);
  mbd_freq_factors(s, c.s_d, c.q_d, st);
  mbd_at_sym(s, c.q_d, c.at_d, st);
  gpuErrchk(hipMemsetAsync(c.pol_d, 0, sizeof(double) * n3f, st[0]));
  for (int k2 = 3; k2 <= n_order; k2++) {
    const double* x_in = (k2 == 3) ? c.at_d : c.x_d;
    mbd_series_step(s, k2, c.s_d, x_in, c.z_d, c.y_d, c.pol_d, c.diag_d, st);
    double* t = c.x_d;
    c.x_d = c.y_d;
    c.y_d = t;
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
  size_t n = s.n_atoms;
  c.da_atom.resize(3 * n);
  for (size_t p = 0; p < n; p++) {
    c.da_atom[p] = da[c.own[p]];
    c.da_atom[n + p] = dO[c.own[p]];
    c.da_atom[2 * n + p] = dr0[c.own[p]];
  }
  up(c.da_d, c.da_atom.data(), 3 * n, st);
  gpuErrchk(hipMemsetAsync(c.dtlr_d, 0, sizeof(double) * 9 * n, st[0]));
  mbd_freq_derivs(s, c.da_d, c.do_d, c.t_d, st);
  mbd_force_kernel<<<c.n_blocks, MBD_FORCE_TPB, 0, st[0]>>>(s, c.k, c3, c.da_d, c.do_d, c.dr0_d, c.s_d, c.q_d, c.t_d, c.part_d,
                                                            c.dtlr_d);
  mbd_force_reduce(s, c.n_blocks, c.part_d, c.out_d, st);
  c.out.resize(8 * (size_t) s.n_freq);
  gpuErrchk(hipMemcpyAsync(c.out.data(), c.out_d, sizeof(double) * 8 * s.n_freq, hipMemcpyDeviceToHost, st[0]));
  gpuErrchk(hipMemcpyAsync(dt_lr_sym, c.dtlr_d, sizeof(double) * 9 * n, hipMemcpyDeviceToHost, st[0]));
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
