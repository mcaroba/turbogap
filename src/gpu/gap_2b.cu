// The two-body descriptor and the core (repulsive) potential: every
// descriptor of each in one pass over the pairs.
//
// Flat maps written once through tg_parallel_for and dispatched to CUDA or
// Kokkos by src/gpu/gpu_backend.h.
#include "gap_gpu.h"
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gpu_scan.h"

#define tpb 64

TG_INLINE_FUNCTION double gpu_spline(int nx, double r, double* x, double* y, double* y2, double rcut, double yp1, double ypn) {
  int j;
  double s = 0.0;
  if (r < rcut) {
    if (r < x[0]) {
      s = y[0] + (r - x[0]) * yp1;
    } else if (r > x[nx - 1]) {
      s = y[nx - 1] + (r - x[nx - 1]) * ypn;
    } else {
      for (j = 0; j < nx - 1; j++)
        if (r < x[j + 1])
          break;
      double h = x[j + 1] - x[j];
      // Written as multiplies, not pow(). nvcc does NOT fold pow(x, 2) with a
      // literal exponent -- the PTX for this file carried six real calls to
      // __internal_accurate_pow, which is log2 + exp2 plus correction terms.
      double h26 = h * h / 6.0;
      double A = (x[j + 1] - r) / h;
      double B = 1.0 - A;
      double C = (A * A * A - A) * h26;
      double D = (B * B * B - B) * h26;
      s = A * y[j] + B * y[j + 1] + C * y2[j] + D * y2[j + 1];
    }
  }
  return s;
}

TG_INLINE_FUNCTION double gpu_spline_der(int nx, double r, double* x, double* y, double* y2, double rcut, double yp1, double ypn) {
  int j;
  double ds = 0.0;
  if (r < rcut) {
    if (r < x[0]) {
      ds = yp1;
    } else if (r > x[nx - 1]) {
      ds = ypn;
    } else {
      for (j = 0; j < nx - 1; j++)
        if (r < x[j + 1])
          break;
      double h = x[j + 1] - x[j];
      double h6 = h / 6.0;
      double A = (x[j + 1] - r) / h;
      double B = 1.0 - A;
      double dAdx = -1.0 / h;
      double dBdx = -dAdx;
      double dCdx = (1.0 - 3.0 * A * A) * h6;
      double dDdx = (3.0 * B * B - 1.0) * h6;
      ds = dAdx * y[j] + dBdx * y[j + 1] + dCdx * y2[j] + dDdx * y2[j + 1];
    }
  }
  return ds;
}

// Each site's first pair, counted from site i_beg: one scan, where every thread
// used to sum n_neigh over all the sites before its own, O(sites^2) a call.
static int* site_first_pair(int i_beg, int i_end, int* n_neigh_d, hipStream_t* stream) {
  int n = i_end - i_beg + 1;
  int* first_d;
  gpuErrchk(hipMallocAsync((void**) &first_d, sizeof(int) * n, stream[0]));
  gpuErrchk(hipMemcpyAsync(first_d, n_neigh_d + i_beg - 1, sizeof(int) * n, hipMemcpyDeviceToDevice, stream[0]));
  inclusiveScan(first_d, n, stream);
  tg_parallel_for("turbogap_site_first_pair", n, stream, TG_LAMBDA(const int t) { first_d[t] -= n_neigh_d[i_beg - 1 + t]; }, tpb);
  return first_d;
}

// The pair terms run a thread per pair: a site has tens of neighbours, so a
// thread per site leaves most of a large device idle. Each pair kernel writes
// the pair's energy and force on its site, ep and fp, and sum_to_sites adds
// them up per site in list order, so the totals do not depend on the schedule.

// Each pair's site, as an offset from i_beg.
static int* pair_sites(int i_beg, int n_sites, int n_pairs, int* first_d, int* n_neigh_d, hipStream_t* stream) {
  int* site_d;
  gpuErrchk(hipMallocAsync((void**) &site_d, sizeof(int) * n_pairs, stream[0]));
  tg_parallel_for(
      "turbogap_pair_sites", n_sites, stream,
      TG_LAMBDA(const int t) {
        for (int j = 0; j < n_neigh_d[i_beg - 1 + t]; j++)
          site_d[first_d[t] + j] = t;
      },
      tpb);
  return site_d;
}

// A site's energy, force and virial from its pairs; its first pair is itself.
static void sum_to_sites(int i_beg, int n_sites, int* first_d, int* n_neigh_d, double* ep_d, double* fp_d, double* xyz_d,
                         bool do_forces, double* energies_d, double* forces_d, double* virial_d, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_sum_to_sites", n_sites, stream,
      TG_LAMBDA(const int t) {
        int i_site = i_beg - 1 + t;
        double e = 0.0;
        double f[3] = {0.0, 0.0, 0.0};
        double v[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int j = 1; j < n_neigh_d[i_site]; j++) {
          int k = first_d[t] + j;
          e += ep_d[k];
          if (do_forces) {
            double fl[3] = {fp_d[3 * k], fp_d[3 * k + 1], fp_d[3 * k + 2]};
            for (int i = 0; i < 3; i++)
              f[i] += fl[i];
            for (int k2 = 0; k2 < 3; k2++)
              for (int k1 = 0; k1 < 3; k1++)
                v[3 * k2 + k1] += -0.5 * (fl[k1] * xyz_d[3 * k + k2] + fl[k2] * xyz_d[3 * k + k1]);
          }
        }
        energies_d[i_site] += e;
        if (do_forces) {
          for (int i = 0; i < 3; i++)
            forces_d[3 * i_site + i] += f[i];
          for (int k1 = 0; k1 < 9; k1++)
            TG_ATOMIC_ADD(&virial_d[k1], 0.5 * v[k1]);
        }
      },
      tpb);
}

// Every 2b descriptor in one pass. desc_int holds, per descriptor, the two
// species, n_sparse and the offset of its sparse points in alphas, cutoff and
// Qs; desc_real holds rcut, delta and sigma.
extern "C" void gpu_get_2b_all(int i_beg, int i_end, int n_pairs, int n_desc, int* desc_int_d, double* desc_real_d,
                               double* alphas_d, double* cutoff_d, double* Qs_d, double buffer, int* n_neigh_d, int* species_d,
                               int* neighbor_species_d, double* rjs_d, double* xyz_d, bool do_forces, double* energies_d,
                               double* forces_d, double* virial_d, hipStream_t* stream) {
  if (i_end < i_beg || n_pairs == 0 || n_desc == 0)
    return;
  int n_sites = i_end - i_beg + 1;
  int* first_d = site_first_pair(i_beg, i_end, n_neigh_d, stream);
  int* site_d = pair_sites(i_beg, n_sites, n_pairs, first_d, n_neigh_d, stream);
  double* ep_d;
  double* fp_d;
  gpuErrchk(hipMallocAsync((void**) &ep_d, sizeof(double) * n_pairs, stream[0]));
  gpuErrchk(hipMallocAsync((void**) &fp_d, sizeof(double) * 3 * n_pairs, stream[0]));
  tg_parallel_for(
      "turbogap_get_2b_pairs", n_pairs, stream,
      TG_LAMBDA(const int k) {
        int t = site_d[k];
        double e = 0.0;
        double c = 0.0;
        if (k != first_d[t]) {
          const double pi = acos(-1.0);
          int si = species_d[i_beg - 1 + t];
          int sj = neighbor_species_d[k];
          double r = rjs_d[k];
          for (int d = 0; d < n_desc; d++) {
            const int* di = desc_int_d + 4 * d;
            if (!((si == di[0] && sj == di[1]) || (si == di[1] && sj == di[0])))
              continue;
            double rcut = desc_real_d[3 * d];
            if (r >= rcut)
              continue;
            double delta = desc_real_d[3 * d + 1];
            double sigma = desc_real_d[3 * d + 2];
            double fcut = 1.0;
            double dfcut = 0.0;
            if (r >= rcut - buffer) {
              fcut = (cos(pi * (r - rcut + buffer) / buffer) + 1.0) / 2.0;
              dfcut = pi / 2.0 / buffer * sin(pi * (r - rcut + buffer) / buffer);
            }
            double sigma2 = sigma * sigma;
            double delta2 = delta * delta;
            double es = 0.0;
            double gs = 0.0;
            for (int s = di[3]; s < di[3] + di[2]; s++) {
              double dq = r - Qs_d[s];
              double tmp = delta2 * alphas_d[s] * cutoff_d[s] * exp(-0.5 * dq * dq / sigma2);
              es += tmp;
              gs += tmp * (dq / sigma2 * fcut + dfcut);
            }
            e += es * fcut;
            c += -2.0 * gs / r;
          }
        }
        ep_d[k] = e;
        if (do_forces)
          for (int i = 0; i < 3; i++)
            fp_d[3 * k + i] = c * xyz_d[3 * k + i];
      },
      tpb);
  sum_to_sites(i_beg, n_sites, first_d, n_neigh_d, ep_d, fp_d, xyz_d, do_forces, energies_d, forces_d, virial_d, stream);
  gpuErrchk(hipFreeAsync(fp_d, stream[0]));
  gpuErrchk(hipFreeAsync(ep_d, stream[0]));
  gpuErrchk(hipFreeAsync(site_d, stream[0]));
  gpuErrchk(hipFreeAsync(first_d, stream[0]));
}

// Every core-potential table in one pass. desc_int holds, per table, the two
// species, its length and its offset in x, V and dVdx2; desc_real holds yp1,
// ypn and the table's last x, its cutoff.
extern "C" void gpu_get_core_pot_all(int i_beg, int i_end, int n_pairs, int n_desc, int* desc_int_d, double* desc_real_d,
                                     double* x_d, double* V_d, double* dVdx2_d, int* n_neigh_d, int* species_d,
                                     int* neighbor_species_d, double* rjs_d, double* xyz_d, bool do_forces, double* energies_d,
                                     double* forces_d, double* virial_d, hipStream_t* stream) {
  if (i_end < i_beg || n_pairs == 0 || n_desc == 0)
    return;
  int n_sites = i_end - i_beg + 1;
  int* first_d = site_first_pair(i_beg, i_end, n_neigh_d, stream);
  int* site_d = pair_sites(i_beg, n_sites, n_pairs, first_d, n_neigh_d, stream);
  double* ep_d;
  double* fp_d;
  gpuErrchk(hipMallocAsync((void**) &ep_d, sizeof(double) * n_pairs, stream[0]));
  gpuErrchk(hipMallocAsync((void**) &fp_d, sizeof(double) * 3 * n_pairs, stream[0]));
  tg_parallel_for(
      "turbogap_get_core_pot_pairs", n_pairs, stream,
      TG_LAMBDA(const int k) {
        int t = site_d[k];
        double e = 0.0;
        double c = 0.0;
        if (k != first_d[t]) {
          int si = species_d[i_beg - 1 + t];
          int sj = neighbor_species_d[k];
          double r = rjs_d[k];
          for (int d = 0; d < n_desc; d++) {
            const int* di = desc_int_d + 4 * d;
            if (!((si == di[0] && sj == di[1]) || (si == di[1] && sj == di[0])))
              continue;
            double yp1 = desc_real_d[3 * d];
            double ypn = desc_real_d[3 * d + 1];
            double rcut = desc_real_d[3 * d + 2];
            if (r >= rcut)
              continue;
            int n = di[2];
            int off = di[3];
            e += 0.5 * gpu_spline(n, r, x_d + off, V_d + off, dVdx2_d + off, rcut, yp1, ypn);
            if (do_forces)
              c += gpu_spline_der(n, r, x_d + off, V_d + off, dVdx2_d + off, rcut, yp1, ypn) / r;
          }
        }
        ep_d[k] = e;
        if (do_forces)
          for (int i = 0; i < 3; i++)
            fp_d[3 * k + i] = c * xyz_d[3 * k + i];
      },
      tpb);
  sum_to_sites(i_beg, n_sites, first_d, n_neigh_d, ep_d, fp_d, xyz_d, do_forces, energies_d, forces_d, virial_d, stream);
  gpuErrchk(hipFreeAsync(fp_d, stream[0]));
  gpuErrchk(hipFreeAsync(ep_d, stream[0]));
  gpuErrchk(hipFreeAsync(site_d, stream[0]));
  gpuErrchk(hipFreeAsync(first_d, stream[0]));
}
