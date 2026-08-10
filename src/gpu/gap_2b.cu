// The two-body descriptor and the core (repulsive) potential, both of which
// evaluate a cubic spline over pair distances.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64

__device__ double gpu_spline(int nx, double r, double* x, double* y, double* y2, double rcut, double yp1, double ypn) {
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
      double h26 = pow(h, 2) / 6.0;
      double A = (x[j + 1] - r) / h;
      double B = 1.0 - A;
      double C = (pow(A, 3) - A) * h26;
      double D = (pow(B, 3) - B) * h26;
      s = A * y[j] + B * y[j + 1] + C * y2[j] + D * y2[j + 1];
    }
  }
  return s;
}

__device__ double gpu_spline_der(int nx, double r, double* x, double* y, double* y2, double rcut, double yp1, double ypn) {
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
      double dCdx = (1.0 - 3.0 * pow(A, 2)) * h6;
      double dDdx = (3.0 * pow(B, 2) - 1.0) * h6;
      ds = dAdx * y[j] + dBdx * y[j + 1] + dCdx * y2[j] + dDdx * y2[j + 1];
    }
  }
  return ds;
}

__global__ void kernel_get_2b(int i_beg, int i_end, int n_sparse, double* energies_d, double e0, int* n_neigh_d, bool do_forces,
                              double* forces_d, double* virial_d, double* rjs_d, double rcut, int* species_d,
                              int* neighbor_species_d, int sp1, int sp2, double buffer, double delta, double* cutoff_d,
                              double* Qs_d, double sigma, double* alphas_d, double* xyz_d) {
  int i_site = i_beg - 1 + threadIdx.x + blockIdx.x * blockDim.x;
  double forces_loc[3], virial_loc[9], energies_loc;
  int i, j, k, s, k1, k2;
  double rjs_k, fcut, dfcut, pi, sigma2, delta2, tmp;
  if (i_site < i_end) {
    energies_loc = e0;
    if (do_forces) {
      for (i = 0; i < 3; i++)
        forces_loc[i] = 0.0;
      for (i = 0; i < 9; i++)
        virial_loc[i] = 0.0;
    }
    if (species_d[i_site] == sp1 || species_d[i_site] == sp2) {
      pi = acos(-1.0);
      k = 0;
      for (i = i_beg - 1; i < i_site; i++)
        k += n_neigh_d[i];
      for (j = 2; j <= n_neigh_d[i_site]; j++) {
        k += 1;
        if (!((species_d[i_site] == sp1 && neighbor_species_d[k] == sp2) ||
              (species_d[i_site] == sp2 && neighbor_species_d[k] == sp1)))
          continue;
        rjs_k = rjs_d[k];
        if (rjs_k < rcut) {
          fcut = (rjs_k < rcut - buffer) ? 1.0 : (cos(pi * (rjs_k - rcut + buffer) / buffer) + 1.0) / 2.0;
          if (do_forces)
            dfcut = (rjs_k < rcut - buffer) ? 0.0 : pi / 2.0 / buffer * sin(pi * (rjs_k - rcut + buffer) / buffer);
          sigma2 = sigma * sigma;
          delta2 = delta * delta;
          for (s = 0; s < n_sparse; s++) {
            tmp = delta2 * alphas_d[s] * cutoff_d[s] * exp(-0.5 * pow(rjs_k - Qs_d[s], 2) / sigma2);
            energies_loc += tmp * fcut;
            if (do_forces) {
              for (i = 0; i < 3; i++) {
                forces_loc[i] = -2.0 * tmp * xyz_d[3 * k + i] / rjs_k * ((rjs_k - Qs_d[s]) / sigma2 * fcut + dfcut);
                forces_d[3 * i_site + i] += forces_loc[i];
              }
              for (k2 = 0; k2 < 3; k2++)
                for (k1 = 0; k1 < 3; k1++)
                  virial_loc[3 * k2 + k1] += -0.5 * (forces_loc[k1] * xyz_d[3 * k + k2] + forces_loc[k2] * xyz_d[3 * k + k1]);
            }
          }
        }
      }
    }
    // Accumulate, like the force and virial writes above: the caller runs
    // this kernel once per descriptor over the same buffers and reads the
    // totals back after the last one. Assigning here kept only the
    // contribution of the final descriptor.
    energies_d[i_site] += energies_loc;
    for (k1 = 0; k1 < 9; k1++)
      atomicAdd(&virial_d[k1], 0.5 * virial_loc[k1]);
  }
}

extern "C" void gpu_get_2b_forces_energies(int i_beg, int i_end, int n_sparse, double* energies_d, double e0, int* n_neigh_d,
                                           bool do_forces, double* forces_d, double* virial_d, double* rjs_d, double rcut,
                                           int* species_d, int* neighbor_species_d, int sp1, int sp2, double buffer, double delta,
                                           double* cutoff_d, double* Qs_d, double sigma, double* alphas_d, double* xyz_d,
                                           hipStream_t* stream) {
  dim3 nblocks = dim3((i_end - i_beg + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_get_2b<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, n_sparse, energies_d, e0, n_neigh_d, do_forces, forces_d,
                                                     virial_d, rjs_d, rcut, species_d, neighbor_species_d, sp1, sp2, buffer, delta,
                                                     cutoff_d, Qs_d, sigma, alphas_d, xyz_d);
  // There was a hipStreamSynchronize here, commented "temporary, to measure
  // timings". It blocked the host on every launch, so the caller's loop over
  // the 2b descriptors could not have overlapped even with a stream each. The
  // results are read back with a stream sync after that loop
  // (add_2b_contribution), which is where the ordering guarantee belongs.
  //
  // It did make time%gap_2b measure kernel execution rather than launch. That
  // accounting moved with it: the "2b" NVTX range and timing bucket now wrap
  // the whole routine, whose final sync still encloses the kernels.
}

__global__ void kernel_get_core_pot(int i_beg, int i_end, bool do_forces, int* species_d, int sp1, int sp2, int* n_neigh_d,
                                    int* neighbor_species_d, double* rjs_d, int n_sparse, double* x_d, double* V_d, double* dVdx2_d,
                                    double yp1, double ypn, double* xyz_d, double* forces_d, double* virial_d, double* energies_d) {
  int i_site = i_beg - 1 + threadIdx.x + blockIdx.x * blockDim.x;
  double forces_loc[3], virial_loc[9], energies_loc;
  int i, j, k, k1, k2;
  double rjs_k, rcut, Vint, d_Vint;

  if (i_site < i_end) {
    energies_loc = 0.0;
    if (do_forces) {
      for (i = 0; i < 3; i++)
        forces_loc[i] = 0.0;
      for (i = 0; i < 9; i++)
        virial_loc[i] = 0.0;
    }
    if (species_d[i_site] == sp1 || species_d[i_site] == sp2) {
      k = 0;
      for (i = i_beg - 1; i < i_site; i++)
        k += n_neigh_d[i];
      for (j = 1; j < n_neigh_d[i_site]; j++) {
        k += 1;
        if (!((species_d[i_site] == sp1 && neighbor_species_d[k] == sp2) ||
              (species_d[i_site] == sp2 && neighbor_species_d[k] == sp1)))
          continue;
        rjs_k = rjs_d[k];
        rcut = x_d[0];
        for (i = 1; i < n_sparse; i++)
          if (rcut <= x_d[i])
            rcut = x_d[i];

        //	printf("GPU_CORE_POT: rjs_k = %lf, k = %d, rcut = %lf, i_site = %d \n",rjs_k, k, rcut,  i_site);
        if (rjs_k < rcut) {
          energies_loc += 0.5 * gpu_spline(n_sparse, rjs_k, x_d, V_d, dVdx2_d, rcut, yp1, ypn);

          if (do_forces) {
            d_Vint = gpu_spline_der(n_sparse, rjs_k, x_d, V_d, dVdx2_d, rcut, yp1, ypn);
            for (i = 0; i < 3; i++) {
              forces_loc[i] = d_Vint * xyz_d[3 * k + i] / rjs_k;
              forces_d[3 * i_site + i] += forces_loc[i];
            }
            for (k2 = 0; k2 < 3; k2++)
              for (k1 = 0; k1 < 3; k1++)
                virial_loc[3 * k2 + k1] += -0.5 * (forces_loc[k1] * xyz_d[3 * k + k2] + forces_loc[k2] * xyz_d[3 * k + k1]);
          }
        }
      }
    }
    // Accumulate, like the force and virial writes above: the caller runs
    // this kernel once per descriptor over the same buffers and reads the
    // totals back after the last one. Assigning here kept only the
    // contribution of the final descriptor.
    energies_d[i_site] += energies_loc;
    for (k1 = 0; k1 < 9; k1++)
      atomicAdd(&virial_d[k1], 0.5 * virial_loc[k1]);
  }
}

extern "C" void gpu_get_core_pot_energy_and_forces(int i_beg, int i_end, bool do_forces, int* species_d, int sp1, int sp2,
                                                   int* n_neigh_d, int* neighbor_species_d, double* rjs_d, int n_sparse,
                                                   double* x_d, double* V_d, double* dVdx2_d, double yp1, double ypn, double* xyz_d,
                                                   double* forces_d, double* virial_d, double* energies_d, hipStream_t* stream) {
  dim3 nblocks = dim3((i_end - i_beg + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  kernel_get_core_pot<<<nblocks, nthreads, 0, stream[0]>>>(i_beg, i_end, do_forces, species_d, sp1, sp2, n_neigh_d,
                                                           neighbor_species_d, rjs_d, n_sparse, x_d, V_d, dVdx2_d, yp1, ypn, xyz_d,
                                                           forces_d, virial_d, energies_d);
}
