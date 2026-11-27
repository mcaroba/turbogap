// The MAD side of the GPU backend: the pair distribution function, the
// structure factor, and the shifted-force electrostatics that shares their
// pair-compaction machinery.

#ifndef TURBOGAP_MAD_GPU_H
#define TURBOGAP_MAD_GPU_H

#include "gpu_backend.h"
#include "gpu_common.h"

// The windowed pdf and its forces (mad_pdf.cu, mad_xrd.cu) evaluate a pair's
// Gaussian at the samples the host evaluates it at -- the index arithmetic of
// pair_kde_row, 0-based -- so the device sums the same terms as the host.

// The channel of the unordered species pair (a, b), species 1-based, in the
// order the Fortran visits them: (1,1), (1,2), ..., (1,n), (2,2), ...
TG_INLINE_FUNCTION int pdf_channel(int a, int b, int n_species) {
  if (a > b) {
    const int t = a;
    a = b;
    b = t;
  }
  return (a - 1) * n_species - ((a - 1) * (a - 2)) / 2 + (b - a);
}

// First and last sample of a pair's window, 0-based.
TG_INLINE_FUNCTION int pdf_window_lo(double r, double r_min, double dr, double kde_sigma, double kde_cut) {
  const int lo = (int) ((r - kde_cut * kde_sigma - r_min) / dr) - 1;
  return lo > 0 ? lo : 0;
}

TG_INLINE_FUNCTION int pdf_window_hi(double r, double r_min, double dr, double kde_sigma, double kde_cut, int n_samples) {
  const int hi = (int) ((r + kde_cut * kde_sigma - r_min) / dr) + 1;
  return hi < n_samples - 1 ? hi : n_samples - 1;
}

extern "C" {

// ---- mad_pdf.cu
void gpu_get_pair_distribution_nk(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, int* n_neigh,
                                  int* neighbor_species, int* species, double* rjs, double* xyz, double r_min, double r_max,
                                  double r_cut, double buffer, int* nk_out_d, int* nk_flags_d, int* nk_flags_sum_d, int species_1,
                                  int species_2, hipStream_t* stream);
void gpu_set_pair_distribution_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, double* rjs,
                                       double* xyz, int* k_index_d, int* j2_index_d, double* rjs_index_d, double* xyz_k_d,
                                       int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_set_pair_distribution_k_index_only(int n_pairs, int* k_index_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_set_pair_distribution_j2_only(int n_pairs, int n_sites0, int* neighbors_list_d, int* j2_index_d, int* nk_sum_flags_d,
                                       hipStream_t* stream);
void gpu_set_pair_distribution_rjs_only(int n_pairs, double* rjs, double* rjs_index_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_set_pair_distribution_xyz_only(int n_pairs, double* xyz, double* xyz_index_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_get_pair_distribution_and_ders(double* pair_distribution_d, double* pair_distribution_der_d, int n_k, int n_samples,
                                        double kde_sigma, double* x_d, double* dV_d, double* rjs_d, double pdf_factor,
                                        double der_factor, hipStream_t* stream);
void gpu_get_pair_distribution_der_only(double* pair_distribution_der_d, int n_k, int n_samples, double kde_sigma, double* x_d,
                                        double* dV_d, double* rjs_d, double pdf_factor, double der_factor, hipStream_t* stream);
void gpu_get_pair_distribution_only(double* pair_distribution_d, int n_k, int n_samples, double kde_sigma, double* x_d,
                                    double* dV_d, double* rjs_d, double pdf_factor, double der_factor, hipStream_t* stream);
void gpu_get_pair_distribution_only_falloc(double* pair_distribution_d, double* pdf_to_reduce, int n_k, int n_samples,
                                           double kde_sigma, double* x_d, double* dV_d, double* rjs_d, double pdf_factor,
                                           double der_factor, hipStream_t* stream);
void gpu_pdf_tag_pairs(int n_sites_b, int n_pairs_b, int* n_neigh_d, int* species_d, int* neighbor_species_d, double* rjs_d,
                       int n_species, double r_lo, double r_hi, int* chan_d, hipStream_t* stream);
void gpu_pdf_compact_pairs(int n_pairs_b, int n_sites0, int n_channels, int* chan_d, int* neighbors_list_d, double* rjs_d,
                           double* xyz_d, double* r_out, double* xyz_out, int* j2_out, int* chan_out, int* off_h, int n_sites_b,
                           int* n_neigh_d, int* seg_out, hipStream_t* stream);
void* gpu_pdf_cells_build(int n, double* pos_h, int* species_h, double lx, double ly, double lz, double rcut, hipStream_t* stream);
void gpu_pdf_cells_count(void* h, int site_beg, int n_site, int* counts_h, hipStream_t* stream);
void gpu_pdf_cells_fill(void* h, int site_beg, int n_site, int* n_neigh_h, int n_pairs, void** n_neigh_out, void** species_out,
                        void** nl_out, void** nspec_out, void** rjs_out, void** xyz_out, hipStream_t* stream);
void gpu_pdf_cells_free(void* h, hipStream_t* stream);
void gpu_pdf_window_kde(int n_channels, int* off_h, int n_samples, double r_min, double dr, double kde_sigma, double kde_cut,
                        double* x_d, double* dV_d, double* r_d, double pdf_factor, double* pdf_out_d, hipStream_t* stream);

// ---- mad_xrd.cu
void gpu_xrd_pair_forces(int nk, int n_samples, double* w_d, double* pdf_der_d, double c_factor, double* xyz_k_d, double* fi_d,
                         hipStream_t* stream);
void gpu_xrd_weights(int n_samples_sf, int n_samples, double* S_d, double* residual_d, double* f_d, double* w_d,
                     hipStream_t* stream);
void gpu_exp_force_virial_collection(int n_k, int n_sites, double3* forces0, double energy_scale, double* fi, int* j2_list,
                                     double* virial, double3* xyz, hipStream_t* stream);
void gpu_xrd_window_weights(int n_samples_sf, int n_samples, int channel, double* S_d, double* residual_d, double* f_d,
                            double* dV_d, double fac, double* w_d, hipStream_t* stream);
void gpu_xrd_window_forces(int n_k, int n_sites, int n_samples, double r_min, double dr, double kde_sigma, double kde_cut,
                           double* x_d, double* w_d, double* r_d, double* xyz_d, int* chan_d, int n_channels, int n_sites_b,
                           int i_off, int* seg_d, double* forces_d, double* virial_d, hipStream_t* stream);

// ---- mad_electrostatics.cu
void gpu_get_electrostatics_nk(int i_beg, int i_end, int n_pairs, int* n_neigh, int* n_neigh_index_d, double* rjs, double* xyz,
                               double r_cut, int* nk_out_d, int* nk_flags_d, int* nk_flags_sum_d, hipStream_t* stream);
void gpu_set_electrostatics_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, double* rjs, double* xyz,
                                    double* charges_d, double* neighbor_charges_index_d, double* charge_gradients_d,
                                    double* charge_gradients_index_d, int* k_index_d, int* j2_index_d, double* rjs_index_d,
                                    double* xyz_k_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_get_electrostatics_energies(const int i_beg, const int nk_max, double* energies_d, double* forces_d, double* virial_d,
                                     double* phi_sum_d, int* j2_index_d, const int n_sites, const int this_n_sites,
                                     const int this_n_pairs, int* n_neigh_index_d, double* charges_d, double* charge_gradients_d,
                                     double* neighbor_charges_index_d, double* rjs_index_d, double* xyz_index_d, const double alpha,
                                     const double rcut, const double rcut_in, const double rcut_width, const double B0_rcut,
                                     const double B0_rcut_der, const bool do_damping, const bool do_forces, hipStream_t* stream);

} // extern "C"

#endif // TURBOGAP_MAD_GPU_H
