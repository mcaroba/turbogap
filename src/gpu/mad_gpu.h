// The MAD side of the GPU backend: the pair distribution function, the
// structure factor, and the shifted-force electrostatics that shares their
// pair-compaction machinery.

#ifndef TURBOGAP_MAD_GPU_H
#define TURBOGAP_MAD_GPU_H

#include "gpu_common.h"

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

// ---- mad_xrd.cu
void gpu_set_Gk(int nk, int n_samples, int* k_index_d, double* Gk_d, double* pair_distribution_partial_der_d, double c_factor,
                hipStream_t* stream);
void gpu_get_Gka(int i, int n_k, int n_samples, double* Gka_d, double* Gk_d, double* xyz_k_d, hipStream_t* stream);
void gpu_get_Gka_inplace(int i, int n_k, int n_samples, double* Gk_d, double* xyz_k_d, hipStream_t* stream);
void gpu_hadamard_vec_mat_product(int n_samples_sf, int n_k, double* all_scattering_factors_d, double* dermat_d,
                                  hipStream_t* stream);
void gpu_get_fi_dgemv(const int i, const int n_samples_sf, const int n_k, double* dermat_d, double* prefactor_d, double* fi_d,
                      hipblasHandle_t handle, hipStream_t* stream);
void gpu_exp_force_virial_collection(int n_k, int n_sites, double3* forces0, double energy_scale, double* fi, int* j2_list,
                                     double* virial, double3* xyz, hipStream_t* stream);

// ---- mad_electrostatics.cu
void gpu_get_electrostatics_nk(int i_beg, int i_end, int n_pairs, int* n_neigh, int* n_neigh_index_d, double* rjs, double* xyz,
                               double r_cut, int* nk_out_d, int* nk_flags_d, int* nk_flags_sum_d, hipStream_t* stream);
void gpu_set_electrostatics_k_index(int i_beg, int i_end, int n_pairs, int n_sites0, int* neighbors_list, double* rjs, double* xyz,
                                    double* charges_d, double* neighbor_charges_index_d, double* charge_gradients_d,
                                    double* charge_gradients_index_d, int* k_index_d, int* j2_index_d, double* rjs_index_d,
                                    double* xyz_k_d, int* nk_sum_flags_d, hipStream_t* stream);
void gpu_get_electrostatics_energies(const int i_beg, const int nk_max, double* energies_d, double* forces_d, double* virial_d,
                                     int* j2_index_d, const int n_sites, const int this_n_sites, const int this_n_pairs,
                                     int* n_neigh_index_d, double* charges_d, double* charge_gradients_d,
                                     double* neighbor_charges_index_d, double* rjs_index_d, double* xyz_index_d, const double alpha,
                                     const double rcut, const double rcut_in, const double rcut_width, const double B0_rcut,
                                     const double B0_rcut_der, const bool do_damping, const bool do_forces, hipStream_t* stream);

} // extern "C"

#endif // TURBOGAP_MAD_GPU_H
