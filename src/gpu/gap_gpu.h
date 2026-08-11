// The GAP side of the GPU backend: 2b, 3b and soap_turbo.
//
// Grouped by descriptor, in the order a prediction uses them --
//   radial + angular -> cnk -> descriptor -> forces
// with 2b, 3b and the prediction arithmetic alongside.

#ifndef TURBOGAP_GAP_GPU_H
#define TURBOGAP_GAP_GPU_H

#include "gpu_common.h"

extern "C" {

// ---- gap_predict.cu
void gpu_kernels_pow(double* a, double* b, double zeta, int size, hipStream_t* stream);
void gpu_axpc(double* a, double dccc, double e0, int size, hipStream_t* stream);
void cuda_matvect_kernels(double* kernels_d, double* alphas_d, int n_sites, int n_sparse, hipStream_t* stream);
void cuda_matvect_qs(double* qs_d, double* qs_copy_d, double* alphas_d, int n_soap, int n_sparse, hipStream_t* stream);

// ---- gap_soap_radial.cu
void gpu_get_radial_exp_coeff_poly3gauss(double* radial_exp_coeff_d, double* radial_exp_coeff_der_d, int* i_beg_d, int* i_end_d,
                                         double* global_scaling_d, int size_radial_exp_coeff_one, int size_radial_exp_coeff_two,
                                         int n_species, bool c_do_derivatives, int bintybint, double* rcut_hard_d, int* k2_i_site_d,
                                         int* k_2start_d, hipStream_t* stream);
void gpu_radial_poly3gauss(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d, int n_sites,
                           int* n_neigh_d, int n_max, int n_temp, bool do_derivatives, double* exp_coeff_d, double* exp_coeff_der_d,
                           double* rcut_soft_d, double* atom_sigma_d, double* exp_coeff_temp1_d, double* exp_coeff_temp2_d,
                           double* exp_coeff_der_temp_d, int* i_beg, int* i_end, double* atom_sigma_scaling_d, int mode,
                           int radial_enhancement, double* amplitude_scaling_d, int* alpha_max_d, double* nf_d, int n_temp_der,
                           double* W_d, hipStream_t* stream);
void gpu_radial_poly3(int n_atom_pairs, int n_species, bool* mask_d, double* rjs_d, double* rcut_hard_d, int n_sites,
                      int* n_neigh_d, int n_max, int n_temp, bool do_derivatives, double* exp_coeff_d, double* exp_coeff_der_d,
                      double* rcut_soft_d, double* atom_sigma_d, double* exp_coeff_temp1_d, double* exp_coeff_temp2_d,
                      double* exp_coeff_der_temp_d, int* i_beg, int* i_end, double* atom_sigma_scaling_d, int mode,
                      int radial_enhancement, double* amplitude_scaling_d, int* alpha_max_d, double* nf_d, int n_temp_der,
                      double* W_d, bool* do_central_d, double* central_weight_d, hipStream_t* stream);

// ---- gap_soap_angular.cu
void gpu_get_cnk(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d, hipDoubleComplex* cnk_d, int* n_neigh_d,
                 int* k2_start_d, int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max, int bintybint,
                 double* atom_sigma_r_d, double* atom_sigma_t_d, double* rcut_hard_d, double* central_weight_d, int* species_d,
                 int* i_beg_d, int* i_end_d, int radial_enhancement, int* species_multiplicity_d, double* W_d, double* S_d,
                 int size_species_1, hipStream_t* stream);
void gpu_get_plm_array_global(double* plm_array_global_d, int n_atom_pairs, int kmax, int lmax, double* thetas_d,
                              hipStream_t* stream);
void gpu_get_exp_coeff_array(hipDoubleComplex* eimphi_global_d, double* rjs_d, double* phis_d, double* thetas_d, bool* mask_d,
                             double* atom_sigma_in_d, double* atom_sigma_scaling_d, double rcut, int n_atom_pairs, int n_species,
                             int lmax, int kmax, double* prefl_global_d, double* plm_array_global_d, double* plm_array_global_der_d,
                             double* prefl_global_der_d, double* preflm_d, hipDoubleComplex* exp_coeff_d, bool c_do_derivatives,
                             hipDoubleComplex* eimphi_rad_der_global_d, hipDoubleComplex* eimphi_azi_der_global_d,
                             double* plm_array_div_sin, double* plm_array_der_mul_sin, hipDoubleComplex* exp_coeff_rad_der_d,
                             hipDoubleComplex* exp_coeff_azi_der_d, hipDoubleComplex* exp_coeff_pol_der_d, hipStream_t* stream);

// ---- gap_soap_descriptor.cu
void gpu_get_sqrt_dot_p(double* sqrt_dot_d, double* soap_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                        bool* skip_soap_component_d, int n_sites, int n_soap, int n_max, int l_max, hipStream_t* stream);
void gpu_get_soap_der(double* soap_d, double* sqrt_dot_d, double3* soap_cart_der_d, double* soap_rad_der_d, double* soap_azi_der_d,
                      double* soap_pol_der_d, double* thetas_d, double* phis_d, double* rjs_d, double* multiplicity_array_d,
                      hipDoubleComplex* cnk_d, hipDoubleComplex* cnk_rad_der_d, hipDoubleComplex* cnk_azi_der_d,
                      hipDoubleComplex* cnk_pol_der_d, int* n_neigh_d, int* i_k2_start_d, int* k2_i_site_d, int* k3_index_d,
                      bool* skip_soap_component_d, int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max,
                      int maxneigh, hipStream_t* stream);
void gpu_soap_normalize(double* soap_d, double* sqrt_dot_d, int n_soap, int n_sites, hipStream_t* stream);
void gpu_get_derivatives(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d, double* radial_exp_coeff_der_d,
                         hipDoubleComplex* angular_exp_coeff_rad_der_d, hipDoubleComplex* angular_exp_coeff_azi_der_d,
                         hipDoubleComplex* angular_exp_coeff_pol_der_d, hipDoubleComplex* cnk_rad_der_d,
                         hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, double* rjs_d, double rcut_max,
                         int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max, int l_max, hipStream_t* stream);

// ---- gap_soap_forces.cu
void gpu_final_soap_forces_virial(int n_sites, double* Qss_d, int n_soap, int* l_index_d, int* j2_index_d, double3* soap_der_d,
                                  double3* xyz_d, double* virial_d, int n_sites0, double* forces_d, int n_pairs,
                                  hipStream_t* stream);
void gpu_local_property_derivatives(int n_sites, double* Qss_d, int n_soap, int* l_index_d, double3* soap_der_d,
                                    double* local_property_cart_der_d, int n_pairs, hipStream_t* stream);
void gpu_soap_dipole(int n_sites, double* Qss_d, int n_soap, int* beg_index_d, double3* soap_der_d, double* dipoles_d,
                     hipStream_t* stream);

// ---- gap_2b.cu
void gpu_get_2b_forces_energies(int i_beg, int i_end, int n_sparse, double* energies_d, double e0, int* n_neigh_d, bool do_forces,
                                double* forces_d, double* virial_d, double* rjs_d, double rcut, int* species_d,
                                int* neighbor_species_d, int sp1, int sp2, double buffer, double delta, double* cutoff_d,
                                double* Qs_d, double sigma, double* alphas_d, double* xyz_d, hipStream_t* stream);
void gpu_get_core_pot_energy_and_forces(int i_beg, int i_end, bool do_forces, int* species_d, int sp1, int sp2, int* n_neigh_d,
                                        int* neighbor_species_d, double* rjs_d, int n_sparse, double* x_d, double* V_d,
                                        double* dVdx2_d, double yp1, double ypn, double* xyz_d, double* forces_d, double* virial_d,
                                        double* energies_d, hipStream_t* stream);

// ---- gap_3b.cc
void gpu_3b(const int n_sparse, const int n_sites, const int n_atom_pairs, const int n_sites0, const int sp0, const int sp1,
            const int sp2, const double* alpha, const double delta, const double e0, const double* cutoff, const hipStream_t* s,
            const double* rjs, const double* xyz, const int* n_neigh, const int* species, const int* neighbors_list,
            const int* neighbor_species, const bool do_forces, const double rcut, const double buffer, const double* sigma,
            const double* qs, const char* str, const int i_beg,
            const int i_end, //not sure is really needed though. need to check with nsites.
            double* energy, double* forces, double* virials, const int* kappas_array_d);

} // extern "C"

#endif // TURBOGAP_GAP_GPU_H
