
! module load gcc/10.3.0 openblas openmpi cuda
! rm cuda_wrappers.o ; nvcc -arch=sm_70 -c src/cuda_wrappers.cu ; make clean; make -B
!  srun  --time=00:10:00 --partition=gputest --account=project_2000634 --nodes=1 --ntasks-per-node=4  --cpus-per-task=1 --gres=gpu:a100:4  ../turbogap_dev/bin/turbogap md
! gnupot
!plot "< paste trajectory_out.xyz_64k_oiriginal trajectory_out.xyz | awk 'NR>2{print $5,$13}'", x
MODULE F_B_C

   INTERFACE
      subroutine gpu_malloc_async(a_d, n, gpu_stream) bind(C, name="gpu_malloc_async")
         use iso_c_binding
         implicit none
         type(c_ptr) :: a_d
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
      end subroutine

      subroutine gpu_memset_async(a_d, valuetoset, n, gpu_stream) bind(C, name="cuda_memset_async")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
         integer(c_int), value :: valuetoset
      end subroutine

      subroutine gpu_malloc_all_blocking(a_d, n) bind(C, name="cuda_malloc_all_blocking")
         use iso_c_binding
         implicit none
         type(c_ptr) :: a_d
         integer(c_size_t), value :: n
      end subroutine

      subroutine gpu_device_reset() bind(C, name="cuda_device_reset")
         use iso_c_binding
         implicit none
      end subroutine
      subroutine gpu_free(a_d) bind(C, name="cuda_free")
         use iso_c_binding
         implicit none
         type(c_ptr) :: a_d
      end subroutine

      subroutine gpu_free_async(a_d, gpu_stream) bind(C, name="cuda_free_async")
         use iso_c_binding
         implicit none
         type(c_ptr) :: a_d
         type(c_ptr) :: gpu_stream
      end subroutine

      subroutine cpy_htod(a, a_d, n, gpu_stream) bind(C, name="cuda_cpy_htod")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: a
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
      end subroutine

      subroutine cpy_htod_blocking(a, a_d, n) bind(C, name="cuda_cpy_htod_blocking")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: a
         integer(c_size_t), value :: n
      end subroutine

      subroutine cpy_dtod(b_d, a_d, n, gpu_stream) bind(C, name="cuda_cpy_dtod")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: b_d
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
      end subroutine

      subroutine cpy_dtoh(a_d, a, n, gpu_stream) bind(C, name="cuda_cpy_dtoh")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: a
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
      end subroutine

      subroutine cpy_dtoh_event(a_d, a, n, gpu_stream) bind(C, name="cuda_cpy_dtoh_event")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: a
         type(c_ptr) :: gpu_stream
         integer(c_size_t), value :: n
      end subroutine

      subroutine cpy_dtoh_blocking(a_d, a, n) bind(C, name="cuda_cpy_dtoh_blocking")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         type(c_ptr), value :: a
         integer(c_size_t), value :: n
      end subroutine
      subroutine gpu_vector_fill_curand(a_d, n, c) bind(C, name="GPU_fill_rand")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: a_d
         integer(c_int), value :: n
         integer(c_int), value :: c
      end subroutine

      subroutine create_cublas_handle(cubhandle, gpu_stream) bind(C, name="create_cublas_handle")
         use iso_c_binding
         implicit none
         type(c_ptr) :: cubhandle
         type(c_ptr) :: gpu_stream
      end subroutine

      subroutine destroy_cublas_handle(cubhandle, gpu_stream) bind(C, name="destroy_cublas_handle")
         use iso_c_binding
         implicit none
         type(c_ptr) :: cubhandle
         type(c_ptr) :: gpu_stream
      end subroutine

      subroutine gpu_dgemm_n_n(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, cubhandle) bind(C, name="gpu_dgemm_n_n")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: cubhandle
         integer(c_int), value :: m
         integer(c_int), value :: n
         integer(c_int), value :: k
         integer(c_int), value :: lda
         integer(c_int), value :: ldb
         integer(c_int), value :: ldc
         real(c_double), value :: alpha
         real(c_double), value :: beta
         type(c_ptr), value :: A
         type(c_ptr), value :: B
         type(c_ptr), value :: C
         ! real(c_double), device, dimension(lda,*) :: A
         ! real(c_double), device, dimension(ldb,*) :: B
         ! real(c_double), device, dimension(ldc,*) :: C
      end subroutine

      subroutine gpu_blas_mmul_t_n(cubhandle, Qs_d, soap_d, kernels_d, &
                                   n_sparse, n_soap, n_sites) bind(C, name="gpu_blas_mmul_t_n")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: cubhandle
         type(c_ptr), value :: kernels_d
         type(c_ptr), value :: Qs_d
         type(c_ptr), value :: soap_d
         integer(c_int), value :: n_sparse
         integer(c_int), value :: n_soap
         integer(c_int), value :: n_sites
      end subroutine

      subroutine gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                                   n_soap, n_sites, cdelta) bind(C, name="gpu_blas_mmul_n_t")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: cubhandle
         type(c_ptr), value :: kernels_der_d
         type(c_ptr), value :: Qs_copy_d
         type(c_ptr), value :: Qss_d
         integer(c_int), value :: n_sparse
         integer(c_int), value :: n_soap
         integer(c_int), value :: n_sites
         real(c_double), value :: cdelta
      end subroutine

      subroutine gpu_blas_mvmul_n(cubhandle, A, B, C, nAx, nAy) bind(C, name="gpu_blas_mvmul_n")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: cubhandle
         type(c_ptr), value :: A
         type(c_ptr), value :: B
         type(c_ptr), value :: C
         integer(c_int), value :: nAx
         integer(c_int), value :: nAy
      end subroutine

      subroutine gpu_kernels_pow(A, B, zeta, N, gpu_stream) bind(C, name="gpu_kernels_pow")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: A
         type(c_ptr), value :: B
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: N
         real(c_double), value :: zeta
      end subroutine

      subroutine gpu_axpe(A, dccc, e0, N, gpu_stream) bind(C, name="gpu_axpc")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: A
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: N
         real(c_double), value :: dccc
         real(c_double), value :: e0
      end subroutine

      subroutine gpu_set_device(my_rank) bind(C, name="cuda_set_device")
         use iso_c_binding
         integer(c_int), value :: my_rank
      end subroutine

      subroutine gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse, gpu_stream) bind(C, name="cuda_matvect_kernels")
         use iso_c_binding
         type(c_ptr), value :: kernels_der_d
         type(c_ptr), value :: alphas_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_sparse
      end subroutine

      subroutine gpu_final_soap_forces_virial(n_sites, &
                                              Qss_d, n_soap, l_index_d, j2_index_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              n_sites0, &
                                              forces_d, &
                                              n_pairs, gpu_stream) &
         bind(C, name="gpu_final_soap_forces_virial")
         use iso_c_binding
         type(c_ptr), value :: Qss_d
         type(c_ptr), value :: l_index_d
         type(c_ptr), value :: forces_d
         type(c_ptr), value :: soap_der_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: virial_d
         type(c_ptr), value :: j2_index_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_soap
         integer(c_int), value :: n_sites0
         integer(c_int), value :: n_pairs
      end subroutine

      subroutine gpu_local_property_derivatives(n_sites, &
                                                Qss_d, n_soap, l_index_d, &
                                                soap_der_d, &
                                                local_property_cart_der_d, &
                                                n_pairs, gpu_stream) &
         bind(C, name="gpu_local_property_derivatives")
         use iso_c_binding
         type(c_ptr), value :: Qss_d
         type(c_ptr), value :: l_index_d
         type(c_ptr), value :: local_property_cart_der_d
         type(c_ptr), value :: soap_der_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_soap
         integer(c_int), value :: n_pairs
      end subroutine

      subroutine gpu_get_sqrt_dot_p(sqrt_dot_d, soap_d, multiplicity_array_d, &
                                    cnk_d, skip_soap_component_d, &
                                    n_sites, n_soap, n_max, l_max, gpu_stream) &
         bind(C, name="gpu_get_sqrt_dot_p")
         use iso_c_binding
         type(c_ptr), value :: cnk_d
         type(c_ptr), value :: skip_soap_component_d
         type(c_ptr) :: gpu_stream
         type(c_ptr), value :: sqrt_dot_d
         type(c_ptr), value :: soap_d
         type(c_ptr), value :: multiplicity_array_d
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_soap
         integer(c_int), value :: n_max
         integer(c_int), value :: l_max
      end subroutine

      subroutine gpu_get_soap_der(soap_d, sqrt_dot_d, soap_cart_der_d, &
                                  soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, &
                                  thetas_d, phis_d, rjs_d, &
                                  multiplicity_array_d, &
                                  cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                  n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d, &
                                  n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh, gpu_stream) &
         bind(C, name="gpu_get_soap_der")
         use iso_c_binding
         type(c_ptr), value :: sqrt_dot_d
         type(c_ptr), value :: soap_d
         type(c_ptr), value :: soap_cart_der_d
         type(c_ptr) :: gpu_stream
         type(c_ptr), value :: soap_rad_der_d
         type(c_ptr), value :: soap_azi_der_d
         type(c_ptr), value :: soap_pol_der_d
         type(c_ptr), value :: thetas_d
         type(c_ptr), value :: phis_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: i_k2_start_d
         type(c_ptr), value :: k2_i_site_d
         type(c_ptr), value :: k3_index_d
         type(c_ptr), value :: skip_soap_component_d
         type(c_ptr), value :: cnk_d
         type(c_ptr), value :: cnk_rad_der_d
         type(c_ptr), value :: cnk_azi_der_d
         type(c_ptr), value :: cnk_pol_der_d
         type(c_ptr), value :: multiplicity_array_d
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_soap
         integer(c_int), value :: k_max
         integer(c_int), value :: n_max
         integer(c_int), value :: l_max
         integer(c_int), value :: maxneigh
      end subroutine

      subroutine gpu_soap_normalize(soap_d, sqrt_dot_d, n_soap, n_sites, gpu_stream) &
         bind(C, name="gpu_soap_normalize")
         use iso_c_binding
         type(c_ptr), value :: sqrt_dot_d
         type(c_ptr), value :: soap_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_soap
      end subroutine

      subroutine gpu_get_derivatives(radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d, &
                                     angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d, &
                                     cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                     rjs_d, &
                                     rcut_max, &
                                     n_atom_pairs, n_sites, n_soap, k_max, n_max, l_max, gpu_stream) &
         bind(C, name="gpu_get_derivatives")
         use iso_c_binding
         type(c_ptr), value :: rjs_d
         type(c_ptr) :: gpu_stream
         type(c_ptr), value :: radial_exp_coeff_d
         type(c_ptr), value :: angular_exp_coeff_d
         type(c_ptr), value :: radial_exp_coeff_der_d
         type(c_ptr), value :: angular_exp_coeff_rad_der_d
         type(c_ptr), value :: angular_exp_coeff_azi_der_d
         type(c_ptr), value :: angular_exp_coeff_pol_der_d
         type(c_ptr), value :: cnk_rad_der_d
         type(c_ptr), value :: cnk_azi_der_d
         type(c_ptr), value :: cnk_pol_der_d
         real(c_double), value :: rcut_max
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_soap
         integer(c_int), value :: k_max
         integer(c_int), value :: n_max
         integer(c_int), value :: l_max
      end subroutine

      subroutine gpu_get_cnk(radial_exp_coeff_d, angular_exp_coeff_d, &
                             cnk_d, &
                             n_neigh_d, k2_start_d, &
                             n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, &
                             bintybint, &
                             atom_sigma_r_d, atom_sigma_t_d, rcut_hard_d, &
                             central_weight_d, species_d, i_beg_d, i_end_d, &
                             radial_enhancement, species_multiplicity_d, &
                             W_d, S_d, size_species_1, gpu_stream) &
         bind(C, name="gpu_get_cnk")
         use iso_c_binding
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: k2_start_d
         type(c_ptr) :: gpu_stream
         type(c_ptr), value :: radial_exp_coeff_d
         type(c_ptr), value :: angular_exp_coeff_d
         type(c_ptr), value :: cnk_d
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_soap
         integer(c_int), value :: k_max
         integer(c_int), value :: n_max
         integer(c_int), value :: l_max
         integer(c_int), value :: bintybint
         integer(c_int), value :: radial_enhancement
         integer(c_int), value :: size_species_1
         type(c_ptr), value :: rcut_hard_d
         type(c_ptr), value :: species_multiplicity_d
         type(c_ptr), value :: atom_sigma_r_d
         type(c_ptr), value :: atom_sigma_t_d
         type(c_ptr), value :: central_weight_d
         type(c_ptr), value :: species_d
         type(c_ptr), value :: i_beg_d
         type(c_ptr), value :: i_end_d
         type(c_ptr), value :: W_d
         type(c_ptr), value :: S_d
      end subroutine

      subroutine gpu_get_plm_array_global(plm_array_global_d, n_atom_pairs, kmax, &
                                          lmax, thetas_d, gpu_stream) &
         bind(C, name="gpu_get_plm_array_global")
         use iso_c_binding
         type(c_ptr), value :: plm_array_global_d
         type(c_ptr), value :: thetas_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: kmax
         integer(c_int), value :: lmax
      end subroutine

      subroutine gpu_get_exp_coeff_array(eimphi_global_d, thetas_d, &
                                         rjs_d, phis_d, &
                                         mask_d, atom_simga_in_d, atom_simga_scaling_d, &
                                         rcut, n_atom_piars, n_spceies, lmax, kmax, &
                                         prefl_array_global_d, plm_array_global_d, &
                                         plm_array_der_global_d, &
                                         prefl_array_global_der_d, &
                                         preflm_d, exp_coeff_d, &
                                         c_do_derivatives, &
                                         eimphi_rad_der_global_d, eimphi_azi_der_global_d, &
                                         plm_array_div_sin_d, plm_array_der_mul_sin_d, &
                                         exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d, &
                                         gpu_stream) &
         bind(C, name="gpu_get_exp_coeff_array")
         use iso_c_binding
         type(c_ptr), value :: eimphi_global_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: phis_d
         type(c_ptr), value :: mask_d
         type(c_ptr), value :: exp_coeff_d
         type(c_ptr) :: gpu_stream
         type(c_ptr), value :: atom_simga_in_d
         type(c_ptr), value :: atom_simga_scaling_d
         type(c_ptr), value :: prefl_array_global_der_d
         integer(c_int), value :: n_atom_piars
         integer(c_int), value :: n_spceies
         integer(c_int), value :: lmax
         integer(c_int), value :: kmax
         real(c_double), value :: rcut
         type(c_ptr), value :: eimphi_rad_der_global_d
         type(c_ptr), value :: eimphi_azi_der_global_d
         type(c_ptr), value :: thetas_d
         type(c_ptr), value :: preflm_d
         type(c_ptr), value :: prefl_array_global_d
         type(c_ptr), value :: plm_array_global_d
         type(c_ptr), value :: plm_array_der_global_d
         logical(c_bool), value :: c_do_derivatives
         type(c_ptr), value :: exp_coeff_rad_der_d
         type(c_ptr), value :: exp_coeff_azi_der_d
         type(c_ptr), value :: exp_coeff_pol_der_d
         type(c_ptr), value :: plm_array_div_sin_d
         type(c_ptr), value :: plm_array_der_mul_sin_d
      end subroutine

      subroutine gpu_get_radial_exp_coeff_poly3gauss(radial_exp_coeff_d, radial_exp_coeff_der_d, &
                                                     i_beg_d, i_end_d, &
                                                     global_scaling_d, &
                                                     size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species, &
                                                     c_do_derivatives, bintybint, &
                                                     rcut_hard_d, &
                                                     k2_i_site_d, k_2start_d, &
                                                     gpu_stream) &
         bind(C, name="gpu_get_radial_exp_coeff_poly3gauss")
         use iso_c_binding
         type(c_ptr), value :: radial_exp_coeff_d
         type(c_ptr), value :: radial_exp_coeff_der_d
         type(c_ptr), value :: i_beg_d
         type(c_ptr), value :: i_end_d
         type(c_ptr), value :: global_scaling_d
         type(c_ptr), value :: rcut_hard_d
         type(c_ptr), value :: k2_i_site_d
         type(c_ptr), value :: k_2start_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: size_radial_exp_coeff_one
         integer(c_int), value :: size_radial_exp_coeff_two
         integer(c_int), value :: n_species
         integer(c_int), value :: bintybint
         logical(c_bool), value :: c_do_derivatives
      end subroutine

      subroutine gpu_radial_poly3gauss(n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, &
                                       ntemp, do_derivatives, exp_coeff_d, exp_coeff_der_d, rcut_soft_d, atom_sigma_d, &
                                       exp_coeff_temp1_d, exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end, &
                                       atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, &
                                       nf_d, ntemp_der_d, W_d, gpu_stream) &
         bind(C, name="gpu_radial_poly3gauss")
         use iso_c_binding
         type(c_ptr), value :: mask_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: rcut_hard_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: exp_coeff_d
         type(c_ptr), value :: exp_coeff_der_d
         type(c_ptr), value :: rcut_soft_d
         type(c_ptr), value :: atom_sigma_d
         type(c_ptr), value :: exp_coeff_temp1_d
         type(c_ptr), value :: exp_coeff_temp2_d
         type(c_ptr), value :: exp_coeff_der_temp_d
         type(c_ptr), value :: i_beg
         type(c_ptr), value :: i_end
         type(c_ptr), value :: atom_sigma_scaling_d
         type(c_ptr), value :: amplitude_scaling_d
         type(c_ptr), value :: alpha_max_d
         type(c_ptr), value :: nf_d
         type(c_ptr), value :: W_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_species
         integer(c_int), value :: n_sites
         integer(c_int), value :: ntemp
         integer(c_int), value :: n_max
         integer(c_int), value :: mode
         integer(c_int), value :: radial_enhancement
         integer(c_int), value :: ntemp_der_d
         logical(c_bool), value :: do_derivatives
      end subroutine

      subroutine gpu_radial_poly3(n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, &
                                  ntemp, do_derivatives, exp_coeff_d, exp_coeff_der_d, rcut_soft_d, atom_sigma_d, &
                                  exp_coeff_temp1_d, exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end, &
                                  atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, &
                                  nf_d, ntemp_der_d, W_d, do_central_d, central_weight_d, gpu_stream) &
         bind(C, name="gpu_radial_poly3")
         use iso_c_binding
         type(c_ptr), value :: mask_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: rcut_hard_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: exp_coeff_d
         type(c_ptr), value :: exp_coeff_der_d
         type(c_ptr), value :: rcut_soft_d
         type(c_ptr), value :: atom_sigma_d
         type(c_ptr), value :: exp_coeff_temp1_d
         type(c_ptr), value :: exp_coeff_temp2_d
         type(c_ptr), value :: exp_coeff_der_temp_d
         type(c_ptr), value :: i_beg
         type(c_ptr), value :: i_end
         type(c_ptr), value :: atom_sigma_scaling_d
         type(c_ptr), value :: amplitude_scaling_d
         type(c_ptr), value :: alpha_max_d
         type(c_ptr), value :: nf_d
         type(c_ptr), value :: W_d
         type(c_ptr), value :: do_central_d
         type(c_ptr), value :: central_weight_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_species
         integer(c_int), value :: n_sites
         integer(c_int), value :: ntemp
         integer(c_int), value :: n_max
         integer(c_int), value :: mode
         integer(c_int), value :: radial_enhancement
         integer(c_int), value :: ntemp_der_d
         logical(c_bool), value :: do_derivatives
      end subroutine

      subroutine gpu_get_2b_forces_energies(i_beg, i_end, n_sparse, energies_d, e0, n_neigh_d, do_forces, forces_d, virial_d, &
                                            rjs_d, rcut, species_d, neighbor_species_d, sp1, sp2, buffer, delta, cutoff_d, &
                                            Qs_d, sigma, alphas_d, xyz_d, gpu_stream) &
         bind(C, name="gpu_get_2b_forces_energies")
         use iso_c_binding
         type(c_ptr), value :: energies_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: forces_d
         type(c_ptr), value :: virial_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: species_d
         type(c_ptr), value :: neighbor_species_d
         type(c_ptr), value :: cutoff_d
         type(c_ptr), value :: Qs_d
         type(c_ptr), value :: alphas_d
         type(c_ptr), value :: xyz_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_sparse
         integer(c_int), value :: sp1
         integer(c_int), value :: sp2
         logical(c_bool), value :: do_forces
         real(c_double), value :: e0
         real(c_double), value :: rcut
         real(c_double), value :: buffer
         real(c_double), value :: delta
         real(c_double), value :: sigma
      end subroutine

      subroutine gpu_get_core_pot_energy_and_forces(i_beg, i_end, do_forces, species_d, sp1, sp2, n_neigh_d, neighbor_species_d, &
                                                    rjs_d, n_sparse, x_d, V_d, dVdx2_d, yp1, ypn, xyz_d, forces_d, virial_d, &
                                                    energies_d, gpu_stream) &
         bind(C, name="gpu_get_core_pot_energy_and_forces")
         use iso_c_binding
         type(c_ptr), value :: species_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: neighbor_species_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: x_d
         type(c_ptr), value :: V_d
         type(c_ptr), value :: dVdx2_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: forces_d
         type(c_ptr), value :: virial_d
         type(c_ptr), value :: energies_d
         type(c_ptr) :: gpu_stream
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: sp1
         integer(c_int), value :: sp2
         integer(c_int), value :: n_sparse
         logical(c_bool), value :: do_forces
         real(c_double), value :: yp1
         real(c_double), value :: ypn
      end subroutine

      subroutine get_orthonormalization_matrix_poly3_gpu(alpha_max, S, W, handle, stream) bind(C,name="get_orthonormalization_matrix_poly3")
         use iso_c_binding
         type(c_ptr), value :: S
         type(c_ptr), value :: W
         type(c_ptr) :: handle
         type(c_ptr) :: stream
         integer(c_int), value :: alpha_max
      end subroutine

      subroutine get_orthonormalization_matrix_poly3gauss_gpu(alpha_max, atom_sigma_in, rcut_hard_in, S, W,handle,stream) bind(C,name="get_orthonormalization_matrix_poly3gauss")
         use iso_c_binding
         type(c_ptr), value :: S
         type(c_ptr), value :: W
         type(c_ptr) :: handle
         type(c_ptr) :: stream
         integer(c_int), value :: alpha_max
         real(c_double), value :: atom_sigma_in
         real(c_double), value :: rcut_hard_in
      end subroutine

      subroutine orthonormalization_copy_to_global_matrix(src, dest,src_rowsize, dest_start, dest_rowsize, stream) bind(C,name="copy_to_global_matrix")
         use iso_c_binding
         type(c_ptr), value :: src
         type(c_ptr), value :: dest
         type(c_ptr) :: stream
         integer(c_int), value :: src_rowsize
         integer(c_int), value :: dest_start
         integer(c_int), value :: dest_rowsize
      end subroutine

      subroutine gpu_3b(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, e0, cutoff, stream, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, do_forces, rcut,buffer, sigma,qs,str,i_beg, i_end, energy_d, forces_d, virials_d, kappas_array_d) bind(C,name="gpu_3b")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: alpha
         type(c_ptr), value :: cutoff
         type(c_ptr) :: stream
         integer(c_int), value :: n_sparse
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_atom_pairs
         integer(c_int), value :: n_sites0
         integer(c_int), value :: sp0
         integer(c_int), value :: sp1
         integer(c_int), value :: sp2
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         real(c_double), value :: delta
         real(c_double), value :: e0
         real(c_double), value :: rcut
         real(c_double), value :: buffer
         type(c_ptr), value :: rjs
         type(c_ptr), value :: xyz
         type(c_ptr), value :: n_neigh
         type(c_ptr), value :: species
         type(c_ptr), value :: neighbors_list
         type(c_ptr), value :: neighbor_species
         type(c_ptr), value :: sigma
         type(c_ptr), value :: qs
         type(c_ptr), value :: energy_d
         type(c_ptr), value :: forces_d
         type(c_ptr), value :: virials_d
         type(c_ptr), value :: kappas_array_d
         logical(c_bool), value :: do_forces
         character(kind=c_char), dimension(*) :: str
      end subroutine

      subroutine gpu_create_kappas(kappas_array_d, n_neigh_host, stream, n_sites) bind(C, name="create_kappas_cwrap")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: kappas_array_d
         type(c_ptr), value :: n_neigh_host
         type(c_ptr) :: stream
         integer(c_int), value :: n_sites
      end subroutine

      subroutine gpu_setup_3bresult_arrays(energy_3b_d, forces_3b_d, virials_3b_d,stream,do_forces, n_sites, n_sites0) bind(C,name="setup_3bresult_arrays_cwrap")
         use iso_c_binding
         implicit none
         type(c_ptr) :: energy_3b_d
         type(c_ptr) :: forces_3b_d
         type(c_ptr) :: virials_3b_d
         type(c_ptr) :: stream
         logical(c_bool), value :: do_forces
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_sites0
      end subroutine

      subroutine gpu_cleanup_3bresult_arrays(energy_3b_d, forces_3b_d, virials_3b_d,energy_3b_h,forces_3b_h,virials_3b_h,stream,do_forces, n_sites, n_sites0) bind(C,name="cleanup_3bresult_arrays_cwrap")
         use iso_c_binding
         implicit none
         type(c_ptr) :: energy_3b_d
         type(c_ptr) :: forces_3b_d
         type(c_ptr) :: virials_3b_d
         type(c_ptr) :: energy_3b_h
         type(c_ptr) :: forces_3b_h
         type(c_ptr) :: virials_3b_h
         type(c_ptr) :: stream
         logical(c_bool), value :: do_forces
         integer(c_int), value :: n_sites
         integer(c_int), value :: n_sites0
      end subroutine

      subroutine gpu_device_sync() bind(C, name="gpu_device_sync")
         use iso_c_binding
         implicit none
      end subroutine gpu_device_sync

      subroutine gpu_stream_sync(stream) bind(C, name="gpu_stream_sync")
         use iso_c_binding
         implicit none
         type(c_ptr) :: stream
      end subroutine gpu_stream_sync

      subroutine gpu_meminfo() bind(C, name="gpu_meminfo")
         use iso_c_binding
         implicit none
      end subroutine

!     The device memory ledger in cuda_wrappers.cu. gpu_meminfo above prints the
!     device's own free/total and nothing else; these also report what THIS rank
!     has asked for, which is the quantity the batching keywords control.
      subroutine gpu_mem_report(label) bind(C, name="gpu_mem_report")
         use iso_c_binding
         implicit none
         character(kind=c_char), dimension(*) :: label
      end subroutine

!     Bytes. Every argument is written; pass all four.
      subroutine gpu_mem_stats(dev_free, dev_total, live, peak) bind(C, name="gpu_mem_stats")
         use iso_c_binding
         implicit none
         integer(c_size_t) :: dev_free
         integer(c_size_t) :: dev_total
         integer(c_size_t) :: live
         integer(c_size_t) :: peak
      end subroutine

!     Tell the C side what the input asked for, so an out-of-memory abort can
!     name the keyword that fixes it and its current value.
      subroutine gpu_mem_set_hint(max_gb, n_batches) bind(C, name="gpu_mem_set_hint")
         use iso_c_binding
         implicit none
         real(c_double), value :: max_gb
         integer(c_int), value :: n_batches
      end subroutine

!     Visible devices. cuda_set_device assigns them round-robin, so this is what
!     says how many ranks share one card's memory.
      integer(c_int) function gpu_device_count() bind(C, name="gpu_device_count")
         use iso_c_binding
         implicit none
      end function

      subroutine gpu_check_error() bind(C, name="gpu_check_error")
         use iso_c_binding
         implicit none
      end subroutine

      subroutine gpu_inclusive_scan_int(size, array_d, stream) bind(C, name="gpu_inclusive_scan_int")
         use iso_c_binding
         integer(c_int), value :: size
         type(c_ptr), value :: array_d
         type(c_ptr) :: stream
      end subroutine

      subroutine gpu_get_electrostatics_nk(i_beg, i_end, n_pairs, &
                                           n_neigh_d, n_neigh_index_d, rjs_d, xyz_d, r_cut, &
                                           nk_out_d, nk_flags_d, nk_flags_sum_d, stream) bind(C, name="gpu_get_electrostatics_nk")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: nk_flags_d
         type(c_ptr), value :: nk_flags_sum_d
         type(c_ptr), value :: n_neigh_index_d
         real(c_double), value :: r_cut
         type(c_ptr), value :: nk_out_d
         type(c_ptr) :: stream
      end subroutine gpu_get_electrostatics_nk

      subroutine gpu_set_electrostatics_k_index(i_beg, i_end, n_pairs, n_sites0, neighbors_list_d, &
                                                rjs_d, xyz_d, &
                                                charges_d, &
                                                neighbor_charges_index_d, &
                                                charge_gradients_d, charge_gradients_index_d, &
                                                k_index_d, j2_index_d, rjs_index_d, xyz_k_d, &
                                                nk_sum_flags_d, stream) bind(C, name="gpu_set_electrostatics_k_index")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_sites0
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: neighbors_list_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: charge_gradients_d
         type(c_ptr), value :: charge_gradients_index_d
         type(c_ptr), value :: k_index_d
         type(c_ptr), value :: j2_index_d
         type(c_ptr), value :: rjs_index_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr), value :: charges_d
         type(c_ptr), value :: neighbor_charges_index_d
         type(c_ptr) :: stream
      end subroutine gpu_set_electrostatics_k_index

      subroutine gpu_get_electrostatics_energies( &
         i_beg, &
         nk_max, &
         energies_d, &
         forces_d, &
         virial_d, &
         j2_index_d, &
         n_sites, &
         this_n_sites, &
         this_n_pairs, &
         n_neigh_index_d, &
         charges_d, &
         charge_gradients_d, &
         neighbor_charges_index_d, &
         rjs_index_d, &
         xyz_k_d, &
         dsf_alpha, &
         rcut, &
         rcut_in, &
         rcut_width, &
         pair_energy_rcut, &
         pair_energy_rcut_der, &
         do_cosine_damping, &
         do_forces, &
         stream) bind(C, name="gpu_get_electrostatics_energies")
         use iso_c_binding
         implicit none
         integer(c_int), value :: this_n_sites
         integer(c_int), value :: n_sites
         integer(c_int), value :: this_n_pairs
         integer(c_int), value :: i_beg
         integer(c_int), value :: nk_max
         type(c_ptr), value :: charges_d
         type(c_ptr), value :: neighbor_charges_index_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr), value :: charge_gradients_d
         type(c_ptr), value :: rjs_index_d
         type(c_ptr), value :: energies_d
         type(c_ptr), value :: n_neigh_index_d
         type(c_ptr), value :: forces_d
         type(c_ptr), value :: j2_index_d
         type(c_ptr), value :: virial_d !
         logical(c_bool), value :: do_forces
         logical(c_bool), value :: do_cosine_damping
         type(c_ptr) :: stream
         real(c_double), value :: dsf_alpha
         real(c_double), value :: pair_energy_rcut
         real(c_double), value :: pair_energy_rcut_der
         real(c_double), value :: rcut
         real(c_double), value :: rcut_in
         real(c_double), value :: rcut_width
      end subroutine gpu_get_electrostatics_energies

      subroutine gpu_get_pair_distribution_nk(i_beg, i_end, n_pairs, n_sites0, neighbors_list_d, &
                                              n_neigh_d, neighbor_species_d, &
                                              species_d, rjs_d, xyz_d, r_min, r_max, r_cut, &
                                              buffer, &
                    nk_out_d, nk_flags_d, nk_flags_sum_d, species_1, species_2, stream) bind(C, name="gpu_get_pair_distribution_nk")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_sites0
         integer(c_int), value :: species_1
         integer(c_int), value :: species_2
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: neighbors_list_d
         type(c_ptr), value :: n_neigh_d
         type(c_ptr), value :: neighbor_species_d
         type(c_ptr), value :: species_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: nk_flags_d
         type(c_ptr), value :: nk_flags_sum_d
         real(c_double), value :: r_min
         real(c_double), value :: r_max
         real(c_double), value :: buffer
         real(c_double), value :: r_cut
         type(c_ptr), value :: nk_out_d
         type(c_ptr) :: stream
      end subroutine gpu_get_pair_distribution_nk

      subroutine gpu_set_pair_distribution_k_index(i_beg, i_end, n_pairs, n_sites0, neighbors_list_d, &
           rjs_d, xyz_d, k_index_d, j2_index_d, rjs_index_d,  xyz_k_d, nk_flags_d,  nk_sum_flags_d, stream)  bind(C,name="gpu_set_pair_distribution_k_index")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_sites0
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: neighbors_list_d
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: nk_flags_d
         type(c_ptr), value :: k_index_d
         type(c_ptr), value :: j2_index_d
         type(c_ptr), value :: rjs_index_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr) :: stream
      end subroutine gpu_set_pair_distribution_k_index

      subroutine gpu_setup_matrix_forces(i_beg, i_end, n_pairs, n_sites0, neighbors_list_d, &
              xyz_all_d, k_index_d, j2_index_d, xyz_k_d, nk_flags_d, nk_flags_sum_d, stream) bind(C, name="gpu_setup_matrix_forces")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i_beg
         integer(c_int), value :: i_end
         integer(c_int), value :: n_sites0
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: neighbors_list_d
         type(c_ptr), value :: xyz_all_d
         type(c_ptr), value :: nk_flags_d
         type(c_ptr), value :: nk_flags_sum_d
         type(c_ptr), value :: k_index_d
         type(c_ptr), value :: j2_index_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr) :: stream
      end subroutine gpu_setup_matrix_forces

      subroutine gpu_exp_force_virial_collection(n_k, forces0, energy_scale, fi, &
                                                 j2_list, virial, xyz, stream) bind(C, name="gpu_exp_force_virial_collection")

         use iso_c_binding
         implicit none
         integer(c_int), value :: n_k
         real(c_double), value :: energy_scale
         type(c_ptr), value :: forces0
         type(c_ptr), value :: xyz
         type(c_ptr), value :: virial
         type(c_ptr), value :: j2_list
         type(c_ptr), value :: fi
         type(c_ptr) :: stream
      end subroutine gpu_exp_force_virial_collection

      subroutine gpu_get_fi_dgemv(i, n_samples_sf, n_k, dermat_d,&
           & prefactor_d, fi_d, cublas_handle, stream) bind(C, name&
           &="gpu_get_fi_dgemv")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i
         integer(c_int), value :: n_samples_sf
         integer(c_int), value :: n_k
         type(c_ptr), value :: dermat_d
         type(c_ptr), value :: prefactor_d
         type(c_ptr), value :: fi_d
         type(c_ptr), value :: cublas_handle
         type(c_ptr), value :: stream
      end subroutine gpu_get_fi_dgemv

      subroutine gpu_hadamard_vec_mat_product(n_samples_sf, n_k,&
           & all_scattering_factors_d, dermat_d, stream) bind(C&
           &, name="gpu_hadamard_vec_mat_product")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_samples_sf
         integer(c_int), value :: n_k
         type(c_ptr), value :: dermat_d
         type(c_ptr), value :: all_scattering_factors_d
         type(c_ptr) :: stream
      end subroutine gpu_hadamard_vec_mat_product

      subroutine gpu_get_Gka(i, n_k, n_samples, Gka_d, Gk_d, xyz_k_d, stream) bind(C&
           &, name="gpu_get_Gka")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i
         integer(c_int), value :: n_samples
         integer(c_int), value :: n_k
         type(c_ptr), value :: Gka_d
         type(c_ptr), value :: Gk_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr) :: stream
      end subroutine gpu_get_Gka

      subroutine gpu_get_Gka_inplace(i, n_k, n_samples, Gk_d, xyz_k_d, stream) bind(C&
           &, name="gpu_get_Gka_inplace")
         use iso_c_binding
         implicit none
         integer(c_int), value :: i
         integer(c_int), value :: n_samples
         integer(c_int), value :: n_k
         type(c_ptr), value :: Gk_d
         type(c_ptr), value :: xyz_k_d
         type(c_ptr) :: stream
      end subroutine gpu_get_Gka_inplace

      subroutine gpu_print_pointer_int(p) bind(C&
           &, name="gpu_print_pointer_int")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: p
      end subroutine gpu_print_pointer_int

      subroutine gpu_print_pointer_double(p) bind(C&
           &, name="gpu_print_pointer_double")
         use iso_c_binding
         implicit none
         type(c_ptr), value :: p
      end subroutine gpu_print_pointer_double

      subroutine gpu_get_pair_distribution_and_ders(pair_distribution_d, pair_distribution_der_d, &
           n_k, n_samples, kde_sigma, x_d, dV_d, rjs_d, pdf_factor, der_factor, stream) bind(C&
           &, name="gpu_get_pair_distribution_and_ders")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_k
         integer(c_int), value :: n_samples
         real(c_double), value :: kde_sigma
         real(c_double), value :: pdf_factor
         real(c_double), value :: der_factor
         type(c_ptr), value :: pair_distribution_d
         type(c_ptr), value :: pair_distribution_der_d
         type(c_ptr), value :: x_d
         type(c_ptr), value :: dV_d
         type(c_ptr), value :: rjs_d
         type(c_ptr) :: stream
      end subroutine gpu_get_pair_distribution_and_ders

      subroutine gpu_get_pair_distribution_only(pair_distribution_d, &
           n_k, n_samples, kde_sigma, x_d, dV_d, rjs_d, pdf_factor, der_factor, stream) bind(C&
           &, name="gpu_get_pair_distribution_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_k
         integer(c_int), value :: n_samples
         real(c_double), value :: kde_sigma
         real(c_double), value :: pdf_factor
         real(c_double), value :: der_factor
         type(c_ptr), value :: pair_distribution_d
         type(c_ptr), value :: x_d
         type(c_ptr), value :: dV_d
         type(c_ptr), value :: rjs_d
         type(c_ptr) :: stream
      end subroutine gpu_get_pair_distribution_only

      subroutine gpu_get_pair_distribution_only_falloc( &
          pair_distribution_d, &
          pdf_to_reduce, &
           n_k, &
           n_samples, &
           kde_sigma, &
           x_d, &
           dV_d, &
           rjs_d, &
           pdf_factor, &
           der_factor, &
           stream) bind(C&
           &, name="gpu_get_pair_distribution_only_falloc")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_k
         integer(c_int), value :: n_samples
         real(c_double), value :: kde_sigma
         real(c_double), value :: pdf_factor
         real(c_double), value :: der_factor
         type(c_ptr), value :: pair_distribution_d
         type(c_ptr), value :: x_d
         type(c_ptr), value :: dV_d
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: pdf_to_reduce
         type(c_ptr) :: stream
      end subroutine gpu_get_pair_distribution_only_falloc

      subroutine gpu_get_pair_distribution_der_only(pair_distribution_der_d, &
           n_k, n_samples, kde_sigma, x_d, dV_d, rjs_d, pdf_factor, der_factor, stream) bind(C&
           &, name="gpu_get_pair_distribution_der_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_k
         integer(c_int), value :: n_samples
         real(c_double), value :: kde_sigma
         real(c_double), value :: pdf_factor
         real(c_double), value :: der_factor
         type(c_ptr), value :: pair_distribution_der_d
         type(c_ptr), value :: x_d
         type(c_ptr), value :: dV_d
         type(c_ptr), value :: rjs_d
         type(c_ptr) :: stream
      end subroutine gpu_get_pair_distribution_der_only

 subroutine gpu_set_Gk(nk, n_samples, k_index_d, Gk_d, pair_distribution_partial_der_d, c_factor, stream) bind(C, name="gpu_set_Gk")
         use iso_c_binding
         implicit none
         integer(c_int), value :: nk
         integer(c_int), value :: n_samples
         real(c_double), value :: c_factor
         type(c_ptr), value :: k_index_d
         type(c_ptr), value :: Gk_d
         type(c_ptr), value :: pair_distribution_partial_der_d
         type(c_ptr) :: stream
      end subroutine gpu_set_Gk

      subroutine  gpu_set_pair_distribution_k_index_only(n_pairs, k_index_d, nk_sum_flags_d, stream )bind(C,name="gpu_set_pair_distribution_k_index_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: k_index_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr) :: stream
      end subroutine gpu_set_pair_distribution_k_index_only

      subroutine  gpu_set_pair_distribution_j2_only(n_pairs, n_sites, neighbors_list_d, j2_d, nk_sum_flags_d, stream )bind(C,name="gpu_set_pair_distribution_j2_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_pairs
         integer(c_int), value :: n_sites
         type(c_ptr), value :: j2_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr), value :: neighbors_list_d
         type(c_ptr) :: stream
      end subroutine gpu_set_pair_distribution_j2_only

      subroutine  gpu_set_pair_distribution_rjs_only(n_pairs, rjs, rjs_d, nk_sum_flags_d, stream )bind(C,name="gpu_set_pair_distribution_rjs_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: rjs_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr), value :: rjs
         type(c_ptr) :: stream
      end subroutine gpu_set_pair_distribution_rjs_only

      subroutine  gpu_set_pair_distribution_xyz_only(n_pairs, xyz, xyz_d, nk_sum_flags_d, stream )bind(C,name="gpu_set_pair_distribution_xyz_only")
         use iso_c_binding
         implicit none
         integer(c_int), value :: n_pairs
         type(c_ptr), value :: xyz_d
         type(c_ptr), value :: nk_sum_flags_d
         type(c_ptr), value :: xyz
         type(c_ptr) :: stream
      end subroutine gpu_set_pair_distribution_xyz_only

      function host_alloc(s) bind(c, name='host_alloc')
         use iso_c_binding
         implicit none
         type(c_ptr) :: host_alloc
         integer(c_size_t), value :: s
      end function host_alloc

      subroutine host_free(p) bind(c, name='host_free')
         use iso_c_binding
         implicit none
         type(c_ptr), intent(in), value :: p
      end subroutine host_free

      subroutine cpy_htoh_pinned(s, d, size) bind(c, name='cpy_htoh_pinned')
         use iso_c_binding
         implicit none
         integer(c_size_t), value :: size
         type(c_ptr), value :: s
         type(c_ptr), value :: d
      end subroutine cpy_htoh_pinned

      !
!      subroutine gpu_2b(n_sparse, n_sites, sp1, sp2, alpha, delta, cutoff, stream, rjs, xyz, n_neigh, species, neighbor_species, do_forces, rcut,buffer, sigma,qs,n_neigh_host) bind(C,name="gpu_2b")
!        use iso_c_binding
!        implicit none
!        type(c_ptr),value :: alpha, cutoff
!        type(c_ptr) :: stream
!        integer(c_int),value :: n_sparse, n_sites, sp1, sp2
!        real(c_double),value :: delta, rcut,buffer,sigma
!        type(c_ptr),value :: rjs,xyz,n_neigh,species,neighbor_species,qs,n_neigh_host
!        logical(c_bool), value :: do_forces
!      end subroutine

   END INTERFACE
END MODULE F_B_C

! subroutine gpu_malloc_double(a_d,n) bind(C,name="cuda_malloc_double")
!   use iso_c_binding
!   implicit none
!   type(c_ptr) :: a_d
!   integer(c_int),value :: n
! end subroutine

! subroutine one_for_all(soap, kernels, kernels_copy, Qs, energies, &
!       delta, zeta, e0, &
!       n_sites, n_soap, n_sparse, &
!       size_kernels, size_soap, size_Qs, size_alphas, size_energies) &
!       bind(C, name="wrappers_all")
!   use iso_c_binding
!   implicit none
!   type(c_ptr), value :: soap, kernels,kernels_copy, Qs, energies
!   real(c_double), value :: delta, zeta, e0
!   integer(c_int), value :: n_sites, n_soap, n_sparse
!   integer(c_int), value :: size_kernels, size_soap, size_Qs, size_alphas, size_energies
! end subroutine

! subroutine gpu_soap_energies_forces_virial(n_neigh_d,n_sites,maxnn, &
!                                         Qss_d,n_soap,neighbors_beg_d, &
!                                         soap_der_d, &
!                                         xyz_d, virial_d, &
!                                         neighbors_list_d, n_sites0, forces_d, &
!                                         cuhandle, kernels_der_d, Qs_copy_d, &
!                                         n_sparse, cdelta_force, &
!                                         alphas_d, &
!                                         kernels_d, mzetam, size_kernels, &
!                                         do_forces, &
!                                         energies_d, cdelta_ene, e0,  size_energies, &
!                                         Qs_d, size_Qs, &
!                                         kernels_copy_d, &
!                                         zeta, &
!                                         soap_d ) &
!             bind(C,name="gpu_soap_energies_forces_virial")
!   use iso_c_binding
!   type(c_ptr), value :: Qss_d, n_neigh_d,neighbors_beg_d, neighbors_list_d
!   type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
!   integer(c_int), value :: n_sites,n_soap,maxnn, n_sites0, size_energies
!   type(c_ptr), value :: cuhandle,kernels_der_d, Qs_copy_d, alphas_d, kernels_d
!   type(c_ptr), value :: energies_d, Qs_d, kernels_copy_d, soap_d
!   integer(c_int),value :: n_sparse, size_kernels, size_Qs
!   logical, value:: do_forces
!   real(c_double), value :: cdelta_ene, cdelta_force, mzetam, e0, zeta
! end subroutine

! subroutine cpy_double_htod(a,a_d,n) bind(C,name="cuda_cpy_double_htod")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_double_complex_htod(a,a_d,n) bind(C,name="cuda_cpy_double_complex_htod")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_int_htod(a,a_d,n) bind(C,name="cuda_cpy_int_htod")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_bool_htod(a,a_d,n) bind(C,name="cuda_cpy_bool_htod")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_double_dtod(b_d,a_d,n) bind(C,name="cuda_cpy_double_dtod")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,b_d
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_bool_dtoh(a_d,a,n) bind(C,name="cuda_cpy_bool_dtoh")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_double_dtoh(a_d,a,n) bind(C,name="cuda_cpy_double_dtoh")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine cpy_double_complex_dtoh(a_d,a,n) bind(C,name="cuda_cpy_double_complex_dtoh")
!   use iso_c_binding
!   implicit none
!   type(c_ptr),value :: a_d,a
!   integer(c_int),value :: n
! end subroutine

! subroutine gpu_malloc_double_complex(a_d,n) bind(C,name="cuda_malloc_double_complex")
!   use iso_c_binding
!   implicit none
!     type(c_ptr) :: a_d
!     integer(c_int),value :: n
! end subroutine

! subroutine gpu_malloc_int(a_d,n) bind(C,name="cuda_malloc_int")
!   use iso_c_binding
!   implicit none
!     type(c_ptr) :: a_d
!     integer(c_int),value :: n
! end subroutine

! subroutine gpu_malloc_bool(a_d,n) bind(C,name="cuda_malloc_bool")
!   use iso_c_binding
!   implicit none
!     type(c_ptr) :: a_d
!     integer(c_int),value :: n
! end subroutine

! gpu_var_mod / gpu_var_int_mod / gpu_var_double_mod were removed here.
!
! They wrapped a device pointer in a class with type-bound galloc/gfree, but no
! instance was ever created: every declaration of gpu_var_int_class or
! gpu_var_double_class in the tree was commented out, and the three modules were
! `use`d by turbogap.f90 and exp_interface.f90 without a single reference. Device
! buffers are allocated directly with gpu_malloc_async / gpu_free_async.
!
! The class also carried a real defect that could not be fixed without changing
! its interface: malloc() took a stream and dealloc() did not, so the latter had
! to call the blocking gpu_free on a pointer from the async pool. Deleting the
! class removes the last such mismatch in the tree -- see
! tools/gpu_check_alloc_pairs.py.

! module gpu_variables was removed here too. It held gpu_neigh and gpu_exp,
! containers of the gpu_var classes above, plus allocate/deallocate/check
! helpers for them. Nothing `use`d the module and nothing declared either type;
! the gpu_neigh_storage_type that electrostatics, exp_interface and gpu_context
! do use is a different type with a similar name.
