
! module load gcc/10.3.0 openblas openmpi cuda
! rm cuda_wrappers.o ; nvcc -arch=sm_70 -c src/cuda_wrappers.cu ; make clean; make -B
!  srun  --time=00:10:00 --partition=gputest --account=project_2000634 --nodes=1 --ntasks-per-node=4  --cpus-per-task=1 --gres=gpu:a100:4  ../turbogap_dev/bin/turbogap md
! gnupot
!plot "< paste trajectory_out.xyz_64k_oiriginal trajectory_out.xyz | awk 'NR>2{print $5,$13}'", x
MODULE F_B_C
    INTERFACE
      subroutine gpu_malloc_all(a_d,n,gpu_stream) bind(C,name="cuda_malloc_all")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d,gpu_stream
        integer(c_size_t),value :: n
      end subroutine

      subroutine gpu_malloc_all_blocking(a_d,n) bind(C,name="cuda_malloc_all_blocking")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
        integer(c_size_t),value :: n
      end subroutine
      
      subroutine gpu_device_reset() bind(C,name="cuda_device_reset")
        use iso_c_binding
        implicit none
      end subroutine
      subroutine gpu_free(a_d) bind(C,name="cuda_free")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
      end subroutine

      subroutine gpu_free_async(a_d,gpu_stream) bind(C,name="cuda_free_async")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
        type(c_ptr) :: gpu_stream
      end subroutine

      subroutine cpy_htod(a,a_d,n, gpu_stream) bind(C,name="cuda_cpy_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        type(c_ptr) :: gpu_stream
        integer(c_size_t),value :: n
      end subroutine

      subroutine cpy_htod_blocking(a,a_d,n) bind(C,name="cuda_cpy_htod_blocking")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_size_t),value :: n
      end subroutine

      subroutine cpy_dtod(b_d,a_d,n, gpu_stream) bind(C,name="cuda_cpy_dtod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,b_d
        type(c_ptr) :: gpu_stream
        integer(c_size_t),value :: n
      end subroutine

      subroutine cpy_dtoh(a_d,a,n, gpu_stream) bind(C,name="cuda_cpy_dtoh")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        type(c_ptr) :: gpu_stream
        integer(c_size_t),value :: n
      end subroutine
      

      subroutine cpy_dtoh_blocking(a_d,a,n) bind(C,name="cuda_cpy_dtoh_blocking")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_size_t),value :: n
      end subroutine
      subroutine gpu_vector_fill_curand(a_d,n,c) bind(C,name="GPU_fill_rand")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d
        integer(c_int),value :: n,c
      end subroutine

      subroutine create_cublas_handle(cubhandle, gpu_stream)bind(C,name="create_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle, gpu_stream
      end subroutine

      subroutine destroy_cublas_handle(cubhandle, gpu_stream)bind(C,name="destroy_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle,gpu_stream
      end subroutine

      subroutine gpu_blas_mmul_t_n(cubhandle, Qs_d, soap_d, kernels_d, &
                         n_sparse, n_soap, n_sites)  bind(C,name="gpu_blas_mmul_t_n")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,kernels_d,Qs_d,soap_d
        integer(c_int),value :: n_sparse, n_soap, n_sites
      end subroutine

      subroutine gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                                  n_soap, n_sites, cdelta) bind(C,name="gpu_blas_mmul_n_t")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,kernels_der_d, Qs_copy_d, Qss_d
        integer(c_int),value :: n_sparse,n_soap, n_sites
        real(c_double), value :: cdelta
      end subroutine

      subroutine gpu_blas_mvmul_n(cubhandle, A, B, C, nAx, nAy) bind(C,name="gpu_blas_mvmul_n")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,A,B,C
        integer(c_int),value :: nAx,nAy
      end subroutine

      subroutine gpu_kernels_pow( A,B,  zeta, N,gpu_stream) bind(C,name="gpu_kernels_pow")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A,B
        type(c_ptr) :: gpu_stream
        integer(c_int),value :: N
        real(c_double), value :: zeta
      end subroutine

      subroutine gpu_axpe( A,  dccc,e0, N, gpu_stream) bind(C,name="gpu_axpc")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A
        type(c_ptr) :: gpu_stream
        integer(c_int),value :: N
        real(c_double), value :: dccc,e0
      end subroutine

      subroutine gpu_set_device(my_rank) bind(C, name="cuda_set_device")
        use iso_c_binding
        integer(c_int), value :: my_rank
      end subroutine

      subroutine gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse,gpu_stream) bind(C,name="cuda_matvect_kernels")
        use iso_c_binding
        type(c_ptr), value :: kernels_der_d, alphas_d
        type(c_ptr) :: gpu_stream
        integer(c_int), value :: n_sites,n_sparse
      end subroutine

      subroutine gpu_final_soap_forces_virial(n_sites, &
                                              Qss_d,n_soap, l_index_d, j2_index_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              n_sites0, &
                                              forces_d, &
                                              n_pairs, gpu_stream) &
                  bind(C,name="gpu_final_soap_forces_virial")
        use iso_c_binding
        type(c_ptr), value :: Qss_d, l_index_d
        type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
        type(c_ptr), value ::  j2_index_d
        type(c_ptr) :: gpu_stream
        integer(c_int), value :: n_sites,n_soap, n_sites0, n_pairs
      end subroutine

      subroutine gpu_get_sqrt_dot_p(sqrt_dot_d, soap_d, multiplicity_array_d, &
                                    cnk_d, skip_soap_component_d,  &
                                    n_sites, n_soap, n_max,l_max, gpu_stream) &
                  bind(C,name="gpu_get_sqrt_dot_p")
        use iso_c_binding
        type(c_ptr), value :: cnk_d, skip_soap_component_d
        type(c_ptr) :: gpu_stream
        type(c_ptr), value :: sqrt_dot_d, soap_d, multiplicity_array_d
        integer(c_int),value :: n_sites, n_soap, n_max,l_max
      end subroutine

      subroutine gpu_get_soap_der(soap_d, sqrt_dot_d, soap_cart_der_d, &
                                  soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, &
                                  thetas_d, phis_d, rjs_d, &
                                  multiplicity_array_d, &
                                  cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                  n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d, &
                                  n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh, gpu_stream)  &
                  bind(C,name="gpu_get_soap_der")
        use iso_c_binding
        type(c_ptr), value :: sqrt_dot_d, soap_d, soap_cart_der_d
        type(c_ptr) :: gpu_stream
        type(c_ptr), value :: soap_rad_der_d, soap_azi_der_d, soap_pol_der_d
        type(c_ptr), value :: thetas_d, phis_d, rjs_d
        type(c_ptr), value :: n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d
        type(c_ptr), value :: cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d
        type(c_ptr), value :: multiplicity_array_d
        integer(c_int),value :: n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh
      end subroutine

      subroutine gpu_soap_normalize(soap_d, sqrt_dot_d, n_soap, n_sites, gpu_stream)  &
                  bind(C,name="gpu_soap_normalize")
        use iso_c_binding
        type(c_ptr), value :: sqrt_dot_d, soap_d
        type(c_ptr) :: gpu_stream
        integer(c_int),value :: n_sites, n_soap
      end subroutine

      subroutine gpu_get_derivatives(radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d, &
                                  angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d, &
                                  cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                  rjs_d, &
                                  rcut_max, &
                                  n_atom_pairs, n_sites, n_soap, k_max, n_max, l_max, gpu_stream) &
                                  bind(C,name="gpu_get_derivatives")
        use iso_c_binding
        type(c_ptr), value :: rjs_d
        type(c_ptr) :: gpu_stream
        type(c_ptr), value :: radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d
        type(c_ptr), value :: angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d
        type(c_ptr), value :: cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d
        real(c_double), value :: rcut_max
        integer(c_int),value :: n_atom_pairs, n_sites, n_soap, k_max, n_max, l_max
      end subroutine


      subroutine gpu_get_cnk(radial_exp_coeff_d, angular_exp_coeff_d, &
                             cnk_d, &
                             n_neigh_d, k2_start_d, &
                             n_sites,  n_atom_pairs,n_soap, k_max, n_max, l_max, &
                             bintybint, &
                             atom_sigma_r_d, atom_sigma_t_d, rcut_hard_d, &
                             central_weight_d, species_d, i_beg_d, i_end_d, &
                             radial_enhancement, species_multiplicity_d, &
                             W_d, S_d, size_species_1, gpu_stream) &
                             bind(C,name="gpu_get_cnk")
        use iso_c_binding
        type(c_ptr), value :: n_neigh_d, k2_start_d
        type(c_ptr) :: gpu_stream
        type(c_ptr), value :: radial_exp_coeff_d, angular_exp_coeff_d
        type(c_ptr), value :: cnk_d
        integer(c_int),value :: n_sites,  n_atom_pairs,n_soap, k_max, n_max, l_max
        integer(c_int), value :: bintybint, radial_enhancement, size_species_1
        type(c_ptr), value :: rcut_hard_d, species_multiplicity_d
        type(c_ptr), value :: atom_sigma_r_d, atom_sigma_t_d
        type(c_ptr), value :: central_weight_d, species_d, i_beg_d, i_end_d
        type(c_ptr), value :: W_d, S_d
      end subroutine

      subroutine gpu_get_plm_array_global(plm_array_global_d, n_atom_pairs, kmax, &
                                      lmax, thetas_d, gpu_stream) &
                                      bind(C,name="gpu_get_plm_array_global")
      use iso_c_binding
      type(c_ptr), value :: plm_array_global_d, thetas_d
      type(c_ptr) :: gpu_stream
      integer(c_int),value :: n_atom_pairs, kmax, lmax
      end subroutine

      subroutine gpu_get_exp_coeff_array(eimphi_global_d, thetas_d, &
                                      rjs_d, phis_d, &
                                      mask_d, atom_simga_in_d, atom_simga_scaling_d, & 
                                      rcut, n_atom_piars, n_spceies, lmax, kmax, &
                                      prefl_array_global_d, plm_array_global_d, &
                                      plm_array_der_global_d, &
                                      prefl_array_global_der_d, &
                                      preflm_d,exp_coeff_d, &
                                      c_do_derivatives, &
                                      eimphi_rad_der_global_d, eimphi_azi_der_global_d, &
                                      plm_array_div_sin_d, plm_array_der_mul_sin_d, &
                                      exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d, &
                                      gpu_stream) & 
                                      bind(C,name="gpu_get_exp_coeff_array")
      use iso_c_binding
      type(c_ptr), value :: eimphi_global_d, rjs_d, phis_d, mask_d, exp_coeff_d
      type(c_ptr) :: gpu_stream
      type(c_ptr), value :: atom_simga_in_d, atom_simga_scaling_d, prefl_array_global_der_d
      integer(c_int), value ::  n_atom_piars, n_spceies, lmax, kmax
      real(c_double), value :: rcut
      type(c_ptr), value :: eimphi_rad_der_global_d, eimphi_azi_der_global_d, thetas_d
      type(c_ptr), value :: preflm_d, prefl_array_global_d, plm_array_global_d, plm_array_der_global_d
      logical(c_bool), value :: c_do_derivatives
      type(c_ptr), value :: exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d
      type(c_ptr), value :: plm_array_div_sin_d, plm_array_der_mul_sin_d
      end subroutine

      subroutine gpu_get_radial_exp_coeff_poly3gauss(radial_exp_coeff_d, radial_exp_coeff_der_d, &
                                      i_beg_d, i_end_d,&
                                      global_scaling_d, &
                                      size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species, &
                                      c_do_derivatives, bintybint, & 
                                      rcut_hard_d, gpu_stream) &
                                      bind(C,name="gpu_get_radial_exp_coeff_poly3gauss")
      use iso_c_binding
      type(c_ptr), value :: radial_exp_coeff_d, radial_exp_coeff_der_d, i_beg_d, i_end_d, global_scaling_d
      type(c_ptr), value :: rcut_hard_d
      type(c_ptr) :: gpu_stream
      integer(c_int), value :: size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species
      integer(c_int), value :: bintybint
      logical(c_bool), value :: c_do_derivatives
      end subroutine

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
