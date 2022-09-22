
! module load gcc/10.3.0 openblas openmpi cuda
! rm cuda_wrappers.o ; nvcc -arch=sm_70 -c src/cuda_wrappers.cu ; make clean; make -B
!  srun  --time=00:10:00 --partition=gputest --account=project_2000634 --nodes=1 --ntasks-per-node=4  --cpus-per-task=1 --gres=gpu:a100:4  ../turbogap_dev/bin/turbogap md
! gnupot
!plot "< paste trajectory_out.xyz_64k_oiriginal trajectory_out.xyz | awk 'NR>2{print $5,$13}'", x
MODULE F_B_C
    INTERFACE
      subroutine gpu_malloc_double(a_d,n) bind(C,name="cuda_malloc_double")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
        integer(c_int),value :: n
      end subroutine

      subroutine gpu_malloc_double_complex(a_d,n) bind(C,name="cuda_malloc_double_complex")
        use iso_c_binding
        implicit none
          type(c_ptr) :: a_d
          integer(c_int),value :: n
      end subroutine

      subroutine gpu_malloc_int(a_d,n) bind(C,name="cuda_malloc_int")
        use iso_c_binding
        implicit none
          type(c_ptr) :: a_d
          integer(c_int),value :: n
      end subroutine

      subroutine gpu_malloc_bool(a_d,n) bind(C,name="cuda_malloc_bool")
        use iso_c_binding
        implicit none
          type(c_ptr) :: a_d
          integer(c_int),value :: n
      end subroutine

      subroutine gpu_free(a_d) bind(C,name="cuda_free")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
      end subroutine

      subroutine cpy_double_htod(a,a_d,n) bind(C,name="cuda_cpy_double_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine


      subroutine cpy_double_complex_htod(a,a_d,n) bind(C,name="cuda_cpy_double_complex_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_int_htod(a,a_d,n) bind(C,name="cuda_cpy_int_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_bool_htod(a,a_d,n) bind(C,name="cuda_cpy_bool_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_double_dtod(b_d,a_d,n) bind(C,name="cuda_cpy_double_dtod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,b_d
        integer(c_int),value :: n
      end subroutine
      
      subroutine cpy_bool_dtoh(a_d,a,n) bind(C,name="cuda_cpy_bool_dtoh")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_double_dtoh(a_d,a,n) bind(C,name="cuda_cpy_double_dtoh")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine gpu_vector_fill_curand(a_d,n,c) bind(C,name="GPU_fill_rand")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d
        integer(c_int),value :: n,c
      end subroutine

      subroutine create_cublas_handle(cubhandle)bind(C,name="create_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle
      end subroutine

      subroutine destroy_cublas_handle(cubhandle)bind(C,name="destroy_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle
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

      subroutine gpu_kernels_pow( A,B,  zeta, N) bind(C,name="gpu_kernels_pow")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A,B
        integer(c_int),value :: N
        real(c_double), value :: zeta
      end subroutine

      subroutine gpu_axpe( A,  dccc,e0, N) bind(C,name="gpu_axpc")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A
        integer(c_int),value :: N
        real(c_double), value :: dccc,e0
      end subroutine

      subroutine one_for_all(soap, kernels, kernels_copy, Qs, energies, &
            delta, zeta, e0, &
            n_sites, n_soap, n_sparse, &
            size_kernels, size_soap, size_Qs, size_alphas, size_energies) &
            bind(C, name="wrappers_all")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: soap, kernels,kernels_copy, Qs, energies
        real(c_double), value :: delta, zeta, e0
        integer(c_int), value :: n_sites, n_soap, n_sparse
        integer(c_int), value :: size_kernels, size_soap, size_Qs, size_alphas, size_energies
      end subroutine

      subroutine gpu_set_device(my_rank) bind(C, name="cuda_set_device")
        use iso_c_binding
        integer(c_int), value :: my_rank
      end subroutine

      subroutine gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse) bind(C,name="cuda_matvect_kernels")
        use iso_c_binding
        type(c_ptr), value :: kernels_der_d, alphas_d
        integer(c_int), value :: n_sites,n_sparse
      end subroutine

      subroutine gpu_final_soap_forces_virial(n_neigh_d,n_sites,maxnn, &
                                              Qss_d,n_soap,neighbors_beg_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              neighbors_list_d, n_sites0, forces_d) &
                  bind(C,name="gpu_final_soap_forces_virial")
        use iso_c_binding
        type(c_ptr), value :: Qss_d, n_neigh_d,neighbors_beg_d, neighbors_list_d
        type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
        integer(c_int), value :: n_sites,n_soap,maxnn, n_sites0
      end subroutine

      subroutine gpu_soap_energies_forces_virial(n_neigh_d,n_sites,maxnn, &
                                              Qss_d,n_soap,neighbors_beg_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              neighbors_list_d, n_sites0, forces_d, &
                                              cuhandle, kernels_der_d, Qs_copy_d, &
                                              n_sparse, cdelta_force, &
                                              alphas_d, &
                                              kernels_d, mzetam, size_kernels, &
                                              do_forces, &
                                              energies_d, cdelta_ene, e0,  size_energies, &
                                              Qs_d, size_Qs, &
                                              kernels_copy_d, &
                                              zeta, &
                                              soap_d ) &
                  bind(C,name="gpu_soap_energies_forces_virial")
        use iso_c_binding
        type(c_ptr), value :: Qss_d, n_neigh_d,neighbors_beg_d, neighbors_list_d
        type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
        integer(c_int), value :: n_sites,n_soap,maxnn, n_sites0, size_energies
        type(c_ptr), value :: cuhandle,kernels_der_d, Qs_copy_d, alphas_d, kernels_d
        type(c_ptr), value :: energies_d, Qs_d, kernels_copy_d, soap_d
        integer(c_int),value :: n_sparse, size_kernels, size_Qs
        logical, value:: do_forces
        real(c_double), value :: cdelta_ene, cdelta_force, mzetam, e0, zeta
      end subroutine

      subroutine gpu_get_sqrt_dot_p(sqrt_dot_d, soap_d, multiplicity_array_d, &
                                    cnk_d, skip_soap_component_d,  &
                                    n_sites, n_soap, n_max,l_max) &
                  bind(C,name="gpu_get_sqrt_dot_p")
        use iso_c_binding
        type(c_ptr), value :: cnk_d, skip_soap_component_d
        type(c_ptr), value :: sqrt_dot_d, soap_d, multiplicity_array_d
        integer(c_int),value :: n_sites, n_soap, n_max,l_max
      end subroutine

      subroutine gpu_get_soap_der(soap_d, sqrt_dot_d, soap_cart_der_d, &
                                  soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, &
                                  thetas_d, phis_d, rjs_d, &
                                  multiplicity_array_d, &
                                  cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                  n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d, &
                                  n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh)  &
                  bind(C,name="gpu_get_soap_der")
        use iso_c_binding
        type(c_ptr), value :: sqrt_dot_d, soap_d, soap_cart_der_d
        type(c_ptr), value :: soap_rad_der_d, soap_azi_der_d, soap_pol_der_d
        type(c_ptr), value :: thetas_d, phis_d, rjs_d
        type(c_ptr), value :: n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d
        type(c_ptr), value :: cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d
        type(c_ptr), value :: multiplicity_array_d
        integer(c_int),value :: n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh
      end subroutine

      subroutine gpu_soap_normalize(soap_d, sqrt_dot_d, n_soap, n_sites)  &
                  bind(C,name="gpu_soap_normalize")
        use iso_c_binding
        type(c_ptr), value :: sqrt_dot_d, soap_d
        integer(c_int),value :: n_sites, n_soap
      end subroutine

    END INTERFACE
  END MODULE F_B_C
