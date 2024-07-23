# 1 "src/soap_turbo/src/soap_turbo.f90"
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2021, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# 26
module soap_turbo_desc

  use soap_turbo_radial
  use soap_turbo_angular
  use F_B_C
  use iso_c_binding
  use mpi

  contains

!**************************************************************************
  subroutine get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                      thetas, phis, alpha_max_d, alpha_max, l_max, rcut_hard_d, rcut_hard, rcut_soft_d, nf_d, global_scaling_d, &
                      atom_sigma_r_d, atom_sigma_r, &
                      atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
                      amplitude_scaling_d, radial_enhancement, central_weight_d, central_weight, basis, scaling_mode, do_timing, &
                      do_derivatives, compress_soap, compress_soap_indices, soap, soap_cart_der, time_get_soap, &
                      soap_d,  soap_cart_der_d, n_neigh_d, k2_i_site_d,cublas_handle , gpu_stream)

  implicit none

!-------------------
! Input variables
  real(c_double), intent(in), target :: rjs(:), thetas(:), phis(:)
! real(c_double), intent(in), target :: amplitude_scaling(:), atom_sigma_r_scaling(:), atom_sigma_t(:), atom_sigma_t_scaling(:)
! real(c_double), intent(in), target :: amplitude_scaling(:), atom_sigma_t(:), atom_sigma_t_scaling(:)
! real(c_double), intent(in), target :: central_weight(:), atom_sigma_r(:), global_scaling(:)
  real(c_double), intent(in), target :: central_weight(:), atom_sigma_r(:)
! real(c_double), intent(in), target :: nf(:), rcut_hard(:), rcut_soft(:)
  real(c_double), intent(in), target :: rcut_hard(:)

  integer(c_int), intent(in), target :: n_species, radial_enhancement, species(:,:), species_multiplicity(:)
  integer, intent(in) :: n_sites,l_max, n_atom_pairs, alpha_max(:), compress_soap_indices(:)
  integer(c_int), intent(in), target ::  n_neigh(:)
  logical, intent(in) :: do_derivatives, do_timing, compress_soap
  logical, intent(in) :: mask(:,:)
  logical(c_bool), allocatable, target :: new_mask(:,:)
  integer :: i_sp

  character(*), intent(in) :: basis, scaling_mode

! Output variables
  real(c_double), intent(inout), target :: soap(:,:), soap_cart_der(:,:,:)
  type(c_ptr), intent(inout) :: cublas_handle, gpu_stream
!-------------------


!-------------------
! Internal variables
  complex(c_double_complex), allocatable, target :: angular_exp_coeff(:,:), cnk(:,:,:)
  complex(c_double_complex), allocatable, target :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:)
  complex(c_double_complex), allocatable, target :: cnk_azi_der(:,:,:), cnk_rad_der(:,:,:), cnk_pol_der(:,:,:)
  complex(c_double_complex), allocatable, target :: angular_exp_coeff_pol_der(:,:)
  complex(c_double_complex), allocatable, target :: eimphi(:), prefm(:), eimphi_rad_der(:)

  real(c_double), allocatable,target, save :: W(:,:), S(:,:), multiplicity_array(:) 
  real(c_double), allocatable,target ::W_check(:,:), S_check(:,:)
  real(c_double), allocatable,target :: soap_rad_der(:,:), sqrt_dot_p(:), soap_azi_der(:,:)
  real*8, allocatable :: W_temp(:,:), S_temp(:,:)
  real(c_double), allocatable,target :: radial_exp_coeff(:,:), soap_pol_der(:,:)
  real*8, allocatable, target :: preflm(:), plm_array(:), prefl(:), fact_array(:), prefl_rad_der(:)
  real(c_double), allocatable,target :: radial_exp_coeff_der(:,:)
  real*8 :: amplitude, multiplicity, pi, rcut_max
  real*8 :: radial_time, angular_time, coeff_time, time3, total_time, soap_time, time1, time2, compress_time, &
            memory_time, basis_time

  integer(c_int), allocatable, target :: i_beg(:), i_end(:)
  integer(c_int), allocatable, target :: k2_i_site(:),k3_index(:),i_k2_start(:), k2_start(:)
  integer, save :: n_max_prev
  integer :: k_max, n_max
  integer :: i, counter, j, k, n_soap, k2, k3, n, l, m, np, counter2
  logical(c_bool), allocatable, target :: skip_soap_component(:,:,:)
  logical, allocatable, target :: do_central(:)
  logical, save :: recompute_basis = .true.
  integer(c_int) :: rank, ierr, n_multiplicity
  type(c_ptr) :: i_beg_d, i_end_d, preflm_d, plm_array_d, eimphi_d, prefm_d
  type(c_ptr) :: prefl_d, eimphi_rad_der_d, prefl_rad_der_d
  type(c_ptr) :: multiplicity_array_d, cnk_d, sqrt_dot_p_d
  type(c_ptr) :: skip_soap_component_d
  real*8 :: ttt(2)=0.d0
  real*8, intent(inout) :: time_get_soap
  integer :: maxneigh
  type(c_ptr) ::  thetas_d,phis_d,rjs_d
  type(c_ptr) ::  k2_i_site_d
  type(c_ptr) :: k3_index_d, i_k2_start_d, n_neigh_d
  type(c_ptr) :: soap_rad_der_d, soap_azi_der_d, soap_pol_der_d
  type(c_ptr) :: cnk_rad_der_d, cnk_azi_der_d,  cnk_pol_der_d
  type(c_ptr) :: radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d
  type(c_ptr) :: radial_exp_coeff_temp1_d, radial_exp_coeff_temp2_d, radial_exp_coeff_der_temp_d
  type(c_ptr) :: angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d
  type(c_ptr) :: angular_exp_coeff_pol_der_d
  type(c_ptr), intent(inout) :: soap_cart_der_d, soap_d, nf_d, rcut_hard_d, rcut_soft_d, global_scaling_d
  type(c_ptr), intent(inout) :: atom_sigma_r_d, atom_sigma_r_scaling_d,atom_sigma_t_d,atom_sigma_t_scaling_d
  type(c_ptr), intent(inout) :: amplitude_scaling_d, alpha_max_d, central_weight_d
  type(c_ptr) :: k2_start_d
  integer(c_size_t) :: st_soap, st_soap_cart_der, st_soap_rap_der
  integer(c_size_t) :: st_n_atom_pairs_int, st_n_sites_int, st_skip_component
  integer(c_size_t) :: st_n_atom_pairs_double, st_n_sites_double, st_n_multiplicity_double
  integer(c_size_t) :: st_cnk,st_cnk_rap_der, st_rad_exp_coeff_der_double, st_ang_exp_coeff_rad_der
  integer(c_size_t) :: st_mask_d, st_size_rcut_hard, st_size_atom_sigma_r,st_size_W, st_size_S
  integer(c_size_t) :: st_size_i_beg,st_size_global_scaling, st_size_i_end, st_size_angular_exp_coeff
  integer(c_size_t) :: st_size_atom_sigma_t, st_size_atom_sigma_t_scaling, st_size_species
  integer(c_size_t) :: st_size_species_multiplicity,st_size_central_weight, st_size_preflm,st_do_central_d
  integer(c_size_t) :: size_tmp_var
  
  integer(c_int) :: bintybint, size_central_weight
  type(c_ptr), save :: W_d, S_d
  logical, save :: W_S_initialized = .false.
  type(c_ptr) :: W_d_work, S_d_work, species_multiplicity_d, species_d
! type(c_ptr) :: central_weight_d, atom_sigma_r_d, atom_sigma_t_d, atom_sigma_t_scaling_d, mask_d, do_central_d
  type(c_ptr) :: mask_d, do_central_d
  integer(c_int) :: size_1_species, size_2_species, size_species_multiplicity, size_i_beg, size_i_end
  integer(c_int) :: size_rcut_hard, size_atom_sigma_t, size_atom_sigma_r, size_2_W, size_2_S,size_1_W, size_1_S
  integer(c_int) :: size_atom_sigma_t_scaling, size_global_scaling, size_alphamax
  integer(c_int) :: tmp_cint_val
! type(c_ptr) :: global_scaling_d, amplitude_scaling_d, alpha_max_d, nf_d
! type(c_ptr) :: amplitude_scaling_d, alpha_max_d
! type(c_ptr) :: alpha_max_d
  logical(c_bool) :: c_do_derivatives
  integer :: kij, ntemp, ntemp_der, mode
  integer(c_size_t) :: ntemp_d, ntemp_der_d
  integer, allocatable, target  :: k_idx(:)
  !real(c_double),allocatable, target :: amplitude_save(:),  amplitude_der_save(:)
  !type(c_ptr) :: amplitude_save_d,  amplitude_der_save_d 
  type(c_ptr) ::  all_exp_coeff_temp2_d, all_exp_coeff_temp1_d
  integer(c_size_t) :: st_size_amplitudes, st_size_exp_coeff_temp,st_size_alphamax
!-------------------


  st_soap=sizeof(soap) !n_soap*n_sites*sizeof(soap(1,1))
  call gpu_malloc_all(soap_d, st_soap, gpu_stream)

  st_soap_cart_der= sizeof(soap_cart_der) !  3*n_soap*n_atom_pairs*sizeof(soap_cart_der(1,1,1))
  call gpu_malloc_all(soap_cart_der_d, st_soap_cart_der, gpu_stream) !call gpu_malloc_all_blocking(soap_cart_der_d, st_soap_cart_der) !call gpu_malloc_all(soap_cart_der_d, st_soap_cart_der, gpu_stream)
!  stop
  if( basis == "poly3gauss" )then
    bintybint=1000
  else if( basis == "poly3" )then
    bintybint=2000
  end if

    c_do_derivatives=logical( .false., kind=c_bool ) 
    if(do_derivatives) then 
    c_do_derivatives=logical( .true., kind=c_bool )
    endif

!Total 24.1 s        
 ! call cpu_time(ttt(1))

  if( do_timing )then
    call cpu_time(time3)
  end if

!-------------------
! Constants
  pi = dacos(-1.d0)
!-------------------



  if( do_timing )then
    call cpu_time(time1)
  end if

! Check if we need to expand the central atom
  allocate( do_central(1:n_species) )
  do i = 1, n_species
    if( central_weight(i) /= 0.d0 )then
      do_central(i) = .true.
    else
      do_central(i) = .false.
    end if
  end do


  rcut_max = 0.d0
  do i = 1, n_species
    if( rcut_hard(i) > rcut_max )then
      rcut_max = rcut_hard(i)
    end if
  end do


! This assigns begin and end indices to the components of the radial basis
  allocate( i_beg(1:n_species) )
  allocate( i_end(1:n_species) ) !integer arrays
  i_beg(1) = 1
  i_end(1) = alpha_max(1)
  do i = 2, n_species
    i_beg(i) = i_end(i-1) + 1
    i_end(i) = i_beg(i) + alpha_max(i) - 1
  end do


  
  
  !23.8
  !call cpu_time(ttt(1))
  ! write(*,*)
  ! write(*,*)
  ! write(*,*)
  ! write(*,*) "Firsty firsty"

! Do some array allocation for the SOAP expansion part
  k_max = 1 + l_max*(l_max+1)/2 + l_max
  allocate( preflm(1:k_max) )
  allocate( plm_array(1:k_max) ) 
  allocate( eimphi(1:k_max) )
  allocate( prefl(0:l_max) ) 
  allocate( prefm(0:l_max) )
  allocate( fact_array(1:l_max) ) !real*8
  call get_preflm(preflm, l_max)
  st_size_preflm=k_max*sizeof(preflm(1))
  call gpu_malloc_all(preflm_d,st_size_preflm, gpu_stream)
  call cpy_htod( c_loc(preflm), preflm_d, st_size_preflm, gpu_stream)
  if( do_derivatives )then
    allocate( prefl_rad_der(0:l_max) ) ! real*8
    allocate( eimphi_rad_der(1:k_max) ) ! complex*16
  end if

  if( do_timing )then
    call cpu_time(time2)
    memory_time = time2 - time1
    time1 = time2
  end if

  ! 23.9 s
  ! call cpu_time(ttt(1))

! This is to build the radial basis
  n_max = 0
  do i = 1, n_species
    n_max = n_max + alpha_max(i)
  end do
  if( n_max_prev /= n_max )then
    n_max_prev = n_max
    recompute_basis = .true.
  end if
  !W and S array are actually never deallocated in the original code! and it use weird "save" semantic to keep allocation between iteration, but never takes care to free the memory.
  !will try to have them on the device only as they don't seem to be used on the host side now.

  if( recompute_basis )then
    if( W_S_initialized )then
      call gpu_free_async(W_d,gpu_stream)
      call gpu_free_async(S_d,gpu_stream)
    end if
    if( allocated(W) .or. allocated(S) )then
      deallocate(W, S)
    end if
    allocate( W(1:n_max, 1:n_max) )
    allocate( S(1:n_max, 1:n_max) )
    W = 0.d0
    S = 0.d0

    size_tmp_var = n_max*n_max*c_double
    call gpu_malloc_all(W_d,size_tmp_var, gpu_stream)
    call gpu_malloc_all(S_d,size_tmp_var, gpu_stream)
    tmp_cint_val = 0
!    write(0,*) tmp_cint_val
    call gpu_memset_async(S_d,tmp_cint_val,size_tmp_var, gpu_stream)
    call gpu_memset_async(W_d,0,size_tmp_var, gpu_stream)
    size_tmp_var = maxval(alpha_max)*maxval(alpha_max)*c_double
    call gpu_malloc_all(W_d_work,size_tmp_var, gpu_stream)
    call gpu_malloc_all(S_d_work,size_tmp_var, gpu_stream)
!   This is done per species. Each loop iteration modifies the slice of the W and S matrices that
!   corresponds to the species in question. The radial basis functions for species A are always
!   assumed orthogonal to the basis functions for species B. W and S are therefore block diagonal.
    do i = 1, n_species
!     We pass these temp arrays with the right size because the Lapack/Blas routines internally fail
!     if the memory is not contiguous
      allocate( S_temp(1:alpha_max(i), 1:alpha_max(i)) )
      allocate( W_temp(1:alpha_max(i), 1:alpha_max(i)) )
     
     S_temp = 0.d0
     W_temp = 0.d0
      if( basis == "poly3gauss" )then
        call get_orthonormalization_matrix_poly3gauss(alpha_max(i), atom_sigma_r(i), rcut_hard(i), S_temp, W_temp)
        !call get_orthonormalization_matrix_poly3gauss_gpu(alpha_max(i), atom_sigma_r(i), rcut_hard(i), S_d_work, W_d_work, cublas_handle, gpu_stream)
      else if( basis == "poly3" )then
        call get_orthonormalization_matrix_poly3(alpha_max(i), S_temp, W_temp)
       ! call get_orthonormalization_matrix_poly3_gpu(alpha_max(i), c_loc(S_temp), c_loc(W_temp), cublas_handle, gpu_stream)
      end if
      S(i_beg(i):i_end(i), i_beg(i):i_end(i)) = S_temp
      W(i_beg(i):i_end(i), i_beg(i):i_end(i)) = W_temp
      !not sure here. -1 should be needed for the fortran to c start indexes.
!        write(0,*) "calling copy"
      !call orthonormalization_copy_to_global_matrix(S_d_work,S_d,alpha_max(i),i_beg(i)-1,n_max,gpu_stream)
!        write(0,*) "calling done"
      !call orthonormalization_copy_to_global_matrix(W_d_work,W_d,alpha_max(i),i_beg(i)-1,n_max,gpu_stream)
!        write(0,*) "second calling done"
!      deallocate( S_temp, W_temp )
    end do
    call gpu_free_async(S_d_work,gpu_stream)
    call gpu_free_async(W_d_work,gpu_stream)
  end if


  ttt(1)=MPI_Wtime() 
  size_tmp_var = n_max*n_max*c_double
  call cpy_htod_blocking(c_loc(S),S_d,size_tmp_var)
  call cpy_htod_blocking(c_loc(W),W_d,size_tmp_var)
  ! call cpy_htod(c_loc(S),S_d,size_tmp_var,gpu_stream)
  ! call cpy_htod(c_loc(W),W_d,size_tmp_var,gpu_stream)
  ! allocate( W_check(1:n_max, 1:n_max) )
  ! allocate( S_check(1:n_max, 1:n_max) )
  ! size_tmp_var = n_max*n_max*c_double
  ! call cpy_dtoh_blocking(S_d,c_loc(S_check), size_tmp_var)
  ! call cpy_dtoh_blocking(W_d,c_loc(W_check), size_tmp_var)
  ! do j=1,n_max
  !   do i=1,n_max
  !     write(*,*) i,j, S(i,j)-S_check(i,j), W(i,j)-W_check(i,j)
  !     enddo
  ! enddo

  ! stop

  write(*,*) ' after copying S and W'
  call gpu_device_sync()
  stop 
  
  if( do_timing )then
    call cpu_time(time2)
    basis_time = time2 - time1
  end if

  ! 23.7
  ! call cpu_time(ttt(1))

! This is for the expansion coefficients and the soap vectors
  allocate( radial_exp_coeff(1:n_max, 1:n_atom_pairs) )
  allocate( angular_exp_coeff(1:k_max, 1:n_atom_pairs) )
  radial_exp_coeff=0.d0
  angular_exp_coeff = 0.d0
  allocate( cnk( 1:k_max, 1:n_max, 1:n_sites) )
  cnk = 0.d0
! Handle SOAP compression here
  allocate( skip_soap_component(0:l_max, 1:n_max, 1:n_max) )
  skip_soap_component = .false.
  if( compress_soap )then
    if( do_timing )then
      call cpu_time(time1)
    end if
    n_soap = size(compress_soap_indices)
    skip_soap_component = .true.
    counter = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          counter = counter + 1
          do i = 1, n_soap
            if( compress_soap_indices(i) == counter )then
              skip_soap_component(l, np, n) = .false.
              exit
            end if
          end do
        end do
      end do
    end do
    if( do_timing )then
      call cpu_time(time2)
      compress_time = time2 - time1
    end if
  else
    n_soap = n_max*(n_max+1)/2 * (l_max+1)
  end if

  if( do_timing )then
    call cpu_time(time1)
  end if

  allocate( sqrt_dot_p(1:n_sites) )

  st_soap_rap_der=n_soap*n_atom_pairs*sizeof(soap_rad_der(1,1))
  st_n_atom_pairs_int=n_atom_pairs*sizeof(k3_index(1))
  st_n_sites_int=n_sites*sizeof(n_neigh(1))
  st_n_sites_double=n_sites*sizeof(sqrt_dot_p(1))
  st_n_atom_pairs_double=n_atom_pairs*sizeof(thetas(1))
  st_skip_component=sizeof(skip_soap_component)

  sqrt_dot_p = 0.d0
  if( do_derivatives )then
    allocate( radial_exp_coeff_der(1:n_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_rad_der(1:k_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_azi_der(1:k_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_pol_der(1:k_max, 1:n_atom_pairs) )
    allocate( cnk_rad_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( cnk_azi_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( cnk_pol_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( soap_rad_der(1:n_soap, 1:n_atom_pairs) )
    allocate( soap_azi_der(1:n_soap, 1:n_atom_pairs) )
    allocate( soap_pol_der(1:n_soap, 1:n_atom_pairs) )
    radial_exp_coeff_der = 0.d0
    angular_exp_coeff_rad_der = 0.d0
    angular_exp_coeff_azi_der = 0.d0
    angular_exp_coeff_pol_der = 0.d0
    cnk_rad_der = 0.d0
    soap_rad_der = 0.d0
    cnk_azi_der = 0.d0
    soap_azi_der = 0.d0
    cnk_pol_der = 0.d0
    soap_pol_der = 0.d0
  else
!   We need this dummy variable defined here. Note the decreased range for the second index
    allocate( radial_exp_coeff_der(1:n_max, 1:1) )
  end if

  !write(*,*) n_sites, n_atom_pairs
  if( do_timing )then
    call cpu_time(time2)
    memory_time = memory_time + time2 - time1
    time1 = time2
  end if
  
  ! 18.6
  !call cpu_time(ttt(1))

  ! write(*,*)
  ! write(*,*) 
  ! write(*,*)
  ! write(*,*) alpha_max
  ! stop

  allocate(k_idx(1:n_sites))
  k_idx(1) = 0

  do i = 2, n_sites
     k_idx(i) = k_idx(i-1) + n_neigh(i-1)
  end do
  
  allocate(new_mask(1:n_atom_pairs,1:n_species))
  new_mask= logical( .false., kind=c_bool ) !.false.
  do i_sp=1,n_species
    do j=1,n_atom_pairs
      if(mask(j,i_sp)) then
        new_mask(j,i_sp)= logical( .true., kind=c_bool ) !.true.
      endif
    enddo
  enddo
  st_mask_d= n_atom_pairs*n_species*sizeof(new_mask(1,1))
  call gpu_malloc_all(mask_d, st_mask_d, gpu_stream)
  call cpy_htod(c_loc(new_mask), mask_d, st_mask_d, gpu_stream)

  st_do_central_d= n_species*sizeof(do_central(1))
  call gpu_malloc_all(do_central_d, st_do_central_d, gpu_stream)
  call cpy_htod(c_loc(do_central), do_central_d, st_do_central_d, gpu_stream)

  size_rcut_hard=size(rcut_hard,1)
  st_size_rcut_hard=size_rcut_hard*sizeof(rcut_hard(1))
! call gpu_malloc_all(rcut_hard_d,st_size_rcut_hard, gpu_stream)
! call cpy_htod(c_loc(rcut_hard),rcut_hard_d,st_size_rcut_hard, gpu_stream)
! call gpu_malloc_all(rcut_soft_d,st_size_rcut_hard, gpu_stream)
! call cpy_htod(c_loc(rcut_soft),rcut_soft_d,st_size_rcut_hard, gpu_stream)
! size_rcut_hard=size(rcut_hard,1)
! call gpu_malloc_all(nf_d,st_size_rcut_hard, gpu_stream)
! call cpy_htod(c_loc(nf),nf_d,st_size_rcut_hard, gpu_stream)

  size_atom_sigma_r=size(atom_sigma_r,1)
  st_size_atom_sigma_r=size_atom_sigma_r*sizeof(atom_sigma_r(1))
! call gpu_malloc_all(atom_sigma_r_d,st_size_atom_sigma_r, gpu_stream)
! call cpy_htod(c_loc(atom_sigma_r),atom_sigma_r_d,st_size_atom_sigma_r, gpu_stream)
! call gpu_malloc_all(atom_sigma_r_scaling_d,st_size_atom_sigma_r, gpu_stream)
! call cpy_htod(c_loc(atom_sigma_r_scaling),atom_sigma_r_scaling_d,st_size_atom_sigma_r, gpu_stream)
! call gpu_malloc_all(amplitude_scaling_d,st_size_atom_sigma_r, gpu_stream)
! call cpy_htod(c_loc(amplitude_scaling),amplitude_scaling_d,st_size_atom_sigma_r, gpu_stream)
! size_alphamax=size(alpha_max,1)
! st_size_alphamax=size_alphamax*sizeof(alpha_max(1))
! call gpu_malloc_all(alpha_max_d,st_size_alphamax, gpu_stream)
! call cpy_htod(c_loc(alpha_max),alpha_max_d,st_size_alphamax, gpu_stream)

  ntemp = maxval(alpha_max)
  ntemp_der = ntemp
  if (do_derivatives) ntemp = ntemp+2
  ntemp_d=ntemp*n_atom_pairs*sizeof(radial_exp_coeff_der(1,1))
  ntemp_der_d=ntemp_der*n_atom_pairs*sizeof(radial_exp_coeff_der(1,1))
  st_rad_exp_coeff_der_double=n_max*n_atom_pairs*sizeof(radial_exp_coeff_der(1,1))
  call gpu_malloc_all(radial_exp_coeff_d, st_rad_exp_coeff_der_double, gpu_stream)
  call gpu_memset_async(radial_exp_coeff_d,0, st_rad_exp_coeff_der_double, gpu_stream)
  !call cpy_htod(c_loc(radial_exp_coeff),radial_exp_coeff_d, st_rad_exp_coeff_der_double, gpu_stream)
  call gpu_malloc_all(radial_exp_coeff_temp1_d, ntemp_d, gpu_stream)
  call gpu_malloc_all(radial_exp_coeff_temp2_d, ntemp_d, gpu_stream)
  call gpu_malloc_all(radial_exp_coeff_der_d, st_rad_exp_coeff_der_double, gpu_stream)
  call gpu_memset_async(radial_exp_coeff_der_d,0, st_rad_exp_coeff_der_double, gpu_stream)
  !call cpy_htod(c_loc(radial_exp_coeff_der),radial_exp_coeff_der_d, st_rad_exp_coeff_der_double, gpu_stream)
  call gpu_malloc_all(radial_exp_coeff_der_temp_d, ntemp_der_d , gpu_stream) 

  size_i_beg=size(i_beg,1)
  st_size_i_beg=size_i_beg*sizeof(i_beg(1))
  call gpu_malloc_all(i_beg_d, st_size_i_beg, gpu_stream)
  call cpy_htod(c_loc(i_beg),i_beg_d,st_size_i_beg, gpu_stream)

  size_i_end=size(i_end,1)
  st_size_i_end=size_i_end*sizeof(i_end(1))
  call gpu_malloc_all(i_end_d, st_size_i_end, gpu_stream)
  call cpy_htod(c_loc(i_end),i_end_d,st_size_i_end, gpu_stream)

!  size_1_W=size(W,1)
!  size_2_W=size(W,2)
!  st_size_W=size_1_W*size_2_W*sizeof(W(1,1))
!  call gpu_malloc_all(W_d,st_size_W, gpu_stream)
!  call cpy_htod(c_loc(W),W_d,st_size_W, gpu_stream)

! size_central_weight=size(central_weight,1)
! st_size_central_weight= size_central_weight*sizeof(central_weight(1))
! call gpu_malloc_all(central_weight_d,st_size_central_weight,gpu_stream)
! call cpy_htod(c_loc(central_weight),central_weight_d,st_size_central_weight, gpu_stream)

  call gpu_malloc_all(rjs_d,st_n_atom_pairs_double, gpu_stream)
  call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)

  call gpu_malloc_all(n_neigh_d,st_n_sites_int,gpu_stream)
  call cpy_htod(c_loc(n_neigh),n_neigh_d, st_n_sites_int,gpu_stream)

  mode = 0
  if( scaling_mode == "polynomial" ) mode =1
  if( basis == "poly3gauss" )then
    call gpu_radial_poly3gauss(n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, ntemp,&
                               c_do_derivatives, radial_exp_coeff_d, radial_exp_coeff_der_d, rcut_soft_d, atom_sigma_r_d, &
                               radial_exp_coeff_temp1_d, radial_exp_coeff_temp2_d, radial_exp_coeff_der_temp_d, i_beg_d, &
                               i_end_d, atom_sigma_r_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, &
                               nf_d, ntemp_der, W_d, gpu_stream) 
!     call get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
!                                                       rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
!                                                       amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
!                                                       scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
!                                                       radial_exp_coeff(i_beg(i):i_end(i), :), &
!                                                       radial_exp_coeff_der(i_beg(i):i_end(i), :), k_idx(1:n_sites))
  else if( basis == "poly3" )then
    call gpu_radial_poly3(n_atom_pairs, n_species, mask_d, rjs_d, rcut_hard_d, n_sites, n_neigh_d, n_max, ntemp,&
                          c_do_derivatives, radial_exp_coeff_d, radial_exp_coeff_der_d, rcut_soft_d, atom_sigma_r_d, &
                          radial_exp_coeff_temp1_d, radial_exp_coeff_temp2_d, radial_exp_coeff_der_temp_d, i_beg_d, &
                          i_end_d, atom_sigma_r_scaling_d, mode, radial_enhancement, amplitude_scaling_d, alpha_max_d, &
                          nf_d, ntemp_der, W_d, do_central_d, central_weight_d, gpu_stream) 
!     call get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
!                                                  rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
!                                                  amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
!                                                  scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
!                                                  do_central(i), central_weight(i), &
!                                                  radial_exp_coeff(i_beg(i):i_end(i), :), &
!                                                  radial_exp_coeff_der(i_beg(i):i_end(i), :) )
  end if
  
  ! do kij=1,n_atom_pairs
  ! do n=1,n_max
  ! if(isnan(radial_exp_coeff(n,kij))) then
  ! write(*,*)
  ! write(*,*) "Nan allert!",n,kij
  ! !stop
  ! endif
  ! enddo
  ! enddo
  ! ttt(1)=MPI_Wtime() 
  allocate(k2_i_site(1:n_atom_pairs))
  !write(*,*) "N atom pairs", n_atom_pairs
  allocate(k2_start(1:n_sites))
  
  k2 = 0
  do i = 1, n_sites
   k2_start(i)=k2
    do j = 1, n_neigh(i)
      k2 = k2 + 1
      k2_i_site(k2)=i
    enddo
  enddo
  call gpu_malloc_all(k2_start_d,st_n_sites_int,gpu_stream)
  call cpy_htod(c_loc(k2_start), k2_start_d, st_n_sites_int, gpu_stream)


  call gpu_malloc_all(k2_i_site_d,st_n_atom_pairs_int,gpu_stream)
  call cpy_htod(c_loc(k2_i_site),k2_i_site_d, st_n_atom_pairs_int, gpu_stream)
  
  !st_size_amplitudes=n_sites*
  !st_size_exp_coeff_temp=size(k_idx,1)*
  ! write(*,*)
  ! write(*,*) i_beg
  ! write(*,*) i_end
  ! write(*,*)
  ! call gpu_malloc_all(amplitude_save_d, st_size_amplitudes, gpu_stream)
  ! call gpu_malloc_all(amplitude_der_save_d, st_size_amplitudes, gpu_stream)

!  st_rad_exp_coeff_der_double=n_max*n_atom_pairs*sizeof(radial_exp_coeff_der(1,1))
!  call gpu_malloc_all(radial_exp_coeff_der_d, st_rad_exp_coeff_der_double , gpu_stream) !call gpu_malloc_all(radial_exp_coeff_der_d, st_rad_exp_coeff_der_double, gpu_stream)
!  call cpy_htod(c_loc(radial_exp_coeff_der),radial_exp_coeff_der_d, &
!               st_rad_exp_coeff_der_double, gpu_stream)


!  call gpu_malloc_all(radial_exp_coeff_d, st_rad_exp_coeff_der_double, gpu_stream)
!  call cpy_htod(c_loc(radial_exp_coeff),radial_exp_coeff_d, & 
!                    st_rad_exp_coeff_der_double, gpu_stream)
  
! size_rcut_hard=size(rcut_hard,1)
! st_size_rcut_hard=size_rcut_hard*sizeof(rcut_hard(1))
! call gpu_malloc_all(rcut_hard_d,st_size_rcut_hard, gpu_stream)
! call cpy_htod(c_loc(rcut_hard),rcut_hard_d,st_size_rcut_hard, gpu_stream)

! size_i_beg=size(i_beg,1)
! st_size_i_beg=size_i_beg*sizeof(i_beg(1))
! call gpu_malloc_all(i_beg_d, st_size_i_beg, gpu_stream)
! call cpy_htod(c_loc(i_beg),i_beg_d,st_size_i_beg, gpu_stream)

! size_i_end=size(i_end,1)
! st_size_i_end=size_i_end*sizeof(i_end(1))
! call cpy_htod(c_loc(i_end),i_end_d,st_size_i_end, gpu_stream)
! call cpy_htod(c_loc(i_end),i_end_d,st_size_i_end, gpu_stream)
  
! size_global_scaling=size(global_scaling,1)
! st_size_global_scaling=size_global_scaling*sizeof(global_scaling(1))
! call gpu_malloc_all(global_scaling_d, st_size_global_scaling, gpu_stream)
! call cpy_htod(c_loc(global_scaling), global_scaling_d, st_size_global_scaling, gpu_stream)

  
! allocate(new_mask(1:n_atom_pairs,1:n_species))
! new_mask= logical( .false., kind=c_bool ) !.false.
! do i_sp=1,n_species
!   do j=1,n_atom_pairs
!     if(mask(j,i_sp)) then
!       new_mask(j,i_sp)= logical( .true., kind=c_bool ) !.true.
!     endif
!   enddo
! enddo
! st_mask_d= n_atom_pairs*n_species*sizeof(new_mask(1,1))
! call gpu_malloc_all(mask_d, st_mask_d, gpu_stream)
! call cpy_htod(c_loc(new_mask), mask_d, st_mask_d, gpu_stream)


! size_atom_sigma_r=size(atom_sigma_r,1)
! st_size_atom_sigma_r=size_atom_sigma_r*sizeof(atom_sigma_r(1))
! call gpu_malloc_all(atom_sigma_r_d,st_size_atom_sigma_r, gpu_stream)
! call cpy_htod(c_loc(atom_sigma_r),atom_sigma_r_d,st_size_atom_sigma_r, gpu_stream)

! size_1_W=size(W,1)
! size_2_W=size(W,2)
! st_size_W=size_1_W*size_2_W*sizeof(W(1,1))
! call gpu_malloc_all(W_d,st_size_W, gpu_stream)
! call cpy_htod(c_loc(W),W_d,st_size_W, gpu_stream)

!  size_1_S=size(S,1)
!  size_2_S=size(S,2)
!  st_size_S=size_1_S*size_2_S*sizeof(S(1,1))
!  call gpu_malloc_all(S_d,st_size_S, gpu_stream)
!  call cpy_htod(c_loc(S),S_d,st_size_S, gpu_stream)
  !  write(*,*) "RoRo ",global_scaling(:), rcut_hard(:)
  if( do_timing )then
    call cpu_time(time2)
    radial_time = time2 - time1
    time1 = time2
  end if
!  stop


!   call gpu_malloc_all(rjs_d,st_n_atom_pairs_double, gpu_stream)
!   call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)

    call gpu_malloc_all(phis_d,st_n_atom_pairs_double, gpu_stream)
    call cpy_htod(c_loc(phis),phis_d, st_n_atom_pairs_double, gpu_stream) 
    call gpu_malloc_all(thetas_d,st_n_atom_pairs_double, gpu_stream)
    call cpy_htod(c_loc(thetas),thetas_d, st_n_atom_pairs_double, gpu_stream)

  st_size_angular_exp_coeff=k_max*n_atom_pairs*sizeof(angular_exp_coeff(1,1))
  call gpu_malloc_all(angular_exp_coeff_d,st_size_angular_exp_coeff,gpu_stream) 
  
  st_ang_exp_coeff_rad_der=k_max*n_atom_pairs*sizeof(angular_exp_coeff_rad_der(1,1))
  call gpu_malloc_all(angular_exp_coeff_rad_der_d,st_ang_exp_coeff_rad_der,gpu_stream)
  call gpu_malloc_all(angular_exp_coeff_azi_der_d,st_ang_exp_coeff_rad_der,gpu_stream)
  call gpu_malloc_all(angular_exp_coeff_pol_der_d,st_ang_exp_coeff_rad_der,gpu_stream)
  
! size_atom_sigma_t=size(atom_sigma_t,1)
! st_size_atom_sigma_t=size_atom_sigma_t*sizeof(atom_sigma_t(1))
! call gpu_malloc_all(atom_sigma_t_d,st_size_atom_sigma_t,gpu_stream)
! call cpy_htod(c_loc(atom_sigma_t),atom_sigma_t_d,st_size_atom_sigma_t, gpu_stream)
  
! size_atom_sigma_t_scaling=size(atom_sigma_t_scaling,1)
! st_size_atom_sigma_t_scaling=size_atom_sigma_t_scaling*sizeof(atom_sigma_t_scaling(1))
! call gpu_malloc_all(atom_sigma_t_scaling_d,st_size_atom_sigma_t_scaling,gpu_stream)
! call cpy_htod(c_loc(atom_sigma_t_scaling),atom_sigma_t_scaling_d, & 
!                                              st_size_atom_sigma_t, gpu_stream)

  call gpu_get_radial_exp_coeff_poly3gauss(radial_exp_coeff_d, radial_exp_coeff_der_d, &
                                i_beg_d, i_end_d, &
                                global_scaling_d, &
                                n_max, n_atom_pairs, n_species, &
                                c_do_derivatives, bintybint, &
                                rcut_hard_d, &
                                k2_i_site_d, k2_start_d, &
                                gpu_stream)
! For the angular expansion the masking works differently, since we do not have a species-augmented basis as in the
! radial expansion part.
! call get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, atom_sigma_t, atom_sigma_t_scaling, &
  call get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, &
                                          rcut_max, l_max, eimphi, preflm, plm_array, prefl, prefm, &
                                          fact_array, mask, n_species, eimphi_rad_der, &
                                          do_derivatives, prefl_rad_der, angular_exp_coeff, angular_exp_coeff_rad_der, &
                                          angular_exp_coeff_azi_der, angular_exp_coeff_pol_der, &
                                          angular_exp_coeff_d, angular_exp_coeff_rad_der_d, &
                                          angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d, n_atom_pairs, &
                                          thetas_d, preflm_d, rjs_d, phis_d, mask_d,  &
                                          atom_sigma_t_d, atom_sigma_t_scaling_d,gpu_stream)
  if( do_timing )then
    call cpu_time(time2)
    angular_time = time2 - time1
  end if
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"

! 14 s
!  call cpu_time(ttt(1)

!! For debugging (only gfortran)
!  if( .false. )then
!    write(*,*) "# Gaussian centered at ", rjs(4)
!    write(*,"(A)",advance="no") "rho(x) = "
!    do i = 1, alpha_max
!      write(*,"(A,I0,A,E16.8,A)",advance="no") "p", i, "(x) *", radial_exp_coeff(i,4), "+"
!    end do
!    write(*,*) "0."
!  end if




!  write(*,*) "Obtaining full expansion coefficients..."
!  write(*,*)
  if( do_timing )then
    call cpu_time(time1)
  end if
 

  st_cnk=k_max*n_max*n_sites*sizeof(cnk(1,1,1))

  call gpu_malloc_all(cnk_d,st_cnk,gpu_stream) ! call gpu_malloc_double_complex(cnk_d, k_max*n_max*n_sites)
   
  

! call gpu_malloc_all(n_neigh_d,st_n_sites_int,gpu_stream)
! call cpy_htod(c_loc(n_neigh),n_neigh_d, st_n_sites_int,gpu_stream)

  size_1_species=size(species,1)
  size_2_species=size(species,2)
  st_size_species=size_1_species*size_2_species*sizeof(species(1,1))
  call gpu_malloc_all(species_d, st_size_species,gpu_stream)
  call cpy_htod(c_loc(species),species_d,st_size_species, gpu_stream)

  size_species_multiplicity=size(species_multiplicity,1)
  st_size_species_multiplicity=size_species_multiplicity*sizeof(species_multiplicity(1))
  call gpu_malloc_all(species_multiplicity_d,st_size_species_multiplicity,gpu_stream)
  call cpy_htod(c_loc(species_multiplicity), species_multiplicity_d, &
                                     st_size_species_multiplicity,gpu_stream)
  
! size_central_weight=size(central_weight,1)
! st_size_central_weight= size_central_weight*sizeof(central_weight(1))
! call gpu_malloc_all(central_weight_d,st_size_central_weight,gpu_stream)
! call cpy_htod(c_loc(central_weight),central_weight_d,st_size_central_weight, gpu_stream)

  call gpu_get_cnk(radial_exp_coeff_d, angular_exp_coeff_d, &
                             cnk_d, &
                             n_neigh_d, k2_start_d, &
                             n_sites,  n_atom_pairs,n_soap, k_max, n_max, l_max, &
                             bintybint, & 
                             atom_sigma_r_d, atom_sigma_t_d, rcut_hard_d, &
                             central_weight_d, species_d, i_beg_d, i_end_d, &
                             radial_enhancement, species_multiplicity_d, &
                             W_d, S_d, size_1_species, gpu_stream)


    !call cpy_htod(c_loc(soap), soap_d, st_soap, gpu_stream)
  call gpu_memset_async(soap_d, 0, st_soap, gpu_stream)
! 12.6 s
  !call cpu_time(ttt(1)) 
  ! ttt(1)=MPI_Wtime()

! Do derivatives
  if( do_derivatives )then

  pi = dacos(-1.d0)
    
    st_cnk_rap_der=k_max*n_max*n_atom_pairs*sizeof(cnk_rad_der(1,1,1))
    call gpu_malloc_all(cnk_rad_der_d,st_cnk_rap_der,gpu_stream)
    call gpu_malloc_all(cnk_azi_der_d,st_cnk_rap_der,gpu_stream)
    call gpu_malloc_all(cnk_pol_der_d,st_cnk_rap_der,gpu_stream)
    ! call gpu_malloc_double_complex(cnk_rad_der_d,k_max*n_max*n_atom_pairs)
    ! call gpu_malloc_double_complex(cnk_azi_der_d,k_max*n_max*n_atom_pairs)
    ! call gpu_malloc_double_complex(cnk_pol_der_d,k_max*n_max*n_atom_pairs)

    call gpu_get_derivatives(radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d, &
                                  angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d, &
                                  cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                                  rjs_d, &
                                  rcut_max, &
                                  n_atom_pairs, n_sites,  n_soap, k_max, n_max, l_max, gpu_stream) 


  end if

  if( do_timing )then
    call cpu_time(time2)
    coeff_time = time2 - time1
  end if




! 7.5 s 

!  write(*,*) "Building SOAP vectors..."

  if( do_timing )then
    call cpu_time(time1)
  end if

! Create the multiplicity array
  if( recompute_basis )then
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          if( skip_soap_component(l, np, n) )cycle
          do m = 0, l
            counter2 = counter2 + 1
          end do
        end do
      end do
    end do
    if( allocated(multiplicity_array) )then
      deallocate(multiplicity_array)
    end if
    allocate( multiplicity_array(1:counter2) )

    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          if( skip_soap_component(l, np, n) )cycle
          do m = 0, l
            counter2 = counter2 + 1
            multiplicity = 1.d0
            if( n /= np )then
              multiplicity = multiplicity * dsqrt(2.d0)
            end if
            if( m > 0 )then
              multiplicity = multiplicity * 2.d0
            end if
            multiplicity_array(counter2) = multiplicity
          end do
        end do
      end do
    end do
  end if
  recompute_basis = .false.

  n_multiplicity=size(multiplicity_array,1)
  st_n_multiplicity_double=n_multiplicity*sizeof(multiplicity_array(1))
  
  call gpu_malloc_all(multiplicity_array_d, st_n_multiplicity_double,gpu_stream)
  call cpy_htod(c_loc( multiplicity_array), multiplicity_array_d, st_n_multiplicity_double, gpu_stream)
  
  call gpu_malloc_all(sqrt_dot_p_d, st_n_sites_double,gpu_stream)
  call gpu_malloc_all(skip_soap_component_d, st_skip_component,gpu_stream)
  call cpy_htod(c_loc(skip_soap_component), skip_soap_component_d, st_skip_component,gpu_stream)
  
  
  call gpu_get_sqrt_dot_p(sqrt_dot_p_d, soap_d, multiplicity_array_d, &
                                    cnk_d, skip_soap_component_d,  &
                                    n_sites, n_soap, n_max,l_max, gpu_stream)


  if( do_derivatives )then
!   Derivatives of the SOAP descriptor in spherical coordinates
!****************************
! Uncomment for detailed timing check

! call cpu_time(time1)
!****************************
    
  allocate(k3_index(1:n_atom_pairs))
  allocate(i_k2_start(1:n_sites))

    k2 = 0
    maxneigh=0
    do i = 1, n_sites
      ! maxneigh=maxneigh+n_neigh(i)
      if(n_neigh(i)>maxneigh) maxneigh=n_neigh(i)
      do j = 1, n_neigh(i)
        k2 = k2 + 1
!       Transform to Cartesian
        if( j == 1 )then
          k3_index(k2) = k2
          i_k2_start(i)=k2
          else
          k3_index(k2)=k3_index(k2-1)
        end if
      end do
    end do
    


    call gpu_malloc_all(i_k2_start_d,st_n_sites_int,gpu_stream)
    call cpy_htod(c_loc(i_k2_start),i_k2_start_d,st_n_sites_int, gpu_stream)

    call gpu_malloc_all(k3_index_d,st_n_atom_pairs_int,gpu_stream)
    call cpy_htod(c_loc(k3_index),k3_index_d, st_n_atom_pairs_int, gpu_stream)
    !call cpy_htod(c_loc(soap_cart_der),soap_cart_der_d, st_soap_cart_der)
    call gpu_malloc_all(soap_rad_der_d, st_soap_rap_der,gpu_stream)
    call gpu_malloc_all(soap_azi_der_d, st_soap_rap_der,gpu_stream)
    call gpu_malloc_all(soap_pol_der_d, st_soap_rap_der,gpu_stream)

    ! call cpy_htod(c_loc(soap_rad_der),soap_rad_der_d, st_soap_rap_der)
    ! call cpy_htod(c_loc(soap_azi_der),soap_azi_der_d, st_soap_rap_der)
    ! call cpy_htod(c_loc(soap_pol_der),soap_pol_der_d, st_soap_rap_der)

    call gpu_get_soap_der(soap_d, sqrt_dot_p_d, soap_cart_der_d, &
                          soap_rad_der_d, soap_azi_der_d, soap_pol_der_d, &
                          thetas_d,phis_d,rjs_d, & 
                          multiplicity_array_d, &
                          cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d, &
                          n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, skip_soap_component_d, &
                          n_sites, n_atom_pairs, n_soap, k_max, n_max, l_max, maxneigh, gpu_stream) 
    
    !call cpy_dtoh(soap_cart_der_d,c_loc(soap_cart_der), st_soap_cart_der)

  end if
! Now we normalize the soap vectors:
  call gpu_soap_normalize(soap_d, sqrt_dot_p_d, n_soap, n_sites, gpu_stream)
  !call cpy_dtoh(soap_d, c_loc(soap), st_soap)


! This is for debugging
 if( .false. )then
    
    call cpy_dtoh(soap_rad_der_d, c_loc(soap_rad_der),st_soap, gpu_stream)
    call cpy_dtoh(soap_azi_der_d, c_loc(soap_azi_der),st_soap, gpu_stream)
    call cpy_dtoh(soap_pol_der_d, c_loc(soap_pol_der),st_soap, gpu_stream)
    call cpy_dtoh(sqrt_dot_p_d, c_loc(sqrt_dot_p), st_soap, gpu_stream) ! blocking or asynchronous?

   open(unit=10, file="gpu_soap_desc.dat", status="unknown")

   do i = 1, 23 
    write(10, *) soap(1:n_soap, i)
   enddo
   open(unit=10, file="gpu_soap_rad_der.dat", status="unknown", position="append")
   do i = 1, 23 
   write(10, *) soap_rad_der(1:n_soap, i)
   enddo
   close(10)
   open(unit=10, file="gpu_soap_azi_der.dat", status="unknown", position="append")
   do i = 1, 23 
   write(10, *) soap_azi_der(1:n_soap, i)
   enddo
   close(10)
   open(unit=10, file="gpu_soap_pol_der.dat", status="unknown", position="append")
   do i = 1, 23 
   write(10, *) soap_pol_der(1:n_soap, i)
   enddo
   close(10)
   open(unit=10, file="gpu_soap_cart_der.dat", status="unknown")
   do i = 1, 23  !431+1, 431+23
     write(10, *) soap_cart_der(1, 1:n_soap, i)
     write(10, *) soap_cart_der(2, 1:n_soap, i)
     write(10, *) soap_cart_der(3, 1:n_soap, i)
     write(10, *)
   end do
   close(10)
   open(unit=10, file="gpu_norm.dat", status="unknown", position="append")
   write(10, *) sqrt_dot_p(317)
   close(10)
  stop

  end if
!  call cpu_time(time2)
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"


  deallocate( eimphi, preflm, plm_array, prefl, prefm, fact_array, radial_exp_coeff, angular_exp_coeff, cnk, &
              i_beg, i_end, do_central, sqrt_dot_p, k2_i_site)
  
  if( do_derivatives )then
    deallocate( radial_exp_coeff_der, angular_exp_coeff_rad_der, soap_rad_der, soap_azi_der, soap_pol_der, &
                cnk_rad_der, cnk_azi_der, cnk_pol_der, eimphi_rad_der, angular_exp_coeff_azi_der, &
                prefl_rad_der, angular_exp_coeff_pol_der,k3_index)
  call gpu_free_async(soap_rad_der_d,gpu_stream)
  call gpu_free_async(soap_azi_der_d,gpu_stream)
  call gpu_free_async(soap_pol_der_d,gpu_stream)
  call gpu_free_async(cnk_rad_der_d,gpu_stream)
  call gpu_free_async(cnk_azi_der_d,gpu_stream)
  call gpu_free_async(cnk_pol_der_d,gpu_stream)
  call gpu_free_async(k3_index_d,gpu_stream)
  call gpu_free_async(i_k2_start_d,gpu_stream)
  call gpu_free_async(radial_exp_coeff_d,gpu_stream)!call gpu_free_async(radial_exp_coeff_der_d,gpu_stream)
  call gpu_free_async(radial_exp_coeff_temp1_d,gpu_stream)!call gpu_free_async(radial_exp_coeff_der_d,gpu_stream)
  call gpu_free_async(radial_exp_coeff_temp2_d,gpu_stream)!call gpu_free_async(radial_exp_coeff_der_d,gpu_stream)
  call gpu_free_async(angular_exp_coeff_d,gpu_stream)
  call gpu_free_async(angular_exp_coeff_rad_der_d,gpu_stream)
  call gpu_free_async(angular_exp_coeff_azi_der_d,gpu_stream)
  call gpu_free_async(angular_exp_coeff_pol_der_d,gpu_stream)

  else
    deallocate( radial_exp_coeff_der )
    call gpu_free_async(radial_exp_coeff_der_d,gpu_stream)
    call gpu_free_async(radial_exp_coeff_der_temp_d,gpu_stream)
  end if

  if( do_timing )then
    call cpu_time(time2)
    soap_time = time2-time1
    total_time = time2 - time3
    write(*,*)'                                       |'
    write(*,*)'SOAP timings:                          |'
    write(*,*)'                                       |'
    write(*,'(A, F8.3, A)') '  *) Radial expansion: ', radial_time, ' seconds |'
    write(*,'(A, F6.3, A)') '  *) Radial basis build: ', basis_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Angular expansion: ', angular_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Expansion coeffs.: ', coeff_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) SOAP vector build: ', soap_time, ' seconds |'
    if( compress_soap )then
      write(*,'(A, F8.3, A)') '  *) SOAP compression: ', compress_time, ' seconds |'
    end if
    write(*,'(A, F7.3, A)') '  *) Memory allocation: ', memory_time, ' seconds |'
    write(*,'(A, F19.3, A)') '  *) Total: ', total_time, ' seconds |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if
  call gpu_free_async(thetas_d,gpu_stream)
  call gpu_free_async(phis_d,gpu_stream)
  call gpu_free_async(rjs_d,gpu_stream)
  call gpu_free_async(mask_d,gpu_stream)
! call gpu_free_async(atom_sigma_t_d,gpu_stream)
! call gpu_free_async(atom_sigma_t_scaling_d,gpu_stream)
! call gpu_free_async(atom_sigma_r_d,gpu_stream)
  call gpu_free_async(k2_start_d,gpu_stream)
  !call gpu_free_async(n_neigh_d,gpu_stream)
  call gpu_free_async(species_d,gpu_stream)
  call gpu_free_async(species_multiplicity_d,gpu_stream)
! call gpu_free_async(rcut_hard_d,gpu_stream)
!!!!!!!!ATTENTION: as the code is now, W_d and S_d are never free'd, only in realloc case they are freed and reallocate'd issue is due to the "save" value of W and S arrays that are needed across different calls of the same function. maybe this can be fixed when reworking the memory
  !call gpu_free_async(W_d,gpu_stream)
  !call gpu_free_async(S_d,gpu_stream)
  !call gpu_free_async(k2_i_site_d,gpu_stream)
  call gpu_free_async(preflm_d,gpu_stream)
  call gpu_free_async(i_beg_d,gpu_stream)
  call gpu_free_async(i_end_d,gpu_stream)
  call gpu_free_async(multiplicity_array_d,gpu_stream)
  call gpu_free_async(skip_soap_component_d,gpu_stream)
  call gpu_free_async(sqrt_dot_p_d,gpu_stream)
! call gpu_free_async(central_weight_d,gpu_stream)
!  call gpu_free_async(global_scaling_d,gpu_stream)
  !call gpu_free_async(cnk_d,gpu_stream)
  call gpu_free(cnk_d)
  
  !call cpu_time(ttt(2))
  ttt(2)=MPI_Wtime()
  time_get_soap=time_get_soap+ttt(2)-ttt(1)
  !stop  
  end subroutine get_soap
!**************************************************************************











!**************************************************************************
  subroutine get_derivatives(radial_exp_coeff, angular_exp_coeff, radial_exp_coeff_der, &
                             angular_exp_coeff_rad_der, angular_exp_coeff_azi_der, angular_exp_coeff_pol_der, &
                             n_sites, n_max, l_max, n_neigh, rjs, cutoff, cnk_rad_der, cnk_azi_der, cnk_pol_der )

    implicit none

    real*8, intent(in) :: radial_exp_coeff(:,:), radial_exp_coeff_der(:,:), rjs(:), cutoff
    complex*16, intent(in) :: angular_exp_coeff(:,:)
    complex*16 :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:), angular_exp_coeff_pol_der(:,:)
    complex*16, intent(out) :: cnk_rad_der(:,:,:), cnk_azi_der(:,:,:), cnk_pol_der(:,:,:)
    integer, intent(in) :: n_sites, n_max, l_max, n_neigh(:)
    integer :: k2, i, j, n, m, k, l
    real*8 :: pi

    pi = dacos(-1.d0)


    k2 = 0
    do i = 1, n_sites
!     We could skip j = 1, which is the central atom
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        if( rjs(k2) < cutoff )then
          do n = 1, n_max
            do l = 0, l_max
              do m = 0, l
                k = 1 + l*(l+1)/2 + m
!               Radial derivatives:
                cnk_rad_der(k, n, k2) = 4.d0*pi * ( angular_exp_coeff(k, k2) * radial_exp_coeff_der(n, k2) + &
                                         angular_exp_coeff_rad_der(k, k2) * radial_exp_coeff(n, k2) )
!               Azimuthal angle derivatives:
                cnk_azi_der(k, n, k2) = 4.d0*pi * angular_exp_coeff_azi_der(k, k2) * radial_exp_coeff(n, k2)
!               Polar angle derivatives:
                cnk_pol_der(k, n, k2) = 4.d0*pi * angular_exp_coeff_pol_der(k, k2) * radial_exp_coeff(n, k2)
              end do
            end do
          end do
        end if
      end do
    end do

    return
  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine assign_species_multiplicity(max_species_multiplicity, species_types_soap, xyz_species, &
                                         xyz_species_supercell, n_species_soap, all_atoms, which_atom, &
                                         indices, n_neigh, n_atom_pairs, neighbors_list, mask_species, &
                                         species, species_multiplicity, species_multiplicity_supercell )
!                                         species, species_multiplicity )

    implicit none

!   Input variables
    integer, intent(in) :: max_species_multiplicity, n_species_soap, which_atom, indices(1:3), &
                           n_neigh(:), n_atom_pairs, neighbors_list(:)
    character*8, intent(in) :: species_types_soap(:), xyz_species(:), xyz_species_supercell(:)
    logical, intent(in) :: all_atoms

!   Output variables
    integer, allocatable, intent(out) :: species(:,:), species_multiplicity(:), &
                                         species_multiplicity_supercell(:)
!    integer, allocatable :: species_multiplicity_supercell(:)
    logical, allocatable, intent(out) :: mask_species(:,:)

!   Internal variables
    integer, allocatable :: species_supercell(:,:)
    integer :: n_sites, n_sites_supercell, i, j, k, i2, j2, k2, ijunk, counter


    n_sites = size(xyz_species)
    n_sites_supercell = size(xyz_species_supercell)

    allocate( species(1:max_species_multiplicity, 1:n_sites) )
    allocate( species_multiplicity(1:n_sites) )
    species = 0
    species_multiplicity = 0

    do i = 1, n_sites
      do j = 1, n_species_soap
        if( xyz_species(i) == species_types_soap(j) )then
          species_multiplicity(i) = species_multiplicity(i) + 1
          species(species_multiplicity(i), i) = j
        end if
      end do
    end do

    if( n_sites_supercell > n_sites )then
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = 0
      species_multiplicity_supercell = 0
!      counter = 0
!      do i2 = 1, indices(1)
!        do j2 = 1, indices(2)
!          do k2 = 1, indices(3)
!            do i = 1, n_sites
!              counter = counter + 1
!              species_supercell(:, counter) = species(:, i)
!            end do
!          end do
!        end do
!      end do
      do i = 1, n_sites_supercell
        do j = 1, n_species_soap
          if( xyz_species_supercell(i) == species_types_soap(j) )then
            species_multiplicity_supercell(i) = species_multiplicity_supercell(i) + 1
            species_supercell(species_multiplicity_supercell(i), i) = j
          end if
        end do
      end do
    else
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = species
      species_multiplicity_supercell = species_multiplicity
    end if

!   This is perhaps not the most efficient way to select only one atom, fix in the future <----- FIX THIS
    if( .not. all_atoms )then
      n_sites = 1
      deallocate( species )
      allocate( species(1:max_species_multiplicity, 1:n_sites) )
      species(1:max_species_multiplicity, 1) = species_supercell(1:max_species_multiplicity, which_atom)
      species_supercell(1:max_species_multiplicity, which_atom) = species_supercell(1:max_species_multiplicity, 1)
      species_supercell(1:max_species_multiplicity, 1) = species(1:max_species_multiplicity, 1)

      ijunk = species_multiplicity(which_atom)
      deallocate( species_multiplicity )
      allocate( species_multiplicity(1:n_sites) )
      species_multiplicity(1) = ijunk
    end if

    allocate( mask_species(1:n_atom_pairs, 1:n_species_soap) )
    mask_species = .false.

    k2 = 0
    do i = 1, n_sites
      do k = 1, n_neigh(i)
        k2 = k2 + 1
        j = neighbors_list(k2)
!        if( k == 1 )then
!          do i2 = 1, species_multiplicity_supercell(i)
!            mask_species(k2, species_supercell(i2, j)) = .true.
!          end do
!        else
!          do i2 = 1, species_multiplicity_supercell(i)
!            mask_species(k2, species_supercell(i2, j)) = .true.
!          end do
!        end if
          do i2 = 1, species_multiplicity_supercell(i)
            if( species_supercell(i2, j) > 0 )then
              mask_species(k2, species_supercell(i2, j)) = .true.
            end if
          end do
      end do
    end do

!    deallocate( species_supercell, species_multiplicity_supercell )
    deallocate( species_supercell )


  end subroutine
!**************************************************************************


end module soap_turbo_desc
