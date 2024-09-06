! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_interface.f90, is copyright (c) 2019-2023,
! HND X   Miguel A. Caro and Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
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

module exp_interface
  use types
  use read_files
#ifdef _MPIF90
  use mpi
#endif
  use exp_utils
  use soap_turbo_functions
  use F_B_C

  use gpu_var_mod  
  use gpu_var_int_mod
  use gpu_var_double_mod  
  
  use iso_c_binding
contains

  ! This module implements the interfaces for the gradient of experimental functions

  ! subroutine galloc_int_1d(gpu_var, array)
  !   implicit none
  !   type(gpu_variable) :: gpu_var
  !   integer, allocatable :: array(:)
  !   integer :: n
  !   integer(c_size_t) :: size 
    
  !   if( gpu_var%allocated == .false. )then
  !      n = size(array, 1)
  !      size = n * c_int
  !      gpu_var%st 

       
  !   else
  !      print *, "Trying to reallocate gpu_variable!"
  !   end if
    
  ! end subroutine galloc
  
  subroutine get_time(time)
    implicit none
    real*8 :: time

#ifdef _MPIF90
    time = MPI_Wtime()
#else
    call cpu_time(_time)
#endif 
  end subroutine get_time
  

  
  subroutine get_write_condition( do_mc, do_md, mc_istep, md_istep, write_xyz, write_condition)
    implicit none
    logical, intent(in) :: do_mc, do_md
    integer, intent(in) :: mc_istep, md_istep, write_xyz
    logical, intent(out) :: write_condition

    if (do_mc)then
       write_condition = mc_istep > -1 .and. (modulo( mc_istep, write_xyz) ==  0)
    elseif (do_md)then
       write_condition = md_istep > -1 .and. (modulo( md_istep, write_xyz) ==  0)
    else
       write_condition = .true.
    end if
  end subroutine get_write_condition

  subroutine get_overwrite_condition( do_mc, do_md, mc_istep, md_istep, write_xyz, write_condition)
    implicit none
    logical, intent(in) :: do_mc, do_md
    integer, intent(in) :: mc_istep, md_istep, write_xyz
    logical, intent(out) :: write_condition

    if (do_mc)then
       write_condition = mc_istep ==  0
    elseif (do_md)then
       write_condition =  md_istep ==  0
    else
       write_condition = .true.
    end if
  end subroutine get_overwrite_condition


  subroutine write_partial_exp(do_mc, do_md, mc_istep, md_istep,&
       & write_xyz, has_partials, n_species, n_samples, n_dim_partial&
       &, x, y, partials, species_types, name )
    implicit none
    logical, intent(in) :: do_mc, do_md, has_partials
    integer, intent(in) :: mc_istep, md_istep, write_xyz, n_species, n_samples, n_dim_partial
    real*8, intent(in) :: x(:), y(:), partials(:,:)
    character*8, allocatable :: species_types(:)
    character(len = *), intent(in) :: name
    integer :: j, k, n_dim_idx
    character*1024 :: filename
    logical :: write_condition, overwrite_condition

    call get_overwrite_condition( do_mc, do_md&
         &, mc_istep, md_istep, write_xyz,&
         & overwrite_condition)

    if (has_partials)then
       n_dim_idx = 1
       outer: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle

             write(filename,'(A)')&
                  & name // '_' // trim(species_types(j)) // '_' // trim(species_types(k)) //&
                  & "_prediction.dat"
             call write_exp_datan(x(1:n_samples),&
                  & partials(1:n_samples, n_dim_idx),&
                  & overwrite_condition, filename, name)

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial )then
                exit outer
             end if

          end do
       end do outer
    end if

    write(filename,'(A)')&
         & name // "_total.dat"
    call write_exp_datan(x(1:n_samples),&
         &y(1:n_samples),&
         & overwrite_condition, filename, name)
  end subroutine write_partial_exp




  subroutine preprocess_exp_data(params, x, y, label, n_sites, V, input, output, exp)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, intent(in), allocatable :: x(:)
    real*8, intent(in) :: V
    real*8, intent(inout), allocatable :: y(:)
    integer, intent(in) :: n_sites
    character*1024, intent(in) :: label
    real*8, parameter :: pi = acos(-1.0)
    real*8 :: mag, dx, rho
    logical, intent(in) :: exp
    character*32, intent(inout) :: output
    character*1024, intent(inout) :: input

    dx = x(2) - x(1)
    rho = dfloat(n_sites) / V

    if     ( trim(label) == "xps" )then
       ! calculate the magnitude and normalize
       mag = sqrt(dot_product(y, y) )
       y = y / mag
       output = "xps"

    elseif ( trim(label) == "pair_distribution" )then
       output = params%pair_distribution_output
       if ( trim(params%pair_distribution_output) == "D(r)" .and. .not. ( exp .and. trim(input) == "D(r)" ) )then
          ! D(r) = 4pi rho * r * G(r)
          ! G(r) = total pair distribution function
          y = 4.d0 * pi * rho * x * (y - 1.d0)

       end if
    elseif ( trim(label) == "xrd" )then
       output = params%xrd_output
       ! mag = sqrt(dot_product(y, y)) * dx
       ! y = y / mag

       if ( trim(params%xrd_output) == "q*i(q)" .and. params%q_units &
            &== "q" )then
          if ( exp .and. ( trim(input) == "i(q)" .or. trim(input) == "F(q)") )then
             y = x * (y - 1.d0)
          end if
       end if
    elseif ( trim(label) == "nd" )then
       output = params%nd_output
       ! mag = sqrt(dot_product(y, y)) * dx
       ! y = y / mag

       if ( trim(params%nd_output) == "q*i(q)" .and. params%q_units &
            &== "q" )then
          if ( exp .and. ( trim(input) == "i(q)" .or. trim(input) == "F(q)") )then
             y = x * (y - 1.d0)
          end if
       end if

    end if
  end subroutine preprocess_exp_data

  subroutine gpu_copy_pdf( n_samples, pdf_d, pdf,  gpu_stream )
    implicit none
    integer :: n_samples
    real*8, intent(inout), target :: pdf(1:n_samples)
    type(c_ptr) :: pdf_d, gpu_stream
    integer(c_size_t) :: size

    size = n_samples * c_double

    call cpy_dtoh( pdf_d, c_loc(pdf), size, gpu_stream)
    
  end subroutine gpu_copy_pdf

  ! subroutine gpu_copy_pdf_der( n_samples, pdf_d, pdf,  gpu_stream )
  !   implicit none
  !   integer :: n_samples
  !   real*8, intent(inout), target :: pdf(1:n_samples)
  !   type(c_ptr) :: pdf_d, gpu_stream
  !   integer(c_size_t) :: size

  !   size = n_samples * c_double

  !   call cpy_dtoh( pdf_d, c_loc(pdf), size, gpu_stream)
    
  ! end subroutine gpu_copy_pdf_der

  subroutine estimate_device_memory_usage( n_sites, n_pairs, nk, n_samples, n_samples_sf, total, standard )
    implicit none
    integer :: n_sites, nk, n_samples, n_samples_sf, n_pairs
    real*8, intent(inout) :: total
    real*8 :: total_exp=0.d0, total_standard=0.d0
    real*8 :: nk_int, nk_float, Gk, dermat, sf, fi, pref, xyz, forces, neigh_list, to_gb
    logical :: standard
    
    to_gb = 1 / dfloat(1024**3)

    neigh_list = dfloat(n_pairs) * 4.d0 * to_gb    
    
    forces   = dfloat(n_sites) * 8.d0 * to_gb    
    nk_int   = dfloat(nk) * 4.d0 * to_gb
    nk_float = dfloat(nk) * 8.d0 * to_gb
    sf       = dfloat(n_samples * n_samples_sf) * 8.d0 * to_gb

    pref = dfloat( n_samples_sf ) * 8.d0 * to_gb
    Gk = n_samples * nk_float
    dermat = n_samples_sf * nk_float
    fi = 3 * nk_float
    xyz  = 3 * nk_float 


    if (standard)then 
       write(*,'(A)') "\n > Estimating memory for normal allocations \n"
       total_standard = total_standard + neigh_list
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device neigh_list_d  ", neigh_list , " Gb\n"

       total_standard = total_standard + neigh_list
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device neigh_spec_d  ", neigh_list , " Gb\n"

       total_standard = total_standard + neigh_list*2.d0 * 3.d0
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device xyz_d  ", neigh_list*2.d0 * 3.d0 , " Gb\n"

       total_standard = total_standard + neigh_list*2.d0 
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device rjs_d  ", neigh_list*2.d0  , " Gb\n"

       total_standard = total_standard + neigh_list*2.d0 
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device rjs_d  ", neigh_list*2.d0  , " Gb\n"

       total_standard = total_standard + forces/3.d0
       write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device species_d  ", forces/3.d0 , " Gb\n"    
    end if
    
    

    
    
    write(*,'(A)') "\n > Estimating memory xrd and pdf calculation\n"
    total_exp = total_exp + nk_int
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device k_index_d  ", nk_int , " Gb\n"

    total_exp = total_exp + 3.d0 * Gk
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device      Gk_d  ", Gk , " Gb\n"
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device     Gka_d  ", Gk , " Gb\n"

    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device par_pdf_d  ", Gk , " Gb\n"
    
    total_exp = total_exp + dermat
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device  dermat_d  ", dermat , " Gb\n"    

    
    total_exp = total_exp + fi
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device      fi_d  ", fi , " Gb\n"    

    total_exp = total_exp + xyz
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device     xyz_d  ", xyz , " Gb\n"    


    total_exp = total_exp + forces
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device   forces_d  ", forces , " Gb\n"    
    

    
    total_exp = total_exp + pref
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device    pref_d  ", pref , " Gb\n"    

    total_exp = total_exp + pref
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device  scat_f_d  ", pref , " Gb\n"    
    
    
    total_exp = total_exp + sf
    write(*, '(A,1X,F7.4,1X,A)') "Gb/core: device  sinc_f_d  ", sf , " Gb\n"


    total = total + total_exp + total_standard
    write(*, '(A,1X,F7.4,1X,A)') "\nTotal device memory usage in block:", total_standard+total_exp , " Gb\n"
    write(*, '(A,1X,F7.4,1X,A)') "\nTotal device memory usage:", total , " Gb\n"
    
  end subroutine estimate_device_memory_usage


  subroutine get_n_atoms_of_species(n_atoms_of_species, n_sites, species, n_species)
    implicit none
    integer, intent(in) :: n_sites, n_species
    integer, allocatable, intent(in) :: species(:)
    real*8, allocatable, intent(out) :: n_atoms_of_species(:)
    integer :: j, i2
    
    if(allocated(n_atoms_of_species)) deallocate(n_atoms_of_species)
    allocate(n_atoms_of_species(1:n_species))

    do j = 1, n_species
       n_atoms_of_species(j) = 0.d0
       do i2 = 1, n_sites
          if ( species(i2) == j)then
             n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
          end if
       end do
    end do
  end subroutine get_n_atoms_of_species


  subroutine setup_batched_pair_distribution( r_min, r_max, r_cut, n_samples, x, dV )
    implicit none
    real*8, intent(in) :: r_min, r_max, r_cut
    integer, intent(in) :: n_samples
    real*8, allocatable, intent(out) :: x(:), dV(:)

    
    if( allocated( x )) deallocate( x )
    if( allocated( dV ) ) deallocate( dV )

    allocate( x( 1:n_samples ) ) 
    allocate( dV( 1:n_samples ) )   

    call setup_pdf_arrays( r_min, r_max, r_cut, n_samples, x, dV)

  end subroutine setup_batched_pair_distribution


  subroutine gpu_malloc_neighbors(gpu_neigh, n_sites, n_pairs,&
       & n_neigh, species, neighbor_species, neighbors_list,  rjs, &
       & xyz, gpu_stream, rank)
    implicit none
    type( gpu_neigh_storage_type ) :: gpu_neigh
    integer, intent(in):: n_sites, n_pairs, rank
    integer, intent(in), target :: n_neigh(:), species(:), neighbor_species(:), neighbors_list(:)
    real*8, intent(in), target :: rjs(:), xyz(:,:)

    integer, allocatable, target :: n_neigh_temp(:), species_temp(:), neighbor_species_temp(:), neighbors_list_temp(:)
    real*8,  allocatable, target :: rjs_temp(:), xyz_temp(:,:)

    
    integer( c_size_t ) :: st_n_sites_int, st_n_atom_pairs_int,  st_n_atom_pairs_double
    type(c_ptr) :: gpu_stream
    integer :: i

    st_n_sites_int = n_sites * c_int
    gpu_neigh % st_n_neigh_d = st_n_sites_int
    gpu_neigh % st_species_d = st_n_sites_int    
    call gpu_malloc_all(gpu_neigh % n_neigh_d,st_n_sites_int,gpu_stream)
    call cpy_htod(c_loc(n_neigh), gpu_neigh % n_neigh_d, st_n_sites_int,gpu_stream)
    call gpu_malloc_all(gpu_neigh % species_d,st_n_sites_int,gpu_stream)
    call cpy_htod(c_loc(species), gpu_neigh % species_d, st_n_sites_int,gpu_stream)

    st_n_atom_pairs_int = n_pairs * c_int
    gpu_neigh % st_neighbor_species_d = st_n_atom_pairs_int
    gpu_neigh % st_neighbors_list_d    = st_n_atom_pairs_int
    call gpu_malloc_all(gpu_neigh % neighbor_species_d,st_n_atom_pairs_int,gpu_stream)
    call cpy_htod(c_loc(neighbor_species),gpu_neigh % neighbor_species_d, st_n_atom_pairs_int,gpu_stream)
    call gpu_malloc_all(gpu_neigh % neighbors_list_d,st_n_atom_pairs_int,gpu_stream)
    call cpy_htod(c_loc(neighbors_list),gpu_neigh % neighbors_list_d, st_n_atom_pairs_int,gpu_stream)
           
    st_n_atom_pairs_double = n_pairs * c_double
    gpu_neigh % st_rjs_d    =     st_n_atom_pairs_double
    gpu_neigh % st_xyz_d    = 3 * st_n_atom_pairs_double
    
    call gpu_malloc_all(gpu_neigh % rjs_d,st_n_atom_pairs_double,gpu_stream)
    call cpy_htod(c_loc(rjs),gpu_neigh % rjs_d, st_n_atom_pairs_double,gpu_stream)
    call gpu_malloc_all(gpu_neigh % xyz_d,3*st_n_atom_pairs_double,gpu_stream)
    call cpy_htod(c_loc(xyz),gpu_neigh % xyz_d,3*st_n_atom_pairs_double,gpu_stream)


    print *, "-- Rank ", rank, " ", " malloc neighbors: n_sites_temp = ", n_sites, " n_pairs_temp = ", n_pairs 

    
    ! allocate( n_neigh_temp( 1:n_sites ) )
    ! n_neigh_temp = n_neigh
    ! allocate( species_temp( 1:n_sites ) )
    ! species_temp = species 
    
    ! allocate( neighbor_species_temp( 1:n_pairs ) )
    ! neighbor_species_temp = neighbor_species
    ! allocate( neighbors_list_temp( 1:n_pairs ) )
    ! neighbors_list_temp = neighbors_list
    ! allocate( rjs_temp( 1:n_pairs ) )
    ! rjs_temp = rjs
    ! allocate( xyz_temp(1:3, 1:n_pairs ) )        
    ! xyz_temp = xyz

    
    ! st_n_sites_int = n_sites * c_int
    ! gpu_neigh % st_n_neigh_d = st_n_sites_int
    ! gpu_neigh % st_species_d = st_n_sites_int    
    ! call gpu_malloc_all(gpu_neigh % n_neigh_d,st_n_sites_int,gpu_stream)
    ! call cpy_htod(c_loc(n_neigh_temp), gpu_neigh % n_neigh_d, st_n_sites_int,gpu_stream)
    ! call gpu_malloc_all(gpu_neigh % species_d,st_n_sites_int,gpu_stream)
    ! call cpy_htod(c_loc(species_temp), gpu_neigh % species_d, st_n_sites_int,gpu_stream)

    ! st_n_atom_pairs_int = n_pairs * c_int
    ! gpu_neigh % st_neighbor_species_d = st_n_atom_pairs_int
    ! gpu_neigh % st_neighbors_list_d    = st_n_atom_pairs_int
    ! call gpu_malloc_all(gpu_neigh % neighbor_species_d,st_n_atom_pairs_int,gpu_stream)
    ! call cpy_htod(c_loc(neighbor_species_temp),gpu_neigh % neighbor_species_d, st_n_atom_pairs_int,gpu_stream)
    ! call gpu_malloc_all(gpu_neigh % neighbors_list_d,st_n_atom_pairs_int,gpu_stream)
    ! call cpy_htod(c_loc(neighbors_list_temp),gpu_neigh % neighbors_list_d, st_n_atom_pairs_int,gpu_stream)
           
    ! st_n_atom_pairs_double = n_pairs * c_double
    ! gpu_neigh % st_rjs_d    =     st_n_atom_pairs_double
    ! gpu_neigh % st_xyz_d    = 3 * st_n_atom_pairs_double
    
    ! call gpu_malloc_all(gpu_neigh % rjs_d,st_n_atom_pairs_double,gpu_stream)
    ! call cpy_htod(c_loc(rjs_temp),gpu_neigh % rjs_d, st_n_atom_pairs_double,gpu_stream)
    ! call gpu_malloc_all(gpu_neigh % xyz_d,3*st_n_atom_pairs_double,gpu_stream)
    ! call cpy_htod(c_loc(xyz_temp),gpu_neigh % xyz_d,3*st_n_atom_pairs_double,gpu_stream)

    ! call gpu_stream_sync( gpu_stream )
    
    
!    deallocate( n_neigh_temp, species_temp, neighbor_species_temp, neighbors_list_temp, rjs_temp, xyz_temp)

    
  end subroutine gpu_malloc_neighbors

  
  subroutine gpu_free_neighbors(gpu_neigh, gpu_stream)
    implicit none
    type( gpu_neigh_storage_type ) :: gpu_neigh
    type(c_ptr) :: gpu_stream
    
    call gpu_free_async(gpu_neigh % n_neigh_d, gpu_stream)
    call gpu_free_async(gpu_neigh % species_d, gpu_stream)
    call gpu_free_async(gpu_neigh % neighbor_species_d, gpu_stream)
    call gpu_free_async(gpu_neigh % neighbors_list_d, gpu_stream)
    call gpu_free_async(gpu_neigh % rjs_d, gpu_stream)
    call gpu_free_async(gpu_neigh % xyz_d, gpu_stream)
    call gpu_stream_sync( gpu_stream )
  end subroutine gpu_free_neighbors


    
    
  

  subroutine collect_batched_pair_distribution( n_batches, gpu_host, n_dim_partial, n_samples, &
       pair_distribution_partial, n_species, n_atoms_of_species, v_uc)
    implicit none
    integer, intent(in) :: n_dim_partial, n_samples, n_batches, n_species 
    type( gpu_host_batch_storage_type ), intent(in), allocatable :: gpu_host(:)
    real*8, allocatable, intent( out ) :: pair_distribution_partial(:,:)
    real*8, allocatable, intent(in) :: n_atoms_of_species(:)
    real*8, intent(in) :: v_uc
    real*8, allocatable :: pair_distribution_partial_temp(:,:), factors(:)
    real*8 :: f 
    integer :: i, j, k , n_dim_idx, ierr


    allocate( factors( 1:n_dim_partial ) )

    n_dim_idx = 1
    outer: do j = 1, n_species
       do k = 1, n_species

          if( j > k ) cycle

          if( j == k ) f = v_uc  /  n_atoms_of_species(j) /  n_atoms_of_species(k)
          if( j /= k ) f = v_uc  /  n_atoms_of_species(j) /  n_atoms_of_species(k) / 2.d0

          factors(n_dim_idx) = f
          
          n_dim_idx = n_dim_idx + 1
          if( n_dim_idx > n_dim_partial ) exit outer
       end do
    end do outer
    
    
    allocate( pair_distribution_partial( 1:n_samples, 1:n_dim_partial ) )    
    pair_distribution_partial = 0.d0    


    do j = 1, n_dim_partial
       f = factors(j)
       do i = 1, n_batches
          
          pair_distribution_partial(1:n_samples, j) = pair_distribution_partial(1:n_samples, j) &
               + gpu_host( i ) % host( j ) % pair_distribution_partial_h( 1:n_samples ) * f
          
       end do
    end do
    
    deallocate(factors)
    
#ifdef _MPIF90
    allocate( pair_distribution_partial_temp( 1:n_samples, 1:n_dim_partial ) )
    pair_distribution_partial_temp = 0.d0
    
    call mpi_reduce(pair_distribution_partial,&
         & pair_distribution_partial_temp, n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
         & 0, MPI_COMM_WORLD, ierr)
    
    pair_distribution_partial =  pair_distribution_partial_temp
    deallocate( pair_distribution_partial_temp )


    call mpi_bcast(pair_distribution_partial, n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
         & MPI_COMM_WORLD, ierr)
    
!    deallocate( pair_distribution_partial_temp )
#endif     

    
  end subroutine collect_batched_pair_distribution


  subroutine collect_batched_forces( n_batches, gpu_host, n_dim_partial, &
       forces_out, virial_out,  n_sites)
    implicit none
    integer, intent(in) :: n_batches, n_sites, n_dim_partial
    type( gpu_host_batch_storage_type ), intent(in), allocatable :: gpu_host(:)
    real*8, allocatable :: forces_out(:,:)
    real*8 :: virial_out(1:3,1:3)
    integer :: i, j, k , n_dim_idx

    ! allocate( forces_out( 1:3, 1:n_sites))
    ! forces_out = 0.d0    
    ! print *, " size forces_out  ", size( forces_out,1 ), " ", size( forces_out,2 )          

    do i = 1, n_batches
       do j = 1, n_dim_partial

!          print *, " size forces_h   i ", i, " j ", j , " ", size( gpu_host( i ) % host( j ) % forces_h,1 ), " ", size( gpu_host( i ) % host( j ) % forces_h,2 )

          forces_out = forces_out + gpu_host( i ) % host( j ) % forces_h
          virial_out = virial_out + gpu_host( i ) % host( j ) % virial_h          

!          print *, "virial h ", i, j, gpu_host( i ) % host( j ) % virial_h 
          
       end do
    end do
    
  end subroutine collect_batched_forces
  
  subroutine free_host_batches(gpu_host, n_batches, n_dim_partial)
    implicit none
    type( gpu_host_batch_storage_type ), allocatable :: gpu_host(:)
    integer, intent(in) :: n_batches, n_dim_partial
    integer :: i, j

    do i = 1, n_batches
       do j = 1, n_dim_partial

          if( allocated( gpu_host( i ) % host( j ) % xyz_k_h ) )&
               deallocate(gpu_host( i ) % host( j ) % xyz_k_h)

          if( allocated( gpu_host( i ) % host( j ) % pair_distribution_partial_h ) )&
               deallocate(gpu_host( i ) % host( j ) % pair_distribution_partial_h)

          if( allocated( gpu_host( i ) % host( j ) % pair_distribution_partial_der_h ) )&
               deallocate(gpu_host( i ) % host( j ) % pair_distribution_partial_der_h)

          if( allocated( gpu_host( i ) % host( j ) % forces_h ) )&
               deallocate(gpu_host( i ) % host( j ) % forces_h)

          if( allocated( gpu_host( i ) % host( j ) % rjs_index_h ) )&
               deallocate(gpu_host( i ) % host( j ) % rjs_index_h)

          if( allocated( gpu_host( i ) % host( j ) % k_index_h ) )&
               deallocate(gpu_host( i ) % host( j ) % k_index_h)
          
       end do
       deallocate( gpu_host( i ) % host )
    end do
    deallocate( gpu_host )
    
    
  end subroutine free_host_batches


    subroutine free_exp_batches(gpu_exp,  n_batches)
    implicit none
    type( gpu_storage_type ), allocatable :: gpu_exp(:)
    integer, intent(in) :: n_batches
    integer :: i, j

    do i = 1, n_batches

       if( allocated( gpu_exp(i) % nk )) deallocate( gpu_exp(i) % nk )

       if( allocated( gpu_exp(i) % nk_d )) deallocate( gpu_exp(i) % nk_d )

       if( allocated( gpu_exp(i) % k_index_d )) deallocate( gpu_exp(i) % k_index_d )

       if( allocated( gpu_exp(i) % j2_index_d )) deallocate( gpu_exp(i) % j2_index_d )

       if( allocated( gpu_exp(i) % xyz_k_d )) deallocate( gpu_exp(i) % xyz_k_d )

       if( allocated( gpu_exp(i) % pair_distribution_partial_d )) deallocate( gpu_exp(i) % pair_distribution_partial_d )

       if( allocated( gpu_exp(i) % pair_distribution_partial_der_d )) deallocate( gpu_exp(i) % pair_distribution_partial_der_d )

       if( allocated( gpu_exp(i) % st_nk_d )) deallocate( gpu_exp(i) % st_nk_d )

       if( allocated( gpu_exp(i) % st_k_index_d )) deallocate( gpu_exp(i) % st_k_index_d )

       if( allocated( gpu_exp(i) % st_j2_index_d )) deallocate( gpu_exp(i) % st_j2_index_d )

       if( allocated( gpu_exp(i) % rjs_index_d )) deallocate( gpu_exp(i) % rjs_index_d )       

       if( allocated( gpu_exp(i) % st_pair_distribution_partial_d )) deallocate( gpu_exp(i) % st_pair_distribution_partial_d )

       if( allocated( gpu_exp(i) % st_pair_distribution_partial_der_d )) deallocate( gpu_exp(i) % st_pair_distribution_partial_der_d )

       if( allocated( gpu_exp(i) % nk_flags_d)) deallocate( gpu_exp(i) % nk_flags_d )
       if( allocated( gpu_exp(i) % nk_flags_sum_d)) deallocate( gpu_exp(i) % nk_flags_sum_d )                     
       


       

    end do
    deallocate( gpu_exp )
    
  end subroutine free_exp_batches

  
  
  subroutine calculate_batched_pair_distribution(gpu_exp,  gpu_host,&
       & gpu_neigh, x, dV, n_atoms_of_species, n_species, n_sites,&
       & i_beg, i_end, j_beg, j_end, n_samples,  r_min, r_max, r_cut, kde_sigma, &
       & gpu_stream, x_d, dV_d, v_uc, rank)
    implicit none
    ! Input variables 
    integer, intent(in) :: n_species, i_beg, i_end, j_beg, j_end, n_sites, n_samples, rank
    real*8, allocatable, intent(in) :: n_atoms_of_species(:), x(:), dV(:)
    real*8, intent(in) :: r_min, r_max, r_cut, kde_sigma, v_uc
    type( gpu_storage_type ), intent( inout ) :: gpu_exp 
    type( gpu_host_batch_storage_type ), intent(inout),  target :: gpu_host
    type( gpu_neigh_storage_type), intent( in ) :: gpu_neigh

    
    
    ! Local variables
    real*8, parameter :: pi = acos(-1.0)    
    integer :: n_dim_partial, n_dim_idx
    type( c_ptr ), intent(in) :: x_d, dV_d
    integer( c_size_t ) :: st_x_d
    integer, target :: nk_temp(1)    
    type( c_ptr ) :: nk_flags_d, nk_flags_sum_d
    integer( c_size_t ) :: st_nk_flags, st_nk_temp
    type( c_ptr ) :: rjs_index_d, pdf_to_reduce_d
    integer( c_size_t ) :: st_rjs_index_d, st_k_index_d, st_pdf_to_reduce_d

    real*8 :: pdf_factor, der_factor = 0.d0, f 
    integer :: i, j, k, l 

    type(c_ptr) :: gpu_stream

    real*8, allocatable, target :: x_check(:)

    
    n_dim_partial = n_species * ( n_species + 1 ) / 2 

    
    allocate( gpu_host % host( 1:n_dim_partial ) )
    
    
    ! allocate( nk_flags_d(1:n_dim_partial) )
    ! allocate( nk_flags_sum_d(1:n_dim_partial) )        
    
    allocate( gpu_exp % nk(1:n_dim_partial) )    
    allocate( gpu_exp % nk_d(1:n_dim_partial) )
    allocate( gpu_exp % k_index_d(1:n_dim_partial) )
    allocate( gpu_exp % j2_index_d(1:n_dim_partial) )
!   allocate( gpu_exp % rjs_index_d(1:n_dim_partial) )    
    allocate( gpu_exp % xyz_k_d(1:n_dim_partial) )    
    allocate( gpu_exp % pair_distribution_partial_d(1:n_dim_partial) )
    allocate( gpu_exp % nk_flags_sum_d(1:n_dim_partial) )
    allocate( gpu_exp % nk_flags_d(1:n_dim_partial) )
    allocate( gpu_exp % rjs_index_d(1:n_dim_partial) )        
    
!   allocate( gpu_exp % pair_distribution_partial_der_d(1:n_dim_partial) )
    allocate( gpu_exp % st_nk_d(1:n_dim_partial) )
    allocate( gpu_exp % st_k_index_d(1:n_dim_partial) )
    allocate( gpu_exp % st_j2_index_d(1:n_dim_partial) )
!   allocate( gpu_exp % st_rjs(1:n_dim_partial) )
    allocate( gpu_exp % st_pair_distribution_partial_d(1:n_dim_partial) )
!   allocate( gpu_exp % st_pair_distribution_partial_der_d(1:n_dim_partial) )

    
    ! st_x_d = n_samples * c_double

    ! call gpu_malloc_all(x_d,      st_x_d, gpu_stream)             
    ! call cpy_htod( c_loc( x ), x_d, st_x_d, gpu_stream )

    ! call gpu_malloc_all(dV_d,      st_x_d, gpu_stream)             
    ! call cpy_htod( c_loc( dV ), dV_d, st_x_d, gpu_stream )

    
    n_dim_idx = 1
    !    call gpu_stream_sync( gpu_stream )
    outer1: do j = 1, n_species
       do k = 1, n_species

          if (j > k) cycle ! We have already calculated the pair correlation function!

          st_nk_temp = 1*c_int
          call gpu_malloc_all(gpu_exp % nk_d(n_dim_idx), st_nk_temp, gpu_stream)                          
          st_nk_flags = (j_end - j_beg + 1)  * c_int
          call gpu_malloc_all(gpu_exp % nk_flags_d(n_dim_idx),      st_nk_flags, gpu_stream)             
          call gpu_memset_async(gpu_exp % nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
          call gpu_malloc_all(gpu_exp % nk_flags_sum_d(n_dim_idx),      st_nk_flags, gpu_stream)                          
!          call gpu_stream_sync(gpu_stream)

          ! print *, "-- Rank ", rank, " ",   "int(gpu_exp % nk_d(n_dim_idx))"
          ! call gpu_print_pointer_int(gpu_exp % nk_d(n_dim_idx))
          ! print *, "-- Rank ", rank, " ",   "int(gpu_exp % nk_flags_d(n_dim_idx))"
          ! call gpu_print_pointer_int(gpu_exp % nk_flags_d(n_dim_idx))
          ! print *, "-- Rank ", rank, " ",   "int(gpu_exp % nk_flags_sum_d(n_dim_idx))"
          ! call gpu_print_pointer_int(gpu_exp % nk_flags_sum_d(n_dim_idx))
          ! print *, "-- Rank ", rank, " ",   "int(gpu_neigh % neighbors_list_d  )"
          ! call gpu_print_pointer_int(gpu_neigh % neighbors_list_d  )
          ! print *, "-- Rank ", rank, " ",   "int(gpu_neigh % n_neigh_d         )"
          ! call gpu_print_pointer_int(gpu_neigh % n_neigh_d         )
          ! print *, "-- Rank ", rank, " ",   "int(gpu_neigh % neighbor_species_d)"
          ! call gpu_print_pointer_int(gpu_neigh % neighbor_species_d)
          ! print *, "-- Rank ", rank, " ",   "int(gpu_neigh % species_d         )"
          ! call gpu_print_pointer_int(gpu_neigh % species_d         )
          ! print *, "-- Rank ", rank, " ",   "double(gpu_neigh % rjs_d             )"
          ! call gpu_print_pointer_double(gpu_neigh % rjs_d             )
          ! print *, "-- Rank ", rank, " ",   "double(gpu_neigh % xyz_d             )"
          ! call gpu_print_pointer_double(gpu_neigh % xyz_d             )
          
          
          call gpu_get_pair_distribution_nk(1, i_end - i_beg + 1, j_end - j_beg + 1, n_sites, & ! i_beg, i_end, j_end, n_sites, &
               gpu_neigh % neighbors_list_d,&
               gpu_neigh % n_neigh_d, &
               gpu_neigh % neighbor_species_d,&
               gpu_neigh % species_d,&
               gpu_neigh % rjs_d, &
               gpu_neigh % xyz_d, &             
               r_min, r_max, r_cut, 6.d0*kde_sigma,&
               gpu_exp % nk_d(n_dim_idx), &
               gpu_exp % nk_flags_d(n_dim_idx), &
               gpu_exp % nk_flags_sum_d(n_dim_idx), &
               j, k, gpu_stream)

          !             call gpu_device_sync()             

          call gpu_free_async(gpu_exp % nk_flags_d(n_dim_idx), gpu_stream)

          ! Now copy the value of nk from the gpu
!          print *, "out of pdf nk kernel "
          st_nk_temp = 1*c_int
          call cpy_dtoh(gpu_exp % nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)
          !                call gpu_stream_sync(gpu_stream)

          gpu_exp % nk(n_dim_idx) = nk_temp(1)
          print *, " -- Rank ", rank, " nk batched = ", gpu_exp % nk(n_dim_idx)

          ! Now we create temporary arrays for the k indices

          st_rjs_index_d = gpu_exp % nk(n_dim_idx) * c_double
          call gpu_malloc_all(gpu_exp % rjs_index_d(n_dim_idx), st_rjs_index_d, gpu_stream)             
          call gpu_memset_async(gpu_exp % rjs_index_d(n_dim_idx), 0, st_rjs_index_d, gpu_stream)             


          !          call gpu_set_pair_distribution_rjs_only(j_end, gpu_exp % rjs_d, rjs_index_d, nk_flags_sum_d, gpu_stream )

          
          gpu_exp % st_k_index_d(n_dim_idx) = gpu_exp % nk(n_dim_idx) * c_int
          call gpu_malloc_all(gpu_exp % k_index_d(n_dim_idx),      gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             
          call gpu_memset_async(gpu_exp % k_index_d(n_dim_idx), 0, gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             

          call gpu_malloc_all(gpu_exp % j2_index_d(n_dim_idx),      gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             
          call gpu_memset_async(gpu_exp % j2_index_d(n_dim_idx), 0, gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             

          ! st_rjs_index_d = gpu_exp % nk(n_dim_idx) * c_double
          ! call gpu_malloc_all(rjs_index_d, st_rjs_index_d, gpu_stream)             
          ! call gpu_memset_async(rjs_index_d, 0, st_rjs_index_d, gpu_stream)             


          ! call gpu_set_pair_distribution_rjs_only(j_end, gpu_exp % rjs_d, rjs_index_d, nk_flags_sum_d, gpu_stream )

          call gpu_malloc_all(gpu_exp % xyz_k_d(n_dim_idx), 3*st_rjs_index_d, gpu_stream)             
          call gpu_memset_async(gpu_exp % xyz_k_d(n_dim_idx), 0, 3*st_rjs_index_d, gpu_stream)             
!          call gpu_stream_sync(gpu_stream)

          !          call gpu_meminfo()

          ! print *, " "
          ! print *, "-- Rank ", rank, " set k ",   "double(gpu_neigh % rjs_d             )"
          ! call gpu_print_pointer_double(gpu_neigh % rjs_d             )
          ! print *, "-- Rank ", rank, " set k ",   "double(gpu_neigh % xyz_d             )"
          ! call gpu_print_pointer_double(gpu_neigh % xyz_d             )
          ! print *, "-- Rank ", rank, " set k ",   "int(gpu_exp % nk_flags_d(n_dim_idx))"
          ! call gpu_print_pointer_int(gpu_exp % nk_flags_d(n_dim_idx))
          ! print *, "-- Rank ", rank, " set k ",   "int(gpu_exp % nk_flags_sum_d(n_dim_idx))"
          ! call gpu_print_pointer_int(gpu_exp % nk_flags_sum_d(n_dim_idx))

          
          call gpu_set_pair_distribution_k_index(1, i_end - i_beg + 1, j_end - j_beg + 1, n_sites,& ! i_beg, i_end, j_end, n_sites,&
               gpu_neigh % neighbors_list_d,&
               gpu_neigh % rjs_d, &
               gpu_neigh % xyz_d, &
               gpu_exp % k_index_d(n_dim_idx),&
               gpu_exp % j2_index_d(n_dim_idx),&
               gpu_exp % rjs_index_d(n_dim_idx), &
               gpu_exp % xyz_k_d(n_dim_idx), &
               gpu_exp % nk_flags_d(n_dim_idx), gpu_exp % nk_flags_sum_d(n_dim_idx),&
               gpu_stream)
          !call gpu_device_sync()             
!          call gpu_stream_sync(gpu_stream)
          
!          print *, " >> Set batch arrays on gpu "
          call gpu_free_async(gpu_exp % nk_flags_sum_d(n_dim_idx), gpu_stream)                          


!          print *, " >> storing k_index_d "
          allocate( gpu_host % host( n_dim_idx ) % k_index_h(1:gpu_exp % nk(n_dim_idx) ) )

          call cpy_dtoh_event(&
               gpu_exp % k_index_d(n_dim_idx), &
               c_loc( gpu_host % host( n_dim_idx ) % k_index_h ), &
               gpu_exp % st_k_index_d(n_dim_idx), &
               gpu_stream)
          call gpu_free_async( gpu_exp % k_index_d(n_dim_idx), gpu_stream )
 !         call gpu_stream_sync( gpu_stream )
          
!          print *, " >> storing j2_index_d "
          allocate( gpu_host % host( n_dim_idx ) % j2_index_h(1:gpu_exp % nk(n_dim_idx) ) )

          call cpy_dtoh_event(&
               gpu_exp % j2_index_d(n_dim_idx), &
               c_loc( gpu_host % host( n_dim_idx ) % j2_index_h ), &
               gpu_exp % st_k_index_d(n_dim_idx), &
               gpu_stream)
          call gpu_free_async( gpu_exp % j2_index_d(n_dim_idx), gpu_stream )

!          print *, " >> storing rjs_index_d "
          allocate( gpu_host % host( n_dim_idx ) % rjs_index_h( 1:gpu_exp % nk(n_dim_idx) ) )

          call cpy_dtoh_event(&
               gpu_exp % rjs_index_d(n_dim_idx), &
               c_loc( gpu_host % host( n_dim_idx ) % rjs_index_h ), &
               st_rjs_index_d, &
               gpu_stream)

!          call gpu_stream_sync( gpu_stream )
          ! do i = 1, gpu_exp % nk( n_dim_idx )
          !    print *,  " rjk check: i ", i, " ", gpu_host % host( n_dim_idx ) % rjs_index_h(i)
          ! end do
          
!          print *, " >> storing xyz_k_d "
          allocate( gpu_host % host( n_dim_idx ) % xyz_k_h(1:3, 1:gpu_exp % nk(n_dim_idx) ) )
          call cpy_dtoh_event(&
               gpu_exp % xyz_k_d(n_dim_idx), &
               c_loc( gpu_host % host( n_dim_idx ) % xyz_k_h ), &
               3 * st_rjs_index_d, &
               gpu_stream)
          call gpu_free_async( gpu_exp % xyz_k_d(n_dim_idx), gpu_stream )

          

          gpu_exp % st_pair_distribution_partial_d(n_dim_idx) = n_samples * c_double 
          call gpu_malloc_all(gpu_exp % pair_distribution_partial_d(n_dim_idx), &
               gpu_exp % st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             
          call gpu_memset_async(gpu_exp % pair_distribution_partial_d(n_dim_idx), 0, &
               gpu_exp % st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             


          st_pdf_to_reduce_d = gpu_exp % nk(n_dim_idx) * n_samples * c_double 
          call gpu_malloc_all(pdf_to_reduce_d, st_pdf_to_reduce_d, gpu_stream)             
          call gpu_memset_async(pdf_to_reduce_d, 0, st_pdf_to_reduce_d, gpu_stream)             
          
          
          pdf_factor =  ( ( r_max - r_min) / dfloat(n_samples) ) / ( sqrt( 2.d0 * pi) * kde_sigma)

!          print *, " >> pdf factor = ", pdf_factor
!          der_factor = 0.d0             

   !        if ( j == k ) f = 1.d0
!           if ( j /= k ) f = 2.d0    
    
!           der_factor = v_uc /  n_atoms_of_species(j) / n_atoms_of_species(k) / f

!           print *, " >> der factor = ", der_factor
!           print *, " >> pdf*der factor = ", pdf_factor * der_factor          
          
!            pdf_factor = pdf_factor * der_factor
! ! !
!          call gpu_stream_sync(gpu_stream)
          !          print *, " >> Getting pdf batch"

!          call gpu_meminfo()
          
          call gpu_get_pair_distribution_only_falloc(&
               gpu_exp % pair_distribution_partial_d(n_dim_idx), pdf_to_reduce_d,&
               gpu_exp % nk(n_dim_idx), &
               n_samples, &
               kde_sigma, &
               x_d, dV_d,&
               gpu_exp % rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)


          !--- check x ---! 
          ! allocate( x_check(1:n_samples) )
          ! st_x_d = n_samples * c_double 
          ! call cpy_dtoh_event(&
          !      x_d, &
          !      c_loc(x_check), &
          !      st_x_d, &
          !      gpu_stream)

          ! do l = 1, n_samples
          !    print *, "x_d, l, ", l, " x_d = ", x_check(l)
          ! end do


          ! call cpy_dtoh_event(&
          !      dV_d, &
          !      c_loc(x_check), &
          !      st_x_d, &
          !      gpu_stream)

          ! do l = 1, n_samples
          !    print *, "dV_d, l, ", l, " dV_d = ", x_check(l)
          ! end do
          
          ! deallocate(x_check)


          
          call gpu_free_async( pdf_to_reduce_d, gpu_stream )
!          print *, " >> freeing rjs_index_d "
          call gpu_free_async(gpu_exp % rjs_index_d(n_dim_idx), gpu_stream)
          
!          print *, " >> storing pair_distribution_partial_d "
          allocate( gpu_host % host( n_dim_idx ) % pair_distribution_partial_h( 1:n_samples ) )           
          call cpy_dtoh_event(&
               gpu_exp % pair_distribution_partial_d(n_dim_idx), &
               c_loc( gpu_host % host( n_dim_idx ) % pair_distribution_partial_h ), &
               gpu_exp % st_pair_distribution_partial_d(n_dim_idx), &
               gpu_stream)


          ! do l = 1, n_samples
          !    if ( mod( l, 1 ) == 0)then
          !       print *, " pdf ", n_dim_idx, " l ", l , " ", gpu_host % host( n_dim_idx ) % pair_distribution_partial_h(l)
          !    end if
          ! end do
          
          
          call gpu_free_async( gpu_exp % pair_distribution_partial_d(n_dim_idx), gpu_stream )


                
          
!          call gpu_stream_sync( gpu_stream )
          ! call gpu_copy_pdf( n_samples, gpu_exp % pair_distribution_partial_d(n_dim_idx), &
          !      gpu_host % pair_distribution_partial_h(n_dim)pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
          !      gpu_stream )



          
          n_dim_idx = n_dim_idx + 1

!          call gpu_meminfo()                          

          if ( n_dim_idx > n_dim_partial )then
             exit outer1
          end if

       end do
    end do outer1


    
    
  end subroutine calculate_batched_pair_distribution



  subroutine calculate_batched_pair_distribution_der(gpu_exp,  gpu_host,&
       &  x, dV, n_atoms_of_species, n_species, n_sites,&
       & i_beg, i_end, j_beg, j_end, n_samples,  r_min, r_max, r_cut, kde_sigma, &
       & gpu_stream, x_d, dV_d, j, k, n_dim_idx, v_uc)
    implicit none
    ! Input variables 
    integer, intent(in) :: n_species, i_beg, i_end, j_beg, j_end, n_sites, n_samples, j, k, n_dim_idx
    real*8, allocatable, intent(in) :: n_atoms_of_species(:), x(:), dV(:)
    real*8, intent(in) :: r_min, r_max, r_cut, kde_sigma, v_uc
    type( gpu_storage_type ), intent( inout ) :: gpu_exp 
    type( gpu_host_batch_storage_type ), intent(inout),  target :: gpu_host

    ! Local variables
    real*8, parameter :: pi = acos(-1.0)    
    integer :: n_dim_partial
    type( c_ptr ), intent(in) :: x_d, dV_d
    integer( c_size_t ) :: st_x_d
    integer, target :: nk_temp(1)    
    type( c_ptr ) :: nk_flags_d, nk_sum_flags_d
    integer( c_size_t ) :: st_nk_flags, st_nk_temp
    type( c_ptr ) :: rjs_index_d
    integer( c_size_t ) :: st_rjs_index_d

    real*8 :: pdf_factor, der_factor = 0.d0, f
    integer :: i

    type(c_ptr) :: gpu_stream
    
    n_dim_partial = n_species * ( n_species + 1 ) / 2 

    if( .not. allocated( gpu_exp % st_pair_distribution_partial_der_d ) ) allocate( gpu_exp % st_pair_distribution_partial_der_d(1:n_dim_partial) )
    if( .not. allocated( gpu_exp % pair_distribution_partial_der_d ) ) allocate( gpu_exp % pair_distribution_partial_der_d(1:n_dim_partial) )

    st_rjs_index_d = gpu_exp % nk(n_dim_idx) * c_double
    call gpu_malloc_all(rjs_index_d, st_rjs_index_d, gpu_stream)             
    call gpu_memset_async(rjs_index_d, 0, st_rjs_index_d, gpu_stream)             


    call cpy_htod(&
         c_loc( gpu_host % host( n_dim_idx ) % rjs_index_h ), &
         rjs_index_d, &
         st_rjs_index_d, &
         gpu_stream)


    
    gpu_exp % st_pair_distribution_partial_der_d(n_dim_idx) = n_samples * gpu_exp % nk(n_dim_idx) * c_double              
    call gpu_malloc_all(gpu_exp %pair_distribution_partial_der_d(n_dim_idx), gpu_exp %st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)             
    call gpu_memset_async(gpu_exp %pair_distribution_partial_der_d(n_dim_idx), 0, gpu_exp %st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

    
    pdf_factor =  ( ( r_max - r_min) / &
         dfloat(n_samples) ) / ( sqrt( 2.d0 * pi) * kde_sigma)

    if ( j == k ) f = 1.d0
    if ( j /= k ) f = 2.d0    
    
    der_factor = v_uc /  n_atoms_of_species(j) / &
         & n_atoms_of_species(k) / f
    
!    print *, " starting pdf der only "
    call gpu_get_pair_distribution_der_only(&
         gpu_exp % pair_distribution_partial_der_d(n_dim_idx),&
         gpu_exp % nk(n_dim_idx), &
         n_samples, &
         kde_sigma, &
         x_d, dV_d,&
         rjs_index_d, pdf_factor, der_factor, gpu_stream)
    
!    print *, " finished pdf der only "
    call gpu_free_async(rjs_index_d, gpu_stream)


          
    ! call gpu_free_async(x_d, gpu_stream)
    ! call gpu_free_async(dV_d)                    

    
    
  end subroutine calculate_batched_pair_distribution_der


  subroutine setup_gpu_xrd_forces( gpu_exp, gpu_host, n_dim_idx, gpu_stream )
    implicit none
    integer, intent(in) :: n_dim_idx 
    type( gpu_storage_type ), intent( inout ) :: gpu_exp 
    type( gpu_host_batch_storage_type ), intent(inout), target :: gpu_host
    integer(c_size_t) :: st_rjs_index_d
    type(c_ptr) :: gpu_stream
    ! copy the xyz, j2 and k_index_d arrays 

    st_rjs_index_d = gpu_exp % nk( n_dim_idx ) * c_double 
    call gpu_malloc_all(gpu_exp % k_index_d(n_dim_idx),      gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             
    call gpu_malloc_all(gpu_exp % j2_index_d(n_dim_idx),      gpu_exp % st_k_index_d(n_dim_idx), gpu_stream)             
    call gpu_malloc_all(gpu_exp % xyz_k_d(n_dim_idx), 3*st_rjs_index_d, gpu_stream)             


    call cpy_htod(&
         c_loc( gpu_host % host( n_dim_idx ) % k_index_h ), &
         gpu_exp % k_index_d(n_dim_idx), &
         gpu_exp % st_k_index_d(n_dim_idx), &
         gpu_stream)


    call cpy_htod(&
         c_loc( gpu_host % host( n_dim_idx ) % j2_index_h ), &
         gpu_exp % j2_index_d(n_dim_idx), &
         gpu_exp % st_k_index_d(n_dim_idx), &
         gpu_stream)

!    st_rjs_index_d = gpu_exp % nk( n_dim_idx ) * c_double 

    call cpy_htod(&
         c_loc( gpu_host % host( n_dim_idx ) % xyz_k_h ), &
         gpu_exp % xyz_k_d(n_dim_idx), &
         3 * st_rjs_index_d, &
         gpu_stream)
          
    
  end subroutine setup_gpu_xrd_forces


  subroutine free_gpu_xrd_forces( gpu_exp, gpu_host, n_dim_idx, gpu_stream )
    implicit none
    integer, intent(in) :: n_dim_idx 
    type( gpu_storage_type ), intent( inout ) :: gpu_exp 
    type( gpu_host_batch_storage_type ), intent(inout), target :: gpu_host
    integer(c_size_t) :: st_rjs_index_d
    type(c_ptr) :: gpu_stream    
    ! copy the xyz, j2 and k_index_d arrays 

    call gpu_free_async(gpu_exp % k_index_d(n_dim_idx), gpu_stream)             
    call gpu_free_async(gpu_exp % j2_index_d(n_dim_idx), gpu_stream)             
    call gpu_free_async(gpu_exp % xyz_k_d(n_dim_idx), gpu_stream)
    call gpu_free_async(gpu_exp % pair_distribution_partial_der_d(n_dim_idx), gpu_stream)                 
    
  end subroutine free_gpu_xrd_forces
    
    
  
  !--- GPU PAIR DISTRIBUTION FUNCTIONS ---!
  subroutine gpu_calculate_pair_distribution(n_dim_partial_out,  params, x_pair_distribution&
       &, y_pair_distribution, y_pair_distribution_temp,&
       & pair_distribution_partial, pair_distribution_partial_temp, &
       & n_species, species_types,  n_atoms_of_species, n_sites, a_box, b_box, c_box,&
       & indices, md_istep, mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs, xyz, &
       & neighbors_list, n_neigh, neighbor_species, species, rank,&
       & do_derivatives, pair_distribution_der, pair_distribution_partial_der, &
       & nk, pair_distribution_d, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
       & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d,&
       & n_neigh_d, species_d,  neighbor_species_d, neighbor_list_d, rjs_d, xyz_d, species_types_d, cublas_handle, gpu_stream,&
       gpu_host_storage, gpu_low_memory,&
       & pair_distribution_partial_temp_der, energies_pair_distribution, forces_pair_distribution, virial)
    implicit none
    type(input_parameters), intent(inout) :: params
    integer, intent(out) :: n_dim_partial_out
    real*8, allocatable, intent(out), target :: x_pair_distribution(:),&
         & y_pair_distribution(:), pair_distribution_partial(:,:),&
         & n_atoms_of_species(:), pair_distribution_partial_temp(:,:),&
         & y_pair_distribution_temp(:), pair_distribution_der(:,:),&
         & pair_distribution_partial_der(:,:,:), &
         & pair_distribution_partial_temp_der(:,:,:), energies_pair_distribution(:), forces_pair_distribution(:,:)
    character*8, allocatable, intent(in) :: species_types(:)
    real*8,  intent(in), allocatable, target :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable, target :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(inout) :: virial(1:3,1:3)
    real*8 :: v_uc, f, pdf_factor, der_factor, total_memory_usage=0.d0
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end
    integer, intent(in) :: indices(1:3), md_istep, mc_istep,  rank
    integer, intent(inout) :: ierr
    real*8, allocatable, target :: factors(:), pair_distribution_der_temp(:), dV(:), &
         pdf_gpu_check(:), rjs_temp(:), ders_temp(:,:)
    integer, allocatable, target :: ks_temp(:), ksd_temp(:)
    integer :: i, j, k, l, i2, n_dim_partial, n_dim_idx
    logical, intent(in) :: do_derivatives
    real*8, parameter :: pi = acos(-1.0)
    logical :: write_condition, overwrite_condition
    character*1024 :: filename
    type(c_ptr) :: cublas_handle, gpu_stream


    type(c_ptr) :: n_neigh_d, species_d, neighbor_species_d, neighbor_list_d, rjs_d, xyz_d, species_types_d, x_d, dV_d
    type(c_ptr) :: pair_distribution_d
    type(c_ptr), allocatable :: nk_d(:), nk_flags_d(:), nk_flags_sum_d(:), k_index_d(:), j2_index_d(:), rjs_index_d(:), xyz_k_d(:), pair_distribution_partial_d(:), pair_distribution_partial_der_d(:)
    integer(c_size_t), allocatable :: st_nk_d(:), st_k_index_d(:), st_j2_index_d(:), st_rjs(:), st_pair_distribution_partial_d(:), st_pair_distribution_partial_der_d(:)    
    integer(c_size_t) :: st_nk_flags, st_nk_temp, st_x_d
    integer, allocatable :: nk(:), k_index_single(:)
    integer, target :: nk_temp(1)
    
    type( gpu_host_storage_type ), intent(inout),  allocatable, target :: gpu_host_storage(:)
    logical, intent(in) :: gpu_low_memory
    integer(c_size_t) :: st_n_sites_int, st_n_atom_pairs_int, st_n_atom_pairs_double, st_species_types_d
    
!            print *, ""
!            print *, " >> Allocating GPU arrays for Exp Calculation << "           
!            print *, ""           
!            st_n_sites_int = n_sites*sizeof(n_neigh(1)) 
!            call gpu_malloc_all(n_neigh_d,st_n_sites_int,gpu_stream)
!            call cpy_htod(c_loc(n_neigh),n_neigh_d, st_n_sites_int,gpu_stream)
!            call gpu_malloc_all(species_d,st_n_sites_int,gpu_stream)
!            call cpy_htod(c_loc(species),species_d, st_n_sites_int,gpu_stream)
!            st_n_atom_pairs_int = j_end * sizeof(neighbor_species(1))

!            call gpu_malloc_all(neighbor_species_d,st_n_atom_pairs_int,gpu_stream)
!            call cpy_htod(c_loc(neighbor_species),neighbor_species_d, st_n_atom_pairs_int,gpu_stream)
!            call gpu_malloc_all(neighbors_list_d,st_n_atom_pairs_int,gpu_stream)
!            print *, " -- n_pairs for neighbor_list = ", j_end           
!            call cpy_htod(c_loc(neighbors_list),neighbors_list_d, st_n_atom_pairs_int,gpu_stream)
           
!            st_n_atom_pairs_double = j_end*sizeof(rjs(1))
!            call gpu_malloc_all(rjs_d,st_n_atom_pairs_double,gpu_stream)
!            call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)
!            call gpu_malloc_all(xyz_d,3*st_n_atom_pairs_double,gpu_stream)
!            call cpy_htod(c_loc(xyz),xyz_d,3*st_n_atom_pairs_double,gpu_stream)

!            st_species_types_d = n_species * c_int
!            call gpu_malloc_all(species_types_d,st_species_types_d,gpu_stream)
!            call cpy_htod(c_loc(species_types),species_types_d,st_species_types_d,gpu_stream)
! !           call gpu_device_sync()

    
    
    ! Seeing if my gpu helper mod actually helps
    ! type( gpu_var_double_class ) :: xc, dVc    
    ! type( gpu_var_double_class ), allocatable :: pdf_partial_d(:)

    ! type( gpu_var_int_class ),    allocatable :: nk_d(:), nk_flags_d(:), nk_flags_sum_d(:), k_index_d(:), j2_index_d(:)
    ! type( gpu_var_double_class ), allocatable :: rjs_index_d(:), xyz_k_d(:), pair_distribution_partial_d(:), pair_distribution_partial_der_d(:)

    
    
    ! Things that are allocated here:
    ! Always:
    !  > x_pair_distribution
    !  > y_pair_distribution
    ! if pair_distribution_partial == .true.
    !  > pair_distribution_partial( n_samples, n_spec * (n_spec + 1)/2 )
    !  if do_derivatives == .true.
    !    > pair_distribution_partial_der( n_samples, n_spec * (n_spec + 1)/2, j_beg : j_end )


    ! first allocate the necessary arrays for the
    ! calculation of the pair correlation function
    if (allocated( x_pair_distribution)) deallocate(x_pair_distribution)
    if (allocated( y_pair_distribution)) deallocate(y_pair_distribution)

    allocate( x_pair_distribution( 1: params%pair_distribution_n_samples) )
    allocate( y_pair_distribution( 1: params%pair_distribution_n_samples) )

    allocate( dV( 1: params%pair_distribution_n_samples) )    
    
    if (params%n_exp > 0)then
       do i = 1, params%n_exp
          if ( trim( params%exp_data(i)%label ) == 'pair_distribution' )then
             x_pair_distribution = params%exp_data(i)%x
          end if
       end do
    end if

    n_dim_partial_out = 0

    if (params%pair_distribution_partial)then
       n_dim_partial = n_species * ( n_species + 1 ) / 2
       n_dim_partial_out = n_dim_partial       
       allocate(factors( 1:n_dim_partial ))

       n_dim_idx = 1
       outer: do i = 1, n_species
          do j = 1, n_species
             if (i > j) cycle

             if (i /= j)then
                factors(n_dim_idx) = 2.d0
             else
                factors(n_dim_idx) = 1.d0
             end if

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial )then
                exit outer
             end if

          end do
       end do outer


       if (.not. allocated(pair_distribution_partial))then   !deallocate(pair_distribution_partial)
          allocate( pair_distribution_partial(1:params%pair_distribution_n_samples,&
               & 1 : n_dim_partial) )
       end if

       pair_distribution_partial = 0.d0


       if (params%do_forces .and. params%exp_forces)then
          allocate( pair_distribution_partial_der(1:params%pair_distribution_n_samples,&
            & 1 : n_dim_partial, j_beg : j_end  ))
          pair_distribution_partial_der = 0.d0

          if (rank == 0 .and. md_istep == 0) write(*, '(A,1X,F7.4,1X,A)') "Gb/core: partial pdfder = ", dfloat(params&
               &%pair_distribution_n_samples * n_dim_partial * j_end)&
               & * 8.d0 / (dfloat(1024*1024*1024)), " Gb  |"
          if (rank == 0 .and. md_istep == 0) write(*,*)'                                       |'

       end if
    else
       if (params%do_forces .and. params%exp_forces)then
          allocate( pair_distribution_partial_der(1:params%pair_distribution_n_samples, 1:1, &
            &  j_beg : j_end  ))
          pair_distribution_partial_der = 0.d0
       end if


    end if

    if(allocated(n_atoms_of_species)) deallocate(n_atoms_of_species)
    allocate(n_atoms_of_species(1:n_species))

    do j = 1, n_species
       n_atoms_of_species(j) = 0.d0
       do i2 = 1, n_sites
          if ( species(i2) == j)then
             n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
          end if
       end do
    end do


#ifdef _MPIF90
    if (params%pair_distribution_partial)then
       allocate( pair_distribution_partial_temp(1:params%pair_distribution_n_samples, 1 : n_dim_partial) )

       pair_distribution_partial_temp = 0.0d0
    end if

    allocate( y_pair_distribution_temp( 1: params%pair_distribution_n_samples) )
    y_pair_distribution_temp = 0.d0

#endif
    v_uc = dot_product( cross_product(a_box,&
         & b_box), c_box ) / (&
         & dfloat(indices(1)*indices(2)&
         &*indices(3)) )


    !#####################################################################!
    !###---   Calculating the partial pair distribution functions   ---###!
    !#####################################################################!
!!    print *, " - Allocating pdf gpu pointers -"
    allocate( nk_d(1:n_dim_partial) )
    allocate( nk_flags_d(1:n_dim_partial) )
    allocate( nk_flags_sum_d(1:n_dim_partial) )        
    allocate( k_index_d(1:n_dim_partial) )
    allocate( j2_index_d(1:n_dim_partial) )
    allocate( rjs_index_d(1:n_dim_partial) )    
    allocate( xyz_k_d(1:n_dim_partial) )    
    allocate( pair_distribution_partial_d(1:n_dim_partial) )
    allocate( pair_distribution_partial_der_d(1:n_dim_partial) )
    allocate( st_nk_d(1:n_dim_partial) )
    allocate( st_k_index_d(1:n_dim_partial) )
    allocate( st_j2_index_d(1:n_dim_partial) )
    allocate( st_rjs(1:n_dim_partial) )
    allocate( st_pair_distribution_partial_d(1:n_dim_partial) )
    allocate( st_pair_distribution_partial_der_d(1:n_dim_partial) )
    allocate( nk(1:n_dim_partial) )    
    
    st_x_d = params%pair_distribution_n_samples * c_double

    call setup_pdf_arrays( params%r_range_min, params%r_range_max,&
         & params%pair_distribution_rcut, params&
         &%pair_distribution_n_samples, x_pair_distribution, dV)

    call gpu_meminfo()    
    call gpu_malloc_all(x_d,      st_x_d, gpu_stream)             
    call cpy_htod( c_loc( x_pair_distribution ), x_d, st_x_d, gpu_stream )

    call gpu_malloc_all(dV_d,      st_x_d, gpu_stream)             
    call cpy_htod( c_loc( dV ), dV_d, st_x_d, gpu_stream )
    
    if ( params%pair_distribution_partial )then
       if( gpu_low_memory )then

          ! Allocate the host storage arrays

          allocate( gpu_host_storage( 1:n_dim_partial ) )
          
          
          n_dim_idx = 1
          outera: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle ! We have already calculated the pair correlation function!

                ! Note that with the calculation of the derivatives here,
                ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
                ! factor, which allows for some freeing of memory

                ! Get all nk_flags 
                call gpu_meminfo()

!                print *, "pdfpairs:  j ", j, " k ",  k
                call gpu_stream_sync(gpu_stream)
                
                st_nk_temp = 1*c_int
                call gpu_malloc_all(nk_d(n_dim_idx),st_nk_temp, gpu_stream)                          
                st_nk_flags = j_end * c_int
                call gpu_malloc_all(nk_flags_d(n_dim_idx),      st_nk_flags, gpu_stream)             
                call gpu_memset_async(nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
                call gpu_malloc_all(nk_flags_sum_d(n_dim_idx),      st_nk_flags, gpu_stream)                          
!                call gpu_stream_sync(gpu_stream)


                ! Here I need the neighbor_list_d, n_neigh_d, neighbor_species_d, species_d, rjs_d

                st_n_sites_int = n_sites*sizeof(n_neigh(1)) 
                call gpu_malloc_all(n_neigh_d,st_n_sites_int,gpu_stream)
                call cpy_htod(c_loc(n_neigh),n_neigh_d, st_n_sites_int,gpu_stream)
                call gpu_malloc_all(species_d,st_n_sites_int,gpu_stream)
                call cpy_htod(c_loc(species),species_d, st_n_sites_int,gpu_stream)
                st_n_atom_pairs_int = j_end * sizeof(neighbor_species(1))

                call gpu_malloc_all(neighbor_species_d,st_n_atom_pairs_int,gpu_stream)
                call cpy_htod(c_loc(neighbor_species),neighbor_species_d, st_n_atom_pairs_int,gpu_stream)
                call gpu_malloc_all(neighbor_list_d,st_n_atom_pairs_int,gpu_stream)
                call cpy_htod(c_loc(neighbors_list),neighbor_list_d, st_n_atom_pairs_int,gpu_stream)

                st_n_atom_pairs_double = j_end*sizeof(rjs(1))
                call gpu_malloc_all(rjs_d,st_n_atom_pairs_double,gpu_stream)
                call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)
           ! call gpu_malloc_all(xyz_d,3*st_n_atom_pairs_double,gpu_stream)
           ! call cpy_htod(c_loc(xyz),xyz_d,3*st_n_atom_pairs_double,gpu_stream)

           ! st_species_types_d = n_species * c_int
           ! call gpu_malloc_all(species_types_d,st_species_types_d,gpu_stream)
           ! call cpy_htod(c_loc(species_types),species_types_d,st_species_types_d,gpu_stream)

                
                ! I don't actually need xyz_d here at all!
                call gpu_get_pair_distribution_nk(i_beg, i_end, j_end, n_sites, neighbor_list_d,&
                     n_neigh_d, neighbor_species_d, species_d,&
                     & rjs_d, xyz_d, params%r_range_min, params%r_range_max, params%pair_distribution_rcut, 6.d0&
                     &*params%pair_distribution_kde_sigma,&
                     & nk_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), j, k, gpu_stream)


                call gpu_free_async(n_neigh_d, gpu_stream)
                call gpu_free_async(species_d, gpu_stream)
                call gpu_free_async(neighbor_species_d, gpu_stream)
                call gpu_free_async(neighbor_list_d, gpu_stream)
                call gpu_free_async(rjs_d, gpu_stream)
                
                !                call gpu_free_async(species_types_actual_d,gpu_stream)
                
                
                call gpu_free(nk_flags_d(n_dim_idx))
!                call gpu_free_async(nk_flags_d(n_dim_idx), gpu_stream)
                
!                print *, "out of pdf nk kernel "
                st_nk_temp = 1*c_int
                call cpy_dtoh(nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)
!                call gpu_stream_sync(gpu_stream)

                nk(n_dim_idx) = nk_temp(1)
!                print *, " nk = ", nk(n_dim_idx)


                ! call estimate_device_memory_usage( n_sites, nk(n_dim_idx), params%pair_distribution_n_samples,&
                !      params%structure_factor_n_samples, total_memory_usage )


                ! Now we create temporary arrays for the k indices
                ! ---------------------------------------------------
                ! -------------------- Setting k --------------------
                ! ---------------------------------------------------                
                st_k_index_d(n_dim_idx) = nk(n_dim_idx) * c_int
                call gpu_malloc_all(k_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(k_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)             

                call gpu_set_pair_distribution_k_index_only(j_end, k_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream )

                allocate( gpu_host_storage( n_dim_idx ) % k_index_h( 1:nk(n_dim_idx) ) )
                call cpy_dtoh( k_index_d(n_dim_idx), c_loc( gpu_host_storage( n_dim_idx ) % k_index_h ), st_k_index_d(n_dim_idx), gpu_stream )
                call gpu_free( k_index_d(n_dim_idx) )
!                call gpu_stream_sync( gpu_stream )


                ! ----------------------------------------------------
                ! -------------------- Setting j2 --------------------
                ! ----------------------------------------------------               
                call gpu_malloc_all(j2_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(j2_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)             

                call gpu_malloc_all(neighbor_list_d,st_n_atom_pairs_int,gpu_stream)
                call cpy_htod(c_loc(neighbors_list),neighbor_list_d, st_n_atom_pairs_int,gpu_stream)

                call gpu_set_pair_distribution_j2_only(j_end, n_sites, neighbor_list_d, j2_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream )

                call gpu_free( neighbor_list_d )
                allocate( gpu_host_storage( n_dim_idx ) % j2_index_h( 1:nk(n_dim_idx) ) )

                call cpy_dtoh( j2_index_d(n_dim_idx), c_loc( gpu_host_storage( n_dim_idx ) % j2_index_h ), st_k_index_d(n_dim_idx),&
                     gpu_stream )
                call gpu_free( j2_index_d(n_dim_idx) )
!                call gpu_stream_sync( gpu_stream )


                ! -----------------------------------------------------
                ! -------------------- Setting xyz --------------------
                ! -----------------------------------------------------
                call gpu_malloc_all(xyz_d,3*st_n_atom_pairs_double,gpu_stream)
                call cpy_htod(c_loc(xyz),xyz_d,3*st_n_atom_pairs_double,gpu_stream)
                
                st_rjs(n_dim_idx) = nk(n_dim_idx) * c_double                
                call gpu_malloc_all(xyz_k_d(n_dim_idx), 3*st_rjs(n_dim_idx), gpu_stream)             
                call gpu_memset_async(xyz_k_d(n_dim_idx), 0, 3*st_rjs(n_dim_idx), gpu_stream)             

                call gpu_set_pair_distribution_xyz_only(j_end, xyz_d, xyz_k_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream )

                allocate( gpu_host_storage( n_dim_idx ) % xyz_k_h(1:3, 1:nk(n_dim_idx) ) )
                call cpy_dtoh( xyz_k_d(n_dim_idx), c_loc( gpu_host_storage( n_dim_idx ) % xyz_k_h ), 3*st_rjs(n_dim_idx), gpu_stream)
                call gpu_free_async(xyz_d, gpu_stream)
                call gpu_stream_sync(gpu_stream)
                
                


                ! -----------------------------------------------------
                ! -------------------- Setting rjs --------------------
                ! -----------------------------------------------------
                st_n_atom_pairs_double = j_end*sizeof(rjs(1))
                call gpu_malloc_all(rjs_d,st_n_atom_pairs_double,gpu_stream)
                call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)
                
                st_rjs(n_dim_idx) = nk(n_dim_idx) * c_double                                
                call gpu_malloc_all(rjs_index_d(n_dim_idx), st_rjs(n_dim_idx), gpu_stream)             
                call gpu_memset_async(rjs_index_d(n_dim_idx), 0, st_rjs(n_dim_idx), gpu_stream)             

                call gpu_set_pair_distribution_rjs_only(j_end, rjs_d, rjs_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream )

                call gpu_free_async(rjs_d, gpu_stream)
                
                ! No need for device storage here 

!                call gpu_stream_sync(gpu_stream)


                call gpu_meminfo()             
                call gpu_free_async(nk_flags_sum_d(n_dim_idx), gpu_stream)                          
                call gpu_stream_sync(gpu_stream)
                
!                print *, "Setting pdf device arrays kernel "
                st_pair_distribution_partial_d(n_dim_idx) = params%pair_distribution_n_samples * c_double 
                call gpu_malloc_all(pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(pair_distribution_partial_d(n_dim_idx), 0, st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             

                pdf_factor =  ( ( params%r_range_max - params%r_range_min) / &
                     dfloat(params%pair_distribution_n_samples) ) / ( sqrt( 2.d0 * pi) * params%pair_distribution_kde_sigma)

                der_factor = v_uc /  n_atoms_of_species(j) / &
                     & n_atoms_of_species(k) / factors(n_dim_idx)             

                call gpu_stream_sync(gpu_stream)
                
!                print *, " starting pdf only "
                call gpu_get_pair_distribution_only(&
                     pair_distribution_partial_d(n_dim_idx),&
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d,&
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                call gpu_copy_pdf( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                     pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                     gpu_stream )

                call gpu_free_async( pair_distribution_partial_d(n_dim_idx), gpu_stream )
                

                call gpu_stream_sync(gpu_stream)
                
                st_pair_distribution_partial_der_d(n_dim_idx) = params%pair_distribution_n_samples * nk(n_dim_idx) * c_double              
                call gpu_malloc_all(pair_distribution_partial_der_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(pair_distribution_partial_der_d(n_dim_idx), 0, st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

!                print *, " starting pdf der only "
                call gpu_get_pair_distribution_der_only(&
                     pair_distribution_partial_der_d(n_dim_idx),&
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d,&
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                
!                call gpu_stream_sync(gpu_stream)


                ! call gpu_get_pair_distribution_and_ders(&
                !      pair_distribution_partial_d(n_dim_idx),&
                !      pair_distribution_partial_der_d(n_dim_idx),&
                !      nk(n_dim_idx), &
                !      params%pair_distribution_n_samples, &
                !      params%pair_distribution_kde_sigma, &
                !      x_d, dV_d,&
                !      rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                ! We can deallocate the rjs as we don't need them any more
                call gpu_free_async(rjs_index_d(n_dim_idx), gpu_stream)

                call gpu_stream_sync(gpu_stream)
                
                allocate( gpu_host_storage( n_dim_idx ) % pair_distribution_partial_der_h( 1:params%pair_distribution_n_samples, 1:nk(n_dim_idx) ) )
                call cpy_dtoh( pair_distribution_partial_der_d(n_dim_idx), c_loc( gpu_host_storage( n_dim_idx ) % pair_distribution_partial_der_h ), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream )
                call gpu_free_async( pair_distribution_partial_der_d(n_dim_idx), gpu_stream )

                

                ! call gpu_copy_pdf_der( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                !      pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                !      gpu_stream )


                
                
!                call gpu_stream_sync(gpu_stream)

                n_dim_idx = n_dim_idx + 1

                call gpu_meminfo()                          

                if ( n_dim_idx > n_dim_partial )then
                   exit outera
                end if

             end do
          end do outera

       else
          n_dim_idx = 1
          call gpu_stream_sync( gpu_stream )
          outer1: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle ! We have already calculated the pair correlation function!

                ! Note that with the calculation of the derivatives here,
                ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
                ! factor, which allows for some freeing of memory

                ! Get all nk_flags 
                call gpu_meminfo()

!                print *, "pdfpairs:  j ", j, " k ",  k
                st_nk_temp = 1*c_int
                call gpu_malloc_all(nk_d(n_dim_idx),st_nk_temp, gpu_stream)                          
                st_nk_flags = j_end * c_int
                call gpu_malloc_all(nk_flags_d(n_dim_idx),      st_nk_flags, gpu_stream)             
                call gpu_memset_async(nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
                call gpu_malloc_all(nk_flags_sum_d(n_dim_idx),      st_nk_flags, gpu_stream)                          
!                call gpu_stream_sync(gpu_stream)

                call gpu_get_pair_distribution_nk(i_beg, i_end, j_end, n_sites, neighbor_list_d,&
                     n_neigh_d, neighbor_species_d, species_d,&
                     & rjs_d, xyz_d, params%r_range_min, params%r_range_max, params%pair_distribution_rcut, 6.d0&
                     &*params%pair_distribution_kde_sigma,&
                     & nk_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), j, k, gpu_stream)
                !             call gpu_device_sync()             

                call gpu_free_async(nk_flags_d(n_dim_idx), gpu_stream)
                
                ! Now copy the value of nk from the gpu
!                print *, "out of pdf nk kernel "
                st_nk_temp = 1*c_int
                call cpy_dtoh(nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)
!                call gpu_stream_sync(gpu_stream)

                nk(n_dim_idx) = nk_temp(1)
!                print *, " nk = ", nk(n_dim_idx)



                call estimate_device_memory_usage( n_sites, 0, nk(n_dim_idx), params%pair_distribution_n_samples,&
                     params%structure_factor_n_samples, total_memory_usage, .false. )


                ! Now we create temporary arrays for the k indices

                st_k_index_d(n_dim_idx) = nk(n_dim_idx) * c_int
                call gpu_malloc_all(k_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(k_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)             

                call gpu_malloc_all(j2_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(j2_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)             

                st_rjs(n_dim_idx) = nk(n_dim_idx) * c_double
                call gpu_malloc_all(rjs_index_d(n_dim_idx), st_rjs(n_dim_idx), gpu_stream)             
                call gpu_memset_async(rjs_index_d(n_dim_idx), 0, st_rjs(n_dim_idx), gpu_stream)             


                call gpu_malloc_all(xyz_k_d(n_dim_idx), 3*st_rjs(n_dim_idx), gpu_stream)             
                call gpu_memset_async(xyz_k_d(n_dim_idx), 0, 3*st_rjs(n_dim_idx), gpu_stream)             
 !               call gpu_stream_sync(gpu_stream)

                call gpu_meminfo()             
                call gpu_set_pair_distribution_k_index(i_beg, i_end, j_end, n_sites, neighbor_list_d,&
                     & rjs_d, xyz_d, k_index_d(n_dim_idx), j2_index_d(n_dim_idx),&
                     & rjs_index_d(n_dim_idx), xyz_k_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx),&
                     & gpu_stream)
                !             call gpu_device_sync()             


                call gpu_free_async(nk_flags_sum_d(n_dim_idx), gpu_stream)                          


                !             call gpu_meminfo()


                !--- CALCULATING THE PAIR DISTRIBUTION FUNCTION ---!
                !             print *, " Allocating pdf arrays   "             
                st_pair_distribution_partial_d(n_dim_idx) = params%pair_distribution_n_samples * c_double 
                call gpu_malloc_all(pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(pair_distribution_partial_d(n_dim_idx), 0, st_pair_distribution_partial_d(n_dim_idx), gpu_stream)             

                st_pair_distribution_partial_der_d(n_dim_idx) = params%pair_distribution_n_samples * nk(n_dim_idx) * c_double              
                call gpu_malloc_all(pair_distribution_partial_der_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)             
                call gpu_memset_async(pair_distribution_partial_der_d(n_dim_idx), 0, st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

                pdf_factor =  ( ( params%r_range_max - params%r_range_min) / &
                     dfloat(params%pair_distribution_n_samples) ) / ( sqrt( 2.d0 * pi) * params%pair_distribution_kde_sigma)

                der_factor = v_uc /  n_atoms_of_species(j) / &
                     & n_atoms_of_species(k) / factors(n_dim_idx)             

!                call gpu_stream_sync(gpu_stream)


                call gpu_get_pair_distribution_and_ders(&
                     pair_distribution_partial_d(n_dim_idx),&
                     pair_distribution_partial_der_d(n_dim_idx),&
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d,&
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)


                ! call get_pair_distribution( n_sites,&
                !      & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end)&
                !      &, neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
                !      & xyz(1:3,j_beg:j_end), params %r_range_min, params&
                !      &%r_range_max, params%pair_distribution_n_samples,&
                !      & x_pair_distribution,&
                !      & pair_distribution_partial(1:params&
                !      &%pair_distribution_n_samples, n_dim_idx), params &
                !      &%pair_distribution_rcut, .false., params&
                !      &%pair_distribution_partial, j, k, params&
                !      &%pair_distribution_kde_sigma, dfloat(n_sites)/v_uc&
                !      &,  params%exp_forces,&
                !      & pair_distribution_partial_der, n_dim_idx, j_beg,&
                !      & j_end )




                ! -------- PDF CHECK ----------!
                !--- CHECKING THAT IT WORKS ---!
                ! print *, "checking pdf"
                ! allocate(pdf_gpu_check(1:params%pair_distribution_n_samples))

                ! call cpy_dtoh( pair_distribution_partial_d(n_dim_idx), c_loc(pdf_gpu_check), &
                !      st_pair_distribution_partial_d(n_dim_idx), gpu_stream )
                ! call gpu_device_sync()

                ! do l = 1, params%pair_distribution_n_samples
                !    print *, " pdfcheck, l = ", l, " gpu ", pdf_gpu_check(l), " cpu ", pair_distribution_partial(l, n_dim_idx) 
                ! end do

                ! deallocate(pdf_gpu_check)

                !------------------------------!
                !-------- PDF DER CHECK -------!

                ! allocate( ks_temp(1:nk(n_dim_idx)) )
                ! call cpy_dtoh( k_index_d(n_dim_idx), c_loc(ks_temp), &
                !      st_k_index_d(n_dim_idx), gpu_stream )             
                ! allocate(ders_temp(1:params%pair_distribution_n_samples, 1 : nk(n_dim_idx) ))
                ! call cpy_dtoh( pair_distribution_partial_der_d(n_dim_idx), c_loc(ders_temp), &
                !      st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream )
                ! call gpu_device_sync()
                ! do l = 1, params%pair_distribution_n_samples
                !    do i = 1, nk(n_dim_idx)
                !       i2 = ks_temp(i)+1

                !       if ( modulo( i , 10000 ) == 0 )then
                !          print *, " pdfdercheck, l = ", l, "i ", i, " i2 ", i2,  " gpu ", ders_temp(l,i), " cpu ", pair_distribution_partial_der(l, n_dim_idx, i2  )
                !       end if
                !    end do
                ! end do
                ! deallocate(ders_temp)
                ! deallocate(ks_temp)             

                ! call cpy_dtoh( pair_distribution_partial_d(n_dim_idx), c_loc(pair_distribution_partial(), &
                !      st_pair_distribution_partial_d(n_dim_idx), gpu_stream )             

                ! will need to do a separate function to copy the pair distribution
                !             call gpu_device_sync()
                call gpu_copy_pdf( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                     pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                     gpu_stream )


                ! call gpu_copy_pdf_der( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                !      pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                !      gpu_stream )


                ! We can deallocate the rjs as we don't need them any more
                call gpu_free_async(rjs_index_d(n_dim_idx), gpu_stream)
!                call gpu_stream_sync(gpu_stream)

                n_dim_idx = n_dim_idx + 1

                call gpu_meminfo()                          

                if ( n_dim_idx > n_dim_partial )then
                   exit outer1
                end if

             end do
          end do outer1
       end if
       call gpu_free_async(x_d, gpu_stream)
       call gpu_free(dV_d)                    
       
       deallocate( rjs_index_d, st_rjs )
       deallocate( nk_flags_d )
       deallocate( nk_flags_sum_d )        
       call gpu_meminfo()             
       
    else
       call get_pair_distribution( n_sites, &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end),&
            & params%r_range_min, params%r_range_max, params &
            &%pair_distribution_n_samples, x_pair_distribution,&
            & y_pair_distribution, params &
            &%pair_distribution_rcut, .false., .false., 1, 1,&
            & params%pair_distribution_kde_sigma, dfloat(n_sites)&
            &/v_uc, params%do_forces .and. params%exp_forces, pair_distribution_partial_der, 1, &
            & j_beg, j_end)
    end if

    deallocate(dV)    

    ! --- MPI communication is here  ---

    if ( params%pair_distribution_partial )then
#ifdef _MPIF90
       call mpi_reduce(pair_distribution_partial,&
            & pair_distribution_partial_temp, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
            & 0, MPI_COMM_WORLD, ierr)

       ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! Note, we have only so far divided by 4 pi r^2 dr
       ! Therefore, we must scale by the density

       pair_distribution_partial =  pair_distribution_partial_temp
       deallocate( pair_distribution_partial_temp )


       call mpi_bcast(pair_distribution_partial, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)


       ! Now, we have the derivatives of the partial pair distribution
       ! function with respect to the atom pairs in that rank
       !
       ! We can keep them in the rank and calculate forces


       ! call mpi_reduce(pair_distribution_partial_der,&
       !      & pair_distribution_partial_der_temp, params&
       !      &%pair_distribution_n_samples * n_species *&
       !      & n_species * 3 * n_pairs_tot, MPI_DOUBLE_PRECISION, MPI_SUM,&
       !      & 0, MPI_COMM_WORLD, ierr)

       ! ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! ! Note, we have only so far divided by 4 pi r^2 dr
       ! ! Therefore, we must scale by the density

       ! pair_distribution_partial_der =  pair_distribution_partial_der_temp
       ! deallocate( pair_distribution_partial_der_temp )



#endif


       if ( params%valid_pdf )then
          allocate(energies_pair_distribution(1:n_sites))
          energies_pair_distribution = 0.d0

          if (params%do_forces .and. params%exp_forces)then
             allocate(forces_pair_distribution(1:3,1:n_sites))
             forces_pair_distribution = 0.d0
          end if

       end if


       !####################################!
       !###---   Accumulate the PDF   ---###!
       !####################################!


       y_pair_distribution = 0.d0
       n_dim_idx = 1
       outer2: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle

             pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) =&
                  & pair_distribution_partial(1:params&
                  &%pair_distribution_n_samples, n_dim_idx) * v_uc &
                  &  /  n_atoms_of_species(j) /  n_atoms_of_species(k) / factors(n_dim_idx) !real(n_sites)

             if (params%do_forces .and. params%exp_forces) then
                pair_distribution_partial_der(1:params&
                     &%pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) =  pair_distribution_partial_der(1:params&
                     & %pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) * v_uc /  n_atoms_of_species(j) / &
                     & n_atoms_of_species(k) / factors(n_dim_idx)!real(n_sites)

             end if

             y_pair_distribution(1:params%pair_distribution_n_samples) = &
                  & y_pair_distribution(1:params%pair_distribution_n_samples)  +  &
                  &  factors(n_dim_idx) * (n_atoms_of_species(j) * n_atoms_of_species(k)) * &
                  & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) &
                  &  /  dfloat(n_sites) / dfloat(n_sites)

             n_dim_idx = n_dim_idx + 1

             if ( n_dim_idx > n_dim_partial )then
                exit outer2
             end if

          end do
       end do outer2

       ! --- Preprocess the pair distribution according to the output --- !
       if ( trim( params%pair_distribution_output ) == "D(r)" )then
          y_pair_distribution = 4.d0 * pi * ( dfloat(n_sites) / v_uc ) * x_pair_distribution * (y_pair_distribution - 1.d0)
       end if



       !######################################!
       !###---   Calculate the forces   ---###!
       !######################################!

       if ( params%valid_pdf .and. allocated( params%exp_energy_scales ))then

          call get_energy_scale( params%do_md, params%do_mc,&
               & md_istep, params%md_nsteps, mc_istep, params&
               &%mc_nsteps, params &
               &%exp_energy_scales_initial(params%pdf_idx), params &
               &%exp_energy_scales_final(params%pdf_idx), params &
               &%exp_energy_scales(params%pdf_idx) )

          call get_exp_energies(params%exp_energy_scales(params&
               &%pdf_idx), params%exp_data(params%pdf_idx)%y&
               &, y_pair_distribution,&
               & params%pair_distribution_n_samples, n_sites,&
               & energies_pair_distribution(i_beg:i_end))



          if (params%do_forces .and.  params%exp_forces )then
             ! forces_pair_distribution = 0.d0
             ! allocate( pair_distribution_der_temp( 1:params%pair_distribution_n_samples ) )
             ! pair_distribution_der_temp = 0.d0

             n_dim_idx = 1
             outerforces: do j = 1, n_species
                do k = 1, n_species

                   if (j > k) cycle

                   call get_pair_distribution_forces(  n_sites, params%exp_energy_scales(params%pdf_idx),&
                        & params%exp_data(params%pdf_idx)%x, params%exp_data(params%pdf_idx)%y,&
                        & forces_pair_distribution, virial,&
                        & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                        & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                        &%r_range_min, params%r_range_max, params&
                        &%pair_distribution_n_samples,&
                        & y_pair_distribution(1:params&
                        &%pair_distribution_n_samples), params%pair_distribution_rcut&
                        &, j, k, pair_distribution_partial_der(1:params &
                        &%pair_distribution_n_samples, n_dim_idx,&
                        & j_beg:j_end), params%pair_distribution_partial,&
                        & params%pair_distribution_kde_sigma,&
                        & ( (n_atoms_of_species(j) *&
                        & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                        & dfloat(n_sites) ), ( dfloat(n_sites) / v_uc ), params%pair_distribution_output)

                   n_dim_idx = n_dim_idx + 1

                   if ( n_dim_idx > n_dim_partial )then
                      exit outerforces
                   end if

                end do
             end do outerforces
          end if

          ! do i = 1, params%pair_distribution_n_samples
          !    print *,  "ppd ", 100, pair_distribution_der_temp( 1 ), pair_distribution_partial_der( 1, 2, 100:102 )
          ! end do

          ! open(unit=1234, file="grad", status="unknown")
          ! do i = 100, 110
          !    write(1234,  '(A,1X,I8,1X,F20.8)'), "dg_dr_0^1 ", i, pair_distribution_der_temp( i )
          ! end do
          ! close(unit=1234)

          ! deallocate( pair_distribution_der_temp)
       end if

       !##################################################################!
       !###---   If not doing partial pair distribution functions   ---###!
       !##################################################################!


    else
#ifdef _MPIF90
       call mpi_reduce(y_pair_distribution,&
            & y_pair_distribution_temp, params&
            &%pair_distribution_n_samples,&
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
            & MPI_COMM_WORLD, ierr)

       y_pair_distribution =  y_pair_distribution_temp
       deallocate( y_pair_distribution_temp )

       call mpi_bcast(y_pair_distribution, params&
            &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)

       ! if ( do_derivatives .and. params%exp_forces .and. allocated( params%exp_energy_scales ) )then
       !    allocate( forces_pair_distribution_temp(1:3, 1:n_sites) )
       !    forces_pair_distribution_temp = 0.d0

       !    call mpi_reduce(forces_pair_distribution,&
       !         & forces_pair_distribution_temp, 3 * n_sites,&
       !         & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
       !         & MPI_COMM_WORLD, ierr)

       !    forces_pair_distribution = forces_pair_distribution_temp
       !    deallocate( forces_pair_distribution_temp )

       !    call mpi_bcast(forces_pair_distribution, 3*n_sites, MPI_DOUBLE_PRECISION, 0,&
       !         & MPI_COMM_WORLD, ierr)

       ! end if

#endif
       y_pair_distribution =y_pair_distribution * &
            & dot_product( cross_product(a_box, b_box),&
            & c_box ) / (dfloat(indices(1)*indices(2)&
            &*indices(3))) / dfloat(n_sites) / dfloat(n_sites)

    end if


    if (params%pair_distribution_partial .and. allocated(factors)) deallocate(factors)


    ! Write out the partial pair distribution functions
    call get_write_condition( params%do_mc, params%do_md&
         &, mc_istep, md_istep, params%write_xyz,&
         & write_condition)


    if (rank == 0 .and. params%write_pair_distribution .and. write_condition) then
       ! call write_partial_exp(params%do_mc, params%do_md, mc_istep, md_istep,&
       !      & params%write_xyz, params%pair_distribution_partial,&
       !      & n_species, params%pair_distribution_n_samples,&
       !      & n_dim_partial , x_pair_distribution(1:params&
       !      &%pair_distribution_n_samples), y_pair_distribution(1:params&
       !      &%pair_distribution_n_samples), pair_distribution_partial(1:params &
       !      &%pair_distribution_n_samples, 1:n_dim_partial),&
       !      & species_types , 'pair_distribution')

       call get_overwrite_condition( params%do_mc, params%do_md ,&
            & mc_istep, md_istep, params%write_xyz,&
            & overwrite_condition)


       if (params%pair_distribution_partial)then
          n_dim_idx = 1
          outer3: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle

                write(filename,'(A)')&
                     & 'pair_distribution_' // trim(params&
                     &%species_types(j)) // '_' // trim(params&
                     &%species_types(k)) //&
                     & "_prediction.dat"
                call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                     & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                     & overwrite_condition, filename, 'pair_distribution')

                n_dim_idx = n_dim_idx + 1
                if ( n_dim_idx > n_dim_partial )then
                   exit outer3
                end if

             end do
          end do outer3
       end if

       write(filename,'(A)')&
            & "pair_distribution_total.dat"
       call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
            &y_pair_distribution(1:params%pair_distribution_n_samples),&
            & overwrite_condition, filename, "pair_distribution  output: " // trim( params&
            &%pair_distribution_output ))

    end if

  end subroutine gpu_calculate_pair_distribution
  


  
  subroutine calculate_pair_distribution( params, x_pair_distribution&
       &, y_pair_distribution, y_pair_distribution_temp,&
       & pair_distribution_partial, pair_distribution_partial_temp, &
       & n_species, species_types,  n_atoms_of_species, n_sites, a_box, b_box, c_box,&
       & indices, md_istep, mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs, xyz, &
       & neighbors_list, n_neigh, neighbor_species, species, rank,&
       & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
       & pair_distribution_partial_temp_der, energies_pair_distribution, forces_pair_distribution, virial)
    implicit none
    type(input_parameters), intent(inout) :: params
    real*8, allocatable, intent(out) :: x_pair_distribution(:),&
         & y_pair_distribution(:), pair_distribution_partial(:,:),&
         & n_atoms_of_species(:), pair_distribution_partial_temp(:,:),&
         & y_pair_distribution_temp(:), pair_distribution_der(:,:),&
         & pair_distribution_partial_der(:,:,:), &
         & pair_distribution_partial_temp_der(:,:,:), energies_pair_distribution(:), forces_pair_distribution(:,:)
    character*8, allocatable, intent(in) :: species_types(:)
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(inout) :: virial(1:3,1:3)
    real*8 :: v_uc, f
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end
    integer, intent(in) :: indices(1:3), md_istep, mc_istep,  rank
    integer, intent(inout) :: ierr
    real*8, allocatable :: factors(:), pair_distribution_der_temp(:)
    integer :: i, j, k, l, i2, n_dim_partial, n_dim_idx
    logical, intent(in) :: do_derivatives
    real*8, parameter :: pi = acos(-1.0)
    logical :: write_condition, overwrite_condition
    character*1024 :: filename

    ! Things that are allocated here:
    ! Always:
    !  > x_pair_distribution
    !  > y_pair_distribution
    ! if pair_distribution_partial == .true.
    !  > pair_distribution_partial( n_samples, n_spec * (n_spec + 1)/2 )
    !  if do_derivatives == .true.
    !    > pair_distribution_partial_der( n_samples, n_spec * (n_spec + 1)/2, j_beg : j_end )


    ! first allocate the necessary arrays for the
    ! calculation of the pair correlation function
    if (allocated( x_pair_distribution)) deallocate(x_pair_distribution)
    if (allocated( y_pair_distribution)) deallocate(y_pair_distribution)

    allocate( x_pair_distribution( 1: params%pair_distribution_n_samples) )
    allocate( y_pair_distribution( 1: params%pair_distribution_n_samples) )

    if (params%n_exp > 0)then
       do i = 1, params%n_exp
          if ( trim( params%exp_data(i)%label ) == 'pair_distribution' )then
             x_pair_distribution = params%exp_data(i)%x
          end if
       end do
    end if



    if (params%pair_distribution_partial)then
       n_dim_partial = n_species * ( n_species + 1 ) / 2
       allocate(factors( 1:n_dim_partial ))

       n_dim_idx = 1
       outer: do i = 1, n_species
          do j = 1, n_species
             if (i > j) cycle

             if (i /= j)then
                factors(n_dim_idx) = 2.d0
             else
                factors(n_dim_idx) = 1.d0
             end if

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial )then
                exit outer
             end if

          end do
       end do outer


       if (.not. allocated(pair_distribution_partial))then   !deallocate(pair_distribution_partial)
          allocate( pair_distribution_partial(1:params%pair_distribution_n_samples,&
               & 1 : n_dim_partial) )
       end if

       pair_distribution_partial = 0.d0


       if (params%do_forces .and. params%exp_forces)then
          allocate( pair_distribution_partial_der(1:params%pair_distribution_n_samples,&
            & 1 : n_dim_partial, j_beg : j_end  ))
          pair_distribution_partial_der = 0.d0

          if (rank == 0 .and. md_istep == 0) write(*, '(A,1X,F7.4,1X,A)') "Gb/core: partial pdfder = ", dfloat(params&
               &%pair_distribution_n_samples * n_dim_partial * j_end)&
               & * 8.d0 / (dfloat(1024*1024*1024)), " Gb  |"
          if (rank == 0 .and. md_istep == 0) write(*,*)'                                       |'

       end if
    else
       if (params%do_forces .and. params%exp_forces)then
          allocate( pair_distribution_partial_der(1:params%pair_distribution_n_samples, 1:1, &
            &  j_beg : j_end  ))
          pair_distribution_partial_der = 0.d0
       end if


    end if

    if(allocated(n_atoms_of_species)) deallocate(n_atoms_of_species)
    allocate(n_atoms_of_species(1:n_species))

    do j = 1, n_species
       n_atoms_of_species(j) = 0.d0
       do i2 = 1, n_sites
          if ( species(i2) == j)then
             n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
          end if
       end do
    end do


#ifdef _MPIF90
    if (params%pair_distribution_partial)then
       allocate( pair_distribution_partial_temp(1:params%pair_distribution_n_samples, 1 : n_dim_partial) )

       pair_distribution_partial_temp = 0.0d0
    end if

    allocate( y_pair_distribution_temp( 1: params%pair_distribution_n_samples) )
    y_pair_distribution_temp = 0.d0

#endif
    v_uc = dot_product( cross_product(a_box,&
         & b_box), c_box ) / (&
         & dfloat(indices(1)*indices(2)&
         &*indices(3)) )


    !#####################################################################!
    !###---   Calculating the partial pair distribution functions   ---###!
    !#####################################################################!

    if ( params%pair_distribution_partial )then
       n_dim_idx = 1
       outer1: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle ! We have already calculated the pair correlation function!

             ! Note that with the calculation of the derivatives here,
             ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
             ! factor, which allows for some freeing of memory

             call get_pair_distribution( n_sites, &
                  & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                  & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                  &%r_range_min, params%r_range_max, params%pair_distribution_n_samples, x_pair_distribution,&
                  & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), params&
                  &%pair_distribution_rcut, .false.,&
                  & params%pair_distribution_partial, j, k,&
                  & params%pair_distribution_kde_sigma,&
                  & dfloat(n_sites)/v_uc,  params%exp_forces,&
                  & pair_distribution_partial_der, n_dim_idx, &
                  & j_beg, j_end )

             n_dim_idx = n_dim_idx + 1

             if ( n_dim_idx > n_dim_partial )then
                exit outer1
             end if

          end do
       end do outer1


    else
       call get_pair_distribution( n_sites, &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end),&
            & params%r_range_min, params%r_range_max, params &
            &%pair_distribution_n_samples, x_pair_distribution,&
            & y_pair_distribution, params &
            &%pair_distribution_rcut, .false., .false., 1, 1,&
            & params%pair_distribution_kde_sigma, dfloat(n_sites)&
            &/v_uc, params%do_forces .and. params%exp_forces, pair_distribution_partial_der, 1, &
            & j_beg, j_end)
    end if


    ! --- MPI communication is here  ---

    if ( params%pair_distribution_partial )then
#ifdef _MPIF90
       call mpi_reduce(pair_distribution_partial,&
            & pair_distribution_partial_temp, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
            & 0, MPI_COMM_WORLD, ierr)

       ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! Note, we have only so far divided by 4 pi r^2 dr
       ! Therefore, we must scale by the density

       pair_distribution_partial =  pair_distribution_partial_temp
       deallocate( pair_distribution_partial_temp )


       call mpi_bcast(pair_distribution_partial, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)


       ! Now, we have the derivatives of the partial pair distribution
       ! function with respect to the atom pairs in that rank
       !
       ! We can keep them in the rank and calculate forces


       ! call mpi_reduce(pair_distribution_partial_der,&
       !      & pair_distribution_partial_der_temp, params&
       !      &%pair_distribution_n_samples * n_species *&
       !      & n_species * 3 * n_pairs_tot, MPI_DOUBLE_PRECISION, MPI_SUM,&
       !      & 0, MPI_COMM_WORLD, ierr)

       ! ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! ! Note, we have only so far divided by 4 pi r^2 dr
       ! ! Therefore, we must scale by the density

       ! pair_distribution_partial_der =  pair_distribution_partial_der_temp
       ! deallocate( pair_distribution_partial_der_temp )



#endif


       if ( params%valid_pdf )then
          allocate(energies_pair_distribution(1:n_sites))
          energies_pair_distribution = 0.d0

          if (params%do_forces .and. params%exp_forces)then
             allocate(forces_pair_distribution(1:3,1:n_sites))
             forces_pair_distribution = 0.d0
          end if

       end if


       !####################################!
       !###---   Accumulate the PDF   ---###!
       !####################################!


       y_pair_distribution = 0.d0
       n_dim_idx = 1
       outer2: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle

             pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) =&
                  & pair_distribution_partial(1:params&
                  &%pair_distribution_n_samples, n_dim_idx) * v_uc &
                  &  /  n_atoms_of_species(j) /  n_atoms_of_species(k) / factors(n_dim_idx) !real(n_sites)

             if (params%do_forces .and. params%exp_forces) then
                pair_distribution_partial_der(1:params&
                     &%pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) =  pair_distribution_partial_der(1:params&
                     & %pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) * v_uc /  n_atoms_of_species(j) / &
                     & n_atoms_of_species(k) / factors(n_dim_idx)!real(n_sites)

             end if

             y_pair_distribution(1:params%pair_distribution_n_samples) = &
                  & y_pair_distribution(1:params%pair_distribution_n_samples)  +  &
                  &  factors(n_dim_idx) * (n_atoms_of_species(j) * n_atoms_of_species(k)) * &
                  & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) &
                  &  /  dfloat(n_sites) / dfloat(n_sites)

             n_dim_idx = n_dim_idx + 1

             if ( n_dim_idx > n_dim_partial )then
                exit outer2
             end if

          end do
       end do outer2

       ! --- Preprocess the pair distribution according to the output --- !
       if ( trim( params%pair_distribution_output ) == "D(r)" )then
          y_pair_distribution = 4.d0 * pi * ( dfloat(n_sites) / v_uc ) * x_pair_distribution * (y_pair_distribution - 1.d0)
       end if



       !######################################!
       !###---   Calculate the forces   ---###!
       !######################################!

       if ( params%valid_pdf .and. allocated( params%exp_energy_scales ))then

          call get_energy_scale( params%do_md, params%do_mc,&
               & md_istep, params%md_nsteps, mc_istep, params&
               &%mc_nsteps, params &
               &%exp_energy_scales_initial(params%pdf_idx), params &
               &%exp_energy_scales_final(params%pdf_idx), params &
               &%exp_energy_scales(params%pdf_idx) )

          call get_exp_energies(params%exp_energy_scales(params&
               &%pdf_idx), params%exp_data(params%pdf_idx)%y&
               &, y_pair_distribution,&
               & params%pair_distribution_n_samples, n_sites,&
               & energies_pair_distribution(i_beg:i_end))



          if (params%do_forces .and.  params%exp_forces )then
             ! forces_pair_distribution = 0.d0
             ! allocate( pair_distribution_der_temp( 1:params%pair_distribution_n_samples ) )
             ! pair_distribution_der_temp = 0.d0

             n_dim_idx = 1
             outerforces: do j = 1, n_species
                do k = 1, n_species

                   if (j > k) cycle

                   call get_pair_distribution_forces(  n_sites, params%exp_energy_scales(params%pdf_idx),&
                        & params%exp_data(params%pdf_idx)%x, params%exp_data(params%pdf_idx)%y,&
                        & forces_pair_distribution, virial,&
                        & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                        & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                        &%r_range_min, params%r_range_max, params&
                        &%pair_distribution_n_samples,&
                        & y_pair_distribution(1:params&
                        &%pair_distribution_n_samples), params%pair_distribution_rcut&
                        &, j, k, pair_distribution_partial_der(1:params &
                        &%pair_distribution_n_samples, n_dim_idx,&
                        & j_beg:j_end), params%pair_distribution_partial,&
                        & params%pair_distribution_kde_sigma,&
                        & ( (n_atoms_of_species(j) *&
                        & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                        & dfloat(n_sites) ), ( dfloat(n_sites) / v_uc ), params%pair_distribution_output)

                   n_dim_idx = n_dim_idx + 1

                   if ( n_dim_idx > n_dim_partial )then
                      exit outerforces
                   end if

                end do
             end do outerforces
          end if

          ! do i = 1, params%pair_distribution_n_samples
          !    print *,  "ppd ", 100, pair_distribution_der_temp( 1 ), pair_distribution_partial_der( 1, 2, 100:102 )
          ! end do

          ! open(unit=1234, file="grad", status="unknown")
          ! do i = 100, 110
          !    write(1234,  '(A,1X,I8,1X,F20.8)'), "dg_dr_0^1 ", i, pair_distribution_der_temp( i )
          ! end do
          ! close(unit=1234)

          ! deallocate( pair_distribution_der_temp)
       end if

       !##################################################################!
       !###---   If not doing partial pair distribution functions   ---###!
       !##################################################################!


    else
#ifdef _MPIF90
       call mpi_reduce(y_pair_distribution,&
            & y_pair_distribution_temp, params&
            &%pair_distribution_n_samples,&
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
            & MPI_COMM_WORLD, ierr)

       y_pair_distribution =  y_pair_distribution_temp
       deallocate( y_pair_distribution_temp )

       call mpi_bcast(y_pair_distribution, params&
            &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)

       ! if ( do_derivatives .and. params%exp_forces .and. allocated( params%exp_energy_scales ) )then
       !    allocate( forces_pair_distribution_temp(1:3, 1:n_sites) )
       !    forces_pair_distribution_temp = 0.d0

       !    call mpi_reduce(forces_pair_distribution,&
       !         & forces_pair_distribution_temp, 3 * n_sites,&
       !         & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
       !         & MPI_COMM_WORLD, ierr)

       !    forces_pair_distribution = forces_pair_distribution_temp
       !    deallocate( forces_pair_distribution_temp )

       !    call mpi_bcast(forces_pair_distribution, 3*n_sites, MPI_DOUBLE_PRECISION, 0,&
       !         & MPI_COMM_WORLD, ierr)

       ! end if

#endif
       y_pair_distribution =y_pair_distribution * &
            & dot_product( cross_product(a_box, b_box),&
            & c_box ) / (dfloat(indices(1)*indices(2)&
            &*indices(3))) / dfloat(n_sites) / dfloat(n_sites)

    end if


    if (params%pair_distribution_partial .and. allocated(factors)) deallocate(factors)


    ! Write out the partial pair distribution functions
    call get_write_condition( params%do_mc, params%do_md&
         &, mc_istep, md_istep, params%write_xyz,&
         & write_condition)


    if (rank == 0 .and. params%write_pair_distribution .and. write_condition) then
       ! call write_partial_exp(params%do_mc, params%do_md, mc_istep, md_istep,&
       !      & params%write_xyz, params%pair_distribution_partial,&
       !      & n_species, params%pair_distribution_n_samples,&
       !      & n_dim_partial , x_pair_distribution(1:params&
       !      &%pair_distribution_n_samples), y_pair_distribution(1:params&
       !      &%pair_distribution_n_samples), pair_distribution_partial(1:params &
       !      &%pair_distribution_n_samples, 1:n_dim_partial),&
       !      & species_types , 'pair_distribution')

       call get_overwrite_condition( params%do_mc, params%do_md ,&
            & mc_istep, md_istep, params%write_xyz,&
            & overwrite_condition)


       if (params%pair_distribution_partial)then
          n_dim_idx = 1
          outer3: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle

                write(filename,'(A)')&
                     & 'pair_distribution_' // trim(params&
                     &%species_types(j)) // '_' // trim(params&
                     &%species_types(k)) //&
                     & "_prediction.dat"
                call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                     & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                     & overwrite_condition, filename, 'pair_distribution')

                n_dim_idx = n_dim_idx + 1
                if ( n_dim_idx > n_dim_partial )then
                   exit outer3
                end if

             end do
          end do outer3
       end if

       write(filename,'(A)')&
            & "pair_distribution_total.dat"
       call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
            &y_pair_distribution(1:params%pair_distribution_n_samples),&
            & overwrite_condition, filename, "pair_distribution  output: " // trim( params&
            &%pair_distribution_output ))

    end if

  end subroutine calculate_pair_distribution


  subroutine finalize_pair_distribution(params, x_pair_distribution&
       &, y_pair_distribution, y_pair_distribution_temp,&
       & pair_distribution_partial, pair_distribution_partial_temp,&
       & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
       & pair_distribution_partial_temp_der,  n_atoms_of_species, rank)
    implicit none
    type(input_parameters), intent(in) :: params
    integer, intent(in) :: rank
    real*8, allocatable, intent(inout) :: x_pair_distribution(:),&
         & y_pair_distribution(:), pair_distribution_partial(:,:),&
         & n_atoms_of_species(:), pair_distribution_partial_temp(:,:),&
         & y_pair_distribution_temp(:), pair_distribution_der(:,:),&
         & pair_distribution_partial_der(:,:,:), &
         & pair_distribution_partial_temp_der(:,:,:)
    logical, intent(in) :: do_derivatives

    ! Naive finalization, include the logic of how things are actually
    ! allocated above rather then allocating and deallocating

    if ( allocated( x_pair_distribution )               ) deallocate(x_pair_distribution)
    if ( allocated( y_pair_distribution )               ) deallocate(y_pair_distribution)
    if ( allocated( y_pair_distribution_temp )          ) deallocate(y_pair_distribution_temp)
    if ( allocated( pair_distribution_partial )         ) deallocate(pair_distribution_partial)
    if ( allocated( pair_distribution_partial_temp )    ) deallocate(pair_distribution_partial_temp)
    if ( allocated( pair_distribution_der )             ) deallocate(pair_distribution_der)
    if ( allocated( pair_distribution_partial_der )     ) deallocate(pair_distribution_partial_der)
    if ( allocated( pair_distribution_partial_temp_der )) deallocate(pair_distribution_partial_temp_der)
    !    if ( allocated( forces_pair_distribution )          ) deallocate(forces_pair_distribution)
    if ( allocated( n_atoms_of_species )          ) deallocate(n_atoms_of_species)


  end subroutine finalize_pair_distribution



  subroutine calculate_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
       & y_structure_factor, y_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & x_pair_distribution, y_pair_distribution, &
       & pair_distribution_partial, n_species, species_types, n_atoms_of_species,&
       & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep,  i_beg,&
       & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
       & neighbor_species, species, rank , q_beg, q_end, ntasks,&
       & sinc_factor_matrix, do_derivatives, &
       & nk, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
       & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d, &
       pair_distribution_partial_der,&
       & energies_sf, forces_sf, virial_sf, use_matrix_forces, cublas_handle, gpu_stream, gpu_host_storage, gpu_low_memory)
    implicit none
    type(input_parameters), intent(inout) :: params
    real*8, allocatable, intent(out) :: x_structure_factor(:), x_structure_factor_temp(:), &
         & y_structure_factor(:), structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:),&
         & y_structure_factor_temp(:)
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:),&
         & x_pair_distribution(:), y_pair_distribution(:),&
         & pair_distribution_partial(:,:), n_atoms_of_species(:), pair_distribution_partial_der(:,:,:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:), nk(:)
    character*8, allocatable, intent(in) :: species_types(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, allocatable, intent(inout) :: sinc_factor_matrix(:,:), energies_sf(:), forces_sf(:,:)
    real*8, intent(inout) :: virial_sf(1:3,1:3)
    real*8, allocatable :: sinc_factor_matrix_temp(:,:), temp_pdf(:,:)
    real*8 :: v_uc
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end, ntasks
    integer, intent(out) :: q_beg, q_end
    integer, intent(in) :: indices(1:3), md_istep, mc_istep, rank
    integer, intent(out) :: ierr
    integer :: i, j, k, l, i2, n_dim_partial, n_dim_idx, n, m
    real*8 :: dq, f, cabh, delta
    real*8, parameter :: pi = acos(-1.0)
    character*1024 :: filename
    logical :: overwrite_condition, write_condition
    logical, intent(in) :: do_derivatives, use_matrix_forces

    type(c_ptr) :: cublas_handle, gpu_stream
    type(c_ptr), allocatable :: nk_d(:), nk_flags_d(:), nk_flags_sum_d(:), k_index_d(:), j2_index_d(:), rjs_index_d(:), xyz_k_d(:), pair_distribution_partial_d(:), pair_distribution_partial_der_d(:)
    integer(c_size_t), allocatable :: st_nk_d(:), st_k_index_d(:), st_j2_index_d(:), st_rjs(:), st_pair_distribution_partial_d(:), st_pair_distribution_partial_der_d(:)    
    type( gpu_host_storage_type ), intent(inout),  allocatable, target :: gpu_host_storage(:)
    logical, intent(in) :: gpu_low_memory
    
    
    v_uc = dot_product( cross_product(a_box,&
            & b_box), c_box ) / (&
            & dfloat(indices(1)*indices(2)&
            &*indices(3)) )

    if (params%structure_factor_from_pdf) then
       q_beg = 1
       q_end = params%structure_factor_n_samples
#ifdef _MPIF90
       ! We need to do an integral for each q value
       ! This can be split among the processes

       ! Split each q the integrals among each of the
       ! processes, just like we do with the atoms,
       ! and then collect accordingly


       if( rank < mod( params%structure_factor_n_samples, ntasks ) )then
          q_beg = 1 + rank*(params%structure_factor_n_samples / ntasks + 1)
       else
          q_beg = 1 + mod(params%structure_factor_n_samples, ntasks)*(params&
               &%structure_factor_n_samples / ntasks + 1) &
               &+ (rank - mod(params%structure_factor_n_samples, ntasks))*(params&
               &%structure_factor_n_samples / ntasks)
       end if
       if( rank < mod( params%structure_factor_n_samples, ntasks ) )then
          q_end = (rank+1)*(params%structure_factor_n_samples / ntasks + 1)
       else
          q_end = q_beg + params%structure_factor_n_samples/ntasks - 1
       end if


#endif
       if (params%structure_factor_from_pdf .and. params%pair_distribution_partial)then

          n_dim_partial = n_species * ( n_species + 1 ) / 2

          if (allocated(structure_factor_partial))  deallocate(structure_factor_partial)
          allocate( structure_factor_partial(1:params&
               &%structure_factor_n_samples, 1 : n_dim_partial) )
          structure_factor_partial = 0.d0
       end if

    end if

    if (allocated(y_structure_factor))deallocate(y_structure_factor)
    allocate(y_structure_factor(1:params%structure_factor_n_samples))
    y_structure_factor = 0.d0

    if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
#ifdef _MPIF90
       allocate( structure_factor_partial_temp(1:params&
            &%structure_factor_n_samples, 1 : n_dim_partial) )

       structure_factor_partial_temp = 0.0d0
#endif

    else

#ifdef _MPIF90
       allocate( y_structure_factor_temp(1:params&
            &%structure_factor_n_samples) )

       y_structure_factor_temp = 0.0d0
#endif

    end if

    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    call linspace( x_structure_factor, params&
         &%q_range_min, params%q_range_max, params&
         &%structure_factor_n_samples, dq )

    call linspace( x_structure_factor_temp, params&
         &%q_range_min, params%q_range_max, params&
         &%structure_factor_n_samples, dq )

    if( trim(params%q_units) == "xrd" .or. params%q_units == "twotheta")then
       ! assume that theta is given for the Q range
       do i = 1, params%structure_factor_n_samples
          ! This gives s, but Q  = 2 pi * s = 4pi sin(theta) / lambda
          x_structure_factor(i) = 2.d0 * sin( pi *&
               & x_structure_factor(i) /180.d0 / 2.d0 ) /&
               & params%xrd_wavelength
       end do
    elseif (trim(params%q_units) == "saxs" .or. params%q_units == "q")then
       do i = 1, params%structure_factor_n_samples
          x_structure_factor(i) =  x_structure_factor(i)  / 2.d0 / pi
       end do
    end if

    ! Here the units of q range are that called Q (1/A) in
    ! the literature, small q in the literature are =
    ! 2sin(theta)/lambda = Q/2/pi

    if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then

       v_uc = dot_product( cross_product(a_box,&
            & b_box), c_box ) / (&
            & dfloat(indices(1)*indices(2)&
            &*indices(3)) )


       if ( .not. params%structure_factor_matrix ) then
          call get_partial_structure_factor(q_beg, q_end, &
               & pair_distribution_partial,&
               & x_structure_factor(1:params%structure_factor_n_samples) ,&
               & x_pair_distribution, params%pair_distribution_rcut,&
               & params %pair_distribution_n_samples, params&
               &%structure_factor_n_samples, n_species, n_dim_partial,&
               & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, &
               & params%structure_factor_window,&
               & structure_factor_partial)



#ifdef _MPIF90
          call mpi_reduce(structure_factor_partial,&
               & structure_factor_partial_temp, params&
               &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
               & MPI_COMM_WORLD, ierr)
          structure_factor_partial = structure_factor_partial_temp
          deallocate(structure_factor_partial_temp)

          call mpi_bcast(structure_factor_partial, params&
               &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
               & MPI_COMM_WORLD, ierr)

#endif
       else
          if ( .not. allocated(sinc_factor_matrix) )then
             call get_sinc_factor_matrix( q_beg, q_end,&
                  & x_structure_factor, x_pair_distribution,&
                  & params%pair_distribution_rcut, params&
                  &%pair_distribution_n_samples, params&
                  &%structure_factor_n_samples, params%structure_factor_window,&
                  & sinc_factor_matrix )

#ifdef _MPIF90
             allocate( sinc_factor_matrix_temp(  1:params%structure_factor_n_samples, 1:params&
                  &%pair_distribution_n_samples ) )
             sinc_factor_matrix_temp = 0.d0

             call mpi_reduce(sinc_factor_matrix, sinc_factor_matrix_temp&
                  &, params%structure_factor_n_samples * params&
                  &%pair_distribution_n_samples ,&
                  & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
                  & ierr)

             sinc_factor_matrix = sinc_factor_matrix_temp

             deallocate(sinc_factor_matrix_temp)
             call mpi_bcast(sinc_factor_matrix, params&
                  &%structure_factor_n_samples * params&
                  &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
                  & MPI_COMM_WORLD, ierr)

#endif
          end if

          ! Now make a blas call to perform the matrix multiplication necessary to obtain the structure factor.


          n = params%structure_factor_n_samples
          m = params%pair_distribution_n_samples
          k = n_dim_partial

          ! [sinc_factor_matrix] = n_samples_sf * n_samples_pc  ( N x M )
          ! [g_ab]               = n_samples_pc * n_dim_partial ( M x K )
          ! [S_ab]               = n_samples_sf * n_dim_partial ( N x K )
          ! S_ab = [sinc_factor_matrix] x [g_ab - 1]
          !      = ( N x M ) . ( M x K ) -> ( N x K )
          ! alpha * (A * B) + beta * C
          ! A === sinc_factor_partial
          ! B === [g_ab - 1]
          ! C === S_ab
          !   TRANS_A, TRANS_B  N_ROWS_A  N_COLS_B  N_COLS_A ( == N_ROWS_B ) alpha,
          !                                                   A,  first_dim_A,    B,      first_dim_B,  beta,  C, first_dim_C

          call dgemm("N", "N", n, k, m, 1.d0, sinc_factor_matrix, n,&
               & pair_distribution_partial - 1.d0, m, 0.d0,&
               & structure_factor_partial, n)

          !##################################!
          !###---   Do derivatives!    ---###!
          !##################################!

          ! Should do this smartly with memory allocation



          ! All processes do this calculation, so they have the whole of the partial structure factors to work with.
          n_dim_idx = 1
          outersfm: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle

                if (j == k) f = 1.d0
                if (j /= k) f = 0.d0

                cabh = ( (n_atoms_of_species(j) / dfloat( n_sites )) &
                     &*  (n_atoms_of_species(k) / dfloat( n_sites )) &
                     & )**(0.5)

                structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx) = f +&
                     & 4.d0 * pi * cabh * (dfloat(n_sites)/v_uc) * &
                     & structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx)

                n_dim_idx = n_dim_idx + 1
                if (n_dim_idx > n_dim_partial) exit outersfm
             end do
          end do outersfm

       end if


! #ifdef _MPIF90
!           call mpi_reduce(structure_factor_partial,&
!                & structure_factor_partial_temp, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
!                & MPI_COMM_WORLD, ierr)
!           structure_factor_partial = structure_factor_partial_temp
!           deallocate(structure_factor_partial_temp)

!           call mpi_bcast(structure_factor_partial, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
!                & MPI_COMM_WORLD, ierr)

! #endif

       ! using dq as a temp variable
       y_structure_factor = 0.d0
       dq = 0.d0
       n_dim_idx = 1
       outer2: do j = 1, n_species
          dq = dq + (n_atoms_of_species(j) / dfloat(n_sites))
          do k = 1, n_species

             if (j > k) cycle

             if (j == k) f = 1.d0
             if (j /= k) f = 2.d0

             if (j == k) delta = 1.d0
             if (j /= k) delta = 0.d0

             y_structure_factor(1:params%structure_factor_n_samples) = &
                  & y_structure_factor(1:params%structure_factor_n_samples)  +  &
                  & f * ( n_atoms_of_species(j) * n_atoms_of_species(k) )**0.5 * &
                  &  (structure_factor_partial(1:params%structure_factor_n_samples, n_dim_idx) &
                  &  - delta)/ dfloat(n_sites) !/ dfloat(n_sites)

             n_dim_idx = n_dim_idx + 1

             if (n_dim_idx > n_dim_partial) exit outer2

          end do
       end do outer2

       y_structure_factor = y_structure_factor + 1.d0
       ! --- Preprocess the structure factor according to the output --- !
       if ( trim( params%sf_output ) == "q*i(q)" )then
          y_structure_factor = x_structure_factor * (y_structure_factor - 1.d0)
       end if



       if ( params%valid_sf) then

          allocate(energies_sf(1:n_sites))
          energies_sf = 0.d0


          if ( params%valid_sf .and. allocated( params%exp_energy_scales ))then

             call get_energy_scale( params%do_md, params%do_mc,&
                  & md_istep, params%md_nsteps, mc_istep, params&
                  &%mc_nsteps, params%exp_energy_scales_initial(params%sf_idx), &
                  & params%exp_energy_scales_final(params%sf_idx), &
                  & params%exp_energy_scales(params%sf_idx) )

             call get_exp_energies(params%exp_energy_scales(params&
                  &%sf_idx), params%exp_data(params%sf_idx)%y&
                  &, y_structure_factor,&
                  & params%structure_factor_n_samples, n_sites,&
                  & energies_sf(i_beg:i_end))



             if (params%do_forces .and. params%exp_forces)then

                allocate(forces_sf(1:3,1:n_sites))
                forces_sf = 0.d0
                virial_sf = 0.d0                

                n_dim_idx = 1
                outerf: do j = 1, n_species
                   do k = 1, n_species

                      if (j > k) cycle

                      if (j == k) f = 1.d0
                      if (j /= k) f = 2.d0

                      if (use_matrix_forces)then
                         call get_structure_factor_forces_matrix(i_beg, i_end,  n_sites, params%exp_energy_scales(params%sf_idx),&
                              & params%exp_data(params%sf_idx)%x,  params%exp_data(params%sf_idx)%y,&
                              & forces_sf, virial_sf,&
                              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                              & neighbor_species(j_beg:j_end), species_types, species(i_beg:i_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                              &%r_range_min, params%r_range_max, params&
                              &%pair_distribution_n_samples, params&
                              &%structure_factor_n_samples, n_species,&
                              & x_structure_factor(1:params&
                              &%structure_factor_n_samples), y_structure_factor(1:params&
                              &%structure_factor_n_samples), params%pair_distribution_rcut&
                              &, j, k, pair_distribution_partial_der(1:params &
                              &%pair_distribution_n_samples, n_dim_idx,&
                              & j_beg:j_end), params%pair_distribution_partial,&
                              & params%pair_distribution_kde_sigma,&
                              & 4.d0 * pi * f * ( (n_atoms_of_species(j) *&
                              & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                              & dfloat(n_sites) ) * ( dfloat(n_sites) /&
                              & v_uc ), sinc_factor_matrix, n_dim_idx,&
                              & .false., params%xrd_output,&
                              & n_atoms_of_species, .false., cublas_handle, gpu_stream, &
                & nk(n_dim_idx), nk_d(n_dim_idx), k_index_d(n_dim_idx), j2_index_d(n_dim_idx), xyz_k_d(n_dim_idx), pair_distribution_partial_d(n_dim_idx), pair_distribution_partial_der_d(n_dim_idx), &
                & st_nk_d(n_dim_idx), st_k_index_d(n_dim_idx), st_j2_index_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), &
                              gpu_host_storage(n_dim_idx), gpu_low_memory)
                      else
                         call get_structure_factor_forces(  n_sites, params%exp_energy_scales(params%sf_idx),&
                              & params%exp_data(params%sf_idx)%x,  params%exp_data(params%sf_idx)%y,&
                              & forces_sf, virial_sf,&
                              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                              & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                              &%r_range_min, params%r_range_max, params&
                              &%pair_distribution_n_samples, params&
                              &%structure_factor_n_samples, n_species,&
                              & x_structure_factor(1:params&
                              &%structure_factor_n_samples), y_structure_factor(1:params&
                              &%structure_factor_n_samples), params%pair_distribution_rcut&
                              &, j, k, pair_distribution_partial_der(1:params &
                              &%pair_distribution_n_samples, n_dim_idx,&
                              & j_beg:j_end), params%pair_distribution_partial,&
                              & params%pair_distribution_kde_sigma,&
                              & 4.d0 * pi * f * ( (n_atoms_of_species(j) *&
                              & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                              & dfloat(n_sites) ) * ( dfloat(n_sites) /&
                              & v_uc ), sinc_factor_matrix, n_dim_idx,&
                              & .false., params%xrd_output,&
                              & n_atoms_of_species, .false. )
                      end if


                      n_dim_idx = n_dim_idx + 1

                      if (n_dim_idx > n_dim_partial) exit outerf

                   end do
                end do outerf
             end if
          end if
       end if


       ! y_structure_factor(1:params%structure_factor_n_samples)&
       !      & = y_structure_factor(1:params&
       !      &%structure_factor_n_samples) !/ dq


    elseif (params%structure_factor_from_pdf)then


       call get_structure_factor_from_pdf(q_beg, q_end, &
            & y_structure_factor(1:params%structure_factor_n_samples) ,&
            & y_pair_distribution,&
            & x_structure_factor(1:params%structure_factor_n_samples) ,&
            & x_pair_distribution, params%pair_distribution_rcut,&
            & params %pair_distribution_n_samples, params&
            &%structure_factor_n_samples, n_species,&
            & n_atoms_of_species, n_sites, dfloat(n_sites) / v_uc,  params%structure_factor_window)

    else

       call get_structure_factor_explicit( n_sites, &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
            & params%structure_factor_n_samples,&
            & x_structure_factor(1:params&
            &%structure_factor_n_samples), y_structure_factor(1:params&
            &%structure_factor_n_samples),params%r_range_min, params%r_range_max,&
            & params%pair_distribution_rcut, .false., 1&
            &, 1, params%structure_factor_window )
    end if


    if (.not. (params%structure_factor_from_pdf .and. params%pair_distribution_partial))then
#ifdef _MPIF90

       call mpi_reduce(y_structure_factor,&
            & y_structure_factor_temp, params&
            &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
            & MPI_COMM_WORLD, ierr)


       y_structure_factor = y_structure_factor_temp
       deallocate(y_structure_factor_temp)
#endif
       y_structure_factor =  y_structure_factor + 1.d0

#ifdef _MPIF90
       call mpi_bcast(y_structure_factor, params &
            &%structure_factor_n_samples,&
            & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

    end if
    ! Write out the partial structure functions
    call get_write_condition( params%do_mc, params%do_md&
         &, mc_istep, md_istep, params%write_xyz,&
         & write_condition)

    if (rank == 0 .and. params%write_structure_factor .and. write_condition) then

       call get_overwrite_condition( params%do_mc, params%do_md&
            &, mc_istep, md_istep, params%write_xyz,&
            & overwrite_condition)


       if (params%structure_factor_from_pdf .and. params%pair_distribution_partial)then
          n_dim_idx = 1
          outer3: do j = 1, n_species
             do k = 1, n_species

                if ( j > k ) cycle

                ! write with the temp data
                write(filename,'(A)')&
                     & 'structure_factor_' // trim(species_types(j)) // '_' // trim(species_types(k)) //&
                     & "_prediction.dat"
                call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples)&
                     &,&
                     & structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx),&
                     & overwrite_condition, filename, "structure_factor: units of "// trim(params%q_units) )

                n_dim_idx = n_dim_idx + 1
                if ( n_dim_idx > n_dim_partial ) exit outer3

             end do
          end do outer3
       end if

       write(filename,'(A)')&
            & 'structure_factor_total.dat'
       call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples),&
            & y_structure_factor(1:params&
            &%structure_factor_n_samples),&
            & overwrite_condition, filename, "structure_factor: units&
            & of "// trim(params%q_units) // " output: " // trim( params&
            &%xrd_output ))


    end if
  end subroutine calculate_structure_factor


  subroutine finalize_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
       & y_structure_factor, y_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & x_pair_distribution, y_pair_distribution, &
       & pair_distribution_partial, sinc_factor_matrix)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(inout) :: x_structure_factor(:), x_structure_factor_temp(:), &
         & y_structure_factor(:), structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:),&
         & y_structure_factor_temp(:)
    real*8,  intent(inout), allocatable :: x_pair_distribution(:), y_pair_distribution(:),&
         & pair_distribution_partial(:,:)
    real*8, allocatable, intent(inout) :: sinc_factor_matrix(:,:)

    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    if (allocated(y_structure_factor)) deallocate(y_structure_factor)
    if (allocated(structure_factor_partial)) deallocate(structure_factor_partial)
    if (allocated(structure_factor_partial_temp)) deallocate(structure_factor_partial_temp)
    if (allocated(y_structure_factor_temp)) deallocate(y_structure_factor_temp)
    if (allocated(x_pair_distribution)) deallocate(x_pair_distribution)
    if (allocated(y_pair_distribution)) deallocate(y_pair_distribution)
    if (allocated(pair_distribution_partial)) deallocate(pair_distribution_partial)
    if (allocated(sinc_factor_matrix)) deallocate(sinc_factor_matrix)

  end subroutine finalize_structure_factor




  subroutine calculate_xrd( params, x_xrd, x_xrd_temp,&
       & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & n_species, species_types,  n_atoms_of_species,&
       & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
       & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
       & neighbor_species, species, rank , q_beg, q_end, ntasks,&
       & sinc_factor_matrix, do_derivatives, &
       & nk, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
       & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d, &
        pair_distribution_partial_der,&
       & energies_xrd, forces_xrd, virial_xrd, neutron, use_matrix_forces, cublas_handle, gpu_stream, gpu_host_storage, gpu_low_memory )
    implicit none
    type(input_parameters), intent(inout) :: params
    real*8, allocatable, intent(out) :: x_xrd(:), x_xrd_temp(:), &
         & y_xrd(:), y_xrd_temp(:), energies_xrd(:), forces_xrd(:,:)
    real*8, allocatable, intent(in) :: structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:), x_structure_factor(:)&
         &, x_structure_factor_temp(:), sinc_factor_matrix(:,:),&
         & pair_distribution_partial_der(:,:,:)
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:),&
         & n_atoms_of_species(:)
    real*8, intent(out) :: virial_xrd(1:3,1:3)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:), nk(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    character*8, allocatable, intent(in) :: species_types(:)
    real*8 :: v_uc
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end, ntasks
    integer, intent(out) :: q_beg, q_end
    integer, intent(in) :: indices(1:3), md_istep, mc_istep, rank
    integer, intent(out) :: ierr
    real*8, allocatable :: y_sub(:)
    integer :: i, j, k, l, i2, n_dim_idx, n_dim_partial
    real*8 :: dq, f
    real*8, parameter :: pi = acos(-1.0)
    character*1024 :: filename
    logical :: write_condition, overwrite_condition, valid_xrd
    logical, intent(in) :: do_derivatives, neutron, use_matrix_forces
    integer :: xrd_idx
    character*32 :: xrd_output

    type(c_ptr) :: cublas_handle, gpu_stream
    type(c_ptr), allocatable :: nk_d(:), nk_flags_d(:), nk_flags_sum_d(:), k_index_d(:), j2_index_d(:), rjs_index_d(:), xyz_k_d(:), pair_distribution_partial_d(:), pair_distribution_partial_der_d(:)
    integer(c_size_t), allocatable :: st_nk_d(:), st_k_index_d(:), st_j2_index_d(:), st_rjs(:), st_pair_distribution_partial_d(:), st_pair_distribution_partial_der_d(:)    
    type( gpu_host_storage_type ), intent(inout),  allocatable, target :: gpu_host_storage(:)
    logical, intent(in) :: gpu_low_memory

    
    if ( neutron ) xrd_idx = params%nd_idx
    if ( .not. neutron ) xrd_idx = params%xrd_idx

    if ( neutron ) xrd_output = params%nd_output
    if ( .not. neutron ) xrd_output = params%xrd_output

    if ( neutron) valid_xrd = params%valid_nd
    if (.not. neutron) valid_xrd = params%valid_xrd


    v_uc = dot_product( cross_product(a_box,&
         & b_box), c_box ) / (&
         & dfloat(indices(1)*indices(2)&
         &*indices(3)) )

    n_dim_partial = n_species * (n_species + 1 ) / 2

    ! Get the XRD from the partial structure factors!
    if (allocated(x_xrd)) deallocate(x_xrd)
    allocate(x_xrd(1:params%structure_factor_n_samples))
    if (allocated(x_xrd_temp)) deallocate(x_xrd_temp)
    allocate(x_xrd_temp(1:params%structure_factor_n_samples))

    x_xrd = x_structure_factor
    x_xrd_temp = x_structure_factor_temp


    if (allocated(y_xrd)) deallocate(y_xrd)
    allocate(y_xrd(1:params%structure_factor_n_samples))
    y_xrd = 0.d0

    ! allocate for the structure factor parameters, so we don't have to look through the horrible list



#ifdef _MPIF90
    allocate( y_xrd_temp(1:params&
         &%structure_factor_n_samples ))

    y_xrd_temp = 0.0d0
#endif

    ! Have the same range for the xrd as the structure factor

    allocate(y_sub(1:params&
         &%structure_factor_n_samples))

    if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
       call get_xrd_from_partial_structure_factors(q_beg, q_end, &
            & structure_factor_partial(1:params&
            &%structure_factor_n_samples,1:n_dim_partial), n_species, species_types&
            &, species, params%xrd_wavelength, params &
            &%xrd_damping, params%xrd_alpha, params &
            &%xrd_method, params%xrd_iwasa, xrd_output, x_xrd(1:params&
            &%structure_factor_n_samples), y_xrd(1:params&
            &%structure_factor_n_samples),&
            & n_atoms_of_species, y_sub, neutron)

    else

       call get_xrd_explicit( n_sites, species_types, n_species,  &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
            & params%structure_factor_n_samples,&
            & x_xrd(1:params&
            &%structure_factor_n_samples), y_xrd(1:params&
            &%structure_factor_n_samples),params%r_range_min, params%r_range_max,&
            & params%pair_distribution_rcut, .false., 1&
            &, 1, params%structure_factor_window )


    end if

    !###################################################################################!
    !###---   Can calculate the Structure factors related to XRD / Neutron here   ---###!
    !###################################################################################!


    ! if (allocated(sf_parameters) ) deallocate( sf_parameters )
    ! allocate( sf_parameters(1:9,1:n_species))
    ! sf_parameters = 0.d0
    ! do i = 1, n_species
    !    call get_scattering_factor_params(params%species_types(i), sf_parameters(1:9,i))
    ! end do



    ! do j = 1, params%structure_factor_n_samples
    !    wfac = 0.d0
    !    do k = 1, n_species
    !       call get_scattering_factor(wfac_temp, sf_parameters(1:9,k), x_xrd(j)/2.d0 )
    !       wfac = wfac + wfac_temp*wfac_temp *n_atoms_of_species(k)  / dfloat(n_sites)
    !    end do
    !    print *, wfac
    !    y_xrd(j) = y_xrd(j) / (wfac)
    ! end do



#ifdef _MPIF90

    call mpi_reduce(y_xrd,&
         & y_xrd_temp, params&
         &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & MPI_COMM_WORLD, ierr)
    y_xrd = y_xrd_temp
    deallocate(y_xrd_temp)

    call mpi_bcast(y_xrd, params&
         &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, 0,&
         & MPI_COMM_WORLD, ierr)

#endif

    ! Already preprocessed the xrd !
    ! --- Preprocess the structure factor according to the output --- !
    ! if ( trim( params%xrd_output ) == "q*i(q)" )then
    !    y_xrd = 2.d0 * pi * x_xrd * ( y_xrd )
    ! end if


    if ( allocated(sinc_factor_matrix) )then
       if (valid_xrd) then

          allocate(energies_xrd(1:n_sites))
          energies_xrd = 0.d0

          if ( valid_xrd .and. allocated( params%exp_energy_scales ))then
             call get_energy_scale( params%do_md, params%do_mc,&
                  & md_istep, params%md_nsteps, mc_istep, params &
                  &%mc_nsteps, params&
                  &%exp_energy_scales_initial(xrd_idx), params&
                  &%exp_energy_scales_final(xrd_idx), params&
                  &%exp_energy_scales(xrd_idx) )

             call get_exp_energies(params%exp_energy_scales(xrd_idx), params%exp_data(xrd_idx)%y&
                  &, y_xrd,&
                  & params%structure_factor_n_samples, n_sites,&
                  & energies_xrd(i_beg:i_end))



             if (params%do_forces .and. params%exp_forces)then

                allocate(forces_xrd(1:3,1:n_sites))
                forces_xrd = 0.d0
                virial_xrd = 0.d0                
                n_dim_idx = 1
                outerf: do j = 1, n_species
                   do k = 1, n_species

                      if (j > k) cycle

                      if (j == k) f = 1.d0
                      if (j /= k) f = 2.d0

                      if (use_matrix_forces) then

                         call get_structure_factor_forces_matrix( i_beg, i_end, n_sites, params%exp_energy_scales(xrd_idx),&
                              & x_xrd, params%exp_data(xrd_idx)%y,&
                              & forces_xrd, virial_xrd,&
                              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                              & neighbor_species(j_beg:j_end), species_types, species(i_beg:i_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                              &%r_range_min, params%r_range_max, params&
                              &%pair_distribution_n_samples, params&
                              &%structure_factor_n_samples, n_species,&
                              & x_xrd(1:params&
                              &%structure_factor_n_samples), y_xrd(1:params&
                              &%structure_factor_n_samples), params%pair_distribution_rcut&
                              &, j, k, pair_distribution_partial_der(1:params &
                              &%pair_distribution_n_samples, n_dim_idx,&
                              & j_beg:j_end), params%pair_distribution_partial,&
                              & params%pair_distribution_kde_sigma,&
                              & 4.d0 * pi * f * ( (n_atoms_of_species(j) *&
                              & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                              & dfloat(n_sites) ) * ( dfloat(n_sites) /&
                              & v_uc ), sinc_factor_matrix, n_dim_idx,&
                              & .true., xrd_output, n_atoms_of_species, neutron, cublas_handle, gpu_stream,&
                & nk(n_dim_idx), nk_d(n_dim_idx), k_index_d(n_dim_idx), j2_index_d(n_dim_idx), xyz_k_d(n_dim_idx), pair_distribution_partial_d(n_dim_idx), pair_distribution_partial_der_d(n_dim_idx), &
                & st_nk_d(n_dim_idx), st_k_index_d(n_dim_idx), st_j2_index_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), &
                              gpu_host_storage(n_dim_idx), gpu_low_memory)
                      else
                         call get_structure_factor_forces(  n_sites, params%exp_energy_scales(xrd_idx),&
                              & x_xrd, params%exp_data(xrd_idx)%y,&
                              & forces_xrd, virial_xrd,&
                              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                              & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                              &%r_range_min, params%r_range_max, params&
                              &%pair_distribution_n_samples, params&
                              &%structure_factor_n_samples, n_species,&
                              & x_xrd(1:params&
                              &%structure_factor_n_samples), y_xrd(1:params&
                              &%structure_factor_n_samples), params%pair_distribution_rcut&
                              &, j, k, pair_distribution_partial_der(1:params &
                              &%pair_distribution_n_samples, n_dim_idx,&
                              & j_beg:j_end), params%pair_distribution_partial,&
                              & params%pair_distribution_kde_sigma,&
                              & 4.d0 * pi * f * ( (n_atoms_of_species(j) *&
                              & n_atoms_of_species(k)) /  dfloat(n_sites) /&
                              & dfloat(n_sites) ) * ( dfloat(n_sites) /&
                              & v_uc ), sinc_factor_matrix, n_dim_idx,&
                              & .true., xrd_output, n_atoms_of_species, neutron)
                      end if


                      n_dim_idx = n_dim_idx + 1

                      if (n_dim_idx > n_dim_partial) exit outerf

                   end do
                end do outerf
             end if
          end if
       end if
    end if




    call get_write_condition( params%do_mc, params%do_md&
         &, mc_istep, md_istep, params%write_xyz,&
         & write_condition)

    if (rank == 0 .and. ( params%write_xrd .or.  params%write_nd)  .and. write_condition) then

       call get_overwrite_condition( params%do_mc, params%do_md&
         &, mc_istep, md_istep, params%write_xyz,&
         & overwrite_condition)


       if (.not. neutron) then
          write(filename,'(A)')&
               & 'xrd_prediction.dat'
          call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
               & y_xrd(1:params&
               &%structure_factor_n_samples),&
               & overwrite_condition, filename, "xrd: units of "//&
               & trim(params%q_units) // " output: " // trim( xrd_output ))
       else
          write(filename,'(A)')&
               & 'nd_prediction.dat'
          call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
               & y_xrd(1:params&
               &%structure_factor_n_samples),&
               & overwrite_condition, filename, "nd: units of "//&
               & trim(params%q_units) // " output: " // trim( xrd_output ))
       end if

    end if

    ! if ( trim(params%xrd_output) == "q*i(q)" .or. trim(params%xrd_output) == "q*F(q)")then
    !    ! output q * i(q) === q * F_x(q)
    !    do l = q_beg, q_end
    !       y_xrd(l) = x_xrd_temp(l) * ( y_xrd(l) - y_sub(l) )
    !    end do

    !    elseif( trim(output) == "F(q)" .or. trim(output) == "i(q)")
    !       ! do nothing,
    !       ! Output the total scattering functon, i(q) === F_x(q)

    !    end if
    if (allocated(y_sub)) deallocate(y_sub)



  end subroutine calculate_xrd


  subroutine finalize_xrd( params, x_xrd, x_xrd_temp,&
       & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(inout) :: x_xrd(:), x_xrd_temp(:), &
         & y_xrd(:), y_xrd_temp(:)
    real*8, allocatable, intent(inout) :: structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:), x_structure_factor(:), x_structure_factor_temp(:)

    if (allocated(x_xrd)) deallocate(x_xrd)
    if (allocated(x_xrd_temp)) deallocate(x_xrd_temp)
    if (allocated(y_xrd)) deallocate(y_xrd)
    if (allocated(y_xrd_temp)) deallocate(y_xrd_temp)
    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    if (allocated(structure_factor_partial)) deallocate(structure_factor_partial)
    if (allocated(structure_factor_partial_temp)) deallocate(structure_factor_partial_temp)
  end subroutine finalize_xrd





  

end module exp_interface
