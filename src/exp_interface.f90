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
contains

  ! This module implements the interfaces for the gradient of experimental functions

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


  subroutine get_pdf_sf_xrd_explicitly_kde( v_uc, n_sites0, n_species, species, species_types, &
       & neighbors_list, n_neigh, neighbor_species, rjs, xyz, r_cut, &
       & r_min, r_max, n_samples,  &!       & q_min, q_max, n_samples_sf, x_sf,  &
       & pair_distribution, pair_distribution_der, kde_sigma, do_derivatives, rank )
       ! & do_forces_pdf, do_forces_sf, do_forces_xrd, &
       ! &    forces_pdf,    forces_sf,    forces_xrd )
    implicit none
    ! Input Variables
    real*8,  intent(in) :: rjs(:), kde_sigma, xyz(:,:), v_uc,  r_min, r_max, r_cut
    integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:), species(:)
    integer, intent(in) :: n_sites0, n_samples, n_species, rank
    logical, intent(in) :: do_derivatives !, do_forces_pdf, do_forces_sf, do_forces_xrd
    character*8, intent(in), allocatable :: species_types(:)
    ! Output Variables
    real*8,  intent(out), allocatable :: pair_distribution(:,:)
    real*8,  intent(out), allocatable :: pair_distribution_der(:,:,:)!, forces_pdf(:,:)

    ! Internal  Variables
    integer :: n_sites, n_pairs, s, species_i, species_j
    integer :: i, j, k, l, i2, j2, n_dim_partial, n_dim_idx, ierr
    real*8  :: r, gauss, c
    real*8, allocatable :: bin_edges(:), dV(:), factors(:), pair_distribution_partial_temp(:,:),&
         & pdf(:), sf(:), xrd(:), forces_sf(:,:), forces_xrd(:,:), n_atoms_of_species(:), &
         prefactor_pdf(:), prefactor_sf(:), prefactor_xrd(:), x_pdf(:)

    integer, allocatable :: species_to_ndim(:,:), species1_partial(:), species2_partial(:)
    character*1024 :: filename
    ! Parameters
    real*8, parameter :: pi = acos(-1.0)



    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)

    allocate(n_atoms_of_species(1:n_species))
    n_atoms_of_species = 0.d0
    do i = 1, size( species, 1 )
       s = species(i)
       n_atoms_of_species( s ) = n_atoms_of_species( s ) + 1
    end do


    allocate(bin_edges(1:n_samples+1))
    allocate(dV(1:n_samples))
    allocate(x_pdf(1:n_samples))


    n_dim_partial = n_species * ( n_species + 1 ) / 2
    allocate(factors( 1:n_dim_partial ))
    allocate(species_to_ndim( 1:n_species, 1:n_species ))

    allocate(species1_partial( 1:n_dim_partial ))
    allocate(species2_partial( 1:n_dim_partial ))

    allocate( pair_distribution(1:n_samples, 1:n_dim_partial) )
    pair_distribution = 0.d0

    if (do_derivatives)then
       allocate( pair_distribution_der(1:n_pairs, 1:n_samples, 1:n_dim_partial) )
       pair_distribution_der = 0.d0
    end if


    n_dim_idx = 1
    outer: do i = 1, n_species
       do j = 1, n_species
          if (i > j) cycle

          if (i /= j)then
             factors(n_dim_idx) = 2.d0
          else
             factors(n_dim_idx) = 1.d0
          end if

          species_to_ndim(i,j) = n_dim_idx
          species_to_ndim(j,i) = n_dim_idx

          species1_partial(n_dim_idx) = i
          species2_partial(n_dim_idx) = j

          n_dim_idx = n_dim_idx + 1

          if ( n_dim_idx > n_dim_partial ) exit outer

       end do
    end do outer


    x_pdf = 0.d0
    bin_edges = 0.d0
    pair_distribution = 0.d0


    do i = 1, n_samples + 1
       bin_edges(i) = r_min  +  ( real( (i-1) ) /  real(n_samples) ) * (r_max - r_min)
    end do

    do i = 1, n_samples
       x_pdf(i) = (bin_edges(i) + bin_edges(i+1)) / 2.d0
       dV(i) = 4.d0 * pi * (bin_edges(i)**2 * (bin_edges(i+1) - bin_edges(i)) )
    end do


    k = 0
    do i = 1, n_sites
       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
       species_i = neighbor_species(k+1)
       k = k + 1
       do j = 2, n_neigh(i)
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
          species_j = neighbor_species(k)
          r = rjs(k) ! atom pair distance

          if ( r > r_cut .or. r < r_min .or. r > r_max + kde_sigma*6.d0  ) cycle

          n_dim_idx = species_to_ndim(species_i, species_j)

          do l = 1,n_samples
             gauss = exp( -( (x_pdf(l) - r) / kde_sigma )**2 / 2.d0 )
             pair_distribution(l, n_dim_idx) = pair_distribution(l, n_dim_idx) + gauss

             ! Construct the initial derivatives here without the ri-rj term here
             if (do_derivatives)then
                pair_distribution_der(k, l, n_dim_idx) = gauss * ( (x_pdf(l) - r) / kde_sigma**2 ) / r
             end if


          end do
       end do
    end do


    pair_distribution = pair_distribution * ( ( r_max - r_min) / dfloat(n_samples) ) &
         & / ( sqrt( 2.d0 * pi) * kde_sigma)

    if (do_derivatives) pair_distribution_der = pair_distribution_der * ( ( r_max - r_min) / dfloat(n_samples) ) &
         & / ( sqrt( 2.d0 * pi) * kde_sigma)


    do i = 1, n_dim_partial
       pair_distribution( 1:n_samples, i) =  pair_distribution(&
            & 1:n_samples, i) * v_uc / n_atoms_of_species(&
            & species1_partial(i) ) / n_atoms_of_species(&
            & species2_partial(i) ) / factors(i)

       pair_distribution( 1:n_samples, i) =  pair_distribution( 1:n_samples, i) / dV

       if (do_derivatives)then
          pair_distribution_der(1:n_pairs,  1:n_samples, i) =  pair_distribution_der(&
               & 1:n_pairs, 1:n_samples, i) * v_uc / n_atoms_of_species(&
               & species1_partial(i) ) / n_atoms_of_species(&
               & species2_partial(i) ) / factors(i)

          do l = 1, n_pairs
             pair_distribution_der(l, 1:n_samples, i) =  pair_distribution_der( l, 1:n_samples, i) / dV
          end do
       end if

    end do

#ifdef _MPIF90
    allocate(pair_distribution_partial_temp( 1:n_samples, 1:n_dim_partial))
    pair_distribution_partial_temp = 0.d0

    call mpi_reduce(pair_distribution,&
         & pair_distribution_partial_temp, n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
         & 0, MPI_COMM_WORLD, ierr)

    pair_distribution =  pair_distribution_partial_temp
    deallocate( pair_distribution_partial_temp )

    call mpi_bcast(pair_distribution, n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
         & MPI_COMM_WORLD, ierr)
#endif

    allocate(pdf(1:n_samples))
    pdf = 0.d0

    do i = 1, n_dim_partial
       c = (n_atoms_of_species( species1_partial(i) ) *&
            & n_atoms_of_species( species1_partial(i) ) ) /&
            & dfloat(n_sites0) / dfloat(n_sites0)

       pdf(1:n_samples) = pdf(1:n_samples) + c * factors(i) * pair_distribution( 1:n_samples, i)
    end do


    ! Write out the partial pair distribution functions
    if (rank == 0 ) then
       n_dim_idx = 1
       outer3: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle

             write(filename,'(A)')&
                  & 'tpair_distribution_' // trim( &
                  & species_types(j)) // '_' // trim( &
                  & species_types(k)) //&
                  & "_prediction.dat"
             call write_exp_datan(x_pdf(1:n_samples),&
                  & pair_distribution(1:n_samples, n_dim_idx),&
                  & .true., filename, 'pair_distribution')

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial )then
                exit outer3
             end if

          end do
       end do outer3

       write(filename,'(A)')&
            & "tpair_distribution_total.dat"
       call write_exp_datan(x_pdf(1:n_samples),&
            &pdf(1:n_samples),&
            & .true., filename, "pair_distribution")

    end if



    ! Now we have the pdf, we can construct the derivatives for the
    ! forces pdf only and we can construct the structure factor / xrd
    ! spectrum

    ! if (do_derivatives .and. ( do_forces_pdf .or. do_forces_sf .or. do_forces_xrd ))then

    !    if (do_forces_pdf) allocate(prefactor_pdf(1:n_samples))
    !    if (do_forces_sf ) allocate(prefactor_sf(1:n_samples_sf))
    !    if (do_forces_xrd) allocate(prefactor_xrd(1:n_samples_sf))

    !    k = 0
    !    do i = 1, n_sites
    !       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
    !       species_i = neighbor_species(k+1)
    !       do j = 1, n_neigh(i)
    !          k = k + 1
    !          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
    !          species_j = neighbor_species(k)
    !          r = rjs(k) ! atom pair distance

    !          if ( r > r_cut .or. r < r_min .or. r > r_max + kde_sigma*6.d0  ) cycle

    !          n_dim_idx = species_to_ndim(species_i, species_j)

    !          if ( .not. all( xyz( 1:3, k ) == 0.d0 ) )then
    !             ! Actual derivative of the pair distribution function given here, no normalisation needed I think

    !             if ( do_forces_pdf )then
    !                call get_this_exp_force(k, xyz(1:3,k), n_samples, n_dim_idx, pair_distribution_der, energy_scale, f,  this_force)
    !                forces_pdf(1:3, j2) = forces_pdf(1:3, j2) + this_force(1:3)
    !             end if

    !             if ( do_forces_sf )then
    !                call get_this_exp_force(k, xyz(1:3,k), n_samples_sf, n_dim_idx, pair_distribution_der, energy_scale, f,  this_force)
    !                forces_sf(1:3, j2) = forces_sf(1:3, j2) + this_force(1:3)
    !             end if

    !             if ( do_forces_xrd )then
    !                call get_this_exp_force(k, xyz(1:3,k), n_samples_sf, n_dim_idx, pair_distribution_der, energy_scale, f,  this_force)
    !                forces_xrd(1:3, j2) = forces_xrd(1:3, j2) + this_force(1:3)
    !             end if

    !          end if
    !       end do
    !    end do
    ! end if


    deallocate(species_to_ndim)
    deallocate(species1_partial)
    deallocate(species2_partial)
    deallocate(factors)
    deallocate(dV)
    deallocate(x_pdf)
    deallocate(bin_edges)
    deallocate(n_atoms_of_species)


  end subroutine get_pdf_sf_xrd_explicitly_kde


  ! subroutine get_this_exp_force(k, xyz, n_samples, n_dim_idx, pair_distribution_der, energy_scale, f, this_force)
  !   implicit none
  !   integer, intent(in) :: k, n_dim_idx, n_samples
  !   real, intent(in), allocatable  :: prefactor(:), pair_distribution_der(:,:,:)
  !   real, intent(in) :: rij(1:3), f, energy_scale
  !   real, intent(out) :: this_force(1:3)

  !   this_force(1) = dot_product( - 2.d0 * rij( 1 ) *&
  !        & pair_distribution_der(k, 1:n_samples,  n_dim_idx),&
  !        & prefactor(1:n_samples))

  !   this_force(2) = dot_product( - 2.d0 * rij( 2 ) *&
  !        & pair_distribution_der(k, 1:n_samples, n_dim_idx),&
  !        & prefactor(1:n_samples))

  !   this_force(3) = dot_product( - 2.d0 * rij( 3 ) *&
  !        & pair_distribution_der(k, 1:n_samples,  n_dim_idx),&
  !        & prefactor(1:n_samples))

  !   this_force(1:3) =  - f * energy_scale  * this_force(1:3)

  ! end subroutine get_this_exp_force


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
       & sinc_factor_matrix, do_derivatives, pair_distribution_partial_der,&
       & energies_sf, forces_sf, virial_sf, use_matrix_forces)
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
         &, neighbor_species(:), species(:)
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

                n_dim_idx = 1
                outerf: do j = 1, n_species
                   do k = 1, n_species

                      if (j > k) cycle

                      if (j == k) f = 1.d0
                      if (j /= k) f = 2.d0

                      if (use_matrix_forces)then
                         call get_structure_factor_forces_matrix(  n_sites, params%exp_energy_scales(params%sf_idx),&
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
       & sinc_factor_matrix, do_derivatives, pair_distribution_partial_der,&
       & energies_xrd, forces_xrd, virial_xrd, neutron, use_matrix_forces )
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
         &, neighbor_species(:), species(:)
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
                  & md_istep, params%md_nsteps, mc_istep, params&
                  &%mc_nsteps, params%exp_energy_scales_initial(xrd_idx), &
                  & params%exp_energy_scales_final(xrd_idx), &
                  & params%exp_energy_scales(xrd_idx) )

             call get_exp_energies(params%exp_energy_scales(xrd_idx), params%exp_data(xrd_idx)%y&
                  &, y_xrd,&
                  & params%structure_factor_n_samples, n_sites,&
                  & energies_xrd(i_beg:i_end))



             if (params%do_forces .and. params%exp_forces)then

                allocate(forces_xrd(1:3,1:n_sites))
                forces_xrd = 0.d0

                n_dim_idx = 1
                outerf: do j = 1, n_species
                   do k = 1, n_species

                      if (j > k) cycle

                      if (j == k) f = 1.d0
                      if (j /= k) f = 2.d0

                      if (use_matrix_forces) then

                         call get_structure_factor_forces_matrix(  n_sites, params%exp_energy_scales(xrd_idx),&
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
