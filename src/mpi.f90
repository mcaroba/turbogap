! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mpi.f90, is copyright (c) 2019-2022, Miguel A. Caro
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

module mpi_helper


  use types


  contains


  subroutine allocate_soap_turbo_hypers(n_soap_turbo, n_species, n_sparse, dim,&! compress_P_nonzero, &
       local_properties_n_sparse, local_properties_dim, has_local_properties, n_local_properties, &
       compress_soap, desc)

!   Input variables
    integer, intent(in) :: n_soap_turbo, n_species(:), n_sparse(:), dim(:), local_properties_n_sparse(:), &
                            n_local_properties(:), local_properties_dim(:) ! compress_P_nonzero(:),
    logical, intent(in) :: compress_soap(:), has_local_properties(:)

!   Output_variables
    type(soap_turbo), allocatable, intent(out) :: desc(:)

!   Internal variables

    integer :: i, n_sp, d, cPnz, counter

    allocate( desc(1:n_soap_turbo) )
    counter = 0
    do i = 1, n_soap_turbo
       n_sp = n_species(i)
       desc(i)%n_species = n_sp
       allocate( desc(i)%nf(1:n_sp) )
       allocate( desc(i)%rcut_hard(1:n_sp) )
       allocate( desc(i)%rcut_soft(1:n_sp) )
       allocate( desc(i)%atom_sigma_r(1:n_sp) )
       allocate( desc(i)%atom_sigma_t(1:n_sp) )
       allocate( desc(i)%atom_sigma_r_scaling(1:n_sp) )
       allocate( desc(i)%atom_sigma_t_scaling(1:n_sp) )
       allocate( desc(i)%amplitude_scaling(1:n_sp) )
       allocate( desc(i)%central_weight(1:n_sp) )
       allocate( desc(i)%global_scaling(1:n_sp) )
       allocate( desc(i)%alpha_max(1:n_sp) )
       allocate( desc(i)%species_types(1:n_sp) )
       n_sp = n_sparse(i)
       desc(i)%n_sparse = n_sp
       d = dim(i)
       desc(i)%dim = d
!       cPnz = compress_P_nonzero(i)
!       desc(i)%compress_P_nonzero = cPnz
       allocate( desc(i)%alphas(1:n_sp) )
       allocate( desc(i)%Qs(1:d, 1:n_sp) )
       !     Currently the cutoff does not get allocated
       !      allocate( desc(i)%cutoff(1:n_sp) )
      if( compress_soap(i) )then
        allocate( desc(i)%compress_soap_indices(1:d) )
      end if
       
       ! if( compress_soap(i) )then
       !    allocate( desc(i)%compress_P_el(1:cPnz) )
       !    allocate( desc(i)%compress_P_i(1:cPnz) )
       !    allocate( desc(i)%compress_P_j(1:cPnz) )
       ! end if
       if( has_local_properties(i) )then
          desc(i)%n_local_properties = n_local_properties(i)
          allocate( desc(i)%local_property_models(1:n_local_properties(i)))
          ! actually this isn't that simple as one has to broadcast and allocate the local_property models
          if ( desc(i)%n_local_properties > 0) then
             do j = 1, desc(i)%n_local_properties
                counter = counter + 1
                n_sp = local_properties_n_sparse(counter)
                d = local_properties_dim(counter)
                desc(i)%local_property_models(j)%n_sparse = n_sp
                desc(i)%local_property_models(j)%dim = d
                allocate( desc(i)%local_property_models(j)%alphas(1:n_sp) )
                ! Assuming same dimension as from vdw implementation
                ! will just change d to the actual one read in
                allocate( desc(i)%local_property_models(j)%Qs(1:d, 1:n_sp) )
             end do
          end if

       end if
    end do

  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine allocate_distance_2b_hypers(n_distance_2b, n_sparse, desc)

!   Input variables
    integer, intent(in) :: n_distance_2b, n_sparse(:)

!   Output_variables
    type(distance_2b), allocatable, intent(out) :: desc(:)

!   Internal variables
    integer :: i, n_sp

    allocate( desc(1:n_distance_2b) )

    do i = 1, n_distance_2b
      n_sp = n_sparse(i)
      desc(i)%n_sparse = n_sp
      allocate( desc(i)%alphas(1:n_sp) )
      allocate( desc(i)%Qs(1:n_sp, 1:1) )
      allocate( desc(i)%cutoff(1:n_sp) )
    end do

  end subroutine
!**************************************************************************





!**************************************************************************
  subroutine allocate_angle_3b_hypers(n_angle_3b, n_sparse, desc)

!   Input variables
    integer, intent(in) :: n_angle_3b, n_sparse(:)

!   Output_variables
    type(angle_3b), allocatable, intent(out) :: desc(:)

!   Internal variables
    integer :: i, n_sp

    allocate( desc(1:n_angle_3b) )

    do i = 1, n_angle_3b
      n_sp = n_sparse(i)
      desc(i)%n_sparse = n_sp
      allocate( desc(i)%alphas(1:n_sp) )
      allocate( desc(i)%Qs(1:n_sp, 1:3) )
      allocate( desc(i)%cutoff(1:n_sp) )
    end do

  end subroutine
!**************************************************************************






!**************************************************************************
  subroutine allocate_core_pot_hypers(n_core_pot, n, desc)

!   Input variables
    integer, intent(in) :: n_core_pot, n(:)

!   Output_variables
    type(core_pot), allocatable, intent(out) :: desc(:)

!   Internal variables
    integer :: i, n_sp

    allocate( desc(1:n_core_pot) )

    do i = 1, n_core_pot
      n_sp = n(i)
      desc(i)%n = n_sp
      allocate( desc(i)%x(1:n_sp) )
      allocate( desc(i)%V(1:n_sp) )
      allocate( desc(i)%dVdx2(1:n_sp) )
    end do

  end subroutine
!**************************************************************************


end module
