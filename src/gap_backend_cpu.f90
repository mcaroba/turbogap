! Copyright (c) 2020-2026 by Albert Bartok and Miguel Caro
!
! The two-body, core-potential and three-body contributions to the energy,
! forces and virial: the CPU implementation of the gap_backend seam.
!
! There are two files defining `module gap_backend` -- this one and
! gap_backend_gpu.f90 on the GPU branch -- with the same three public
! procedures and the same argument lists. The Makefile compiles exactly one.
! That is the mechanism by which a single turbogap.f90 can eventually serve
! both branches: the driver calls these names and never learns which
! implementation it got.
!
! The interface is deliberately physics only. No device pointer, stream or
! cuBLAS handle appears in it; those belong inside the GPU implementation.
! The GPU branch's current internal procedures take ten such buffers from the
! driver by host association, and moving that ownership inside is the work
! this interface defines rather than the work it does.
!
! Lifted verbatim from turbogap.f90, where the three loops ran inline between
! the "Loop through distance_2b descriptors" banner and the MPI reduce block.
module gap_backend

  use types
  use gap
  use timing

  implicit none

  private
  public :: add_2b_contribution, add_core_pot_contribution, add_3b_contribution

contains

!**************************************************************************
! Accumulate the two-body energies, forces and virial.
  subroutine add_2b_contribution( n_distance_2b, distance_2b_hypers, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       energies_2b, forces_2b, virial_2b, time_2b )

    implicit none

    integer, intent(in) :: n_distance_2b
    type(distance_2b), allocatable, intent(inout) :: distance_2b_hypers(:)
    type(input_parameters), intent(inout) :: params
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: n_neigh(:), species(:), neighbor_species(:)
    integer, intent(in) :: i_beg, i_end, j_beg, j_end
!   Scratch owned by the driver and reused across the three calls, rather than
!   local, so the extraction does not add an allocation per descriptor.
    real*8, intent(inout), allocatable :: this_energies(:), this_forces(:,:)
    real*8, intent(inout) :: this_virial(1:3,1:3)
    real*8, intent(inout), allocatable :: energies_2b(:), forces_2b(:,:)
    real*8, intent(inout) :: virial_2b(1:3,1:3), time_2b(1:3)

    integer :: i

    do i = 1, n_distance_2b
       call get_time(time_2b(1))
       this_energies = 0.d0
       if( params%do_forces )then
          this_forces = 0.d0
          this_virial = 0.d0
       end if
       call get_2b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), distance_2b_hypers(i)%alphas, &
            distance_2b_hypers(i)%cutoff, &
            distance_2b_hypers(i)%rcut, 0.5d0, distance_2b_hypers(i)%delta, &
            distance_2b_hypers(i)%sigma, 0.d0, distance_2b_hypers(i)%Qs(:,1), &
            n_neigh(i_beg:i_end), params%do_forces, params%do_timing, &
            species(i_beg:i_end), neighbor_species(j_beg:j_end), &
            distance_2b_hypers(i)%species1, distance_2b_hypers(i)%species2, &
            params%species_types, this_energies(i_beg:i_end), this_forces(1:3, i_beg:i_end), &
            this_virial )
       energies_2b = energies_2b + this_energies
       if( params%do_forces )then
          forces_2b = forces_2b + this_forces
          virial_2b = virial_2b + this_virial
       end if
       call get_time(time_2b(2))
       time_2b(3) = time_2b(3) + time_2b(2) - time_2b(1)
    end do

  end subroutine add_2b_contribution

!**************************************************************************
! Accumulate the core-potential energies, forces and virial.
  subroutine add_core_pot_contribution( n_core_pot, core_pot_hypers, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       energies_core_pot, forces_core_pot, virial_core_pot, time_core_pot )

    implicit none

    integer, intent(in) :: n_core_pot
    type(core_pot), allocatable, intent(inout) :: core_pot_hypers(:)
    type(input_parameters), intent(inout) :: params
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: n_neigh(:), species(:), neighbor_species(:)
    integer, intent(in) :: i_beg, i_end, j_beg, j_end
!   Scratch owned by the driver and reused across the three calls, rather than
!   local, so the extraction does not add an allocation per descriptor.
    real*8, intent(inout), allocatable :: this_energies(:), this_forces(:,:)
    real*8, intent(inout) :: this_virial(1:3,1:3)
    real*8, intent(inout), allocatable :: energies_core_pot(:), forces_core_pot(:,:)
    real*8, intent(inout) :: virial_core_pot(1:3,1:3), time_core_pot(1:3)

    integer :: i

    do i = 1, n_core_pot
       call get_time(time_core_pot(1))
       this_energies = 0.d0
       if( params%do_forces )then
          this_forces = 0.d0
          this_virial = 0.d0
       end if
       call get_core_pot_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), &
            core_pot_hypers(i)%x, core_pot_hypers(i)%V, &
            core_pot_hypers(i)%yp1, core_pot_hypers(i)%ypn, &
            core_pot_hypers(i)%dVdx2, n_neigh(i_beg:i_end), params%do_forces, &
            params%do_timing, species(i_beg:i_end), neighbor_species(j_beg:j_end), &
            core_pot_hypers(i)%species1, core_pot_hypers(i)%species2, &
            params%species_types, this_energies(i_beg:i_end), this_forces(1:3, i_beg:i_end), &
            this_virial )
       energies_core_pot = energies_core_pot + this_energies
       if( params%do_forces )then
          forces_core_pot = forces_core_pot + this_forces
          virial_core_pot = virial_core_pot + this_virial
       end if
       call get_time(time_core_pot(2))
       time_core_pot(3) = time_core_pot(3) + time_core_pot(2) - time_core_pot(1)
    end do

  end subroutine add_core_pot_contribution

!**************************************************************************
! Accumulate the three-body energies, forces and virial.
  subroutine add_3b_contribution( n_angle_3b, angle_3b_hypers, neighbors_list, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       energies_3b, forces_3b, virial_3b, time_3b )

    implicit none

    integer, intent(in) :: n_angle_3b
    type(angle_3b), allocatable, intent(inout) :: angle_3b_hypers(:)
    integer, intent(in), allocatable :: neighbors_list(:)
    type(input_parameters), intent(inout) :: params
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: n_neigh(:), species(:), neighbor_species(:)
    integer, intent(in) :: i_beg, i_end, j_beg, j_end
!   Scratch owned by the driver and reused across the three calls, rather than
!   local, so the extraction does not add an allocation per descriptor.
    real*8, intent(inout), allocatable :: this_energies(:), this_forces(:,:)
    real*8, intent(inout) :: this_virial(1:3,1:3)
    real*8, intent(inout), allocatable :: energies_3b(:), forces_3b(:,:)
    real*8, intent(inout) :: virial_3b(1:3,1:3), time_3b(1:3)

    integer :: i

    do i = 1, n_angle_3b
       call get_time(time_3b(1))
       this_energies = 0.d0
       if( params%do_forces )then
          this_forces = 0.d0
          this_virial = 0.d0
       end if
       call get_3b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), angle_3b_hypers(i)%alphas, &
            angle_3b_hypers(i)%cutoff, &
            angle_3b_hypers(i)%rcut, 0.5d0, angle_3b_hypers(i)%delta, &
            angle_3b_hypers(i)%sigma, 0.d0, angle_3b_hypers(i)%Qs, n_neigh(i_beg:i_end), &
            neighbors_list(j_beg:j_end), &
            params%do_forces, params%do_timing, angle_3b_hypers(i)%kernel_type, &
            species(i_beg:i_end), neighbor_species(j_beg:j_end), angle_3b_hypers(i)%species_center, &
            angle_3b_hypers(i)%species1, angle_3b_hypers(i)%species2, params%species_types, &
            this_energies(i_beg:i_end), this_forces, this_virial)
       energies_3b = energies_3b + this_energies
       if( params%do_forces )then
          forces_3b = forces_3b + this_forces
          virial_3b = virial_3b + this_virial
       end if
       call get_time(time_3b(2))
       time_3b(3) = time_3b(3) + time_3b(2) - time_3b(1)
    end do

  end subroutine add_3b_contribution

end module gap_backend
