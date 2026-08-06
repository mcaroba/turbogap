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

   use kinds

   use types
   use gap
   use timing

   implicit none

   private
   public :: gap_backend_begin
   public :: gap_backend_end
   public :: add_2b_contribution
   public :: add_core_pot_contribution
   public :: add_3b_contribution

contains

!**************************************************************************
! Bracket the three contribution calls.
!
! Nothing to do on the CPU: the arrays the kernels read are already in host
! memory. These exist because the GPU implementation needs somewhere to upload
! the neighbour data once for all three calls rather than three times, and that
! has to happen without the driver ever holding a device pointer. An empty pair
! here is the price of keeping the interface physics-only.
   subroutine gap_backend_begin(params, rjs, xyz, n_neigh, species, neighbor_species, &
                                neighbors_list, i_beg, i_end, j_beg, j_end)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
   end subroutine gap_backend_begin

   subroutine gap_backend_end()
      implicit none
   end subroutine gap_backend_end

!**************************************************************************
! Accumulate the two-body energies, forces and virial.
   subroutine add_2b_contribution(n_distance_2b, distance_2b_hypers, &
                                  params, rjs, xyz, n_neigh, species, neighbor_species, &
                                  i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                  energies_2b, forces_2b, virial_2b, time)

      implicit none

!   ---- Input: describes the system and this descriptor set ----
      integer, intent(in) :: n_distance_2b
      type(distance_2b), allocatable, intent(inout) :: distance_2b_hypers(:)
      type(input_parameters), intent(inout) :: params
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end

!   ---- Output: accumulated into by this routine ----
      real(dp), intent(inout), allocatable :: energies_2b(:)
      real(dp), intent(inout), allocatable :: forces_2b(:, :)
      real(dp), intent(inout) :: virial_2b(1:3, 1:3)
      type(times_t), intent(inout) :: time

!   ---- Scratch: owned by the driver and reused across the three calls,
!        rather than local, so this does not add an allocation per descriptor ----
      real(dp), intent(inout), allocatable :: this_energies(:)
      real(dp), intent(inout), allocatable :: this_forces(:, :)
      real(dp), intent(inout) :: this_virial(1:3, 1:3)

!   ---- Internal ----
      integer :: i

      do i = 1, n_distance_2b
         call time_start(time%gap_2b)
         this_energies = 0.d0
         if (params%do_forces) then
            this_forces = 0.d0
            this_virial = 0.d0
         end if
         call get_2b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), distance_2b_hypers(i)%alphas, &
                                       distance_2b_hypers(i)%cutoff, &
                                       distance_2b_hypers(i)%rcut, 0.5d0, distance_2b_hypers(i)%delta, &
                                       distance_2b_hypers(i)%sigma, 0.d0, distance_2b_hypers(i)%Qs(:, 1), &
                                       n_neigh(i_beg:i_end), params%do_forces, params%do_timing, &
                                       species(i_beg:i_end), neighbor_species(j_beg:j_end), &
                                       distance_2b_hypers(i)%species1, distance_2b_hypers(i)%species2, &
                                       params%species_types, this_energies(i_beg:i_end), this_forces(1:3, i_beg:i_end), &
                                       this_virial)
         energies_2b = energies_2b + this_energies
         if (params%do_forces) then
            forces_2b = forces_2b + this_forces
            virial_2b = virial_2b + this_virial
         end if
         call time_end(time%gap_2b)
      end do

   end subroutine add_2b_contribution

!**************************************************************************
! Accumulate the core-potential energies, forces and virial.
   subroutine add_core_pot_contribution(n_core_pot, core_pot_hypers, &
                                        params, rjs, xyz, n_neigh, species, neighbor_species, &
                                        i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                        energies_core_pot, forces_core_pot, virial_core_pot, time)

      implicit none

!   ---- Input: describes the system and this descriptor set ----
      integer, intent(in) :: n_core_pot
      type(core_pot), allocatable, intent(inout) :: core_pot_hypers(:)
      type(input_parameters), intent(inout) :: params
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end

!   ---- Output: accumulated into by this routine ----
      real(dp), intent(inout), allocatable :: energies_core_pot(:)
      real(dp), intent(inout), allocatable :: forces_core_pot(:, :)
      real(dp), intent(inout) :: virial_core_pot(1:3, 1:3)
      type(times_t), intent(inout) :: time

!   ---- Scratch: owned by the driver and reused across the three calls,
!        rather than local, so this does not add an allocation per descriptor ----
      real(dp), intent(inout), allocatable :: this_energies(:)
      real(dp), intent(inout), allocatable :: this_forces(:, :)
      real(dp), intent(inout) :: this_virial(1:3, 1:3)

!   ---- Internal ----
      integer :: i

      do i = 1, n_core_pot
         call time_start(time%gap_core_pot)
         this_energies = 0.d0
         if (params%do_forces) then
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
                                             this_virial)
         energies_core_pot = energies_core_pot + this_energies
         if (params%do_forces) then
            forces_core_pot = forces_core_pot + this_forces
            virial_core_pot = virial_core_pot + this_virial
         end if
         call time_end(time%gap_core_pot)
      end do

   end subroutine add_core_pot_contribution

!**************************************************************************
! Accumulate the three-body energies, forces and virial.
   subroutine add_3b_contribution(n_angle_3b, angle_3b_hypers, neighbors_list, &
                                  params, rjs, xyz, n_neigh, species, neighbor_species, &
                                  i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                  energies_3b, forces_3b, virial_3b, time)

      implicit none

!   ---- Input: describes the system and this descriptor set ----
      integer, intent(in) :: n_angle_3b
      type(angle_3b), allocatable, intent(inout) :: angle_3b_hypers(:)
      type(input_parameters), intent(inout) :: params
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end

!   ---- Output: accumulated into by this routine ----
      real(dp), intent(inout), allocatable :: energies_3b(:)
      real(dp), intent(inout), allocatable :: forces_3b(:, :)
      real(dp), intent(inout) :: virial_3b(1:3, 1:3)
      type(times_t), intent(inout) :: time

!   ---- Scratch: owned by the driver and reused across the three calls,
!        rather than local, so this does not add an allocation per descriptor ----
      real(dp), intent(inout), allocatable :: this_energies(:)
      real(dp), intent(inout), allocatable :: this_forces(:, :)
      real(dp), intent(inout) :: this_virial(1:3, 1:3)

!   ---- Internal ----
      integer :: i

      do i = 1, n_angle_3b
         call time_start(time%gap_3b)
         this_energies = 0.d0
         if (params%do_forces) then
            this_forces = 0.d0
            this_virial = 0.d0
         end if
         call get_3b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), angle_3b_hypers(i)%alphas, &
                                       angle_3b_hypers(i)%cutoff, &
                                       angle_3b_hypers(i)%rcut, 0.5d0, angle_3b_hypers(i)%delta, &
                                       angle_3b_hypers(i)%sigma, 0.d0, angle_3b_hypers(i)%Qs, n_neigh(i_beg:i_end), &
                                       neighbors_list(j_beg:j_end), &
                                       params%do_forces, params%do_timing, angle_3b_hypers(i)%kernel_type, &
                                       species(i_beg:i_end), neighbor_species(j_beg:j_end), angle_3b_hypers(i)%species_center, &
                                       angle_3b_hypers(i)%species1, angle_3b_hypers(i)%species2, params%species_types, &
                                       this_energies(i_beg:i_end), this_forces, this_virial)
         energies_3b = energies_3b + this_energies
         if (params%do_forces) then
            forces_3b = forces_3b + this_forces
            virial_3b = virial_3b + this_virial
         end if
         call time_end(time%gap_3b)
      end do

   end subroutine add_3b_contribution

end module gap_backend
