! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_domain.f90, is copyright (c) 2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
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

!  How the sites are split over ranks, and the neighbour list each rank holds
!  for its share. The kernels see the split only as i_beg:i_end over sites and
!  j_beg:j_end over this rank's pairs.
module turbogap_domain

   use kinds, only: dp
   use turbogap_comm, only: comm_t, comm_bcast
   use types, only: input_parameters
   use timing, only: times_t, time_start, time_end
   use turbogap_structure, only: state_t

   implicit none

   private
   public :: domain_sync_state

   type, public :: neighbors_t
      real(dp), allocatable :: rjs(:)
      real(dp), allocatable :: thetas(:)
      real(dp), allocatable :: phis(:)
      real(dp), allocatable :: xyz(:, :)
      integer, allocatable :: n_neigh(:)
      integer, allocatable :: n_neigh_local(:)
      integer, allocatable :: neighbors_list(:)
      integer, allocatable :: neighbor_species(:)
      integer :: n_atom_pairs
      integer :: n_atom_pairs_total
      logical :: rebuild_neighbors_list = .true.
   end type neighbors_t

   type, public :: domain_t
      integer :: i_beg
      integer :: i_end
      integer :: j_beg
      integer :: j_end
      logical, allocatable :: do_list(:)
      integer, allocatable :: site_in_rank(:)
      integer, allocatable :: this_site_in_rank(:)
      integer, allocatable :: n_atom_pairs_by_rank(:)
      integer :: n_atom_pairs_by_rank_prev = 0
   end type domain_t

contains

!  Give every rank the structure rank 0 holds. Replicated data: the whole of it,
!  sized from rank 0's arrays.
   subroutine domain_sync_state(dom, comm, state, params, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(times_t), intent(inout) :: time
      integer :: n_pos
      integer :: n_sp
      integer :: n_sp_sc

      if (comm%rank == 0) then
         n_pos = size(state%positions, 2)
         n_sp = size(state%xyz_species, 1)
         n_sp_sc = size(state%xyz_species_supercell, 1)
      end if
      call time_start(time%mpi)
      call comm_bcast(comm, n_pos)
      call comm_bcast(comm, n_sp)
      call comm_bcast(comm, n_sp_sc)
      call comm_bcast(comm, state%n_sites)
      call time_end(time%mpi)

      if (comm%rank /= 0) then
         if (allocated(state%positions)) deallocate (state%positions)
         allocate (state%positions(1:3, n_pos))
         if (params%do_md .or. params%do_nested_sampling .or. params%do_mc) then
            if (allocated(state%velocities)) deallocate (state%velocities)
            allocate (state%velocities(1:3, n_pos))
            if (allocated(state%masses)) deallocate (state%masses)
            allocate (state%masses(1:n_sp))
         end if
         if (allocated(state%xyz_species)) deallocate (state%xyz_species)
         allocate (state%xyz_species(1:n_sp))
         if (allocated(state%species)) deallocate (state%species)
         allocate (state%species(1:n_sp))
         if (allocated(state%xyz_species_supercell)) deallocate (state%xyz_species_supercell)
         allocate (state%xyz_species_supercell(1:n_sp_sc))
         if (allocated(state%species_supercell)) deallocate (state%species_supercell)
         allocate (state%species_supercell(1:n_sp_sc))
         if (allocated(state%fix_atom)) deallocate (state%fix_atom)
         allocate (state%fix_atom(1:3, 1:n_sp))
      end if
      call time_start(time%mpi_positions)
      call comm_bcast(comm, state%positions, 3*n_pos)
      if (params%do_md .or. params%do_nested_sampling .or. params%do_mc .or. params%mc_hamiltonian) then
         call comm_bcast(comm, state%velocities, 3*n_pos)
         call comm_bcast(comm, state%masses, n_sp)
         call comm_bcast(comm, state%fix_atom, 3*n_sp)
      end if
      call comm_bcast(comm, state%xyz_species, 8*n_sp)
      call comm_bcast(comm, state%xyz_species_supercell, 8*n_sp_sc)
      call comm_bcast(comm, state%species, n_sp)
      call comm_bcast(comm, state%species_supercell, n_sp_sc)
      call comm_bcast(comm, state%indices, 3)
      call comm_bcast(comm, state%a_box, 3)
      call comm_bcast(comm, state%b_box, 3)
      call comm_bcast(comm, state%c_box, 3)
      call time_end(time%mpi_positions)
   end subroutine domain_sync_state

end module turbogap_domain
