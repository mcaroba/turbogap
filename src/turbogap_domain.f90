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
   use turbogap_comm, only: comm_t, comm_bcast, comm_allgather, comm_sum_to_root
   use types, only: input_parameters, any_has_local_properties
   use timing, only: times_t, time_start, time_end
   use turbogap_structure, only: state_t
   use turbogap_setup, only: model_t
   use turbogap_loop, only: loop_t
   use turbogap_results, only: results_t
   use read_files, only: read_xyz
   use neighbors, only: build_neighbors_list

   implicit none

   private
   public :: domain_sync_state
   public :: domain_build
   public :: domain_complete_sites

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

!  Split the sites over the ranks and build each rank's neighbour list for its
!  share. Replicated data: a contiguous block of sites per rank, every rank
!  holding every position. The supercell is refreshed first because the box may
!  have crossed the cutoff sphere since the last build.
   subroutine domain_build(dom, nl, comm, state, params, model, loop, mc_file, time)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(inout) :: params
      type(model_t), intent(in) :: model
      type(loop_t), intent(inout) :: loop
      character(len=*), intent(in) :: mc_file
      type(times_t), intent(inout) :: time
      integer :: i

      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. loop%md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
      end if
      call time_start(time%neigh)
      !   Parallel neighbors list build
      call comm_bcast(comm, nl%rebuild_neighbors_list)

      !   If we're using a box rescaling algorithm or a barostat, then the box size can
      !   become smaller or bigger than the cutoff sphere. If that happens, and the current
      !   situation is different from before, then we need to figure out if we need to
      !   construct a supercell (i.e., the box was bigger than the cutoff sphere and now
      !   is smaller -> makes computations slower) or default back to the primitive unit cell
      !   (i.e., the box was smaller and now is bigger -> makes computations faster).
      !   We only need to check if rebuild_neighbors_list = .true.
      if (nl%rebuild_neighbors_list .and. params%do_mc .and. loop%mc_istep > 0) then
         call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                       model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                       state%positions, params%do_md, state%velocities, params%masses_types, state%masses, state%xyz_species, &
                       state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                       state%b_box, state%c_box, &
                       state%n_sites, .true., state%fix_atom, params%t_beg, &
                       params%write_array_property(6), .false., params%randomize_velocities)
      else if (nl%rebuild_neighbors_list) then
         call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                       model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                       state%positions, params%do_md, state%velocities, params%masses_types, state%masses, state%xyz_species, &
                       state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                       state%b_box, state%c_box, &
                       state%n_sites, .true., state%fix_atom, params%t_beg, params%write_array_property(6), &
                       .false., params%randomize_velocities)

      end if

      !   Overlapping domain decomposition with subcommunicators goes here <------------------- TO DO

      !   This is some trivial MPI parallelization to make sure the code works fine
      if (comm%rank < mod(state%n_sites, comm%size)) then
         dom%i_beg = 1 + comm%rank*(state%n_sites/comm%size + 1)
      else
         dom%i_beg = 1 + mod(state%n_sites, comm%size)*(state%n_sites/comm%size + 1) + (comm%rank - mod(state%n_sites, &
                                                                                               comm%size))*(state%n_sites/comm%size)
      end if
      if (comm%rank < mod(state%n_sites, comm%size)) then
         dom%i_end = (comm%rank + 1)*(state%n_sites/comm%size + 1)
      else
         dom%i_end = dom%i_beg + state%n_sites/comm%size - 1
      end if

      dom%do_list = .false.
      dom%do_list(dom%i_beg:dom%i_end) = .true.

      call build_neighbors_list(state%positions, state%a_box, state%b_box, state%c_box, params%do_timing, &
                                state%species_supercell, model%rcut_max, nl%n_atom_pairs, nl%rjs, &
                                nl%thetas, nl%phis, nl%xyz, nl%n_neigh_local, nl%neighbors_list, nl%neighbor_species, &
                                state%n_sites, state%indices, &
                                nl%rebuild_neighbors_list, dom%do_list, comm%rank)
      if (nl%rebuild_neighbors_list) then
         !     Get total number of atom pairs
         call comm_allgather(comm, nl%n_atom_pairs, dom%n_atom_pairs_by_rank)
         nl%n_atom_pairs_total = sum(dom%n_atom_pairs_by_rank)
         nl%n_atom_pairs = nl%n_atom_pairs_total

         !     Get number of neighbors
         if (.not. allocated(nl%n_neigh)) allocate (nl%n_neigh(1:state%n_sites))
         call comm_sum_to_root(comm, nl%n_neigh_local, nl%n_neigh, state%n_sites)
         call comm_bcast(comm, nl%n_neigh, state%n_sites)

         dom%j_beg = 1
         dom%j_end = dom%n_atom_pairs_by_rank(comm%rank + 1)
      end if
!   Store by which rank each site is being handled
      if (allocated(dom%site_in_rank)) then
         if (size(dom%site_in_rank) /= state%n_sites) then
            deallocate (dom%site_in_rank, dom%this_site_in_rank)
         end if
      end if
      if (.not. allocated(dom%site_in_rank)) then
         allocate (dom%site_in_rank(1:state%n_sites))
         allocate (dom%this_site_in_rank(1:state%n_sites))
      end if
      dom%site_in_rank = 0
      dom%this_site_in_rank = 0
      do i = dom%i_beg, dom%i_end
         dom%this_site_in_rank(i) = comm%rank
      end do
      call comm_sum_to_root(comm, dom%this_site_in_rank, dom%site_in_rank, state%n_sites)
      call comm_bcast(comm, dom%site_in_rank, state%n_sites)
      call time_end(time%neigh)
   end subroutine domain_build

!  Make every rank's copy of the per-site properties whole. Replicated data:
!  each rank computed its own sites and left zeros elsewhere, so a sum is the
!  whole reduction; rank 0 sums and broadcasts.
   subroutine domain_complete_sites(dom, comm, res, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), intent(inout) :: res
      type(state_t), intent(in) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time

      if (any_has_local_properties(model%soap_turbo_hypers)) then
         call time_start(time%mpi)
         call comm_sum_to_root(comm, res%local_properties, res%this_local_properties, state%n_sites*params%n_local_properties)
         res%local_properties = res%this_local_properties
         call comm_bcast(comm, res%local_properties, state%n_sites*params%n_local_properties)

         call time_end(time%mpi)
      end if

!     Each rank owns a slice of the sites, so its local_dipoles is zero
!     everywhere else and a plain sum is the whole reduction.
      if (params%do_dipole) then
         call time_start(time%mpi)
         call comm_sum_to_root(comm, res%local_dipoles, res%this_local_dipoles, 3*state%n_sites)
         res%local_dipoles = res%this_local_dipoles
         call comm_bcast(comm, res%local_dipoles, 3*state%n_sites)

         call comm_sum_to_root(comm, res%energies_dipole, res%this_energies_dipole, state%n_sites)
         res%energies_dipole = res%this_energies_dipole
         call comm_bcast(comm, res%energies_dipole, state%n_sites)
         call time_end(time%mpi)
      end if
   end subroutine domain_complete_sites

end module turbogap_domain
