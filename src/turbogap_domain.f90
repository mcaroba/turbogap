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
   use types, only: input_parameters, any_has_local_properties, contribution_ref
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
   public :: domain_complete_e0
   public :: domain_complete_contributions

!   The eleven additive contribution families that are reduced together after the
!   descriptor loop. Their predicates used to be written out three times -- to
!   count the slots, to pack them and to unpack them -- and evaluated
!   independently each time. Two copies disagreeing shifts counter2 and
!   silently attributes one term's energies to another. That is the same shape
!   as the ts+mbd predicate defect (KNOWN_ISSUES.md #1), so it is killed the
!   same way: evaluated once, into contrib_on, and only read thereafter.
   integer, parameter :: C_SOAP = 1
   integer, parameter :: C_VDW = 2
   integer, parameter :: C_ESTAT = 3
   integer, parameter :: C_LP = 4
   integer, parameter :: C_PDF = 5
   integer, parameter :: C_SF = 6
   integer, parameter :: C_XRD = 7
   integer, parameter :: C_ND = 8
   integer, parameter :: C_2B = 9
   integer, parameter :: C_CP = 10
   integer, parameter :: C_3B = 11
   integer, parameter :: N_CONTRIB = 11

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

!  Complete the e0 baseline each rank wrote for its own sites. Replicated data:
!  a sum to rank 0, whose copy is the one the output reads.
   subroutine domain_complete_e0(dom, comm, res, state, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), intent(inout) :: res
      type(state_t), intent(in) :: state
      type(times_t), intent(inout) :: time

      call time_start(time%mpi_ef)
      call comm_sum_to_root(comm, res%energies, res%this_energies, state%n_sites)
      call time_end(time%mpi_ef)
      res%energies = res%this_energies
   end subroutine domain_complete_e0

!  Complete every contribution family: each rank holds partial sums for the
!  sites it owns and the neighbours they touch. Replicated data: one packed
!  reduce to rank 0, whose totals are the ones the rest of the step reads.
   subroutine domain_complete_contributions(dom, comm, res, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), target, intent(inout) :: res
      type(state_t), intent(in) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time
      real(dp), allocatable :: all_energies(:, :)
      real(dp), allocatable :: all_forces(:, :, :)
      real(dp), allocatable :: all_virial(:, :, :)
      real(dp), allocatable :: all_this_energies(:, :)
      real(dp), allocatable :: all_this_forces(:, :, :)
      real(dp), allocatable :: all_this_virial(:, :, :)
      logical :: contrib_on(1:N_CONTRIB)
      type(contribution_ref) :: contrib(1:N_CONTRIB)
      integer :: n_active
      integer :: i_contrib
      integer :: counter2

      call time_start(time%mpi_ef)
!       One evaluation of the eleven predicates, and one list built from them.
!       The pack and unpack walks below read only that list, so they cannot
!       disagree about which slot belongs to which family -- the failure mode
!       this replaces was three independent copies of these conditions, where
!       any two disagreeing shifts the slot numbering and silently attributes
!       one family's energies and forces to another.
      contrib_on(C_SOAP) = (model%n_soap_turbo > 0)
      contrib_on(C_VDW) = allocated(res%this_energies_vdw)
      contrib_on(C_ESTAT) = allocated(res%this_energies_estat)
      contrib_on(C_LP) = allocated(res%this_energies_lp)
      contrib_on(C_PDF) = allocated(res%this_energies_pdf) .and. params%valid_pdf
      contrib_on(C_SF) = allocated(res%this_energies_sf) .and. params%valid_sf
      contrib_on(C_XRD) = allocated(res%this_energies_xrd) .and. params%valid_xrd
      contrib_on(C_ND) = allocated(res%this_energies_nd) .and. params%valid_nd
      contrib_on(C_2B) = (model%n_distance_2b > 0)
      contrib_on(C_CP) = (model%n_core_pot > 0)
      contrib_on(C_3B) = (model%n_angle_3b > 0)

      n_active = 0
      if (contrib_on(C_SOAP)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%energies_soap
         contrib(n_active)%e_dst => res%energies_soap
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%forces_soap
            contrib(n_active)%v_src => res%virial_soap
            contrib(n_active)%f_dst => res%forces_soap
            contrib(n_active)%v_dst => res%virial_soap
         end if
      end if
      if (contrib_on(C_VDW)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_vdw
         contrib(n_active)%e_dst => res%energies_vdw
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_vdw
            contrib(n_active)%v_src => res%this_virial_vdw
            contrib(n_active)%f_dst => res%forces_vdw
            contrib(n_active)%v_dst => res%virial_vdw
         end if
      end if
      if (contrib_on(C_ESTAT)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_estat
         contrib(n_active)%e_dst => res%energies_estat
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_estat
            contrib(n_active)%v_src => res%this_virial_estat
            contrib(n_active)%f_dst => res%forces_estat
            contrib(n_active)%v_dst => res%virial_estat
         end if
      end if
      if (contrib_on(C_LP)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_lp
         contrib(n_active)%e_dst => res%energies_lp
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_lp
            contrib(n_active)%v_src => res%this_virial_lp
            contrib(n_active)%f_dst => res%forces_lp
            contrib(n_active)%v_dst => res%virial_lp
         end if
      end if
      if (contrib_on(C_PDF)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_pdf
         contrib(n_active)%e_dst => res%energies_pdf
         contrib(n_active)%forces = params%do_forces .and. params%exp_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_pdf
            contrib(n_active)%v_src => res%this_virial_pdf
            contrib(n_active)%f_dst => res%forces_pdf
            contrib(n_active)%v_dst => res%virial_pdf
         end if
      end if
      if (contrib_on(C_SF)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_sf
         contrib(n_active)%e_dst => res%energies_sf
         contrib(n_active)%forces = params%do_forces .and. params%exp_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_sf
            contrib(n_active)%v_src => res%this_virial_sf
            contrib(n_active)%f_dst => res%forces_sf
            contrib(n_active)%v_dst => res%virial_sf
         end if
      end if
      if (contrib_on(C_XRD)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_xrd
         contrib(n_active)%e_dst => res%energies_xrd
         contrib(n_active)%forces = params%do_forces .and. params%exp_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_xrd
            contrib(n_active)%v_src => res%this_virial_xrd
            contrib(n_active)%f_dst => res%forces_xrd
            contrib(n_active)%v_dst => res%virial_xrd
         end if
      end if
      if (contrib_on(C_ND)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%this_energies_nd
         contrib(n_active)%e_dst => res%energies_nd
         contrib(n_active)%forces = params%do_forces .and. params%exp_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%this_forces_nd
            contrib(n_active)%v_src => res%this_virial_nd
            contrib(n_active)%f_dst => res%forces_nd
            contrib(n_active)%v_dst => res%virial_nd
         end if
      end if
      if (contrib_on(C_2B)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%energies_2b
         contrib(n_active)%e_dst => res%energies_2b
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%forces_2b
            contrib(n_active)%v_src => res%virial_2b
            contrib(n_active)%f_dst => res%forces_2b
            contrib(n_active)%v_dst => res%virial_2b
         end if
      end if
      if (contrib_on(C_CP)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%energies_core_pot
         contrib(n_active)%e_dst => res%energies_core_pot
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%forces_core_pot
            contrib(n_active)%v_src => res%virial_core_pot
            contrib(n_active)%f_dst => res%forces_core_pot
            contrib(n_active)%v_dst => res%virial_core_pot
         end if
      end if
      if (contrib_on(C_3B)) then
         n_active = n_active + 1
         contrib(n_active)%e_src => res%energies_3b
         contrib(n_active)%e_dst => res%energies_3b
         contrib(n_active)%forces = params%do_forces
         if (contrib(n_active)%forces) then
            contrib(n_active)%f_src => res%forces_3b
            contrib(n_active)%v_src => res%virial_3b
            contrib(n_active)%f_dst => res%forces_3b
            contrib(n_active)%v_dst => res%virial_3b
         end if
      end if

      counter2 = n_active

      allocate (all_energies(1:state%n_sites, 1:counter2))
      allocate (all_this_energies(1:state%n_sites, 1:counter2))
      if (params%do_forces) then
         allocate (all_forces(1:3, 1:state%n_sites, 1:counter2))
         allocate (all_this_forces(1:3, 1:state%n_sites, 1:counter2))
         allocate (all_virial(1:3, 1:3, 1:counter2))
         allocate (all_this_virial(1:3, 1:3, 1:counter2))
      end if

!       Pack. A family owns a slot whenever it is active, but only contributes
!       forces when it carries them -- the exp-spectra families additionally
!       need exp_forces. Their slot must still be cleared: all_forces is
!       allocated and never zeroed, and mpi_reduce below reads the whole array
!       regardless of who wrote what into it.
      do i_contrib = 1, n_active
         all_energies(1:state%n_sites, i_contrib) = contrib(i_contrib)%e_src(1:state%n_sites)
         if (contrib(i_contrib)%forces) then
            all_forces(1:3, 1:state%n_sites, i_contrib) = contrib(i_contrib)%f_src(1:3, 1:state%n_sites)
            all_virial(1:3, 1:3, i_contrib) = contrib(i_contrib)%v_src(1:3, 1:3)
         else if (params%do_forces) then
            all_forces(1:3, 1:state%n_sites, i_contrib) = 0.d0
            all_virial(1:3, 1:3, i_contrib) = 0.d0
         end if
      end do

      !       Here we communicate
      call comm_sum_to_root(comm, all_energies, all_this_energies, state%n_sites*counter2)
      if (params%do_forces) then
         call comm_sum_to_root(comm, all_forces, all_this_forces, 3*state%n_sites*counter2)
         call comm_sum_to_root(comm, all_virial, all_this_virial, 9*counter2)
      end if

!       Unpack. For the six families packed from a this_ array this is where
!       the reduced result lands in the un-prefixed one.
      do i_contrib = 1, n_active
         contrib(i_contrib)%e_dst(1:state%n_sites) = all_this_energies(1:state%n_sites, i_contrib)
         if (contrib(i_contrib)%forces) then
            contrib(i_contrib)%f_dst(1:3, 1:state%n_sites) = all_this_forces(1:3, 1:state%n_sites, i_contrib)
            contrib(i_contrib)%v_dst(1:3, 1:3) = all_this_virial(1:3, 1:3, i_contrib)
         end if
      end do

!       Release the this_ arrays now that their contents have been unpacked.
!       Kept explicit rather than folded into the loop: an allocatable cannot
!       be deallocated through a pointer, and this_local_virial_vdw_diag has no
!       counterpart in the list.
      if (contrib_on(C_VDW)) then
         deallocate (res%this_energies_vdw)
         if (params%do_forces) deallocate (res%this_forces_vdw, res%this_local_virial_vdw_diag)
      end if
      if (contrib_on(C_ESTAT)) then
         deallocate (res%this_energies_estat)
         if (params%do_forces) deallocate (res%this_forces_estat)
      end if
      if (contrib_on(C_LP)) then
         deallocate (res%this_energies_lp)
         if (params%do_forces) deallocate (res%this_forces_lp)
      end if
      if (contrib_on(C_PDF)) then
         deallocate (res%this_energies_pdf)
         if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_pdf)
      end if
      if (contrib_on(C_SF)) then
         deallocate (res%this_energies_sf)
         if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_sf)
      end if
      if (contrib_on(C_XRD)) then
         deallocate (res%this_energies_xrd)
         if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_xrd)
      end if
      if (contrib_on(C_ND)) then
         deallocate (res%this_energies_nd)
         if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_nd)
      end if

      !       Clean up
      deallocate (all_energies, all_this_energies)
      if (params%do_forces) then
         deallocate (all_forces, all_this_forces, all_virial, all_this_virial)
      end if

      call time_end(time%mpi_ef)
   end subroutine domain_complete_contributions

end module turbogap_domain
