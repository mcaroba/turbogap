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
   use turbogap_comm, only: comm_t, comm_bcast, comm_allgather, comm_sum_to_root, comm_with_mpi
   use types, only: input_parameters, any_has_local_properties, contribution_ref, perform_t
   use timing, only: times_t, time_start, time_end
   use turbogap_structure, only: state_t
   use turbogap_setup, only: model_t
   use turbogap_loop, only: loop_t
   use turbogap_results, only: results_t, results_resize_sites
   use domain_legacy, only: legacy_dd_t, legacy_init, legacy_free, legacy_check_widths, legacy_distribute, &
                            legacy_reverse, legacy_gather
   use error, only: turbogap_abort
   use read_files, only: read_xyz
   use neighbors, only: build_neighbors_list

   implicit none

   private
   public :: domain_init
   public :: domain_configure
   public :: domain_free
   public :: domain_sync_state
   public :: domain_sync_after_md
   public :: domain_sync_after_ipi
   public :: domain_end_step
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

!   dd_type: everything on every rank, or the legacy spatial domains.
   integer, parameter, public :: DD_NONE = 0
   integer, parameter, public :: DD_LEGACY = 1

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
      integer :: dd_type = DD_NONE
      type(legacy_dd_t) :: legacy
!     Legacy, rank 0: the whole structure, held here while state has the
!     domain's for the evaluation.
      integer :: g_n_sites = 0
      real(dp), allocatable :: g_positions(:, :)
      integer, allocatable :: g_species(:)
      character*8, allocatable :: g_xyz_species(:)
   end type domain_t

contains

!  Replicated data: every rank holds every atom.
   subroutine domain_init(dom, comm)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm

      allocate (dom%n_atom_pairs_by_rank(1:comm%size))
   end subroutine domain_init

!  Choose the decomposition, once the input has been read.
   subroutine domain_configure(dom, comm, params, model, perform, mode)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(perform_t), intent(in) :: perform
      character(len=*), intent(in) :: mode

      select case (trim(params%dd_type))
      case ("none")
         dom%dd_type = DD_NONE
      case ("legacy")
         dom%dd_type = DD_LEGACY
         call legacy_refuse(comm, params, model, perform, mode)
         call legacy_init(dom%legacy, comm, params%dd_grid, params%dd_grid_affinity)
      case default
         if (comm%rank == 0) write (*, '(A,A,A)') 'ERROR: dd_type is "', trim(params%dd_type), &
            '"; it must be none or legacy'
         call turbogap_abort()
      end select
   end subroutine domain_configure

!  What the legacy decomposition cannot compute correctly is refused, not run;
!  docs/LEGACY_DD.md says why for each.
   subroutine legacy_refuse(comm, params, model, perform, mode)
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(perform_t), intent(in) :: perform
      character(len=*), intent(in) :: mode
      character(len=64) :: why

      why = ""
#ifdef _GPU
      why = "the device build"
#endif
      if (.not. comm_with_mpi) why = "a build without MPI"
      if (.not. params%do_prediction) why = "a run without do_prediction"
      if (params%write_soap .or. params%write_derivatives) why = "write_soap or write_derivatives"
      if (params%do_mc) why = "mc"
      if (params%do_nested_sampling) why = "nested sampling"
      if (trim(mode) == "ipi") why = "i-PI"
      if (trim(params%vdw_type) /= "none") why = "van der Waals"
      if (model%valid_estat_charges .or. trim(params%estat_method) /= "none") why = "electrostatics"
      if (params%do_exp .or. params%do_xps .or. model%valid_xps .or. perform%pdf .or. perform%sf .or. &
          perform%xrd .or. perform%nd) why = "experimental observables"
      if (params%do_dipole .or. params%valid_ir) why = "dipoles or MAD-IR"
      if (why /= "") then
         if (comm%rank == 0) write (*, '(A,A)') 'ERROR: dd_type legacy does not support ', trim(why)
         call turbogap_abort()
      end if
   end subroutine legacy_refuse

!  Give every rank the structure rank 0 holds. Replicated data: the whole of it,
!  sized from rank 0's arrays. Sending only what an MC move changed was tried,
!  keyed by a broadcast flag array, and segfaulted; hence everything, every time.
   subroutine domain_sync_state(dom, comm, state, params, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(times_t), intent(inout) :: time
      integer :: n_pos
      integer :: n_sp
      integer :: n_sp_sc

!     Legacy: domain_build hands each rank its domain.
      if (dom%dd_type == DD_LEGACY) return
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

!  After the integrator. Replicated data: rank 0 moved the atoms, so its
!  positions and velocities, and whether the neighbour list went stale, go to
!  every rank.
   subroutine domain_sync_after_md(dom, comm, state, nl, params, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(neighbors_t), intent(inout) :: nl
      type(input_parameters), intent(in) :: params
      type(times_t), intent(inout) :: time
      integer :: n_pos

!     Legacy: the halos move with the atoms, so every step rebuilds.
      if (dom%dd_type == DD_LEGACY) then
         nl%rebuild_neighbors_list = .true.
         return
      end if
      call comm_bcast(comm, nl%rebuild_neighbors_list)
      if (params%do_md) then
         call time_start(time%mpi_positions)
         n_pos = size(state%positions, 2)
         call comm_bcast(comm, state%positions, 3*n_pos)
         call comm_bcast(comm, state%velocities, 3*n_pos)
         call time_end(time%mpi_positions)
      end if
   end subroutine domain_sync_after_md

!  After an i-PI exchange: rank 0 holds i-PI's next configuration, or the
!  instruction to stop.
   subroutine domain_sync_after_ipi(dom, comm, state, nl, loop)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(neighbors_t), intent(inout) :: nl
      type(loop_t), intent(inout) :: loop
      integer :: n_pos

      call comm_bcast(comm, loop%exit_loop)
      call comm_bcast(comm, nl%rebuild_neighbors_list)
      if (.not. loop%exit_loop) then
         n_pos = size(state%positions, 2)
         call comm_bcast(comm, state%positions, 3*n_pos)
         call comm_bcast(comm, state%a_box, 3)
         call comm_bcast(comm, state%b_box, 3)
         call comm_bcast(comm, state%c_box, 3)
      end if
   end subroutine domain_sync_after_ipi

!  End of a pass: the site mask for the next one, whether the neighbour list
!  is to be rebuilt, and the arrays a finished MD or MC run no longer needs.
   subroutine domain_end_step(dom, nl, comm, state, params, loop)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(loop_t), intent(in) :: loop

      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. loop%md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
      end if
      !   Parallel neighbors list build
      call comm_bcast(comm, nl%rebuild_neighbors_list)

      if (nl%rebuild_neighbors_list) then
         deallocate (nl%rjs, nl%xyz, nl%thetas, nl%phis, nl%neighbor_species)
         deallocate (nl%neighbors_list, nl%n_neigh)
         deallocate (nl%n_neigh_local)
      end if
      if ((params%do_nested_sampling .and. .not. params%do_mc) .and. &
          (params%do_md .and. (loop%md_istep == params%md_nsteps .or. loop%exit_loop))) then
         deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, state%species_supercell, &
                     dom%do_list)
         if (allocated(state%velocities)) deallocate (state%velocities)
      end if
      if (params%do_mc .and. params%do_md) then
         if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop)) then
            deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, &
                        state%species_supercell, dom%do_list)
            if (allocated(state%velocities)) deallocate (state%velocities)
         end if
      end if

      if ((params%do_md .and. .not. params%do_mc) .and. &
          (loop%md_istep == params%md_nsteps .or. loop%exit_loop) .and. comm%rank == 0) then
         deallocate (state%positions_prev, state%forces_prev)
      end if
      if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop) .and. comm%rank == 0) then
         if (allocated(state%forces_prev)) deallocate (state%forces_prev)
         if (allocated(state%positions_prev)) deallocate (state%positions_prev)
      end if
      dom%n_atom_pairs_by_rank_prev = dom%n_atom_pairs_by_rank(comm%rank + 1)
   end subroutine domain_end_step

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

      if (dom%dd_type == DD_LEGACY) then
         call legacy_build(dom, nl, comm, state, params, model, loop, time)
         return
      end if
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

!  Legacy: rank 0 keeps the whole structure aside and every rank gets its
!  domain's, owned sites first, then halo. A domain's ranks split the owned
!  sites as the replicated mode splits all of them.
   subroutine legacy_build(dom, nl, comm, state, params, model, loop, time)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(inout) :: params
      type(model_t), intent(in) :: model
      type(loop_t), intent(inout) :: loop
      type(times_t), intent(inout) :: time
      integer, allocatable :: pairs(:)
      integer :: n_global
      integer :: n_owned
      integer :: n_ranks
      integer :: r

      call time_start(time%neigh)
      if (comm%rank == 0) then
!        The replicated build's own supercell check, for a box that has shrunk.
         call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                       model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                       state%positions, params%do_md, state%velocities, params%masses_types, state%masses, state%xyz_species, &
                       state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                       state%b_box, state%c_box, &
                       state%n_sites, .true., state%fix_atom, params%t_beg, params%write_array_property(6), &
                       .false., params%randomize_velocities)
         n_global = state%n_sites
      end if
      call comm_bcast(comm, n_global)
      call comm_bcast(comm, state%indices, 3)
      call comm_bcast(comm, state%a_box, 3)
      call comm_bcast(comm, state%b_box, 3)
      call comm_bcast(comm, state%c_box, 3)
      if (any(state%indices /= 1)) then
         if (comm%rank == 0) write (*, '(A)') 'ERROR: dd_type legacy needs a cell at least 2 rcut_max across; '// &
            'this one needs a supercell'
         call turbogap_abort()
      end if
      call legacy_check_widths(dom%legacy, state%a_box, state%b_box, state%c_box, model%rcut_max, comm%rank)

      if (comm%rank == 0) then
         dom%g_n_sites = state%n_sites
         call move_alloc(state%positions, dom%g_positions)
         call move_alloc(state%species, dom%g_species)
         call move_alloc(state%xyz_species, dom%g_xyz_species)
      else if (.not. allocated(dom%g_positions)) then
         allocate (dom%g_positions(1:3, 0), dom%g_species(0), dom%g_xyz_species(0))
      end if
      call legacy_distribute(dom%legacy, comm, n_global, dom%g_positions, dom%g_species, dom%g_xyz_species, &
                             state%a_box, state%b_box, state%c_box, model%rcut_max, &
                             state%positions, state%species, state%xyz_species)
      state%n_sites = dom%legacy%n_local
      state%species_supercell = state%species
      state%xyz_species_supercell = state%xyz_species

      n_owned = dom%legacy%n_owned
      n_ranks = dom%legacy%local%size
      r = dom%legacy%local%rank
      if (r < mod(n_owned, n_ranks)) then
         dom%i_beg = 1 + r*(n_owned/n_ranks + 1)
         dom%i_end = dom%i_beg + n_owned/n_ranks
      else
         dom%i_beg = 1 + mod(n_owned, n_ranks)*(n_owned/n_ranks + 1) + (r - mod(n_owned, n_ranks))*(n_owned/n_ranks)
         dom%i_end = dom%i_beg + n_owned/n_ranks - 1
      end if
      if (allocated(dom%do_list)) deallocate (dom%do_list)
      allocate (dom%do_list(1:state%n_sites))
      dom%do_list = .false.
      dom%do_list(dom%i_beg:dom%i_end) = .true.

      nl%rebuild_neighbors_list = .true.
      if (allocated(nl%rjs)) deallocate (nl%rjs, nl%xyz, nl%thetas, nl%phis, nl%neighbor_species)
      if (allocated(nl%neighbors_list)) deallocate (nl%neighbors_list)
      if (allocated(nl%n_neigh)) deallocate (nl%n_neigh)
      if (allocated(nl%n_neigh_local)) deallocate (nl%n_neigh_local)
      call build_neighbors_list(state%positions, state%a_box, state%b_box, state%c_box, params%do_timing, &
                                state%species_supercell, model%rcut_max, nl%n_atom_pairs, nl%rjs, &
                                nl%thetas, nl%phis, nl%xyz, nl%n_neigh_local, nl%neighbors_list, nl%neighbor_species, &
                                state%n_sites, state%indices, .true., dom%do_list, comm%rank)

      allocate (pairs(1:n_ranks))
      call comm_allgather(dom%legacy%local, nl%n_atom_pairs, pairs)
      dom%n_atom_pairs_by_rank = 0
      dom%n_atom_pairs_by_rank(comm%rank + 1) = nl%n_atom_pairs
      nl%n_atom_pairs_total = sum(pairs)
      nl%n_atom_pairs = nl%n_atom_pairs_total
      allocate (nl%n_neigh(1:state%n_sites))
      call comm_sum_to_root(dom%legacy%local, nl%n_neigh_local, nl%n_neigh, state%n_sites)
      call comm_bcast(dom%legacy%local, nl%n_neigh, state%n_sites)
      dom%j_beg = 1
      dom%j_end = pairs(r + 1)

!     Only van der Waals reads these, and it is refused.
      if (allocated(dom%site_in_rank)) deallocate (dom%site_in_rank, dom%this_site_in_rank)
      allocate (dom%site_in_rank(1:state%n_sites))
      allocate (dom%this_site_in_rank(1:state%n_sites))
      dom%site_in_rank = 0
      dom%this_site_in_rank = 0
      call time_end(time%neigh)
   end subroutine legacy_build

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

!     Legacy: summed inside the domain only; rank 0 gets them in the gather.
      if (dom%dd_type == DD_LEGACY) then
         if (any_has_local_properties(model%soap_turbo_hypers)) then
            call time_start(time%mpi)
            call comm_sum_to_root(dom%legacy%local, res%local_properties, res%this_local_properties, &
                                  state%n_sites*params%n_local_properties)
            if (dom%legacy%local%rank == 0) res%local_properties = res%this_local_properties
            call time_end(time%mpi)
         end if
         return
      end if
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
      if (dom%dd_type == DD_LEGACY) then
         call comm_sum_to_root(dom%legacy%local, res%energies, res%this_energies, state%n_sites)
         if (dom%legacy%local%rank == 0) res%energies = res%this_energies
         call time_end(time%mpi_ef)
         return
      end if
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
      type(state_t), intent(inout) :: state
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

      if (dom%dd_type == DD_LEGACY) then
         call legacy_complete(dom, comm, res, state, params, model, time)
         return
      end if
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

!  Legacy: what landed on halo sites goes home, the owned sites' results come
!  to rank 0 in global order, and rank 0 gets its whole structure back with
!  result arrays to match. The four families are the only ones legacy_refuse
!  lets through.
   subroutine legacy_complete(dom, comm, res, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), intent(inout) :: res
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time
      real(dp), allocatable :: part(:, :)
      real(dp), allocatable :: total(:, :)
      real(dp), allocatable :: owned(:, :)
      real(dp), allocatable :: out(:, :)
      real(dp) :: vir_part(1:9, 1:4)
      real(dp) :: vir_dom(1:9, 1:4)
      real(dp) :: vir_all(1:9, 1:4)
      logical :: fam_on(1:4)
      integer :: n_lp
      integer :: per_fam
      integer :: nv
      integer :: n
      integer :: k
      integer :: col
      integer :: i

      call time_start(time%mpi_ef)
      if (allocated(res%this_energies_vdw) .or. allocated(res%this_energies_estat) .or. &
          allocated(res%this_energies_lp) .or. allocated(res%this_energies_pdf) .or. &
          allocated(res%this_energies_sf) .or. allocated(res%this_energies_xrd) .or. &
          allocated(res%this_energies_nd)) then
         if (comm%rank == 0) write (*, '(A)') 'ERROR: dd_type legacy met a term it cannot reduce'
         call turbogap_abort()
      end if
      fam_on = [model%n_soap_turbo > 0, model%n_distance_2b > 0, model%n_core_pot > 0, model%n_angle_3b > 0]
      per_fam = 1
      if (params%do_forces) per_fam = 4
      nv = per_fam*4
      n_lp = 0
      if (any_has_local_properties(model%soap_turbo_hypers)) n_lp = params%n_local_properties
      n = state%n_sites

      allocate (part(1:nv, 1:max(n, 1)))
      allocate (total(1:nv, 1:max(n, 1)))
      part = 0.d0
      total = 0.d0
      vir_part = 0.d0
      do k = 1, 4
         if (.not. fam_on(k)) cycle
         col = (k - 1)*per_fam + 1
         select case (k)
         case (1)
            part(col, 1:n) = res%energies_soap(1:n)
            if (params%do_forces) part(col + 1:col + 3, 1:n) = res%forces_soap(1:3, 1:n)
            if (params%do_forces) vir_part(1:9, k) = reshape(res%virial_soap, [9])
         case (2)
            part(col, 1:n) = res%energies_2b(1:n)
            if (params%do_forces) part(col + 1:col + 3, 1:n) = res%forces_2b(1:3, 1:n)
            if (params%do_forces) vir_part(1:9, k) = reshape(res%virial_2b, [9])
         case (3)
            part(col, 1:n) = res%energies_core_pot(1:n)
            if (params%do_forces) part(col + 1:col + 3, 1:n) = res%forces_core_pot(1:3, 1:n)
            if (params%do_forces) vir_part(1:9, k) = reshape(res%virial_core_pot, [9])
         case (4)
            part(col, 1:n) = res%energies_3b(1:n)
            if (params%do_forces) part(col + 1:col + 3, 1:n) = res%forces_3b(1:3, 1:n)
            if (params%do_forces) vir_part(1:9, k) = reshape(res%virial_3b, [9])
         end select
      end do

      call comm_sum_to_root(dom%legacy%local, part, total, nv*n)
      call comm_sum_to_root(dom%legacy%local, vir_part, vir_dom, 36)

      if (dom%legacy%local%rank == 0) then
         call legacy_reverse(dom%legacy, total)
         call comm_sum_to_root(dom%legacy%roots, vir_dom, vir_all, 36)
         allocate (owned(1:1 + nv + n_lp, 1:max(dom%legacy%n_owned, 1)))
         do i = 1, dom%legacy%n_owned
            owned(1, i) = res%energies(i)
            owned(2:nv + 1, i) = total(1:nv, i)
            if (n_lp > 0) owned(nv + 2:nv + 1 + n_lp, i) = res%local_properties(i, 1:n_lp)
         end do
      else
         allocate (owned(1:1 + nv + n_lp, 1))
      end if

      if (comm%rank == 0) then
         allocate (out(1:1 + nv + n_lp, 1:dom%g_n_sites))
      else
         allocate (out(1:1 + nv + n_lp, 1))
      end if
      call legacy_gather(dom%legacy, comm, owned, out)

      if (comm%rank == 0) then
         state%n_sites = dom%g_n_sites
         call move_alloc(dom%g_positions, state%positions)
         call move_alloc(dom%g_species, state%species)
         call move_alloc(dom%g_xyz_species, state%xyz_species)
         state%species_supercell = state%species
         state%xyz_species_supercell = state%xyz_species
         n = state%n_sites
         call results_resize_sites(res, n)
         res%energies(1:n) = out(1, 1:n)
         do k = 1, 4
            if (.not. fam_on(k)) cycle
            col = 1 + (k - 1)*per_fam + 1
            select case (k)
            case (1)
               res%energies_soap(1:n) = out(col, 1:n)
               if (params%do_forces) res%forces_soap(1:3, 1:n) = out(col + 1:col + 3, 1:n)
               if (params%do_forces) res%virial_soap = reshape(vir_all(1:9, k), [3, 3])
            case (2)
               res%energies_2b(1:n) = out(col, 1:n)
               if (params%do_forces) res%forces_2b(1:3, 1:n) = out(col + 1:col + 3, 1:n)
               if (params%do_forces) res%virial_2b = reshape(vir_all(1:9, k), [3, 3])
            case (3)
               res%energies_core_pot(1:n) = out(col, 1:n)
               if (params%do_forces) res%forces_core_pot(1:3, 1:n) = out(col + 1:col + 3, 1:n)
               if (params%do_forces) res%virial_core_pot = reshape(vir_all(1:9, k), [3, 3])
            case (4)
               res%energies_3b(1:n) = out(col, 1:n)
               if (params%do_forces) res%forces_3b(1:3, 1:n) = out(col + 1:col + 3, 1:n)
               if (params%do_forces) res%virial_3b = reshape(vir_all(1:9, k), [3, 3])
            end select
         end do
         if (n_lp > 0) then
            do i = 1, n_lp
               res%local_properties(1:n, i) = out(nv + 1 + i, 1:n)
            end do
         end if
      end if
      call time_end(time%mpi_ef)
   end subroutine legacy_complete

   subroutine domain_free(dom)
      type(domain_t), intent(inout) :: dom

      deallocate (dom%n_atom_pairs_by_rank)
      if (allocated(dom%do_list)) deallocate (dom%do_list)
      if (dom%dd_type == DD_LEGACY) call legacy_free(dom%legacy)
      if (allocated(dom%g_positions)) deallocate (dom%g_positions, dom%g_species, dom%g_xyz_species)
   end subroutine domain_free

end module turbogap_domain
