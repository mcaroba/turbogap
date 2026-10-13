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
   use turbogap_comm, only: comm_t, comm_bcast, comm_allgather, comm_sum_to_root, comm_with_mpi, comm_gather
   use types, only: input_parameters, any_has_local_properties, contribution_ref
   use timing, only: times_t, time_start, time_end
   use turbogap_structure, only: state_t
   use turbogap_setup, only: model_t
   use turbogap_loop, only: loop_t
   use turbogap_results, only: results_t, results_resize_sites, results_bcast_sites
   use domain_legacy, only: legacy_dd_t, legacy_init, legacy_free, legacy_check_widths, legacy_distribute, &
                            legacy_reverse, legacy_gather, legacy_refresh, legacy_auto_grid, legacy_make_grid, &
                            legacy_forward, legacy_migrate, legacy_rebuild_owned, REC_DYN
   use error, only: turbogap_abort
   use printing, only: print_parameters
   use read_files, only: read_xyz
   use neighbors, only: build_neighbors_list

   implicit none

   private
   public :: domain_init
   public :: domain_configure
   public :: domain_report
   public :: domain_total_dipole
   public :: domain_n_total
   public :: domain_site_values
   public :: domain_site_ids
   public :: domain_exp_view
   public :: domain_localize
   public :: domain_globalize
   public :: domain_integrating

   interface domain_globalize
      module procedure globalize_r2
      module procedure globalize_r3
   end interface domain_globalize
   public :: domain_decomposed
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

!   The whole structure on every rank, split the replicated way, for the terms
!   that are a property of all of it (pair distribution, structure factor, XRD,
!   ND). Built under dd_type legacy only when one of them is on.
   type, public :: exp_view_t
      type(state_t) :: state
      type(neighbors_t) :: nl
      integer :: i_beg = 1
      integer :: i_end = 0
      integer :: j_beg = 1
      integer :: j_end = 0
      logical :: active = .false.
   end type exp_view_t

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
!     Legacy: this rank's local species, kept for the steps between rebuilds;
!     rank 0's state holds the whole structure then, and so does every rank's
!     after a replicate_after step.
      integer, allocatable :: l_species(:)
      character*8, allocatable :: l_xyz_species(:)
      type(exp_view_t) :: view
!     Legacy: give every rank the whole structure and results after the
!     evaluation, for the code that runs on every rank afterwards and expects
!     them (MAD-IR, nested sampling).
      logical :: replicate_after = .false.
!     dd_integrate = domain: each domain integrates its own atoms, which
!     migrate between domains at rebuilds; rank 0 holds the structure only
!     long enough to write it out.
      logical :: integrate = .false.
   end type domain_t

contains

!  Replicated data: every rank holds every atom.
   subroutine domain_init(dom, comm)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm

      allocate (dom%n_atom_pairs_by_rank(1:comm%size))
   end subroutine domain_init

!  Choose the decomposition, once the input has been read.
   subroutine domain_configure(dom, comm, params)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      character(len=64) :: why

      select case (trim(params%dd_type))
      case ("none")
         dom%dd_type = DD_NONE
      case ("legacy")
         why = legacy_unsupported(params)
         if (why /= "") then
            if (comm%rank == 0) write (*, '(A,A,A)') ' dd_type legacy: ', trim(why), &
               ' is not decomposed; this run evaluates replicated, as dd_type none'
            dom%dd_type = DD_NONE
            return
         end if
         dom%dd_type = DD_LEGACY
         dom%replicate_after = params%valid_ir .or. params%do_ir .or. params%do_nested_sampling
         if (trim(params%dd_integrate) == "domain") then
            why = domain_integration_unsupported(params)
            if (why == "") then
               dom%integrate = .true.
            else if (comm%rank == 0) then
               write (*, '(A,A,A)') ' dd_integrate domain: ', trim(why), &
                  ' is integrated on rank 0; the evaluation stays on the domains'
            end if
         else if (trim(params%dd_integrate) /= "rank0") then
            if (comm%rank == 0) write (*, '(A,A,A)') 'ERROR: dd_integrate is "', trim(params%dd_integrate), &
               '"; it must be rank0 or domain'
            call turbogap_abort()
         end if
         call legacy_init(dom%legacy, comm, params%dd_grid, params%dd_grid_affinity, params%dd_balance)
      case default
         if (comm%rank == 0) write (*, '(A,A,A)') 'ERROR: dd_type is "', trim(params%dd_type), &
            '"; it must be none or legacy'
         call turbogap_abort()
      end select
   end subroutine domain_configure

!  The first feature of the run that the domains do not compute yet, or blank.
!  Such a run evaluates replicated instead; docs/LEGACY_DD.md lists them.
   function legacy_unsupported(params) result(why)
      type(input_parameters), intent(in) :: params
      character(len=64) :: why

      why = ""
      if (.not. comm_with_mpi) why = "a build without MPI"
      if (.not. params%do_prediction) why = "a run without do_prediction"
      if (trim(params%vdw_type) == "mbd" .or. trim(params%vdw_type) == "ts+mbd") why = "MBD van der Waals"
   end function legacy_unsupported

!  The first part of the run that the domain integrator does not do yet, or
!  blank. Such a run integrates on rank 0, which is always correct.
   function domain_integration_unsupported(params) result(why)
      type(input_parameters), intent(in) :: params
      character(len=64) :: why

      why = ""
      if (.not. params%do_md) why = "a run without MD"
      if (params%optimize /= "vv") why = "optimize = "//trim(params%optimize)
      if (params%thermostat /= "none" .and. params%thermostat /= "berendsen" .and. params%thermostat /= "bussi") &
         why = "thermostat = "//trim(params%thermostat)
      if (params%barostat /= "none" .or. params%scale_box) why = "a barostat or box scaling"
      if (params%variable_time_step .or. params%adaptive_time) why = "a variable time step"
      if (params%electronic_stopping .or. params%nonadiabatic_processes) why = "electronic stopping"
      if (params%do_mc .or. params%do_nested_sampling) why = "MC or nested sampling"
      if (params%valid_ir .or. params%do_ir) why = "MAD-IR"
      if (params%do_pair_distribution .or. params%do_structure_factor .or. params%do_xrd .or. params%do_nd) &
         why = "a whole-structure observable"
   end function domain_integration_unsupported

!  Is this run integrating on the domains?
   pure logical function domain_integrating(dom) result(yes)
      type(domain_t), intent(in) :: dom

      yes = (dom%dd_type == DD_LEGACY .and. dom%integrate)
   end function domain_integrating

!  Give every rank the structure rank 0 holds. Replicated data: the whole of it,
!  sized from rank 0's arrays. Sending only what an MC move changed was tried,
!  keyed by a broadcast flag array, and segfaulted; hence everything, every time.
   subroutine domain_sync_state(dom, comm, state, params, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(times_t), intent(inout) :: time

!     Legacy: domain_build hands each rank its domain.
      if (dom%dd_type == DD_LEGACY) return
      call replicate_state(comm, state, params, time)
   end subroutine domain_sync_state

!  Rank 0's structure to every rank, as the replicated mode keeps it.
   subroutine replicate_state(comm, state, params, time)
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
   end subroutine replicate_state

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

!     Legacy: whether to rebuild, decided by rank 0's skin test as in the
!     replicated mode. The positions go out in domain_build.
      if (dom%dd_type == DD_LEGACY) then
         call comm_bcast(comm, nl%rebuild_neighbors_list)
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
!     Legacy: rank 0 holds i-PI's configuration, and domain_build sends it out.
      if (dom%dd_type == DD_LEGACY) return
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
         if (dom%dd_type == DD_LEGACY) return
!        The cell cannot be decomposed: this step and the rest go replicated.
         call domain_sync_state(dom, comm, state, params, time)
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
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
      integer :: grid(1:3)
      integer :: n_global

      call time_start(time%neigh)
      dom%view%active = .false.
      call comm_bcast(comm, nl%rebuild_neighbors_list)
      if (dom%integrate .and. dom%legacy%n_builds > 0) then
         call legacy_domain_build(dom, nl, comm, state, params, model, time)
         call time_end(time%neigh)
         return
      end if
      if (params%do_md .and. dom%legacy%n_builds > 0 .and. .not. nl%rebuild_neighbors_list) then
         call legacy_reuse(dom, nl, comm, state, params, model, time)
         call time_end(time%neigh)
         return
      end if
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
         call legacy_fall_back(dom, nl, comm, 'the cell needs a supercell to hold one cutoff sphere')
         call time_end(time%neigh)
         return
      end if
      if (.not. dom%legacy%grid_ready) then
         grid = legacy_auto_grid(dom%legacy%n_domains, state%a_box, state%b_box, state%c_box, model%rcut_max)
         if (all(grid == 0)) then
            call legacy_fall_back(dom, nl, comm, 'no grid of that many domains keeps each one rcut_max wide')
            call time_end(time%neigh)
            return
         end if
         call legacy_make_grid(dom%legacy, comm, grid)
         if (comm%rank == 0) call print_parameters("dd_grid (chosen)", grid)
      end if
      if (.not. legacy_check_widths(dom%legacy, state%a_box, state%b_box, state%c_box, model%rcut_max)) then
         call legacy_fall_back(dom, nl, comm, 'a domain of dd_grid is narrower than rcut_max')
         call time_end(time%neigh)
         return
      end if

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
                             state%positions, state%species, state%xyz_species, &
                             time%dd_scatter, time%dd_halo, time%dd_share)
      state%n_sites = dom%legacy%n_local
      state%species_supercell = state%species
      state%xyz_species_supercell = state%xyz_species
      dom%l_species = state%species
      dom%l_xyz_species = state%xyz_species

      if (dom%integrate) call legacy_first_dynamics(dom, comm, state)
      call legacy_lists(dom, nl, comm, state, params, model)
      call time_end(time%neigh)
   end subroutine legacy_build

!  The owned sites split over the domain's ranks, and their neighbour lists
!  built from scratch.
   subroutine legacy_lists(dom, nl, comm, state, params, model)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      integer, allocatable :: pairs(:)
      integer :: n_owned
      integer :: n_ranks
      integer :: r

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
   end subroutine legacy_lists

!  dd_integrate = domain, first build: every local site gets its velocity, mass
!  and fixed coordinates from rank 0's whole structure, by global id. The
!  previous positions and forces are set by the first velocity-Verlet step.
!  Rank 0 keeps nothing aside after this; the domains hold the structure.
   subroutine legacy_first_dynamics(dom, comm, state)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      real(dp), allocatable :: vel(:, :)
      real(dp), allocatable :: mass(:)
      logical, allocatable :: fix(:, :)
      integer :: n
      integer :: n_local

      n = dom%legacy%n_global
      n_local = dom%legacy%n_local
      allocate (vel(1:3, 1:n), mass(1:n), fix(1:3, 1:n))
      if (comm%rank == 0) then
         vel = state%velocities(1:3, 1:n)
         mass = state%masses(1:n)
         fix = state%fix_atom(1:3, 1:n)
      end if
      call comm_bcast(comm, vel, 3*n)
      call comm_bcast(comm, mass, n)
      call comm_bcast(comm, fix, 3*n)
      call set_local_dynamics(state, n_local)
      state%velocities(1:3, 1:n_local) = vel(1:3, dom%legacy%ids)
      state%masses(1:n_local) = mass(dom%legacy%ids)
      state%fix_atom(1:3, 1:n_local) = fix(1:3, dom%legacy%ids)
      state%positions_prev(1:3, 1:n_local) = state%positions(1:3, 1:n_local)
      state%forces_prev = 0.d0
      state%positions_diff = 0.d0
      if (allocated(dom%g_positions)) deallocate (dom%g_positions, dom%g_species, dom%g_xyz_species)
   end subroutine legacy_first_dynamics

!  The per-site dynamic arrays sized for n local sites.
   subroutine set_local_dynamics(state, n)
      type(state_t), intent(inout) :: state
      integer, intent(in) :: n

      if (allocated(state%velocities)) deallocate (state%velocities)
      if (allocated(state%masses)) deallocate (state%masses)
      if (allocated(state%fix_atom)) deallocate (state%fix_atom)
      if (allocated(state%positions_prev)) deallocate (state%positions_prev)
      if (allocated(state%forces_prev)) deallocate (state%forces_prev)
      if (allocated(state%positions_diff)) deallocate (state%positions_diff)
      allocate (state%velocities(1:3, 1:max(n, 1)), state%masses(1:max(n, 1)), state%fix_atom(1:3, 1:max(n, 1)))
      allocate (state%positions_prev(1:3, 1:max(n, 1)), state%forces_prev(1:3, 1:max(n, 1)))
      allocate (state%positions_diff(1:3, 1:max(n, 1)))
      state%velocities = 0.d0
      state%masses = 0.d0
      state%fix_atom = .false.
      state%positions_prev = 0.d0
      state%forces_prev = 0.d0
      state%positions_diff = 0.d0
   end subroutine set_local_dynamics

!  dd_integrate = domain, after the first build. Between rebuilds the owned
!  positions, moved by the domain's own integrator, are copied to their halo
!  copies and the lists re-measured. At a rebuild the owned atoms migrate,
!  with everything the integrator keeps for them, the halos are rebuilt and
!  the lists built afresh.
   subroutine legacy_domain_build(dom, nl, comm, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time
      real(dp), allocatable :: t(:, :)
      real(dp), allocatable :: dyn(:, :)
      integer :: n_owned
      integer :: n_local

      n_local = dom%legacy%n_local
      if (.not. nl%rebuild_neighbors_list) then
         allocate (t(1:3, 1:max(n_local, 1)))
         t = 0.d0
         t(1:3, 1:n_local) = state%positions(1:3, 1:n_local)
         if (dom%legacy%local%rank == 0) then
            call time_start(time%dd_halo, "dd_halo")
            call legacy_forward(dom%legacy, t)
            call time_end(time%dd_halo, "dd_halo")
         end if
         call comm_bcast(dom%legacy%local, t, 3*n_local)
         state%positions(1:3, 1:n_local) = t(1:3, 1:n_local)
         dom%legacy%n_refreshes = dom%legacy%n_refreshes + 1
         call build_neighbors_list(state%positions, state%a_box, state%b_box, state%c_box, params%do_timing, &
                                   state%species_supercell, model%rcut_max, nl%n_atom_pairs, nl%rjs, &
                                   nl%thetas, nl%phis, nl%xyz, nl%n_neigh_local, nl%neighbors_list, nl%neighbor_species, &
                                   state%n_sites, state%indices, .false., dom%do_list, comm%rank)
         return
      end if

      n_owned = dom%legacy%n_owned
      allocate (dyn(1:REC_DYN, 1:max(n_owned, 1)))
      if (dom%legacy%local%rank == 0) then
         call pack_dynamics(dom, state, dyn)
         call time_start(time%dd_scatter, "dd_scatter")
         call legacy_migrate(dom%legacy, dyn, state%a_box, state%b_box, state%c_box, model%rcut_max)
         call time_end(time%dd_scatter, "dd_scatter")
      end if
      call legacy_rebuild_owned(dom%legacy, dyn, state%a_box, state%b_box, state%c_box, model%rcut_max, &
                                state%positions, state%species, state%xyz_species, time%dd_halo, time%dd_share)
      state%n_sites = dom%legacy%n_local
      state%species_supercell = state%species
      state%xyz_species_supercell = state%xyz_species
      dom%l_species = state%species
      dom%l_xyz_species = state%xyz_species
      call set_local_dynamics(state, dom%legacy%n_local)
      if (dom%legacy%local%rank == 0) call unpack_dynamics(dom, dyn, state)
      call legacy_lists(dom, nl, comm, state, params, model)
   end subroutine legacy_domain_build

!  Owned atom i as a migration record: the REC_ATOM fields, then what the
!  integrator keeps for it.
   subroutine pack_dynamics(dom, state, dyn)
      type(domain_t), intent(in) :: dom
      type(state_t), intent(in) :: state
      real(dp), intent(inout) :: dyn(:, :)
      integer :: i

      do i = 1, dom%legacy%n_owned
         dyn(1, i) = dble(dom%legacy%ids(i))
         dyn(2, i) = dble(state%species(i))
         dyn(3:5, i) = state%positions(1:3, i)
         dyn(6, i) = transfer(state%xyz_species(i), 0.d0)
         dyn(7:9, i) = state%velocities(1:3, i)
         dyn(10, i) = state%masses(i)
         dyn(11:13, i) = merge(1.d0, 0.d0, state%fix_atom(1:3, i))
         dyn(14:16, i) = state%positions_prev(1:3, i)
         dyn(17:19, i) = state%forces_prev(1:3, i)
         dyn(20:22, i) = state%positions_diff(1:3, i)
      end do
   end subroutine pack_dynamics

   subroutine unpack_dynamics(dom, dyn, state)
      type(domain_t), intent(in) :: dom
      real(dp), intent(in) :: dyn(:, :)
      type(state_t), intent(inout) :: state
      integer :: i

      do i = 1, dom%legacy%n_owned
         state%velocities(1:3, i) = dyn(7:9, i)
         state%masses(i) = dyn(10, i)
         state%fix_atom(1:3, i) = dyn(11:13, i) > 0.5d0
         state%positions_prev(1:3, i) = dyn(14:16, i)
         state%forces_prev(1:3, i) = dyn(17:19, i)
         state%positions_diff(1:3, i) = dyn(20:22, i)
      end do
   end subroutine unpack_dynamics

!  The cell cannot be split as asked: say why, and evaluate replicated from
!  here on. Rank 0 still holds the whole structure, so nothing is lost.
   subroutine legacy_fall_back(dom, nl, comm, why)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      character(len=*), intent(in) :: why

      if (comm%rank == 0) write (*, '(A,A,A)') ' dd_type legacy: ', why, '; evaluating replicated from here on'
      dom%dd_type = DD_NONE
      dom%integrate = .false.
      nl%rebuild_neighbors_list = .true.
      if (allocated(nl%rjs)) deallocate (nl%rjs, nl%xyz, nl%thetas, nl%phis, nl%neighbor_species)
      if (allocated(nl%neighbors_list)) deallocate (nl%neighbors_list)
      if (allocated(nl%n_neigh)) deallocate (nl%n_neigh)
      if (allocated(nl%n_neigh_local)) deallocate (nl%n_neigh_local)
   end subroutine legacy_fall_back

!  Legacy, between rebuilds: the domains, halos and neighbour lists stay; the
!  positions are refreshed and the pair geometry recomputed, as the replicated
!  mode does when the skin test passes.
   subroutine legacy_reuse(dom, nl, comm, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(neighbors_t), intent(inout) :: nl
      type(comm_t), intent(in) :: comm
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time

      call comm_bcast(comm, state%a_box, 3)
      call comm_bcast(comm, state%b_box, 3)
      call comm_bcast(comm, state%c_box, 3)
      if (comm%rank == 0) then
         dom%g_n_sites = state%n_sites
         call move_alloc(state%positions, dom%g_positions)
         call move_alloc(state%species, dom%g_species)
         call move_alloc(state%xyz_species, dom%g_xyz_species)
      end if
      state%species = dom%l_species
      state%xyz_species = dom%l_xyz_species
      state%species_supercell = state%species
      state%xyz_species_supercell = state%xyz_species
      call legacy_refresh(dom%legacy, comm, dom%g_positions, state%positions, time%dd_scatter, time%dd_share)
      state%n_sites = dom%legacy%n_local
      dom%do_list = .false.
      dom%do_list(dom%i_beg:dom%i_end) = .true.
      call build_neighbors_list(state%positions, state%a_box, state%b_box, state%c_box, params%do_timing, &
                                state%species_supercell, model%rcut_max, nl%n_atom_pairs, nl%rjs, &
                                nl%thetas, nl%phis, nl%xyz, nl%n_neigh_local, nl%neighbors_list, nl%neighbor_species, &
                                state%n_sites, state%indices, .false., dom%do_list, comm%rank)
   end subroutine legacy_reuse

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

!     Legacy: summed inside the domain, then copied from owners to their halo
!     copies, since vdW, electrostatics and XPS read a neighbour's property.
!     Rank 0 gets the owned values in the gather.
      if (dom%dd_type == DD_LEGACY) then
         call time_start(time%mpi)
         if (any_has_local_properties(model%soap_turbo_hypers)) then
            call comm_sum_to_root(dom%legacy%local, res%local_properties, res%this_local_properties, &
                                  state%n_sites*params%n_local_properties)
            if (dom%legacy%local%rank == 0) then
               res%local_properties = res%this_local_properties
               call legacy_share_sites(dom, res%local_properties)
            end if
            call comm_bcast(dom%legacy%local, res%local_properties, state%n_sites*params%n_local_properties)
         end if
         if (params%do_dipole) then
            call comm_sum_to_root(dom%legacy%local, res%local_dipoles, res%this_local_dipoles, 3*state%n_sites)
            if (dom%legacy%local%rank == 0) res%local_dipoles = res%this_local_dipoles
            call comm_sum_to_root(dom%legacy%local, res%energies_dipole, res%this_energies_dipole, state%n_sites)
            if (dom%legacy%local%rank == 0) res%energies_dipole = res%this_energies_dipole
         end if
         call time_end(time%mpi)
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
      integer :: family(1:N_CONTRIB)
      integer :: n_active
      integer :: i_contrib
      integer :: counter2

      if (dom%dd_type == DD_LEGACY) then
         call legacy_complete(dom, comm, res, state, params, model, time)
         return
      end if
      call time_start(time%mpi_ef)
      call contributions_list(res, params, model, contrib_on, contrib, n_active, family)

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

      call contributions_release(res, params, contrib_on)

      !       Clean up
      deallocate (all_energies, all_this_energies)
      if (params%do_forces) then
         deallocate (all_forces, all_this_forces, all_virial, all_this_virial)
      end if

      call time_end(time%mpi_ef)
   end subroutine domain_complete_contributions

!  The families active this step, as source and destination arrays. Built from
!  one evaluation of the eleven predicates, so the pack and unpack walks cannot
!  disagree about which slot belongs to which family.
   subroutine contributions_list(res, params, model, contrib_on, contrib, n_active, family)
      type(results_t), target, intent(inout) :: res
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      logical, intent(out) :: contrib_on(1:N_CONTRIB)
      type(contribution_ref), intent(out) :: contrib(1:N_CONTRIB)
      integer, intent(out) :: n_active
!     Which family each slot holds.
      integer, intent(out) :: family(1:N_CONTRIB)

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
         family(n_active) = C_SOAP
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
         family(n_active) = C_VDW
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
         family(n_active) = C_ESTAT
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
         family(n_active) = C_LP
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
         family(n_active) = C_PDF
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
         family(n_active) = C_SF
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
         family(n_active) = C_XRD
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
         family(n_active) = C_ND
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
         family(n_active) = C_2B
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
         family(n_active) = C_CP
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
         family(n_active) = C_3B
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

   end subroutine contributions_list

!  The this_ arrays, once their contents are in the un-prefixed ones.
   subroutine contributions_release(res, params, contrib_on)
      type(results_t), intent(inout) :: res
      type(input_parameters), intent(in) :: params
      logical, intent(in) :: contrib_on(1:N_CONTRIB)

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

   end subroutine contributions_release

!  Legacy, roots: copy per-site values, sites on the first index, from owners
!  to halo copies.
   subroutine legacy_share_sites(dom, x)
      type(domain_t), intent(in) :: dom
      real(dp), intent(inout) :: x(:, :)
      real(dp), allocatable :: t(:, :)

      t = transpose(x)
      call legacy_forward(dom%legacy, t)
      x = transpose(t)
   end subroutine legacy_share_sites

!  Atoms in the whole system, whatever this rank holds.
   pure integer function domain_n_total(dom, state) result(n)
      type(domain_t), intent(in) :: dom
      type(state_t), intent(in) :: state

      n = state%n_sites
      if (dom%dd_type == DD_LEGACY) n = dom%legacy%n_global
   end function domain_n_total

!  A per-site value of every atom in the system, in global order, on every
!  rank. Replicated: the sites here already are the system. Legacy: the owners
!  send theirs to rank 0, which broadcasts the whole array.
   subroutine domain_site_values(dom, comm, state, x, all)
      type(domain_t), intent(in) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(in) :: state
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: all(:)
      real(dp), allocatable :: owned(:, :)
      real(dp), allocatable :: out(:, :)
      integer :: n

      if (dom%dd_type /= DD_LEGACY) then
         all = x(1:state%n_sites)
         return
      end if
      n = dom%legacy%n_global
      allocate (owned(1:1, 1:max(dom%legacy%n_owned, 1)))
      owned(1, 1:dom%legacy%n_owned) = x(1:dom%legacy%n_owned)
      if (comm%rank == 0) then
         allocate (out(1:1, 1:n))
      else
         allocate (out(1:1, 1))
      end if
      call legacy_gather(dom%legacy, comm, owned, out)
      allocate (all(1:n))
      if (comm%rank == 0) all = out(1, 1:n)
      call comm_bcast(comm, all, n)
   end subroutine domain_site_values

!  The global index of each site this rank holds.
   subroutine domain_site_ids(dom, state, ids)
      type(domain_t), intent(in) :: dom
      type(state_t), intent(in) :: state
      integer, allocatable, intent(out) :: ids(:)
      integer :: i

      if (dom%dd_type == DD_LEGACY) then
         ids = dom%legacy%ids
      else
         allocate (ids(1:state%n_sites))
         ids = [(i, i=1, state%n_sites)]
      end if
   end subroutine domain_site_ids

!  Legacy: the whole structure on every rank and a replicated neighbour list
!  for it, for the experimental observables. Rank 0 holds the structure aside
!  during the evaluation; everyone else gets it from there.
   subroutine domain_exp_view(dom, comm, state, params, model)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(state_t), intent(in) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      logical, allocatable :: do_list(:)
      integer, allocatable :: pairs(:)
      integer :: n

      associate (v => dom%view)
         v%active = .true.
         n = dom%legacy%n_global
         v%state%n_sites = n
         v%state%a_box = state%a_box
         v%state%b_box = state%b_box
         v%state%c_box = state%c_box
         v%state%indices = state%indices
         v%state%v_uc = state%v_uc
         if (allocated(v%state%positions)) deallocate (v%state%positions, v%state%species, v%state%species_supercell)
         allocate (v%state%positions(1:3, 1:n), v%state%species(1:n))
         if (comm%rank == 0) then
            v%state%positions = dom%g_positions(1:3, 1:n)
            v%state%species = dom%g_species(1:n)
         end if
         call comm_bcast(comm, v%state%positions, 3*n)
         call comm_bcast(comm, v%state%species, n)
         v%state%species_supercell = v%state%species

         if (comm%rank < mod(n, comm%size)) then
            v%i_beg = 1 + comm%rank*(n/comm%size + 1)
            v%i_end = v%i_beg + n/comm%size
         else
            v%i_beg = 1 + mod(n, comm%size)*(n/comm%size + 1) + (comm%rank - mod(n, comm%size))*(n/comm%size)
            v%i_end = v%i_beg + n/comm%size - 1
         end if
         allocate (do_list(1:n))
         do_list = .false.
         do_list(v%i_beg:v%i_end) = .true.

         if (allocated(v%nl%rjs)) deallocate (v%nl%rjs, v%nl%xyz, v%nl%thetas, v%nl%phis, v%nl%neighbor_species)
         if (allocated(v%nl%neighbors_list)) deallocate (v%nl%neighbors_list)
         if (allocated(v%nl%n_neigh)) deallocate (v%nl%n_neigh)
         if (allocated(v%nl%n_neigh_local)) deallocate (v%nl%n_neigh_local)
         call build_neighbors_list(v%state%positions, v%state%a_box, v%state%b_box, v%state%c_box, params%do_timing, &
                                   v%state%species_supercell, model%rcut_max, v%nl%n_atom_pairs, v%nl%rjs, &
                                   v%nl%thetas, v%nl%phis, v%nl%xyz, v%nl%n_neigh_local, v%nl%neighbors_list, &
                                   v%nl%neighbor_species, n, v%state%indices, .true., do_list, comm%rank)
         allocate (pairs(1:comm%size))
         call comm_allgather(comm, v%nl%n_atom_pairs, pairs)
         v%j_beg = 1
         v%j_end = pairs(comm%rank + 1)
         v%nl%n_atom_pairs_total = sum(pairs)
         v%nl%n_atom_pairs = v%nl%n_atom_pairs_total
         allocate (v%nl%n_neigh(1:n))
         call comm_sum_to_root(comm, v%nl%n_neigh_local, v%nl%n_neigh, n)
         call comm_bcast(comm, v%nl%n_neigh, n)
      end associate
   end subroutine domain_exp_view

!  Legacy: reorder a per-site array, sites on the second index, from global
!  order into this rank's local order, for code that reads it by local site
!  during the evaluation. x must have at least n_local columns.
   subroutine domain_localize(dom, x)
      type(domain_t), intent(in) :: dom
      real(dp), intent(inout) :: x(:, :)
      real(dp), allocatable :: t(:, :)

      if (dom%dd_type /= DD_LEGACY) return
      t = x(:, dom%legacy%ids)
      x(:, 1:dom%legacy%n_local) = t
   end subroutine domain_localize

!  Legacy: the reverse, for a per-site sum each rank accumulated by local site
!  during the evaluation. On return rank 0 holds the complete sums in global
!  order and every other rank zeros, so a following all-reduce, as the
!  replicated mode does, completes it everywhere.
   subroutine globalize_r2(dom, comm, x)
      type(domain_t), intent(in) :: dom
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: x(:, :)
      real(dp), allocatable :: part(:, :)
      real(dp), allocatable :: total(:, :)
      real(dp), allocatable :: out(:, :)
      integer :: m
      integer :: n

      if (dom%dd_type /= DD_LEGACY) return
      m = size(x, 1)
      n = dom%legacy%n_local
      part = x(:, 1:max(n, 1))
      if (n == 0) part = 0.d0
      allocate (total(1:m, 1:max(n, 1)))
      total = 0.d0
      call comm_sum_to_root(dom%legacy%local, part, total, m*n)
      if (dom%legacy%local%rank == 0) call legacy_reverse(dom%legacy, total)
      if (comm%rank == 0) then
         allocate (out(1:m, 1:dom%legacy%n_global))
      else
         allocate (out(1:m, 1))
      end if
      call legacy_gather(dom%legacy, comm, total(:, 1:max(dom%legacy%n_owned, 1)), out)
      x = 0.d0
      if (comm%rank == 0) x(:, 1:dom%legacy%n_global) = out
   end subroutine globalize_r2

!  The same for a tensor per site: (a, b, site).
   subroutine globalize_r3(dom, comm, x)
      type(domain_t), intent(in) :: dom
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: x(:, :, :)
      real(dp), allocatable :: flat(:, :)

      if (dom%dd_type /= DD_LEGACY) return
      flat = reshape(x, [size(x, 1)*size(x, 2), size(x, 3)])
      call globalize_r2(dom, comm, flat)
      x = reshape(flat, shape(x))
   end subroutine globalize_r3

!  Does this rank hold only its domain's atoms during the evaluation?
   pure logical function domain_decomposed(dom) result(yes)
      type(domain_t), intent(in) :: dom

      yes = (dom%dd_type == DD_LEGACY)
   end function domain_decomposed

!  The total dipole: a sum over every site, each counted once.
   subroutine domain_total_dipole(dom, comm, res, state)
      type(domain_t), intent(in) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), intent(inout) :: res
      type(state_t), intent(in) :: state
      real(dp) :: mine(1:3)
      real(dp) :: all(1:3)
      integer :: n

      if (dom%dd_type /= DD_LEGACY) then
         res%dipole(1) = sum(res%local_dipoles(1, 1:state%n_sites))
         res%dipole(2) = sum(res%local_dipoles(2, 1:state%n_sites))
         res%dipole(3) = sum(res%local_dipoles(3, 1:state%n_sites))
         return
      end if
      mine = 0.d0
      if (dom%legacy%local%rank == 0) then
         n = dom%legacy%n_owned
         mine(1) = sum(res%local_dipoles(1, 1:n))
         mine(2) = sum(res%local_dipoles(2, 1:n))
         mine(3) = sum(res%local_dipoles(3, 1:n))
      end if
      call comm_sum_to_root(comm, mine, all, 3)
      if (comm%rank == 0) res%dipole = all
      call comm_bcast(comm, res%dipole, 3)
   end subroutine domain_total_dipole

!  Legacy: what landed on halo sites goes home, the owned sites' results come
!  to rank 0 in global order, and rank 0 gets its whole structure back with
!  result arrays to match. Every active family goes the same way, from the
!  same list the replicated reduction uses.
   subroutine legacy_complete(dom, comm, res, state, params, model, time)
      type(domain_t), intent(inout) :: dom
      type(comm_t), intent(in) :: comm
      type(results_t), target, intent(inout) :: res
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(times_t), intent(inout) :: time
      real(dp), allocatable :: part(:, :)
      real(dp), allocatable :: total(:, :)
      real(dp), allocatable :: owned(:, :)
      real(dp), allocatable :: out(:, :)
      real(dp) :: vir_part(1:9, 1:N_CONTRIB)
      real(dp) :: vir_dom(1:9, 1:N_CONTRIB)
      real(dp) :: vir_all(1:9, 1:N_CONTRIB)
      logical :: contrib_on(1:N_CONTRIB)
      type(contribution_ref) :: contrib(1:N_CONTRIB)
      integer :: family(1:N_CONTRIB)
      logical :: on_view(1:N_CONTRIB)
      real(dp), allocatable :: gpart(:, :)
      real(dp), allocatable :: gtot(:, :)
      integer :: n_global
      integer :: n_active
      integer :: n_lp
      integer :: n_dip
      integer :: per_fam
      integer :: nv
      integer :: n_out
      integer :: n
      integer :: k
      integer :: col
      integer :: i

      call time_start(time%mpi_ef)
      call contributions_list(res, params, model, contrib_on, contrib, n_active, family)
      per_fam = 1
      if (params%do_forces) per_fam = 4
      nv = per_fam*n_active
      n_lp = 0
      if (any_has_local_properties(model%soap_turbo_hypers)) n_lp = params%n_local_properties
      n_dip = 0
      if (params%do_dipole) n_dip = 4
      n = state%n_sites

!     The families computed on the whole-structure view are indexed by global
!     atom, not by local site: they skip the halo return and the gather, and are
!     summed over all ranks straight into global order.
      on_view = .false.
      do k = 1, n_active
         on_view(k) = dom%view%active .and. any(family(k) == [C_PDF, C_SF, C_XRD, C_ND])
      end do
      n_global = dom%legacy%n_global
      if (any(on_view)) then
         allocate (gpart(1:max(nv, 1), 1:n_global))
         allocate (gtot(1:max(nv, 1), 1:n_global))
         gpart = 0.d0
         gtot = 0.d0
         do k = 1, n_active
            if (.not. on_view(k)) cycle
            col = (k - 1)*per_fam + 1
            gpart(col, 1:n_global) = contrib(k)%e_src(1:n_global)
            if (contrib(k)%forces) gpart(col + 1:col + 3, 1:n_global) = contrib(k)%f_src(1:3, 1:n_global)
         end do
         call comm_sum_to_root(comm, gpart, gtot, size(gpart))
      end if

!     Families in columns: energy, then the three force components.
      allocate (part(1:max(nv, 1), 1:max(n, 1)))
      allocate (total(1:max(nv, 1), 1:max(n, 1)))
      part = 0.d0
      total = 0.d0
      vir_part = 0.d0
      do k = 1, n_active
         col = (k - 1)*per_fam + 1
         if (contrib(k)%forces) vir_part(1:9, k) = reshape(contrib(k)%v_src, [9])
         if (on_view(k)) cycle
         part(col, 1:n) = contrib(k)%e_src(1:n)
         if (contrib(k)%forces) part(col + 1:col + 3, 1:n) = contrib(k)%f_src(1:3, 1:n)
      end do

      call comm_sum_to_root(dom%legacy%local, part, total, size(part))
      call comm_sum_to_root(dom%legacy%local, vir_part, vir_dom, 9*N_CONTRIB)

      n_out = 1 + nv + n_lp + n_dip
      if (dom%legacy%local%rank == 0) then
         call time_start(time%dd_reverse, "dd_reverse")
         call legacy_reverse(dom%legacy, total)
         call time_end(time%dd_reverse, "dd_reverse")
         call comm_sum_to_root(dom%legacy%roots, vir_dom, vir_all, 9*N_CONTRIB)
         allocate (owned(1:n_out, 1:max(dom%legacy%n_owned, 1)))
         do i = 1, dom%legacy%n_owned
            owned(1, i) = res%energies(i)
            owned(2:nv + 1, i) = total(1:nv, i)
            if (n_lp > 0) owned(nv + 2:nv + 1 + n_lp, i) = res%local_properties(i, 1:n_lp)
            if (n_dip > 0) then
               owned(nv + n_lp + 2:nv + n_lp + 4, i) = res%local_dipoles(1:3, i)
               owned(nv + n_lp + 5, i) = res%energies_dipole(i)
            end if
         end do
      else
         allocate (owned(1:n_out, 1))
      end if

!     dd_integrate = domain: the results stay where the atoms are. The owned
!     sites get their completed sums for the domain's own integrator, and
!     everyone the whole virial; nothing is gathered.
      if (dom%integrate) then
         call comm_bcast(comm, vir_all, 9*N_CONTRIB)
         if (dom%legacy%local%rank == 0) then
            n = dom%legacy%n_owned
            do k = 1, n_active
               col = (k - 1)*per_fam + 1
               contrib(k)%e_dst(1:n) = total(col, 1:n)
               if (contrib(k)%forces) then
                  contrib(k)%f_dst(1:3, 1:n) = total(col + 1:col + 3, 1:n)
                  contrib(k)%v_dst = reshape(vir_all(1:9, k), [3, 3])
               end if
            end do
         end if
         call contributions_release(res, params, contrib_on)
         call time_end(time%mpi_ef)
         return
      end if

      if (comm%rank == 0) then
         allocate (out(1:n_out, 1:dom%g_n_sites))
      else
         allocate (out(1:n_out, 1))
      end if
      call time_start(time%dd_gather, "dd_gather")
      call legacy_gather(dom%legacy, comm, owned, out)
      call time_end(time%dd_gather, "dd_gather")

      if (comm%rank /= 0) call contributions_release(res, params, contrib_on)
      if (comm%rank == 0) then
         state%n_sites = dom%g_n_sites
         call move_alloc(dom%g_positions, state%positions)
         call move_alloc(dom%g_species, state%species)
         call move_alloc(dom%g_xyz_species, state%xyz_species)
         state%species_supercell = state%species
         state%xyz_species_supercell = state%xyz_species
         n = state%n_sites
         call results_resize_sites(res, n)
!        The list again, since the destinations were just reallocated. The
!        resize kept every this_ array allocated, so its predicates still hold.
         call contributions_list(res, params, model, contrib_on, contrib, n_active, family)
         res%energies(1:n) = out(1, 1:n)
         do k = 1, n_active
            col = 1 + (k - 1)*per_fam + 1
            if (on_view(k)) then
               contrib(k)%e_dst(1:n) = gtot(col - 1, 1:n)
               if (contrib(k)%forces) contrib(k)%f_dst(1:3, 1:n) = gtot(col:col + 2, 1:n)
            else
               contrib(k)%e_dst(1:n) = out(col, 1:n)
               if (contrib(k)%forces) contrib(k)%f_dst(1:3, 1:n) = out(col + 1:col + 3, 1:n)
            end if
            if (contrib(k)%forces) contrib(k)%v_dst = reshape(vir_all(1:9, k), [3, 3])
         end do
         do i = 1, n_lp
            res%local_properties(1:n, i) = out(nv + 1 + i, 1:n)
         end do
         if (n_dip > 0) then
            res%local_dipoles(1:3, 1:n) = out(nv + n_lp + 2:nv + n_lp + 4, 1:n)
            res%energies_dipole(1:n) = out(nv + n_lp + 5, 1:n)
         end if
         call contributions_release(res, params, contrib_on)
      end if
      call time_end(time%mpi_ef)
      if (dom%replicate_after) then
         call comm_bcast(comm, n)
         if (comm%rank /= 0) call results_resize_sites(res, n)
         call results_bcast_sites(res, comm)
         call replicate_state(comm, state, params, time)
      end if
   end subroutine legacy_complete

!  Legacy: one line per rank at the end of the run, for load balance. Times
!  are this rank's totals; atoms and pairs are from the last step.
   subroutine domain_report(dom, comm, time)
      type(domain_t), intent(in) :: dom
      type(comm_t), intent(in) :: comm
      type(times_t), intent(in) :: time
      integer, parameter :: N_COL = 13
      real(dp) :: mine(1:N_COL)
      real(dp), allocatable :: all(:)
      integer :: r
      integer :: k

      if (dom%dd_type /= DD_LEGACY) return
      mine = [dble(dom%legacy%local%rank), dble(dom%legacy%n_owned), dble(dom%legacy%n_local), &
              dble(dom%i_end - dom%i_beg + 1), dble(dom%j_end - dom%j_beg + 1), time%gap(3), time%neigh(3), &
              time%dd_scatter(3), time%dd_halo(3), time%dd_share(3), time%dd_reverse(3), time%dd_gather(3), &
              dble(dom%legacy%n_migrated)]
      allocate (all(1:N_COL*comm%size))
      call comm_gather(comm, mine, N_COL, all)
      if (comm%rank /= 0) return
      write (*, '(A,I0,A,I0,A,A)') ' dd_type legacy: ', dom%legacy%n_builds, ' rebuilds, ', dom%legacy%n_refreshes, &
         ' steps reusing the domains', trim(merge(', integrated on the domains', '                           ', &
                                                  dom%integrate))
      if (dom%integrate) write (*, '(A,I0)') ' dd_type legacy: atoms migrated between domains: ', &
         nint(sum(all(13::N_COL)))
      write (*, '(A)') ' dd_type legacy, per rank (seconds over the run; sites and pairs of the last step):'
      write (*, '(A)') '  rank  dom.rank  owned   local  slice     pairs      gap    neigh  scatter     halo    share'// &
         '  reverse   gather'
      do r = 0, comm%size - 1
         k = r*N_COL
         write (*, '(I6,I10,I7,I8,I7,I10,7F9.3)') r, nint(all(k + 1)), nint(all(k + 2)), nint(all(k + 3)), &
            nint(all(k + 4)), nint(all(k + 5)), all(k + 6:k + 12)
      end do
   end subroutine domain_report

   subroutine domain_free(dom)
      type(domain_t), intent(inout) :: dom

      deallocate (dom%n_atom_pairs_by_rank)
      if (allocated(dom%do_list)) deallocate (dom%do_list)
      call legacy_free(dom%legacy)
      if (allocated(dom%g_positions)) deallocate (dom%g_positions, dom%g_species, dom%g_xyz_species)
      if (allocated(dom%l_species)) deallocate (dom%l_species, dom%l_xyz_species)
   end subroutine domain_free

end module turbogap_domain
