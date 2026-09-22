! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_structure.f90, is copyright (c) 2026, Miguel A. Caro
! HND X   and Tigany Zarrouk
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

!  The atomic structure the driver steps: what read_xyz produces and the MD
!  integrator advances. Rank 0 writes it; the other ranks hold a copy.
module turbogap_structure

   use kinds, only: dp
   use turbogap_comm, only: comm_t, comm_bcast
   use types, only: input_parameters
   use timing, only: times_t, time_start, time_end
   use turbogap_setup, only: model_t
   use turbogap_loop, only: loop_t
   use read_files, only: read_xyz

   implicit none

   private
   public :: structure_acquire

   type, public :: state_t
      integer :: n_sites
      real(dp), allocatable :: positions(:, :)
      real(dp), allocatable :: positions_prev(:, :)
      real(dp), allocatable :: positions_diff(:, :)
      real(dp), allocatable :: forces_prev(:, :)
      real(dp), allocatable :: velocities(:, :)
      real(dp), allocatable :: masses(:)
      integer, allocatable :: species(:)
      integer, allocatable :: species_supercell(:)
      character*8, allocatable :: xyz_species(:)
      character*8, allocatable :: xyz_species_supercell(:)
      logical, allocatable :: fix_atom(:, :)
      real(dp) :: a_box(1:3)
      real(dp) :: b_box(1:3)
      real(dp) :: c_box(1:3)
      integer :: indices(1:3)
!     Volume of the primitive cell, a_box/indices(1) etc.
      real(dp) :: v_uc
!     The time= tag of the frame just read, and whether it was there at all.
      real(dp) :: frame_time = 0.d0
      logical :: has_frame_time = .false.
   end type state_t

contains

!  Rank 0 reads the structure this step works on: the input file on the first
!  MD step or for each frame outside MD, the MC trial file after a move. The
!  other ranks learn whether there is another frame to come; the structure
!  itself reaches them through domain_sync_state.
   subroutine structure_acquire(state, rebuild, loop, params, model, comm, mc_file, time)
      type(state_t), intent(inout) :: state
      logical, intent(inout) :: rebuild
      type(loop_t), intent(inout) :: loop
      type(input_parameters), intent(inout) :: params
      type(model_t), intent(in) :: model
      type(comm_t), intent(in) :: comm
      character(len=*), intent(in) :: mc_file
      type(times_t), intent(inout) :: time

      !   This chunk of code does all the reading/neighbor builds etc for each snapshot
      !   or MD step
      !   Read in XYZ file and build neighbors lists

      if ((params%do_md .and. loop%md_istep == 0)) then
         call time_start(time%read_xyz)
         if (comm%rank == 0) then
            if (loop%mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites,.not. params%mc_write_xyz, state%fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz, params%randomize_velocities)
               rebuild = .true.

            else if (.not. params%do_nested_sampling .or. loop%mc_istep == 0) then
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .false., state%fix_atom, params%t_beg, &
                             params%write_array_property(6), .false., params%randomize_velocities)

            end if

            ! call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
            !               n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
            !               positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
            !               xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
            !               n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), .true. )
            !     Only rank 0 handles these variables
            !      allocate( positions_prev(1:3, 1:size(positions,2)) )
            !      allocate( positions_diff(1:3, 1:size(positions,2)) )
            if (.not. allocated(state%forces_prev)) allocate (state%forces_prev(1:3, 1:state%n_sites))
            if (.not. allocated(state%positions_prev)) allocate (state%positions_prev(1:3, 1:state%n_sites))
            if (.not. allocated(state%positions_diff)) allocate (state%positions_diff(1:3, 1:state%n_sites))
            state%positions_diff = 0.d0
            rebuild = .true.
         end if
         call time_end(time%read_xyz)
         !     If we're doing MD, we don't read beyond the first snapshot in the XYZ file
         loop%repeat_xyz = .false.
         !     At the moment, we can't do prediction if the unit cell doesn't fit a whole cutoff sphere
      else if (.not. params%do_md) then
         call time_start(time%read_xyz)
         if (comm%rank == 0) then
            if (loop%mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites,.not. params%mc_write_xyz, state%fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz, params%randomize_velocities)
               rebuild = .true.
            else
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .false., state%fix_atom, params%t_beg, params%write_array_property(6), &
                             .false., params%randomize_velocities, state%frame_time, state%has_frame_time)
            end if
         end if
         call time_end(time%read_xyz)
         call time_start(time%mpi)
         call comm_bcast(comm, loop%repeat_xyz)
!        The frame's time label. Every rank pushes the same dipole into the
!        same buffer -- the trajectory is the ensemble and it is replicated,
!        not distributed -- so every rank needs the same time with it.
         call comm_bcast(comm, state%frame_time)
         call comm_bcast(comm, state%has_frame_time)
         call time_end(time%mpi)
         rebuild = .true.
      end if
   end subroutine structure_acquire

end module turbogap_structure
