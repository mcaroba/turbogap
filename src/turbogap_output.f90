! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_output.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  What the driver itself prints: the banner and the options summary.
module turbogap_output

   use kinds, only: dp
   use turbogap_comm, only: comm_t, comm_with_mpi
   use types, only: input_parameters, perform_t
   use timing, only: times_t, get_time, sum_times
   use md, only: wrap_pbc
   use xyz_module, only: get_xyz_energy_string, write_extxyz
   use turbogap_setup, only: model_t
   use turbogap_structure, only: state_t
   use turbogap_domain, only: domain_t
   use turbogap_results, only: results_t
   use turbogap_loop, only: loop_t
   use turbogap_md, only: dynamics_t
   use turbogap_ir, only: ir_run_t, ir_report

   implicit none

   private
   public :: print_banner
   public :: print_options
   public :: print_single_point_energies
   public :: write_debug_forces
   public :: write_single_point
   public :: print_nothing_to_do
   public :: print_timing_report

contains

   subroutine print_banner(comm)
      type(comm_t), intent(in) :: comm

      if (comm%rank == 0) then
         write (*, *) '_________________________________________________________________ '
         write (*, *) '                             _                                   \'
         write (*, *) ' ___________            __   \\ /\        _____     ___   _____  |'
         write (*, *) '/____  ____/           / / /\|*\|*\/\    / ___ \   /   | |  _  \ |'
         write (*, *) '    / / __  __  __    / /  \********/   / /  /_/  / /| | | / | | |'
         write (*, *) '   / / / / / / / /_  / /__  \**__**/   / / ____  / / | | | |_/ / |'
         write (*, *) '  / / / / / / / __/ / ___ \ /*/  \*\  / / /_  / / /__| | |  __/  |'
         write (*, *) ' / / / /_/ / / /   / /__/ / \ \__/ / / /___/ / / ____  | | |     |'
         write (*, *) '/_/_/_____/_/_/___/______/___\____/__\______/_/_/____|_|_|_|____ |'
         write (*, *) '_____________________________________________________________  / |'
         write (*, *) '*************************************************************|/  |'
         write (*, *) '                  Welcome to the TurboGAP code                   |'
         write (*, *) '                         Maintained by                           |'
         write (*, *) '                                                                 |'
         write (*, *) '               Miguel A. Caro and Tigany Zarrouk                 |'
         write (*, *) '                       mcaroba@gmail.com                         |'
         write (*, *) '                      miguel.caro@aalto.fi                       |'
         write (*, *) '                                                                 |'
         write (*, *) '          Department of Chemistry and Materials Science          |'
         write (*, *) '                     Aalto University, Finland                   |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '====================>>>>>  turbogap.fi  <<<<<====================|'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) 'Contributors (code and methodology) in chronological order:      |'
         write (*, *) '                                                                 |'
         write (*, *) 'Miguel A. Caro, Patricia Hernández-León, Suresh Kondati          |'
         write (*, *) 'Natarajan, Albert P. Bartók, Eelis V. Mielonen, Heikki Muhli,    |'
         write (*, *) 'Mikhail Kuklin, Gábor Csányi, Jan Kloppenburg, Richard Jana,     |'
         write (*, *) 'Tigany Zarrouk, Cristian V. Achim                                |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '                     Last updated: Jun. 2026                     |'
         write (*, *) '                                        _________________________/'
         write (*, *) '.......................................|'
         if (comm_with_mpi) then
            write (*, *) '                                       |'
            write (*, *) 'Running TurboGAP with MPI support:     |'
            write (*, *) '                                       |'
            write (*, '(A,I6,A)') ' Running TurboGAP on ', comm%size, ' MPI tasks   |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         else
            write (*, *) '                                       |'
            write (*, *) 'Running the serial version of TurboGAP |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if
   end subroutine print_banner

   subroutine print_options(comm, params, model)
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      character*8 :: i_char
      integer :: i

      if (comm%rank == 0) then
         ! Print out chosen options:
         write (*, *) '                                       |'
         write (*, '(1X,A)') 'You specified the following options:   |'
         write (*, *) '                                       |'
         write (*, *) '---------------------------------      |'
         if (len(trim(params%atoms_file)) > 20) then
            write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '...   |'
         else
            write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '      |'
         end if
         write (*, *) '---------------------------------      |'
         write (i_char, '(I8)') model%n_species
         write (*, '(1X,A,A8,A)') 'No. of species   = ', adjustl(i_char), '            |'
         do i = 1, model%n_species
            write (i_char, '(I8)') i
           write (*, '(1X,A,A2,A,A8,A)') '  *) Species #', adjustl(i_char), ' =       ', adjustr(params%species_types(i)), '      |'
         end do
         write (*, *) '---------------------------------      |'
         write (*, '(1X,A,F15.4,A)') 'rcut_max = ', model%rcut_max, ' Angst.      |'
         write (*, *) '---------------------------------      |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end subroutine print_options

!  The per-term energies of a single-point prediction.
   subroutine print_single_point_energies(comm, res, params, model, perform)
      type(comm_t), intent(in) :: comm
      type(results_t), intent(in) :: res
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(perform_t), intent(in) :: perform

      if (comm%rank == 0) then
         write (*, *) '                                       |'
         write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
         write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
         write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
         write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
         write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
         write (*, '(A,1X,F21.8,1X,A)') ' estat energy:', sum(res%energies_estat), 'eV |'
         write (*, '(A,1X,F22.8,1X,A)') ' Exp. energy:', sum(res%energies_exp), 'eV |'
         if (model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(res%energies_lp), 'eV |'
         if (perform%pdf)&
              & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
              & sum(res%energies_pdf), 'eV |'
         if (perform%sf)&
              & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
              & sum(res%energies_sf), 'eV |'
         if (perform%xrd)&
              & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
              & sum(res%energies_xrd), 'eV |'
         if (perform%nd)&
              & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
              & sum(res%energies_nd), 'eV |'

         write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(res%energies), 'eV |'

         write (*, *) '                                       |'
         write (*, *) 'Energy & forces in "trajectory_out.xyz"|'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end subroutine print_single_point_energies

!  The force and virial dumps the print_*_forces keywords ask for.
   subroutine write_debug_forces(comm, res, state, dom, params, model, perform, loop)
      type(comm_t), intent(in) :: comm
      type(results_t), intent(in) :: res
      type(state_t), intent(in) :: state
      type(domain_t), intent(in) :: dom
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(perform_t), intent(in) :: perform
      type(loop_t), intent(in) :: loop
      character*1024 :: temp_string
      character*1024 :: temp_string2
      integer :: i
      integer :: j

      if (comm%rank == 0 .and. params%print_vdw_forces) then
         print *, "> Virial ESTAT "
         do i = 1, 3
            do j = 1, 3
               print *, " i, ", i, " j ", j, " ", res%virial_estat(i, j)
            end do
         end do

         print *, "> Virial soap "
         do i = 1, 3
            do j = 1, 3
               print *, " i, ", i, " j ", j, " ", res%virial_soap(i, j)
            end do
         end do

         print *, "> Virial 2b "
         do i = 1, 3
            do j = 1, 3
               print *, " i, ", i, " j ", j, " ", res%virial_2b(i, j)
            end do
         end do

         print *, "> Virial 3b "
         do i = 1, 3
            do j = 1, 3
               print *, " i, ", i, " j ", j, " ", res%virial_3b(i, j)
            end do
         end do

         print *, "> Virial core_pot "
         do i = 1, 3
            do j = 1, 3
               print *, " i, ", i, " j ", j, " ", res%virial_core_pot(i, j)
            end do
         end do

         if (perform%xrd_forces) then
            print *, "> Virial xrd "
            do i = 1, 3
               do j = 1, 3
                  print *, " i, ", i, " j ", j, " ", res%virial_xrd(i, j)
               end do
            end do
            temp_string = ""
            temp_string2 = ""
            write (temp_string, "(I8)") loop%md_istep
            write (temp_string2, "(A)") "forces_xrd_"//trim(adjustl(temp_string))
            open (unit=90, file=temp_string2, status="unknown")
            do i = 1, state%n_sites
               write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                  res%forces_xrd(1, i), res%forces_xrd(2, i), res%forces_xrd(3, i)
            end do
            close (90)

         end if

      end if

      if (params%print_vdw_forces) then
         open (unit=90, file="forces_vdw", status="unknown")
         do i = 1, state%n_sites
            write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
               res%forces_vdw(1, i), res%forces_vdw(2, i), res%forces_vdw(3, i)
         end do
         close (90)

      end if

      if (comm%rank == 0 .and. params%print_estat_forces) then
         open (unit=90, file="forces_estat", status="unknown")
         do i = 1, state%n_sites
            write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
               res%forces_estat(1, i), res%forces_estat(2, i), res%forces_estat(3, i)
         end do
         close (90)

         open (unit=90, file="charge_gradients_estat", status="unknown")
         do i = 1, dom%n_atom_pairs_by_rank(comm%rank + 1)
            write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
               res%local_properties_cart_der(1, i, model%charge_lp_index), &
               res%local_properties_cart_der(2, i, model%charge_lp_index), &
               res%local_properties_cart_der(3, i, model%charge_lp_index)
         end do
         close (90)

      end if
   end subroutine write_debug_forces

!  A single-point prediction's frame of trajectory_out.xyz.
   subroutine write_single_point(comm, res, state, dyn, params, model, loop)
      type(comm_t), intent(in) :: comm
      type(results_t), intent(in) :: res
      type(state_t), intent(inout) :: state
      type(dynamics_t), intent(in) :: dyn
      type(input_parameters), intent(inout) :: params
      type(model_t), intent(in) :: model
      type(loop_t), intent(in) :: loop
      character*1024 :: string

      if (comm%rank == 0) then
         !       Write energy and forces if we're just doing static predictions
         !       The masses should be divided by 103.6426965268d0 to have amu units, but
         !       since masses is not allocated for single point calculations, it would
         !       likely lead to a segfault
         call wrap_pbc(state%positions(1:3, 1:state%n_sites), state%a_box&
              &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
              & state%c_box/dfloat(state%indices(3)))
         call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
              & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
              &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
              & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
              & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
              & params%do_dipole, res%dipole, res%energies_dipole)

         call write_extxyz(state%n_sites, -loop%n_xyz, dyn%md_time, dyn%time_step,&
              & dyn%instant_temp, dyn%instant_pressure, state%a_box&
              &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
              & state%c_box/dfloat(state%indices(3)), res%virial, state%xyz_species,&
              & state%positions(1:3, 1:state%n_sites), state%velocities, res%forces,&
              & res%energies(1:state%n_sites), state%masses, params&
              &%write_property, params%write_array_property,&
              & params%write_local_properties, model%local_property_labels, res%local_properties, &
              & state%fix_atom, "trajectory_out.xyz", string, .false.,&
              & params%do_dipole, res%local_dipoles(1:3, 1:state%n_sites))

      end if
   end subroutine write_single_point

   subroutine print_nothing_to_do(comm)
      type(comm_t), intent(in) :: comm

      if (comm%rank == 0) then
         !     Do nothing
         write (*, *) '                                       |'
         write (*, *) 'You didn''t ask me to do anything!      |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end subroutine print_nothing_to_do

!  The timing report at the end of the run, on rank 0.
   subroutine print_timing_report(comm, params, model, loop, ir, do_electrostatics, time3, time)
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(loop_t), intent(in) :: loop
      type(ir_run_t), intent(inout) :: ir
      logical, intent(in) :: do_electrostatics
      real(dp), intent(in) :: time3
      type(times_t), intent(inout) :: time
      real(dp) :: time2

      if (params%do_md .or. params%do_prediction .or. params%do_mc) then
         call get_time(time2)
         if (comm%rank == 0) then
            if (params%do_md .and. .not. params%do_nested_sampling) then
               write (*, *) '                                       |'
               write (*, '(I8,A,F13.3,A)') loop%md_istep, ' MD steps:', time2 - time3, ' seconds |'
            end if
            if (params%do_mc) then
               write (*, *)
               write (*, *) '                                       |'
               write (*, '(I8,A,F13.3,A)') loop%mc_istep, ' MC steps:', time2 - time3, ' seconds |'
            end if

            write (*, *) '                                       |'
            write (*, '(A,F13.3,A)') ' *          Setup:', time%setup(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     - input+pot.:', time%read_input(3), ' seconds |'
            if (comm_with_mpi) then
               write (*, '(A,F13.3,A)') '     -  MPI setup:', time%mpi_setup(3), ' seconds |'
            end if
            write (*, '(A,F13.3,A)') ' * Read XYZ files:', time%read_xyz(3), ' seconds |'
            write (*, '(A,F13.3,A)') ' * Neighbor lists:', time%neigh(3), ' seconds |'
            write (*, '(A,F13.3,A)') ' *  GAP desc/pred:', time%gap(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     - soap_turbo:', time%soap(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -         2b:', time%gap_2b(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -         3b:', time%gap_3b(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -   core_pot:', time%gap_core_pot(3), ' seconds |'
!       vdw is a parent, not one of the GAP children above: compute_vdw runs
!       outside the time%gap region and sum_times adds it in its own right.
!       Printing it indented under GAP said otherwise.
            if (params%vdw_type /= "none") then
               write (*, '(A,F13.3,A)') ' *            vdw:', time%vdw(3), ' seconds |'
            end if
            if (model%valid_xps .or. params%do_pair_distribution .or. params&
                 &%do_structure_factor .or. params%do_xrd .or. params%do_nd) write (*, '(A&
                 &,F13.3,A)') ' *  Exp. pred.   :', time%pdf(3) + time%sf(3) + time%xrd(3) + time%nd(3), ' seconds&
                 & |'
            if (model%valid_xps) write (*, '(A,F13.3,A)') '     -        xps:',&
                 & time%xps(3), ' seconds |'
            if (params%do_pair_distribution) write (*, '(A,F13.3,A)') '     -        pdf:', time%pdf(3), ' seconds |'
            if (params%do_structure_factor) write (*, '(A,F13.3,A)') '     -         sf:', time%sf(3), ' seconds |'
            if (params%do_xrd) write (*, '(A,F13.3,A)') '     -        xrd:', time%xrd(3), ' seconds |'
            if (params%do_nd) write (*, '(A,F13.3,A)') '     -         nd:', time%nd(3), ' seconds |'

            call ir_report(ir, params, time)

            if (do_electrostatics) then
               write (*, '(A,F13.3,A)') ' * Electrostatics:', time%estat(3), ' seconds |'
            end if
            if (params%do_md) then
               write (*, '(A,F13.3,A)') ' *  MD algorithms:', time%md(3), ' seconds |'
            end if
            if (params%do_mc) then
               write (*, '(A,F13.3,A)') ' *  MC algorithms:', time%mc(3), ' seconds |'
            end if

            if (comm_with_mpi) then
               write (*, '(A,F13.3,A)') ' *  MPI comms.   :', time%mpi(3) + time%mpi_positions(3) + time%mpi_ef(3), ' seconds |'
               write (*, '(A,F13.3,A)') '     -  pos & vel:', time%mpi_positions(3), ' seconds |'
               write (*, '(A,F13.3,A)') '     - E & F brc.:', time%mpi_ef(3), ' seconds |'
               write (*, '(A,F13.3,A)') '     -  MPI misc.:', time%mpi(3), ' seconds |'
            end if
!       Miscellaneous is what the parent buckets do not account for.  It used
!       to be written out here as one long subtraction, which is how it came to
!       subtract time%gap and the mpi_ef reduce nested inside it and print a
!       negative number.  sum_times owns the list now (src/timing.f90), so the
!       set summed here and the set declared as parents there cannot disagree.
!
!       Accounted-for is printed beside it so the arithmetic is visible: the
!       three numbers below have to add up, and a reader can see at a glance how
!       much of the run the buckets actually name.  A large Miscellaneous is a
!       statement that something real is not being measured -- which is how the
!       setup bucket above came to exist.
            time%total(3) = time2 - time3
            write (*, *) '                                       |'
            write (*, '(A,F13.3,A)') ' *  Accounted for:', sum_times(time), ' seconds |'
            write (*, '(A,F13.3,A)') ' *  Miscellaneous:', time%total(3) - sum_times(time), ' seconds |'
            write (*, '(A,F13.3,A)') ' *     Total time:', time%total(3), ' seconds |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if
   end subroutine print_timing_report

end module turbogap_output
