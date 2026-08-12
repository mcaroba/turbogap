! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_md.f90, is copyright (c) 2019-2026, Miguel A. Caro and
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

!  One MD step: the thermostat and barostat, the velocity-Verlet update, the
!  box scaling, the electronic-stopping and electron-phonon coupling, the
!  trajectory and thermo output, and the broadcast of the new positions.
!
!  Lifted verbatim from turbogap.f90. The block carried its own
!  `#ifdef _MPIF90 / IF (rank == 0)` guard and the position broadcast that
!  follows it, so both moved too and the driver is left with a single call.
module turbogap_md

   use kinds
   use types
   use md
   use bussi
   use xyz_module
   use timing
   use mpi_helper
   use adaptive_time
   use electronic_stopping
   use eph_beta
   use eph_fdm
   use eph_electronic_stopping
#ifdef _MPIF90
   use mpi
#endif

   implicit none

   private
   public :: compute_md

contains

!**************************************************************************
   subroutine compute_md(params, rank, ierr, n_sites, n_species, md_istep, md_time, &
                         time_step, positions, positions_prev, positions_diff, velocities, forces, &
                         forces_prev, masses, masses_types, xyz, xyz_species, a_box, b_box, c_box, indices, &
                         v_uc, virial, energy, energy_prev, energies, energies_soap, energies_2b, &
                         energies_3b, energies_core_pot, energies_vdw, energies_lp, energies_exp, &
                         energies_pdf, energies_sf, energies_xrd, energies_nd, local_properties, &
                         local_property_labels, instant_temp, instant_pressure, instant_pressure_prev, &
                         e_kin, e_kinetic, kb, evpera3tobar, fix_atom, exit_loop, rebuild_neighbors_list, &
                         i_image, i_nested, n_pos, nrows, filename, string, allelstopdata, ephbeta, ephfdm, &
                         ephlsc, time, cum_eel, gd_istep, target_temp, time_step_prev, &
                         dipole, local_dipoles, energies_dipole)

      implicit none

      type(input_parameters), intent(inout) :: params
      integer, intent(in) :: rank
      integer, intent(inout) :: ierr
      integer, intent(inout) :: n_sites
      integer, intent(inout) :: n_species
      integer, intent(inout) :: md_istep
      real(dp), intent(inout) :: md_time
      real(dp), intent(inout) :: time_step
      real(dp), allocatable, intent(inout) :: positions(:, :)
      real(dp), allocatable, intent(inout) :: positions_prev(:, :)
      real(dp), allocatable, intent(inout) :: positions_diff(:, :)
      real(dp), allocatable, intent(inout) :: velocities(:, :)
      real(dp), allocatable, intent(inout) :: forces(:, :)
      real(dp), allocatable, intent(inout) :: forces_prev(:, :)
      real(dp), allocatable, intent(in) :: masses(:)
      real(dp), allocatable, intent(in) :: masses_types(:)
      real(dp), allocatable, intent(in) :: xyz(:, :)
      character*8, allocatable, intent(in) :: xyz_species(:)
      real(dp), intent(inout) :: a_box(1:3)
      real(dp), intent(inout) :: b_box(1:3)
      real(dp), intent(inout) :: c_box(1:3)
      integer, intent(inout) :: indices(1:3)
      real(dp), intent(in) :: v_uc
      real(dp), intent(in) :: virial(1:3, 1:3)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: energy_prev
      real(dp), allocatable, intent(in) :: energies(:)
      real(dp), allocatable, intent(inout) :: energies_soap(:)
      real(dp), allocatable, intent(inout) :: energies_2b(:)
      real(dp), allocatable, intent(inout) :: energies_3b(:)
      real(dp), allocatable, intent(inout) :: energies_core_pot(:)
      real(dp), allocatable, intent(inout) :: energies_vdw(:)
      real(dp), allocatable, intent(inout) :: energies_lp(:)
      real(dp), allocatable, intent(inout) :: energies_exp(:)
      real(dp), allocatable, intent(inout) :: energies_pdf(:)
      real(dp), allocatable, intent(inout) :: energies_sf(:)
      real(dp), allocatable, intent(inout) :: energies_xrd(:)
      real(dp), allocatable, intent(inout) :: energies_nd(:)
      real(dp), allocatable, intent(in) :: local_properties(:, :)
!     Dipole model output, written alongside the trajectory. energies_dipole is
!     the model's fictitious scalar and is reported separately from energy=.
      real(dp), intent(in) :: dipole(1:3)
      real(dp), allocatable, intent(in) :: local_dipoles(:, :)
      real(dp), allocatable, intent(in) :: energies_dipole(:)
      character*1024, allocatable, intent(in) :: local_property_labels(:)
      real(dp), intent(inout) :: instant_temp
      real(dp), intent(inout) :: instant_pressure
      real(dp), intent(in) :: instant_pressure_prev
      real(dp), intent(in) :: e_kin
      real(dp), intent(inout) :: e_kinetic
      real(dp), intent(in) :: kb
      real(dp), intent(in) :: evpera3tobar
      logical, allocatable, intent(in) :: fix_atom(:, :)
      logical, intent(inout) :: exit_loop
      logical, intent(inout) :: rebuild_neighbors_list
      integer, intent(in) :: i_image
      integer, intent(in) :: i_nested
      integer, intent(inout) :: n_pos
      integer, intent(in) :: nrows
      character*1024, intent(inout) :: filename
      character*1024, intent(inout) :: string
      real(dp), allocatable, intent(in) :: allelstopdata(:)
      type(EPH_Beta_class), intent(inout) :: ephbeta
      type(EPH_FDM_class), intent(inout) :: ephfdm
      type(EPH_LangevinSpatialCorrelation_class), intent(inout) :: ephlsc
      type(times_t), intent(inout) :: time
      real(dp), intent(inout) :: cum_eel
!     Counts steps taken by the gd-box relaxation. Only the convergence test
!     below reads it; the optimizer keeps its own state.
      integer, intent(inout) :: gd_istep
      real(dp), intent(inout) :: target_temp
      real(dp), intent(inout) :: time_step_prev

      character*64 :: cjunk
      real(dp) :: instant_pressure_tensor(1:3, 1:3)
      real(dp) :: lv(1:3, 1:3)
!     Loop scratch. Written before read here; the driver rewrites i and j in
!     do-loops after the call and never reads i2, j2 or k2 again.
      integer :: i, i2, j, j2, k2

      !**************************************************************************
      !   Do MD stuff here
#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         if (params%do_md .and. md_istep > -1) then
            call time_start(time%md)
            !     Define the time_step and md_time prior to possible scaling (see variable_time_step below)
            if (md_istep > 0) then
               md_time = md_time + time_step
            else
               md_time = 0.d0
               time_step = params%md_step
            end if
            !     We wrap the positions and remoce CM velocity
            call wrap_pbc(positions(1:3, 1:n_sites), a_box&
                 &/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box&
                 &/dfloat(indices(3)))
            call remove_cm_vel(velocities(1:3, 1:n_sites),&
                 & masses(1:n_sites))

            !     First we check if this is a variable time step simulation
            if (params%variable_time_step) then
               call variable_time_step(md_istep == 0, velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
                                       params%target_pos_step, params%tau_dt, params%md_step, time_step)
            end if

           !! ------- option for radiation cascade simulation with electronic stopping

            if (params%electronic_stopping) then
               call electron_stopping_velocity_dependent(md_istep, n_species, params%eel_cut, params%eel_freq_out, &
                                                         velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
                                                   params%masses_types, time_step, md_time, nrows, allelstopdata, cum_EEL, 'forces')
            end if

           !! -----------------------------------        ******** until here for electronic stopping

           !! ------- option for electronic stopping based on eph model

            if (params%nonadiabatic_processes) then
               call ephlsc%eph_LangevinForces(velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), &
                                              masses(1:n_sites), params%masses_types, md_istep, time_step, md_time, &
                                              positions(1:3, 1:n_sites), n_species, ephbeta, ephfdm)
            end if

          !! -----------------------------------        ******** until here for electronic stopping basd on eph model

          !! ------- option for doing simulation with adaptive time step

            if (params%adaptive_time) then
               if (MOD(md_istep, params%adapt_tstep_interval) == 0) then
                  call variable_time_step_adaptive(md_istep == 0, velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), &
                                                   masses(1:n_sites), params%adapt_tmin, params%adapt_tmax, params%adapt_xmax, &
                                                   params%adapt_emax, params%md_step, time_step)
               end if
            end if

          !! ----------------------------------        ******** until here for adaptive time

            !     This takes care of NVE
            !     Velocity Verlet takes positions for t, positions_prev for t-dt, and velocities for t-dt and returns everything
            !     dt later. forces are taken at t, and forces_prev at t-dt. forces is left unchanged by the routine, and
            !     forces_prev is returned as equal to forces (both arrays contain the same information on return)
            if (params%optimize == "vv") then
               call velocity_verlet(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                                forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), time_step, time_step_prev, &
                                    md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                                    fix_atom(1:3, 1:n_sites))
            else if (params%optimize == "gd") then
               call gradient_descent(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                                     forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), &
                                     params%max_opt_step, md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), &
                                     c_box/dfloat(indices(3)), fix_atom(1:3, 1:n_sites), energy)
            else if (params%optimize == "gd-box" .or. params%optimize == "gd-box-ortho") then
               !       Positions and lattice descend together, in the
               !       preconditioned variables of Gubler et al. (2023). The
               !       old scheme alternated -- relax positions to convergence,
               !       then the box to convergence, then back -- and each half
               !       undid part of the other's work. The lattice update that
               !       used to live in the box-rescaling section below is part
               !       of this one call now.
               positions_prev(1:3, 1:n_sites) = positions(1:3, 1:n_sites)
               forces_prev(1:3, 1:n_sites) = forces(1:3, 1:n_sites)
               call gradient_descent_positions_and_lattice(positions(1:3, 1:n_sites), &
                                                           velocities(1:3, 1:n_sites), &
                                                           forces(1:3, 1:n_sites), virial(1:3, 1:3), &
                                                           energy, a_box, b_box, c_box, indices, &
                                                           fix_atom(1:3, 1:n_sites), &
                                                           params%gd_box_weight, params%max_opt_step, &
                                                           params%optimize == "gd-box-ortho", &
                                                           md_istep == 0)
               !       Restarts at zero for each relaxation, so that the
               !       gd_istep > 1 guard on the convergence test below means
               !       "this relaxation has taken a step or two" even under mc,
               !       where compute_md is re-entered once per accepted move.
               if (md_istep == 0) gd_istep = 0
               gd_istep = gd_istep + 1
            else
               !       If nothing happens we still update these variables
               positions_prev(1:3, 1:n_sites) = positions(1:3, 1:n_sites)
               forces_prev(1:3, 1:n_sites) = forces(1:3, 1:n_sites)
            end if

        !! ------- option for radiation cascade simulation with electronic stopping

            if (params%electronic_stopping) then
               call electron_stopping_velocity_dependent(md_istep, n_species, params%eel_cut, params%eel_freq_out, &
                                                         velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
                                                   params%masses_types, time_step, md_time, nrows, allelstopdata, cum_EEL, 'energy')
            end if

        !! -----------------------------------                ******** until here for electronic stopping

        !! ------- option for electronic stopping based on eph model

            if (params%nonadiabatic_processes) then
               call ephlsc%eph_LangevinEnergyDissipation(md_istep, md_time, velocities(1:3, 1:n_sites), &
                                                         positions(1:3, 1:n_sites), time_step, ephfdm)
            end if

        !! -----------------------------------                ******** until here for electronic stopping basd on eph model

            !     Compute kinetic energy from current velocities. Because Velocity Verlet
            !     works with the velocities at t-dt (except for the first time step) we
            !     have to compute the velocities after call Verlet
            E_kinetic = 0.d0
            do i = 1, n_sites
               E_kinetic = E_kinetic + 0.5d0*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
            end do
            instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic

            !     Instant pressure in bar
            instant_pressure = (kB*dfloat(n_sites - 1)*instant_temp&
                 &+ (virial(1, 1) + virial(2, 2) + virial(3, 3))/3.d0)&
                 &/v_uc*eVperA3tobar
            instant_pressure_tensor(1:3, 1:3) = virial(1:3, 1:3)/v_uc&
                                               &*eVperA3tobar
            do i = 1, 3
               instant_pressure_tensor(i, i) =&
                    & instant_pressure_tensor(i, i) + (kB&
                    &*dfloat(n_sites - 1)*instant_temp)/v_uc*eVperA3tobar
            end do

            !     Here we write thermodynamic information -> THIS NEEDS CLEAN UP AND IMPROVEMENT
            if (md_istep == 0 .and. .not. params%do_nested_sampling) then
               open (unit=10, file="thermo.log", status="unknown")
               write (10, "(A,A)") "#     Step             Time      Temperature                E_kin                     E_pot", &
                  "             Pressure"
            else if (md_istep == 0 .and. i_nested == 1) then
               open (unit=10, file="thermo.log", status="unknown")
               write (10, "(A,A)") "#     Step             Time      Temperature                E_kin                     E_pot", &
                  "             Pressure"
            else
               open (unit=10, file="thermo.log", status="old", position="append")
            end if
            if (.not. params%do_mc .and. (md_istep == 0 .or. md_istep == params%md_nsteps &
                                          .or. modulo(md_istep, params%write_thermo) == 0)) then
               !       Organize this better so that the user can have more freedom about what gets printed to thermo.log
               !       There should also be a header preceded by # specifying what gets printed
               if (params%do_exp) then
                  write (10, "(I10, 1X, F16.4, 1X, F16.4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8)", advance="no") &
                     md_istep, md_time, instant_temp, E_kinetic, sum(energies), sum(energies_exp), instant_pressure
               else
                  write (10, "(I10, 1X, F16.4, 1X, F16.4, 1X, F20.8, 1X, F20.8, 1X, F20.8)", advance="no") &
                     md_istep, md_time, instant_temp, E_kinetic, sum(energies), instant_pressure
               end if

               if (params%write_lv) then
                  write (10, "(1X, 9F20.8)", advance="no") a_box(1:3)/dfloat(indices(1)), &
                     b_box(1:3)/dfloat(indices(2)), &
                     c_box(1:3)/dfloat(indices(3))
               end if
               !       Further printouts should go here
               !       <<HERE>>
               !
               !       This is to make the pointer advance
               write (10, *)
            end if
            close (10)
            !
            !     Check if we have converged a relaxation calculation
            !     Check if we have converged a relaxation calculation
            if (params%do_md .and. params%optimize == "gd" .and. md_istep > 0 .and. &
                abs(energy - energy_prev) < params%e_tol*dfloat(n_sites) .and. &
                maxval(abs(forces)) < params%f_tol .and. rank == 0) then
               exit_loop = .true.
               if (params%do_mc) exit_loop = .false.
               !     THIS CONDITION ON INSTANT PRESSURE WILL NEED TO BE FINE TUNED, TO ACCOUNT FOR ARBITRARY TARGET PRESSURES
               !     BUT ALSO TO ACCOMMODATE NON-TRICLINIC TARGET BOX SHAPES, WHERE IT MIGHT NOT BE POSSIBLE TO CONVERGE THE
               !     TOTAL PRESSURE BELOW A CERTAIN MINIMUM (DUE TO THE BOX SHAPE CONSTRAINTS)
            else if (params%do_md .and. (params%optimize == "gd-box"&
                 & .or. params%optimize == "gd-box-ortho") .and.&
                 & gd_istep > 1 .and. abs(energy - energy_prev) < params&
                 &%e_tol*dfloat(n_sites) .and. abs(instant_pressure -&
                 & instant_pressure_prev) < params%p_tol .and.&
                 & maxval(abs(forces)) < params%f_tol .and. rank == 0&
                 & ) then
               exit_loop = .true.
               if (params%do_mc) exit_loop = .false.
            end if

            !     We write out the trajectory file. We write positions_prev which is the one for which we have computed
            !     the properties. positions_prev and velocities are synchronous
            if ((md_istep == 0 .and. .not. params%do_nested_sampling) .or. &
                (md_istep == params%md_nsteps .and. .not. params%do_nested_sampling) &
                .or. (modulo(md_istep, params%write_xyz) == 0 .and. .not. params%do_nested_sampling) .or. &
                exit_loop) then
               call wrap_pbc(positions_prev(1:3, 1:n_sites), a_box&
                    &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                    & c_box/dfloat(indices(3)))
               call get_xyz_energy_string(energies_soap, energies_2b,&
                    & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                    &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                    & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                    & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                    & params%do_dipole, dipole, energies_dipole)

               call write_extxyz(n_sites, md_istep, md_time, time_step,&
                    & instant_temp, instant_pressure, a_box&
                    &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                    & c_box/dfloat(indices(3)), virial, xyz_species,&
                    & positions_prev(1:3, 1:n_sites), velocities,&
                    & forces, energies(1:n_sites), masses, params&
                    &%write_property, params%write_array_property,&
                    & params%write_local_properties,&
                    & local_property_labels, local_properties,&
                    & fix_atom(1:3, 1:n_sites), "trajectory_out.xyz", string, &
                    & md_istep == 0, params%do_dipole, local_dipoles(1:3, 1:n_sites))
            else if (md_istep == params%md_nsteps .and. params%do_nested_sampling) then
               write (cjunk, '(I8)') i_image
               write (filename, '(A,A,A)') "walkers/", trim(adjustl(cjunk)), ".xyz"
               call wrap_pbc(positions_prev(1:3, 1:n_sites), &
                             a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
               call get_xyz_energy_string(energies_soap, energies_2b,&
                    & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                    &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                    & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                    & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                    & params%do_dipole, dipole, energies_dipole)

               call write_extxyz(n_sites, md_istep, md_time, time_step, instant_temp, instant_pressure, &
                    a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                    virial, xyz_species, &
                    positions_prev(1:3, 1:n_sites), velocities, &
                    forces, energies(1:n_sites), masses, &
                    params%write_property, params%write_array_property&
                    &, params%write_local_properties,&
                    & local_property_labels, local_properties,&
                    & fix_atom(1:3, 1:n_sites), filename, string, .true.,&
                    & params%do_dipole, local_dipoles(1:3, 1:n_sites))

            end if
            !
            !     If there are pressure/box rescaling operations they happen here
            if (params%scale_box) then
               call box_scaling(positions(1:3, 1:n_sites), a_box(1:3), b_box(1:3), c_box(1:3), &
                                indices, md_istep, params%md_nsteps, params%box_scaling_factor)
            else if (params%barostat == "berendsen") then
               lv(1:3, 1) = a_box(1:3)
               lv(1:3, 2) = b_box(1:3)
               lv(1:3, 3) = c_box(1:3)
               call berendsen_barostat(lv(1:3, 1:3), &
                                       params%p_beg + (params%p_end - params%p_beg)*dfloat(md_istep + 1)/float(params%md_nsteps), &
                                       instant_pressure_tensor, params%barostat_sym, params%tau_p, params%gamma_p, time_step)
               a_box(1:3) = lv(1:3, 1)
               b_box(1:3) = lv(1:3, 2)
               c_box(1:3) = lv(1:3, 3)
               call berendsen_barostat(positions(1:3, 1:n_sites), &
                                       params%p_beg + (params%p_end - params%p_beg)*dfloat(md_istep + 1)/float(params%md_nsteps), &
                                       instant_pressure_tensor, params%barostat_sym, params%tau_p, params%gamma_p, time_step)
            end if
            !     gd-box has no separate box-rescaling stage any more: the
            !     lattice moved with the positions in the single
            !     gradient_descent_positions_and_lattice call above.
            !     If there are thermostating operations they happen here
            if (params%thermostat == "berendsen") then
               call get_target_temp(params%t_beg, params%t_end,&
                    & md_istep, params%md_nsteps, params%n_t_hold, &
                    & params%t_hold, target_temp)
               call berendsen_thermostat(velocities(1:3, 1:n_sites), &
                                         target_temp, &
                                         instant_temp, params%tau_t, time_step)
            else if (params%thermostat == "bussi") then
               call get_target_temp(params%t_beg, params%t_end,&
                    & md_istep, params%md_nsteps, params%n_t_hold, &
                    & params%t_hold, target_temp)
               if (E_kinetic > 0.0d0) then
                  velocities(1:3, 1:n_sites) = velocities(1:3, 1:n_sites) &
                                               *dsqrt(resamplekin(E_kinetic, target_temp, &
                                                                  3*n_sites - 3, params%tau_t, time_step)/E_kinetic)
               else
                  velocities(1:3, 1:n_sites) = 0.0d0
               end if
            end if
            !     Check what's the maximum atomic displacement since last neighbors build
            positions_diff = positions_diff + positions(1:3, 1:n_sites) - positions_prev(1:3, 1:n_sites)
            rebuild_neighbors_list = .false.
            !--------
            ! CHECK THIS OUT and fix it at some point
            ! Here we set the neighbors list rebuild to always true if the supercell and the primitive unit cell are not
            ! the same. This is because of how atoms get wrapped around the PBC during MD (they get wrapped around the
            ! primitive unit cell) making the neighbors lists obsolete. This is an issue with wrapping, not with the neighbors
            ! lists. A possible solution would be to wrap around the supercell, instead of the unit cell, to maintain the
            ! internal consistency of the positions(:,:) array, and then do a wrapping around the primitive unit cell for
            ! printing the XYZ coordinates only (i.e., keeping the other wrapping convention internally for positions(:,:))
            if (any(indices > 1)) then
               rebuild_neighbors_list = .true.
            end if
            !--------
            do i = 1, n_sites
               if (positions_diff(1, i)**2 + positions_diff(2, i)**2 + positions_diff(3, i)**2 > params%neighbors_buffer/2.d0) then
                  rebuild_neighbors_list = .true.
                  positions_diff = 0.d0
                  exit
               end if
            end do
            !       We make sure the atoms in the supercell have the same positions and velocities as in the unit cell
            j = 0
            do i2 = 1, indices(1)
               do j2 = 1, indices(2)
                  do k2 = 1, indices(3)
                     do i = 1, n_sites
                        j = j + 1
                        if (j > n_sites) then
                           positions(1:3, j) = positions(1:3, i) + dfloat(i2 - 1)/dfloat(indices(1))*a_box &
                                               + dfloat(j2 - 1)/dfloat(indices(2))*b_box &
                                               + dfloat(k2 - 1)/dfloat(indices(3))*c_box
                           velocities(1:3, j) = velocities(1:3, i)
                        end if
                     end do
                  end do
               end do
            end do
            call time_end(time%md)
         end if
#ifdef _MPIF90
      END IF
      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
      !   Make sure all ranks have correct positions and velocities
#ifdef _MPIF90
      if (params%do_md) then
         call time_start(time%mpi_positions)
         n_pos = size(positions, 2)
         call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(velocities, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call time_end(time%mpi_positions)
      end if
#endif

   end subroutine compute_md
!**************************************************************************

end module turbogap_md
