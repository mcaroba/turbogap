! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_sampling.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  Nested sampling and Monte Carlo: the pool of images they keep, and what each
!  carries from one step to the next.
module turbogap_sampling

   use kinds, only: dp
   use types, only: image, input_parameters, perform_t, from_properties_to_image, from_image_to_properties
   use md, only: remove_cm_vel, get_ns_unbiased_volume_proposal, volume_preserving_strain_transformation, &
                 randomize_velocities, wrap_pbc
   use mc, only: get_accessible_volume, get_mc_acceptance, perform_mc_step
   use xyz_module, only: get_xyz_energy_string, write_extxyz
   use read_files, only: read_xyz
   use timing, only: times_t, time_start, time_end
   use turbogap_setup, only: model_t
   use soap_turbo_functions, only: cross_product
   use turbogap_comm, only: comm_t, comm_bcast
   use turbogap_results, only: results_t
   use turbogap_domain, only: neighbors_t
   use turbogap_structure, only: state_t
   use turbogap_loop, only: loop_t, creturn
   use turbogap_md, only: dynamics_t

   implicit none

   private
   public :: sampling_init
   public :: mc_prepare_step
   public :: nested_step
   public :: mc_step

   type, public :: sampling_t
!     The image pool: the nested-sampling walkers, or MC's current and trial
!     configurations.
      type(image), allocatable :: images(:)
      type(image), allocatable :: images_temp(:)
      integer :: i_image
      integer :: i_current_image = 1
      integer :: i_trial_image = 2
!     Nested sampling
      integer :: i_nested
      integer :: i_max
      real(dp) :: e_max
      real(dp) :: rand
      real(dp) :: rand_scale(1:6)
!     Monte Carlo
      character*32 :: mc_move = "none"
      character*1024 :: mc_file = "mc_trial.xyz"
      integer, allocatable :: mc_id(:)
      integer :: mc_mu_id = 1
      integer, allocatable :: n_mc_species(:)
      integer, allocatable :: n_mc_species_prev(:)
      integer, allocatable :: species_idx(:)
      logical :: do_mc_relax = .false.
!  Whether the trial about to be tested came out of an MD or relaxation burst,
!  captured before params%do_md is cleared just below.
      logical :: trial_came_from_md = .false.
!  Molecule bookkeeping for grand-canonical moves that exchange whole molecules
!  (mc_molecule_files). mc_mol_id tags each atom with the inserted copy it
!  belongs to, mc_mol_mu with the chemical potential that copy came from, and
!  both are zero for a free atom. mc_mol_next hands out the tags. Rank 0 only:
!  the energy evaluation never reads any of it, so none of it is broadcast.
      integer, allocatable :: mc_mol_id(:)
      integer, allocatable :: mc_mol_mu(:)
      integer :: mc_mol_next = 0
      integer :: temp_md_nsteps
      real(dp) :: v_uc_prev
      real(dp) :: v_a_uc
      real(dp) :: v_a_uc_prev
      real(dp) :: ranf
      real(dp) :: disp(1:3)
      real(dp) :: d_disp
      real(dp) :: p_accept
      real(dp) :: virial_prev(1:3, 1:3)
   end type sampling_t

contains

!  Hamiltonian MC: before each move rank 0 redraws the velocities at t_beg,
!  keeping the kinetic energy of the current state after the first move.
   subroutine mc_prepare_step(smp, dyn, state, params, loop, comm)
      type(sampling_t), intent(inout) :: smp
      type(dynamics_t), intent(inout) :: dyn
      type(state_t), intent(inout) :: state
      type(input_parameters), intent(in) :: params
      type(loop_t), intent(in) :: loop
      type(comm_t), intent(in) :: comm
      integer :: i

      if (comm%rank == 0) then
         if (params%do_mc .and. (smp%mc_move /= "md" .or. loop%md_istep == 0) .and. params%mc_hamiltonian) then
            if (loop%mc_istep > 0) dyn%E_kinetic_prev = dyn%E_kinetic
            call random_number(state%velocities)
            call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), state%masses(1:state%n_sites))
            dyn%E_kinetic = 0.d0
            do i = 1, state%n_sites
               dyn%E_kinetic = dyn%E_kinetic + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), state%velocities(1:3, i))
            end do
            dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
            state%velocities = state%velocities*dsqrt(params%t_beg/dyn%instant_temp)
            if (loop%mc_istep > 0) then
               dyn%E_kinetic = dyn%E_kinetic_prev
               dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
               ! Reversing as we want it to be at the instant temp and not at t_beg
               state%velocities = state%velocities*dsqrt(dyn%instant_temp/params%t_beg)

               do i = 1, state%n_sites
                  dyn%E_kinetic = dyn%E_kinetic + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), &
                                                                                    state%velocities(1:3, i))
               end do
            else
               dyn%E_kinetic = dyn%E_kinetic*params%t_beg/dyn%instant_temp
            end if
         end if
      end if
   end subroutine mc_prepare_step

   subroutine sampling_init(smp)
      type(sampling_t), intent(inout) :: smp

      smp%i_nested = 0
      smp%i_image = 0
   end subroutine sampling_init

!  Nested sampling: collect the initial walkers, then replace the highest-
!  enthalpy walker with an MD-decorrelated clone of another one.
   subroutine nested_step(smp, state, res, dyn, nl, params, loop, comm)
      type(sampling_t), intent(inout) :: smp
      type(state_t), intent(inout) :: state
      type(results_t), intent(inout) :: res
      type(dynamics_t), intent(inout) :: dyn
      type(neighbors_t), intent(inout) :: nl
      type(input_parameters), intent(inout) :: params
      type(loop_t), intent(inout) :: loop
      type(comm_t), intent(in) :: comm
      real(dp) :: u
      integer :: i

      !   This runs at the beginning to read in the initial images
      if (params%do_nested_sampling .and. loop%n_xyz > smp%i_image .and. .not. params%do_md) then
         smp%i_image = smp%i_image + 1
         if (.not. allocated(smp%images)) then
            allocate (smp%images(1:smp%i_image))
         else
            allocate (smp%images_temp(1:smp%i_image))
            smp%images_temp(1:smp%i_image - 1) = smp%images(1:smp%i_image - 1)
            deallocate (smp%images)
            allocate (smp%images(1:smp%i_image))
            smp%images = smp%images_temp
            deallocate (smp%images_temp)
         end if
         !     Save initial pool of structures
         state%velocities = 0.d0
         call from_properties_to_image(smp%images(smp%i_image), state%positions, state%velocities, state%masses, &
                                       res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                       res%energy_exp, dyn%E_kinetic, &
                                       state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                       state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                       res%local_dipoles, res%energies_dipole, res%dipole)
      end if

      !   This handles the nested sampling iterations after all images have
      !   been read and their energies computed
      if (params%do_nested_sampling .and. .not. loop%repeat_xyz) then
         if (smp%i_nested == 0) then
            loop%md_istep = -1
            params%write_xyz = params%md_nsteps
            params%do_md = .true.
            if (comm%rank == 0) then
               write (*, *) '                                       |'
               write (*, *) 'Running nested sampling algorithm with |'
               write (*, '(1X,I6,A)') loop%n_xyz, ' walkers.                        |'
               write (*, *) '                                       |'
               write (*, *) 'Target pressure in nested sampling:    |'
               write (*, '(A,ES15.7,A)') ' P = ', params%p_nested, ' bar.               |'
               write (*, *) '                                       |'
               write (*, *) '[P = 0 means total energy, rather than |'
               write (*, *) 'total enthalphy, simulation]           |'
            end if
         end if
         !     At the end of the MD/MC moves we add the image to the pool if its energy has decreased
         if (loop%md_istep == params%md_nsteps) then
            loop%md_istep = -1
            state%velocities = 0.d0
            !       Unit cell volume
            state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                                     state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
            !       We check enthalpy, not internal energy (they are the same for P = 0)
            if (res%energy + dyn%E_kinetic + params%p_nested/dyn%eVperA3tobar*state%v_uc < smp%e_max) then
               call from_properties_to_image(smp%images(smp%i_image), state%positions, state%velocities, state%masses, &
                                             res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                             res%energy_exp, dyn%E_kinetic, &
                                             state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                             state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                             res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)
            end if
         end if
         !     This selects the highest energy image from the pool
         if (loop%md_istep == -1 .and. smp%i_nested < params%n_nested) then
            smp%i_nested = smp%i_nested + 1
            nl%rebuild_neighbors_list = .true.
            smp%i_max = 0
            smp%e_max = -1.d100
            do i = 1, loop%n_xyz
               state%v_uc = dot_product(cross_product(smp%images(i)%a_box, smp%images(i)%b_box), smp%images(i)%c_box)/ &
                            (dfloat(smp%images(i)%indices(1)*smp%images(i)%indices(2)*smp%images(i)%indices(3)))
               !         We check enthalpy, not potential energy (they are the same for P = 0)
               if (smp%images(i)%energy + smp%images(i)%e_kin + params%p_nested/dyn%eVperA3tobar*state%v_uc > smp%e_max) then
                  smp%e_max = smp%images(i)%energy + smp%images(i)%e_kin + params%p_nested/dyn%eVperA3tobar*state%v_uc
                  smp%i_max = i
               end if
            end do
            smp%i_image = smp%i_max
            deallocate (state%positions, state%velocities, state%masses, res%forces, state%species, &
                        state%species_supercell, state%fix_atom, state%xyz_species, state%xyz_species_supercell)
            !       Make a copy of a randonmly chosen image which is not i_image
            if (loop%n_xyz == 1) then
               i = smp%i_image
            else
               i = smp%i_image
!              random_number, which random_seed reaches; irand was seeded from the clock.
               do while (i == smp%i_image)
                  call random_number(u)
                  i = int(u*loop%n_xyz) + 1
               end do
            end if
            if (comm%rank == 0) then
               loop%counter = 1
               write (*, *) '                                       |'
               write (*, '(A,I8,A,I8,A)') "Nested sampling iter.:", smp%i_nested, "/", params%n_nested, " |"
               write (*, '(A,I8,A)') " - Highest enthalpy walker:    ", smp%i_image, " |"
               write (*, '(A,I8,A)') " - Walker selected for cloning:", i, " |"
               write (*, '(A,F15.7,A)') " - Max. enthalpy: ", smp%e_max, " eV |"
            end if
            call from_image_to_properties(smp%images(i), state%positions, state%velocities, state%masses, &
                                          res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                          res%energy_exp, dyn%E_kinetic, &
                                          state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                          state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                          res%local_dipoles, res%energies_dipole, res%dipole)
            state%v_uc = dot_product(cross_product(smp%images(i)%a_box, smp%images(i)%b_box), smp%images(i)%c_box)/ &
                         (dfloat(smp%images(i)%indices(1)*smp%images(i)%indices(2)*smp%images(i)%indices(3)))
            !       This only gets triggered if we are doing box rescaling, i.e., if the target nested sampling pressure (*not* the
            !       actual pressure for the atomic configuration) is > 0
!!!!!!!!!!!!!!!!!!!!!!!!!! Temporary hack
            if (params%scale_box_nested) then
               params%scale_box = .true.
               call random_number(smp%rand_scale)
!!!!!!!!!!!!!!! The size of the scaling should also decrease as we reach convergence (otherwise all trial moves will be rejected)
!!!!!!!!!!!!!!! Finally, there should be a limit for the acceptable aspect ratio of the simulation box
               smp%rand_scale = 2.d0*(smp%rand_scale - 0.5d0)*params%nested_max_strain
               params%box_scaling_factor = reshape([1.d0 + smp%rand_scale(1), smp%rand_scale(6)/2.d0, smp%rand_scale(5)/2.d0, &
                                                    smp%rand_scale(6)/2.d0, 1.d0 + smp%rand_scale(2), smp%rand_scale(4)/2.d0, &
                                                    smp%rand_scale(5)/2.d0, smp%rand_scale(4)/2.d0, 1.d0 + smp%rand_scale(3)], &
                                                   [3, 3])
               ! Make the transformation volume-preserving
               call volume_preserving_strain_transformation(state%a_box, state%b_box, state%c_box, params%box_scaling_factor)
               ! Volume scaling
               call get_ns_unbiased_volume_proposal(1.d0 - params%nested_max_volume_change, &
                                                    1.d0 + params%nested_max_volume_change, state%n_sites, smp%rand)
               params%box_scaling_factor = params%box_scaling_factor*(smp%rand)**(1.d0/3.d0)
               ! Each MPI process has a different set of random numbers so we need to broadcast
               call comm_bcast(comm, params%box_scaling_factor, 9)
            end if
            !       This is the so-called total enthalpy Hamiltonian Montecarlo approach (with physical masses)
            !       We do not need to broadcast the velocities here since they get broadcasted later on; otherwise
            !       we would have to do it since each MPI rank may see a different random number
            call random_number(state%velocities)
            call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), state%masses(1:state%n_sites))
            dyn%e_kin = 0.d0
            do i = 1, state%n_sites
               dyn%e_kin = dyn%e_kin + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), state%velocities(1:3, i))
            end do
            call random_number(smp%rand)
            state%velocities = state%velocities/sqrt(dyn%e_kin)*sqrt(smp%rand*(smp%e_max - res%energy - &
                                                                               params%p_nested/dyn%eVperA3tobar*state%v_uc))
         else if (smp%i_nested == params%n_nested) then
            loop%exit_loop = .true.
         end if
      end if
   end subroutine nested_step

!  One Monte Carlo step on rank 0: accept or reject the trial just evaluated,
!  log it, write the current configuration, and propose the next move.
   subroutine mc_step(smp, state, res, dyn, nl, model, params, perform, loop, comm, time)
      type(sampling_t), intent(inout) :: smp
      type(state_t), intent(inout) :: state
      type(results_t), intent(inout) :: res
      type(dynamics_t), intent(inout) :: dyn
      type(neighbors_t), intent(inout) :: nl
      type(model_t), intent(in) :: model
      type(input_parameters), intent(inout) :: params
      type(perform_t), intent(in) :: perform
      type(loop_t), intent(inout) :: loop
      type(comm_t), intent(in) :: comm
      type(times_t), intent(inout) :: time
      character*1024 :: string
      character*1024 :: temp_string
      character*1024 :: temp_string2
      integer :: i
      integer :: j

      if (comm%rank == 0) then

         if (params%do_mc) then
            if (loop%mc_istep == params%mc_nsteps) then
               loop%exit_loop = .true.
            else
               loop%exit_loop = .false.
            end if

            if (.not. loop%exit_loop .and. ( &
                (loop%md_istep == -1) .or. &
                (params%do_md .and. ( &
                 (loop%md_istep == params%md_nsteps) .or. &
                 ((abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites)) .and. (maxval(abs(res%forces)) < &
                                                                                                  params%f_tol)) &
                 )))) then
               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               call time_start(time%mc)

               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               if (loop%mc_istep > 0) then
                  !       Evaluate the conditions for acceptance
                  !       > We have the mc conditions in mc.f90
                  !       > We care about comparing e_store to the energy of the new configuration based on the mc_movw

                  ! Reset the parameters for md / relaxation
                  smp%trial_came_from_md = params%do_md
                  if (params%do_md) then
                     loop%md_istep = -1
                     params%do_md = .false.
                     smp%do_mc_relax = .false.
                     ! Assume that the number of steps has already been set.
                  end if

                  if (.not. params%mc_hamiltonian) dyn%E_kinetic = 0.d0

!                 The trial configuration has to be the one `energy` belongs to.
!                 An "md" move, or a relaxation after any other move, leaves
!                 `positions` one integrator step *past* the last force
!                 evaluation: compute_md advances them after the energy was
!                 computed, and stashes the configuration it was computed at in
!                 positions_prev. md.f90 says as much -- "velocities and
!                 positions_prev are synchronous, positions is dt ahead of
!                 velocities" -- and compute_md writes positions_prev to the
!                 trajectory for exactly this reason.
!
!                 Storing `positions` here accepted or rejected x_(n+1) on the
!                 strength of E(x_n), and wrote a frame to mc_all.xyz whose
!                 energy was not the energy of its own coordinates: ~1 eV out on
!                 512 atoms after a 0.5 fs velocity-Verlet burst. Rewind by the
!                 one step, which also pairs the stored positions with the
!                 stored velocities.
                  if (smp%trial_came_from_md) then
                     state%positions(1:3, 1:state%n_sites) = state%positions_prev(1:3, 1:state%n_sites)
                  end if

                  call from_properties_to_image(smp%images(smp%i_trial_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, dyn%E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

                  if (params%verb > 50) write (*, *) '.......................................|'
                  if (params%verb > 50) write (*, '(A,1X,I0)') ' MC Iteration:', loop%mc_istep
                  if (params%verb > 50) write (*, '(A,1X,A)') '    Move type:', smp%mc_move

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_prev:', smp%images(smp%i_current_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_prev:', smp%images(smp%i_current_image)&
                       &%energy + smp%images(smp%i_current_image)%e_kin

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_new:', smp%images(smp%i_trial_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_new :', smp%images(smp%i_trial_image)%energy &
                       &+ smp%images(smp%i_trial_image)%e_kin

                  state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                                           state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))

                  if (params%accessible_volume) then
                     call get_accessible_volume(state%v_uc, smp%v_a_uc, state%species, params%radii)
                     if (params%verb > 50) write (*, '(A,F12.6,A,F12.6&
                          &,1X,A)') ' V_acc new: ', smp%v_a_uc, ' A^3&
                          & V_acc old ', smp%v_a_uc_prev, 'A^3 |'
                  else
                     smp%v_a_uc = state%v_uc
                  end if

                  call get_mc_acceptance(smp%mc_move, smp%p_accept, &
                       res%energy + dyn%E_kinetic, &
                       smp%images(smp%i_current_image)%energy + smp%images(smp%i_current_image)%e_kin, &
                       params%t_beg, smp%mc_mu_id, &
                       params%mc_mu, smp%n_mc_species, state%v_uc, smp%v_uc_prev,&
                       & smp%v_a_uc, smp%v_a_uc_prev, params%mc_exchange_mass, &
                       & params%mc_exchange_e0, params%mc_mu_reference, &
                       & params%p_beg, state%n_sites)

!                 call get_mc_acceptance(mc_move, p_accept, &
!                      energy + E_kinetic, &
!                      images(i_current_image)%energy + images(i_current_image)%e_kin, &
!                      params%t_beg, &
!                      params%mc_mu(mc_mu_id), n_mc_species(mc_mu_id), v_uc, v_uc_prev,&
!                      & v_a_uc, v_a_uc_prev, params&
!                      &%masses_types(mc_id(mc_mu_id)), params%p_beg)

                  call random_number(smp%ranf)

                  if (smp%mc_move == "insertion") smp%n_mc_species(smp%mc_mu_id) = smp%n_mc_species(smp%mc_mu_id) + 1
                  if (smp%mc_move == "removal") smp%n_mc_species(smp%mc_mu_id) = smp%n_mc_species(smp%mc_mu_id) - 1

                  !    ACCEPT OR REJECT
                  if (params%verb > 50) write (*, '(A,1X,A,1X,A,L4,1X&
                       &,A,ES12.6,1X,A,1X,ES12.6)') 'Is ',&
                       & trim(smp%mc_move), 'accepted?', smp%p_accept >&
                       & smp%ranf, ' p_accept =', smp%p_accept, ' ranf = ',&
                       & smp%ranf

                  if (loop%mc_istep == 1) then
                     open (unit=200, file="mc.log", status="unknown")
                     if (res%energy_exp > 0.d0) then
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current             E_exp_trial&
                             &          E_exp_current  N_tot_trial &
                             & N_mc_species_trial'
                     else
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current &
                             &          N_tot_trial  N_mc_species_trial'
                     end if

                  end if
                  if (loop%mc_istep > 1) then
                     open (unit=200, file="mc.log", status="old", position="append")
                  end if

                  ! collect the strings for the species etc
                  temp_string = ""
                  temp_string2 = ""

                  do i = 1, params%n_mc_mu
                     temp_string = ""
                     write (temp_string, "(A,1X,I8)") trim(params%mc_species(i)), smp%n_mc_species(i)
                     temp_string2 = trim(temp_string2)//" "//trim(temp_string)
                  end do

                  if (res%energy_exp > 0.d0) then

                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                          loop%mc_istep, trim(adjustl(smp%mc_move)), smp%p_accept > smp%ranf, res%energy + dyn%E_kinetic, &
                          smp%images(smp%i_current_image)%energy +&
                          & smp%images(smp%i_current_image)%e_kin, res%energy_exp,&
                          & smp%images(smp%i_current_image)%energy_exp,&
                          & smp%images(smp%i_trial_image)%n_sites,&
                          & trim(temp_string2)
                  else
                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                        loop%mc_istep, trim(adjustl(smp%mc_move)), smp%p_accept > smp%ranf, res%energy + dyn%E_kinetic, &
                        smp%images(smp%i_current_image)%energy + smp%images(smp%i_current_image)%e_kin, &
                        smp%images(smp%i_trial_image)%n_sites, trim(temp_string2)

                  end if

                  if (loop%mc_istep >= 1) close (200)

                  if (smp%p_accept > smp%ranf) then
                     !             Accept
                     ! Set variables
                     loop%n_sites_prev = state%n_sites
                     smp%v_uc_prev = state%v_uc
                     smp%v_a_uc_prev = smp%v_a_uc
                     smp%virial_prev = res%virial
                     !   Assigning the default image with the accepted one
                     smp%images(smp%i_current_image) = smp%images(smp%i_trial_image)

                     if (params%n_mc_mu > 0) then
                        smp%n_mc_species_prev = smp%n_mc_species
                     end if

                  end if
                  if (state%n_sites > 1) then
                     dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
                     dyn%instant_pressure = (dyn%kB*dfloat(state%n_sites - 1)*dyn%instant_temp&
                          &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                          &/state%v_uc*dyn%eVperA3tobar
                  else
                     dyn%instant_temp = 0.0d0
                     dyn%instant_pressure = 0.0d0
                  end if

                  if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') '&
                          & Writing mc_current.xyz and&
                          & mc_all.xyz '
                     call wrap_pbc(smp%images(smp%i_current_image)&
                          &%positions(1:3,&
                          & 1:smp%images(smp%i_current_image)%n_sites),&
                          & smp%images(smp%i_current_image)%a_box&
                          &/dfloat(state%indices(1)),&
                          & smp%images(smp%i_current_image)%b_box&
                          &/dfloat(state%indices(2)),&
                          & smp%images(smp%i_current_image)%c_box&
                          &/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf,&
                          & params%valid_xrd, params%valid_nd,&
                          & params%do_pair_distribution, params&
                          &%do_structure_factor, params%do_xrd,&
                          & params%do_nd, string, params%do_dipole,&
                          & smp%images(smp%i_current_image)%dipole,&
                          & smp%images(smp%i_current_image)%energies_dipole)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 0, 1.0d0, 0.d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 1, 1.0d0, 0.d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & smp%images(smp%i_current_image)%local_properties,&
                          & smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .false., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                  end if

                  !          Add acceptance to the log file else dont
                  call time_end(time%mc)

               else ! if (mc_istep == 0)
                  smp%temp_md_nsteps = params%md_nsteps
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, *) 'Starting MC, using parameters:         |'
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nsteps     = ', params%mc_nsteps, '     &
                       &        |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_types    = ', params%n_mc_types, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_types:                              |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '     ', params%mc_types(i), '|'
                  end do
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_accept_ratio:                       |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,F12.8,1X&
                          &,A)') '   ', params%mc_acceptance(i), '    &
                          &                  |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_swaps    = ', params%n_mc_swaps, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_swaps:  &
                       &                            |'
                  do i = 1, 2*params%n_mc_swaps
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '   ', params%mc_swaps(i), '              &
                          &        |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_move_max   = ', params%mc_move_max, &
                       & 'A   |'

                  do i = 1, params%n_mc_mu
                     write (*, '(1X,A,1X,F17.8,1X,A)') 'mc_mu         = ', params%mc_mu(1), 'eV  |'
                     write (*, '(1X,A,1X,A,1X,A)') 'mc_species    = ', trim(params%mc_species(i)), '                    |'
                  end do

                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_min_dist   = ', params%mc_min_dist, &
                       & 'A   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_lnvol_max  = ', params%mc_lnvol_max,&
                       & '    |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_write_xyz  = ', params%mc_write_xyz, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_relax      = ', params%mc_relax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nrelax     = ', params%mc_nrelax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_relax_opt  = ', params%mc_relax_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_hybrid_opt = ', params%mc_hybrid_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_optimize_exp = ', params%mc_optimize_exp&
                       &, '  |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_hamiltonian = ', params%mc_hamiltonian, '&
                       &  |'

                  if (params%verb > 50) write (*, *) '                                       |'
                  ! t_beg must

                  if (.not. allocated(smp%images) .and. .not. params%do_nested_sampling) then
                     allocate (smp%images(1:2))
                  else if (.not. allocated(smp%images) .and. params%do_nested_sampling) then
                     allocate (smp%images(1:2*smp%i_image))
                  end if

                  if (.not. allocated(smp%mc_mol_id)) then
                     allocate (smp%mc_mol_id(1:state%n_sites), smp%mc_mol_mu(1:state%n_sites))
                     smp%mc_mol_id = 0
                     smp%mc_mol_mu = 0
                  end if

                  if (.not. allocated(smp%mc_id) .and. params%n_mc_mu > 0) then
                     allocate (smp%mc_id(1:params%n_mc_mu))
                     allocate (smp%n_mc_species(1:params%n_mc_mu))
                     allocate (smp%n_mc_species_prev(1:params%n_mc_mu))

                     smp%mc_id = 1
                     smp%n_mc_species = 0

                     !    get the mc species types

                     do j = 1, params%n_mc_mu
                        do i = 1, model%n_species
                           if (params%species_types(i) == params%mc_species(j)) then
                              smp%mc_id(j) = i
                           end if
                        end do
                     end do
                  end if

                  !       Now use the image construct to store this as the image to compare to
                  call from_properties_to_image(smp%images(smp%i_current_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, dyn%E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

                  dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
                  dyn%instant_pressure = (dyn%kB*dfloat(state%n_sites - 1)*dyn%instant_temp&
                       &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                       &/state%v_uc*dyn%eVperA3tobar

                  if ((loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') ' Writing mc_current.xyz and mc_all.xyz '
                     call wrap_pbc(smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                                   smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                                   smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                                   smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                          & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                          & params%do_dipole, smp%images(smp%i_current_image)%dipole,&
                          & smp%images(smp%i_current_image)%energies_dipole)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 0, 1.0d0, 0.0d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 1, 1.0d0, 0.0d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     smp%v_uc_prev = dot_product(cross_product(state%a_box, state%b_box), &
                                                 state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
                     if (params%accessible_volume) then
                        call get_accessible_volume(smp%v_uc_prev, smp%v_a_uc_prev, state%species, params%radii)
                     else
                        smp%v_a_uc_prev = smp%v_uc_prev
                     end if
                  end if

               end if

               !  Now start the mc logic: first, use the stored images properties
               call from_image_to_properties(smp%images(smp%i_current_image), state%positions, state%velocities, state%masses, &
                                             res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                             res%energy_exp, dyn%E_kinetic, &
                                             state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                             state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                             res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

               call perform_mc_step(&
                    & state%positions, state%species, state%xyz_species, state%masses, state%fix_atom,&
                    & state%velocities, state%positions_prev, state%positions_diff, smp%disp, smp%d_disp, &
                       params%n_local_properties,&
                    & params%mc_acceptance, params%mc_mu_acceptance, res%local_properties, &
                    smp%images(smp%i_current_image)%local_properties, res%energies,&
                    & res%forces, state%forces_prev, state%n_sites, params%n_mc_mu, smp%mc_mu_id, smp%n_mc_species,&
                    & smp%mc_move, params%mc_species,&
                    & params%mc_move_max, params%mc_min_dist, params%mc_max_dist, params%mc_max_insertion_trials, &
                    params%mc_lnvol_max, params%mc_types, params%masses_types, smp%species_idx,&
                    & smp%images(smp%i_current_image)%positions,&
                    & smp%images(smp%i_current_image)%species,&
                    & smp%images(smp%i_current_image)%xyz_species,&
                    & smp%images(smp%i_current_image)%fix_atom,&
                    & smp%images(smp%i_current_image)%masses, state%a_box(1:3), state%b_box(1:3),&
                    & state%c_box(1:3), state%indices, params%do_md, params%mc_relax,&
                    & loop%md_istep, smp%mc_id, dyn%E_kinetic, dyn%instant_temp, params%t_beg,&
                    & params%n_mc_swaps, params%mc_swaps, params%mc_swaps_id, &
                    & params%species_types, params%mc_hamiltonian,&
                    & params%n_mc_relax_after, params&
                    &%mc_relax_after, smp%do_mc_relax, params%verb, &
                    params%mc_n_planes, params%mc_planes, params%mc_max_dist_to_planes, &
                    params%mc_planes_restrict_to_polyhedron, &
                    params%mc_molecules, smp%mc_mol_id, smp%mc_mol_mu, &
                    smp%images(smp%i_current_image)%mc_mol_id, smp%images(smp%i_current_image)%mc_mol_mu, smp%mc_mol_next)

               nl%rebuild_neighbors_list = .true.

               ! NOTE: the species_supercell and xyz_species_supercell are
               ! not commensurate with the new image as these have not been
               ! calculated. If reading from an outputted xyz file, then it
               ! should be okay but really the new atoms should be added to
               ! the supercell in the usual way, but for convenience, one has
               ! not done that.

               if (params%mc_relax .and. smp%do_mc_relax) then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_relax_opt
                  params%md_nsteps = params%mc_nrelax

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, dyn%E_kinetic, state%masses, dyn%instant_temp, &
                                            params%t_beg, &
                                            params%velocity_distribution)

                  if (params%mc_hamiltonian) dyn%E_kinetic_prev = dyn%E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if
               ! If doing md, don't relax
               if (smp%mc_move == 'md') then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_hybrid_opt
                  params%md_nsteps = smp%temp_md_nsteps

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, dyn%E_kinetic, state%masses, dyn%instant_temp, &
                                            params%t_beg, &
                                            params%velocity_distribution)
                  if (params%mc_hamiltonian) dyn%E_kinetic_prev = dyn%E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if

               if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                    modulo(loop%mc_istep, params%write_xyz) == 0)) then

                  call wrap_pbc(state%positions(1:3, 1:state%n_sites), &
                                state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                                state%c_box/dfloat(state%indices(3)))
                  call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                       & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                       &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                       & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                       & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                  call write_extxyz(state%n_sites, 0, 1.0d0, 0.0d0, dyn%instant_temp, dyn%instant_pressure, &
                       state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                          state%c_box/dfloat(state%indices(3)), &
                       res%virial, state%xyz_species, &
                       state%positions(1:3, 1:state%n_sites), state%velocities, &
                       res%forces, res%energies(1:state%n_sites), state%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & model%local_property_labels, res%local_properties&
                       &, state%fix_atom, smp%mc_file, string, .true.)
               end if
               ! As we have moved/added/removed, we must check the supercell and  broadcast the results

               call read_xyz(smp%mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .true., state%fix_atom, params%t_beg, &
                             params%write_array_property(6), .true., params%randomize_velocities)

            else
               if (smp%mc_move == 'md') then
                  if (params%print_progress .and. loop%md_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%md_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%md_istep == params%md_nsteps - 1 .or. &
                           (abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites) .and. &
                            maxval(abs(res%forces)) < params%f_tol) .and. loop%md_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%md_istep < params%md_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%md_istep + 1)/params%md_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%md_istep + 1)/params%md_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%mc_hamiltonian) then
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step: H =&
                          & T + V = ", res%energy + dyn%E_kinetic, ",&
                          & iteration ", loop%md_istep, "/", params&
                          &%md_nsteps
                  else
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step:&
                          & energy = ", res%energy, ", iteration ",&
                          & loop%md_istep, "/", params%md_nsteps
                  end if

                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') '&
                       & core_pot energy:', sum(res%energies_core_pot),&
                       & 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') '&
                       & vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               else
                  if (params%print_progress .and. loop%mc_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%mc_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%mc_istep == params%mc_nsteps - 1 .and. loop%mc_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%mc_istep < params%mc_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%mc_istep + 1)/params%mc_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%mc_istep + 1)/params%mc_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%verb > 50 .and. smp%do_mc_relax) write (*, '(1X,A,1X,F20.8,1X,A&
                       &,1X,I8,1X,A,1X,I8)') "MC Relax md step: energy &
                       &= ", res%energy, ", iteration ", loop%md_istep, "/",&
                       & params%mc_nrelax
                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               end if
            end if
         end if

      end if
   end subroutine mc_step

end module turbogap_sampling
