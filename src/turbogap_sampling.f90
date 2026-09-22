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
   use types, only: image, input_parameters, from_properties_to_image, from_image_to_properties
   use md, only: remove_cm_vel, get_ns_unbiased_volume_proposal, volume_preserving_strain_transformation
   use soap_turbo_functions, only: cross_product
   use turbogap_comm, only: comm_t, comm_bcast
   use turbogap_results, only: results_t
   use turbogap_domain, only: neighbors_t
   use turbogap_structure, only: state_t
   use turbogap_loop, only: loop_t
   use turbogap_md, only: dynamics_t

   implicit none

   private
   public :: mc_prepare_step
   public :: nested_step

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
               do while (i == smp%i_image)
                  i = mod(irand(), loop%n_xyz) + 1
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

end module turbogap_sampling
