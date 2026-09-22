! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_results.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  What one evaluation of the potential produces: energies, forces and virials
!  per contribution family and in total, and the per-site properties computed
!  on the way. The this_ arrays hold one rank's partial sums before the
!  reduction. Declare the variable with target: the contribution table and the
!  local-property pointers point into it.
module turbogap_results

   use kinds, only: dp
   use types, only: input_parameters, perform_t, any_has_local_properties
   use turbogap_structure, only: state_t
   use turbogap_setup, only: model_t
   use turbogap_loop, only: loop_t
   use exp_utils, only: exp_dissimilarity, exp_dissim_ref

   implicit none

   private
   public :: results_free
   public :: results_prepare

   type, public :: results_t
      real(dp) :: energy
      real(dp) :: energy_prev
      real(dp) :: energy_exp
      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: forces(:, :)
      real(dp) :: virial(1:3, 1:3)
      real(dp), allocatable :: this_energies(:)
      real(dp), allocatable :: this_forces(:, :)
      real(dp) :: this_virial(1:3, 1:3)
      real(dp), allocatable :: energies_exp(:)

      real(dp), allocatable :: energies_soap(:)
      real(dp), allocatable :: forces_soap(:, :)
      real(dp) :: virial_soap(1:3, 1:3)
      real(dp), allocatable :: energies_2b(:)
      real(dp), allocatable :: forces_2b(:, :)
      real(dp) :: virial_2b(1:3, 1:3)
      real(dp), allocatable :: energies_3b(:)
      real(dp), allocatable :: forces_3b(:, :)
      real(dp) :: virial_3b(1:3, 1:3)
      real(dp), allocatable :: energies_core_pot(:)
      real(dp), allocatable :: forces_core_pot(:, :)
      real(dp) :: virial_core_pot(1:3, 1:3)

      real(dp), allocatable :: energies_vdw(:)
      real(dp), allocatable :: forces_vdw(:, :)
      real(dp) :: virial_vdw(1:3, 1:3)
      real(dp), allocatable :: this_energies_vdw(:)
      real(dp), allocatable :: this_forces_vdw(:, :)
      real(dp) :: this_virial_vdw(1:3, 1:3)
      real(dp), allocatable :: energies_estat(:)
      real(dp), allocatable :: forces_estat(:, :)
      real(dp) :: virial_estat(1:3, 1:3)
      real(dp), allocatable :: this_energies_estat(:)
      real(dp), allocatable :: this_forces_estat(:, :)
      real(dp) :: this_virial_estat(1:3, 1:3)
      real(dp), allocatable :: energies_lp(:)
      real(dp), allocatable :: forces_lp(:, :)
      real(dp) :: virial_lp(1:3, 1:3)
      real(dp), allocatable :: this_energies_lp(:)
      real(dp), allocatable :: this_forces_lp(:, :)
      real(dp) :: this_virial_lp(1:3, 1:3)
      real(dp), allocatable :: energies_pdf(:)
      real(dp), allocatable :: forces_pdf(:, :)
      real(dp) :: virial_pdf(1:3, 1:3)
      real(dp), allocatable :: this_energies_pdf(:)
      real(dp), allocatable :: this_forces_pdf(:, :)
      real(dp) :: this_virial_pdf(1:3, 1:3)
      real(dp), allocatable :: energies_sf(:)
      real(dp), allocatable :: forces_sf(:, :)
      real(dp) :: virial_sf(1:3, 1:3)
      real(dp), allocatable :: this_energies_sf(:)
      real(dp), allocatable :: this_forces_sf(:, :)
      real(dp) :: this_virial_sf(1:3, 1:3)
      real(dp), allocatable :: energies_xrd(:)
      real(dp), allocatable :: forces_xrd(:, :)
      real(dp) :: virial_xrd(1:3, 1:3)
      real(dp), allocatable :: this_energies_xrd(:)
      real(dp), allocatable :: this_forces_xrd(:, :)
      real(dp) :: this_virial_xrd(1:3, 1:3)
      real(dp), allocatable :: energies_nd(:)
      real(dp), allocatable :: forces_nd(:, :)
      real(dp) :: virial_nd(1:3, 1:3)
      real(dp), allocatable :: this_energies_nd(:)
      real(dp), allocatable :: this_forces_nd(:, :)
      real(dp) :: this_virial_nd(1:3, 1:3)

      real(dp), allocatable :: mbd_ts_scaling(:)
      real(dp), allocatable :: this_mbd_ts_scaling(:)
      real(dp), allocatable :: local_virial_vdw_diag(:, :)
      real(dp), allocatable :: local_virial_vdw_diag_corr(:, :)
      real(dp), allocatable :: this_local_virial_vdw_diag(:, :)
      real(dp), allocatable :: energies_vdw_corr(:)
      real(dp), allocatable :: forces_vdw_corr(:, :)
      logical :: update_mbd_ts_scaling = .true.

      real(dp), allocatable :: local_properties(:, :)
      real(dp), allocatable :: local_properties_cart_der(:, :, :)
      real(dp), allocatable :: this_local_properties(:, :)
      real(dp), allocatable :: this_local_properties_cart_der(:, :, :)
      real(dp), pointer :: this_local_properties_pt(:, :)
      real(dp), pointer :: this_local_properties_cart_der_pt(:, :, :)

!     The dipole model's fictitious scalar is reported as energy_dipole and is
!     deliberately absent from energies; the model has no forces or virial.
      real(dp), allocatable :: local_dipoles(:, :)
      real(dp), allocatable :: this_local_dipoles(:, :)
      real(dp), allocatable :: energies_dipole(:)
      real(dp), allocatable :: this_energies_dipole(:)
      real(dp) :: dipole(1:3)
   end type results_t

contains

!  Size the result arrays for this structure, reallocating only when the site
!  count changed (or under MC, which can change it silently), and zero what this
!  evaluation accumulates into. RESIZED says whether the arrays were rebuilt.
   subroutine results_prepare(res, state, params, model, perform, loop, n_pairs, n_pairs_prev, resized)
      type(results_t), target, intent(inout) :: res
      type(state_t), intent(in) :: state
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      type(perform_t), intent(in) :: perform
      type(loop_t), intent(in) :: loop
      integer, intent(in) :: n_pairs
      integer, intent(in) :: n_pairs_prev
      logical, intent(out) :: resized

      resized = (state%n_sites /= loop%n_sites_prev .or. params%do_mc)

      !     We only need to reallocate the arrays if the number of sites changes
      ! REMOVE TRUE FROM IF STATEMENT
      if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
         if (allocated(res%energies)) deallocate (res%energies, &
                                                  res%energies_soap, &
                                                  res%energies_2b, &
                                                  res%energies_3b, &
                                                  res%energies_core_pot, &
                                                  res%this_energies, &
                                                  res%energies_vdw, &
                                                  res%energies_vdw_corr, &
                                                  res%mbd_ts_scaling, &
                                                  res%this_forces, &
                                                  res%energies_lp, &
                                                  res%energies_exp, &
                                                  res%energies_estat, &
                                                  res%this_mbd_ts_scaling)
         allocate (res%energies(1:state%n_sites))
         allocate (res%this_energies(1:state%n_sites))
         allocate (res%energies_soap(1:state%n_sites))
         allocate (res%energies_2b(1:state%n_sites))
         allocate (res%energies_3b(1:state%n_sites))
         allocate (res%energies_core_pot(1:state%n_sites))
         allocate (res%energies_vdw(1:state%n_sites))
         allocate (res%energies_vdw_corr(1:state%n_sites))
         allocate (res%energies_lp(1:state%n_sites))
         allocate (res%energies_estat(1:state%n_sites))
         allocate (res%energies_exp(1:state%n_sites))
!          We do this allocations for van der Waals corrections
         allocate (res%mbd_ts_scaling(1:state%n_sites))
         allocate (res%this_mbd_ts_scaling(1:state%n_sites))
!          Allocated whether or not a dipole model is loaded: they are passed to
!          get_gap_soap unconditionally, and 4 doubles per atom is not worth a
!          second code path.
         if (allocated(res%local_dipoles)) deallocate (res%local_dipoles, res%this_local_dipoles, &
                                                       res%energies_dipole, res%this_energies_dipole)
         allocate (res%local_dipoles(1:3, 1:state%n_sites))
         allocate (res%this_local_dipoles(1:3, 1:state%n_sites))
         allocate (res%energies_dipole(1:state%n_sites))
         allocate (res%this_energies_dipole(1:state%n_sites))

         if (perform%pdf) then
            if (allocated(res%energies_pdf)) deallocate (res%energies_pdf)
            allocate (res%energies_pdf(1:state%n_sites))
         end if

         if (perform%sf) then
            if (allocated(res%energies_sf)) deallocate (res%energies_sf)
            allocate (res%energies_sf(1:state%n_sites))
         end if

         if (perform%xrd) then
            if (allocated(res%energies_xrd)) deallocate (res%energies_xrd)
            allocate (res%energies_xrd(1:state%n_sites))
         end if

         if (perform%nd) then
            if (allocated(res%energies_nd)) deallocate (res%energies_nd)
            allocate (res%energies_nd(1:state%n_sites))
         end if

         !       This needs to be allocated even if no force prediction is needed:
         allocate (res%this_forces(1:3, 1:state%n_sites))
      end if
      res%energies = 0.d0
      res%energies_soap = 0.d0
      res%energies_2b = 0.d0
      res%energies_3b = 0.d0
      res%energies_core_pot = 0.d0
      res%energies_vdw = 0.d0
      res%energies_estat = 0.d0
      res%energies_lp = 0.d0
      res%energies_exp = 0.d0
!        The dissimilarity accumulators belong to the same step as energies_exp
!        and are zeroed with it. This is the only point that knows a new
!        evaluation has begun; get_exp_energies is called once per observable
!        and mad_ir separately again, so neither can reset them itself.
      exp_dissimilarity = 0.d0
      exp_dissim_ref = 0.d0
      res%local_dipoles = 0.d0
      res%energies_dipole = 0.d0
      res%dipole = 0.d0

      if (perform%pdf) res%energies_pdf = 0.d0
      if (perform%sf) res%energies_sf = 0.d0
      if (perform%xrd) res%energies_xrd = 0.d0
      if (perform%nd) res%energies_nd = 0.d0

      ! Adding allocation of local properties

      ! Now one could use pointers such that hirshfeld_v(:) acts as an alias for local_properties(vdw_index,:)...
      if (any_has_local_properties(model%soap_turbo_hypers)) then
         if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
            if (allocated(res%local_properties)) then
               nullify (res%this_local_properties_pt)
               deallocate (res%this_local_properties, res%local_properties)
               if (params%do_forces) then
                  nullify (res%this_local_properties_cart_der_pt)
                  deallocate (res%this_local_properties_cart_der, res%local_properties_cart_der)
               end if
            end if
            allocate (res%local_properties(1:state%n_sites, 1:params%n_local_properties))
            allocate (res%this_local_properties(1:state%n_sites, 1:params%n_local_properties))
            res%this_local_properties_pt => res%this_local_properties

            !         I don't remember why this needs a pointer <----------------------------------------- CHECK

         end if
         res%local_properties = 0.d0

         if (params%do_forces) then
            if (n_pairs /= n_pairs_prev) then
               if (allocated(res%local_properties_cart_der)) deallocate (res%local_properties_cart_der, &
                                                                         res%this_local_properties_cart_der)
               allocate (res%local_properties_cart_der(1:3, 1:n_pairs, 1:params%n_local_properties))
               allocate (res%this_local_properties_cart_der(1:3, 1:n_pairs, &
                                                            1:params%n_local_properties))
            end if
            if (.not. allocated(res%local_properties_cart_der)) then
               allocate (res%local_properties_cart_der(1:3, 1:n_pairs, 1:params%n_local_properties))
               allocate (res%this_local_properties_cart_der(1:3, 1:n_pairs, &
                                                            1:params%n_local_properties))
            end if

            res%local_properties_cart_der = 0.d0
            res%this_local_properties_cart_der_pt =>&
                 & res%this_local_properties_cart_der(1:3,&
                 & 1:n_pairs, 1:params&
                 &%n_local_properties)
         end if
      end if

      ! Now go through the soap turbo hypers, and see if any are vdw or
      ! otherwise, if vdw, one can have pointers to point to the data
      ! structures such that it makes things clearer. One needs to check
      ! that this allocation still works iwth if(allocated(hirsh_v))
      ! statements

      if (params%do_forces) then
         if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
            if (allocated(res%forces)) deallocate (res%forces, res%forces_soap, res%forces_2b, res%forces_3b, &
               res%forces_core_pot, res%forces_vdw,&
                 & res%forces_lp, res%forces_estat, res%local_virial_vdw_diag, res%local_virial_vdw_diag_corr)
            allocate (res%forces(1:3, 1:state%n_sites))
            allocate (res%forces_soap(1:3, 1:state%n_sites))
            allocate (res%forces_2b(1:3, 1:state%n_sites))
            allocate (res%forces_3b(1:3, 1:state%n_sites))
            allocate (res%forces_core_pot(1:3, 1:state%n_sites))
            allocate (res%forces_vdw(1:3, 1:state%n_sites))
            if (allocated(res%forces_vdw_corr)) deallocate (res%forces_vdw_corr)
            allocate (res%forces_vdw_corr(1:3, 1:state%n_sites))
            allocate (res%forces_lp(1:3, 1:state%n_sites))
            allocate (res%forces_estat(1:3, 1:state%n_sites))
            allocate (res%local_virial_vdw_diag_corr(1:3, 1:state%n_sites))
            allocate (res%local_virial_vdw_diag(1:3, 1:state%n_sites))

            if (perform%pdf_forces) then
               if (allocated(res%forces_pdf)) deallocate (res%forces_pdf)
               allocate (res%forces_pdf(1:3, 1:state%n_sites))
            end if

            if (perform%sf_forces) then
               if (allocated(res%forces_sf)) deallocate (res%forces_sf)
               allocate (res%forces_sf(1:3, 1:state%n_sites))
            end if

            if (perform%xrd_forces) then
               if (allocated(res%forces_xrd)) deallocate (res%forces_xrd)
               allocate (res%forces_xrd(1:3, 1:state%n_sites))
            end if

            if (perform%nd_forces) then
               if (allocated(res%forces_nd)) deallocate (res%forces_nd)
               allocate (res%forces_nd(1:3, 1:state%n_sites))
            end if

         end if
         res%forces = 0.d0
         res%forces_soap = 0.d0
         res%forces_2b = 0.d0
         res%forces_3b = 0.d0
         res%forces_core_pot = 0.d0
         res%forces_vdw = 0.d0
         res%forces_estat = 0.d0
         res%forces_lp = 0.d0
         res%virial = 0.d0
         res%virial_soap = 0.d0
         res%virial_2b = 0.d0
         res%virial_3b = 0.d0
         res%virial_core_pot = 0.d0
         res%virial_vdw = 0.d0
         res%virial_estat = 0.d0
         res%virial_lp = 0.d0
         res%local_virial_vdw_diag = 0.d0
         if (perform%pdf_forces) then
            res%forces_pdf = 0.d0
            res%virial_pdf = 0.d0
            res%this_virial_pdf = 0.d0
         end if

         if (perform%sf_forces) then
            res%forces_sf = 0.d0
            res%virial_sf = 0.d0
            res%this_virial_sf = 0.d0
         end if

         if (perform%xrd_forces) then
            res%forces_xrd = 0.d0
            res%virial_xrd = 0.d0
            res%this_virial_xrd = 0.d0
         end if

         if (perform%nd_forces) then
            res%forces_nd = 0.d0
            res%virial_nd = 0.d0
            res%this_virial_nd = 0.d0
         end if
      end if
   end subroutine results_prepare

   subroutine results_free(res)
      type(results_t), intent(inout) :: res

      if (allocated(res%energies)) deallocate (res%energies)
      if (allocated(res%local_dipoles)) deallocate (res%local_dipoles, res%this_local_dipoles)
      if (allocated(res%energies_dipole)) deallocate (res%energies_dipole, res%this_energies_dipole)
      if (allocated(res%energies_soap)) deallocate (res%energies_soap)
      if (allocated(res%energies_2b)) deallocate (res%energies_2b)
      if (allocated(res%energies_3b)) deallocate (res%energies_3b)
      if (allocated(res%energies_core_pot)) deallocate (res%energies_core_pot)
      if (allocated(res%energies_vdw)) deallocate (res%energies_vdw)
      if (allocated(res%energies_exp)) deallocate (res%energies_exp)
      if (allocated(res%energies_lp)) deallocate (res%energies_lp)
      if (allocated(res%energies_pdf)) deallocate (res%energies_pdf)
      if (allocated(res%energies_sf)) deallocate (res%energies_sf)
      if (allocated(res%energies_xrd)) deallocate (res%energies_xrd)
      if (allocated(res%energies_nd)) deallocate (res%energies_nd)
      if (allocated(res%this_energies)) deallocate (res%this_energies)
      if (allocated(res%this_energies_vdw)) deallocate (res%this_energies_vdw)
      if (allocated(res%this_energies_lp)) deallocate (res%this_energies_lp)
      if (allocated(res%this_energies_pdf)) deallocate (res%this_energies_pdf)
      if (allocated(res%this_energies_sf)) deallocate (res%this_energies_sf)
      if (allocated(res%this_energies_xrd)) deallocate (res%this_energies_xrd)
      if (allocated(res%this_energies_nd)) deallocate (res%this_energies_nd)
      if (allocated(res%forces)) deallocate (res%forces)
      if (allocated(res%forces_soap)) deallocate (res%forces_soap)
      if (allocated(res%forces_2b)) deallocate (res%forces_2b)
      if (allocated(res%forces_3b)) deallocate (res%forces_3b)
      if (allocated(res%forces_core_pot)) deallocate (res%forces_core_pot)
      if (allocated(res%forces_vdw)) deallocate (res%forces_vdw)
      if (allocated(res%forces_lp)) deallocate (res%forces_lp)
      if (allocated(res%forces_pdf)) deallocate (res%forces_pdf)
      if (allocated(res%forces_sf)) deallocate (res%forces_sf)
      if (allocated(res%forces_xrd)) deallocate (res%forces_xrd)
      if (allocated(res%forces_nd)) deallocate (res%forces_nd)
      if (allocated(res%this_forces)) deallocate (res%this_forces)
      if (allocated(res%this_forces_vdw)) deallocate (res%this_forces_vdw)
      if (allocated(res%this_forces_lp)) deallocate (res%this_forces_lp)
      if (allocated(res%this_forces_pdf)) deallocate (res%this_forces_pdf)
      if (allocated(res%this_forces_sf)) deallocate (res%this_forces_sf)
      if (allocated(res%this_forces_xrd)) deallocate (res%this_forces_xrd)
      if (allocated(res%this_forces_nd)) deallocate (res%this_forces_nd)
      if (allocated(res%local_properties)) deallocate (res%local_properties)
      if (allocated(res%local_properties_cart_der)) deallocate (res%local_properties_cart_der)
      if (allocated(res%this_local_properties)) deallocate (res%this_local_properties)
      if (allocated(res%this_local_properties_cart_der)) deallocate (res%this_local_properties_cart_der)
   end subroutine results_free

end module turbogap_results
