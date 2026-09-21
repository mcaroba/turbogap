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

   implicit none

   private

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

end module turbogap_results
