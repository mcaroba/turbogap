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

   implicit none

   private

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

end module turbogap_structure
