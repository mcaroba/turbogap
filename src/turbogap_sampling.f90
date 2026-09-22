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
   use types, only: image

   implicit none

   private

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

end module turbogap_sampling
