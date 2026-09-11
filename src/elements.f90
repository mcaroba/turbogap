! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, elements.f90, is copyright (c) 2019-2026, Miguel A. Caro
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

! Per-element radii, and the symbol -> Z lookup that reaches them.
!
! Two tables, because the two questions are different and the wrong answer to
! either is silent:
!
!   COVALENT radii decide whether two atoms are BONDED. r_i + r_j is a bond
!   length, so H-C comes out at 0.31 + 0.76 = 1.07 A against a real 1.09 A.
!   This is what topology.f90 uses.
!
!   VAN DER WAALS radii decide whether two atoms OVERLAP -- whether a trial
!   insertion has landed on top of something. They are roughly twice the
!   covalent ones (H 1.20, C 1.70, so 2.90 A), and using them for bonding makes
!   every atom in a condensed phase bonded to every neighbour, which collapses
!   the whole structure into one connected component and matches nothing.
!
! Covalent radii: Cordero et al., Dalton Trans. 2008, 2832 -- the low-spin
! values where that paper gives two. Van der Waals radii: Bondi, J. Phys. Chem.
! 68, 441 (1964), extended by Alvarez, Dalton Trans. 2013, 8617 for the
! elements Bondi does not cover.
!
! Untabulated elements return a negative radius rather than zero or a guess. A
! zero would silently make an atom bond to nothing; a guess would silently make
! it bond to the wrong things. Every caller has to say what it does about that.
module elements

   use kinds

   implicit none

   private
   public :: symbol_to_z, covalent_radius, vdw_radius, element_is_known
   public :: Z_MAX

!  Through curium. Past that the radii are extrapolations nobody has measured,
!  and a potential that needs them does not exist.
   integer, parameter :: Z_MAX = 96

   character(len=2), parameter :: SYMBOLS(1:Z_MAX) = [character(len=2) :: &
                                                      "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", &
                                                      "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", &
                                                      "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
                                                      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", &
                                                      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
                                                      "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
                                                      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
                                                      "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
                                                      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
                                                      "Pa", "U", "Np", "Pu", "Am", "Cm"]

   real(dp), parameter :: COVALENT(1:Z_MAX) = [ &
                          0.31d0, 0.28d0, 1.28d0, 0.96d0, 0.84d0, 0.76d0, 0.71d0, 0.66d0, 0.57d0, 0.58d0, &
                          1.66d0, 1.41d0, 1.21d0, 1.11d0, 1.07d0, 1.05d0, 1.02d0, 1.06d0, 2.03d0, 1.76d0, &
                          1.70d0, 1.60d0, 1.53d0, 1.39d0, 1.39d0, 1.32d0, 1.26d0, 1.24d0, 1.32d0, 1.22d0, &
                          1.22d0, 1.20d0, 1.19d0, 1.20d0, 1.20d0, 1.16d0, 2.20d0, 1.95d0, 1.90d0, 1.75d0, &
                          1.64d0, 1.54d0, 1.47d0, 1.46d0, 1.42d0, 1.39d0, 1.45d0, 1.44d0, 1.42d0, 1.39d0, &
                          1.39d0, 1.38d0, 1.39d0, 1.40d0, 2.44d0, 2.15d0, 2.07d0, 2.04d0, 2.03d0, 2.01d0, &
                          1.99d0, 1.98d0, 1.98d0, 1.96d0, 1.94d0, 1.92d0, 1.92d0, 1.89d0, 1.90d0, 1.87d0, &
                          1.87d0, 1.75d0, 1.70d0, 1.62d0, 1.51d0, 1.44d0, 1.41d0, 1.36d0, 1.36d0, 1.32d0, &
                          1.45d0, 1.46d0, 1.48d0, 1.40d0, 1.50d0, 1.50d0, 2.60d0, 2.21d0, 2.15d0, 2.06d0, &
                          2.00d0, 1.96d0, 1.90d0, 1.87d0, 1.80d0, 1.69d0]

   real(dp), parameter :: VDW(1:Z_MAX) = [ &
                          1.20d0, 1.43d0, 2.12d0, 1.98d0, 1.91d0, 1.77d0, 1.66d0, 1.50d0, 1.46d0, 1.58d0, &
                          2.50d0, 2.51d0, 2.25d0, 2.19d0, 1.90d0, 1.89d0, 1.82d0, 1.83d0, 2.73d0, 2.62d0, &
                          2.58d0, 2.46d0, 2.42d0, 2.45d0, 2.45d0, 2.44d0, 2.40d0, 2.40d0, 2.38d0, 2.39d0, &
                          2.32d0, 2.29d0, 1.88d0, 1.82d0, 1.86d0, 2.25d0, 3.21d0, 2.84d0, 2.75d0, 2.52d0, &
                          2.56d0, 2.45d0, 2.44d0, 2.46d0, 2.44d0, 2.15d0, 2.53d0, 2.49d0, 2.43d0, 2.42d0, &
                          2.47d0, 1.99d0, 2.04d0, 2.06d0, 3.48d0, 3.03d0, 2.98d0, 2.88d0, 2.92d0, 2.95d0, &
                          2.90d0, 2.87d0, 2.83d0, 2.79d0, 2.87d0, 2.81d0, 2.83d0, 2.79d0, 2.80d0, 2.74d0, &
                          2.63d0, 2.53d0, 2.57d0, 2.49d0, 2.48d0, 2.41d0, 2.29d0, 2.32d0, 2.45d0, 2.47d0, &
                          2.60d0, 2.54d0, 2.50d0, 2.50d0, 2.50d0, 2.50d0, 2.50d0, 2.80d0, 2.93d0, 2.88d0, &
                          2.71d0, 2.82d0, 2.81d0, 2.83d0, 3.05d0, 3.38d0]

contains

!  Atomic number of a chemical symbol, or 0 if it is not one.
!
!  Case is normalised rather than required: an xyz file may say "SI" or "si",
!  and refusing those would be a parse error about nothing.
   function symbol_to_z(symbol) result(z)

      implicit none

      character(len=*), intent(in) :: symbol
      integer :: z
      character(len=2) :: name
      integer :: i

      z = 0
      name = adjustl(trim(symbol))
      if (len_trim(name) == 0) return

      i = iachar(name(1:1))
      if (i >= iachar("a") .and. i <= iachar("z")) name(1:1) = achar(i - 32)
      i = iachar(name(2:2))
      if (i >= iachar("A") .and. i <= iachar("Z")) name(2:2) = achar(i + 32)

      do i = 1, Z_MAX
         if (name == SYMBOLS(i)) then
            z = i
            return
         end if
      end do

   end function symbol_to_z

   function element_is_known(z) result(known)

      implicit none

      integer, intent(in) :: z
      logical :: known

      known = (z >= 1 .and. z <= Z_MAX)

   end function element_is_known

!  Covalent radius in Angstrom, or -1 for an element with no tabulated value.
   function covalent_radius(z) result(r)

      implicit none

      integer, intent(in) :: z
      real(dp) :: r

      r = -1.d0
      if (element_is_known(z)) r = COVALENT(z)

   end function covalent_radius

!  Van der Waals radius in Angstrom, or -1 for an element with no tabulated value.
   function vdw_radius(z) result(r)

      implicit none

      integer, intent(in) :: z
      real(dp) :: r

      r = -1.d0
      if (element_is_known(z)) r = VDW(z)

   end function vdw_radius

end module elements
