! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, fingerprint.f90, is copyright (c) 2026, Tigany Zarrouk
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

!  A 64-bit fingerprint of an array's bits and shape, so ranks can tell whether
!  their copies agree without sending them. Two 32-bit murmur-style lanes with
!  31-bit multipliers: every product stays below 2**63, so nothing overflows.
module fingerprint

   use kinds, only: dp
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64

   implicit none

   private
   public :: fingerprint_of

   integer(int64), parameter :: MASK32 = 4294967295_int64
   integer(int64), parameter :: M1 = 1540483477_int64
   integer(int64), parameter :: M2 = 668265261_int64
   integer(int64), parameter :: UNALLOCATED = -1_int64

   interface fingerprint_of
      module procedure fingerprint_r1
      module procedure fingerprint_r2
      module procedure fingerprint_i1
      module procedure fingerprint_l2
      module procedure fingerprint_c1
   end interface fingerprint_of

contains

   function fingerprint_r1(x) result(h)
      real(dp), allocatable, intent(in) :: x(:)
      integer(int64) :: h

      h = UNALLOCATED
      if (allocated(x)) h = hash_words(transfer(x, [0_int32]), [size(x), 1])
   end function fingerprint_r1

   function fingerprint_r2(x) result(h)
      real(dp), allocatable, intent(in) :: x(:, :)
      integer(int64) :: h

      h = UNALLOCATED
      if (allocated(x)) h = hash_words(transfer(x, [0_int32]), shape(x))
   end function fingerprint_r2

   function fingerprint_i1(x) result(h)
      integer, allocatable, intent(in) :: x(:)
      integer(int64) :: h

      h = UNALLOCATED
      if (allocated(x)) h = hash_words(transfer(x, [0_int32]), [size(x), 1])
   end function fingerprint_i1

   function fingerprint_l2(x) result(h)
      logical, allocatable, intent(in) :: x(:, :)
      integer(int64) :: h

      h = UNALLOCATED
      if (allocated(x)) h = hash_words(transfer(x, [0_int32]), shape(x))
   end function fingerprint_l2

!  Bytes padded to whole words, so any string length works.
   function fingerprint_c1(x) result(h)
      character(len=*), allocatable, intent(in) :: x(:)
      integer(int64) :: h
      integer(int8), allocatable :: bytes(:)
      integer :: n_bytes

      h = UNALLOCATED
      if (.not. allocated(x)) return
      n_bytes = len(x)*size(x)
      allocate (bytes(1:4*((n_bytes + 3)/4)))
      bytes = 0_int8
      if (n_bytes > 0) bytes(1:n_bytes) = transfer(x, [0_int8])
      h = hash_words(transfer(bytes, [0_int32]), [size(x), len(x)])
   end function fingerprint_c1

   function hash_words(words, dims) result(h)
      integer(int32), intent(in) :: words(:)
      integer, intent(in) :: dims(1:2)
      integer(int64) :: h
      integer(int64) :: h1
      integer(int64) :: h2
      integer :: i

      h1 = iand(ieor(2166136261_int64, int(dims(1), int64)), MASK32)
      h2 = iand(ieor(3266489917_int64, int(dims(2), int64)), MASK32)
      h1 = lane_step(h1, int(dims(2), int64), M1)
      h2 = lane_step(h2, int(dims(1), int64), M2)
      do i = 1, size(words)
         h1 = lane_step(h1, int(words(i), int64), M1)
         h2 = lane_step(h2, int(words(i), int64), M2)
      end do
      h = ior(ishft(lane_final(h1, M1), 32), lane_final(h2, M2))
   end function hash_words

   pure elemental function lane_step(h, w, m) result(hn)
      integer(int64), intent(in) :: h
      integer(int64), intent(in) :: w
      integer(int64), intent(in) :: m
      integer(int64) :: hn
      integer(int64) :: k

      k = iand(w, MASK32)
      k = iand(k*m, MASK32)
      k = ieor(k, ishft(k, -24))
      k = iand(k*m, MASK32)
      hn = ieor(iand(h*m, MASK32), k)
   end function lane_step

   pure elemental function lane_final(h, m) result(hf)
      integer(int64), intent(in) :: h
      integer(int64), intent(in) :: m
      integer(int64) :: hf

      hf = ieor(h, ishft(h, -13))
      hf = iand(hf*m, MASK32)
      hf = ieor(hf, ishft(hf, -15))
   end function lane_final

end module fingerprint
