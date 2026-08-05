! Numeric kind parameters for TurboGAP.
!
! The tree declared storage as real*8, complex*16 and integer*8. Those forms
! are not standard Fortran, they are not portable, and they say nothing about
! what the number is for. dp is taken from iso_fortran_env's real64, so it is
! exactly the same 64-bit representation those declarations already had: this
! is a rename, not a change of precision, and no result moves.
!
! src/soap_turbo is a git submodule with its own upstream, so it must not
! depend on this module. Those files declare the same parameters locally
! instead, which keeps the submodule self-contained.
!
! src/third_party is left alone as well: the NSWC code there already defines
! its own dp through constants_nswc, and re-declaring it here would collide.
module kinds

   use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64

   implicit none

   public

! Real kinds
   integer, parameter :: sp = real32
   integer, parameter :: dp = real64

! Integer kinds
   integer, parameter :: i32 = int32
   integer, parameter :: i64 = int64

end module kinds
