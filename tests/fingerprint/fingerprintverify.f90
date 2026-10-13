!  Drives src/fingerprint.f90. Prints reference values for run.sh to check
!  against an independent Python computation, then property checks.
program fingerprintverify

   use kinds, only: dp
   use fingerprint, only: fingerprint_of
   use, intrinsic :: iso_fortran_env, only: int32, int64

   implicit none

   real(dp), allocatable :: a(:, :)
   real(dp), allocatable :: b(:, :)
   real(dp), allocatable :: none(:, :)
   real(dp), allocatable :: empty(:, :)
   integer, allocatable :: sp(:)
   character*8, allocatable :: names(:)
   logical, allocatable :: fix(:, :)
   integer(int64) :: h0
   integer(int32), allocatable :: words(:)
   integer :: i
   integer :: j
   integer :: bit
   integer :: n_fail

   n_fail = 0
   allocate (a(1:3, 1:5))
   do j = 1, 5
      do i = 1, 3
         a(i, j) = real(i, dp)*0.1_dp - real(j, dp)*1.7_dp
      end do
   end do
   sp = [1, 2, 2, 1, -7]
   names = [character(len=8) :: "C", "O", "Si"]
   allocate (fix(1:3, 1:2))
   fix = .false.
   fix(2, 1) = .true.

!  Reference lines: kind, the input words, the fingerprint.
   words = transfer(a, [0_int32])
   write (*, '(A,2(1X,I0),*(1X,I0))') 'REF r2', shape(a), words
   write (*, '(A,1X,I0)') 'HASH', fingerprint_of(a)
   words = transfer(sp, [0_int32])
   write (*, '(A,2(1X,I0),*(1X,I0))') 'REF i1', size(sp), 1, words
   write (*, '(A,1X,I0)') 'HASH', fingerprint_of(sp)

   h0 = fingerprint_of(a)
   b = a
   if (fingerprint_of(b) /= h0) call fail('a copy changed the fingerprint')

!  Every single-bit flip of every word.
   words = transfer(a, [0_int32])
   do i = 1, size(words)
      do bit = 0, 31
         words(i) = ieor(words(i), ishft(1_int32, bit))
         b = reshape(transfer(words, a), shape(a))
         if (fingerprint_of(b) == h0) call fail('a bit flip kept the fingerprint')
         words(i) = ieor(words(i), ishft(1_int32, bit))
      end do
   end do

   b = a
   b(:, 2) = a(:, 4)
   b(:, 4) = a(:, 2)
   if (fingerprint_of(b) == h0) call fail('swapping two atoms kept the fingerprint')

   b = reshape(a, [5, 3])
   if (fingerprint_of(b) == h0) call fail('a reshape kept the fingerprint')

   allocate (empty(1:3, 1:0))
   if (fingerprint_of(none) == fingerprint_of(empty)) call fail('unallocated looks like empty')

   h0 = fingerprint_of(names)
   names(3) = "Sj"
   if (fingerprint_of(names) == h0) call fail('a changed species name kept the fingerprint')
   h0 = fingerprint_of(fix)
   fix(3, 2) = .true.
   if (fingerprint_of(fix) == h0) call fail('a changed fix_atom kept the fingerprint')

   if (n_fail == 0) then
      write (*, '(A)') 'PROPERTIES OK'
   else
      write (*, '(A,I0)') 'PROPERTY FAILURES ', n_fail
   end if

contains

   subroutine fail(what)
      character(len=*), intent(in) :: what

      write (*, '(A,A)') 'FAIL: ', what
      n_fail = n_fail + 1
   end subroutine fail

end program fingerprintverify
