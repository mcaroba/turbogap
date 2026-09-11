module neighbors_skin

!  The Verlet skin: deciding when a neighbours list built with a padded cutoff
!  has stopped being valid.
!
!  The list is built out to rcut + buffer, so it holds every pair inside rcut
!  now, plus every pair that could come inside rcut later. Two atoms approach
!  each other by at most the sum of their individual displacements, so the list
!  stays complete for exactly as long as
!
!     d_1 + d_2 < buffer
!
!  with d_1 and d_2 the two largest displacements since the list was built.
!  Larger displacements than that can hide a pair that was outside the padded
!  cutoff at build time and is inside rcut now, and no test on the list itself
!  can find it -- the pair is simply absent.
!
!  Two things this gets right that the inline test it replaces did not.
!
!  The comparison is between two lengths. The inline test compared a SQUARED
!  displacement against buffer/2, which fires at |d| = sqrt(buffer/2). For any
!  buffer below 2 A that is later than the bound above, so the list went stale
!  before it was rebuilt and pairs were silently dropped. That is why a nonzero
!  neighbors_buffer could not be used, and why the shipped default is 0.
!
!  The accumulator is a displacement, not a coordinate difference. Positions are
!  wrapped into the primitive cell every MD step, so an atom crossing a boundary
!  registers a full box length in positions - positions_prev. Under the old test
!  that forced a rebuild the list did not need; under the bound above it would
!  do the same. The minimum image convention removes it. The displacement of one
!  atom in one step is far below half a cell, so rounding the fractional offset
!  to the nearest integer is exact here -- no image search is needed.
!
!  buffer <= 0 means "rebuild every step" and is the shipped default, so a deck
!  that does not set the keyword behaves exactly as before.

   use kinds, only: dp
   use error, only: turbogap_abort

   implicit none
   private
   public :: skin_accumulate
   public :: skin_needs_rebuild
   public :: skin_two_largest

contains

!  Add this step's displacement to the running total since the last build,
!  under the minimum image convention of the primitive cell. cell holds the
!  lattice vectors as its columns.
   subroutine skin_accumulate(positions, positions_prev, cell, positions_diff)

      implicit none

      real(dp), intent(in) :: positions(:, :)
      real(dp), intent(in) :: positions_prev(:, :)
      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(inout) :: positions_diff(:, :)
      real(dp) :: cell_inv(1:3, 1:3)
      real(dp) :: dr(1:3)
      real(dp) :: s(1:3)
      integer :: i
      integer :: n

      n = size(positions_diff, 2)
      if (size(positions, 2) < n .or. size(positions_prev, 2) < n) then
         write (*, *) "ERROR: skin_accumulate was given fewer positions than accumulated displacements <-- ERROR"
         call turbogap_abort()
      end if
      call invert_cell(cell, cell_inv)
      do i = 1, n
         dr(1:3) = positions(1:3, i) - positions_prev(1:3, i)
         s(1:3) = matmul(cell_inv, dr)
         s(1:3) = s(1:3) - dnint(s(1:3))
         positions_diff(1:3, i) = positions_diff(1:3, i) + matmul(cell, s)
      end do

   end subroutine skin_accumulate

!  Is the list built when positions_diff was last zeroed still complete?
   pure function skin_needs_rebuild(positions_diff, buffer) result(rebuild)

      implicit none

      real(dp), intent(in) :: positions_diff(:, :)
      real(dp), intent(in) :: buffer
      logical :: rebuild
      real(dp) :: d1
      real(dp) :: d2

      if (buffer <= 0.d0) then
         rebuild = .true.
         return
      end if
      call skin_two_largest(positions_diff, d1, d2)
      rebuild = (d1 + d2 >= buffer)

   end function skin_needs_rebuild

!  The two largest displacement magnitudes. Two and not one: the bound is on
!  how fast a PAIR closes, and the pair is made of two atoms that may both be
!  moving. Selected in the square to keep the square roots out of the loop.
   pure subroutine skin_two_largest(positions_diff, d1, d2)

      implicit none

      real(dp), intent(in) :: positions_diff(:, :)
      real(dp), intent(out) :: d1
      real(dp), intent(out) :: d2
      real(dp) :: dsq
      real(dp) :: dsq1
      real(dp) :: dsq2
      integer :: i

      dsq1 = 0.d0
      dsq2 = 0.d0
      do i = 1, size(positions_diff, 2)
         dsq = positions_diff(1, i)**2 + positions_diff(2, i)**2 + positions_diff(3, i)**2
         if (dsq > dsq1) then
            dsq2 = dsq1
            dsq1 = dsq
         else if (dsq > dsq2) then
            dsq2 = dsq
         end if
      end do
      d1 = dsqrt(dsq1)
      d2 = dsqrt(dsq2)

   end subroutine skin_two_largest

!  Analytic inverse of the 3x3 cell. Kept here rather than taken from misc so
!  that the module depends on kinds and nothing else, which is what lets the
!  verifier compile it on its own.
   subroutine invert_cell(cell, cell_inv)

      implicit none

      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(out) :: cell_inv(1:3, 1:3)
      real(dp) :: det

      cell_inv(1, 1) = cell(2, 2)*cell(3, 3) - cell(2, 3)*cell(3, 2)
      cell_inv(1, 2) = cell(1, 3)*cell(3, 2) - cell(1, 2)*cell(3, 3)
      cell_inv(1, 3) = cell(1, 2)*cell(2, 3) - cell(1, 3)*cell(2, 2)
      cell_inv(2, 1) = cell(2, 3)*cell(3, 1) - cell(2, 1)*cell(3, 3)
      cell_inv(2, 2) = cell(1, 1)*cell(3, 3) - cell(1, 3)*cell(3, 1)
      cell_inv(2, 3) = cell(1, 3)*cell(2, 1) - cell(1, 1)*cell(2, 3)
      cell_inv(3, 1) = cell(2, 1)*cell(3, 2) - cell(2, 2)*cell(3, 1)
      cell_inv(3, 2) = cell(1, 2)*cell(3, 1) - cell(1, 1)*cell(3, 2)
      cell_inv(3, 3) = cell(1, 1)*cell(2, 2) - cell(1, 2)*cell(2, 1)
      det = cell(1, 1)*cell_inv(1, 1) + cell(1, 2)*cell_inv(2, 1) + cell(1, 3)*cell_inv(3, 1)
      cell_inv = cell_inv/det

   end subroutine invert_cell

end module neighbors_skin
