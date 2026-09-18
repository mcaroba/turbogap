! Does the shortcut in get_distance return what the 27-image search returns?
!
! The non-orthorhombic branch brute-forces the 27 images around the rounded
! fractional offset. It now skips them whenever the rounded image is closer
! than half the smallest perpendicular width of the cell, which makes that
! image provably the nearest one. This checks the two agree bit for bit.
!
! The reference is the same routine with the shortcut switched off, compiled
! from a copy of the source into a second binary: what a change claiming to be
! bit-exact must reproduce is itself without the change. Both binaries run
! this driver and dump every distance they return; run.sh compares the dumps.
!
! The cells are chosen for how much of the shortcut they use -- hexagonal,
! where nearly every pair takes it, down to a thin sheared cell where almost
! none do. The driver prints that fraction per cell, because a cell where it
! is zero would make the comparison vacuous.
program minimageverify

   use kinds, only: dp
   use neighbors, only: get_distance

   implicit none

   integer, parameter :: n_cases = 6
   integer, parameter :: n_pairs = 5000
   real(dp) :: cell(1:3, 1:3, 1:n_cases)
   character(len=28) :: label(1:n_cases)
   integer :: c

!  The glassy carbon MAD cell, which is what this was written for.
   cell(1:3, 1, 1) = [29.52d0, 0.d0, 0.d0]
   cell(1:3, 2, 1) = [-9.84d0, 17.0433799464d0, 0.d0]
   cell(1:3, 3, 1) = [0.d0, 0.d0, 40.2d0]
   label(1) = "hexagonal, MAD"

   cell(1:3, 1, 2) = [20.d0, 0.d0, 0.d0]
   cell(1:3, 2, 2) = [0.5d0, 20.d0, 0.d0]
   cell(1:3, 3, 2) = [0.d0, 0.5d0, 20.d0]
   label(2) = "barely sheared"

   cell(1:3, 1, 3) = [20.d0, 0.d0, 0.d0]
   cell(1:3, 2, 3) = [3.d0, 20.d0, 0.d0]
   cell(1:3, 3, 3) = [2.d0, 1.5d0, 20.d0]
   label(3) = "sheared"

   cell(1:3, 1, 4) = [18.d0, 0.d0, 0.d0]
   cell(1:3, 2, 4) = [7.d0, 16.d0, 0.d0]
   cell(1:3, 3, 4) = [5.d0, 4.d0, 20.d0]
   label(4) = "strongly sheared"

!  Perpendicular width 2 A on one axis against a 10 A edge: the shortcut is
!  almost never available, so this is the fallback under test.
   cell(1:3, 1, 5) = [10.d0, 0.d0, 0.d0]
   cell(1:3, 2, 5) = [9.5d0, 2.d0, 0.d0]
   cell(1:3, 3, 5) = [0.d0, 0.d0, 10.d0]
   label(5) = "thin and sheared"

   cell(1:3, 1, 6) = [12.d0, 0.d0, 0.d0]
   cell(1:3, 2, 6) = [11.d0, 3.d0, 0.d0]
   cell(1:3, 3, 6) = [6.d0, 5.d0, 4.d0]
   label(6) = "skewed on every axis"

   open (unit=10, file="dump.txt", status="replace", action="write")
   do c = 1, n_cases
      call one_case(cell(1:3, 1:3, c), label(c))
   end do
   close (10)

contains

!  Positions from a fixed linear congruential generator rather than
!  random_number, whose sequence is not specified and differs between
!  compilers. A test that cannot be reproduced on another machine is not one.
   subroutine fill_positions(cell_in, n, positions)

      real(dp), intent(in) :: cell_in(1:3, 1:3)
      integer, intent(in) :: n
      real(dp), intent(out) :: positions(1:3, 1:n)
      integer(kind=8) :: state
      integer :: i
      integer :: k
      real(dp) :: s(1:3)

      state = 20260918_8
      do i = 1, n
         do k = 1, 3
            state = modulo(6364136223846793005_8*state + 1442695040888963407_8, 9223372036854775807_8)
            s(k) = dble(modulo(state/65536_8, 1000000_8))/1000000.d0
         end do
         positions(1:3, i) = s(1)*cell_in(1:3, 1) + s(2)*cell_in(1:3, 2) + s(3)*cell_in(1:3, 3)
      end do

   end subroutine fill_positions

   function cross(u, v) result(w)

      real(dp), intent(in) :: u(1:3)
      real(dp), intent(in) :: v(1:3)
      real(dp) :: w(1:3)

      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)

   end function cross

!  Distance between the two cell planes normal to each lattice direction.
   function widths(cell_in) result(w)

      real(dp), intent(in) :: cell_in(1:3, 1:3)
      real(dp) :: w(1:3)
      real(dp) :: bxc(1:3)
      real(dp) :: cxa(1:3)
      real(dp) :: axb(1:3)
      real(dp) :: vol

      bxc = cross(cell_in(1:3, 2), cell_in(1:3, 3))
      cxa = cross(cell_in(1:3, 3), cell_in(1:3, 1))
      axb = cross(cell_in(1:3, 1), cell_in(1:3, 2))
      vol = dot_product(cell_in(1:3, 1), bxc)
      w(1) = dabs(vol)/dsqrt(dot_product(bxc, bxc))
      w(2) = dabs(vol)/dsqrt(dot_product(cxa, cxa))
      w(3) = dabs(vol)/dsqrt(dot_product(axb, axb))

   end function widths

   subroutine one_case(cell_in, label_in)

      real(dp), intent(in) :: cell_in(1:3, 1:3)
      character(len=*), intent(in) :: label_in
      real(dp), allocatable :: pos(:, :)
      real(dp) :: dist(1:3)
      real(dp) :: d
      real(dp) :: w(1:3)
      real(dp) :: near
      integer :: i_shift(1:3)
      integer :: p
      integer :: n_near

      allocate (pos(1:3, 1:2*n_pairs))
      call fill_positions(cell_in, 2*n_pairs, pos)
      w = widths(cell_in)
      near = 0.5d0*minval(w)
      n_near = 0
      do p = 1, n_pairs
         call get_distance(pos(1:3, 2*p - 1), pos(1:3, 2*p), cell_in(1:3, 1), cell_in(1:3, 2), &
                           cell_in(1:3, 3), [.true., .true., .true.], dist, d, i_shift)
         write (10, '(A,1X,I0,4ES26.17)') trim(label_in), p, dist(1), dist(2), dist(3), d
         if (d < near) n_near = n_near + 1
      end do
      write (*, '(A,A28,A,I5,A,I5,A,F5.1,A)') "  ", label_in, ": ", n_near, " of ", n_pairs, &
         " pairs inside half the smallest width (", 100.d0*dfloat(n_near)/dfloat(n_pairs), " %)"
      deallocate (pos)

   end subroutine one_case

end program minimageverify
