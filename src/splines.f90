! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, splines.f90, is copyright (c) 2019-2021, Miguel A. Caro
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


module splines


  contains


!**************************************************************************
  function spline(x, y, y2, yp1, ypn, r, rcut) result(s)

!   This function computes a cubic spline on all the points defined in
!   array r(:), for those r(i) < rcut, but constructing a cubic spline
!   interpolation. The tabulated data is given in arrays x(:) and y(:),
!   which should have the same size, assuming y(i) = y[x(i)].
!   y2(:) is the (precomputed) array of second derivatives of y(:) at the x(:)
!   points. yp1 and ypn give the constant derivatives for values to
!   the left and right of the tabulated domain. This code is an
!   implementation of the cubic spline algorithm outlined in Numerical
!   Recipes, but the code has been written from scratch from the
!   equations.
!   The algorithm is as follows for interpolating s(i) at r(i):
!
!   if r(i) > rcut: don't do anything
!
!   if r(i) < x(1): s(i) = y(1) + (r(i) - x(1)) * yp1
!
!   if r(i) > x(n): s(i) = y(n) + (r(i) - x(n)) * ypn
!
!   else:
!
!     find the index j such that x(j) <= r(i) < x(j+1)
!     h = x(j+1) - x(j) [it must be > 0.d0!!!]
!     A = (x(j+1) - r(i)) / h
!     B = 1 - A
!     C = (A^3 - A) * h^2 / 6
!     D = (B^3 - B) * h^2 / 6
!     s(i) = A*y(j) + B*y(j+1) + C*y2(j) + D*y2(j+1)

    implicit none

    real*8, intent(in) :: x(:), y(:), y2(:), r(:), yp1, ypn, rcut
    real*8, dimension(1:size(r)) :: s

    real*8 :: h, h26, A, B, C, D
    integer :: j, n, m, i

    n = size(x)
    m = size(r)

    s = 0.d0

    do i = 1, m
      if( r(i) < rcut )then
        if( r(i) < x(1) )then
          s(i) = y(1) + (r(i) - x(1)) * yp1
        else if( r(i) > x(n) )then
          s(i) = y(n) + (r(i) - x(n)) * ypn
        else
          do j = 1, n-1
            if( r(i) < x(j+1) )then
              exit
            end if
          end do

          h = x(j+1) - x(j)
          h26 = h**2 / 6.d0
          if( h == 0.d0 )then
            write(*, *) "spline: h=0 -> check your tabulated values!"
            stop
          end if

          A = (x(j+1) - r(i)) / h
          B = 1.d0 - A
          C = (A**3 - A) * h26
          D = (B**3 - B) * h26
          s(i) = A*y(j) + B*y(j+1) + C*y2(j) + D*y2(j+1)
        end if
      end if
    end do

  end function

  function spline_der(x, y, y2, yp1, ypn, r, rcut) result(ds)

!   This just gives the derivative of spline()

    implicit none

    real*8, intent(in) :: x(:), y(:), y2(:), r(:), yp1, ypn, rcut
    real*8, dimension(1:size(r)) :: ds

    real*8 :: h, h6, dAdx, dBdx, dCdx, dDdx, A, B
    integer :: j, n, m, i

    n = size(x)
    m = size(r)

    ds = 0.d0

    do i = 1, m
      if( r(i) < rcut )then
        if( r(i) < x(1) )then
          ds(i) = yp1
        else if( r(i) > x(n) )then
          ds(i) = ypn
        else
          do j = 1, n-1
            if( r(i) < x(j+1) )then
              exit
            end if
          end do

          h = x(j+1) - x(j)
          h6 = h / 6.d0
          if( h == 0.d0 )then
            write(*, *) "spline: h=0 -> check your tabulated values!"
            stop
          end if

          A = (x(j+1) - r(i)) / h
          B = 1.d0 - A
          dAdx = -1.d0 / h
          dBdx = -dAdx
          dCdx = (1.d0 - 3.d0*A**2) * h6
          dDdx = (3.d0*B**2 - 1.d0) * h6

          ds(i) = dAdx*y(j) + dBdx*y(j+1) + dCdx*y2(j) + dDdx*y2(j+1)
        end if
      end if
    end do

  end function
!**************************************************************************


end module
