module misc

  contains

  subroutine integrate(method, y, x, integral)

    implicit none

!   Input variables
    real*8, intent(in) :: y(:), x(:)
    character(len=*), intent(in) :: method
!   Output variables
    real*8, intent(out) :: integral

!   Internal variables
    integer :: n

    if( size(x) /= size(y) )then
      write(*,*) "Integrate subroutine: x and y should have the same dimension!"
      stop
    end if

    n = size(x)

    if( method == "trapezoidal" )then
!     Implements trapezoidal rule
    else if( method == "simpson" )then
!     Implements Simpson's rule
    end if

  end subroutine

end module misc
