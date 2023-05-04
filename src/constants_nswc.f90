!Support file for qxgs, in public domain.

module constants_nswc
! Contains the NSWC functions IPMPAR, SPMPAR, DPMPAR, EPSLN, DEPSLN,
! EXPARG & DXPARG
!-----------------------------------------------------------------------
!     WRITTEN using F90 intrinsics by
!        Alan Miller
!        CSIRO Mathematical & Information Sciences
!        CLAYTON, VICTORIA, AUSTRALIA 3169
!     Latest revision - 1 February 1997
!-----------------------------------------------------------------------

integer, parameter     :: dp = SELECTED_real_KIND(15, 60)

CONTAINS

function ipmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE integer MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN integer
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  integerS.

!     ASSUME integerS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND real FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

integer, intent(in) :: i
integer             :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0D0)
  CASE( 9)
    fn_val = MINEXPONENT(1.0D0)
  CASE(10)
    fn_val = MAXEXPONENT(1.0D0)
  CASE DEFAULT
    RETURN
end SELECT

RETURN
end function ipmpar



function spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN integer HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

integer, intent(in) :: i
real                :: fn_val

! Local variable
real                :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
end SELECT

RETURN
end function spmpar



function dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE real MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN integer HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     real ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

integer, intent(in) :: i
real (dp)           :: fn_val

! Local variable
real (dp)    :: one = 1._dp

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
end SELECT

RETURN
end function dpmpar


function epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
real                :: fn_val

! Local variable
real                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
end function epsln


function exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
integer, intent(in) :: l
real                :: fn_val

! Local variable
real                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
end IF
RETURN
end function exparg


function depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
real (dp)           :: fn_val

! Local variable
real (dp)    :: one = 1._dp

fn_val = LOG( EPSILON(one) )
RETURN
end function depsln


function dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
integer, intent(in) :: l
real (dp)           :: fn_val

! Local variable
real (dp)    :: one = 1._dp

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
end IF
RETURN
end function dxparg

end module constants_nswc
