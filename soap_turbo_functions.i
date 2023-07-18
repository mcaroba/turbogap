# 1 "src/soap_turbo/src/soap_turbo_functions.f90"
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2021, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
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

# 26
module soap_turbo_functions

  contains

!**************************************************************************
  function plm_double(l, m, x)
!   Associated Legendre polynomial Plm(x).
!   l is > 0, m can take values from 0 to l and x is within the [-1, 1] domain
    implicit none
    integer :: l, m, i, l2
    real*8 :: plm_double, x, dfact_arg, pl2m, pmm, pmp1m, sq
    if(l < 0 .or. m < 0 .or. m > l .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
    pmm = 1.d0
    if(m > 0) then
      sq = dsqrt( (1.d0 - x**2) )
      dfact_arg = -1.d0
      do i = 1, m
        dfact_arg = dfact_arg + 2.d0
        pmm =  -pmm * dfact_arg * sq
      end do
    end if
    if(l == m)then
      plm_double = pmm
    else
      pmp1m = x * (2.d0*dfloat(m) + 1.d0) * pmm
      if(l == m+1) then
        plm_double = pmp1m
      else
        do l2 = m+2, l
          pl2m = ( x*dfloat(2*l2-1)*pmp1m - dfloat(l2+m-1)*pmm ) / dfloat(l2-m)
          pmm = pmp1m
          pmp1m = pl2m
        end do
        plm_double = pl2m
      end if
    end if
  end function
!**************************************************************************







!**************************************************************************
  function plm_single(l, m, x)
!   Associated Legendre polynomial Plm(x).
!   l is > 0, m can take values from 0 to l and x is within the [-1, 1] domain
    implicit none
    integer :: l, m, i, l2
    real*4 :: plm_single, x, dfact_arg, pl2m, pmm, pmp1m, sq
    if(l < 0 .or. m < 0 .or. m > l .or. abs(x) > 1.)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
    pmm = 1.
    if(m > 0) then
      sq = sqrt( (1. - x**2) )
      dfact_arg = -1.
      do i = 1, m
        dfact_arg = dfact_arg + 2.
        pmm =  -pmm * dfact_arg * sq
      end do
    end if
    if(l == m)then
      plm_single = pmm
    else
      pmp1m = x * (2.*float(m) + 1.) * pmm
      if(l == m+1) then
        plm_single = pmp1m
      else
        do l2 = m+2, l
          pl2m = ( x*float(2*l2-1)*pmp1m - float(l2+m-1)*pmm ) / float(l2-m)
          pmm = pmp1m
          pmp1m = pl2m
        end do
        plm_single = pl2m
      end if
    end if
  end function
!**************************************************************************








!**************************************************************************
  function ylm_double(l, m, theta, phi)
    implicit none
    real*8 :: plm, theta, phi, pref, pi, ylm_r, ylm_i, fact1, fact2
    complex*16 :: ylm_double
    integer :: l, m, i, sgn
    pi = dacos(-1.d0)
!   Factorials for |m|
    fact1 = 1.d0
    do i = 1, l-abs(m)
      fact1 = fact1 * dfloat(i)
    end do
    fact2 = 1.d0
    do i = 1, l+abs(m)
      fact2 = fact2 * dfloat(i)
    end do
!   Sign
    sgn = 1
    do i = 1, m
      sgn = -1 * sgn
    end do
!   Prefactor
    pref = dsqrt( dfloat(2*l+1)/(4.d0*pi) * fact1 / fact2 )
!   This is pl|m| 
    plm = plm_double(l, abs(m), dcos(theta))
!   Real part
    ylm_r = pref * plm * dcos( dfloat(abs(m))*phi )
!   Imaginary part
    ylm_i = pref * plm * dsin( dfloat(abs(m))*phi )
!   For m >= 0
    if( m >= 0 )then
      ylm_double = ylm_r * (1.d0, 0.d0) + ylm_i * (0.d0, 1.d0)
!   For m < 0
    else
      ylm_double = dfloat(sgn) * ( ylm_r * (1.d0, 0.d0) - ylm_i * (0.d0, 1.d0) )
    end if
  end function
!**************************************************************************







!**************************************************************************
  subroutine get_eimphi(eimphi, prefl, prefm, lmax, phi, rj2s2)
    implicit none
    complex*16 :: eimphi(:)
    real*8 :: phi, mf, rj2s2, pref
    integer :: lmax, l, m, k
!    logical, save :: init = .true.
!    real*8, allocatable, save :: prefl(:)
!    complex*16, allocatable, save :: prefm(:)
    real*8 :: prefl(0:)
    complex*16 :: prefm(0:)

!    if( init )then
!      allocate( prefl(0:lmax) )
!      allocate( prefm(0:lmax) )
!      init = .false.
!    end if
    do l = 0, lmax
      prefl(l) = ilexp_double(l, rj2s2)
      mf = dfloat(l)
      prefm(l) = dcos(mf*phi)*(1.d0, 0.d0) + dsin(mf*phi)*(0.d0, 1.d0)
    end do
    k = 1
    do l = 0, lmax
      pref = prefl(l)
      do m = 0, l
        eimphi(k) = pref * prefm(m)
        k = k + 1
      end do
    end do
  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine get_preflm(preflm, lmax)
    implicit none
    real*8 :: preflm(:)
    real*8 :: ppi, fact1, fact2
    integer :: lmax, l, m, k, kmax, i
    kmax = 1 + lmax*(lmax+1)/2 + lmax
    k = 1
    do l = 0, lmax
      ppi = dsqrt(dfloat(2*l+1)/4.d0/dacos(-1.d0))
!     This calculates the factorial of l
      fact2 = 1.d0
      do i = 1, l
        fact2 = fact2 * dfloat(i)
      end do
      fact1 = fact2
      do m = 0, l
        if( m > 0 )then
!         This is the factorial of (l-m)
          fact1 = fact1 / dfloat(l+1-m)
!         This is the factorial of (l+m)
          fact2 = fact2 * dfloat(l+m)
        end if
        preflm(k) = ppi * dsqrt(fact1/fact2)
        k = k + 1
      end do
    end do
  end subroutine
!**************************************************************************






!**************************************************************************
! This function is the product of the modified spherical Bessel function o
! the first kind i_l(x^2) times a decaying exponential exp(-x^2):
!
! ilexp(l, x) = i_l(x^2) * exp(-x^2)
!
! This function is based on the recursion relation:
!
! ilexp(l, x) = ilexp(l-2, x) - (2l-1) / x^2 * ilexp(l-1, x)
!
! We need the two first functions:
!
! ilexp(0, x) = (1 - exp(-2 x^2)) / 2 / x^2
! ilexp(1, x) = (x^2 - 1 + exp(-2 x^2)*(x^2+1)) / 2 / x^4
!
! For x close to zero we use the Taylor expansion of the
! exponential to avoid the singularity:
!
! ilexp(0, x) = 1 - x^2
! ilexp(l, x) = x^2l / (2l+1)!!
!
! E.g.:
! ilexp(1,x) = x^2 / 3
! ilexp(2,x) = x^4 / 3 / 5
! ilexp(3,x) = x^6 / 3 / 5 / 7
! etc.
!
! With these settings the numerical noise stays below 1.d-7, at least
! up till l = 19
!
  function ilexp_double(l, x)
    implicit none

    real*8 :: ilexp_double, x, x2, x4, fl, flm1, flm2, xcut = 1.d-7, fact
    integer :: l, i

    if( l < 0 )then
      write(*,*) "Bad argument (l<0) for function ilexp_double!"
      ilexp_double = 0.d0
      ilexp_double = 1.d0/ilexp_double
      return
    end if

    x2 = x**2
    x4 = x**4

!   Approximation between zero and xcut
    if( l == 0 .and. x < xcut )then
      ilexp_double = 1.d0 - x2
      return
    end if
!   Semifactorial calculation
    if( l > 0 )then
      fact = 1.d0
      do i = 1, l
        fact = fact * (2.d0*dfloat(i) + 1.d0)
      end do
    end if
    if( l == 1 .and. x2/1000.d0 < xcut )then
      ilexp_double = (x2 - x4)/ fact
      return
    else if( l > 1 .and. x2**l / fact * dfloat(l) < xcut )then
      ilexp_double = x2**l / fact
      return
    end if

!   Full calculation. This is numerically unstable for small x
    flm2 = dabs( (1.d0 - dexp(-2.d0 * x2)) / 2.d0 / x2 )
    flm1 = dabs( (x2 - 1.d0 + dexp(-2.d0 * x2)*(x2+1.d0)) / 2.d0 / x4 )
    if( l == 0 )then
      ilexp_double = flm2
    else if( l == 1 )then
      ilexp_double = flm1
    else
      do i = 2, l
        fl = dabs( flm2 - (2.d0*dfloat(i) - 1.d0)/x2 * flm1 )
        flm2 = flm1
        flm1 = fl
      end do
      ilexp_double = fl
    end if
  end function
!**************************************************************************






!**************************************************************************
  function cross_product(u, v)

    implicit none

    real*8, intent(in) :: u(1:3), v(1:3)
    real*8 :: cross_product(1:3)

    cross_product(1) = u(2)*v(3) - u(3)*v(2)
    cross_product(2) = u(3)*v(1) - u(1)*v(3)
    cross_product(3) = u(1)*v(2) - u(2)*v(1)

  end function cross_product
!**************************************************************************



end module soap_turbo_functions
