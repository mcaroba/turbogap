# 1 "src/soap_turbo/src/soap_turbo_angular.f90"
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
module soap_turbo_angular
  use F_B_C
  use iso_c_binding

  contains

!**************************************************************************
! Returns an array with the Associated Legendre polynomials Plm(x), where
! Pll(x) is the highest order in the series, together with, all the lower
! order ones. l is > 0, m can take values from 0 to l and x is within the
! [-1, 1] domain plm_array is a 1D array with modified index lm -> k, where
! k = 1 + l*(l+1)/2 + m
  subroutine get_plm_array(plm_array, lmax, x)

    implicit none
    integer :: l, m, lmax, k
    real*8 :: plm_array(:), x
    if(lmax < 0 .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
!   We need these 6 polynomials to initialize the recursion series
!   P_00(x), k = 1
    plm_array(1) = 1.d0
    if( lmax > 0 )then
!     P_10(x), k = 2
      plm_array(2) = x
!     P_11(x), k = 3
      plm_array(3) = - dsqrt(1.d0-x**2)
      if( lmax > 1 )then
!       P_20(x), k = 4
        plm_array(4) = 1.5d0*x**2 - 0.5d0
!       P_21(x), k = 5
        plm_array(5) = -3.d0 * x * dsqrt(1.d0-x**2)
!       P_22(x), k = 6
        plm_array(6) = 3.d0 - 3.d0*x**2
      else
        return
      end if
    else
      return
    end if
    if( lmax == 2 )then
      return
    else
      do l = 3, lmax
!       First we need to obtain Pl0, Pl1, up to m=l-2 with the recursion formula on l:
!       Plm(x) = ( (2l-1)*x*Pl-1m(x) - (l-1+m)*Pl-2m(x) ) / (l-m)
        do m = 0, l-2
          k = 1 + l*(l+1)/2 + m
          plm_array(k) = (dfloat(2*l-1)*x*plm_array(k-l) - dfloat(l-1+m)*plm_array(k-2*l+1)) / dfloat(l-m)
        end do
!       Now we get Pll-1 and Pll with the recursion formulas
!       Pll-1(x) = x*(2l-1)*Pl-1l-1(x)
!       Pll(x) = -(2l-1)*sqrt(1-x^2)*Pl-1l-1(x)
        k = k + 1
        plm_array(k) = x * dfloat(2*l-1) * plm_array(k-l)
        k = k + 1
        plm_array(k) = - dfloat(2*l-1) * dsqrt(1.d0 - x**2) * plm_array(k-l-1)
      end do
    end if
  return
  end subroutine
!**************************************************************************






!**************************************************************************
! This subroutine computes arrays of associated Legendre polynomials of the
! first kind which have been suitably modified to avoid the singularity at
! the poles (theta = 0, pi). For the derivative wrt the azimuthal angle phi
! we return -m * Plm(cos(t)) / sin(t). For m = 0, the coefficients are
! zero.
! For the derivative wrt the azimuthal angle theta, we return:
! sin(t) * d(Plm(cos(t)))/d(cos(t)).
! These tricks remove all the singularities in the gradient in spherical
! coordinates except for that at r = 0.
!
! The subroutine takes as input the Plm, which have been computed up to
! lmax+1, since that is needed in one of the recursions 
  subroutine get_plm_array_der(plm_array, lmax, x, plm_array_div_sin, plm_array_der_mul_sin)

    implicit none
    integer :: l, m, lmax, k, k_lp1_mp1, k_lp1_mm1, k_l_mp1, k_l_mm1, k_temp
    real*8, intent(in) :: plm_array(:), x
    real*8, intent(out) :: plm_array_div_sin(:), plm_array_der_mul_sin(:)
    real*8 :: part1, part2

    if(lmax < 0 .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if

!   This is for
!   sin(t) * d(Plm(cos(t)))/d(cos(t)) = 1/2 * [ (l+m)*(l-m+1)*P_l^{m-1}(cost(t))
!                                               - P_l^{m+1}(cost(t)) ]
    plm_array_der_mul_sin = 0.d0
    do l = 0, lmax
      do m = 0, l
        k = 1 + l*(l+1)/2 + m
        k_l_mp1 = k + 1
        k_l_mm1 = k - 1
!       If m = 0 then we are asking for P_l^{-1}, which is not defined. We need
!       to rewrite in terms of P_l^1:
        if( m == 0 )then
!         P_0^1 = 0
          if( l == 0 )then
            part1 = 0.d0
!         P_l^{-1} = - (l-1)!/(l+1)! * P_l^1
          else
            k_temp = 1 + l*(l+1)/2 + 1
            part1 = -0.5d0 * plm_array(k_temp)
          end if
        else
          part1 = 0.5d0 * dfloat(l+m) * dfloat(l-m+1) * plm_array(k_l_mm1)
        end if
        if( m == l )then
          part2 = 0.d0
        else
          part2 = -0.5d0 * plm_array(k_l_mp1)
        end if
        plm_array_der_mul_sin(k) = part1 + part2
      end do
    end do

!   This is for
!   -m * Plm(cos(t))/sin(t) = 1/2 * [ (l-m+1)*(l-m+2)*P_{l+1}^{m-1}(cost(t))
!                                      - P_{l+1}^{m+1}(cost(t)) ]
    plm_array_div_sin = 0.d0
    do l = 0, lmax
!     Note m starts at 1 here
      do m = 1, l
        k = 1 + l*(l+1)/2 + m
        k_lp1_mp1 = 1 + (l+1)*(l+2)/2 + m + 1
        k_lp1_mm1 = 1 + (l+1)*(l+2)/2 + m - 1
!       If m = 0 then we are asking for P_{l+1}^{-1}, which is not defined. We would need
!       to rewrite in terms of P_{l+1}^1, except that for m = 0 these coefficients are
!       zero anyway.
        part1 = 0.5d0 * dfloat(l-m+1) * dfloat(l-m+2) * plm_array(k_lp1_mm1)
        part2 = 0.5d0 * plm_array(k_lp1_mp1)
        plm_array_div_sin(k) = part1 + part2
      end do
    end do
  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the complex conjugate of all the e^{imphi}.
! It also calls the ilexp function and multiplies its result times the
! complex exponential, returning the final prefactors which depend on l
! and m. It also lets you retrieve the product of the complex exponential
! times the radial derivative of the ilexp function.
!
  subroutine get_eimphi_conjg(eimphi, prefl, prefm, fact_array, lmax, phi, rj, atom_sigma, &
                              scaling, do_derivatives, prefl_rad_der, eimphi_rad_der, eimphi_azi_der, &
                              prefl_new, prefl_der_new )
    implicit none
    real*8, intent(in) :: rj, atom_sigma, scaling
    complex*16 :: eimphi(:), eimphi_rad_der(:), eimphi_azi_der(:)
    real*8 :: phi, rjbysigma, pref, cosm2, sinm2, cosm1, sinm1, cos0, sin0, cosphi2
    integer :: lmax, l, m, k
    real*8 :: prefl_der_new(0:), prefl_new(0:)
    real*8 :: prefl(0:), fact_array(:), prefl_rad_der(0:), pref_rad_der
    complex*16 :: prefm(0:)
    logical, intent(in) :: do_derivatives

    rjbysigma = rj/atom_sigma

!   This is fast
    call get_ilexp(prefl, fact_array, lmax, rjbysigma)
    !write(*,*) prefl
    prefl=prefl_new
    !write(*,*) prefl_new
    !stop
!   Complex exponential using Euler's formula and Chebyshev recursion
    cosm2 = dcos(phi)
    cosphi2 = 2.d0 * cosm2
    sinm2 = -dsin(phi)
    cosm1 = 1.d0
    sinm1 = 0.d0
    prefm(0) = (1.d0, 0.d0)
    do l = 1, lmax
      cos0 = cosphi2 * cosm1 - cosm2
      sin0 = cosphi2 * sinm1 - sinm2
      cosm2 = cosm1
      sinm2 = sinm1
      cosm1 = cos0
      sinm1 = sin0
      prefm(l) = cos0*(1.d0, 0.d0) - sin0*(0.d0, 1.d0)
    end do
    k = 1
    do l = 0, lmax
      pref = prefl(l)
      do m = 0, l
!       This is somewhat slow but not critical
        eimphi(k) = pref * prefm(m)
        k = k + 1
      end do
    end do
    if( do_derivatives )then
!     Get ilexp derivatives
      call get_ilexp_der(prefl, lmax, rj, atom_sigma, scaling, prefl_rad_der)
      prefl_rad_der=prefl_der_new
      k = 1
      do l = 0, lmax
        pref_rad_der = prefl_rad_der(l)
        do m = 0, l
!         This is somewhat slow but not critical
!         Radial derivative:
          eimphi_rad_der(k) = pref_rad_der * prefm(m)
          k = k + 1
        end do
      end do
!     Azimuthal angle derivative:
!     In the expression for the derivatives, there is a -i*m term multiplying the rest,
!     which we omit here, only retaining i, since we include -m as part of the Plm_array
      eimphi_azi_der = eimphi * dcmplx(0.d0, 1.d0)
    end if
    return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! For a description see definition of the function ilexp in functions.f90
!
  subroutine get_ilexp(ilexp_array, fact_array, lmax, x)
    implicit none

    real*8 :: ilexp_array(0:), x, x2, x4, fl, flm1, flm2, xcut = 1.d-7, fact
    integer :: l, i, lmax
    real*8 :: fact_array(:)

    if( lmax < 0 )then
      write(*,*) "Bad argument (l<0) for function ilexp_double!"
      return
    end if
    x2 = x**2
    x4 = x**4

!   Semifactorial calculation
    if( lmax > 0 )then
      fact = 1.d0
      do i = 1, lmax
        fact = fact * (2.d0*dfloat(i) + 1.d0)
        fact_array(i) = fact
      end do
    end if
  !write(*,*) fact_array
  !stop
!   Full calculation. This is numerically unstable for small x, that's why
!   we have cases below

    if( x > 0.0d0 ) then
       flm2 = dabs( (1.d0 - dexp(-2.d0 * x2)) / 2.d0 / x2 )
       flm1 = dabs( (x2 - 1.d0 + dexp(-2.d0 * x2)*(x2+1.d0)) / 2.d0 / x4 )
    else
       flm2 = 1.0d0
       flm1 = 0.0d0
    endif

    do l = 0, lmax
      if( l == 0 )then
        if( x < xcut )then
          ilexp_array(0) = 1.d0 - x2
        else
          ilexp_array(0) = flm2
        end if
      else if( l == 1 )then
        if( x2/1000.d0 < xcut )then
          ilexp_array(1) = (x2 - x4)/ fact_array(1)
        else
          ilexp_array(1) = flm1
        end if
      else
        if( x2**l / fact_array(l) * dfloat(l) < xcut )then
          fl = x2**l / fact_array(l)
        else
          fl = dabs( flm2 - (2.d0*dfloat(l) - 1.d0)/x2 * flm1 )
        end if
        flm2 = flm1
        flm1 = fl
        ilexp_array(l) = fl
      end if
    end do
  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the radial derivatives of ilexp(rj/srj; l=0,lmax)
! If rj = 0. then the derivative is zero. The atom_sigma here is the angular
! scaled atom sigma and the scaling is the angular atom sigma scaling
!
  subroutine get_ilexp_der(ilexp_array, lmax, rj, atom_sigma, scaling, ilexp_array_der)
    implicit none

    real*8, intent(in) :: ilexp_array(0:), rj, atom_sigma, scaling
    integer, intent(in) :: lmax
    real*8, intent(inout) :: ilexp_array_der(0:)
    integer :: l
    real*8 :: coeff1, coeff2

    coeff1 = 2.d0*rj/atom_sigma**2
    coeff2 = 1.d0 - scaling*rj/atom_sigma

    ilexp_array_der = 0.d0

!   Set derivative to zero if rj is (approximately) zero
    if( rj < 1.d-5 )then
      return
    end if

    ilexp_array_der(0) = coeff1 * (ilexp_array(1) - ilexp_array(0))
    do l = 1, lmax
      ilexp_array_der(l) = (-coeff1 - dfloat(2*l+2)/rj) * ilexp_array(l) + coeff1*ilexp_array(l-1)
    end do

    ilexp_array_der = ilexp_array_der * coeff2

  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine gets the angular expansion coefficients.
!
! subroutine get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, atom_sigma_in, &
!                                               atom_sigma_scaling, rcut, lmax, &
  subroutine get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, &
                                                rcut, lmax, &
                                                eimphi, preflm, plm_array, prefl, prefm, fact_array, &
                                                mask, n_species, eimphi_rad_der, do_derivatives, &
                                                prefl_rad_der, exp_coeff, exp_coeff_rad_der, &
                                                exp_coeff_azi_der, exp_coeff_pol_der, &
                                                exp_coeff_d, exp_coeff_rad_der_d, &
                                                exp_coeff_azi_der_d, exp_coeff_pol_der_d, n_atom_pairs, &
                                                thetas_d, preflm_d, rjs_d, phis_d, mask_d, &
                                                atom_sigma_in_d, atom_sigma_scaling_d,gpu_stream)

    implicit none

    integer(c_int), intent(in) :: n_species, n_atom_pairs
    integer (c_int):: lmax, kmax, n_neigh(:), n_sites, i, j, k, kmax_der, i_sp, k_int,l,m,lmpo
    complex*16, intent(out), target :: exp_coeff(:,:), exp_coeff_rad_der(:,:), exp_coeff_azi_der(:,:), exp_coeff_pol_der(:,:)
!   real*8 :: thetas(:), phis(:), atom_sigma_in(:), atom_sigma, atom_sigma_scaling(:), rjs(:), x, theta, phi, rj
    real*8 :: thetas(:), phis(:),  atom_sigma, rjs(:), x, theta, phi, rj
    real(c_double), intent(in) :: rcut
    real*8 :: amplitude
!   I should probably allocate and save these variables internally to minimize the number of variables that
!   need to be passed to this subroutine
    complex*16 :: eimphi(:), eimphi_rad_der(:)
    real*8 :: preflm(:), plm_array(:)
    real*8 :: prefl(0:), fact_array(:), prefl_rad_der(0:)
    complex*16 :: prefm(0:)
    logical, intent(in) :: mask(:,:)
    logical, intent(in) ::  do_derivatives
    real*8, allocatable :: plm_array_der(:), plm_array_div_sin(:), plm_array_der_mul_sin(:)
    real(c_double), allocatable,target :: prefl_array_global(:,:), prefl_array_global_der(:,:)
    real(c_double), allocatable,target :: plm_array_global(:,:), plm_array_der_global(:,:)
    complex*16, allocatable :: eimphi_azi_der(:)
    complex*16, allocatable, target  :: eimphi_global(:,:), eimphi_rad_der_global(:,:)
    complex*16, allocatable, target  :: tr_exp_coeff(:,:), eimphi_azi_der_global(:,:)
    !complex*16, allocatable, target :: tr_exp_coeff_rad_der(:,:), &
    !                tr_exp_coeff_azi_der(:,:), tr_exp_coeff_pol_der (:,:)
    type(c_ptr) :: preflm_d, prefl_array_global_d
    type(c_ptr) :: thetas_d, plm_array_global_d, plm_array_der_global_d,eimphi_global_d
    type(c_ptr) :: exp_coeff_d, exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d
    type(c_ptr) :: rjs_d, phis_d, mask_d
    type(c_ptr) ::  atom_sigma_in_d, atom_sigma_scaling_d
    logical(c_bool) :: c_do_derivatives
    type(c_ptr) :: eimphi_rad_der_global_d, eimphi_azi_der_global_d
    type(c_ptr) :: prefl_array_global_der_d
    real(c_double), allocatable, target :: plm_array_div_sin_global(:,:), plm_array_der_mul_sin_global(:,:)
    type(c_ptr) :: plm_array_div_sin_d, plm_array_der_mul_sin_d
    integer(c_size_t) :: st_prefl_array_global, st_plm_array_global,  st_eimphi_global, st_plm_array_der_global
    type(c_ptr), intent(inout) :: gpu_stream
    
    c_do_derivatives=logical( .false., kind=c_bool ) 
    if(do_derivatives) then 
    c_do_derivatives=logical( .true., kind=c_bool )
    endif

    kmax = 1 + lmax*(lmax+1)/2 + lmax

    exp_coeff = 0.d0
    if( do_derivatives )then
      kmax_der = 1 + (lmax+1)*(lmax+2)/2 + lmax+1
    end if
    
    st_prefl_array_global= (lmax+1)*n_atom_pairs*sizeof(rcut) !c_double !sizeof(rcut)
    st_plm_array_global= kmax*n_atom_pairs*sizeof(rcut) !c_double !sizeof(rcut)
    st_eimphi_global = kmax*n_atom_pairs*sizeof(exp_coeff(1,1)) !c_double_complex !sizeof(exp_coeff(1,1))
    call gpu_malloc_all(prefl_array_global_d, st_prefl_array_global,gpu_stream) ! call gpu_malloc_double(prefl_array_global_d, (lmax+1)*n_atom_pairs)
    call gpu_malloc_all(plm_array_global_d, st_plm_array_global,gpu_stream) ! call gpu_malloc_double(plm_array_global_d, kmax*n_atom_pairs)
    call gpu_malloc_all(eimphi_global_d, st_eimphi_global,gpu_stream) !call gpu_malloc_double_complex(eimphi_global_d, kmax*n_atom_pairs)

    if(do_derivatives) then
    call gpu_malloc_all(prefl_array_global_der_d, st_prefl_array_global,gpu_stream) !call gpu_malloc_double(prefl_array_global_der_d, (lmax+1)*n_atom_pairs)

    call gpu_malloc_all(eimphi_rad_der_global_d, st_eimphi_global ,gpu_stream) ! call gpu_malloc_double_complex(eimphi_rad_der_global_d, kmax*n_atom_pairs)
    call gpu_malloc_all(eimphi_azi_der_global_d, st_eimphi_global,gpu_stream) ! call gpu_malloc_double_complex(eimphi_azi_der_global_d, kmax*n_atom_pairs)
    
    call gpu_malloc_all(plm_array_div_sin_d, st_plm_array_global,gpu_stream) !call gpu_malloc_double(plm_array_div_sin_d, kmax*n_atom_pairs)
    call gpu_malloc_all(plm_array_der_mul_sin_d, st_plm_array_global,gpu_stream) !call gpu_malloc_double(plm_array_der_mul_sin_d, kmax*n_atom_pairs)

    endif

    call  gpu_get_plm_array_global(plm_array_global_d, n_atom_pairs, kmax, &
                              lmax, thetas_d, gpu_stream) 
    
    
    if(do_derivatives) then
    lmpo=lmax+1 
    st_plm_array_der_global= kmax_der*n_atom_pairs*sizeof(rcut)
    call gpu_malloc_all(plm_array_der_global_d, st_plm_array_der_global, gpu_stream) !call gpu_malloc_double(plm_array_der_global_d, kmax_der*n_atom_pairs)
    call gpu_get_plm_array_global(plm_array_der_global_d, n_atom_pairs, kmax_der, &
                              lmpo, thetas_d, gpu_stream )
    
    endif

    call gpu_get_exp_coeff_array(eimphi_global_d, &
                                  rjs_d, phis_d, thetas_d, &
                                  mask_d, atom_sigma_in_d, atom_sigma_scaling_d, & 
                                  rcut, n_atom_pairs, n_species, lmax, kmax, &
                                  prefl_array_global_d,plm_array_global_d, &
                                  plm_array_der_global_d, &
                                  prefl_array_global_der_d, &
                                  preflm_d,exp_coeff_d, &
                                  c_do_derivatives, &
                                  eimphi_rad_der_global_d, eimphi_azi_der_global_d, &
                                  plm_array_div_sin_d, plm_array_der_mul_sin_d, &
                                  exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d, &
                                  gpu_stream) 



    if( do_derivatives )then
      call gpu_free_async(plm_array_der_mul_sin_d,gpu_stream)
      call gpu_free_async(plm_array_div_sin_d,gpu_stream)
      call gpu_free_async(plm_array_der_global_d,gpu_stream)
      call gpu_free_async(prefl_array_global_der_d,gpu_stream)
      call gpu_free_async(eimphi_azi_der_global_d,gpu_stream)
      call gpu_free_async(eimphi_rad_der_global_d,gpu_stream)
    end if
    call gpu_free_async(plm_array_global_d,gpu_stream)
    call gpu_free_async(eimphi_global_d,gpu_stream)
    call gpu_free_async(prefl_array_global_d,gpu_stream)
  return
  end subroutine get_angular_expansion_coefficients
!**************************************************************************




end module soap_turbo_angular
