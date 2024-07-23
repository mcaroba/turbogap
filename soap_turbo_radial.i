# 1 "src/soap_turbo/src/soap_turbo_radial.f90"
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


# 27
module soap_turbo_radial

  use soap_turbo_functions
  use iso_c_binding

  contains

!**************************************************************************
!
! This function returns the normalization coefficients for the polynomial
! basis functions used to construct the orthonormal basis.
!
  function N_a(rcut, a)
    implicit none

    integer :: a, b
    real*8 :: rcut, N_a

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    N_a = dsqrt( rcut / dfloat(b) )
!   *********************************************

  return
  end function
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and the saparable Gaussian representation for the
! atomic sites.
!
  subroutine get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                     rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                     amplitude_scaling, nf, W, scaling_mode, mask, &
                                                     radial_enhancement, do_derivatives, do_central, &
                                                     central_weight, exp_coeff, exp_coeff_der)
!   Expansion coefficients using the polynomial basis with smooth filter
!
!   Apparently, with this basis the W matrix becomes complex for 12 and higher
!   basis functions, so alpha_max should not be higher than 11.
!   ADD ERROR MESSAGE ABOUT THIS IN THE CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in,  rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
    real(c_double),intent(in) :: rcut_hard_in
!
    integer :: n, i, j, k, alpha_max_der
    real*8 :: I_n, I_np1, I_np2, pi, sq2, rj, N_n, N_np1, N_np2, C1, C2, dr, s2, sf2
    real*8 :: W(:,:)
    real*8 :: atom_sigma_f, rj_f
!   Results will be stored in exp_coeff, which is an array of dimension
!   (alpha_max, n_atom_pairs)
    real(c_double),intent(inout), target :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_temp1(:), exp_coeff_temp2(:), exp_coeff_der_temp(:)
    logical, save :: print_basis = .false.
    real*8 :: denom, der_sjf_rj, der_rjf_rj, amplitude_der, pref_f, der_pref_f

!   If the user requests derivatives, we need to get the expansion coefficients up to
!   alpha_max + 2
    if( do_derivatives )then
      alpha_max_der = alpha_max + 2
    else
      alpha_max_der = alpha_max
    end if
    
    allocate( exp_coeff_temp1(1:alpha_max_der) )
    allocate( exp_coeff_temp2(1:alpha_max_der) )
    allocate( exp_coeff_der_temp(1:alpha_max) )


!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

    pi = dacos(-1.d0)
    sq2 = dsqrt(2.d0)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if
!
!   **************** New basis ******************
!   Redefine all the distances by dividing them by rcut_hard
!   We do this to avoid numerical instability when the value
!   of alpha is too high
!
!    rcut_soft = rcut_soft_in
!    rcut_hard = rcut_hard_in
!    atom_sigma = atom_sigma_in
!    dr = rcut_hard - rcut_soft
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
    dr = 1.d0 - rcut_soft_in/rcut_hard_in
!   *********************************************
!


    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
!       Check if we need to do the central atom
        if( j == 1 .and. .not. do_central )then
          cycle
        end if
        if( rjs_in(k) < rcut_hard_in .and. mask(k) )then
          pref_f = 0.d0
          exp_coeff_temp1 = 0.d0
          exp_coeff_temp2 = 0.d0
          exp_coeff_der_temp = 0.d0
!   **************** New basis ******************
!          rj = rjs_in(k)
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
!         We leave this here because in the future atom_sigma will be rj dependent
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitude = 1.d0 / atom_sigma_scaled
              amplitude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              amplitude = 0.d0
              amplitude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                amplitude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              else
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                amplitude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( j == 1 )then
            amplitude = central_weight * amplitude
            amplitude_der = central_weight * amplitude_der
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitude_der = amplitude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            amplitude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            amplitude = amplitude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            amplitude_der = amplitude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            amplitude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            amplitude = amplitude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if       
!         We have the recursion series starting at n = 0, which means alpha = -2
!         However, we only need to save the expansion coefficients for alpha >= 1
!         This is I_-1
          I_n = 0.d0
          N_n = 1.d0
!         This is I_0
          N_np1 = N_a(rcut_hard, -2)
          I_np1 = dsqrt(pi/2.d0) * atom_sigma_scaled * ( derf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - &
                                                  derf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1
!         Speed up the computation of these coefficients
          if( rcut_hard_in == rcut_soft_in )then
            C1 = 0.d0
          else
            C1 = s2 / dr * dexp(-0.5d0 * (rcut_soft - rj)**2 / s2)
          end if
          C2 = s2 / rcut_hard * dexp(-0.5d0 * rj**2 / s2)
          do n = -1, alpha_max_der
            C1 = C1 * dr
            C2 = C2 * rcut_hard
            N_np2 = N_a(rcut_hard, n)
!           This is I_alpha
            I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                    - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 &
                    + C1 / N_np2 &
                    - C2 / N_np2
            if(n > 0)then
              exp_coeff_temp1(n) = I_np2
            end if
            N_n = N_np1
            N_np1 = N_np2
            I_n = I_np1
            I_np1 = I_np2
          end do
!         Compute the contribution to the derivative for this part (excludes the amplitude)
          if( do_derivatives )then
            do n = 1, alpha_max
              exp_coeff_der_temp(n) = (rj - rcut_hard)/s2 * ( atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled &
                                  - 1.d0 ) * exp_coeff_temp1(n) + &
                                 rcut_hard*N_a(rcut_hard, n+1)/s2/N_a(rcut_hard, n) * ( 2.d0 * atom_sigma_scaling &
                                  * (rj - rcut_hard) / atom_sigma_scaled - 1.d0 ) * exp_coeff_temp1(n+1) + &
                                 atom_sigma_scaling*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_scaled**3/ &
                                  N_a(rcut_hard, n) * exp_coeff_temp1(n+2)
            end do
          end if
!         If the atom is less than 4 standard deviations away from the soft cutoff, we add
!         also this correction to the integrals. This corresponds to a Gaussian filter. We
!         integrate between rcut_soft and rcut_hard in this case
!
!         This explicit ".true." or ".false." logical statement is there for debugging. For
!         regular code use it can be set to .false.
          if( .false. .or. (rcut_soft - rj) < 4.d0*atom_sigma_scaled )then
            atom_sigma_f = atom_sigma_scaled * dr / nf / dsqrt(s2 + dr**2/nf**2)
            rj_f = (s2 * rcut_soft + dr**2/nf**2*rj) / (s2 + dr**2/nf**2)
!           We leave this here because in the future atom_sigma will be rj dependent
            sf2 = atom_sigma_f**2
!           The products of two Gaussians is a Gaussian, but we need to add a prefactor
            pref_f = dexp( -0.5d0 * (rcut_soft-rj)**2 / ( s2 + dr**2/nf**2 ) )
!           We have the recursion series starting at n = 0, which means alpha = -2
!           However, we only need to save the expansion coefficients for alpha >= 1
!           This is I_-1
            I_n = 0.d0
            N_n = 1.d0
!           This is I_0
            N_np1 = N_a(rcut_hard, -2)
            I_np1 = dsqrt(pi/2.d0) * atom_sigma_f * ( derf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - &
                                                      derf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1
!           Speed up the computation of these coefficients
            C2 = sf2 / dr * dexp(-0.5d0 * (rcut_soft - rj_f)**2 / sf2)
            do n = -1, alpha_max_der
              C2 = C2 * dr
              N_np2 = N_a(rcut_hard, n)
!             This is I_alpha
              I_np2 = sf2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                      - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 &
                      - C2 / N_np2
              if(n > 0)then
                exp_coeff_temp2(n) = I_np2
              end if
              N_n = N_np1
              N_np1 = N_np2
              I_n = I_np1
              I_np1 = I_np2
            end do
!           Compute the contribution to the derivative for this part (excludes the amplitude)
            if( do_derivatives )then
              denom = s2 + dr**2/nf**2
              der_pref_f = pref_f * ( (rcut_soft - rj) / denom + (rcut_soft - rj)**2 / denom**2 * &
                                      atom_sigma_scaled * atom_sigma_scaling )
              der_rjf_rj = (2.d0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + dr**2/nf**2) / denom &
                           - (s2*rcut_soft + dr**2/nf**2 * rj) * 2.d0 * atom_sigma_scaled * &
                             atom_sigma_scaling / denom**2
              der_sjf_rj = atom_sigma_scaling * dr/nf / dsqrt(denom) * (1.d0 - atom_sigma_scaled**2/denom)
              do n = 1, alpha_max
                exp_coeff_der_temp(n) = exp_coeff_der_temp(n) + &
                                      pref_f * ( &
                                      (rj_f - rcut_hard)/sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f &
                                        - der_rjf_rj ) * exp_coeff_temp2(n) + &
                                      rcut_hard*N_a(rcut_hard, n+1)/sf2/N_a(rcut_hard, n) * ( 2.d0 * der_sjf_rj &
                                        * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj ) * exp_coeff_temp2(n+1) + &
                                      der_sjf_rj*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_f**3/ &
                                      N_a(rcut_hard, n) * exp_coeff_temp2(n+2) ) + &
                                      der_pref_f * &
                                      exp_coeff_temp2(n)
              end do
            end if
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          if( do_derivatives )then
            exp_coeff_der_temp(1:alpha_max) = amplitude * exp_coeff_der_temp(1:alpha_max) + amplitude_der * &
                                              (exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max))
            exp_coeff_der(1:alpha_max, k) = matmul( W, exp_coeff_der_temp(1:alpha_max) )
          end if
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max) )
        end if
      end do
    end do




!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    ! exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    ! if( do_derivatives )then
    !   exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    ! end if
!   *********************************************

!   This is for debugging
    if( .false. )then
      open(10, file="radial_expansion_coefficients.dat", status="unknown", position="append")
      write(10,*) exp_coeff(1:alpha_max, 2)
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_derivatives.dat", status="unknown", position="append")
        write(10,*) exp_coeff_der(1:alpha_max, 2)
        close(10)
      end if
     end if

    deallocate( exp_coeff_temp1, exp_coeff_temp2, exp_coeff_der_temp )

  return
  end subroutine
!**************************************************************************










!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set augmented with a Gaussian function at the origin
! and the saparable Gaussian representation for the atomic sites.
!
  subroutine get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                          rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                          amplitude_scaling, nf, W, scaling_mode, mask, &
                                                          radial_enhancement, do_derivatives, exp_coeff, &
                                                          exp_coeff_der, k_idx)
!   Expansion coefficients using the polynomial basis with smooth filter plus a Gaussian centered at the origin
!
!   Apparently, with this basis the W matrix becomes complex for 13 and higher
!   basis functions, so alpha_max should not be higher than 12.
!   ADD ERROR MESSAGE ABOUT THIS IN THE CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   TRY OUT: Check for very small numbers (that could lead to negative numbers) and then truncate them
!   to zero
!
    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in,rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, ampli_tude
    logical, intent(in) :: mask(:), do_derivatives
    character(*), intent(in) :: scaling_mode
    real(c_double), intent(in) ::  rcut_hard_in
!
    integer :: n, i, j, k, alpha_max_der, kij
    real*8 :: I_n, I_np1, I_np2, pi, sq2, rj, N_n, N_np1, N_np2, C1, C2, dr, s2, sf2
    real*8 :: W(:,:)
    real*8 :: atom_sigma_f, rj_f, N_gauss
!   Results will be stored in exp_coeff, which is an array of dimension
!   (n_sites, n_neigh_max, alpha_max)
    real(c_double),intent(inout), target :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_temp1(:), exp_coeff_temp2(:), exp_coeff_der_temp(:)
    logical, save :: print_basis = .false.
    real*8 :: denom, der_sjf_rj, der_rjf_rj, ampli_tude_der, pref_f, der_pref_f, sigma_star
    integer(c_int), intent(in), target :: k_idx(:)
    !real(c_double),allocatable, target :: amplitude_save(:),  amplitude_der_save(:)
    !integer :: n_atom_pairs
    
    
!   If the user requests derivatives, we need to get the expansion coefficients up to
!   alpha_max - 1 + 2. The "-1" is there because the Gaussian basis at the origin does not
!   participate in the calculation of the derivatives for the polynomial basis functions
    if( do_derivatives )then
      alpha_max_der = alpha_max + 2
    else
      alpha_max_der = alpha_max
    end if
    allocate( exp_coeff_temp1(1:alpha_max_der) )
    allocate( exp_coeff_temp2(1:alpha_max_der) )
    allocate( exp_coeff_der_temp(1:alpha_max) )
    
!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
!    if( .false. .and. print_basis )then
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      write(*,*) "g(x,s) = exp(-0.5*x**2/s**2) * sqrt(2./s) / pi**0.25 "
      do j = 1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max-1
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,"(A,F16.10,A,E16.8,A)",advance="no") "g(x,", atom_sigma_in, ") *", W(j,alpha_max), "+"
        write(*,*) "0."
      end do
    end if

    pi = dacos(-1.d0)
    sq2 = dsqrt(2.d0)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if
!
!   **************** New basis ******************
!   Redefine all the distances by dividing them by rcut_hard
!   We do this to avoid numerical instability when the value
!   of alpha is too high
!
!    rcut_soft = rcut_soft_in
!    rcut_hard = rcut_hard_in
!    atom_sigma = atom_sigma_in
!    dr = rcut_hard - rcut_soft
!    N_gauss = dsqrt(2.d0/atom_sigma_in) / pi**0.25
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
    dr = 1.d0 - rcut_soft_in/rcut_hard_in
    N_gauss = dsqrt(2.d0/atom_sigma) / pi**0.25d0
!   *********************************************
!
    k = 0
    do i = 1, n_sites
       k = k_idx(i)
      do j = 1, n_neigh(i)
         k = k + 1
!       IMPORTANT: for this basis, we skip i itself, which is neighbor number 1
        if( j == 1 )then
          cycle
        end if
        ! if(rjs_in(k) <= rcut_hard_in) then
        ! if(mask(k).eqv..false.) then
        ! write(*,*) "Mask is false at ", k
        ! endif 
        ! endif 

        if( rjs_in(k) <= rcut_hard_in .and. mask(k) )then
          pref_f = 0.d0
          exp_coeff_temp1 = 0.d0
          exp_coeff_temp2 = 0.d0
          exp_coeff_der_temp = 0.d0
!   **************** New basis ******************
!          rj = rjs_in(k)
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
!         We leave this here because in the future atom_sigma will be rj dependent
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              ampli_tude = 1.d0 / atom_sigma_scaled
              ampli_tude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              ampli_tude = 0.d0
              ampli_tude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                ampli_tude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                ampli_tude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * ampli_tude
              else
                ampli_tude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                ampli_tude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * ampli_tude
              end if
            end if
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            ampli_tude_der = ampli_tude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            ampli_tude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            ampli_tude = ampli_tude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            ampli_tude_der = ampli_tude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            ampli_tude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            ampli_tude = ampli_tude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if
          !amplitude_save(k)=ampli_tude
          !amplitude_der_save(k)=ampli_tude_der
!         We have the recursion series starting at n = 0, which means alpha = -2
!         However, we only need to save the expansion coefficients for alpha >= 1
!         This is I_-1
          I_n = 0.d0
          N_n = 1.d0
!         This is I_0
          N_np1 = N_a(rcut_hard, -2)
          I_np1 = dsqrt(pi/2.d0) * atom_sigma_scaled * ( derf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - &
                                                  derf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1
!         Speed up the computation of these coefficients
          if( rcut_hard_in == rcut_soft_in )then
            C1 = 0.d0
          else
            C1 = s2 / dr * dexp(-0.5d0 * (rcut_soft - rj)**2 / s2)
          end if
          C2 = s2 / rcut_hard * dexp(-0.5d0 * rj**2 / s2)
!         This is different wrt the regular polynomial basis, we only go up to alpha_max-1
          do n = -1, alpha_max_der-1
            C1 = C1 * dr
            C2 = C2 * rcut_hard
            N_np2 = N_a(rcut_hard, n)
!           This is I_alpha
            I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                    - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 &
                    + C1 / N_np2 &
                    - C2 / N_np2
            if(n > 0)then
              exp_coeff_temp1(n) = I_np2
            end if
            N_n = N_np1
            N_np1 = N_np2
            I_n = I_np1
            I_np1 = I_np2
          end do
!         Compute the contribution to the derivative for this part (excludes the amplitude)
          if( do_derivatives )then
            do n = 1, alpha_max-1
              exp_coeff_der_temp(n) = (rj - rcut_hard)/s2 * ( atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled &
                                  - 1.d0 ) * exp_coeff_temp1(n) + &
                                 rcut_hard*N_a(rcut_hard, n+1)/s2/N_a(rcut_hard, n) * ( 2.d0 * atom_sigma_scaling &
                                  * (rj - rcut_hard) / atom_sigma_scaled - 1.d0 ) * exp_coeff_temp1(n+1) + &
                                 atom_sigma_scaling*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_scaled**3/ &
                                  N_a(rcut_hard, n) * exp_coeff_temp1(n+2)
            end do
          end if
!         If the atom is less than 4 standard deviations away from the soft cutoff, we add
!         also this correction to the integrals. This corresponds to a Gaussian filter. We
!         integrate between rcut_soft and rcut_hard in this case
!
!         This explicit ".true." or ".false." logical statement is there for debugging. For
!         regular code use it can be set to .false.
          if( .false. .or. (rcut_soft - rj) < 4.d0*atom_sigma_scaled )then
            atom_sigma_f = atom_sigma_scaled * dr / nf / dsqrt(s2 + dr**2/nf**2)
            rj_f = (s2 * rcut_soft + dr**2/nf**2*rj) / (s2 + dr**2/nf**2)
!           We leave this here because in the future atom_sigma will be rj dependent
            sf2 = atom_sigma_f**2
!           The products of two Gaussians is a Gaussian, but we need to add a prefactor
            pref_f = dexp( -0.5d0 * (rcut_soft-rj)**2 / ( s2 + dr**2/nf**2 ) )
!           We have the recursion series starting at n = 0, which means alpha = -2
!           However, we only need to save the expansion coefficients for alpha >= 1
!           This is I_-1
            I_n = 0.d0
            N_n = 1.d0
!           This is I_0
            N_np1 = N_a(rcut_hard, -2)
            I_np1 = dsqrt(pi/2.d0) * atom_sigma_f * ( derf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - &
                                                      derf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1
!           Speed up the computation of these coefficients
            C2 = sf2 / dr * dexp(-0.5d0 * (rcut_soft - rj_f)**2 / sf2)
!           This is different wrt the regular polynomial basis, we only go up to alpha_max-1
            do n = -1, alpha_max_der-1
              C2 = C2 * dr
              N_np2 = N_a(rcut_hard, n)
!             This is I_alpha
              I_np2 = sf2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                      - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 &
                      - C2 / N_np2
              if(n > 0)then
                exp_coeff_temp2(n) = exp_coeff_temp2(n) + I_np2
              end if
              N_n = N_np1
              N_np1 = N_np2
              I_n = I_np1
              I_np1 = I_np2
            end do
!           Compute the contribution to the derivative for this part (excludes the amplitude)
            if( do_derivatives )then
              denom = s2 + dr**2/nf**2
              der_pref_f = pref_f * ( (rcut_soft - rj) / denom + (rcut_soft - rj)**2 / denom**2 * &
                                      atom_sigma_scaled * atom_sigma_scaling )
              der_rjf_rj = (2.d0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + dr**2/nf**2) / denom &
                           - (s2*rcut_soft + dr**2/nf**2 * rj) * 2.d0 * atom_sigma_scaled * &
                             atom_sigma_scaling / denom**2
              der_sjf_rj = atom_sigma_scaling * dr/nf / dsqrt(denom) * (1.d0 - atom_sigma_scaled**2/denom)
              do n = 1, alpha_max-1
                exp_coeff_der_temp(n) = exp_coeff_der_temp(n) + &
                                      pref_f * ( &
                                      (rj_f - rcut_hard)/sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f &
                                        - der_rjf_rj ) * exp_coeff_temp2(n) + &
                                      rcut_hard*N_a(rcut_hard, n+1)/sf2/N_a(rcut_hard, n) * ( 2.d0 * der_sjf_rj &
                                        * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj ) * exp_coeff_temp2(n+1) + &
                                      der_sjf_rj*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_f**3/ &
                                      N_a(rcut_hard, n) * exp_coeff_temp2(n+2) ) + &
                                      der_pref_f * &
                                      exp_coeff_temp2(n)
              end do
            end if
          end if
!         Now we obtain the overlap integral between the atomic density and the Gaussian centered at the origin
!         We assume that at the soft cutoff the Gaussian basis function is approx. zero, and we only
!         compute the overlap coefficient if the atom is close to the origin
!         Note that the sigma for the Gaussian basis function and the atom's sigma are not the same if there is
!         sigma scaling
!         We reset these to zero (temp2 will stay zero since it's the filter one, and we do not apply filter to
!         the overlap with the central Gaussian. This should be a good approximation if rcut_soft >= 4 * atom_sigma
          exp_coeff_temp1(alpha_max) = 0.d0
          exp_coeff_temp2(alpha_max) = 0.d0
          if( .false. .or. rj < 4.d0*(atom_sigma+atom_sigma_scaled) )then
              sigma_star = dsqrt(atom_sigma**2 + s2)
              exp_coeff_temp1(alpha_max) = dexp(- 0.5d0 * rj**2 / sigma_star**2 ) * dsqrt(pi/2.d0) * &
                                        atom_sigma_scaled*atom_sigma / sigma_star * ( 1.d0 &
                                        + derf(atom_sigma/atom_sigma_scaled*rj/sq2/sigma_star) )* N_gauss
            if( do_derivatives )then
              exp_coeff_der_temp(alpha_max) = ( rj**2 * atom_sigma_scaling / atom_sigma_scaled**3 - rj/sigma_star**2 + &
                                                atom_sigma_scaling*rj**2*atom_sigma**4/atom_sigma_scaled**3/sigma_star**4 + &
                                                atom_sigma_scaling*atom_sigma**2/atom_sigma_scaled/sigma_star**2 - &
                                                2.d0*rj**2*atom_sigma_scaling*atom_sigma**2/atom_sigma_scaled**3/sigma_star**2 &
                                              ) * exp_coeff_temp1(alpha_max) + &
                                              (1./s2 - 2.d0*rj*atom_sigma_scaling/atom_sigma_scaled**3) * s2 * atom_sigma**2 / &
                                                sigma_star**2 * dsqrt(2.d0/atom_sigma) / pi**0.25d0 * &
                                                dexp(-0.5d0 * rj**2 / sigma_star**2 * (1.d0 + atom_sigma**2 / s2) ) + &
                                              dsqrt(2.d0/atom_sigma) / pi**0.25d0 * dexp(-0.5d0 * rj**2 / sigma_star**2 * &
                                                (1.d0 + atom_sigma**2 / s2) ) * atom_sigma_scaling / atom_sigma_scaled * &
                                                rj*atom_sigma**4/sigma_star**4
            end if
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          if( do_derivatives )then
            exp_coeff_der_temp(1:alpha_max) = ampli_tude * exp_coeff_der_temp(1:alpha_max) + ampli_tude_der * &
                                              (exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max))
            exp_coeff_der(1:alpha_max, k) = matmul( W, exp_coeff_der_temp(1:alpha_max) )
          end if
          exp_coeff(1:alpha_max, k) = ampli_tude * matmul( W, exp_coeff_temp1(1:alpha_max) &
          + pref_f * exp_coeff_temp2(1:alpha_max) )
          !write(*,*) "alpha_max", alpha_max
          !stop
        end if
        ! do n=1,alpha_max
        !   if(isnan(exp_coeff(n,k)).or.isnan(exp_coeff_temp1(n)).or.isnan(exp_coeff_temp2(n)).or.isnan(pref_f)) then
        !     write(*,*)
        !     write(*,*) "Nan allert! in radial subroutine",n,k,i,j
        !     write(*,*) pref_f , exp_coeff(n,k)
        !     stop
        !   endif
        !  enddo
      !   if(k.eq.19855) then
      !     open(unit=666, position=append)
      !     write(666,*) "Test"
      !     write(666,*) exp_coeff(3,k), exp_coeff(6,k), pref_f,exp_coeff_temp1(3),exp_coeff_temp1(6),exp_coeff_temp2(3),exp_coeff_temp2(6)
      !     close(666)
      !   endif
       end do
    end do
!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    ! exp_coeff = exp_coeff*dsqrt(rcut_hard_in)
    ! if( do_derivatives )then
    !   exp_coeff_der = exp_coeff_der/ dsqrt(rcut_hard_in)
    ! end if
!   ***********************************
    
  
    ! k=0
    ! do i = 1, n_sites
    !   do j = 1, n_neigh(i)
    !     k = k + 1
    !     do n=1,alpha_max
    !       if(isnan(exp_coeff(n,k))) then
    !         write(*,*)
    !         write(*,*) "Nan allert! in radial subroutine",n,k,i
    !         write(*,*) pref_f, exp_coeff(n,k)
    !       endif
    !     enddo
    !   enddo
    ! enddo
!   This is for debugging
    if( .false. )then
      open(10, file="coefficients.dat", status="unknown", position="append")
      write(10,*) exp_coeff(1:alpha_max, 1)
      close(10)
      if( do_derivatives )then
        open(10, file="derivatives.dat", status="unknown", position="append")
        write(10,*) exp_coeff_der(1:alpha_max, 1)
        close(10)
      end if
    end if

    deallocate( exp_coeff_temp1, exp_coeff_temp2, exp_coeff_der_temp )
    !deallocate(amplitude_der_save,amplitude_save)

  return
  end subroutine
!**************************************************************************










!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3(alpha_max, S, W)

    implicit none

    integer :: alpha_max, i, j, info
    real*8, intent(inout) :: W(:,:), S(:,:)
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:)
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

    do i = 1, alpha_max
      S(i,i) = 1.d0
      do j = i+1, alpha_max
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
    S = matmul(U,S)
    S = matmul(S,VT)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv )

  return
  end subroutine
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set augmented with a central Gaussian function. It is
! also needed to transform between original basis and orthonormal basis.
! It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3gauss(alpha_max, atom_sigma_in, rcut_hard_in, S, W)

    implicit none

    integer :: alpha_max, i, j, info, n
    real*8, intent(inout) :: W(:,:), S(:,:)
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:)
    real*8 :: s2, I_n, N_n, N_np1, I_np1, N_np2, I_np2, C2, sq2, pi, atom_sigma, rcut_hard
    real*8, intent(in) :: rcut_hard_in, atom_sigma_in
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

!   These are the overlap integrals for the polynomial functions
    do i = 1, alpha_max-1
      S(i,i) = 1.d0
      do j = i+1, alpha_max-1
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

!   These are the overlap integrals between the Gaussian and the polynomials
!   See derivation of radial expansion coefficients to understand this code
!   **** New basis ****
    atom_sigma = atom_sigma_in/rcut_hard_in
    rcut_hard = 1.d0
!   *******************
    s2 = atom_sigma**2
    sq2 = dsqrt(2.d0)
    pi = dacos(-1.d0)
    I_n = 0.d0
    N_n = 1.d0
    N_np1 = N_a(rcut_hard, -2)
    I_np1 = dsqrt(pi/2.d0) * atom_sigma * derf( rcut_hard/sq2/atom_sigma ) / N_np1
    C2 = s2 / rcut_hard
    do n = -1, alpha_max-1
      C2 = C2 * rcut_hard
      N_np2 = N_a(rcut_hard, n)
      I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
              + N_np1 * rcut_hard / N_np2 * I_np1 &
              - C2 / N_np2
      if(n > 0)then
!       Include the normalization factor of the Gaussian
        S(alpha_max, n) = I_np2 * sq2 / dsqrt(atom_sigma) / pi**0.25d0
        S(n, alpha_max) = S(alpha_max, n)
      end if
      N_n = N_np1
      N_np1 = N_np2
      I_n = I_np1
      I_np1 = I_np2
    end do
    S(alpha_max, alpha_max) = 1.d0

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
    S = matmul(U,S)
    S = matmul(S,VT)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv )

  return
  end subroutine
!**************************************************************************




end module soap_turbo_radial
