! Two isolated first/second-derivative checks, upstream of anything the
! central-atom Hessian assembles:
!
!   A) the radial expansion coefficients: exp_coeff_der and exp_coeff_der2
!      against central differences of exp_coeff in rj.
!
!   B) F_l(r) = amplitude(r) * ilexp(l, r/sigma), the radial part of the
!      angular coefficient, and its first two derivatives -- the chain rule
!      through z = (r/sigma)^2 with a sigma that itself depends on r.
!
! Both are pure functions of one variable, so a converging h-scan is a
! complete test of them.
!
program rdcheck
   use soap_turbo_radial
   use soap_turbo_angular
   use soap_turbo_functions
   implicit none
   integer, parameter :: dp = kind(0.d0)
   integer, parameter :: amax = 8, lmax = 8
   integer :: nr, i, l, ih
   real(dp) :: rcut_h, rcut_s, sig_r, sig_r_s, ampsc, nfv, h
   real(dp), allocatable :: rjs(:), ec(:, :), ecd(:, :), ecd2(:, :)
   real(dp), allocatable :: ecp(:, :), ecdp(:, :), ecm(:, :), ecdm(:, :)
   real(dp), allocatable :: W(:, :), S(:, :)
   logical, allocatable :: mask(:)
   integer, allocatable :: nn(:)
   real(dp) :: w1, w2, fd, e1, e2

   legacy_filter_seed = .false.
   rcut_h = 4.5d0
   rcut_s = 4.0d0
   sig_r = 0.5d0
   sig_r_s = 0.05d0
   ampsc = 1.d0
   nfv = 4.d0

   nr = 40
   allocate (rjs(nr), nn(1), mask(nr))
   allocate (ec(amax, nr), ecd(amax, nr), ecd2(amax, nr))
   allocate (ecp(amax, nr), ecdp(amax, nr), ecm(amax, nr), ecdm(amax, nr))
   allocate (W(amax, amax), S(amax, amax))
   call get_orthonormalization_matrix_poly3gauss(amax, sig_r, rcut_h, S, W)
   nn(1) = nr
   mask = .true.
   do i = 1, nr
      rjs(i) = 0.15d0 + 4.2d0*dble(i - 1)/dble(nr - 1)
   end do

   write (*, '(A)') "A) radial expansion coefficients, poly3gauss, alpha_max = 8"
   write (*, '(A)') "   h          max rel err d/drj      max rel err d2/drj2"
   do ih = 1, 5
      h = 1.d-2*10.d0**(-(ih - 1))
      call radial(rjs, ec, ecd, ecd2)
      call radial(rjs + h, ecp, ecdp, ecd2)
      call radial(rjs - h, ecm, ecdm, ecd2)
      call radial(rjs, ec, ecd, ecd2)
      w1 = 0.d0; w2 = 0.d0
      e1 = maxval(abs(ecd)); e2 = maxval(abs(ecd2))
      do i = 1, nr
         do l = 1, amax
            fd = (ecp(l, i) - ecm(l, i))/(2.d0*h)
            w1 = max(w1, abs(fd - ecd(l, i))/e1)
            fd = (ecdp(l, i) - ecdm(l, i))/(2.d0*h)
            w2 = max(w2, abs(fd - ecd2(l, i))/e2)
         end do
      end do
      write (*, '(3X,ES8.1,2(6X,ES12.4))') h, w1, w2
   end do

   write (*, *)
   write (*, '(A)') "B) F_l(r) = amplitude(r)*ilexp(l, r/sigma(r)), l = 0..8"
   write (*, '(A)') "   h          max rel err dF/dr       max rel err d2F/dr2"
   do ih = 1, 5
      h = 1.d-2*10.d0**(-(ih - 1))
      call fcheck(h, w1, w2)
      write (*, '(3X,ES8.1,2(6X,ES12.4))') h, w1, w2
   end do

contains

   subroutine radial(r, a, b, c)
      real(dp), intent(in) :: r(:)
      real(dp), intent(out) :: a(:, :), b(:, :), c(:, :)
      a = 0.d0; b = 0.d0; c = 0.d0
      call get_radial_expansion_coefficients_poly3gauss(1, nn, r, amax, rcut_s, rcut_h, sig_r, sig_r_s, &
                                                        ampsc, nfv, W, "polynomial", mask, 1, .true., a, b, c)
   end subroutine radial

! F and its two derivatives, exactly as pair_geometry in soap_turbo.f90 forms
! them, checked against central differences of the value.
   subroutine fcheck(h, e1, e2)
      real(dp), intent(in) :: h
      real(dp), intent(out) :: e1, e2
      real(dp) :: f0(0:lmax), f1(0:lmax), f2(0:lmax)
      real(dp) :: g0(0:lmax), g1(0:lmax), g2(0:lmax)
      real(dp) :: m0(0:lmax), m1(0:lmax), m2(0:lmax)
      real(dp) :: r, sc, fd
      integer :: ii, ll
      e1 = 0.d0; e2 = 0.d0
      do ii = 1, 30
         r = 0.3d0 + 4.0d0*dble(ii - 1)/29.d0
         call fval(r, f0, f1, f2)
         call fval(r + h, g0, g1, g2)
         call fval(r - h, m0, m1, m2)
         sc = maxval(abs(f1)); if (sc == 0.d0) sc = 1.d0
         do ll = 0, lmax
            fd = (g0(ll) - m0(ll))/(2.d0*h)
            e1 = max(e1, abs(fd - f1(ll))/sc)
         end do
         sc = maxval(abs(f2)); if (sc == 0.d0) sc = 1.d0
         do ll = 0, lmax
            fd = (g1(ll) - m1(ll))/(2.d0*h)
            e2 = max(e2, abs(fd - f2(ll))/sc)
         end do
      end do
   end subroutine fcheck

   subroutine fval(r, fa, fb, fc)
      real(dp), intent(in) :: r
      real(dp), intent(out) :: fa(0:), fb(0:), fc(0:)
      real(dp) :: sig, ss, amp, ampd, ampdd, x, z, invz, zp, zpp, rl
      real(dp) :: il(0:lmax + 2), dz(0:lmax), d2z(0:lmax), fact(lmax + 2)
      integer :: ll
      ss = 0.1d0                     ! atom_sigma_t_scaling
      sig = 0.2d0 + ss*r              ! atom_sigma_t + scaling*r
      amp = 4.5d0**2/sig**2
      ampd = -2.d0*ss/sig*amp
      ampdd = 6.d0*ss**2/sig**2*amp
      x = r/sig
      z = x*x
      invz = 1.d0/z
      call get_semifactorial_array(fact, lmax + 2)
      call get_ilexp(il, fact, lmax + 2, x)
      do ll = 0, lmax
         rl = dble(ll)
         dz(ll) = il(ll + 1) + (rl*invz - 1.d0)*il(ll)
         d2z(ll) = il(ll + 2) + (dble(2*ll + 1)*invz - 2.d0)*il(ll + 1) &
                   + (rl*(rl - 1.d0)*invz*invz - 2.d0*rl*invz + 1.d0)*il(ll)
      end do
      zp = 2.d0*r/sig**2 - 2.d0*r**2*ss/sig**3
      zpp = 2.d0/sig**2 - 8.d0*r*ss/sig**3 + 6.d0*r**2*ss**2/sig**4
      do ll = 0, lmax
         fa(ll) = amp*il(ll)
         fb(ll) = ampd*il(ll) + amp*(dz(ll)*zp)
         fc(ll) = ampdd*il(ll) + 2.d0*ampd*(dz(ll)*zp) + amp*(d2z(ll)*zp**2 + dz(ll)*zpp)
      end do
   end subroutine fval

end program rdcheck
