! Independent check of the poly3operator radial expansion coefficients.
!
! The code under test evaluates
!
!    c_alpha(rj) = amplitude(rj) * Int_0^1 rho(r; rj) (1-r)**(alpha+2) sqrt(2 alpha + 5) dr
!
! in closed form, by integrating the piecewise-cubic atomic density rho against
! the poly3 basis by parts. This program evaluates the same integral by
! Gauss-Legendre quadrature in real128 instead, and compares.
!
! The quadrature is not an approximation. rho is cubic on each of its pieces and
! the filter that multiplies it above rcut_soft is cubic too, so the integrand
! is a polynomial of degree at most alpha + 8 on every panel between its
! breakpoints; 16-node Gauss-Legendre is exact through degree 31. Splitting the
! range at rj - width, rj, rcut_soft and rj + width therefore gives the exact
! integral to real128 rounding, and the difference from the double-precision
! closed form is that form's own conditioning, nothing else.
!
! Input is radial_exp_coeff_dump.dat, written by either build under
! TURBOGAP_DUMP_RADIAL. The dump carries the post-W, post-global-scaling
! coefficients, so the reference applies the same orthonormalisation matrix --
! which is checked here against the analytic overlap matrix rather than trusted.
!
! Usage:
!   poly3operatorverify <dump> <alpha_max> <rcut_hard> <rcut_soft> <atom_sigma>
!                       <amplitude_scaling> <radial_enhancement> <central_weight>
!                       [tol]
program poly3operatorverify

   use soap_turbo_radial, only: get_orthonormalization_matrix_poly3_tabulated

   implicit none

   integer, parameter :: dp = selected_real_kind(15, 307)
   integer, parameter :: qp = selected_real_kind(30, 291)
   integer, parameter :: ngl = 16

   real(qp) :: xgl(1:ngl)
   real(qp) :: wgl(1:ngl)

   character(len=512) :: dump_file
   character(len=64) :: arg
   integer :: alpha_max
   integer :: radial_enhancement
   real(dp) :: rcut_hard
   real(dp) :: rcut_soft
   real(dp) :: atom_sigma
   real(dp) :: amplitude_scaling
   real(dp) :: central_weight
   real(dp) :: tol

   real(dp), allocatable :: W(:, :)
   real(dp), allocatable :: S(:, :)
   real(qp), allocatable :: Wq(:, :)
   real(qp), allocatable :: ref(:)
   real(qp), allocatable :: ref_p(:)
   real(qp), allocatable :: ref_m(:)

   integer :: unit_in
   integer :: iostatus
   integer :: k2
   integer :: n
   integer :: alpha
   integer :: n_checked
   integer :: n_lines
   real(dp) :: rj
   real(dp) :: val
   real(dp) :: der
   real(dp) :: rj_prev
   real(dp) :: err
   real(dp) :: worst_val
   real(dp) :: worst_der
   real(dp) :: worst_val_rj
   real(dp) :: worst_der_rj
   real(dp) :: scale_val
   real(dp) :: scale_der
   real(dp) :: h
   integer, parameter :: n_scan = 3
   integer :: n_scan_have
   real(dp) :: scan_rj(1:n_scan)
   real(dp) :: scan_der(1:n_scan)
   integer :: scan_alpha(1:n_scan)
   character(len=256) :: line
   logical :: ok

   call gauss_legendre(xgl, wgl)

   if (command_argument_count() < 8) then
      write (*, *) "usage: poly3operatorverify <dump> <alpha_max> <rcut_hard> <rcut_soft> &
                   &<atom_sigma> <amplitude_scaling> <radial_enhancement> <central_weight> [tol]"
      stop 1
   end if
   call get_command_argument(1, dump_file)
   call get_command_argument(2, arg); read (arg, *) alpha_max
   call get_command_argument(3, arg); read (arg, *) rcut_hard
   call get_command_argument(4, arg); read (arg, *) rcut_soft
   call get_command_argument(5, arg); read (arg, *) atom_sigma
   call get_command_argument(6, arg); read (arg, *) amplitude_scaling
   call get_command_argument(7, arg); read (arg, *) radial_enhancement
   call get_command_argument(8, arg); read (arg, *) central_weight
   tol = 1.0e-9_dp
   if (command_argument_count() >= 9) then
      call get_command_argument(9, arg); read (arg, *) tol
   end if

   allocate (W(1:alpha_max, 1:alpha_max))
   allocate (S(1:alpha_max, 1:alpha_max))
   allocate (Wq(1:alpha_max, 1:alpha_max))
   allocate (ref(1:alpha_max))
   allocate (ref_p(1:alpha_max))
   allocate (ref_m(1:alpha_max))
   W = 0.0_dp
   S = 0.0_dp
   call get_orthonormalization_matrix_poly3_tabulated(alpha_max, S, W)
   Wq = real(W, qp)

   ok = .true.
   write (*, '(A)') "== check 1: the tabulated basis matrices, against the analytic overlap =="
   call check_basis_matrices(S, W, alpha_max, ok)
   if (.not. ok) then
      write (*, '(A)') "RESULT: FAIL"
      stop 1
   end if

   write (*, '(A)') ""
   write (*, '(A)') "== check 2: coefficients and their radial derivatives, against quadrature =="
   worst_val = 0.0_dp
   worst_der = 0.0_dp
   worst_val_rj = 0.0_dp
   worst_der_rj = 0.0_dp
   scale_val = 0.0_dp
   scale_der = 0.0_dp
   n_checked = 0
   n_lines = 0
   n_scan_have = 0
   rj_prev = -1.0_dp
   h = 1.0e-10_dp

   open (newunit=unit_in, file=trim(dump_file), status="old", action="read", iostat=iostatus)
   if (iostatus /= 0) then
      write (*, '(A,A)') "cannot open ", trim(dump_file)
      stop 1
   end if
   do
      read (unit_in, '(A)', iostat=iostatus) line
      if (iostatus /= 0) exit
      if (line(1:1) == "#") cycle
      read (line, *, iostat=iostatus) k2, n, rj, val, der
      if (iostatus /= 0) cycle
      n_lines = n_lines + 1
      if (val == 0.0_dp .and. der == 0.0_dp) cycle
      alpha = mod(n - 1, alpha_max) + 1

!     The reference is the same for every alpha at a given rj, so evaluate the
!     whole vector once and reuse it across the alpha_max lines that follow.
      if (rj /= rj_prev) then
         call reference_vector(real(rj, qp), ref)
         call reference_vector(real(rj, qp) + real(h, qp), ref_p)
         call reference_vector(real(rj, qp) - real(h, qp), ref_m)
         rj_prev = rj
      end if

      scale_val = max(scale_val, abs(val))
      scale_der = max(scale_der, abs(der))
      err = abs(val - real(ref(alpha), dp))
      if (err > worst_val) then
         worst_val = err
         worst_val_rj = rj
      end if
!     Skip the finite difference within h of a kink in rj, where the analytic
!     derivative is right and the central difference straddles two branches.
      if (near_kink(rj, 1.0e-6_dp)) cycle
      err = abs(der - real((ref_p(alpha) - ref_m(alpha))/(2.0_qp*real(h, qp)), dp))
      if (err > worst_der) then
         worst_der = err
         worst_der_rj = rj
      end if
      n_checked = n_checked + 1
!     Keep a few well-separated samples for the h-scan below.
      if (n_scan_have < n_scan .and. alpha == 1) then
         if (n_scan_have == 0) then
            if (rj > 0.25_dp*rcut_hard .and. rj < 0.40_dp*rcut_hard) call keep_scan(rj, der, alpha)
         else if (n_scan_have == 1) then
            if (rj > 0.60_dp*rcut_hard .and. rj < 0.75_dp*rcut_hard) call keep_scan(rj, der, alpha)
         else
            if (rj > 0.90_dp*rcut_hard .and. rj < 0.97_dp*rcut_hard) call keep_scan(rj, der, alpha)
         end if
      end if
   end do
   close (unit_in)

   write (*, '(A,I0,A,I0,A)') "  ", n_checked, " nonzero coefficients checked (", n_lines, " dump lines)"
   write (*, '(A,ES12.4,A,ES12.4,A,ES12.4)') "  value      maxabsdiff=", worst_val, &
      "  max|ref|=", scale_val, "  rel=", worst_val/max(scale_val, tiny(1.0_dp))
   write (*, '(A,F8.4)') "                                            worst at rj=", worst_val_rj
   write (*, '(A,ES12.4,A,ES12.4,A,ES12.4)') "  derivative maxabsdiff=", worst_der, &
      "  max|ref|=", scale_der, "  rel=", worst_der/max(scale_der, tiny(1.0_dp))
   write (*, '(A,F8.4)') "                                            worst at rj=", worst_der_rj

   ok = (worst_val/max(scale_val, tiny(1.0_dp)) <= tol) .and. &
        (worst_der/max(scale_der, tiny(1.0_dp)) <= tol)

   write (*, '(A)') ""
   write (*, '(A)') "== check 3: h-scan of the derivative, error must fall ~4x per halving =="
   if (n_scan_have == 0) then
      write (*, '(A)') "  SKIP: the dump held no sample away from a kink in rj"
      ok = .false.
   end if
   do k2 = 1, n_scan_have
      call h_scan(scan_rj(k2), scan_der(k2), scan_alpha(k2), ok)
   end do

   write (*, '(A)') ""
   if (ok) then
      write (*, '(A)') "RESULT: PASS"
   else
      write (*, '(A)') "RESULT: FAIL"
      stop 1
   end if

contains

!  True when rj sits within eps of a point where a limit of integration
!  switches branch, and the coefficient is only C0 in rj.
   function near_kink(rj, eps) result(is_near)
      real(dp), intent(in) :: rj
      real(dp), intent(in) :: eps
      logical :: is_near
      real(dp) :: width

      width = 2.0_dp*sqrt(2.0_dp*log(2.0_dp))*atom_sigma
      is_near = abs(rj) < eps .or. &
                abs(rj - rcut_soft) < eps .or. &
                abs(rj - width - rcut_soft) < eps .or. &
                abs(rj + width - rcut_soft) < eps .or. &
                abs(rj - width) < eps .or. &
                abs(rj + width - rcut_hard) < eps .or. &
                abs(rj - rcut_hard) < eps
   end function near_kink

!  The full coefficient vector at one rj: quadrature, amplitude, W, and the
!  sqrt(rcut_hard) that the change of variable leaves behind.
   subroutine reference_vector(rj_in, c)
      real(qp), intent(in) :: rj_in
      real(qp), intent(out) :: c(:)
      real(qp) :: raw(1:alpha_max)
      real(qp) :: amp
      integer :: a
      integer :: b

      call raw_coefficients(rj_in, raw)
      amp = amplitude(rj_in)
      do a = 1, alpha_max
         c(a) = 0.0_qp
         do b = 1, alpha_max
            c(a) = c(a) + Wq(a, b)*amp*raw(b)
         end do
         c(a) = c(a)*sqrt(real(rcut_hard, qp))
      end do
   end subroutine reference_vector

!  Int_0^1 rho(r) (1-r)**(alpha+2) sqrt(2 alpha + 5) dr, panel by panel, with
!  everything in units of rcut_hard.
   subroutine raw_coefficients(rj_in, raw)
      real(qp), intent(in) :: rj_in
      real(qp), intent(out) :: raw(:)
      real(qp) :: rj
      real(qp) :: width
      real(qp) :: rcs
      real(qp) :: fw
      real(qp) :: brk(1:4)
      real(qp) :: lo
      real(qp) :: hi
      real(qp) :: mid
      real(qp) :: half
      real(qp) :: r
      real(qp) :: rho
      integer :: nbrk
      integer :: i
      integer :: g
      integer :: a

      rj = rj_in/real(rcut_hard, qp)
      rcs = real(rcut_soft, qp)/real(rcut_hard, qp)
      width = 2.0_qp*sqrt(2.0_qp*log(2.0_qp))*real(atom_sigma, qp)/real(rcut_hard, qp)
      fw = 2.0_qp*sqrt(2.0_qp*log(2.0_qp))*(1.0_qp - rcs)

      raw = 0.0_qp
      lo = max(0.0_qp, rj - width)
      hi = min(1.0_qp, rj + width)
      if (hi <= lo) return

      nbrk = 0
      call push(brk, nbrk, lo)
      call push(brk, nbrk, min(max(rj, lo), hi))
      if (rcs > lo .and. rcs < hi) call push(brk, nbrk, rcs)
      call push(brk, nbrk, hi)
      call sort_ascending(brk, nbrk)

      do i = 1, nbrk - 1
         if (brk(i + 1) <= brk(i)) cycle
         mid = 0.5_qp*(brk(i) + brk(i + 1))
         half = 0.5_qp*(brk(i + 1) - brk(i))
         do g = 1, ngl
            r = mid + half*xgl(g)
            rho = bump(r, rj, width)
            if (r > rcs) rho = rho*filter(r, rcs, fw)
            do a = 1, alpha_max
               raw(a) = raw(a) + half*wgl(g)*rho*(1.0_qp - r)**(a + 2)*sqrt(real(2*a + 5, qp))
            end do
         end do
      end do
   end subroutine raw_coefficients

!  The piecewise-cubic atomic density: unit height at rj, zero and C1 at the
!  edges of [rj - width, rj + width].
   function bump(r, rj, width) result(f)
      real(qp), intent(in) :: r
      real(qp), intent(in) :: rj
      real(qp), intent(in) :: width
      real(qp) :: f
      real(qp) :: x

      x = (r - rj)/width
      if (x < -1.0_qp .or. x > 1.0_qp) then
         f = 0.0_qp
      else if (x <= 0.0_qp) then
         f = 1.0_qp - 3.0_qp*x*x - 2.0_qp*x*x*x
      else
         f = 1.0_qp - 3.0_qp*x*x + 2.0_qp*x*x*x
      end if
   end function bump

!  The cutoff filter that multiplies the density above rcut_soft.
   function filter(r, rcs, fw) result(f)
      real(qp), intent(in) :: r
      real(qp), intent(in) :: rcs
      real(qp), intent(in) :: fw
      real(qp) :: f
      real(qp) :: y

      y = (r - rcs)/fw
      f = 1.0_qp - 3.0_qp*y*y + 2.0_qp*y*y*y
   end function filter

   function amplitude(rj_in) result(amp)
      real(qp), intent(in) :: rj_in
      real(qp) :: amp
      real(qp) :: rj
      real(qp) :: sigma
      real(qp) :: pi

      pi = acos(-1.0_qp)
      rj = rj_in/real(rcut_hard, qp)
      sigma = real(atom_sigma, qp)/real(rcut_hard, qp)
      if (amplitude_scaling == 0.0_dp) then
         amp = 1.0_qp/sigma
      else
         amp = (1.0_qp/sigma)*(1.0_qp + 2.0_qp*rj**3 - 3.0_qp*rj**2)**real(amplitude_scaling, qp)
      end if
      if (rj == 0.0_qp) amp = amp*real(central_weight, qp)
      if (radial_enhancement == 1) then
         amp = amp*(rj + sqrt(2.0_qp/pi)*sigma)
      else if (radial_enhancement == 2) then
         amp = amp*(rj*rj + sigma*sigma + sqrt(8.0_qp/pi)*sigma*rj)
      end if
   end function amplitude

!  S(a,b) = Int_0^1 (1-r)**(a+2) (1-r)**(b+2) sqrt(2a+5) sqrt(2b+5) dr,
!  and W must be its inverse square root.
   subroutine check_basis_matrices(S_in, W_in, am, all_ok)
      real(dp), intent(in) :: S_in(:, :)
      real(dp), intent(in) :: W_in(:, :)
      integer, intent(in) :: am
      logical, intent(inout) :: all_ok
      real(qp) :: Sa(1:am, 1:am)
      real(qp) :: T(1:am, 1:am)
      real(qp) :: e
      real(qp) :: worst_s
      real(qp) :: worst_i
      integer :: a
      integer :: b
      integer :: c

      worst_s = 0.0_qp
      do a = 1, am
         do b = 1, am
            Sa(a, b) = sqrt(real((2*a + 5)*(2*b + 5), qp))/real(a + b + 5, qp)
            worst_s = max(worst_s, abs(Sa(a, b) - real(S_in(a, b), qp)))
         end do
      end do

      T = 0.0_qp
      do a = 1, am
         do b = 1, am
            do c = 1, am
               T(a, b) = T(a, b) + real(W_in(a, c), qp)*Sa(c, b)
            end do
         end do
      end do
      worst_i = 0.0_qp
      do a = 1, am
         do b = 1, am
            e = 0.0_qp
            do c = 1, am
               e = e + T(a, c)*real(W_in(c, b), qp)
            end do
            if (a == b) e = e - 1.0_qp
            worst_i = max(worst_i, abs(e))
         end do
      end do

      write (*, '(A,ES12.4)') "  |S_tabulated - S_analytic|_max = ", real(worst_s, dp)
      write (*, '(A,ES12.4)') "  |W S W - I|_max                = ", real(worst_i, dp)
      if (worst_s > 1.0e-13_qp .or. worst_i > 1.0e-8_qp) then
         write (*, '(A)') "  the tabulated matrices are not S and S**(-1/2) for this alpha_max"
         all_ok = .false.
      end if
   end subroutine check_basis_matrices

!  The analytic derivative under test, against central differences of the
!  quadrature reference at shrinking h. Second-order convergence is the
!  assertion: a derivative wrong by a constant factor passes any single-h
!  tolerance, but cannot make this ratio approach 4.
   subroutine h_scan(rj, dana, alpha_in, all_ok)
      real(dp), intent(in) :: rj
      real(dp), intent(in) :: dana
      integer, intent(in) :: alpha_in
      logical, intent(inout) :: all_ok
      integer, parameter :: nh = 5
      real(qp) :: cp(1:alpha_max)
      real(qp) :: cm(1:alpha_max)
      real(dp) :: hh
      real(dp) :: e(1:nh)
      real(dp) :: ratio
      integer :: i

      do i = 1, nh
         hh = 1.0e-2_dp/real(2**(i - 1), dp)
         call reference_vector(real(rj, qp) + real(hh, qp), cp)
         call reference_vector(real(rj, qp) - real(hh, qp), cm)
         e(i) = abs(real((cp(alpha_in) - cm(alpha_in))/(2.0_qp*real(hh, qp)), dp) - dana)
      end do

      write (*, '(A,F8.4,A,ES11.3,A)', advance="no") "  rj=", rj, "  analytic d/dr=", dana, "  ratios:"
      do i = 2, nh
         ratio = e(i - 1)/max(e(i), tiny(1.0_dp))
         write (*, '(1X,F7.3)', advance="no") ratio
         if (ratio < 3.6_dp .or. ratio > 4.4_dp) all_ok = .false.
      end do
      write (*, '(A,ES10.2)') "   final |FD-analytic|=", e(nh)
   end subroutine h_scan

   subroutine keep_scan(rj, der, alpha_in)
      real(dp), intent(in) :: rj
      real(dp), intent(in) :: der
      integer, intent(in) :: alpha_in

      n_scan_have = n_scan_have + 1
      scan_rj(n_scan_have) = rj
      scan_der(n_scan_have) = der
      scan_alpha(n_scan_have) = alpha_in
   end subroutine keep_scan

   subroutine push(v, nv, x)
      real(qp), intent(inout) :: v(:)
      integer, intent(inout) :: nv
      real(qp), intent(in) :: x

      nv = nv + 1
      v(nv) = x
   end subroutine push

   subroutine sort_ascending(v, nv)
      real(qp), intent(inout) :: v(:)
      integer, intent(in) :: nv
      integer :: i
      integer :: j
      real(qp) :: t

      do i = 1, nv - 1
         do j = i + 1, nv
            if (v(j) < v(i)) then
               t = v(i)
               v(i) = v(j)
               v(j) = t
            end if
         end do
      end do
   end subroutine sort_ascending

!  Gauss-Legendre nodes and weights on [-1,1], by Newton on P_n.
   subroutine gauss_legendre(x, w)
      real(qp), intent(out) :: x(:)
      real(qp), intent(out) :: w(:)
      integer :: n
      integer :: i
      integer :: j
      integer :: it
      real(qp) :: pi
      real(qp) :: z
      real(qp) :: p0
      real(qp) :: p1
      real(qp) :: p2
      real(qp) :: dp1

      n = size(x)
      pi = acos(-1.0_qp)
      do i = 1, (n + 1)/2
         z = cos(pi*(real(i, qp) - 0.25_qp)/(real(n, qp) + 0.5_qp))
         do it = 1, 100
            p1 = 1.0_qp
            p2 = 0.0_qp
            do j = 1, n
               p0 = p2
               p2 = p1
               p1 = ((2.0_qp*real(j, qp) - 1.0_qp)*z*p2 - (real(j, qp) - 1.0_qp)*p0)/real(j, qp)
            end do
            dp1 = real(n, qp)*(z*p1 - p2)/(z*z - 1.0_qp)
            z = z - p1/dp1
            if (abs(p1/dp1) < 1.0e-32_qp) exit
         end do
         x(i) = -z
         x(n + 1 - i) = z
         w(i) = 2.0_qp/((1.0_qp - z*z)*dp1*dp1)
         w(n + 1 - i) = w(i)
      end do
   end subroutine gauss_legendre

end program poly3operatorverify
