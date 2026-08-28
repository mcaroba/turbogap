! Validation of the envelope-targeted resonator bank, ir_bias_mode = aux.
!
! Every check here compares against something that is NOT this code: an
! analytic solution of the driven damped oscillator, an analytic dipole
! Jacobian, the linearised control loop's own period, or finite differences.
!
! Checks 1-3 are the ones that matter most, because they are the ones that
! separate a working scheme from the one this module replaced:
!
!   1  the controller moves R to R_target at all, over decades of target
!   2  the loop does not RING -- and does ring, with period 2 pi tau, when the
!      proportional gain is removed, which pins the loop's transfer function
!      rather than only its fixed point
!   3  the controller does not move the resonant FREQUENCY, which a bias
!      potential in the amplitude would, and which would corrupt the observable
!
! Run from tests/ir_aux/run.sh.

program iraux_verify

   use kinds, only: dp
   use mad_ir, only: mad_ir_type, mad_ir_init, mad_ir_free, mad_ir_push, &
                     mad_ir_ready, mad_ir_spectrum, CM_PER_INV_FS
   use ir_auxiliary_dynamics, only: ir_aux_type, ir_aux_init, ir_aux_free, &
                                    ir_aux_advance, ir_aux_evaluate, ir_aux_forces, &
                                    ir_aux_calibrate, ir_aux_save, ir_aux_load, &
                                    ir_aux_stability, ir_aux_escale_max

   implicit none

   integer :: n_fail

   n_fail = 0

   write (*, *) "=================================================================="
   write (*, *) " ir_aux: the envelope-targeted resonator bank"
   write (*, *) "=================================================================="

   call check_units(n_fail)
   call check_free_resonator(n_fail)
   call check_tracking(n_fail)
   call check_ringing(n_fail)
   call check_frequency(n_fail)
   call check_linear_response(n_fail)
   call check_forces(n_fail)
   call check_windup(n_fail)
   call check_calibration(n_fail)
   call check_stability(n_fail)
   call check_restart(n_fail)

   write (*, *) "------------------------------------------------------------------"
   if (n_fail == 0) then
      write (*, *) "ir_aux: ALL CHECKS PASSED"
   else
      write (*, '(A,I0,A)') " ir_aux: ", n_fail, " CHECK(S) FAILED"
      stop 1
   end if

contains

!**************************************************************************
!  A bank of one mode, built directly rather than calibrated, so the
!  controller and the propagator can be exercised in isolation.
   subroutine one_mode(this, nu_cm, mass_amu, tau, kappa, g, r_target, eta0, eta_frac, dt)
      type(ir_aux_type), intent(inout) :: this
      real(dp), intent(in) :: nu_cm
      real(dp), intent(in) :: mass_amu
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: kappa
      real(dp), intent(in) :: g
      real(dp), intent(in) :: r_target
      real(dp), intent(in) :: eta0
      real(dp), intent(in) :: eta_frac
      real(dp), intent(in) :: dt
      real(dp) :: two_pi
      real(dp) :: amu

      two_pi = 2.0_dp*dacos(-1.0_dp)
      amu = 103.6426965268_dp

      call ir_aux_free(this)
      this%n_modes = 1
      this%n_sites = 1
      this%dt = dt
      this%nu_power = 2.0_dp
      this%eta_max_frac = eta_frac
      this%normalisation = 1.0_dp
      allocate (this%kmap(1), this%nu(1), this%omega(1), this%eff_mass(1))
      allocate (this%gamma_k(1), this%tau(1), this%kappa(1), this%g_k(1))
      allocate (this%X(3, 1), this%P(3, 1), this%eta(1), this%R(1), this%R_target(1))
      allocate (this%I_exp(1), this%wgt(1), this%I_calc(1))
      allocate (this%muted(1))
      this%muted(1) = .false.
      this%n_live = 1
      this%kmap(1) = 1
      this%nu(1) = nu_cm
      this%omega(1) = two_pi*nu_cm/CM_PER_INV_FS
      this%eff_mass(1) = mass_amu*amu
      this%gamma_k(1) = eta0
      this%tau(1) = tau
      this%kappa(1) = kappa
      this%g_k(1) = g
      this%R_target(1) = r_target
      this%I_exp(1) = 1.0_dp
      this%wgt(1) = 1.0_dp
      this%I_calc(1) = 0.0_dp
      this%eta(1) = eta0
      this%X(1:3, 1) = 0.0_dp
      this%P(1:3, 1) = 0.0_dp
      this%active = .true.
      this%calibrated = .true.
   end subroutine one_mode

!**************************************************************************
!  Seed the single mode at amplitude r with a chosen phase.
   subroutine seed(this, r, phase)
      type(ir_aux_type), intent(inout) :: this
      real(dp), intent(in) :: r
      real(dp), intent(in) :: phase
      this%X(1, 1) = r*dcos(phase)
      this%X(2, 1) = 0.0_dp
      this%X(3, 1) = 0.0_dp
      this%P(1, 1) = this%eff_mass(1)*this%omega(1)*r*dsin(phase)
      this%P(2, 1) = 0.0_dp
      this%P(3, 1) = 0.0_dp
      this%R(1) = r
   end subroutine seed

   real(dp) function amp(this)
      type(ir_aux_type), intent(in) :: this
      real(dp) :: inv
      inv = 1.0_dp/(this%eff_mass(1)*this%omega(1))**2
      amp = dsqrt(this%X(1, 1)**2 + this%X(2, 1)**2 + this%X(3, 1)**2 &
                  + (this%P(1, 1)**2 + this%P(2, 1)**2 + this%P(3, 1)**2)*inv)
   end function amp

   subroutine verdict(name, ok, n_fail)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      integer, intent(inout) :: n_fail
      if (ok) then
         write (*, '(A,A)') "    PASS  ", name
      else
         write (*, '(A,A)') "    FAIL  ", name
         n_fail = n_fail + 1
      end if
   end subroutine verdict

!**************************************************************************
!  1. UNITS. A resonator seeded from a 1000 cm^-1 grid entry must oscillate
!  with period CM_PER_INV_FS/1000 = 33.356 fs. Treating the wavenumber as a
!  rate directly is off by 2 pi/3.34e4; multiplying by CM_PER_INV_FS instead
!  of dividing is off by ten orders of magnitude. Both were once in the file
!  this replaced, and neither shows up in any check of the controller.
   subroutine check_units(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: expect
      real(dp) :: got
      logical :: ok

      write (*, *) ""
      write (*, *) " 1. UNITS: cm^-1 -> rad/fs"

      call one_mode(b, 1000.0_dp, 100.0_dp, 1.0e12_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.1_dp)
      expect = CM_PER_INV_FS/1000.0_dp
      got = 2.0_dp*dacos(-1.0_dp)/b%omega(1)
      write (*, '(A,F14.6,A,F14.6,A)') "      period ", got, " fs, expected ", expect, " fs"
      ok = dabs(got - expect) < 1.0e-6_dp*expect
      call verdict("1000 cm^-1 is a 33.356 fs period", ok, n_fail)
      call ir_aux_free(b)
   end subroutine check_units

!**************************************************************************
!  2. THE FREE RESONATOR. With no drive, no damping and no controller, the
!  exact propagator must conserve R to round-off over many periods. This is
!  the propagator's own baseline: every later check reads R, so a propagator
!  that leaked amplitude would make all of them meaningless.
   subroutine check_free_resonator(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: r0
      real(dp) :: r1
      integer :: i
      logical :: ok

      write (*, *) ""
      write (*, *) " 2. FREE RESONATOR: R conserved with no drive, damping or control"

!     R_target = R so the controller is idle; kappa = 0 and tau huge as well.
      call one_mode(b, 1000.0_dp, 100.0_dp, 1.0e12_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      call seed(b, 1.0_dp, 0.0_dp)
      r0 = amp(b)
      do i = 1, 200000
         call ir_aux_advance(b, [0.0_dp, 0.0_dp, 0.0_dp])
      end do
      r1 = amp(b)
      write (*, '(A,ES14.6,A,ES14.6,A,ES10.2)') "      R0 = ", r0, "  R1 = ", r1, &
         "  drift = ", dabs(r1 - r0)/r0
      ok = dabs(r1 - r0)/r0 < 1.0e-10_dp
      call verdict("R conserved to 1e-10 over 200000 steps", ok, n_fail)
      call ir_aux_free(b)
   end subroutine check_free_resonator

!**************************************************************************
!  3. TRACKING. Driven on resonance, the controller must bring R to R_target
!  over decades of target. This is the check whose absence let the first
!  version of the resonator bias in mad_ir_xl.f90 push the wrong way, and the
!  one that separates this scheme from the bias potential it replaced -- a
!  potential in R returns the SAME amplitude for every target.
   subroutine check_tracking(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: targets(4)
      real(dp) :: got(4)
      real(dp) :: err
      real(dp) :: worst
      integer :: j
      logical :: ok

      write (*, *) ""
      write (*, *) " 3. TRACKING: R -> R_target, driven on resonance"

      targets = [0.5_dp, 2.0_dp, 8.0_dp, 32.0_dp]
      worst = 0.0_dp
      write (*, *) "      R_target        R_final     rel. error"
      do j = 1, 4
         call one_mode(b, 1000.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 5.0e-3_dp, &
                       targets(j), 0.0_dp, 0.5_dp, 0.5_dp)
         call seed(b, 1.0_dp, 0.0_dp)
         call drive_on_resonance(b, 200000)
         got(j) = amp(b)
         err = dabs(got(j) - targets(j))/targets(j)
         if (err > worst) worst = err
         write (*, '(A,F12.4,F15.5,ES15.3)') "  ", targets(j), got(j), err
         call ir_aux_free(b)
      end do
      write (*, '(A,ES10.2)') "      worst relative error = ", worst
      ok = worst < 0.02_dp
      call verdict("tracks every target to better than 2%", ok, n_fail)
   end subroutine check_tracking

!**************************************************************************
   subroutine drive_on_resonance(b, nstep)
      type(ir_aux_type), intent(inout) :: b
      integer, intent(in) :: nstep
      real(dp) :: t
      real(dp) :: m(3)
      integer :: i
      do i = 1, nstep
         t = dfloat(i - 1)*b%dt
         m(1) = dcos(b%omega(1)*t)
         m(2) = 0.0_dp
         m(3) = 0.0_dp
         call ir_aux_advance(b, m)
      end do
   end subroutine drive_on_resonance

!**************************************************************************
!  4. THE CONTROL LOOP. The linearised envelope obeys
!  nu'' + kappa nu' + nu/tau^2 = 0, so with kappa = 0 it is marginally stable
!  and rings with period exactly 2 pi tau, and at kappa = 2/tau it is
!  critically damped. Both halves are asserted: that the ringing IS there
!  without the proportional term pins the loop's transfer function, not just
!  its fixed point, and it is what would catch a future change that silently
!  detunes the controller.
   subroutine check_ringing(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: tau
      real(dp) :: band
      real(dp) :: period
      real(dp) :: expect
      real(dp) :: gains(3)
      integer :: j
      logical :: ok

      write (*, *) ""
      write (*, *) " 4. CONTROL LOOP: rings at 2 pi tau without the proportional term,"
      write (*, *) "                  and does not ring with it"

      tau = 500.0_dp

!     (a) kappa = 0: a pure integral controller. Must ring, at 2 pi tau.
      call one_mode(b, 1000.0_dp, 100.0_dp, tau, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 5.0_dp, 0.5_dp)
      call seed(b, 4.8_dp, 0.0_dp)
      call envelope_stats(b, 60000, band, period)
      expect = 2.0_dp*dacos(-1.0_dp)*tau
      write (*, '(A,F12.1,A,F12.1,A)') "      kappa = 0: ring period ", period, &
         " fs, predicted 2*pi*tau = ", expect, " fs"
      write (*, '(A,F10.4)') "                 late-time band / R_target = ", band
      ok = (band > 0.10_dp) .and. (dabs(period - expect) < 0.10_dp*expect)
      call verdict("pure integral control rings with period 2*pi*tau", ok, n_fail)
      call ir_aux_free(b)

!     (b) kappa*tau/2 = 0.5, 1, 2: must not ring.
      gains = [0.5_dp, 1.0_dp, 2.0_dp]
      do j = 1, 3
         call one_mode(b, 1000.0_dp, 100.0_dp, tau, gains(j)*2.0_dp/tau, 0.0_dp, &
                       4.0_dp, 0.0_dp, 5.0_dp, 0.5_dp)
         call seed(b, 4.8_dp, 0.0_dp)
         call envelope_stats(b, 60000, band, period)
         write (*, '(A,F6.2,A,F10.5)') "      kappa*tau/2 = ", gains(j), &
            ":  late-time band / R_target = ", band
         ok = band < 0.02_dp
         call verdict("no ringing at this gain", ok, n_fail)
         call ir_aux_free(b)
      end do
   end subroutine check_ringing

!**************************************************************************
!  Late-time spread of the envelope about its target, and the period of its
!  oscillation, measured from crossings of R = R_target over the second half
!  of the run.
   subroutine envelope_stats(b, nstep, band, period)
      type(ir_aux_type), intent(inout) :: b
      integer, intent(in) :: nstep
      real(dp), intent(out) :: band
      real(dp), intent(out) :: period
      real(dp) :: r
      real(dp) :: r_lo
      real(dp) :: r_hi
      real(dp) :: prev
      real(dp) :: t_first
      real(dp) :: t_last
      real(dp) :: t
      integer :: i
      integer :: n_cross

      r_lo = 1.0e30_dp
      r_hi = -1.0e30_dp
      n_cross = 0
      t_first = -1.0_dp
      t_last = -1.0_dp
      prev = amp(b) - b%R_target(1)

      do i = 1, nstep
         call ir_aux_advance(b, [0.0_dp, 0.0_dp, 0.0_dp])
         t = dfloat(i)*b%dt
         r = amp(b)
!        Crossings over the whole run, so a heavily damped case still gives a
!        period if it has one; the band only over the second half, so the
!        initial approach is not counted as ringing.
         if (prev < 0.0_dp .and. (r - b%R_target(1)) >= 0.0_dp) then
            n_cross = n_cross + 1
            if (t_first < 0.0_dp) t_first = t
            t_last = t
         end if
         prev = r - b%R_target(1)
         if (i > nstep/2) then
            if (r < r_lo) r_lo = r
            if (r > r_hi) r_hi = r
         end if
      end do

      band = (r_hi - r_lo)/b%R_target(1)
      if (n_cross > 1) then
         period = (t_last - t_first)/dfloat(n_cross - 1)
      else
         period = 0.0_dp
      end if
   end subroutine envelope_stats

!**************************************************************************
!  5. FREQUENCY PRESERVATION. The controller must not move the resonant
!  frequency: a channel that has been detuned is no longer reporting the
!  intensity at its own w_k, and the linear-response calibration then refers
!  to the wrong frequency. This is exactly what a bias POTENTIAL in R does and
!  the radial controller does not -- applying one scalar to both X and P is a
!  radial flow, and in the complex quadrature z = X + iP/(mu w) it multiplies
!  |z| while leaving arg(z) alone.
   subroutine check_frequency(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: p_free
      real(dp) :: p_ctrl
      real(dp) :: expect
      logical :: ok

      write (*, *) ""
      write (*, *) " 5. FREQUENCY: the controller rescales R without detuning w_k"

      expect = CM_PER_INV_FS/1000.0_dp

!     Free.
      call one_mode(b, 1000.0_dp, 100.0_dp, 1.0e12_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      call seed(b, 1.0_dp, 0.0_dp)
      p_free = zero_cross_period(b, 200000)
      call ir_aux_free(b)

!     With the controller working hard: R starts a factor of 8 from target, so
!     the radial factor is far from 1 throughout.
      call one_mode(b, 1000.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 0.0_dp, &
                    8.0_dp, 0.0_dp, 5.0_dp, 0.5_dp)
      call seed(b, 1.0_dp, 0.0_dp)
      p_ctrl = zero_cross_period(b, 200000)
      call ir_aux_free(b)

      write (*, '(A,F14.6,A)') "      free period       ", p_free, " fs"
      write (*, '(A,F14.6,A)') "      under control     ", p_ctrl, " fs"
      write (*, '(A,F14.6,A)') "      analytic          ", expect, " fs"
      write (*, '(A,ES10.2)') "      relative shift    ", dabs(p_ctrl - expect)/expect
      ok = dabs(p_ctrl - expect)/expect < 1.0e-4_dp
      call verdict("no frequency shift beyond 1e-4 under active control", ok, n_fail)
   end subroutine check_frequency

!**************************************************************************
   real(dp) function zero_cross_period(b, nstep)
      type(ir_aux_type), intent(inout) :: b
      integer, intent(in) :: nstep
      real(dp) :: prev
      real(dp) :: t_first
      real(dp) :: t_last
      real(dp) :: frac
      real(dp) :: cur
      integer :: i
      integer :: n

      n = 0
      t_first = -1.0_dp
      t_last = -1.0_dp
      prev = b%X(1, 1)
      do i = 1, nstep
         call ir_aux_advance(b, [0.0_dp, 0.0_dp, 0.0_dp])
         cur = b%X(1, 1)
         if (prev < 0.0_dp .and. cur >= 0.0_dp) then
!           Linear interpolation of the crossing, so the answer is not limited
!           to the step size: with dt = 0.5 fs and a 33 fs period, the
!           uninterpolated period quantises at 1.5%, far above the 1e-4 the
!           check is trying to resolve.
            frac = -prev/(cur - prev)
            n = n + 1
            if (t_first < 0.0_dp) t_first = (dfloat(i - 1) + frac)*b%dt
            t_last = (dfloat(i - 1) + frac)*b%dt
         end if
         prev = cur
      end do
      if (n > 1) then
         zero_cross_period = (t_last - t_first)/dfloat(n - 1)
      else
         zero_cross_period = 0.0_dp
      end if
   end function zero_cross_period

!**************************************************************************
!  6. LINEAR RESPONSE. With the controller off and eta held at Gamma, the bank
!  is a plain damped driven oscillator, whose steady-state amplitude under
!  M(t) = A cos(w t) is ANALYTIC:
!
!     R = g A / (mu Gamma w)
!
!  This is the normalisation the whole g_k calibration rests on -- (10) is
!  derived from exactly this response -- so it is checked against the closed
!  form rather than against another run of the same code.
   subroutine check_linear_response(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: gam
      real(dp) :: g
      real(dp) :: expect
      real(dp) :: got
      real(dp) :: err
      logical :: ok

      write (*, *) ""
      write (*, *) " 6. LINEAR RESPONSE: steady state against R = g A/(mu Gamma w)"

      gam = 2.0_dp*dacos(-1.0_dp)*10.0_dp/CM_PER_INV_FS
      g = 1.0e-2_dp

!     kappa = 0 and tau enormous, so neither half of the controller acts and
!     eta stays at Gamma.
      call one_mode(b, 1000.0_dp, 100.0_dp, 1.0e12_dp, 0.0_dp, g, 1.0_dp, gam, 0.5_dp, 0.2_dp)
      call seed(b, 0.0_dp, 0.0_dp)
      call drive_on_resonance(b, 400000)
      got = amp(b)
      expect = g/(b%eff_mass(1)*gam*b%omega(1))
      err = dabs(got - expect)/expect
      write (*, '(A,ES14.6)') "      integrated steady-state R = ", got
      write (*, '(A,ES14.6)') "      analytic  g A/(mu Gam w)  = ", expect
      write (*, '(A,ES10.2)') "      relative error            = ", err
      ok = err < 1.0e-3_dp
      call verdict("steady state matches the closed-form response", ok, n_fail)
      call ir_aux_free(b)
   end subroutine check_linear_response

!**************************************************************************
!  7. THE BACK-REACTION FORCE, by h-scan against finite differences of the
!  coupling energy for an ANALYTIC dipole model whose Jacobian is written out
!  by hand. The error must fall as h^2 -- a factor of 4 per halving -- and
!  then turn round on round-off. The ratio is asserted, not a tolerance: a
!  gradient wrong by a constant factor passes any single-h check, which is how
!  mad_ir_xl's first version got a factor of 1.5 wrong and kept it.
   subroutine check_forces(n_fail)
      integer, intent(inout) :: n_fail
      integer, parameter :: na = 4
      type(ir_aux_type) :: b
      real(dp) :: r(3, na)
      real(dp) :: j_an(3, 3, na)
      real(dp) :: f(3, na)
      real(dp) :: h
      real(dp) :: e_p
      real(dp) :: e_m
      real(dp) :: f_fd
      real(dp) :: err
      real(dp) :: err_prev
      real(dp) :: ratio
      real(dp) :: escale
      real(dp) :: worst
      integer :: i
      integer :: ia
      integer :: ib
      integer :: is
      integer :: n_good
      logical :: ok

      write (*, *) ""
      write (*, *) " 7. BACK-REACTION FORCE: h-scan against finite differences"

      escale = 0.7_dp
      call one_mode(b, 1000.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 3.0e-2_dp, &
                    2.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      b%n_sites = na
      call seed(b, 1.7_dp, 0.9_dp)
!     Give the three components genuinely different values, so a transposed
!     contraction cannot pass.
      b%X(2, 1) = 0.4_dp*b%X(1, 1)
      b%X(3, 1) = -0.9_dp*b%X(1, 1)

      r(1:3, 1) = [0.31_dp, -0.72_dp, 0.55_dp]
      r(1:3, 2) = [1.13_dp, 0.24_dp, -0.41_dp]
      r(1:3, 3) = [-0.87_dp, 0.66_dp, 1.02_dp]
      r(1:3, 4) = [0.05_dp, -1.31_dp, -0.28_dp]

      call dipole_jacobian(r, na, j_an)
      f = 0.0_dp
      call ir_aux_forces(b, escale, j_an, f)

      write (*, *) "            h        max|F_fd - F|      ratio"
      worst = 0.0_dp
      n_good = 0
      err_prev = -1.0_dp
      h = 1.0e-2_dp
      do is = 1, 9
         err = 0.0_dp
         do ia = 1, na
            do ib = 1, 3
               r(ib, ia) = r(ib, ia) + h
               e_p = coupling_energy(b, r, na, escale)
               r(ib, ia) = r(ib, ia) - 2.0_dp*h
               e_m = coupling_energy(b, r, na, escale)
               r(ib, ia) = r(ib, ia) + h
!              F = -dE/dr
               f_fd = -(e_p - e_m)/(2.0_dp*h)
               err = max(err, dabs(f_fd - f(ib, ia)))
            end do
         end do
         if (err_prev > 0.0_dp) then
            ratio = err_prev/max(err, 1.0e-300_dp)
            write (*, '(A,ES12.3,ES18.6,F12.3)') "  ", h, err, ratio
!           h^2 convergence, before round-off takes over.
            if (is <= 5 .and. ratio > 3.0_dp .and. ratio < 5.0_dp) n_good = n_good + 1
         else
            write (*, '(A,ES12.3,ES18.6)') "  ", h, err
         end if
         err_prev = err
         h = 0.5_dp*h
      end do

      ok = n_good >= 3
      call verdict("error falls as h^2 over at least 3 halvings", ok, n_fail)
      call ir_aux_free(b)
   end subroutine check_forces

!**************************************************************************
!  An analytic per-atom dipole with a non-symmetric, position-dependent
!  Jacobian, so the contraction in ir_aux_forces is tested and not merely
!  reproduced.
!     m_i = ( sin(x) + 0.3 y ,  y^2 + 0.2 z ,  x z )
   subroutine dipole_jacobian(r, na, j_an)
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: na
      real(dp), intent(out) :: j_an(3, 3, na)
      integer :: i

      j_an = 0.0_dp
      do i = 1, na
         j_an(1, 1, i) = dcos(r(1, i))
         j_an(1, 2, i) = 0.3_dp
         j_an(1, 3, i) = 0.0_dp
         j_an(2, 1, i) = 0.0_dp
         j_an(2, 2, i) = 2.0_dp*r(2, i)
         j_an(2, 3, i) = 0.2_dp
         j_an(3, 1, i) = r(3, i)
         j_an(3, 2, i) = 0.0_dp
         j_an(3, 3, i) = r(1, i)
      end do
   end subroutine dipole_jacobian

!**************************************************************************
!  E = -escale sum_k g_k X_k . M,  M = sum_i m_i(r_i). The generator of the
!  force ir_aux_forces adds.
   real(dp) function coupling_energy(b, r, na, escale)
      type(ir_aux_type), intent(in) :: b
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: na
      real(dp), intent(in) :: escale
      real(dp) :: m(3)
      real(dp) :: e_eff(3)
      integer :: i
      integer :: k

      m = 0.0_dp
      do i = 1, na
         m(1) = m(1) + dsin(r(1, i)) + 0.3_dp*r(2, i)
         m(2) = m(2) + r(2, i)**2 + 0.2_dp*r(3, i)
         m(3) = m(3) + r(1, i)*r(3, i)
      end do
      e_eff = 0.0_dp
      do k = 1, b%n_modes
         e_eff(1:3) = e_eff(1:3) + b%g_k(k)*b%X(1:3, k)
      end do
      coupling_energy = -escale*(e_eff(1)*m(1) + e_eff(2)*m(2) + e_eff(3)*m(3))
   end function coupling_energy

!**************************************************************************
!  8. ANTI-WINDUP. An unreachable target makes the integral in etadot grow
!  without bound. The clamp must hold eta inside its fraction of w_k, and the
!  channel must still be resonating afterwards -- past |eta| = 2 w it would be
!  overdamped, and would report an intensity for a frequency it cannot see.
   subroutine check_windup(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      real(dp) :: cap
      real(dp) :: per
      real(dp) :: expect
      logical :: ok

      write (*, *) ""
      write (*, *) " 8. ANTI-WINDUP: an unreachable target does not run eta away"

!     No drive at all and a target of 100: the controller can never get there,
!     so the integral winds up for the whole run.
      call one_mode(b, 1000.0_dp, 100.0_dp, 50.0_dp, 0.0_dp, 0.0_dp, 100.0_dp, 0.0_dp, 0.1_dp, 0.5_dp)
      call seed(b, 1.0_dp, 0.0_dp)
      call drive_on_resonance(b, 200000)
      cap = b%eta_max_frac*b%omega(1)
      write (*, '(A,ES14.6,A,ES14.6)') "      eta = ", b%eta(1), "   clamp = ", cap
      ok = dabs(b%eta(1)) <= cap*(1.0_dp + 1.0e-12_dp)
      call verdict("eta stays inside the clamp", ok, n_fail)

      expect = CM_PER_INV_FS/1000.0_dp
      per = zero_cross_period(b, 100000)
      write (*, '(A,F12.6,A,F12.6,A)') "      period still ", per, " fs (free ", expect, " fs)"
!     A clamp of 0.1 w permits a second-order shift of sqrt(w^2-(eta/2)^2),
!     i.e. 0.125%, so the channel is still tuned to what it reports.
      ok = dabs(per - expect)/expect < 5.0e-3_dp
      call verdict("the channel is still resonating near w_k", ok, n_fail)
      call ir_aux_free(b)
   end subroutine check_windup

!**************************************************************************
!  9. CALIBRATION, through the path the driver actually uses: a mad_ir
!  ensemble filled with a known dipole, then ir_aux_init and ir_aux_calibrate.
!  Asserts the identity (10) rearranges to <R^2> = R_target^2 -- which is what
!  "calibrated" means -- and that no channel came out NaN, which is the failure
!  mode a zero in S_MM produces and which would poison every force through the
!  shared E_eff.
   subroutine check_calibration(n_fail)
      integer, intent(inout) :: n_fail
      integer, parameter :: nf = 40
      integer, parameter :: nlag = 256
      type(mad_ir_type) :: par
      type(ir_aux_type) :: b
      real(dp) :: nu(nf)
      real(dp) :: iexp(nf)
      real(dp) :: wgt(nf)
      real(dp) :: mu(3)
      real(dp) :: dt
      real(dp) :: w0
      real(dp) :: t
      real(dp) :: pred
      real(dp) :: worst
      real(dp) :: s_mm
      character(len=512) :: msg
      integer :: k
      integer :: i
      integer :: n_bad
      logical :: ok

      write (*, *) ""
      write (*, *) " 9. CALIBRATION: the driver's path, mad_ir -> ir_aux_calibrate"

      dt = 2.0_dp
      do k = 1, nf
         nu(k) = 200.0_dp*dfloat(k)
         iexp(k) = 1.0_dp + 0.5_dp*dsin(0.3_dp*dfloat(k))
         wgt(k) = 1.0_dp
      end do

      call mad_ir_init(par, dt, nlag, 2*nlag, nu, iexp, wgt, .true., 2.0_dp, &
                       "hann", .true., .true., .true., .false.)

!     A broadband dipole, so every channel has some power and none is muted
!     for want of it.
      w0 = 2.0_dp*dacos(-1.0_dp)*1400.0_dp/CM_PER_INV_FS
      do i = 1, 3*nlag
         t = dfloat(i - 1)*dt
         mu(1) = dcos(w0*t) + 0.5_dp*dcos(2.3_dp*w0*t) + 0.3_dp*dcos(0.41_dp*w0*t)
         mu(2) = 0.7_dp*dsin(1.7_dp*w0*t) + 0.2_dp*dcos(0.7_dp*w0*t)
         mu(3) = 0.4_dp*dcos(3.1_dp*w0*t + 0.6_dp) + 0.6_dp*dsin(0.23_dp*w0*t)
         call mad_ir_push(par, mu)
      end do
      write (*, '(A,L2)') "      mad_ir_ready = ", mad_ir_ready(par)

      call ir_aux_init(b, par, 8, dt, 100.0_dp, 10.0_dp, 500.0_dp, -1.0_dp, 0.1_dp, ok, msg)
      write (*, *) "      ", trim(msg)
      call verdict("ir_aux_init succeeded", ok, n_fail)

      call ir_aux_calibrate(b, par, 300.0_dp, 1.0e-3_dp, ok, msg)
      write (*, *) "      ", trim(msg)
      call verdict("ir_aux_calibrate succeeded", ok, n_fail)

!     (10): g_k = mu w sqrt(Gam R_t^2/S_MM)  =>  g^2 S_MM/(mu^2 w^2 Gam) = R_t^2.
      worst = 0.0_dp
      n_bad = 0
      do k = 1, b%n_modes
         if (b%g_k(k) == 0.0_dp) cycle
         s_mm = par%I_calc(b%kmap(k))/b%nu(k)**b%nu_power
         pred = b%g_k(k)**2*s_mm/(b%eff_mass(k)**2*b%omega(k)**2*b%gamma_k(k))
         if (b%R_target(k) > 0.0_dp) then
            worst = max(worst, dabs(dsqrt(dabs(pred)) - b%R_target(k))/b%R_target(k))
         end if
         if (b%g_k(k) /= b%g_k(k)) n_bad = n_bad + 1
         if (b%R_target(k) /= b%R_target(k)) n_bad = n_bad + 1
      end do
      write (*, '(A,ES10.2)') "      worst |sqrt(<R^2>_LR) - R_target|/R_target = ", worst
      write (*, '(A,I0)') "      NaN couplings or targets: ", n_bad
      call verdict("the calibration identity holds to 1e-10", worst < 1.0e-10_dp, n_fail)
      call verdict("no NaN in any coupling or target", n_bad == 0, n_fail)
!     R_MAX = 1 is the scale the strongest band is put at -- the convention
!     mad_ir_notes.md section 4 intends. It is asserted rather than merely
!     printed because the back-reaction goes as R_MAX SQUARED, so a change here
!     silently rescales every force on the atoms by its square.
      write (*, '(A,F10.5)') "      max R_target = ", maxval(b%R_target)
      call verdict("the strongest band sits at R_MAX = 1", &
                   dabs(maxval(b%R_target) - 1.0_dp) < 1.0e-10_dp, n_fail)

      call ir_aux_free(b)
      call mad_ir_free(par)
   end subroutine check_calibration

!**************************************************************************
!  9b. THE STABILITY BOUND. The coupling -g_k X_k.s is bilinear, hence unbounded
!  below, and is held only by the resonator spring and the signal's own
!  stiffness. The bound (S3) is checked three ways: against an INDEPENDENT
!  evaluation of the formula written out here from the primitive quantities,
!  for the linearity in energy_scale that (S3) demands, and for the exact
!  inverse relation (S4) -- escale_max must be the scale at which Lambda is 1.
!
!  This is the check that would have saved a water box: the calibration matches
!  amplitudes and never asks whether the coupling it produces is below
!  threshold, so nothing else in this suite can catch it.
   subroutine check_stability(n_fail)
      integer, intent(inout) :: n_fail
      integer, parameter :: nf = 40
      integer, parameter :: nlag = 256
      type(mad_ir_type) :: par
      type(ir_aux_type) :: b
      real(dp), parameter :: KB = 8.6173303e-5_dp
      real(dp) :: nu(nf)
      real(dp) :: iexp(nf)
      real(dp) :: wgt(nf)
      real(dp) :: mu(3)
      real(dp) :: dt
      real(dp) :: w0
      real(dp) :: t
      real(dp) :: temp
      real(dp) :: escale
      real(dp) :: sum_g2
      real(dp) :: lam_ref
      real(dp) :: lam_1
      real(dp) :: lam_2
      real(dp) :: lam_at_max
      character(len=512) :: msg
      integer :: k
      integer :: i
      logical :: ok

      write (*, *) ""
      write (*, *) " 9b. STABILITY BOUND: Lambda = gamma sum(g^2/mu w^2) C(0)/(3 kB T)"

      dt = 2.0_dp
      temp = 300.0_dp
      escale = 1.0e-3_dp
      do k = 1, nf
         nu(k) = 200.0_dp*dfloat(k)
         iexp(k) = 1.0_dp + 0.5_dp*dsin(0.3_dp*dfloat(k))
         wgt(k) = 1.0_dp
      end do
      call mad_ir_init(par, dt, nlag, 2*nlag, nu, iexp, wgt, .true., 2.0_dp, &
                       "hann", .true., .true., .true., .false.)
      w0 = 2.0_dp*dacos(-1.0_dp)*1400.0_dp/CM_PER_INV_FS
      do i = 1, 3*nlag
         t = dfloat(i - 1)*dt
         mu(1) = dcos(w0*t) + 0.5_dp*dcos(2.3_dp*w0*t) + 0.3_dp*dcos(0.41_dp*w0*t)
         mu(2) = 0.7_dp*dsin(1.7_dp*w0*t) + 0.2_dp*dcos(0.7_dp*w0*t)
         mu(3) = 0.4_dp*dcos(3.1_dp*w0*t + 0.6_dp) + 0.6_dp*dsin(0.23_dp*w0*t)
         call mad_ir_push(par, mu)
      end do

      call ir_aux_init(b, par, 8, dt, 100.0_dp, 10.0_dp, 500.0_dp, -1.0_dp, 0.1_dp, ok, msg)
      call ir_aux_calibrate(b, par, temp, escale, ok, msg)
      call verdict("calibration with the bound succeeded", ok, n_fail)

!     An independent evaluation of (S3) from the primitives.
      sum_g2 = 0.0_dp
      do k = 1, b%n_modes
         if (b%muted(k)) cycle
         sum_g2 = sum_g2 + b%g_k(k)**2/(b%eff_mass(k)*b%omega(k)**2)
      end do
      lam_ref = escale/dfloat(b%n_live)*sum_g2*par%acf(0)/(3.0_dp*KB*temp)
      write (*, '(A,ES14.6)') "      Lambda (module)      = ", b%stab
      write (*, '(A,ES14.6)') "      Lambda (independent) = ", lam_ref
      call verdict("Lambda matches an independent evaluation", &
                   dabs(b%stab - lam_ref) <= 1.0e-12_dp*max(1.0_dp, dabs(lam_ref)), n_fail)

!     (S3) is linear in the back-reaction scale.
      lam_1 = ir_aux_stability(b, escale)
      lam_2 = ir_aux_stability(b, 2.0_dp*escale)
      write (*, '(A,F12.6)') "      Lambda(2s)/Lambda(s) = ", lam_2/lam_1
      call verdict("Lambda is linear in exp_energy_scales", &
                   dabs(lam_2/lam_1 - 2.0_dp) < 1.0e-12_dp, n_fail)

!     (S4) inverts (S3) exactly: at escale_max the number is 1.
      lam_at_max = ir_aux_stability(b, ir_aux_escale_max(b))
      write (*, '(A,ES14.6)') "      escale_max           = ", ir_aux_escale_max(b)
      write (*, '(A,F14.10)') "      Lambda(escale_max)   = ", lam_at_max
      call verdict("escale_max is exactly where Lambda = 1", &
                   dabs(lam_at_max - 1.0_dp) < 1.0e-10_dp, n_fail)

      call ir_aux_free(b)
      call mad_ir_free(par)
   end subroutine check_stability

!**************************************************************************
!  10. RESTART. A round trip must return X, P and eta bit for bit -- eta
!  especially, because it is an integrator and dropping it throws away
!  everything the controller had learned. And a file whose frequency grid
!  disagrees must be REFUSED rather than adopted.
   subroutine check_restart(n_fail)
      integer, intent(inout) :: n_fail
      type(ir_aux_type) :: b
      type(ir_aux_type) :: c
      character(len=512) :: msg
      real(dp) :: d
      logical :: ok

      write (*, *) ""
      write (*, *) "10. RESTART: round trip, and refusal of a mismatched grid"

      call one_mode(b, 1000.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 3.0e-2_dp, &
                    2.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      call seed(b, 1.3_dp, 0.4_dp)
      call drive_on_resonance(b, 5000)

      call ir_aux_save(b, "ir_aux_roundtrip.dat", ok, msg)
      write (*, *) "      ", trim(msg)
      call verdict("save succeeded", ok, n_fail)

      call one_mode(c, 1000.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 0.0_dp, &
                    1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      call ir_aux_load(c, "ir_aux_roundtrip.dat", ok, msg)
      write (*, *) "      ", trim(msg)
      call verdict("load succeeded", ok, n_fail)

      d = maxval(dabs(c%X - b%X)) + maxval(dabs(c%P - b%P)) &
          + dabs(c%eta(1) - b%eta(1)) + dabs(c%g_k(1) - b%g_k(1)) &
          + dabs(c%R_target(1) - b%R_target(1))
      write (*, '(A,ES10.2)') "      total round-trip difference = ", d
      call verdict("X, P, eta, g_k and R_target survive the round trip", d == 0.0_dp, n_fail)
      call ir_aux_free(c)

!     A grid that has moved means the modes mean something else.
      call one_mode(c, 1200.0_dp, 100.0_dp, 500.0_dp, 2.0_dp/500.0_dp, 0.0_dp, &
                    1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp)
      call ir_aux_load(c, "ir_aux_roundtrip.dat", ok, msg)
      write (*, *) "      ", trim(msg)
      call verdict("a mismatched frequency grid is refused",.not. ok, n_fail)
      call ir_aux_free(c)

      call ir_aux_free(b)
   end subroutine check_restart

end program iraux_verify
