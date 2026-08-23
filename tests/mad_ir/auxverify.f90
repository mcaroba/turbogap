! The MAD IR bias in AUXILIARY-VARIABLE form: ir_acf_mode = exponential.
!
! In this mode the running dipole autocorrelation is not recomputed from the
! stored buffer each step. It is carried as a set of auxiliary degrees of
! freedom, one per lag, integrated alongside the atoms:
!
!   s_n'(t) = -(1/tau_mem) s_n(t) + (1/tau_mem) [ d(t) . d(t - tau_n) ]
!
! which makes the bias a functional of an exponentially decaying window of the
! trajectory rather than of a hard one -- the Markovian embedding of a
! generalized Langevin bias, with the target spectrum in the role of the bath.
!
! What that changes, and therefore what has to be tested, is the GRADIENT. The
! newest dipole reaches the loss through exactly one filter update rather than
! through every pair in the buffer, so none of the block estimator's gradient
! applies and the whole path is new code. In increasing order of how much they
! can hide:
!
!   1. The filter is the ODE it claims to be. Fed a constant correlation, s_n
!      must follow c(1 - (1-alpha)^k) with alpha = 1 - exp(-dt/tau_mem) -- the
!      exact discretisation, not the Euler step dt/tau_mem, which is a
!      different filter that additionally goes unstable for dt > tau_mem.
!
!   2. Every lag charges up together. The claim that the fill transient is a
!      pure overall factor -- which is what lets match_scale absorb it instead
!      of a bias correction -- is true only if every s_n has had the same number
!      of updates. Fed a signal whose correlation is the same at every lag, all
!      of them must be EQUAL, not merely similar.
!
!   3. THE GRADIENT, against central differences of the loss in the newest
!      dipole. This is the quantity the force is built from and the only one
!      that has to be right to the last digit. Run as an h-scan, over every
!      combination of the mean subtraction and the fitted scale, because both
!      enter the gradient and both are things that are wrong silently.
!
!      The h-scan is also the only test that catches the factor of two at zero
!      lag. s_0 correlates the newest dipole with itself and so depends on it
!      twice over; a gradient written as one line for every lag misses that,
!      and since win(0)C(0) is the largest single term in the transform the
!      resulting force is a few per cent wrong everywhere -- large enough to
!      matter and far too small to notice in a spectrum.
!
!   4. Restart. The auxiliary variables are state; a file written under a
!      different tau_mem or a different estimator describes a different
!      observable and must be refused rather than adopted.
!
!   5. That the two modes are estimating the same thing: on stationary data
!      they must put the peak in the same place, or the shared transform is not
!      in fact shared.
!
program auxverify

   use kinds
   use mad_ir

   implicit none

   integer :: n_fail
   real(dp) :: escale

   n_fail = 0
!  exp_energy_scales for this observable, as in irverify.
   escale = 3.7d0

   call test_filter_is_the_ode()
   call test_lags_charge_together()
   call test_gradient()
   call test_restart()
   call test_modes_agree()
   call test_readiness()

   write (*, *)
   if (n_fail == 0) then
      write (*, '(A)') "==> auxverify: all checks passed"
   else
      write (*, '(A,I0,A)') "==> auxverify: ", n_fail, " CHECK(S) FAILED"
      stop 1
   end if

contains

   subroutine check(name, got, want, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got, want, tol
      real(dp) :: err
      err = dabs(got - want)
      if (dabs(want) > 1.d0) err = err/dabs(want)
      if (err <= tol .and. .not. (got /= got)) then
         write (*, '(A,A,A,ES13.6,A,ES13.6,A,ES9.2,A)') "   ok   ", name, &
            "  got ", got, "  want ", want, "  (err ", err, ")"
      else
         write (*, '(A,A,A,ES13.6,A,ES13.6,A,ES9.2,A,ES9.2,A)') "   FAIL ", name, &
            "  got ", got, "  want ", want, "  (err ", err, " > ", tol, ")"
         n_fail = n_fail + 1
      end if
   end subroutine check

   subroutine check_true(name, got)
      character(len=*), intent(in) :: name
      logical, intent(in) :: got
      if (got) then
         write (*, '(A,A)') "   ok   ", name
      else
         write (*, '(A,A)') "   FAIL ", name
         n_fail = n_fail + 1
      end if
   end subroutine check_true

!  A small but honest configuration: a frequency grid a real fit would use,
!  and an experiment with a band in it so the loss has a gradient worth
!  differentiating.
   subroutine build(st, dt, n_lag, n_window, tau_mem, subtract_mean, match_scale)
      type(mad_ir_type), intent(inout) :: st
      real(dp), intent(in) :: dt, tau_mem
      integer, intent(in) :: n_lag, n_window
      logical, intent(in) :: subtract_mean, match_scale
      real(dp), allocatable :: nu(:), Iexp(:), wgt(:)
      integer :: n_freq, k
      n_freq = 120
      allocate (nu(n_freq), Iexp(n_freq), wgt(n_freq))
      do k = 1, n_freq
         nu(k) = 3500.d0*dfloat(k)/dfloat(n_freq)
!        A band near 1500 cm^-1, so the residual is not a constant and the
!        gradient is not dominated by one frequency.
         Iexp(k) = 1.d3*dexp(-((nu(k) - 1500.d0)/220.d0)**2)
         wgt(k) = 1.d0
      end do
      call mad_ir_init(st, dt, n_lag, n_window, nu, Iexp, wgt, match_scale, &
                       2.d0, "hann", subtract_mean, .true., .true., .false., &
                       "exponential", tau_mem)
      deallocate (nu, Iexp, wgt)
   end subroutine build

!  A dipole trace with structure at two frequencies, an offset, and a
!  deterministic wobble. Deterministic on purpose: an h-scan has to be
!  repeatable.
   subroutine dipole_at(i, dt, mu)
      integer, intent(in) :: i
      real(dp), intent(in) :: dt
      real(dp), intent(out) :: mu(1:3)
      real(dp) :: t, two_pi, w1, w2
      two_pi = 2.d0*dacos(-1.d0)
      t = dfloat(i)*dt
      w1 = two_pi*1500.d0/33356.40952d0
      w2 = two_pi*800.d0/33356.40952d0
      mu(1) = 0.60d0*dcos(w1*t) + 0.20d0*dcos(w2*t + 0.7d0) + 0.35d0
      mu(2) = 0.45d0*dsin(w1*t + 0.3d0) + 0.15d0*dcos(w2*t) - 0.20d0
      mu(3) = 0.30d0*dcos(w1*t + 1.1d0) + 0.25d0*dsin(w2*t + 0.4d0) + 0.10d0
   end subroutine dipole_at

!  ------------------------------------------------- 1. the filter is the ODE

   subroutine test_filter_is_the_ode()
      type(mad_ir_type) :: st
      real(dp) :: dt, tau, mu(1:3), c, alpha, want
      integer :: n_lag, n_window, i, k
      write (*, '(A)') "1. the filter is the exact discretisation of the memory ODE"

      dt = 2.d0
      tau = 40.d0
      n_lag = 6
      n_window = 20
!     Mean subtraction off, so that a constant dipole gives a constant
!     correlation and the filter's own response is what is being measured.
      call build(st, dt, n_lag, n_window, tau, .false., .true.)

      alpha = 1.d0 - dexp(-dt/tau)
      call check("alpha = 1 - exp(-dt/tau_mem)    ", st%alpha, alpha, 1.d-15)
      call check_true("and is not the Euler step dt/tau", dabs(st%alpha - dt/tau) > 1.d-3)

      mu = [0.3d0, -0.5d0, 0.2d0]
      c = dot_product(mu, mu)
!     Push enough frames that the filter has started, then count updates.
      do i = 1, n_window
         call mad_ir_push(st, mu)
      end do
      k = st%n_ema
      want = c*(1.d0 - (1.d0 - alpha)**k)
      call check("s_0 after k updates             ", st%acf(0), want, 1.d-13)
      call check("s_n at the longest lag          ", st%acf(n_lag), want, 1.d-13)
      write (*, '(A,I0,A,F8.5)') "        k = ", k, "   fill fraction = ", &
         1.d0 - (1.d0 - alpha)**k

      call mad_ir_free(st)
      write (*, *)
   end subroutine test_filter_is_the_ode

!  ------------------------------------------ 2. every lag charges up together

   subroutine test_lags_charge_together()
      type(mad_ir_type) :: st
      real(dp) :: dt, tau, mu(1:3), spread
      integer :: n_lag, n_window, i, t
      write (*, '(A)') "2. every lag has had the same number of updates"

      dt = 2.d0
      tau = 60.d0
      n_lag = 8
      n_window = 24
      call build(st, dt, n_lag, n_window, tau, .false., .true.)

      mu = [0.4d0, 0.1d0, -0.3d0]
!     A constant dipole makes d(0).d(n) the same at every lag, so any spread
!     between the s_n can only come from their having been updated a different
!     number of times.
      do i = 1, n_window
         call mad_ir_push(st, mu)
      end do
      spread = 0.d0
      do t = 0, n_lag
         spread = max(spread, dabs(st%acf(t) - st%acf(0)))
      end do
      call check("spread across lags              ", spread, 0.d0, 1.d-300)
      call check_true("and they are not all zero       ", dabs(st%acf(0)) > 1.d-6)

!     Nothing may have been filtered before the longest lag existed, or the
!     lags would not be able to charge together.
      call mad_ir_free(st)
      call build(st, dt, n_lag, n_window, tau, .false., .true.)
      do i = 1, n_lag
         call mad_ir_push(st, mu)
      end do
      call check("no update before every lag exists", dfloat(st%n_ema), 0.d0, 1.d-300)
      call mad_ir_push(st, mu)
      call check("first update at n_lag+1 frames  ", dfloat(st%n_ema), 1.d0, 1.d-300)

      call mad_ir_free(st)
      write (*, *)
   end subroutine test_lags_charge_together

!  --------------------------------------------------- 3. the gradient, by h-scan

   subroutine test_gradient()
      integer :: i, j
      write (*, '(A)') "3. lambda = dL/dmu(newest) against central differences"
      write (*, '(A)') "        subtract_mean  match_scale        h        worst rel. error"
      do i = 0, 1
         do j = 0, 1
            call gradient_hscan(i == 1, j == 1)
         end do
      end do
      write (*, *)
   end subroutine test_gradient

   subroutine gradient_hscan(subtract_mean, match_scale)
      logical, intent(in) :: subtract_mean, match_scale
      type(mad_ir_type) :: st, snap
      real(dp) :: dt, tau, mu(1:3), munew(1:3), mup(1:3)
      real(dp) :: loss, lambda(1:3), lp, lm, fd, h, worst, best
      real(dp) :: lam_ref(1:3)
      integer :: n_lag, n_window, i, c, ih
      logical :: ok

      dt = 2.d0
      tau = 50.d0
      n_lag = 10
      n_window = 30

      call build(st, dt, n_lag, n_window, tau, subtract_mean, match_scale)
!     Fill past readiness, leaving the LAST frame to be pushed inside the scan
!     so that the same filter update is the one being differentiated.
      do i = 1, n_window + 12
         call dipole_at(i, dt, mu)
         call mad_ir_push(st, mu)
      end do
      call dipole_at(n_window + 13, dt, munew)

!     Snapshot before the final push. Intrinsic assignment of a derived type
!     with allocatable components is a deep copy, which is what makes it
!     possible to re-run one filter update against a perturbed dipole.
      snap = st

      st = snap
      call mad_ir_push(st, munew)
      call mad_ir_evaluate(st, escale, loss, lambda)
      lam_ref = lambda
      ok = mad_ir_ready(st)
      if (.not. ok) then
         write (*, '(A)') "   FAIL the ensemble is not ready; the scan would be vacuous"
         n_fail = n_fail + 1
         call mad_ir_free(st)
         return
      end if

      best = huge(1.d0)
      do ih = 1, 6
         h = 1.d-2/(4.d0**(ih - 1))
         worst = 0.d0
         do c = 1, 3
            mup = munew
            mup(c) = munew(c) + h
            st = snap
            call mad_ir_push(st, mup)
            call mad_ir_evaluate(st, escale, lp, lambda)

            mup = munew
            mup(c) = munew(c) - h
            st = snap
            call mad_ir_push(st, mup)
            call mad_ir_evaluate(st, escale, lm, lambda)

            fd = (lp - lm)/(2.d0*h)
            worst = max(worst, dabs(fd - lam_ref(c))/max(1.d-30, dabs(lam_ref(c))))
         end do
         write (*, '(A,L8,L13,ES14.3,ES18.4)') "     ", subtract_mean, match_scale, h, worst
         best = min(best, worst)
      end do

      call check("      best over the h-scan      ", best, 0.d0, 1.d-7)

      call mad_ir_free(st)
      call mad_ir_free(snap)
   end subroutine gradient_hscan

!  ------------------------------------------------------------- 4. restart

   subroutine test_restart()
      type(mad_ir_type) :: st, st2
      real(dp) :: dt, tau, mu(1:3), worst
      integer :: n_lag, n_window, i, t, iu
      logical :: ok
      character(len=512) :: msg
      write (*, '(A)') "4. restart of the auxiliary variables"

      dt = 2.d0
      tau = 50.d0
      n_lag = 8
      n_window = 24

      call build(st, dt, n_lag, n_window, tau, .true., .true.)
      do i = 1, n_window + 5
         call dipole_at(i, dt, mu)
         call mad_ir_push(st, mu)
      end do
      call mad_ir_save(st, "aux_restart_test.dat", ok, msg)
      call check_true("saved                           ", ok)

      call build(st2, dt, n_lag, n_window, tau, .true., .true.)
      call mad_ir_load(st2, "aux_restart_test.dat", ok, msg)
      call check_true("loaded                          ", ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg)
      worst = 0.d0
      do t = 0, n_lag
         worst = max(worst, dabs(st2%acf(t) - st%acf(t)))
      end do
      call check("auxiliary variables round-trip  ", worst, 0.d0, 1.d-14)
      call check("update count round-trips        ", dfloat(st2%n_ema), dfloat(st%n_ema), 1.d-300)
      call check("running mean round-trips        ", &
                 maxval(dabs(st2%mu_bar_raw - st%mu_bar_raw)), 0.d0, 1.d-14)
      call mad_ir_free(st2)

!     A different memory constant is a different observable.
      call build(st2, dt, n_lag, n_window, 2.d0*tau, .true., .true.)
      call mad_ir_load(st2, "aux_restart_test.dat", ok, msg)
      call check_true("different tau_mem refused       ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(92, len_trim(msg))))
      call mad_ir_free(st2)

!     So is the block estimator's state.
      call build_block(st2, dt, n_lag, n_window)
      call mad_ir_load(st2, "aux_restart_test.dat", ok, msg)
      call check_true("block run refuses an aux file   ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(92, len_trim(msg))))
      call mad_ir_free(st2)

      call mad_ir_free(st)
      open (newunit=iu, file="aux_restart_test.dat", status="old"); close (iu, status="delete")
      write (*, *)
   end subroutine test_restart

   subroutine build_block(st, dt, n_lag, n_window)
      type(mad_ir_type), intent(inout) :: st
      real(dp), intent(in) :: dt
      integer, intent(in) :: n_lag, n_window
      real(dp), allocatable :: nu(:), Iexp(:), wgt(:)
      integer :: n_freq, k
      n_freq = 120
      allocate (nu(n_freq), Iexp(n_freq), wgt(n_freq))
      do k = 1, n_freq
         nu(k) = 3500.d0*dfloat(k)/dfloat(n_freq)
         Iexp(k) = 1.d3*dexp(-((nu(k) - 1500.d0)/220.d0)**2)
         wgt(k) = 1.d0
      end do
      call mad_ir_init(st, dt, n_lag, n_window, nu, Iexp, wgt, .true., 2.d0, &
                       "hann", .true., .true., .true., .false.)
      deallocate (nu, Iexp, wgt)
   end subroutine build_block

!  ------------------------------------------- 5. the two modes agree on a peak

   subroutine test_modes_agree()
      type(mad_ir_type) :: sa, sb
      real(dp) :: dt, mu(1:3), peak_a, peak_b
      integer :: n_lag, n_window, i
      write (*, '(A)') "5. block and exponential put the peak in the same place"

      dt = 2.d0
      n_lag = 40
      n_window = 160
!     tau_mem long compared with the buffer, so the exponential filter is
!     averaging over roughly what the block estimator averages over and the two
!     are estimating the same correlation.
      call build(sa, dt, n_lag, n_window, 400.d0, .true., .true.)
      call build_block(sb, dt, n_lag, n_window)

      do i = 1, n_window
         call dipole_at(i, dt, mu)
         call mad_ir_push(sa, mu)
         call mad_ir_push(sb, mu)
      end do
      call mad_ir_spectrum(sa)
      call mad_ir_spectrum(sb)
      peak_a = sa%nu(maxloc(sa%I_calc, 1))
      peak_b = sb%nu(maxloc(sb%I_calc, 1))
      write (*, '(A,F9.2,A,F9.2,A)') "        exponential ", peak_a, &
         " cm^-1     block ", peak_b, " cm^-1"
!     The drive is at 1500 cm^-1; both must find it, and the grid spacing is
!     3500/120 = 29 cm^-1 so agreement means the same bin or its neighbour.
      call check("peaks agree                     ", peak_a, peak_b, 60.d0)
      call check("and it is the driven frequency  ", peak_a, 1500.d0, 60.d0)

      call mad_ir_free(sa)
      call mad_ir_free(sb)
      write (*, *)
   end subroutine test_modes_agree

!  ------------------------------------------------------------ 6. readiness

   subroutine test_readiness()
      type(mad_ir_type) :: st
      real(dp) :: dt, mu(1:3), loss, lambda(1:3)
      integer :: n_lag, n_window, i
      logical :: seen_ready
      write (*, '(A)') "6. no bias until the filter has run"

      dt = 2.d0
      n_lag = 10
      n_window = 30
      call build(st, dt, n_lag, n_window, 50.d0, .true., .true.)

      seen_ready = .false.
      do i = 1, n_window - 1
         call dipole_at(i, dt, mu)
         call mad_ir_push(st, mu)
         if (mad_ir_ready(st)) seen_ready = .true.
         call mad_ir_evaluate(st, escale, loss, lambda)
         if (maxval(dabs(lambda)) > 0.d0) then
            write (*, '(A,I0)') "   FAIL a bias was produced at frame ", i
            n_fail = n_fail + 1
            exit
         end if
      end do
      call check_true("withheld while filling          ",.not. seen_ready)

      do i = n_window, n_window + 3
         call dipole_at(i, dt, mu)
         call mad_ir_push(st, mu)
      end do
      call mad_ir_evaluate(st, escale, loss, lambda)
      call check_true("applied once the buffer is full ", mad_ir_ready(st))
      call check_true("and lambda is non-zero          ", maxval(dabs(lambda)) > 0.d0)
      write (*, '(A,ES12.4,A,3ES12.4)') "        loss = ", loss, "   lambda = ", lambda

      call mad_ir_free(st)
      write (*, *)
   end subroutine test_readiness

end program auxverify
