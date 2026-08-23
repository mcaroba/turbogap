! The MAD IR bias in AUXILIARY-VARIABLE form: ir_bias_mode = xl.
!
! Here the spectrum is not transformed out of a stored trajectory. A bank of
! damped resonators, a quadrature pair per fitted frequency, is integrated
! alongside the atoms and driven by the dipole:
!
!   xddot + gamma xdot + w_k^2 x = m(t)
!   yddot + gamma ydot + w_k^2 y = mdot(t) / w_k
!
! and the experiment is compared with R_k^2 = |x_k|^2 + |y_k|^2. Nothing is
! shared with the autocorrelation path except the fitted grid, so none of
! irverify's or auxverify's coverage carries over.
!
! What has to be true, in increasing order of how quietly it can be wrong:
!
!   1. THE PROPAGATOR IS EXACT. The bank is advanced by the closed-form
!      solution of the damped oscillator over one interval with the drive held
!      constant, not by a Verlet step. Against a constant drive -- for which
!      the zero-order hold is not an approximation at all -- it must reproduce
!      the analytic transient to machine precision at any step size, including
!      step sizes where Verlet has long since gone unstable.
!
!   2. THE RESPONSE IS THE ONE THE PREFACTOR ASSUMES: driven on resonance by a
!      tone of amplitude M0, a resonator settles at M0/(gamma w). Every
!      intensity the bank reports is proportional to its square.
!
!   3. THE QUADRATURE PAIR IS A QUADRATURE PAIR. x^2 alone oscillates at 2w, so
!      x^2 + y^2 must be flat where x^2 ripples by order unity. This is also
!      the test that catches a first-order backward difference for mdot: that
!      form estimates the derivative at t - h/2 rather than at t, and half a
!      step of phase error puts the pair out of quadrature by w h / 2, which
!      shows up here as a ripple of exactly that size and nowhere else.
!
!   4. THE PREFACTOR CARRIES EXACTLY THE RIGHT POWER OF w. I = [nu^p gamma w^2]
!      R^2, and the w^2 is there to cancel the 1/(gamma w)^2 of the response.
!      Omit it and the predicted spectrum is TILTED by w^4 across the fitted
!      range. It still looks like a spectrum, it still has bands in the right
!      places, and no single fitted scale can absorb it. Two tones of equal
!      amplitude at different frequencies, fitted with nu_power = 0, must give
!      equal intensities; that is a factor of nine apart at 1000 and 3000
!      cm^-1 if the w^2 is missing.
!
!   5. COHERENT IS NOT INCOHERENT. Two sites in antiphase have zero total
!      dipole and therefore no infrared absorption, but each has a perfectly
!      good local one. The coherent bank must see the cancellation and the
!      incoherent one must not. If they agree here, the sum is being taken in
!      the wrong place and every spectrum the bank produces has lost the
!      interference between sites.
!
!   6. THE WEIGHT IS THE GRADIENT OF THE LOSS. This is the test that matters,
!      and it is the reason the bias is built as a gradient at all.
!
!      The first version of this scheme put the restraint on the auxiliary
!      coordinates, as an extended Lagrangian in R^2, and it biased the wrong
!      way: dU/dx is proportional to x, so it is not a force on the resonator
!      but a shift in its spring constant, and softening a resonator DETUNES it
!      from the drive it was measuring. A band that was too weak got weaker --
!      by two orders of magnitude, monotonically -- and the resonator stopped
!      reporting the intensity at w_k while it happened. No choice of sign or
!      magnitude repairs that, because the response has a maximum at zero shift
!      rather than being monotone in it.
!
!      The sign of a spectral bias cannot be reasoned about from the phase of a
!      filtered signal; the paragraph above is what that looks like when it
!      goes wrong. A gradient is the only construction that is right by
!      definition, so the weight is one, and this is the h-scan that says so:
!      central differences of the loss in the driving dipole, over both
!      amplitude modes and both settings of the fitted scale. It must fall as
!      h^2 and then turn on round-off. Anything flat in h is a systematic
!      error, whatever the spectra look like.
!
!   7. Warm-up and restart: the bank applies no weight while charging, and a
!      bank written under a different filter is refused rather than adopted.
!
program xlverify

   use kinds
   use mad_ir
   use mad_ir_xl

   implicit none

   integer :: n_fail
   real(dp) :: escale
   real(dp), parameter :: TWOPI = 6.283185307179586d0

   n_fail = 0
   escale = 3.7d0

   call test_propagator()
   call test_resonant_response()
   call test_quadrature_envelope()
   call test_prefactor_power()
   call test_coherent_vs_incoherent()
   call test_gradient()
   call test_warmup()
   call test_restart()

   write (*, *)
   if (n_fail == 0) then
      write (*, '(A)') "==> xlverify: all checks passed"
   else
      write (*, '(A,I0,A)') "==> xlverify: ", n_fail, " CHECK(S) FAILED"
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

!  ------------------------------------------------------------------
!  A parent observable holding nothing but a fitted grid. The bank takes nu,
!  I_exp, wgt, nu_power, match_scale and match_offset from it and uses no other
!  part of it, so the autocorrelation sizing here is whatever keeps mad_ir_init
!  happy.
!
!  nu_lo/nu_hi are given explicitly because several tests need the grid to
!  contain a particular frequency exactly: a resonator tuned a quarter of a
!  linewidth off a tone is a different amplitude, and the test would be
!  measuring the grid rather than the response.
!
!  band /= 0 puts a Gaussian band in the experiment so that the residual is not
!  a constant and the gradient is not dominated by one frequency.
!
   subroutine build_parent(st, nu_lo, nu_hi, n_freq, nu_power, match_scale, &
                           Iexp_amp, band)
      type(mad_ir_type), intent(inout) :: st
      real(dp), intent(in) :: nu_lo, nu_hi, nu_power, Iexp_amp, band
      integer, intent(in) :: n_freq
      logical, intent(in) :: match_scale
      real(dp), allocatable :: nu(:), Iexp(:), wgt(:)
      integer :: k
      allocate (nu(n_freq), Iexp(n_freq), wgt(n_freq))
      do k = 1, n_freq
         if (n_freq == 1) then
            nu(k) = nu_lo
         else
            nu(k) = nu_lo + (nu_hi - nu_lo)*dfloat(k - 1)/dfloat(n_freq - 1)
         end if
         if (band > 0.d0) then
            Iexp(k) = Iexp_amp*dexp(-((nu(k) - band)/220.d0)**2)
         else
            Iexp(k) = Iexp_amp
         end if
         wgt(k) = 1.d0
      end do
      call mad_ir_init(st, 1.d0, 8, 32, nu, Iexp, wgt, match_scale, nu_power, &
                       "hann", .true., .true., .true., .false.)
   end subroutine build_parent

!  Wavenumber -> angular frequency, the one conversion the whole file turns on.
   real(dp) function omega_of(nu)
      real(dp), intent(in) :: nu
      omega_of = TWOPI*nu/CM_PER_INV_FS
   end function omega_of

!  ==================================================================
!  1. The propagator, against the analytic step response.
!  ==================================================================
!
!  A resonator released from rest under a constant force F has
!
!    x(t) = (F/w^2) [ 1 - e^{-lam t} ( cos(Om t) + (lam/Om) sin(Om t) ) ]
!
!  with lam = gamma/2 and Om = sqrt(w^2 - lam^2). A constant force is exactly
!  what the zero-order hold represents, so this is not an accuracy test with a
!  tolerance to argue about -- the integrator either reproduces it to round-off
!  or it is not the propagator it claims to be.
!
!  Run at w*h = 3.8, past where velocity Verlet diverges (its stability limit
!  is w*h = 2), because unconditional stability is the reason for choosing a
!  closed form over a Verlet step in the first place.
!
   subroutine test_propagator()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, lam, Om, h, F, xa, t, worst, err
      real(dp), allocatable :: mu(:, :)
      integer :: i, nstep
      write (*, '(A)') "1. the propagator is the closed form, not an integrator"

      nu0 = 3000.d0
      w0 = omega_of(nu0)
      h = 3.8d0/w0
      tau = 4000.d0
      call build_parent(parent, nu0, nu0, 1, 2.d0, .false., 1.d0, 0.d0)

      call mad_ir_xl_init(mad_ir_xl_state, parent, 1, 1, h, tau, .true., 0.d0)
      allocate (mu(1:3, 1:1))
      mu = 0.d0
      mu(1, 1) = 0.5d0
      F = mu(1, 1)

      lam = 0.5d0*mad_ir_xl_state%gamma
      Om = dsqrt(w0**2 - lam**2)
      write (*, '(A,F8.4,A,ES12.5,A)') "        w*h = ", w0*h, &
         "   (Verlet is unstable above 2)   gamma = ", mad_ir_xl_state%gamma, " 1/fs"

      worst = 0.d0
!     Ten memory times, so that the stationary check below is actually at the
!     stationary value; the step-response comparison is exact at every step and
!     does not care how many there are.
      nstep = nint(10.d0*tau/h)
      do i = 1, nstep
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
         t = dfloat(i)*h
         xa = (F/w0**2)*(1.d0 - dexp(-lam*t)*(dcos(Om*t) + (lam/Om)*dsin(Om*t)))
         err = dabs(mad_ir_xl_state%x(1, 1, 1) - xa)/max(1.d-300, dabs(F/w0**2))
         if (err > worst) worst = err
      end do
      call check("step response, worst rel err    ", worst, 0.d0, 1.d-12)
!     After ten memory times the transient envelope is exp(-t/tau) = 4.5e-5, so
!     this is a check that the fixed point is F/w^2 and not a check of the
!     propagator -- the exact comparison above is that, to 1e-13.
      call check("stationary displacement F/w^2   ", mad_ir_xl_state%x(1, 1, 1), &
                 F/w0**2, 1.d-4)

!     And the sensitivity the gradient uses, dx(t+h)/dF, against a difference
!     of the propagator itself. sens is twice it, from equation (6).
      call check("sens = 2 dx/dF                  ", &
                 mad_ir_xl_state%sens(1), 2.d0*dxdF_numeric(w0, lam, h), 1.d-9)

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_propagator

!  One interval of the propagator from rest, differentiated in F by hand.
   real(dp) function dxdF_numeric(w0, lam, h)
      real(dp), intent(in) :: w0, lam, h
      real(dp) :: Om, F, x1, x2, d
      Om = dsqrt(w0**2 - lam**2)
      d = 1.d-6
      F = 1.d0
      x1 = one_step_from_rest(F + d, w0, lam, Om, h)
      x2 = one_step_from_rest(F - d, w0, lam, Om, h)
      dxdF_numeric = (x1 - x2)/(2.d0*d)
   end function dxdF_numeric

   real(dp) function one_step_from_rest(F, w0, lam, Om, h)
      real(dp), intent(in) :: F, w0, lam, Om, h
      real(dp) :: xp, u0
      xp = F/w0**2
      u0 = -xp
      one_step_from_rest = xp + dexp(-lam*h)*(u0*dcos(Om*h) + lam*u0*dsin(Om*h)/Om)
   end function one_step_from_rest

!  ==================================================================
!  2. Resonant response: |x| = M0 / (gamma w).
!  ==================================================================
!
!  Drive a resonator tuned to w0 with a tone of amplitude M0 at w0 and let it
!  settle. This is the identity the whole R^2 -> intensity conversion rests on.
!
!  The zero-order hold makes the tone a staircase, which is the tone times
!  sinc(w h / 2) plus images; h is taken small enough (w*h ~ 0.13) that the
!  correction is under a per cent, and the test is given the correction rather
!  than pretending it is absent.
!
   subroutine test_resonant_response()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, h, M0, t, amp, want, sinc
      real(dp), allocatable :: mu(:, :)
      integer :: i, nstep
      write (*, '(A)') "2. the resonant response is M0 / (gamma w)"

      nu0 = 2000.d0
      w0 = omega_of(nu0)
      h = 0.35d0
      tau = 2000.d0
      M0 = 0.4d0
      call build_parent(parent, nu0, nu0, 1, 2.d0, .false., 1.d0, 0.d0)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 1, 1, h, tau, .true., 0.d0)
      allocate (mu(1:3, 1:1))

!     Ten memory times: the transient is exp(-t/tau) and this leaves it at 5e-5.
      nstep = nint(10.d0*tau/h)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         mu = 0.d0
         mu(1, 1) = M0*dcos(w0*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do

      amp = dsqrt(mad_ir_xl_state%x(1, 1, 1)**2 + mad_ir_xl_state%y(1, 1, 1)**2)
      want = M0/(mad_ir_xl_state%gamma*w0)
      sinc = dsin(0.5d0*w0*h)/(0.5d0*w0*h)
      write (*, '(A,F10.6,A,F10.6)') "        zero-order-hold sinc(wh/2) = ", sinc, &
         "   w*h = ", w0*h
      call check("envelope amplitude              ", amp, want*sinc, 5.d-3)

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_resonant_response

!  ==================================================================
!  3. x^2 + y^2 is an envelope; x^2 is not.
!  ==================================================================
!
!  The peak-to-mean ripple of x^2 over one drive period is 1 by construction --
!  it is a squared cosine. The pair's ripple measures how nearly y is in
!  quadrature with x, and so directly measures the phase error in the mdot
!  estimate: a first-order backward difference leaves w h / 2 of it, which at
!  w h = 0.16 is a 6% ripple that looks like noise.
!
   subroutine test_quadrature_envelope()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, h, M0, t
      real(dp) :: x2min, x2max, r2min, r2max, x2, r2, rip_x, rip_r
      real(dp), allocatable :: mu(:, :)
      integer :: i, nstep, nper
      write (*, '(A)') "3. the quadrature pair gives an envelope, x alone does not"

      nu0 = 2500.d0
      w0 = omega_of(nu0)
      h = 0.25d0
      tau = 1500.d0
      M0 = 0.3d0
      call build_parent(parent, nu0, nu0, 1, 2.d0, .false., 1.d0, 0.d0)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 1, 1, h, tau, .true., 0.d0)
      allocate (mu(1:3, 1:1))

      nstep = nint(10.d0*tau/h)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         mu = 0.d0
         mu(1, 1) = M0*dcos(w0*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do

      nper = max(4, nint(TWOPI/(w0*h)))
      x2min = huge(1.d0); x2max = -huge(1.d0)
      r2min = huge(1.d0); r2max = -huge(1.d0)
      do i = 1, nper
         t = dfloat(nstep + i - 1)*h
         mu = 0.d0
         mu(1, 1) = M0*dcos(w0*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
         x2 = mad_ir_xl_state%x(1, 1, 1)**2
         r2 = x2 + mad_ir_xl_state%y(1, 1, 1)**2
         x2min = min(x2min, x2); x2max = max(x2max, x2)
         r2min = min(r2min, r2); r2max = max(r2max, r2)
      end do
      rip_x = (x2max - x2min)/(x2max + x2min)
      rip_r = (r2max - r2min)/(r2max + r2min)
      write (*, '(A,F10.6,A,F10.6,A,F10.6)') "        ripple x^2 = ", rip_x, &
         "   ripple x^2+y^2 = ", rip_r, "   w*h/2 = ", 0.5d0*w0*h
      call check_true("x^2 ripples by order unity      ", rip_x > 0.8d0)
      call check_true("x^2+y^2 is flat to 1%           ", rip_r < 1.d-2)
      call check_true("and well below the wh/2 a two-  ", rip_r < 0.25d0*w0*h)
      write (*, '(A)') "        point mdot would leave"

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_quadrature_envelope

!  ==================================================================
!  4. THE PREFACTOR. Two tones, equal amplitude, nu_power = 0.
!  ==================================================================
!
   subroutine test_prefactor_power()
      type(mad_ir_type) :: parent
      real(dp) :: nu1, nu2, w1, w2, tau, h, M0, t, I1, I2, e, s1, s2
      real(dp), allocatable :: mu(:, :)
      integer :: i, nstep
      write (*, '(A)') "4. the prefactor carries exactly w^2 (the tilt test)"

      nu1 = 1000.d0
      nu2 = 3000.d0
      w1 = omega_of(nu1)
      w2 = omega_of(nu2)
      h = 0.3d0
      tau = 1500.d0
      M0 = 0.3d0

      call build_parent(parent, nu1, nu2, 2, 0.d0, .false., 1.d0, 0.d0)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 1, 2, h, tau, .true., 0.d0)
      call check("mode 1 sits on the first tone   ", mad_ir_xl_state%nu(1), nu1, 1.d-12)
      call check("mode 2 sits on the second tone  ", mad_ir_xl_state%nu(2), nu2, 1.d-12)

      allocate (mu(1:3, 1:1))
      nstep = nint(10.d0*tau/h)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         mu = 0.d0
!        Two tones, equal amplitude, on different Cartesian components so that
!        neither resonator sees any part of the other's drive.
         mu(1, 1) = M0*dcos(w1*t)
         mu(2, 1) = M0*dcos(w2*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do

      call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
      I1 = mad_ir_xl_state%I_calc(1)
      I2 = mad_ir_xl_state%I_calc(2)
      s1 = dsin(0.5d0*w1*h)/(0.5d0*w1*h)
      s2 = dsin(0.5d0*w2*h)/(0.5d0*w2*h)
      write (*, '(A,ES13.6,A,ES13.6)') "        I(1000) = ", I1, "   I(3000) = ", I2
!     sinc(wh/2) is the leading zero-order-hold correction, not the whole of
!     it -- the staircase also folds images back in -- so 2% is left for the
!     rest. The failure this catches is a factor of 3 or 9, not 2%.
      call check("I(3000)/I(1000), sinc-corrected ", (I2/I1)*(s1/s2)**2, 1.d0, 3.d-2)
!     A line of dipole amplitude M0 has I = nu^p M0^2 / gamma; with p = 0 that
!     is a number this test can state outright rather than compare to itself.
      call check("absolute line intensity, mode 1 ", I1, &
                 M0**2/mad_ir_xl_state%gamma*s1**2, 1.d-2)

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_prefactor_power

!  ==================================================================
!  5. Coherent sees the cancellation; incoherent does not.
!  ==================================================================
!
   subroutine test_coherent_vs_incoherent()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, h, M0, t, e, I_coh, I_inc
      real(dp), allocatable :: mu(:, :)
      integer :: i, nstep
      write (*, '(A)') "5. coherent cancels antiphase sites, incoherent does not"

      nu0 = 2000.d0
      w0 = omega_of(nu0)
      h = 0.3d0
      tau = 1200.d0
      M0 = 0.35d0
      nstep = nint(10.d0*tau/h)
      call build_parent(parent, nu0, nu0, 1, 2.d0, .false., 1.d0, 0.d0)
      allocate (mu(1:3, 1:2))

      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 1, h, tau, .true., 0.d0)
      call check("coherent needs one bank replica ", &
                 dfloat(mad_ir_xl_state%n_bank), 1.d0, 0.d0)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         mu = 0.d0
         mu(1, 1) = M0*dcos(w0*t)
         mu(1, 2) = -M0*dcos(w0*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do
      call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
      I_coh = mad_ir_xl_state%I_calc(1)
      call mad_ir_xl_free(mad_ir_xl_state)

      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 1, h, tau, .false., 0.d0)
      call check("incoherent needs one per site   ", &
                 dfloat(mad_ir_xl_state%n_bank), 2.d0, 0.d0)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         mu = 0.d0
         mu(1, 1) = M0*dcos(w0*t)
         mu(1, 2) = -M0*dcos(w0*t)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do
      call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
      I_inc = mad_ir_xl_state%I_calc(1)

      write (*, '(A,ES13.6,A,ES13.6)') "        I_coherent = ", I_coh, &
         "   I_incoherent = ", I_inc
      call check_true("coherent cancels to < 1e-20 of  ", I_coh < 1.d-20*I_inc)
      call check_true("incoherent keeps both sites     ", I_inc > 0.d0)

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_coherent_vs_incoherent

!  ==================================================================
!  6. THE GRADIENT, against central differences of the loss.
!  ==================================================================
!
!  The claim is exact: the weight is -dU/dm of the CURRENT dipole, where U is
!  the loss evaluated after that dipole has been folded into the bank. So drive
!  a bank into a stationary state, freeze it, and then for each Cartesian
!  component perturb the next dipole by +-h, advance, evaluate, and compare the
!  central difference of the loss with -W.
!
!  Freezing is a whole-derived-type assignment, which deep-copies the
!  allocatable components: the perturbed runs must each start from the same
!  bank, and re-driving it from rest for each of them would take longer than
!  the rest of the file put together.
!
!  Run over both amplitude modes and both settings of the fitted scale, because
!  the scale enters dLdI and the envelope-theorem argument that it contributes
!  no extra term is exactly the kind of claim an h-scan is for.
!
   subroutine test_gradient()
      write (*, '(A)') "6. the weight is -dU/dm, against central differences"
      call gradient_case(.true., .true., "coherent  , scale on ")
      call gradient_case(.true., .false., "coherent  , scale off")
      call gradient_case(.false., .true., "incoherent, scale on ")
      call gradient_case(.false., .false., "incoherent, scale off")
      write (*, *)
   end subroutine test_gradient

   subroutine gradient_case(coherent, match_scale, label)
      logical, intent(in) :: coherent, match_scale
      character(len=*), intent(in) :: label
      type(mad_ir_type) :: parent
      type(mad_ir_xl_type) :: frozen
      real(dp) :: nu0, w0, tau, h, M0, t, e, hd, up, dn, num, ana, err, best
      real(dp), allocatable :: mu(:, :), w_out(:, :), mu_next(:, :)
      integer :: i, nstep, isite, a, ih
      integer, parameter :: n_site_test = 3

      nu0 = 2000.d0
      w0 = omega_of(nu0)
      h = 1.5d0
      tau = 400.d0
      M0 = 0.3d0

!     A grid with a band in it, so the residual varies across the modes and the
!     gradient is a sum of terms that do not all have the same sign.
      call build_parent(parent, 800.d0, 3200.d0, 40, 2.d0, match_scale, 5.d2, 1600.d0)
      call mad_ir_xl_init(mad_ir_xl_state, parent, n_site_test, 12, h, tau, &
                          coherent, 0.d0)
      allocate (mu(1:3, 1:n_site_test), mu_next(1:3, 1:n_site_test))
      allocate (w_out(1:3, 1:n_site_test))

      nstep = nint(8.d0*tau/h)
      do i = 1, nstep
         t = dfloat(i - 1)*h
         call signal_at(t, w0, M0, mu)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do
      frozen = mad_ir_xl_state

!     The reference weight: advance once more with mu_next, evaluate, and take
!     the weight from the resulting state. That is exactly the driver's order.
      t = dfloat(nstep)*h
      call signal_at(t, w0, M0, mu_next)
      call mad_ir_xl_advance(mad_ir_xl_state, mu_next)
      call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
      call mad_ir_xl_weights(mad_ir_xl_state, w_out)

!     Site 2 and the y component: an arbitrary but fixed choice, so the number
!     printed is the same one every run. In coherent mode every site has the
!     same weight and the choice does not matter; in incoherent mode it does.
      isite = 2
      a = 2
      ana = -w_out(a, isite)

      best = huge(1.d0)
      write (*, '(A,A,A,ES14.6)') "        ", label, "   analytic dU/dm = ", ana
      do ih = 1, 7
         hd = 1.d-2/(4.d0**(ih - 1))
         mad_ir_xl_state = frozen
         mu(1:3, 1:n_site_test) = mu_next(1:3, 1:n_site_test)
         mu(a, isite) = mu(a, isite) + hd
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
         call mad_ir_xl_evaluate(mad_ir_xl_state, escale, up)

         mad_ir_xl_state = frozen
         mu(1:3, 1:n_site_test) = mu_next(1:3, 1:n_site_test)
         mu(a, isite) = mu(a, isite) - hd
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
         call mad_ir_xl_evaluate(mad_ir_xl_state, escale, dn)

         num = (up - dn)/(2.d0*hd)
         err = dabs(num - ana)/max(1.d-300, dabs(ana))
         write (*, '(A,ES9.2,A,ES14.6,A,ES9.2)') "          h = ", hd, &
            "   numeric = ", num, "   rel err = ", err
         if (err < best) best = err
      end do
      call check_true("   "//label//": h-scan reaches 1e-8", best < 1.d-8)

      deallocate (mu, mu_next, w_out)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_free(frozen)
      call mad_ir_free(parent)
   end subroutine gradient_case

!  A three-site dipole signal with structure: different amplitudes, different
!  phases, two frequencies, and a static offset, so that no accidental symmetry
!  makes a wrong gradient look right.
   subroutine signal_at(t, w0, M0, mu)
      real(dp), intent(in) :: t, w0, M0
      real(dp), intent(out) :: mu(:, :)
      integer :: i
      real(dp) :: ph
      do i = 1, size(mu, 2)
         ph = 0.7d0*dfloat(i)
         mu(1, i) = M0*(1.d0 + 0.3d0*dfloat(i))*dcos(w0*t + ph)
         mu(2, i) = M0*0.6d0*dsin(1.4d0*w0*t + ph) + 0.05d0*dfloat(i)
         mu(3, i) = M0*0.25d0*dcos(0.6d0*w0*t - ph)
      end do
   end subroutine signal_at

!  ==================================================================
!  7. No weight while the bank is charging.
!  ==================================================================
!
   subroutine test_warmup()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, h, M0, t, e
      real(dp), allocatable :: mu(:, :), w_out(:, :)
      integer :: i
      logical :: any_nonzero
      write (*, '(A)') "7. no weight while charging"

      nu0 = 2000.d0
      w0 = omega_of(nu0)
      h = 1.d0
      tau = 40.d0
      M0 = 0.3d0
      call build_parent(parent, 1500.d0, 2500.d0, 6, 2.d0, .false., 1.d0, 0.d0)
!     warm_factor 2 with tau/h = 40 gives n_warm = 80 advances.
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 6, h, tau, .true., 2.d0)
      call check("n_warm from warm_factor*tau/dt  ", dfloat(mad_ir_xl_state%n_warm), &
                 80.d0, 1.d-12)

      allocate (mu(1:3, 1:2), w_out(1:3, 1:2))
      any_nonzero = .false.
      do i = 1, mad_ir_xl_state%n_warm - 1
         t = dfloat(i - 1)*h
         call signal_at(t, w0, M0, mu)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
         call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
         call mad_ir_xl_weights(mad_ir_xl_state, w_out)
         if (maxval(dabs(w_out)) > 0.d0) any_nonzero = .true.
      end do
      call check_true("weight is zero while charging   ",.not. any_nonzero)
      call check_true("and the bank is not ready       ",.not. mad_ir_xl_ready(mad_ir_xl_state))
      call check("and the loss is withheld too    ", e, 0.d0, 0.d0)

      call mad_ir_xl_advance(mad_ir_xl_state, mu)
      call mad_ir_xl_evaluate(mad_ir_xl_state, escale, e)
      call mad_ir_xl_weights(mad_ir_xl_state, w_out)
      call check_true("ready once charged              ", mad_ir_xl_ready(mad_ir_xl_state))
      call check_true("and the weight is non-zero      ", maxval(dabs(w_out)) > 0.d0)
!     Coherent: every site carries the same weight, which is the honest form of
!     the statement that localisation buys nothing against a global target.
      call check("coherent weights are site-free  ", w_out(1, 1), w_out(1, 2), 0.d0)

      deallocate (mu, w_out)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      write (*, *)
   end subroutine test_warmup

!  ==================================================================
!  8. Restart: round trip, and every refusal.
!  ==================================================================
!
   subroutine test_restart()
      type(mad_ir_type) :: parent
      real(dp) :: nu0, w0, tau, h, M0, t, x_saved
      real(dp), allocatable :: mu(:, :)
      integer :: i
      logical :: ok, resumed
      character(len=512) :: msg
      write (*, '(A)') "8. the bank restarts, and refuses a different filter"

      nu0 = 2000.d0
      w0 = omega_of(nu0)
      h = 1.d0
      tau = 60.d0
      M0 = 0.3d0
      call build_parent(parent, 1500.d0, 2500.d0, 5, 2.d0, .false., 1.d0, 0.d0)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, h, tau, .false., 1.d0)
      allocate (mu(1:3, 1:2))
      do i = 1, 200
         t = dfloat(i - 1)*h
         call signal_at(t, w0, M0, mu)
         call mad_ir_xl_advance(mad_ir_xl_state, mu)
      end do
      x_saved = mad_ir_xl_state%x(1, 3, 2)
      call mad_ir_xl_save(mad_ir_xl_state, "xlverify_bank.dat", ok, msg)
      call check_true("saved                           ", ok)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, h, tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("resumed                         ", resumed)
      call check("x came back bit for bit         ", mad_ir_xl_state%x(1, 3, 2), &
                 x_saved, 0.d0)
      call check("and so did the advance count    ", dfloat(mad_ir_xl_state%n_steps), &
                 200.d0, 0.d0)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, h, 2.d0*tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("refuses a different tau_mem     ",.not. resumed)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, 2.d0*h, tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("refuses a different interval    ",.not. resumed)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, h, tau, .true., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("refuses the other amplitude     ",.not. resumed)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 3, 5, h, tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("refuses a different site count  ",.not. resumed)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 3, h, tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "xlverify_bank.dat", resumed, msg)
      call check_true("refuses a different mode count  ",.not. resumed)

      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_xl_init(mad_ir_xl_state, parent, 2, 5, h, tau, .false., 1.d0)
      call mad_ir_xl_load(mad_ir_xl_state, "no_such_bank.dat", resumed, msg)
      call check_true("a missing file is not fatal     ",.not. resumed)

      deallocate (mu)
      call mad_ir_xl_free(mad_ir_xl_state)
      call mad_ir_free(parent)
      open (unit=91, file="xlverify_bank.dat", status="old")
      close (91, status="delete")
      write (*, *)
   end subroutine test_restart

end program xlverify
