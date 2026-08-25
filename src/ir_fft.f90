! IR spectra by the Wiener-Khinchin route: a Fortran translation of the
! GPUMD/TNEP pipeline in TNEP/spectroscopy.py (compute_dipole_acf and
! compute_ir_spectrum), following Xu et al., J. Chem. Theory Comput. 20, 3273
! (2024).
!
! This is a SECOND IR estimator, not a replacement for mad_ir.f90. Both take
! the same physics -- the spectrum is a cosine transform of a dipole
! autocorrelation -- but they differ in four ways that matter, and the point of
! having both is to be able to see whether those differences move the answer:
!
!   1. THE CORRELATION IS COMPUTED BY FFT, over the whole trajectory at once,
!      rather than lag by lag over a rolling buffer. Same estimator (see
!      "biased" below), O(T log T) instead of O(T * n_lag).
!   2. THE LONGEST LAG IS A FRACTION OF THE TRAJECTORY (acf_ratio, default 0.1)
!      rather than a resolution asked for up front. That is the GPUMD
!      convention and it makes the resolution an OUTPUT of the run length.
!   3. THE QUANTUM CORRECTION IS EXPLICIT AND DEFAULTS TO "harmonic",
!      w (1 - exp(-hbar w / kT)), rather than mad_ir's fixed nu^p with p = 2.
!      These agree only for hbar w << kT. See the block below; it is the single
!      largest difference between the two spectra on water.
!   4. THE RESULT IS SMOOTHED AND PEAK-NORMALISED, so what comes out is a
!      lineshape to overlay on an experiment rather than an absolute
!      alpha(nu) n(nu).
!
! The module is deliberately free of I/O, MPI and of any TurboGAP type: it
! takes a dipole trajectory and gives back a spectrum, so it can be driven from
! a test program and checked against the Python line by line. ir_fft_loss adds
! the adjoint of the whole pipeline, which is what makes it usable as a MAD
! bias.
!
!==========================================================================
! THE PIPELINE, IN THE ORDER IT RUNS
!==========================================================================
!
! Input: mu(1:3, 1:T), CHRONOLOGICAL -- mu(:,1) is the oldest frame and
! mu(:,T) the newest -- sampled every dt fs.
!
!   (1) mu <- mu - <mu>_t                      subtract the mean of each
!                                              Cartesian component over time
!
!   (2) C(tau) = (1/T) sum_{a=0}^{T-1-tau} mu(a) . mu(a+tau)      tau = 0..T-1
!
!       computed as irfft(|rfft(mu, n)|^2)[0:T] summed over x,y,z, with
!       n >= 2T-1 so the circular correlation the FFT computes is the linear
!       one that is wanted. Note the divisor is the CONSTANT T, not T - tau:
!       this is the BIASED estimator, and it is not a mistake. C(tau) with
!       divisor T is a periodogram, hence a positive semi-definite sequence,
!       hence a transform that cannot go negative; with divisor T - tau the
!       long-lag end is amplified exactly where it is least well determined.
!       GPUMD, Xu et al. and mad_ir's ir_estimator = biased all make the same
!       choice for the same reason.
!
!   (3) L = int(T * acf_ratio),  keep C(0..L-1)
!
!       The longest lag sets the frequency resolution and NOTHING ELSE does:
!
!           d(nu) = 33356.40952 / (L * dt)    cm^-1
!
!       A 100 ps run at acf_ratio = 0.1 resolves 3.3 cm^-1; the same run at
!       acf_ratio = 0.01 resolves 33 cm^-1 no matter how long it was. The
!       trade is variance: C(tau) at lag tau averages T - tau products, so the
!       last lag kept is averaged over T(1 - acf_ratio) of them and the
!       inflation over lag zero is capped at 1/(1 - acf_ratio).
!
!   (4) a(tau) = C(tau) * w(tau) * kron(tau),  kron(0) = 1, kron(tau>0) = 2
!
!       The Kronecker factor turns the one-sided sum into the two-sided
!       integral of equation (1) of mad_ir.f90's header: C is even, so
!       INT_{-inf}^{inf} = C(0) + 2 sum_{tau>0}. w is the lag window; see
!       ir_fft_window below for why "none" is not an option so much as a
!       different, worse, window.
!
!   (5) M(k) = sum_{tau=0}^{L-1} a(tau) cos(2 pi k tau / N2),   N2 = 2L - 1
!
!       the real part of the length-N2 DFT of a zero-padded to N2. This is the
!       power spectrum, M(w) >= 0 by construction of (2).
!
!   (6) nu_k = k * 33356.40952 / (N2 * dt)     cm^-1,  k = 0 .. L-1
!
!   (7) I(k) = P(nu_k) * M(k)                  quantum correction, below
!
!   (8) keep nu_k <= max_freq_cm, THEN smooth, THEN peak-normalise
!
!       The order matters: smoothing after the cut means the edge treatment
!       sees the cut edge, and normalising last means the peak is the peak of
!       what is returned. The Python does it in this order and so does this.
!
!==========================================================================
! THE QUANTUM CORRECTION, WHICH IS WHERE THE TWO ESTIMATORS PART COMPANY
!==========================================================================
!
! Linear response gives the absorption as
!
!   alpha(w) n(w) ~ w (1 - exp(-beta hbar w)) INT dt exp(-i w t) <dM(0).dM(t)>_qm
!
! and a classical trajectory cannot produce <.>_qm. mad_ir.f90 takes the
! harmonic quantum correction factor, which cancels (1 - exp(-beta hbar w))
! exactly and leaves w^2 C_cl(w); that is ir_nu_power = 2, and it is right.
!
! This module keeps the SAME physics but writes it the other way round, as
! GPUMD does: apply w (1 - exp(-beta hbar w)) to the classical lineshape
! directly. The two are the same correction seen from opposite sides --
!
!   w^2 C_cl   vs   w (1 - exp(-x)) C_cl,     x = hbar w / kT
!
! -- and they agree only in the limit x -> 0, where (1 - exp(-x)) -> x. At
! 300 K, kT = 25.85 meV and hbar c = 1.23984e-4 eV cm, so x = 1 at 208.5
! cm^-1. Below that the two forms agree; above it they diverge without bound,
! because (1 - exp(-x)) saturates at 1 while x does not:
!
!     nu (cm^-1)      x       w^2 form     w(1-e^-x) form     ratio
!        200        0.959      1.00            1.00           1.00
!       1000        4.80       1.00            0.208          4.8
!       1650        7.91       1.00            0.126          7.9
!       3400       16.3        1.00            0.0614        16.3
!
! (columns scaled to the w^2 form at each frequency). So the choice is not
! cosmetic: relative to "classical", "harmonic" suppresses the O-H stretch by
! a factor 16 against the librational region. Neither is more correct than the
! other -- they are different QCFs applied to the same classical ACF, and
! Ramirez et al. (JCP 121, 3973 (2004)) compared them -- but they answer to
! different experimental quantities, and mixing them up is the fastest way to
! conclude that a perfectly good dipole model has the wrong band intensities.
!
! quantum_correction = "classical" reproduces mad_ir with ir_nu_power = 2 up
! to the grid, the smoothing and the normalisation, and that is how the two
! implementations are cross-checked against each other.
!
!==========================================================================
! DEVIATIONS FROM spectroscopy.py, ALL DELIBERATE, ALL HERE
!==========================================================================
!
!   * BLACKMAN IS THE HALF WINDOW, NOT numpy's. The Python calls np.blackman(L)
!     and applies it to lags 0..L-1. np.blackman is SYMMETRIC: it is ~0 at
!     index 0, rises to 1 at L/2, and falls to ~0 at L-1. Applied to a lag
!     sequence that is a taper which DELETES C(0) -- the largest and best
!     determined term, the one carrying the total intensity -- and gives the
!     most weight to lags around L/2. That is not a lag window, and a spectrum
!     computed with it is not the spectrum. This module uses the half window
!
!         w(tau) = 0.42 + 0.5 cos(pi tau / L) + 0.08 cos(2 pi tau / L)
!
!     which is 1 at tau = 0 and 0 at tau = L, i.e. the right half of
!     np.blackman, and is what "Blackman lag window" means everywhere else.
!     window = "hann" and "none" agree with the Python exactly.
!
!   * THE FFT LENGTH IS THE NEXT POWER OF TWO >= 2T, where the Python asks
!     scipy for next_fast_len(2T-1). Any n >= 2T-1 gives the identical linear
!     correlation -- the padding only has to be long enough that the wrap-round
!     of the circular correlation cannot reach the lags being kept -- so this
!     is a difference in speed, not in result. It buys a radix-2 FFT with no
!     external dependency.
!
!   * 1/c IS 33356.40952 cm/fs (the module constant shared with mad_ir.f90),
!     against the Python's 2.99792458e-5 cm/fs for c. The two differ in the
!     12th significant figure.
!
!   * PEAK NORMALISATION IS OPTIONAL AND OFF FOR THE BIAS. Dividing by
!     max|I| is a discontinuous function of the trajectory: the derivative
!     jumps whenever the tallest bin changes, and under a bias that is
!     deliberately moving band heights around it changes often. The prediction
!     path normalises (and reports the divisor); the MAD path does not, and
!     fits a scale against the experiment instead, which is the same freedom
!     applied smoothly.
!
!==========================================================================
! COST
!==========================================================================
!
! Per call, with T frames, L = acf_ratio*T lags kept and K bins below
! max_freq_cm:
!
!   ACF          3 forward + 3 inverse FFTs of length 2^ceil(log2 2T)
!   transform    O(K * L)   -- direct, because N2 = 2L-1 is odd and a length-N2
!                             FFT would need Bluestein for no gain: K is a few
!                             hundred where L is thousands, and the bins above
!                             max_freq_cm are thrown away anyway, so evaluating
!                             only the ones that survive is cheaper than
!                             transforming all of them.
!   adjoint      O(K * L) again, plus O(T) for the prefix sums
!
! On 6000 frames with acf_ratio = 0.2 and max_freq_cm = 4000 that is L = 1200,
! K = 288, so ~7e5 multiply-adds for the transform: microseconds. It is the
! FFT length, not the transform, that grows with the trajectory.
!
module ir_fft

   use kinds

   implicit none

   private
   public :: ir_fft_config_type, ir_fft_result_type
   public :: ir_fft_check_config, ir_fft_free
   public :: ir_fft_autocorrelation, ir_fft_spectrum, ir_fft_loss
   public :: ir_fft_next_pow2, ir_fft_transform, ir_fft_resolution
   public :: ir_fft_n_lag, ir_fft_lag_window
   public :: IR_FFT_CM_PER_INV_FS, IR_FFT_HBAR_C_EV_CM, IR_FFT_KB_EV

!  Wavenumber in cm^-1 of a frequency of 1/fs. Same value as mad_ir's
!  CM_PER_INV_FS; duplicated rather than used from there so that this module
!  depends on nothing but kinds and can be driven from a bare test program.
   real(dp), parameter :: IR_FFT_CM_PER_INV_FS = 33356.40952d0
!  hbar*c in eV cm, so that hbar*w in eV is IR_FFT_HBAR_C_EV_CM * nu[cm^-1].
   real(dp), parameter :: IR_FFT_HBAR_C_EV_CM = 1.23984d-4
!  Boltzmann's constant in eV/K.
   real(dp), parameter :: IR_FFT_KB_EV = 8.617333d-5

   type :: ir_fft_config_type
!     Interval between stored frames, in fs. This is md_step*ir_stride for a
!     run, or the spacing of the time= labels for a trajectory read from disk.
      real(dp) :: dt_fs = 1.d0
!     Lag window: "hann" (default), "blackman" (the half window; see the
!     header), "none" (a boxcar, which rings), and mad_ir's "bartlett",
!     "welch" and "lorch" for comparison against that estimator.
      character(len=32) :: window = "hann"
!     Highest wavenumber returned. Bins above it are dropped BEFORE smoothing.
      real(dp) :: max_freq_cm = 4000.d0
!     Fraction of the trajectory kept as lags. Sets the resolution; see (3).
      real(dp) :: acf_ratio = 0.1d0
!     Smoothing strength. Gaussian FWHM in bins, or box width in bins. 0 or 1
!     disables it.
      integer :: smooth_k = 10
      character(len=32) :: smooth_kind = "gaussian"
!     Temperature for the harmonic correction, in K.
      real(dp) :: temperature = 300.d0
!     "harmonic" (default), "classical" (alias "quadratic"), "linear", "none".
      character(len=32) :: quantum_correction = "harmonic"
!     Bins below this are excluded from the power spectrum's peak normaliser.
!     M(w) has a large peak at w = 0 that would otherwise crush every
!     vibrational feature to a per cent of full scale. 0 keeps the DC bin.
      real(dp) :: power_dc_cutoff_cm = 100.d0
!     Correlate mu - <mu> rather than mu. Linear response wants the
!     fluctuation; the Python always does this and there is no good reason not
!     to, but the switch exists so the leakage it prevents can be demonstrated.
      logical :: subtract_mean = .true.
!     Peak-normalise the returned intensity and power. Off for the MAD bias;
!     see the header.
      logical :: normalise = .true.
   end type ir_fft_config_type

   type :: ir_fft_result_type
      integer :: n_freq = 0        ! bins returned (after the max_freq cut and
      !                              after box smoothing, which shortens)
!     Bins after the max_freq cut but BEFORE smoothing. Equal to n_freq except
!     under the box smoother, which shortens by smooth_k - 1. Carried rather
!     than inferred: the smoother is skipped when there are fewer bins than its
!     width, so n_freq alone does not say whether it ran, and an adjoint that
!     guesses wrong there is wrong silently.
      integer :: n_freq_raw = 0
      integer :: n_lag = 0         ! L, the lags kept
      integer :: n_frames = 0      ! T
      real(dp) :: dt_fs = 0.d0
      real(dp) :: resolution = 0.d0   ! cm^-1, = CM/(L dt)
      real(dp) :: nyquist = 0.d0      ! cm^-1, = CM/(2 dt)
      real(dp) :: d_nu = 0.d0         ! bin spacing, = CM/(N2 dt)
      real(dp) :: peak = 0.d0         ! divisor used on the intensity
      real(dp) :: peak_power = 0.d0   ! divisor used on the power
      real(dp) :: mu_mean(1:3) = 0.d0 ! the mean that was subtracted
      real(dp), allocatable :: freq(:)          ! (1:n_freq) cm^-1
      real(dp), allocatable :: intensity(:)     ! (1:n_freq) normalised if asked
      real(dp), allocatable :: intensity_raw(:) ! (1:n_freq) before normalising
      real(dp), allocatable :: power(:)         ! (1:n_freq) M(w), normalised
      real(dp), allocatable :: power_raw(:)     ! (1:n_freq) M(w)
      real(dp), allocatable :: acf(:)           ! (0:n_lag-1) C(tau), unwindowed
   end type ir_fft_result_type

contains

!**************************************************************************
!
! Smallest power of two that is at least n. Used for the FFT length; see the
! header on why any length >= 2T-1 is equivalent.
!
   function ir_fft_next_pow2(n) result(m)

      implicit none

      integer, intent(in) :: n
      integer :: m

      m = 1
      do while (m < n)
         m = 2*m
      end do

   end function ir_fft_next_pow2

!**************************************************************************
!
! In-place radix-2 Cooley-Tukey FFT, n a power of two.
!
!   isign = -1   X(k) = sum_j x(j) exp(-2 pi i j k / n)      forward, as numpy
!   isign = +1   X(k) = sum_j x(j) exp(+2 pi i j k / n)      inverse, UNSCALED
!
! The caller divides by n after an inverse transform. Written out rather than
! taken from a library because the only two transforms this module needs are a
! forward and an inverse of the same power-of-two length, and a dependency for
! that is not worth the build complexity.
!
   subroutine ir_fft_transform(z, n, isign)

      implicit none

      integer, intent(in) :: n
      integer, intent(in) :: isign
      complex(dp), intent(inout) :: z(0:n - 1)
      integer :: i, j, k, m, mmax, istep
      real(dp) :: theta, pi
      complex(dp) :: w, temp

      if (n < 2) return

!     Bit-reversal permutation.
      j = 0
      do i = 0, n - 2
         if (i < j) then
            temp = z(i)
            z(i) = z(j)
            z(j) = temp
         end if
         m = n/2
         do while (m >= 1 .and. j >= m)
            j = j - m
            m = m/2
         end do
         j = j + m
      end do

!     Danielson-Lanczos. At each stage mmax is the half-length of the
!     butterflies being combined, so the twiddle is exp(isign i pi m / mmax).
      pi = dacos(-1.d0)
      mmax = 1
      do while (n > mmax)
         istep = 2*mmax
         do m = 0, mmax - 1
            theta = dfloat(isign)*pi*dfloat(m)/dfloat(mmax)
            w = dcmplx(dcos(theta), dsin(theta))
            do i = m, n - 1, istep
               k = i + mmax
               temp = w*z(k)
               z(k) = z(i) - temp
               z(i) = z(i) + temp
            end do
         end do
         mmax = istep
      end do

   end subroutine ir_fft_transform

!**************************************************************************
!
! The dipole autocorrelation, by Wiener-Khinchin.
!
!   C(tau) = (1/T) sum_{a=0}^{T-1-tau} d(a) . d(a+tau),   d = mu - <mu>
!
! summed (not averaged) over the three Cartesian components, and divided by
! the constant T -- the biased estimator, for the reason in the header.
!
! mu is chronological. mu_mean returns the mean that was subtracted, which the
! caller wants for the gradient and for diagnostics.
!
! The zero padding is what makes the circular correlation the FFT computes
! equal to the linear one: with n >= 2T the wrap-round term at lag tau
! involves products of d(a) with d(a + tau - n), and a + tau - n < 0 for every
! a < T, so it is a product with a zero.
!
   subroutine ir_fft_autocorrelation(mu, n_frames, subtract_mean, acf, mu_mean)

      implicit none

      integer, intent(in) :: n_frames
      real(dp), intent(in) :: mu(1:3, 1:n_frames)
      logical, intent(in) :: subtract_mean
      real(dp), intent(out) :: acf(0:n_frames - 1)
      real(dp), intent(out) :: mu_mean(1:3)
      complex(dp), allocatable :: z(:)
      real(dp) :: mean, re, im
      integer :: n, t, d

      acf = 0.d0
      mu_mean = 0.d0
      if (n_frames < 1) return

      n = ir_fft_next_pow2(2*n_frames)
      allocate (z(0:n - 1))

      do d = 1, 3
         mean = 0.d0
         if (subtract_mean) then
            mean = sum(mu(d, 1:n_frames))/dfloat(n_frames)
         end if
         mu_mean(d) = mean

         z = dcmplx(0.d0, 0.d0)
         do t = 0, n_frames - 1
            z(t) = dcmplx(mu(d, t + 1) - mean, 0.d0)
         end do

         call ir_fft_transform(z, n, -1)
!        |Z|^2. Real and even, so its inverse transform is real.
         do t = 0, n - 1
            re = real(z(t), kind=dp)
            im = aimag(z(t))
            z(t) = dcmplx(re*re + im*im, 0.d0)
         end do
         call ir_fft_transform(z, n, +1)

         do t = 0, n_frames - 1
            acf(t) = acf(t) + real(z(t), kind=dp)/dfloat(n)
         end do
      end do

      acf = acf/dfloat(n_frames)

      deallocate (z)

   end subroutine ir_fft_autocorrelation

!**************************************************************************
!
! The lag window w(tau), tau = 0 .. L-1, over a half-length of L.
!
! Every one of these is 1 at tau = 0 and falls to (about) 0 at tau = L. A lag
! window that is not 1 at tau = 0 discards the total intensity; see the header
! on np.blackman.
!
! "hann" and "none" reproduce spectroscopy.py exactly. "bartlett", "welch" and
! "lorch" are mad_ir.f90's set, offered so that the two estimators can be run
! with the same taper and the remaining difference attributed to something
! else.
!
   subroutine ir_fft_lag_window(kind, L, w)

      implicit none

      integer, intent(in) :: L
      character(len=*), intent(in) :: kind
      real(dp), intent(out) :: w(0:L - 1)
      real(dp) :: pi, x
      integer :: t

      pi = dacos(-1.d0)

      select case (trim(kind))
      case ("none", "boxcar")
         w = 1.d0
      case ("blackman")
!        The HALF Blackman, not np.blackman; see the header.
         do t = 0, L - 1
            x = pi*dfloat(t)/dfloat(L)
            w(t) = 0.42d0 + 0.5d0*dcos(x) + 0.08d0*dcos(2.d0*x)
         end do
      case ("bartlett")
         do t = 0, L - 1
            w(t) = 1.d0 - dfloat(t)/dfloat(L)
         end do
      case ("welch")
         do t = 0, L - 1
            w(t) = 1.d0 - (dfloat(t)/dfloat(L))**2
         end do
      case ("lorch")
         w(0) = 1.d0
         do t = 1, L - 1
            x = pi*dfloat(t)/dfloat(L)
            w(t) = dsin(x)/x
         end do
      case default
!        Hann: 0.5*(cos(pi tau / L) + 1), exactly the Python's expression.
         do t = 0, L - 1
            w(t) = 0.5d0*(dcos(pi*dfloat(t)/dfloat(L)) + 1.d0)
         end do
      end select

   end subroutine ir_fft_lag_window

!**************************************************************************
!
! How many lags acf_ratio buys out of T frames, and what resolution that is.
! Truncating rather than rounding matches the Python's int().
!
   function ir_fft_n_lag(n_frames, acf_ratio) result(L)

      implicit none

      integer, intent(in) :: n_frames
      real(dp), intent(in) :: acf_ratio
      integer :: L

      L = int(dfloat(n_frames)*acf_ratio)
      if (L < 1) L = 1
      if (L > n_frames) L = n_frames

   end function ir_fft_n_lag

   function ir_fft_resolution(n_lag, dt_fs) result(dnu)

      implicit none

      integer, intent(in) :: n_lag
      real(dp), intent(in) :: dt_fs
      real(dp) :: dnu

      dnu = IR_FFT_CM_PER_INV_FS/(dfloat(max(n_lag, 1))*dt_fs)

   end function ir_fft_resolution

!**************************************************************************
!
! The quantum-correction prefactor P(nu), by which the power spectrum M is
! multiplied to give the absorption lineshape. See the header for what the
! choice does.
!
   subroutine ir_fft_prefactor(kind, temperature, nu, n, pref, ok, msg)

      implicit none

      integer, intent(in) :: n
      character(len=*), intent(in) :: kind
      real(dp), intent(in) :: temperature
      real(dp), intent(in) :: nu(1:n)
      real(dp), intent(out) :: pref(1:n)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: kT, x
      integer :: k

      ok = .true.
      msg = ""

      select case (trim(kind))
      case ("classical", "quadratic")
         pref(1:n) = nu(1:n)**2
      case ("linear")
         pref(1:n) = nu(1:n)
      case ("none")
         pref(1:n) = 1.d0
      case ("harmonic")
         kT = IR_FFT_KB_EV*temperature
         if (kT < 1.d-30) kT = 1.d-30
         do k = 1, n
            x = IR_FFT_HBAR_C_EV_CM*nu(k)/kT
!           1 - exp(-x) -> x as x -> 0, so this is continuous at nu = 0 and
!           reduces to the classical nu^2 form there. No guard is needed: at
!           x = 0 the expression is exactly 0.
            pref(k) = nu(k)*(1.d0 - dexp(-x))
         end do
      case default
         ok = .false.
         msg = "ir_fft: quantum_correction must be harmonic, classical "// &
               "(alias quadratic), linear or none"
      end select

   end subroutine ir_fft_prefactor

!**************************************************************************
!
! Validate a configuration before anything is allocated, so that a bad input
! is a message rather than a crash three routines down.
!
   subroutine ir_fft_check_config(cfg, n_frames, ok, msg)

      implicit none

      type(ir_fft_config_type), intent(in) :: cfg
      integer, intent(in) :: n_frames
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: pref(1:1), nu1(1:1)
      logical :: ok2
      character(len=512) :: msg2
      integer :: L

      ok = .false.
      msg = ""

      if (cfg%dt_fs <= 0.d0) then
         msg = "ir_fft: the sampling interval must be positive"
         return
      end if
      if (cfg%acf_ratio <= 0.d0 .or. cfg%acf_ratio > 1.d0) then
         msg = "ir_fft: ir_fft_acf_ratio must be in (0, 1]"
         return
      end if
      if (n_frames < 2) then
         msg = "ir_fft: at least two frames are needed"
         return
      end if
      L = ir_fft_n_lag(n_frames, cfg%acf_ratio)
      if (L < 2) then
         write (msg, '(A,I0,A,F8.4,A)') &
            "ir_fft: ", n_frames, " frames at ir_fft_acf_ratio = ", cfg%acf_ratio, &
            " leaves fewer than two lags; lengthen the run or raise the ratio"
         return
      end if
      if (cfg%max_freq_cm <= 0.d0) then
         msg = "ir_fft: ir_nu_max must be positive"
         return
      end if
!     Above Nyquist nothing is represented: power there does not go missing,
!     it folds back down into the range being reported. Same check, same
!     reason, as mad_ir_size_window.
      if (cfg%max_freq_cm > IR_FFT_CM_PER_INV_FS/(2.d0*cfg%dt_fs)) then
         write (msg, '(A,F10.1,A,F10.1,A)') &
            "ir_fft: a sampling interval of ", cfg%dt_fs, &
            " fs only reaches ", IR_FFT_CM_PER_INV_FS/(2.d0*cfg%dt_fs), &
            " cm^-1; lower ir_nu_max or sample more often"
         return
      end if
      select case (trim(cfg%smooth_kind))
      case ("gaussian", "box")
      case default
         msg = "ir_fft: ir_fft_smooth_kind must be gaussian or box"
         return
      end select
      select case (trim(cfg%window))
      case ("hann", "blackman", "none", "boxcar", "bartlett", "welch", "lorch")
      case default
         msg = "ir_fft: ir_window must be one of hann blackman none "// &
               "bartlett welch lorch"
         return
      end select
      if (trim(cfg%quantum_correction) == "harmonic" .and. cfg%temperature <= 0.d0) then
         msg = "ir_fft: the harmonic quantum correction needs a positive "// &
               "temperature; set ir_fft_temperature"
         return
      end if
      nu1(1) = 1.d0
      call ir_fft_prefactor(cfg%quantum_correction, max(cfg%temperature, 1.d0), &
                            nu1, 1, pref, ok2, msg2)
      if (.not. ok2) then
         msg = trim(msg2)
         return
      end if

      ok = .true.

   end subroutine ir_fft_check_config

!**************************************************************************
!
! Gaussian smoothing, matching scipy.ndimage.gaussian_filter1d with
! mode = "nearest" and the default truncate = 4.
!
!   sigma  = smooth_k / 2.355        FWHM in bins -> standard deviation
!   radius = int(4 sigma + 0.5)
!   w(i)   ~ exp(-i^2 / 2 sigma^2),  i = -radius .. radius, normalised to 1
!   y(j)   = sum_i w(i) x(clamp(j+i))
!
! The kernel is symmetric, so correlation and convolution coincide and the
! adjoint below only has to undo the edge clamping.
!
   subroutine ir_fft_smooth_gaussian(x, n, smooth_k, y)

      implicit none

      integer, intent(in) :: n, smooth_k
      real(dp), intent(in) :: x(1:n)
      real(dp), intent(out) :: y(1:n)
      real(dp), allocatable :: wk(:)
      real(dp) :: sigma, s
      integer :: radius, i, j, m

      sigma = dfloat(smooth_k)/2.355d0
      radius = int(4.d0*sigma + 0.5d0)
      if (radius < 1) then
         y = x
         return
      end if

      allocate (wk(-radius:radius))
      s = 0.d0
      do i = -radius, radius
         wk(i) = dexp(-0.5d0*(dfloat(i)/sigma)**2)
         s = s + wk(i)
      end do
      wk = wk/s

      do j = 1, n
         y(j) = 0.d0
         do i = -radius, radius
            m = j + i
            if (m < 1) m = 1
            if (m > n) m = n
            y(j) = y(j) + wk(i)*x(m)
         end do
      end do

      deallocate (wk)

   end subroutine ir_fft_smooth_gaussian

!  Adjoint of the above: given dL/dy, accumulate dL/dx.
   subroutine ir_fft_smooth_gaussian_adj(gy, n, smooth_k, gx)

      implicit none

      integer, intent(in) :: n, smooth_k
      real(dp), intent(in) :: gy(1:n)
      real(dp), intent(inout) :: gx(1:n)
      real(dp), allocatable :: wk(:)
      real(dp) :: sigma, s
      integer :: radius, i, j, m

      sigma = dfloat(smooth_k)/2.355d0
      radius = int(4.d0*sigma + 0.5d0)
      if (radius < 1) then
         gx = gx + gy
         return
      end if

      allocate (wk(-radius:radius))
      s = 0.d0
      do i = -radius, radius
         wk(i) = dexp(-0.5d0*(dfloat(i)/sigma)**2)
         s = s + wk(i)
      end do
      wk = wk/s

      do j = 1, n
         do i = -radius, radius
            m = j + i
            if (m < 1) m = 1
            if (m > n) m = n
            gx(m) = gx(m) + wk(i)*gy(j)
         end do
      end do

      deallocate (wk)

   end subroutine ir_fft_smooth_gaussian_adj

!**************************************************************************
!
! Box smoothing: a moving average of width k, in numpy's mode = "valid". The
! output is SHORTER by k-1 and its first bin is centred k-1 bins in half a
! width, which is why the frequency axis moves with it:
!
!   y(j)  = (1/k) sum_{m=0}^{k-1} x(j+m),      j = 1 .. n-k+1
!   nu(j) = nu_1 + (k-1)/2 d_nu + (j-1) d_nu
!
! This is GPUMD's smoother. It has sidelobes -- a boxcar in the frequency
! domain is a sinc in the lag domain -- so "gaussian" is the default here even
! though the Python offers both.
!
   subroutine ir_fft_smooth_box(x, n, smooth_k, y, n_out)

      implicit none

      integer, intent(in) :: n, smooth_k
      real(dp), intent(in) :: x(1:n)
      real(dp), intent(out) :: y(1:n)
      integer, intent(out) :: n_out
      real(dp) :: acc
      integer :: j, m

      n_out = n - smooth_k + 1
      if (n_out < 1) then
         n_out = n
         y(1:n) = x(1:n)
         return
      end if
      do j = 1, n_out
         acc = 0.d0
         do m = 0, smooth_k - 1
            acc = acc + x(j + m)
         end do
         y(j) = acc/dfloat(smooth_k)
      end do

   end subroutine ir_fft_smooth_box

   subroutine ir_fft_smooth_box_adj(gy, n_out, n, smooth_k, gx)

      implicit none

      integer, intent(in) :: n_out, n, smooth_k
      real(dp), intent(in) :: gy(1:n_out)
      real(dp), intent(inout) :: gx(1:n)
      integer :: j, m

      if (n_out == n) then
         gx(1:n) = gx(1:n) + gy(1:n)
         return
      end if
      do j = 1, n_out
         do m = 0, smooth_k - 1
            gx(j + m) = gx(j + m) + gy(j)/dfloat(smooth_k)
         end do
      end do

   end subroutine ir_fft_smooth_box_adj

!**************************************************************************
!
! Release a result.
!
   subroutine ir_fft_free(res)

      implicit none

      type(ir_fft_result_type), intent(inout) :: res

      if (allocated(res%freq)) deallocate (res%freq)
      if (allocated(res%intensity)) deallocate (res%intensity)
      if (allocated(res%intensity_raw)) deallocate (res%intensity_raw)
      if (allocated(res%power)) deallocate (res%power)
      if (allocated(res%power_raw)) deallocate (res%power_raw)
      if (allocated(res%acf)) deallocate (res%acf)
      res%n_freq = 0
      res%n_freq_raw = 0
      res%n_lag = 0
      res%n_frames = 0

   end subroutine ir_fft_free

!**************************************************************************
!
! THE PIPELINE. Steps (1) to (8) of the header, in that order.
!
! mu(1:3, 1:n_frames) is chronological. Everything else comes from cfg, and
! res is filled from scratch (any previous contents are released first).
!
   subroutine ir_fft_spectrum(mu, n_frames, cfg, res, ok, msg)

      implicit none

      integer, intent(in) :: n_frames
      real(dp), intent(in) :: mu(1:3, 1:n_frames)
      type(ir_fft_config_type), intent(in) :: cfg
      type(ir_fft_result_type), intent(inout) :: res
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: acf_full(:), w(:), a(:), M(:), nu(:), pref(:)
      real(dp), allocatable :: inten(:), pw(:), tmp(:)
      real(dp) :: mu_mean(1:3), two_pi, theta, acc, dnu, peak, peakp
      integer :: L, N2, k, t, k_max, n_out, n_box
      logical :: ok2
      character(len=512) :: msg2

      ok = .false.
      msg = ""

      call ir_fft_check_config(cfg, n_frames, ok2, msg2)
      if (.not. ok2) then
         msg = trim(msg2)
         return
      end if

      call ir_fft_free(res)

      L = ir_fft_n_lag(n_frames, cfg%acf_ratio)
      N2 = 2*L - 1
      dnu = IR_FFT_CM_PER_INV_FS/(dfloat(N2)*cfg%dt_fs)

!     (1)+(2) the mean-subtracted autocorrelation, all T lags.
      allocate (acf_full(0:n_frames - 1))
      call ir_fft_autocorrelation(mu, n_frames, cfg%subtract_mean, acf_full, mu_mean)

!     (3) keep the first L.
      allocate (res%acf(0:L - 1))
      res%acf(0:L - 1) = acf_full(0:L - 1)
      deallocate (acf_full)

!     (4) taper and Kronecker-double.
      allocate (w(0:L - 1), a(0:L - 1))
      call ir_fft_lag_window(cfg%window, L, w)
      a(0) = res%acf(0)*w(0)
      do t = 1, L - 1
         a(t) = 2.d0*res%acf(t)*w(t)
      end do

!     (6) the grid, and (8)'s cut applied first so that the transform is only
!     evaluated where it will be kept. k_max is the last bin at or below
!     max_freq_cm; bins run k = 0 .. L-1 in the full Python grid.
!     Found by scanning with the same expression that builds nu below, rather
!     than by int(max_freq_cm/dnu): the two differ in the last bit at a bin
!     that lands exactly on max_freq_cm, and numpy's freq_cm <= max_freq_cm
!     keeps that bin. Getting this wrong is a silent one-bin disagreement with
!     the Python at some grids and not others.
      k_max = -1
      do k = 0, L - 1
         if (dfloat(k)*dnu <= cfg%max_freq_cm) then
            k_max = k
         else
            exit
         end if
      end do
      if (k_max > L - 1) k_max = L - 1
      if (k_max < 1) then
         deallocate (w, a)
         call ir_fft_free(res)
         write (msg, '(A,F10.3,A,F10.1,A)') &
            "ir_fft: the bin spacing is ", dnu, " cm^-1, so no bin lands at or "// &
            "below ir_nu_max = ", cfg%max_freq_cm, "; keep more lags"
         return
      end if
      n_out = k_max + 1

      allocate (nu(1:n_out), M(1:n_out), pref(1:n_out))
      do k = 1, n_out
         nu(k) = dfloat(k - 1)*dnu
      end do

!     (5) the cosine transform.
      two_pi = 2.d0*dacos(-1.d0)
      do k = 1, n_out
         acc = 0.d0
         do t = 0, L - 1
            theta = two_pi*dfloat(k - 1)*dfloat(t)/dfloat(N2)
            acc = acc + a(t)*dcos(theta)
         end do
         M(k) = acc
      end do

!     (7) the quantum correction.
      call ir_fft_prefactor(cfg%quantum_correction, cfg%temperature, nu, n_out, &
                            pref, ok2, msg2)
      if (.not. ok2) then
         deallocate (w, a, nu, M, pref)
         msg = trim(msg2)
         return
      end if

      allocate (inten(1:n_out), pw(1:n_out))
      do k = 1, n_out
         inten(k) = pref(k)*M(k)
         pw(k) = M(k)
      end do

!     (8) smooth. Gaussian keeps the length; box shortens it and shifts the
!     axis, exactly as np.convolve(mode="valid") does.
      n_box = n_out
      if (cfg%smooth_k > 1 .and. n_out > cfg%smooth_k) then
         allocate (tmp(1:n_out))
         if (trim(cfg%smooth_kind) == "box") then
            call ir_fft_smooth_box(inten, n_out, cfg%smooth_k, tmp, n_box)
            inten(1:n_box) = tmp(1:n_box)
            call ir_fft_smooth_box(pw, n_out, cfg%smooth_k, tmp, n_box)
            pw(1:n_box) = tmp(1:n_box)
            do k = 1, n_box
               nu(k) = dfloat(k - 1)*dnu + 0.5d0*dfloat(cfg%smooth_k - 1)*dnu
            end do
         else
            call ir_fft_smooth_gaussian(inten, n_out, cfg%smooth_k, tmp)
            inten(1:n_out) = tmp(1:n_out)
            call ir_fft_smooth_gaussian(pw, n_out, cfg%smooth_k, tmp)
            pw(1:n_out) = tmp(1:n_out)
         end if
         deallocate (tmp)
      end if

      res%n_freq = n_box
      res%n_freq_raw = n_out
      res%n_lag = L
      res%n_frames = n_frames
      res%dt_fs = cfg%dt_fs
      res%d_nu = dnu
      res%resolution = ir_fft_resolution(L, cfg%dt_fs)
      res%nyquist = IR_FFT_CM_PER_INV_FS/(2.d0*cfg%dt_fs)
      res%mu_mean = mu_mean

      allocate (res%freq(1:n_box), res%intensity(1:n_box), res%intensity_raw(1:n_box))
      allocate (res%power(1:n_box), res%power_raw(1:n_box))
      res%freq(1:n_box) = nu(1:n_box)
      res%intensity_raw(1:n_box) = inten(1:n_box)
      res%power_raw(1:n_box) = pw(1:n_box)

!     Peak normalisation. The intensity is normalised over everything
!     returned; the power over the vibrational range only, because M(w) has a
!     DC peak that has nothing to do with a vibration and would otherwise set
!     the scale for the whole plot.
      peak = 0.d0
      do k = 1, n_box
         peak = max(peak, dabs(inten(k)))
      end do
      peakp = 0.d0
      do k = 1, n_box
         if (nu(k) >= cfg%power_dc_cutoff_cm) peakp = max(peakp, dabs(pw(k)))
      end do
!     If the cut-off excluded everything -- a very coarse grid, or a cut-off
!     above max_freq_cm -- fall back to the full range rather than divide by
!     zero and return a spectrum of NaNs.
      if (peakp <= 0.d0) then
         do k = 1, n_box
            peakp = max(peakp, dabs(pw(k)))
         end do
      end if
      res%peak = peak
      res%peak_power = peakp

      if (cfg%normalise .and. peak > 0.d0) then
         res%intensity(1:n_box) = inten(1:n_box)/peak
      else
         res%intensity(1:n_box) = inten(1:n_box)
      end if
      if (cfg%normalise .and. peakp > 0.d0) then
         res%power(1:n_box) = pw(1:n_box)/peakp
      else
         res%power(1:n_box) = pw(1:n_box)
      end if

      deallocate (w, a, nu, M, pref, inten, pw)
      ok = .true.

   end subroutine ir_fft_spectrum

!**************************************************************************
!
! THE MAD BIAS: the mismatch with an experiment, and its exact gradient with
! respect to the NEWEST dipole in the buffer.
!
!   I_fit(j) = interp( I(nu), nu_exp(j) )                    linear
!   Lo       = 1/2 e_s sum_j wgt_j ( s I_fit(j) + b - I_exp(j) )^2
!   lambda   = dLo / dmu(:, n_frames)
!
! with s (and b, when match_offset) fitted by weighted least squares, exactly
! as mad_ir_evaluate does. Both sit at their own optimum, so by the envelope
! theorem neither contributes a term to dLo/dI.
!
! WHY THE NEWEST FRAME AND NOTHING ELSE. The bias is a force applied now, and
! the only configuration a force can act on is the current one. mu of every
! other frame in the buffer is a number already written down: its gradient
! with respect to today's positions is zero. This is the same restriction the
! ACF bias operates under and it is not an approximation -- it is what
! "biasing a trajectory" means.
!
! THE CHAIN, backwards through the pipeline of the header:
!
!   dLo/dI_fit(j) = e_s wgt_j s ( s I_fit(j) + b - I_exp(j) )
!   dLo/dI(k)     = sum_j dLo/dI_fit(j) * (interp weight of k in j)
!   dLo/dIs(k)    = smoothing adjoint of the above
!   dLo/dM(k)     = pref(k) dLo/dIs(k)
!   dLo/da(tau)   = sum_k dLo/dM(k) cos(2 pi (k-1) tau / N2)
!   dLo/dC(tau)   = dLo/da(tau) w(tau) kron(tau)
!   lambda_c      = sum_tau dLo/dC(tau) dC(tau)/dmu_c(newest)
!
! and the last factor is where the mean subtraction earns its keep. Writing
! d(a) = mu(a) - (1/T) sum_b mu(b) with a the AGE (0 = newest), so that
! dd(a)/dmu_c(0) = delta_{a,0} - 1/T,
!
!   dC(0)/dmu_c   = [ 2 d_c(0) - P_c(0)/T ] / T
!   dC(tau)/dmu_c = [   d_c(tau) - P_c(tau)/T ] / T          tau > 0
!
!   P_c(tau) = sum_{a=0}^{T-1-tau} [ d_c(a) + d_c(a+tau) ]
!            = S_c(T-1-tau) - S_c(tau-1),   S_c(k) = sum_{a=0}^{k} d_c(a)
!
! so one prefix sum gives every P_c(tau). P_c(0) = 2 S_c(T-1) = 0 identically,
! which is worth keeping as a check on the prefix sum rather than special-
! casing away. These are the same expressions as mad_ir_evaluate's, and they
! must be: the two estimators compute the same C(tau) by different routes.
!
! NOT NORMALISED. cfg%normalise is ignored here and treated as .false.; see
! the header on why dividing by max|I| is not a thing to differentiate.
!
   subroutine ir_fft_loss(mu, n_frames, cfg, nu_exp, I_exp, wgt, n_exp, &
                          match_scale, match_offset, energy_scale, &
                          energy, lambda, I_fit, scale, offset, &
                          dissim, dissim_ref, ok, msg)

      implicit none

      integer, intent(in) :: n_frames
      integer, intent(in) :: n_exp
      real(dp), intent(in) :: mu(1:3, 1:n_frames)
      type(ir_fft_config_type), intent(in) :: cfg
      real(dp), intent(in) :: nu_exp(1:n_exp), I_exp(1:n_exp), wgt(1:n_exp)
      logical, intent(in) :: match_scale, match_offset
      real(dp), intent(in) :: energy_scale
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: lambda(1:3)
      real(dp), intent(out) :: I_fit(1:n_exp)
      real(dp), intent(out) :: scale, offset
      real(dp), intent(out) :: dissim, dissim_ref
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      type(ir_fft_config_type) :: cfg_raw
      type(ir_fft_result_type) :: res
      real(dp), allocatable :: gI(:), gIs(:), gM(:), ga(:), gC(:), pref(:)
      real(dp), allocatable :: w(:), d(:, :), S(:, :)
      real(dp) :: two_pi, theta, x, tt, resid, c, fmean
      real(dp) :: s11, s10, s00, s1y, s0y, det
      real(dp) :: P(1:3), dCdmu(1:3)
      integer :: L, N2, n_out, n_box, k, t, j, kk, a, m
      logical :: ok2
      character(len=512) :: msg2

      ok = .false.
      msg = ""
      energy = 0.d0
      lambda = 0.d0
      I_fit = 0.d0
      scale = 1.d0
      offset = 0.d0
      dissim = 0.d0
      dissim_ref = 0.d0

!     The forward pass, unnormalised. Everything the adjoint needs that is not
!     cheap to recompute comes back in res.
      cfg_raw = cfg
      cfg_raw%normalise = .false.
      call ir_fft_spectrum(mu, n_frames, cfg_raw, res, ok2, msg2)
      if (.not. ok2) then
         msg = trim(msg2)
         return
      end if

      L = res%n_lag
      N2 = 2*L - 1
      n_box = res%n_freq
!     The length before smoothing, from the forward pass rather than inferred
!     from n_box: the box smoother shortens, but only when it ran at all, and
!     it is skipped when there are fewer bins than its width.
      n_out = res%n_freq_raw

!     ---- interpolate onto the experimental grid ----------------------
!     The bins are uniform, so the bracket is arithmetic. Outside the grid the
!     value is clamped to the end bin, and so is its gradient: an experimental
!     point the run cannot resolve should pull on the nearest thing it can,
!     not on nothing.
      allocate (gI(1:n_box))
      gI = 0.d0
      do j = 1, n_exp
         x = (nu_exp(j) - res%freq(1))/res%d_nu
         if (x < 0.d0) then
            kk = -1
         else
            kk = int(x)
         end if
         if (kk < 0) then
            I_fit(j) = res%intensity(1)
         else if (kk >= n_box - 1) then
            I_fit(j) = res%intensity(n_box)
         else
            tt = x - dfloat(kk)
            I_fit(j) = (1.d0 - tt)*res%intensity(kk + 1) + tt*res%intensity(kk + 2)
         end if
      end do

!     ---- fit the scale (and offset) by weighted least squares --------
!     Identical to mad_ir_evaluate: without s the loss compares two things in
!     different units and its gradient means nothing; b exists for
!     experimental files whose baseline has been subtracted to a hard zero.
      scale = 1.d0
      offset = 0.d0
      if (match_scale) then
         s11 = 0.d0; s10 = 0.d0; s00 = 0.d0; s1y = 0.d0; s0y = 0.d0
         do j = 1, n_exp
            s11 = s11 + wgt(j)*I_fit(j)*I_fit(j)
            s10 = s10 + wgt(j)*I_fit(j)
            s00 = s00 + wgt(j)
            s1y = s1y + wgt(j)*I_fit(j)*I_exp(j)
            s0y = s0y + wgt(j)*I_exp(j)
         end do
         if (match_offset) then
            det = s11*s00 - s10*s10
            if (dabs(det) > 1.d-300) then
               scale = (s1y*s00 - s0y*s10)/det
               offset = (s11*s0y - s10*s1y)/det
            end if
         else
            if (s11 > 1.d-300) scale = s1y/s11
         end if
      end if

!     ---- the loss, and the diagnostics that are not the loss ---------
      do j = 1, n_exp
         resid = scale*I_fit(j) + offset - I_exp(j)
         dissim = dissim + wgt(j)*resid*resid
         dissim_ref = dissim_ref + wgt(j)*I_exp(j)*I_exp(j)
      end do
      energy = 0.5d0*energy_scale*dissim

      if (.not. (dabs(energy_scale) > 0.d0)) then
!        A reference run: the spectrum is wanted, the force is not, and
!        computing an adjoint that will be multiplied by zero is waste.
         call ir_fft_free(res)
         deallocate (gI)
         ok = .true.
         return
      end if

!     ---- adjoint: experimental grid -> spectrum bins -----------------
      do j = 1, n_exp
         resid = scale*I_fit(j) + offset - I_exp(j)
         c = energy_scale*wgt(j)*scale*resid
         x = (nu_exp(j) - res%freq(1))/res%d_nu
         if (x < 0.d0) then
            kk = -1
         else
            kk = int(x)
         end if
         if (kk < 0) then
            gI(1) = gI(1) + c
         else if (kk >= n_box - 1) then
            gI(n_box) = gI(n_box) + c
         else
            tt = x - dfloat(kk)
            gI(kk + 1) = gI(kk + 1) + (1.d0 - tt)*c
            gI(kk + 2) = gI(kk + 2) + tt*c
         end if
      end do

!     ---- adjoint: smoothing ------------------------------------------
      allocate (gIs(1:n_out))
      gIs = 0.d0
!     The same guard the forward pass used, in the same form, so the two cannot
!     drift apart.
      if (cfg%smooth_k > 1 .and. n_out > cfg%smooth_k) then
         if (trim(cfg%smooth_kind) == "box") then
            call ir_fft_smooth_box_adj(gI, n_box, n_out, cfg%smooth_k, gIs)
         else
            call ir_fft_smooth_gaussian_adj(gI, n_out, cfg%smooth_k, gIs)
         end if
      else
         gIs(1:n_out) = gI(1:n_box)
      end if

!     ---- adjoint: the quantum-correction prefactor -------------------
      allocate (pref(1:n_out), gM(1:n_out))
      do k = 1, n_out
         pref(k) = dfloat(k - 1)*res%d_nu
      end do
      call ir_fft_prefactor(cfg%quantum_correction, cfg%temperature, pref, n_out, &
                            gM, ok2, msg2)
      if (.not. ok2) then
         call ir_fft_free(res)
         deallocate (gI, gIs, pref, gM)
         msg = trim(msg2)
         return
      end if
      pref(1:n_out) = gM(1:n_out)
      do k = 1, n_out
         gM(k) = pref(k)*gIs(k)
      end do

!     ---- adjoint: the cosine transform -------------------------------
      allocate (ga(0:L - 1), gC(0:L - 1), w(0:L - 1))
      two_pi = 2.d0*dacos(-1.d0)
      do t = 0, L - 1
         x = 0.d0
         do k = 1, n_out
            theta = two_pi*dfloat(k - 1)*dfloat(t)/dfloat(N2)
            x = x + gM(k)*dcos(theta)
         end do
         ga(t) = x
      end do

!     ---- adjoint: taper and Kronecker doubling -----------------------
      call ir_fft_lag_window(cfg%window, L, w)
      gC(0) = ga(0)*w(0)
      do t = 1, L - 1
         gC(t) = 2.d0*ga(t)*w(t)
      end do

!     ---- adjoint: C(tau) -> the newest dipole ------------------------
!     d is indexed by AGE: d(:,0) is the newest frame, which is mu(:,n_frames).
      allocate (d(1:3, 0:n_frames - 1), S(1:3, 0:n_frames - 1))
      do a = 0, n_frames - 1
         d(1:3, a) = mu(1:3, n_frames - a) - res%mu_mean(1:3)
      end do
      S(1:3, 0) = d(1:3, 0)
      do a = 1, n_frames - 1
         S(1:3, a) = S(1:3, a - 1) + d(1:3, a)
      end do

!     fmean is 1/T when the mean was subtracted and 0 when it was not: with no
!     mean there is nothing for mu(newest) to reach the other frames through,
!     the P terms are absent, and the expression reduces to the plain one.
      if (cfg%subtract_mean) then
         fmean = 1.d0/dfloat(n_frames)
      else
         fmean = 0.d0
      end if

      lambda = 0.d0
      do t = 0, L - 1
         if (t == 0) then
!           P_c(0) = 2 S_c(T-1), which is 0 identically once the mean has been
!           subtracted. Kept as an expression rather than special-cased away
!           because it is then a check on the prefix sum rather than an
!           assumption about it.
            P(1:3) = 2.d0*S(1:3, n_frames - 1)
            dCdmu(1:3) = (2.d0*d(1:3, 0) - P(1:3)*fmean)/dfloat(n_frames)
         else
            m = n_frames - 1 - t
            if (m >= 0) then
               P(1:3) = S(1:3, m)
            else
               P(1:3) = 0.d0
            end if
            if (t - 1 >= 0) P(1:3) = P(1:3) - S(1:3, t - 1)
            dCdmu(1:3) = (d(1:3, t) - P(1:3)*fmean)/dfloat(n_frames)
         end if
         lambda(1:3) = lambda(1:3) + gC(t)*dCdmu(1:3)
      end do

      call ir_fft_free(res)
      deallocate (gI, gIs, gM, pref, ga, gC, w, d, S)
      ok = .true.

   end subroutine ir_fft_loss

end module ir_fft
