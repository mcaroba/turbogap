! Experimental IR forces for MAD, and IR spectrum prediction.
!
! The observable is an IR spectrum obtained from the autocorrelation function of
! the total dipole along a trajectory, and the MAD force is the gradient of the
! mismatch between that spectrum and an experimental one.
!
!==========================================================================
! WHERE THE FORMULA COMES FROM
!==========================================================================
!
! A weak IR field couples to the cell through H' = -M.E(t), where M is the
! total dipole. Linear response (Gordon 1965; McQuarrie ch. 21) gives the
! absorption coefficient as a Fourier transform of a dipole autocorrelation:
!
!   alpha(w) n(w) = (2 pi w / 3 hbar c V eps0) (1 - exp(-beta hbar w))
!                   INT dt exp(-i w t) <dM(0).dM(t)>_qm                    (1)
!
! Note dM = M - <M>. It is the FLUCTUATION that correlates, not the dipole.
! This is not a convention that can be dropped: a static dipole does not
! absorb, and the mean drops out of the derivation for that reason. See
! ir_subtract_mean below.
!
! (1) is quantum. A classical trajectory cannot produce <.>_qm, so the standard
! bridge is a quantum correction factor. With the HARMONIC QCF,
!
!   C_qm(w) ~= [beta hbar w / (1 - exp(-beta hbar w))] C_cl(w)             (2)
!
! the (1 - exp(-beta hbar w)) in (1) cancels exactly and one power of w is left:
!
!   alpha(w) n(w) ~= (2 pi beta w^2 / 3 c V eps0)
!                    INT dt exp(-i w t) <dM(0).dM(t)>_cl                   (3)
!
! So: nu^2 times the cosine transform of the classical dipole ACF. That is
! ir_nu_power = 2, and Ramirez et al. (JCP 121, 3973 (2004)) compared the
! available QCFs for IR lineshapes and found the harmonic one best. Leave it
! at 2 unless you know why you are changing it.
!
! Because C is even, INT_{-inf}^{inf} C cos = 2 INT_0^inf C cos, which is where
! the "C(0) + 2*sum" shape of the discrete transform below comes from.
!
! WHAT THE INTENSITY IS. (3) produces alpha*n, the absorption coefficient times
! the refractive index. That is what the experimental file must also be, and it
! is NOT the same as alpha alone or as the extinction coefficient k: for water n
! runs from ~1.45 in the librational region to ~1.16 across the O-H stretch, so
! the choice moves relative band intensities by tens of percent.
!
!==========================================================================
! WHAT A LAG IS, AND WHY THE LONGEST ONE IS THE PARAMETER THAT MATTERS
!==========================================================================
!
! This is the single most misunderstood part of the method, so it is worth
! setting out from the beginning.
!
! We have a sequence of stored total dipoles, one every dt fs:
!
!     mu(0), mu(1), mu(2), ... mu(N-1)          N = n_window
!
! A LAG, tau, is a separation in time measured in stored frames. Lag tau = 5 at
! dt = 2 fs means "10 fs apart". The autocorrelation at lag tau asks a single
! question:
!
!     if I know the dipole now, how much does that tell me about the dipole
!     tau frames from now?
!
! and answers it by averaging the dot product over every pair in the buffer
! that is tau apart:
!
!     C(tau) = (1/M) SUM_a mu(a) . mu(a+tau)
!
! There are N - tau such pairs. That fact drives everything below.
!
! WHY THE TRANSFORM OF C IS A SPECTRUM. If the dipole contains an oscillation
! at frequency w, then mu(a) and mu(a+tau) stay in phase whenever w*tau*dt is a
! multiple of 2 pi and go out of phase in between, so C(tau) itself oscillates
! at w. Every vibration in the system plants a cosine of its own frequency in
! C(tau), and the cosine transform reads them back out. C(tau) is the spectrum,
! written in the time domain. This is the Wiener-Khinchin theorem.
!
! CONSEQUENCE 1: THE LONGEST LAG SETS THE RESOLUTION, NOT THE RUN LENGTH.
! To tell two frequencies apart you must watch long enough for them to drift
! out of phase with each other. Two cosines differing by dnu drift apart in a
! time 1/(c dnu), so the longest lag you keep must be at least that long:
!
!     dnu = CM_PER_INV_FS / (n_lag * dt)
!
! A 40 ps run truncated at n_lag*dt = 800 fs has 40 cm^-1 resolution, not the
! 0.8 cm^-1 its length suggests. Everything past the longest lag is discarded
! before the transform ever sees it. If you want finer resolution you must
! raise n_lag -- and then, by consequence 2, lengthen the run to support it.
!
! CONSEQUENCE 2: LONG LAGS ARE NOISY, AND UNAVOIDABLY SO.
! C(tau) is an average over N - tau pairs. At tau = 0 that is N pairs; at
! tau = N-1 it is one. The variance of the estimate therefore grows without
! limit as tau approaches N, and the part of C(tau) that carries the fine
! structure is the part that is least well determined. This is why n_window is
! lag_factor * n_lag rather than n_lag: with lag_factor = 2 the longest lag is
! still averaged over N/2 pairs and the variance inflation is capped at 2x.
!
! CONSEQUENCE 3: TRUNCATING C(tau) IS ITSELF A SPECTRAL OPERATION.
! Keeping lags 0..n_lag and discarding the rest is multiplication of C by a
! boxcar, hence CONVOLUTION of the spectrum with that boxcar's transform. That
! transform is a sinc, and a sinc has large oscillating side lobes -- so a hard
! truncation makes every band ring, and an absorption spectrum that rings goes
! NEGATIVE. That is what the lag window is for, and it is the subject of the
! next block.
!
!==========================================================================
! THE LAG WINDOW: WHY HANN, AND WHY THE XRD-STYLE SINC IS ALSO OFFERED
!==========================================================================
!
! ir_window chooses w(tau), applied to C(tau) before transforming. Setting it
! to "none" is NOT "no window" -- it is a boxcar window, and the sinc ringing
! is its spectral kernel. There is no way to avoid choosing one; the only
! question is which.
!
! The trade is always the same: a taper that falls to zero smoothly has small
! side lobes (little ringing, little leakage between bands) but a wider main
! lobe (coarser effective resolution). Numbers, for a lag window of half-length
! L = 417, kernel measured directly by ana/formalism_windows.py:
!
!   ir_window   FWHM/dnu   peak side lobe   most negative   kernel >= 0 ?
!   ---------   --------   --------------   -------------   -------------
!   none          0.60        -13.3 dB         -21.7 %          no
!   welch         0.80        -21.3 dB          -8.6 %          no
!   lorch         0.87        -26.4 dB          -4.8 %          no
!   bartlett      0.89        -26.5 dB          -0.0 %          YES
!   hann          1.00        -31.5 dB          -2.7 %          no
!
!   dnu = CM_PER_INV_FS/(L dt), the nominal resolution. "most negative" is the
!   deepest excursion of the kernel below zero, i.e. how hard it rings.
!
! and the leakage that actually matters here -- a band 40x stronger sitting
! 15-25 resolution elements away, which is the O-H stretch seen from the
! transparency minimum near 2460 cm^-1:
!
!   none      worst 41.84 %   typical 23.31 %   of the weak feature
!   bartlett  worst  1.70 %   typical  0.69 %
!   welch     worst  1.36 %   typical  0.60 %
!   lorch     worst  0.77 %   typical  0.34 %
!   hann      worst  0.04 %   typical  0.02 %
!
! WHY HANN IS THE DEFAULT. Its far-field leakage is 13 dB below the next best
! option, because its kernel decays as 1/f^3 where the others go as 1/f or
! 1/f^2. For a spectrum whose bands span more than an order of magnitude --
! water's O-H stretch is ~40x the transparency window at 2500 cm^-1 -- that is
! what decides whether the weak region is a measurement or is spillover from
! the strong band next to it. With no window at all, the stretch band alone
! puts up to 42% of the weak feature's own size into that window; with Hann,
! 0.04%. Hann is also the only one whose main lobe FWHM equals the advertised
! dnu exactly, so mad_ir_size_window's promise and the delivered resolution
! agree with no fudge factor.
!
! IF YOU NEED GUARANTEED POSITIVITY, USE BARTLETT. The triangular lag window is
! the one case where the kernel is the Fejer kernel, which is non-negative
! everywhere (it is |D_L|^2/L, a squared magnitude). A non-negative kernel
! convolved with a non-negative spectrum cannot produce a negative one, so
! ir_window = "bartlett" with ir_estimator = "biased" makes a non-negative
! spectrum a THEOREM rather than an empirical observation.
!
! Note that the estimator alone is NOT enough for that: the periodogram
! argument applies before windowing, and measured on one trajectory with every
! other fix on, ir_window = "none" still gave 124 negative points of 400 while
! hann, lorch and bartlett all gave zero. The price of bartlett is 5 dB more
! leakage than Hann. Hann's residual ringing is small enough (-2.7% kernel dip,
! and no negatives observed) that it is the better default, but if something
! downstream will divide by the spectrum or take its log, take the theorem.
!
! WHY NOT JUST TRANSFORM, AS THE XRD PATH DOES. The XRD/structure-factor code
! in exp_utils.f90 applies sinc(pi r / r_cut) to g(r) before transforming --
! the Lorch modification function -- under structure_factor_window. That is the
! same operation in the conjugate pair (r,q) that this is in (tau,nu), and it
! is offered here as ir_window = "lorch" for exactly that reason: if you want
! the IR path to make the same approximation the XRD path makes, you can.
!
! Lorch is a reasonable choice and is the diffraction community's standard, but
! it is not the better one here. It buys 13% narrower bands for 5 dB more
! leakage and 1.8x more ringing. Diffraction favours it because g(r) is sharply
! peaked, resolution is at a premium, and the dynamic range is modest; an IR
! spectrum is the opposite case on all three counts.
!
! "none" is provided so the ringing can be seen rather than argued about. Run
! the same trajectory with ir_window = none and ir_window = hann and compare
! the number of negative intensities. It is not a supported production setting.
!
!==========================================================================
! WHAT SETS THE SIZE OF THE ENSEMBLE
!==========================================================================
!
! Three separate numbers, and they are set by three different requirements. Get
! them confused and the spectrum is quietly wrong rather than obviously wrong.
!
!   dt        the interval between STORED configurations, not the MD timestep.
!             It sets the Nyquist limit: nothing above
!             CM_PER_INV_FS/(2 dt) cm^-1 can be represented, and power above it
!             folds back down into the range you are fitting. So dt is fixed by
!             the HIGHEST wavenumber wanted.
!
!   n_lag     the longest correlation lag entering the cosine transform. It,
!             not the ensemble length, sets the frequency RESOLUTION:
!             d(nu) = CM_PER_INV_FS/(n_lag dt). Wanting to resolve two bands
!             4 cm^-1 apart with 1 fs sampling means n_lag ~ 8300.
!
!   n_window  how many configurations are kept. C(tau) is averaged over
!             n_window - tau pairs, so the ensemble has to be longer than the
!             longest lag or the tail of the ACF is built from a handful of
!             samples and is pure noise. n_window = lag_factor * n_lag with
!             lag_factor >= 2 is the usual compromise.
!
! mad_ir_size_window turns (dt, resolution, highest wavenumber) into n_lag and
! n_window and refuses combinations that cannot work.
!
!==========================================================================
! THE ESTIMATOR, AND THE THREE KNOBS THAT CONTROL IT
!==========================================================================
!
! ir_subtract_mean (default .true.)
!   Correlate mu - <mu> rather than mu, as equation (1) requires. The defence
!   for not doing so -- "a constant only puts a delta at nu = 0, which nu^2
!   kills" -- is false. Writing mu = mubar + d,
!
!       C(tau) = |mubar|^2 + cross(tau) + C_d(tau)
!
!   only the FIRST term is constant. cross(tau) is a partial sum of a mean-zero
!   series divided by its own length, i.e. a random walk in tau, whose amplitude
!   GROWS as the lag lengthens and the averaging shortens. Measured on a 2 fs
!   liquid-water buffer (834 frames, 100 molecules): |mubar|^2 was 86% of C(0),
!   and past ~300 fs of lag the rms cross term was 17x the rms of the real
!   correlation. It is broadband noise deposited preferentially into the fine
!   structure, not a nu = 0 artefact. Turning this off reproduces the old
!   behaviour and is for reproducing old results only.
!
! ir_estimator (default "biased")
!   "unbiased" divides C(tau) by N - tau, the number of pairs actually averaged.
!   "biased" divides by N always. The unbiased one is unbiased and is also not
!   a positive semi-definite sequence, so its cosine transform can be negative
!   -- and an absorption spectrum cannot be. The biased one IS a positive
!   semi-definite sequence: its transform is a periodogram (Percival & Walden
!   ch. 6; Numerical Recipes 13.4), at the price of scaling C(tau) by
!   (1 - tau/N). Since the lag window already tapers far harder than that,
!   nothing is lost. This is what simulation practice does and it is the
!   default here.
!
!   BE PRECISE ABOUT WHAT THIS GUARANTEES. The periodogram argument covers the
!   UNWINDOWED estimate. Windowing convolves the spectrum with the window's
!   kernel, so a kernel that goes negative can still drive the result negative
!   -- and measurably does: with ir_estimator = biased and every other fix on,
!   ir_window = "none" still produced 124 negative points out of 400, because
!   the boxcar kernel dips to -21.7%. The guarantee is only unconditional when
!   the kernel is non-negative too, which of the shapes offered here means
!   ir_window = "bartlett" (Fejer). Hann's kernel dips to -2.7% and in practice
!   gives no negatives at all (0 of 400 on the same trajectory), but that is an
!   observation, not a theorem.
!
! ir_taper_partial (default .true.)
!   Only matters while the buffer is still filling, i.e. during a do_ir
!   prediction run. The window is built for a half-length of n_lag, but early on
!   only N-1 < n_lag lags exist, so C would be truncated at a point where the
!   explicit window is still near 1 -- a hard boxcar cut, with all the ringing
!   of consequence 3 above. With this on, the taper is rebuilt over the lag
!   range actually available.
!
!   ITS EFFECT DEPENDS ON ir_estimator, WHICH IS NOT OBVIOUS AND IS WORTH
!   STATING. Dividing by N rather than N - tau multiplies C(tau) by
!   (1 - tau/N), and that IS a triangular taper -- the biased estimator applies
!   an implicit Bartlett lag window on top of the explicit one. In the
!   partial-fill regime the truncation sits at tau = N-1, where the implicit
!   taper has fallen to 1/N, so it has already gone to zero by itself and no
!   ringing appears. Measured on one buffer with the window still sized for
!   n_lag = 417 and only N = 51 frames stored, i.e. the explicit window still at
!   0.965 where C is cut:
!
!       ir_estimator = unbiased : 160 negative points of 400
!       ir_estimator = biased   :   0
!
!   So with the default biased estimator this switch does not change the
!   ringing. What it still buys is an honest header -- the resolution quoted is
!   the one the transformed lags actually support -- and correct behaviour for
!   anyone who selects the unbiased estimator, where the ringing is severe.
!   The 157-of-400 first block that motivated this fix came from the original
!   code, which had the unbiased estimator AND no partial taper.
!
!==========================================================================
! FITTING TO THE EXPERIMENT
!==========================================================================
!
! ir_match_scale (default .true.)
!   The transform's output is in arbitrary units, so a scale must be fitted or
!   the loss compares two things on different axes and its gradient is
!   meaningless.
!
! ir_match_offset (default .false.)
!   Fit an additive baseline b alongside the scale s, i.e. compare
!   s*I_calc + b with I_exp. Needed when the experimental file has had a
!   baseline removed such that its minimum is exactly zero: the model spectrum
!   is strictly positive, so without an offset the residual in the transparency
!   window can never vanish, and a quadratic loss responds by shrinking the
!   whole prediction -- fighting the band intensity the bias is trying to build.
!   Both s and b sit at their own least-squares optimum, so by the envelope
!   theorem dL/dI picks up no extra term from either.
!
!   The two-parameter solve is unconstrained and will return a NEGATIVE s when
!   the predicted shape does not resemble the experiment, which on this
!   potential it does (s = -9.7e-9, b = +2.70 measured). dL/dI carries a factor
!   of s, so that reverses the bias. mad_ir_evaluate detects it and falls back
!   to the scale-only fit; see the guard there.
!
! ir_weight_by_spacing (default .true.)
!   The fit grid IS the experimental grid (interpolating the experiment would
!   invent structure between its points and then fit to it), but experimental
!   grids are rarely uniform. The Downing & Williams water data shipped with the
!   test has 23.9 cm^-1 spacing across the O-H stretch and ~55 cm^-1 elsewhere,
!   so 45% of an unweighted loss falls in the top quarter of the range -- a
!   weighting chosen by whoever digitised the paper, not by the physics. Setting
!   wgt_k proportional to the local grid spacing turns the sum into an integral
!   over frequency and removes that dependence. Weights are normalised to mean
!   1 so that exp_energy_scales keeps its magnitude.
!
!==========================================================================
! WHERE THE FORCE ACTS, AND WHAT IT IS
!==========================================================================
!
! The spectrum is a functional of the whole stored trajectory, but only the
! newest configuration can still be moved. So the loss is differentiated with
! respect to mu of the newest frame alone,
!
!     lambda_a = dL/dmu_a(newest),
!
! and the force on atom j is
!
!     f_jb = - sum_a lambda_a  dmu_a/dr_jb,
!
! with dmu/dr from the SOAP dipole machinery (gap.f90's accumulate_dmu_dr).
! The older frames enter only through their stored dipole values, which are
! constants here.
!
! READ THIS BEFORE INTERPRETING A BIASED RUN. Unwinding the chain rule,
!
!     lambda = SUM_tau S(tau) mu(t - tau) / M_tau
!
! so the force is a CONVOLUTION of the past dipole trajectory with a kernel
! S(tau). That is the definition of a filtered, coloured external field, and
! transforming S for a water run against the Downing target puts its peak at
! ~3510 cm^-1 -- next to the 3416 cm^-1 experimental band the model is short of.
! The bias is an optical pumping field tuned to the deficient band. Three things
! follow, and none of them are bugs:
!
!   1. It does net work on the system with no dissipative counterpart, so there
!      is no fluctuation-dissipation balance and the run heats.
!   2. A global thermostat (Bussi) removes that energy democratically over
!      3N-3 degrees of freedom, so the steady state is hot where the bias pumps
!      and cold everywhere else. A temperature quoted from it is not a
!      temperature in the usual sense.
!   3. Biasing an ENSEMBLE AVERAGE (pdf, sf, xrd) has a maximum-entropy
!      justification: a least-biased reweighted ensemble exists. No such
!      theorem covers a TIME CORRELATION FUNCTION. There is no reweighting that
!      turns a biased trajectory back into a physical one.
!
! So: this is not the gradient of any potential, it does not conserve energy,
! and it is not the gradient of anything at all, because the target is not a
! function of the state. Report the output of a biased run as a driven
! non-equilibrium trajectory, not as a corrected ensemble. Keep
! exp_energy_scales small and check what else moved.
!
! One further reporting note. turbogap.f90 folds the IR mismatch into the total
! energy when exp_energies is set, but the force above is the frozen-past
! PARTIAL gradient -- the stored dipoles are functions of earlier positions and
! are treated as constants. The reported total energy is therefore not the
! generator of the reported forces, and its drift is NOT an integration
! accuracy diagnostic in a biased run. Use exp_energies = .false. if you want
! the conservation check to keep its usual meaning.
!
! Because lambda is a single 3-vector for the whole cell, and dmu/dr is
! accumulated per atom rather than per pair, the entire per-step cost beyond
! the descriptor Hessian is one cosine transform and a 9*n_atoms contraction.
!
module mad_ir

   use kinds

   implicit none

   private
   public :: mad_ir_type, mad_ir_size_window, mad_ir_init, mad_ir_free
   public :: mad_ir_push, mad_ir_ready, mad_ir_evaluate, mad_ir_forces
   public :: mad_ir_spectrum
   public :: mad_ir_select_range, mad_ir_save, mad_ir_load
   public :: mad_ir_setup, mad_ir_setup_predict, mad_ir_state
   public :: mad_ir_dmu_dr, mad_ir_collect
   public :: mad_ir_need_dmu
   public :: mad_ir_append_spectrum, mad_ir_write_exp_spectrum
   public :: mad_ir_write_spectrum, mad_ir_build_window
   public :: CM_PER_INV_FS

!  Wavenumber in cm^-1 of a frequency of 1/fs: 1/(c) with c in cm/fs.
   real(dp), parameter :: CM_PER_INV_FS = 33356.40952d0

!  ---------------------------------------------------------------------------
!  Run-wide state.
!
!  It lives here rather than being threaded through get_gap_soap's already very
!  long argument list. mad_ir_dmu_dr is accumulated across descriptor batches
!  and across MPI ranks during the force pass, and only contracted with lambda
!  once the total dipole is known -- which is why lambda does not have to exist
!  before the descriptors are evaluated.
!  ---------------------------------------------------------------------------
   real(dp), allocatable, save :: mad_ir_dmu_dr(:, :, :)   ! (3, 3, n_atoms)
   logical, save :: mad_ir_collect = .false.               ! gate in gap_interface
!  Is dmu/dr wanted at all? The spectrum needs only the dipole, which the
!  dipole model produces anyway; dmu/dr is needed solely to turn lambda into a
!  force. A run with exp_forces off is a reference trajectory -- the spectrum
!  is predicted and written but nothing is biased -- and it should cost what
!  plain MD costs, not plain MD plus a descriptor Hessian per step. Set once
!  from exp_forces; gap_interface reads it alongside mad_ir_collect.
   logical, save :: mad_ir_need_dmu = .true.

   type :: mad_ir_type
      logical :: active = .false.
!     sizing
      integer :: n_window = 0        ! stored configurations
      integer :: n_lag = 0           ! longest lag used in the transform
      integer :: n_stored = 0        ! how many have been pushed (caps at n_window)
      integer :: head = 0            ! circular slot holding the newest frame
      real(dp) :: dt = 0.d0          ! fs between stored frames
!     target spectrum
      integer :: n_freq = 0
      real(dp) :: nu_power = 2.d0    ! I(nu) = nu**nu_power * FT[C]
      logical :: match_scale = .true.
      real(dp) :: scale = 1.d0       ! the scale fitted on the last evaluate
      logical :: match_offset = .false.
      real(dp) :: offset = 0.d0      ! the baseline fitted on the last evaluate
!     estimator choices; see the header block for what each one costs and buys
      logical :: subtract_mean = .true.    ! correlate mu - <mu>, as (1) requires
      logical :: biased_estimator = .true. ! divide by N, not N - tau
      logical :: taper_partial = .true.    ! rebuild the taper while filling
      character(len=32) :: window_kind = "hann"
!     n_lag is what the run was SIZED for; n_lag_used is how many lags the
!     buffer can currently supply, and they differ only while a do_ir ensemble
!     is still filling. n_lag_win records the half-length win() is currently
!     built for, so the taper is rebuilt when it changes and not otherwise.
      integer :: n_lag_used = 0
      integer :: n_lag_win = -1
      real(dp) :: mu_mean(1:3) = 0.d0
!     arrays
      real(dp), allocatable :: mu_hist(:, :)   ! (1:3, 1:n_window)
      real(dp), allocatable :: nu(:)           ! (1:n_freq), cm^-1
      real(dp), allocatable :: I_exp(:)        ! (1:n_freq)
      real(dp), allocatable :: wgt(:)          ! (1:n_freq)
      real(dp), allocatable :: win(:)          ! (0:n_lag), lag window
      real(dp), allocatable :: acf(:)          ! (0:n_lag)
      real(dp), allocatable :: I_calc(:)       ! (1:n_freq)
   end type mad_ir_type

   type(mad_ir_type), save :: mad_ir_state

contains

!**************************************************************************
!
! Turn the spectral requirements into buffer sizes, and refuse the ones that
! cannot be met. dt is in fs and is the interval between stored frames, which
! is the MD timestep times whatever stride the caller uses.
!
   subroutine mad_ir_size_window(dt, nu_res, nu_max, lag_factor, n_lag, n_window, ok, msg)

      implicit none

      real(dp), intent(in) :: dt
      real(dp), intent(in) :: nu_res
      real(dp), intent(in) :: nu_max
      integer, intent(in) :: lag_factor
      integer, intent(out) :: n_lag
      integer, intent(out) :: n_window
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: nyquist

      ok = .false.
      n_lag = 0
      n_window = 0
      msg = ""

      if (dt <= 0.d0) then
         msg = "mad_ir: the sampling interval must be positive"
         return
      end if
      if (nu_res <= 0.d0) then
         msg = "mad_ir: the requested resolution must be positive"
         return
      end if
      if (lag_factor < 1) then
         msg = "mad_ir: lag_factor must be at least 1"
         return
      end if

!     Nyquist. Anything above this folds back into the fitted range instead of
!     being absent from it, so it is an error rather than a warning.
      nyquist = CM_PER_INV_FS/(2.d0*dt)
      if (nu_max > nyquist) then
         write (msg, '(A,F10.1,A,F10.1,A)') &
            "mad_ir: the sampling interval only reaches ", nyquist, &
            " cm^-1 but ", nu_max, " cm^-1 was asked for; sample more often"
         return
      end if

!     Resolution is set by the longest lag, not by the ensemble length.
      n_lag = ceiling(CM_PER_INV_FS/(nu_res*dt))
      if (n_lag < 1) n_lag = 1
      n_window = lag_factor*n_lag
      ok = .true.

   end subroutine mad_ir_size_window

!**************************************************************************
!
! Build the lag window over a half-length of L, i.e. w(0) = 1 falling to
! w(L) = 0, and zero beyond L so that lags outside the range contribute
! nothing. See the header block for what each shape costs and buys.
!
! Split out from mad_ir_init because L is not fixed for the life of the run:
! while a do_ir ensemble is filling, only n_stored - 1 lags exist, and tapering
! over n_lag instead would cut C(tau) off at a point where the window is still
! near 1 -- a boxcar, with all of its ringing. Rebuilt only when L changes, so
! a bias run (whose buffer is always full) builds it once.
!
   subroutine mad_ir_build_window(this, L)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      integer, intent(in) :: L
      real(dp) :: pi, x
      integer :: t

      pi = dacos(-1.d0)
      this%win = 0.d0

      select case (trim(this%window_kind))
      case ("none")
!        A boxcar. Not "no window": its kernel is a sinc, side lobes at
!        -13.3 dB, and it rings to -21.7%. Provided for comparison.
         do t = 0, L
            this%win(t) = 1.d0
         end do
      case ("lorch")
!        sinc(pi t / L). The Lorch modification function, i.e. the same taper
!        exp_utils.f90 applies to g(r) under structure_factor_window.
         this%win(0) = 1.d0
         do t = 1, L
            x = pi*dfloat(t)/dfloat(L)
            this%win(t) = dsin(x)/x
         end do
      case ("bartlett")
!        Triangular. The one shape whose kernel (Fejer) is non-negative
!        everywhere: it is |D_L|^2/L, a squared magnitude. With the biased
!        estimator -- itself a positive semi-definite sequence -- a
!        non-negative spectrum is then a theorem rather than an observation,
!        because a non-negative kernel convolved with a non-negative spectrum
!        cannot go negative.
         do t = 0, L
            this%win(t) = 1.d0 - dfloat(t)/dfloat(L)
         end do
      case ("welch")
         do t = 0, L
            this%win(t) = 1.d0 - (dfloat(t)/dfloat(L))**2
         end do
      case default
!        Hann. Lowest far-field leakage of the set and the only one whose main
!        lobe FWHM equals the advertised resolution exactly.
         do t = 0, L
            this%win(t) = 0.5d0*(1.d0 + dcos(pi*dfloat(t)/dfloat(L)))
         end do
      end select

      this%n_lag_win = L

   end subroutine mad_ir_build_window

!**************************************************************************
!
! window_kind: see mad_ir_build_window and the header. A lag window is not
! cosmetic here -- C(tau) at tau near n_lag is averaged over very few pairs,
! and transforming it unwindowed puts that noise straight into the spectrum
! and hence into the force.
!
   subroutine mad_ir_init(this, dt, n_lag, n_window, nu, I_exp, wgt, match_scale, nu_power, &
                          window_kind, subtract_mean, biased_estimator, taper_partial, &
                          match_offset)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      real(dp), intent(in) :: dt
      integer, intent(in) :: n_lag
      integer, intent(in) :: n_window
      real(dp), intent(in) :: nu(:)
      real(dp), intent(in) :: I_exp(:)
      real(dp), intent(in) :: wgt(:)
      logical, intent(in) :: match_scale
      real(dp), intent(in) :: nu_power
      character(len=*), intent(in) :: window_kind
      logical, intent(in) :: subtract_mean, biased_estimator, taper_partial
      logical, intent(in) :: match_offset

      call mad_ir_free(this)

      this%dt = dt
      this%n_lag = n_lag
      this%n_window = n_window
      this%n_stored = 0
      this%head = 0
      this%n_freq = size(nu)
      this%match_scale = match_scale
      this%nu_power = nu_power
      this%scale = 1.d0
      this%offset = 0.d0
      this%match_offset = match_offset .and. match_scale
      this%subtract_mean = subtract_mean
      this%biased_estimator = biased_estimator
      this%taper_partial = taper_partial
      this%window_kind = window_kind
      this%n_lag_used = 0
      this%n_lag_win = -1
      this%mu_mean = 0.d0

      allocate (this%mu_hist(1:3, 1:n_window))
      allocate (this%nu(1:this%n_freq), this%I_exp(1:this%n_freq), this%wgt(1:this%n_freq))
      allocate (this%I_calc(1:this%n_freq))
      allocate (this%win(0:n_lag), this%acf(0:n_lag))

      this%mu_hist = 0.d0
      this%nu = nu
      this%I_exp = I_exp
      this%wgt = wgt
      this%I_calc = 0.d0
      this%acf = 0.d0

!     Built for the full half-length here; mad_ir_spectrum rebuilds it for a
!     shorter one while the ensemble is still filling, if taper_partial is set.
      call mad_ir_build_window(this, n_lag)

      this%active = .true.

   end subroutine mad_ir_init

!**************************************************************************
   subroutine mad_ir_free(this)
      implicit none
      type(mad_ir_type), intent(inout) :: this
      if (allocated(this%mu_hist)) deallocate (this%mu_hist)
      if (allocated(this%nu)) deallocate (this%nu)
      if (allocated(this%I_exp)) deallocate (this%I_exp)
      if (allocated(this%wgt)) deallocate (this%wgt)
      if (allocated(this%win)) deallocate (this%win)
      if (allocated(this%acf)) deallocate (this%acf)
      if (allocated(this%I_calc)) deallocate (this%I_calc)
      this%active = .false.
      this%n_stored = 0
      this%head = 0
      this%n_lag_win = -1
   end subroutine mad_ir_free

!**************************************************************************
!
! Store the current total dipole. The buffer is circular, so the newest frame
! is at head and the frame of age a is at head - a wrapped into range.
!
   subroutine mad_ir_push(this, mu)
      implicit none
      type(mad_ir_type), intent(inout) :: this
      real(dp), intent(in) :: mu(1:3)
      this%head = this%head + 1
      if (this%head > this%n_window) this%head = 1
      this%mu_hist(1:3, this%head) = mu(1:3)
      if (this%n_stored < this%n_window) this%n_stored = this%n_stored + 1
   end subroutine mad_ir_push

!**************************************************************************
!
! No spectrum until the ensemble is full. Producing one from a partly filled
! buffer would give a resolution that changes from step to step, and a force
! that reflects the fill state rather than the structure.
!
   logical function mad_ir_ready(this)
      implicit none
      type(mad_ir_type), intent(in) :: this
      mad_ir_ready = this%active .and. (this%n_stored >= this%n_window)
   end function mad_ir_ready

!**************************************************************************
!
! The autocorrelation of the stored dipoles and its cosine transform. This is
! the whole of the prediction; everything mad_ir_evaluate does beyond it --
! the fitted scale and offset, the loss, lambda -- exists only because there is
! an experiment to compare against.
!
!   d(a)    = mu(a) - <mu>                              a = age, 0 is newest
!   C(tau)  = 1/M_tau sum_a d(a) . d(a+tau)             M_tau = N or N - tau
!   I(nu_k) = nu_k^p dt [ win(0)C(0) + 2 sum_tau win(tau) C(tau) cos(th_k tau) ]
!
! Three estimator choices enter here and each is documented at length in the
! header block:
!
!   subtract_mean     whether d = mu - <mu> or just mu. Equation (1) of the
!                     header wants the fluctuation. Off reproduces the old
!                     behaviour, in which a static offset that was 86% of C(0)
!                     on a real water buffer leaked broadband noise into the
!                     long-lag end.
!   biased_estimator  M_tau = N (guaranteed non-negative transform, since the
!                     result is a periodogram) or M_tau = N - tau (unbiased,
!                     but not a positive semi-definite sequence).
!   taper_partial     whether the lag window is rebuilt for the lags actually
!                     available while the ensemble is still filling, or left at
!                     the full half-length so that C is cut off with a boxcar.
!
! Split out so that do_ir prediction runs, which have no experiment and want
! the spectrum only at the end, share the arithmetic with the bias rather than
! reimplementing it. Two cosine transforms that must agree is the shape of
! defect this file's header warns about.
!
   subroutine mad_ir_spectrum(this)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      real(dp) :: two_pi, theta, acc, denom
      real(dp) :: mub(1:3)
      integer :: t, a, k, m, m2, N, L

      if (.not. this%active) return
      if (this%n_stored < 1) return

!     A partly filled buffer averages over what it actually holds, not over
!     n_window: in a prediction run the ensemble is the whole trajectory and
!     is full only at the last step, and dividing by n_window before then
!     would scale C(tau) by the fill fraction.
      N = this%n_stored
      two_pi = 2.d0*dacos(-1.d0)

!     ---- the mean, which is what the correlation is taken about ------
      mub = 0.d0
      if (this%subtract_mean) then
         do a = 0, N - 1
            m = this%head - a
            if (m < 1) m = m + this%n_window
            mub(1:3) = mub(1:3) + this%mu_hist(1:3, m)
         end do
         mub(1:3) = mub(1:3)/dfloat(N)
      end if
      this%mu_mean(1:3) = mub(1:3)

!     ---- how many lags exist, and the taper that fits them -----------
      L = min(this%n_lag, N - 1)
      if (L < 1) L = 1
      this%n_lag_used = L
      if (this%taper_partial) then
         if (this%n_lag_win /= L) call mad_ir_build_window(this, L)
      else
         if (this%n_lag_win /= this%n_lag) call mad_ir_build_window(this, this%n_lag)
      end if

      do t = 0, L
         acc = 0.d0
         do a = 0, N - 1 - t
            m = this%head - a
            if (m < 1) m = m + this%n_window
            m2 = this%head - a - t
            if (m2 < 1) m2 = m2 + this%n_window
            acc = acc + (this%mu_hist(1, m) - mub(1))*(this%mu_hist(1, m2) - mub(1)) &
                  + (this%mu_hist(2, m) - mub(2))*(this%mu_hist(2, m2) - mub(2)) &
                  + (this%mu_hist(3, m) - mub(3))*(this%mu_hist(3, m2) - mub(3))
         end do
         if (this%biased_estimator) then
            denom = dfloat(N)
         else
            denom = dfloat(N - t)
         end if
         this%acf(t) = acc/denom
      end do
!     Lags the buffer cannot reach yet contribute nothing rather than stale
!     values from a previous call.
      do t = L + 1, this%n_lag
         this%acf(t) = 0.d0
      end do

      do k = 1, this%n_freq
         acc = this%win(0)*this%acf(0)
         do t = 1, L
            theta = two_pi*this%nu(k)*dfloat(t)*this%dt/CM_PER_INV_FS
            acc = acc + 2.d0*this%win(t)*this%acf(t)*dcos(theta)
         end do
         this%I_calc(k) = (this%nu(k)**this%nu_power)*this%dt*acc
      end do

   end subroutine mad_ir_spectrum

!**************************************************************************
!
! Spectrum, loss, and the sensitivity of the loss to the newest dipole.
!
!   d(a)    = mu(a) - <mu>                              a = age, 0 is newest
!   C(tau)  = 1/M_tau sum_a d(a) . d(a+tau)             M_tau = N or N - tau
!   I(nu_k) = nu_k^p dt [ win(0)C(0) + 2 sum_tau win(tau) C(tau) cos(th_k tau) ]
!   L       = 1/2 energy_scale sum_k wgt_k ( s I_k + b - I_exp_k )^2
!
! with s (and b, if match_offset) fitted by weighted least squares. The
! computed spectrum is in arbitrary units, so without s the loss compares two
! things on different scales and its gradient is meaningless; b is there for
! experimental files whose baseline has been subtracted to a hard zero, where
! a strictly positive model can never reach the transparency minimum and the
! quadratic loss responds by shrinking the whole prediction. Both sit at their
! own optimum, so by the envelope theorem dL/dI_k picks up no extra term from
! either.
!
! THE GRADIENT. mu of the newest frame (age 0) enters C(tau) twice over: once
! directly, and once through the mean that is subtracted from every frame.
! Writing d(a) = mu(a) - (1/N) sum_b mu(b), so that dd(a)/dmu_c(0) =
! delta_{a,0} - 1/N,
!
!   dC(0)/dmu_c   = [ 2 d_c(0) - P_c(0)/N ] / M_0
!   dC(tau)/dmu_c = [   d_c(tau) - P_c(tau)/N ] / M_tau       tau > 0
!
!   P_c(tau) = sum_{a=0}^{N-1-tau} [ d_c(a) + d_c(a+tau) ]
!            = T_c(N-1-tau) - T_c(tau-1),   T_c(k) = sum_{a=0}^{k} d_c(a)
!
! so one prefix sum over the buffer gives every P_c(tau) and the whole thing
! is still a single pass. P_c(0) = 2 T_c(N-1) = 0 identically, which is a
! useful check on the prefix sum. With subtract_mean off the P terms are absent
! and this reduces to the plain expression.
!
! The correction is not large -- P_c is a partial sum of a mean-zero series, so
! it is O(sqrt(N)) against the O(N) leading term. Measured on a real water
! buffer it is 0.30% of |lambda| at N = 180, and it shrinks as 1/sqrt(N) from
! there. It is kept because it is two lines and it makes the force exactly the
! gradient of the loss: ana/check_gradient.py verifies that against finite
! differences for every combination of the switches, and an approximate
! gradient would pass every spectrum-shaped test while still being wrong.
!
   subroutine mad_ir_evaluate(this, energy_scale, energy, lambda)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      real(dp), intent(in) :: energy_scale
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: lambda(1:3)
      real(dp), allocatable :: dLdI(:)
      real(dp), allocatable :: S(:)
      real(dp), allocatable :: Tpre(:, :)
      real(dp) :: two_pi, theta, acc, s_fit, b_fit, pref, c, denom
      real(dp) :: swi2, swi, sw, swie, swe, det
      real(dp) :: mub(1:3), dnew(1:3), dlag, Pc
      integer :: t, k, m, N, L, c1, a

      energy = 0.d0
      lambda = 0.d0
      if (.not. mad_ir_ready(this)) return

      two_pi = 2.d0*dacos(-1.d0)

!     ---- autocorrelation, then spectrum -----------------------------
!     This also sets n_lag_used, the taper, and mu_mean, all of which the
!     gradient below has to agree with. Recomputing any of them here would be
!     the classic way for the force to stop being the gradient of the loss.
      call mad_ir_spectrum(this)
      N = this%n_stored
      L = this%n_lag_used
      mub(1:3) = this%mu_mean(1:3)

!     ---- overall scale, and baseline if asked for --------------------
      s_fit = 1.d0
      b_fit = 0.d0
      if (this%match_scale) then
         swi2 = 0.d0
         swi = 0.d0
         sw = 0.d0
         swie = 0.d0
         swe = 0.d0
         do k = 1, this%n_freq
            swi2 = swi2 + this%wgt(k)*this%I_calc(k)**2
            swi = swi + this%wgt(k)*this%I_calc(k)
            sw = sw + this%wgt(k)
            swie = swie + this%wgt(k)*this%I_calc(k)*this%I_exp(k)
            swe = swe + this%wgt(k)*this%I_exp(k)
         end do
         if (this%match_offset) then
!           [ swi2 swi ] [ s ]   [ swie ]
!           [ swi  sw  ] [ b ] = [ swe  ]
            det = swi2*sw - swi*swi
            if (dabs(det) > 1.d-300) then
               s_fit = (swie*sw - swe*swi)/det
               b_fit = (swi2*swe - swi*swie)/det
            else if (swi2 > 0.d0) then
               s_fit = swie/swi2
            end if
!           A NEGATIVE SCALE IS NOT A FIT, IT IS AN INVERSION. The two-parameter
!           solve is unconstrained, and on a prediction whose shape does not
!           resemble the experiment it will happily return s < 0 -- putting the
!           model's peaks where the experiment's troughs are and absorbing the
!           rest into b. Measured on this potential, whose dipole model produces
!           no vibrational bands, the offset fit returned s = -9.7e-9 with
!           b = +2.70.
!
!           That is not merely a poor fit. dL/dI_k carries a factor of s, so a
!           negative scale REVERSES the bias: it would drive the model away from
!           having the bands it is being asked to grow. Fall back to the
!           scale-only fit, which cannot change sign, and leave the baseline at
!           zero. The header records the offset actually used, so a run that
!           took this path is identifiable after the fact.
            if (s_fit <= 0.d0) then
               b_fit = 0.d0
               if (swi2 > 0.d0) then
                  s_fit = swie/swi2
               else
                  s_fit = 1.d0
               end if
            end if
         else
            if (swi2 > 0.d0) s_fit = swie/swi2
         end if
      end if
      this%scale = s_fit
      this%offset = b_fit

!     ---- energy and its sensitivity to the spectrum ------------------
!     The same form every other MAD observable uses (exp_utils'
!     get_exp_energies): E = 1/2 * energy_scale * sum (y_pred - y_exp)^2, so
!     that exp_energy_scales means the same thing here as it does for a pair
!     distribution or a structure factor, and the force really is the gradient
!     of the spectral difference rather than of something proportional to it.
      allocate (dLdI(1:this%n_freq))
      do k = 1, this%n_freq
         c = s_fit*this%I_calc(k) + b_fit - this%I_exp(k)
         energy = energy + 0.5d0*energy_scale*this%wgt(k)*c**2
         dLdI(k) = energy_scale*this%wgt(k)*s_fit*c
      end do

!     ---- back through the transform, to dL/dC(tau) -------------------
      allocate (S(0:this%n_lag))
      S = 0.d0
      do k = 1, this%n_freq
         pref = dLdI(k)*(this%nu(k)**this%nu_power)*this%dt
         S(0) = S(0) + pref*this%win(0)
         do t = 1, L
            theta = two_pi*this%nu(k)*dfloat(t)*this%dt/CM_PER_INV_FS
            S(t) = S(t) + pref*2.d0*this%win(t)*dcos(theta)
         end do
      end do

!     ---- prefix sums of the fluctuation, for the mean-subtraction term
!     Tpre(:,k) = sum_{a=0}^{k} d(a). Only needed when the mean is subtracted;
!     otherwise the newest dipole enters C only directly.
      if (this%subtract_mean) then
         allocate (Tpre(1:3, -1:N - 1))
         Tpre(1:3, -1) = 0.d0
         do a = 0, N - 1
            m = this%head - a
            if (m < 1) m = m + this%n_window
            Tpre(1:3, a) = Tpre(1:3, a - 1) + (this%mu_hist(1:3, m) - mub(1:3))
         end do
      end if

!     ---- and through the autocorrelation, to dL/dmu(newest) ----------
      m = this%head
      dnew(1:3) = this%mu_hist(1:3, m) - mub(1:3)
      do c1 = 1, 3
!        tau = 0 is N pairs under either estimator, and P(0) = 2*T(N-1) = 0
!        identically, so there is no correction term at zero lag.
         acc = 2.d0*S(0)*dnew(c1)/dfloat(N)
         do t = 1, L
            m = this%head - t
            if (m < 1) m = m + this%n_window
            dlag = this%mu_hist(c1, m) - mub(c1)
            if (this%biased_estimator) then
               denom = dfloat(N)
            else
               denom = dfloat(N - t)
            end if
            if (this%subtract_mean) then
               Pc = Tpre(c1, N - 1 - t) - Tpre(c1, t - 1)
               acc = acc + S(t)*(dlag - Pc/dfloat(N))/denom
            else
               acc = acc + S(t)*dlag/denom
            end if
         end do
         lambda(c1) = acc
      end do

      if (allocated(Tpre)) deallocate (Tpre)
      deallocate (dLdI, S)

   end subroutine mad_ir_evaluate

!**************************************************************************
!
! Everything a driver needs to start a MAD IR run: size the ensemble from the
! spectral requirements, read the experiment, allocate the per-atom gradient
! accumulator, and pick up a previous history buffer if there is one.
!
! dt_md is the MD timestep in fs and stride the number of steps between stored
! configurations; the sampling interval is their product, and that is what the
! Nyquist check is applied to.
!
! resumed says whether a restart buffer was adopted. A refused or missing one
! is not an error -- the run simply starts filling a fresh ensemble -- but the
! caller should say so, because it means no bias for the next n_window
! samples.
!
   subroutine mad_ir_setup(dt_md, stride, nu_res, nu_min, nu_max, lag_factor, &
                           nu_in, I_in, restart_file, match_scale, nu_power, &
                           window_kind, subtract_mean, biased_estimator, &
                           taper_partial, match_offset, weight_by_spacing, &
                           n_atoms, ok, resumed, msg)

      implicit none

      real(dp), intent(in) :: dt_md
      integer, intent(in) :: stride
      real(dp), intent(in) :: nu_res, nu_min, nu_max
      integer, intent(in) :: lag_factor
      real(dp), intent(in) :: nu_in(:), I_in(:)
      character(len=*), intent(in) :: restart_file, window_kind
      logical, intent(in) :: match_scale
      logical, intent(in) :: subtract_mean, biased_estimator, taper_partial
      logical, intent(in) :: match_offset, weight_by_spacing
      real(dp), intent(in) :: nu_power
      integer, intent(in) :: n_atoms
      logical, intent(out) :: ok, resumed
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: nu(:), I_exp(:), wgt(:)
      real(dp) :: dt
      integer :: n_lag, n_window
      logical :: ok2
      character(len=512) :: msg2

      ok = .false.
      resumed = .false.
      msg = ""

      if (stride < 1) then
         msg = "mad_ir: mad_ir_stride must be at least 1"
         return
      end if
      dt = dt_md*dfloat(stride)

      call mad_ir_size_window(dt, nu_res, nu_max, lag_factor, n_lag, n_window, ok2, msg)
      if (.not. ok2) return

      call mad_ir_select_range(nu_in, I_in, nu_min, nu_max, weight_by_spacing, &
                               nu, I_exp, wgt, ok2, msg)
      if (.not. ok2) return

      call mad_ir_init(mad_ir_state, dt, n_lag, n_window, nu, I_exp, wgt, &
                       match_scale, nu_power, window_kind, subtract_mean, &
                       biased_estimator, taper_partial, match_offset)

      if (allocated(mad_ir_dmu_dr)) deallocate (mad_ir_dmu_dr)
      allocate (mad_ir_dmu_dr(1:3, 1:3, 1:n_atoms))
      mad_ir_dmu_dr = 0.d0

      if (trim(restart_file) /= "none") then
         call mad_ir_load(mad_ir_state, restart_file, resumed, msg2)
         if (.not. resumed) msg = trim(msg2)
      end if

      ok = .true.

   end subroutine mad_ir_setup

!**************************************************************************
!
! Set up a PREDICTION run: no experiment, no bias, one spectrum at the end.
!
! The sizing question is the opposite way round from the bias. There, the
! ensemble is a rolling window whose length follows from the resolution asked
! for, and the run has to be long enough to fill it. Here the trajectory is
! the ensemble -- "run for this long, then tell me the spectrum" -- so the
! length is given and the resolution follows from it:
!
!   n_window = n_frames                      every frame the run will produce
!   n_lag    = n_window / lag_factor         so C(tau) at the longest lag is
!                                            still averaged over n_lag pairs
!   d(nu)    = CM_PER_INV_FS / (n_lag dt)    what that buys
!
! A resolution asked for explicitly is honoured when the run is long enough
! for it and refused with the required step count when it is not, because the
! alternative -- quietly giving a coarser spectrum than the input asked for --
! is the failure mode this whole file is organised against.
!
   subroutine mad_ir_setup_predict(dt_md, stride, n_frames, nu_res, nu_min, nu_max, &
                                   lag_factor, n_samples, nu_power, window_kind, &
                                   subtract_mean, biased_estimator, taper_partial, &
                                   n_atoms, ok, msg)

      implicit none

      real(dp), intent(in) :: dt_md
      integer, intent(in) :: stride, n_frames, lag_factor, n_samples, n_atoms
      real(dp), intent(in) :: nu_res, nu_min, nu_max, nu_power
      character(len=*), intent(in) :: window_kind
      logical, intent(in) :: subtract_mean, biased_estimator, taper_partial
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: nu(:), I_exp(:), wgt(:)
      real(dp) :: dt, nyquist, dnu
      integer :: n_lag, n_window, n_freq, k, n_lag_want

      ok = .false.
      msg = ""

      if (stride < 1) then
         msg = "mad_ir: ir_stride must be at least 1"
         return
      end if
      if (lag_factor < 1) then
         msg = "mad_ir: ir_lag_factor must be at least 1"
         return
      end if
      dt = dt_md*dfloat(stride)
      if (dt <= 0.d0) then
         msg = "mad_ir: the sampling interval must be positive"
         return
      end if

!     The Nyquist check is the same one the bias applies, and for the same
!     reason: power above the limit does not go missing, it folds back down
!     into the range being reported.
      nyquist = CM_PER_INV_FS/(2.d0*dt)
      if (nu_max > nyquist) then
         write (msg, '(A,F10.1,A,F10.1,A)') &
            "mad_ir: a sampling interval of md_step*ir_stride only reaches ", nyquist, &
            " cm^-1 but ir_nu_max is ", nu_max, " cm^-1; sample more often"
         return
      end if
      if (nu_max <= nu_min) then
         msg = "mad_ir: ir_nu_max must be greater than ir_nu_min"
         return
      end if

      n_window = n_frames
      if (n_window < 2*lag_factor) then
         write (msg, '(A,I0,A)') &
            "mad_ir: a do_ir run needs at least ", 2*lag_factor, &
            " sampled frames; raise md_nsteps or lower ir_stride"
         return
      end if
      n_lag = n_window/lag_factor

      if (nu_res > 0.d0) then
         n_lag_want = ceiling(CM_PER_INV_FS/(nu_res*dt))
         if (n_lag_want > n_lag) then
            write (msg, '(A,F8.2,A,I0,A,I0,A)') &
               "mad_ir: ir_resolution = ", nu_res, &
               " cm^-1 needs a longest lag of ", n_lag_want, &
               ", so md_nsteps must be at least ", n_lag_want*lag_factor*stride, &
               "; drop ir_resolution to take whatever the run length gives"
            return
         end if
         n_lag = n_lag_want
      end if
      if (n_lag < 1) n_lag = 1

!     The grid. One point per resolution element is the honest default -- a
!     finer one interpolates the lag window and invites the eye to read
!     structure that the transform cannot carry -- but a finer grid does make
!     the curve easier to look at, so n_samples overrides it.
      dnu = CM_PER_INV_FS/(dfloat(n_lag)*dt)
      if (n_samples > 1) then
         n_freq = n_samples
      else
         n_freq = int((nu_max - nu_min)/dnu) + 1
         if (n_freq < 2) n_freq = 2
      end if

      allocate (nu(1:n_freq), I_exp(1:n_freq), wgt(1:n_freq))
      do k = 1, n_freq
         nu(k) = nu_min + (nu_max - nu_min)*dfloat(k - 1)/dfloat(n_freq - 1)
      end do
!     There is no experiment. I_exp and wgt exist because the type carries
!     them; match_scale is off so nothing is fitted against them and the
!     intensity stays in the arbitrary units the transform produces.
      I_exp = 0.d0
      wgt = 1.d0

      call mad_ir_init(mad_ir_state, dt, n_lag, n_window, nu, I_exp, wgt, &
                       .false., nu_power, window_kind, subtract_mean, &
                       biased_estimator, taper_partial, .false.)

!     No bias means no dmu/dr, but the accumulator is allocated anyway so that
!     the gate in gap_interface is the only thing deciding whether it is
!     filled, rather than allocation state.
      if (allocated(mad_ir_dmu_dr)) deallocate (mad_ir_dmu_dr)
      allocate (mad_ir_dmu_dr(1:3, 1:3, 1:n_atoms))
      mad_ir_dmu_dr = 0.d0

      deallocate (nu, I_exp, wgt)
      ok = .true.

   end subroutine mad_ir_setup_predict

!**************************************************************************
!
! The spectrum, with the provenance needed to read it.
!
! A dipole-autocorrelation spectrum is trustworthy over a band, not
! everywhere, and the band is set by the sampling rather than by the physics:
!
!   above CM_PER_INV_FS/(2 dt), the Nyquist limit, nothing is represented, and
!   what was there has folded back down into the range below. The top of the
!   range is therefore not merely empty, it is contaminated, which is why the
!   header names a margin as well as the limit.
!
!   below CM_PER_INV_FS/(n_lag dt), one resolution element, nothing is
!   resolved: two bands closer than that are one band here, and a feature
!   narrower than that is the lag window.
!
!   and across all of it, C(tau) at lag tau is an average over n_frames - tau
!   pairs, so the long-lag end -- which is what carries the fine structure --
!   is the noisiest part of the transform however long the run was.
!
! Writing those numbers into the file rather than into a log is deliberate:
! the file outlives the run that made it, and a spectrum without its sampling
! interval cannot be checked by anyone who did not launch the job.
!
! with_exp adds the experimental column and the fitted scale, i.e. the file
! says what the bias was comparing. A do_ir run has neither, and printing a
! column of zeros for the experiment would invite it to be plotted.
!
   subroutine mad_ir_write_spectrum(this, fname, with_exp, n_md_steps, md_step_fs)
      implicit none
      type(mad_ir_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(in) :: with_exp
      integer, intent(in) :: n_md_steps
      real(dp), intent(in) :: md_step_fs
      integer :: iu, ios, k, n_lag_here
      real(dp) :: nyquist, dnu, alias_edge, trust_lo, trust_hi

      if (.not. this%active) return
      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) return

      nyquist = CM_PER_INV_FS/(2.d0*this%dt)
!     The lags ACTUALLY transformed, not the ones the run was sized for. While
!     a do_ir ensemble is filling these differ, and quoting the sized-for
!     resolution then would promise a sharpness the file does not have.
      n_lag_here = this%n_lag_used
      if (n_lag_here < 1) n_lag_here = this%n_lag
      dnu = CM_PER_INV_FS/(dfloat(n_lag_here)*this%dt)
      alias_edge = 0.9d0*nyquist
!     The band quoted is the one this FILE can be read over, so it is the
!     sampling limits intersected with the grid actually written. Quoting a
!     Nyquist limit eight times the highest wavenumber in the file would be
!     true and useless.
      trust_lo = max(dnu, this%nu(1))
      trust_hi = min(alias_edge, this%nu(this%n_freq))

      write (iu, '(A)') "#  TurboGAP IR spectrum"
      write (iu, '(A)') "#  Fourier transform of the total-dipole autocorrelation function,"
      write (iu, '(A)') "#  I(nu) = nu^p * dt * FT[C(tau)], intensity in arbitrary units."
      write (iu, '(A)') "#"
      write (iu, '(A,F14.4,A)') "#  MD timestep                  = ", md_step_fs, " fs"
      write (iu, '(A,I14)') "#  MD steps run                 = ", n_md_steps
      write (iu, '(A,F14.4,A)') "#  sampling interval dt         = ", this%dt, " fs"
      write (iu, '(A,I14)') "#  frames in the ensemble       = ", this%n_stored
      write (iu, '(A,I14)') "#  longest lag transformed      = ", n_lag_here
      write (iu, '(A,I14)') "#  longest lag sized for        = ", this%n_lag
      write (iu, '(A,F14.4)') "#  nu power p                   = ", this%nu_power
      write (iu, '(A,A14)') "#  lag window                   = ", trim(this%window_kind)
      if (this%subtract_mean) then
         write (iu, '(A)') "#  correlating                  =  mu - <mu>"
      else
         write (iu, '(A)') "#  correlating                  =  mu  (mean NOT subtracted)"
      end if
      if (this%biased_estimator) then
         write (iu, '(A)') "#  ACF estimator                =  biased, 1/N"
      else
         write (iu, '(A)') "#  ACF estimator                =  unbiased, 1/(N-tau)"
      end if
      if (with_exp) then
         write (iu, '(A,ES14.6)') "#  fitted scale                 = ", this%scale
         if (this%match_offset) then
            write (iu, '(A,ES14.6)') "#  fitted offset                = ", this%offset
         end if
      end if
      write (iu, '(A)') "#"
      if (trust_hi > trust_lo) then
         write (iu, '(A,F12.2,A,F12.2,A)') "#  TRUSTWORTHY RANGE:", trust_lo, "  to", trust_hi, "  cm^-1"
      else
         write (iu, '(A)') "#  TRUSTWORTHY RANGE: none of it."
         write (iu, '(A)') "#    The grid written lies entirely outside the band the sampling"
         write (iu, '(A)') "#    supports. Raise the run length (lowers the resolution figure"
         write (iu, '(A)') "#    below) or lower md_step*ir_stride (raises the Nyquist limit)."
      end if
      write (iu, '(A)') "#"
      write (iu, '(A,F12.2,A)') "#    Nyquist limit          = ", nyquist, " cm^-1 = 1/(2 dt)."
      write (iu, '(A)') "#      Nothing above it is represented. Power that was there is not"
      write (iu, '(A)') "#      absent from this file, it has folded back down into the range"
      write (iu, '(A)') "#      below, so read the last 10% up to the Nyquist limit -- above"
      write (iu, '(A)') "#      the upper bound quoted above -- as possibly aliased, not signal."
      write (iu, '(A,F12.2,A)') "#    resolution             = ", dnu, " cm^-1 = 1/(n_lag dt)."
      write (iu, '(A)') "#      Two bands closer together than this are one band here, and"
      write (iu, '(A)') "#      structure narrower than it is the lag window, not the spectrum."
      write (iu, '(A)') "#      Below one resolution element nothing is resolved at all."
      write (iu, '(A)') "#"
      write (iu, '(A,F12.2,A,F12.2,A)') "#    grid written           = ", this%nu(1), " to", &
         this%nu(this%n_freq), " cm^-1; the range above is"
      write (iu, '(A)') "#      that grid intersected with [resolution, 0.9*Nyquist]."
      write (iu, '(A)') "#"
      write (iu, '(A)') "#    C(tau) is averaged over n_frames - tau pairs, so the long-lag"
      write (iu, '(A)') "#    end is the noisiest part of the transform: fine structure near"
      write (iu, '(A)') "#    the resolution limit needs a longer run, not a finer grid."
      write (iu, '(A)') "#"
      if (with_exp) then
         write (iu, '(A)') "#  wavenumber_cm-1        predicted           experimental"
         do k = 1, this%n_freq
            write (iu, '(1X,ES20.10,1X,ES20.10,1X,ES20.10)') &
               this%nu(k), this%scale*this%I_calc(k) + this%offset, this%I_exp(k)
         end do
      else
         write (iu, '(A)') "#  wavenumber_cm-1        intensity"
         do k = 1, this%n_freq
            write (iu, '(1X,ES20.10,1X,ES20.10)') this%nu(k), this%I_calc(k)
         end do
      end if
      close (iu)
   end subroutine mad_ir_write_spectrum

!**************************************************************************
!
! The predicted spectrum appended as one more block, so that the whole
! trajectory of predictions survives the run rather than only its last frame.
!
! The file format is the one every other MAD observable's prediction uses
! (read_files' write_exp_datan, as written for xrd_prediction.dat): a "# label"
! line, then one two-column block per write, blocks separated by a blank line.
! That is gnuplot's `index` format, so `plot "ir_prediction.dat" index n` picks
! out a frame and `every :::n::n` a range, with no post-processing.
!
! It is a separate routine from write_exp_datan rather than a call to it
! because mad_ir sits below read_files in the module graph -- read_files uses
! types, types would then have to know about the spectrum -- and duplicating
! twelve lines of formatting is cheaper than that cycle.
!
! `overwrite` cannot be derived from the step number the way it is for the
! per-frame observables: the first block appears whenever the ensemble first
! fills, which is not step zero, and opening with status="old" before anything
! has been written is a runtime error rather than an empty file. The caller
! owns that flag and flips it after the first successful write.
!
   subroutine mad_ir_append_spectrum(this, fname, overwrite, istep, time_fs)
      implicit none
      type(mad_ir_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(in) :: overwrite
      integer, intent(in) :: istep
      real(dp), intent(in) :: time_fs
      integer :: iu, ios, k

      if (.not. this%active) return

      if (overwrite) then
         open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
         if (ios /= 0) return
         write (iu, '(A)') "#  ir: wavenumber_cm-1  predicted_intensity (scaled to the experiment)"
      else
         open (newunit=iu, file=trim(fname), status="old", position="append", &
               action="write", iostat=ios)
         if (ios /= 0) return
         write (iu, '(A)') " "
      end if

!     A header per block: gnuplot ignores it as a comment, and without it the
!     blocks are only identifiable by counting them.
      write (iu, '(A,I10,A,ES16.8,A,ES16.8,A,ES16.8,A,I8)') "# step ", istep, &
         "   time_fs ", time_fs, "   scale ", this%scale, "   offset ", this%offset, &
         "   n_lag ", this%n_lag_used
      do k = 1, this%n_freq
         write (iu, '(1X,ES20.10,1X,ES20.10)') this%nu(k), &
            this%scale*this%I_calc(k) + this%offset
      end do
      close (iu)
   end subroutine mad_ir_append_spectrum

!**************************************************************************
!
! The experimental spectrum as the fit actually sees it: restricted to
! [nu_min, nu_max] and on its own grid. The file named by exp_data_files is
! the raw one, so plotting the prediction against it is plotting against
! points that were never fitted. Written once, mirroring the "<label>_exp.dat"
! that the per-frame observables get from turbogap.f90.
!
   subroutine mad_ir_write_exp_spectrum(this, fname)
      implicit none
      type(mad_ir_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      integer :: iu, ios, k
      if (.not. this%active) return
      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) return
      write (iu, '(A)') "#  ir: wavenumber_cm-1  experimental_intensity  fit_weight"
      do k = 1, this%n_freq
         write (iu, '(1X,ES20.10,1X,ES20.10,1X,ES20.10)') this%nu(k), this%I_exp(k), this%wgt(k)
      end do
      close (iu)
   end subroutine mad_ir_write_exp_spectrum

!**************************************************************************
!
! Keep the experimental points inside [nu_min, nu_max].
!
! The data itself arrives already read: read_exp_data pulled it in when
! exp_data_files was parsed, exactly as it does for every other observable, so
! there is no second reader here and no second file format to keep in step.
!
! The fit grid IS the experimental grid, restricted. Interpolating the
! experiment onto a grid of our own choosing would invent structure between its
! points and then fit to it.
!
   subroutine mad_ir_select_range(nu_in, I_in, nu_min, nu_max, weight_by_spacing, &
                                  nu, I_exp, wgt, ok, msg)

      implicit none

      real(dp), intent(in) :: nu_in(:)
      real(dp), intent(in) :: I_in(:)
      real(dp), intent(in) :: nu_min
      real(dp), intent(in) :: nu_max
      logical, intent(in) :: weight_by_spacing
      real(dp), allocatable, intent(out) :: nu(:)
      real(dp), allocatable, intent(out) :: I_exp(:)
      real(dp), allocatable, intent(out) :: wgt(:)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: wsum
      integer :: n, k, j

      ok = .false.
      msg = ""
      n = 0
      do k = 1, size(nu_in)
         if (nu_in(k) >= nu_min .and. nu_in(k) <= nu_max) n = n + 1
      end do
      if (n < 2) then
         write (msg, '(A,F9.1,A,F9.1,A)') &
            "mad_ir: the experimental spectrum has fewer than two points between ", &
            nu_min, " and ", nu_max, " cm^-1"
         return
      end if
      allocate (nu(1:n), I_exp(1:n), wgt(1:n))
      j = 0
      do k = 1, size(nu_in)
         if (nu_in(k) >= nu_min .and. nu_in(k) <= nu_max) then
            j = j + 1
            nu(j) = nu_in(k)
            I_exp(j) = I_in(k)
         end if
      end do
!     The loss is a SUM over experimental points, and an experimental grid is
!     rarely uniform: the Downing & Williams water data shipped with the test
!     is sampled at 23.9 cm^-1 across the O-H stretch and ~55 cm^-1 elsewhere,
!     so 45% of an unweighted loss falls in the top quarter of the range -- a
!     weighting chosen by whoever digitised the paper. Trapezoidal weights make
!     it an integral over frequency instead, which is what was meant.
!
!     Normalised to mean 1 so that exp_energy_scales keeps the magnitude it has
!     everywhere else; without that the loss would pick up a factor of the
!     frequency range and every scale in every old input would silently change
!     meaning.
      if (weight_by_spacing .and. n > 2) then
         wgt(1) = nu(2) - nu(1)
         wgt(n) = nu(n) - nu(n - 1)
         do k = 2, n - 1
            wgt(k) = 0.5d0*(nu(k + 1) - nu(k - 1))
         end do
         wsum = sum(wgt(1:n))
         if (wsum > 0.d0) then
            wgt = wgt*dfloat(n)/wsum
         else
            wgt = 1.d0
         end if
      else
         wgt = 1.d0
      end if
      ok = .true.

   end subroutine mad_ir_select_range

!**************************************************************************
!
! Restart. The history buffer is the expensive thing in a MAD IR run -- a
! 2085-lag window at 2 fs sampling is 4 ps of trajectory -- so losing it means
! 4 ps with no force applied while it refills. Small enough to write as text,
! which also makes it inspectable.
!
   subroutine mad_ir_save(this, fname, ok, msg)

      implicit none

      type(mad_ir_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: iu, ios, k

      ok = .false.
      msg = ""
      if (.not. this%active) then
         msg = "mad_ir: nothing to save"
         return
      end if
      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         msg = "mad_ir: cannot write "//trim(fname)
         return
      end if
      write (iu, '(A)') "MAD_IR_RESTART 1"
      write (iu, '(ES24.16)') this%dt
      write (iu, '(3(1X,I0))') this%n_lag, this%n_window, this%n_stored
      write (iu, '(I0)') this%head
      do k = 1, this%n_window
         write (iu, '(3ES24.16)') this%mu_hist(1, k), this%mu_hist(2, k), this%mu_hist(3, k)
      end do
      close (iu)
      ok = .true.

   end subroutine mad_ir_save

!**************************************************************************
!
! Load a history buffer into an already-initialised state.
!
! The sizing is checked rather than adopted. A buffer written with a different
! dt or n_lag describes a different spectrum, and silently continuing with it
! would change the observable mid-run in a way nothing downstream could
! notice. Refusing is the only safe answer; the caller can then choose to start
! a fresh buffer.
!
   subroutine mad_ir_load(this, fname, ok, msg)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: iu, ios, k, n_lag_in, n_window_in, n_stored_in, head_in, ver
      real(dp) :: dt_in
      character(len=64) :: tag

      ok = .false.
      msg = ""
      if (.not. this%active) then
         msg = "mad_ir: load called before init"
         return
      end if
      open (newunit=iu, file=trim(fname), status="old", action="read", iostat=ios)
      if (ios /= 0) then
         msg = "mad_ir: no restart file "//trim(fname)
         return
      end if
      read (iu, *, iostat=ios) tag, ver
      if (ios /= 0 .or. trim(tag) /= "MAD_IR_RESTART") then
         close (iu)
         msg = "mad_ir: "//trim(fname)//" is not a MAD IR restart file"
         return
      end if
      read (iu, *, iostat=ios) dt_in
      if (ios == 0) read (iu, *, iostat=ios) n_lag_in, n_window_in, n_stored_in
      if (ios == 0) read (iu, *, iostat=ios) head_in
      if (ios /= 0) then
         close (iu)
         msg = "mad_ir: "//trim(fname)//" is truncated"
         return
      end if

      if (n_lag_in /= this%n_lag .or. n_window_in /= this%n_window .or. &
          dabs(dt_in - this%dt) > 1.d-10*max(1.d0, dabs(this%dt))) then
         close (iu)
         write (msg, '(A,I0,A,I0,A,F10.4,A,I0,A,I0,A,F10.4)') &
            "mad_ir: restart was written for n_lag=", n_lag_in, " n_window=", n_window_in, &
            " dt=", dt_in, " but this run wants n_lag=", this%n_lag, &
            " n_window=", this%n_window, " dt=", this%dt
         return
      end if

      do k = 1, this%n_window
         read (iu, *, iostat=ios) this%mu_hist(1, k), this%mu_hist(2, k), this%mu_hist(3, k)
         if (ios /= 0) then
            close (iu)
            this%n_stored = 0
            this%head = 0
            msg = "mad_ir: "//trim(fname)//" ended early; the buffer is unusable"
            return
         end if
      end do
      close (iu)
      this%n_stored = n_stored_in
      this%head = head_in
      ok = .true.

   end subroutine mad_ir_load

!**************************************************************************
!
! f_jb = - dE/dr_jb = - sum_a lambda_a dmu_a/dr_jb, with lambda = dE/dmu of the
! newest configuration. E is the MAD energy, so this is the gradient of the
! spectral difference and nothing else.
!
! dmu_dr(a,b,j) is what gap.f90's accumulate_dmu_dr builds. The forces are
! ADDED to whatever is already in forces, since a MAD run carries a real
! potential alongside the bias.
!
   subroutine mad_ir_forces(lambda, dmu_dr, forces)
      implicit none
      real(dp), intent(in) :: lambda(1:3)
      real(dp), intent(in) :: dmu_dr(:, :, :)
      real(dp), intent(inout) :: forces(:, :)
      integer :: j, b, a, n_atoms
      real(dp) :: acc
      n_atoms = size(dmu_dr, 3)
      do j = 1, n_atoms
         do b = 1, 3
            acc = 0.d0
            do a = 1, 3
               acc = acc + lambda(a)*dmu_dr(a, b, j)
            end do
            forces(b, j) = forces(b, j) - acc
         end do
      end do
   end subroutine mad_ir_forces

end module mad_ir
