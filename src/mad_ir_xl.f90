! Auxiliary-variable MAD IR: a resonator bank carrying the spectrum in place of
! a stored trajectory.
!
! Selected with ir_bias_mode = "xl". The observable is still an infrared
! spectrum and the target is still an experimental one, but nothing here forms
! an autocorrelation, stores a trajectory, or differentiates backwards through
! time. The spectrum is carried by extra degrees of freedom integrated
! alongside the atoms, and the bias is the gradient of the mismatch with
! respect to the dipoles that drive them.
!
! THE SCHEME
!
! Attach to the dipole a pair of damped resonators for every fitted frequency
! w_k, driven by the dipole and by its time derivative:
!
!   xddot + gamma xdot + w_k^2 x = m(t)                                 (1)
!   yddot + gamma ydot + w_k^2 y = mdot(t) / w_k                        (2)
!
! Equation (1) is a Lorentzian bandpass filter written as an equation of
! motion: the bank IS the Fourier transform, integrated rather than
! transformed. Equation (2) is its quadrature partner, scaled so that the two
! have equal stationary variance, so that
!
!   R_k^2 = |x_k|^2 + |y_k|^2                                           (3)
!
! is the slowly varying ENVELOPE of the response at w_k rather than something
! oscillating at 2 w_k. The experiment is compared with R_k^2, and the bias is
!
!   U = 1/2 energy_scale sum_k wgt_k ( s I_k + b - I_exp_k )^2,
!   I_k = pref_k R_k^2                                                  (4)
!
! -- the same form, the same fitted scale and the same meaning of
! exp_energy_scales as every other MAD observable.
!
! This is the extended-Lagrangian idea in the only form that turns out to work.
! What it is not is set out below, because the difference is the whole content
! of the file.
!
! WHY THE RESTRAINT IS NOT A POTENTIAL IN THE AUXILIARY COORDINATES
!
! The natural way to write this scheme is to make U a term in an extended
! Lagrangian alongside the resonators, coupling them to the atoms by g x.m(q)
! and letting U act on x:
!
!   L = L_MD + sum_k [ 1/2 xdot^2 - 1/2 w_k^2 x^2 + g x_k . m(q) ] - U({R})
!
! It is an appealing picture -- a tuning fork bolted to every molecule, with
! the experiment pushing on the forks -- and it does not work. The reason is
! worth stating precisely, because it is not obvious and it is invisible in any
! test that does not check the SIGN of the response.
!
! dU/dx_k = 2 dLdI_k pref_k x_k is proportional to x itself. It is therefore
! not a force on the resonator at all: it is a shift in its spring constant,
!
!   w_k^2  ->  w_k^2 + 2 dLdI_k pref_k
!
! and a band that is too weak gives dLdI_k < 0, which SOFTENS the resonator. A
! softened resonator is no longer tuned to w_k. It is detuned from the very
! frequency whose intensity it was measuring, and a driven oscillator detuned
! from its drive responds LESS, not more:
!
!   |x| = F / sqrt( (w_eff^2 - w_d^2)^2 + gamma^2 w_d^2 )
!
! is maximal at w_eff = w_d and falls off either side. So the restraint that
! was supposed to raise the band lowers it. Measured, in
! tests/mad_ir/xlverify.f90 as it was first written: a target four times the
! prediction drove R^2 down by two orders of magnitude, monotonically, in both
! the coherent and the incoherent form.
!
! Two things are wrong at once, and the second is worse than the first. The
! bias moves the band the wrong way; and a detuned resonator is no longer
! reporting the intensity at w_k, so the OBSERVABLE is corrupted by the
! restraint acting on it. There is no choice of sign or magnitude that repairs
! either, because the response is not monotone in the shift: it has a maximum
! at zero shift.
!
! What survives is the part of the picture that was doing the work. The
! resonators are a filter bank -- an analogue spectrometer running alongside
! the dynamics -- and the bias belongs on the ATOMS, as the gradient of the
! mismatch with respect to the dipole that drives the bank. The bank itself
! runs free, at its own frequencies, undisturbed by what is being asked of the
! spectrum.
!
! THE GRADIENT
!
! The dipole of the current configuration enters the bank through exactly one
! advance, and the propagator below is a closed form, so its sensitivity is
! available in closed form too. Over one interval h with the drive held
! constant at F,
!
!   x(t+h) = F/w^2 + e^{-lam h} [ u0 cos(Om h) + (v0 + lam u0) sin(Om h)/Om ]
!
! with u0 = x(t) - F/w^2, so
!
!   P_k = dx(t+h)/dF = (1/w_k^2) [ 1 - e^{-lam h} ( cos(Om h)
!                                                 + lam sin(Om h)/Om ) ]   (5)
!
! and the y drive is mdot/w_k, so its sensitivity is P_k times whatever
! d(mdot)/dm(t) the difference formula gives. THAT COEFFICIENT IS NOT 1/h. The
! three-point backward form used below, mdot = (3m(t) - 4m(t-h) + m(t-2h))/(2h),
! differentiates to 3/(2h), and using 1/h instead leaves the gradient wrong by
! a factor of order 1.5 -- large, constant, and completely invisible in any
! spectrum. It was found by the h-scan below coming out flat, which is what an
! h-scan is for. Therefore
!
!   dR_k^2(t+h)/dm_a(t) = 2 P_k [ x_{k,a}(t+h) + (3/(2 h w_k)) y_{k,a}(t+h) ] (6)
!
! and the weight the descriptor pass contracts against the ML dipole gradient
! is the negative gradient of the bias energy,
!
!   W_a = - sum_k dLdI_k pref_k * 2 P_k [ x_{k,a} + 3 y_{k,a}/(2 h w_k) ]  (7)
!   f_jb = sum_a W_a d m_a / d r_jb                                        (8)
!
! (7) is EXACT, not a model of the gradient: tests/mad_ir/xlverify.f90 checks
! it against central differences of the loss in the driving dipole and the
! error falls as h^2 to round-off. That matters more than it might seem. The
! sign of a spectral bias is not something that can be reasoned about from the
! phase of a filtered signal -- the argument above is an example of getting it
! confidently wrong -- and a gradient is the only construction that is right by
! definition.
!
! ONE FRAME OF LAG, AND WHY IT IS THE PRICE OF THE WHOLE APPROACH. (7) needs
! the bank state AFTER the current dipole has been folded in, so it cannot be
! formed until the descriptor pass that produced that dipole is over. The
! weight contracted at stored frame n is therefore the one formed at the end of
! frame n-1: the exact gradient of the bias energy with respect to the dipoles
! of frame n-1, contracted against the dipole gradient at frame n's positions.
!
! The alternative is the ACF bias's: keep the whole (3,3,n_atoms) tensor from
! the pass and contract it afterwards with a weight belonging to the current
! frame. That is exact in the positions and costs nine numbers per atom plus a
! 9*n_atoms all-reduce; this costs three, contracts inside the pass, and is a
! third of the work in the descriptor term. One stored frame of lag is far
! below the bias's own time resolution -- the bank has bandwidth gamma and
! cannot respond to anything faster than 1/gamma = tau/2 anyway -- so it is the
! right trade here and the wrong one there.
!
! FOUR THINGS THAT ARE WRONG QUIETLY
!
! 1. AN UNDAMPED RESONATOR HAS NO AMPLITUDE TO MEASURE. Driven on resonance it
!    grows without limit, so R_k never reaches a value that can be compared
!    with anything. The bank must be damped, and then gamma is not a nuisance
!    parameter but THE RESOLUTION: (1) is a Lorentzian of full width gamma, so
!
!      d(nu) = CM_PER_INV_FS gamma / (2 pi) = CM_PER_INV_FS / (pi tau)
!
!    with tau = 2/gamma the memory time (ir_xl_tau_mem). This is the same
!    quantity that n_lag*dt sets for the block estimator, reached from the other
!    side: 4 cm^-1 needs about 2650 fs either way.
!
! 2. R^2 IS NOT THE INTENSITY, AND THE FACTOR BETWEEN THEM DEPENDS ON w. The
!    stationary variance of (1) driven by a signal of two-sided power spectral
!    density S, smooth across the filter width, is
!
!      <x^2> = S(w_k) INT dw/2pi |H_k(w)|^2 = S(w_k) / (2 gamma w_k^2)
!
!    exactly, because INT dw/2pi 1/((w_k^2-w^2)^2 + gamma^2 w^2) = 1/(2 gamma
!    w_k^2). With the quadrature partner carrying an equal share and the IR
!    convention I(nu) = nu^p S(w),
!
!      I_k = [ nu_k^p gamma w_k^2 ] R_k^2                                (9)
!
!    and that bracket is pref. Leaving it out does not rescale the spectrum --
!    it TILTS it by w^4 across the fitted range, which no single fitted scale
!    can absorb, and which looks like a spectrum rather than like a bug.
!
! 3. A FIRST-ORDER BACKWARD DIFFERENCE BREAKS THE QUADRATURE. mdot is estimated
!    from stored dipoles, and (m(t) - m(t-h))/h is the derivative at t - h/2,
!    not at t. Half a step of phase error between the two drives puts x and y
!    out of quadrature by w h / 2, which is a ripple of that size in R^2 and
!    half of it in the amplitude: 6% and 3% at w h = 0.13, measured. It looks
!    like noise and it is a systematic error. The three-point form
!    (3m(t) - 4m(t-h) + m(t-2h))/(2h) estimates the derivative AT t and takes
!    the ripple to O((w h)^2). Two previous frames are stored for it.
!
! 4. COHERENT VERSUS INCOHERENT: sum_i |x_i|^2 IS NOT AN IR SPECTRUM.
!    Absorption comes from the correlation of the TOTAL dipole, and the cross
!    terms between sites are transition-dipole coupling -- for water, the
!    reason the O-H stretch has the shape it has -- not noise. Summing the
!    squares site by site throws them away.
!
!      coherent (default)  one bank, driven by M = sum_i m_i
!      incoherent          one bank per site, driven by m_i, R^2 summed after
!
!    Because (1) is linear, the coherent bank driven by the total dipole is
!    exactly the sum of the per-site banks, so the coherent case needs no
!    per-site storage at all: n_modes resonators, not n_modes * n_atoms. The
!    weight (7) is then the same for every site, which is worth being explicit
!    about -- THE LOCALISATION BUYS NOTHING UNLESS THE TARGET IS LOCALISED.
!    Per-site amplitudes are only meaningful against a per-site target, which
!    is Improvement A below and does not exist yet. The incoherent mode is
!    provided because it is the half of that idea that can be built today.
!
! WHAT IS NOT IMPLEMENTED, AND WHERE IT WOULD GO
!
! IMPROVEMENT A, environment-aware targets. The experiment is one macroscopic
! curve; an incoherent bank is per site. Matching the sum against it lets the
! bias overdrive a few sites and leave others dead, and the fix is a per-site
! target -- weighting site i's contribution to the w_k restraint by how much
! its SOAP environment resembles the structure responsible for that band. That
! needs a mapping from a local environment to a spectral weight, which does not
! exist and is not something to invent silently. Where it would attach: I_exp
! and wgt become (n_modes, n_sites), and the incoherent branch of
! mad_ir_xl_evaluate accumulates the loss per site instead of summing R^2 over
! sites first. Nothing else in this file would change.
!
! COST
!
! Coherent: 96 * n_modes bytes, which is nothing, so ir_xl_n_modes can keep the
! whole experimental grid. Incoherent: 96 * n_modes * n_atoms bytes per rank,
! 46 MB at 64 modes and 2500 sites, which is why ir_xl_n_modes subsamples. The
! subsample is a SUBSET of the experimental grid, never an interpolation onto a
! grid of our choosing, for the reason mad_ir.f90 gives: interpolating the
! experiment invents structure between its points and then fits to it. There is
! no point asking for more modes than (nu_max - nu_min)/d(nu) either way -- the
! bank cannot resolve bins narrower than its own bandwidth.
!
! Related: mad_ir.f90 (the ACF bias, and the experimental grid this shares),
! gle.f90 (the same Markovian-embedding device applied to a thermostat).
module mad_ir_xl

   use kinds
   use mad_ir, only: mad_ir_type, CM_PER_INV_FS

   implicit none

   private
   public :: mad_ir_xl_type, mad_ir_xl_state
   public :: mad_ir_xl_init, mad_ir_xl_free, mad_ir_xl_setup
   public :: mad_ir_xl_evaluate, mad_ir_xl_advance, mad_ir_xl_weights
   public :: mad_ir_xl_ready, mad_ir_xl_resolution, mad_ir_xl_memory_bytes
   public :: mad_ir_xl_check, mad_ir_xl_pick_modes
   public :: mad_ir_xl_save, mad_ir_xl_load, mad_ir_xl_write_spectrum
   public :: mad_ir_xl_site_w, mad_ir_xl_active, mad_ir_xl_collect
   public :: mad_ir_xl_force
   public :: MAD_IR_XL_RESTART_VERSION

!  The restart format's own version, independent of mad_ir's. A bank is not a
!  history buffer and the two files are never interchangeable.
   integer, parameter :: MAD_IR_XL_RESTART_VERSION = 1

!  Run-wide state, for the same reason mad_ir_dmu_dr is run-wide: it is read
!  inside the descriptor pass, several batches deep in an argument list that is
!  already too long.
!
!  mad_ir_xl_site_w(1:3, 1:n_atoms) is the weight of equation (7) -- already
!  negated, so it is the FORCE weight and the descriptor pass adds
!  sum_a W_a dm_a/dr_jb straight into the force. Unlike the ACF bias's lambda
!  it exists before the pass, which is what lets it be contracted inside it.
!
!  mad_ir_xl_force(1:3, 1:n_atoms) is where the pass leaves the result. It is a
!  separate accumulator from the run's own forces array because gap_interface
!  is several batches and one MPI reduction away from the caller that owns
!  those forces, and threading them down would mean adding an argument to a
!  chain that already carries thirty.
   real(dp), allocatable, save :: mad_ir_xl_site_w(:, :)
   real(dp), allocatable, save :: mad_ir_xl_force(:, :)
   logical, save :: mad_ir_xl_active = .false.
!  Set once per step by the caller, exactly as mad_ir_collect is: is this a
!  step on which a stored frame is taken, and is a force wanted from it.
   logical, save :: mad_ir_xl_collect = .false.

   type :: mad_ir_xl_type
      logical :: active = .false.
!     sizing
      integer :: n_modes = 0
      integer :: n_sites = 0         ! dipole sites in the system
!     BANK REPLICAS. 1 in coherent mode -- one bank driven by the total dipole
!     is exactly the sum of the per-site banks, because (1) is linear -- and
!     n_sites in incoherent mode. Everything below is written against n_bank so
!     that the two modes are one code path.
      integer :: n_bank = 0
      real(dp) :: dt = 0.d0          ! fs between advances (md_step * ir_stride)
!     filter
      real(dp) :: tau_mem = 0.d0     ! fs
      real(dp) :: gamma = 0.d0       ! 1/fs, = 2/tau_mem
      logical :: coherent = .true.
!     target matching, same meaning as in mad_ir_type
      real(dp) :: nu_power = 2.d0
      logical :: match_scale = .true.
      logical :: match_offset = .false.
      real(dp) :: scale = 1.d0
      real(dp) :: offset = 0.d0
      real(dp) :: dissim = 0.d0
      real(dp) :: dissim_ref = 0.d0
!     progress and diagnostics
      integer :: n_steps = 0         ! advances applied
      integer :: n_warm = 0          ! advances to hold the bias off for
      integer :: n_prev = 0          ! previous frames held, for mdot (max 2)
!     Did the advance just performed actually use an mdot? Until two previous
!     frames exist the y bank is not driven at all, so the newest dipole does
!     not reach it and its share of the gradient is zero. n_prev alone cannot
!     answer this: it is bumped by the same advance, so a step that ran with
!     mdot = 0 leaves n_prev at 2 exactly like the step after it.
      logical :: mdot_live = .false.
      real(dp) :: w_rms = 0.d0       ! RMS of the weight, a proxy for force size
!     mode grid: a subset of the parent's fitted grid
      integer, allocatable :: kmap(:)      ! (n_modes) -> index in parent nu
      real(dp), allocatable :: nu(:)       ! (n_modes) cm^-1
      real(dp), allocatable :: w(:)        ! (n_modes) rad/fs
      real(dp), allocatable :: pref(:)     ! (n_modes) R^2 -> I, equation (9)
      real(dp), allocatable :: sens(:)     ! (n_modes) 2 P_k, equation (5)
!     The y bank's share of the same sensitivity: 3/(2 h w_k), the derivative
!     of the three-point mdot with respect to the newest dipole. Kept beside
!     sens rather than written out at the point of use, because the difference
!     formula in mad_ir_xl_advance and its derivative here have to agree and
!     the only way to keep them agreeing is to have one place each.
      real(dp), allocatable :: yfac(:)     ! (n_modes)
      real(dp), allocatable :: I_exp(:)    ! (n_modes)
      real(dp), allocatable :: wgt(:)      ! (n_modes)
      real(dp), allocatable :: I_calc(:)   ! (n_modes)
      real(dp), allocatable :: R2(:)       ! (n_modes)
      real(dp), allocatable :: dLdI(:)     ! (n_modes), set by evaluate
!     the bank itself, (3, n_modes, n_bank)
      real(dp), allocatable :: x(:, :, :)
      real(dp), allocatable :: xd(:, :, :)
      real(dp), allocatable :: y(:, :, :)
      real(dp), allocatable :: yd(:, :, :)
!     two previous drives, for the three-point mdot; (3, n_bank)
      real(dp), allocatable :: mu_p1(:, :)
      real(dp), allocatable :: mu_p2(:, :)
   end type mad_ir_xl_type

   type(mad_ir_xl_type), save :: mad_ir_xl_state

contains

!
! The resolution the bank actually has, in cm^-1.
!
! A Lorentzian of full width gamma in angular frequency is a line of width
! gamma/(2 pi c) in wavenumber. This is the number to compare with
! CM_PER_INV_FS/(n_lag*dt) for the block estimator; they are the same quantity,
! and a run should not be asked for a target grid much finer than either.
   real(dp) function mad_ir_xl_resolution(tau_mem)
      implicit none
      real(dp), intent(in) :: tau_mem
      if (tau_mem <= 0.d0) then
         mad_ir_xl_resolution = huge(1.d0)
      else
         mad_ir_xl_resolution = CM_PER_INV_FS/(dacos(-1.d0)*tau_mem)
      end if
   end function mad_ir_xl_resolution

!
! Bytes the bank occupies, per MPI rank: four (3, n_modes, n_bank) arrays.
   real(dp) function mad_ir_xl_memory_bytes(n_modes, n_bank)
      implicit none
      integer, intent(in) :: n_modes, n_bank
      mad_ir_xl_memory_bytes = 4.d0*3.d0*8.d0*dfloat(n_modes)*dfloat(n_bank)
   end function mad_ir_xl_memory_bytes

!
! Is this a combination the bank can be run with? Refuse rather than silently
! repair, and say which requirement failed.
!
! THE UNDERDAMPED CONDITION IS THE ONE THAT MATTERS. gamma < 2 w_k has to hold
! at the LOWEST fitted frequency, because gamma is common to the whole bank. A
! resonator with gamma >= 2 w is overdamped: its response has no peak, and "the
! amplitude at w_k" is not a spectral estimate of anything. In wavenumbers,
!
!   tau_mem > CM_PER_INV_FS / (2 pi nu_min)
!
! which is 53 fs at 100 cm^-1, so it bites only for a bank taken down towards
! zero frequency. It is checked anyway: the failure mode is a spectrum that
! looks smooth and means nothing.
   subroutine mad_ir_xl_check(tau_mem, n_modes, nu_min, dt, max_mem, n_bank, ok, msg)

      implicit none

      real(dp), intent(in) :: tau_mem, nu_min, dt, max_mem
      integer, intent(in) :: n_modes, n_bank
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: tau_floor, bytes

      ok = .false.
      msg = ""

      if (tau_mem <= 0.d0) then
         msg = "mad_ir_xl: ir_bias_mode = xl needs a positive ir_xl_tau_mem"
         return
      end if

!     A memory shorter than the sampling interval cannot be resolved by a
!     filter advanced once per interval. Same requirement, and same reason, as
!     ir_tau_mem's.
      if (tau_mem <= dt) then
         write (msg, '(A,F0.4,A,F0.4,A)') &
            "mad_ir_xl: ir_xl_tau_mem (", tau_mem, " fs) must exceed the interval "// &
            "between stored frames, md_step*ir_stride = ", dt, " fs"
         return
      end if

      if (nu_min <= 0.d0) then
         msg = "mad_ir_xl: the fitted range must start above zero wavenumber"
         return
      end if
      tau_floor = CM_PER_INV_FS/(2.d0*dacos(-1.d0)*nu_min)
      if (tau_mem <= tau_floor) then
         write (msg, '(A,F0.2,A,F0.2,A,F0.1,A)') &
            "mad_ir_xl: ir_xl_tau_mem (", tau_mem, " fs) is too short for the "// &
            "lowest fitted wavenumber: the bank is overdamped below ", tau_floor, &
            " fs at ", nu_min, " cm^-1"
         return
      end if

      if (n_modes < 1) then
         msg = "mad_ir_xl: ir_xl_n_modes must be at least 1"
         return
      end if

      bytes = mad_ir_xl_memory_bytes(n_modes, n_bank)
      if (max_mem > 0.d0 .and. bytes > max_mem) then
         write (msg, '(A,F0.2,A,F0.2,A)') &
            "mad_ir_xl: the resonator bank needs ", bytes/1.048576d6, &
            " MB per rank, over the ir_xl_max_memory limit of ", &
            max_mem/1.048576d6, " MB. Lower ir_xl_n_modes, use "// &
            "ir_xl_amplitude = coherent, or raise the limit."
         return
      end if

      ok = .true.

   end subroutine mad_ir_xl_check

!
! Choose which of the parent's fitted frequencies get resonators.
!
! Evenly spaced IN INDEX over the parent grid, endpoints included, and always a
! strict subset of it -- the experimental points themselves, never a
! resampling. n_modes_ask at or above n_freq keeps every point.
!
! Evenly spaced in index rather than in wavenumber because the parent grid is
! already whatever spacing the experiment came on, and weight_by_spacing has
! already accounted for it; respacing here would double-count that.
   subroutine mad_ir_xl_pick_modes(n_freq, n_modes_ask, kmap, n_modes)

      implicit none

      integer, intent(in) :: n_freq, n_modes_ask
      integer, allocatable, intent(out) :: kmap(:)
      integer, intent(out) :: n_modes
      integer :: m, j, jprev

      if (n_modes_ask >= n_freq .or. n_modes_ask <= 0) then
         n_modes = n_freq
         allocate (kmap(1:n_modes))
         do m = 1, n_modes
            kmap(m) = m
         end do
         return
      end if

      if (n_modes_ask == 1) then
         n_modes = 1
         allocate (kmap(1:1))
         kmap(1) = (n_freq + 1)/2
         return
      end if

!     Two passes: place, then drop duplicates. Rounding can land two requests
!     on the same experimental point when the grid is short, and a bank with
!     two identical resonators double-counts that point in the loss.
      allocate (kmap(1:n_modes_ask))
      n_modes = 0
      jprev = 0
      do m = 1, n_modes_ask
         j = 1 + nint(dfloat((m - 1)*(n_freq - 1))/dfloat(n_modes_ask - 1))
         if (j < 1) j = 1
         if (j > n_freq) j = n_freq
         if (j > jprev) then
            n_modes = n_modes + 1
            kmap(n_modes) = j
            jprev = j
         end if
      end do

   end subroutine mad_ir_xl_pick_modes

!
! Allocate the bank and precompute everything that depends only on the grid.
!
! parent supplies the fitted grid: nu, I_exp, wgt, nu_power, match_scale and
! match_offset all come from mad_ir_state, which read_exp_data has already
! filled and mad_ir_select_range has already restricted to [nu_min, nu_max].
! Duplicating any of that here would be a second reader of the same file, which
! mad_ir.f90 is explicit about not wanting.
   subroutine mad_ir_xl_init(this, parent, n_sites, n_modes_ask, dt, tau_mem, &
                             coherent, warm_factor)

      implicit none

      type(mad_ir_xl_type), intent(inout) :: this
      type(mad_ir_type), intent(in) :: parent
      integer, intent(in) :: n_sites, n_modes_ask
      real(dp), intent(in) :: dt, tau_mem
      logical, intent(in) :: coherent
      real(dp), intent(in) :: warm_factor
      integer :: m, k
      real(dp) :: two_pi, lam, Om, w2

      call mad_ir_xl_free(this)

      two_pi = 2.d0*dacos(-1.d0)

      call mad_ir_xl_pick_modes(parent%n_freq, n_modes_ask, this%kmap, this%n_modes)

      this%n_sites = n_sites
      this%coherent = coherent
      if (coherent) then
         this%n_bank = 1
      else
         this%n_bank = n_sites
      end if
      this%dt = dt
      this%tau_mem = tau_mem
      this%gamma = 2.d0/tau_mem
      this%nu_power = parent%nu_power
      this%match_scale = parent%match_scale
      this%match_offset = parent%match_offset

      allocate (this%nu(1:this%n_modes), this%w(1:this%n_modes))
      allocate (this%pref(1:this%n_modes), this%sens(1:this%n_modes))
      allocate (this%yfac(1:this%n_modes))
      allocate (this%I_exp(1:this%n_modes), this%wgt(1:this%n_modes))
      allocate (this%I_calc(1:this%n_modes), this%R2(1:this%n_modes))
      allocate (this%dLdI(1:this%n_modes))

      lam = 0.5d0*this%gamma
      do m = 1, this%n_modes
         k = this%kmap(m)
         this%nu(m) = parent%nu(k)
         this%I_exp(m) = parent%I_exp(k)
         this%wgt(m) = parent%wgt(k)
!        rad/fs from cm^-1
         this%w(m) = two_pi*this%nu(m)/CM_PER_INV_FS
         w2 = this%w(m)**2
!        equation (9)
         this%pref(m) = (this%nu(m)**this%nu_power)*this%gamma*w2
!        equation (5), times the 2 of equation (6): the sensitivity of the
!        resonator to one interval of drive. Precomputed because it depends
!        only on the mode and the step, and because writing it out at every use
!        is how the propagator and its derivative drift apart.
         Om = dsqrt(max(1.d-300, w2 - lam*lam))
         this%sens(m) = 2.d0*(1.d0/w2)*(1.d0 - dexp(-lam*dt)* &
                                        (dcos(Om*dt) + lam*dsin(Om*dt)/Om))
!        d(mdot)/dm(t) = 3/(2h) for the three-point backward difference, over
!        the w_k that scales the y drive.
         this%yfac(m) = 1.5d0/(dt*this%w(m))
      end do

      allocate (this%x(1:3, 1:this%n_modes, 1:this%n_bank))
      allocate (this%xd(1:3, 1:this%n_modes, 1:this%n_bank))
      allocate (this%y(1:3, 1:this%n_modes, 1:this%n_bank))
      allocate (this%yd(1:3, 1:this%n_modes, 1:this%n_bank))
      allocate (this%mu_p1(1:3, 1:this%n_bank), this%mu_p2(1:3, 1:this%n_bank))

      this%x = 0.d0
      this%xd = 0.d0
      this%y = 0.d0
      this%yd = 0.d0
      this%mu_p1 = 0.d0
      this%mu_p2 = 0.d0
      this%I_calc = 0.d0
      this%R2 = 0.d0
      this%dLdI = 0.d0

      this%n_steps = 0
      this%n_prev = 0
      this%mdot_live = .false.
      this%scale = 1.d0
      this%offset = 0.d0

!     THE BANK STARTS AT REST AND HAS TO CHARGE UP. Every resonator approaches
!     its stationary amplitude with the envelope 1 - exp(-gamma t/2), which is
!     the SAME at every frequency because gamma is -- so the fill deficit is a
!     pure overall factor and ir_match_scale absorbs it exactly, for the same
!     reason it absorbs the exponential ACF estimator's. What it does not
!     absorb is the transient ringing at w_k left by the initial condition,
!     which is not common across the bank. Holding the bias off for a few
!     memory times lets that die.
      this%n_warm = max(2, nint(warm_factor*tau_mem/dt))

      this%active = .true.

   end subroutine mad_ir_xl_init

   subroutine mad_ir_xl_free(this)
      implicit none
      type(mad_ir_xl_type), intent(inout) :: this
      if (allocated(this%kmap)) deallocate (this%kmap)
      if (allocated(this%nu)) deallocate (this%nu)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%pref)) deallocate (this%pref)
      if (allocated(this%sens)) deallocate (this%sens)
      if (allocated(this%yfac)) deallocate (this%yfac)
      if (allocated(this%I_exp)) deallocate (this%I_exp)
      if (allocated(this%wgt)) deallocate (this%wgt)
      if (allocated(this%I_calc)) deallocate (this%I_calc)
      if (allocated(this%R2)) deallocate (this%R2)
      if (allocated(this%dLdI)) deallocate (this%dLdI)
      if (allocated(this%x)) deallocate (this%x)
      if (allocated(this%xd)) deallocate (this%xd)
      if (allocated(this%y)) deallocate (this%y)
      if (allocated(this%yd)) deallocate (this%yd)
      if (allocated(this%mu_p1)) deallocate (this%mu_p1)
      if (allocated(this%mu_p2)) deallocate (this%mu_p2)
      this%active = .false.
      this%n_modes = 0
      this%n_sites = 0
      this%n_bank = 0
      this%n_steps = 0
      this%n_prev = 0
   end subroutine mad_ir_xl_free

!
! Has the bank charged enough for its amplitudes to mean anything?
   logical function mad_ir_xl_ready(this)
      implicit none
      type(mad_ir_xl_type), intent(in) :: this
      mad_ir_xl_ready = this%active .and. (this%n_steps >= this%n_warm)
   end function mad_ir_xl_ready

!
! Advance the bank by one stored frame, driven by the dipoles of the current
! configuration.
!
! THE INTEGRATOR IS EXACT FOR THE LINEAR PART. Over one interval h the drive is
! held constant -- a zero-order hold, the standard impulse-invariant
! discretisation -- and the homogeneous part is propagated by its closed form
! rather than by a Verlet step:
!
!   lam = gamma/2,  Om = sqrt(w^2 - lam^2),  E = exp(-lam h)
!   xp  = F / w^2
!   u0  = x0 - xp
!   x1  = xp + E [ u0 cos(Om h) + (v0 + lam u0) sin(Om h)/Om ]
!   v1  =      E [ v0 cos(Om h) - (w^2 u0 + lam v0) sin(Om h)/Om ]
!
! This is unconditionally stable at any h, which matters: w h at 4000 cm^-1 and
! a 5 fs sampling interval is 3.8, and velocity Verlet has been unstable since
! w h = 2. It is also the same choice mad_ir.f90 makes for alpha -- the exact
! coefficient of the ODE over the interval, not its Euler approximation -- for
! the same reason. The sensitivity used by the gradient, equation (5), is the
! derivative of exactly this expression.
!
! ORDERING. This is called AFTER the descriptor pass of the frame it is given,
! and mad_ir_xl_evaluate and mad_ir_xl_weights follow it. The weight the next
! frame contracts is therefore formed from a bank that already knows this
! frame's dipole; see the header's note on the one frame of lag.
   subroutine mad_ir_xl_advance(this, mu_site)

      implicit none

      type(mad_ir_xl_type), intent(inout) :: this
      real(dp), intent(in) :: mu_site(:, :)   ! (3, n_sites), local dipoles
      real(dp) :: lam, h, w2, Om, Ec, cs, snOm
      real(dp) :: Fx(1:3), mdot(1:3)
      real(dp) :: xp, u0, v0, invw2
      real(dp), allocatable :: drive(:, :)
      integer :: m, i, a

      if (.not. this%active) return

      h = this%dt
      lam = 0.5d0*this%gamma

!     ---- the drive -------------------------------------------------------
!     In coherent mode the single bank is driven by the total dipole, which is
!     exactly what the sum of the per-site banks would give because (1) is
!     linear. Forming it here rather than asking the caller for it keeps the
!     two modes on one interface.
      allocate (drive(1:3, 1:this%n_bank))
      if (this%coherent) then
         do a = 1, 3
            drive(a, 1) = sum(mu_site(a, 1:this%n_sites))
         end do
      else
         drive(1:3, 1:this%n_bank) = mu_site(1:3, 1:this%n_bank)
      end if

      do m = 1, this%n_modes

         w2 = this%w(m)**2
         invw2 = 1.d0/w2
         Om = dsqrt(w2 - lam*lam)
         Ec = dexp(-lam*h)
         cs = dcos(Om*h)
         snOm = dsin(Om*h)/Om

         do i = 1, this%n_bank

!           ---- mdot, estimated AT t and not at t - h/2 ------------------
!           The three-point backward form. With the two-point one the y bank's
!           drive represents a different instant from the x bank's, and half a
!           step of phase error puts them out of quadrature by w h / 2 -- a
!           systematic ripple in R^2 of that size, which looks like noise. Two
!           previous frames have to exist for it, and until they do the y bank
!           is simply not driven: an approximate derivative at the start of a
!           run that is about to be discarded as warm-up is not worth having.
            if (this%n_prev >= 2) then
               mdot(1:3) = (3.d0*drive(1:3, i) - 4.d0*this%mu_p1(1:3, i) &
                            + this%mu_p2(1:3, i))/(2.d0*h)
            else
               mdot = 0.d0
            end if

            Fx(1:3) = drive(1:3, i)
            do a = 1, 3
               xp = Fx(a)*invw2
               u0 = this%x(a, m, i) - xp
               v0 = this%xd(a, m, i)
               this%x(a, m, i) = xp + Ec*(u0*cs + (v0 + lam*u0)*snOm)
               this%xd(a, m, i) = Ec*(v0*cs - (w2*u0 + lam*v0)*snOm)

               xp = (mdot(a)/this%w(m))*invw2
               u0 = this%y(a, m, i) - xp
               v0 = this%yd(a, m, i)
               this%y(a, m, i) = xp + Ec*(u0*cs + (v0 + lam*u0)*snOm)
               this%yd(a, m, i) = Ec*(v0*cs - (w2*u0 + lam*v0)*snOm)
            end do

         end do
      end do

      this%mdot_live = (this%n_prev >= 2)
      this%mu_p2(1:3, 1:this%n_bank) = this%mu_p1(1:3, 1:this%n_bank)
      this%mu_p1(1:3, 1:this%n_bank) = drive(1:3, 1:this%n_bank)
      if (this%n_prev < 2) this%n_prev = this%n_prev + 1
      this%n_steps = this%n_steps + 1
      deallocate (drive)

   end subroutine mad_ir_xl_advance

!
! Amplitudes, predicted spectrum, loss, and dU/dI.
!
! Called once per stored frame, immediately after mad_ir_xl_advance, so that
! the amplitudes it reads are the ones that already know the current dipole.
!
! energy is the bias energy, with the same 1/2 and the same energy_scale that
! get_exp_energies produces for every other observable. Nothing here touches
! the atoms; dLdI is what mad_ir_xl_weights turns into a force.
   subroutine mad_ir_xl_evaluate(this, energy_scale, energy)

      implicit none

      type(mad_ir_xl_type), intent(inout) :: this
      real(dp), intent(in) :: energy_scale
      real(dp), intent(out) :: energy
      real(dp) :: swi2, swi, sw, swie, swe, det, s_fit, b_fit, c, acc
      integer :: m, i, a

      energy = 0.d0
      this%dLdI = 0.d0
      this%dissim = 0.d0
      this%dissim_ref = 0.d0
      if (.not. this%active) return

!     ---- amplitudes -------------------------------------------------
!     One bank in coherent mode, so the sum over i is a sum over one thing and
!     the distinction lives entirely in what the bank was driven by.
      do m = 1, this%n_modes
         acc = 0.d0
         do i = 1, this%n_bank
            do a = 1, 3
               acc = acc + this%x(a, m, i)**2 + this%y(a, m, i)**2
            end do
         end do
         this%R2(m) = acc
         this%I_calc(m) = this%pref(m)*acc
      end do

      if (.not. mad_ir_xl_ready(this)) return

!     ---- overall scale, and baseline if asked for --------------------
!     Identical in form, and in the negative-scale guard, to
!     mad_ir_evaluate's: a negative scale reverses the bias, driving the model
!     away from the bands it is being asked to grow, so the two-parameter solve
!     falls back to the scale-only one whenever it returns s <= 0.
      s_fit = 1.d0
      b_fit = 0.d0
      if (this%match_scale) then
         swi2 = 0.d0; swi = 0.d0; sw = 0.d0; swie = 0.d0; swe = 0.d0
         do m = 1, this%n_modes
            swi2 = swi2 + this%wgt(m)*this%I_calc(m)**2
            swi = swi + this%wgt(m)*this%I_calc(m)
            sw = sw + this%wgt(m)
            swie = swie + this%wgt(m)*this%I_calc(m)*this%I_exp(m)
            swe = swe + this%wgt(m)*this%I_exp(m)
         end do
         if (this%match_offset) then
            det = swi2*sw - swi*swi
            if (dabs(det) > 1.d-300) then
               s_fit = (swie*sw - swe*swi)/det
               b_fit = (swi2*swe - swi*swie)/det
            else if (swi2 > 0.d0) then
               s_fit = swie/swi2
            end if
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

!     ---- loss, and its sensitivity to the predicted spectrum ---------
      do m = 1, this%n_modes
         c = s_fit*this%I_calc(m) + b_fit - this%I_exp(m)
         energy = energy + 0.5d0*energy_scale*this%wgt(m)*c**2
         this%dLdI(m) = energy_scale*this%wgt(m)*s_fit*c
         this%dissim = this%dissim + this%wgt(m)*c**2
         this%dissim_ref = this%dissim_ref + this%wgt(m)*this%I_exp(m)**2
      end do

   end subroutine mad_ir_xl_evaluate

!
! The per-site weight the descriptor pass contracts against d m_i / d r_j:
! equation (7), already negated, so it is the force weight and the pass adds.
!
!   W_{i,a} = - sum_k dLdI_k pref_k sens_k [ x_{k,i,a} + y_{k,a}/(h w_k) ]
!
! In coherent mode the bank has one replica and the weight is the same for
! every site -- which is the honest statement that localisation buys nothing
! against a macroscopic target. In incoherent mode each site gets its own.
!
! Zero while the bank is charging, so that a warming step applies exactly no
! force rather than a transient one whose size depends on when the run started.
! The gate is here, once, rather than in every consumer.
   subroutine mad_ir_xl_weights(this, w_out)

      implicit none

      type(mad_ir_xl_type), intent(inout) :: this
      real(dp), intent(out) :: w_out(:, :)     ! (3, n_sites)
      real(dp) :: cf, s2, yf, acc(1:3)
      integer :: i, i2, m, a

      w_out = 0.d0
      this%w_rms = 0.d0
      if (.not. mad_ir_xl_ready(this)) return

      do i = 1, this%n_bank
         acc = 0.d0
         do m = 1, this%n_modes
            cf = this%dLdI(m)*this%pref(m)*this%sens(m)
!           The y bank contributes only if the advance that just ran actually
!           had an mdot to drive it with.
            if (this%mdot_live) then
               yf = this%yfac(m)
            else
               yf = 0.d0
            end if
            do a = 1, 3
               acc(a) = acc(a) + cf*(this%x(a, m, i) + yf*this%y(a, m, i))
            end do
         end do
         if (this%coherent) then
!           One bank, and dX/dm_i is the same P for every site because the
!           drive is the SUM of the local dipoles. So is the weight.
            do i2 = 1, this%n_sites
               w_out(1:3, i2) = -acc(1:3)
            end do
         else
            w_out(1:3, i) = -acc(1:3)
         end if
      end do

      s2 = 0.d0

      do i = 1, this%n_sites
         do a = 1, 3
            s2 = s2 + w_out(a, i)**2
         end do
      end do
      if (this%n_sites > 0) this%w_rms = dsqrt(s2/dfloat(this%n_sites))

   end subroutine mad_ir_xl_weights

!
! Set the bank up from the parent observable and, if there is one, a restart.
!
! resumed says whether a saved bank was adopted. A refused or missing one is
! not fatal: the run charges a fresh bank, which costs n_warm stored frames of
! unbiased dynamics, and the caller should say so.
   subroutine mad_ir_xl_setup(parent, n_sites, n_modes_ask, dt_md, stride, &
                              tau_mem, coherent, warm_factor, max_mem, &
                              restart_file, ok, resumed, msg)

      implicit none

      type(mad_ir_type), intent(in) :: parent
      integer, intent(in) :: n_sites, n_modes_ask, stride
      real(dp), intent(in) :: dt_md, tau_mem, warm_factor, max_mem
      logical, intent(in) :: coherent
      character(len=*), intent(in) :: restart_file
      logical, intent(out) :: ok, resumed
      character(len=*), intent(out) :: msg
      real(dp) :: dt, nu_lo
      integer :: n_modes_eff, n_bank_eff
      integer, allocatable :: kmap_try(:)
      character(len=512) :: msg2

      ok = .false.
      resumed = .false.
      msg = ""

      if (.not. parent%active) then
         msg = "mad_ir_xl: the ir observable must be set up before its bank"
         return
      end if
      if (parent%n_freq < 1) then
         msg = "mad_ir_xl: the fitted grid is empty"
         return
      end if
      if (n_sites < 1) then
         msg = "mad_ir_xl: no dipole sites to drive the bank"
         return
      end if
      if (stride < 1) then
         msg = "mad_ir_xl: ir_stride must be at least 1"
         return
      end if

      dt = dt_md*dfloat(stride)
      nu_lo = minval(parent%nu(1:parent%n_freq))

!     Size first, so that the memory check sees the number of modes that will
!     actually be allocated rather than the number asked for.
      call mad_ir_xl_pick_modes(parent%n_freq, n_modes_ask, kmap_try, n_modes_eff)
      deallocate (kmap_try)
      if (coherent) then
         n_bank_eff = 1
      else
         n_bank_eff = n_sites
      end if

      call mad_ir_xl_check(tau_mem, n_modes_eff, nu_lo, dt, max_mem, n_bank_eff, ok, msg)
      if (.not. ok) return

      call mad_ir_xl_init(mad_ir_xl_state, parent, n_sites, n_modes_ask, dt, &
                          tau_mem, coherent, warm_factor)

      if (allocated(mad_ir_xl_site_w)) deallocate (mad_ir_xl_site_w)
      allocate (mad_ir_xl_site_w(1:3, 1:n_sites))
      mad_ir_xl_site_w = 0.d0
      if (allocated(mad_ir_xl_force)) deallocate (mad_ir_xl_force)
      allocate (mad_ir_xl_force(1:3, 1:n_sites))
      mad_ir_xl_force = 0.d0

      if (trim(restart_file) /= "none") then
         call mad_ir_xl_load(mad_ir_xl_state, restart_file, resumed, msg2)
         if (.not. resumed) msg = trim(msg2)
      end if

      ok = .true.

   end subroutine mad_ir_xl_setup

!
! Persist the bank.
!
! Unformatted, and carrying every parameter the state depends on, for the
! reason mad_ir_save gives: x is a filtered dipole with a particular tau_mem,
! sampling interval and mode grid baked into it, and adopting it under
! different ones would change the observable mid-run in a way nothing
! downstream could notice. mad_ir_xl_load refuses each mismatch by name.
   subroutine mad_ir_xl_save(this, fname, ok, msg)

      implicit none

      type(mad_ir_xl_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: u, ios

      ok = .false.
      msg = ""
      if (.not. this%active) then
         msg = "mad_ir_xl_save: nothing to save"
         return
      end if

      open (newunit=u, file=trim(fname), form="unformatted", status="replace", &
            action="write", iostat=ios)
      if (ios /= 0) then
         msg = "mad_ir_xl_save: cannot open "//trim(fname)
         return
      end if

      write (u) MAD_IR_XL_RESTART_VERSION
      write (u) this%n_modes, this%n_bank, this%n_sites
      write (u) this%dt, this%tau_mem, this%nu_power
      write (u) this%coherent
      write (u) this%n_steps, this%n_prev, this%mdot_live
      write (u) this%nu(1:this%n_modes)
      write (u) this%x
      write (u) this%xd
      write (u) this%y
      write (u) this%yd
      write (u) this%mu_p1
      write (u) this%mu_p2
      close (u)

      ok = .true.

   end subroutine mad_ir_xl_save

   subroutine mad_ir_xl_load(this, fname, ok, msg)

      implicit none

      type(mad_ir_xl_type), intent(inout) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: u, ios, ver, nm, nb, ns, nst, npv
      logical :: mlive
      real(dp) :: dtf, tauf, powf
      logical :: coh, exists
      real(dp), allocatable :: nuf(:)

      ok = .false.
      msg = ""

      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
         msg = "mad_ir_xl: no bank restart file "//trim(fname)//"; charging a fresh bank"
         return
      end if

      open (newunit=u, file=trim(fname), form="unformatted", status="old", &
            action="read", iostat=ios)
      if (ios /= 0) then
         msg = "mad_ir_xl: cannot read "//trim(fname)//"; charging a fresh bank"
         return
      end if

      read (u, iostat=ios) ver
      if (ios /= 0 .or. ver /= MAD_IR_XL_RESTART_VERSION) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" is not a bank restart of this version; "// &
               "charging a fresh bank"
         return
      end if

      read (u, iostat=ios) nm, nb, ns
      if (ios /= 0) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" is truncated; charging a fresh bank"
         return
      end if
      if (nm /= this%n_modes .or. nb /= this%n_bank .or. ns /= this%n_sites) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" was written for a different bank size; "// &
               "charging a fresh bank"
         return
      end if

      read (u, iostat=ios) dtf, tauf, powf
      if (ios /= 0) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" is truncated; charging a fresh bank"
         return
      end if
!     A bank filtered with a different memory or sampling interval is a
!     different object wearing the same name.
      if (dabs(dtf - this%dt) > 1.d-9*max(1.d0, dabs(this%dt)) .or. &
          dabs(tauf - this%tau_mem) > 1.d-9*max(1.d0, dabs(this%tau_mem)) .or. &
          dabs(powf - this%nu_power) > 1.d-9*max(1.d0, dabs(this%nu_power))) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" was written with a different "// &
               "ir_stride, ir_xl_tau_mem or ir_nu_power; charging a fresh bank"
         return
      end if

      read (u, iostat=ios) coh
      if (ios /= 0) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" is truncated; charging a fresh bank"
         return
      end if
!     .neqv. binds looser than .and., so this comparison gets its own statement
!     rather than being folded into the read's iostat test.
      if (coh .neqv. this%coherent) then
         close (u)
         msg = "mad_ir_xl: "//trim(fname)//" was written with the other "// &
               "ir_xl_amplitude; charging a fresh bank"
         return
      end if

      allocate (nuf(1:nm))
      read (u, iostat=ios) nst, npv, mlive
      if (ios == 0) read (u, iostat=ios) nuf
      if (ios /= 0) then
         close (u)
         deallocate (nuf)
         msg = "mad_ir_xl: "//trim(fname)//" is truncated; charging a fresh bank"
         return
      end if
      if (any(dabs(nuf - this%nu(1:nm)) > 1.d-6*max(1.d0, maxval(dabs(this%nu))))) then
         close (u)
         deallocate (nuf)
         msg = "mad_ir_xl: "//trim(fname)//" was written on a different mode grid; "// &
               "charging a fresh bank"
         return
      end if
      deallocate (nuf)

      read (u, iostat=ios) this%x
      if (ios == 0) read (u, iostat=ios) this%xd
      if (ios == 0) read (u, iostat=ios) this%y
      if (ios == 0) read (u, iostat=ios) this%yd
      if (ios == 0) read (u, iostat=ios) this%mu_p1
      if (ios == 0) read (u, iostat=ios) this%mu_p2
      close (u)
      if (ios /= 0) then
         this%x = 0.d0; this%xd = 0.d0; this%y = 0.d0; this%yd = 0.d0
         this%mu_p1 = 0.d0; this%mu_p2 = 0.d0
         msg = "mad_ir_xl: "//trim(fname)//" is truncated; charging a fresh bank"
         return
      end if

      this%n_steps = nst
      this%n_prev = npv
      this%mdot_live = mlive
!     n_warm is NOT adopted: it is derived from this run's tau_mem and dt,
!     which the checks above have already established are the file's too.

      ok = .true.

   end subroutine mad_ir_xl_load

!
! The predicted spectrum on the mode grid, with the experiment beside it.
!
! The fitted scale and offset are applied, so the two columns are directly
! comparable; the raw R^2 is kept as a fourth column because it is the quantity
! the bias acts on, and a band that is not growing is diagnosed there rather
! than in the scaled spectrum.
   subroutine mad_ir_xl_write_spectrum(this, fname, with_exp)

      implicit none

      type(mad_ir_xl_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(in) :: with_exp
      integer :: u, ios, m

      if (.not. this%active) return

      open (newunit=u, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) return

      write (u, '(A)') "# auxiliary-variable (resonator bank) IR prediction"
      if (this%coherent) then
         write (u, '(A,I0,A,I0,A)') "# modes: ", this%n_modes, "   sites: ", &
            this%n_sites, "   amplitude: coherent"
      else
         write (u, '(A,I0,A,I0,A)') "# modes: ", this%n_modes, "   sites: ", &
            this%n_sites, "   amplitude: incoherent"
      end if
      write (u, '(A,F0.4,A,F0.4,A)') "# tau_mem: ", this%tau_mem, " fs   resolution: ", &
         mad_ir_xl_resolution(this%tau_mem), " cm^-1"
      write (u, '(A,ES16.8,A,ES16.8)') "# fitted scale: ", this%scale, &
         "   offset: ", this%offset
      write (u, '(A,I0,A,I0)') "# advances: ", this%n_steps, "   warm-up: ", this%n_warm
      if (with_exp) then
         write (u, '(A)') "#      nu (cm^-1)        I_pred (scaled)          I_exp                 R^2"
      else
         write (u, '(A)') "#      nu (cm^-1)        I_pred (scaled)                                R^2"
      end if
      do m = 1, this%n_modes
         if (with_exp) then
            write (u, '(4ES22.10)') this%nu(m), this%scale*this%I_calc(m) + this%offset, &
               this%I_exp(m), this%R2(m)
         else
            write (u, '(2ES22.10,22X,ES22.10)') this%nu(m), &
               this%scale*this%I_calc(m) + this%offset, this%R2(m)
         end if
      end do
      close (u)

   end subroutine mad_ir_xl_write_spectrum

end module mad_ir_xl
