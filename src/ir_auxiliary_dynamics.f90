! Envelope-targeted auxiliary dynamics for MAD IR: a driven resonator bank whose
! own amplitude is held at the experimental one by feedback, and which pushes the
! atoms only through the dipole coupling.
!
! Selected with ir_bias_mode = "aux".
!
!==========================================================================
! THE SCHEME
!==========================================================================
!
! To every fitted frequency w_k attach a 3-vector resonator (X_k, P_k) with
! fictitious mass mu_k, driven by the total dipole M(q):
!
!   Xdot_k = P_k/mu_k
!   Pdot_k = -mu_k w_k^2 X_k + g_k M(q) - eta_k P_k                      (1)
!
! The quantity the spectrum is read from is the phase-space amplitude
!
!   R_k^2 = |X_k|^2 + |P_k|^2/(mu_k^2 w_k^2)                             (2)
!
! which is the slowly varying ENVELOPE of the response at w_k: for a free
! resonator it is exactly constant, where |X_k|^2 alone would oscillate at
! 2 w_k. R_k is to (2) what mad_ir_xl's |x|^2 + |y|^2 is to its quadrature
! pair, reached by using the canonical momentum as the quadrature partner
! instead of integrating a second bank.
!
! The atoms feel the bank ONLY through the coupling term -g_k X_k.M(q):
!
!   F_i = F_MLP,i + escale sum_k g_k (grad_i m_i)^T X_k                  (3)
!
! which is local in i even though X_k couples to the total dipole -- the
! back-reaction needs only atom i's own dipole Jacobian, so (3) is O(N).
! escale is exp_energy_scales(ir_idx), ramped as for every other observable.
! The energy reported alongside (3) is -escale sum_k g_k X_k.M, which IS the
! generator of (3); unlike the ACF bias there is no frozen-past caveat here.
!
!==========================================================================
! WHY THE RESTRAINT IS FEEDBACK AND NOT A POTENTIAL
!==========================================================================
!
! The obvious way to hold R_k at a target is to add a bias potential
! U(R_k) = 1/2 gamma (R_k - R_target)^2 to the extended Hamiltonian and let
! Hamilton's equations do the rest. That is what an earlier draft of this file
! did, and it CANNOT work. The proof is one line in action-angle variables.
!
! R_k^2 is proportional to the action J_k of the resonator, so U(R_k) = U(J_k),
! and therefore
!
!   Jdot_k    = -dU/dphi_k = 0        IDENTICALLY                        (4)
!   phidot_k  = w_k + U'(J_k)
!
! A function of the action alone generates shifts in the ANGLE and never in the
! action. And the action is exactly what the spectrum reads out. Equivalently,
! in the complex quadrature z_k = X_k + i P_k/(mu_k w_k), for which
! |z_k| = R_k exactly, the canonical equations with U(R) give
!
!   zdot_k = -i w_k beta_k z_k ,     beta_k = 1 + U'(R_k)/(mu_k w_k^2 R_k)  (5)
!
! with beta_k REAL: |z_k| is exactly conserved and the bias is a pure frequency
! modulation w_k -> w_k beta_k. Measured, in tests/ir_aux/iraux_verify.f90
! check 1: an undriven resonator started at R = 1 stays at R = 1.000000 for
! targets of 0.01, 1 and 100.
!
! With the drive on, the only channel left is that frequency shift, and it
! points the wrong way in BOTH directions. Since U' = gamma(R - R_t),
!
!   w_eff^2 = w_k^2 + gamma (R_k - R_t)/(mu_k R_k)
!
! so a band that is too weak SOFTENS and a band that is too strong STIFFENS,
! and each channel is tuned to its own drive by construction. The driven
! response
!
!   |X| = (g/mu)|M| / sqrt((w_eff^2 - w_d^2)^2 + Gamma^2 w_d^2)
!
! is PEAKED at zero detuning, so either sign of the error reduces the response.
! No choice of sign or magnitude of gamma repairs this, because the response is
! not monotone in the shift: it has a maximum at zero. Driven on resonance the
! canonical form returns R = 1.0 for targets 0.5, 2, 8 and 32 alike -- zero
! discriminating power. This is the same result mad_ir_xl.f90 records having
! measured for its own quadratic-in-x bias; (4) shows it was never specific to
! that form, but holds for ANY U that is a function of R alone.
!
! An amplitude cannot be steered by a potential. It can be steered by friction.
!
!==========================================================================
! THE CONTROLLER
!==========================================================================
!
! eta_k in (1) is a per-mode friction that is free to go negative, integrating
! the envelope error, plus a proportional term that damps the loop:
!
!   etadot_k = (1/tau_k^2) (R_k^2/R_target_k^2 - 1)                      (6)
!
!   plus, applied to X_k and P_k with the SAME scalar,
!     c_k = -kappa_k (R_k - R_target_k)/R_k                              (7)
!
! (7) is a radial flow in the (X, P) plane. Applying one scalar to both
! components gives, exactly,
!
!   zdot_k = (c_k - i w_k) z_k    =>   dR_k/dt = c_k R_k                 (8)
!
! -- amplitude rescaled, frequency EXACTLY unchanged, no O(c^2) shift at all.
! It is the precise complement of (5): a potential in R changes the frequency
! and never the amplitude, a radial friction changes the amplitude and never
! the frequency. Amplitude is what the spectrum reads out, and the frequency
! must not move or channel k stops reporting the intensity at w_k and the
! calibration below refers to the wrong frequency.
!
! Together the envelope obeys
!
!   d(R_k^2)/dt = 2 g_k (M.P_k)/(mu_k^2 w_k^2)   injected power
!                 - 2 eta_k |P_k|^2/(mu_k^2 w_k^2)   integral feedback
!                 + 2 c_k R_k^2                      proportional feedback
!
! WHY BOTH TERMS. eta alone is a pure integral controller and has no damping.
! Cycle-averaging (<|P|^2> = mu^2 w^2 R^2/2, so <d(R^2)/dt> = -eta R^2 when
! undriven) and linearising R^2 = R_t^2 (1 + nu) gives
!
!   nudot = -eta,  etadot = nu/tau^2      =>      nu'' = -nu/tau^2
!
! which is MARGINALLY STABLE: it rings forever, with period exactly 2 pi tau.
! Measured (check 2): 1573, 3148, 6293 and 12600 fs against 2 pi tau of 1571,
! 3142, 6283 and 12566. Undriven from R = 1 towards a target of 4, R swung
! between 1.02 and 8.27 for 40 ps with no decay. This is the classical
! Nose-Hoover ringing problem, and it is worse here than in a thermostat for
! two reasons: R_k IS the observable, so the ringing goes straight into the
! reported spectrum with a 12-19 ps period for tau of a few ps and with each
! channel at its own phase; and the only damping the driven case supplies comes
! from the drive's own amplitude dependence, which vanishes for weakly driven
! channels -- exactly the channels with R_k < R_t that the bias exists for.
!
! Adding (7) makes the loop
!
!   nu'' + kappa nudot + nu/tau^2 = 0                                    (9)
!
! CRITICALLY DAMPED AT kappa_k = 2/tau_k, which is the default. Measured, as
! the late-time envelope band about a target of 4 for tau = 500 fs:
!
!   kappa*tau/2 = 0.00   3.3156 .. 4.6834    rings forever
!   kappa*tau/2 = 0.12   3.9028 .. 4.0651    ringing
!   kappa*tau/2 = 0.50   3.9999 .. 4.0003    no ring
!   kappa*tau/2 = 1.00   4.0000 .. 4.0000    critical
!   kappa*tau/2 = 10.0   3.9989 .. 3.9999    no ring
!
! And the integral term is kept because the proportional one alone has a
! steady-state droop of Pi/kappa against an injected power Pi: measured
! +0.260, +0.062 and +0.013 for kappa = 0.005, 0.02 and 0.08. With both, the
! steady-state error is exactly zero and the approach is monotone.
!
! H_ext is NOT conserved under (6)-(7), and cannot be: holding a mode at a
! target amplitude means pumping it. ir_aux_energy_pumped carries the running
! total for thermo.log, the way the GLE thermostat's ledger does.
!
!==========================================================================
! CALIBRATION
!==========================================================================
!
! Run unbiased first. Once mad_ir's ensemble is full its autocorrelation gives
! the unconstrained dipole power spectral density, and linear response fixes
! the couplings. For (1) with eta -> Gamma_k and no controller,
!
!   X_k(w) = chi_k(w) M(w),   chi_k(w) = (g_k/mu_k)/(w_k^2 - w^2 + i w Gamma_k)
!
!   <R_k^2> = INT dw/2pi |chi_k(w)|^2 (1 + w^2/w_k^2) S_MM(w)
!           ~= (g_k^2/(mu_k^2 w_k^2 Gamma_k)) S_MM(w_k)     sharp Gamma
!
! so that R_k ~ R_target when the trajectory already matches experiment:
!
!   g_k = mu_k w_k sqrt( Gamma_k <R_target_k^2> / S_MM(w_k) )           (10)
!
! S_MM(w_k) comes straight from mad_ir: its I_calc(k) is nu_k^p dt FT[win*acf],
! so S_MM at channel k is I_calc(k)/nu_k^p. No second transform is needed and
! no frequency-indexed autocorrelation has to be invented -- mad_ir's acf is
! indexed by LAG, and dividing out the nu^p prefactor is the whole conversion.
!
! R_target itself is the experiment on the same footing, scaled so that the
! strongest fitted band sits at R_MAX:
!
!   R_target_k^2 = norm I_exp(k)/nu_k^p,   norm = R_MAX^2/max_k(I_exp/nu^p)
!
!==========================================================================
! UNITS
!==========================================================================
!
! Everything internal is in rad/fs. The experimental grid, ir_damping and the
! output are wavenumbers in cm^-1, and the conversion is
!
!   w [rad/fs] = 2 pi nu [cm^-1] / CM_PER_INV_FS
!
! with CM_PER_INV_FS = 33356.40952 imported from mad_ir rather than redefined.
! Getting this wrong is not a small error: treating nu as rad/fs directly is
! off by 2 pi/3.34e4 ~ 1.9e-4 in every frequency, and multiplying by
! CM_PER_INV_FS instead of dividing puts a 10 cm^-1 bandwidth at 3.3e5 fs^-1
! rather than 1.9e-3. Check 6 exists because both were once in this file.
!
!==========================================================================
! INTEGRATOR
!==========================================================================
!
! Strang split, so the bank is stable for any dt the MD can take. The linear
! part of (1) -- harmonic, eta damping, drive held constant over the step -- is
! propagated in CLOSED FORM, as in mad_ir_xl:
!
!   xp = g M/(mu w^2),  u0 = X - xp,  v0 = P/mu,  lam = eta/2,
!   Om = sqrt(w^2 - lam^2)
!   X' = xp + e^{-lam h} (u0 cos(Om h) + (v0 + lam u0) sin(Om h)/Om)
!   v' =      e^{-lam h} (v0 cos(Om h) - (w^2 u0 + lam v0) sin(Om h)/Om)
!
! and the controller (7) is a half-step multiplicative rescaling on either
! side of it, which is exact because (8) is exactly multiplicative in R.
!
! eta is clamped to |eta| <= eta_max_frac * w_k. Without it the integral in (6)
! winds up without bound when the target is unreachable, and at |eta| > 2 w_k
! the channel goes overdamped and stops resonating at all -- it would then
! report an intensity for a frequency it is no longer sensitive to.

module ir_auxiliary_dynamics

   use kinds, only: dp
   use mad_ir, only: mad_ir_type, mad_ir_spectrum, CM_PER_INV_FS

   implicit none

   private
   public :: ir_aux_type, ir_aux_state
   public :: ir_aux_init, ir_aux_free, ir_aux_setup
   public :: ir_aux_advance, ir_aux_evaluate, ir_aux_forces
   public :: ir_aux_calibrate, ir_aux_calibrated, ir_aux_ready
   public :: ir_aux_active
   public :: ir_aux_save, ir_aux_load, ir_aux_write_spectrum
   public :: ir_aux_bank_energy, ir_aux_energy_pumped
   public :: IR_AUX_RESTART_VERSION

!  amu -> eV fs^2 / A^2, so that mu w^2 X is an eV/A force
   real(dp), parameter :: AMU = 103.6426965268_dp
   real(dp), parameter :: DEFAULT_EFF_MASS_AMU = 100.0_dp
   real(dp), parameter :: DEFAULT_DAMPING_CM = 10.0_dp
   real(dp), parameter :: DEFAULT_TAU_FS = 500.0_dp
   real(dp), parameter :: DEFAULT_ETA_MAX_FRAC = 0.1_dp
!  The amplitude the strongest fitted band is scaled to. Only the ratio
!  R_k/R_target matters to the controller, so this is a units choice -- but NOT
!  a free one: the back-reaction (4) goes as g_k X_k with g_k proportional to
!  R_target, so the force on the atoms scales as R_MAX SQUARED. 1 is the value
!  mad_ir_notes.md section 4 intends ("so that R_k ~ 1 when the MD trajectory
!  matches the experimental intensity"); an earlier draft used 3, which is a
!  factor of 9 on every force and helped blow a water box to 1e8 K.
   real(dp), parameter :: R_MAX = 1.0_dp
!  Channels below this are dropped: w = 0 has no resonator, and nu^p in the
!  calibration divides by it. ir_nu_min defaults to 0, so this is reachable.
   real(dp), parameter :: NU_FLOOR_CM = 1.0_dp

   integer, parameter :: IR_AUX_RESTART_VERSION = 1

   type :: ir_aux_type
      logical :: active = .false.
      logical :: calibrated = .false.
      integer :: n_modes = 0
      integer :: n_sites = 0
      integer :: n_steps = 0
      integer :: n_dropped = 0
!     Modes that carry a target and a coupling. The back-reaction is
!     normalised by this; see ir_aux_forces.
      integer :: n_live = 0
      real(dp) :: dt = 0.0_dp
      real(dp) :: eta_max_frac = DEFAULT_ETA_MAX_FRAC
      real(dp) :: nu_power = 2.0_dp
      real(dp) :: normalisation = 0.0_dp
      real(dp) :: dissim = 0.0_dp
      real(dp) :: dissim_ref = 0.0_dp
      real(dp) :: e_pump = 0.0_dp
      integer, allocatable :: kmap(:)
      real(dp), allocatable :: nu(:)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: eff_mass(:)
      real(dp), allocatable :: gamma_k(:)
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: kappa(:)
      real(dp), allocatable :: g_k(:)
      real(dp), allocatable :: X(:, :)
      real(dp), allocatable :: P(:, :)
      real(dp), allocatable :: eta(:)
      real(dp), allocatable :: R(:)
      real(dp), allocatable :: R_target(:)
      real(dp), allocatable :: I_exp(:)
      real(dp), allocatable :: wgt(:)
      real(dp), allocatable :: I_calc(:)
!     A channel the experiment or the trajectory says nothing about: no
!     target to hold and no coupling to hold it with. Carried explicitly
!     rather than inferred from g_k == 0, because "uncoupled" and "nothing
!     to do" are different statements and conflating them once already made
!     three checks in tests/ir_aux pass by propagating nothing at all.
      logical, allocatable :: muted(:)
   end type ir_aux_type

   type(ir_aux_type), save :: ir_aux_state
   logical, save :: ir_aux_active = .false.

contains

!**************************************************************************
!  True once the bank has been calibrated against a full ensemble. Before
!  that it applies no force: the couplings are not known yet.
   logical function ir_aux_calibrated(this)
      implicit none
      type(ir_aux_type), intent(in) :: this
      ir_aux_calibrated = this%active .and. this%calibrated
   end function ir_aux_calibrated

!**************************************************************************
   logical function ir_aux_ready(this)
      implicit none
      type(ir_aux_type), intent(in) :: this
      ir_aux_ready = ir_aux_calibrated(this) .and. (this%n_steps >= 1)
   end function ir_aux_ready

!**************************************************************************
!  The bank's own mechanical energy, sum_k [ P^2/2mu + 1/2 mu w^2 X^2 ], which
!  for this parametrisation is sum_k 1/2 mu_k w_k^2 R_k^2 exactly.
   real(dp) function ir_aux_bank_energy(this)
      implicit none
      type(ir_aux_type), intent(in) :: this
      integer :: m
      ir_aux_bank_energy = 0.0_dp
      if (.not. this%active) return
      do m = 1, this%n_modes
         ir_aux_bank_energy = ir_aux_bank_energy &
                              + 0.5_dp*this%eff_mass(m)*this%omega(m)**2*this%R(m)**2
      end do
   end function ir_aux_bank_energy

!**************************************************************************
!  Running total of the energy the controller has put into the bank (or taken
!  out of it). H_ext is not conserved by construction, so this is the ledger
!  entry that says by how much.
   real(dp) function ir_aux_energy_pumped(this)
      implicit none
      type(ir_aux_type), intent(in) :: this
      ir_aux_energy_pumped = this%e_pump
   end function ir_aux_energy_pumped

!**************************************************************************
   subroutine ir_aux_free(this)
      implicit none
      type(ir_aux_type), intent(inout) :: this

      if (allocated(this%kmap)) deallocate (this%kmap)
      if (allocated(this%nu)) deallocate (this%nu)
      if (allocated(this%omega)) deallocate (this%omega)
      if (allocated(this%eff_mass)) deallocate (this%eff_mass)
      if (allocated(this%gamma_k)) deallocate (this%gamma_k)
      if (allocated(this%tau)) deallocate (this%tau)
      if (allocated(this%kappa)) deallocate (this%kappa)
      if (allocated(this%g_k)) deallocate (this%g_k)
      if (allocated(this%X)) deallocate (this%X)
      if (allocated(this%P)) deallocate (this%P)
      if (allocated(this%eta)) deallocate (this%eta)
      if (allocated(this%R)) deallocate (this%R)
      if (allocated(this%R_target)) deallocate (this%R_target)
      if (allocated(this%I_exp)) deallocate (this%I_exp)
      if (allocated(this%wgt)) deallocate (this%wgt)
      if (allocated(this%I_calc)) deallocate (this%I_calc)
      if (allocated(this%muted)) deallocate (this%muted)

      this%active = .false.
      this%calibrated = .false.
      this%n_modes = 0
      this%n_steps = 0
      this%n_dropped = 0
      this%e_pump = 0.0_dp

   end subroutine ir_aux_free

!**************************************************************************
!  Allocate the bank on the frequencies of parent that can carry a resonator.
!  No experiment is read here and no coupling is set: that is ir_aux_calibrate,
!  which needs a full ensemble and therefore cannot run at setup time.
   subroutine ir_aux_init(this, parent, n_sites, dt, eff_mass_amu, damping_cm, tau_fs, &
                          kappa_in, eta_max_frac, ok, msg)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      type(mad_ir_type), intent(in) :: parent
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: eff_mass_amu
      real(dp), intent(in) :: damping_cm
      real(dp), intent(in) :: tau_fs
      real(dp), intent(in) :: kappa_in
      real(dp), intent(in) :: eta_max_frac
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg

      real(dp) :: two_pi
      real(dp) :: mass_use
      real(dp) :: damp_use
      real(dp) :: tau_use
      integer :: k
      integer :: m
      integer :: n_keep

      call ir_aux_free(this)

      ok = .false.
      msg = ""
      two_pi = 2.0_dp*dacos(-1.0_dp)

      if (n_sites < 1) then
         msg = "ir_aux: n_sites is less than one"
         return
      end if
      if (parent%n_freq < 1 .or. .not. allocated(parent%nu)) then
         msg = "ir_aux: the parent mad_ir state has no frequency grid"
         return
      end if
      if (dt <= 0.0_dp) then
         msg = "ir_aux: the stored-frame interval is not positive"
         return
      end if

      mass_use = eff_mass_amu
      if (mass_use <= 0.0_dp) mass_use = DEFAULT_EFF_MASS_AMU
      damp_use = damping_cm
      if (damp_use <= 0.0_dp) damp_use = DEFAULT_DAMPING_CM
      tau_use = tau_fs
      if (tau_use <= 0.0_dp) tau_use = DEFAULT_TAU_FS

!     Which channels can carry a resonator at all. w = 0 cannot, and the
!     calibration divides by nu^p, so the floor is not cosmetic.
      n_keep = 0
      do k = 1, parent%n_freq
         if (parent%nu(k) > NU_FLOOR_CM) n_keep = n_keep + 1
      end do
      if (n_keep < 1) then
         msg = "ir_aux: no fitted frequency above the floor; check ir_nu_min and ir_nu_max"
         return
      end if

      this%n_modes = n_keep
      this%n_dropped = parent%n_freq - n_keep
      this%n_sites = n_sites
      this%dt = dt
      this%nu_power = parent%nu_power
      this%eta_max_frac = eta_max_frac
      if (this%eta_max_frac <= 0.0_dp) this%eta_max_frac = DEFAULT_ETA_MAX_FRAC

      allocate (this%kmap(1:this%n_modes))
      allocate (this%nu(1:this%n_modes))
      allocate (this%omega(1:this%n_modes))
      allocate (this%eff_mass(1:this%n_modes))
      allocate (this%gamma_k(1:this%n_modes))
      allocate (this%tau(1:this%n_modes))
      allocate (this%kappa(1:this%n_modes))
      allocate (this%g_k(1:this%n_modes))
      allocate (this%X(1:3, 1:this%n_modes))
      allocate (this%P(1:3, 1:this%n_modes))
      allocate (this%eta(1:this%n_modes))
      allocate (this%R(1:this%n_modes))
      allocate (this%R_target(1:this%n_modes))
      allocate (this%I_exp(1:this%n_modes))
      allocate (this%wgt(1:this%n_modes))
      allocate (this%I_calc(1:this%n_modes))
      allocate (this%muted(1:this%n_modes))

      m = 0
      do k = 1, parent%n_freq
         if (parent%nu(k) <= NU_FLOOR_CM) cycle
         m = m + 1
         this%kmap(m) = k
         this%nu(m) = parent%nu(k)
!        cm^-1 -> rad/fs. The single most error-prone line in the module.
         this%omega(m) = two_pi*parent%nu(k)/CM_PER_INV_FS
         this%eff_mass(m) = mass_use*AMU
         this%gamma_k(m) = two_pi*damp_use/CM_PER_INV_FS
         this%tau(m) = tau_use
         if (kappa_in > 0.0_dp) then
            this%kappa(m) = kappa_in
         else
!           Critical damping of (9).
            this%kappa(m) = 2.0_dp/tau_use
         end if
         this%I_exp(m) = parent%I_exp(k)
         this%wgt(m) = parent%wgt(k)
      end do

      this%g_k = 0.0_dp
      this%X = 0.0_dp
      this%P = 0.0_dp
      this%eta = 0.0_dp
      this%R = 0.0_dp
      this%R_target = 0.0_dp
      this%I_calc = 0.0_dp
      this%muted = .false.
      this%e_pump = 0.0_dp
      this%n_steps = 0
      this%n_live = 0
      this%calibrated = .false.
      this%active = .true.

      ok = .true.
      write (msg, '(A,I0,A,I0,A)') "ir_aux: bank of ", this%n_modes, " modes (", &
         this%n_dropped, " channels below the frequency floor dropped)"

   end subroutine ir_aux_init

!**************************************************************************
!  Set the couplings and the targets from a full ensemble, and seed the bank at
!  its target amplitude. Call once, when mad_ir_ready(parent) first holds.
!
!  The seeding is deterministic -- no RNG -- for two reasons: the bank is
!  replicated on every MPI rank and advanced identically, so a per-rank draw
!  would make the ranks disagree about the forces; and by construction of (10)
!  the target amplitude IS the expected steady-state amplitude, so starting
!  there means there is no transient to wait out.
   subroutine ir_aux_calibrate(this, parent, ok, msg)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      type(mad_ir_type), intent(inout) :: parent
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg

      real(dp) :: two_pi
      real(dp) :: s_mm
      real(dp) :: i_scaled
      real(dp) :: max_scaled
      real(dp) :: rk_sq
      real(dp) :: phase
      real(dp) :: golden
      real(dp) :: e_eff_scale
      integer :: m
      integer :: k
      integer :: n_mute
      integer :: a

      ok = .false.
      msg = ""
      two_pi = 2.0_dp*dacos(-1.0_dp)

      if (.not. this%active) then
         msg = "ir_aux: calibrate called on an inactive bank"
         return
      end if

!     Make sure I_calc belongs to the ensemble as it stands now: the transform
!     is what S_MM is read from, and under ir_bias_mode = aux nothing else has
!     asked mad_ir for it this step.
      call mad_ir_spectrum(parent)

!     The experiment, on the same nu^p footing as the prediction, scaled so the
!     strongest fitted band sits at R_MAX.
      max_scaled = 0.0_dp
      do m = 1, this%n_modes
         i_scaled = this%I_exp(m)/this%nu(m)**this%nu_power
         if (i_scaled > max_scaled) max_scaled = i_scaled
      end do
      if (max_scaled <= 0.0_dp) then
         msg = "ir_aux: the experimental spectrum has no positive intensity in the fitted range"
         return
      end if
      this%normalisation = R_MAX**2/max_scaled

      golden = 0.5_dp*(dsqrt(5.0_dp) - 1.0_dp)
      n_mute = 0

      do m = 1, this%n_modes

         k = this%kmap(m)

         i_scaled = this%I_exp(m)/this%nu(m)**this%nu_power
!        Experimental noise can put a fitted point below zero. Such a channel
!        has no target to hold, so it is muted rather than driven to sqrt(-x).
         if (i_scaled <= 0.0_dp) then
            this%R_target(m) = 0.0_dp
            this%g_k(m) = 0.0_dp
            this%X(1:3, m) = 0.0_dp
            this%P(1:3, m) = 0.0_dp
            this%muted(m) = .true.
            n_mute = n_mute + 1
            cycle
         end if

         rk_sq = this%normalisation*i_scaled
         this%R_target(m) = dsqrt(rk_sq)

!        S_MM at this channel. I_calc = nu^p dt FT[win*acf], so dividing out
!        nu^p is the whole conversion; mad_ir's acf is indexed by lag and never
!        by frequency, and no second transform is needed.
         s_mm = parent%I_calc(k)/this%nu(m)**this%nu_power

!        A channel the trajectory puts no power into cannot be calibrated: (10)
!        would divide by zero, and a NaN in one g_k poisons every force through
!        the shared E_eff. Mute it and say how many.
         if (s_mm <= 0.0_dp) then
!           Without a coupling there is nothing to drive the channel up to its
!           target, so the target goes too: leaving it set would have the
!           controller chase an amplitude it can never reach and wind eta into
!           its clamp for the whole run.
            this%R_target(m) = 0.0_dp
            this%g_k(m) = 0.0_dp
            this%X(1:3, m) = 0.0_dp
            this%P(1:3, m) = 0.0_dp
            this%muted(m) = .true.
            n_mute = n_mute + 1
            cycle
         end if

!        (10): g_k = mu_k w_k sqrt(Gamma_k <R_target^2> / S_MM(w_k))
         this%g_k(m) = this%eff_mass(m)*this%omega(m) &
                       *dsqrt(this%gamma_k(m)*rk_sq/s_mm)
         this%muted(m) = .false.

!        Seed at the target amplitude, with a deterministic low-discrepancy
!        phase per mode and the three components in quadrature thirds, so the
!        bank neither starts in lockstep nor needs a random number.
         phase = two_pi*(dfloat(m)*golden - dint(dfloat(m)*golden))
         do a = 1, 3
            this%X(a, m) = this%R_target(m)/dsqrt(3.0_dp) &
                           *dcos(phase + two_pi*dfloat(a - 1)/3.0_dp)
            this%P(a, m) = this%eff_mass(m)*this%omega(m)*this%R_target(m)/dsqrt(3.0_dp) &
                           *dsin(phase + two_pi*dfloat(a - 1)/3.0_dp)
         end do

      end do

      this%eta = 0.0_dp
      this%e_pump = 0.0_dp
      this%n_steps = 0
      this%n_live = this%n_modes - n_mute
      call ir_aux_amplitude(this)
      this%calibrated = .true.

      ok = .true.
!     The scale of the field the atoms are about to feel, reported at once.
!     |E_eff| times a dipole gradient of order one electron is roughly the
!     back-reaction in eV/A, and an atom in liquid water carries a few eV/A of
!     real force. A number far above that here means the run will not survive
!     the first biased step, and the whole point of printing it is to see that
!     now rather than nine hundred steps later.
      e_eff_scale = 0.0_dp
      do m = 1, this%n_modes
         if (this%muted(m)) cycle
         e_eff_scale = e_eff_scale + this%g_k(m)*this%R_target(m)
      end do
      e_eff_scale = e_eff_scale/dfloat(max(1, this%n_live))
      write (msg, '(A,I0,A,I0,A,ES11.3,A)') "ir_aux: calibrated ", this%n_live, &
         " of ", this%n_modes, " modes; |E_eff| ~ ", e_eff_scale, " eV/A per unit dipole gradient"

   end subroutine ir_aux_calibrate

!**************************************************************************
!  R_k from (2). Kept in one place: it is the observable, and three routines
!  need it.
   subroutine ir_aux_amplitude(this)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      real(dp) :: inv
      integer :: m

      do m = 1, this%n_modes
         inv = 1.0_dp/(this%eff_mass(m)*this%omega(m))**2
         this%R(m) = dsqrt(this%X(1, m)**2 + this%X(2, m)**2 + this%X(3, m)**2 &
                           + (this%P(1, m)**2 + this%P(2, m)**2 + this%P(3, m)**2)*inv)
      end do

   end subroutine ir_aux_amplitude

!**************************************************************************
!  One step of (1) with (6) and (7), Strang split: half a step of the radial
!  controller, an exact step of the linear resonator with the drive held
!  constant, half a step of the controller again, then the integral update.
   subroutine ir_aux_advance(this, dipole)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      real(dp), intent(in) :: dipole(1:3)

      real(dp) :: h
      real(dp) :: half
      real(dp) :: lam
      real(dp) :: om_d
      real(dp) :: w2
      real(dp) :: ec
      real(dp) :: cs
      real(dp) :: sn
      real(dp) :: xp
      real(dp) :: u0
      real(dp) :: v0
      real(dp) :: scal
      real(dp) :: e_before
      real(dp) :: e_after
      real(dp) :: work_drive
      real(dp) :: eta_cap
      real(dp) :: x_old(1:3)
      integer :: m
      integer :: a

      if (.not. ir_aux_calibrated(this)) return

      h = this%dt
      half = 0.5_dp*h

      call ir_aux_amplitude(this)

      do m = 1, this%n_modes

         if (this%muted(m)) cycle

         w2 = this%omega(m)**2
         e_before = 0.5_dp*this%eff_mass(m)*w2*this%R(m)**2
         x_old(1:3) = this%X(1:3, m)

!        ---- half step of the radial controller (7) -----------------------
!        Exactly multiplicative in R, so one scalar on both X and P rescales
!        the amplitude and leaves the phase and the frequency untouched.
         scal = ir_aux_radial_factor(this, m, half)
         this%X(1:3, m) = scal*this%X(1:3, m)
         this%P(1:3, m) = scal*this%P(1:3, m)

!        ---- exact step of the linear resonator --------------------------
         lam = 0.5_dp*this%eta(m)
         om_d = w2 - lam*lam
!        The eta clamp keeps this positive; the guard is for a restart file
!        that disagrees.
         if (om_d <= 0.0_dp) then
            lam = 0.0_dp
            om_d = w2
         end if
         om_d = dsqrt(om_d)
         ec = dexp(-lam*h)
         cs = dcos(om_d*h)
         sn = dsin(om_d*h)/om_d

         do a = 1, 3
            xp = this%g_k(m)*dipole(a)/(this%eff_mass(m)*w2)
            u0 = this%X(a, m) - xp
            v0 = this%P(a, m)/this%eff_mass(m)
            this%X(a, m) = xp + ec*(u0*cs + (v0 + lam*u0)*sn)
            this%P(a, m) = this%eff_mass(m)*ec*(v0*cs - (w2*u0 + lam*v0)*sn)
         end do

!        ---- second half step of the controller --------------------------
         call ir_aux_amplitude_one(this, m)
         scal = ir_aux_radial_factor(this, m, half)
         this%X(1:3, m) = scal*this%X(1:3, m)
         this%P(1:3, m) = scal*this%P(1:3, m)

         call ir_aux_amplitude_one(this, m)

!        ---- the integral update (6) -------------------------------------
         if (this%R_target(m) > 0.0_dp) then
            this%eta(m) = this%eta(m) &
                          + h/this%tau(m)**2*((this%R(m)/this%R_target(m))**2 - 1.0_dp)
!           Anti-windup. Beyond |eta| = 2 w the channel stops resonating and
!           would report an intensity for a frequency it no longer senses.
            eta_cap = this%eta_max_frac*this%omega(m)
            if (this%eta(m) > eta_cap) this%eta(m) = eta_cap
            if (this%eta(m) < -eta_cap) this%eta(m) = -eta_cap
         end if

!        ---- the ledger --------------------------------------------------
!        Everything the step changed, less the work the dipole did through the
!        coupling, is what eta and the controller injected.
         e_after = 0.5_dp*this%eff_mass(m)*w2*this%R(m)**2
         work_drive = 0.0_dp
         do a = 1, 3
            work_drive = work_drive + this%g_k(m)*dipole(a)*(this%X(a, m) - x_old(a))
         end do
         this%e_pump = this%e_pump + (e_after - e_before) - work_drive

      end do

      this%n_steps = this%n_steps + 1

   end subroutine ir_aux_advance

!**************************************************************************
   subroutine ir_aux_amplitude_one(this, m)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      integer, intent(in) :: m
      real(dp) :: inv

      inv = 1.0_dp/(this%eff_mass(m)*this%omega(m))**2
      this%R(m) = dsqrt(this%X(1, m)**2 + this%X(2, m)**2 + this%X(3, m)**2 &
                        + (this%P(1, m)**2 + this%P(2, m)**2 + this%P(3, m)**2)*inv)

   end subroutine ir_aux_amplitude_one

!**************************************************************************
!  exp(c_k dt) for the radial flow (7)-(8). Returned as a factor rather than a
!  rate because that is what makes the half step exact.
   real(dp) function ir_aux_radial_factor(this, m, dt_use)

      implicit none

      type(ir_aux_type), intent(in) :: this
      integer, intent(in) :: m
      real(dp), intent(in) :: dt_use
      real(dp) :: c

      ir_aux_radial_factor = 1.0_dp
      if (this%R_target(m) <= 0.0_dp) return
      if (this%R(m) <= 0.0_dp) return

      c = -this%kappa(m)*(this%R(m) - this%R_target(m))/this%R(m)
      ir_aux_radial_factor = dexp(c*dt_use)

   end function ir_aux_radial_factor

!**************************************************************************
!  The coupling energy -escale sum_k g_k X_k.M, which is the generator of the
!  force ir_aux_forces adds, and the dissimilarity for thermo.log.
   subroutine ir_aux_evaluate(this, energy_scale, dipole, energy)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      real(dp), intent(in) :: energy_scale
      real(dp), intent(in) :: dipole(1:3)
      real(dp), intent(out) :: energy

      real(dp) :: acc
      real(dp) :: num
      real(dp) :: den
      real(dp) :: d
      integer :: m

      energy = 0.0_dp
      this%dissim = 0.0_dp
      this%dissim_ref = 0.0_dp
      if (.not. ir_aux_calibrated(this)) return

      call ir_aux_amplitude(this)

      acc = 0.0_dp
      num = 0.0_dp
      den = 0.0_dp
      do m = 1, this%n_modes
         if (this%muted(m)) cycle
         acc = acc + this%g_k(m)*(this%X(1, m)*dipole(1) &
                                  + this%X(2, m)*dipole(2) + this%X(3, m)*dipole(3))
!        Back to intensity units, so ir_aux_spectrum.dat is on the same axes as
!        ir_exp.dat: R^2 = norm I/nu^p inverts to I = nu^p R^2/norm.
         if (this%normalisation > 0.0_dp) then
            this%I_calc(m) = this%nu(m)**this%nu_power*this%R(m)**2/this%normalisation
         else
            this%I_calc(m) = 0.0_dp
         end if
         if (this%R_target(m) > 0.0_dp) then
            d = this%R(m)**2 - this%R_target(m)**2
            num = num + this%wgt(m)*d*d
            den = den + this%wgt(m)*this%R_target(m)**4
         end if
      end do

!     The same 1/n_live as the force, so the reported energy remains exactly
!     the generator of the reported force.
      energy = -energy_scale*acc/dfloat(max(1, this%n_live))

      if (den > 0.0_dp) then
         this%dissim = dsqrt(num/den)
         this%dissim_ref = 1.0_dp
      end if

   end subroutine ir_aux_evaluate

!**************************************************************************
!  F_jb += (escale/n_live) sum_a (sum_k g_k X_ka) dmu_a/dr_jb.
!
!  dmu_dr(a, b, j) = d mu_a / d r_jb, exactly as accumulate_dmu_dr leaves it.
!  The sign is + : the coupling energy is -g X.M, so the force is +g J^T X.
!  mad_ir_forces subtracts because its lambda is dL/dmu; here the vector
!  sum_k g_k X_k plays the part of -lambda.
!
!  WHY THE 1/n_live. Every mode is driven by the SAME total dipole, so the X_k
!  are strongly correlated in direction and sum_k g_k X_k grows like n_live,
!  not like sqrt(n_live). The bank is an OVER-COMPLETE measurement of one
!  signal: each mode measures M in its own band, and for that each needs its
!  full calibrated g_k, but the total back-action must not count the same
!  dipole once per channel. Two consequences, and the second is the reason this
!  is not merely cosmetic:
!
!    - Without it the force is ~n_live times too large. Measured: a 64-molecule
!      water box at 300 K with 101 fitted channels reached 1e8 K within ten
!      steps of the bias switching on, with a coupling energy of -115 eV on a
!      -950 eV system in the FIRST biased step.
!    - Without it exp_energy_scales would silently mean something different for
!      every choice of ir_nu_min/ir_nu_max, because changing the fitted range
!      changes n_live. A bias whose strength depends on the plotting range is
!      not a bias anyone can calibrate.
!
!  The price is that the drive and the back-reaction no longer use the same
!  coupling, so the coupling term is not strictly Hamiltonian even before the
!  controller is added. Given that the controller has already given up
!  conservation this is not a further loss, but it is a real asymmetry and it is
!  stated rather than hidden. The principled alternative is to make the bank a
!  genuine partition of unity -- Gamma_k tied to the grid spacing, so that
!  sum_k chi_k(w) ~ 1 and the sum is a true reconstruction of M rather than
!  n_live copies of it -- which is the better long-term answer.
   subroutine ir_aux_forces(this, energy_scale, dmu_dr, forces)

      implicit none

      type(ir_aux_type), intent(in) :: this
      real(dp), intent(in) :: energy_scale
      real(dp), intent(in) :: dmu_dr(:, :, :)
      real(dp), intent(inout) :: forces(:, :)

      real(dp) :: e_eff(1:3)
      real(dp) :: acc
      integer :: m
      integer :: j
      integer :: a
      integer :: b
      integer :: n_atoms

      if (.not. ir_aux_calibrated(this)) return

      n_atoms = size(dmu_dr, 3)

      e_eff = 0.0_dp
      do m = 1, this%n_modes
         if (this%muted(m)) cycle
         e_eff(1) = e_eff(1) + this%g_k(m)*this%X(1, m)
         e_eff(2) = e_eff(2) + this%g_k(m)*this%X(2, m)
         e_eff(3) = e_eff(3) + this%g_k(m)*this%X(3, m)
      end do
      e_eff = energy_scale*e_eff/dfloat(max(1, this%n_live))

      !$omp parallel do private(j, b, a, acc) schedule(static)
      do j = 1, n_atoms
         do b = 1, 3
            acc = 0.0_dp
            do a = 1, 3
               acc = acc + e_eff(a)*dmu_dr(a, b, j)
            end do
            forces(b, j) = forces(b, j) + acc
         end do
      end do
      !$omp end parallel do

   end subroutine ir_aux_forces

!**************************************************************************
!  Lazy first-use setup, for the same reason mad_ir's and the GLE bath's are
!  lazy: n_sites is not known at input-reading time.
   subroutine ir_aux_setup(parent, n_sites, dt_md, stride, eff_mass_amu, damping_cm, &
                           tau_fs, kappa_in, eta_max_frac, restart_file, ok, resumed, msg)

      implicit none

      type(mad_ir_type), intent(in) :: parent
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: dt_md
      integer, intent(in) :: stride
      real(dp), intent(in) :: eff_mass_amu
      real(dp), intent(in) :: damping_cm
      real(dp), intent(in) :: tau_fs
      real(dp), intent(in) :: kappa_in
      real(dp), intent(in) :: eta_max_frac
      character(len=*), intent(in) :: restart_file
      logical, intent(out) :: ok
      logical, intent(out) :: resumed
      character(len=*), intent(out) :: msg

      real(dp) :: dt_use
      logical :: lok
      character(len=512) :: lmsg

      resumed = .false.
      dt_use = dt_md*dfloat(max(1, stride))

      call ir_aux_init(ir_aux_state, parent, n_sites, dt_use, eff_mass_amu, &
                       damping_cm, tau_fs, kappa_in, eta_max_frac, ok, msg)
      if (.not. ok) return

      if (len_trim(restart_file) > 0 .and. trim(restart_file) /= "none") then
         call ir_aux_load(ir_aux_state, restart_file, lok, lmsg)
         if (lok) then
            resumed = .true.
            msg = trim(msg)//"; "//trim(lmsg)
         end if
!        A load failure is not fatal: the run calibrates a fresh bank once the
!        ensemble fills, which is what it would have done without the file.
      end if

      ir_aux_active = .true.

   end subroutine ir_aux_setup

!**************************************************************************
   subroutine ir_aux_save(this, fname, ok, msg)

      implicit none

      type(ir_aux_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg

      integer :: iu
      integer :: ios
      integer :: m

      ok = .false.
      msg = ""
      if (.not. this%active) then
         msg = "ir_aux_save: nothing to write, the bank is not active"
         return
      end if

      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         msg = "ir_aux_save: cannot open "//trim(fname)
         return
      end if

!     The file says what it is. A loader that identifies a file by its length
!     will one day adopt another observable's state.
      write (iu, '(A,1X,I0)') "IR_AUX_RESTART", IR_AUX_RESTART_VERSION
      write (iu, '(A)') "#  n_modes  n_steps  calibrated  normalisation  e_pump"
      write (iu, '(2(I10,1X),L2,1X,2(ES24.16,1X))') &
         this%n_modes, this%n_steps, this%calibrated, this%normalisation, this%e_pump
      write (iu, '(A)') "#  nu_cm-1  X(1:3)  P(1:3)  eta  eff_mass  omega  g_k  R_target"
      do m = 1, this%n_modes
         write (iu, '(12(ES24.16,1X))') &
            this%nu(m), &
            this%X(1, m), this%X(2, m), this%X(3, m), &
            this%P(1, m), this%P(2, m), this%P(3, m), &
            this%eta(m), &
            this%eff_mass(m), &
            this%omega(m), &
            this%g_k(m), &
            this%R_target(m)
      end do
      close (iu)

      ok = .true.
      write (msg, '(A,A,A,I0,A)') "ir_aux_save: wrote ", trim(fname), " (", this%n_modes, " modes)"

   end subroutine ir_aux_save

!**************************************************************************
!  Refuse a file that does not describe this bank rather than adopt it. The
!  frequencies are the identity: a grid that has moved means the modes mean
!  something else, and silently carrying X and P across would be worse than
!  starting fresh.
   subroutine ir_aux_load(this, fname, ok, msg)

      implicit none

      type(ir_aux_type), intent(inout) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg

      character(len=64) :: magic
      character(len=512) :: line
      integer :: iu
      integer :: ios
      integer :: ver
      integer :: n_modes_f
      integer :: n_steps_f
      integer :: m
      logical :: cal_f
      logical :: there
      real(dp) :: norm_f
      real(dp) :: pump_f
      real(dp) :: nu_f
      real(dp) :: xf(1:3)
      real(dp) :: pf(1:3)
      real(dp) :: eta_f
      real(dp) :: mass_f
      real(dp) :: om_f
      real(dp) :: g_f
      real(dp) :: rt_f

      ok = .false.
      msg = ""

      inquire (file=trim(fname), exist=there)
      if (.not. there) then
         msg = "ir_aux_load: "//trim(fname)//" does not exist; starting a fresh bank"
         return
      end if

      open (newunit=iu, file=trim(fname), status="old", action="read", iostat=ios)
      if (ios /= 0) then
         msg = "ir_aux_load: cannot open "//trim(fname)
         return
      end if

      read (iu, *, iostat=ios) magic, ver
      if (ios /= 0 .or. trim(magic) /= "IR_AUX_RESTART") then
         close (iu)
         msg = "ir_aux_load: "//trim(fname)//" is not an ir_aux restart file"
         return
      end if
      if (ver /= IR_AUX_RESTART_VERSION) then
         close (iu)
         write (msg, '(A,I0,A,I0)') "ir_aux_load: restart version ", ver, &
            " but this build writes ", IR_AUX_RESTART_VERSION
         return
      end if

      read (iu, '(A)', iostat=ios) line
      read (iu, *, iostat=ios) n_modes_f, n_steps_f, cal_f, norm_f, pump_f
      if (ios /= 0) then
         close (iu)
         msg = "ir_aux_load: the header of "//trim(fname)//" is unreadable"
         return
      end if
      if (n_modes_f /= this%n_modes) then
         close (iu)
         write (msg, '(A,I0,A,I0,A)') "ir_aux_load: the file holds ", n_modes_f, &
            " modes but this run has ", this%n_modes, "; starting a fresh bank"
         return
      end if

      read (iu, '(A)', iostat=ios) line
      do m = 1, n_modes_f
         read (iu, *, iostat=ios) nu_f, xf(1), xf(2), xf(3), pf(1), pf(2), pf(3), &
            eta_f, mass_f, om_f, g_f, rt_f
         if (ios /= 0) then
            close (iu)
            write (msg, '(A,I0,A)') "ir_aux_load: truncated at mode ", m, "; starting a fresh bank"
            return
         end if
         if (dabs(nu_f - this%nu(m)) > 1.0e-6_dp*max(1.0_dp, dabs(this%nu(m)))) then
            close (iu)
            write (msg, '(A,I0,A)') "ir_aux_load: frequency grid disagrees at mode ", m, &
               "; starting a fresh bank"
            return
         end if
         this%X(1:3, m) = xf(1:3)
         this%P(1:3, m) = pf(1:3)
         this%eta(m) = eta_f
         this%eff_mass(m) = mass_f
         this%omega(m) = om_f
         this%g_k(m) = g_f
         this%R_target(m) = rt_f
      end do
      close (iu)

      this%n_steps = n_steps_f
      this%normalisation = norm_f
      this%e_pump = pump_f
      this%calibrated = cal_f
!     calibrate() leaves a muted channel with BOTH g_k and R_target at
!     zero and a live one with both positive, so the flag is recoverable
!     and does not need a column of its own.
      this%n_live = 0
      do m = 1, this%n_modes
         this%muted(m) = (this%g_k(m) == 0.0_dp) .or. (this%R_target(m) <= 0.0_dp)
         if (.not. this%muted(m)) this%n_live = this%n_live + 1
      end do
      call ir_aux_amplitude(this)

      ok = .true.
      write (msg, '(A,A)') "ir_aux_load: resumed from ", trim(fname)

   end subroutine ir_aux_load

!**************************************************************************
   subroutine ir_aux_write_spectrum(this, fname, with_exp)

      implicit none

      type(ir_aux_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(in) :: with_exp

      integer :: iu
      integer :: ios
      integer :: m

      if (.not. this%active) return

      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) return

      if (with_exp) then
         write (iu, '(A)') "#  ir_aux: wavenumber_cm-1   predicted   experimental   R   R_target   eta"
      else
         write (iu, '(A)') "#  ir_aux: wavenumber_cm-1   predicted   R   eta"
      end if
      do m = 1, this%n_modes
         if (with_exp) then
            write (iu, '(6(ES20.10,1X))') this%nu(m), this%I_calc(m), this%I_exp(m), &
               this%R(m), this%R_target(m), this%eta(m)
         else
            write (iu, '(4(ES20.10,1X))') this%nu(m), this%I_calc(m), this%R(m), this%eta(m)
         end if
      end do
      close (iu)

   end subroutine ir_aux_write_spectrum

end module ir_auxiliary_dynamics
