! Experimental IR forces for MAD.
!
! The observable is an IR spectrum obtained from the autocorrelation function of
! the total dipole along a trajectory, and the MAD force is the gradient of the
! mismatch between that spectrum and an experimental one.
!
! WHAT SETS THE SIZE OF THE ENSEMBLE
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
! WHERE THE FORCE ACTS
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
! constants here. One consequence worth being explicit about: this is not the
! gradient of any potential, so it does not conserve energy and should be
! thought of as a bias, not a force field.
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
   public :: mad_ir_select_range, mad_ir_save, mad_ir_load
   public :: mad_ir_setup, mad_ir_state, mad_ir_dmu_dr, mad_ir_collect
   public :: mad_ir_write_calc_spectrum
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
! window_kind: "hann" (default) or "none". A lag window is not cosmetic here --
! C(tau) at tau near n_lag is averaged over very few pairs, and transforming it
! unwindowed puts that noise straight into the spectrum and hence into the
! force.
!
   subroutine mad_ir_init(this, dt, n_lag, n_window, nu, I_exp, wgt, match_scale, nu_power, window_kind)

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
      real(dp) :: pi
      integer :: t

      call mad_ir_free(this)

      pi = dacos(-1.d0)
      this%dt = dt
      this%n_lag = n_lag
      this%n_window = n_window
      this%n_stored = 0
      this%head = 0
      this%n_freq = size(nu)
      this%match_scale = match_scale
      this%nu_power = nu_power
      this%scale = 1.d0

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

      if (trim(window_kind) == "none") then
         this%win = 1.d0
      else
         do t = 0, n_lag
            this%win(t) = 0.5d0*(1.d0 + dcos(pi*dfloat(t)/dfloat(n_lag)))
         end do
      end if

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
! Spectrum, loss, and the sensitivity of the loss to the newest dipole.
!
!   C(tau)  = 1/(N-tau) sum_a mu(a) . mu(a+tau)          a = age, 0 is newest
!   I(nu_k) = nu_k^p dt [ win(0)C(0) + 2 sum_tau win(tau) C(tau) cos(th_k tau) ]
!   L       = sum_k wgt_k ( s I_k - I_exp_k )^2
!
! with s fitted by least squares when match_scale is on -- the computed
! spectrum is in arbitrary units, so without that the loss compares two things
! on different scales and its gradient is meaningless. s sits at its own
! optimum, so by the envelope theorem dL/dI_k picks up no extra term from it.
!
! Then, since mu of the newest frame (age 0) appears in C(0) twice and in
! C(tau>0) once,
!
!   dC(0)/dmu_c   = 2 mu_c(0)/N,   dC(tau)/dmu_c = mu_c(tau)/(N-tau)
!
! so lambda falls out as a single pass over the stored dipoles.
!
   subroutine mad_ir_evaluate(this, energy_scale, energy, lambda)

      implicit none

      type(mad_ir_type), intent(inout) :: this
      real(dp), intent(in) :: energy_scale
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: lambda(1:3)
      real(dp), allocatable :: dLdI(:)
      real(dp), allocatable :: S(:)
      real(dp) :: two_pi, theta, acc, num, den, s_fit, pref, c
      integer :: t, a, k, m, m2, N, c1

      energy = 0.d0
      lambda = 0.d0
      if (.not. mad_ir_ready(this)) return

      N = this%n_window
      two_pi = 2.d0*dacos(-1.d0)

!     ---- autocorrelation of the stored dipoles ----------------------
      do t = 0, this%n_lag
         acc = 0.d0
         do a = 0, N - 1 - t
            m = this%head - a
            if (m < 1) m = m + N
            m2 = this%head - a - t
            if (m2 < 1) m2 = m2 + N
            acc = acc + this%mu_hist(1, m)*this%mu_hist(1, m2) &
                  + this%mu_hist(2, m)*this%mu_hist(2, m2) &
                  + this%mu_hist(3, m)*this%mu_hist(3, m2)
         end do
         this%acf(t) = acc/dfloat(N - t)
      end do

!     ---- spectrum ---------------------------------------------------
      do k = 1, this%n_freq
         acc = this%win(0)*this%acf(0)
         do t = 1, this%n_lag
            theta = two_pi*this%nu(k)*dfloat(t)*this%dt/CM_PER_INV_FS
            acc = acc + 2.d0*this%win(t)*this%acf(t)*dcos(theta)
         end do
         this%I_calc(k) = (this%nu(k)**this%nu_power)*this%dt*acc
      end do

!     ---- overall scale ----------------------------------------------
      if (this%match_scale) then
         num = 0.d0
         den = 0.d0
         do k = 1, this%n_freq
            num = num + this%wgt(k)*this%I_calc(k)*this%I_exp(k)
            den = den + this%wgt(k)*this%I_calc(k)**2
         end do
         if (den > 0.d0) then
            s_fit = num/den
         else
            s_fit = 1.d0
         end if
      else
         s_fit = 1.d0
      end if
      this%scale = s_fit

!     ---- energy and its sensitivity to the spectrum ------------------
!     The same form every other MAD observable uses (exp_utils'
!     get_exp_energies): E = 1/2 * energy_scale * sum (y_pred - y_exp)^2, so
!     that exp_energy_scales means the same thing here as it does for a pair
!     distribution or a structure factor, and the force really is the gradient
!     of the spectral difference rather than of something proportional to it.
      allocate (dLdI(1:this%n_freq))
      do k = 1, this%n_freq
         c = s_fit*this%I_calc(k) - this%I_exp(k)
         energy = energy + 0.5d0*energy_scale*this%wgt(k)*c**2
         dLdI(k) = energy_scale*this%wgt(k)*s_fit*c
      end do

!     ---- back through the transform, to dL/dC(tau) -------------------
      allocate (S(0:this%n_lag))
      S = 0.d0
      do k = 1, this%n_freq
         pref = dLdI(k)*(this%nu(k)**this%nu_power)*this%dt
         S(0) = S(0) + pref*this%win(0)
         do t = 1, this%n_lag
            theta = two_pi*this%nu(k)*dfloat(t)*this%dt/CM_PER_INV_FS
            S(t) = S(t) + pref*2.d0*this%win(t)*dcos(theta)
         end do
      end do

!     ---- and through the autocorrelation, to dL/dmu(newest) ----------
      do c1 = 1, 3
         acc = 2.d0*S(0)*this%mu_hist(c1, this%head)/dfloat(N)
         do t = 1, this%n_lag
            m = this%head - t
            if (m < 1) m = m + N
            acc = acc + S(t)*this%mu_hist(c1, m)/dfloat(N - t)
         end do
         lambda(c1) = acc
      end do

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
                           window_kind, n_atoms, ok, resumed, msg)

      implicit none

      real(dp), intent(in) :: dt_md
      integer, intent(in) :: stride
      real(dp), intent(in) :: nu_res, nu_min, nu_max
      integer, intent(in) :: lag_factor
      real(dp), intent(in) :: nu_in(:), I_in(:)
      character(len=*), intent(in) :: restart_file, window_kind
      logical, intent(in) :: match_scale
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

      call mad_ir_select_range(nu_in, I_in, nu_min, nu_max, nu, I_exp, wgt, ok2, msg)
      if (.not. ok2) return

      call mad_ir_init(mad_ir_state, dt, n_lag, n_window, nu, I_exp, wgt, &
                       match_scale, nu_power, window_kind)

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
! The computed spectrum next to the experimental one, for watching a fit.
! Columns: wavenumber, computed (scaled as the loss sees it), experimental.
!
   subroutine mad_ir_write_calc_spectrum(this, fname)
      implicit none
      type(mad_ir_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      integer :: iu, ios, k
      if (.not. this%active) return
      open (newunit=iu, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) return
      write (iu, '(A)') "# wavenumber_cm-1   computed_scaled   experimental"
      do k = 1, this%n_freq
         write (iu, '(3ES20.10)') this%nu(k), this%scale*this%I_calc(k), this%I_exp(k)
      end do
      close (iu)
   end subroutine mad_ir_write_calc_spectrum

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
   subroutine mad_ir_select_range(nu_in, I_in, nu_min, nu_max, nu, I_exp, wgt, ok, msg)

      implicit none

      real(dp), intent(in) :: nu_in(:)
      real(dp), intent(in) :: I_in(:)
      real(dp), intent(in) :: nu_min
      real(dp), intent(in) :: nu_max
      real(dp), allocatable, intent(out) :: nu(:)
      real(dp), allocatable, intent(out) :: I_exp(:)
      real(dp), allocatable, intent(out) :: wgt(:)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
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
      wgt = 1.d0
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
