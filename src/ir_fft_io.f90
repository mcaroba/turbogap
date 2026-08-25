! Everything ir_fft.f90 deliberately does not do: hold a trajectory, work out
! its timestep, write files, and adapt mad_ir's ring buffer to the array
! ir_fft_loss wants.
!
! ir_fft.f90 is pure numerics so that it can be driven from a bare test
! program and compared with the Python line by line. This module is where it
! meets TurboGAP.
!
!==========================================================================
! THE TWO WAYS IN
!==========================================================================
!
! PREDICTION (turbogap predict, do_ir = .true.). The frames are already on
! disk: an extended-xyz trajectory, one configuration per frame, each carrying
! its time in a time= tag on the comment line -- which is exactly what
! TurboGAP's own write_extxyz emits, so a trajectory_out.xyz from a previous
! run is a valid input with no preparation at all. The dipole model is
! evaluated on each frame as it is read, ir_fft_frames_push accumulates
! (time, mu), and at the end of the file ir_fft_frames_finish transforms the
! lot.
!
! THE TIMESTEP IS READ, NOT ASSUMED, and that is the point of the time= tag.
! A trajectory written every write_xyz steps of an md_step run has a frame
! interval of write_xyz*md_step, and getting that wrong scales the whole
! wavenumber axis by the ratio -- a spectrum that is wrong by a factor of two
! in frequency and looks completely plausible. ir_fft_frames_dt takes the
! interval from the labels and REFUSES a trajectory whose spacing is not
! uniform, because an unevenly sampled series has no Fourier transform of the
! kind being computed and the failure is otherwise silent.
!
! MD (turbogap md, ir_bias_mode = fft). The ensemble is mad_ir's rolling
! buffer, already filled by mad_ir_push for every estimator; ir_fft_md_unroll
! copies it out in chronological order and ir_fft_loss returns the gradient
! with respect to the newest dipole. That gradient is dL/dmu of the current
! configuration, which is precisely what mad_ir_forces contracts with dmu/dr,
! so the force path downstream is the ACF bias's, unchanged.
!
!==========================================================================
! WHAT IS WRITTEN
!==========================================================================
!
!   ir_fft_spectrum.dat   nu, intensity, power, and the raw (un-normalised)
!                         columns of each, with the sizing in the header. When
!                         there is an experiment it also carries the fitted
!                         scale and the experimental curve on the same grid.
!   ir_fft_dipoles.dat    time, mux, muy, muz -- the trajectory the spectrum
!                         was computed from.
!
! ir_fft_dipoles.dat exists so that the run can be checked against the Python
! it was translated from without rerunning anything: load those three columns
! and call compute_ir_spectrum on them. That is a stronger test than any
! internal one, and it costs 4 columns of text.
!
! ir_spectrum.dat is NOT touched. It stays what it has always been -- the
! block ACF estimator of mad_ir.f90 -- so that a run with ir_bias_mode = fft
! carries two independent estimates of the same spectrum, exactly as an xl run
! does.
!
module ir_fft_io

   use kinds
   use ir_fft

   implicit none

   private
   public :: ir_fft_frames_type, ir_fft_frames
   public :: ir_fft_frames_reset, ir_fft_frames_push, ir_fft_frames_dt
   public :: ir_fft_frames_finish, ir_fft_write_dipoles, ir_fft_write_spectrum
   public :: ir_fft_md_unroll
   public :: ir_fft_config_from_params

!  A growing (time, dipole) buffer. Doubling, so pushing n frames is O(n)
!  copies in total and the trajectory length does not have to be known in
!  advance -- which it is not, because the frame count of an xyz file is
!  discovered by reaching the end of it.
   type :: ir_fft_frames_type
      integer :: n = 0
      integer :: cap = 0
!     .false. as soon as ONE frame arrives without a time= tag. A partly
!     labelled trajectory is not usable: the interval cannot be checked, and
!     silently falling back to a nominal dt for the unlabelled ones would
!     produce a spectrum from a time axis that is partly invented.
      logical :: have_times = .true.
      integer :: n_missing = 0
      real(dp), allocatable :: mu(:, :)   ! (1:3, 1:cap)
      real(dp), allocatable :: t(:)       ! (1:cap), fs
   end type ir_fft_frames_type

   type(ir_fft_frames_type), save :: ir_fft_frames

contains

   subroutine ir_fft_frames_reset(this)

      implicit none

      type(ir_fft_frames_type), intent(inout) :: this

      if (allocated(this%mu)) deallocate (this%mu)
      if (allocated(this%t)) deallocate (this%t)
      this%n = 0
      this%cap = 0
      this%have_times = .true.
      this%n_missing = 0

   end subroutine ir_fft_frames_reset

!**************************************************************************
!
! Append one frame. have_time says whether the comment line carried a time=
! tag; when it did not, t is stored as a placeholder and have_times latches
! false so that ir_fft_frames_dt can say which frames were unlabelled.
!
   subroutine ir_fft_frames_push(this, mu, time_fs, have_time)

      implicit none

      type(ir_fft_frames_type), intent(inout) :: this
      real(dp), intent(in) :: mu(1:3)
      real(dp), intent(in) :: time_fs
      logical, intent(in) :: have_time
      real(dp), allocatable :: tmp_mu(:, :), tmp_t(:)
      integer :: new_cap

      if (this%n >= this%cap) then
         new_cap = max(2*this%cap, 1024)
         allocate (tmp_mu(1:3, 1:new_cap), tmp_t(1:new_cap))
         if (this%n > 0) then
            tmp_mu(1:3, 1:this%n) = this%mu(1:3, 1:this%n)
            tmp_t(1:this%n) = this%t(1:this%n)
         end if
         if (allocated(this%mu)) deallocate (this%mu)
         if (allocated(this%t)) deallocate (this%t)
         call move_alloc(tmp_mu, this%mu)
         call move_alloc(tmp_t, this%t)
         this%cap = new_cap
      end if

      this%n = this%n + 1
      this%mu(1:3, this%n) = mu(1:3)
      this%t(this%n) = time_fs
      if (.not. have_time) then
         this%have_times = .false.
         this%n_missing = this%n_missing + 1
      end if

   end subroutine ir_fft_frames_push

!**************************************************************************
!
! The interval between frames, from the labels, with the uniformity check
! that makes taking it from the labels worth doing.
!
!   dt = (t(n) - t(1)) / (n - 1)
!
! and every consecutive spacing must agree with that to within tol (relative).
! The default tolerance is loose -- 1e-3, i.e. a tenth of a per cent -- because
! the labels are written with finite precision (write_extxyz uses F16.6, so a
! 0.5 fs step is exact but a 1/3 fs step is not) and the check is meant to
! catch a trajectory that is missing frames or was concatenated from two runs,
! not to police the last digit.
!
! dt_fallback is used, with a warning left in msg, only when NO frame carried
! a label. That is the case of a hand-built xyz, and refusing it outright
! would be unhelpful; refusing a PARTLY labelled one is not, because there the
! file is telling us something inconsistent about itself.
!
   subroutine ir_fft_frames_dt(this, dt_fallback, tol, dt, ok, msg)

      implicit none

      type(ir_fft_frames_type), intent(in) :: this
      real(dp), intent(in) :: dt_fallback
      real(dp), intent(in) :: tol
      real(dp), intent(out) :: dt
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: gap, worst, rel
      integer :: i, i_worst

      ok = .false.
      dt = 0.d0
      msg = ""

      if (this%n < 2) then
         msg = "ir_fft: a spectrum needs at least two frames"
         return
      end if

      if (.not. this%have_times) then
         if (this%n_missing < this%n) then
            write (msg, '(A,I0,A,I0,A)') &
               "ir_fft: ", this%n_missing, " of ", this%n, &
               " frames carry no time= tag. A partly labelled trajectory "// &
               "cannot be checked for uniform spacing; either label them all "// &
               "or none, and set ir_frame_dt for the unlabelled case"
            return
         end if
         if (dt_fallback <= 0.d0) then
            msg = "ir_fft: no frame carries a time= tag and ir_frame_dt was "// &
                  "not set, so there is no time axis. Give ir_frame_dt the "// &
                  "interval between frames in fs"
            return
         end if
         dt = dt_fallback
         write (msg, '(A,F12.6,A)') &
            "no time= tags in the trajectory; using ir_frame_dt = ", dt, &
            " fs. The wavenumber axis is only as right as that number."
         ok = .true.
         return
      end if

      dt = (this%t(this%n) - this%t(1))/dfloat(this%n - 1)
      if (dt <= 0.d0) then
         write (msg, '(A,F16.6,A,F16.6,A)') &
            "ir_fft: the time labels do not increase: t(1) = ", this%t(1), &
            " and t(n) = ", this%t(this%n), &
            ". Is the trajectory in order?"
         return
      end if

      worst = 0.d0
      i_worst = 1
      do i = 2, this%n
         gap = this%t(i) - this%t(i - 1)
         rel = dabs(gap - dt)/dt
         if (rel > worst) then
            worst = rel
            i_worst = i
         end if
      end do

      if (worst > tol) then
         write (msg, '(A,I0,A,F16.6,A,F16.6,A)') &
            "ir_fft: the frames are not evenly spaced. Frame ", i_worst, &
            " is ", this%t(i_worst) - this%t(i_worst - 1), &
            " fs after its predecessor where the mean spacing is ", dt, &
            " fs. An unevenly sampled series has no spectrum of this kind; "// &
            "trim the trajectory or raise ir_frame_dt_tol if the labels are "// &
            "merely imprecise."
         return
      end if

      ok = .true.

   end subroutine ir_fft_frames_dt

!**************************************************************************
!
! Copy mad_ir's circular buffer out in chronological order, oldest first, so
! that mu(:, n) is the newest frame -- the one ir_fft_loss differentiates with
! respect to.
!
! head is the slot holding the newest frame and ages run backwards from it,
! which is mad_ir_push's convention; n_stored caps at n_window.
!
   subroutine ir_fft_md_unroll(mu_hist, n_window, n_stored, head, mu, n)

      implicit none

      integer, intent(in) :: n_window, n_stored, head
      real(dp), intent(in) :: mu_hist(1:3, 1:n_window)
      real(dp), intent(out) :: mu(1:3, 1:n_stored)
      integer, intent(out) :: n
      integer :: j, a, m

      n = n_stored
      do j = 1, n_stored
!        age of the j-th chronological frame: the newest (j = n) has age 0.
         a = n_stored - j
         m = head - a
         do while (m < 1)
            m = m + n_window
         end do
         mu(1:3, j) = mu_hist(1:3, m)
      end do

   end subroutine ir_fft_md_unroll

!**************************************************************************
!
! Assemble a config from the parameters the input file carries. Kept here
! rather than in turbogap.f90 so that the mapping from keyword to field is in
! one place and the test program can use the same one.
!
! temperature <= 0 means "take the run's target temperature", which is what
! t_beg is; the harmonic correction needs a number and the run already knows
! one, so making the user repeat it is an invitation to have them disagree.
!
   subroutine ir_fft_config_from_params(dt_fs, window, acf_ratio, max_freq_cm, &
                                        smooth_k, smooth_kind, quantum_correction, &
                                        temperature, t_beg, power_dc_cutoff_cm, &
                                        subtract_mean, normalise, cfg)

      implicit none

      real(dp), intent(in) :: dt_fs, acf_ratio, max_freq_cm, temperature, t_beg
      real(dp), intent(in) :: power_dc_cutoff_cm
      character(len=*), intent(in) :: window, smooth_kind, quantum_correction
      integer, intent(in) :: smooth_k
      logical, intent(in) :: subtract_mean, normalise
      type(ir_fft_config_type), intent(out) :: cfg

      cfg%dt_fs = dt_fs
      cfg%window = window
      cfg%acf_ratio = acf_ratio
      cfg%max_freq_cm = max_freq_cm
      cfg%smooth_k = smooth_k
      cfg%smooth_kind = smooth_kind
      cfg%quantum_correction = quantum_correction
      if (temperature > 0.d0) then
         cfg%temperature = temperature
      else
         cfg%temperature = t_beg
      end if
      cfg%power_dc_cutoff_cm = power_dc_cutoff_cm
      cfg%subtract_mean = subtract_mean
      cfg%normalise = normalise

   end subroutine ir_fft_config_from_params

!**************************************************************************
!
! The dipole trajectory the spectrum was computed from.
!
   subroutine ir_fft_write_dipoles(this, fname)

      implicit none

      type(ir_fft_frames_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      integer :: u, i

      open (newunit=u, file=trim(fname), status="replace", action="write")
      write (u, '(A)') "# The total dipole per frame, as predicted by the dipole model."
      write (u, '(A)') "# Columns 2-4 are what ir_fft_spectrum.dat was computed from, so"
      write (u, '(A)') "# TNEP/spectroscopy.py compute_ir_spectrum() on them reproduces it."
      write (u, '(A)') "#"
      write (u, '(A)') "#     time_fs                mu_x                 mu_y                 mu_z"
      do i = 1, this%n
         write (u, '(F16.4,3ES22.13)') this%t(i), this%mu(1:3, i)
      end do
      close (u)

   end subroutine ir_fft_write_dipoles

!**************************************************************************
!
! The spectrum, with enough provenance in the header to read it.
!
! Everything in that header is a limit of the calculation rather than a
! property of the sample, and every one of them has been mistaken for physics
! at some point:
!
!   nyquist      above it nothing is represented; power there FOLDS BACK down
!                into the range being reported rather than going missing.
!   resolution   set by the longest lag kept, i.e. by acf_ratio and the run
!                length, and NOT by the bin spacing. The bins are twice as
!                fine as the resolution by construction (N2 = 2L-1), so a
!                feature narrower than `resolution` is the lag window, not a
!                mode.
!   smoothing    broadens further, by roughly smooth_k bins.
!   peak         the divisor. The intensity column is dimensionless; multiply
!                by it to get back what the transform produced.
!
! With an experiment present the fitted scale and the experimental curve
! interpolated onto the same grid are written too, so the file is
! self-contained for plotting.
!
   subroutine ir_fft_write_spectrum(res, cfg, fname, nu_exp, I_exp, n_exp, &
                                    has_exp, scale, offset, dissim, dissim_ref, &
                                    label)

      implicit none

      type(ir_fft_result_type), intent(in) :: res
      type(ir_fft_config_type), intent(in) :: cfg
      character(len=*), intent(in) :: fname
      integer, intent(in) :: n_exp
      real(dp), intent(in) :: nu_exp(1:n_exp), I_exp(1:n_exp)
      logical, intent(in) :: has_exp
      real(dp), intent(in) :: scale, offset, dissim, dissim_ref
      character(len=*), intent(in) :: label
      real(dp) :: e_here, x, tt, rel
      integer :: u, k, j, kk

      open (newunit=u, file=trim(fname), status="replace", action="write")
      write (u, '(A)') "# IR spectrum, FFT estimator (ir_fft.f90), a translation of"
      write (u, '(A)') "# TNEP/spectroscopy.py -- GPUMD / Xu et al. JCTC 20, 3273 (2024)."
      if (len_trim(label) > 0) write (u, '(A,A)') "# ", trim(label)
      write (u, '(A)') "#"
      write (u, '(A,I0)') "#   frames               : ", res%n_frames
      write (u, '(A,F14.6,A)') "#   frame interval       : ", res%dt_fs, " fs"
      write (u, '(A,I0,A,F12.4,A)') "#   lags kept            : ", res%n_lag, &
         "   (acf_ratio = ", cfg%acf_ratio, ")"
      write (u, '(A,F14.4,A)') "#   resolution           : ", res%resolution, &
         " cm^-1   <- set by the longest lag, not by the bin spacing"
      write (u, '(A,F14.4,A)') "#   bin spacing          : ", res%d_nu, " cm^-1"
      write (u, '(A,F14.1,A)') "#   Nyquist              : ", res%nyquist, &
         " cm^-1   <- power above this folds back down into the range below"
      write (u, '(A,A)') "#   lag window           : ", trim(cfg%window)
      write (u, '(A,I0,A,A,A)') "#   smoothing            : ", cfg%smooth_k, &
         " bins, ", trim(cfg%smooth_kind), &
         "   <- broadens beyond the resolution above"
      write (u, '(A,A)') "#   quantum correction   : ", trim(cfg%quantum_correction)
      if (trim(cfg%quantum_correction) == "harmonic") then
         write (u, '(A,F10.2,A)') "#   temperature          : ", cfg%temperature, " K"
      end if
      write (u, '(A,3ES16.7)') "#   mean dipole removed  : ", res%mu_mean(1:3)
!     Say whether the divisor was APPLIED, not just what it was. The MAD path
!     turns normalisation off (see the ir_fft.f90 header on why max|I| is not a
!     thing to differentiate), and a header reporting a divisor next to a column
!     that has not been divided is exactly the kind of half-truth that gets
!     plotted without being read.
      if (cfg%normalise) then
         write (u, '(A,ES20.10)') "#   intensity divided by : ", res%peak
         write (u, '(A,ES20.10)') "#   power divided by     : ", res%peak_power
      else
         write (u, '(A)') "#   NOT peak-normalised. The intensity and power columns are the raw"
         write (u, '(A)') "#   transform; the fitted scale below is the comparison with experiment."
         write (u, '(A,ES20.10)') "#   (peak would have been: ", res%peak
         write (u, '(A,ES20.10)') "#    power peak would have been: ", res%peak_power
      end if
      if (has_exp) then
         write (u, '(A)') "#"
         write (u, '(A,ES16.7)') "#   fitted scale         : ", scale
         write (u, '(A,ES16.7)') "#   fitted offset        : ", offset
         write (u, '(A,ES16.7)') "#   dissimilarity        : ", dissim
         rel = 0.d0
         if (dissim_ref > 0.d0) rel = dsqrt(dissim/dissim_ref)
         write (u, '(A,ES16.7,A)') "#   relative mismatch    : ", rel, &
            "   = sqrt(dissim / sum w I_exp^2)"
      end if
      write (u, '(A)') "#"
      if (has_exp) then
         write (u, '(A)') "#           nu        intensity            power    "// &
            "intensity_raw        power_raw   fitted_vs_exp     experiment"
      else
         write (u, '(A)') "#           nu        intensity            power    "// &
            "intensity_raw        power_raw"
      end if

      do k = 1, res%n_freq
         if (has_exp) then
!           The experiment interpolated onto the spectrum's grid, and the
!           fitted prediction beside it, so the two columns can be plotted
!           against each other without any further arithmetic.
            e_here = 0.d0
            if (n_exp >= 2) then
               if (res%freq(k) <= nu_exp(1)) then
                  e_here = I_exp(1)
               else if (res%freq(k) >= nu_exp(n_exp)) then
                  e_here = I_exp(n_exp)
               else
                  kk = 1
                  do j = 1, n_exp - 1
                     if (res%freq(k) >= nu_exp(j) .and. res%freq(k) <= nu_exp(j + 1)) then
                        kk = j
                        exit
                     end if
                  end do
                  x = nu_exp(kk + 1) - nu_exp(kk)
                  tt = 0.d0
                  if (x > 0.d0) tt = (res%freq(k) - nu_exp(kk))/x
                  e_here = (1.d0 - tt)*I_exp(kk) + tt*I_exp(kk + 1)
               end if
            end if
            write (u, '(7ES17.8)') res%freq(k), res%intensity(k), res%power(k), &
               res%intensity_raw(k), res%power_raw(k), &
               scale*res%intensity_raw(k) + offset, e_here
         else
            write (u, '(5ES17.8)') res%freq(k), res%intensity(k), res%power(k), &
               res%intensity_raw(k), res%power_raw(k)
         end if
      end do
      close (u)

   end subroutine ir_fft_write_spectrum

!**************************************************************************
!
! End of a prediction run: work out the timestep, transform, compare with the
! experiment if there is one, and write both files.
!
! Called once, on rank 0, after the last frame. Everything it needs is in the
! frame buffer and in cfg apart from the experiment, which is optional.
!
   subroutine ir_fft_frames_finish(this, cfg_in, dt_fallback, dt_tol, &
                                   nu_exp, I_exp, wgt, n_exp, has_exp, &
                                   match_scale, match_offset, &
                                   spectrum_file, dipole_file, write_files, &
                                   write_dipoles, &
                                   res, dt_used, scale, offset, dissim, dissim_ref, &
                                   ok, msg)

      implicit none

      type(ir_fft_frames_type), intent(in) :: this
      type(ir_fft_config_type), intent(in) :: cfg_in
      real(dp), intent(in) :: dt_fallback, dt_tol
      integer, intent(in) :: n_exp
      real(dp), intent(in) :: nu_exp(1:n_exp), I_exp(1:n_exp), wgt(1:n_exp)
      logical, intent(in) :: has_exp, match_scale, match_offset
      character(len=*), intent(in) :: spectrum_file, dipole_file
!     write_files gates the I/O on rank 0. EVERY rank runs the transform --
!     the ensemble is replicated and the arithmetic is identical, so there is
!     nothing to reduce and nothing to broadcast -- but exactly one of them may
!     open the output, or they race on the same path and the loser truncates
!     the winner's file.
      logical, intent(in) :: write_files
      logical, intent(in) :: write_dipoles
      type(ir_fft_result_type), intent(inout) :: res
      real(dp), intent(out) :: dt_used, scale, offset, dissim, dissim_ref
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      type(ir_fft_config_type) :: cfg
      real(dp), allocatable :: I_fit(:)
      real(dp) :: energy, lambda(1:3)
      logical :: ok2
      character(len=512) :: msg2, msg_dt

      ok = .false.
      scale = 1.d0
      offset = 0.d0
      dissim = 0.d0
      dissim_ref = 0.d0
      dt_used = 0.d0
      msg = ""

      call ir_fft_frames_dt(this, dt_fallback, dt_tol, dt_used, ok2, msg_dt)
      if (.not. ok2) then
         msg = trim(msg_dt)
         return
      end if

      cfg = cfg_in
      cfg%dt_fs = dt_used

      call ir_fft_spectrum(this%mu, this%n, cfg, res, ok2, msg2)
      if (.not. ok2) then
         msg = trim(msg2)
         return
      end if

!     With an experiment, run the loss too. energy_scale is zero -- this is a
!     prediction, nothing is being biased -- which makes ir_fft_loss skip the
!     adjoint and return only the fit and the mismatch, which is all that is
!     wanted here.
      if (has_exp .and. n_exp >= 2) then
         allocate (I_fit(1:n_exp))
         call ir_fft_loss(this%mu, this%n, cfg, nu_exp, I_exp, wgt, n_exp, &
                          match_scale, match_offset, 0.d0, energy, lambda, &
                          I_fit, scale, offset, dissim, dissim_ref, ok2, msg2)
         deallocate (I_fit)
         if (.not. ok2) then
            msg = trim(msg2)
            return
         end if
      end if

      if (write_files) then
         call ir_fft_write_spectrum(res, cfg, spectrum_file, nu_exp, I_exp, n_exp, &
                                    has_exp .and. n_exp >= 2, scale, offset, &
                                    dissim, dissim_ref, trim(msg_dt))
         if (write_dipoles) call ir_fft_write_dipoles(this, dipole_file)
      end if

      msg = trim(msg_dt)
      ok = .true.

   end subroutine ir_fft_frames_finish

end module ir_fft_io
