! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_ir.f90, is copyright (c) 2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!  The driver's side of the IR spectrum: which route a run takes, the state of
!  the bias across steps, and what it cost. The estimators themselves live in
!  mad_ir, mad_ir_xl, ir_auxiliary_dynamics and ir_fft.
module turbogap_ir

   use kinds, only: dp
   use types, only: input_parameters
   use timing, only: times_t, time_start, time_end, get_time
   use error, only: turbogap_abort
   use turbogap_comm, only: comm_t, comm_sum_all
   use turbogap_structure, only: state_t
   use turbogap_results, only: results_t
   use turbogap_loop, only: loop_t
   use exp_utils, only: exp_dissimilarity, exp_dissim_ref, get_energy_scale
   use exp_interface, only: get_write_condition
   use mad_ir, only: CM_PER_INV_FS, mad_ir_state, mad_ir_collect, mad_ir_need_dmu, mad_ir_dmu_dr, &
                     mad_ir_setup, mad_ir_setup_predict, mad_ir_push, mad_ir_ready, mad_ir_evaluate, &
                     mad_ir_forces, mad_ir_save, mad_ir_spectrum, mad_ir_write_spectrum, &
                     mad_ir_append_spectrum, mad_ir_write_exp_spectrum, mad_ir_select_range
   use mad_ir_xl, only: mad_ir_xl_state, mad_ir_xl_active, mad_ir_xl_collect, mad_ir_xl_force, &
                        mad_ir_xl_site_w, mad_ir_xl_setup, mad_ir_xl_advance, mad_ir_xl_evaluate, &
                        mad_ir_xl_weights, mad_ir_xl_ready, mad_ir_xl_save, mad_ir_xl_write_spectrum, &
                        mad_ir_xl_resolution, mad_ir_xl_memory_bytes
   use ir_auxiliary_dynamics, only: ir_aux_state, ir_aux_active, ir_aux_setup, ir_aux_advance, &
                                    ir_aux_evaluate, ir_aux_forces, ir_aux_calibrate, ir_aux_calibrated, &
                                    ir_aux_save, ir_aux_write_spectrum, ir_aux_escale_max
   use ir_fft, only: ir_fft_config_type, ir_fft_result_type, ir_fft_loss, ir_fft_spectrum, ir_fft_free
   use ir_fft_io, only: ir_fft_frames, ir_fft_config_from_params, ir_fft_frames_reset, ir_fft_frames_push, &
                        ir_fft_frames_finish, ir_fft_md_unroll, ir_fft_write_spectrum

   implicit none

   private
   public :: ir_init
   public :: ir_step_begin
   public :: ir_before_evaluate
   public :: ir_push_frame
   public :: ir_after_forces
   public :: ir_step_end
   public :: ir_finish
   public :: ir_report

   type, public :: ir_run_t
!  MAD IR bias. lambda is dL/dmu of the newest configuration; mad_ir_applied
!  says whether the ensemble was full enough for a force to have been added.
      real(dp) :: mad_ir_lambda(1:3) = 0.d0
      real(dp) :: mad_ir_energy = 0.d0
      real(dp) :: mad_ir_scale = 0.d0
      logical :: mad_ir_applied = .false.
      logical :: mad_ir_ok
      logical :: mad_ir_resumed
!  The envelope-targeted bank, ir_bias_mode = aux. Calibrated once, when the
!  ACF ensemble first fills, so the flag says whether that has happened yet.
      logical :: ir_aux_ok
      logical :: ir_aux_resumed
      character(len=512) :: ir_aux_msg
!  The largest back-reaction scale the run will reach, and the temperature the
!  ensemble was collected at: the two inputs to the stability bound.
      real(dp) :: ir_aux_escale_top = 0.d0
      real(dp) :: ir_aux_temp = 0.d0
!  mad_ir_evaluate is called under aux purely for the INDEPENDENT dissimilarity
!  it leaves behind, so its energy is discarded here.
      real(dp) :: ir_aux_acf_energy = 0.d0
!  The SECOND constraint. Lambda < 1 stops the bilinear runaway, but the
!  controller still pumps the bank and the bank still pumps the atoms, and a
!  thermostat of time constant tau_t absorbs a steady power P only at the cost
!  of a standing temperature offset
!
!     dT = 2 P tau_t / (3 N kB)
!
!  which is what a bias that is stable but too strong looks like: bounded, and
!  a thousand degrees hot. Measured here from the pumped-energy ledger rather
!  than predicted, because P depends on the dynamics in a way Lambda does not.
      real(dp) :: ir_aux_power = 0.d0
      real(dp) :: ir_aux_work = 0.d0
      real(dp) :: ir_aux_pump_prev = 0.d0
      real(dp) :: ir_aux_pump_time = 0.d0
      real(dp) :: ir_aux_dT = 0.d0
      logical :: ir_aux_warned_hot = .false.
      character(len=512) :: mad_ir_msg
!  Has anything been appended to ir_prediction.dat yet? The first block cannot
!  be identified by its step number the way the per-frame observables' can --
!  it appears whenever the ensemble first fills, which is not step zero -- and
!  appending to a file that does not exist is a runtime error, so the flag is
!  carried rather than derived.
      logical :: mad_ir_wrote_prediction = .false.
!  The extended-Lagrangian bias (ir_bias_mode = "xl"). mad_ir_xl_ok/msg mirror
!  the ACF ones; mad_ir_xl_resumed says whether a saved resonator bank was
!  adopted, which is the difference between biasing from the first stored frame
!  and charging a fresh bank for ir_xl_warm_factor memory times first.
      logical :: mad_ir_xl_ok = .false.
      logical :: mad_ir_xl_resumed = .false.
      character(len=512) :: mad_ir_xl_msg
!  What the bias costs. The ensemble fills partway into the run, so the same
!  run measures both sides of the question: mad_ir_t_pre accumulates the
!  wall-clock of the steps before the first spectrum and mad_ir_t_post that of
!  the steps after, and their per-step means are directly comparable because
!  nothing else about the step changes at that boundary.
      real(dp) :: mad_ir_t_first = -1.d0
      integer :: mad_ir_step_first = -1
      real(dp) :: mad_ir_t_pre = 0.d0
      real(dp) :: mad_ir_t_post = 0.d0
      integer :: mad_ir_n_pre = 0
      integer :: mad_ir_n_post = 0
      real(dp) :: mad_ir_step_beg = 0.d0
      real(dp) :: mad_ir_t_now = 0.d0
      real(dp) :: mad_ir_rate_pre
      real(dp) :: mad_ir_rate_post
      real(dp) :: mad_ir_res_ask
      logical :: mad_ir_have_spectrum = .false.
!  The FFT estimator (ir_bias_mode = fft, and the only estimator `turbogap
!  predict` can use). ir_from_traj distinguishes the two ways in: reading a
!  trajectory off disk frame by frame, versus riding along on an MD run and
!  transforming mad_ir's rolling buffer.
      logical :: ir_fft_active = .false.
      logical :: ir_from_traj = .false.
      logical :: ir_fft_ok
      character(len=1024) :: ir_fft_msg
      type(ir_fft_config_type) :: ir_fft_cfg
      type(ir_fft_result_type) :: ir_fft_res
      real(dp) :: ir_fft_dt_used = 0.d0
      real(dp) :: ir_fft_scale_fit = 1.d0
      real(dp) :: ir_fft_offset_fit = 0.d0
      real(dp) :: ir_fft_dissim = 0.d0
      real(dp) :: ir_fft_dissim_ref = 0.d0
      real(dp), allocatable :: ir_fft_mu_chron(:, :)
      real(dp), allocatable :: ir_nu_exp(:)
      real(dp), allocatable :: ir_I_exp(:)
      real(dp), allocatable :: ir_wgt_exp(:)
      real(dp), allocatable :: ir_fft_I_fit(:)
      integer :: ir_fft_n_chron = 0
   end type ir_run_t

contains

!  Decide which IR route the run takes and set up what it needs before the
!  first frame.
   subroutine ir_init(ir, params, comm)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(inout) :: params
      type(comm_t), intent(in) :: comm

!  WHICH IR ROUTE. There are two, and they are not variations of one thing.
!
!  ir_from_traj: `turbogap predict` with do_ir. The configurations are already
!  on disk and are read one at a time; the ensemble is the file, its length is
!  discovered by reaching the end of it, and the sampling interval comes from
!  the time= tags rather than from md_step. Nothing is biased -- there is no
!  dynamics to bias -- so this is prediction only, and it uses ir_fft.f90
!  because that is the estimator whose sizing is an output rather than an
!  input. mad_ir's ring buffer is not set up at all: it wants n_window before
!  the first frame, and in this mode nobody knows it.
!
!  Otherwise: MD or MC, where mad_ir sizes and fills a rolling buffer as
!  before and ir_bias_mode picks which estimator transforms it.
      ir%ir_from_traj = params%do_ir .and. .not. params%do_md .and. .not. params%do_mc
      ir%ir_fft_active = ir%ir_from_traj .or. &
                         ((params%valid_ir .or. params%do_ir) .and. &
                          trim(params%ir_bias_mode) == "fft")

!  One config for both routes; only dt_fs differs, and it is filled in at the
!  point of use because in one route it comes from the file and in the other
!  from md_step*ir_stride. normalise is off for MD: under a bias the spectrum
!  is fitted against the experiment with a scale, and dividing by max|I| as
!  well would be a second, discontinuous, normalisation of the same freedom.
      if (ir%ir_fft_active) then
         call ir_fft_config_from_params(1.d0, params%ir_window, params%ir_fft_acf_ratio, &
                                        params%ir_nu_max, params%ir_fft_smooth_k, &
                                        params%ir_fft_smooth_kind, &
                                        params%ir_fft_quantum_correction, &
                                        params%ir_fft_temperature, params%t_beg, &
                                        params%ir_fft_power_dc_cutoff, &
                                        params%ir_subtract_mean, ir%ir_from_traj, ir%ir_fft_cfg)
      end if

      if (ir%ir_from_traj) then
!     Without do_prediction the descriptor pass never runs, so no dipole is
!     ever formed and the frame buffer stays empty. That surfaces much later as
!     "a spectrum needs at least two frames", which is true and unhelpful.
         if (.not. params%do_prediction) then
            write (*, *) "ERROR: do_ir in predict mode needs do_prediction = .true."
            write (*, *) "       Without it no descriptor is evaluated and no dipole exists."
            stop 1
         end if
         if (.not. params%do_dipole) then
            write (*, *) "ERROR: do_ir in predict mode needs a dipole model. Add"
            write (*, *) "       dipole_model = .true. to one of the soap_turbo blocks."
            stop 1
         end if
         if (trim(params%ir_bias_mode) /= "fft" .and. trim(params%ir_bias_mode) /= "acf") then
            write (*, *) "ERROR: ir_bias_mode = ", trim(params%ir_bias_mode), &
               " has no meaning in predict mode."
            write (*, *) "       A trajectory read from disk is transformed by the fft"
            write (*, *) "       estimator; leave ir_bias_mode unset or set it to fft."
            stop 1
         end if
         call ir_fft_frames_reset(ir_fft_frames)
         if (comm%rank == 0) then
            write (*, *) '                                       |'
            write (*, *) 'IR prediction from a trajectory:       |'
            write (*, '(A,A)') '  *) estimator:         fft (ir_fft.f90)  |'
            write (*, '(A,A20,A)') '  *) lag window:     ', trim(params%ir_window), '  |'
            write (*, '(A,F12.4,A)') '  *) acf_ratio:         ', params%ir_fft_acf_ratio, '         |'
            write (*, '(A,A20,A)') '  *) quantum corr.:  ', trim(params%ir_fft_quantum_correction), '  |'
            write (*, '(A,I12,A)') '  *) smoothing bins:    ', params%ir_fft_smooth_k, '         |'
            write (*, '(A,F12.1,A)') '  *) nu_max:            ', params%ir_nu_max, ' cm^-1   |'
            write (*, *) '  *) sizing follows the file; see       |'
            write (*, *) '     ir_fft_spectrum.dat at the end.    |'
            write (*, *) '                                       |'
         end if
      end if
   end subroutine ir_init

   subroutine ir_step_begin(ir, params)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(in) :: params

!     One stamp per iteration, so that the cost of a step with the IR bias
!     active can be compared with the cost of one without. Closed at the
!     bottom of the loop.
      if (params%valid_ir) call get_time(ir%mad_ir_step_beg)
   end subroutine ir_step_begin

!  Before the descriptors: set up on first use, and decide whether this step
!  contributes a configuration, since get_soap must know before it builds.
   subroutine ir_before_evaluate(ir, params, state, res, loop, comm)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(inout) :: params
      type(state_t), intent(in) :: state
      type(results_t), intent(inout) :: res
      type(loop_t), intent(in) :: loop
      type(comm_t), intent(in) :: comm

!           Decide here, before any descriptor is evaluated, whether this step
!           contributes a configuration to the IR ensemble: get_soap has to be
!           told to produce second derivatives before it builds anything, and
!           gap_interface reads mad_ir_collect to do that.
      if ((params%valid_ir .or. params%do_ir) .and. .not. ir%ir_from_traj) then
!              Set up on first use: n_sites is known by now, and doing it here
!              rather than in the setup phase keeps the sizing next to the
!              place that consumes it.
         if (.not. mad_ir_state%active) then
            if (.not. params%do_dipole) then
               write (*, *) "ERROR: an ir observable needs a dipole model. Add"
               write (*, *) "       dipole_model = .true. to one of the soap_turbo blocks."
               stop
            end if
            if (params%soap_radial_legacy_filter .and. params%valid_ir) then
               write (*, *) "ERROR: an ir observable needs the descriptor second derivatives,"
               write (*, *) "       which require soap_radial_legacy_filter = .false."
               stop
            end if
            if (params%valid_ir) then
               call mad_ir_setup(params%md_step, params%ir_stride, params%ir_resolution, &
                                 params%ir_nu_min, params%ir_nu_max, params%ir_lag_factor, &
                                 params%exp_data(params%ir_idx)%data(1, :), &
                                 params%exp_data(params%ir_idx)%data(2, :), params%ir_restart_file, &
                                 params%ir_match_scale, params%ir_nu_power, &
                                 params%ir_window, params%ir_subtract_mean, &
                                 trim(params%ir_estimator) /= "unbiased", &
                                 params%ir_taper_partial, params%ir_match_offset, &
                                 params%ir_weight_by_spacing, &
                                 state%n_sites, ir%mad_ir_ok, ir%mad_ir_resumed, ir%mad_ir_msg, &
                                 params%ir_acf_mode, params%ir_tau_mem)
            else
!                    Prediction: the ensemble is the whole trajectory, so its
!                    length is md_nsteps/ir_stride + 1 -- every step for which
!                    modulo(md_istep, ir_stride) is zero, counting step zero.
               ir%mad_ir_resumed = .false.
!                    A negative resolution means "whatever the run gives"; the
!                    default value of ir_resolution cannot be told from a
!                    chosen one by its value, hence the flag.
               if (params%ir_resolution_set) then
                  ir%mad_ir_res_ask = params%ir_resolution
               else
                  ir%mad_ir_res_ask = -1.d0
               end if
               call mad_ir_setup_predict(params%md_step, params%ir_stride, &
                                         params%md_nsteps/params%ir_stride + 1, &
                                         ir%mad_ir_res_ask, params%ir_nu_min, &
                                         params%ir_nu_max, params%ir_lag_factor, &
                                         params%ir_n_samples, params%ir_nu_power, &
                                         params%ir_window, params%ir_subtract_mean, &
                                         trim(params%ir_estimator) /= "unbiased", &
                                         params%ir_taper_partial, &
                                         state%n_sites, ir%mad_ir_ok, ir%mad_ir_msg, &
                                         params%ir_acf_mode, params%ir_tau_mem)
            end if
            if (.not. ir%mad_ir_ok) then
               write (*, *) "ERROR: ", trim(ir%mad_ir_msg)
               stop
            end if
            if (comm%rank == 0) then
               write (*, *) '                                       |'
               if (params%valid_ir) then
                  write (*, *) 'MAD IR bias:                           |'
               else
                  write (*, *) 'IR prediction:                         |'
               end if
               write (*, '(A,F12.4,A)') '  *) sampling interval: ', mad_ir_state%dt, ' fs      |'
               write (*, '(A,I12,A)') '  *) longest lag:       ', mad_ir_state%n_lag, '         |'
               write (*, '(A,I12,A)') '  *) ensemble size:     ', mad_ir_state%n_window, '         |'
               write (*, '(A,F12.4,A)') '  *) resolution:        ', &
                  CM_PER_INV_FS/(dfloat(mad_ir_state%n_lag)*mad_ir_state%dt), ' cm^-1   |'
               write (*, '(A,F12.1,A)') '  *) Nyquist:           ', &
                  CM_PER_INV_FS/(2.d0*mad_ir_state%dt), ' cm^-1   |'
               write (*, '(A,I12,A)') '  *) spectrum points:   ', mad_ir_state%n_freq, '         |'
               if (.not. params%valid_ir) then
                  write (*, *) '  *) no experiment and no bias; the   |'
                  write (*, *) '     spectrum is written at the end.  |'
               else if (ir%mad_ir_resumed) then
                  write (*, '(A,I12,A)') '  *) resumed, frames:   ', mad_ir_state%n_stored, '         |'
               else
                  write (*, *) '  *) fresh ensemble; no bias is       |'
                  write (*, *) '     applied until it is full.        |'
                  if (len_trim(ir%mad_ir_msg) > 0) write (*, *) '     ', trim(ir%mad_ir_msg)
               end if
               write (*, *) '.......................................|'
            end if
!                 The resonator bank, if this is an extended-Lagrangian run. It
!                 is set up after the ACF observable rather than instead of it
!                 because it takes the fitted grid from mad_ir_state: the
!                 experiment is read once, by read_exp_data, and a second
!                 reader of the same file is exactly the sort of divergence
!                 that shows up months later as a spectrum that does not match
!                 the one the header quoted.
            if (trim(params%ir_bias_mode) == "xl") then
               if (.not. params%valid_ir) then
                  write (*, *) "ERROR: ir_bias_mode = xl needs an experimental spectrum."
                  write (*, *) "       Add ir to exp_labels with a file in exp_data_files,"
                  write (*, *) "       or leave ir_bias_mode at acf for a prediction run."
                  stop
               end if
               call mad_ir_xl_setup(mad_ir_state, state%n_sites, params%ir_xl_n_modes, &
                                    params%md_step, params%ir_stride, &
                                    params%ir_xl_tau_mem, &
                                    trim(params%ir_xl_amplitude) == "coherent", &
                                    params%ir_xl_warm_factor, &
                                    params%ir_xl_max_memory*1.048576d6, &
                                    params%ir_xl_restart_file, ir%mad_ir_xl_ok, &
                                    ir%mad_ir_xl_resumed, ir%mad_ir_xl_msg)
               if (.not. ir%mad_ir_xl_ok) then
                  write (*, *) "ERROR: ", trim(ir%mad_ir_xl_msg)
                  stop
               end if
               mad_ir_xl_active = .true.
               if (comm%rank == 0) then
                  write (*, *) 'MAD IR, extended Lagrangian:           |'
                  write (*, '(A,I12,A)') '  *) resonator modes:   ', &
                     mad_ir_xl_state%n_modes, '         |'
                  write (*, '(A,F12.4,A)') '  *) memory time:       ', &
                     mad_ir_xl_state%tau_mem, ' fs      |'
                  write (*, '(A,F12.4,A)') '  *) bank resolution:   ', &
                     mad_ir_xl_resolution(mad_ir_xl_state%tau_mem), ' cm^-1   |'
                  write (*, '(A,F12.3,A)') '  *) bank memory/rank:  ', &
                     mad_ir_xl_memory_bytes(mad_ir_xl_state%n_modes, &
                                            mad_ir_xl_state%n_bank)/1.048576d6, &
                     ' MB      |'
                  if (mad_ir_xl_state%coherent) then
                     write (*, *) '  *) amplitude:            coherent   |'
                  else
                     write (*, *) '  *) amplitude:          incoherent   |'
                  end if
                  if (ir%mad_ir_xl_resumed) then
                     write (*, '(A,I12,A)') '  *) resumed, advances: ', &
                        mad_ir_xl_state%n_steps, '         |'
                  else
                     write (*, '(A,I12,A)') '  *) charging, advances:', &
                        mad_ir_xl_state%n_warm, '         |'
                     if (len_trim(ir%mad_ir_xl_msg) > 0) write (*, *) '     ', trim(ir%mad_ir_xl_msg)
                  end if
                  write (*, *) '  *) the ACF ensemble above is still  |'
                  write (*, *) '     filled, but only as a check: it  |'
                  write (*, *) '     produces no force in this mode.  |'
                  write (*, *) '.......................................|'
               end if
            end if
!                 The envelope-targeted bank. Set up here for the same reason
!                 the extended-Lagrangian one is: it takes the fitted grid from
!                 mad_ir_state, so the experiment is still read exactly once.
!                 What it does NOT do here is calibrate -- the couplings come
!                 from the autocorrelation by linear response, and there is no
!                 autocorrelation until the ensemble fills. That happens in the
!                 bias block below, on the first step mad_ir_ready holds.
            if (trim(params%ir_bias_mode) == "aux") then
               if (.not. params%valid_ir) then
                  write (*, *) "ERROR: ir_bias_mode = aux needs an experimental spectrum."
                  write (*, *) "       Add ir to exp_labels with a file in exp_data_files,"
                  write (*, *) "       or leave ir_bias_mode at acf for a prediction run."
                  stop
               end if
!                    No virial is formed from the coupling, for the reason every
!                    other IR bias here gives: a wrong stress is worse than a
!                    missing one. Under a barostat that stops being a missing
!                    diagnostic and becomes a wrong cell, so it is refused rather
!                    than left to be discovered in the density.
               if (trim(params%barostat) /= "none") then
                  write (*, *) "ERROR: ir_bias_mode = aux does not form a virial, so the cell"
                  write (*, *) "       would relax against an incomplete stress. Run it at fixed"
                  write (*, *) "       volume, or equilibrate the cell first with barostat on and"
                  write (*, *) "       the bias off."
                  call turbogap_abort()
               end if
               call ir_aux_setup(mad_ir_state, state%n_sites, params%md_step, params%ir_stride, &
                                 params%ir_aux_eff_mass, params%ir_aux_damping, &
                                 params%ir_aux_tau, params%ir_aux_gain, &
                                 params%ir_aux_eta_max, params%ir_aux_restart_file, &
                                 ir%ir_aux_ok, ir%ir_aux_resumed, ir%ir_aux_msg)
               if (.not. ir%ir_aux_ok) then
                  write (*, *) "ERROR: ", trim(ir%ir_aux_msg)
                  stop
               end if
               if (comm%rank == 0) then
                  write (*, *) 'MAD IR, envelope-targeted bank:        |'
                  write (*, '(A,I12,A)') '  *) resonator modes:   ', &
                     ir_aux_state%n_modes, '         |'
                  write (*, '(A,F12.4,A)') '  *) fictitious mass:   ', &
                     params%ir_aux_eff_mass, ' amu     |'
                  write (*, '(A,F12.4,A)') '  *) bank resolution:   ', &
                     params%ir_aux_damping, ' cm^-1   |'
                  write (*, '(A,F12.4,A)') '  *) controller time:   ', &
                     ir_aux_state%tau(1), ' fs      |'
                  write (*, '(A,ES12.4,A)') '  *) controller gain:   ', &
                     ir_aux_state%kappa(1), ' 1/fs    |'
                  write (*, '(A,F12.4,A)') '  *) critical gain:     ', &
                     2.d0/ir_aux_state%tau(1), ' 1/fs    |'
                  if (ir%ir_aux_resumed) then
                     write (*, '(A,I12,A)') '  *) resumed, advances: ', &
                        ir_aux_state%n_steps, '         |'
                  else
                     write (*, *) '  *) uncalibrated: the bank waits    |'
                     write (*, *) '     for the ACF ensemble to fill   |'
                     write (*, *) '     and applies no force until it  |'
                     write (*, *) '     does.                          |'
                  end if
                  if (len_trim(ir%ir_aux_msg) > 0) write (*, *) '     ', trim(ir%ir_aux_msg)
                  write (*, *) '.......................................|'
               end if
            end if
         end if
!              Only a biased run needs dmu/dr; see mad_ir_need_dmu. Set every
!              step rather than once, because it costs nothing and there is no
!              earlier point at which params is known to be final.
!
!              Under the extended Lagrangian nothing wants the (3,3,n_atoms)
!              tensor at all: the weight is already known, so the descriptor
!              pass contracts it and leaves a force behind instead. The two
!              are mutually exclusive and mad_ir_need_dmu is what keeps them so.
         mad_ir_need_dmu = params%valid_ir .and. params%exp_forces &
                           .and. .not. mad_ir_xl_active
         mad_ir_collect = (loop%md_istep >= 0) .and. (modulo(loop%md_istep, params%ir_stride) == 0)
         if (mad_ir_collect .and. mad_ir_need_dmu) mad_ir_dmu_dr = 0.d0
         mad_ir_xl_collect = mad_ir_xl_active .and. mad_ir_collect &
                             .and. params%valid_ir .and. params%exp_forces
!              mad_ir_xl_site_w ALREADY HOLDS THE WEIGHT. It was formed at the
!              end of the previous stored frame, from a bank that had just been
!              advanced with that frame's dipole, and it is the exact gradient
!              of the bias energy with respect to those dipoles. That is why it
!              can be contracted inside the descriptor pass instead of leaving
!              a (3,3,n_atoms) tensor behind for a weight that does not exist
!              yet; the price is one stored frame of lag in the POSITIONS the
!              gradient is contracted at, which is far below the bank's own
!              time resolution. It is identically zero while the bank charges.
         if (mad_ir_xl_collect) mad_ir_xl_force = 0.d0
      end if
   end subroutine ir_before_evaluate

   subroutine ir_push_frame(ir, res, state)
      type(ir_run_t), intent(inout) :: ir
      type(results_t), intent(in) :: res
      type(state_t), intent(in) :: state

      if (ir%ir_from_traj) then
         call ir_fft_frames_push(ir_fft_frames, res%dipole, state%frame_time, state%has_frame_time)
      end if
   end subroutine ir_push_frame

!  After the forces: push this frame's dipole, apply whichever bias the run
!  asked for, and write the restart and spectrum files on their schedule.
   subroutine ir_after_forces(ir, params, state, res, loop, comm, md_time, time3, write_condition, time)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(inout) :: params
      type(state_t), intent(inout) :: state
      type(results_t), intent(inout) :: res
      type(loop_t), intent(in) :: loop
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: md_time
      real(dp), intent(in) :: time3
      logical, intent(inout) :: write_condition
      type(times_t), intent(inout) :: time

      if ((params%valid_ir .or. params%do_ir) .and. mad_ir_collect) then
!              Nothing to reduce when no force is being formed from it.
         if (mad_ir_need_dmu) then
            call time_start(time%mpi)
!                 every rank holds a slice, so a plain sum is the whole reduction;
!                 all-reduce so each rank can form the same lambda and bias its own
!                 atoms without a second broadcast of the forces
            call comm_sum_all(comm, mad_ir_dmu_dr, 9*state%n_sites)
            call time_end(time%mpi)
         end if
!              The extended Lagrangian reduces a FORCE rather than a tensor, so
!              a third of the traffic. It is still an all-reduce and not a
!              reduce: the bank itself is replicated and advanced identically on
!              every rank, driven by local_dipoles, which is already broadcast,
!              so every rank must end the step holding the same forces.
         if (mad_ir_xl_collect) then
            call time_start(time%mpi)
            call comm_sum_all(comm, mad_ir_xl_force, 3*state%n_sites)
            call time_end(time%mpi)
         end if
         call time_start(time%ir)
!              The ACF ensemble is filled under BOTH biases. Under the extended
!              Lagrangian it produces no force -- mad_ir_evaluate is never
!              called -- but it costs three doubles a frame and it is the only
!              independent estimate of the same spectrum the bank is claiming,
!              so ir_spectrum.dat stays meaningful as a cross-check.
         call mad_ir_push(mad_ir_state, res%dipole)
!              A prediction run has no experiment to compare against, so it
!              does none of what follows: it only accumulates, and transforms
!              once when the file is written.
         if (mad_ir_xl_active) then
!                 THE RESONATOR BANK (ir_bias_mode = xl).
!
!                 Four things happen here in an order that is not
!                 interchangeable: advance the bank with this frame's dipoles,
!                 evaluate the spectrum and the loss from the state that
!                 produces, add the force the descriptor pass already left, and
!                 form the weight the NEXT frame's pass will contract.
!
!                 The force added here was built with the weight formed at the
!                 end of the PREVIOUS stored frame. That is the one frame of lag
!                 the scheme trades for contracting inside the descriptor pass;
!                 mad_ir_xl.f90's header has why it is the right trade here and
!                 the wrong one for the ACF bias.
            call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                  loop%mc_istep, params%mc_nsteps, &
                                  params%exp_energy_scales_initial(params%ir_idx), &
                                  params%exp_energy_scales_final(params%ir_idx), ir%mad_ir_scale)
            call time_start(time%ir_predict)
!                 Advance BEFORE evaluate: the loss and the gradient both belong
!                 to a bank that already knows this frame's dipole. Evaluating
!                 first would report the previous frame's spectrum and, worse,
!                 would make the weight the gradient of a loss the current
!                 dipole had not yet entered -- which is identically zero, not
!                 merely inaccurate.
            call mad_ir_xl_advance(mad_ir_xl_state, res%local_dipoles(1:3, 1:state%n_sites))
            call mad_ir_xl_evaluate(mad_ir_xl_state, ir%mad_ir_scale, ir%mad_ir_energy)
            call time_end(time%ir_predict)
            res%energies_exp = res%energies_exp + ir%mad_ir_energy/dfloat(state%n_sites)
            exp_dissimilarity = exp_dissimilarity + mad_ir_xl_state%dissim
            exp_dissim_ref = exp_dissim_ref + mad_ir_xl_state%dissim_ref
            if (params%exp_energies) then
               res%energies = res%energies + ir%mad_ir_energy/dfloat(state%n_sites)
               res%energy = sum(res%energies)
            end if
            res%energy_exp = sum(res%energies_exp)
!                 mad_ir_xl_force is whatever the descriptor pass left, which is
!                 identically zero while the bank is charging because the weight
!                 it was contracted with was. There is no readiness test here for
!                 that reason: the gate is in mad_ir_xl_weights, where it can be
!                 applied once instead of in every consumer.
            if (params%exp_forces .and. mad_ir_xl_collect) then
               call time_start(time%ir_forces)
               res%forces(1:3, 1:state%n_sites) = res%forces(1:3, 1:state%n_sites) &
                                                  + mad_ir_xl_force(1:3, 1:state%n_sites)
!                    And the weight the NEXT stored frame's descriptor pass
!                    will contract, from the bank as it now stands.
               call mad_ir_xl_weights(mad_ir_xl_state, mad_ir_xl_site_w)
               call time_end(time%ir_forces)
            end if
            if (mad_ir_xl_ready(mad_ir_xl_state) .and. .not. ir%mad_ir_applied) then
               call get_time(ir%mad_ir_t_now)
               ir%mad_ir_t_first = ir%mad_ir_t_now - time3
               ir%mad_ir_step_first = loop%md_istep
            end if
            ir%mad_ir_applied = mad_ir_xl_ready(mad_ir_xl_state)
            if (.not. ir%mad_ir_applied) ir%mad_ir_energy = 0.d0
         else if (ir_aux_active) then
!                 THE ENVELOPE-TARGETED BANK (ir_bias_mode = aux).
!
!                 The bank is a filter bank whose OWN amplitude is held at the
!                 experimental one by feedback friction, and the atoms feel it
!                 only through the dipole coupling. So unlike the other three
!                 modes there is no loss to differentiate here: the force is
!                 escale * sum_k g_k J_i^T X_k, whose generator is the coupling
!                 energy ir_aux_evaluate reports. ir_auxiliary_dynamics.f90's
!                 header has why the restraint cannot be a potential instead.
!
!                 Calibration first, and once. It needs S_MM(w_k) from the
!                 autocorrelation, so it cannot happen at setup time; this is
!                 the first step on which the ensemble is full. Until then the
!                 branch does nothing at all and the run is plain MD, which is
!                 exactly the unbiased trajectory the calibration wants.
            if (.not. ir_aux_calibrated(ir_aux_state) .and. mad_ir_ready(mad_ir_state)) then
               call time_start(time%ir_predict)
!                    The bound is tested against the WORST case the ramp will
!                    reach, and at the LOWER of the two temperatures, since a
!                    colder signal is a softer one and softens the threshold.
               ir%ir_aux_escale_top = max(params%exp_energy_scales_initial(params%ir_idx), &
                                          params%exp_energy_scales_final(params%ir_idx))
               ir%ir_aux_temp = min(params%t_beg, params%t_end)
               if (ir%ir_aux_temp <= 0.d0) ir%ir_aux_temp = max(params%t_beg, params%t_end)
               call ir_aux_calibrate(ir_aux_state, mad_ir_state, ir%ir_aux_temp, &
                                     ir%ir_aux_escale_top, ir%ir_aux_ok, ir%ir_aux_msg)
               call time_end(time%ir_predict)
               if (comm%rank == 0) then
                  if (ir%ir_aux_ok) then
                     write (*, *) 'MAD IR: ', trim(ir%ir_aux_msg)
                  else
                     write (*, *) 'WARNING: ', trim(ir%ir_aux_msg)
                  end if
               end if
               if (.not. ir%ir_aux_ok) then
                  write (*, *) "ERROR: ir_bias_mode = aux could not calibrate the bank."
                  stop
               end if
!                    THE STABILITY BOUND. The coupling -g_k X_k.s is bilinear and
!                    therefore unbounded below; it is held only by the resonator
!                    spring and by the stiffness of the signal against the physical
!                    potential. Past Lambda = 1 the combined quadratic form is
!                    indefinite and the pair runs away exponentially -- measured,
!                    for 64 H2O at the defaults, as 1e8 K within ten steps.
!
!                    The linear-response calibration knows nothing about this: it
!                    matches amplitudes, and whether the resulting coupling is
!                    below threshold is a separate question it never asks. So the
!                    check is here, it is made before the first biased step, and
!                    it is fatal -- a run above threshold does not produce a worse
!                    answer, it produces no answer at all.
               if (ir_aux_state%stab >= 1.d0) then
                  if (comm%rank == 0) then
                     write (*, *) "ERROR: ir_bias_mode = aux is above its stability threshold."
                     write (*, '(A,ES12.4)') "        stability number Lambda = ", ir_aux_state%stab
                     write (*, *) "        Lambda must be below 1; the bilinear coupling runs away above it."
                     write (*, '(A,ES12.4)') "        largest usable exp_energy_scales = ", &
                        ir_aux_escale_max(ir_aux_state)
                     write (*, '(A,ES12.4)') "        this run asked for               = ", ir%ir_aux_escale_top
                     write (*, *) "        Lower exp_energy_scales, or raise ir_aux_damping"
                     write (*, *) "        (which lowers every g_k), and try again."
                  end if
                  call turbogap_abort()
               else if (ir_aux_state%stab >= 0.5d0 .and. comm%rank == 0) then
                  write (*, '(A,ES10.2,A)') " WARNING: ir_aux stability number ", &
                     ir_aux_state%stab, " is above 0.5; the bias is close to runaway."
               end if
               call get_time(ir%mad_ir_t_now)
               ir%mad_ir_t_first = ir%mad_ir_t_now - time3
               ir%mad_ir_step_first = loop%md_istep
            end if
            if (ir_aux_calibrated(ir_aux_state)) then
               call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                     loop%mc_istep, params%mc_nsteps, &
                                     params%exp_energy_scales_initial(params%ir_idx), &
                                     params%exp_energy_scales_final(params%ir_idx), ir%mad_ir_scale)
               call time_start(time%ir_predict)
!                    Advance before evaluate, for the same reason the extended
!                    Lagrangian does: the energy and the force both belong to a
!                    bank that already knows this frame's dipole.
               call ir_aux_advance(ir_aux_state, res%dipole)
               call ir_aux_evaluate(ir_aux_state, ir%mad_ir_scale, res%dipole, ir%mad_ir_energy)
               call time_end(time%ir_predict)
               res%energies_exp = res%energies_exp + ir%mad_ir_energy/dfloat(state%n_sites)
!                    THE FIGURE OF MERIT IS THE ACF SPECTRUM, NOT THE BANK'S.
!
!                    The controller drives R_k onto R_target by construction, so
!                    the bank's own mismatch falls to zero whether or not the ATOMS
!                    have learned anything -- it measures the controller, not the
!                    physics. The independent estimate is the one mad_ir forms from
!                    the stored dipoles, which no part of this mode steers, and that
!                    is what goes in thermo.log.
!
!                    Called with a zero energy scale: mad_ir_evaluate sets dissim
!                    and dissim_ref regardless of it, so this buys the honest number
!                    and contributes no energy and no force.
               call time_start(time%ir_predict)
               call mad_ir_evaluate(mad_ir_state, 0.d0, ir%ir_aux_acf_energy, ir%mad_ir_lambda)
               call time_end(time%ir_predict)
               exp_dissimilarity = exp_dissimilarity + mad_ir_state%dissim
               exp_dissim_ref = exp_dissim_ref + mad_ir_state%dissim_ref
               if (params%exp_energies) then
                  res%energies = res%energies + ir%mad_ir_energy/dfloat(state%n_sites)
                  res%energy = sum(res%energies)
               end if
               res%energy_exp = sum(res%energies_exp)
               if (params%exp_forces) then
                  call time_start(time%ir_forces)
                  call ir_aux_forces(ir_aux_state, ir%mad_ir_scale, mad_ir_dmu_dr, res%forces, &
                                     state%velocities(1:3, 1:state%n_sites), ir%ir_aux_power)
                  call time_end(time%ir_forces)
!                       Integrated with the stored-frame interval, since that is
!                       how often the force is refreshed.
                  ir%ir_aux_work = ir%ir_aux_work &
                                   + ir%ir_aux_power*params%md_step*dfloat(params%ir_stride)
               end if
               ir%mad_ir_applied = .true.
!                    ---- the thermal-fidelity check -------------------------
!                    Over a window of stored frames, how much energy did the
!                    controller put in, and what standing temperature offset
!                    does the thermostat therefore have to hold against?
!                    50 fs of biased dynamics is enough to average the pump
!                    rate over many resonator periods (the fastest fitted band
!                    is ~8 fs) while still reporting inside a short run.
               if (md_time - ir%ir_aux_pump_time > 50.d0) then
                  if (ir%ir_aux_pump_time > 0.d0 .and. params%tau_t > 0.d0) then
!                          The work the BIAS FORCE did on the atoms, which is
!                          the channel that actually heats: the controller's own
!                          injection into the bank is a different and, for a bank
!                          far off target, much smaller number.
                     ir%ir_aux_dT = 2.d0*(ir%ir_aux_work - ir%ir_aux_pump_prev) &
                                    /(md_time - ir%ir_aux_pump_time)*params%tau_t &
                                    /(3.d0*dfloat(state%n_sites)*8.6173303d-5)
                     if (comm%rank == 0 .and. .not. ir%ir_aux_warned_hot .and. &
                         dabs(ir%ir_aux_dT) > 0.1d0*max(1.d0, params%t_beg)) then
                        write (*, '(A,F10.1,A)') " WARNING: ir_aux is pumping hard enough for a standing", &
                           ir%ir_aux_dT, " K offset."
                        write (*, *) "          The bias is below its stability bound but above what the"
                        write (*, *) "          thermostat can absorb quietly. Lower exp_energy_scales."
                        ir%ir_aux_warned_hot = .true.
                     end if
                  end if
                  ir%ir_aux_pump_prev = ir%ir_aux_work
                  ir%ir_aux_pump_time = md_time
               end if
            else
               ir%mad_ir_energy = 0.d0
               ir%mad_ir_applied = .false.
            end if
         else if (params%valid_ir .and. trim(params%ir_bias_mode) == "fft" &
                  .and. mad_ir_ready(mad_ir_state)) then
!                 THE FFT ESTIMATOR (ir_bias_mode = fft).
!
!                 The ensemble is mad_ir's rolling buffer -- mad_ir_push filled
!                 it above, as it does under every mode -- so all that differs
!                 from the ACF branch below is which routine turns that buffer
!                 into a loss and a lambda. The buffer is circular and
!                 ir_fft_loss wants a chronological array, so it is unrolled
!                 first: 3*n_window doubles copied per stored frame, against a
!                 transform that is already O(n_lag * n_bins).
!
!                 lambda is dL/dmu of the NEWEST configuration, the same
!                 quantity mad_ir_evaluate returns, so mad_ir_forces contracts
!                 it with the same dmu/dr and the force path is unchanged. That
!                 is the whole reason this fits in as a branch rather than as a
!                 second bias: the two estimators disagree about the spectrum
!                 and agree exactly about what a bias on a dipole is.
            call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                  loop%mc_istep, params%mc_nsteps, &
                                  params%exp_energy_scales_initial(params%ir_idx), &
                                  params%exp_energy_scales_final(params%ir_idx), ir%mad_ir_scale)
            call time_start(time%ir_predict)
            if (allocated(ir%ir_fft_mu_chron)) then
               if (size(ir%ir_fft_mu_chron, 2) /= mad_ir_state%n_stored) &
                  deallocate (ir%ir_fft_mu_chron)
            end if
            if (.not. allocated(ir%ir_fft_mu_chron)) &
               allocate (ir%ir_fft_mu_chron(1:3, 1:mad_ir_state%n_stored))
            if (.not. allocated(ir%ir_fft_I_fit)) &
               allocate (ir%ir_fft_I_fit(1:mad_ir_state%n_freq))
            call ir_fft_md_unroll(mad_ir_state%mu_hist, mad_ir_state%n_window, &
                                  mad_ir_state%n_stored, mad_ir_state%head, &
                                  ir%ir_fft_mu_chron, ir%ir_fft_n_chron)
            ir%ir_fft_cfg%dt_fs = mad_ir_state%dt
            call ir_fft_loss(ir%ir_fft_mu_chron, ir%ir_fft_n_chron, ir%ir_fft_cfg, &
                             mad_ir_state%nu, mad_ir_state%I_exp, mad_ir_state%wgt, &
                             mad_ir_state%n_freq, params%ir_match_scale, &
                             params%ir_match_offset, ir%mad_ir_scale, &
                             ir%mad_ir_energy, ir%mad_ir_lambda, ir%ir_fft_I_fit, &
                             ir%ir_fft_scale_fit, ir%ir_fft_offset_fit, &
                             ir%ir_fft_dissim, ir%ir_fft_dissim_ref, ir%ir_fft_ok, ir%ir_fft_msg)
            call time_end(time%ir_predict)
            if (.not. ir%ir_fft_ok) then
               write (*, *) "ERROR: ", trim(ir%ir_fft_msg)
               stop
            end if
            res%energies_exp = res%energies_exp + ir%mad_ir_energy/dfloat(state%n_sites)
            exp_dissimilarity = exp_dissimilarity + ir%ir_fft_dissim
            exp_dissim_ref = exp_dissim_ref + ir%ir_fft_dissim_ref
            if (params%exp_energies) then
               res%energies = res%energies + ir%mad_ir_energy/dfloat(state%n_sites)
               res%energy = sum(res%energies)
            end if
            res%energy_exp = sum(res%energies_exp)
            if (params%exp_forces) then
               call time_start(time%ir_forces)
               call mad_ir_forces(ir%mad_ir_lambda, mad_ir_dmu_dr, res%forces)
               call time_end(time%ir_forces)
            end if
            if (.not. ir%mad_ir_applied) then
               call get_time(ir%mad_ir_t_now)
               ir%mad_ir_t_first = ir%mad_ir_t_now - time3
               ir%mad_ir_step_first = loop%md_istep
            end if
            ir%mad_ir_applied = .true.
         else if (params%valid_ir .and. mad_ir_ready(mad_ir_state)) then
!                 The weight is exp_energy_scales, ramped over the run exactly
!                 as every other MAD observable's is.
            call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                  loop%mc_istep, params%mc_nsteps, &
                                  params%exp_energy_scales_initial(params%ir_idx), &
                                  params%exp_energy_scales_final(params%ir_idx), ir%mad_ir_scale)
            call time_start(time%ir_predict)
            call mad_ir_evaluate(mad_ir_state, ir%mad_ir_scale, ir%mad_ir_energy, ir%mad_ir_lambda)
            call time_end(time%ir_predict)
!                 The mismatch is an energy like the others, spread over the
!                 sites; the force is its gradient, and only if exp_forces.
            res%energies_exp = res%energies_exp + ir%mad_ir_energy/dfloat(state%n_sites)
!                 IR does not go through get_exp_energies -- it has its own
!                 fitted scale and offset -- so it contributes to the shared
!                 dissimilarity accumulator here instead. Same definition:
!                 the residual sum of squares with no energy scale on it.
            exp_dissimilarity = exp_dissimilarity + mad_ir_state%dissim
            exp_dissim_ref = exp_dissim_ref + mad_ir_state%dissim_ref
!                 The sum of energies_exp into energies happened above, before
!                 the dipole of this configuration existed. Folding the IR term
!                 in there is not possible -- it needs the forces pass -- so it
!                 is folded in here instead. Without this the mismatch never
!                 reached the reported total energy at all: energies_exp is
!                 zeroed at the top of the next step.
            if (params%exp_energies) then
               res%energies = res%energies + ir%mad_ir_energy/dfloat(state%n_sites)
               res%energy = sum(res%energies)
            end if
            res%energy_exp = sum(res%energies_exp)
            if (params%exp_forces) then
               call time_start(time%ir_forces)
               call mad_ir_forces(ir%mad_ir_lambda, mad_ir_dmu_dr, res%forces)
               call time_end(time%ir_forces)
            end if
            if (.not. ir%mad_ir_applied) then
!                    The first spectrum of the run: how long the ensemble took
!                    to fill, measured from the same origin as the total.
               call get_time(ir%mad_ir_t_now)
               ir%mad_ir_t_first = ir%mad_ir_t_now - time3
               ir%mad_ir_step_first = loop%md_istep
            end if
            ir%mad_ir_applied = .true.
         else
            ir%mad_ir_energy = 0.d0
            ir%mad_ir_applied = .false.
         end if
!              Persist the ensemble alongside the trajectory. Losing it costs
!              n_window samples of unbiased dynamics on the next restart.
         call time_start(time%ir_io)
         if (comm%rank == 0 .and. params%write_xyz > 0 .and. params%valid_ir) then
            if (modulo(loop%md_istep, params%write_xyz) == 0 .or. loop%md_istep == params%md_nsteps) then
               if (trim(params%ir_restart_file) /= "none") then
                  call mad_ir_save(mad_ir_state, params%ir_restart_file, ir%mad_ir_ok, ir%mad_ir_msg)
                  if (.not. ir%mad_ir_ok) write (*, *) "WARNING: ", trim(ir%mad_ir_msg)
               end if
!                    The bank, in its own file. Losing it costs the warm-up
!                    again, and it is the larger of the two by orders of
!                    magnitude, which is why it is a separate write that can be
!                    turned off on its own.
               if (mad_ir_xl_active .and. trim(params%ir_xl_restart_file) /= "none") then
                  call mad_ir_xl_save(mad_ir_xl_state, params%ir_xl_restart_file, &
                                      ir%mad_ir_xl_ok, ir%mad_ir_xl_msg)
                  if (.not. ir%mad_ir_xl_ok) write (*, *) "WARNING: ", trim(ir%mad_ir_xl_msg)
               end if
!                    And the envelope-targeted bank. X, P and eta are state
!                    in the same sense the velocities are, and eta especially:
!                    it is an integrator, so dropping it throws away everything
!                    the controller had learned about the mismatch.
               if (ir_aux_active .and. ir_aux_calibrated(ir_aux_state) .and. &
                   trim(params%ir_aux_restart_file) /= "none") then
                  call ir_aux_save(ir_aux_state, params%ir_aux_restart_file, &
                                   ir%ir_aux_ok, ir%ir_aux_msg)
                  if (.not. ir%ir_aux_ok) write (*, *) "WARNING: ", trim(ir%ir_aux_msg)
               end if
            end if
         end if
         call time_end(time%ir_io)
!              The spectrum, on the same schedule every other observable's
!              prediction is written on. ir_spectrum.dat is the current one
!              with its sampling limits in the header; ir_prediction.dat
!              accumulates one block per write and so covers the trajectory.
!
!              Two conditions, not one, because the two modes become ready at
!              different times: a fit has a spectrum once the rolling window
!              is full, a prediction has one as soon as there are two frames
!              to correlate -- and the prediction's last write, at the final
!              step, is the one the run is for.
!              get_write_condition takes modulo(md_istep, write_xyz), so it
!              cannot be asked anything when write_xyz is zero -- which is its
!              default, and the natural setting for a prediction run that
!              wants one spectrum and no trajectory.
         if (params%write_xyz > 0) then
            call get_write_condition(params%do_mc, params%do_md, &
                                     loop%mc_istep, loop%md_istep, params%write_xyz, write_condition)
         else
            write_condition = .false.
         end if
         if (params%do_ir) then
            ir%mad_ir_have_spectrum = mad_ir_state%n_stored > 1
!                 The last COLLECTED step, not the last step: with an
!                 ir_stride that does not divide md_nsteps the two differ, and
!                 the final frame is the one the whole run was for.
            write_condition = write_condition .or. &
                              (loop%md_istep > params%md_nsteps - params%ir_stride)
         else
            ir%mad_ir_have_spectrum = ir%mad_ir_applied
         end if
         if (comm%rank == 0 .and. params%write_ir .and. ir%mad_ir_have_spectrum &
             .and. write_condition) then
!                 A fit already transformed this step's ensemble on its way to
!                 the loss; a prediction has not, and this is the only place
!                 that asks for it. Timed as predict rather than io, and
!                 outside the io region, so the two buckets stay disjoint and
!                 "i/o" means the filesystem.
            if (params%do_ir) then
               call time_start(time%ir_predict)
               call mad_ir_spectrum(mad_ir_state)
               call time_end(time%ir_predict)
            end if
            call time_start(time%ir_io)
!                 THE FFT ESTIMATOR'S OWN SPECTRUM. ir_spectrum.dat below is
!                 always the block ACF -- that is what it has always meant and
!                 changing it would silently rewrite the meaning of every
!                 existing plotting script -- so under ir_bias_mode = fft the
!                 quantity actually being biased has to be written somewhere
!                 else, or it cannot be looked at at all. Same argument as
!                 ir_xl_spectrum.dat.
!
!                 The transform is redone here rather than cached from the
!                 bias: ir_fft_loss frees its intermediates, and this happens
!                 on write_xyz steps rather than every step, so recomputing is
!                 cheaper than keeping a copy alive across the whole run.
            if (ir%ir_fft_active .and. .not. ir%ir_from_traj &
                .and. mad_ir_state%n_stored > 1) then
               if (allocated(ir%ir_fft_mu_chron)) then
                  if (size(ir%ir_fft_mu_chron, 2) /= mad_ir_state%n_stored) &
                     deallocate (ir%ir_fft_mu_chron)
               end if
               if (.not. allocated(ir%ir_fft_mu_chron)) &
                  allocate (ir%ir_fft_mu_chron(1:3, 1:mad_ir_state%n_stored))
               call ir_fft_md_unroll(mad_ir_state%mu_hist, mad_ir_state%n_window, &
                                     mad_ir_state%n_stored, mad_ir_state%head, &
                                     ir%ir_fft_mu_chron, ir%ir_fft_n_chron)
               ir%ir_fft_cfg%dt_fs = mad_ir_state%dt
               call ir_fft_spectrum(ir%ir_fft_mu_chron, ir%ir_fft_n_chron, ir%ir_fft_cfg, &
                                    ir%ir_fft_res, ir%ir_fft_ok, ir%ir_fft_msg)
               if (ir%ir_fft_ok) then
                  call ir_fft_write_spectrum(ir%ir_fft_res, ir%ir_fft_cfg, &
                                             "ir_fft_spectrum.dat", &
                                             mad_ir_state%nu, mad_ir_state%I_exp, &
                                             mad_ir_state%n_freq, params%valid_ir, &
                                             ir%ir_fft_scale_fit, ir%ir_fft_offset_fit, &
                                             ir%ir_fft_dissim, ir%ir_fft_dissim_ref, &
                                             "from the rolling MD ensemble")
                  call ir_fft_free(ir%ir_fft_res)
               else
                  write (*, *) "WARNING: ir_fft_spectrum.dat not written: ", &
                     trim(ir%ir_fft_msg)
               end if
            end if
            call mad_ir_write_spectrum(mad_ir_state, "ir_spectrum.dat", &
                                       params%valid_ir, loop%md_istep, params%md_step)
            call mad_ir_append_spectrum(mad_ir_state, "ir_prediction.dat", &
                                        .not. ir%mad_ir_wrote_prediction, loop%md_istep, &
                                        dfloat(loop%md_istep)*params%md_step)
            if (.not. ir%mad_ir_wrote_prediction) then
               if (params%valid_ir) then
!                       "<label>_exp.dat" is the convention, but for label
!                       "ir" that is a name a user may well have given the
!                       file in exp_data_files -- tests/mad_ir/md_run.sh does
!                       exactly that -- and writing it would destroy the input
!                       midway through the run. Reading the same path back on
!                       the next restart would then fit the prediction to
!                       itself.
                  if (trim(params%exp_data(params%ir_idx)%file_data) == "ir_exp.dat") then
                     write (*, *) "WARNING: not writing ir_exp.dat; it is the file named"
                     write (*, *) "         in exp_data_files and would be overwritten."
                  else
                     call mad_ir_write_exp_spectrum(mad_ir_state, "ir_exp.dat")
                  end if
               end if
               ir%mad_ir_wrote_prediction = .true.
            end if
!                 The bank's own spectrum, beside the ACF one. The two are
!                 independent estimates of the same quantity from the same
!                 trajectory -- one a windowed transform of a stored ensemble,
!                 the other the standing amplitude of a filter bank -- and
!                 whether they agree is the single most useful check there is
!                 that the bank is measuring what it claims to.
            if (mad_ir_xl_active) then
               call mad_ir_xl_write_spectrum(mad_ir_xl_state, "ir_xl_spectrum.dat", &
                                             params%valid_ir)
            end if
            if (ir_aux_active .and. ir_aux_calibrated(ir_aux_state)) then
               call ir_aux_write_spectrum(ir_aux_state, "ir_aux_spectrum.dat", &
                                          params%valid_ir)
            end if
            call time_end(time%ir_io)
         end if
         call time_end(time%ir)
      end if
   end subroutine ir_after_forces

   subroutine ir_step_end(ir, params, loop)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(in) :: params
      type(loop_t), intent(in) :: loop

!     Close the per-step stopwatch and charge it to whichever side of the
!     boundary this step fell on. mad_ir_applied was set for THIS step in the
!     force block above, so the two accumulators separate exactly at the step
!     the ensemble filled.
      if (params%valid_ir .and. params%do_md .and. loop%md_istep >= 0) then
         call get_time(ir%mad_ir_t_now)
         if (ir%mad_ir_applied) then
            ir%mad_ir_t_post = ir%mad_ir_t_post + (ir%mad_ir_t_now - ir%mad_ir_step_beg)
            ir%mad_ir_n_post = ir%mad_ir_n_post + 1
         else
            ir%mad_ir_t_pre = ir%mad_ir_t_pre + (ir%mad_ir_t_now - ir%mad_ir_step_beg)
            ir%mad_ir_n_pre = ir%mad_ir_n_pre + 1
         end if
      end if
   end subroutine ir_step_end

!  IR from a trajectory: the transform, once the whole file has been read.
   subroutine ir_finish(ir, params, comm, time)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(in) :: params
      type(comm_t), intent(in) :: comm
      type(times_t), intent(inout) :: time

!
!  IR PREDICTION FROM A TRAJECTORY: the transform, now that the file has been
!  read to the end and its length is known.
!
!  This is the whole of the do_ir predict path. Everything before it only
!  accumulated (time, dipole) pairs; the resolution, the grid and the sizing
!  all follow from the number of pairs, which is why none of it could happen
!  earlier.
!
!  Rank 0 writes, but every rank ran the transform on an identical buffer, so
!  there is nothing to reduce and no rank can be holding a different answer.
      if (ir%ir_from_traj) then
         call time_start(time%ir_predict)
         if (params%valid_ir) then
!        There is an experiment: restrict it to [ir_nu_min, ir_nu_max] and
!        weight it exactly as the MAD bias would, so the mismatch printed here
!        is the same number a biased run would be minimising.
            call mad_ir_select_range(params%exp_data(params%ir_idx)%data(1, :), &
                                     params%exp_data(params%ir_idx)%data(2, :), &
                                     params%ir_nu_min, params%ir_nu_max, &
                                     params%ir_weight_by_spacing, &
                                     ir%ir_nu_exp, ir%ir_I_exp, ir%ir_wgt_exp, ir%mad_ir_ok, ir%mad_ir_msg)
            if (.not. ir%mad_ir_ok) then
               write (*, *) "ERROR: ", trim(ir%mad_ir_msg)
               stop
            end if
         else
            allocate (ir%ir_nu_exp(1:1), ir%ir_I_exp(1:1), ir%ir_wgt_exp(1:1))
            ir%ir_nu_exp = 0.d0; ir%ir_I_exp = 0.d0; ir%ir_wgt_exp = 1.d0
         end if

         call ir_fft_frames_finish(ir_fft_frames, ir%ir_fft_cfg, params%ir_frame_dt, &
                                   params%ir_frame_dt_tol, ir%ir_nu_exp, ir%ir_I_exp, &
                                   ir%ir_wgt_exp, size(ir%ir_nu_exp), params%valid_ir, &
                                   params%ir_match_scale, params%ir_match_offset, &
                                   "ir_fft_spectrum.dat", "ir_fft_dipoles.dat", &
                                   comm%rank == 0, params%ir_fft_write_dipoles, &
                                   ir%ir_fft_res, ir%ir_fft_dt_used, ir%ir_fft_scale_fit, &
                                   ir%ir_fft_offset_fit, ir%ir_fft_dissim, ir%ir_fft_dissim_ref, &
                                   ir%ir_fft_ok, ir%ir_fft_msg)
         call time_end(time%ir_predict)

         if (.not. ir%ir_fft_ok) then
            if (comm%rank == 0) then
               write (*, *) ""
               write (*, *) "ERROR: ", trim(ir%ir_fft_msg)
            end if
!        A NONZERO exit, unlike the bare `stop` used elsewhere in this file.
!        This path is driven by scripts -- post-processing a directory of
!        trajectories is the obvious use -- and the whole design here is
!        organised against failing silently. Exiting 0 with no spectrum
!        written is precisely that failure wearing a message.
            stop 1
         end if

         if (comm%rank == 0) then
            write (*, *) '                                       |'
            write (*, *) 'IR spectrum from the trajectory:       |'
            write (*, '(A,I12,A)') '  *) frames read:       ', ir_fft_frames%n, '         |'
            write (*, '(A,F12.4,A)') '  *) frame interval:    ', ir%ir_fft_dt_used, ' fs      |'
            write (*, '(A,I12,A)') '  *) lags kept:         ', ir%ir_fft_res%n_lag, '         |'
            write (*, '(A,F12.4,A)') '  *) resolution:        ', ir%ir_fft_res%resolution, ' cm^-1   |'
            write (*, '(A,F12.4,A)') '  *) bin spacing:       ', ir%ir_fft_res%d_nu, ' cm^-1   |'
            write (*, '(A,F12.1,A)') '  *) Nyquist:           ', ir%ir_fft_res%nyquist, ' cm^-1   |'
            write (*, '(A,I12,A)') '  *) bins written:      ', ir%ir_fft_res%n_freq, '         |'
            if (params%valid_ir .and. ir%ir_fft_dissim_ref > 0.d0) then
               write (*, '(A,F12.6,A)') '  *) rel. mismatch:     ', &
                  dsqrt(ir%ir_fft_dissim/ir%ir_fft_dissim_ref), '         |'
            end if
            if (len_trim(ir%ir_fft_msg) > 0) then
               write (*, *) '  *) ', trim(ir%ir_fft_msg)
            end if
            write (*, *) '                                       |'
         end if

         call ir_fft_free(ir%ir_fft_res)
         call ir_fft_frames_reset(ir_fft_frames)
         deallocate (ir%ir_nu_exp, ir%ir_I_exp, ir%ir_wgt_exp)
      end if
   end subroutine ir_finish

!  The IR lines of the timing report, on rank 0.
   subroutine ir_report(ir, params, time)
      type(ir_run_t), intent(inout) :: ir
      type(input_parameters), intent(in) :: params
      type(times_t), intent(in) :: time

!       The MAD IR bias, and what it costs.
!
!       Three separate numbers, because they answer three separate questions
!       and quoting one for another is how a bias gets blamed for a cost it
!       does not carry:
!
!         - the buckets: where the time inside the bias goes. predict is the
!           autocorrelation and the cosine transform, and it grows with
!           n_lag*n_freq, not with the number of atoms. forces is the
!           contraction with dmu/dr and grows with n_atoms.
!
!         - time to the first spectrum: the ensemble has to fill before any
!           spectrum exists, so this is n_window*ir_stride steps of ordinary
!           MD and is the latency before the fit can begin at all.
!
!         - the per-step rates either side of that point. Note what does NOT
!           change across it: the descriptor second derivatives that produce
!           dmu/dr are built on every collected step from step zero, so their
!           cost is already inside the "no bias" rate. The ratio below is
!           therefore the cost of the spectrum and its gradient alone, and the
!           full price of asking for IR at all is larger -- compare the soap
!           bucket here against a run without the observable.
      if (params%valid_ir .or. params%do_ir) then
         if (params%valid_ir) then
            write (*, '(A,F13.3,A)') ' *  MAD IR bias  :', time%ir(3), ' seconds |'
         else
            write (*, '(A,F13.3,A)') ' *  IR prediction:', time%ir(3), ' seconds |'
         end if
         write (*, '(A,F13.3,A)') '     -    predict:', time%ir_predict(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -     forces:', time%ir_forces(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -        i/o:', time%ir_io(3), ' seconds |'
!          Only a bias has a moment at which the spectrum arrives and the
!          cost changes: a prediction accumulates from step zero and
!          transforms at the end, so there is no before and after to compare.
         if (params%valid_ir) then
            if (ir%mad_ir_step_first >= 0) then
               write (*, '(A,I13,A)') '   1st spectrum @:', ir%mad_ir_step_first, ' step    |'
               write (*, '(A,F13.3,A)') '   1st spectrum @:', ir%mad_ir_t_first, ' seconds |'
            else
               write (*, *) '   no spectrum: ensemble never filled  |'
            end if
         end if
         if (ir%mad_ir_n_pre > 0 .and. ir%mad_ir_n_post > 0) then
            ir%mad_ir_rate_pre = ir%mad_ir_t_pre/dfloat(ir%mad_ir_n_pre)
            ir%mad_ir_rate_post = ir%mad_ir_t_post/dfloat(ir%mad_ir_n_post)
            write (*, '(A,F13.5,A)') '  s/step unbiased:', ir%mad_ir_rate_pre, ' seconds |'
            write (*, '(A,F13.5,A)') '  s/step   biased:', ir%mad_ir_rate_post, ' seconds |'
            if (ir%mad_ir_rate_pre > 0.d0) then
               write (*, '(A,F13.3,A)') '  bias slowdown  :', &
                  ir%mad_ir_rate_post/ir%mad_ir_rate_pre, ' x       |'
            end if
         end if
!          HOW HARD THE BIAS IS PULLING. The weight is dU/dm of the bias
!          energy, so its RMS is the size of the force per unit dipole
!          gradient; it is the number to watch when deciding whether
!          exp_energy_scales is doing anything or doing too much. The fitted
!          scale is beside it because the two move together: a scale that
!          drifts by orders of magnitude means the bank and the experiment are
!          not on comparable footings and the weight is not interpretable.
         if (mad_ir_xl_active) then
            write (*, *) '                                       |'
            write (*, *) ' *  XL resonators:                      |'
            write (*, '(A,I13,A)') '     -   advances:', mad_ir_xl_state%n_steps, '         |'
            write (*, '(A,ES13.4,A)') '     - rms weight:', mad_ir_xl_state%w_rms, '         |'
            write (*, '(A,ES13.4,A)') '     - fit scale :', mad_ir_xl_state%scale, '         |'
         end if
      end if
   end subroutine ir_report

end module turbogap_ir
