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
   use ir_fft, only: ir_fft_config_type, ir_fft_result_type

   implicit none

   private

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

end module turbogap_ir
