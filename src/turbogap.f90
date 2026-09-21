! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap.f90, is copyright (c) 2019-2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X   Uttiyoarnab saha
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

program turbogap

   use kinds

   use timing
   use neighbors
   use soap_turbo_desc
   use gap
   use read_files
   use md
   use adaptive_time                        ! for adaptive time simulation (TurboGAP will use these five modules for radiation cascades)
   use electronic_stopping                ! for electronic stopping correction in radiation cascades
   use eph_fdm                                ! for T - dependent parameters - elec. stop. - eph model
   use eph_beta                                ! for the atomic electronic densities  - elec. stop. - eph model
   use eph_electronic_stopping                ! for electronic stopping based in radiation cascades on the eph model
   use mc
   use gap_interface
   use types
   use vdw
   use electrostatics, only: compute_coulomb_direct, compute_coulomb_dsf, compute_coulomb_lamichhane
   use turbogap_setup
   use turbogap_structure, only: state_t
   use turbogap_domain, only: domain_t, neighbors_t, domain_sync_state
   use turbogap_results, only: results_t
   use turbogap_loop, only: loop_t
   use turbogap_exp
   use turbogap_md
   use ipi_driver, only: ipi_driver_open, ipi_driver_exchange, ipi_driver_close
   use gap_backend
   use gpu_context
   use turbogap_vdw
   use turbogap_estat
   use exp_utils
   use exp_interface
   use soap_turbo_functions
   use mad_ir
   use mad_ir_xl
   use ir_auxiliary_dynamics, only: ir_aux_state, ir_aux_active, ir_aux_setup, &
                                    ir_aux_advance, ir_aux_evaluate, ir_aux_forces, &
                                    ir_aux_calibrate, ir_aux_calibrated, ir_aux_save, &
                                    ir_aux_write_spectrum, ir_aux_bank_energy, &
                                    ir_aux_energy_pumped, ir_aux_stability, &
                                    ir_aux_escale_max
   use ir_fft
   use ir_fft_io
   use turbogap_comm, only: comm_t, comm_init, comm_finalize, comm_bcast, comm_sum_to_root, &
                            comm_sum_all, comm_allgather, comm_with_mpi
   use bussi
   use xyz_module
   use keyword_help
#ifdef _GPU
   use F_B_C
   use iso_c_binding
#endif

   implicit none

   ! Variable definitions
   real(dp), allocatable :: soap(:, :)
   real(dp), allocatable :: soap_cart_der(:, :, :)
   real(dp) :: v_uc_prev
   real(dp) :: v_a_uc
   real(dp) :: v_a_uc_prev
   real(dp) :: eVperA3tobar = 1602176.6208d0
   real(dp) :: ranf
   real(dp) :: disp(1:3)
   real(dp) :: d_disp
   real(dp) :: p_accept
   real(dp) :: virial_prev(1:3, 1:3)
   real(dp), allocatable :: masses_types(:)

!  MAD IR bias. lambda is dL/dmu of the newest configuration; mad_ir_applied
!  says whether the ensemble was full enough for a force to have been added.
   real(dp) :: mad_ir_lambda(1:3) = 0.d0
   real(dp) :: mad_ir_energy = 0.d0
   real(dp) :: mad_ir_scale = 0.d0
   logical :: mad_ir_applied = .false.
   logical :: mad_ir_ok, mad_ir_resumed
!  The envelope-targeted bank, ir_bias_mode = aux. Calibrated once, when the
!  ACF ensemble first fills, so the flag says whether that has happened yet.
   logical :: ir_aux_ok, ir_aux_resumed
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
   logical :: mad_ir_xl_ok = .false., mad_ir_xl_resumed = .false.
   character(len=512) :: mad_ir_xl_msg
!  What the bias costs. The ensemble fills partway into the run, so the same
!  run measures both sides of the question: mad_ir_t_pre accumulates the
!  wall-clock of the steps before the first spectrum and mad_ir_t_post that of
!  the steps after, and their per-step means are directly comparable because
!  nothing else about the step changes at that boundary.
   real(dp) :: mad_ir_t_first = -1.d0
   integer :: mad_ir_step_first = -1
   real(dp) :: mad_ir_t_pre = 0.d0, mad_ir_t_post = 0.d0
   integer :: mad_ir_n_pre = 0, mad_ir_n_post = 0
   real(dp) :: mad_ir_step_beg = 0.d0, mad_ir_t_now = 0.d0
   real(dp) :: mad_ir_rate_pre, mad_ir_rate_post
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
   real(dp) :: ir_fft_scale_fit = 1.d0, ir_fft_offset_fit = 0.d0
   real(dp) :: ir_fft_dissim = 0.d0, ir_fft_dissim_ref = 0.d0
   real(dp), allocatable :: ir_fft_mu_chron(:, :)
   real(dp), allocatable :: ir_nu_exp(:), ir_I_exp(:), ir_wgt_exp(:)
   real(dp), allocatable :: ir_fft_I_fit(:)
   integer :: ir_fft_n_chron = 0

   real(dp), allocatable :: all_energies(:, :)
   real(dp), allocatable :: all_forces(:, :, :)
   real(dp), allocatable :: all_virial(:, :, :)
   real(dp), allocatable :: all_this_energies(:, :)
   real(dp), allocatable :: all_this_forces(:, :, :)
   real(dp), allocatable :: all_this_virial(:, :, :)
   real(dp) :: instant_temp
   real(dp) :: kB = 8.6173303d-5
   real(dp) :: E_kinetic = 0.d0
   real(dp) :: E_kinetic_prev
   real(dp) :: time1
   real(dp) :: time2
   real(dp) :: time3
!   Every wall-clock bucket lives in one times_t (src/timing.f90), so the
!   extracted modules take a single argument instead of thirteen and the two
!   branches' signatures agree. time_step and time_step_prev below are the MD
!   integration step in fs, not timers, and deliberately stay separate.
   type(times_t) :: time
   real(dp) :: instant_pressure
   real(dp) :: time_step
   real(dp) :: md_time
   real(dp) :: instant_pressure_prev
   integer, allocatable :: mc_id(:)
   integer :: gd_istep = 0
   logical :: write_condition = .false.
   logical :: overwrite_condition = .false.
   character*1 :: creturn = achar(13)

  !! these decalarations are for time step and electronic stopping by different methods
   real(dp) :: time_step_prev
   integer :: nrows
   real(dp) :: cum_EEL = 0.0d0
   real(dp), allocatable :: allelstopdata(:)
   type(EPH_Beta_class) :: ephbeta
   type(EPH_FDM_class) :: ephfdm
   type(EPH_LangevinSpatialCorrelation_class) :: ephlsc

   ! Clean up these variables after code refactoring !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer, allocatable :: alpha_max(:)
   integer, allocatable :: der_neighbors(:)
   integer, allocatable :: der_neighbors_list(:)
   integer, allocatable :: i_beg_list(:)
   integer, allocatable :: i_end_list(:)
   integer, allocatable :: j_beg_list(:)
   integer, allocatable :: j_end_list(:)
   integer, allocatable :: species_idx(:)
   integer, allocatable :: n_mc_species(:)
   integer, allocatable :: n_mc_species_prev(:)
   integer :: i
   integer :: j
   integer :: k
   integer :: i2
   integer :: n_soap
   integer :: k2
   integer :: n_sites_this
   integer :: ierr
   integer :: rank
   integer :: ntasks
   type(comm_t) :: comm
   type(state_t) :: state
   type(domain_t) :: dom
   type(neighbors_t) :: nl
   type(model_t), target :: model
   type(results_t), target :: res
   type(loop_t) :: loop
   integer :: n_sp
   integer :: n_pos
   integer :: this_i_beg
   integer :: this_i_end
   integer :: this_j_beg
   integer :: this_j_end
   integer :: this_n_sites_mpi
   integer :: n_lp_count = 0

   integer :: l_max
   integer :: n_max
   integer :: central_species = 0
   integer :: iostatus
   integer :: counter2

!   The eleven additive contribution families that are reduced together after the
!   descriptor loop. Their predicates used to be written out three times -- to
!   count the slots, to pack them and to unpack them -- and evaluated
!   independently each time. Two copies disagreeing shifts counter2 and
!   silently attributes one term's energies to another. That is the same shape
!   as the ts+mbd predicate defect (KNOWN_ISSUES.md #1), so it is killed the
!   same way: evaluated once, into contrib_on, and only read thereafter.
   integer, parameter :: C_SOAP = 1
   integer, parameter :: C_VDW = 2
   integer, parameter :: C_ESTAT = 3
   integer, parameter :: C_LP = 4
   integer, parameter :: C_PDF = 5
   integer, parameter :: C_SF = 6
   integer, parameter :: C_XRD = 7
   integer, parameter :: C_ND = 8
   integer, parameter :: C_2B = 9
   integer, parameter :: C_CP = 10
   integer, parameter :: C_3B = 11
   integer, parameter :: N_CONTRIB = 11
   logical :: contrib_on(1:N_CONTRIB)
   type(contribution_ref) :: contrib(1:N_CONTRIB)
   type(perform_t) :: perform
   integer :: n_active
   integer :: i_contrib
   integer :: which_atom = 0
   integer :: n_omp = 1
   integer :: radial_enhancement = 0
   integer :: mc_mu_id = 1
   logical :: do_mc_relax = .false.
!  Whether the trial about to be tested came out of an MD or relaxation burst,
!  captured before params%do_md is cleared just below.
   logical :: trial_came_from_md = .false.
!  Molecule bookkeeping for grand-canonical moves that exchange whole molecules
!  (mc_molecule_files). mc_mol_id tags each atom with the inserted copy it
!  belongs to, mc_mol_mu with the chemical potential that copy came from, and
!  both are zero for a free atom. mc_mol_next hands out the tags. Rank 0 only:
!  the energy evaluation never reads any of it, so none of it is broadcast.
   integer, allocatable :: mc_mol_id(:)
   integer, allocatable :: mc_mol_mu(:)
   integer :: mc_mol_next = 0

   character*1024 :: filename
   character*1024 :: mc_file = "mc_trial.xyz"
   character*1024 :: string
   character*1024 :: temp_string
   character*1024 :: temp_string2
   character*8 :: i_char

   ! This is the mode in which we run TurboGAP
   character*16 :: mode = "none"
   character*16 :: help_topic = ""
   character*32 :: mc_move = "none"
   character*32 :: exp_output = "none"

   ! Here we store the input parameters
   type(input_parameters) :: params

   integer :: temp_md_nsteps

   logical :: do_electrostatics = .true.
! Persistent ts+mbd correction state, owned by turbogap_vdw
   type(vdw_state) :: vdw_ws

   ! Nested sampling
   real(dp) :: e_max
   real(dp) :: e_kin
   real(dp) :: rand
   real(dp) :: rand_scale(1:6)
   real(dp) :: target_temp
   integer :: i_nested
   integer :: i_max
   integer :: i_image
   integer :: i_current_image = 1
   integer :: i_trial_image = 2
   type(image), allocatable :: images(:)
   type(image), allocatable :: images_temp(:)
   character*32 :: implemented_exp_observables(1:5)
#ifdef _GPU
   integer :: omp_task
   integer(c_size_t) :: st_size_nf
   type(c_ptr) :: alphas_d
   type(c_ptr) :: qs_d
   type(c_ptr) :: nf_d
   type(c_ptr) :: rcut_hard_d
   type(c_ptr) :: rcut_soft_d
   type(c_ptr) :: global_scaling_d
   type(c_ptr) :: atom_sigma_r_d
   type(c_ptr) :: atom_sigma_r_scaling_d
   type(c_ptr) :: atom_sigma_t_d
   type(c_ptr) :: atom_sigma_t_scaling_d
   type(c_ptr) :: amplitude_scaling_d
   type(c_ptr) :: alpha_max_d
   type(c_ptr) :: central_weight_d
   integer :: n_sparse
   integer :: dim
   integer :: n_pairs_temp
   integer :: n_sites_temp
   character*8, allocatable, target :: species_types_actual(:)
#endif

!  --help answers from the generated keyword reference and exits. It is the
!  first thing the program does because it must work with no input file, no
!  GPU and no MPI: everything below this point assumes at least one of those.
   call get_command_argument(1, mode)
   if (mode == "--help" .or. mode == "-h" .or. mode == "help") then
      call get_command_argument(2, help_topic)
!     Validated against the SAME list the error message prints, which
!     keyword_help.f90 generates from tools/keyword_docs.py. It used to be a
!     hardcoded chain of comparisons beside a generated message, and the two
!     drifted the moment a mode was added: --help ipi was rejected by a message
!     that listed ipi as valid. Slashes on both sides so that a topic cannot
!     match a substring of another.
      if (len_trim(help_topic) > 0 .and. &
          index("/"//trim(keyword_help_topics())//"/", "/"//trim(help_topic)//"/") == 0) then
         write (*, '(A)') 'ERROR: unknown help topic "'//trim(help_topic)// &
            '". turbogap --help ['//trim(keyword_help_topics())//']'
         stop 1
      end if
      call print_keyword_help(help_topic)
      stop
   end if
   mode = "none"

   implemented_exp_observables(1) = "xps"
   implemented_exp_observables(2) = "xrd"
   implemented_exp_observables(3) = "saxs"
   implemented_exp_observables(4) = "pair_distribution"
   implemented_exp_observables(5) = "structure_factor"

!  Bring the device context up. Empty on this branch (src/gpu_context.f90);
!  on the GPU branch the same two names create the streams and cuBLAS handles.
   call time_start(time%create_streams)
   call gpu_context_init(params, rank, n_omp)
   call gap_backend_init()
   call time_end(time%create_streams)

   ! Start recording the time
   call get_time(time1)
   time3 = time1
!  Everything before the first pass of the main loop: the input file, the
!  potential files, and the allocation and broadcast that follow them. Without
!  this bucket the pre-loop cost fell into Miscellaneous, which is why a run
!  whose real work took 1.2 s reported 0.4 s "miscellaneous" and a 31 s run
!  reported the same 0.4 s -- a constant, and therefore obviously a setup cost,
!  but not one the report could name.
   call time_start(time%setup)
   ! Start random seed
   call srand(int(time1*1000))

   call comm_init(comm)
   rank = comm%rank
   ntasks = comm%size
   allocate (dom%n_atom_pairs_by_rank(1:ntasks))

   ! Read the mode. It should be "soap", "predict" or "md"
   call get_command_argument(1, mode)
   if (mode == "" .or. mode == "none") then
      write (*, *) "ERROR: you need to run 'turbogap md', 'turbogap mc', 'turbogap predict'"
      write (*, *) "       or 'turbogap ipi' (forces for an i-PI server; see ipi_address)"
      write (*, *) "       'turbogap --help [predict|md|mc|soap|gap]' lists the keywords"
      stop
      ! THIS SHOULD BE FIXED, IN CASE THE USER JUST WANT TO OUTPUT THE SOAP DESCRIPTORS
      mode = "soap"
   end if

   ! Prints some welcome message and reads in the input file
   if (rank == 0) then
      write (*, *) '_________________________________________________________________ '
      write (*, *) '                             _                                   \'
      write (*, *) ' ___________            __   \\ /\        _____     ___   _____  |'
      write (*, *) '/____  ____/           / / /\|*\|*\/\    / ___ \   /   | |  _  \ |'
      write (*, *) '    / / __  __  __    / /  \********/   / /  /_/  / /| | | / | | |'
      write (*, *) '   / / / / / / / /_  / /__  \**__**/   / / ____  / / | | | |_/ / |'
      write (*, *) '  / / / / / / / __/ / ___ \ /*/  \*\  / / /_  / / /__| | |  __/  |'
      write (*, *) ' / / / /_/ / / /   / /__/ / \ \__/ / / /___/ / / ____  | | |     |'
      write (*, *) '/_/_/_____/_/_/___/______/___\____/__\______/_/_/____|_|_|_|____ |'
      write (*, *) '_____________________________________________________________  / |'
      write (*, *) '*************************************************************|/  |'
      write (*, *) '                  Welcome to the TurboGAP code                   |'
      write (*, *) '                         Maintained by                           |'
      write (*, *) '                                                                 |'
      write (*, *) '               Miguel A. Caro and Tigany Zarrouk                 |'
      write (*, *) '                       mcaroba@gmail.com                         |'
      write (*, *) '                      miguel.caro@aalto.fi                       |'
      write (*, *) '                                                                 |'
      write (*, *) '          Department of Chemistry and Materials Science          |'
      write (*, *) '                     Aalto University, Finland                   |'
      write (*, *) '                                                                 |'
      write (*, *) '.................................................................|'
      write (*, *) '                                                                 |'
      write (*, *) '====================>>>>>  turbogap.fi  <<<<<====================|'
      write (*, *) '                                                                 |'
      write (*, *) '.................................................................|'
      write (*, *) '                                                                 |'
      write (*, *) 'Contributors (code and methodology) in chronological order:      |'
      write (*, *) '                                                                 |'
      write (*, *) 'Miguel A. Caro, Patricia Hernández-León, Suresh Kondati          |'
      write (*, *) 'Natarajan, Albert P. Bartók, Eelis V. Mielonen, Heikki Muhli,    |'
      write (*, *) 'Mikhail Kuklin, Gábor Csányi, Jan Kloppenburg, Richard Jana,     |'
      write (*, *) 'Tigany Zarrouk, Cristian V. Achim                                |'
      write (*, *) '                                                                 |'
      write (*, *) '.................................................................|'
      write (*, *) '                                                                 |'
      write (*, *) '                     Last updated: Jun. 2026                     |'
      write (*, *) '                                        _________________________/'
      write (*, *) '.......................................|'
      if (comm_with_mpi) then
         write (*, *) '                                       |'
         write (*, *) 'Running TurboGAP with MPI support:     |'
         write (*, *) '                                       |'
         write (*, '(A,I6,A)') ' Running TurboGAP on ', ntasks, ' MPI tasks   |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      else
         write (*, *) '                                       |'
         write (*, *) 'Running the serial version of TurboGAP |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end if

   ! Read input file and other files
   call read_input_and_gap_files(mode, rank, ntasks, params, &
                                 model%soap_turbo_hypers, model%distance_2b_hypers, model%angle_3b_hypers, model%core_pot_hypers, &
                                 model%n_soap_turbo, model%n_distance_2b, model%n_angle_3b, model%n_core_pot, model%n_species, &
                                 model%rcut_max, &
                                 model%valid_xps, model%xps_idx, model%vdw_lp_index, model%core_be_lp_index, &
                                 model%valid_estat_charges, model%charge_lp_index, &
                                 model%local_property_labels, model%local_property_indexes, model%n_local_properties_mpi, &
                                 model%has_local_properties_mpi, model%local_properties_n_sparse_mpi_soap_turbo, &
                                 model%local_properties_dim_mpi_soap_turbo, nrows, allelstopdata, &
                                 ephbeta, ephfdm, ephlsc, time)

!  The host memory budget, which has to sit exactly here.
!
!  After read_input_and_gap_files, because it reads mem_fraction and writes
!  max_Gbytes_per_process -- placed next to gpu_context_init it would run before
!  the input existed and size the loop from defaults whatever the input said.
!
!  ntasks is passed for the case where MPI cannot say how the ranks are laid
!  out; the routine prefers to ask MPI which ranks share a node, because that is
!  the set that shares the memory it is dividing.
!
!  The GPU branch calls the same name in the same place, where it budgets from
!  the device instead of from the node.
   call gpu_memory_budget_init(params, rank, ntasks)

   ! <----------------------------------------------------------------------------------------------- Finish printouts
   if (rank == 0) then
      ! Print out chosen options:
      write (*, *) '                                       |'
      write (*, '(1X,A)') 'You specified the following options:   |'
      write (*, *) '                                       |'
      write (*, *) '---------------------------------      |'
      if (len(trim(params%atoms_file)) > 20) then
         write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '...   |'
      else
         write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '      |'
      end if
      write (*, *) '---------------------------------      |'
      write (i_char, '(I8)') model%n_species
      write (*, '(1X,A,A8,A)') 'No. of species   = ', adjustl(i_char), '            |'
      do i = 1, model%n_species
         write (i_char, '(I8)') i
         write (*, '(1X,A,A2,A,A8,A)') '  *) Species #', adjustl(i_char), ' =       ', adjustr(params%species_types(i)), '      |'
      end do
      write (*, *) '---------------------------------      |'
      write (*, '(1X,A,F15.4,A)') 'rcut_max = ', model%rcut_max, ' Angst.      |'
      write (*, *) '---------------------------------      |'
      write (*, *) '                                       |'
      write (*, *) '.......................................|'
   end if

   ! Print progress bar and initialize timers

   model%xps_idx = params%xps_idx
   loop%md_istep = -1
   loop%mc_istep = -1
   loop%n_xyz = 0
   i_nested = 0
   i_image = 0

   if (params%do_md) then
      if (rank == 0) then
         write (*, *) '                                       |'
         write (*, *) 'Doing molecular dynamics...            |'
         if (params%print_progress .and. loop%md_istep > 0) then
            write (*, *) '                                       |'
            write (*, *) 'Progress:                              |'
            write (*, *) '                                       |'
            write (*, '(1X,A)', advance='no') '[                                    ] |'
         end if
      end if
      loop%update_bar = params%md_nsteps/36
      if (loop%update_bar < 1) then
         loop%update_bar = 1
      end if
      loop%counter = 1
   end if

   ! This checks if we need to do the SOAP calculation more than once, if there are several concatenated
   ! structures in the xyz file provided or we're doing molecular dynamics

!   The exp-observable decisions, evaluated once.  Every input is a params
!   field or valid_xps, none of which changes inside the main loop.
!
!   This closes a defect.  The allocation guards asked do_X .and. valid_X, the
!   zeroing guards asked do_X .and. exp_forces .and. valid_X, and the force
!   accumulation asked only exp_forces .and. valid_X -- so a deck supplying an
!   experimental dataset for an observable it had not switched on, with
!   exp_forces set, accumulated forces_X and virial_X that the allocation
!   guard had skipped.  do_X and valid_X are independent: valid_X is set from
!   a label in the experimental data file, do_X is its own input keyword.
!   Same shape as the electrostatics guard and as has_vdw against
!   has_local_properties.
   perform%pdf = params%do_pair_distribution .and. params%valid_pdf
   perform%sf = params%do_structure_factor .and. params%valid_sf
   perform%xrd = params%do_xrd .and. params%valid_xrd
   perform%nd = params%do_nd .and. params%valid_nd

   perform%pdf_forces = perform%pdf .and. params%exp_forces
   perform%sf_forces = perform%sf .and. params%exp_forces
   perform%xrd_forces = perform%xrd .and. params%exp_forces
   perform%nd_forces = perform%nd .and. params%exp_forces
   perform%xps_forces = model%valid_xps .and. params%exp_forces

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
   ir_from_traj = params%do_ir .and. .not. params%do_md .and. .not. params%do_mc
   ir_fft_active = ir_from_traj .or. &
                   ((params%valid_ir .or. params%do_ir) .and. &
                    trim(params%ir_bias_mode) == "fft")

!  One config for both routes; only dt_fs differs, and it is filled in at the
!  point of use because in one route it comes from the file and in the other
!  from md_step*ir_stride. normalise is off for MD: under a bias the spectrum
!  is fitted against the experiment with a scale, and dividing by max|I| as
!  well would be a second, discontinuous, normalisation of the same freedom.
   if (ir_fft_active) then
      call ir_fft_config_from_params(1.d0, params%ir_window, params%ir_fft_acf_ratio, &
                                     params%ir_nu_max, params%ir_fft_smooth_k, &
                                     params%ir_fft_smooth_kind, &
                                     params%ir_fft_quantum_correction, &
                                     params%ir_fft_temperature, params%t_beg, &
                                     params%ir_fft_power_dc_cutoff, &
                                     params%ir_subtract_mean, ir_from_traj, ir_fft_cfg)
   end if

   if (ir_from_traj) then
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
      if (rank == 0) then
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

   call time_end(time%setup)

!  Connect before the first force call, so that a missing or unstarted i-PI
!  server is reported now rather than after the first GAP evaluation.
   if (mode == "ipi") call ipi_driver_open(params%ipi_address, rank)

   do while (loop%repeat_xyz .or. (params%do_md .and. loop%md_istep < params%md_nsteps) &
             .or. (params%do_mc .and. loop%mc_istep < params%mc_nsteps))
      loop%exit_loop = .false.

!     One stamp per iteration, so that the cost of a step with the IR bias
!     active can be compared with the cost of one without. Closed at the
!     bottom of the loop.
      if (params%valid_ir) call get_time(mad_ir_step_beg)

      if (params%do_mc) then
         loop%mc_istep = loop%mc_istep + 1
         ! Undo if the step is md related
         if (loop%md_istep > -1) loop%mc_istep = loop%mc_istep - 1

      end if

      if (params%do_md) then
         loop%md_istep = loop%md_istep + 1
      else
         loop%n_xyz = loop%n_xyz + 1
      end if

      !   Update progress bar
      !
      !   md_nsteps = 0 is a legitimate input -- one configuration, forces, no
      !   dynamics -- and the bar divides by it. Integer division by zero is a
      !   SIGFPE, so the run died on the first step with a backtrace and no
      !   message rather than producing the single frame it was asked for.
      loop%bar_frac = 0
      if (params%md_nsteps > 0) loop%bar_frac = 36*loop%md_istep/params%md_nsteps
      if (loop%bar_frac > 36) loop%bar_frac = 36
      if (params%print_progress .and. loop%counter == loop%update_bar .and. (.not. params%do_mc)) then
         if (rank == 0) then
            do j = 1, 36 + 3
               write (*, "(A)", advance="no") creturn
            end do
            write (*, "(1X,A)", advance="no") "["
            do i = 1, loop%bar_frac
               write (*, "(A)", advance="no") "."
            end do
            do i = loop%bar_frac + 1, 36
               write (*, "(A)", advance="no") " "
            end do
            write (*, "(A)", advance="no") "] |"
            if (loop%md_istep == params%md_nsteps) then
               write (*, *)
            end if
         end if
         loop%counter = 1
      else
         loop%counter = loop%counter + 1
      end if

      !   This chunk of code does all the reading/neighbor builds etc for each snapshot
      !   or MD step
      !   Read in XYZ file and build neighbors lists

      if ((params%do_md .and. loop%md_istep == 0)) then
         call time_start(time%read_xyz)
         if (rank == 0) then
            if (loop%mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites,.not. params%mc_write_xyz, state%fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz, params%randomize_velocities)
               nl%rebuild_neighbors_list = .true.

            else if (.not. params%do_nested_sampling .or. loop%mc_istep == 0) then
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .false., state%fix_atom, params%t_beg, &
                             params%write_array_property(6), .false., params%randomize_velocities)

            end if

            ! call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
            !               n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
            !               positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
            !               xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
            !               n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), .true. )
            !     Only rank 0 handles these variables
            !      allocate( positions_prev(1:3, 1:size(positions,2)) )
            !      allocate( positions_diff(1:3, 1:size(positions,2)) )
            if (.not. allocated(state%forces_prev)) allocate (state%forces_prev(1:3, 1:state%n_sites))
            if (.not. allocated(state%positions_prev)) allocate (state%positions_prev(1:3, 1:state%n_sites))
            if (.not. allocated(state%positions_diff)) allocate (state%positions_diff(1:3, 1:state%n_sites))
            state%positions_diff = 0.d0
            nl%rebuild_neighbors_list = .true.
         end if
         call time_end(time%read_xyz)
         !     If we're doing MD, we don't read beyond the first snapshot in the XYZ file
         loop%repeat_xyz = .false.
         !     At the moment, we can't do prediction if the unit cell doesn't fit a whole cutoff sphere
         if (rank == 0) then
            !     CLEAN THIS UP <------------------------------------------------------------------- LOOK HERE
            !      if( size(positions,2) /= n_sites )then
            if (.false.) then
               write (*, *) "Sorry, at the moment TurboGAP can't do MD for unit cells smaller than ", &
                  "a cutoff sphere <-- ERROR"
               call comm_finalize(comm)
               stop
            end if
         end if
      else if (.not. params%do_md) then
         call time_start(time%read_xyz)
         if (rank == 0) then
            if (loop%mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites,.not. params%mc_write_xyz, state%fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz, params%randomize_velocities)
               nl%rebuild_neighbors_list = .true.
            else
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .false., state%fix_atom, params%t_beg, params%write_array_property(6), &
                             .false., params%randomize_velocities, state%frame_time, state%has_frame_time)
            end if
         end if
         call time_end(time%read_xyz)
         call time_start(time%mpi)
         call comm_bcast(comm, loop%repeat_xyz)
!        The frame's time label. Every rank pushes the same dipole into the
!        same buffer -- the trajectory is the ensemble and it is replicated,
!        not distributed -- so every rank needs the same time with it.
         call comm_bcast(comm, state%frame_time)
         call comm_bcast(comm, state%has_frame_time)
         call time_end(time%mpi)
         nl%rebuild_neighbors_list = .true.
      end if
      !   Broadcast the info in the XYZ file: positions, velocities, masses, xyz_species, xyz_species_supercell,
      !   species, species_supercell, indices, a_box, b_box, c_box and n_sites. I should put this into a module!!!!!!!

      if (rank == 0) then
         if (params%randomize_velocities .and. loop%md_istep == 0) then
            call randomize_velocities(state%velocities, state%n_sites, E_kinetic, state%masses, instant_temp, params%t_beg, &
                                      params%velocity_distribution)
         end if
         if (params%do_mc .and. (mc_move /= "md" .or. loop%md_istep == 0) .and. params%mc_hamiltonian) then
            if (loop%mc_istep > 0) E_kinetic_prev = E_kinetic
            call random_number(state%velocities)
            call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), state%masses(1:state%n_sites))
            E_kinetic = 0.d0
            do i = 1, state%n_sites
               E_kinetic = E_kinetic + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), state%velocities(1:3, i))
            end do
            instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/kB*E_kinetic
            state%velocities = state%velocities*dsqrt(params%t_beg/instant_temp)
            if (loop%mc_istep > 0) then
               E_kinetic = E_kinetic_prev
               instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/kB*E_kinetic
               ! Reversing as we want it to be at the instant temp and not at t_beg
               state%velocities = state%velocities*dsqrt(instant_temp/params%t_beg)

               do i = 1, state%n_sites
                  E_kinetic = E_kinetic + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), state%velocities(1:3, i))
               end do
            else
               E_kinetic = E_kinetic*params%t_beg/instant_temp
            end if
         end if

      end if
      call domain_sync_state(dom, comm, state, params, time)
      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. loop%md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
      end if
      call time_start(time%neigh)
      !   Parallel neighbors list build
      call comm_bcast(comm, nl%rebuild_neighbors_list)

      !   If we're using a box rescaling algorithm or a barostat, then the box size can
      !   become smaller or bigger than the cutoff sphere. If that happens, and the current
      !   situation is different from before, then we need to figure out if we need to
      !   construct a supercell (i.e., the box was bigger than the cutoff sphere and now
      !   is smaller -> makes computations slower) or default back to the primitive unit cell
      !   (i.e., the box was smaller and now is bigger -> makes computations faster).
      !   We only need to check if rebuild_neighbors_list = .true.
      if (nl%rebuild_neighbors_list .and. params%do_mc .and. loop%mc_istep > 0) then
         call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                       model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                       state%positions, params%do_md, state%velocities, params%masses_types, state%masses, state%xyz_species, &
                       state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                       state%b_box, state%c_box, &
                       state%n_sites, .true., state%fix_atom, params%t_beg, &
                       params%write_array_property(6), .false., params%randomize_velocities)
      else if (nl%rebuild_neighbors_list) then
         call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                       model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                       state%positions, params%do_md, state%velocities, params%masses_types, state%masses, state%xyz_species, &
                       state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                       state%b_box, state%c_box, &
                       state%n_sites, .true., state%fix_atom, params%t_beg, params%write_array_property(6), &
                       .false., params%randomize_velocities)

      end if

      !   Overlapping domain decomposition with subcommunicators goes here <------------------- TO DO

      !   This is some trivial MPI parallelization to make sure the code works fine
      if (rank < mod(state%n_sites, ntasks)) then
         dom%i_beg = 1 + rank*(state%n_sites/ntasks + 1)
      else
         dom%i_beg = 1 + mod(state%n_sites, ntasks)*(state%n_sites/ntasks + 1) + (rank - mod(state%n_sites, &
                                                                                             ntasks))*(state%n_sites/ntasks)
      end if
      if (rank < mod(state%n_sites, ntasks)) then
         dom%i_end = (rank + 1)*(state%n_sites/ntasks + 1)
      else
         dom%i_end = dom%i_beg + state%n_sites/ntasks - 1
      end if

      dom%do_list = .false.
      dom%do_list(dom%i_beg:dom%i_end) = .true.

      call build_neighbors_list(state%positions, state%a_box, state%b_box, state%c_box, params%do_timing, &
                                state%species_supercell, model%rcut_max, nl%n_atom_pairs, nl%rjs, &
                                nl%thetas, nl%phis, nl%xyz, nl%n_neigh_local, nl%neighbors_list, nl%neighbor_species, &
                                state%n_sites, state%indices, &
                                nl%rebuild_neighbors_list, dom%do_list, rank)
      if (nl%rebuild_neighbors_list) then
         !     Get total number of atom pairs
         call comm_allgather(comm, nl%n_atom_pairs, dom%n_atom_pairs_by_rank)
         nl%n_atom_pairs_total = sum(dom%n_atom_pairs_by_rank)
         nl%n_atom_pairs = nl%n_atom_pairs_total

         !     Get number of neighbors
         if (.not. allocated(nl%n_neigh)) allocate (nl%n_neigh(1:state%n_sites))
         call comm_sum_to_root(comm, nl%n_neigh_local, nl%n_neigh, state%n_sites)
         call comm_bcast(comm, nl%n_neigh, state%n_sites)

         dom%j_beg = 1
         dom%j_end = dom%n_atom_pairs_by_rank(rank + 1)
      end if
!   Store by which rank each site is being handled
      if (allocated(dom%site_in_rank)) then
         if (size(dom%site_in_rank) /= state%n_sites) then
            deallocate (dom%site_in_rank, dom%this_site_in_rank)
         end if
      end if
      if (.not. allocated(dom%site_in_rank)) then
         allocate (dom%site_in_rank(1:state%n_sites))
         allocate (dom%this_site_in_rank(1:state%n_sites))
      end if
      dom%site_in_rank = 0
      dom%this_site_in_rank = 0
      do i = dom%i_beg, dom%i_end
         dom%this_site_in_rank(i) = rank
      end do
      call comm_sum_to_root(comm, dom%this_site_in_rank, dom%site_in_rank, state%n_sites)
      call comm_bcast(comm, dom%site_in_rank, state%n_sites)
      !   Compute the volume of the "primitive" unit cell
      state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                               state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
      call time_end(time%neigh)

      !   If we are doing prediction, we run this chunk of code
      if (params%do_prediction .or. params%write_soap .or. params%write_derivatives) then

         !     We only need to reallocate the arrays if the number of sites changes
         ! REMOVE TRUE FROM IF STATEMENT
         if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
            if (allocated(res%energies)) deallocate (res%energies, &
                                                     res%energies_soap, &
                                                     res%energies_2b, &
                                                     res%energies_3b, &
                                                     res%energies_core_pot, &
                                                     res%this_energies, &
                                                     res%energies_vdw, &
                                                     res%energies_vdw_corr, &
                                                     res%mbd_ts_scaling, &
                                                     res%this_forces, &
                                                     res%energies_lp, &
                                                     res%energies_exp, &
                                                     res%energies_estat, &
                                                     res%this_mbd_ts_scaling)
            allocate (res%energies(1:state%n_sites))
            allocate (res%this_energies(1:state%n_sites))
            allocate (res%energies_soap(1:state%n_sites))
            allocate (res%energies_2b(1:state%n_sites))
            allocate (res%energies_3b(1:state%n_sites))
            allocate (res%energies_core_pot(1:state%n_sites))
            allocate (res%energies_vdw(1:state%n_sites))
            allocate (res%energies_vdw_corr(1:state%n_sites))
            allocate (res%energies_lp(1:state%n_sites))
            allocate (res%energies_estat(1:state%n_sites))
            allocate (res%energies_exp(1:state%n_sites))
!          We do this allocations for van der Waals corrections
            allocate (res%mbd_ts_scaling(1:state%n_sites))
            allocate (res%this_mbd_ts_scaling(1:state%n_sites))
!          Allocated whether or not a dipole model is loaded: they are passed to
!          get_gap_soap unconditionally, and 4 doubles per atom is not worth a
!          second code path.
            if (allocated(res%local_dipoles)) deallocate (res%local_dipoles, res%this_local_dipoles, &
                                                          res%energies_dipole, res%this_energies_dipole)
            allocate (res%local_dipoles(1:3, 1:state%n_sites))
            allocate (res%this_local_dipoles(1:3, 1:state%n_sites))
            allocate (res%energies_dipole(1:state%n_sites))
            allocate (res%this_energies_dipole(1:state%n_sites))

! Read in file for ts+mbd van der Waals mode if it exists
! Initialise the TS scaling factors. md_istep <= 0 rather than == 0 because
! predict and mc never advance md_istep past -1, and without this they reach
! get_ts_energy_and_forces with mbd_ts_scaling never having been set. For MD
! this is still exactly the first step, so the MD path is unchanged.
            if (params%vdw_type == "ts+mbd" .and. loop%md_istep <= 0) then
               if (rank == 0) then
                  open (unit=30, file="mbd_ts_scaling.dat", status="old", iostat=iostatus)
                  if (iostatus == 0) then
                     write (*, *) '                                       |'
                     write (*, *) '.......................................|'
                     write (*, *) '                                       |'
                     write (*, *) 'Reading TS scaling factors from file   |'
                     write (*, *) 'mbd_ts_scaling.dat                     |'
                     write (*, *) '                                       |'
                     write (*, *) '.......................................|'
                     write (*, *) '                                       |'
                     do i = 1, state%n_sites
                        read (30, *) res%mbd_ts_scaling(i)
                     end do
                     res%update_mbd_ts_scaling = .false.
                  else
                     res%mbd_ts_scaling = 1.d0
                  end if
                  close (30)
                  res%this_mbd_ts_scaling = res%mbd_ts_scaling
               end if
               call comm_bcast(comm, res%this_mbd_ts_scaling, state%n_sites)
            end if

            if (perform%pdf) then
               if (allocated(res%energies_pdf)) deallocate (res%energies_pdf)
               allocate (res%energies_pdf(1:state%n_sites))
            end if

            if (perform%sf) then
               if (allocated(res%energies_sf)) deallocate (res%energies_sf)
               allocate (res%energies_sf(1:state%n_sites))
            end if

            if (perform%xrd) then
               if (allocated(res%energies_xrd)) deallocate (res%energies_xrd)
               allocate (res%energies_xrd(1:state%n_sites))
            end if

            if (perform%nd) then
               if (allocated(res%energies_nd)) deallocate (res%energies_nd)
               allocate (res%energies_nd(1:state%n_sites))
            end if

            !       This needs to be allocated even if no force prediction is needed:
            allocate (res%this_forces(1:3, 1:state%n_sites))
         end if
         res%energies = 0.d0
         res%energies_soap = 0.d0
         res%energies_2b = 0.d0
         res%energies_3b = 0.d0
         res%energies_core_pot = 0.d0
         res%energies_vdw = 0.d0
         res%energies_estat = 0.d0
         res%energies_lp = 0.d0
         res%energies_exp = 0.d0
!        The dissimilarity accumulators belong to the same step as energies_exp
!        and are zeroed with it. This is the only point that knows a new
!        evaluation has begun; get_exp_energies is called once per observable
!        and mad_ir separately again, so neither can reset them itself.
         exp_dissimilarity = 0.d0
         exp_dissim_ref = 0.d0
         res%local_dipoles = 0.d0
         res%energies_dipole = 0.d0
         res%dipole = 0.d0

         if (perform%pdf) res%energies_pdf = 0.d0
         if (perform%sf) res%energies_sf = 0.d0
         if (perform%xrd) res%energies_xrd = 0.d0
         if (perform%nd) res%energies_nd = 0.d0

         ! Adding allocation of local properties

         ! Now one could use pointers such that hirshfeld_v(:) acts as an alias for local_properties(vdw_index,:)...
         if (any_has_local_properties(model%soap_turbo_hypers)) then
            if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
               if (allocated(res%local_properties)) then
                  nullify (res%this_local_properties_pt)
                  deallocate (res%this_local_properties, res%local_properties)
                  if (params%do_forces) then
                     nullify (res%this_local_properties_cart_der_pt)
                     deallocate (res%this_local_properties_cart_der, res%local_properties_cart_der)
                  end if
               end if
               allocate (res%local_properties(1:state%n_sites, 1:params%n_local_properties))
               allocate (res%this_local_properties(1:state%n_sites, 1:params%n_local_properties))
               res%this_local_properties_pt => res%this_local_properties

               !         I don't remember why this needs a pointer <----------------------------------------- CHECK

            end if
            res%local_properties = 0.d0

            if (params%do_forces) then
               if (dom%n_atom_pairs_by_rank(rank + 1) /= dom%n_atom_pairs_by_rank_prev) then
                  if (allocated(res%local_properties_cart_der)) deallocate (res%local_properties_cart_der, &
                                                                            res%this_local_properties_cart_der)
                  allocate (res%local_properties_cart_der(1:3, 1:dom%n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
                  allocate (res%this_local_properties_cart_der(1:3, 1:dom%n_atom_pairs_by_rank(rank + 1), &
                                                               1:params%n_local_properties))
               end if
               if (.not. allocated(res%local_properties_cart_der)) then
                  allocate (res%local_properties_cart_der(1:3, 1:dom%n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
                  allocate (res%this_local_properties_cart_der(1:3, 1:dom%n_atom_pairs_by_rank(rank + 1), &
                                                               1:params%n_local_properties))
               end if

               res%local_properties_cart_der = 0.d0
               res%this_local_properties_cart_der_pt =>&
                    & res%this_local_properties_cart_der(1:3,&
                    & 1:dom%n_atom_pairs_by_rank(rank + 1), 1:params&
                    &%n_local_properties)
            end if
         end if

         ! Now go through the soap turbo hypers, and see if any are vdw or
         ! otherwise, if vdw, one can have pointers to point to the data
         ! structures such that it makes things clearer. One needs to check
         ! that this allocation still works iwth if(allocated(hirsh_v))
         ! statements

         if (params%do_forces) then
            if (state%n_sites /= loop%n_sites_prev .or. params%do_mc) then
               if (allocated(res%forces)) deallocate (res%forces, res%forces_soap, res%forces_2b, res%forces_3b, &
                  res%forces_core_pot, res%forces_vdw,&
                    & res%forces_lp, res%forces_estat, res%local_virial_vdw_diag, res%local_virial_vdw_diag_corr)
               allocate (res%forces(1:3, 1:state%n_sites))
               allocate (res%forces_soap(1:3, 1:state%n_sites))
               allocate (res%forces_2b(1:3, 1:state%n_sites))
               allocate (res%forces_3b(1:3, 1:state%n_sites))
               allocate (res%forces_core_pot(1:3, 1:state%n_sites))
               allocate (res%forces_vdw(1:3, 1:state%n_sites))
               if (allocated(res%forces_vdw_corr)) deallocate (res%forces_vdw_corr)
               allocate (res%forces_vdw_corr(1:3, 1:state%n_sites))
               allocate (res%forces_lp(1:3, 1:state%n_sites))
               allocate (res%forces_estat(1:3, 1:state%n_sites))
               allocate (res%local_virial_vdw_diag_corr(1:3, 1:state%n_sites))
               allocate (res%local_virial_vdw_diag(1:3, 1:state%n_sites))

               if (perform%pdf_forces) then
                  if (allocated(res%forces_pdf)) deallocate (res%forces_pdf)
                  allocate (res%forces_pdf(1:3, 1:state%n_sites))
               end if

               if (perform%sf_forces) then
                  if (allocated(res%forces_sf)) deallocate (res%forces_sf)
                  allocate (res%forces_sf(1:3, 1:state%n_sites))
               end if

               if (perform%xrd_forces) then
                  if (allocated(res%forces_xrd)) deallocate (res%forces_xrd)
                  allocate (res%forces_xrd(1:3, 1:state%n_sites))
               end if

               if (perform%nd_forces) then
                  if (allocated(res%forces_nd)) deallocate (res%forces_nd)
                  allocate (res%forces_nd(1:3, 1:state%n_sites))
               end if

            end if
!           Decide here, before any descriptor is evaluated, whether this step
!           contributes a configuration to the IR ensemble: get_soap has to be
!           told to produce second derivatives before it builds anything, and
!           gap_interface reads mad_ir_collect to do that.
            if ((params%valid_ir .or. params%do_ir) .and. .not. ir_from_traj) then
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
                                       state%n_sites, mad_ir_ok, mad_ir_resumed, mad_ir_msg, &
                                       params%ir_acf_mode, params%ir_tau_mem)
                  else
!                    Prediction: the ensemble is the whole trajectory, so its
!                    length is md_nsteps/ir_stride + 1 -- every step for which
!                    modulo(md_istep, ir_stride) is zero, counting step zero.
                     mad_ir_resumed = .false.
!                    A negative resolution means "whatever the run gives"; the
!                    default value of ir_resolution cannot be told from a
!                    chosen one by its value, hence the flag.
                     if (params%ir_resolution_set) then
                        mad_ir_res_ask = params%ir_resolution
                     else
                        mad_ir_res_ask = -1.d0
                     end if
                     call mad_ir_setup_predict(params%md_step, params%ir_stride, &
                                               params%md_nsteps/params%ir_stride + 1, &
                                               mad_ir_res_ask, params%ir_nu_min, &
                                               params%ir_nu_max, params%ir_lag_factor, &
                                               params%ir_n_samples, params%ir_nu_power, &
                                               params%ir_window, params%ir_subtract_mean, &
                                               trim(params%ir_estimator) /= "unbiased", &
                                               params%ir_taper_partial, &
                                               state%n_sites, mad_ir_ok, mad_ir_msg, &
                                               params%ir_acf_mode, params%ir_tau_mem)
                  end if
                  if (.not. mad_ir_ok) then
                     write (*, *) "ERROR: ", trim(mad_ir_msg)
                     stop
                  end if
                  if (rank == 0) then
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
                     else if (mad_ir_resumed) then
                        write (*, '(A,I12,A)') '  *) resumed, frames:   ', mad_ir_state%n_stored, '         |'
                     else
                        write (*, *) '  *) fresh ensemble; no bias is       |'
                        write (*, *) '     applied until it is full.        |'
                        if (len_trim(mad_ir_msg) > 0) write (*, *) '     ', trim(mad_ir_msg)
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
                                          params%ir_xl_restart_file, mad_ir_xl_ok, &
                                          mad_ir_xl_resumed, mad_ir_xl_msg)
                     if (.not. mad_ir_xl_ok) then
                        write (*, *) "ERROR: ", trim(mad_ir_xl_msg)
                        stop
                     end if
                     mad_ir_xl_active = .true.
                     if (rank == 0) then
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
                        if (mad_ir_xl_resumed) then
                           write (*, '(A,I12,A)') '  *) resumed, advances: ', &
                              mad_ir_xl_state%n_steps, '         |'
                        else
                           write (*, '(A,I12,A)') '  *) charging, advances:', &
                              mad_ir_xl_state%n_warm, '         |'
                           if (len_trim(mad_ir_xl_msg) > 0) write (*, *) '     ', trim(mad_ir_xl_msg)
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
                                       ir_aux_ok, ir_aux_resumed, ir_aux_msg)
                     if (.not. ir_aux_ok) then
                        write (*, *) "ERROR: ", trim(ir_aux_msg)
                        stop
                     end if
                     if (rank == 0) then
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
                        if (ir_aux_resumed) then
                           write (*, '(A,I12,A)') '  *) resumed, advances: ', &
                              ir_aux_state%n_steps, '         |'
                        else
                           write (*, *) '  *) uncalibrated: the bank waits    |'
                           write (*, *) '     for the ACF ensemble to fill   |'
                           write (*, *) '     and applies no force until it  |'
                           write (*, *) '     does.                          |'
                        end if
                        if (len_trim(ir_aux_msg) > 0) write (*, *) '     ', trim(ir_aux_msg)
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
            res%forces = 0.d0
            res%forces_soap = 0.d0
            res%forces_2b = 0.d0
            res%forces_3b = 0.d0
            res%forces_core_pot = 0.d0
            res%forces_vdw = 0.d0
            res%forces_estat = 0.d0
            res%forces_lp = 0.d0
            res%virial = 0.d0
            res%virial_soap = 0.d0
            res%virial_2b = 0.d0
            res%virial_3b = 0.d0
            res%virial_core_pot = 0.d0
            res%virial_vdw = 0.d0
            res%virial_estat = 0.d0
            res%virial_lp = 0.d0
            res%local_virial_vdw_diag = 0.d0
            if (perform%pdf_forces) then
               res%forces_pdf = 0.d0
               res%virial_pdf = 0.d0
#ifdef _MPIF90
               res%this_virial_pdf = 0.d0
#endif
            end if

            if (perform%sf_forces) then
               res%forces_sf = 0.d0
               res%virial_sf = 0.d0
#ifdef _MPIF90
               res%this_virial_sf = 0.d0
#endif
            end if

            if (perform%xrd_forces) then
               res%forces_xrd = 0.d0
               res%virial_xrd = 0.d0
#ifdef _MPIF90
               res%this_virial_xrd = 0.d0
#endif
            end if

            if (perform%nd_forces) then
               res%forces_nd = 0.d0
               res%virial_nd = 0.d0
#ifdef _MPIF90
               res%this_virial_nd = 0.d0
#endif
            end if
         end if

         if (params%do_prediction) then
            !       Assign the e0 to each atom according to its species
            do i = dom%i_beg, dom%i_end
               do j = 1, model%n_species
                  if (state%xyz_species(i) == params%species_types(j)) then
                     res%energies(i) = params%e0(j)
                  end if
               end do
            end do
         end if
         !     Collect all energies
         call time_start(time%mpi_ef)
         call comm_sum_to_root(comm, res%energies, res%this_energies, state%n_sites)
         call time_end(time%mpi_ef)
         res%energies = res%this_energies

         !     Loop through soap_turbo descriptors - we always call this routine, even if we don't want to do prediction
         n_lp_count = 0 ! This counts the local properties
         call time_start(time%gap)
         do i = 1, model%n_soap_turbo
            call time_start(time%soap)
            !       Compute number of pairs for this SOAP. SOAP has in general a different cutoff than overall max
            !       cutoff, so the number of pairs may be a lot smaller for the SOAP subset.
            !       This subroutine splits the load optimally so as to not use more memory per MPI process than available.
            !       TurboGAP does not check how much memory is available, it just relies on heuristics and a user provided
            !       max_Gbytes_per_process (default = 1.d0)
#ifdef _GPU
            if (params%n_batches > 0) then
               call get_number_of_atom_pairs_batches(params%n_batches, nl%n_neigh(dom%i_beg:dom%i_end), &
                                                     nl%rjs(dom%j_beg:dom%j_end), model%soap_turbo_hypers(i)%rcut_max, &
                                                     model%soap_turbo_hypers(i)%l_max, &
                                                     model%soap_turbo_hypers(i)%n_max, &
                                                     model%soap_turbo_hypers(i)%dim, &
                                                     model%soap_turbo_hypers(i)%n_species, &
                                                     params%max_Gbytes_per_process, i_beg_list, &
                                                     i_end_list, j_beg_list, j_end_list)
            else
               call get_number_of_atom_pairs(nl%n_neigh(dom%i_beg:dom%i_end), nl%rjs(dom%j_beg:dom%j_end), &
                                             model%soap_turbo_hypers(i)%rcut_max, &
                                             model%soap_turbo_hypers(i)%l_max, &
                                             model%soap_turbo_hypers(i)%n_max, &
                                             model%soap_turbo_hypers(i)%dim, &
                                             model%soap_turbo_hypers(i)%n_species, &
                                             params%max_Gbytes_per_process, i_beg_list, &
                                             i_end_list, j_beg_list, j_end_list)
            end if

            n_sp = model%soap_turbo_hypers(i)%n_species

            st_size_nf = n_sp*sizeof(model%soap_turbo_hypers(i)%nf(1))
            call gpu_malloc_async(nf_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%nf), nf_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(rcut_hard_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%rcut_hard), rcut_hard_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(rcut_soft_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%rcut_soft), rcut_soft_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(global_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%global_scaling), global_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_r_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_r), atom_sigma_r_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_r_scaling), atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_t_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_t), atom_sigma_t_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_t_scaling), atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(amplitude_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%amplitude_scaling), amplitude_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(central_weight_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%central_weight), central_weight_d, st_size_nf, gpu_stream)
            st_size_nf = n_sp*sizeof(model%soap_turbo_hypers(i)%alpha_max(1))
            call gpu_malloc_async(alpha_max_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(model%soap_turbo_hypers(i)%alpha_max), alpha_max_d, st_size_nf, gpu_stream)
            n_sparse = model%soap_turbo_hypers(i)%n_sparse
            dim = model%soap_turbo_hypers(i)%dim
            call soap_backend_begin(model%soap_turbo_hypers(i))

            if (model%soap_turbo_hypers(i)%has_local_properties) then
               ! Allocate gpu memory
               do j = 1, model%soap_turbo_hypers(i)%n_local_properties
                  model%soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas = &
                     model%soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                     sizeof(model%soap_turbo_hypers(i)%local_property_models(j)%alphas(1))
                  call gpu_malloc_async(model%soap_turbo_hypers(i)%local_property_models(j)%alphas_d, &
                                        model%soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas, gpu_stream)
                  call cpy_htod(c_loc(model%soap_turbo_hypers(i)&
                    &%local_property_models(j)%alphas), &
                    & model%soap_turbo_hypers(i)%local_property_models(j)&
                    &%alphas_d, model%soap_turbo_hypers(i)&
                    &%local_property_models(j)%st_size_alphas,&
                    & gpu_stream)

                  model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs = &
                     model%soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                     model%soap_turbo_hypers(i)%local_property_models(j)%dim* &
                     sizeof(model%soap_turbo_hypers(i)%local_property_models(j)%Qs(1, 1))

                  call gpu_malloc_async(model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                                        model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, gpu_stream)
                  call cpy_htod(c_loc(model%soap_turbo_hypers(i)%local_property_models(j)%Qs), &
                                model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                                model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, &
                                gpu_stream)

               end do

            end if
#else
            call get_number_of_atom_pairs(nl%n_neigh(dom%i_beg:dom%i_end), nl%rjs(dom%j_beg:dom%j_end), &
                                          model%soap_turbo_hypers(i)%rcut_max, &
                                          model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%n_max, &
                                          model%soap_turbo_hypers(i)%dim, model%soap_turbo_hypers(i)%n_species, &
                                          params%max_Gbytes_per_process, i_beg_list, i_end_list, j_beg_list, j_end_list)
#endif

            do j = 1, size(i_beg_list)
               this_i_beg = dom%i_beg - 1 + i_beg_list(j)
               this_i_end = dom%i_beg - 1 + i_end_list(j)
               this_j_beg = dom%j_beg - 1 + j_beg_list(j)
               this_j_end = dom%j_beg - 1 + j_end_list(j)
               this_n_sites_mpi = this_i_end - this_i_beg + 1
               res%this_energies = 0.d0
               if (params%do_forces) then
                  res%this_forces = 0.d0
                  res%this_virial = 0.d0
               end if
               if (model%soap_turbo_hypers(i)%is_dipole_model) then
                  res%this_local_dipoles = 0.d0
                  res%this_energies_dipole = 0.d0
               end if
               if (model%soap_turbo_hypers(i)%has_local_properties) then
                  res%this_local_properties = 0.d0
                  if (params%do_forces) then
                     res%this_local_properties_cart_der = 0.d0
                     !             I don't remember why this needs a pointer <----------------------------------------- CHECK
                     nullify (res%this_local_properties_cart_der_pt)
                     res%this_local_properties_cart_der_pt =>&
                          & res%this_local_properties_cart_der(1:3,&
                          & this_j_beg:this_j_end, 1:params&
                          &%n_local_properties)
                  end if
               end if

#ifdef _GPU
               call get_gap_soap( &
                  n_sparse, state%n_sites, this_n_sites_mpi, nl%n_neigh(this_i_beg:this_i_end), &
                  nl%neighbors_list(this_j_beg:this_j_end), model%soap_turbo_hypers(i)%n_species, &
                  model%soap_turbo_hypers(i)%species_types, nl%rjs(this_j_beg:this_j_end), nl%thetas(this_j_beg:this_j_end), &
                  nl%phis(this_j_beg:this_j_end), nl%xyz(1:3, this_j_beg:this_j_end), alpha_max_d, &
                  model%soap_turbo_hypers(i)%alpha_max, model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%dim, &
                  rcut_hard_d, &
                  model%soap_turbo_hypers(i)%rcut_hard, rcut_soft_d, nf_d, global_scaling_d, atom_sigma_r_d, &
                  model%soap_turbo_hypers(i)%atom_sigma_r, atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
                  amplitude_scaling_d, model%soap_turbo_hypers(i)%radial_enhancement, central_weight_d, &
                  model%soap_turbo_hypers(i)%central_weight, model%soap_turbo_hypers(i)%basis, &
                  model%soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
                  params%do_prediction, params%write_soap, params%write_derivatives, &
                  model%soap_turbo_hypers(i)%compress_soap, model%soap_turbo_hypers(i)%compress_soap_indices, &
                  model%soap_turbo_hypers(i)%delta, model%soap_turbo_hypers(i)%zeta, model%soap_turbo_hypers(i)%central_species, &
                  state%xyz_species(this_i_beg:this_i_end), state%xyz_species_supercell, params%all_atoms, &
                  params%which_atom, state%indices, soap, soap_cart_der, der_neighbors, der_neighbors_list, &
                  model%soap_turbo_hypers(i)%has_local_properties, model%soap_turbo_hypers(i)%n_local_properties, &
                  model%soap_turbo_hypers(i)%local_property_models, n_lp_count, res%energies_soap, res%forces_soap, &
                  res%this_local_properties_pt, res%this_local_properties_cart_der_pt, model%local_property_indexes, &
                  res%this_virial, &
                  time%soap_lin(3), time%get_soap(3), model%soap_turbo_hypers(i)%W_d, model%soap_turbo_hypers(i)%S_d, &
                  model%soap_turbo_hypers(i)%multiplicity_array_d, model%soap_turbo_hypers(i)%st_W_d, &
                  model%soap_turbo_hypers(i)%st_S_d, model%soap_turbo_hypers(i)%st_multiplicity_array_d, &
                  model%soap_turbo_hypers(i)%recompute_basis, time%local_prop, &
                  model%soap_turbo_hypers(i)%is_dipole_model, res%local_dipoles, res%energies_dipole)
#else
               call soap_backend_begin(model%soap_turbo_hypers(i))
               call get_gap_soap(state%n_sites, this_n_sites_mpi, nl%n_neigh(this_i_beg:this_i_end), &
                  nl%neighbors_list(this_j_beg:this_j_end), &
                    model%soap_turbo_hypers(i)%n_species, model%soap_turbo_hypers(i)%species_types, &
                    nl%rjs(this_j_beg:this_j_end), nl%thetas(this_j_beg:this_j_end), nl%phis(this_j_beg:this_j_end), &
                    nl%xyz(1:3, this_j_beg:this_j_end), &
                    model%soap_turbo_hypers(i)%alpha_max, &
                    model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%dim, model%soap_turbo_hypers(i)%rcut_hard, &
                    model%soap_turbo_hypers(i)%rcut_soft, model%soap_turbo_hypers(i)%nf, &
                       model%soap_turbo_hypers(i)%global_scaling, &
                    model%soap_turbo_hypers(i)%atom_sigma_r, model%soap_turbo_hypers(i)%atom_sigma_r_scaling, &
                    model%soap_turbo_hypers(i)%atom_sigma_t, model%soap_turbo_hypers(i)%atom_sigma_t_scaling, &
                    model%soap_turbo_hypers(i)%amplitude_scaling, model%soap_turbo_hypers(i)%radial_enhancement, &
                    model%soap_turbo_hypers(i)%central_weight, model%soap_turbo_hypers(i)%basis, &
                    model%soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
                    params%do_prediction, params%write_soap, params%write_derivatives, &
                    model%soap_turbo_hypers(i)%compress_soap, model%soap_turbo_hypers(i)%compress_P_nonzero, &
                    model%soap_turbo_hypers(i)%compress_P_i, model%soap_turbo_hypers(i)%compress_P_j, &
                    model%soap_turbo_hypers(i)%compress_P_el, &
                    model%soap_turbo_hypers(i)%delta, model%soap_turbo_hypers(i)%zeta, model%soap_turbo_hypers(i)%central_species, &
                    state%xyz_species(this_i_beg:this_i_end), state%xyz_species_supercell, &
                    params%all_atoms, params%which_atom, state%indices, soap, soap_cart_der, &
                    der_neighbors, der_neighbors_list, &
                    & model%soap_turbo_hypers(i)%has_local_properties,&
                    & model%soap_turbo_hypers(i)%n_local_properties,&
                    & model%soap_turbo_hypers(i)%local_property_models,&
                    & res%this_energies, res%this_forces, res%this_local_properties_pt,&
                    & res%this_local_properties_cart_der_pt,&
                    & model%local_property_indexes, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                    & res%this_virial, n_lp_count, model%soap_turbo_hypers(i)%is_dipole_model, &
                    & res%this_local_dipoles, res%this_energies_dipole)

               call soap_backend_end()
#endif

               ! We can have a pointer to specific parts of this_local_properties array to then

!              A dipole descriptor leaves this_energies and this_forces at the
!              zero they were set to above -- get_gap_soap never writes them --
!              so its fictitious energy stays out of energies_soap and its
!              gradient out of forces_soap. It is carried separately.
               res%energies_soap = res%energies_soap + res%this_energies

               if (model%soap_turbo_hypers(i)%is_dipole_model) then
                  res%local_dipoles = res%local_dipoles + res%this_local_dipoles
                  res%energies_dipole = res%energies_dipole + res%this_energies_dipole
               end if

               if (model%soap_turbo_hypers(i)%has_local_properties) then

                  res%local_properties(:, :) = res%local_properties(:, :) + res%this_local_properties(:, :)
                  if (any(model%soap_turbo_hypers(i)&
                       &%local_property_models(:)%do_derivatives) &
                       & .and. params%do_derivatives) then
                     res%local_properties_cart_der(:, :, :) =&
                          & res%local_properties_cart_der(:, :, :) +&
                          & res%this_local_properties_cart_der(:, :, :)
                  end if

               end if
               if (params%do_forces) then
                  res%forces_soap = res%forces_soap + res%this_forces
                  res%virial_soap = res%virial_soap + res%this_virial
               end if
            end do
            n_lp_count = n_lp_count + model%soap_turbo_hypers(i)%n_local_properties

#ifdef _GPU
            call gpu_free_async(nf_d, gpu_stream)
            call gpu_free_async(rcut_hard_d, gpu_stream)
            call gpu_free_async(rcut_soft_d, gpu_stream)
            call gpu_free_async(global_scaling_d, gpu_stream)
            call gpu_free_async(atom_sigma_r_d, gpu_stream)
            call gpu_free_async(atom_sigma_r_scaling_d, gpu_stream)
            call gpu_free_async(atom_sigma_t_d, gpu_stream)
            call gpu_free_async(atom_sigma_t_scaling_d, gpu_stream)
            call gpu_free_async(amplitude_scaling_d, gpu_stream)
            call gpu_free_async(alpha_max_d, gpu_stream)
            call gpu_free_async(central_weight_d, gpu_stream)

            if (model%soap_turbo_hypers(i)%has_local_properties) then
               do j = 1, model%soap_turbo_hypers(i)%n_local_properties
                  call gpu_free_async(model%soap_turbo_hypers(i)%local_property_models(j)%alphas_d, gpu_stream)
                  call gpu_free_async(model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, gpu_stream)
               end do
            end if

            call soap_backend_end()

            call get_time(time%soap_solo(2))
#endif
            deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)
#ifdef _GPU
            time%soap_solo(3) = time%soap_solo(3) + time%soap_solo(2) - time%soap_solo(1)
#endif

            ! THIS WON'T WORK! THE SOAP AND SOAP DERIVATIVES NEED TO BE COLLECTED FROM ALL RANKS <--------------------- FIX THIS!!!!
            ! AT THE MOMENT I'M MAKING THE CODE PRINT AN ERROR MESSAGE AND STOP EXECUTION IF THE USER TRIES TO WRITE OUT THESE
            ! FILES WITH MORE THAN ONE MPI TASK
            if (rank == 0) then
               !       Write out stuff - THIS SHOULD PROBABLY BE PUT IN A MODULE
               if (model%n_soap_turbo == 1) then
                  i_char = ""
               else
                  write (i_char, '(I7)') i
                  i_char = "_"//adjustl(i_char)
               end if
               !       Write the SOAP vectors - NOT THE OPTIMAL STRATEGY IN TERMS OF DISK SPACE SINCE SOME ATOMS HAVE SOAP = 0
               if (params%write_soap) then
                  if (loop%n_xyz == 1 .or. loop%md_istep == 0) then
                     open (unit=10, file="soap"//trim(i_char)//".dat", status="unknown")
                  else
                     open (unit=10, file="soap"//trim(i_char)//".dat", status="old", position="append")
                  end if
                  if (.not. params%do_md .or. &
                      (params%do_md .and. (loop%md_istep == 0 .or. loop%md_istep == params%md_nsteps .or. &
                                           modulo(loop%md_istep, params%write_xyz) == 0))) then
                     n_sites_this = size(soap, 2)
                     n_soap = size(soap, 1)
                     write (10, *) n_sites_this, n_soap
                     do i2 = 1, n_sites_this
                        write (10, '(*(ES24.15))') soap(1:n_soap, i2)
                     end do
                  end if
                  close (10)
               end if
               if (allocated(soap)) deallocate (soap)

               !       Optionally, write out the derivatives (might take a lot of disk space)
               if ((params%do_derivatives .or. params%do_derivatives_fd) .and. params%write_derivatives) then
                  if (loop%n_xyz == 1 .or. loop%md_istep == 0) then
                     open (unit=10, file="soap_der"//trim(i_char)//".dat", status="unknown")
                  else
                     open (unit=10, file="soap_der"//trim(i_char)//".dat", status="old", position="append")
                  end if
                  if (.not. params%do_md .or. &
                      (params%do_md .and. (loop%md_istep == 0 .or. loop%md_istep == params%md_nsteps .or. &
                                           modulo(loop%md_istep, params%write_xyz) == 0))) then
                     !           Note, this n_sites is not the same as the total number of sites, it's just the total number
                     !           of sites that have a derivative, since the first neighbor of each site is itself, the site
                     !           ID can always be retrieved from there. Note also that the sites are not necessarily given in
                     !           order
                     n_sites_this = size(der_neighbors, 1)
                     n_soap = size(soap_cart_der, 2)
                     nl%n_atom_pairs = size(der_neighbors_list, 1)
                     write (10, *) state%n_sites, n_soap, nl%n_atom_pairs
                     k = 1
                     k2 = 0
                     do i2 = 1, n_sites_this
                        write (10, *) der_neighbors_list(k), der_neighbors(i2), der_neighbors_list(k:k + der_neighbors(i2) - 1)
                        k = k + der_neighbors(i)
                        do j = 1, der_neighbors(i)
                           k2 = k2 + 1
                           write (10, '(*(ES24.15))') soap_cart_der(1, 1:n_soap, k2)
                           write (10, '(*(ES24.15))') soap_cart_der(2, 1:n_soap, k2)
                           write (10, '(*(ES24.15))') soap_cart_der(3, 1:n_soap, k2)
                        end do
                     end do
                  end if
                  close (10)
               end if
               if (params%write_derivatives) then
                  deallocate (soap_cart_der, der_neighbors, der_neighbors_list)
               end if
            end if

            call time_end(time%soap)

         end do
         call time_end(time%gap)

         if (any_has_local_properties(model%soap_turbo_hypers)) then
            call time_start(time%mpi)
            call comm_sum_to_root(comm, res%local_properties, res%this_local_properties, state%n_sites*params%n_local_properties)
            res%local_properties = res%this_local_properties
            call comm_bcast(comm, res%local_properties, state%n_sites*params%n_local_properties)

            call time_end(time%mpi)
         end if

!        Each rank owns a slice of the sites, so its local_dipoles is zero
!        everywhere else and a plain sum is the whole reduction.
         if (params%do_dipole) then
            call time_start(time%mpi)
            call comm_sum_to_root(comm, res%local_dipoles, res%this_local_dipoles, 3*state%n_sites)
            res%local_dipoles = res%this_local_dipoles
            call comm_bcast(comm, res%local_dipoles, 3*state%n_sites)

            call comm_sum_to_root(comm, res%energies_dipole, res%this_energies_dipole, state%n_sites)
            res%energies_dipole = res%this_energies_dipole
            call comm_bcast(comm, res%energies_dipole, state%n_sites)
            call time_end(time%mpi)
         end if

         if (params%do_dipole) then
            res%dipole(1) = sum(res%local_dipoles(1, 1:state%n_sites))
            res%dipole(2) = sum(res%local_dipoles(2, 1:state%n_sites))
            res%dipole(3) = sum(res%local_dipoles(3, 1:state%n_sites))
         end if

!        IR PREDICTION FROM A TRAJECTORY. This frame's total dipole joins the
!        ensemble, with the time its comment line claimed. Nothing is
!        transformed yet: the file's length is not known until the end of it,
!        and the resolution follows from that length.
!
!        Every rank keeps the same buffer. local_dipoles was all-reduced and
!        broadcast just above, so the sums agree bit for bit, and having the
!        ensemble replicated means the final transform needs no communication.
!        Three doubles a frame; a 100 ps trajectory at 1 fs is 2.4 MB.
         if (ir_from_traj) then
            call ir_fft_frames_push(ir_fft_frames, res%dipole, state%frame_time, state%has_frame_time)
         end if

         !     Compute vdW energies and forces

!        Compute ELECTROSTATIC energies and forces
!
!        Ported from the GPU branch. That branch additionally routes the gsf
!        method through a batched device implementation when params%gpu_batched
!        is set; here gsf always takes the compute_coulomb_lamichhane path,
!        which is what the GPU branch itself falls back to.
!        valid_estat_charges is part of the guard, not an afterthought: without it a
!        deck that asks for electrostatics against a GAP with no atomic_charge local
!        property indexes local_properties with an uninitialised charge_lp_index and
!        segfaults. Same shape as the has_vdw/has_local_properties defect.
!        Moved to src/turbogap_estat.f90. The #ifdef is here, at the one call,
!        rather than inside three continued argument lists where nothing
!        Fortran-aware could parse it.
#ifdef _GPU
#ifdef _MPIF90
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, state%species, nl%neighbor_species, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, n_omp, &
                            res%this_energies_estat, res%this_forces_estat, res%this_virial_estat, time)
#else
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, state%species, nl%neighbor_species, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, n_omp, &
                            res%energies_estat, res%forces_estat, res%virial_estat, time)
#endif
#else
#ifdef _MPIF90
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, &
                            res%this_energies_estat, res%this_forces_estat, res%this_virial_estat, time)
#else
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, &
                            res%energies_estat, res%forces_estat, res%virial_estat, time)
#endif
#endif

         call compute_vdw(params, any_has_vdw(model%soap_turbo_hypers), state%n_sites, &
                          nl%n_neigh, nl%neighbors_list, nl%neighbor_species, nl%rjs, nl%xyz, &
                          res%local_properties, res%local_properties_cart_der, model%vdw_lp_index, &
                          dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, dom%n_atom_pairs_by_rank, dom%site_in_rank, &
                          state%indices, rank, ntasks, loop%md_istep, vdw_ws, &
                          res%energies_vdw, res%forces_vdw, res%virial_vdw, res%local_virial_vdw_diag, &
                          res%this_energies_vdw, res%this_forces_vdw, res%this_virial_vdw, &
                          res%this_local_virial_vdw_diag, res%energies_vdw_corr, res%forces_vdw_corr, &
                          res%local_virial_vdw_diag_corr, res%mbd_ts_scaling, res%this_mbd_ts_scaling, &
                          res%update_mbd_ts_scaling, time)

         !--- EXPERIMENTAL SPECTRUM CALCULATION AND FORCES ---!

         ! --- Changing the implementation:
         !     > All experimental prediction should be done here
         !     > do_exp is the variable which says whether calculation should be done
         !     > experimental_forces = .true. will add forces to the calculation

         !###---   Compute Experimental Data Interpolation   ---###!

         if (params%do_exp) then
            do i = 1, params%n_exp
               ! If we want to compute the experimental interpolation, we do it now.

               call get_write_condition(params%do_mc, params%do_md&
                    &, loop%mc_istep, loop%md_istep, params%write_xyz,&
                    & write_condition)

               if (params%exp_data(i)%compute_exp) then
                  if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
                  if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
                  call calculate_exp_interpolation(params%exp_data(i)&
                       &%x, params%exp_data(i)%y, params%exp_data(i)&
                       &%n_samples, params%exp_data(i)%data)

!                 The weights live on the same grid as the experiment, so they
!                 are built here, from the same x, every time it is rebuilt.
                  call build_exp_weights(params%exp_data(i)%x, params%exp_data(i)%w, &
                                         params%exp_data(i)%weights_data, &
                                         params%exp_data(i)%n_weights, &
                                         params%exp_data(i)%data, &
                                         params%exp_data(i)%data_weights, &
                                         params%exp_data(i)%n_data_weights, &
                                         params%exp_data(i)%n_data, &
                                         trim(params%exp_data(i)%file_data_weights))

                  call preprocess_exp_data(params, params%exp_data(i)%x,&
                       & params%exp_data(i)%y, params%exp_data(i)%label,&
                       & state%n_sites, dot_product(cross_product(state%a_box,&
                       & state%b_box), state%c_box)/(dfloat(state%indices(1)*state%indices(2) &
                       &*state%indices(3))), params%exp_data(i)%input, exp_output, .true.)

                  if (params%write_exp .and. .not. params&
                       &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition) then

                     call get_overwrite_condition(params%do_mc,&
                          & params%do_md, loop%mc_istep, loop%md_istep, params&
                          &%write_xyz, overwrite_condition)

                     call write_exp_data(params%exp_data(i)%x, params&
                          &%exp_data(i)%y, overwrite_condition,&
                          & trim(params%exp_data(i)%label)//&
                          & "_exp.dat", params%exp_data(i)%label)
                  end if

               end if

               if (params%exp_data(i)%compute_exp .and. .not. params&
                    &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition) then

                  if (params%write_exp) then
                     write (filename, '(A)')&
                          & trim(params%exp_data(i)%label)//"_exp_fit.dat"

                     call get_overwrite_condition(params%do_mc,&
                          & params%do_md, loop%mc_istep, loop%md_istep, params&
                          &%write_xyz, overwrite_condition)

                     call write_exp_data(params%exp_data(i)%x, params%exp_data(i)%y,&
                          & overwrite_condition, trim(params&
                          &%exp_data(i)%label)//"_exp_fit.dat",&
                          & trim(params%exp_data(i)%label)//" : output = "&
                          & //trim(exp_output))

                  end if

               end if

               params%exp_data(i)%wrote_exp = .true.
               params%exp_data(i)%compute_exp = .true.

            end do
         end if

         !###---   XPS Forces and Spectra Prediction   ---###!

         !     Compute core_electron_be energies and forces
         !
         ! Moved to src/turbogap_exp.f90. The #ifdef is here, at the one call,
         ! rather than inside a continued argument list where nothing
         ! Fortran-aware could parse it.
#ifdef _MPIF90
         call compute_exp_xps(params, state%n_sites, loop%n_xyz, nl%xyz, nl%neighbors_list, nl%n_neigh, &
                              res%local_properties, res%local_properties_cart_der, model%soap_turbo_hypers, &
                              state%a_box, state%b_box, state%c_box, state%indices, dom%i_beg, dom%i_end, dom%j_beg, &
                              dom%j_end, rank, &
                              loop%md_istep, loop%mc_istep, model%valid_xps, model%xps_idx, model%core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              res%this_energies_lp, res%this_forces_lp, res%this_virial_lp, time)
#else
         call compute_exp_xps(params, state%n_sites, loop%n_xyz, nl%xyz, nl%neighbors_list, nl%n_neigh, &
                              res%local_properties, res%local_properties_cart_der, model%soap_turbo_hypers, &
                              state%a_box, state%b_box, state%c_box, state%indices, dom%i_beg, dom%i_end, dom%j_beg, &
                              dom%j_end, rank, &
                              loop%md_istep, loop%mc_istep, model%valid_xps, model%xps_idx, model%core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              res%energies_lp, res%forces_lp, res%virial_lp, time)
#endif

         !###---   (Partial) Pair distribution functions and XRD   ---###!
         !
         ! Moved to src/turbogap_exp.f90. The #ifdef below is the whole reason
         ! it is here rather than inside: the exp_interface routines take the
         ! this_-prefixed arrays under MPI and the plain ones otherwise. Choosing
         ! once, at the call, is what let four preprocessor-interrupted argument
         ! lists disappear from the moved code.
#ifdef _GPU
#ifdef _MPIF90
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, &
                                  dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, res%this_energies_pdf, &
                                  res%this_forces_pdf, res%this_virial_pdf, res%this_energies_sf, &
                                  res%this_forces_sf, res%this_virial_sf, res%this_energies_xrd, res%this_forces_xrd, &
                                  res%this_virial_xrd, res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, time, &
                                  i_beg_list, i_end_list, j_beg_list, &
                                  j_end_list, n_omp, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                                  n_sites_temp, n_pairs_temp, write_condition, overwrite_condition, &
                                  temp_string, species_types_actual, state%v_uc)
#else
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, &
                                  dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, res%energies_pdf, res%forces_pdf, &
                                  res%virial_pdf, res%energies_sf, res%forces_sf, &
                                  res%virial_sf, res%energies_xrd, res%forces_xrd, res%virial_xrd, res%energies_nd, res%forces_nd, &
                                  res%virial_nd, time, i_beg_list, &
                                  i_end_list, j_beg_list, j_end_list, n_omp, omp_task, this_i_beg, this_i_end, &
                                  this_j_beg, this_j_end, n_sites_temp, n_pairs_temp, write_condition, &
                                  overwrite_condition, temp_string, species_types_actual, state%v_uc)
#endif
#else
#ifdef _MPIF90
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, &
                                  res%this_energies_pdf, res%this_forces_pdf, res%this_virial_pdf, &
                                  res%this_energies_sf, res%this_forces_sf, res%this_virial_sf, &
                                  res%this_energies_xrd, res%this_forces_xrd, res%this_virial_xrd, &
                                  res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, &
                                  time)
#else
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, &
                                  res%energies_pdf, res%forces_pdf, res%virial_pdf, &
                                  res%energies_sf, res%forces_sf, res%virial_sf, &
                                  res%energies_xrd, res%forces_xrd, res%virial_xrd, &
                                  res%energies_nd, res%forces_nd, res%virial_nd, &
                                  time)
#endif
#endif

         if (params%do_prediction) then
            !       Two-body, core-potential and three-body contributions, via the
            !       gap_backend seam. The CPU implementation is in
            !       src/gap_backend_cpu.f90; the GPU branch provides the same three
            !       names from src/gap_backend_gpu.f90 and the Makefile picks one.
            call time_start(time%gap)

            call gap_backend_begin(params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                   nl%neighbors_list, dom%i_beg, dom%i_end, dom%j_beg, dom%j_end)

            call add_2b_contribution(model%n_distance_2b, model%distance_2b_hypers, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%energies_2b, res%forces_2b, res%virial_2b, time)

            call add_core_pot_contribution(model%n_core_pot, model%core_pot_hypers, &
                                           params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                           dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                           res%this_virial, &
                                           res%energies_core_pot, res%forces_core_pot, res%virial_core_pot, time)

#ifdef _GPU
            call add_3b_contribution(model%n_angle_3b, model%angle_3b_hypers, nl%neighbors_list, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%forces, res%energies_3b, res%forces_3b, res%virial_3b, time)
#else
            call add_3b_contribution(model%n_angle_3b, model%angle_3b_hypers, nl%neighbors_list, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%energies_3b, res%forces_3b, res%virial_3b, time)
#endif

            call gap_backend_end()

            call time_end(time%gap)
            !       Communicate all energies and forces here for all
            !       terms
#ifdef _MPIF90
            call time_start(time%mpi_ef)
!       One evaluation of the eleven predicates, and one list built from them.
!       The pack and unpack walks below read only that list, so they cannot
!       disagree about which slot belongs to which family -- the failure mode
!       this replaces was three independent copies of these conditions, where
!       any two disagreeing shifts the slot numbering and silently attributes
!       one family's energies and forces to another.
            contrib_on(C_SOAP) = (model%n_soap_turbo > 0)
            contrib_on(C_VDW) = allocated(res%this_energies_vdw)
            contrib_on(C_ESTAT) = allocated(res%this_energies_estat)
            contrib_on(C_LP) = allocated(res%this_energies_lp)
            contrib_on(C_PDF) = allocated(res%this_energies_pdf) .and. params%valid_pdf
            contrib_on(C_SF) = allocated(res%this_energies_sf) .and. params%valid_sf
            contrib_on(C_XRD) = allocated(res%this_energies_xrd) .and. params%valid_xrd
            contrib_on(C_ND) = allocated(res%this_energies_nd) .and. params%valid_nd
            contrib_on(C_2B) = (model%n_distance_2b > 0)
            contrib_on(C_CP) = (model%n_core_pot > 0)
            contrib_on(C_3B) = (model%n_angle_3b > 0)

            n_active = 0
            if (contrib_on(C_SOAP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%energies_soap
               contrib(n_active)%e_dst => res%energies_soap
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%forces_soap
                  contrib(n_active)%v_src => res%virial_soap
                  contrib(n_active)%f_dst => res%forces_soap
                  contrib(n_active)%v_dst => res%virial_soap
               end if
            end if
            if (contrib_on(C_VDW)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_vdw
               contrib(n_active)%e_dst => res%energies_vdw
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_vdw
                  contrib(n_active)%v_src => res%this_virial_vdw
                  contrib(n_active)%f_dst => res%forces_vdw
                  contrib(n_active)%v_dst => res%virial_vdw
               end if
            end if
            if (contrib_on(C_ESTAT)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_estat
               contrib(n_active)%e_dst => res%energies_estat
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_estat
                  contrib(n_active)%v_src => res%this_virial_estat
                  contrib(n_active)%f_dst => res%forces_estat
                  contrib(n_active)%v_dst => res%virial_estat
               end if
            end if
            if (contrib_on(C_LP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_lp
               contrib(n_active)%e_dst => res%energies_lp
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_lp
                  contrib(n_active)%v_src => res%this_virial_lp
                  contrib(n_active)%f_dst => res%forces_lp
                  contrib(n_active)%v_dst => res%virial_lp
               end if
            end if
            if (contrib_on(C_PDF)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_pdf
               contrib(n_active)%e_dst => res%energies_pdf
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_pdf
                  contrib(n_active)%v_src => res%this_virial_pdf
                  contrib(n_active)%f_dst => res%forces_pdf
                  contrib(n_active)%v_dst => res%virial_pdf
               end if
            end if
            if (contrib_on(C_SF)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_sf
               contrib(n_active)%e_dst => res%energies_sf
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_sf
                  contrib(n_active)%v_src => res%this_virial_sf
                  contrib(n_active)%f_dst => res%forces_sf
                  contrib(n_active)%v_dst => res%virial_sf
               end if
            end if
            if (contrib_on(C_XRD)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_xrd
               contrib(n_active)%e_dst => res%energies_xrd
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_xrd
                  contrib(n_active)%v_src => res%this_virial_xrd
                  contrib(n_active)%f_dst => res%forces_xrd
                  contrib(n_active)%v_dst => res%virial_xrd
               end if
            end if
            if (contrib_on(C_ND)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%this_energies_nd
               contrib(n_active)%e_dst => res%energies_nd
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%this_forces_nd
                  contrib(n_active)%v_src => res%this_virial_nd
                  contrib(n_active)%f_dst => res%forces_nd
                  contrib(n_active)%v_dst => res%virial_nd
               end if
            end if
            if (contrib_on(C_2B)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%energies_2b
               contrib(n_active)%e_dst => res%energies_2b
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%forces_2b
                  contrib(n_active)%v_src => res%virial_2b
                  contrib(n_active)%f_dst => res%forces_2b
                  contrib(n_active)%v_dst => res%virial_2b
               end if
            end if
            if (contrib_on(C_CP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%energies_core_pot
               contrib(n_active)%e_dst => res%energies_core_pot
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%forces_core_pot
                  contrib(n_active)%v_src => res%virial_core_pot
                  contrib(n_active)%f_dst => res%forces_core_pot
                  contrib(n_active)%v_dst => res%virial_core_pot
               end if
            end if
            if (contrib_on(C_3B)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => res%energies_3b
               contrib(n_active)%e_dst => res%energies_3b
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => res%forces_3b
                  contrib(n_active)%v_src => res%virial_3b
                  contrib(n_active)%f_dst => res%forces_3b
                  contrib(n_active)%v_dst => res%virial_3b
               end if
            end if

            counter2 = n_active

            allocate (all_energies(1:state%n_sites, 1:counter2))
            allocate (all_this_energies(1:state%n_sites, 1:counter2))
            if (params%do_forces) then
               allocate (all_forces(1:3, 1:state%n_sites, 1:counter2))
               allocate (all_this_forces(1:3, 1:state%n_sites, 1:counter2))
               allocate (all_virial(1:3, 1:3, 1:counter2))
               allocate (all_this_virial(1:3, 1:3, 1:counter2))
            end if

!       Pack. A family owns a slot whenever it is active, but only contributes
!       forces when it carries them -- the exp-spectra families additionally
!       need exp_forces. Their slot must still be cleared: all_forces is
!       allocated and never zeroed, and mpi_reduce below reads the whole array
!       regardless of who wrote what into it.
            do i_contrib = 1, n_active
               all_energies(1:state%n_sites, i_contrib) = contrib(i_contrib)%e_src(1:state%n_sites)
               if (contrib(i_contrib)%forces) then
                  all_forces(1:3, 1:state%n_sites, i_contrib) = contrib(i_contrib)%f_src(1:3, 1:state%n_sites)
                  all_virial(1:3, 1:3, i_contrib) = contrib(i_contrib)%v_src(1:3, 1:3)
               else if (params%do_forces) then
                  all_forces(1:3, 1:state%n_sites, i_contrib) = 0.d0
                  all_virial(1:3, 1:3, i_contrib) = 0.d0
               end if
            end do

            !       Here we communicate
            call comm_sum_to_root(comm, all_energies, all_this_energies, state%n_sites*counter2)
            if (params%do_forces) then
               call comm_sum_to_root(comm, all_forces, all_this_forces, 3*state%n_sites*counter2)
               call comm_sum_to_root(comm, all_virial, all_this_virial, 9*counter2)
            end if

!       Unpack. For the six families packed from a this_ array this is where
!       the reduced result lands in the un-prefixed one.
            do i_contrib = 1, n_active
               contrib(i_contrib)%e_dst(1:state%n_sites) = all_this_energies(1:state%n_sites, i_contrib)
               if (contrib(i_contrib)%forces) then
                  contrib(i_contrib)%f_dst(1:3, 1:state%n_sites) = all_this_forces(1:3, 1:state%n_sites, i_contrib)
                  contrib(i_contrib)%v_dst(1:3, 1:3) = all_this_virial(1:3, 1:3, i_contrib)
               end if
            end do

!       Release the this_ arrays now that their contents have been unpacked.
!       Kept explicit rather than folded into the loop: an allocatable cannot
!       be deallocated through a pointer, and this_local_virial_vdw_diag has no
!       counterpart in the list.
            if (contrib_on(C_VDW)) then
               deallocate (res%this_energies_vdw)
               if (params%do_forces) deallocate (res%this_forces_vdw, res%this_local_virial_vdw_diag)
            end if
            if (contrib_on(C_ESTAT)) then
               deallocate (res%this_energies_estat)
               if (params%do_forces) deallocate (res%this_forces_estat)
            end if
            if (contrib_on(C_LP)) then
               deallocate (res%this_energies_lp)
               if (params%do_forces) deallocate (res%this_forces_lp)
            end if
            if (contrib_on(C_PDF)) then
               deallocate (res%this_energies_pdf)
               if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_pdf)
            end if
            if (contrib_on(C_SF)) then
               deallocate (res%this_energies_sf)
               if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_sf)
            end if
            if (contrib_on(C_XRD)) then
               deallocate (res%this_energies_xrd)
               if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_xrd)
            end if
            if (contrib_on(C_ND)) then
               deallocate (res%this_energies_nd)
               if (params%do_forces .and. params%exp_forces) deallocate (res%this_forces_nd)
            end if

            !       Clean up
            deallocate (all_energies, all_this_energies)
            if (params%do_forces) then
               deallocate (all_forces, all_this_forces, all_virial, all_this_virial)
            end if

            call time_end(time%mpi_ef)
#endif

            !       Add up all the energy terms
            res%energies = res%energies + res%energies_soap + res%energies_2b +&
                 & res%energies_3b + res%energies_core_pot + res%energies_vdw + res%energies_estat !+energies_lp

            if (model%valid_xps) res%energies_exp = res%energies_exp + res%energies_lp
            if (perform%pdf) res%energies_exp = res%energies_exp + res%energies_pdf
            if (perform%sf) res%energies_exp = res%energies_exp + res%energies_sf
            if (perform%xrd) res%energies_exp = res%energies_exp + res%energies_xrd
            if (perform%nd) res%energies_exp = res%energies_exp + res%energies_nd

            if (params%exp_energies) res%energies = res%energies + res%energies_exp

            res%energy_prev = res%energy
            instant_pressure_prev = instant_pressure
            res%energy = sum(res%energies)
            res%energy_exp = sum(res%energies_exp)

         end if

         if (.not. params%do_md .and. .not. params%do_mc) then
            if (rank == 0) then
               write (*, *) '                                       |'
               write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
               write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
               write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
               write (*, '(A,1X,F21.8,1X,A)') ' estat energy:', sum(res%energies_estat), 'eV |'
               write (*, '(A,1X,F22.8,1X,A)') ' Exp. energy:', sum(res%energies_exp), 'eV |'
               if (model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(res%energies_lp), 'eV |'
               if (perform%pdf)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                    & sum(res%energies_pdf), 'eV |'
               if (perform%sf)&
                    & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                    & sum(res%energies_sf), 'eV |'
               if (perform%xrd)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                    & sum(res%energies_xrd), 'eV |'
               if (perform%nd)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                    & sum(res%energies_nd), 'eV |'

               if (.not. params%do_mc .or. (params%do_mc .and. loop%mc_istep <= 1)) then
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(res%energies), 'eV |'
               else
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(images(i_trial_image)%energies), 'eV |'
               end if

               if (.not. params%do_mc) then
                  write (*, *) '                                       |'
                  write (*, *) 'Energy & forces in "trajectory_out.xyz"|'
                  write (*, *) '                                       |'
                  write (*, *) '.......................................|'
               else if (loop%mc_istep == 0) then
                  write (*, *) '                                       |'
                  write (*, *) ' MC configs in "mc_current.xyz" and    |'
                  write (*, *) '               "mc_trial.xyz"          |'
                  write (*, *) '               "mc_all.xyz"            |'
                  write (*, *) '.......................................|'
               end if
            end if
         end if

         if (params%do_forces) then
            res%forces = res%forces_soap + res%forces_2b + res%forces_3b + res%forces_core_pot + res%forces_vdw
            res%virial = res%virial_soap + res%virial_2b + res%virial_3b + res%virial_core_pot + res%virial_vdw

!           MAD IR bias. The dipole of this configuration joins the ensemble,
!           the spectrum is compared with the experiment, and the gradient of
!           the mismatch with respect to THIS configuration is added to the
!           forces. Nothing is applied until the ensemble is full, because a
!           partly filled one has a resolution that changes step to step.
!
!           No virial: the bias is a function of the dipole, not of the cell,
!           and a stress from it would be wrong rather than merely missing.
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
                                        params%exp_energy_scales_final(params%ir_idx), mad_ir_scale)
                  call time_start(time%ir_predict)
!                 Advance BEFORE evaluate: the loss and the gradient both belong
!                 to a bank that already knows this frame's dipole. Evaluating
!                 first would report the previous frame's spectrum and, worse,
!                 would make the weight the gradient of a loss the current
!                 dipole had not yet entered -- which is identically zero, not
!                 merely inaccurate.
                  call mad_ir_xl_advance(mad_ir_xl_state, res%local_dipoles(1:3, 1:state%n_sites))
                  call mad_ir_xl_evaluate(mad_ir_xl_state, mad_ir_scale, mad_ir_energy)
                  call time_end(time%ir_predict)
                  res%energies_exp = res%energies_exp + mad_ir_energy/dfloat(state%n_sites)
                  exp_dissimilarity = exp_dissimilarity + mad_ir_xl_state%dissim
                  exp_dissim_ref = exp_dissim_ref + mad_ir_xl_state%dissim_ref
                  if (params%exp_energies) then
                     res%energies = res%energies + mad_ir_energy/dfloat(state%n_sites)
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
                  if (mad_ir_xl_ready(mad_ir_xl_state) .and. .not. mad_ir_applied) then
                     call get_time(mad_ir_t_now)
                     mad_ir_t_first = mad_ir_t_now - time3
                     mad_ir_step_first = loop%md_istep
                  end if
                  mad_ir_applied = mad_ir_xl_ready(mad_ir_xl_state)
                  if (.not. mad_ir_applied) mad_ir_energy = 0.d0
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
                     ir_aux_escale_top = max(params%exp_energy_scales_initial(params%ir_idx), &
                                             params%exp_energy_scales_final(params%ir_idx))
                     ir_aux_temp = min(params%t_beg, params%t_end)
                     if (ir_aux_temp <= 0.d0) ir_aux_temp = max(params%t_beg, params%t_end)
                     call ir_aux_calibrate(ir_aux_state, mad_ir_state, ir_aux_temp, &
                                           ir_aux_escale_top, ir_aux_ok, ir_aux_msg)
                     call time_end(time%ir_predict)
                     if (rank == 0) then
                        if (ir_aux_ok) then
                           write (*, *) 'MAD IR: ', trim(ir_aux_msg)
                        else
                           write (*, *) 'WARNING: ', trim(ir_aux_msg)
                        end if
                     end if
                     if (.not. ir_aux_ok) then
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
                        if (rank == 0) then
                           write (*, *) "ERROR: ir_bias_mode = aux is above its stability threshold."
                           write (*, '(A,ES12.4)') "        stability number Lambda = ", ir_aux_state%stab
                           write (*, *) "        Lambda must be below 1; the bilinear coupling runs away above it."
                           write (*, '(A,ES12.4)') "        largest usable exp_energy_scales = ", &
                              ir_aux_escale_max(ir_aux_state)
                           write (*, '(A,ES12.4)') "        this run asked for               = ", ir_aux_escale_top
                           write (*, *) "        Lower exp_energy_scales, or raise ir_aux_damping"
                           write (*, *) "        (which lowers every g_k), and try again."
                        end if
                        call turbogap_abort()
                     else if (ir_aux_state%stab >= 0.5d0 .and. rank == 0) then
                        write (*, '(A,ES10.2,A)') " WARNING: ir_aux stability number ", &
                           ir_aux_state%stab, " is above 0.5; the bias is close to runaway."
                     end if
                     call get_time(mad_ir_t_now)
                     mad_ir_t_first = mad_ir_t_now - time3
                     mad_ir_step_first = loop%md_istep
                  end if
                  if (ir_aux_calibrated(ir_aux_state)) then
                     call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                           loop%mc_istep, params%mc_nsteps, &
                                           params%exp_energy_scales_initial(params%ir_idx), &
                                           params%exp_energy_scales_final(params%ir_idx), mad_ir_scale)
                     call time_start(time%ir_predict)
!                    Advance before evaluate, for the same reason the extended
!                    Lagrangian does: the energy and the force both belong to a
!                    bank that already knows this frame's dipole.
                     call ir_aux_advance(ir_aux_state, res%dipole)
                     call ir_aux_evaluate(ir_aux_state, mad_ir_scale, res%dipole, mad_ir_energy)
                     call time_end(time%ir_predict)
                     res%energies_exp = res%energies_exp + mad_ir_energy/dfloat(state%n_sites)
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
                     call mad_ir_evaluate(mad_ir_state, 0.d0, ir_aux_acf_energy, mad_ir_lambda)
                     call time_end(time%ir_predict)
                     exp_dissimilarity = exp_dissimilarity + mad_ir_state%dissim
                     exp_dissim_ref = exp_dissim_ref + mad_ir_state%dissim_ref
                     if (params%exp_energies) then
                        res%energies = res%energies + mad_ir_energy/dfloat(state%n_sites)
                        res%energy = sum(res%energies)
                     end if
                     res%energy_exp = sum(res%energies_exp)
                     if (params%exp_forces) then
                        call time_start(time%ir_forces)
                        call ir_aux_forces(ir_aux_state, mad_ir_scale, mad_ir_dmu_dr, res%forces, &
                                           state%velocities(1:3, 1:state%n_sites), ir_aux_power)
                        call time_end(time%ir_forces)
!                       Integrated with the stored-frame interval, since that is
!                       how often the force is refreshed.
                        ir_aux_work = ir_aux_work &
                                      + ir_aux_power*params%md_step*dfloat(params%ir_stride)
                     end if
                     mad_ir_applied = .true.
!                    ---- the thermal-fidelity check -------------------------
!                    Over a window of stored frames, how much energy did the
!                    controller put in, and what standing temperature offset
!                    does the thermostat therefore have to hold against?
!                    50 fs of biased dynamics is enough to average the pump
!                    rate over many resonator periods (the fastest fitted band
!                    is ~8 fs) while still reporting inside a short run.
                     if (md_time - ir_aux_pump_time > 50.d0) then
                        if (ir_aux_pump_time > 0.d0 .and. params%tau_t > 0.d0) then
!                          The work the BIAS FORCE did on the atoms, which is
!                          the channel that actually heats: the controller's own
!                          injection into the bank is a different and, for a bank
!                          far off target, much smaller number.
                           ir_aux_dT = 2.d0*(ir_aux_work - ir_aux_pump_prev) &
                                       /(md_time - ir_aux_pump_time)*params%tau_t &
                                       /(3.d0*dfloat(state%n_sites)*8.6173303d-5)
                           if (rank == 0 .and. .not. ir_aux_warned_hot .and. &
                               dabs(ir_aux_dT) > 0.1d0*max(1.d0, params%t_beg)) then
                              write (*, '(A,F10.1,A)') " WARNING: ir_aux is pumping hard enough for a standing", &
                                 ir_aux_dT, " K offset."
                              write (*, *) "          The bias is below its stability bound but above what the"
                              write (*, *) "          thermostat can absorb quietly. Lower exp_energy_scales."
                              ir_aux_warned_hot = .true.
                           end if
                        end if
                        ir_aux_pump_prev = ir_aux_work
                        ir_aux_pump_time = md_time
                     end if
                  else
                     mad_ir_energy = 0.d0
                     mad_ir_applied = .false.
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
                                        params%exp_energy_scales_final(params%ir_idx), mad_ir_scale)
                  call time_start(time%ir_predict)
                  if (allocated(ir_fft_mu_chron)) then
                     if (size(ir_fft_mu_chron, 2) /= mad_ir_state%n_stored) &
                        deallocate (ir_fft_mu_chron)
                  end if
                  if (.not. allocated(ir_fft_mu_chron)) &
                     allocate (ir_fft_mu_chron(1:3, 1:mad_ir_state%n_stored))
                  if (.not. allocated(ir_fft_I_fit)) &
                     allocate (ir_fft_I_fit(1:mad_ir_state%n_freq))
                  call ir_fft_md_unroll(mad_ir_state%mu_hist, mad_ir_state%n_window, &
                                        mad_ir_state%n_stored, mad_ir_state%head, &
                                        ir_fft_mu_chron, ir_fft_n_chron)
                  ir_fft_cfg%dt_fs = mad_ir_state%dt
                  call ir_fft_loss(ir_fft_mu_chron, ir_fft_n_chron, ir_fft_cfg, &
                                   mad_ir_state%nu, mad_ir_state%I_exp, mad_ir_state%wgt, &
                                   mad_ir_state%n_freq, params%ir_match_scale, &
                                   params%ir_match_offset, mad_ir_scale, &
                                   mad_ir_energy, mad_ir_lambda, ir_fft_I_fit, &
                                   ir_fft_scale_fit, ir_fft_offset_fit, &
                                   ir_fft_dissim, ir_fft_dissim_ref, ir_fft_ok, ir_fft_msg)
                  call time_end(time%ir_predict)
                  if (.not. ir_fft_ok) then
                     write (*, *) "ERROR: ", trim(ir_fft_msg)
                     stop
                  end if
                  res%energies_exp = res%energies_exp + mad_ir_energy/dfloat(state%n_sites)
                  exp_dissimilarity = exp_dissimilarity + ir_fft_dissim
                  exp_dissim_ref = exp_dissim_ref + ir_fft_dissim_ref
                  if (params%exp_energies) then
                     res%energies = res%energies + mad_ir_energy/dfloat(state%n_sites)
                     res%energy = sum(res%energies)
                  end if
                  res%energy_exp = sum(res%energies_exp)
                  if (params%exp_forces) then
                     call time_start(time%ir_forces)
                     call mad_ir_forces(mad_ir_lambda, mad_ir_dmu_dr, res%forces)
                     call time_end(time%ir_forces)
                  end if
                  if (.not. mad_ir_applied) then
                     call get_time(mad_ir_t_now)
                     mad_ir_t_first = mad_ir_t_now - time3
                     mad_ir_step_first = loop%md_istep
                  end if
                  mad_ir_applied = .true.
               else if (params%valid_ir .and. mad_ir_ready(mad_ir_state)) then
!                 The weight is exp_energy_scales, ramped over the run exactly
!                 as every other MAD observable's is.
                  call get_energy_scale(params%do_md, params%do_mc, loop%md_istep, params%md_nsteps, &
                                        loop%mc_istep, params%mc_nsteps, &
                                        params%exp_energy_scales_initial(params%ir_idx), &
                                        params%exp_energy_scales_final(params%ir_idx), mad_ir_scale)
                  call time_start(time%ir_predict)
                  call mad_ir_evaluate(mad_ir_state, mad_ir_scale, mad_ir_energy, mad_ir_lambda)
                  call time_end(time%ir_predict)
!                 The mismatch is an energy like the others, spread over the
!                 sites; the force is its gradient, and only if exp_forces.
                  res%energies_exp = res%energies_exp + mad_ir_energy/dfloat(state%n_sites)
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
                     res%energies = res%energies + mad_ir_energy/dfloat(state%n_sites)
                     res%energy = sum(res%energies)
                  end if
                  res%energy_exp = sum(res%energies_exp)
                  if (params%exp_forces) then
                     call time_start(time%ir_forces)
                     call mad_ir_forces(mad_ir_lambda, mad_ir_dmu_dr, res%forces)
                     call time_end(time%ir_forces)
                  end if
                  if (.not. mad_ir_applied) then
!                    The first spectrum of the run: how long the ensemble took
!                    to fill, measured from the same origin as the total.
                     call get_time(mad_ir_t_now)
                     mad_ir_t_first = mad_ir_t_now - time3
                     mad_ir_step_first = loop%md_istep
                  end if
                  mad_ir_applied = .true.
               else
                  mad_ir_energy = 0.d0
                  mad_ir_applied = .false.
               end if
!              Persist the ensemble alongside the trajectory. Losing it costs
!              n_window samples of unbiased dynamics on the next restart.
               call time_start(time%ir_io)
               if (rank == 0 .and. params%write_xyz > 0 .and. params%valid_ir) then
                  if (modulo(loop%md_istep, params%write_xyz) == 0 .or. loop%md_istep == params%md_nsteps) then
                     if (trim(params%ir_restart_file) /= "none") then
                        call mad_ir_save(mad_ir_state, params%ir_restart_file, mad_ir_ok, mad_ir_msg)
                        if (.not. mad_ir_ok) write (*, *) "WARNING: ", trim(mad_ir_msg)
                     end if
!                    The bank, in its own file. Losing it costs the warm-up
!                    again, and it is the larger of the two by orders of
!                    magnitude, which is why it is a separate write that can be
!                    turned off on its own.
                     if (mad_ir_xl_active .and. trim(params%ir_xl_restart_file) /= "none") then
                        call mad_ir_xl_save(mad_ir_xl_state, params%ir_xl_restart_file, &
                                            mad_ir_xl_ok, mad_ir_xl_msg)
                        if (.not. mad_ir_xl_ok) write (*, *) "WARNING: ", trim(mad_ir_xl_msg)
                     end if
!                    And the envelope-targeted bank. X, P and eta are state
!                    in the same sense the velocities are, and eta especially:
!                    it is an integrator, so dropping it throws away everything
!                    the controller had learned about the mismatch.
                     if (ir_aux_active .and. ir_aux_calibrated(ir_aux_state) .and. &
                         trim(params%ir_aux_restart_file) /= "none") then
                        call ir_aux_save(ir_aux_state, params%ir_aux_restart_file, &
                                         ir_aux_ok, ir_aux_msg)
                        if (.not. ir_aux_ok) write (*, *) "WARNING: ", trim(ir_aux_msg)
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
                  mad_ir_have_spectrum = mad_ir_state%n_stored > 1
!                 The last COLLECTED step, not the last step: with an
!                 ir_stride that does not divide md_nsteps the two differ, and
!                 the final frame is the one the whole run was for.
                  write_condition = write_condition .or. &
                                    (loop%md_istep > params%md_nsteps - params%ir_stride)
               else
                  mad_ir_have_spectrum = mad_ir_applied
               end if
               if (rank == 0 .and. params%write_ir .and. mad_ir_have_spectrum &
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
                  if (ir_fft_active .and. .not. ir_from_traj &
                      .and. mad_ir_state%n_stored > 1) then
                     if (allocated(ir_fft_mu_chron)) then
                        if (size(ir_fft_mu_chron, 2) /= mad_ir_state%n_stored) &
                           deallocate (ir_fft_mu_chron)
                     end if
                     if (.not. allocated(ir_fft_mu_chron)) &
                        allocate (ir_fft_mu_chron(1:3, 1:mad_ir_state%n_stored))
                     call ir_fft_md_unroll(mad_ir_state%mu_hist, mad_ir_state%n_window, &
                                           mad_ir_state%n_stored, mad_ir_state%head, &
                                           ir_fft_mu_chron, ir_fft_n_chron)
                     ir_fft_cfg%dt_fs = mad_ir_state%dt
                     call ir_fft_spectrum(ir_fft_mu_chron, ir_fft_n_chron, ir_fft_cfg, &
                                          ir_fft_res, ir_fft_ok, ir_fft_msg)
                     if (ir_fft_ok) then
                        call ir_fft_write_spectrum(ir_fft_res, ir_fft_cfg, &
                                                   "ir_fft_spectrum.dat", &
                                                   mad_ir_state%nu, mad_ir_state%I_exp, &
                                                   mad_ir_state%n_freq, params%valid_ir, &
                                                   ir_fft_scale_fit, ir_fft_offset_fit, &
                                                   ir_fft_dissim, ir_fft_dissim_ref, &
                                                   "from the rolling MD ensemble")
                        call ir_fft_free(ir_fft_res)
                     else
                        write (*, *) "WARNING: ir_fft_spectrum.dat not written: ", &
                           trim(ir_fft_msg)
                     end if
                  end if
                  call mad_ir_write_spectrum(mad_ir_state, "ir_spectrum.dat", &
                                             params%valid_ir, loop%md_istep, params%md_step)
                  call mad_ir_append_spectrum(mad_ir_state, "ir_prediction.dat", &
                                              .not. mad_ir_wrote_prediction, loop%md_istep, &
                                              dfloat(loop%md_istep)*params%md_step)
                  if (.not. mad_ir_wrote_prediction) then
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
                     mad_ir_wrote_prediction = .true.
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

            if (model%valid_estat_charges) res%forces = res%forces + res%forces_estat
            if (model%valid_estat_charges) res%virial = res%virial + res%virial_estat

            if (perform%xps_forces) res%forces = res%forces + res%forces_lp
            if (perform%xps_forces) res%virial = res%virial + res%virial_lp

            if (perform%pdf_forces) res%forces = res%forces + res%forces_pdf
            if (perform%pdf_forces) res%virial = res%virial + res%virial_pdf

            if (perform%sf_forces) res%forces = res%forces + res%forces_sf
            if (perform%sf_forces) res%virial = res%virial + res%virial_sf

            if (perform%xrd_forces) res%forces = res%forces + res%forces_xrd
            if (perform%xrd_forces) res%virial = res%virial + res%virial_xrd

            if (perform%nd_forces) res%forces = res%forces + res%forces_nd
            if (perform%nd_forces) res%virial = res%virial + res%virial_nd

            if (rank == 0 .and. params%print_vdw_forces) then
               print *, "> Virial ESTAT "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_estat(i, j)
                  end do
               end do

               print *, "> Virial soap "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_soap(i, j)
                  end do
               end do

               print *, "> Virial 2b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_2b(i, j)
                  end do
               end do

               print *, "> Virial 3b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_3b(i, j)
                  end do
               end do

               print *, "> Virial core_pot "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_core_pot(i, j)
                  end do
               end do

               if (perform%xrd_forces) then
                  print *, "> Virial xrd "
                  do i = 1, 3
                     do j = 1, 3
                        print *, " i, ", i, " j ", j, " ", res%virial_xrd(i, j)
                     end do
                  end do
                  temp_string = ""
                  temp_string2 = ""
                  write (temp_string, "(I8)") loop%md_istep
                  write (temp_string2, "(A)") "forces_xrd_"//trim(adjustl(temp_string))
                  open (unit=90, file=temp_string2, status="unknown")
                  do i = 1, state%n_sites
                     write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                        res%forces_xrd(1, i), res%forces_xrd(2, i), res%forces_xrd(3, i)
                  end do
                  close (90)

               end if

            end if

            if (params%print_vdw_forces) then
               open (unit=90, file="forces_vdw", status="unknown")
               do i = 1, state%n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%forces_vdw(1, i), res%forces_vdw(2, i), res%forces_vdw(3, i)
               end do
               close (90)

            end if

            if (rank == 0 .and. params%print_estat_forces) then
               open (unit=90, file="forces_estat", status="unknown")
               do i = 1, state%n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%forces_estat(1, i), res%forces_estat(2, i), res%forces_estat(3, i)
               end do
               close (90)

               open (unit=90, file="charge_gradients_estat", status="unknown")
               do i = 1, dom%n_atom_pairs_by_rank(rank + 1)
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%local_properties_cart_der(1, i, model%charge_lp_index), &
                     res%local_properties_cart_der(2, i, model%charge_lp_index), &
                     res%local_properties_cart_der(3, i, model%charge_lp_index)
               end do
               close (90)

            end if

         end if
         ! For debugging the virial implementation
         if (rank == 0 .and. .false.) then
            write (*, *) "pressure_soap: ", res%virial_soap/3.d0/state%v_uc
            write (*, *) "pressure_vdw: ", res%virial_vdw/3.d0/state%v_uc
            write (*, *) "pressure_lp: ", res%virial_lp/3.d0/state%v_uc
            write (*, *) "pressure_2b: ", res%virial_2b/3.d0/state%v_uc
            write (*, *) "pressure_3b: ", res%virial_3b/3.d0/state%v_uc
            write (*, *) "pressure_core_pot: ", res%virial_core_pot/3.d0/state%v_uc
         end if
! For debugging the virial implementation
         if (rank == 0 .and. .false.) then
            write (*, *) "pressure_soap: ", res%virial_soap/3.d0/state%v_uc
            write (*, *) "pressure_vdw: ", res%virial_vdw/3.d0/state%v_uc
            do i = 1, 3
               write (*, *) res%virial_vdw(i, :)/state%v_uc
            end do
            write (*, *) "Trace of vdw pressure:", (res%virial_vdw(1, 1) + res%virial_vdw(2, 2) + res%virial_vdw(3, &
                                                                                                                 3))/3.d0/state%v_uc
            write (*, *) "pressure_2b: ", res%virial_2b/3.d0/state%v_uc
            write (*, *) "pressure_3b: ", res%virial_3b/3.d0/state%v_uc
            write (*, *) "pressure_core_pot: ", res%virial_core_pot/3.d0/state%v_uc
            write (*, *) "full vdw forces"
            do i = 1, state%n_sites
               write (*, *) i, res%forces_vdw(1:3, i)
            end do
            write (*, *) "Local virial", res%local_virial_vdw_diag
         end if

         if (params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc) then
            if (rank == 0) then
               !       Write energy and forces if we're just doing static predictions
               !       The masses should be divided by 103.6426965268d0 to have amu units, but
               !       since masses is not allocated for single point calculations, it would
               !       likely lead to a segfault
               call wrap_pbc(state%positions(1:3, 1:state%n_sites), state%a_box&
                    &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
                    & state%c_box/dfloat(state%indices(3)))
               call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                    & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                    &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                    & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                    & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                    & params%do_dipole, res%dipole, res%energies_dipole)

               call write_extxyz(state%n_sites, -loop%n_xyz, md_time, time_step,&
                    & instant_temp, instant_pressure, state%a_box&
                    &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
                    & state%c_box/dfloat(state%indices(3)), res%virial, state%xyz_species,&
                    & state%positions(1:3, 1:state%n_sites), state%velocities, res%forces,&
                    & res%energies(1:state%n_sites), state%masses, params&
                    &%write_property, params%write_array_property,&
                    & params%write_local_properties, model%local_property_labels, res%local_properties, &
                    & state%fix_atom, "trajectory_out.xyz", string, .false.,&
                    & params%do_dipole, res%local_dipoles(1:3, 1:state%n_sites))

            end if
         end if
      else
         if (rank == 0) then
            !     Do nothing
            write (*, *) '                                       |'
            write (*, *) 'You didn''t ask me to do anything!      |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if

      !   Do MD stuff here. Moved to src/turbogap_md.f90; the rank guard and the
      !   position broadcast moved with it.
!     In i-PI mode the integrator is i-PI's, so the forces just computed go
!     out over the socket and the next coordinates come back in. Everything
!     compute_md does AROUND the integration -- the skin accounting, the
!     supercell refresh, the broadcast -- happens inside the exchange.
      if (mode == "ipi") then
         call ipi_driver_exchange(rank, state%n_sites, state%positions, state%positions_prev, state%positions_diff, &
                                  state%velocities, state%a_box, state%b_box, state%c_box, state%indices, params%neighbors_buffer, &
                                  res%forces, res%energy, res%virial, loop%exit_loop, nl%rebuild_neighbors_list)
      else
         call compute_md(params, rank, ierr, state%n_sites, model%n_species, loop%md_istep, md_time, time_step, &
                         state%positions, state%positions_prev, state%positions_diff, state%velocities, res%forces, &
                         state%forces_prev, state%masses, &
                         masses_types, nl%xyz, state%xyz_species, state%a_box, state%b_box, state%c_box, state%indices, &
                         state%v_uc, res%virial, res%energy, &
                         res%energy_prev, res%energies, res%energies_soap, res%energies_2b, res%energies_3b, &
                         res%energies_core_pot, &
                         res%energies_vdw, res%energies_lp, res%energies_exp, res%energies_pdf, res%energies_sf, res%energies_xrd, &
                         res%energies_nd, res%local_properties, model%local_property_labels, instant_temp, &
                         instant_pressure, instant_pressure_prev, e_kin, e_kinetic, kb, evpera3tobar, &
                         state%fix_atom, loop%exit_loop, nl%rebuild_neighbors_list, i_image, i_nested, n_pos, nrows, &
                         filename, string, allelstopdata, ephbeta, ephfdm, ephlsc, time, &
                         cum_eel, gd_istep, &
                         target_temp, time_step_prev, res%dipole, res%local_dipoles, res%energies_dipole)
      end if

      !   Nested sampling
      !   PUT THIS INTO A MODULE!!!!!!!!!!!!!!

      !   This runs at the beginning to read in the initial images
      if (params%do_nested_sampling .and. loop%n_xyz > i_image .and. .not. params%do_md) then
         i_image = i_image + 1
         if (.not. allocated(images)) then
            allocate (images(1:i_image))
         else
            allocate (images_temp(1:i_image))
            images_temp(1:i_image - 1) = images(1:i_image - 1)
            deallocate (images)
            allocate (images(1:i_image))
            images = images_temp
            deallocate (images_temp)
         end if
         !     Save initial pool of structures
         state%velocities = 0.d0
         call from_properties_to_image(images(i_image), state%positions, state%velocities, state%masses, &
                                       res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                       res%energy_exp, E_kinetic, &
                                       state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                       state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                       res%local_dipoles, res%energies_dipole, res%dipole)
      end if

      !   This handles the nested sampling iterations after all images have
      !   been read and their energies computed
      if (params%do_nested_sampling .and. .not. loop%repeat_xyz) then
         if (i_nested == 0) then
            loop%md_istep = -1
            params%write_xyz = params%md_nsteps
            params%do_md = .true.
            if (rank == 0) then
               write (*, *) '                                       |'
               write (*, *) 'Running nested sampling algorithm with |'
               write (*, '(1X,I6,A)') loop%n_xyz, ' walkers.                        |'
               write (*, *) '                                       |'
               write (*, *) 'Target pressure in nested sampling:    |'
               write (*, '(A,ES15.7,A)') ' P = ', params%p_nested, ' bar.               |'
               write (*, *) '                                       |'
               write (*, *) '[P = 0 means total energy, rather than |'
               write (*, *) 'total enthalphy, simulation]           |'
            end if
         end if
         !     At the end of the MD/MC moves we add the image to the pool if its energy has decreased
         if (loop%md_istep == params%md_nsteps) then
            loop%md_istep = -1
            state%velocities = 0.d0
            !       Unit cell volume
            state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                                     state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
            !       We check enthalpy, not internal energy (they are the same for P = 0)
            if (res%energy + E_kinetic + params%p_nested/eVperA3tobar*state%v_uc < e_max) then
               call from_properties_to_image(images(i_image), state%positions, state%velocities, state%masses, &
                                             res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                             res%energy_exp, E_kinetic, &
                                             state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                             state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                             res%local_dipoles, res%energies_dipole, res%dipole, mc_mol_id, mc_mol_mu)
            end if
         end if
         !     This selects the highest energy image from the pool
         if (loop%md_istep == -1 .and. i_nested < params%n_nested) then
            i_nested = i_nested + 1
            nl%rebuild_neighbors_list = .true.
            i_max = 0
            e_max = -1.d100
            do i = 1, loop%n_xyz
               state%v_uc = dot_product(cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box)/ &
                            (dfloat(images(i)%indices(1)*images(i)%indices(2)*images(i)%indices(3)))
               !         We check enthalpy, not potential energy (they are the same for P = 0)
               if (images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*state%v_uc > e_max) then
                  e_max = images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*state%v_uc
                  i_max = i
               end if
            end do
            i_image = i_max
            deallocate (state%positions, state%velocities, state%masses, res%forces, state%species, &
                        state%species_supercell, state%fix_atom, state%xyz_species, state%xyz_species_supercell)
            !       Make a copy of a randonmly chosen image which is not i_image
            if (loop%n_xyz == 1) then
               i = i_image
            else
               i = i_image
               do while (i == i_image)
                  i = mod(irand(), loop%n_xyz) + 1
               end do
            end if
            if (rank == 0) then
               loop%counter = 1
               write (*, *) '                                       |'
               write (*, '(A,I8,A,I8,A)') "Nested sampling iter.:", i_nested, "/", params%n_nested, " |"
               write (*, '(A,I8,A)') " - Highest enthalpy walker:    ", i_image, " |"
               write (*, '(A,I8,A)') " - Walker selected for cloning:", i, " |"
               write (*, '(A,F15.7,A)') " - Max. enthalpy: ", e_max, " eV |"
            end if
            call from_image_to_properties(images(i), state%positions, state%velocities, state%masses, &
                                          res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                          res%energy_exp, E_kinetic, &
                                          state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                          state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                          res%local_dipoles, res%energies_dipole, res%dipole)
            state%v_uc = dot_product(cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box)/ &
                         (dfloat(images(i)%indices(1)*images(i)%indices(2)*images(i)%indices(3)))
            !       This only gets triggered if we are doing box rescaling, i.e., if the target nested sampling pressure (*not* the
            !       actual pressure for the atomic configuration) is > 0
!!!!!!!!!!!!!!!!!!!!!!!!!! Temporary hack
            if (params%scale_box_nested) then
               params%scale_box = .true.
               call random_number(rand_scale)
!!!!!!!!!!!!!!! The size of the scaling should also decrease as we reach convergence (otherwise all trial moves will be rejected)
!!!!!!!!!!!!!!! Finally, there should be a limit for the acceptable aspect ratio of the simulation box
               rand_scale = 2.d0*(rand_scale - 0.5d0)*params%nested_max_strain
               params%box_scaling_factor = reshape([1.d0 + rand_scale(1), rand_scale(6)/2.d0, rand_scale(5)/2.d0, &
                                                    rand_scale(6)/2.d0, 1.d0 + rand_scale(2), rand_scale(4)/2.d0, &
                                                    rand_scale(5)/2.d0, rand_scale(4)/2.d0, 1.d0 + rand_scale(3)], [3, 3])
               ! Make the transformation volume-preserving
               call volume_preserving_strain_transformation(state%a_box, state%b_box, state%c_box, params%box_scaling_factor)
               ! Volume scaling
               call get_ns_unbiased_volume_proposal(1.d0 - params%nested_max_volume_change, &
                                                    1.d0 + params%nested_max_volume_change, state%n_sites, rand)
               params%box_scaling_factor = params%box_scaling_factor*(rand)**(1.d0/3.d0)
               ! Each MPI process has a different set of random numbers so we need to broadcast
               call comm_bcast(comm, params%box_scaling_factor, 9)
            end if
            !       This is the so-called total enthalpy Hamiltonian Montecarlo approach (with physical masses)
            !       We do not need to broadcast the velocities here since they get broadcasted later on; otherwise
            !       we would have to do it since each MPI rank may see a different random number
            call random_number(state%velocities)
            call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), state%masses(1:state%n_sites))
            e_kin = 0.d0
            do i = 1, state%n_sites
               e_kin = e_kin + 0.5d0*state%masses(i)*dot_product(state%velocities(1:3, i), state%velocities(1:3, i))
            end do
            call random_number(rand)
            state%velocities = state%velocities/sqrt(e_kin)*sqrt(rand*(e_max - res%energy - &
                                                                       params%p_nested/eVperA3tobar*state%v_uc))
         else if (i_nested == params%n_nested) then
            loop%exit_loop = .true.
         end if
      end if

      if (rank == 0) then

         if (params%do_mc) then
            if (loop%mc_istep == params%mc_nsteps) then
               loop%exit_loop = .true.
            else
               loop%exit_loop = .false.
            end if

            if (.not. loop%exit_loop .and. ( &
                (loop%md_istep == -1) .or. &
                (params%do_md .and. ( &
                 (loop%md_istep == params%md_nsteps) .or. &
                 ((abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites)) .and. (maxval(abs(res%forces)) < &
                                                                                                  params%f_tol)) &
                 )))) then
               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               call time_start(time%mc)

               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               if (loop%mc_istep > 0) then
                  !       Evaluate the conditions for acceptance
                  !       > We have the mc conditions in mc.f90
                  !       > We care about comparing e_store to the energy of the new configuration based on the mc_movw

                  ! Reset the parameters for md / relaxation
                  trial_came_from_md = params%do_md
                  if (params%do_md) then
                     loop%md_istep = -1
                     params%do_md = .false.
                     do_mc_relax = .false.
                     ! Assume that the number of steps has already been set.
                  end if

                  if (.not. params%mc_hamiltonian) E_kinetic = 0.d0

!                 The trial configuration has to be the one `energy` belongs to.
!                 An "md" move, or a relaxation after any other move, leaves
!                 `positions` one integrator step *past* the last force
!                 evaluation: compute_md advances them after the energy was
!                 computed, and stashes the configuration it was computed at in
!                 positions_prev. md.f90 says as much -- "velocities and
!                 positions_prev are synchronous, positions is dt ahead of
!                 velocities" -- and compute_md writes positions_prev to the
!                 trajectory for exactly this reason.
!
!                 Storing `positions` here accepted or rejected x_(n+1) on the
!                 strength of E(x_n), and wrote a frame to mc_all.xyz whose
!                 energy was not the energy of its own coordinates: ~1 eV out on
!                 512 atoms after a 0.5 fs velocity-Verlet burst. Rewind by the
!                 one step, which also pairs the stored positions with the
!                 stored velocities.
                  if (trial_came_from_md) then
                     state%positions(1:3, 1:state%n_sites) = state%positions_prev(1:3, 1:state%n_sites)
                  end if

                  call from_properties_to_image(images(i_trial_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, mc_mol_id, mc_mol_mu)

                  if (params%verb > 50) write (*, *) '.......................................|'
                  if (params%verb > 50) write (*, '(A,1X,I0)') ' MC Iteration:', loop%mc_istep
                  if (params%verb > 50) write (*, '(A,1X,A)') '    Move type:', mc_move

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_prev:', images(i_current_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_prev:', images(i_current_image)&
                       &%energy + images(i_current_image)%e_kin

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_new:', images(i_trial_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_new :', images(i_trial_image)%energy &
                       &+ images(i_trial_image)%e_kin

                  state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                                           state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))

                  if (params%accessible_volume) then
                     call get_accessible_volume(state%v_uc, v_a_uc, state%species, params%radii)
                     if (params%verb > 50) write (*, '(A,F12.6,A,F12.6&
                          &,1X,A)') ' V_acc new: ', v_a_uc, ' A^3&
                          & V_acc old ', v_a_uc_prev, 'A^3 |'
                  else
                     v_a_uc = state%v_uc
                  end if

                  call get_mc_acceptance(mc_move, p_accept, &
                       res%energy + E_kinetic, &
                       images(i_current_image)%energy + images(i_current_image)%e_kin, &
                       params%t_beg, mc_mu_id, &
                       params%mc_mu, n_mc_species, state%v_uc, v_uc_prev,&
                       & v_a_uc, v_a_uc_prev, params%mc_exchange_mass, &
                       & params%mc_exchange_e0, params%mc_mu_reference, &
                       & params%p_beg, state%n_sites)

!                 call get_mc_acceptance(mc_move, p_accept, &
!                      energy + E_kinetic, &
!                      images(i_current_image)%energy + images(i_current_image)%e_kin, &
!                      params%t_beg, &
!                      params%mc_mu(mc_mu_id), n_mc_species(mc_mu_id), v_uc, v_uc_prev,&
!                      & v_a_uc, v_a_uc_prev, params&
!                      &%masses_types(mc_id(mc_mu_id)), params%p_beg)

                  call random_number(ranf)

                  if (mc_move == "insertion") n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) + 1
                  if (mc_move == "removal") n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) - 1

                  !    ACCEPT OR REJECT
                  if (params%verb > 50) write (*, '(A,1X,A,1X,A,L4,1X&
                       &,A,ES12.6,1X,A,1X,ES12.6)') 'Is ',&
                       & trim(mc_move), 'accepted?', p_accept >&
                       & ranf, ' p_accept =', p_accept, ' ranf = ',&
                       & ranf

                  if (loop%mc_istep == 1) then
                     open (unit=200, file="mc.log", status="unknown")
                     if (res%energy_exp > 0.d0) then
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current             E_exp_trial&
                             &          E_exp_current  N_tot_trial &
                             & N_mc_species_trial'
                     else
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current &
                             &          N_tot_trial  N_mc_species_trial'
                     end if

                  end if
                  if (loop%mc_istep > 1) then
                     open (unit=200, file="mc.log", status="old", position="append")
                  end if

                  ! collect the strings for the species etc
                  temp_string = ""
                  temp_string2 = ""

                  do i = 1, params%n_mc_mu
                     temp_string = ""
                     write (temp_string, "(A,1X,I8)") trim(params%mc_species(i)), n_mc_species(i)
                     temp_string2 = trim(temp_string2)//" "//trim(temp_string)
                  end do

                  if (res%energy_exp > 0.d0) then

                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                          loop%mc_istep, trim(adjustl(mc_move)), p_accept > ranf, res%energy + E_kinetic, &
                          images(i_current_image)%energy +&
                          & images(i_current_image)%e_kin, res%energy_exp,&
                          & images(i_current_image)%energy_exp,&
                          & images(i_trial_image)%n_sites,&
                          & trim(temp_string2)
                  else
                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                        loop%mc_istep, trim(adjustl(mc_move)), p_accept > ranf, res%energy + E_kinetic, &
                        images(i_current_image)%energy + images(i_current_image)%e_kin, &
                        images(i_trial_image)%n_sites, trim(temp_string2)

                  end if

                  if (loop%mc_istep >= 1) close (200)

                  if (p_accept > ranf) then
                     !             Accept
                     ! Set variables
                     loop%n_sites_prev = state%n_sites
                     v_uc_prev = state%v_uc
                     v_a_uc_prev = v_a_uc
                     virial_prev = res%virial
                     !   Assigning the default image with the accepted one
                     images(i_current_image) = images(i_trial_image)

                     if (params%n_mc_mu > 0) then
                        n_mc_species_prev = n_mc_species
                     end if

                  end if
                  if (state%n_sites > 1) then
                     instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/kB*E_kinetic
                     instant_pressure = (kB*dfloat(state%n_sites - 1)*instant_temp&
                          &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                          &/state%v_uc*eVperA3tobar
                  else
                     instant_temp = 0.0d0
                     instant_pressure = 0.0d0
                  end if

                  if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') '&
                          & Writing mc_current.xyz and&
                          & mc_all.xyz '
                     call wrap_pbc(images(i_current_image)&
                          &%positions(1:3,&
                          & 1:images(i_current_image)%n_sites),&
                          & images(i_current_image)%a_box&
                          &/dfloat(state%indices(1)),&
                          & images(i_current_image)%b_box&
                          &/dfloat(state%indices(2)),&
                          & images(i_current_image)%c_box&
                          &/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf,&
                          & params%valid_xrd, params%valid_nd,&
                          & params%do_pair_distribution, params&
                          &%do_structure_factor, params%do_xrd,&
                          & params%do_nd, string, params%do_dipole,&
                          & images(i_current_image)%dipole,&
                          & images(i_current_image)%energies_dipole)

                     call write_extxyz(images(i_current_image)%n_sites, 0, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                          images(i_current_image)%a_box/dfloat(state%indices(1)), &
                          images(i_current_image)%b_box/dfloat(state%indices(2)), &
                          images(i_current_image)%c_box/dfloat(state%indices(3)), &
                          virial_prev, images(i_current_image)%xyz_species, &
                          images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                          images(i_current_image)%velocities, &
                          images(i_current_image)%forces, &
                          images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                          images(i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & images(i_current_image)%local_properties&
                          &, images(i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & images(i_current_image)%local_dipoles)

                     call write_extxyz(images(i_current_image)%n_sites, 1, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                          images(i_current_image)%a_box/dfloat(state%indices(1)), &
                          images(i_current_image)%b_box/dfloat(state%indices(2)), &
                          images(i_current_image)%c_box/dfloat(state%indices(3)), &
                          virial_prev, images(i_current_image)%xyz_species, &
                          images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                          images(i_current_image)%velocities, &
                          images(i_current_image)%forces, &
                          images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                          images(i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & images(i_current_image)%local_properties,&
                          & images(i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .false., &
                          & params%do_dipole,&
                          & images(i_current_image)%local_dipoles)

                  end if

                  !          Add acceptance to the log file else dont
                  call time_end(time%mc)

               else ! if (mc_istep == 0)
                  temp_md_nsteps = params%md_nsteps
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, *) 'Starting MC, using parameters:         |'
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nsteps     = ', params%mc_nsteps, '     &
                       &        |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_types    = ', params%n_mc_types, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_types:                              |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '     ', params%mc_types(i), '|'
                  end do
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_accept_ratio:                       |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,F12.8,1X&
                          &,A)') '   ', params%mc_acceptance(i), '    &
                          &                  |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_swaps    = ', params%n_mc_swaps, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_swaps:  &
                       &                            |'
                  do i = 1, 2*params%n_mc_swaps
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '   ', params%mc_swaps(i), '              &
                          &        |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_move_max   = ', params%mc_move_max, &
                       & 'A   |'

                  do i = 1, params%n_mc_mu
                     write (*, '(1X,A,1X,F17.8,1X,A)') 'mc_mu         = ', params%mc_mu(1), 'eV  |'
                     write (*, '(1X,A,1X,A,1X,A)') 'mc_species    = ', trim(params%mc_species(i)), '                    |'
                  end do

                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_min_dist   = ', params%mc_min_dist, &
                       & 'A   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_lnvol_max  = ', params%mc_lnvol_max,&
                       & '    |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_write_xyz  = ', params%mc_write_xyz, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_relax      = ', params%mc_relax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nrelax     = ', params%mc_nrelax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_relax_opt  = ', params%mc_relax_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_hybrid_opt = ', params%mc_hybrid_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_optimize_exp = ', params%mc_optimize_exp&
                       &, '  |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_hamiltonian = ', params%mc_hamiltonian, '&
                       &  |'

                  if (params%verb > 50) write (*, *) '                                       |'
                  ! t_beg must

                  if (.not. allocated(images) .and. .not. params%do_nested_sampling) then
                     allocate (images(1:2))
                  else if (.not. allocated(images) .and. params%do_nested_sampling) then
                     allocate (images(1:2*i_image))
                  end if

                  if (.not. allocated(mc_mol_id)) then
                     allocate (mc_mol_id(1:state%n_sites), mc_mol_mu(1:state%n_sites))
                     mc_mol_id = 0
                     mc_mol_mu = 0
                  end if

                  if (.not. allocated(mc_id) .and. params%n_mc_mu > 0) then
                     allocate (mc_id(1:params%n_mc_mu))
                     allocate (n_mc_species(1:params%n_mc_mu))
                     allocate (n_mc_species_prev(1:params%n_mc_mu))

                     mc_id = 1
                     n_mc_species = 0

                     !    get the mc species types

                     do j = 1, params%n_mc_mu
                        do i = 1, model%n_species
                           if (params%species_types(i) == params%mc_species(j)) then
                              mc_id(j) = i
                           end if
                        end do
                     end do
                  end if

                  !       Now use the image construct to store this as the image to compare to
                  call from_properties_to_image(images(i_current_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, mc_mol_id, mc_mol_mu)

                  instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/kB*E_kinetic
                  instant_pressure = (kB*dfloat(state%n_sites - 1)*instant_temp&
                       &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                       &/state%v_uc*eVperA3tobar

                  if ((loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') ' Writing mc_current.xyz and mc_all.xyz '
                     call wrap_pbc(images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                                   images(i_current_image)%a_box/dfloat(state%indices(1)), &
                                   images(i_current_image)%b_box/dfloat(state%indices(2)), &
                                   images(i_current_image)%c_box/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                          & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                          & params%do_dipole, images(i_current_image)%dipole,&
                          & images(i_current_image)%energies_dipole)

                     call write_extxyz(images(i_current_image)%n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                          images(i_current_image)%a_box/dfloat(state%indices(1)), &
                          images(i_current_image)%b_box/dfloat(state%indices(2)), &
                          images(i_current_image)%c_box/dfloat(state%indices(3)), &
                          virial_prev, images(i_current_image)%xyz_species, &
                          images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                          images(i_current_image)%velocities, &
                          images(i_current_image)%forces, &
                          images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                          images(i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, images(i_current_image)%local_properties&
                          &, images(i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & images(i_current_image)%local_dipoles)

                     call write_extxyz(images(i_current_image)%n_sites, 1, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                          images(i_current_image)%a_box/dfloat(state%indices(1)), &
                          images(i_current_image)%b_box/dfloat(state%indices(2)), &
                          images(i_current_image)%c_box/dfloat(state%indices(3)), &
                          virial_prev, images(i_current_image)%xyz_species, &
                          images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                          images(i_current_image)%velocities, &
                          images(i_current_image)%forces, &
                          images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                          images(i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, images(i_current_image)%local_properties&
                          &, images(i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .true., &
                          & params%do_dipole,&
                          & images(i_current_image)%local_dipoles)

                     v_uc_prev = dot_product(cross_product(state%a_box, state%b_box), &
                                             state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
                     if (params%accessible_volume) then
                        call get_accessible_volume(v_uc_prev, v_a_uc_prev, state%species, params%radii)
                     else
                        v_a_uc_prev = v_uc_prev
                     end if
                  end if

               end if

               !  Now start the mc logic: first, use the stored images properties
               call from_image_to_properties(images(i_current_image), state%positions, state%velocities, state%masses, &
                                             res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                             res%energy_exp, E_kinetic, &
                                             state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                             state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                             res%local_dipoles, res%energies_dipole, res%dipole, mc_mol_id, mc_mol_mu)

               call perform_mc_step(&
                    & state%positions, state%species, state%xyz_species, state%masses, state%fix_atom,&
                    & state%velocities, state%positions_prev, state%positions_diff, disp, d_disp, params%n_local_properties,&
                    & params%mc_acceptance, params%mc_mu_acceptance, res%local_properties, &
                    images(i_current_image)%local_properties, res%energies,&
                    & res%forces, state%forces_prev, state%n_sites, params%n_mc_mu, mc_mu_id, n_mc_species,&
                    & mc_move, params%mc_species,&
                    & params%mc_move_max, params%mc_min_dist, params%mc_max_dist, params%mc_max_insertion_trials, &
                    params%mc_lnvol_max, params%mc_types, params%masses_types, species_idx,&
                    & images(i_current_image)%positions,&
                    & images(i_current_image)%species,&
                    & images(i_current_image)%xyz_species,&
                    & images(i_current_image)%fix_atom,&
                    & images(i_current_image)%masses, state%a_box(1:3), state%b_box(1:3),&
                    & state%c_box(1:3), state%indices, params%do_md, params%mc_relax,&
                    & loop%md_istep, mc_id, E_kinetic, instant_temp, params%t_beg,&
                    & params%n_mc_swaps, params%mc_swaps, params%mc_swaps_id, &
                    & params%species_types, params%mc_hamiltonian,&
                    & params%n_mc_relax_after, params&
                    &%mc_relax_after, do_mc_relax, params%verb, &
                    params%mc_n_planes, params%mc_planes, params%mc_max_dist_to_planes, &
                    params%mc_planes_restrict_to_polyhedron, &
                    params%mc_molecules, mc_mol_id, mc_mol_mu, &
                    images(i_current_image)%mc_mol_id, images(i_current_image)%mc_mol_mu, mc_mol_next)

               nl%rebuild_neighbors_list = .true.

               ! NOTE: the species_supercell and xyz_species_supercell are
               ! not commensurate with the new image as these have not been
               ! calculated. If reading from an outputted xyz file, then it
               ! should be okay but really the new atoms should be added to
               ! the supercell in the usual way, but for convenience, one has
               ! not done that.

               if (params%mc_relax .and. do_mc_relax) then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_relax_opt
                  params%md_nsteps = params%mc_nrelax

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, E_kinetic, state%masses, instant_temp, params%t_beg, &
                                            params%velocity_distribution)

                  if (params%mc_hamiltonian) E_kinetic_prev = E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if
               ! If doing md, don't relax
               if (mc_move == 'md') then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_hybrid_opt
                  params%md_nsteps = temp_md_nsteps

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, E_kinetic, state%masses, instant_temp, params%t_beg, &
                                            params%velocity_distribution)
                  if (params%mc_hamiltonian) E_kinetic_prev = E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if

               if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                    modulo(loop%mc_istep, params%write_xyz) == 0)) then

                  call wrap_pbc(state%positions(1:3, 1:state%n_sites), &
                                state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                                state%c_box/dfloat(state%indices(3)))
                  call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                       & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                       &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                       & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                       & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                  call write_extxyz(state%n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                       state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                          state%c_box/dfloat(state%indices(3)), &
                       res%virial, state%xyz_species, &
                       state%positions(1:3, 1:state%n_sites), state%velocities, &
                       res%forces, res%energies(1:state%n_sites), state%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & model%local_property_labels, res%local_properties&
                       &, state%fix_atom, mc_file, string, .true.)
               end if
               ! As we have moved/added/removed, we must check the supercell and  broadcast the results

               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .true., state%fix_atom, params%t_beg, &
                             params%write_array_property(6), .true., params%randomize_velocities)

            else
               if (mc_move == 'md') then
                  if (params%print_progress .and. loop%md_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%md_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%md_istep == params%md_nsteps - 1 .or. &
                           (abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites) .and. &
                            maxval(abs(res%forces)) < params%f_tol) .and. loop%md_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%md_istep < params%md_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%md_istep + 1)/params%md_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%md_istep + 1)/params%md_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%mc_hamiltonian) then
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step: H =&
                          & T + V = ", res%energy + E_kinetic, ",&
                          & iteration ", loop%md_istep, "/", params&
                          &%md_nsteps
                  else
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step:&
                          & energy = ", res%energy, ", iteration ",&
                          & loop%md_istep, "/", params%md_nsteps
                  end if

                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') '&
                       & core_pot energy:', sum(res%energies_core_pot),&
                       & 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') '&
                       & vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               else
                  if (params%print_progress .and. loop%mc_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%mc_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%mc_istep == params%mc_nsteps - 1 .and. loop%mc_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%mc_istep < params%mc_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%mc_istep + 1)/params%mc_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%mc_istep + 1)/params%mc_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%verb > 50 .and. do_mc_relax) write (*, '(1X,A,1X,F20.8,1X,A&
                       &,1X,I8,1X,A,1X,I8)') "MC Relax md step: energy &
                       &= ", res%energy, ", iteration ", loop%md_istep, "/",&
                       & params%mc_nrelax
                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               end if
            end if
         end if

      end if

      ! NOTE!! One tried for far far too long to be smart and implement some
      ! sort of conditional broadcasting: having a logical array named
      ! broadcast, which perform_mc_step would then to set values to
      ! true. Specific indexes referenced specific quantities to be
      ! broadcasted, which allowed for the broadcasting amount to be
      ! dependent on the step, e.g. if it were an insertion step then
      ! positions, masses, n_sites, etc would have to be broadcast, whereas
      ! for a simple move only positions had to be broadcasted. This array
      ! would then subsequently be broadcast to all other ranks, thereby
      ! allowing for the minimum number of allocations and
      ! communication. BUT, for some reason, this led to segfaults
      ! (corrupted unsorted chunks or something of that sort).

      ! This doesn't make sense to be as all ranks have the same broadcast
      ! array (as it is broadcasted before) so it seems like it should work
      ! but it does not! Hence, in the following broadcasting, everything is
      ! transmitted.

      ! This can be optimised, so please do if you are smarter than me

      call time_start(time%mpi)
      call comm_bcast(comm, params%do_md)
      call comm_bcast(comm, loop%md_istep)
      call time_end(time%mpi)
      call domain_sync_state(dom, comm, state, params, time)
      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. loop%md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
      end if
      call get_time(time1)
      !   Parallel neighbors list build
      call comm_bcast(comm, nl%rebuild_neighbors_list)

      if (nl%rebuild_neighbors_list) then
         deallocate (nl%rjs, nl%xyz, nl%thetas, nl%phis, nl%neighbor_species)
         deallocate (nl%neighbors_list, nl%n_neigh)
         deallocate (nl%n_neigh_local)
      end if
      if ((params%do_nested_sampling .and. .not. params%do_mc) .and. &
          (params%do_md .and. (loop%md_istep == params%md_nsteps .or. loop%exit_loop))) then
         deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, state%species_supercell, &
                     dom%do_list)
         if (allocated(state%velocities)) deallocate (state%velocities)
      end if
      if (params%do_mc .and. params%do_md) then
         if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop)) then
            deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, &
                        state%species_supercell, dom%do_list)
            if (allocated(state%velocities)) deallocate (state%velocities)
         end if
      end if

      if ((params%do_md .and. .not. params%do_mc) .and. &
          (loop%md_istep == params%md_nsteps .or. loop%exit_loop) .and. rank == 0) then
         deallocate (state%positions_prev, state%forces_prev)
      end if
      if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop) .and. rank == 0) then
         if (allocated(state%forces_prev)) deallocate (state%forces_prev)
         if (allocated(state%positions_prev)) deallocate (state%positions_prev)
      end if

      if (params%exp_forces .and. (loop%md_istep == params%md_nsteps .or.&
           & loop%mc_istep == params%mc_nsteps .or. loop%exit_loop)) then
         do i = 1, params%n_exp
            if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
            if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
            if (allocated(params%exp_data(i)%y_pred)) deallocate (params%exp_data(i)%y_pred)
         end do
      end if

!     Close the per-step stopwatch and charge it to whichever side of the
!     boundary this step fell on. mad_ir_applied was set for THIS step in the
!     force block above, so the two accumulators separate exactly at the step
!     the ensemble filled.
      if (params%valid_ir .and. params%do_md .and. loop%md_istep >= 0) then
         call get_time(mad_ir_t_now)
         if (mad_ir_applied) then
            mad_ir_t_post = mad_ir_t_post + (mad_ir_t_now - mad_ir_step_beg)
            mad_ir_n_post = mad_ir_n_post + 1
         else
            mad_ir_t_pre = mad_ir_t_pre + (mad_ir_t_now - mad_ir_step_beg)
            mad_ir_n_pre = mad_ir_n_pre + 1
         end if
      end if

      if (.not. params%do_mc) loop%n_sites_prev = state%n_sites
      dom%n_atom_pairs_by_rank_prev = dom%n_atom_pairs_by_rank(rank + 1)

      call comm_bcast(comm, loop%exit_loop)
      if (loop%exit_loop) exit
      ! End of loop through structures in the xyz file or MD steps
   end do

!  i-PI has said EXIT, or something else ended the loop. Close the socket
!  before the reports below, so that i-PI sees the driver leave cleanly
!  rather than timing out on a half-open connection.
   if (mode == "ipi") call ipi_driver_close(rank)

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
   if (ir_from_traj) then
      call time_start(time%ir_predict)
      if (params%valid_ir) then
!        There is an experiment: restrict it to [ir_nu_min, ir_nu_max] and
!        weight it exactly as the MAD bias would, so the mismatch printed here
!        is the same number a biased run would be minimising.
         call mad_ir_select_range(params%exp_data(params%ir_idx)%data(1, :), &
                                  params%exp_data(params%ir_idx)%data(2, :), &
                                  params%ir_nu_min, params%ir_nu_max, &
                                  params%ir_weight_by_spacing, &
                                  ir_nu_exp, ir_I_exp, ir_wgt_exp, mad_ir_ok, mad_ir_msg)
         if (.not. mad_ir_ok) then
            write (*, *) "ERROR: ", trim(mad_ir_msg)
            stop
         end if
      else
         allocate (ir_nu_exp(1:1), ir_I_exp(1:1), ir_wgt_exp(1:1))
         ir_nu_exp = 0.d0; ir_I_exp = 0.d0; ir_wgt_exp = 1.d0
      end if

      call ir_fft_frames_finish(ir_fft_frames, ir_fft_cfg, params%ir_frame_dt, &
                                params%ir_frame_dt_tol, ir_nu_exp, ir_I_exp, &
                                ir_wgt_exp, size(ir_nu_exp), params%valid_ir, &
                                params%ir_match_scale, params%ir_match_offset, &
                                "ir_fft_spectrum.dat", "ir_fft_dipoles.dat", &
                                rank == 0, params%ir_fft_write_dipoles, &
                                ir_fft_res, ir_fft_dt_used, ir_fft_scale_fit, &
                                ir_fft_offset_fit, ir_fft_dissim, ir_fft_dissim_ref, &
                                ir_fft_ok, ir_fft_msg)
      call time_end(time%ir_predict)

      if (.not. ir_fft_ok) then
         if (rank == 0) then
            write (*, *) ""
            write (*, *) "ERROR: ", trim(ir_fft_msg)
         end if
!        A NONZERO exit, unlike the bare `stop` used elsewhere in this file.
!        This path is driven by scripts -- post-processing a directory of
!        trajectories is the obvious use -- and the whole design here is
!        organised against failing silently. Exiting 0 with no spectrum
!        written is precisely that failure wearing a message.
         stop 1
      end if

      if (rank == 0) then
         write (*, *) '                                       |'
         write (*, *) 'IR spectrum from the trajectory:       |'
         write (*, '(A,I12,A)') '  *) frames read:       ', ir_fft_frames%n, '         |'
         write (*, '(A,F12.4,A)') '  *) frame interval:    ', ir_fft_dt_used, ' fs      |'
         write (*, '(A,I12,A)') '  *) lags kept:         ', ir_fft_res%n_lag, '         |'
         write (*, '(A,F12.4,A)') '  *) resolution:        ', ir_fft_res%resolution, ' cm^-1   |'
         write (*, '(A,F12.4,A)') '  *) bin spacing:       ', ir_fft_res%d_nu, ' cm^-1   |'
         write (*, '(A,F12.1,A)') '  *) Nyquist:           ', ir_fft_res%nyquist, ' cm^-1   |'
         write (*, '(A,I12,A)') '  *) bins written:      ', ir_fft_res%n_freq, '         |'
         if (params%valid_ir .and. ir_fft_dissim_ref > 0.d0) then
            write (*, '(A,F12.6,A)') '  *) rel. mismatch:     ', &
               dsqrt(ir_fft_dissim/ir_fft_dissim_ref), '         |'
         end if
         if (len_trim(ir_fft_msg) > 0) then
            write (*, *) '  *) ', trim(ir_fft_msg)
         end if
         write (*, *) '                                       |'
      end if

      call ir_fft_free(ir_fft_res)
      call ir_fft_frames_reset(ir_fft_frames)
      deallocate (ir_nu_exp, ir_I_exp, ir_wgt_exp)
   end if

   if (params%do_md .or. params%do_prediction .or. params%do_mc) then
      call get_time(time2)
      if (rank == 0) then
         if (params%do_md .and. .not. params%do_nested_sampling) then
            write (*, *) '                                       |'
            write (*, '(I8,A,F13.3,A)') loop%md_istep, ' MD steps:', time2 - time3, ' seconds |'
         end if
         if (params%do_mc) then
            write (*, *)
            write (*, *) '                                       |'
            write (*, '(I8,A,F13.3,A)') loop%mc_istep, ' MC steps:', time2 - time3, ' seconds |'
         end if

         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *          Setup:', time%setup(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - input+pot.:', time%read_input(3), ' seconds |'
         if (comm_with_mpi) then
            write (*, '(A,F13.3,A)') '     -  MPI setup:', time%mpi_setup(3), ' seconds |'
         end if
         write (*, '(A,F13.3,A)') ' * Read XYZ files:', time%read_xyz(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' * Neighbor lists:', time%neigh(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' *  GAP desc/pred:', time%gap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - soap_turbo:', time%soap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         2b:', time%gap_2b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         3b:', time%gap_3b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -   core_pot:', time%gap_core_pot(3), ' seconds |'
!       vdw is a parent, not one of the GAP children above: compute_vdw runs
!       outside the time%gap region and sum_times adds it in its own right.
!       Printing it indented under GAP said otherwise.
         if (params%vdw_type /= "none") then
            write (*, '(A,F13.3,A)') ' *            vdw:', time%vdw(3), ' seconds |'
         end if
         if (model%valid_xps .or. params%do_pair_distribution .or. params&
              &%do_structure_factor .or. params%do_xrd .or. params%do_nd) write (*, '(A&
              &,F13.3,A)') ' *  Exp. pred.   :', time%pdf(3) + time%sf(3) + time%xrd(3) + time%nd(3), ' seconds&
              & |'
         if (model%valid_xps) write (*, '(A,F13.3,A)') '     -        xps:',&
              & time%xps(3), ' seconds |'
         if (params%do_pair_distribution) write (*, '(A,F13.3,A)') '     -        pdf:', time%pdf(3), ' seconds |'
         if (params%do_structure_factor) write (*, '(A,F13.3,A)') '     -         sf:', time%sf(3), ' seconds |'
         if (params%do_xrd) write (*, '(A,F13.3,A)') '     -        xrd:', time%xrd(3), ' seconds |'
         if (params%do_nd) write (*, '(A,F13.3,A)') '     -         nd:', time%nd(3), ' seconds |'

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
               if (mad_ir_step_first >= 0) then
                  write (*, '(A,I13,A)') '   1st spectrum @:', mad_ir_step_first, ' step    |'
                  write (*, '(A,F13.3,A)') '   1st spectrum @:', mad_ir_t_first, ' seconds |'
               else
                  write (*, *) '   no spectrum: ensemble never filled  |'
               end if
            end if
            if (mad_ir_n_pre > 0 .and. mad_ir_n_post > 0) then
               mad_ir_rate_pre = mad_ir_t_pre/dfloat(mad_ir_n_pre)
               mad_ir_rate_post = mad_ir_t_post/dfloat(mad_ir_n_post)
               write (*, '(A,F13.5,A)') '  s/step unbiased:', mad_ir_rate_pre, ' seconds |'
               write (*, '(A,F13.5,A)') '  s/step   biased:', mad_ir_rate_post, ' seconds |'
               if (mad_ir_rate_pre > 0.d0) then
                  write (*, '(A,F13.3,A)') '  bias slowdown  :', &
                     mad_ir_rate_post/mad_ir_rate_pre, ' x       |'
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

         if (do_electrostatics) then
            write (*, '(A,F13.3,A)') ' * Electrostatics:', time%estat(3), ' seconds |'
         end if
         if (params%do_md) then
            write (*, '(A,F13.3,A)') ' *  MD algorithms:', time%md(3), ' seconds |'
         end if
         if (params%do_mc) then
            write (*, '(A,F13.3,A)') ' *  MC algorithms:', time%mc(3), ' seconds |'
         end if

         if (comm_with_mpi) then
            write (*, '(A,F13.3,A)') ' *  MPI comms.   :', time%mpi(3) + time%mpi_positions(3) + time%mpi_ef(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -  pos & vel:', time%mpi_positions(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     - E & F brc.:', time%mpi_ef(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -  MPI misc.:', time%mpi(3), ' seconds |'
         end if
!       Miscellaneous is what the parent buckets do not account for.  It used
!       to be written out here as one long subtraction, which is how it came to
!       subtract time%gap and the mpi_ef reduce nested inside it and print a
!       negative number.  sum_times owns the list now (src/timing.f90), so the
!       set summed here and the set declared as parents there cannot disagree.
!
!       Accounted-for is printed beside it so the arithmetic is visible: the
!       three numbers below have to add up, and a reader can see at a glance how
!       much of the run the buckets actually name.  A large Miscellaneous is a
!       statement that something real is not being measured -- which is how the
!       setup bucket above came to exist.
         time%total(3) = time2 - time3
         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *  Accounted for:', sum_times(time), ' seconds |'
         write (*, '(A,F13.3,A)') ' *  Miscellaneous:', time%total(3) - sum_times(time), ' seconds |'
         write (*, '(A,F13.3,A)') ' *     Total time:', time%total(3), ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end if

#ifdef _GPU
   do i = 1, model%n_soap_turbo
      if (.not. model%soap_turbo_hypers(i)%recompute_basis) then
         call gpu_free_async(model%soap_turbo_hypers(i)%W_d, gpu_stream)
         call gpu_free_async(model%soap_turbo_hypers(i)%S_d, gpu_stream)
         call gpu_free_async(model%soap_turbo_hypers(i)%multiplicity_array_d, gpu_stream)
      end if
   end do
#endif
   if (allocated(state%fix_atom)) deallocate (state%fix_atom)
   if (allocated(state%positions)) deallocate (state%positions)
   if (allocated(state%velocities)) deallocate (state%velocities)
   if (allocated(state%positions_diff)) deallocate (state%positions_diff)

   if (allocated(res%energies)) deallocate (res%energies)
   if (allocated(res%local_dipoles)) deallocate (res%local_dipoles, res%this_local_dipoles)
   if (allocated(res%energies_dipole)) deallocate (res%energies_dipole, res%this_energies_dipole)
   if (allocated(res%energies_soap)) deallocate (res%energies_soap)
   if (allocated(res%energies_2b)) deallocate (res%energies_2b)
   if (allocated(res%energies_3b)) deallocate (res%energies_3b)
   if (allocated(res%energies_core_pot)) deallocate (res%energies_core_pot)
   if (allocated(res%energies_vdw)) deallocate (res%energies_vdw)
   if (allocated(res%energies_exp)) deallocate (res%energies_exp)
   if (allocated(res%energies_lp)) deallocate (res%energies_lp)
   if (allocated(res%energies_pdf)) deallocate (res%energies_pdf)
   if (allocated(res%energies_sf)) deallocate (res%energies_sf)
   if (allocated(res%energies_xrd)) deallocate (res%energies_xrd)
   if (allocated(res%energies_nd)) deallocate (res%energies_nd)

   if (allocated(res%this_energies)) deallocate (res%this_energies)
   if (allocated(res%this_energies_vdw)) deallocate (res%this_energies_vdw)
   if (allocated(res%this_energies_lp)) deallocate (res%this_energies_lp)
   if (allocated(res%this_energies_pdf)) deallocate (res%this_energies_pdf)
   if (allocated(res%this_energies_sf)) deallocate (res%this_energies_sf)
   if (allocated(res%this_energies_xrd)) deallocate (res%this_energies_xrd)
   if (allocated(res%this_energies_nd)) deallocate (res%this_energies_nd)

   if (allocated(res%forces)) deallocate (res%forces)
   if (allocated(res%forces_soap)) deallocate (res%forces_soap)
   if (allocated(res%forces_2b)) deallocate (res%forces_2b)
   if (allocated(res%forces_3b)) deallocate (res%forces_3b)
   if (allocated(res%forces_core_pot)) deallocate (res%forces_core_pot)
   if (allocated(res%forces_vdw)) deallocate (res%forces_vdw)
   if (allocated(res%forces_lp)) deallocate (res%forces_lp)
   if (allocated(res%forces_pdf)) deallocate (res%forces_pdf)
   if (allocated(res%forces_sf)) deallocate (res%forces_sf)
   if (allocated(res%forces_xrd)) deallocate (res%forces_xrd)
   if (allocated(res%forces_nd)) deallocate (res%forces_nd)

   if (allocated(res%this_forces)) deallocate (res%this_forces)
   if (allocated(res%this_forces_vdw)) deallocate (res%this_forces_vdw)
   if (allocated(res%this_forces_lp)) deallocate (res%this_forces_lp)
   if (allocated(res%this_forces_pdf)) deallocate (res%this_forces_pdf)
   if (allocated(res%this_forces_sf)) deallocate (res%this_forces_sf)
   if (allocated(res%this_forces_xrd)) deallocate (res%this_forces_xrd)
   if (allocated(res%this_forces_nd)) deallocate (res%this_forces_nd)

   if (allocated(res%local_properties)) deallocate (res%local_properties)
   if (allocated(res%local_properties_cart_der)) deallocate (res%local_properties_cart_der)
   if (allocated(res%this_local_properties)) deallocate (res%this_local_properties)
   if (allocated(res%this_local_properties_cart_der)) deallocate (res%this_local_properties_cart_der)

   if (allocated(model%soap_turbo_hypers)) deallocate (model%soap_turbo_hypers)
   if (allocated(model%distance_2b_hypers)) deallocate (model%distance_2b_hypers)
   if (allocated(model%angle_3b_hypers)) deallocate (model%angle_3b_hypers)
   if (allocated(model%core_pot_hypers)) deallocate (model%core_pot_hypers)

   deallocate (dom%n_atom_pairs_by_rank)
   if (allocated(model%n_local_properties_mpi)) deallocate (model%n_local_properties_mpi)
   if (allocated(model%local_properties_n_sparse_mpi_soap_turbo)) deallocate (model%local_properties_n_sparse_mpi_soap_turbo)
   if (allocated(model%local_properties_dim_mpi_soap_turbo)) deallocate (model%local_properties_dim_mpi_soap_turbo)
   if (allocated(model%has_local_properties_mpi)) deallocate (model%has_local_properties_mpi)

   if (allocated(model%local_property_labels)) deallocate (model%local_property_labels)
   if (allocated(model%local_property_indexes)) deallocate (model%local_property_indexes)
   if (allocated(dom%do_list)) deallocate (dom%do_list)
   if (allocated(params%write_local_properties)) deallocate (params%write_local_properties)

   if (params%vdw_type == "ts+mbd") then
      if (rank == 0) then
         open (unit=30, file="mbd_ts_scaling.dat", status="unknown")
         do i = 1, state%n_sites
#ifdef _MPIF90
            write (30, *) res%this_mbd_ts_scaling(i)
#else
            write (30, *) res%mbd_ts_scaling(i)
#endif
         end do
         close (30)
      end if
   end if

   if (rank == 0) then
      write (*, *) '                                       |'
      write (*, *) 'End of execution                       |'
      write (*, *) '_______________________________________/'
   end if

#ifdef _GPU
!  The high-water mark, which is the number that sizes the next run.
!
!  Before gpu_context_finalize, which calls hipDeviceReset and takes the whole
!  context down -- after it there is nothing left to ask. Printed unconditionally
!  and to stderr: it costs one line, and "what did that actually use" is the
!  first question asked after any run that was close to the limit.
   if (rank == 0) call gpu_memory_report("end of run")
#endif
   call comm_finalize(comm)

   call gpu_context_finalize(params, n_omp)

end program turbogap
