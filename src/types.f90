! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, types.f90, is copyright (c) 2019-2022, Miguel A. Caro and
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

module types

   use kinds

   implicit none

   ! GAP+descriptor data structure for SOAP
   type exp_data_container
      character*1024 :: file_data = "none"
      character*1024 :: label
      character*1024 :: input = "default"
      integer :: n_data
      integer :: n_samples = 200
      logical :: compute_similarity = .false.
      logical :: compute_exp = .false.
      logical :: wrote_exp = .false.
      logical :: user_range = .false.
      logical :: compute_forces = .false.
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: y_pred(:)
      real(dp), allocatable :: y_pred_prev(:)
      real(dp) :: similarity
      real(dp) :: range_min = 0.d0
      real(dp) :: range_max = 1.d0
      real(dp) :: mag
   end type exp_data_container

   type local_property_soap_turbo
      real(dp), allocatable :: Qs(:, :)
      real(dp), allocatable :: alphas(:)
      real(dp), allocatable :: cutoff(:)
      real(dp) :: zeta
      real(dp) :: delta
      real(dp) :: V0
      character*1024 :: file_alphas
      character*1024 :: file_desc
      character*1024 :: label
      integer :: n_sparse
      integer :: dim
      logical :: do_derivatives = .false.
      logical :: compute = .true.
   end type local_property_soap_turbo

   type soap_turbo
      real(dp), allocatable :: nf(:)
      real(dp), allocatable :: rcut_hard(:)
      real(dp), allocatable :: rcut_soft(:)
      real(dp), allocatable :: atom_sigma_r(:)
      real(dp), allocatable :: atom_sigma_t(:)
      real(dp), allocatable :: atom_sigma_r_scaling(:)
      real(dp), allocatable :: atom_sigma_t_scaling(:)
      real(dp), allocatable :: amplitude_scaling(:)
      real(dp), allocatable :: central_weight(:)
      real(dp), allocatable :: global_scaling(:)
      real(dp), allocatable :: alphas(:)
      real(dp), allocatable :: Qs(:, :)
      real(dp), allocatable :: cutoff(:)
      real(dp), allocatable :: vdw_Qs(:, :)
      real(dp), allocatable :: vdw_alphas(:)
      real(dp), allocatable :: vdw_cutoff(:)
      real(dp), allocatable :: compress_P_el(:)
      real(dp) :: zeta = 2.d0
      real(dp) :: delta = 1.d0
      real(dp) :: rcut_max
      real(dp) :: vdw_zeta
      real(dp) :: vdw_delta
      real(dp) :: vdw_V0
      integer, allocatable :: alpha_max(:)
      integer, allocatable :: compress_P_i(:)
      integer, allocatable :: compress_P_j(:)
      integer :: n_species
      integer :: central_species = 0
      integer :: dim
      integer :: l_max
      integer :: radial_enhancement = 0
      integer :: n_max
      integer :: n_sparse
      integer :: vdw_n_sparse
      integer :: compress_P_nonzero
      integer :: n_local_properties = 0
      integer :: vdw_index = 0
      integer :: core_electron_be_index = 0
      character*1024 :: file_alphas
      character*1024 :: file_desc
      character*1024 :: file_compress = "none"
      character*1024 :: file_vdw_alphas
      character*1024 :: file_vdw_desc
      character*64 :: basis = "poly3"
      character*64 :: compress_mode = "none"
      character*32 :: scaling_mode = "polynomial"
      character*8, allocatable :: species_types(:)
      logical :: compress_soap = .false.
      logical :: has_vdw = .false.
      logical :: has_charges = .false.
      logical :: has_core_electron_be = .false.
      logical :: has_local_properties = .false.
!     A dipole model is a GAP whose "energy" is a fictitious scalar fitted so that
!     mu_i = dE_i/dr_i (the derivative w.r.t. the central atom's *own* position)
!     reproduces the local dipole. Such a descriptor contributes to mu only: its
!     energy, forces and virial are meaningless and are not accumulated.
      logical :: is_dipole_model = .false.
      type(local_property_soap_turbo), allocatable :: local_property_models(:)
   end type soap_turbo

! GAP+descriptor data structure for distance_2b
   type distance_2b
      real(dp), allocatable :: cutoff(:)
      real(dp), allocatable :: alphas(:)
      real(dp), allocatable :: Qs(:, :)
      real(dp) :: delta = 1.d0
      real(dp) :: sigma = 1.d0
      real(dp) :: rcut
      real(dp) :: buffer = 0.5d0
      integer :: dim = 1
      integer :: n_sparse
      character*1024 :: file_alphas
      character*1024 :: file_desc
      character*8 :: species1
      character*8 :: species2
   end type distance_2b

! GAP+descriptor data structure for distance_2b
   type angle_3b
      real(dp), allocatable :: cutoff(:)
      real(dp), allocatable :: alphas(:)
      real(dp), allocatable :: Qs(:, :)
      real(dp) :: delta = 1.d0
      real(dp) :: sigma(1:3) = 1.d0
      real(dp) :: rcut
      real(dp) :: buffer = 0.5d0
      integer :: dim = 3
      integer :: n_sparse
      character*1024 :: file_alphas
      character*1024 :: file_desc
      character*8 :: species_center
      character*8 :: species1
      character*8 :: species2
      character*3 :: kernel_type = "exp"
   end type angle_3b

! Data structure for core_pot
   type core_pot
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: V(:)
      real(dp), allocatable :: dVdx2(:)
      real(dp) :: yp1
      real(dp) :: ypn
      integer :: n
      character*1024 :: core_pot_file
      character*8 :: species1
      character*8 :: species2
   end type core_pot

! These is the type for the input parameters
   type options_estat
      logical :: damped = .false.
      logical :: tsf = .true.
      logical :: sp = .false.
      logical :: gsf = .false.
      logical :: self_energy_correction = .false.
      logical :: damped_cosine = .false.
   end type options_estat

!  A rigid molecule that grand-canonical moves may insert and remove as one
!  object. Read from an xyz file by read_mc_molecule; positions are stored
!  relative to the centre of mass so that an insertion is a rotation about the
!  origin followed by a translation.
   type mc_molecule
      integer :: n_atoms = 0
      real(dp), allocatable :: positions(:, :)
      real(dp), allocatable :: masses(:)
      integer, allocatable :: species(:)
      character*8, allocatable :: xyz_species(:)
!     Total mass, for the thermal de Broglie wavelength, and the sum of the e0
!     of its atoms, for mc_mu_reference = "e0".
      real(dp) :: total_mass = 0.d0
      real(dp) :: e0_total = 0.d0
!     Largest distance of any atom from the centre of mass. An insertion needs
!     it to know how far the molecule reaches from the point it is placed at.
      real(dp) :: radius = 0.d0
      logical :: is_molecule = .false.
      character*1024 :: file = "none"
   end type mc_molecule

   type input_parameters

!     ==================================================================
!     GENERAL
!     ==================================================================
      character*1024 :: atoms_file
      character*8, allocatable :: species_types(:)
      real(dp), allocatable :: masses_types(:)
      real(dp), allocatable :: e0(:)
      real(dp), allocatable :: radii(:)
      integer :: which_atom = 0
      logical :: all_atoms = .true.
      logical :: do_timing = .false.
      logical :: print_progress = .true.
!     Verbosity for debugging. 0 is quiet.
      integer :: verb = 0
!     Seed for the intrinsic pseudo-random number generator. Zero (the default)
!     leaves the compiler's default sequence untouched; any other value makes
!     runs reproducible, which is what the CPU/GPU regression comparisons rely
!     on when initial velocities have to be randomized.
      integer :: random_seed_value = 0

!     ==================================================================
!     WHAT THE RUN COMPUTES
!     ==================================================================
      logical :: do_prediction = .false.
      logical :: do_forces = .false.
      logical :: do_derivatives = .false.
      logical :: do_derivatives_fd = .false.
!     Which of the two filter recursions the SOAP radial expansion uses.
!     Applied to soap_turbo_radial's legacy_filter_seed in turbogap_setup; the
!     comment at the top of soap_turbo_radial.f90 says what the two mean.
      logical :: soap_radial_legacy_filter = .true.
      logical :: do_md = .false.
      logical :: do_mc = .false.
      logical :: do_nested_sampling = .false.
      logical :: do_exp = .false.
      logical :: do_xps = .false.
!     Set from the .gap file, not the input file: true when at least one
!     soap_turbo descriptor is flagged dipole_model
      logical :: do_dipole = .false.

!     MAD IR. Bias a trajectory towards an experimental infrared spectrum
!     computed from the dipole autocorrelation function; see src/mad_ir.f90.
!     Needs a dipole model, do_derivatives, and soap_radial_legacy_filter off.
!     The observable is switched on the same way as every other one: name "ir"
!     in exp_labels with a data file beside it. Its weight is
!     exp_energy_scales(ir_idx), ramped like the rest, and exp_forces decides
!     whether the gradient of the mismatch reaches the forces.
      logical :: valid_ir = .false.
!     Plain IR prediction, no experiment and no bias: run the trajectory,
!     accumulate the dipole over all of it, transform once at the end. The
!     sizing runs the other way round from the bias -- the run length is given
!     and the resolution follows from it -- and the result lands in
!     ir_prediction.dat with the sampling limits in its header.
      logical :: do_ir = .false.
!     Points on the predicted grid. Zero, the default, is one point per
!     resolution element, which is as much as the transform can carry; a
!     larger number is a smoother curve over the same information.
      integer :: ir_n_samples = 0
!     true when some observable other than IR is asked for, i.e. when the
!     per-frame prediction pipeline actually has work to do
      logical :: do_exp_structural = .false.
      integer :: ir_idx = 0
      logical :: ir_match_scale = .true.
!     Write the predicted spectrum out. Named and gated like every other
!     observable's prediction (write_xrd, write_pair_distribution, ...): on
!     every write_xyz step a block is appended to ir_prediction.dat, so the
!     record covers the whole trajectory instead of only its last frame.
!     ir_spectrum.dat holds the latest frame alone, and ir_exp.dat the
!     experiment restricted to the fitted range.
      logical :: write_ir = .false.
!     The name this switch had before it wrote a trajectory. Kept because
!     inputs use it; it sets write_ir and nothing else.
      logical :: ir_write_spectrum = .false.
      character*1024 :: ir_restart_file = "ir_restart.dat"
      character*32 :: ir_window = "hann"
      integer :: ir_stride = 1
      integer :: ir_lag_factor = 2
      real(dp) :: ir_resolution = 20.d0
!     Did the input name ir_resolution? A prediction run takes whatever
!     resolution the trajectory length gives unless one was asked for, and the
!     default value cannot be told from a chosen one by its value alone.
      logical :: ir_resolution_set = .false.
      real(dp) :: ir_nu_min = 0.d0
      real(dp) :: ir_nu_max = 4000.d0
      real(dp) :: ir_nu_power = 2.d0
!     Estimator choices. The defaults are the correct ones; every switch here
!     exists so the old behaviour can be recovered for reproducing an earlier
!     result, not because the choice is a matter of taste. mad_ir.f90's header
!     block gives the derivation and the measured cost of each.
!
!     ir_subtract_mean      correlate mu - <mu>, as linear response requires.
!                           Off: a static dipole offset (86% of C(0) on a real
!                           water buffer) leaks broadband noise into the
!                           long-lag end of the ACF.
!     ir_estimator          "biased" divides C(tau) by N, giving a positive
!                           semi-definite sequence whose UNWINDOWED transform
!                           cannot be negative; "unbiased" divides by N - tau
!                           and is not, so it can produce negative absorption.
!                           Windowing can still take a biased estimate negative
!                           unless the kernel is non-negative too: pair it with
!                           ir_window = bartlett for an unconditional guarantee.
!     ir_taper_partial      rebuild the lag window over the lags actually
!                           available while a do_ir ensemble is still filling,
!                           so the header quotes the resolution the transformed
!                           lags support. Off with ir_estimator = unbiased, C is
!                           cut off with a boxcar and the early blocks ring hard
!                           enough to go negative; with the default biased
!                           estimator the implicit (1 - tau/N) taper already
!                           covers that case.
!     ir_weight_by_spacing  trapezoidal weights, so the loss is an integral
!                           over frequency rather than a sum over however the
!                           experimental file happened to be sampled.
!     ir_match_offset       fit an additive baseline alongside the scale, for
!                           experimental data whose baseline was subtracted to
!                           a hard zero.
      logical :: ir_subtract_mean = .true.
      character*32 :: ir_estimator = "biased"
      logical :: ir_taper_partial = .true.
      logical :: ir_weight_by_spacing = .true.
      logical :: ir_match_offset = .false.

!     ==================================================================
!     NEIGHBOUR LISTS AND THE CORE POTENTIAL
!     ==================================================================
      real(dp) :: neighbors_buffer = 0.d0
      real(dp) :: core_pot_cutoff = 1.d10
      real(dp) :: core_pot_buffer = 1.d0

!     ==================================================================
!     MEMORY BUDGET AND SOAP BATCHING
!     ==================================================================
      real(dp) :: max_GBytes_per_process = 1.d0
!     Did the input name max_Gbytes_per_process? gpu_memory_budget_init sizes
!     that keyword from the node, and must not do so over a value someone chose
!     deliberately.
      logical :: max_Gbytes_set = .false.
      real(dp) :: mem_fraction = 0.25d0

!     ==================================================================
!     MOLECULAR DYNAMICS
!     ==================================================================
      integer :: md_nsteps = 1
      real(dp) :: md_step = 1.d0
      character*16 :: optimize = "vv"
      logical :: randomize_velocities = .false.
!     Which distribution randomize_velocities draws from: "uniform" (the historical
!     draw, kept as the default so old trajectories reproduce) or "maxwell".
      character*16 :: velocity_distribution = "uniform"

!     Thermostat
      character*32 :: thermostat = "none"
      real(dp) :: t_beg = 300.d0
      real(dp) :: t_end = 300.d0
      real(dp) :: tau_t = 100.d0
      real(dp), allocatable :: t_hold(:)
      integer :: n_t_hold = 0

!     Barostat
      character*32 :: barostat = "none"
      character*32 :: barostat_sym = "isotropic"
      real(dp) :: p_beg = 1.d0
      real(dp) :: p_end = 1.d0
      real(dp) :: tau_p = 1000.d0
      real(dp) :: gamma_p = 1.d0
      logical :: scale_box = .false.
      real(dp) :: box_scaling_factor(3, 3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3, 3])

!     Structural relaxation (optimize = "gd" and friends)
      real(dp) :: e_tol = 1.d-6
      real(dp) :: f_tol = 0.01d0
      real(dp) :: p_tol = 0.01d0
      real(dp) :: max_opt_step = 0.1d0
      real(dp) :: gd_box_weight = 2.d0

!     Variable time step
      logical :: variable_time_step = .false.
      real(dp) :: target_pos_step
      real(dp) :: tau_dt = 100.d0

!     ==================================================================
!     NESTED SAMPLING
!     ==================================================================
      integer :: n_nested = 0
      real(dp) :: p_nested = 0.d0
      real(dp) :: nested_max_strain = 0.d0
      real(dp) :: nested_max_volume_change = 0.d0
      logical :: scale_box_nested = .false.

!     ==================================================================
!     MONTE CARLO
!     ==================================================================
      integer :: mc_nsteps = 1
      character*32, allocatable :: mc_types(:)
      real(dp), allocatable :: mc_acceptance(:)
      integer :: n_mc_types = 0
      real(dp) :: mc_move_max = 1.d0
      real(dp) :: mc_lnvol_max = 0.01d0
      real(dp) :: mc_min_dist = 0.2d0
      real(dp) :: mc_max_dist = 100000000.d0
      integer :: mc_max_insertion_trials = 500
      logical :: mc_write_xyz = .false.
      logical :: mc_hamiltonian = .false.
      logical :: accessible_volume = .false.
      character*16 :: mc_hybrid_opt = "vv"

!     Grand canonical: chemical potentials per species
      character*8, allocatable :: mc_species(:)
      real(dp), allocatable :: mc_mu(:)
      real(dp), allocatable :: mc_mu_acceptance(:)
      integer :: n_mc_mu = 0

!     Swap moves
      character*8, allocatable :: mc_swaps(:)
      integer, allocatable :: mc_swaps_id(:)
      integer :: n_mc_swaps = 0

!     Relaxation interleaved with the MC walk
      logical :: mc_relax = .false.
      integer :: mc_nrelax = 0
      character*16 :: mc_relax_opt = "gd"
!     Rigid molecules the grand-canonical moves exchange, one entry per chemical
!     potential. An entry with is_molecule false is an ordinary single atom.
      type(mc_molecule), allocatable :: mc_molecules(:)
      character*1024, allocatable :: mc_molecule_files(:)
!     What mc_mu is measured from: "absolute" leaves the acceptance ratio as
!     written, "e0" adds the exchanged object's e0 to it, so that mc_mu is
!     quoted relative to the isolated-species reference energy rather than to
!     zero.
      character*16 :: mc_mu_reference = "absolute"
!     Per chemical potential, filled once at setup: the mass that enters the
!     thermal de Broglie wavelength and the e0 that mc_mu_reference = "e0" adds.
!     For a molecule these are summed over its atoms.
      real(dp), allocatable :: mc_exchange_mass(:)
      real(dp), allocatable :: mc_exchange_e0(:)
      character*32, allocatable :: mc_relax_after(:)
      integer :: n_mc_relax_after = 0

!     Restricting moves to a polyhedron
      integer :: mc_n_planes = 0
      real(dp), allocatable :: mc_planes(:)
      real(dp), allocatable :: mc_max_dist_to_planes(:)
      logical :: mc_planes_restrict_to_polyhedron = .false.

!     Reverse Monte Carlo against experimental data
      logical :: mc_optimize_exp = .false.

!     ==================================================================
!     VAN DER WAALS
!     ==================================================================
      character*32 :: vdw_type = "none"
      real(dp) :: vdw_sr = 0.94d0
      real(dp) :: vdw_d = 20.d0
      real(dp) :: vdw_rcut = 25.d0
      real(dp) :: vdw_buffer = 1.d0
      real(dp) :: vdw_rcut_inner = 0.5d0
      real(dp) :: vdw_buffer_inner = 0.5d0
      real(dp) :: vdw_scs_rcut = 5.d0
      real(dp), allocatable :: vdw_c6_ref(:)
      real(dp), allocatable :: vdw_r0_ref(:)
      real(dp), allocatable :: vdw_alpha0_ref(:)
      logical :: vdw_hirsh_grad = .true.
      logical :: print_vdw_forces = .false.

!     Many-body dispersion
      real(dp) :: vdw_mbd_rcut = 15.d0
      real(dp) :: vdw_mbd_rcut2 = 8.d0
      real(dp) :: vdw_2b_rcut = 15.d0
      real(dp) :: vdw_2b_rcut2 = 8.d0
      real(dp) :: vdw_omega_ref = 1.3d0
      real(dp) :: vdw_loc_rcut = 5.d0
      real(dp) :: vdw_d_mbd = 6.d0
      real(dp) :: vdw_sr_mbd = 0.83d0
      integer :: vdw_mbd_norder = 6
      integer :: vdw_mbd_nfreq = 13
      integer :: mbd_correction_freq = 100
      logical :: vdw_mbd_grad = .true.
      logical :: vdw_mbd_cent_appr = .true.
      logical :: vdw_polynomial = .false.
      real(dp) :: poly_cut_xmin = 3.d0
      real(dp) :: poly_cut_xmax = 10.d0
      logical :: do_nnls = .false.

!     ==================================================================
!     ELECTROSTATICS
!     ==================================================================
      character*32 :: estat_method = "none"
      type(options_estat) :: estat_options
      real(dp) :: estat_rcut = 10.d0
      real(dp) :: estat_rcut_inner = 4.0d0
      real(dp) :: estat_inner_width = 1.d0
      real(dp) :: estat_dsf_alpha = -1.d0
      logical :: print_estat_forces = .false.

!     ==================================================================
!     LOCAL PROPERTIES
!     ==================================================================
      integer :: n_local_properties = 0
      character*1024, allocatable :: compute_local_properties(:)
      logical, allocatable :: write_local_properties(:)

!     ==================================================================
!     EXPERIMENTAL DATA (MAD): SHARED
!     ==================================================================
      integer :: n_exp = 0
      type(exp_data_container), allocatable :: exp_data(:)
      character*32 :: exp_similarity_type = "squared_diff"
      logical :: exp_forces = .false.
      logical :: exp_energies = .true.
      real(dp), allocatable :: exp_energy_scales(:)
      real(dp), allocatable :: exp_energy_scales_initial(:)
      real(dp), allocatable :: exp_energy_scales_final(:)
      integer :: n_moments = 0
      logical :: write_exp = .true.

!     ==================================================================
!     XPS
!     ==================================================================
      type(exp_data_container) :: xps
      real(dp) :: xps_sigma = 0.4d0
      real(dp) :: xps_e_min = 280.0
      real(dp) :: xps_e_max = 300.0
      integer :: xps_n_samples = 200

!     ==================================================================
!     PAIR DISTRIBUTION
!     ==================================================================
      logical :: do_pair_distribution = .false.
      logical :: write_pair_distribution = .false.
      logical :: pair_distribution_partial = .true.
      character*32 :: pair_distribution_output = "pdf"
      integer :: pair_distribution_n_samples = 200
      real(dp) :: pair_distribution_rcut = 4.d0
      real(dp) :: pair_distribution_kde_sigma = 0.d0
      real(dp) :: r_range_min = 1.0
      real(dp) :: r_range_max = 5.d0

!     ==================================================================
!     STRUCTURE FACTOR
!     ==================================================================
      logical :: do_structure_factor = .false.
      logical :: write_structure_factor = .false.
      logical :: structure_factor_from_pdf = .true.
      logical :: structure_factor_window = .true.
      logical :: structure_factor_matrix = .true.
      logical :: structure_factor_matrix_forces = .true.
      character*32 :: sf_output = "xrd"
      integer :: structure_factor_n_samples = 200
      real(dp) :: q_range_min = 1.0
      real(dp) :: q_range_max = 5.d0
      character*32 :: q_units = "q"

!     ==================================================================
!     X-RAY AND NEUTRON DIFFRACTION
!     ==================================================================
      logical :: do_xrd = .false.
      logical :: write_xrd = .false.
      character*32 :: xrd_method = "xrd"
      character*32 :: xrd_output = "xrd"
      integer :: xrd_n_samples = 200
      real(dp) :: xrd_wavelength = 1.5405981d0
      real(dp) :: xrd_rcut = 4.d0
      real(dp) :: xrd_alpha = 1.01d0
      real(dp) :: xrd_damping = 0.0d0
      logical :: do_nd = .false.
      logical :: write_nd = .false.
      character*32 :: nd_output = "xrd"
      integer :: nd_n_samples = 200
      real(dp) :: nd_rcut = 4.d0
!     Compute the XRD/ND pattern from the N^2 Debye sum over atomic positions
!     instead of Fourier-transforming the partial pair distributions. The two
!     routes answer the same question by different approximations: the pdf/sf
!     route bins pair distances and truncates at pair_distribution_rcut, the
!     Debye route sums every pair in the cell exactly. Turning this on makes
!     the pdf and sf calculations unnecessary, and read_files switches off
!     whichever of them was enabled only to feed the XRD.
      logical :: xrd_debye = .false.
!     Multiply the predicted Debye pattern by the powder Lorentz-polarization
!     factor. Only the Debye route applies it, where the energy and the
!     gradient stay consistent because it is one multiplicative weight per q.
      logical :: xrd_lorentz_polarization = .false.
!     Degree of polarization of the incident beam in the Lorentz-polarization
!     factor (1 + P cos^2 2theta) / (sin^2 theta cos theta). P = 1 is an
!     unpolarized source; a graphite monochromator at 2theta_M gives
!     P = cos^2(2theta_M).
      real(dp) :: xrd_lp_polarization = 1.d0
!     Below this |sin(theta)| the Lorentz factor 1/(sin^2 theta cos theta) is
!     unusable, so the whole factor is set to zero there rather than blown up.
      real(dp) :: xrd_lp_sin_theta_min = 1.d-3
!     Deprecated, read and ignored.
      logical :: xrd_iwasa = .true.

!     ==================================================================
!     GPU
!     ==================================================================
!     Accepted by the CPU build too, so that one input deck runs on both.
      logical :: gpu_batched = .true.
      logical :: gpu_low_memory = .true.
      integer :: gpu_n_batches = 1
      integer :: n_batches = 0
      real(dp) :: gpu_max_batch_size = 1.d0

!     ==================================================================
!     OUTPUT
!     ==================================================================
      integer :: write_xyz = 0
      integer :: write_thermo = 1
      logical :: write_soap = .false.
      logical :: write_derivatives = .false.
      logical :: write_lv = .false.
      logical :: write_forces = .true.
      logical :: write_velocities = .true.
      logical :: write_virial = .true.
      logical :: write_pressure = .true.
      logical :: write_stress = .true.
      logical :: write_local_energies = .true.
      logical :: write_masses = .false.
      logical :: write_fixes = .true.
!     write_forces and the rest above are the user-facing switches;
!     these two arrays are what xyz.f90 actually reads, filled from them
!     once read_input_file has finished.
      logical :: write_property(1:11) = .true.
      logical :: write_array_property(1:8) = .true.

!     ==================================================================
!     ADAPTIVE TIME STEP
!     ==================================================================
      logical :: adaptive_time = .false.
      integer :: adapt_tstep_interval = 1
      real(dp) :: adapt_tmin = 1.0d-3
      real(dp) :: adapt_tmax = 1.0d0
      real(dp) :: adapt_xmax = 1.0d-2
      real(dp) :: adapt_emax = 1.0d+1

!     ==================================================================
!     RADIATION CASCADE: ELECTRONIC STOPPING
!     ==================================================================
      logical :: electronic_stopping = .false.
      real(dp) :: eel_cut = 1.0d0
      integer :: eel_freq_out = 1
      character*1024 :: estop_filename = 'NULL'

!     ==================================================================
!     NON-ADIABATIC ENERGY EXCHANGE: EPH MODEL
!     ==================================================================
      logical :: nonadiabatic_processes = .false.
      integer :: model_eph = 1
      integer :: eph_fdm_option = 1
      integer :: eph_friction_option = 1
      integer :: eph_random_option = 1
      integer :: eph_md_last_step = 0
      integer :: eph_freq_Tout = 1
      integer :: eph_freq_mesh_Tout = 1
      integer :: eph_fdm_steps = 1
      integer :: eph_gsx = 1
      integer :: eph_gsy = 1
      integer :: eph_gsz = 1
      real(dp) :: eph_rho_e = 1.0
      real(dp) :: eph_C_e = 1.0
      real(dp) :: eph_kappa_e = 1.0
      real(dp) :: eph_Ti_e = 300.0
      real(dp) :: eph_E_prev_time = 0.0d0
      real(dp) :: eph_md_prev_time = 0.0d0
      real(dp), dimension(6) :: eph_box_limits = (/-100.0, 100.0, -100.0, 100.0, -100.0, 100.0/)
      real(dp) :: in_x0 = -100.0
      real(dp) :: in_x1 = 100.0
      real(dp) :: in_y0 = -100.0
      real(dp) :: in_y1 = 100.0
      real(dp) :: in_z0 = -100.0
      real(dp) :: in_z1 = 100.0
      character*128 :: eph_Tinfile = 'NULL'
      character*128 :: eph_Toutfile = 'NULL'
      character*128 :: eph_betafile = 'NULL'

!     ==================================================================
!     INTERNAL: set by read_files from the keywords above, never read
!     directly from the input file
!     ==================================================================
!     Index into exp_data of each observable, or 0 when that observable
!     has no data file
      integer :: xps_idx
      integer :: xrd_idx
      integer :: pdf_idx
      integer :: sf_idx
      integer :: nd_idx
!     An observable is "valid" when it was both requested and backed by data
      logical :: valid_pdf = .false.
      logical :: valid_sf = .false.
      logical :: valid_xrd = .false.
      logical :: valid_nd = .false.
!     Whether the deck asked for the pdf/sf in as many words. The XRD keywords
!     switch both on as a side effect, and xrd_debye has to be able to tell a
!     side effect it may undo from a request it may not.
      logical :: do_pair_distribution_explicit = .false.
      logical :: do_structure_factor_explicit = .false.

   end type input_parameters

! This is a container for atomic images
   type image
      real(dp), allocatable :: positions(:, :)
      real(dp), allocatable :: positions_prev(:, :)
      real(dp), allocatable :: velocities(:, :)
      real(dp), allocatable :: masses(:)
      real(dp), allocatable :: forces(:, :)
      real(dp), allocatable :: forces_prev(:, :)
      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: local_properties(:, :)
!     Dipole model output. It travels with the image because MC accepts or
!     rejects whole configurations, and a rejected move must not leave the
!     dipole of the configuration that was thrown away.
      real(dp), allocatable :: local_dipoles(:, :)
      real(dp), allocatable :: energies_dipole(:)
      real(dp) :: dipole(1:3) = 0.d0
      real(dp) :: a_box(1:3)
      real(dp) :: b_box(1:3)
      real(dp) :: c_box(1:3)
      real(dp) :: energy
      real(dp) :: e_kin
      real(dp) :: energy_exp
      integer, allocatable :: species(:)
      integer, allocatable :: species_supercell(:)
!     Which inserted molecule each atom belongs to: mc_mol_id is a serial that
!     is unique within the run, mc_mol_mu says which chemical potential it was
!     inserted under. Both are zero for an atom that is not part of a molecule.
!     Rank 0 only -- the energy evaluation never needs them, so they are not
!     broadcast.
      integer, allocatable :: mc_mol_id(:)
      integer, allocatable :: mc_mol_mu(:)
      integer :: n_sites
      integer :: indices(1:3)
      logical, allocatable :: fix_atom(:, :)
      character*8, allocatable :: xyz_species(:)
      character*8, allocatable :: xyz_species_supercell(:)
   end type image

!**************************************************************************
!   One additive contribution family as it crosses the MPI reduce block in
!   turbogap.f90: where its energies, forces and virial are packed from, and
!   where the reduced result is unpacked to.
!
!   Six of the ten families are computed into a this_ prefixed array and land
!   in the un-prefixed one, so src and dst differ; for the other four they are
!   the same array. Holding both ends here is what lets the pack and unpack
!   walks collapse to one loop each, so they cannot disagree about which slot
!   belongs to which family.
!
!   The force and virial pointers are associated only when that family
!   actually carries forces -- for the exp-spectra families that additionally
!   requires exp_forces, and their this_forces_ arrays are not allocated
!   otherwise.
!**************************************************************************
!   Per-run decisions, evaluated once and only read thereafter.
!
!   The rule this exists to enforce: a condition written at N sites is free to
!   disagree with itself, and every time it has, in this codebase, it has been
!   a defect -- ts+mbd, the MPI reduce triple walk, has_vdw against
!   has_local_properties, gap_interface counting < and filling <=, the
!   electrostatics guard, the duplicate keyword handlers, time_gap.  Adding a
!   decision here and reading the field is the structural answer.
!
!   Every field below is a function of params and valid_xps alone, none of
!   which is reassigned inside the main loop, so this is filled once before
!   the loop rather than per iteration.
   type perform_t
!     an experimental observable is both requested and backed by data
      logical :: pdf = .false.
      logical :: sf = .false.
      logical :: xrd = .false.
      logical :: nd = .false.
!     ... and its forces are wanted, so its forces_/virial_ arrays exist
      logical :: pdf_forces = .false.
      logical :: sf_forces = .false.
      logical :: xrd_forces = .false.
      logical :: nd_forces = .false.
      logical :: xps_forces = .false.
   end type perform_t
!**************************************************************************

   type contribution_ref
      real(dp), pointer :: e_src(:) => null()
      real(dp), pointer :: f_src(:, :) => null()
      real(dp), pointer :: v_src(:, :) => null()
      real(dp), pointer :: e_dst(:) => null()
      real(dp), pointer :: f_dst(:, :) => null()
      real(dp), pointer :: v_dst(:, :) => null()
      logical :: forces = .false.
   end type contribution_ref
!**************************************************************************

contains

!**************************************************************************
! This provides a way to pass all the individual arrays/variables in the main code to an image container
!**************************************************************************
!  Whether any soap_turbo descriptor sets a given flag.
!
!  A potential file need not contain a soap_turbo block at all -- a pure 2b, 3b
!  or core_pot potential is a legitimate thing to run, and the finite-difference
!  suite builds exactly those. read_gap_file only allocates soap_turbo_hypers
!  when there is at least one block, so the bare `any(soap_turbo_hypers(:)%flag)`
!  these replace read an unallocated array descriptor and segfaulted before the
!  first energy was computed.
   pure function any_has_local_properties(hypers) result(flag)
      type(soap_turbo), allocatable, intent(in) :: hypers(:)
      logical :: flag

      flag = .false.
      if (allocated(hypers)) flag = any(hypers(:)%has_local_properties)

   end function any_has_local_properties

   pure function any_has_vdw(hypers) result(flag)
      type(soap_turbo), allocatable, intent(in) :: hypers(:)
      logical :: flag

      flag = .false.
      if (allocated(hypers)) flag = any(hypers(:)%has_vdw)

   end function any_has_vdw

   pure function any_has_core_electron_be(hypers) result(flag)
      type(soap_turbo), allocatable, intent(in) :: hypers(:)
      logical :: flag

      flag = .false.
      if (allocated(hypers)) flag = any(hypers(:)%has_core_electron_be)

   end function any_has_core_electron_be
!**************************************************************************

! In time I should make the image data type the default way to store these properties!!!!!!!
   subroutine from_properties_to_image(this_image, positions, velocities, masses, &
                                       forces, a_box, b_box, c_box, energy, energies, energy_exp, e_kin, &
                                       species, species_supercell, n_sites, indices, fix_atom, &
                                       xyz_species, xyz_species_supercell, local_properties, &
                                       local_dipoles, energies_dipole, dipole, mc_mol_id, mc_mol_mu)
      implicit none

!   Input variables
      real(dp), intent(in) :: positions(:, :)
      real(dp), intent(in) :: velocities(:, :)
      real(dp), intent(in) :: masses(:)
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: forces(:, :)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: e_kin
      real(dp), intent(in) :: energy_exp
      real(dp), allocatable, intent(in) :: local_properties(:, :)
      real(dp), allocatable, intent(in), optional :: local_dipoles(:, :)
      real(dp), allocatable, intent(in), optional :: energies_dipole(:)
      real(dp), intent(in), optional :: dipole(1:3)
      integer, intent(in) :: species(:)
      integer, intent(in) :: species_supercell(:)
      integer, intent(in) :: n_sites
      integer, intent(in) :: indices(1:3)
      logical, intent(in) :: fix_atom(:, :)
      character*8, intent(in) :: xyz_species(:)
      character*8, intent(in) :: xyz_species_supercell(:)
!   Molecule bookkeeping. Optional because only a grand-canonical run that
!   exchanges whole molecules has any, and it is rank-0 state: the energy
!   evaluation never reads it, so it is not broadcast.
      integer, allocatable, intent(in), optional :: mc_mol_id(:)
      integer, allocatable, intent(in), optional :: mc_mol_mu(:)
!   In/out variables
      type(image), intent(inout) :: this_image
!   Internal variables
      integer :: n
      integer :: n2

      n = size(positions, 2)
      if (allocated(this_image%positions)) deallocate (this_image%positions)
      allocate (this_image%positions(1:3, 1:n))
      this_image%positions = positions

      n = size(velocities, 2)
      if (allocated(this_image%velocities)) deallocate (this_image%velocities)
      allocate (this_image%velocities(1:3, 1:n))
      this_image%velocities = velocities

      n = size(masses, 1)
      if (allocated(this_image%masses)) deallocate (this_image%masses)
      allocate (this_image%masses(1:n))
      this_image%masses = masses

      n = size(energies, 1)
      if (allocated(this_image%energies)) deallocate (this_image%energies)
      allocate (this_image%energies(1:n))
      this_image%energies = energies

      n = size(forces, 2)
      if (allocated(this_image%forces)) deallocate (this_image%forces)
      allocate (this_image%forces(1:3, 1:n))
      this_image%forces = forces

      this_image%a_box = a_box

      this_image%b_box = b_box

      this_image%c_box = c_box

      this_image%energy = energy

      this_image%energy_exp = energy_exp

      this_image%e_kin = e_kin

      n = size(species, 1)
      if (allocated(this_image%species)) deallocate (this_image%species)
      allocate (this_image%species(1:n))
      this_image%species = species

      n = size(species_supercell, 1)
      if (allocated(this_image%species_supercell)) deallocate (this_image%species_supercell)
      allocate (this_image%species_supercell(1:n))
      this_image%species_supercell = species_supercell

      this_image%n_sites = n_sites

      if (present(mc_mol_id)) then
         if (allocated(mc_mol_id)) then
            n = size(mc_mol_id, 1)
            if (allocated(this_image%mc_mol_id)) deallocate (this_image%mc_mol_id)
            allocate (this_image%mc_mol_id(1:n))
            this_image%mc_mol_id = mc_mol_id
         end if
      end if
      if (present(mc_mol_mu)) then
         if (allocated(mc_mol_mu)) then
            n = size(mc_mol_mu, 1)
            if (allocated(this_image%mc_mol_mu)) deallocate (this_image%mc_mol_mu)
            allocate (this_image%mc_mol_mu(1:n))
            this_image%mc_mol_mu = mc_mol_mu
         end if
      end if

      this_image%indices = indices

      n = size(fix_atom, 2)
      if (allocated(this_image%fix_atom)) deallocate (this_image%fix_atom)
      allocate (this_image%fix_atom(1:3, 1:n))
      this_image%fix_atom = fix_atom

      n = size(xyz_species, 1)
      if (allocated(this_image%xyz_species)) deallocate (this_image%xyz_species)
      allocate (this_image%xyz_species(1:n))
      this_image%xyz_species = xyz_species

      n = size(xyz_species_supercell, 1)
      if (allocated(this_image%xyz_species_supercell)) deallocate (this_image%xyz_species_supercell)
      allocate (this_image%xyz_species_supercell(1:n))
      this_image%xyz_species_supercell = xyz_species_supercell

      if (allocated(local_properties)) then
         n = size(local_properties, 1)
         n2 = size(local_properties, 2)
         if (allocated(this_image%local_properties)) deallocate (this_image%local_properties)
         allocate (this_image%local_properties(1:n, 1:n2))
         this_image%local_properties = local_properties
      end if

      if (present(local_dipoles)) then
         if (allocated(local_dipoles)) then
            n = size(local_dipoles, 2)
            if (allocated(this_image%local_dipoles)) deallocate (this_image%local_dipoles)
            allocate (this_image%local_dipoles(1:3, 1:n))
            this_image%local_dipoles = local_dipoles
         end if
      end if

      if (present(energies_dipole)) then
         if (allocated(energies_dipole)) then
            n = size(energies_dipole, 1)
            if (allocated(this_image%energies_dipole)) deallocate (this_image%energies_dipole)
            allocate (this_image%energies_dipole(1:n))
            this_image%energies_dipole = energies_dipole
         end if
      end if

      if (present(dipole)) this_image%dipole = dipole

   end subroutine
!**************************************************************************

!**************************************************************************
   subroutine from_image_to_properties(this_image, positions, velocities, masses, &
                                       forces, a_box, b_box, c_box, energy, energies, energy_exp, e_kin, &
                                       species, species_supercell, n_sites, indices, fix_atom, &
                                       xyz_species, xyz_species_supercell, local_properties, &
                                       local_dipoles, energies_dipole, dipole, mc_mol_id, mc_mol_mu)
      implicit none

!   Input variables
      type(image), intent(in) :: this_image
!   Output variables
      real(dp), allocatable, intent(out) :: positions(:, :)
      real(dp), allocatable, intent(out) :: velocities(:, :)
      real(dp), allocatable, intent(out) :: masses(:)
      real(dp), allocatable, intent(out) :: forces(:, :)
      real(dp), allocatable, intent(out) :: energies(:)
      real(dp), allocatable, intent(out) :: local_properties(:, :)
      real(dp), allocatable, intent(out), optional :: local_dipoles(:, :)
      real(dp), allocatable, intent(out), optional :: energies_dipole(:)
      real(dp), intent(out), optional :: dipole(1:3)
      real(dp), intent(out) :: a_box(1:3)
      real(dp), intent(out) :: b_box(1:3)
      real(dp), intent(out) :: c_box(1:3)
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: e_kin
      real(dp), intent(out) :: energy_exp
      integer, allocatable, intent(out) :: species(:)
      integer, allocatable, intent(out) :: species_supercell(:)
      integer, intent(out) :: n_sites
      integer, intent(out) :: indices(1:3)
      logical, allocatable, intent(out) :: fix_atom(:, :)
      character*8, allocatable, intent(out) :: xyz_species(:)
      character*8, allocatable, intent(out) :: xyz_species_supercell(:)
!   Molecule bookkeeping; see from_properties_to_image.
      integer, allocatable, intent(inout), optional :: mc_mol_id(:)
      integer, allocatable, intent(inout), optional :: mc_mol_mu(:)
!   Internal variables
      integer :: n
      integer :: n2

      n = size(this_image%positions, 2)
      allocate (positions(1:3, 1:n))
      positions = this_image%positions

      n = size(this_image%velocities, 2)
      allocate (velocities(1:3, 1:n))
      velocities = this_image%velocities

      n = size(this_image%masses, 1)
      allocate (masses(1:n))
      masses = this_image%masses

      n = size(this_image%energies, 1)
      allocate (energies(1:n))
      energies = this_image%energies

      n = size(this_image%forces, 2)
      allocate (forces(1:3, 1:n))
      forces = this_image%forces

      a_box = this_image%a_box

      b_box = this_image%b_box

      c_box = this_image%c_box

      energy_exp = this_image%energy_exp

      energy = this_image%energy

      e_kin = this_image%e_kin

      n = size(this_image%species, 1)
      allocate (species(1:n))
      species = this_image%species

      n = size(this_image%species_supercell, 1)
      allocate (species_supercell(1:n))
      species_supercell = this_image%species_supercell

      n_sites = this_image%n_sites

      if (present(mc_mol_id)) then
         if (allocated(this_image%mc_mol_id)) then
            n = size(this_image%mc_mol_id, 1)
            if (allocated(mc_mol_id)) deallocate (mc_mol_id)
            allocate (mc_mol_id(1:n))
            mc_mol_id = this_image%mc_mol_id
         end if
      end if
      if (present(mc_mol_mu)) then
         if (allocated(this_image%mc_mol_mu)) then
            n = size(this_image%mc_mol_mu, 1)
            if (allocated(mc_mol_mu)) deallocate (mc_mol_mu)
            allocate (mc_mol_mu(1:n))
            mc_mol_mu = this_image%mc_mol_mu
         end if
      end if

      indices = this_image%indices

      n = size(this_image%fix_atom, 2)
      allocate (fix_atom(1:3, 1:n))
      fix_atom = this_image%fix_atom

      n = size(this_image%xyz_species, 1)
      allocate (xyz_species(1:n))
      xyz_species = this_image%xyz_species

      n = size(this_image%xyz_species_supercell, 1)
      allocate (xyz_species_supercell(1:n))
      xyz_species_supercell = this_image%xyz_species_supercell

      if (allocated(this_image%local_properties)) then
         n = size(this_image%local_properties, 1)
         n2 = size(this_image%local_properties, 2)
         allocate (local_properties(1:n, 1:n2))
         local_properties = this_image%local_properties
      end if

      if (present(local_dipoles)) then
         if (allocated(this_image%local_dipoles)) then
            n = size(this_image%local_dipoles, 2)
            allocate (local_dipoles(1:3, 1:n))
            local_dipoles = this_image%local_dipoles
         end if
      end if

      if (present(energies_dipole)) then
         if (allocated(this_image%energies_dipole)) then
            n = size(this_image%energies_dipole, 1)
            allocate (energies_dipole(1:n))
            energies_dipole = this_image%energies_dipole
         end if
      end if

      if (present(dipole)) dipole = this_image%dipole

   end subroutine
!**************************************************************************

end module
