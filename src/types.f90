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
   use iso_c_binding

   implicit none

! dp came from here before src/kinds.f90 existed. It is now use-associated from
! kinds and re-exported, so every module that gets dp through  still
! does, and there is one definition rather than two.

   type gpu_host_storage_type
      real(dp), allocatable :: xyz_k_h(:, :)
      real(dp), allocatable :: pair_distribution_partial_h(:)
      real(dp), allocatable :: pair_distribution_partial_der_h(:, :)
      real(dp), allocatable :: forces_h(:, :)
      real(dp), allocatable :: rjs_index_h(:)
      real(dp) :: virial_h(1:3, 1:3)
      integer, allocatable :: k_index_h(:)
      integer, allocatable :: j2_index_h(:)
   end type gpu_host_storage_type

   type gpu_host_batch_storage_type
      type(gpu_host_storage_type), allocatable :: host(:)
   end type gpu_host_batch_storage_type

   ! each one of these will be allocated 1:n_dim_partial
   type gpu_storage_type
      integer, allocatable :: nk(:)
      type(c_ptr), allocatable :: nk_d(:), k_index_d(:), j2_index_d(:), xyz_k_d(:), &
                pair_distribution_partial_d(:), pair_distribution_partial_der_d(:), nk_flags_sum_d(:), nk_flags_d(:), rjs_index_d(:)
      integer(c_size_t), allocatable :: st_nk_d(:)
      integer(c_size_t), allocatable :: st_k_index_d(:)
      integer(c_size_t), allocatable :: st_j2_index_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_der_d(:)
   end type gpu_storage_type

   ! type gpu_batch_storage_type
   !    type( gpu_storage_type ), allocatable :: device(:)
   ! end type gpu_batch_storage_type

   type gpu_neigh_storage_type
      type(c_ptr) :: n_neigh_d, species_d, neighbors_list_d, neighbor_species_d, rjs_d, xyz_d
      integer(c_size_t) :: st_n_neigh_d
      integer(c_size_t) :: st_species_d
      integer(c_size_t) :: st_neighbors_list_d
      integer(c_size_t) :: st_neighbor_species_d
      integer(c_size_t) :: st_rjs_d
      integer(c_size_t) :: st_xyz_d
   end type gpu_neigh_storage_type

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

   type exp_pred_container
      integer             :: n_samples = 200
      logical             :: write = .false.
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp) :: range_min = 0.d0
      real(dp) :: range_max = 1.d0
   end type exp_pred_container

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
      logical :: zero_trunc = .false.
      type(c_ptr)         :: Qs_d, alphas_d
      integer(c_int)      :: n_sparse_d
      integer(c_size_t) :: st_size_alphas
      integer(c_size_t) :: st_size_Qs
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
      real(dp), allocatable :: vdw_cutoff(:)!
      real(dp), allocatable :: compress_P_el(:)
      real(dp) :: zeta = 2.d0
      real(dp) :: delta = 1.d0
      real(dp) :: rcut_max
      real(dp) :: vdw_zeta
      real(dp) :: vdw_delta
      real(dp) :: vdw_V0
      integer, allocatable :: alpha_max(:)
      integer, allocatable :: compress_soap_indices(:)!
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
      integer :: n_local_properties = 0& !  compress_P_nonzero
                 &, vdw_index = 0, core_electron_be_index = 0
      ! NOTE!! We still have the compress_P_i etc here as we
      ! have not merged properly with the newest version of soap_turbo
      ! which contains support for compression without specifying the
      ! compression indices. We can do this later as it takes time
      ! validating the reimplementation of the GPU kernels.

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
!     A dipole model is a GAP whose "energy" is a fictitious scalar fitted so
!     that mu_i = dE_i/dr_i (the derivative w.r.t. the central atom's *own*
!     position) reproduces the local dipole. Such a descriptor contributes to mu
!     only: its energy, forces and virial are meaningless and are not accumulated.
      logical :: is_dipole_model = .false.
      logical :: recompute_basis = .true.
      type(local_property_soap_turbo), allocatable :: local_property_models(:)
      type(c_ptr) :: W_d, S_d, multiplicity_array_d
      integer(c_size_t) :: st_W_d
      integer(c_size_t) :: st_S_d
      integer(c_size_t) :: st_multiplicity_array_d
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

   type options_estat
      logical :: damped = .false.
      logical :: tsf = .true.
      logical :: sp = .false.
      logical :: gsf = .false.
      logical :: self_energy_correction = .false.
      logical :: damped_cosine = .false.
   end type options_estat

! These is the type for the input parameters
   type input_parameters

      ! General Parameters
      character*1024 :: atoms_file
      character*8, allocatable :: species_types(:)
      real(dp), allocatable :: masses_types(:)
      real(dp), allocatable :: e0(:)
      logical :: do_prediction = .false.
      logical :: do_forces = .false.
      logical :: do_derivatives = .false.
      logical :: do_derivatives_fd = .false.
      integer :: which_atom = 0

      real(dp), allocatable :: mc_acceptance(:)
      real(dp), allocatable :: mc_mu_acceptance(:)
      real(dp), allocatable :: mc_mu(:)
      real(dp), allocatable :: radii(:)
      real(dp), allocatable :: t_hold(:)

      ! State Variables / Thermo/barostat parameters
      ! Temperature
      real(dp) :: t_beg = 300.d0
      real(dp) :: t_end = 300.d0
      real(dp) :: tau_t = 100.d0

      ! Pressure
      real(dp) :: p_beg = 1.d0
      real(dp) :: p_end = 1.d0
      real(dp) :: tau_p = 1000.d0
      real(dp) :: gamma_p = 1.d0
      logical :: scale_box = .false.
      real(dp) :: box_scaling_factor(3, 3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3, 3])

      ! General MD parameters
      logical :: do_md = .false.
      integer :: md_nsteps = 1
      real(dp) :: md_step = 1.d0
      character*16 :: optimize = "vv"
      character*32 :: barostat = "none"
      character*32 :: thermostat = "none"
      character*32 :: barostat_sym = "isotropic"
      ! Gradient Descent
      real(dp) :: e_tol = 1.d-6
      real(dp) :: gamma0 = 0.01d0
      real(dp) :: max_opt_step = 0.1d0
      real(dp) :: f_tol = 0.01d0
      real(dp) :: p_tol = 0.01d0
      real(dp) :: max_opt_step_eps = 0.05d0

      ! Neighbors parameters
      real(dp) :: neighbors_buffer = 0.d0
      real(dp) :: core_pot_cutoff = 1.d10
      real(dp) :: core_pot_buffer = 1.d0

      ! SOAP calculation batching
      real(dp) :: max_GBytes_per_process = 1.d0

      ! Variable time step
      real(dp) :: tau_dt = 100.d0
      real(dp) :: target_pos_step
      logical :: variable_time_step = .false.

      ! VDW parameters
      character*32 :: vdw_type = "none"
      real(dp) :: vdw_sr = 0.94d0
      real(dp) :: vdw_d = 20.d0
      real(dp) :: vdw_rcut = 10.d0
      real(dp) :: vdw_buffer = 1.d0
      real(dp) :: vdw_rcut_inner = 0.5d0
      real(dp) :: vdw_buffer_inner = 0.5d0
      real(dp) :: vdw_scs_rcut = 4.d0
      logical :: print_vdw_forces = .false.
      real(dp), allocatable :: vdw_c6_ref(:)
      real(dp), allocatable :: vdw_r0_ref(:)
      real(dp), allocatable :: vdw_alpha0_ref(:)

      ! MBD-VDW parameters
      real(dp) :: poly_cut_xmin = 3.d0
      real(dp) :: poly_cut_xmax = 10.d0
      real(dp) :: vdw_mbd_rcut = 15.d0
      real(dp) :: vdw_mbd_rcut2 = 8.d0
      real(dp) :: vdw_2b_rcut = 15.d0
      real(dp) :: vdw_2b_rcut2 = 8.d0
      real(dp) :: vdw_omega_ref = 1.3d0
      real(dp) :: vdw_loc_rcut = 5.d0
      real(dp) :: vdw_d_mbd = 6.d0
      real(dp) :: vdw_sr_mbd = 0.83d0
      logical :: vdw_hirsh_grad = .true.
      logical :: vdw_polynomial = .false.
      logical :: do_nnls = .false.
      logical :: vdw_mbd_cent_appr = .true.
      integer :: mbd_correction_freq = 100
      integer :: vdw_mbd_norder = 6
      integer :: vdw_mbd_nfreq = 11
      logical :: vdw_mbd_grad = .false.

      ! Monte-Carlo
      logical :: do_mc = .false.
      integer :: mc_nsteps = 1
      integer :: n_mc_types = 0
      integer :: n_mc_mu = 0
      integer :: n_mc_swaps = 0
      integer :: mc_nrelax = 0
      integer :: n_mc_relax_after = 0
      real(dp) :: mc_move_max = 1.d0
      real(dp) :: mc_lnvol_max = 0.01d0
      real(dp) :: mc_min_dist = 0.2d0
      character*8, allocatable :: mc_swaps(:)
      character*8, allocatable :: mc_species(:)
      character*32, allocatable :: mc_types(:)
      character*32, allocatable :: mc_relax_after(:)
      character*16 :: mc_relax_opt = "gd"
      character*16 :: mc_hybrid_opt = "vv"
      logical :: mc_write_xyz = .false.
      logical :: mc_relax = .false.
      logical :: mc_hamiltonian = .false.
      logical :: accessible_volume = .false.
      ! Reverse-Monte Carlo
      logical :: mc_optimize_exp = .false.
      ! MC internal
      integer :: mc_idx = 1
      integer, allocatable :: mc_swaps_id(:)

      ! Nested Sampling
      logical :: do_nested_sampling = .false.
      integer :: n_nested = 0
      real(dp) :: p_nested = 0.d0
      real(dp) :: nested_max_strain = 0.d0
      real(dp) :: nested_max_volume_change = 0.d0
      logical :: scale_box_nested = .false.

      ! XPS parameters
      logical :: do_xps = .false.
      real(dp) :: xps_sigma = 0.4d0
      real(dp) :: xps_e_min = 280.0
      real(dp) :: xps_e_max = 300.0
      integer :: xps_n_samples = 200
      character*32 :: xps_force_type = "similarity"

      ! XRD/ND/PDF parameters
      logical :: pair_distribution_partial = .true.
      logical :: structure_factor_from_pdf = .true.
      logical :: structure_factor_window = .true.
      logical :: write_pair_distribution = .false.
      logical :: write_structure_factor = .false.
      logical :: do_pair_distribution = .false.
      logical :: do_structure_factor = .false.
      logical :: do_xrd = .false.
      logical :: do_nd = .false.
      logical :: write_xrd = .false.
      logical :: write_nd = .false.
      logical :: structure_factor_matrix = .true.
      logical :: structure_factor_matrix_forces = .true.
      real(dp) :: xrd_wavelength = 1.5405981d0
      real(dp) :: xrd_damping = 0.0d0
      real(dp) :: nd_wavelength = 1.5405981d0
      real(dp) :: xrd_alpha = 1.01d0
      real(dp) :: xrd_rcut = 4.d0
      real(dp) :: nd_rcut = 4.d0
      real(dp) :: q_range_min = 1.0
      real(dp) :: q_range_max = 5.d0
      integer :: pair_distribution_n_samples = 200
      integer :: structure_factor_n_samples = 200
      integer :: xrd_n_samples = 200
      integer :: nd_n_samples = 200
      character*32 :: q_units = "q"
      character*32 :: sf_output = "xrd"
      character*32 :: nd_output = "xrd"
      character*32 :: pair_distribution_output = "pdf"
      type(exp_data_container), allocatable :: exp_data(:)
      type(exp_pred_container) :: pair_distribution_params
      type(exp_pred_container) :: structure_factor_params
      type(exp_pred_container) :: xrd_params
      logical :: write_exp = .true.
      logical :: valid_pdf = .false.
      logical :: valid_sf = .false.
      logical :: valid_xrd = .false.
      logical :: valid_nd = .false.
      character*32 :: xrd_method = "xrd"
      character*32 :: xrd_output = "xrd"
      ! Depreceated
      logical :: xrd_iwasa = .true.
      ! Internal
      integer :: xps_idx
      integer :: xrd_idx
      integer :: saxs_idx
      integer :: pdf_idx
      integer :: sf_idx
      integer :: nd_idx

      ! PDF parameters
      real(dp) :: r_range_min = 1.0
      real(dp) :: r_range_max = 5.d0
      real(dp) :: pair_distribution_rcut = 4.d0
      real(dp) :: pair_distribution_kde_sigma = 0.d0

      ! MAD parameters
      logical :: do_exp = .false.
      integer :: n_exp = 0
      character*32 :: exp_similarity_type = "squared_diff"
      character*32 :: estat_method = "none"
      logical :: exp_forces = .false.
      logical :: exp_energies = .true.
      real(dp), allocatable :: exp_energy_scales(:)
      real(dp), allocatable :: exp_energy_scales_initial(:)
      real(dp), allocatable :: exp_energy_scales_final(:)

      ! Electrostatics Parameters
      real(dp) :: estat_rcut = 10.d0
      real(dp) :: estat_dsf_alpha = -1.d0
      real(dp) :: estat_rcut_inner = 4.0d0
      real(dp) :: estat_inner_width = 1.d0

      ! GPU Parameters
      logical :: gpu_batched = .true.
      integer :: gpu_n_batches = 1
      integer :: n_batches = 0
      logical :: gpu_low_memory = .true.
      real(dp) :: gpu_mem_fraction = 0.d0
      real(dp) :: gpu_max_batch_size = 1.d0

      ! I/O parameters
      integer :: write_xyz = 0
      integer :: write_thermo = 1

      ! Verbosity for debugging
      integer :: verb = 0

      ! Internal
      integer :: n_local_properties = 0
      integer :: n_moments = 0
      integer :: n_t_hold = 0
      integer :: n_exp_opt = 0
!   Seed for the intrinsic pseudo-random number generator. Zero (the default)
!   leaves the compiler's default sequence untouched; any other value makes
!   runs reproducible, which is what the CPU/GPU regression comparisons rely
!   on when initial velocities have to be randomized.
      integer :: random_seed_value = 0

      character*1024, allocatable :: compute_local_properties(:)

!     Set from the .gap file, not the input file: true when at least one
!     soap_turbo descriptor is flagged dipole_model
      logical :: do_dipole = .false.

      logical :: write_soap = .false.
      logical :: write_derivatives = .false.
      logical :: do_timing = .false.
      logical :: all_atoms = .true.
      logical :: print_progress = .true.
      logical :: write_lv = .false.
      logical :: write_forces = .true.
      logical :: write_velocities = .true.
      logical :: write_hirshfeld_v = .true.
      logical :: write_virial = .true.
      logical :: write_pressure = .true.
      logical :: write_stress = .true.
      logical :: write_local_energies = .true.
      logical :: write_property(1:11) = .true.
      logical :: write_array_property(1:8) = .true.
      logical :: write_masses = .false.
      logical :: write_fixes = .true.

      logical :: print_lp_forces = .false.
      logical :: print_estat_forces = .false.

      type(options_estat) :: estat_options

      logical, allocatable :: write_local_properties(:)

      ! ------- option for doing simulation with adaptive time step
      logical :: adaptive_time = .false.
      integer :: adapt_tstep_interval = 1
      real(dp) :: adapt_tmin = 1.0d-3
      real(dp) :: adapt_tmax = 1.0d0
      real(dp) :: adapt_xmax = 1.0d-2
      real(dp) :: adapt_emax = 1.0d+1
      ! ----------------------------------------------                ******** until here for adaptive time

      ! ------- option for radiation cascade simulation with electronic stopping
      logical :: electronic_stopping = .false.
      real(dp) :: eel_cut = 1.0d0
      integer :: eel_freq_out = 1
      character*1024 :: estop_filename = 'NULL'
      ! ----------------------------------------------                ******** until here for electronic stopping

      ! ------- option for non-adiabatic processes of energy exchange through EPH model

      logical :: nonadiabatic_processes = .false.
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
      integer :: model_eph = 1
      real(dp) :: eph_rho_e = 1.0
      real(dp) :: eph_C_e = 1.0
      real(dp) :: eph_kappa_e = 1.0
      real(dp) :: eph_Ti_e = 300.0
      real(dp) :: in_x0 = -100.0
      real(dp) :: in_x1 = 100.0
      real(dp) :: in_y0 = -100.0
      real(dp) :: in_y1 = 100.0
      real(dp) :: in_z0 = -100.0
      real(dp) :: in_z1 = 100.0
      real(dp) :: eph_E_prev_time = 0.0d0
      real(dp) :: eph_md_prev_time = 0.0d0
      real(dp), dimension(6) :: eph_box_limits = (/-100.0, 100.0, -100.0, 100.0, -100.0, 100.0/)
      character*128 :: eph_Tinfile = 'NULL'
      character*128 :: eph_Toutfile = 'NULL'
      character*128 :: eph_betafile = 'NULL'

      ! ---------------------------------------------                ******** until here for electronic stopping based on EPH model

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
!   Seven of the eleven families are computed into a this_ prefixed array and
!   land in the un-prefixed one, so src and dst differ; for the other four they
!   are the same array. Holding both ends here is what lets the pack and unpack
!   walks collapse to one loop each, so they cannot disagree about which slot
!   belongs to which family.
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
      logical         :: forces = .false.
   end type contribution_ref
!**************************************************************************

contains

!**************************************************************************
! This provides a way to pass all the individual arrays/variables in the main code to an image container
! In time I should make the image data type the default way to store these properties!!!!!!!
   subroutine from_properties_to_image(this_image, positions, velocities, masses, &
                                       forces, a_box, b_box, c_box, energy, energies, energy_exp, e_kin, &
                                       species, species_supercell, n_sites, indices, fix_atom, &
                                       xyz_species, xyz_species_supercell, local_properties, &
                                       local_dipoles, energies_dipole, dipole)
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
                                       local_dipoles, energies_dipole, dipole)
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
