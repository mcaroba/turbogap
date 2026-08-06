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
   use turbogap_setup
   use turbogap_exp
   use turbogap_md
   use turbogap_estat
   use turbogap_vdw
   use gap_backend
   use neighbors
   use soap_turbo_desc
   use gap
   use read_files
   use md
   !  use adaptive_time                        ! for adaptive time simulation (TurboGAP will use these five modules for radiation cascades)
   !  use electronic_stopping                ! for electronic stopping correction in radiation cascades
   !  use eph_fdm                                ! for T - dependent parameters - elec. stop. - eph model
   !  use eph_beta                                ! for the atomic electronic densities  - elec. stop. - eph model
   !  use eph_electronic_stopping                ! for electronic stopping based in radiation cascades on the eph model
   use mc
   use gap_interface
   use types
   use vdw
   use electrostatics, only: compute_coulomb_direct, compute_coulomb_dsf, compute_coulomb_lamichhane, &
                             calculate_batched_electrostatics
   use exp_utils
   use exp_interface
   use soap_turbo_functions
   use F_B_C
   use iso_c_binding
   use gpu_context

#ifdef _MPIF90
   use mpi
   use mpi_helper
#endif
   use bussi
   use xyz_module

!$ use omp_lib

   implicit none

   !**************************************************************************
   ! Variable definitions
   !
   real(dp), allocatable, target :: rjs(:)
   real(dp), allocatable, target :: thetas(:)
   real(dp), allocatable, target :: phis(:)
   real(dp), allocatable, target :: xyz(:, :)
   real(dp), allocatable :: positions(:, :)
   real(dp), allocatable :: positions_prev(:, :)
   real(dp), allocatable :: soap(:, :)
   real(dp), allocatable :: soap_cart_der(:, :, :)
   real(dp), allocatable :: positions_diff(:, :)
   real(dp), allocatable :: forces_prev(:, :)
   real(dp), allocatable :: frac_positions(:, :)
   real(dp) :: rcut_max
   real(dp) :: a_box(1:3)
   real(dp) :: b_box(1:3)
   real(dp) :: c_box(1:3)
   real(dp) :: max_displacement
   real(dp) :: energy
   real(dp) :: energy_prev
   real(dp), target :: virial(1:3, 1:3)
   real(dp), target :: this_virial(1:3, 1:3)
   real(dp), target :: virial_soap(1:3, 1:3)
   real(dp), target :: virial_2b(1:3, 1:3)
   real(dp), target :: virial_3b(1:3, 1:3)
   real(dp), target :: virial_core_pot(1:3, 1:3)
   real(dp), target :: virial_vdw(1:3, 1:3)
   real(dp), target :: virial_lp(1:3, 1:3)
   real(dp), target :: this_virial_vdw(1:3, 1:3)
   real(dp), target :: this_virial_lp(1:3, 1:3)
   real(dp), target :: virial_pdf(1:3, 1:3)
   real(dp), target :: this_virial_pdf(1:3, 1:3)
   real(dp), target :: v_uc
   real(dp), target :: virial_sf(1:3, 1:3)
   real(dp), target :: this_virial_sf(1:3, 1:3)
   real(dp), target :: virial_xrd(1:3, 1:3)
   real(dp), target :: this_virial_xrd(1:3, 1:3)
   real(dp), target :: virial_nd(1:3, 1:3)
   real(dp), target :: this_virial_nd(1:3, 1:3)
   real(dp), target :: virial_estat(1:3, 1:3)
   real(dp), target :: this_virial_estat(1:3, 1:3)
   real(dp) :: v_uc_prev
   real(dp) :: v_a_uc
   real(dp) :: v_a_uc_prev
   real(dp) :: eVperA3tobar = 1602176.6208d0
   real(dp) :: ranf
   real(dp) :: ranv(1:3)
   real(dp) :: disp(1:3)
   real(dp) :: d_disp
   real(dp) :: e_mc_prev
   real(dp) :: p_accept
   real(dp) :: virial_prev(1:3, 1:3)
   real(dp) :: sim_exp_pred
   real(dp) :: sim_exp_prev
   real(dp) :: sim_exp_pred_der(1:3)
   real(dp), allocatable, target :: energies(:)
   real(dp), allocatable, target :: forces(:, :)
   real(dp), allocatable, target :: energies_soap(:)
   real(dp), allocatable, target :: forces_soap(:, :)
   real(dp), allocatable, target :: this_energies(:)
   real(dp), allocatable, target :: this_forces(:, :)
   real(dp), allocatable, target :: energies_2b(:)
   real(dp), allocatable, target :: forces_2b(:, :)
   real(dp), allocatable, target :: energies_3b(:)
   real(dp), allocatable, target :: forces_3b(:, :)
   real(dp), allocatable, target :: energies_core_pot(:)
   real(dp), allocatable, target :: forces_core_pot(:, :)
   real(dp), allocatable, target :: velocities(:, :)
   real(dp), allocatable, target :: masses_types(:)
   real(dp), allocatable, target :: masses(:)
   real(dp), allocatable, target :: hirshfeld_v_temp(:)
   real(dp), allocatable, target :: masses_temp(:)
   real(dp), allocatable, target :: energies_exp(:)
   !  real(dp), allocatable :: this_hirshfeld_v(:), this_hirshfeld_v_cart_der(:,:)
   !  real(dp), pointer :: this_hirshfeld_v_pt(:), this_hirshfeld_v_cart_der_pt(:,:)

   real(dp), allocatable, target :: local_properties(:, :)
   real(dp), allocatable, target :: local_properties_cart_der(:, :, :)
   ! Have one rank lower for the pointer, such that it just relates to a sub array of the local properties/cart_der
   real(dp), pointer :: local_properties_pt(:)
   real(dp), pointer :: local_properties_cart_der_pt(:, :)
   !  real(dp), pointer :: hirshfeld_v(:), hirshfeld_v_cart_der(:,:)
   real(dp), allocatable, target :: this_local_properties(:, :)
   real(dp), allocatable, target :: this_local_properties_cart_der(:, :, :)
   real(dp), pointer :: this_local_properties_pt(:, :)
   real(dp), pointer :: this_local_properties_cart_der_pt(:, :, :)
   real(dp), allocatable :: y_i_pred_all(:, :)
   real(dp), allocatable :: moments(:)
   real(dp), allocatable :: moments_exp(:)

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
   real(dp) :: charge_sum
   real(dp) :: time1
   real(dp) :: time2
   real(dp) :: time3
!   Every wall-clock bucket lives in one times_t (src/timing.f90), so the
!   extracted modules take a single argument instead of thirteen and the two
!   branches' signatures agree. time_step and time_step_prev below are the MD
!   integration step in fs, not timers, and deliberately stay separate.
   type(times_t) :: time
   real(dp) :: instant_pressure
   real(dp) :: lv(1:3, 1:3)
   real(dp) :: instant_pressure_tensor(1:3, 1:3)
   real(dp) :: time_step
   real(dp) :: md_time
   real(dp) :: t1
   real(dp) :: instant_pressure_prev
   real(dp) :: wfac
   real(dp) :: wfac_temp
   real(dp) :: energy_exp
   integer, allocatable :: displs(:)
   integer, allocatable :: displs2(:)
   integer, allocatable :: counts(:)
   integer, allocatable :: counts2(:)
   integer, allocatable :: in_to_out_pairs(:)
   integer, allocatable :: in_to_out_site(:)
   integer, allocatable :: mc_id(:)
   integer :: update_bar
   integer :: n_sparse
   integer :: idx
   integer :: gd_istep = 0
   integer :: nprop
   integer :: n_pairs_temp
   integer :: n_sites_temp
   logical, allocatable :: do_list(:)
   logical, allocatable :: has_local_properties_mpi(:)
   logical, allocatable :: fix_atom(:, :)
   logical :: rebuild_neighbors_list = .true.
   logical :: exit_loop = .true.
   logical :: gd_box_do_pos = .true.
   logical :: restart_box_optim = .false.
   logical :: valid_xps = .false.
   logical :: valid_vdw = .false.
   logical :: valid_estat_charges = .false.
   logical :: write_condition = .false.
   logical :: overwrite_condition = .false.
   logical :: do_electrostatics = .true.

   character*1 :: creturn = achar(13)

   !--- TODO: Create variables for GOU allocation of variables which are needed! ---!

  !! these decalarations are for time step and electronic stopping by different methods
   real(dp) :: time_step_prev
   integer :: nrows
   real(dp) :: cum_EEL = 0.0d0
   real(dp), allocatable :: allelstopdata(:)
   ! type (EPH_Beta_class) :: ephbeta
   ! type (EPH_FDM_class) :: ephfdm
   ! type (EPH_LangevinSpatialCorrelation_class) :: ephlsc

   ! Clean up these variables after code refactoring !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer, allocatable, target :: n_neigh(:)
   integer, allocatable, target :: neighbors_list(:)
   integer, allocatable, target :: alpha_max(:)
   integer, allocatable, target :: species(:)
   integer, allocatable, target :: species_supercell(:)
   integer, allocatable, target :: neighbor_species(:)
   integer, allocatable, target :: der_neighbors(:)
   integer, allocatable, target :: der_neighbors_list(:)
   integer, allocatable, target :: i_beg_list(:)
   integer, allocatable, target :: i_end_list(:)
   integer, allocatable, target :: j_beg_list(:)
   integer, allocatable, target :: j_end_list(:)
   integer, allocatable, target :: species_idx(:)
   integer, allocatable, target :: n_neigh_out(:)
   integer, allocatable, target :: n_local_properties_mpi(:)
   integer, allocatable, target :: local_property_indexes(:)
   integer, allocatable, target :: n_mc_species(:)
   integer, allocatable, target :: n_mc_species_prev(:)
   integer, allocatable, target :: kappas(:)
   integer :: n_sites
   integer :: i
   integer :: j
   integer :: k
   integer :: i2
   integer :: j2
   integer :: n_soap
   integer :: k2
   integer :: k3
   integer :: l
   integer :: n_sites_this
   integer :: ierr
   integer :: rank
   integer :: ntasks
   integer :: dim
   integer :: n_sp
   integer :: n_pos
   integer :: n_sp_sc
   integer :: this_i_beg
   integer :: this_i_end
   integer :: this_j_beg
   integer :: this_j_end
   integer :: this_n_sites_mpi
   integer :: n_sites_prev = 0
   integer :: n_atom_pairs_by_rank_prev = 0
   integer :: cPnz
   integer :: n_pairs
   integer :: n_all_sites
   integer :: n_sites_out
   integer :: n_local_properties_tot = 0
   integer :: n_lp_count = 0
   integer :: charge_lp_index
   integer :: vdw_lp_index
   integer :: core_be_lp_index
   integer :: xps_idx

   integer :: l_max
   integer :: n_atom_pairs
   integer :: n_max
   integer :: ijunk
   integer :: central_species = 0
   integer :: n_atom_pairs_total
   integer :: iostatus
   integer :: counter = 0
   integer :: counter2

!   The eleven additive contribution families reduced together after the
!   descriptor loop. Their predicates used to be written out three times -- to
!   count the slots, to pack them and to unpack them -- and evaluated
!   independently each time. Two copies disagreeing shifts counter2 and
!   silently attributes one family's energies to another.
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
   integer :: n_species = 1
   integer :: n_xyz
   integer :: indices(1:3)
   integer :: radial_enhancement = 0
   integer :: md_istep
   integer :: mc_istep
   integer :: mc_mu_id = 1
   integer :: n_mc
   character*8, allocatable, target :: species_types_actual(:)
   character*1024, allocatable :: local_property_labels(:)
   character*1024, allocatable :: local_property_labels_temp(:)
   character*1024, allocatable :: local_property_labels_temp2(:)
   logical :: repeat_xyz = .true.
   logical :: overwrite = .false.
   logical :: check_species
   logical :: valid_local_properties = .false.
   logical :: label_in_list
   logical :: do_mc_relax = .false.

   character*1024 :: filename
   character*1024 :: cjunk
   character*1024 :: file_compress_soap
   character*1024 :: file_alphas
   character*1024 :: file_soap
   character*1024 :: file_2b
   character*1024 :: file_alphas_2b
   character*1024 :: file_3b
   character*1024 :: file_alphas_3b
   character*1024 :: file_gap = "none"
   character*1024 :: mc_file = "mc_trial.xyz"
   character*1024 :: string
   character*1024 :: temp_string
   character*1024 :: temp_string2
   character*64 :: keyword
   character*16 :: lattice_string(1:9)
   character*8 :: i_char
   character*8, allocatable :: xyz_species(:)
   character*8, allocatable :: xyz_species_supercell(:)
   character*8, allocatable :: species_type_temp(:)

   character*1 :: keyword_first

   ! This is the mode in which we run TurboGAP
   character*16 :: mode = "none"
   character*32 :: mc_move = "none"
   character*32 :: exp_output = "none"

   ! Here we store the input parameters
   type(input_parameters) :: params

   ! These are the containers for the hyperparameters of descriptors and GAPs
   integer :: n_soap_turbo = 0
   integer :: n_distance_2b = 0
   integer :: n_angle_3b = 0
   integer :: n_core_pot = 0
   integer :: counter_lp_names = 0
   integer :: temp_md_nsteps
   real(dp), parameter :: pi = acos(-1.0)
   type(soap_turbo), allocatable, target :: soap_turbo_hypers(:)
   type(distance_2b), allocatable, target :: distance_2b_hypers(:)
   type(angle_3b), allocatable, target :: angle_3b_hypers(:)
   type(core_pot), allocatable, target :: core_pot_hypers(:)

   !vdw crap
   real(dp), allocatable, target :: v_neigh_vdw(:)
!   Persistent ts+mbd correction state, owned by turbogap_vdw.
   type(vdw_state) :: vdw_ws
   real(dp), allocatable :: energies_vdw_corr(:)
   real(dp), allocatable :: forces_vdw_corr(:, :)
   real(dp), allocatable :: local_virial_vdw_diag_corr(:, :)
   logical :: update_mbd_ts_scaling = .true.
   real(dp), allocatable, target :: energies_vdw(:)
!  The TS scaling factor and the diagonal local virial that master's
!  get_ts_energy_and_forces takes. For plain TS the scaling is 1.d0; it
!  exists so the ts+mbd scheme can reuse the same routine.
   real(dp), allocatable :: mbd_ts_scaling(:), this_mbd_ts_scaling(:)
   real(dp), allocatable :: local_virial_vdw_diag(:, :), this_local_virial_vdw_diag(:, :)
   real(dp), allocatable, target :: forces_vdw(:, :)
   real(dp), allocatable, target :: this_energies_vdw(:)
   real(dp), allocatable, target :: this_forces_vdw(:, :)
   real(dp), allocatable, target :: v_neigh_lp(:)
   real(dp), allocatable, target :: energies_lp(:)
   real(dp), allocatable, target :: forces_lp(:, :)
   real(dp), allocatable, target :: this_energies_lp(:)
   real(dp), allocatable, target :: this_forces_lp(:, :)
   real(dp), allocatable, target :: chg_neigh_estat(:)
   real(dp), allocatable, target :: energies_estat(:)
   real(dp), allocatable, target :: forces_estat(:, :)
   real(dp), allocatable, target :: this_energies_estat(:)
   real(dp), allocatable, target :: this_forces_estat(:, :)
   real(dp), allocatable, target :: energies_pdf(:)
   real(dp), allocatable, target :: forces_pdf(:, :)
   real(dp), allocatable, target :: this_energies_pdf(:)
   real(dp), allocatable, target :: this_forces_pdf(:, :)
   real(dp), allocatable, target :: energies_sf(:)
   real(dp), allocatable, target :: forces_sf(:, :)
   real(dp), allocatable, target :: this_energies_sf(:)
   real(dp), allocatable, target :: this_forces_sf(:, :)
   real(dp), allocatable, target :: energies_xrd(:)
   real(dp), allocatable, target :: forces_xrd(:, :)
   real(dp), allocatable, target :: this_energies_xrd(:)
   real(dp), allocatable, target :: this_forces_xrd(:, :)
   real(dp), allocatable, target :: energies_nd(:)
   real(dp), allocatable, target :: forces_nd(:, :)
   real(dp), allocatable, target :: this_energies_nd(:)
   real(dp), allocatable, target :: this_forces_nd(:, :)
   ! MPI stuff
   real(dp), allocatable, target :: temp_1d(:)
   real(dp), allocatable, target :: temp_1d_bis(:)
   real(dp), allocatable, target :: temp_2d(:, :)
   real(dp), allocatable, target :: structure_factor_partial_der(:, :, :)
   real(dp), allocatable, target :: structure_factor_partial_temp_der(:, :)
   real(dp), allocatable, target :: y_xrd_der(:, :, :)
   real(dp), allocatable, target :: y_xrd_temp_der(:, :, :)
   real(dp), allocatable, target :: y_nd_der(:, :, :)
   real(dp), allocatable, target :: y_nd_temp_der(:, :, :)
   real(dp), allocatable, target :: all_scattering_factors(:)
   integer, allocatable :: temp_1d_int(:)
   integer, allocatable :: n_atom_pairs_by_rank(:)
   integer, allocatable :: site_in_rank(:)
   integer, allocatable :: this_site_in_rank(:)
   integer, allocatable :: n_species_mpi(:)
   integer, allocatable :: n_sparse_mpi_soap_turbo(:)
   integer, allocatable :: dim_mpi(:)
   integer, allocatable :: n_sparse_mpi_distance_2b(:)
   integer, allocatable :: n_sparse_mpi_angle_3b(:)
   integer, allocatable :: n_mpi_core_pot(:)
   integer, allocatable :: local_properties_n_sparse_mpi_soap_turbo(:)
   integer, allocatable :: local_properties_dim_mpi_soap_turbo(:)
   integer, allocatable :: n_neigh_local(:)
   integer, allocatable :: compress_P_nonzero_mpi(:)
   integer :: i_beg
   integer :: i_end
   integer :: n_sites_mpi
   integer :: j_beg
   integer :: j_end
   integer :: size_soap_turbo
   integer :: size_distance_2b
   integer :: size_angle_3b
   integer :: n_nonzero
   logical, allocatable :: compress_soap_mpi(:)

   !--- GPU VARIABLES FOR ALLOCATION ---!

   integer :: n_omp
   integer :: omp_task
   integer :: omp_n_sites
   integer :: n_omp_temp
   integer, allocatable :: i_beg_omp(:)
   integer, allocatable :: i_end_omp(:)
   integer, allocatable :: j_beg_omp(:)
   integer, allocatable :: j_end_omp(:)
   !**************************************************************************
   integer :: sp1
   integer :: sp2
   logical(c_bool) :: c_do_forces
   integer(c_size_t) :: st_n_sites_int
   integer(c_size_t) :: st_n_sites_double
   integer(c_size_t) :: st_n_atom_pairs_int
   integer(c_size_t) :: st_n_atom_pairs_double
   integer(c_size_t) :: st_n_sparse_double
   integer(c_size_t) :: st_virial
   integer(c_size_t) :: st_size_nf
   integer(c_size_t) :: st_size_rcut_hard
   integer(c_size_t) :: st_species_types_actual_d
   type(c_ptr) :: n_neigh_d
   type(c_ptr) :: species_d
   type(c_ptr) :: neighbor_species_d
   type(c_ptr) :: rjs_d
   type(c_ptr) :: alphas_d
   type(c_ptr) :: cutoff_d
   type(c_ptr) :: qs_d
   type(c_ptr) :: xyz_d
   type(c_ptr) :: species_types_actual_d
   type(c_ptr) :: energies_2b_d
   type(c_ptr) :: forces_2b_d
   type(c_ptr) :: virial_2b_d
   type(c_ptr) :: x_d
   type(c_ptr) :: V_d
   type(c_ptr) :: dVdx2_d
   type(c_ptr) :: forces_core_pot_d
   type(c_ptr) :: virial_core_pot_d
   type(c_ptr) :: energies_core_pot_d
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
   !**************************************************************************
   !local variables for 3benergy and forces gpu
   character(kind=c_char, len=4) :: c_name_3b
   integer :: sp0_3b
   integer :: sp1_3b
   integer :: sp2_3b
   integer :: max_np
   type(c_ptr) :: energies_3b_d
   type(c_ptr) :: forces_3b_d
   type(c_ptr) :: virial_3b_d
   type(c_ptr) :: kappas_array_d
   type(c_ptr) :: sigma_d
   type(c_ptr) :: neighbors_list_d
   integer(c_size_t) :: size_maxnp_bytes
   integer(c_size_t) :: size_maxnp_qs_bytes
   integer(c_size_t) :: size_alphas_bytes
   integer(c_size_t) :: size_energy3b
   integer(c_size_t) :: size_forces3b
   integer(c_size_t) :: size_virial3b

   !**************************************************************************

   !--- GPU variables for experimental ---!
   !  type( gpu_exp ) :: gpu_exp_vars
   type(c_ptr) :: pair_distribution_d
   type(c_ptr), allocatable :: rjs_index_d(:)

   real(dp), allocatable, target :: charges_temp(:)
   type(c_ptr) :: charges_d
   integer(c_size_t) :: st_charges_d

   ! Nested sampling
   real(dp) :: e_max
   real(dp) :: e_kin
   real(dp) :: rand
   real(dp) :: rand_scale(1:6)
   real(dp) :: mag
   real(dp) :: n_total_cutoff
   real(dp) :: n_total_cutoff_temp
   real(dp) :: dq
   real(dp) :: target_temp
   real(dp) :: f
   integer :: i_nested
   integer :: i_max
   integer :: i_image
   integer :: i_current_image = 1
   integer :: i_trial_image = 2
   type(image), allocatable :: images(:)
   type(image), allocatable :: images_temp(:)
   type(exp_data_container) :: temp_exp_container

   ! Storage of host arrays which are compatible with gpu implementation
   type(gpu_host_storage_type) :: gpu_host_temp

   ! this is a type of
   ! gpu_batch_storage( i_batch ) % host( i_n_dim_idx ) % pair_distribution_h( 1:n_samples )

   ! integer, parameter :: nstr=7
   ! type(c_ptr) :: virial_d(nstr),tmp_forces0_d(nstr),tmp_energies0_d(nstr)

   !**************************************************************************

   !--- TODO: Add in random seeds for repeatable calculations which rely on random numbers! ---!

   ! call random_number(AA)
   ! call random_number(BB)

   !**************************************************************************
   ! Start recording the time

   !  call cpu_time(time1)
   ! Start random seed
   !call srand(int(time1*1000))
   !**************************************************************************

   !**************************************************************************
   ! MPI stuff
#ifdef _MPIF90
   call mpi_init(ierr)
   call mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

   ! !time1=MPI_Wtime()
   !call get_time( !time1 )

   call get_time(time1)

   time3 = time1
   !  allocate( displs(1:ntasks) )
   !  allocate( displs2(1:ntasks) )
   !  allocate( counts(1:ntasks) )
   !  allocate( counts2(1:ntasks) )

#else
   call get_time(time1)

   time3 = time1

   rank = 0
   ntasks = 1
#endif
   allocate (n_atom_pairs_by_rank(1:ntasks))
   !**************************************************************************

   ! !---------------------------------------------------------------------!
   ! !------   DELETE THIS CODE THIS IS JUST FOR CHECKIGN CUDA-GDB   ------!
   ! !---------------------------------------------------------------------!

   ! if ( rank == 0 ) then
   !    i = 0
   !    print *, "rank ", rank, ": pid ", getpid(), " on %s ready for attach\n.";
   !    do while (0 == i)
   !       call sleep(5)
   !    end do
   ! end if

   !################################################!
   !###---   OPENMP PARALLELIZATION STARTUP   ---###!
   !################################################!

   n_omp = 1
   omp_task = 0

   !    !$omp parallel private(omp_task)
   !    !$ omp_task = omp_get_thread_num()
   !    !$omp end parallel

   !--- Creating GPU communication ---!
   call time_start(time%create_streams)
   call gpu_context_init(params, rank, n_omp)
   call gap_backend_init()
   call time_end(time%create_streams)
   ! print *, " "
   ! print *, " Time to create streams = ", time%create_streams(3), " Seconds"
   !call create_cublas_handle(cublas_handle)

   ! write(*,*) "Starting dummy kernel"

   ! dut1= MPI_Wtime()
   ! do n_ii=1,10

   !      do i_ii=1, 1024
   !      do j_jj=1,1024
   !      CC(i_ii,j_jj)=0.0
   !      do k_ii=1,1024
   !      CC(i_ii,j_jj)=CC(i_ii,j_jj)+AA(i_ii,k_ii)*BB(k_ii,j_jj)
   !      enddo
   !      enddo
   !      enddo
   ! enddo
   ! dut2= MPI_Wtime()
   ! write(*,*) "Ending dummy region"
   ! write(*,*) "Time spent in dummy region", dut2-dut1

   !**************************************************************************
   ! Read the mode. It should be "soap", "predict" or "md"
   !

   call get_command_argument(1, mode)
   if (mode == "" .or. mode == "none") then
      write (*, *) "ERROR: you need to run 'turbogap md' or 'turbogap predict'"
      stop
      ! THIS SHOULD BE FIXED, IN CASE THE USER JUST WANT TO OUTPUT THE SOAP DESCRIPTORS
      mode = "soap"
   end if
   !**************************************************************************

   !**************************************************************************
   ! Prints some welcome message and reads in the input file
   !
#ifdef _MPIF90
   IF (rank == 0) THEN
#endif
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
      write (*, *) '                         Miguel A. Caro                          |'
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
      write (*, *) 'Tigany Zarrouk                                                   |'
      write (*, *) '                                                                 |'
      write (*, *) '.................................................................|'
      write (*, *) '                                                                 |'
      write (*, *) '                     Last updated: June. 2023                     |'
      write (*, *) '                                        _________________________/'
      write (*, *) '.......................................|'
#ifdef _MPIF90
      write (*, *) '                                       |'
      write (*, *) 'Running TurboGAP with MPI+GPU support: |'
      write (*, *) '                                       |'
      write (*, '(A,I6,A)') ' Running TurboGAP on ', ntasks, ' MPI tasks   |'
      write (*, *) '                                       |'
      write (*, *) '.......................................|'
#else
      write (*, *) '                                       |'
      write (*, *) 'Running the serial version of TurboGAP |'
      write (*, *) '                                       |'
      write (*, *) '.......................................|'
#endif
#ifdef _MPIF90
   END IF
#endif
   !**************************************************************************

   !**************************************************************************
   ! Read input file and other files
   !
   call read_input_and_gap_files(mode, rank, ntasks, params, &
                                 soap_turbo_hypers, distance_2b_hypers, angle_3b_hypers, core_pot_hypers, &
                                 n_soap_turbo, n_distance_2b, n_angle_3b, n_core_pot, n_species, rcut_max, &
                                 valid_xps, xps_idx, vdw_lp_index, core_be_lp_index, &
                                 valid_estat_charges, charge_lp_index, &
                                 local_property_labels, local_property_indexes, n_local_properties_mpi, &
                                 has_local_properties_mpi, local_properties_n_sparse_mpi_soap_turbo, &
                                 local_properties_dim_mpi_soap_turbo, time)
   !**************************************************************************
   ! <----------------------------------------------------------------------------------------------- Finish printouts
#ifdef _MPIF90
   IF (rank == 0) THEN
#endif
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
      write (i_char, '(I8)') n_species
      write (*, '(1X,A,A8,A)') 'No. of species   = ', adjustl(i_char), '            |'
      do i = 1, n_species
         write (i_char, '(I8)') i
         write (*, '(1X,A,A2,A,A8,A)') '  *) Species#', adjustl(i_char), ' =       ', adjustr(params%species_types(i)), '      |'
      end do
      write (*, *) '---------------------------------      |'
      write (*, '(1X,A,F15.4,A)') 'rcut_max = ', rcut_max, ' Angst.      |'
      write (*, *) '---------------------------------      |'
      write (*, *) '                                       |'
      write (*, *) '.......................................|'
#ifdef _MPIF90
   END IF
#endif
   !**************************************************************************

   !**************************************************************************
   ! Print progress bar and initialize timers

   xps_idx = params%xps_idx
   md_istep = -1
   mc_istep = -1
   n_xyz = 0
   i_nested = 0
   i_image = 0

   if (params%do_md) then
#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         write (*, *) '                                       |'
         write (*, *) 'Doing molecular dynamics...            |'
         if (params%print_progress .and. md_istep > 0) then
            write (*, *) '                                       |'
            write (*, *) 'Progress:                              |'
            write (*, *) '                                       |'
            write (*, '(1X,A)', advance='no') '[                                    ] |'
         end if
#ifdef _MPIF90
      END IF
#endif
      update_bar = params%md_nsteps/36
      if (update_bar < 1) then
         update_bar = 1
      end if
      counter = 1
   end if
   !**************************************************************************

   !**************************************************************************
   !**************************************************************************
   !**************************************************************************
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
   perform%xps_forces = valid_xps .and. params%exp_forces

   do while (repeat_xyz .or. (params%do_md .and. md_istep < params%md_nsteps) &
             .or. (params%do_mc .and. mc_istep < params%mc_nsteps))
      exit_loop = .false.

      if (params%do_mc) then
         mc_istep = mc_istep + 1
         ! Undo if the step is md related
         if (md_istep > -1) mc_istep = mc_istep - 1

      end if

      if (params%do_md) then
         md_istep = md_istep + 1
      else
         n_xyz = n_xyz + 1
      end if

      !   Update progress bar
      if (params%print_progress .and. counter == update_bar .and. (.not. params%do_mc)) then
#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            do j = 1, 36 + 3
               write (*, "(A)", advance="no") creturn
            end do
            write (*, "(1X,A)", advance="no") "["
            do i = 1, 36*md_istep/params%md_nsteps
               write (*, "(A)", advance="no") "."
            end do
            do i = 36*md_istep/params%md_nsteps + 1, 36
               write (*, "(A)", advance="no") " "
            end do
            write (*, "(A)", advance="no") "] |"
            if (md_istep == params%md_nsteps) then
               write (*, *)
            end if

#ifdef _MPIF90
         END IF
#endif
         counter = 1
      else
         counter = counter + 1
      end if
      !**************************************************************************

      !**************************************************************************
      !   This chunk of code does all the reading/neighbor builds etc for each snapshot
      !   or MD step
      !   Read in XYZ file and build neighbors lists

      if ((params%do_md .and. md_istep == 0)) then

         !time%read_xyz(1) = MPI_wtime()
         call time_start(time%read_xyz)

#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            if (mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                             n_sites,.not. params%mc_write_xyz, fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz)
               rebuild_neighbors_list = .true.

            else if (.not. params%do_nested_sampling .or. mc_istep == 0) then
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                             n_sites, .false., fix_atom, params%t_beg, &
                             params%write_array_property(6), .false.)

            end if

            ! call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
            !               n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
            !               positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
            !               xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
            !               n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), .true. )
            !     Only rank 0 handles these variables
            !      allocate( positions_prev(1:3, 1:size(positions,2)) )
            !      allocate( positions_diff(1:3, 1:size(positions,2)) )
            if (.not. allocated(forces_prev)) allocate (forces_prev(1:3, 1:n_sites))
            if (.not. allocated(positions_prev)) allocate (positions_prev(1:3, 1:n_sites))
            if (.not. allocated(positions_diff)) allocate (positions_diff(1:3, 1:n_sites))
            positions_diff = 0.d0
            rebuild_neighbors_list = .true.
#ifdef _MPIF90
         END IF
#endif

         !time%read_xyz(2) = MPI_wtime()
         call time_end(time%read_xyz)
         !     If we're doing MD, we don't read beyond the first snapshot in the XYZ file
         repeat_xyz = .false.
         !     At the moment, we can't do prediction if the unit cell doesn't fit a whole cutoff sphere
#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            !     CLEAN THIS UP <------------------------------------------------------------------- LOOK HERE
            !      if( size(positions,2) /= n_sites )then
            if (.false.) then
               write (*, *) "Sorry, at the moment TurboGAP can't do MD for unit cells smaller than ", &
                  "a cutoff sphere <-- ERROR"
#ifdef _MPIF90
               call mpi_finalize(ierr)
#endif
               stop
            end if
#ifdef _MPIF90
         END IF
#endif
      else if (.not. params%do_md) then

         !time%read_xyz(1) = MPI_wtime()
         call time_start(time%read_xyz)

#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            if (mc_istep > 0) then
               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                             n_sites,.not. params%mc_write_xyz, fix_atom, params%t_beg, &
                             params%write_array_property(6),.not. params%mc_write_xyz)
               rebuild_neighbors_list = .true.
            else
               call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                             n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), &
                             .false.)
            end if
#ifdef _MPIF90
         END IF
#endif

         !time%read_xyz(2) = MPI_wtime()
         call time_end(time%read_xyz)
#ifdef _MPIF90

         !! time%mpi(1)=MPI_Wtime()
         !        call get_time( ! time%mpi(1) )

         call time_start(time%mpi)

         call mpi_bcast(repeat_xyz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         !! time%mpi(2)=MPI_Wtime()
         !        call get_time( ! time%mpi(2) )

         call time_end(time%mpi)
#endif
         rebuild_neighbors_list = .true.
      end if
      !   Broadcast the info in the XYZ file: positions, velocities, masses, xyz_species, xyz_species_supercell,
      !   species, species_supercell, indices, a_box, b_box, c_box and n_sites. I should put this into a module!!!!!!!

#ifdef _MPIF90
      IF (rank == 0) THEN
         n_pos = size(positions, 2)
         n_sp = size(xyz_species, 1)
         n_sp_sc = size(xyz_species_supercell, 1)
         if (params%do_mc .and. (mc_move /= "md" .or. md_istep == 0) .and. params%mc_hamiltonian) then
            if (mc_istep > 0) E_kinetic_prev = E_kinetic
            call random_number(velocities)
            call remove_cm_vel(velocities(1:3, 1:n_sites), masses(1:n_sites))
            E_kinetic = 0.d0
            do i = 1, n_sites
               E_kinetic = E_kinetic + 0.5d0*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
            end do
            instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic
            velocities = velocities*dsqrt(params%t_beg/instant_temp)
            if (mc_istep > 0) then
               E_kinetic = E_kinetic_prev
               instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic
               ! Reversing as we want it to be at the instant temp and not at t_beg
               velocities = velocities*dsqrt(instant_temp/params%t_beg)
            else
               E_kinetic = E_kinetic*params%t_beg/instant_temp
            end if
         end if

      END IF
       !! time%mpi(1)=MPI_Wtime()
      !call get_time( ! time%mpi(1) )

      call time_start(time%mpi)

      call mpi_bcast(n_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sp_sc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! time%mpi(2)=MPI_Wtime()
      call time_end(time%mpi)

      IF (rank /= 0) THEN
         if (allocated(positions)) deallocate (positions)
         allocate (positions(1:3, n_pos))
         if (params%do_md .or. params%do_nested_sampling .or. params%do_mc) then
            if (allocated(velocities)) deallocate (velocities)
            allocate (velocities(1:3, n_pos))
            !      allocate( masses(n_pos) )
            if (allocated(masses)) deallocate (masses)
            allocate (masses(1:n_sp))
            ! if(allocated( fix_atom ))deallocate( fix_atom )
            ! allocate( fix_atom(1:3, 1:n_sp) )
         end if
         if (allocated(xyz_species)) deallocate (xyz_species)
         allocate (xyz_species(1:n_sp))
         if (allocated(species)) deallocate (species)
         allocate (species(1:n_sp))
         if (allocated(xyz_species_supercell)) deallocate (xyz_species_supercell)
         allocate (xyz_species_supercell(1:n_sp_sc))
         if (allocated(species_supercell)) deallocate (species_supercell)
         allocate (species_supercell(1:n_sp_sc))
         if (allocated(fix_atom)) deallocate (fix_atom)
         allocate (fix_atom(1:3, 1:n_sp))

      END IF

      !time%mpi_positions(1) = MPI_wtime()
      call time_start(time%mpi_positions)

      call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (params%do_md .or. params%do_nested_sampling .or. params%do_mc .or. params%mc_hamiltonian) then
         call mpi_bcast(velocities, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(masses, n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(fix_atom, 3*n_sp, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
      call mpi_bcast(xyz_species, 8*n_sp, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(xyz_species_supercell, 8*n_sp_sc, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(species, n_sp, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(species_supercell, n_sp_sc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      !time%mpi_positions(2) = MPI_wtime()
      call time_end(time%mpi_positions)
#endif
      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(do_list)) deallocate (do_list)
         allocate (do_list(1:n_sites))
         do_list = .true.
      end if
      !

      call time_start(time%neigh)

      ! call get_time( time1 )

      ! !     !     call cpu_time(time1)
      ! call get_time(  )
      ! call get_time( time1 )

#ifdef _MPIF90
      !   Parallel neighbors list build
      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

      !   If we're using a box rescaling algorithm or a barostat, then the box size can
      !   become smaller or bigger than the cutoff sphere. If that happens, and the current
      !   situation is different from before, then we need to figure out if we need to
      !   construct a supercell (i.e., the box was bigger than the cutoff sphere and now
      !   is smaller -> makes computations slower) or default back to the primitive unit cell
      !   (i.e., the box was smaller and now is bigger -> makes computations faster).
      !   We only need to check if rebuild_neighbors_list = .true.
      if (rebuild_neighbors_list .and. params%do_mc .and. mc_istep > 0) then
         call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                       n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                       positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                       xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                       n_sites, .true., fix_atom, params%t_beg, &
                       params%write_array_property(6), .false.)
      else if (rebuild_neighbors_list) then
         call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                       n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                       positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                       xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                       n_sites, .true., fix_atom, params%t_beg, params%write_array_property(6), &
                       .false.)

      end if

#ifdef _MPIF90
      !   Overlapping domain decomposition with subcommunicators goes here <------------------- TO DO

      !   This is some trivial MPI parallelization to make sure the code works fine
      if (rank < mod(n_sites, ntasks)) then
         i_beg = 1 + rank*(n_sites/ntasks + 1)
      else
         i_beg = 1 + mod(n_sites, ntasks)*(n_sites/ntasks + 1) + (rank - mod(n_sites, ntasks))*(n_sites/ntasks)
      end if
      if (rank < mod(n_sites, ntasks)) then
         i_end = (rank + 1)*(n_sites/ntasks + 1)
      else
         i_end = i_beg + n_sites/ntasks - 1
      end if

      do_list = .false.
      do_list(i_beg:i_end) = .true.

      !      if( rebuild_neighbors_list )then
      !         if(allocated( rjs))deallocate( rjs)
      !         if(allocated( xyz))deallocate( xyz)
      !         if(allocated( thetas))deallocate( thetas)
      !         if(allocated( phis))deallocate( phis)
      !         if(allocated( neighbor_species ))deallocate( neighbor_species )
      !         if(allocated( neighbors_list))deallocate( neighbors_list )
      !         if(allocated(n_neigh))deallocate( n_neigh )
      !#ifdef _MPIF90
      !         if(allocated(n_neigh_local))deallocate( n_neigh_local )
      !#endif
      !      end if

      call build_neighbors_list(positions, a_box, b_box, c_box, params%do_timing, &
                                species_supercell, rcut_max, n_atom_pairs, rjs, &
                                thetas, phis, xyz, n_neigh_local, neighbors_list, neighbor_species, n_sites, indices, &
                                rebuild_neighbors_list, do_list, rank)

      if (rebuild_neighbors_list) then
         !     Get total number of atom pairs
         call mpi_allgather(n_atom_pairs, 1, MPI_INTEGER, n_atom_pairs_by_rank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
         n_atom_pairs_total = sum(n_atom_pairs_by_rank)
         n_atom_pairs = n_atom_pairs_total

         !     Get number of neighbors
         if (.not. allocated(n_neigh)) allocate (n_neigh(1:n_sites))
         call mpi_reduce(n_neigh_local, n_neigh, n_sites, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(n_neigh, n_sites, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         j_beg = 1
         j_end = n_atom_pairs_by_rank(rank + 1)
      end if
      !     print *, "-- Rank ", rank , " > set ", " i_beg = ", i_beg, " i_end = ", i_end,  " j_beg = ", j_beg, " j_end = ", j_end
#else

      call build_neighbors_list(positions, a_box, b_box, c_box, params%do_timing, &
                                species_supercell, rcut_max, n_atom_pairs, rjs, &
                                thetas, phis, xyz, n_neigh, neighbors_list, neighbor_species, n_sites, indices, &
                                rebuild_neighbors_list, do_list, rank)
      i_beg = 1
      i_end = n_sites
      n_sites_mpi = n_sites
      j_beg = 1
      j_end = n_atom_pairs
      n_atom_pairs_by_rank(rank + 1) = n_atom_pairs
#endif
!   Store by which rank each site is being handled
      if (allocated(site_in_rank)) then
         if (size(site_in_rank) /= n_sites) then
            deallocate (site_in_rank, this_site_in_rank)
         end if
      end if
      if (.not. allocated(site_in_rank)) then
         allocate (site_in_rank(1:n_sites))
         allocate (this_site_in_rank(1:n_sites))
      end if
      site_in_rank = 0
      this_site_in_rank = 0
      do i = i_beg, i_end
         this_site_in_rank(i) = rank
      end do
#ifdef _MPIF90
      call mpi_reduce(this_site_in_rank, site_in_rank, n_sites, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(site_in_rank, n_sites, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
      site_in_rank = this_site_in_rank
#endif

      !   Compute the volume of the "primitive" unit cell
      v_uc = dot_product(cross_product(a_box, b_box), c_box)/(dfloat(indices(1)*indices(2)*indices(3)))

      !    call cpu_time(time2)
      call get_time(time2)

       !! time2=MPI_Wtime()
      !call get_time( ! time2 )

      call time_end(time%neigh)
      !**************************************************************************

      !**************************************************************************
      !   If we are doing prediction, we run this chunk of code
      if (params%do_prediction .or. params%write_soap .or. params%write_derivatives) then

         !        print *, rank, " Allocating prediction arrays"
         !     We only need to reallocate the arrays if the number of sites changes
         ! REMOVE TRUE FROM IF STATEMENT
         if (n_sites /= n_sites_prev .or. params%do_mc) then
            if (allocated(energies)) deallocate (energies, energies_soap, energies_2b, energies_3b, energies_core_pot, &
                                              this_energies, energies_vdw, energies_estat, this_forces, energies_lp, energies_exp, &
                                                 energies_vdw_corr, mbd_ts_scaling, this_mbd_ts_scaling)
            allocate (energies(1:n_sites))
            allocate (this_energies(1:n_sites))
            allocate (energies_soap(1:n_sites))
            allocate (energies_2b(1:n_sites))
            allocate (energies_3b(1:n_sites))
            allocate (energies_core_pot(1:n_sites))
            allocate (energies_vdw(1:n_sites))
            allocate (energies_estat(1:n_sites))
            allocate (energies_lp(1:n_sites))
            allocate (energies_exp(1:n_sites))
            allocate (energies_vdw_corr(1:n_sites))
!          We do this allocations for van der Waals corrections
            allocate (mbd_ts_scaling(1:n_sites))
            allocate (this_mbd_ts_scaling(1:n_sites))

! Read in file for ts+mbd van der Waals mode if it exists
! Initialise the TS scaling factors. md_istep <= 0 rather than == 0 because
! predict and mc never advance md_istep past -1, and without this they reach
! get_ts_energy_and_forces with mbd_ts_scaling never having been set. For MD
! this is still exactly the first step, so the MD path is unchanged.
            if (params%vdw_type == "ts+mbd" .and. md_istep <= 0) then
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
                     do i = 1, n_sites
                        read (30, *) mbd_ts_scaling(i)
                     end do
                     update_mbd_ts_scaling = .false.
                  else
                     mbd_ts_scaling = 1.d0
                  end if
                  close (30)
                  this_mbd_ts_scaling = mbd_ts_scaling
               end if
#ifdef _MPIF90
               call mpi_bcast(this_mbd_ts_scaling, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
            end if

            if (perform%pdf) then
               if (allocated(energies_pdf)) deallocate (energies_pdf)
               allocate (energies_pdf(1:n_sites))
            end if

            if (perform%sf) then
               if (allocated(energies_sf)) deallocate (energies_sf)
               allocate (energies_sf(1:n_sites))
            end if

            if (perform%xrd) then
               if (allocated(energies_xrd)) deallocate (energies_xrd)
               allocate (energies_xrd(1:n_sites))
            end if

            if (perform%nd) then
               if (allocated(energies_nd)) deallocate (energies_nd)
               allocate (energies_nd(1:n_sites))
            end if

            !       This needs to be allocated even if no force prediction is needed:
            allocate (this_forces(1:3, 1:n_sites))
         end if
         energies = 0.d0
         energies_soap = 0.d0
         energies_2b = 0.d0
         energies_3b = 0.d0
         energies_core_pot = 0.d0
         energies_vdw = 0.d0
         energies_estat = 0.d0
         energies_lp = 0.d0
         energies_exp = 0.d0

         if (perform%pdf) energies_pdf = 0.d0
         if (perform%sf) energies_sf = 0.d0
         if (perform%xrd) energies_xrd = 0.d0
         if (perform%nd) energies_nd = 0.d0

         ! Adding allocation of local properties

         ! Now one could use pointers such that hirshfeld_v(:) acts as an alias for local_properties(vdw_index,:)...
         !        print *, rank, " Associating local properties "
         if (any(soap_turbo_hypers(:)%has_local_properties)) then
            if (n_sites /= n_sites_prev .or. params%do_mc) then
               if (allocated(local_properties)) then
                  nullify (this_local_properties_pt)
                  deallocate (this_local_properties, local_properties)
                  if (params%do_forces) then
                     nullify (this_local_properties_cart_der_pt)
                     deallocate (this_local_properties_cart_der, local_properties_cart_der)
                  end if
               end if
               allocate (local_properties(1:n_sites, 1:params%n_local_properties))
               allocate (this_local_properties(1:n_sites, 1:params%n_local_properties))
               this_local_properties_pt => this_local_properties

               !         I don't remember why this needs a pointer <----------------------------------------- CHECK

            end if
            local_properties = 0.d0

            if (params%do_forces) then
               if (n_atom_pairs_by_rank(rank + 1) /= n_atom_pairs_by_rank_prev) then
                  if (allocated(local_properties_cart_der)) deallocate (local_properties_cart_der, this_local_properties_cart_der)
                  allocate (local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
                  allocate (this_local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
               end if
               if (.not. allocated(local_properties_cart_der)) then
                  allocate (local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
                  allocate (this_local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank + 1), 1:params%n_local_properties))
               end if

               local_properties_cart_der = 0.d0
               this_local_properties_cart_der_pt =>&
                 & this_local_properties_cart_der(1:3,&
                 & 1:n_atom_pairs_by_rank(rank + 1), 1:params&
                 &%n_local_properties)
            end if
         end if

         ! Now go through the soap turbo hypers, and see if any are vdw or
         ! otherwise, if vdw, one can have pointers to point to the data
         ! structures such that it makes things clearer. One needs to check
         ! that this allocation still works iwth if(allocated(hirsh_v))
         ! statements

         if (params%do_forces) then
            if (n_sites /= n_sites_prev .or. params%do_mc) then
               if (allocated(forces)) deallocate (forces, forces_soap, forces_2b, forces_3b, &
                                                  forces_core_pot, forces_vdw, forces_estat, &
                                                  forces_lp)
               allocate (forces(1:3, 1:n_sites))
               allocate (forces_soap(1:3, 1:n_sites))
               allocate (forces_2b(1:3, 1:n_sites))
               allocate (forces_3b(1:3, 1:n_sites))
               allocate (forces_core_pot(1:3, 1:n_sites))
               allocate (forces_vdw(1:3, 1:n_sites))
               if (allocated(forces_vdw_corr)) deallocate (forces_vdw_corr)
               allocate (forces_vdw_corr(1:3, 1:n_sites))
               allocate (local_virial_vdw_diag_corr(1:3, 1:n_sites))
               allocate (local_virial_vdw_diag(1:3, 1:n_sites))
               allocate (forces_estat(1:3, 1:n_sites))
               allocate (forces_lp(1:3, 1:n_sites))

               if (perform%pdf_forces) then
                  if (allocated(forces_pdf)) deallocate (forces_pdf)
                  allocate (forces_pdf(1:3, 1:n_sites))
               end if

               if (perform%sf_forces) then
                  if (allocated(forces_sf)) deallocate (forces_sf)
                  allocate (forces_sf(1:3, 1:n_sites))
               end if

               if (perform%xrd_forces) then
                  if (allocated(forces_xrd)) deallocate (forces_xrd)
                  allocate (forces_xrd(1:3, 1:n_sites))
               end if

               if (perform%nd_forces) then
                  if (allocated(forces_nd)) deallocate (forces_nd)
                  allocate (forces_nd(1:3, 1:n_sites))
               end if

            end if
            forces = 0.d0
            forces_soap = 0.d0
            forces_2b = 0.d0
            forces_3b = 0.d0
            forces_core_pot = 0.d0
            forces_vdw = 0.d0
            forces_estat = 0.d0
            forces_lp = 0.d0
            virial = 0.d0
            virial_soap = 0.d0
            virial_2b = 0.d0
            virial_3b = 0.d0
            virial_core_pot = 0.d0
            virial_vdw = 0.d0
            local_virial_vdw_diag = 0.d0
            virial_estat = 0.d0
            virial_lp = 0.d0

            if (perform%pdf_forces) then
               forces_pdf = 0.d0
               virial_pdf = 0.d0
            end if

            if (perform%sf_forces) then
               forces_sf = 0.d0
               virial_sf = 0.d0
            end if

            if (perform%xrd_forces) then
               forces_xrd = 0.d0
               virial_xrd = 0.d0
            end if

            if (perform%nd_forces) then
               forces_nd = 0.d0
               virial_nd = 0.d0
            end if
         end if

         if (params%do_prediction) then
            !       Assign the e0 to each atom according to its species
            !        do i = 1, n_sites
            do i = i_beg, i_end
               do j = 1, n_species
                  if (xyz_species(i) == params%species_types(j)) then
                     energies(i) = params%e0(j)
                  end if
               end do
            end do
         end if
         !     Collect all energies
#ifdef _MPIF90

         !time%mpi_ef(1) = MPI_wtime()
         call time_start(time%mpi_ef)

         call mpi_reduce(energies, this_energies, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

         !time%mpi_ef(2) = MPI_wtime()
         call time_end(time%mpi_ef)
         energies = this_energies
#endif

         !     Loop through soap_turbo descriptors - we always call this routine, even if we don't want to do prediction
         n_lp_count = 0 ! This counts the local properties

         !!! COMMENTING OUT

         !#################################################################

         !              write(*,*) " > Starting get_gap_soap loop "
         call time_start(time%gap)
         do i = 1, n_soap_turbo

            call time_start(time%soap)

            !       Compute number of pairs for this SOAP. SOAP has in general a different cutoff than overall max
            !       cutoff, so the number of pairs may be a lot smaller for the SOAP subset.
            !       This subroutine splits the load optimally so as to not use more memory per MPI process than available.
            !       TurboGAP does not check how much memory is available, it just relies on heuristics and a user provided
            !       max_Gbytes_per_process (default = 1.d0)
            if (params%n_batches > 0) then
    call get_number_of_atom_pairs_batches(params%n_batches, n_neigh(i_beg:i_end), rjs(j_beg:j_end), soap_turbo_hypers(i)%rcut_max, &
                                                     soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%n_max, &
                                                     params%max_Gbytes_per_process, i_beg_list, i_end_list, j_beg_list, j_end_list)
            else
               call get_number_of_atom_pairs(n_neigh(i_beg:i_end), rjs(j_beg:j_end), soap_turbo_hypers(i)%rcut_max, &
                                             soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%n_max, &
                                             params%max_Gbytes_per_process, i_beg_list, i_end_list, j_beg_list, j_end_list)
            end if

            n_sp = soap_turbo_hypers(i)%n_species

            st_size_nf = n_sp*sizeof(soap_turbo_hypers(i)%nf(1))
            call gpu_malloc_async(nf_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%nf), nf_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(rcut_hard_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%rcut_hard), rcut_hard_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(rcut_soft_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%rcut_soft), rcut_soft_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(global_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%global_scaling), global_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_r_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%atom_sigma_r), atom_sigma_r_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%atom_sigma_r_scaling), atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_t_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%atom_sigma_t), atom_sigma_t_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%atom_sigma_t_scaling), atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(amplitude_scaling_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%amplitude_scaling), amplitude_scaling_d, st_size_nf, gpu_stream)
            call gpu_malloc_async(central_weight_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%central_weight), central_weight_d, st_size_nf, gpu_stream)
            st_size_nf = n_sp*sizeof(soap_turbo_hypers(i)%alpha_max(1))
            call gpu_malloc_async(alpha_max_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%alpha_max), alpha_max_d, st_size_nf, gpu_stream)
            n_sparse = soap_turbo_hypers(i)%n_sparse
            st_size_nf = n_sparse*sizeof(soap_turbo_hypers(i)%nf(1))
            call gpu_malloc_async(alphas_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%alphas), alphas_d, st_size_nf, gpu_stream)
            dim = soap_turbo_hypers(i)%dim
            st_size_nf = n_sparse*dim*sizeof(soap_turbo_hypers(i)%nf(1))
            call gpu_malloc_async(Qs_d, st_size_nf, gpu_stream)
            call cpy_htod(c_loc(soap_turbo_hypers(i)%Qs), Qs_d, st_size_nf, gpu_stream)

            if (soap_turbo_hypers(i)%has_local_properties) then
               ! Allocate gpu memory
               do j = 1, soap_turbo_hypers(i)%n_local_properties
                  soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas = &
                     soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                     sizeof(soap_turbo_hypers(i)%local_property_models(j)%alphas(1))
                  call gpu_malloc_async(soap_turbo_hypers(i)%local_property_models(j)%alphas_d, &
                                        soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas, gpu_stream)
                  call cpy_htod(c_loc(soap_turbo_hypers(i)&
                    &%local_property_models(j)%alphas), &
                    & soap_turbo_hypers(i)%local_property_models(j)&
                    &%alphas_d, soap_turbo_hypers(i)&
                    &%local_property_models(j)%st_size_alphas,&
                    & gpu_stream)

                  soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs = &
                     soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                     soap_turbo_hypers(i)%local_property_models(j)%dim* &
                     sizeof(soap_turbo_hypers(i)%local_property_models(j)%Qs(1, 1))

                  call gpu_malloc_async(soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                                        soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, gpu_stream)
                  call cpy_htod(c_loc(soap_turbo_hypers(i)%local_property_models(j)%Qs), &
                                soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                                soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, &
                                gpu_stream)

               end do

            end if

            do j = 1, size(i_beg_list)
               this_i_beg = i_beg - 1 + i_beg_list(j)
               this_i_end = i_beg - 1 + i_end_list(j)
               this_j_beg = j_beg - 1 + j_beg_list(j)
               this_j_end = j_beg - 1 + j_end_list(j)
               this_n_sites_mpi = this_i_end - this_i_beg + 1
               this_energies = 0.d0
               if (params%do_forces) then
                  this_forces = 0.d0
                  this_virial = 0.d0
               end if
               if (soap_turbo_hypers(i)%has_local_properties) then
                  this_local_properties = 0.d0
                  if (params%do_forces) then
                     this_local_properties_cart_der = 0.d0
                     !             I don't remember why this needs a pointer <----------------------------------------- CHECK
                     nullify (this_local_properties_cart_der_pt)
                     this_local_properties_cart_der_pt =>&
                       & this_local_properties_cart_der(1:3,&
                       & this_j_beg:this_j_end, 1:params&
                       &%n_local_properties)
                  end if
               end if

               call get_gap_soap( &
                  n_sparse, n_sites, this_n_sites_mpi, n_neigh(this_i_beg:this_i_end), &
                  neighbors_list(this_j_beg:this_j_end), soap_turbo_hypers(i)%n_species, &
                  soap_turbo_hypers(i)%species_types, rjs(this_j_beg:this_j_end), thetas(this_j_beg:this_j_end), &
                  phis(this_j_beg:this_j_end), xyz(1:3, this_j_beg:this_j_end), alpha_max_d, &
                  soap_turbo_hypers(i)%alpha_max, soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%dim, rcut_hard_d, &
                  soap_turbo_hypers(i)%rcut_hard, rcut_soft_d, nf_d, global_scaling_d, atom_sigma_r_d, &
                  soap_turbo_hypers(i)%atom_sigma_r, atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
                  amplitude_scaling_d, soap_turbo_hypers(i)%radial_enhancement, central_weight_d, &
                  soap_turbo_hypers(i)%central_weight, soap_turbo_hypers(i)%basis, &
                  soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
                  params%do_prediction, params%write_soap, params%write_derivatives, &
                  soap_turbo_hypers(i)%compress_soap, soap_turbo_hypers(i)%compress_soap_indices, &
                  soap_turbo_hypers(i)%delta, soap_turbo_hypers(i)%zeta, soap_turbo_hypers(i)%central_species, &
                  xyz_species(this_i_beg:this_i_end), xyz_species_supercell, alphas_d, Qs_d, params%all_atoms, &
                  params%which_atom, indices, soap, soap_cart_der, der_neighbors, der_neighbors_list, &
                  soap_turbo_hypers(i)%has_local_properties, soap_turbo_hypers(i)%n_local_properties, &
                  soap_turbo_hypers(i)%local_property_models, n_lp_count, this_energies, this_forces, &
                  this_local_properties_pt, this_local_properties_cart_der_pt, local_property_indexes, this_virial, &
                  time%soap_lin(3), time%get_soap(3), soap_turbo_hypers(i)%W_d, soap_turbo_hypers(i)%S_d, &
                  soap_turbo_hypers(i)%multiplicity_array_d, soap_turbo_hypers(i)%st_W_d, &
                  soap_turbo_hypers(i)%st_S_d, soap_turbo_hypers(i)%st_multiplicity_array_d, &
                  soap_turbo_hypers(i)%recompute_basis, time%local_prop, cublas_handle, gpu_stream)

               energies_soap = energies_soap + this_energies

               if (soap_turbo_hypers(i)%has_local_properties) then

                  local_properties(:, :) = local_properties(:, :) + this_local_properties(:, :)
                  if (any(soap_turbo_hypers(i)&
                    &%local_property_models(:)%do_derivatives) &
                    & .and. params%do_derivatives) then
                     local_properties_cart_der(:, :, :) =&
                       & local_properties_cart_der(:, :, :) +&
                       & this_local_properties_cart_der(:, :, :)
                  end if

                  if (soap_turbo_hypers(i)%has_vdw) then
                     !                    hirshfeld_v => local_properties( :, vdw_lp_index)
                     ! if (any(soap_turbo_hypers(i)&
                     !      &%local_property_models(:)%do_derivatives)&
                     !      & .and. params%do_derivatives)&
                     !      & hirshfeld_v_cart_der =>&
                     !      & local_properties_cart_der( :, :,&
                     !      & vdw_lp_index)
                  end if

               end if
               if (params%do_forces) then
                  forces_soap = forces_soap + this_forces
                  virial_soap = virial_soap + this_virial
               end if
            end do

            n_lp_count = n_lp_count + soap_turbo_hypers(i)%n_local_properties

            !           print *, rank, " >> Freeing gpu memory  "
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
            call gpu_free_async(alphas_d, gpu_stream)

            if (soap_turbo_hypers(i)%has_local_properties) then
               do j = 1, soap_turbo_hypers(i)%n_local_properties
                  call gpu_free_async(soap_turbo_hypers(i)%local_property_models(j)%alphas_d, gpu_stream)
                  call gpu_free_async(soap_turbo_hypers(i)%local_property_models(j)%Qs_d, gpu_stream)
               end do
            end if

            call gpu_free_async(Qs_d, gpu_stream)

           !!time%soap_solo(2 = MPI_wtime()
            !        call get_time( time%soap_solo(2  )

            ! ! time%soap_solo(2)=MPI_Wtime()
            call get_time(time%soap_solo(2))

            deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)

            time%soap_solo(3) = time%soap_solo(3) + time%soap_solo(2) - time%soap_solo(1)

            ! THIS WON'T WORK! THE SOAP AND SOAP DERIVATIVES NEED TO BE COLLECTED FROM ALL RANKS <--------------------- FIX THIS!!!!
            ! AT THE MOMENT I'M MAKING THE CODE PRINT AN ERROR MESSAGE AND STOP EXECUTION IF THE USER TRIES TO WRITE OUT THESE
            ! FILES WITH MORE THAN ONE MPI TASK
#ifdef _MPIF90
            IF (rank == 0) THEN
#endif
               !       Write out stuff - THIS SHOULD PROBABLY BE PUT IN A MODULE
               if (n_soap_turbo == 1) then
                  i_char = ""
               else
                  write (i_char, '(I7)') i
                  i_char = "_"//adjustl(i_char)
               end if
               !       Write the SOAP vectors - NOT THE OPTIMAL STRATEGY IN TERMS OF DISK SPACE SINCE SOME ATOMS HAVE SOAP = 0
               if (params%write_soap) then
                  if (n_xyz == 1 .or. md_istep == 0) then
                     open (unit=10, file="soap"//trim(i_char)//".dat", status="unknown")
                  else
                     open (unit=10, file="soap"//trim(i_char)//".dat", status="old", position="append")
                  end if
                  if (.not. params%do_md .or. &
                      (params%do_md .and. (md_istep == 0 .or. md_istep == params%md_nsteps .or. &
                                           modulo(md_istep, params%write_xyz) == 0))) then
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
                  if (n_xyz == 1 .or. md_istep == 0) then
                     open (unit=10, file="soap_der"//trim(i_char)//".dat", status="unknown")
                  else
                     open (unit=10, file="soap_der"//trim(i_char)//".dat", status="old", position="append")
                  end if
                  if (.not. params%do_md .or. &
                      (params%do_md .and. (md_istep == 0 .or. md_istep == params%md_nsteps .or. &
                                           modulo(md_istep, params%write_xyz) == 0))) then
                     !           Note, this n_sites is not the same as the total number of sites, it's just the total number
                     !           of sites that have a derivative, since the first neighbor of each site is itself, the site
                     !           ID can always be retrieved from there. Note also that the sites are not necessarily given in
                     !           order
                     n_sites_this = size(der_neighbors, 1)
                     n_soap = size(soap_cart_der, 2)
                     n_atom_pairs = size(der_neighbors_list, 1)
                     write (10, *) n_sites, n_soap, n_atom_pairs
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
#ifdef _MPIF90
            END IF
#endif

            !time%soap(2) = MPI_wtime()
            call time_end(time%soap)

         end do
         call time_end(time%gap)

         !#################################################################

#ifdef _MPIF90
         if (any(soap_turbo_hypers(:)%has_local_properties)) then
            ! time%mpi(1)=MPI_Wtime()
            call time_start(time%mpi)

            call mpi_reduce(local_properties, this_local_properties, n_sites*params%n_local_properties,&
              & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
              & ierr)
            !           if( any( soap_turbo_hypers(:)%has_vdw ) )then
            ! call mpi_reduce(hirshfeld_v, this_hirshfeld_v, n_sites,&
            !      & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
            !      & ierr)
            !        if( params%do_forces )then
            !         I'm not sure if this is necessary at all... CHECK
            !          call mpi_reduce(hirshfeld_v_cart_der,
            !          this_hirshfeld_v_cart_der, 3*n_atom_pairs,
            !          MPI_DOUBLE_PRECISION, MPI_SUM, 0,
            !          MPI_COMM_WORLD, ierr)
            !          hirshfeld_v_cart_der = this_hirshfeld_v_cart_der
            !        end if
            !            hirshfeld_v = this_hirshfeld_v
            local_properties = this_local_properties
            !           call mpi_bcast(hirshfeld_v, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(local_properties, n_sites*params%n_local_properties, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

            ! time%mpi(2)=MPI_Wtime()
            call time_end(time%mpi)
         end if
#endif

         !     Compute vdW energies and forces
!   vdW / MBD, via turbogap_vdw.f90, adopted from the CPU branch.  What was
!   inline here handled vdw_type == "ts" only; compute_vdw handles ts, mbd and
!   ts+mbd, and carries the ts+mbd correction state in vdw_ws rather than in
!   driver locals whose lifetime was invisible.  vdw.f90 and misc.f90 are
!   already byte-identical between the trees, so this is a straight transplant.
         call compute_vdw(params, any(soap_turbo_hypers(:)%has_vdw), n_sites, &
                          n_neigh, neighbors_list, neighbor_species, rjs, xyz, &
                          local_properties, local_properties_cart_der, vdw_lp_index, &
                          i_beg, i_end, j_beg, j_end, n_atom_pairs_by_rank, site_in_rank, &
                          indices, rank, ntasks, md_istep, vdw_ws, &
                          energies_vdw, forces_vdw, virial_vdw, local_virial_vdw_diag, &
                          this_energies_vdw, this_forces_vdw, this_virial_vdw, &
                          this_local_virial_vdw_diag, energies_vdw_corr, forces_vdw_corr, &
                          local_virial_vdw_diag_corr, mbd_ts_scaling, this_mbd_ts_scaling, &
                          update_mbd_ts_scaling, time)

         !     Compute ELECTROSTATIC energies and forces
!        valid_estat_charges is part of the guard, not an afterthought: without it a
!        deck that asks for electrostatics against a GAP with no atomic_charge local
!        property indexes local_properties with an uninitialised charge_lp_index and
!        segfaults. Reproduced on the CPU branch, where the same block was ported.
!        Moved to src/turbogap_estat.f90. The #ifdef is here, at the one call,
!        rather than inside four continued argument lists where nothing
!        Fortran-aware could parse it.
#ifdef _MPIF90
         call compute_estat(params, do_electrostatics, valid_estat_charges, charge_lp_index, &
                            n_sites, n_neigh, neighbors_list, species, neighbor_species, rjs, xyz, &
                            local_properties, local_properties_cart_der, &
                            i_beg, i_end, j_beg, j_end, rank, n_omp, &
                            this_energies_estat, this_forces_estat, this_virial_estat, time)
#else
         call compute_estat(params, do_electrostatics, valid_estat_charges, charge_lp_index, &
                            n_sites, n_neigh, neighbors_list, species, neighbor_species, rjs, xyz, &
                            local_properties, local_properties_cart_der, &
                            i_beg, i_end, j_beg, j_end, rank, n_omp, &
                            energies_estat, forces_estat, virial_estat, time)
#endif

         !----------------------------------------------------!
         !--- EXPERIMENTAL SPECTRUM CALCULATION AND FORCES ---!
         !----------------------------------------------------!

         ! --- Changing the implementation:
         !     > All experimental prediction should be done here
         !     > do_exp is the variable which says whether calculation should be done
         !     > experimental_forces = .true. will add forces to the calculation

         !#########################################################!
         !###---   Compute Experimental Data Interpolation   ---###!
         !#########################################################!

         if (params%do_exp) then
            do i = 1, params%n_exp
               ! If we want to compute the experimental interpolation, we do it now.

               call get_write_condition(params%do_mc, params%do_md&
                 &, mc_istep, md_istep, params%write_xyz,&
                 & write_condition)

               if (params%exp_data(i)%compute_exp) then
                  if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
                  if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
                  ! print *, " exp n samples ", params%exp_data(i)%n_samples
                  ! print *, " sf n samples  ", params%structure_factor_n_samples
                  call calculate_exp_interpolation(params%exp_data(i)&
                    &%x, params%exp_data(i)%y, params%exp_data(i)&
                    &%n_samples, params%exp_data(i)%data)

                  call preprocess_exp_data(params, params%exp_data(i)%x,&
                    & params%exp_data(i)%y, params%exp_data(i)%label,&
                    & n_sites, dot_product(cross_product(a_box,&
                    & b_box), c_box)/(dfloat(indices(1)*indices(2) &
                    &*indices(3))), params%exp_data(i)%input, exp_output, .true.)

                  if (params%write_exp .and. .not. params&
                    &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition) then

                     call get_overwrite_condition(params%do_mc,&
                       & params%do_md, mc_istep, md_istep, params&
                       &%write_xyz, overwrite_condition)

                     call write_exp_data(params%exp_data(i)%x, params&
                       &%exp_data(i)%y, overwrite_condition,&
                       & trim(params%exp_data(i)%label)//&
                       & "_exp.dat", params%exp_data(i)%label)
                  end if

               end if

               !          if ( params%exp_data(i)%compute_exp .and. .not.  params&
               !            &%exp_data(i)%wrote_exp .and. rank == 0  .and. write_condition) then

               !            if (params%write_exp) then
               !              write(filename,'(A)')&
               !                & trim(params%exp_data(i)%label) // "_exp_fit.dat"

               !              call get_overwrite_condition( params%do_mc,&
               !                & params%do_md, mc_istep, md_istep, params&
               !                &%write_xyz, overwrite_condition)

               !              call write_exp_data(params%exp_data(i)%x, params%exp_data(i)%y,&
               !                & overwrite_condition, trim(params&
               !                &%exp_data(i)%label) // "_exp_fit.dat",&
               !                & trim(params%exp_data(i)%label) // " : output = "&
               !                & // trim( exp_output ))

               !            end if

               !          end if

               params%exp_data(i)%wrote_exp = .true.
               params%exp_data(i)%compute_exp = .true.

            end do
         end if

         !###################################################!
         !###---   XPS Forces and Spectra Prediction   ---###!
         !###################################################!

         !     Compute core_electron_be energies and forces
         !
         ! Moved to src/turbogap_exp.f90. The #ifdef is here, at the one call,
         ! rather than inside a continued argument list where nothing
         ! Fortran-aware could parse it.
#ifdef _MPIF90
         call compute_exp_xps(params, n_sites, n_xyz, xyz, neighbors_list, n_neigh, &
                              local_properties, local_properties_cart_der, soap_turbo_hypers, &
                              a_box, b_box, c_box, indices, i_beg, i_end, j_beg, j_end, rank, &
                              md_istep, mc_istep, valid_xps, xps_idx, core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              this_energies_lp, this_forces_lp, this_virial_lp, time)
#else
         call compute_exp_xps(params, n_sites, n_xyz, xyz, neighbors_list, n_neigh, &
                              local_properties, local_properties_cart_der, soap_turbo_hypers, &
                              a_box, b_box, c_box, indices, i_beg, i_end, j_beg, j_end, rank, &
                              md_istep, mc_istep, valid_xps, xps_idx, core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              energies_lp, forces_lp, virial_lp, time)
#endif

         !##############################################################!
         !###---   (Partial) Pair distribution functions and XRD   ---###!
         !##############################################################!
         !
         ! Moved to src/turbogap_exp.f90. The #ifdef is here, at the one
         ! call, rather than inside four continued argument lists where
         ! nothing Fortran-aware could parse it.
#ifdef _MPIF90
         call compute_exp_spectra(params, n_sites, species, rjs, xyz, neighbors_list, &
                                  n_neigh, neighbor_species, indices, a_box, b_box, c_box, i_beg, i_end, j_beg, &
                                  j_end, rank, ntasks, ierr, md_istep, mc_istep, this_energies_sf, &
                                  this_forces_sf, this_virial_sf, this_energies_xrd, this_forces_xrd, &
                                  this_virial_xrd, this_energies_nd, this_forces_nd, this_virial_nd, time, &
                                  i_beg_list, i_end_list, j_beg_list, &
                                  j_end_list, n_omp, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                                  n_sites_temp, n_pairs_temp, write_condition, overwrite_condition, &
                                  temp_string, species_types_actual, v_uc)
#else
         call compute_exp_spectra(params, n_sites, species, rjs, xyz, neighbors_list, &
                                  n_neigh, neighbor_species, indices, a_box, b_box, c_box, i_beg, i_end, j_beg, &
                                  j_end, rank, ntasks, ierr, md_istep, mc_istep, energies_sf, forces_sf, &
                                  virial_sf, energies_xrd, forces_xrd, virial_xrd, energies_nd, forces_nd, &
                                  virial_nd, time, i_beg_list, &
                                  i_end_list, j_beg_list, j_end_list, n_omp, omp_task, this_i_beg, this_i_end, &
                                  this_j_beg, this_j_end, n_sites_temp, n_pairs_temp, write_condition, &
                                  overwrite_condition, temp_string, species_types_actual, v_uc)
#endif

         if (params%do_prediction) then
            call time_start(time%gap)

            if (n_core_pot > 0 .or. n_distance_2b > 0 .or. n_angle_3b > 0) then

               ! Two-body, core-potential and three-body contributions, via the
               ! gap_backend seam. The device buffers these used to reach by host
               ! association now live inside src/gap_backend_gpu.f90, so nothing
               ! here names a device pointer and the calls match the CPU branch's
               ! gap_backend_cpu.f90 exactly.
               call gap_backend_begin(params, rjs, xyz, n_neigh, species, neighbor_species, &
                                      neighbors_list, i_beg, i_end, j_beg, j_end)

               call add_2b_contribution(n_distance_2b, distance_2b_hypers, &
                                        params, rjs, xyz, n_neigh, species, neighbor_species, &
                                        i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                        energies_2b, forces_2b, virial_2b, time)

               call add_core_pot_contribution(n_core_pot, core_pot_hypers, &
                                              params, rjs, xyz, n_neigh, species, neighbor_species, &
                                              i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                              energies_core_pot, forces_core_pot, virial_core_pot, time)

               call add_3b_contribution(n_angle_3b, angle_3b_hypers, neighbors_list, &
                                        params, rjs, xyz, n_neigh, species, neighbor_species, &
                                        i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
                                        forces, energies_3b, forces_3b, virial_3b, time)

               call gap_backend_end()

            end if

            call time_end(time%gap)
            !       Communicate all energies and forces here for all
            !       terms
#ifdef _MPIF90

            !time%mpi_ef(1) = MPI_wtime()
            call time_start(time%mpi_ef)

!       One evaluation of the eleven predicates, and one list built from them.
!       The pack and unpack walks below read only that list, so they cannot
!       disagree about which slot belongs to which family.
            contrib_on(C_SOAP) = (n_soap_turbo > 0)
            contrib_on(C_VDW) = allocated(this_energies_vdw)
            contrib_on(C_ESTAT) = allocated(this_energies_estat)
            contrib_on(C_LP) = allocated(this_energies_lp)
            contrib_on(C_PDF) = allocated(this_energies_pdf) .and. params%valid_pdf
            contrib_on(C_SF) = allocated(this_energies_sf) .and. params%valid_sf
            contrib_on(C_XRD) = allocated(this_energies_xrd) .and. params%valid_xrd
            contrib_on(C_ND) = allocated(this_energies_nd) .and. params%valid_nd
            contrib_on(C_2B) = (n_distance_2b > 0)
            contrib_on(C_CP) = (n_core_pot > 0)
            contrib_on(C_3B) = (n_angle_3b > 0)

            n_active = 0
            if (contrib_on(C_SOAP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => energies_soap
               contrib(n_active)%e_dst => energies_soap
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => forces_soap
                  contrib(n_active)%v_src => virial_soap
                  contrib(n_active)%f_dst => forces_soap
                  contrib(n_active)%v_dst => virial_soap
               end if
            end if
            if (contrib_on(C_VDW)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_vdw
               contrib(n_active)%e_dst => energies_vdw
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_vdw
                  contrib(n_active)%v_src => this_virial_vdw
                  contrib(n_active)%f_dst => forces_vdw
                  contrib(n_active)%v_dst => virial_vdw
               end if
            end if
            if (contrib_on(C_ESTAT)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_estat
               contrib(n_active)%e_dst => energies_estat
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_estat
                  contrib(n_active)%v_src => this_virial_estat
                  contrib(n_active)%f_dst => forces_estat
                  contrib(n_active)%v_dst => virial_estat
               end if
            end if
            if (contrib_on(C_LP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_lp
               contrib(n_active)%e_dst => energies_lp
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_lp
                  contrib(n_active)%v_src => this_virial_lp
                  contrib(n_active)%f_dst => forces_lp
                  contrib(n_active)%v_dst => virial_lp
               end if
            end if
            if (contrib_on(C_PDF)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_pdf
               contrib(n_active)%e_dst => energies_pdf
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_pdf
                  contrib(n_active)%v_src => this_virial_pdf
                  contrib(n_active)%f_dst => forces_pdf
                  contrib(n_active)%v_dst => virial_pdf
               end if
            end if
            if (contrib_on(C_SF)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_sf
               contrib(n_active)%e_dst => energies_sf
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_sf
                  contrib(n_active)%v_src => this_virial_sf
                  contrib(n_active)%f_dst => forces_sf
                  contrib(n_active)%v_dst => virial_sf
               end if
            end if
            if (contrib_on(C_XRD)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_xrd
               contrib(n_active)%e_dst => energies_xrd
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_xrd
                  contrib(n_active)%v_src => this_virial_xrd
                  contrib(n_active)%f_dst => forces_xrd
                  contrib(n_active)%v_dst => virial_xrd
               end if
            end if
            if (contrib_on(C_ND)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => this_energies_nd
               contrib(n_active)%e_dst => energies_nd
               contrib(n_active)%forces = params%do_forces .and. params%exp_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => this_forces_nd
                  contrib(n_active)%v_src => this_virial_nd
                  contrib(n_active)%f_dst => forces_nd
                  contrib(n_active)%v_dst => virial_nd
               end if
            end if
            if (contrib_on(C_2B)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => energies_2b
               contrib(n_active)%e_dst => energies_2b
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => forces_2b
                  contrib(n_active)%v_src => virial_2b
                  contrib(n_active)%f_dst => forces_2b
                  contrib(n_active)%v_dst => virial_2b
               end if
            end if
            if (contrib_on(C_CP)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => energies_core_pot
               contrib(n_active)%e_dst => energies_core_pot
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => forces_core_pot
                  contrib(n_active)%v_src => virial_core_pot
                  contrib(n_active)%f_dst => forces_core_pot
                  contrib(n_active)%v_dst => virial_core_pot
               end if
            end if
            if (contrib_on(C_3B)) then
               n_active = n_active + 1
               contrib(n_active)%e_src => energies_3b
               contrib(n_active)%e_dst => energies_3b
               contrib(n_active)%forces = params%do_forces
               if (contrib(n_active)%forces) then
                  contrib(n_active)%f_src => forces_3b
                  contrib(n_active)%v_src => virial_3b
                  contrib(n_active)%f_dst => forces_3b
                  contrib(n_active)%v_dst => virial_3b
               end if
            end if

            counter2 = n_active

            allocate (all_energies(1:n_sites, 1:counter2))
            allocate (all_this_energies(1:n_sites, 1:counter2))
            if (params%do_forces) then
               allocate (all_forces(1:3, 1:n_sites, 1:counter2))
               allocate (all_this_forces(1:3, 1:n_sites, 1:counter2))
               allocate (all_virial(1:3, 1:3, 1:counter2))
               allocate (all_this_virial(1:3, 1:3, 1:counter2))
            end if

!       Pack. A family owns a slot whenever it is active, but only contributes
!       forces when it carries them. Its slot must still be cleared: all_forces
!       is allocated and never zeroed, and mpi_reduce reads the whole array.
            do i_contrib = 1, n_active
               all_energies(1:n_sites, i_contrib) = contrib(i_contrib)%e_src(1:n_sites)
               if (contrib(i_contrib)%forces) then
                  all_forces(1:3, 1:n_sites, i_contrib) = contrib(i_contrib)%f_src(1:3, 1:n_sites)
                  all_virial(1:3, 1:3, i_contrib) = contrib(i_contrib)%v_src(1:3, 1:3)
               else if (params%do_forces) then
                  all_forces(1:3, 1:n_sites, i_contrib) = 0.d0
                  all_virial(1:3, 1:3, i_contrib) = 0.d0
               end if
            end do

            !       Here we communicate
            call mpi_reduce(all_energies, all_this_energies, n_sites&
                 &*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                 & MPI_COMM_WORLD, ierr)
            if (params%do_forces) then
               call mpi_reduce(all_forces, all_this_forces, 3*n_sites&
                    &*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                    & MPI_COMM_WORLD, ierr)
               call mpi_reduce(all_virial, all_this_virial, 9*counter2&
                    &, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                    & MPI_COMM_WORLD, ierr)
            end if

!       Unpack. For the seven families packed from a this_ array this is where
!       the reduced result lands in the un-prefixed one.
            do i_contrib = 1, n_active
               contrib(i_contrib)%e_dst(1:n_sites) = all_this_energies(1:n_sites, i_contrib)
               if (contrib(i_contrib)%forces) then
                  contrib(i_contrib)%f_dst(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, i_contrib)
                  contrib(i_contrib)%v_dst(1:3, 1:3) = all_this_virial(1:3, 1:3, i_contrib)
               end if
            end do

!       Release the this_ arrays now that their contents have been unpacked.
!       Kept explicit: an allocatable cannot be deallocated through a pointer.
!       this_local_virial_vdw_diag has no counterpart in the list, so it is
!       released here by name; compute_vdw allocates it per call.
            if (contrib_on(C_VDW)) then
               deallocate (this_energies_vdw)
               if (params%do_forces) deallocate (this_forces_vdw, this_local_virial_vdw_diag)
            end if
            if (contrib_on(C_ESTAT)) then
               deallocate (this_energies_estat)
               if (params%do_forces) deallocate (this_forces_estat)
            end if
            if (contrib_on(C_LP)) then
               deallocate (this_energies_lp)
               if (params%do_forces) deallocate (this_forces_lp)
            end if
            if (contrib_on(C_PDF)) then
               deallocate (this_energies_pdf)
               if (params%do_forces .and. params%exp_forces) deallocate (this_forces_pdf)
            end if
            if (contrib_on(C_SF)) then
               deallocate (this_energies_sf)
               if (params%do_forces .and. params%exp_forces) deallocate (this_forces_sf)
            end if
            if (contrib_on(C_XRD)) then
               deallocate (this_energies_xrd)
               if (params%do_forces .and. params%exp_forces) deallocate (this_forces_xrd)
            end if
            if (contrib_on(C_ND)) then
               deallocate (this_energies_nd)
               if (params%do_forces .and. params%exp_forces) deallocate (this_forces_nd)
            end if

            !       Clean up
            deallocate (all_energies, all_this_energies)
            if (params%do_forces) then
               deallocate (all_forces, all_this_forces, all_virial, all_this_virial)
            end if

            !time%mpi_ef(2) = MPI_wtime()
            call time_end(time%mpi_ef)
#endif

            !       Add up all the energy terms
            energies = energies + energies_soap + energies_2b +&
              & energies_3b + energies_core_pot + energies_vdw + energies_estat

            if (valid_xps) energies_exp = energies_exp + energies_lp
            if (perform%pdf) energies_exp = energies_exp + energies_pdf
            if (perform%sf) energies_exp = energies_exp + energies_sf
            if (perform%xrd) energies_exp = energies_exp + energies_xrd
            if (perform%nd) energies_exp = energies_exp + energies_nd

            if (params%exp_energies) energies = energies + energies_exp

            energy_prev = energy
            instant_pressure_prev = instant_pressure
            energy = sum(energies)
            energy_exp = sum(energies_exp)

         end if

         if (.not. params%do_md .and. .not. params%do_mc) then
#ifdef _MPIF90
            IF (rank == 0) then
#endif
               write (*, *) '                                       |'
               write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(energies_soap), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(energies_2b), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(energies_3b), 'eV |'
               write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(energies_core_pot), 'eV |'
               write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(energies_vdw), 'eV |'
               write (*, '(A,1X,F21.8,1X,A)') ' estat energy:', sum(energies_estat), 'eV |'
               write (*, '(A,1X,F22.8,1X,A)') ' Exp. energy:', sum(energies_exp), 'eV |'
               if (valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(energies_lp), 'eV |'
               if (perform%pdf)&
                 & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                 & sum(energies_pdf), 'eV |'
               if (perform%sf)&
                 & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                 & sum(energies_sf), 'eV |'
               if (perform%xrd)&
                 & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                 & sum(energies_xrd), 'eV |'
               if (perform%nd)&
                 & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                 & sum(energies_nd), 'eV |'

               if (.not. params%do_mc .or. (params%do_mc .and. mc_istep <= 1)) then
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(energies), 'eV |'
               else
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(images(i_trial_image)%energies), 'eV |'
               end if

               if (.not. params%do_mc) then
                  write (*, *) '                                       |'
                  write (*, *) 'Energy & forces in "trajectory_out.xyz"|'
                  write (*, *) '                                       |'
                  write (*, *) '.......................................|'
               else if (mc_istep == 0) then
                  write (*, *) '                                       |'
                  write (*, *) ' MC configs in "mc_current.xyz" and    |'
                  write (*, *) '               "mc_trial.xyz"          |'
                  write (*, *) '               "mc_all.xyz"            |'
                  write (*, *) '.......................................|'
               end if
#ifdef _MPIF90
            END IF
#endif
         end if

         if (params%do_forces) then
            forces = forces_soap + forces_2b + forces_3b + forces_core_pot + forces_vdw
            virial = virial_soap + virial_2b + virial_3b + virial_core_pot + virial_vdw

            if (valid_estat_charges) forces = forces + forces_estat
            if (valid_estat_charges) virial = virial + virial_estat

            if (perform%xps_forces) forces = forces + forces_lp
            if (perform%xps_forces) virial = virial + virial_lp

            if (perform%pdf_forces) forces = forces + forces_pdf
            if (perform%pdf_forces) virial = virial + virial_pdf

            if (perform%sf_forces) forces = forces + forces_sf
            if (perform%sf_forces) virial = virial + virial_sf

            if (perform%xrd_forces) forces = forces + forces_xrd
            if (perform%xrd_forces) virial = virial + virial_xrd

            if (perform%nd_forces) forces = forces + forces_nd
            if (perform%nd_forces) virial = virial + virial_nd

            if (rank == 0 .and. params%print_vdw_forces) then
               print *, "> Virial ESTAT "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", virial_estat(i, j)
                  end do
               end do

               print *, "> Virial soap "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", virial_soap(i, j)
                  end do
               end do

               print *, "> Virial 2b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", virial_2b(i, j)
                  end do
               end do

               print *, "> Virial 3b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", virial_3b(i, j)
                  end do
               end do

               print *, "> Virial core_pot "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", virial_core_pot(i, j)
                  end do
               end do

               if (perform%xrd_forces) then
                  print *, "> Virial xrd "
                  do i = 1, 3
                     do j = 1, 3
                        print *, " i, ", i, " j ", j, " ", virial_xrd(i, j)
                     end do
                  end do
                  temp_string = ""
                  temp_string2 = ""
                  write (temp_string, "(I8)") md_istep
                  write (temp_string2, "(A)") "forces_xrd_"//trim(adjustl(temp_string))
                  open (unit=90, file=temp_string2, status="unknown")
                  do i = 1, n_sites
                     write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                        forces_xrd(1, i), forces_xrd(2, i), forces_xrd(3, i)
                  end do
                  close (90)

               end if

            end if

            if (params%print_vdw_forces) then
               open (unit=90, file="forces_vdw", status="unknown")
               do i = 1, n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     forces_vdw(1, i), forces_vdw(2, i), forces_vdw(3, i)
               end do
               close (90)

            end if

            if (rank == 0 .and. params%print_estat_forces) then
               open (unit=90, file="forces_estat", status="unknown")
               do i = 1, n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     forces_estat(1, i), forces_estat(2, i), forces_estat(3, i)
               end do
               close (90)

               open (unit=90, file="charge_gradients_estat", status="unknown")
               do i = 1, n_atom_pairs_by_rank(rank + 1)
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     local_properties_cart_der(1, i, charge_lp_index), &
                     local_properties_cart_der(2, i, charge_lp_index), &
                     local_properties_cart_der(3, i, charge_lp_index)
               end do
               close (90)

            end if

         end if
         ! For debugging the virial implementation
         if (rank == 0 .and. .false.) then
            write (*, *) "pressure_soap: ", virial_soap/3.d0/v_uc
            write (*, *) "pressure_vdw: ", virial_vdw/3.d0/v_uc
            write (*, *) "pressure_estat:", virial_estat/3.d0/v_uc
            write (*, *) "pressure_lp: ", virial_lp/3.d0/v_uc
            write (*, *) "pressure_2b: ", virial_2b/3.d0/v_uc
            write (*, *) "pressure_3b: ", virial_3b/3.d0/v_uc
            write (*, *) "pressure_core_pot: ", virial_core_pot/3.d0/v_uc
         end if

         if (params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc) then
#ifdef _MPIF90
            IF (rank == 0) then
#endif
               !       Write energy and forces if we're just doing static predictions
               !       The masses should be divided by 103.6426965268d0 to have amu units, but
               !       since masses is not allocated for single point calculations, it would
               !       likely lead to a segfault
               call wrap_pbc(positions(1:3, 1:n_sites), a_box&
                 &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                 & c_box/dfloat(indices(3)))
               call get_xyz_energy_string(energies_soap, energies_2b,&
                 & energies_3b, energies_core_pot, energies_vdw, energies_estat, energies_exp&
                 &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                 & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                 & params%do_structure_factor, params%do_xrd, params%do_nd, string)

               call write_extxyz(n_sites, -n_xyz, md_time, time_step,&
                 & instant_temp, instant_pressure, a_box&
                 &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                 & c_box/dfloat(indices(3)), virial, xyz_species,&
                 & positions(1:3, 1:n_sites), velocities, forces,&
                 & energies(1:n_sites), masses, params&
                 &%write_property, params%write_array_property,&
                 & params%write_local_properties, local_property_labels, local_properties, &
                 & fix_atom, "trajectory_out.xyz", string, .false.)

#ifdef _MPIF90
            END IF
#endif
         end if
      else
#ifdef _MPIF90
         IF (rank == 0) then
#endif
            !     Do nothing
            write (*, *) '                                       |'
            write (*, *) 'You didn''t ask me to do anything!      |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
#ifdef _MPIF90
         END IF
#endif
      end if
      !**************************************************************************

      !**************************************************************************
      !   Do MD stuff here. Moved to src/turbogap_md.f90; the rank guard and the
      !   position broadcast moved with it.
      call compute_md(params, rank, ierr, n_sites, md_istep, md_time, time_step, positions, &
                      positions_prev, positions_diff, velocities, forces, forces_prev, masses, xyz, &
                      xyz_species, a_box, b_box, c_box, indices, v_uc, virial, energy, energy_prev, &
                      energies, energies_soap, energies_2b, energies_3b, energies_core_pot, &
                      energies_estat, energies_vdw, energies_lp, energies_exp, energies_pdf, energies_sf, &
                      energies_xrd, energies_nd, local_properties, local_property_labels, instant_temp, &
                      instant_pressure, instant_pressure_prev, e_kin, e_kinetic, kb, evpera3tobar, &
                      fix_atom, exit_loop, rebuild_neighbors_list, i_image, i_nested, n_pos, string, &
                      time, gd_box_do_pos, gd_istep, restart_box_optim, &
                      target_temp, time_step_prev)

      !**************************************************************************
      !   Nested sampling
      !   PUT THIS INTO A MODULE!!!!!!!!!!!!!!

      !   This runs at the beginning to read in the initial images
      if (params%do_nested_sampling .and. n_xyz > i_image .and. .not. params%do_md) then
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
         velocities = 0.d0
         call from_properties_to_image(images(i_image), positions, velocities, masses, &
                                       forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                       species, species_supercell, n_sites, indices, fix_atom, &
                                       xyz_species, xyz_species_supercell, local_properties)
      end if

      !   This handles the nested sampling iterations after all images have
      !   been read and their energies computed
      if (params%do_nested_sampling .and. .not. repeat_xyz) then
         if (i_nested == 0) then
            md_istep = -1
            params%write_xyz = params%md_nsteps
            params%do_md = .true.
            if (rank == 0) then
               write (*, *) '                                       |'
               write (*, *) 'Running nested sampling algorithm with |'
               write (*, '(1X,I6,A)') n_xyz, ' walkers.                        |'
               write (*, *) '                                       |'
               write (*, *) 'Target pressure in nested sampling:    |'
               write (*, '(A,ES15.7,A)') ' P = ', params%p_nested, ' bar.               |'
               write (*, *) '                                       |'
               write (*, *) '[P = 0 means total energy, rather than |'
               write (*, *) 'total enthalphy, simulation]           |'
            end if
         end if
         !     At the end of the MD/MC moves we add the image to the pool if its energy has decreased
         if (md_istep == params%md_nsteps) then
            md_istep = -1
            velocities = 0.d0
            !       Unit cell volume
            v_uc = dot_product(cross_product(a_box, b_box), c_box)/(dfloat(indices(1)*indices(2)*indices(3)))
            !       We check enthalpy, not internal energy (they are the same for P = 0)
            if (energy + E_kinetic + params%p_nested/eVperA3tobar*v_uc < e_max) then
               call from_properties_to_image(images(i_image), positions, velocities, masses, &
                                             forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                             species, species_supercell, n_sites, indices, fix_atom, &
                                             xyz_species, xyz_species_supercell, local_properties)
            end if
         end if
         !     This selects the highest energy image from the pool
         if (md_istep == -1 .and. i_nested < params%n_nested) then
            i_nested = i_nested + 1
            rebuild_neighbors_list = .true.
            i_max = 0
            e_max = -1.d100
            do i = 1, n_xyz
               v_uc = dot_product(cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box)/ &
                      (dfloat(images(i)%indices(1)*images(i)%indices(2)*images(i)%indices(3)))
               !         We check enthalpy, not potential energy (they are the same for P = 0)
               if (images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*v_uc > e_max) then
                  e_max = images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*v_uc
                  i_max = i
               end if
            end do
            i_image = i_max
            deallocate (positions, velocities, masses, forces, species, &
                        species_supercell, fix_atom, xyz_species, xyz_species_supercell)
            !       Make a copy of a randonmly chosen image which is not i_image
            if (n_xyz == 1) then
               i = i_image
            else
               i = i_image
               do while (i == i_image)
                  call random_number(ranf)
                  i = 1 + floor((n_xyz)*ranf)
                  ! i2 = random_
                  ! i = mod(irand(), n_xyz) + 1
               end do
            end if
            if (rank == 0) then
               counter = 1
               !          write(*,*)
               write (*, *) '                                       |'
               write (*, '(A,I8,A,I8,A)') "Nested sampling iter.:", i_nested, "/", params%n_nested, " |"
               write (*, '(A,I8,A)') " - Highest enthalpy walker:    ", i_image, " |"
               write (*, '(A,I8,A)') " - Walker selected for cloning:", i, " |"
               write (*, '(A,F15.7,A)') " - Max. enthalpy: ", e_max, " eV |"
            end if
            call from_image_to_properties(images(i), positions, velocities, masses, &
                                          forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                          species, species_supercell, n_sites, indices, fix_atom, &
                                          xyz_species, xyz_species_supercell, local_properties)
            v_uc = dot_product(cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box)/ &
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
               call volume_preserving_strain_transformation(a_box, b_box, c_box, params%box_scaling_factor)
               ! Volume scaling
               call get_ns_unbiased_volume_proposal(1.d0 - params%nested_max_volume_change, &
                                                    1.d0 + params%nested_max_volume_change, n_sites, rand)
               params%box_scaling_factor = params%box_scaling_factor*(rand)**(1.d0/3.d0)
               ! Each MPI process has a different set of random numbers so we need to broadcast
#ifdef _MPIF90
               call mpi_bcast(params%box_scaling_factor, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
            end if
            !       This is the so-called total enthalpy Hamiltonian Montecarlo approach (with physical masses)
            !       We do not need to broadcast the velocities here since they get broadcasted later on; otherwise
            !       we would have to do it since each MPI rank may see a different random number
            call random_number(velocities)
            call remove_cm_vel(velocities(1:3, 1:n_sites), masses(1:n_sites))
            e_kin = 0.d0
            do i = 1, n_sites
               e_kin = e_kin + 0.5d0*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
            end do
            call random_number(rand)
            !        rand = rand * 4.d0/3.d0 - 1.d0/3.d0
            !        velocities = velocities / sqrt(e_kin) * sqrt(e_max - energy - params%p_nested/eVperA3tobar*v_uc + &
            !                                                     1.5d0*real(n_sites-1)*kB*params%t_extra*max(0.d0, rand))
            velocities = velocities/sqrt(e_kin)*sqrt(rand*(e_max - energy - params%p_nested/eVperA3tobar*v_uc))
         else if (i_nested == params%n_nested) then
            exit_loop = .true.
         end if
      end if

      !   NOW THIS IS HANDLED AT THE BEGINNING OF THE CODE WHEN WE CHECK IF THE NUMBER OF SITES HAS CHANGED
      !   Clean up
      !    deallocate( energies, energies_soap, energies_2b, energies_3b, energies_core_pot, this_energies, energies_vdw )
      !    if( params%do_forces )then
      !      deallocate( forces, forces_soap, forces_2b, forces_3b, forces_core_pot, this_forces, forces_vdw )
      !    end if
      !    if( any( soap_turbo_hypers(:)%has_vdw ) )then
      !      nullify( this_hirshfeld_v_pt )
      !      deallocate( this_hirshfeld_v, hirshfeld_v )
      !      if( params%do_forces )then
      !        nullify( this_hirshfeld_v_cart_der_pt )
      !        deallocate( this_hirshfeld_v_cart_der, hirshfeld_v_cart_der )
      !      end if
      !    end if

#ifdef _MPIF90
      IF (rank == 0) THEN
#endif

         if (params%do_mc) then
            if (mc_istep == params%mc_nsteps) then
               exit_loop = .true.
            else
               exit_loop = .false.
            end if

            if (.not. exit_loop .and. ( &
                (md_istep == -1) .or. &
                (params%do_md .and. ( &
                 (md_istep == params%md_nsteps) .or. &
                 ((abs(energy - energy_prev) < params%e_tol*dfloat(n_sites)) .and. (maxval(forces) < params%f_tol)) &
                 )))) then
               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               !time%mc(1) = MPI_wtime()
               call time_start(time%mc)

               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               if (mc_istep > 0) then
                  !       Evaluate the conditions for acceptance
                  !       > We have the mc conditions in mc.f90
                  !       > We care about comparing e_store to the energy of the new configuration based on the mc_movw

                  ! Reset the parameters for md / relaxation
                  if (params%do_md) then
                     md_istep = -1
                     params%do_md = .false.
                     do_mc_relax = .false.
                     ! Assume that the number of steps has already been set.
                  end if

                  if (.not. params%mc_hamiltonian) E_kinetic = 0.d0

                  call from_properties_to_image(images(i_trial_image), positions, velocities, masses, &
                                                forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                                species, species_supercell, n_sites, indices, fix_atom, &
                                                xyz_species, xyz_species_supercell, local_properties)

                  if (params%verb > 50) write (*, *) '.......................................|'
                  if (params%verb > 50) write (*, '(A,1X,I0)') ' MC Iteration:', mc_istep
                  if (params%verb > 50) write (*, '(A,1X,A)') '    Move type:', mc_move
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                    & Etot_prev:', images(i_current_image)&
                    &%energy + images(i_current_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                    & Etot_new :', images(i_trial_image)%energy &
                    &+ images(i_trial_image)%e_kin

                  v_uc = dot_product(cross_product(a_box, b_box), c_box)/(dfloat(indices(1)*indices(2)*indices(3)))

                  if (params%accessible_volume) then
                     call get_accessible_volume(v_uc, v_a_uc, species, params%radii)
                     if (params%verb > 50) write (*, '(A,F12.6,A,F12.6&
                       &,1X,A)') ' V_acc new: ', v_a_uc, ' A^3&
                       & V_acc old ', v_a_uc_prev, 'A^3 |'
                  else
                     v_a_uc = v_uc
                  end if

                  call get_mc_acceptance(mc_move, p_accept, &
                    energy + E_kinetic, &
                    images(i_current_image)%energy + images(i_current_image)%e_kin, &
                    params%t_beg, &
                    params%mc_mu(mc_mu_id), n_mc_species(mc_mu_id), v_uc, v_uc_prev,&
                    & v_a_uc, v_a_uc_prev, params&
                    &%masses_types(mc_id(mc_mu_id)), params%p_beg)

                  call random_number(ranf)

                  if (mc_move == "insertion") n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) + 1
                  if (mc_move == "removal") n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) - 1

                  !    ACCEPT OR REJECT
                  if (params%verb > 50) write (*, '(A,1X,A,1X,A,L4,1X&
                    &,A,ES12.6,1X,A,1X,ES12.6)') 'Is ',&
                    & trim(mc_move), 'accepted?', p_accept >&
                    & ranf, ' p_accept =', p_accept, ' ranf = ',&
                    & ranf

                  if (mc_istep == 1) then
                     open (unit=200, file="mc.log", status="unknown")
                     if (energy_exp > 0.d0) then
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
                  if (mc_istep > 1) then
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

                  if (energy_exp > 0.d0) then

                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                       mc_istep, trim(adjustl(mc_move)), p_accept > ranf, energy + E_kinetic, &
                       images(i_current_image)%energy +&
                       & images(i_current_image)%e_kin, energy_exp,&
                       & images(i_current_image)%energy_exp,&
                       & images(i_trial_image)%n_sites,&
                       & trim(temp_string2)
                  else
                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                        mc_istep, trim(adjustl(mc_move)), p_accept > ranf, energy + E_kinetic, &
                        images(i_current_image)%energy + images(i_current_image)%e_kin, &
                        images(i_trial_image)%n_sites, trim(temp_string2)

                  end if

                  if (mc_istep >= 1) close (200)

                  if (p_accept > ranf) then
                     !             Accept
                     ! Set variables
                     n_sites_prev = n_sites
                     v_uc_prev = v_uc
                     v_a_uc_prev = v_a_uc
                     virial_prev = virial
                     !   Assigning the default image with the accepted one
                     images(i_current_image) = images(i_trial_image)

                     if (params%n_mc_mu > 0) then
                        n_mc_species_prev = n_mc_species
                     end if

                  end if
                  instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic
                  instant_pressure = (kB*dfloat(n_sites - 1)*instant_temp&
                    &+ (virial(1, 1) + virial(2, 2) + virial(3, 3))/3.d0)&
                    &/v_uc*eVperA3tobar

                  if ((params%mc_write_xyz .or. mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                       modulo(mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') '&
                       & Writing mc_current.xyz and&
                       & mc_all.xyz '
                     call wrap_pbc(images(i_current_image)&
                       &%positions(1:3,&
                       & 1:images(i_current_image)%n_sites),&
                       & images(i_current_image)%a_box&
                       &/dfloat(indices(1)),&
                       & images(i_current_image)%b_box&
                       &/dfloat(indices(2)),&
                       & images(i_current_image)%c_box&
                       &/dfloat(indices(3)))
                     call get_xyz_energy_string(energies_soap, energies_2b,&
                       & energies_3b, energies_core_pot, energies_vdw, energies_estat, energies_exp&
                       &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                       & params%valid_pdf, params%valid_sf,&
                       & params%valid_xrd, params%valid_nd,&
                       & params%do_pair_distribution, params&
                       &%do_structure_factor, params%do_xrd,&
                       & params%do_nd, string)

                     call write_extxyz(images(i_current_image)%n_sites, 0, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                       images(i_current_image)%a_box/dfloat(indices(1)), &
                       images(i_current_image)%b_box/dfloat(indices(2)), &
                       images(i_current_image)%c_box/dfloat(indices(3)), &
                       virial_prev, images(i_current_image)%xyz_species, &
                       images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                       images(i_current_image)%velocities, &
                       images(i_current_image)%forces, &
                       images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                       images(i_current_image)%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & local_property_labels,&
                       & images(i_current_image)%local_properties&
                       &, images(i_current_image)%fix_atom,&
                       & "mc_current.xyz", string, .true.)

                     call write_extxyz(images(i_current_image)%n_sites, 1, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                       images(i_current_image)%a_box/dfloat(indices(1)), &
                       images(i_current_image)%b_box/dfloat(indices(2)), &
                       images(i_current_image)%c_box/dfloat(indices(3)), &
                       virial_prev, images(i_current_image)%xyz_species, &
                       images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                       images(i_current_image)%velocities, &
                       images(i_current_image)%forces, &
                       images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                       images(i_current_image)%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & local_property_labels,&
                       & images(i_current_image)%local_properties,&
                       & images(i_current_image)%fix_atom,&
                       & "mc_all.xyz", string, .false.)

                  end if

                  !          Add acceptance to the log file else dont

                  !time%mc(2) = MPI_wtime()
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
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                    &   'mc_reverse    = ', params%mc_reverse, '    &
                    & |'
                  if (params%verb > 50) write (*, '(1X,A,1X,F12.6,1X&
                    &,A)') 'mc_reverse_lambda = ', params&
                    &%mc_reverse_lambda, '|'

                  if (params%verb > 50) write (*, *) '                                       |'
                  ! t_beg must

                  if (.not. allocated(images) .and. .not. params%do_nested_sampling) then
                     allocate (images(1:2))
                  else if (.not. allocated(images) .and. params%do_nested_sampling) then
                     allocate (images(1:2*i_image))
                  end if

                  if (.not. allocated(mc_id) .and. params%n_mc_mu > 0) then
                     allocate (mc_id(1:params%n_mc_mu))
                     allocate (n_mc_species(1:params%n_mc_mu))
                     allocate (n_mc_species_prev(1:params%n_mc_mu))

                     mc_id = 1
                     n_mc_species = 0

                     !    get the mc species types

                     do j = 1, params%n_mc_mu
                        do i = 1, n_species
                           if (params%species_types(i) == params%mc_species(j)) then
                              mc_id(j) = i
                           end if
                        end do
                     end do
                  end if

                  !       Now use the image construct to store this as the image to compare to
                  call from_properties_to_image(images(i_current_image), positions, velocities, masses, &
                                                forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                                species, species_supercell, n_sites, indices, fix_atom, &
                                                xyz_species, xyz_species_supercell, local_properties)

                  instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic
                  instant_pressure = (kB*dfloat(n_sites - 1)*instant_temp&
                    &+ (virial(1, 1) + virial(2, 2) + virial(3, 3))/3.d0)&
                    &/v_uc*eVperA3tobar

                  if ((mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                       modulo(mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') ' Writing mc_current.xyz and mc_all.xyz '
                     call wrap_pbc(images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                                   images(i_current_image)%a_box/dfloat(indices(1)), &
                                   images(i_current_image)%b_box/dfloat(indices(2)), &
                                   images(i_current_image)%c_box/dfloat(indices(3)))
                     call get_xyz_energy_string(energies_soap, energies_2b,&
                       & energies_3b, energies_core_pot, energies_vdw, energies_estat, energies_exp&
                       &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                       & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                       & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                     call write_extxyz(images(i_current_image)%n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                       images(i_current_image)%a_box/dfloat(indices(1)), &
                       images(i_current_image)%b_box/dfloat(indices(2)), &
                       images(i_current_image)%c_box/dfloat(indices(3)), &
                       virial_prev, images(i_current_image)%xyz_species, &
                       images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                       images(i_current_image)%velocities, &
                       images(i_current_image)%forces, &
                       images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                       images(i_current_image)%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & local_property_labels, images(i_current_image)%local_properties&
                       &, images(i_current_image)%fix_atom,&
                       & "mc_current.xyz", string, .true.)

                     call write_extxyz(images(i_current_image)%n_sites, 1, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                       images(i_current_image)%a_box/dfloat(indices(1)), &
                       images(i_current_image)%b_box/dfloat(indices(2)), &
                       images(i_current_image)%c_box/dfloat(indices(3)), &
                       virial_prev, images(i_current_image)%xyz_species, &
                       images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                       images(i_current_image)%velocities, &
                       images(i_current_image)%forces, &
                       images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                       images(i_current_image)%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & local_property_labels, images(i_current_image)%local_properties&
                       &, images(i_current_image)%fix_atom,&
                       & "mc_all.xyz", string, .true.)

                     v_uc_prev = dot_product(cross_product(a_box, b_box), c_box)/(dfloat(indices(1)*indices(2)*indices(3)))
                     if (params%accessible_volume) then
                        call get_accessible_volume(v_uc_prev, v_a_uc_prev, species, params%radii)
                     else
                        v_a_uc_prev = v_uc_prev
                     end if
                  end if

               end if

               !  Now start the mc logic: first, use the stored images properties
               call from_image_to_properties(images(i_current_image), positions, velocities, masses, &
                                             forces, a_box, b_box, c_box, energy, energies, energy_exp, E_kinetic, &
                                             species, species_supercell, n_sites, indices, fix_atom, &
                                             xyz_species, xyz_species_supercell, local_properties)

               call perform_mc_step(&
                 & positions, species, xyz_species, masses, fix_atom,&
                 & velocities, positions_prev, positions_diff, disp, d_disp, params%n_local_properties,&
                 & params%mc_acceptance, params%mc_mu_acceptance, local_properties, &
                 images(i_current_image)%local_properties, energies,&
                 & forces, forces_prev, n_sites, params%n_mc_mu, mc_mu_id, n_mc_species,&
                 & mc_move, params%mc_species,&
                 & params%mc_move_max, params%mc_min_dist, params%mc_lnvol_max, params&
                 &%mc_types, params%masses_types, species_idx,&
                 & images(i_current_image)%positions,&
                 & images(i_current_image)%species,&
                 & images(i_current_image)%xyz_species,&
                 & images(i_current_image)%fix_atom,&
                 & images(i_current_image)%masses, a_box(1:3), b_box(1:3),&
                 & c_box(1:3), indices, params%do_md, params%mc_relax,&
                 & md_istep, mc_id, E_kinetic, instant_temp, params%t_beg,&
                 & params%n_mc_swaps, params%mc_swaps, params%mc_swaps_id, &
                 & params%species_types, params%mc_hamiltonian,&
                 & params%n_mc_relax_after, params&
                 &%mc_relax_after, do_mc_relax, params%verb)

               rebuild_neighbors_list = .true.
               ! end if

               ! NOTE: the species_supercell and xyz_species_supercell are
               ! not commensurate with the new image as these have not been
               ! calculated. If reading from an outputted xyz file, then it
               ! should be okay but really the new atoms should be added to
               ! the supercell in the usual way, but for convenience, one has
               ! not done that.

               if (params%mc_relax .and. do_mc_relax) then
                  ! Set the parameters for relaxatrino
                  md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_relax_opt
                  params%md_nsteps = params%mc_nrelax
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if
               ! If doing md, don't relax
               if (mc_move == 'md') then
                  ! Set the parameters for relaxatrino
                  md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_hybrid_opt
                  params%md_nsteps = temp_md_nsteps
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if

               if ((params%mc_write_xyz .or. mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                    modulo(mc_istep, params%write_xyz) == 0)) then

                  call wrap_pbc(positions(1:3, 1:n_sites), &
                                a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
                  call get_xyz_energy_string(energies_soap, energies_2b,&
                    & energies_3b, energies_core_pot, energies_vdw, energies_estat, energies_exp&
                    &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                    & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                    & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                  call write_extxyz(n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                    a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                    virial, xyz_species, &
                    positions(1:3, 1:n_sites), velocities, &
                    forces, energies(1:n_sites), masses, &
                    params%write_property, params&
                    &%write_array_property, params&
                    &%write_local_properties,&
                    & local_property_labels, local_properties&
                    &, fix_atom, mc_file, string, .true.)
               end if
               ! As we have moved/added/removed, we must check the supercell and  broadcast the results

               call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                             n_sites, .true., fix_atom, params%t_beg, &
                             params%write_array_property(6), .true.)

            else
               if (mc_move == 'md') then
                  if (params%print_progress .and. md_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     update_bar = params%md_nsteps/36
                     if (update_bar < 1) then
                        update_bar = 1
                     end if
                     counter = 1
                  else if (md_istep == params%md_nsteps - 1 .or. &
                           (abs(energy - energy_prev) < params%e_tol*dfloat(n_sites) .and. &
                            maxval(forces) < params%f_tol) .and. md_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. counter == update_bar .and. md_istep < params%md_nsteps - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(md_istep + 1)/params%md_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(md_istep + 1)/params%md_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     counter = 1
                  else
                     counter = counter + 1
                  end if

                  if (params%mc_hamiltonian) then
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                       &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step: H =&
                       & T + V = ", energy + E_kinetic, ",&
                       & iteration ", md_istep, "/", params&
                       &%md_nsteps
                  else
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                       &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step:&
                       & energy = ", energy, ", iteration ",&
                       & md_istep, "/", params%md_nsteps
                  end if

                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') '&
                    & core_pot energy:', sum(energies_core_pot),&
                    & 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') '&
                    & vdw energy:', sum(energies_vdw), 'eV |'
                  if (params%verb > 50 .and. valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                    & sum(energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                    & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                    & sum(energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                    & sum(energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                    & sum(energies_nd), 'eV |'

               else
                  if (params%print_progress .and. mc_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     update_bar = params%mc_nsteps/36
                     if (update_bar < 1) then
                        update_bar = 1
                     end if
                     counter = 1
                  else if (mc_istep == params%mc_nsteps - 1 .and. mc_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. counter == update_bar .and. mc_istep < params%mc_nsteps - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(mc_istep + 1)/params%mc_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(mc_istep + 1)/params%mc_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     counter = 1
                  else
                     counter = counter + 1
                  end if

                  if (params%verb > 50 .and. do_mc_relax) write (*, '(1X,A,1X,F20.8,1X,A&
                    &,1X,I8,1X,A,1X,I8)') "MC Relax md step: energy &
                    &= ", energy, ", iteration ", md_istep, "/",&
                    & params%mc_nrelax
                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(energies_core_pot), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(energies_vdw), 'eV |'
                  if (params%verb > 50 .and. valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                    & sum(energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                    & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                    & sum(energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                    & sum(energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                    & sum(energies_nd), 'eV |'

               end if
            end if
         end if

#ifdef _MPIF90
      END IF
#endif

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

#ifdef _MPIF90
      IF (params%do_mc .and. md_istep == -1 .and. rank == 0) THEN
         n_pos = size(positions, 2)
         n_sp = size(xyz_species, 1)
         n_sp_sc = size(xyz_species_supercell, 1)
      END IF
       !! time%mpi(1)=MPI_Wtime()
      !call get_time( ! time%mpi(1) )

      call time_start(time%mpi)

      call mpi_bcast(n_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sp_sc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(params%do_md, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(md_istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       !! time%mpi(2)=MPI_Wtime()
      !call get_time( ! time%mpi(2) )

      call time_end(time%mpi)
      IF (rank /= 0) THEN !.and. (mc_move == "insertion" .or. mc_move == "removal")
         if (allocated(positions)) deallocate (positions)
         allocate (positions(1:3, n_pos))
         if (params%do_md .or. params%do_nested_sampling .or. params%do_mc) then
            if (allocated(velocities)) deallocate (velocities)
            allocate (velocities(1:3, n_pos))
            if (allocated(masses)) deallocate (masses)
            allocate (masses(1:n_sp))
            if (allocated(fix_atom)) deallocate (fix_atom)
            allocate (fix_atom(1:3, 1:n_sp))

            ! if(allocated(forces_prev))deallocate(forces_prev)
            ! allocate( forces_prev(1:3, 1:n_sites) )
            ! if(allocated(positions_prev))deallocate(positions_prev)
            ! allocate( positions_prev(1:3, 1:n_sites) )
            ! if(allocated(positions_diff))deallocate(positions_diff)
            ! allocate( positions_diff(1:3, 1:n_sites) )
            ! positions_diff = 0.d0

         end if
         if (allocated(xyz_species)) deallocate (xyz_species)
         allocate (xyz_species(1:n_sp))
         if (allocated(species)) deallocate (species)
         allocate (species(1:n_sp))
         if (allocated(xyz_species_supercell)) deallocate (xyz_species_supercell)
         allocate (xyz_species_supercell(1:n_sp_sc))
         if (allocated(species_supercell)) deallocate (species_supercell)
         allocate (species_supercell(1:n_sp_sc))
      END IF

      !time%mpi_positions(1) = MPI_wtime()
      call time_start(time%mpi_positions)

      call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (params%do_md .or. params%do_nested_sampling .or. params%do_mc) then
         call mpi_bcast(velocities, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(masses, n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(fix_atom, 3*n_sp, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
      call mpi_bcast(xyz_species, 8*n_sp, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(xyz_species_supercell, 8*n_sp_sc, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(species, n_sp, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(species_supercell, n_sp_sc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(n_sites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      !time%mpi_positions(2) = MPI_wtime()
      call time_end(time%mpi_positions)
#endif
      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(do_list)) deallocate (do_list)
         allocate (do_list(1:n_sites))
         do_list = .true.
      end if
      !

      !     !     !     !     call cpu_time(time1)
      call get_time(time1)

#ifdef _MPIF90
      !   Parallel neighbors list build
      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

      if (rebuild_neighbors_list) then
         deallocate (rjs, xyz, thetas, phis, neighbor_species)
         deallocate (neighbors_list, n_neigh)
#ifdef _MPIF90
         deallocate (n_neigh_local)
#endif
      end if
      if ((params%do_nested_sampling .and. .not. params%do_mc) .and. &
          (params%do_md .and. (md_istep == params%md_nsteps .or. exit_loop))) then
         deallocate (positions, xyz_species, xyz_species_supercell, species, species_supercell, do_list)
         if (allocated(velocities)) deallocate (velocities)
      end if
      if (params%do_mc .and. params%do_md) then
         if (params%do_mc .and. (mc_istep == params%mc_nsteps .or. exit_loop)) then
            deallocate (positions, xyz_species, xyz_species_supercell, species, species_supercell, do_list)
            if (allocated(velocities)) deallocate (velocities)
         end if
      end if

      if ((params%do_md .and. .not. params%do_mc) .and. &
          (md_istep == params%md_nsteps .or. exit_loop) .and. rank == 0) then
         deallocate (positions_prev, forces_prev)
      end if
      if (params%do_mc .and. (mc_istep == params%mc_nsteps .or. exit_loop) .and. rank == 0) then
         if (allocated(forces_prev)) deallocate (forces_prev)
         if (allocated(positions_prev)) deallocate (positions_prev)
      end if

      if (params%exp_forces .and. (md_istep == params%md_nsteps .or.&
        & mc_istep == params%mc_nsteps .or. exit_loop)) then
         do i = 1, params%n_exp
            if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
            if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
            if (allocated(params%exp_data(i)%y_pred)) deallocate (params%exp_data(i)%y_pred)
         end do
      end if

      if (.not. params%do_mc) n_sites_prev = n_sites
      n_atom_pairs_by_rank_prev = n_atom_pairs_by_rank(rank + 1)

#ifdef _MPIF90
      call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
      if (exit_loop) exit
      ! End of loop through structures in the xyz file or MD steps
   end do

   if (params%do_md .or. params%do_prediction .or. params%do_mc) then
       !! time2=MPI_Wtime()
      !call get_time( ! time2 )

      call get_time(time2)

#ifdef _MPIF90
      IF (rank == 0) then
#endif

         if (params%do_md .and. .not. params%do_nested_sampling) then
            !      write(*,'(A)')'] |'
            !      write(*,*)
            write (*, *) '                                       |'
            !      write(*,'(I8,A,F13.3,A)') params%md_nsteps, ' MD steps:', time2-time3, ' seconds |'
            write (*, '(I8,A,F13.3,A)') md_istep, ' MD steps:', time2 - time3, ' seconds |'
         end if
         if (params%do_mc) then
            !      write(*,'(A)')'] |'
            write (*, *)
            write (*, *) '                                       |'
            !      write(*,'(I8,A,F13.3,A)') params%md_nsteps, ' MD steps:', time2-time3, ' seconds |'
            write (*, '(I8,A,F13.3,A)') mc_istep, ' MC steps:', time2 - time3, ' seconds |'
         end if

         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *     Read input:', time%read_input(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' * Read XYZ files:', time%read_xyz(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' * Neighbor lists:', time%neigh(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' *  GAP desc/pred:', time%gap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - soap_turbo:', time%soap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - local_prop:', time%local_prop(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - lolo__soap:', time%soap_solo(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - get___soap:', time%get_soap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - lin__turbo:', time%soap_lin(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         2b:', time%gap_2b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         3b:', time%gap_3b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -   core_pot:', time%gap_core_pot(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -        vdw:', time%vdw(3), ' seconds |'
         if (valid_xps .or. params%do_pair_distribution .or. params&
           &%do_structure_factor .or. params%do_xrd .or. params%do_nd) write (*, '(A&
           &,F13.3,A)') ' *  Exp. pred.   :', time%pdf(3) + time%sf(3) + time%xrd(3) + time%nd(3), ' seconds&
           & |'
         if (valid_xps) write (*, '(A,F13.3,A)') '     -        xps:',&
           & time%xps(3), ' seconds |'
         if (params%do_pair_distribution) write (*, '(A,F13.3,A)') '     -        pdf:', time%pdf(3), ' seconds |'
         if (params%do_structure_factor) write (*, '(A,F13.3,A)') '     -         sf:', time%sf(3), ' seconds |'
         if (params%do_xrd) write (*, '(A,F13.3,A)') '     -        xrd:', time%xrd(3), ' seconds |'
         if (params%do_nd) write (*, '(A,F13.3,A)') '     -         nd:', time%nd(3), ' seconds |'

         if ((params%estat_method /= "none") .and. params%do_prediction) &
            write (*, '(A,F13.3,A)') '     -      estat:', time%estat(3), ' seconds |'
         if (params%do_md) then
            write (*, '(A,F13.3,A)') ' *  MD algorithms:', time%md(3), ' seconds |'
         end if
         if (params%do_mc) then
            write (*, '(A,F13.3,A)') ' *  MC algorithms:', time%mc(3), ' seconds |'
         end if

#ifdef _MPIF90
         write (*, '(A,F13.3,A)') ' *  MPI comms.   :', time%mpi(3) + time%mpi_positions(3) + time%mpi_ef(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -  pos & vel:', time%mpi_positions(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - E & F brc.:', time%mpi_ef(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -  MPI misc.:', time%mpi(3), ' seconds |'
#endif
!       Miscellaneous is what the parent buckets do not account for.  It used
!       to be written out here as one long subtraction, which is how it came to
!       subtract time%gap and the mpi_ef reduce nested inside it and print a
!       negative number.  sum_times owns the list now (src/timing.f90), so the
!       set summed here and the set declared as parents there cannot disagree.
         time%total(3) = time2 - time3
         write (*, '(A,F13.3,A)') ' *  Miscellaneous:', time%total(3) - sum_times(time), ' seconds |'
         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *     Total time:', time2 - time3, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
#ifdef _MPIF90
      END IF
#endif
   end if

   do i = 1, n_soap_turbo
      if (.not. soap_turbo_hypers(i)%recompute_basis) then
         call gpu_free_async(soap_turbo_hypers(i)%W_d, gpu_stream)
         call gpu_free_async(soap_turbo_hypers(i)%S_d, gpu_stream)
         call gpu_free_async(soap_turbo_hypers(i)%multiplicity_array_d, gpu_stream)
      end if
   end do

   if (allocated(fix_atom)) deallocate (fix_atom)
   if (allocated(positions)) deallocate (positions)
   if (allocated(velocities)) deallocate (velocities)
   if (allocated(positions_diff)) deallocate (positions_diff)

   if (allocated(energies)) deallocate (energies)
   if (allocated(energies_soap)) deallocate (energies_soap)
   if (allocated(energies_2b)) deallocate (energies_2b)
   if (allocated(energies_3b)) deallocate (energies_3b)
   if (allocated(energies_core_pot)) deallocate (energies_core_pot)
   if (allocated(energies_vdw)) deallocate (energies_vdw)
   if (allocated(energies_exp)) deallocate (energies_exp)
   if (allocated(energies_lp)) deallocate (energies_lp)
   if (allocated(energies_pdf)) deallocate (energies_pdf)
   if (allocated(energies_sf)) deallocate (energies_sf)
   if (allocated(energies_xrd)) deallocate (energies_xrd)
   if (allocated(energies_nd)) deallocate (energies_nd)

   if (allocated(this_energies)) deallocate (this_energies)
   if (allocated(this_energies_vdw)) deallocate (this_energies_vdw)
   if (allocated(this_energies_lp)) deallocate (this_energies_lp)
   if (allocated(this_energies_pdf)) deallocate (this_energies_pdf)
   if (allocated(this_energies_sf)) deallocate (this_energies_sf)
   if (allocated(this_energies_xrd)) deallocate (this_energies_xrd)
   if (allocated(this_energies_nd)) deallocate (this_energies_nd)

   if (allocated(forces)) deallocate (forces)
   if (allocated(forces_soap)) deallocate (forces_soap)
   if (allocated(forces_2b)) deallocate (forces_2b)
   if (allocated(forces_3b)) deallocate (forces_3b)
   if (allocated(forces_core_pot)) deallocate (forces_core_pot)
   if (allocated(forces_vdw)) deallocate (forces_vdw)
   if (allocated(forces_lp)) deallocate (forces_lp)
   if (allocated(forces_pdf)) deallocate (forces_pdf)
   if (allocated(forces_sf)) deallocate (forces_sf)
   if (allocated(forces_xrd)) deallocate (forces_xrd)
   if (allocated(forces_nd)) deallocate (forces_nd)

   if (allocated(this_forces)) deallocate (this_forces)
   if (allocated(this_forces_vdw)) deallocate (this_forces_vdw)
   if (allocated(this_forces_lp)) deallocate (this_forces_lp)
   if (allocated(this_forces_pdf)) deallocate (this_forces_pdf)
   if (allocated(this_forces_sf)) deallocate (this_forces_sf)
   if (allocated(this_forces_xrd)) deallocate (this_forces_xrd)
   if (allocated(this_forces_nd)) deallocate (this_forces_nd)

   if (allocated(local_properties)) deallocate (local_properties)
   if (allocated(local_properties_cart_der)) deallocate (local_properties_cart_der)
   if (allocated(this_local_properties)) deallocate (this_local_properties)
   if (allocated(this_local_properties_cart_der)) deallocate (this_local_properties_cart_der)

   if (allocated(soap_turbo_hypers)) deallocate (soap_turbo_hypers)
   if (allocated(distance_2b_hypers)) deallocate (distance_2b_hypers)
   if (allocated(angle_3b_hypers)) deallocate (angle_3b_hypers)
   if (allocated(core_pot_hypers)) deallocate (core_pot_hypers)

   !call destroy_cublas_handle(cublas_handle)

   ! write(*,*) "Starting dummy kernel"

     !! ! dut1= MPI_Wtime()
   !call get_time( ! ! dut1 )

   !  call get_time( ! dut1 )

   ! do n_ii=1,10

   !      do i_ii=1, 1024
   !      do j_jj=1,1024
   !      CC(i_ii,j_jj)=0.0
   !      do k_ii=1,1024
   !      CC(i_ii,j_jj)=CC(i,j)+AA(i_ii,k_ii)*BB(k_ii,j_jj)
   !      enddo
   !      enddo
   !      enddo
   ! enddo
   ! ! ! dut2= MPI_Wtime()
   ! call get_time( ! ! dut2 )

   !call get_time( ! dut2 )

   ! write(*,*) "Ending dummy region"
   ! write(*,*) "Time spent in dummy region", dut2-dut1

   deallocate (n_atom_pairs_by_rank)
   if (allocated(n_local_properties_mpi)) deallocate (n_local_properties_mpi)
   if (allocated(local_properties_n_sparse_mpi_soap_turbo)) deallocate (local_properties_n_sparse_mpi_soap_turbo)
   if (allocated(local_properties_dim_mpi_soap_turbo)) deallocate (local_properties_dim_mpi_soap_turbo)
   if (allocated(has_local_properties_mpi)) deallocate (has_local_properties_mpi)

   if (allocated(local_property_labels)) deallocate (local_property_labels)
   if (allocated(local_property_indexes)) deallocate (local_property_indexes)
   if (allocated(do_list)) deallocate (do_list)
   if (allocated(params%write_local_properties)) deallocate (params%write_local_properties)

#ifdef _MPIF90
   IF (rank == 0) then
#endif
      write (*, *) '                                       |'
      write (*, *) 'End of execution                       |'
      write (*, *) '_______________________________________/'
#ifdef _MPIF90
   END IF
#endif

   !write(*,*) "    - lin__turbo:", time%soap_lin, rank

#ifdef _MPIF90
   call mpi_finalize(ierr)
#endif

   call gpu_context_finalize(params, n_omp)

end program turbogap
