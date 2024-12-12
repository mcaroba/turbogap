! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap.f90, is copyright (c) 2019-2023, Miguel A. Caro and
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

  use neighbors
  use soap_turbo_desc
  use gap
  use read_files
  use md
  use adaptive_time			! for adaptive time simulation (TurboGAP will use these five modules for radiation cascades)
  use electronic_stopping		! for electronic stopping correction in radiation cascades
  use eph_fdm				! for T - dependent parameters - elec. stop. - eph model
  use eph_beta				! for the atomic electronic densities  - elec. stop. - eph model
  use eph_electronic_stopping		! for electronic stopping based in radiation cascades on the eph model
  use mc
  use gap_interface
  use types
  use vdw
  use exp_utils
  use exp_interface
  use soap_turbo_functions
#ifdef _MPIF90
  use mpi
  use mpi_helper
#endif
  use bussi
  use xyz_module

  implicit none


  !**************************************************************************
  ! Variable definitions
  !
  real*8, allocatable :: rjs(:), thetas(:), phis(:), xyz(:,:), sph_temp(:), sph_temp3(:,:)
  real*8, allocatable :: positions(:,:), positions_prev(:,:), soap(:,:), soap_cart_der(:,:,:), &
       positions_diff(:,:), forces_prev(:,:), frac_positions(:,:)
  real*8 :: rcut_max, a_box(1:3), b_box(1:3), c_box(1:3), max_displacement, energy, energy_prev
  real*8 :: virial(1:3, 1:3), this_virial(1:3, 1:3), virial_soap(1:3, 1:3), virial_2b(1:3, 1:3), &
       virial_3b(1:3,1:3), virial_core_pot(1:3, 1:3), virial_vdw(1:3, 1:3), virial_lp(1:3,1:3), &
       this_virial_vdw(1:3, 1:3), this_virial_lp(1:3, 1:3), virial_pdf(1:3,1:3), this_virial_pdf(1:3,1:3), v_uc,&
       & virial_sf(1:3,1:3), this_virial_sf(1:3,1:3), &
       & virial_xrd(1:3,1:3), this_virial_xrd(1:3,1:3), &
       & virial_nd(1:3,1:3), this_virial_nd(1:3,1:3), &
       & v_uc_prev, v_a_uc, v_a_uc_prev, eVperA3tobar =&
       & 1602176.6208d0, ranf, ranv(1:3), disp(1:3), d_disp, &
       & e_mc_prev, p_accept, virial_prev(1:3, 1:3), sim_exp_pred,&
       & sim_exp_prev, sim_exp_pred_der(1:3)
  real*8, allocatable :: energies(:), forces(:,:), energies_soap(:),&
       & forces_soap(:,:), this_energies(:), this_forces(:,:),&
       & energies_2b(:), forces_2b(:,:), energies_3b(:), forces_3b(:&
       &,:), energies_core_pot(:), forces_core_pot(:,:), velocities(:&
       &,:), masses_types(:), masses(:),  hirshfeld_v_temp(:),&
       & masses_temp(:), sinc_factor_matrix(:,:), energies_exp(:)
!  real*8, allocatable, target :: this_hirshfeld_v(:), this_hirshfeld_v_cart_der(:,:)
!  real*8, pointer :: this_hirshfeld_v_pt(:), this_hirshfeld_v_cart_der_pt(:,:)

  real*8, allocatable, target :: local_properties(:,:), local_properties_cart_der(:,:,:)
  ! Have one rank lower for the pointer, such that it just relates to a sub array of the local properties/cart_der
  real*8, pointer :: local_properties_pt(:), local_properties_cart_der_pt(:,:)
!  real*8, pointer :: hirshfeld_v(:), hirshfeld_v_cart_der(:,:)
  real*8, allocatable, target :: this_local_properties(:,:), this_local_properties_cart_der(:,:,:)
  real*8, pointer :: this_local_properties_pt(:,:), this_local_properties_cart_der_pt(:,:,:)
  real*8, allocatable ::  y_i_pred_all(:,:), moments(:), moments_exp(:)


  real*8, allocatable :: all_energies(:,:), all_forces(:,:,:), all_virial(:,:,:)
  real*8, allocatable :: all_this_energies(:,:), all_this_forces(:,:,:), all_this_virial(:,:,:)
  real*8 :: instant_temp, kB = 8.6173303d-5, E_kinetic=0.d0, E_kinetic_prev, time1, time2, time3, time_neigh, &
       time_gap, time_soap(1:3), time_2b(1:3), time_3b(1:3), time_read_input(1:3), time_read_xyz(1:3), &
       time_mpi(1:3) = 0.d0, time_core_pot(1:3), time_vdw(1:3),&
       & time_pdf(1:3), time_sf(1:3), time_xrd(1:3), time_nd(1:3), time_xps(1:3), time_mc(1:3), &
       & instant_pressure, lv(1:3,1:3), time_mpi_positions(1:3) =&
       & 0.d0, time_mpi_ef(1:3) = 0.d0, time_md(3) = 0.d0,&
       & instant_pressure_tensor(1:3, 1:3), time_step, md_time,&
       & instant_pressure_prev, wfac, wfac_temp, energy_exp
  integer, allocatable :: displs(:), displs2(:), counts(:), counts2(:), in_to_out_pairs(:), in_to_out_site(:), mc_id(:)
  integer :: update_bar, n_sparse, idx, gd_istep = 0, nprop
  logical, allocatable :: do_list(:), has_local_properties_mpi(:), fix_atom(:,:)
  logical :: rebuild_neighbors_list = .true., exit_loop = .true.,&
       & gd_box_do_pos = .true., restart_box_optim = .false.,&
       & valid_xps=.false., valid_vdw=.false.,  write_condition=.false., overwrite_condition=.false.
  character*1 :: creturn = achar(13)

  !! these decalarations are for time step and electronic stopping by different methods
  real*8 :: time_step_prev
  integer :: nrows
  real*8 :: cum_EEL = 0.0d0
  real*8, allocatable :: allelstopdata(:)
  type (EPH_Beta_class) :: ephbeta
  type (EPH_FDM_class) :: ephfdm
  type (EPH_LangevinSpatialCorrelation_class) :: ephlsc
  
  ! Clean up these variables after code refactoring !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, allocatable :: n_neigh(:), neighbors_list(:), alpha_max(:), species(:), species_supercell(:), &
       neighbor_species(:), sph_temp_int(:), der_neighbors(:), der_neighbors_list(:), &
       i_beg_list(:), i_end_list(:), j_beg_list(:), j_end_list(:),&
       & species_idx(:), n_neigh_out(:), n_local_properties_mpi(:),&
       & local_property_indexes(:), n_mc_species(:), n_mc_species_prev(:)
  integer :: n_sites, i, j, k, i2, j2, n_soap, k2, k3, l,&
       & n_sites_this, ierr, rank, ntasks, dim, n_sp, n_pos, n_sp_sc,&
       & this_i_beg, this_i_end, this_j_beg, this_j_end,&
       & this_n_sites_mpi, n_sites_prev = 0,&
       & n_atom_pairs_by_rank_prev=0, cPnz, n_pairs, n_all_sites,&
       & n_sites_out, n_local_properties_tot=0, n_lp_count=0,&
       & vdw_lp_index, core_be_lp_index, xps_idx

  integer :: l_max, n_atom_pairs, n_max, ijunk, central_species = 0,&
       & n_atom_pairs_total
  integer :: iostatus, counter = 0, counter2
  integer :: which_atom = 0, n_species = 1, n_species_actual, n_xyz, indices(1:3)
  integer :: radial_enhancement = 0
  integer :: md_istep, mc_istep, mc_mu_id=1, n_mc
  character*8, allocatable :: species_types_actual(:)
  character*1024, allocatable ::  local_property_labels(:), local_property_labels_temp(:), local_property_labels_temp2(:)
  logical :: repeat_xyz = .true., overwrite = .false., check_species,&
       & valid_local_properties=.false., label_in_list, do_mc_relax&
       &=.false.

  character*1024 :: filename, cjunk, file_compress_soap, file_alphas, file_soap, file_2b, file_alphas_2b, &
       file_3b, file_alphas_3b, file_gap = "none", mc_file = "mc_trial.xyz", string, temp_string, temp_string2
  character*64 :: keyword
  character*16 :: lattice_string(1:9)
  character*8 :: i_char
  character*8, allocatable ::  xyz_species(:), xyz_species_supercell(:), &
       species_type_temp(:)

  character*1 :: keyword_first

  ! This is the mode in which we run TurboGAP
  character*16 :: mode = "none"
  character*32 :: mc_move = "none", exp_output="none"

  ! Here we store the input parameters
  type(input_parameters) :: params

  ! These are the containers for the hyperparameters of descriptors and GAPs
  integer :: n_soap_turbo = 0, n_distance_2b = 0, n_angle_3b = 0, n_core_pot = 0, counter_lp_names=0, temp_md_nsteps
  real*8, parameter :: pi = acos(-1.0)
  type(soap_turbo), allocatable :: soap_turbo_hypers(:)
  type(distance_2b), allocatable :: distance_2b_hypers(:)
  type(angle_3b), allocatable :: angle_3b_hypers(:)
  type(core_pot), allocatable :: core_pot_hypers(:)

  !vdw crap
  real*8, allocatable :: v_neigh_vdw(:), energies_vdw(:), forces_vdw(:,:), this_energies_vdw(:), this_forces_vdw(:,:)
  real*8, allocatable :: v_neigh_lp(:), energies_lp(:), forces_lp(:,:), this_energies_lp(:), this_forces_lp(:,:)
  real*8, allocatable :: energies_pdf(:) , forces_pdf(:,:), this_energies_pdf(:), this_forces_pdf(:,:)
  real*8, allocatable :: energies_sf(:) , forces_sf(:,:), this_energies_sf(:), this_forces_sf(:,:)
  real*8, allocatable :: energies_xrd(:), forces_xrd(:,:), this_energies_xrd(:), this_forces_xrd(:,:)
  real*8, allocatable :: energies_nd(:), forces_nd(:,:), this_energies_nd(:), this_forces_nd(:,:)
  ! MPI stuff
  real*8, allocatable :: temp_1d(:), temp_1d_bis(:), temp_2d(:,:),&
       & pair_distribution_partial(:,:), pair_distribution_der(:,:), pair_distribution_partial_der(:,:,:), &
       & pair_distribution_partial_temp(:,:),&
       & pair_distribution_partial_temp_der(:,:,:),&
       & n_atoms_of_species(:), structure_factor_partial(:,:),&
       & structure_factor_partial_temp(:,:), structure_factor_partial_der(:,:,:),&
       & structure_factor_partial_temp_der(:,:), x_pair_distribution(:)&
       &, y_pair_distribution(:), y_pair_distribution_temp(:),&
       & x_structure_factor(:), x_structure_factor_temp(:),&
       & y_structure_factor(:), y_structure_factor_temp(:),&
       x_xrd(:), x_xrd_temp(:), y_xrd(:), y_xrd_temp(:), y_xrd_der(:,:,:), y_xrd_temp_der(:,:,:), &
       x_nd(:), x_nd_temp(:), y_nd(:), y_nd_temp(:), y_nd_der(:,:,:), y_nd_temp_der(:,:,:)
  integer, allocatable :: temp_1d_int(:), n_atom_pairs_by_rank(:),&
       & displ(:)
  integer, allocatable :: n_species_mpi(:), n_sparse_mpi_soap_turbo(:), dim_mpi(:), n_sparse_mpi_distance_2b(:), &
       n_sparse_mpi_angle_3b(:), n_mpi_core_pot(:),&
       & local_properties_n_sparse_mpi_soap_turbo(:),&
       & local_properties_dim_mpi_soap_turbo(:), n_neigh_local(:),&
       & compress_P_nonzero_mpi(:)
  integer :: i_beg, i_end, n_sites_mpi, j_beg, j_end, size_soap_turbo, size_distance_2b, size_angle_3b
  integer :: n_nonzero, q_beg, q_end
  logical, allocatable :: compress_soap_mpi(:)

  ! Nested sampling
  real*8 :: e_max, e_kin, rand, rand_scale(1:6), mag, n_total_cutoff, n_total_cutoff_temp, dq, target_temp
  integer :: i_nested, i_max, i_image, i_current_image=1, i_trial_image=2
  type(image), allocatable :: images(:), images_temp(:)
  type(exp_data_container) :: temp_exp_container
  character*32 :: implemented_exp_observables(1:5)

  real*8, allocatable :: x_xps(:), y_xps(:)

  implemented_exp_observables(1) = "xps"
  implemented_exp_observables(2) = "xrd"
  implemented_exp_observables(3) = "saxs"
  implemented_exp_observables(4) = "pair_distribution"
  implemented_exp_observables(5) = "structure_factor"

  !**************************************************************************

  !**************************************************************************
  ! Start recording the time
  call cpu_time(time1)
  time3 = time1
  ! Start random seed
  call srand(int(time1*1000))
  !**************************************************************************


  !**************************************************************************
  ! MPI stuff
#ifdef _MPIF90
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  !  allocate( displs(1:ntasks) )
  !  allocate( displs2(1:ntasks) )
  !  allocate( counts(1:ntasks) )
  !  allocate( counts2(1:ntasks) )
  ! allocate( displ(1:ntasks) )
#else
  rank = 0
  ntasks = 1
#endif
  allocate( n_atom_pairs_by_rank(1:ntasks) )
  !**************************************************************************



  !**************************************************************************
  ! Read the mode. It should be "soap", "predict" or "md"
  !
  call get_command_argument(1,mode)
  if( mode == "" .or. mode == "none" )then
     write(*,*) "ERROR: you need to run 'turbogap md' or 'turbogap predict'"
     stop
     ! THIS SHOULD BE FIXED, IN CASE THE USER JUST WANT TO OUTPUT THE SOAP DESCRIPTORS
     mode = "soap"
  end if
  !**************************************************************************




  !**************************************************************************
  ! Prints some welcome message and reads in the input file
  !
#ifdef _MPIF90
  IF( rank == 0 )THEN
#endif
  write(*,*)'_________________________________________________________________ '
  write(*,*)'                             _                                   \'
  write(*,*)' ___________            __   \\ /\        _____     ___   _____  |'
  write(*,*)'/____  ____/           / / /\|*\|*\/\    / ___ \   /   | |  _  \ |'
  write(*,*)'    / / __  __  __    / /  \********/   / /  /_/  / /| | | / | | |'
  write(*,*)'   / / / / / / / /_  / /__  \**__**/   / / ____  / / | | | |_/ / |'
  write(*,*)'  / / / / / / / __/ / ___ \ /*/  \*\  / / /_  / / /__| | |  __/  |'
  write(*,*)' / / / /_/ / / /   / /__/ / \ \__/ / / /___/ / / ____  | | |     |'
  write(*,*)'/_/_/_____/_/_/___/______/___\____/__\______/_/_/____|_|_|_|____ |'
  write(*,*)'_____________________________________________________________  / |'
  write(*,*)'*************************************************************|/  |'
  write(*,*)'                  Welcome to the TurboGAP code                   |'
  write(*,*)'                         Maintained by                           |'
  write(*,*)'                                                                 |'
  write(*,*)'                         Miguel A. Caro                          |'
  write(*,*)'                       mcaroba@gmail.com                         |'
  write(*,*)'                      miguel.caro@aalto.fi                       |'
  write(*,*)'                                                                 |'
  write(*,*)'          Department of Chemistry and Materials Science          |'
  write(*,*)'                     Aalto University, Finland                   |'
  write(*,*)'                                                                 |'
  write(*,*)'.................................................................|'
  write(*,*)'                                                                 |'
  write(*,*)'====================>>>>>  turbogap.fi  <<<<<====================|'
  write(*,*)'                                                                 |'
  write(*,*)'.................................................................|'
  write(*,*)'                                                                 |'
  write(*,*)'Contributors (code and methodology) in chronological order:      |'
  write(*,*)'                                                                 |'
  write(*,*)'Miguel A. Caro, Patricia Hernández-León, Suresh Kondati          |'
  write(*,*)'Natarajan, Albert P. Bartók, Eelis V. Mielonen, Heikki Muhli,    |'
  write(*,*)'Mikhail Kuklin, Gábor Csányi, Jan Kloppenburg, Richard Jana,     |'
  write(*,*)'Tigany Zarrouk                                                   |'
  write(*,*)'                                                                 |'
  write(*,*)'.................................................................|'
  write(*,*)'                                                                 |'
  write(*,*)'                     Last updated: Sep. 2024                     |'
  write(*,*)'                                        _________________________/'
  write(*,*)'.......................................|'
#ifdef _MPIF90
     write(*,*)'                                       |'
     write(*,*)'Running TurboGAP with MPI support:     |'
     write(*,*)'                                       |'
     write(*,'(A,I6,A)')' Running TurboGAP on ', ntasks, ' MPI tasks   |'
     write(*,*)'                                       |'
     write(*,*)'.......................................|'
#else
     write(*,*)'                                       |'
     write(*,*)'Running the serial version of TurboGAP |'
     write(*,*)'                                       |'
     write(*,*)'.......................................|'
#endif
#ifdef _MPIF90
  END IF
#endif
  !**************************************************************************







  !**************************************************************************
  ! Read input file and other files
  !
  time_read_input(3) = 0.d0
  call cpu_time(time_read_input(1))
  open(unit=10,file='input',status='old',iostat=iostatus)
  ! Check for existence of input file
#ifdef _MPIF90
  IF( rank == 0 )THEN
#endif
     write(*,*)'                                       |'
     write(*,*)'Checking input file...                 |'
#ifdef _MPIF90
  END IF
#endif
  if(iostatus/=0)then
     close(10)
#ifdef _MPIF90
     IF( rank == 0 )THEN
#endif
        write(*,*)'                                       |'
        write(*,*)'ERROR: input file could not be found   |  <-- ERROR'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
        write(*,*)'                                       |'
        write(*,*)'End of execution                       |'
        write(*,*)'_______________________________________/'
#ifdef _MPIF90
     END IF
     call mpi_finalize(ierr)
#endif
     stop
  end if
  !
  ! First, we look for n_species, which determines how we allocate the species-specific arrays
  do while(iostatus==0)
     read(10, *, iostat=iostatus) keyword
     keyword = trim(keyword)
     if(iostatus/=0)then
        exit
     end if
     keyword_first = keyword(1:1)
     if(keyword_first=='#' .or. keyword_first=='!')then
        continue
     else if(keyword=='n_species')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, n_species
        if( n_species < 1 )then
#ifdef _MPIF90
           IF( rank == 0 )THEN
#endif
              write(*,*)'                                       |'
              write(*,*)'ERROR: n_species must be > 0           |  <-- ERROR'
              write(*,*)'                                       |'
              write(*,*)'.......................................|'
#ifdef _MPIF90
           END IF
           call mpi_finalize(ierr)
#endif
           stop
        end if
     end if
  end do
  ! Let's look for those and other options in the input file
  rewind(10)
  call read_input_file(n_species, mode, params, rank)

!! If electronic stopping is required to be done, then read the stopping data file for once
!! Reading and storing the elctronic stopping data 
  if ( params%electronic_stopping ) then
		if (params%estop_filename == 'NULL') then
			write(*,*) "ERROR: No stopping data file is provided."
			stop
		else
			call read_electronic_stopping_file (n_species,params%species_types,params%estop_filename,nrows,allelstopdata)
		end if
  end if
!! -----------------------------				---- untill here for reading stopping data

  ! TEMPORARY ERROR, FIX THE UNDERLYING ISSUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _MPIF90
  IF( rank == 0 .and. ntasks > 1 .and. (params%write_soap .or. params%write_derivatives) )THEN
     write(*,*)'                                       |'
     write(*,*)'ERROR: writing of SOAP and/or SOAP     |  <-- ERROR'
     write(*,*)'derivatives cannot currently be done in|'
     write(*,*)'parallel. Try using the serial code or |'
     write(*,*)'running the MPI code on a single CPU.  |'
     write(*,*)'                                       |'
     write(*,*)'.......................................|'
     write(*,*)'                                       |'
     write(*,*)'End of execution                       |'
     write(*,*)'_______________________________________/'
  END IF
#endif
  !
  ! Second, we look for pot_file, which contains the GAP difinitions
  rewind(10)
  iostatus = 0
  do while(iostatus==0)
     read(10, *, iostat=iostatus) keyword
     keyword = trim(keyword)
     if(iostatus/=0)then
        exit
     end if
     keyword_first = keyword(1:1)
     if(keyword_first=='#' .or. keyword_first=='!')then
        continue
     else if(keyword=='pot_file')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, file_gap
        exit
     end if
  end do
  close(10)
  ! Now, read file_gap and register each GAP, including the hypers for its descriptor
  if( file_gap /= "none" )then
#ifdef _MPIF90
     IF( rank == 0 )then
#endif
        call read_gap_hypers(file_gap, &
             n_soap_turbo, soap_turbo_hypers, &
             n_distance_2b, distance_2b_hypers, &
             n_angle_3b, angle_3b_hypers, &
             n_core_pot, core_pot_hypers, &
             rcut_max, params%do_prediction, &
             params )
        !   Check if vdw_rcut is bigger
        if( params%vdw_rcut > rcut_max )then
           rcut_max = params%vdw_rcut
        end if
        if( params%xrd_rcut > rcut_max )then
           rcut_max = params%xrd_rcut
        end if
        if( params%nd_rcut > rcut_max )then
           rcut_max = params%nd_rcut
        end if
        if( params%pair_distribution_rcut > rcut_max )then
           rcut_max = params%pair_distribution_rcut
        end if

        !   We increase rcut_max by the neighbors buffer
        rcut_max = rcut_max + params%neighbors_buffer

     ! Check that the local properties we want to compute are valid,
     ! so things are commensurate between input file and .gap model

     ! Need to set the number of local properties, and get an array of labels and sizes for broadcasting

        write(*,*)'                                       |'
        write(*,*)'.......................................|'
        write(*,*)'                                       |'

        call get_irreducible_local_properties(params, n_local_properties_tot, n_soap_turbo, soap_turbo_hypers, &
             local_property_labels, local_property_labels_temp, local_property_labels_temp2, local_property_indexes, &
             valid_vdw, vdw_lp_index, core_be_lp_index, valid_xps, xps_idx )

        if( params%n_local_properties > 0)then
           write(*,*)'                                       |'
           write(*,*)' Irreducible local properties:         |'
           do i = 1, params%n_local_properties
              write(*,'(A41)') trim( local_property_labels(i) ) // ' |'
           end do


           allocate( params%write_local_properties(1:params%n_local_properties) )
           params%write_local_properties = .true.
        end if




#ifdef _MPIF90
     END IF
#endif
     !   THIS CHUNK HERE DISTRIBUTES THE INPUT DATA AMONG ALL THE PROCESSES
     !   Broadcast number of descriptors to other processes
#ifdef _MPIF90
     call cpu_time(time_mpi(1))
     call mpi_bcast(n_soap_turbo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_distance_2b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_angle_3b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_core_pot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_local_properties_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(params%n_local_properties, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     !   Broadcast the maximum cutoff distance
     call mpi_bcast(rcut_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(valid_xps, 1,&
          & MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(valid_vdw, 1,&
          & MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

     !   Processes other than 0 need to allocate the data structures on their own
     call cpu_time(time_mpi(2))
     time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
     allocate( n_species_mpi(1:n_soap_turbo) )
     allocate( n_sparse_mpi_soap_turbo(1:n_soap_turbo) )
     allocate( dim_mpi(1:n_soap_turbo) )
     allocate( n_local_properties_mpi(1:n_soap_turbo))
     if (n_local_properties_tot > 0) then
        allocate( local_properties_n_sparse_mpi_soap_turbo(1:n_local_properties_tot))
        allocate( local_properties_dim_mpi_soap_turbo(1:n_local_properties_tot))
        if (rank /= 0) allocate( local_property_indexes(1:n_local_properties_tot))
        if (rank /= 0) allocate( params%write_local_properties(1:params%n_local_properties))
     end if
     allocate( has_local_properties_mpi(1:n_soap_turbo) )
     allocate( compress_soap_mpi(1:n_soap_turbo) )
     allocate( n_sparse_mpi_distance_2b(1:n_distance_2b) )
     allocate( n_sparse_mpi_angle_3b(1:n_angle_3b) )
     allocate( n_mpi_core_pot(1:n_core_pot) )
     allocate( compress_P_nonzero_mpi(1:n_soap_turbo) )
     IF( rank == 0 )THEN
        n_species_mpi = soap_turbo_hypers(1:n_soap_turbo)%n_species
        n_sparse_mpi_soap_turbo = soap_turbo_hypers(1:n_soap_turbo)%n_sparse
        dim_mpi = soap_turbo_hypers(1:n_soap_turbo)%dim
        compress_soap_mpi = soap_turbo_hypers(1:n_soap_turbo)%compress_soap
        n_sparse_mpi_distance_2b = distance_2b_hypers(1:n_distance_2b)%n_sparse
        n_sparse_mpi_angle_3b = angle_3b_hypers(1:n_angle_3b)%n_sparse
        n_mpi_core_pot = core_pot_hypers(1:n_core_pot)%n
        compress_P_nonzero_mpi = soap_turbo_hypers(1:n_soap_turbo)%compress_P_nonzero

        has_local_properties_mpi = soap_turbo_hypers(1:n_soap_turbo)%has_local_properties
        n_local_properties_mpi = soap_turbo_hypers(1:n_soap_turbo)%n_local_properties


        ! Allocate the arrays have n_sparse and n_data
        if( any( soap_turbo_hypers(:)%has_local_properties ))then
           n_lp_count = 0
           do i = 1, n_soap_turbo
              if (n_local_properties_mpi(i) > 0)then
                 do j = 1, n_local_properties_mpi(i)
                    n_lp_count = n_lp_count + 1
                    local_properties_n_sparse_mpi_soap_turbo(n_lp_count) = &
                         & soap_turbo_hypers(i)%local_property_models(j)&
                         &%n_sparse

                    soap_turbo_hypers(i)%local_property_models(j)%dim = soap_turbo_hypers(i)%dim
                    local_properties_dim_mpi_soap_turbo(n_lp_count) = soap_turbo_hypers(i)%dim

                    ! Changing the dim to that of the descriptor, silly !
                    ! &
                    !      & soap_turbo_hypers(i)%local_property_models(j)&
                    !      &%dim

                 end do
              end if
           end do
        end if

     END IF
     call cpu_time(time_mpi(1))
     call mpi_bcast(n_species_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sparse_mpi_soap_turbo, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(dim_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


     if (n_local_properties_tot > 0)then
        call mpi_bcast(local_properties_n_sparse_mpi_soap_turbo,&
             & n_local_properties_tot, MPI_INTEGER, 0,&
             & MPI_COMM_WORLD, ierr)
        call mpi_bcast(local_properties_dim_mpi_soap_turbo,&
             & n_local_properties_tot, MPI_INTEGER, 0,&
             & MPI_COMM_WORLD, ierr)
        call mpi_bcast(local_property_indexes,&
             & n_local_properties_tot, MPI_INTEGER, 0,&
             & MPI_COMM_WORLD, ierr)
        call mpi_bcast(params%write_local_properties,&
             & params%n_local_properties, MPI_LOGICAL, 0,&
             & MPI_COMM_WORLD, ierr)
     end if

     call mpi_bcast(has_local_properties_mpi, n_soap_turbo,&
          & MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

     call mpi_bcast(n_local_properties_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(compress_soap_mpi, n_soap_turbo, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sparse_mpi_distance_2b, n_distance_2b, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sparse_mpi_angle_3b, n_angle_3b, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_mpi_core_pot, n_core_pot, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(compress_P_nonzero_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call cpu_time(time_mpi(2))
     time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)

     IF( rank /= 0 )THEN
        call allocate_soap_turbo_hypers(n_soap_turbo, n_species_mpi, n_sparse_mpi_soap_turbo, dim_mpi, &
             compress_P_nonzero_mpi,&
             & local_properties_n_sparse_mpi_soap_turbo,&
             & local_properties_dim_mpi_soap_turbo,&
             & has_local_properties_mpi, n_local_properties_mpi,&
             & compress_soap_mpi, soap_turbo_hypers)
        call allocate_distance_2b_hypers(n_distance_2b, n_sparse_mpi_distance_2b, distance_2b_hypers)
        call allocate_angle_3b_hypers(n_angle_3b, n_sparse_mpi_angle_3b, angle_3b_hypers)
        call allocate_core_pot_hypers(n_core_pot, n_mpi_core_pot, core_pot_hypers)
     END IF
     !   Now broadcast the data structures
     !    call mpi_bcast(soap_turbo_hypers, size_soap_turbo, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
     !    call mpi_bcast(distance_2b_hypers, size_distance_2b, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
     !    call mpi_bcast(angle_3b_hypers, size_angle_3b, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
     !   VERY IMPORTANT: the broadcasting above only affects the non-allocatable items in the data structures;
     !   we need to manually broadcast allocatable arrays within the structures. To avoid error, we broadcast
     !   also non allocatable variables. It looks ugly as fuck, but
     !   putting this into a subroutine is even worse because we need to get swifty with the communications. So
     !   my solution is to go the ugly way. All that said, I'll probably put this ugly motherfucker in a module.
     !   I should also make some arrays of the correct type and pass all the variables in the stack (of the same
     !   type) at once via broadcasting, to reduce the total number of MPI calls to the minimum. This will be
     !   done at the module's subroutine's level.
     !   soap_turbo allocatable structures
     call cpu_time(time_mpi(1))
     do i = 1, n_soap_turbo
        n_sp = soap_turbo_hypers(i)%n_species
        call mpi_bcast(soap_turbo_hypers(i)%nf(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%rcut_hard(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%rcut_soft(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%rcut_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%atom_sigma_r(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%atom_sigma_t(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%atom_sigma_r_scaling(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%atom_sigma_t_scaling(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%amplitude_scaling(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%central_weight(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%global_scaling(1:n_sp), n_sp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%alpha_max(1:n_sp), n_sp, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%species_types(1:n_sp), 8*n_sp, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        n_sparse = soap_turbo_hypers(i)%n_sparse
        dim = soap_turbo_hypers(i)%dim
        n_nonzero = soap_turbo_hypers(i)%compress_P_nonzero
        call mpi_bcast(soap_turbo_hypers(i)%alphas(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%Qs(1:dim, 1:n_sparse), n_sparse*dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%delta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%zeta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%basis, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%scaling_mode, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%species_types(1:n_sp), 8*n_sp, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%central_species, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%l_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%n_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%radial_enhancement, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%compress_soap, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if( soap_turbo_hypers(i)%compress_soap )then
           cPnz = soap_turbo_hypers(i)%compress_P_nonzero
           call mpi_bcast(soap_turbo_hypers(i)%compress_P_el(1:cPnz), cPnz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call mpi_bcast(soap_turbo_hypers(i)%compress_P_i(1:cPnz), cPnz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
           call mpi_bcast(soap_turbo_hypers(i)%compress_P_j(1:cPnz), cPnz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end if
        call mpi_bcast(soap_turbo_hypers(i)%has_local_properties, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%has_core_electron_be, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (valid_xps .or. params%do_xps_standalone) call mpi_bcast(core_be_lp_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (valid_vdw) call mpi_bcast(vdw_lp_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call mpi_bcast(soap_turbo_hypers(i)%has_vdw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%n_local_properties, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if( soap_turbo_hypers(i)%has_local_properties )then
           do j = 1, soap_turbo_hypers(i)%n_local_properties
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%n_sparse, 1,&
                   & MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%label, 1024,&
                   & MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%delta, 1,&
                   & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%zeta, 1,&
                   & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%V0, 1,&
                   & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%dim, 1, MPI_INTEGER, 0,&
                   & MPI_COMM_WORLD, ierr)
              n_sparse = soap_turbo_hypers(i)&
                   &%local_property_models(j)%n_sparse
              dim = soap_turbo_hypers(i)%local_property_models(j)%dim
              call mpi_bcast(soap_turbo_hypers(i)&
                   &%local_property_models(j)%alphas(1:n_sparse) ,&
                   & n_sparse, MPI_DOUBLE_PRECISION, 0,&
                   & MPI_COMM_WORLD, ierr)
! Fortran runtime warning: An array temporary was created for
              ! argument 'buffer' of procedure 'mpi_bcast'
              call mpi_bcast(soap_turbo_hypers(i) &
                   &%local_property_models(j)%Qs(1:dim, 1:n_sparse),&
                   & n_sparse*dim, MPI_DOUBLE_PRECISION, 0,&
                   & MPI_COMM_WORLD, ierr)
              call mpi_bcast(soap_turbo_hypers(i)%local_property_models(j)%do_derivatives, &
                      & 1, MPI_LOGICAL, 0,&
                      & MPI_COMM_WORLD, ierr)

              call mpi_bcast(soap_turbo_hypers(i)%local_property_models(j)%compute, &
                      & 1, MPI_LOGICAL, 0,&
                      & MPI_COMM_WORLD, ierr)

           end do
        end if
     end do
     do i = 1, n_distance_2b
        n_sparse = distance_2b_hypers(i)%n_sparse
        call mpi_bcast(distance_2b_hypers(i)%alphas(1:n_sparse),&
             & n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(distance_2b_hypers(i)%cutoff(1:n_sparse),&
             & n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(distance_2b_hypers(i)%Qs(1:n_sparse, 1:1),&
             & n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(distance_2b_hypers(i)%delta, 1,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(distance_2b_hypers(i)%sigma, 1,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(distance_2b_hypers(i)%rcut, 1,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(distance_2b_hypers(i)%buffer, 1,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(distance_2b_hypers(i)%species1, 8,&
             & MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(distance_2b_hypers(i)%species2, 8,&
             & MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     end do
     do i = 1, n_angle_3b
        n_sparse = angle_3b_hypers(i)%n_sparse
        call mpi_bcast(angle_3b_hypers(i)%alphas(1:n_sparse),&
             & n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(angle_3b_hypers(i)%cutoff(1:n_sparse),&
             & n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(angle_3b_hypers(i)%Qs(1:n_sparse, 1:3), 3&
             &*n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,&
             & ierr)
        call mpi_bcast(angle_3b_hypers(i)%delta, 1,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%sigma(1:3), 3,&
             & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%rcut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%buffer, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%species_center, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%species1, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%species2, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(angle_3b_hypers(i)%kernel_type, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     end do
     do i = 1, n_core_pot
        n_sparse = core_pot_hypers(i)%n
        call mpi_bcast(core_pot_hypers(i)%x(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%V(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%dVdx2(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%yp1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%ypn, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%species1, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(core_pot_hypers(i)%species2, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     end do
     call cpu_time(time_mpi(2))
     time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
     !   Clean up
     deallocate( n_species_mpi, n_sparse_mpi_soap_turbo, dim_mpi, compress_soap_mpi, n_sparse_mpi_distance_2b, &
          n_sparse_mpi_angle_3b, n_mpi_core_pot, compress_P_nonzero_mpi, n_local_properties_mpi, has_local_properties_mpi )
     if (allocated(local_properties_dim_mpi_soap_turbo)) deallocate(&
          & local_properties_dim_mpi_soap_turbo,&
          & local_properties_n_sparse_mpi_soap_turbo )
#endif
  else
#ifdef _MPIF90
     IF( rank == 0 )THEN
#endif
        write(*,*)'                                       |'
        write(*,*)'ERROR: you must provide a "pot_file"   |  <-- ERROR'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
#ifdef _MPIF90
     END IF
     call mpi_finalize(ierr)
#endif
     stop
  end if
  call cpu_time(time_read_input(2))
  time_read_input(3) = time_read_input(3) + time_read_input(2) - time_read_input(1)
  !**************************************************************************

!! If electronic stopping based on eph model is to be calculated, these data structures are required to be
!! initialized first.
  if ( params%nonadiabatic_processes ) then
	if ( params%eph_Tinfile /= "NULL" ) then
		call ephfdm%EPH_FDM_input_file(params%eph_Tinfile,params%eph_md_last_step)
	end if
	if ( params%eph_Tinfile == "NULL" ) then
		call ephfdm%EPH_FDM_input_params(params%eph_md_last_step, params%eph_gsx, &
		params%eph_gsy, params%eph_gsz, params%in_x0, params%in_x1, params%in_y0, params%in_y1, &
		params%in_z0, params%in_z1, params%eph_Ti_e, params%eph_C_e, params%eph_rho_e, &
		params%eph_kappa_e, params%eph_fdm_steps)
	end if
	call ephbeta%beta_parameters(params%eph_betafile,n_species)
	call ephlsc%eph_InitialValues (params%eph_friction_option, params%eph_random_option, params%eph_fdm_option, &
			params%eph_Toutfile, params%eph_freq_Tout, params%eph_freq_mesh_Tout, params%model_eph, &
			params%eph_E_prev_time, params%eph_md_prev_time)	
  end if
!! -------------------------				---- untill here for initializing eph model elec. stopping




  !**************************************************************************
  ! <----------------------------------------------------------------------------------------------- Finish printouts
#ifdef _MPIF90
  IF( rank == 0 )THEN
#endif
     ! Print out chosen options:
     write(*,*)'                                       |'
     write(*,'(1X,A)')'You specified the following options:   |'
     write(*,*)'                                       |'
     write(*,*)'---------------------------------      |'
     if( len(trim(params%atoms_file)) > 20 )then
        write(*,'(1X,A,A20,A)')'Atoms file = ', adjustr(trim(params%atoms_file)), '...   |'
     else
        write(*,'(1X,A,A20,A)')'Atoms file = ', adjustr(trim(params%atoms_file)), '      |'
     end if
     write(*,*)'---------------------------------      |'
     write(i_char, '(I8)') n_species
     write(*,'(1X,A,A8,A)')'No. of species   = ', adjustl(i_char), '            |'
     do i = 1, n_species
        write(i_char, '(I8)') i
        write(*,'(1X,A,A2,A,A8,A)')'  *) Species #', adjustl(i_char), ' =       ', adjustr(params%species_types(i)), '      |'
     end do
     write(*,*)'---------------------------------      |'
     write(*,'(1X,A,F15.4,A)')'rcut_max = ', rcut_max, ' Angst.      |'
     write(*,*)'---------------------------------      |'
     write(*,*)'                                       |'
     write(*,*)'.......................................|'
#ifdef _MPIF90
  END IF
#endif
  !**************************************************************************










  !**************************************************************************
  ! Print progress bar and initialize timers
  time_neigh = 0.d0
  time_gap = 0.d0
  time_soap = 0.d0
  time_2b = 0.d0
  time_3b = 0.d0
  time_core_pot = 0.d0
  time_vdw = 0.d0
  time_read_xyz = 0.d0
  time_pdf = 0.d0
  time_sf = 0.d0
  time_mc = 0.d0
  time_xrd = 0.d0
  time_nd = 0.d0
  time_xps = 0.d0

  xps_idx = params%xps_idx
  md_istep = -1
  mc_istep = -1
  n_xyz = 0
  i_nested = 0
  i_image = 0


  if( params%do_md )then
#ifdef _MPIF90
     IF( rank == 0 )THEN
#endif
        write(*,*)'                                       |'
        write(*,*)'Doing molecular dynamics...            |'
        if( params%print_progress .and. md_istep > 0 )then
           write(*,*)'                                       |'
           write(*,*)'Progress:                              |'
           write(*,*)'                                       |'
           write(*,'(1X,A)',advance='no')'[                                    ] |'
        end if
#ifdef _MPIF90
     END IF
#endif
     update_bar = params%md_nsteps/36
     if( update_bar < 1 )then
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

  do while( repeat_xyz .or. ( params%do_md .and. md_istep < params%md_nsteps) &
       .or. ( params%do_mc .and. mc_istep < params%mc_nsteps))
     exit_loop = .false.

     if( params%do_mc )then
        mc_istep = mc_istep + 1
        ! Undo if the step is md related
        if(md_istep > -1) mc_istep = mc_istep - 1

     end if

     if( params%do_md )then
        md_istep = md_istep + 1
     else
        n_xyz = n_xyz + 1
     end if

     !   Update progress bar
     if( params%print_progress .and. counter == update_bar .and. (.not. params%do_mc) )then
#ifdef _MPIF90
        IF( rank == 0 )THEN
#endif
           do j = 1, 36+3
              write(*,"(A)", advance="no") creturn
           end do
           write (*,"(1X,A)",advance="no") "["
           do i = 1, 36*md_istep/params%md_nsteps
              write (*,"(A)",advance="no") "."
           end do
           do i = 36*md_istep/params%md_nsteps+1, 36
              write (*,"(A)",advance="no") " "
           end do
           write (*,"(A)",advance="no") "] |"
           if( md_istep == params%md_nsteps )then
              write(*,*)
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

     if( (params%do_md .and. md_istep == 0) )then
        call cpu_time(time_read_xyz(1))
#ifdef _MPIF90
        IF( rank == 0 )THEN
#endif
           if(mc_istep > 0)then
              call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                   n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                   positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                   xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                   n_sites, .not. params%mc_write_xyz, fix_atom, params%t_beg, &
                   params%write_array_property(6), .not. params%mc_write_xyz)
              rebuild_neighbors_list = .true.

           else if( .not. params%do_nested_sampling .or. mc_istep == 0 )then
              call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                   n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                   positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                   xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                   n_sites, .false., fix_atom, params%t_beg, &
                   params%write_array_property(6), .false. )

           end if

           ! call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
           !               n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
           !               positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
           !               xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
           !               n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), .true. )
           !     Only rank 0 handles these variables
           !      allocate( positions_prev(1:3, 1:size(positions,2)) )
           !      allocate( positions_diff(1:3, 1:size(positions,2)) )
           if(.not. allocated(forces_prev))allocate( forces_prev(1:3, 1:n_sites) )
           if(.not. allocated(positions_prev))allocate( positions_prev(1:3, 1:n_sites) )
           if(.not. allocated(positions_diff))allocate( positions_diff(1:3, 1:n_sites) )
           positions_diff = 0.d0
           rebuild_neighbors_list = .true.
#ifdef _MPIF90
        END IF
#endif
        call cpu_time(time_read_xyz(2))
        time_read_xyz(3) = time_read_xyz(3) + time_read_xyz(2) - time_read_xyz(1)
        !     If we're doing MD, we don't read beyond the first snapshot in the XYZ file
        repeat_xyz = .false.
        !     At the moment, we can't do prediction if the unit cell doesn't fit a whole cutoff sphere
#ifdef _MPIF90
        IF( rank == 0 )THEN
#endif
           !     CLEAN THIS UP <------------------------------------------------------------------- LOOK HERE
           !      if( size(positions,2) /= n_sites )then
           if( .false. )then
              write(*,*) "Sorry, at the moment TurboGAP can't do MD for unit cells smaller than ", &
                   "a cutoff sphere <-- ERROR"
#ifdef _MPIF90
              call mpi_finalize(ierr)
#endif
              stop
           end if
#ifdef _MPIF90
        END IF
#endif
     else if( .not. params%do_md )then
        call cpu_time(time_read_xyz(1))
#ifdef _MPIF90
        IF( rank == 0 )THEN
#endif
           if(mc_istep > 0)then
              call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                   n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                   positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                   xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                   n_sites, .not. params%mc_write_xyz, fix_atom, params%t_beg, &
                   params%write_array_property(6), .not. params%mc_write_xyz )
              rebuild_neighbors_list = .true.
           else
              call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                   n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                   positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                   xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                   n_sites, .false., fix_atom, params%t_beg, params%write_array_property(6), &
                   .false. )
           end if
#ifdef _MPIF90
        END IF
#endif
        call cpu_time(time_read_xyz(2))
        time_read_xyz(3) = time_read_xyz(3) + time_read_xyz(2) - time_read_xyz(1)
#ifdef _MPIF90
        call cpu_time(time_mpi(1))
        call mpi_bcast(repeat_xyz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call cpu_time(time_mpi(2))
        time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
#endif
        rebuild_neighbors_list = .true.
     end if
     !   Broadcast the info in the XYZ file: positions, velocities, masses, xyz_species, xyz_species_supercell,
     !   species, species_supercell, indices, a_box, b_box, c_box and n_sites. I should put this into a module!!!!!!!




#ifdef _MPIF90
     IF( rank == 0 )THEN
        n_pos = size(positions,2)
        n_sp = size(xyz_species,1)
        n_sp_sc = size(xyz_species_supercell,1)
        if ( params%do_mc .and. (mc_move /= "md" .or. md_istep == 0) .and. params%mc_hamiltonian )then
           if(mc_istep > 0) E_kinetic_prev = E_kinetic
           call random_number( velocities )
           call remove_cm_vel(velocities(1:3,1:n_sites), masses(1:n_sites))
           E_kinetic = 0.d0
           do i = 1, n_sites
              E_kinetic = E_kinetic + 0.5d0 * masses(i) * dot_product(velocities(1:3, i), velocities(1:3, i))
           end do
           instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic
           velocities = velocities * dsqrt(params%t_beg/instant_temp)
           if (mc_istep > 0)then
              E_kinetic = E_kinetic_prev
              instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic
              ! Reversing as we want it to be at the instant temp and not at t_beg
              velocities = velocities * dsqrt(instant_temp / params%t_beg)
           else
              E_kinetic = E_kinetic * params%t_beg/instant_temp
           end if
        end if

     END IF
     call cpu_time(time_mpi(1))
     call mpi_bcast(n_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sp_sc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call cpu_time(time_mpi(2))
     time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)

     IF( rank /= 0 )THEN
        if(allocated(positions))deallocate(positions)
        allocate( positions(1:3, n_pos) )
        if( params%do_md .or. params%do_nested_sampling .or. params%do_mc )then
           if(allocated(velocities))deallocate(velocities)
           allocate( velocities(1:3, n_pos) )
           !      allocate( masses(n_pos) )
           if(allocated( masses ))deallocate( masses )
           allocate( masses(1:n_sp) )
           ! if(allocated( fix_atom ))deallocate( fix_atom )
           ! allocate( fix_atom(1:3, 1:n_sp) )
        end if
        if(allocated( xyz_species ))deallocate( xyz_species )
        allocate( xyz_species(1:n_sp) )
        if(allocated( species ))deallocate( species )
        allocate( species(1:n_sp) )
        if(allocated( xyz_species_supercell ))deallocate( xyz_species_supercell )
        allocate( xyz_species_supercell(1:n_sp_sc) )
        if(allocated( species_supercell ))deallocate( species_supercell )
        allocate( species_supercell(1:n_sp_sc) )
        if(allocated( fix_atom ))deallocate( fix_atom )
        allocate( fix_atom(1:3,1:n_sp) )

     END IF
     call cpu_time(time_mpi_positions(1))
     call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     if( params%do_md .or. params%do_nested_sampling .or. params%do_mc .or. params%mc_hamiltonian)then
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
     call cpu_time(time_mpi_positions(2))
     time_mpi_positions(3) = time_mpi_positions(3) + time_mpi_positions(2) - time_mpi_positions(1)
#endif
     !   Now that all ranks know the size of n_sites, we allocate do_list
     if( .not. params%do_md .or. (params%do_md .and. md_istep == 0) .or. &
          ( params%do_mc ) )then
        if( allocated(do_list))deallocate(do_list)
        allocate( do_list(1:n_sites) )
        do_list = .true.
     end if
     !
     call cpu_time(time1)
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
     if( rebuild_neighbors_list .and.  params%do_mc .and. mc_istep > 0 )then
        call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
             n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
             positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
             xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
             n_sites, .true., fix_atom, params%t_beg, &
             params%write_array_property(6), .false.)
     else if( rebuild_neighbors_list )then
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
     if( rank < mod( n_sites, ntasks ) )then
        i_beg = 1 + rank*(n_sites / ntasks + 1)
     else
        i_beg = 1 + mod(n_sites, ntasks)*(n_sites / ntasks + 1) + (rank - mod( n_sites, ntasks))*(n_sites / ntasks)
     end if
     if( rank < mod( n_sites, ntasks ) )then
        i_end = (rank+1)*(n_sites / ntasks + 1)
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
! #ifdef _MPIF90
!         if(allocated(n_neigh_local))deallocate( n_neigh_local )
! #endif
!      end if

     call build_neighbors_list(positions, a_box, b_box, c_box, params%do_timing, &
          species_supercell, rcut_max, n_atom_pairs, rjs, &
          thetas, phis, xyz, n_neigh_local, neighbors_list, neighbor_species, n_sites, indices, &
          rebuild_neighbors_list, do_list, rank )
     if( rebuild_neighbors_list )then
        !     Get total number of atom pairs
        call mpi_allgather(n_atom_pairs, 1, MPI_INTEGER, n_atom_pairs_by_rank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        n_atom_pairs_total = sum(n_atom_pairs_by_rank)
        n_atom_pairs = n_atom_pairs_total

        !     Get number of neighbors
        if( .not. allocated(n_neigh))allocate( n_neigh(1:n_sites) )
        call mpi_reduce(n_neigh_local, n_neigh, n_sites, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(n_neigh, n_sites, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        j_beg = 1
        j_end = n_atom_pairs_by_rank(rank+1)
     end if
#else
     call build_neighbors_list(positions, a_box, b_box, c_box, params%do_timing, &
          species_supercell, rcut_max, n_atom_pairs, rjs, &
          thetas, phis, xyz, n_neigh, neighbors_list, neighbor_species, n_sites, indices, &
          rebuild_neighbors_list, do_list, rank )
     i_beg = 1
     i_end = n_sites
     n_sites_mpi = n_sites
     j_beg = 1
     j_end = n_atom_pairs
     n_atom_pairs_by_rank(rank+1) = n_atom_pairs
#endif
     !   Compute the volume of the "primitive" unit cell
     v_uc = dot_product( cross_product(a_box, b_box), c_box ) / (dfloat(indices(1)*indices(2)*indices(3)))
     call cpu_time(time2)
     time_neigh = time_neigh + time2 - time1
     !**************************************************************************




     !**************************************************************************
     !   If we are doing prediction, we run this chunk of code
     if( params%do_prediction .or. params%write_soap .or. params%write_derivatives)then
        call cpu_time(time1)

        !     We only need to reallocate the arrays if the number of sites changes
        ! REMOVE TRUE FROM IF STATEMENT
        if( n_sites /= n_sites_prev .or. params%do_mc  )then
           if( allocated(energies) )deallocate( energies, energies_soap, energies_2b, energies_3b, energies_core_pot, &
                this_energies, energies_vdw, this_forces, energies_lp, energies_exp  )
           allocate( energies(1:n_sites) )
           allocate( this_energies(1:n_sites) )
           allocate( energies_soap(1:n_sites) )
           allocate( energies_2b(1:n_sites) )
           allocate( energies_3b(1:n_sites) )
           allocate( energies_core_pot(1:n_sites) )
           allocate( energies_vdw(1:n_sites) )
           allocate( energies_lp(1:n_sites) )
           allocate( energies_exp(1:n_sites) )

           if (params%do_pair_distribution .and. params%valid_pdf)then
              if ( allocated( energies_pdf ) ) deallocate(energies_pdf)
              allocate( energies_pdf(1:n_sites) )
           end if

           if (params%do_structure_factor .and. params%valid_sf)then
              if ( allocated( energies_sf ) ) deallocate(energies_sf)
              allocate( energies_sf(1:n_sites) )
           end if

           if (params%do_xrd .and. params%valid_xrd)then
              if ( allocated( energies_xrd ) ) deallocate(energies_xrd)
              allocate( energies_xrd(1:n_sites) )
           end if

           if (params%do_nd .and. params%valid_nd)then
              if ( allocated( energies_nd ) ) deallocate(energies_nd)
              allocate( energies_nd(1:n_sites) )
           end if

           !       This needs to be allocated even if no force prediction is needed:
           allocate( this_forces(1:3, 1:n_sites) )
        end if
        energies = 0.d0
        energies_soap = 0.d0
        energies_2b = 0.d0
        energies_3b = 0.d0
        energies_core_pot = 0.d0
        energies_vdw = 0.d0
        energies_lp = 0.d0
        energies_exp = 0.d0

        if (params%do_pair_distribution .and. params%valid_pdf) energies_pdf = 0.d0
        if (params%do_structure_factor  .and. params%valid_sf)  energies_sf = 0.d0
        if (params%do_xrd .and. params%valid_xrd)               energies_xrd = 0.d0
        if (params%do_nd  .and. params%valid_nd)                energies_nd = 0.d0


        ! Adding allocation of local properties

        ! Now one could use pointers such that hirshfeld_v(:) acts as an alias for local_properties(vdw_index,:)...

        if( any( soap_turbo_hypers(:)%has_local_properties ) )then
           if( n_sites /= n_sites_prev .or.  params%do_mc  )then
              if( allocated(local_properties) )then
                 nullify( this_local_properties_pt )
                 deallocate( this_local_properties, local_properties )
                 if( params%do_forces )then
                    nullify( this_local_properties_cart_der_pt )
                    deallocate( this_local_properties_cart_der, local_properties_cart_der )
                 end if
              end if
              allocate( local_properties(1:n_sites, 1:params%n_local_properties) )
              allocate( this_local_properties(1:n_sites, 1:params%n_local_properties) )
              this_local_properties_pt => this_local_properties

              !         I don't remember why this needs a pointer <----------------------------------------- CHECK

           end if
           local_properties = 0.d0


           if( params%do_forces )then
              if( n_atom_pairs_by_rank(rank+1) /= n_atom_pairs_by_rank_prev  )then
                 if( allocated(local_properties_cart_der) )deallocate( local_properties_cart_der, this_local_properties_cart_der )
                 allocate( local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1), 1:params%n_local_properties) )
                 allocate( this_local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1), 1:params%n_local_properties) )
              end if
              if(.not. allocated(local_properties_cart_der) )then
                 allocate( local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1), 1:params%n_local_properties) )
                 allocate( this_local_properties_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1), 1:params%n_local_properties) )
              end if

              local_properties_cart_der = 0.d0
              this_local_properties_cart_der_pt =>&
                   & this_local_properties_cart_der(1:3,&
                   & 1:n_atom_pairs_by_rank(rank+1), 1:params&
                   &%n_local_properties)
           end if
        end if

        ! Now go through the soap turbo hypers, and see if any are vdw or
        ! otherwise, if vdw, one can have pointers to point to the data
        ! structures such that it makes things clearer. One needs to check
        ! that this allocation still works iwth if(allocated(hirsh_v))
        ! statements



        if( params%do_forces )then
           if( n_sites /= n_sites_prev .or.  params%do_mc )then
              if( allocated(forces) )deallocate( forces, forces_soap, forces_2b, forces_3b, forces_core_pot, forces_vdw,&
                   & forces_lp )
              allocate( forces(1:3, 1:n_sites) )
              allocate( forces_soap(1:3, 1:n_sites) )
              allocate( forces_2b(1:3, 1:n_sites) )
              allocate( forces_3b(1:3, 1:n_sites) )
              allocate( forces_core_pot(1:3, 1:n_sites) )
              allocate( forces_vdw(1:3,1:n_sites) )
              allocate( forces_lp(1:3,1:n_sites) )

              if (params%do_pair_distribution .and. params%exp_forces .and. params%valid_pdf)then
                 if ( allocated( forces_pdf ) ) deallocate(forces_pdf)
                 allocate( forces_pdf(1:3,1:n_sites) )
              end if

              if (params%do_structure_factor .and. params%exp_forces .and. params%valid_sf)then
                 if ( allocated( forces_sf ) ) deallocate(forces_sf)
                 allocate( forces_sf(1:3,1:n_sites) )
              end if

              if (params%do_xrd .and. params%exp_forces .and. params%valid_xrd)then
                 if ( allocated( forces_xrd ) ) deallocate(forces_xrd)
                 allocate( forces_xrd(1:3,1:n_sites) )
              end if

              if (params%do_nd .and. params%exp_forces .and. params%valid_nd)then
                 if ( allocated( forces_nd ) ) deallocate(forces_nd)
                 allocate( forces_nd(1:3,1:n_sites) )
              end if


           end if
           forces = 0.d0
           forces_soap = 0.d0
           forces_2b = 0.d0
           forces_3b = 0.d0
           forces_core_pot = 0.d0
           forces_vdw = 0.d0
           forces_lp = 0.d0
           virial = 0.d0
           virial_soap = 0.d0
           virial_2b = 0.d0
           virial_3b = 0.d0
           virial_core_pot = 0.d0
           virial_vdw = 0.d0
           virial_lp = 0.d0

           if (params%do_pair_distribution .and. params%exp_forces .and. params%valid_pdf)then
              forces_pdf = 0.d0
              virial_pdf = 0.d0
           end if

           if (params%do_structure_factor .and. params%exp_forces .and. params%valid_sf)then
              forces_sf = 0.d0
              virial_sf = 0.d0
           end if

           if (params%do_xrd .and. params%exp_forces .and. params%valid_xrd)then
              forces_xrd = 0.d0
              virial_xrd = 0.d0
           end if

           if (params%do_nd .and. params%exp_forces .and. params%valid_nd)then
              forces_nd = 0.d0
              virial_nd = 0.d0
           end if


        end if

        if( params%do_prediction )then
           !       Assign the e0 to each atom according to its species
           !        do i = 1, n_sites
           do i = i_beg, i_end
              do j = 1, n_species
                 if( xyz_species(i) ==  params%species_types(j) )then
                    energies(i) = params%e0(j)
                 end if
              end do
           end do
        end if
        !     Collect all energies
#ifdef _MPIF90
        call cpu_time(time_mpi_ef(1))
        call mpi_reduce(energies, this_energies, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call cpu_time(time_mpi_ef(2))
        time_mpi_ef(3) = time_mpi_ef(3) + time_mpi_ef(2) - time_mpi_ef(1)
        energies = this_energies
#endif

        !     Loop through soap_turbo descriptors - we always call this routine, even if we don't want to do prediction
        n_lp_count = 0 ! This counts the local properties
        do i = 1, n_soap_turbo
           call cpu_time(time_soap(1))
           !       Compute number of pairs for this SOAP. SOAP has in general a different cutoff than overall max
           !       cutoff, so the number of pairs may be a lot smaller for the SOAP subset.
           !       This subroutine splits the load optimally so as to not use more memory per MPI process than available.
           !       TurboGAP does not check how much memory is available, it just relies on heuristics and a user provided
           !       max_Gbytes_per_process (default = 1.d0)
           call get_number_of_atom_pairs( n_neigh(i_beg:i_end), rjs(j_beg:j_end), soap_turbo_hypers(i)%rcut_max, &
                soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%n_max, &
                params%max_Gbytes_per_process, i_beg_list, i_end_list, j_beg_list, j_end_list )
           do j = 1, size(i_beg_list)
              this_i_beg = i_beg - 1 + i_beg_list(j)
              this_i_end = i_beg - 1 + i_end_list(j)
              this_j_beg = j_beg - 1 + j_beg_list(j)
              this_j_end = j_beg - 1 + j_end_list(j)
              this_n_sites_mpi = this_i_end - this_i_beg + 1
              this_energies = 0.d0
              if( params%do_forces )then
                 this_forces = 0.d0
                 this_virial = 0.d0
              end if
              if( soap_turbo_hypers(i)%has_local_properties )then
                 this_local_properties = 0.d0
                 if( params%do_forces )then
                    this_local_properties_cart_der = 0.d0
                    !             I don't remember why this needs a pointer <----------------------------------------- CHECK
                    nullify(this_local_properties_cart_der_pt)
                    this_local_properties_cart_der_pt =>&
                         & this_local_properties_cart_der(1:3,&
                         & this_j_beg:this_j_end, 1:params&
                         &%n_local_properties)
                 end if
              end if


              ! call get_gap_soap(n_sites, this_n_sites_mpi, n_neigh(this_i_beg:this_i_end), neighbors_list(this_j_beg:this_j_end), &
              !      soap_turbo_hypers(i)%n_species, soap_turbo_hypers(i)%species_types, &
              !      rjs(this_j_beg:this_j_end), thetas(this_j_beg:this_j_end), phis(this_j_beg:this_j_end), &
              !      xyz(1:3, this_j_beg:this_j_end), &
              !      soap_turbo_hypers(i)%alpha_max, &
              !      soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%dim, soap_turbo_hypers(i)%rcut_hard, &
              !      soap_turbo_hypers(i)%rcut_soft, soap_turbo_hypers(i)%nf, soap_turbo_hypers(i)%global_scaling, &
              !      soap_turbo_hypers(i)%atom_sigma_r, soap_turbo_hypers(i)%atom_sigma_r_scaling, &
              !      soap_turbo_hypers(i)%atom_sigma_t, soap_turbo_hypers(i)%atom_sigma_t_scaling, &
              !      soap_turbo_hypers(i)%amplitude_scaling, soap_turbo_hypers(i)%radial_enhancement, &
              !      soap_turbo_hypers(i)%central_weight, soap_turbo_hypers(i)%basis, &
              !      soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
              !      params%do_prediction, params%write_soap, params%write_derivatives, &
              !      soap_turbo_hypers(i)%compress_soap, soap_turbo_hypers(i)%compress_P_nonzero, &
              !      soap_turbo_hypers(i)%compress_P_i, soap_turbo_hypers(i)%compress_P_j, &
              !      soap_turbo_hypers(i)%compress_P_el, &
              !      soap_turbo_hypers(i)%delta, soap_turbo_hypers(i)%zeta, soap_turbo_hypers(i)%central_species, &
              !      xyz_species(this_i_beg:this_i_end), xyz_species_supercell, soap_turbo_hypers(i)%alphas, &
              !      soap_turbo_hypers(i)%Qs, params%all_atoms, params%which_atom, indices, soap, soap_cart_der, &
              !      der_neighbors, der_neighbors_list, &
              !      this_energies, this_forces, &
              !      this_virial, soap_turbo_hypers(i)&
              !      &%has_local_properties, n_pairs, in_to_out_pairs, n_all_sites, in_to_out_site, n_neigh_out, n_sites_out )

              call get_gap_soap(n_sites, this_n_sites_mpi, n_neigh(this_i_beg:this_i_end), neighbors_list(this_j_beg:this_j_end), &
                   soap_turbo_hypers(i)%n_species, soap_turbo_hypers(i)%species_types, &
                   rjs(this_j_beg:this_j_end), thetas(this_j_beg:this_j_end), phis(this_j_beg:this_j_end), &
                   xyz(1:3, this_j_beg:this_j_end), &
                   soap_turbo_hypers(i)%alpha_max, &
                   soap_turbo_hypers(i)%l_max, soap_turbo_hypers(i)%dim, soap_turbo_hypers(i)%rcut_hard, &
                   soap_turbo_hypers(i)%rcut_soft, soap_turbo_hypers(i)%nf, soap_turbo_hypers(i)%global_scaling, &
                   soap_turbo_hypers(i)%atom_sigma_r, soap_turbo_hypers(i)%atom_sigma_r_scaling, &
                   soap_turbo_hypers(i)%atom_sigma_t, soap_turbo_hypers(i)%atom_sigma_t_scaling, &
                   soap_turbo_hypers(i)%amplitude_scaling, soap_turbo_hypers(i)%radial_enhancement, &
                   soap_turbo_hypers(i)%central_weight, soap_turbo_hypers(i)%basis, &
                   soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
                   params%do_prediction, params%write_soap, params%write_derivatives, &
                   soap_turbo_hypers(i)%compress_soap, soap_turbo_hypers(i)%compress_P_nonzero, &
                   soap_turbo_hypers(i)%compress_P_i, soap_turbo_hypers(i)%compress_P_j, &
                   soap_turbo_hypers(i)%compress_P_el, &
                   soap_turbo_hypers(i)%delta, soap_turbo_hypers(i)%zeta, soap_turbo_hypers(i)%central_species, &
                   xyz_species(this_i_beg:this_i_end), xyz_species_supercell, soap_turbo_hypers(i)%alphas, &
                   soap_turbo_hypers(i)%Qs, params%all_atoms, params%which_atom, indices, soap, soap_cart_der, &
                   der_neighbors, der_neighbors_list, &
                   & soap_turbo_hypers(i)%has_local_properties,&
                   & soap_turbo_hypers(i)%n_local_properties,&
                   & soap_turbo_hypers(i)%local_property_models,&
                   & this_energies, this_forces, this_local_properties_pt,&
                   & this_local_properties_cart_der_pt,&
                   & local_property_indexes, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                   & this_virial, n_lp_count )


              ! We can have a pointer to specific parts of this_local_properties array to then


              ! if (soap_turbo_hypers(i)%has_local_properties)then

              !    ! only iterating over the computed properties
              !    this_local_properties = 0.d0
              !    if( params%do_forces )then
              !       this_local_properties_cart_der = 0.d0
              !    end if

              !    do l = 1, soap_turbo_hypers(i)%n_local_properties
              !       ! We can increment the counter for the local properties
              !       i2 = i2 + 1
              !       if (soap_turbo_hypers(i)%local_property_models(l)%compute)then
              !          ! compute the local property!
              !          ! Above the local properties passes are this_local_properties so one must change
              !          !
              !          ! Allocate the pointer

              !          this_local_properties_pt =>&
              !               & this_local_properties(1:n_sites,&
              !               & local_property_indexes(i2))
              !             !             I don't remember why this
              !          !             needs a pointer
              !          !
              !          !<-----------------------------------------
              !          !CHECK
              !          this_local_properties_cart_der_pt =>&
              !               & this_local_properties_cart_der(1:3,&
              !               & this_j_beg:this_j_end,&
              !               & local_property_indexes(i2))

              !          call get_local_properties( soap, &
              !               soap_turbo_hypers(i)%local_property_models(l)%Qs, &
              !               soap_turbo_hypers(i)%local_property_models(l)%alphas, &
              !               soap_turbo_hypers(i)%local_property_models(l)%V0, &
              !               soap_turbo_hypers(i)%local_property_models(l)%delta, &
              !               soap_turbo_hypers(i)%local_property_models(l)%zeta, &
              !               this_local_properties_pt, &
              !               soap_turbo_hypers(i)%local_property_models(l)%do_derivatives, &
              !               soap_cart_der, &
              !               n_neigh_out, &
              !               this_local_properties_cart_der_pt, n_pairs,&
              !               & in_to_out_pairs, n_all_sites,&
              !               & in_to_out_site,  n_sites_out )


              !          nullify(this_local_properties_pt)
              !          nullify(this_local_properties_cart_der_pt)
              !       end if
              !    end do
              !    ! Now deallocate the arrays which were not deallocated in get_gap_soap
              !    deallocate( in_to_out_site, in_to_out_pairs, n_neigh_out )

              !    if( .not. params%write_soap  ) deallocate( soap )

              !    if( params%do_derivatives .and. .not. params%write_derivatives)then
              !       deallocate( soap_cart_der )
              !    end if
              ! end if



              energies_soap = energies_soap + this_energies

              if( soap_turbo_hypers(i)%has_local_properties )then

                 local_properties(:,:) = local_properties(:,:) + this_local_properties(:,:)
                 if(  any(soap_turbo_hypers(i)&
                      &%local_property_models(:)%do_derivatives) &
                      & .and. params%do_derivatives)then
                    local_properties_cart_der(:,:,:) =&
                         & local_properties_cart_der(:,:,:) +&
                         & this_local_properties_cart_der(:,:,:)
                 end if


                 ! if( soap_turbo_hypers(i)%has_vdw )then
                 !    !                    hirshfeld_v => local_properties( :, vdw_lp_index)
                 !    if (any(soap_turbo_hypers(i)&
                 !         &%local_property_models(:)%do_derivatives)&
                 !         & .and. params%do_derivatives)&
                 !         & hirshfeld_v_cart_der =>&
                 !         & local_properties_cart_der( :, :,&
                 !         & vdw_lp_index)
                 ! end if

              end if
              if( params%do_forces )then
                 forces_soap = forces_soap + this_forces
                 virial_soap = virial_soap + this_virial
              end if
           end do
           n_lp_count = n_lp_count + soap_turbo_hypers(i)%n_local_properties

           deallocate( i_beg_list, i_end_list, j_beg_list, j_end_list )


           ! THIS WON'T WORK! THE SOAP AND SOAP DERIVATIVES NEED TO BE COLLECTED FROM ALL RANKS <--------------------- FIX THIS!!!!
           ! AT THE MOMENT I'M MAKING THE CODE PRINT AN ERROR MESSAGE AND STOP EXECUTION IF THE USER TRIES TO WRITE OUT THESE
           ! FILES WITH MORE THAN ONE MPI TASK
#ifdef _MPIF90
           IF( rank == 0 )THEN
#endif
              !       Write out stuff - THIS SHOULD PROBABLY BE PUT IN A MODULE
              if( n_soap_turbo == 1 )then
                 i_char = ""
              else
                 write(i_char, '(I7)') i
                 i_char = "_" // adjustl(i_char)
              end if
              !       Write the SOAP vectors - NOT THE OPTIMAL STRATEGY IN TERMS OF DISK SPACE SINCE SOME ATOMS HAVE SOAP = 0
              if( params%write_soap )then
                 if( n_xyz == 1 .or. md_istep == 0 )then
                    open(unit=10, file="soap" // trim(i_char) // ".dat", status="unknown")
                 else
                    open(unit=10, file="soap" // trim(i_char) // ".dat", status="old", position="append")
                 end if
                 if( .not. params%do_md .or. &
                      (params%do_md .and. (md_istep == 0 .or. md_istep == params%md_nsteps .or. &
                      modulo(md_istep, params%write_xyz) == 0)) )then
                    n_sites_this = size(soap,2)
                    n_soap = size(soap,1)
                    write(10, *) n_sites_this, n_soap
                    do i2 = 1, n_sites_this
                       write(10, '(*(ES24.15))') soap(1:n_soap, i2)
                    end do
                 end if
                 close(10)
              end if
              if(allocated(soap)) deallocate(soap)

              !       Optionally, write out the derivatives (might take a lot of disk space)
              if( ( params%do_derivatives .or. params%do_derivatives_fd ) .and. params%write_derivatives )then
                 if( n_xyz == 1 .or. md_istep == 0 )then
                    open(unit=10, file="soap_der" // trim(i_char) // ".dat", status="unknown")
                 else
                    open(unit=10, file="soap_der" // trim(i_char) // ".dat", status="old", position="append")
                 end if
                 if( .not. params%do_md .or. &
                      (params%do_md .and. (md_istep == 0 .or. md_istep == params%md_nsteps .or. &
                      modulo(md_istep, params%write_xyz) == 0)) )then
                    !           Note, this n_sites is not the same as the total number of sites, it's just the total number
                    !           of sites that have a derivative, since the first neighbor of each site is itself, the site
                    !           ID can always be retrieved from there. Note also that the sites are not necessarily given in
                    !           order
                    n_sites_this = size(der_neighbors,1)
                    n_soap = size(soap_cart_der, 2)
                    n_atom_pairs = size(der_neighbors_list,1)
                    write(10, *) n_sites, n_soap, n_atom_pairs
                    k = 1
                    k2 = 0
                    do i2 = 1, n_sites_this
                       write(10,*) der_neighbors_list(k), der_neighbors(i2), der_neighbors_list(k:k+der_neighbors(i2)-1)
                       k = k + der_neighbors(i)
                       do j = 1, der_neighbors(i)
                          k2 = k2 + 1
                          write(10, '(*(ES24.15))') soap_cart_der(1, 1:n_soap, k2)
                          write(10, '(*(ES24.15))') soap_cart_der(2, 1:n_soap, k2)
                          write(10, '(*(ES24.15))') soap_cart_der(3, 1:n_soap, k2)
                       end do
                    end do
                 end if
                 close(10)
              end if
              if( params%write_derivatives )then
                 deallocate( soap_cart_der, der_neighbors, der_neighbors_list )
              end if
#ifdef _MPIF90
           END IF
#endif

           call cpu_time(time_soap(2))
           time_soap(3) = time_soap(3) + time_soap(2) - time_soap(1)

           call cpu_time(time2)
           time_gap = time_gap + time2 - time1

        end do


#ifdef _MPIF90
        if( any( soap_turbo_hypers(:)%has_local_properties) )then
           call cpu_time(time_mpi(1))
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

           call cpu_time(time_mpi(2))
           time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
        end if
#endif

        !     Compute vdW energies and forces
        if( any( soap_turbo_hypers(:)%has_vdw ) .and.( params%do_prediction ) &
             .and. params%vdw_type == "ts" )then
           call cpu_time(time_vdw(1))
#ifdef _MPIF90
           allocate( this_energies_vdw(1:n_sites) )
           this_energies_vdw = 0.d0
           if( params%do_forces )then
              allocate( this_forces_vdw(1:3,1:n_sites) )
              this_forces_vdw = 0.d0
              this_virial_vdw = 0.d0
           end if
#endif
           allocate(v_neigh_vdw(1:j_end-j_beg+1))
           v_neigh_vdw = 0.d0
           k = 0
           do i = i_beg, i_end
              do j = 1, n_neigh(i)
                 !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
                 j2 = mod(neighbors_list(j_beg + k)-1, n_sites) + 1
                 k = k + 1
                 !                 v_neigh_vdw(k) = hirshfeld_v(j2)
                 v_neigh_vdw(k) = local_properties(j2, vdw_lp_index)
              end do
           end do

           call get_ts_energy_and_forces( local_properties(i_beg:i_end, vdw_lp_index), &
                & local_properties_cart_der(1:3, j_beg:j_end, vdw_lp_index), &
                n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                neighbor_species(j_beg:j_end), &
                params%vdw_rcut, params%vdw_buffer, &
                params%vdw_rcut_inner, params%vdw_buffer_inner, &
                rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
                params%vdw_alpha0_ref, params%do_forces, &
#ifdef _MPIF90
                this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw)
#else
           energies_vdw(i_beg:i_end), forces_vdw, virial_vdw )
#endif
           call cpu_time(time_vdw(2))
           time_vdw(3) = time_vdw(2) - time_vdw(1)

           deallocate(v_neigh_vdw)
        end if


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

        if ( params%do_exp )then
           do i = 1, params%n_exp
              ! If we want to compute the experimental interpolation, we do it now.

              call get_write_condition( params%do_mc, params%do_md&
                   &, mc_istep, md_istep, params%write_xyz,&
                   & write_condition)


              if ( params%exp_data(i)%compute_exp )then
                 if (allocated(params%exp_data(i)%x)) deallocate(params%exp_data(i)%x)
                 if (allocated(params%exp_data(i)%y)) deallocate(params%exp_data(i)%y)
                 call calculate_exp_interpolation(params%exp_data(i)&
                      &%x, params%exp_data(i)%y, params%exp_data(i)&
                      &%n_samples, params%exp_data(i)%data)


                 call preprocess_exp_data(params, params%exp_data(i)%x,&
                      & params%exp_data(i)%y, params%exp_data(i)%label,&
                      & n_sites, dot_product( cross_product(a_box,&
                      & b_box), c_box ) / (dfloat(indices(1)*indices(2) &
                      &*indices(3)) ), params%exp_data(i)%input, exp_output, .true. )

                 if (params%write_exp .and. .not.  params&
                      &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition ) then

                    call get_overwrite_condition( params%do_mc,&
                         & params%do_md, mc_istep, md_istep, params&
                         &%write_xyz, overwrite_condition)

                    call write_exp_data(params%exp_data(i)%x, params&
                         &%exp_data(i)%y, overwrite_condition,&
                         & trim(params%exp_data(i) %label) //&
                         & "_exp.dat", params%exp_data(i) %label  )
                 end if


              end if

              if ( params%exp_data(i)%compute_exp .and. .not.  params&
                   &%exp_data(i)%wrote_exp .and. rank == 0  .and. write_condition) then

                 if (params%write_exp) then
                    write(filename,'(A)')&
                         & trim(params%exp_data(i)%label) // "_exp_fit.dat"

                    call get_overwrite_condition( params%do_mc,&
                         & params%do_md, mc_istep, md_istep, params&
                         &%write_xyz, overwrite_condition)

                    call write_exp_data(params%exp_data(i)%x, params%exp_data(i)%y,&
                         & overwrite_condition, trim(params&
                         &%exp_data(i)%label) // "_exp_fit.dat",&
                         & trim(params%exp_data(i)%label) // " : output = "&
                         & // trim( exp_output ))

                 end if

              end if

              params%exp_data(i)%wrote_exp = .true.
              params%exp_data(i)%compute_exp = .true.

           end do
        end if


        !###################################################!
        !###---   XPS Forces and Spectra Prediction   ---###!
        !###################################################!

        !     Compute core_electron_be energies and forces
        if( any( soap_turbo_hypers(:)%has_core_electron_be ) .and.( params%do_prediction ) &
             .and. valid_xps )then
           call cpu_time(time_xps(1))

#ifdef _MPIF90
           allocate( this_energies_lp(1:n_sites) )
           this_energies_lp = 0.d0
           if( params%do_forces )then
              allocate( this_forces_lp(1:3,1:n_sites) )
              this_forces_lp = 0.d0
              this_virial_lp = 0.d0
           end if
#endif
           allocate(v_neigh_lp(1:j_end-j_beg+1))
           v_neigh_lp = 0.d0
           k = 0
           do i = i_beg, i_end
              do j = 1, n_neigh(i)
                 !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
                 j2 = mod(neighbors_list(j_beg + k)-1, n_sites) + 1
                 k = k + 1
                 !                 v_neigh_lp(k) = hirshfeld_v(j2)
                 v_neigh_lp(k) = local_properties(j2, core_be_lp_index)
              end do
           end do

           call get_xps_spectra(params%exp_data(xps_idx)%data(1,:),&
                & params%exp_data(xps_idx)%data(2,:), params&
                &%xps_sigma, params%exp_data(xps_idx)%n_samples, mag,&
                & params%exp_data(xps_idx)%x, params&
                &%exp_data(xps_idx)%y_pred, y_i_pred_all,&
                & local_properties(1:n_sites, core_be_lp_index),&
                & .true.)


           ! call get_compare_xps_spectra(params%exp_data(xps_idx)%data&
           !      & , local_properties(1:n_sites, core_be_lp_index),&
           !      & params%xps_sigma, params%exp_data(xps_idx) &
           !      &%n_samples, mag, params%exp_data(xps_idx)%similarity&
           !      & , params%exp_data(xps_idx)%x, params &
           !      &%exp_data(xps_idx)%y, params%exp_data(xps_idx) &
           !      &%y_pred, y_i_pred_all, .not. allocated(params &
           !      &%exp_data(xps_idx)%x), params%exp_similarity_type )

           ! print *, params%exp_data(xps_idx)%n_samples, xps_idx
           call get_energy_scale( params%do_md, params%do_mc,&
                & md_istep, params%md_nsteps, mc_istep, params&
                &%mc_nsteps, params &
                &%exp_energy_scales_initial(xps_idx), params &
                &%exp_energy_scales_final(xps_idx), params &
                &%exp_energy_scales(xps_idx) )

           call get_exp_pred_spectra_energies_forces( params&
                &%exp_energy_scales(xps_idx),&
                & local_properties(i_beg:i_end,core_be_lp_index),&
                & local_properties_cart_der(1:3, j_beg:j_end,&
                & core_be_lp_index ), n_neigh(i_beg:i_end),&
                & neighbors_list(j_beg:j_end), params%xps_sigma,&
                & params%exp_data(xps_idx)%n_samples, mag, params&
                &%exp_data(xps_idx)%x, params %exp_data(xps_idx)%y,&
                & params%exp_data(xps_idx) %y_pred,&
                & y_i_pred_all(i_beg:i_end, 1:params &
                &%exp_data(xps_idx)%n_samples), params %do_forces, &
                & xyz(1:3, j_beg:j_end),&
#ifdef _MPIF90
                & this_energies_lp(i_beg:i_end), this_forces_lp, this_virial_lp, params%exp_similarity_type, rank )
#else
           & energies_lp(i_beg:i_end), forces_lp, virial_lp, params%exp_similarity_type, rank )
#endif

           ! if (rank == 0)then
           !    open(unit=11, file="tg_xps.dat", status="unknown")
           !    do i = 1, params%exp_data(xps_idx)%n_samples
           !       write(11, '(1X,F20.8,1X,F20.8)') params%exp_data(xps_idx)%x(i), params%exp_data(xps_idx)%y_pred(i)
           !    end do
           !    close(11)
           ! end if


           call get_write_condition( params%do_mc, params%do_md&
                &, mc_istep, md_istep, params%write_xyz,&
                & write_condition)

           if (rank == 0 .and. params%write_exp .and. write_condition)then

              call get_overwrite_condition( params%do_mc, params%do_md&
                   &, mc_istep, md_istep, params%write_xyz, overwrite_condition)

              call write_exp_datan(params%exp_data(xps_idx)&
                   &%x(1:params%exp_data(xps_idx)%n_samples), params&
                   &%exp_data(xps_idx)%y_pred(1:params&
                   &%exp_data(xps_idx)%n_samples),&
                   & overwrite_condition, "xps_prediction.dat",&
                   & params%exp_data(xps_idx)%label)

              if ( .not.  params%exp_data(xps_idx)%wrote_exp ) then


                 call preprocess_exp_data(params, params%exp_data(xps_idx)%x,&
                      & params%exp_data(xps_idx)%y, params%exp_data(xps_idx)%label,&
                      & n_sites, dot_product( cross_product(a_box,&
                      & b_box), c_box ) / (dfloat(indices(1)*indices(2) &
                      &*indices(3)) ), params%exp_data(xps_idx)%input, exp_output, .true. )

                 call write_exp_datan(params%exp_data(xps_idx)&
                      &%x(1:params%exp_data(xps_idx)%n_samples),&
                      & params%exp_data(xps_idx)%y(1:params&
                      &%exp_data(xps_idx)%n_samples),&
                      & overwrite_condition, "xps_exp.dat" , params&
                      &%exp_data(xps_idx)%label)
                 params%exp_data(xps_idx)%wrote_exp = .true.
              end if

              ! else
              !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y_pred, mc_istep == 0, "xps_prediction.dat" )
              !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y, mc_istep == 0, "xps_exp.dat" )

           end if


           !deallocate( params%exp_data(xps_idx)%y_pred )
           if (allocated(y_i_pred_all)) deallocate(y_i_pred_all)
           ! sim_exp_pred would be an energy if multiplied by some energy scale \gamma * ( 1 - sim )
           ! sim_exp_pred_der would be the array of forces if multiplied by (- \gamma )
           deallocate(v_neigh_lp)


           call cpu_time(time_xps(2))
           time_xps(3) = time_xps(3) + time_xps(2) - time_xps(1)
           !           if (rank == 0) print *, rank, " TIME_XPS = ", time_xps(3)

           else if ( any( soap_turbo_hypers(:)%has_core_electron_be) .and. params%do_xps_standalone )then
             ! Get the linspace of the xps spectrum and then perform the
             ! calculation and write to the prediction file
             !
              if (rank == 0)then
                 call get_xps_spectra_standalone(&
                      & params%xps_e_min,&
                      & params%xps_e_max, &
                      & params%xps_sigma, &
                      & params%xps_n_samples,&
                      & x_xps, &
                      & y_xps, &
                      & local_properties(1:n_sites, core_be_lp_index))

                 call get_overwrite_condition( params%do_mc, params%do_md&
                      &, mc_istep, md_istep, params%write_xyz, overwrite_condition)

                 if (n_xyz > 0)then
                    overwrite_condition = ( n_xyz == 1 )
                 end if

                 call write_exp_datan( x_xps(1:params%xps_n_samples), &
                      & y_xps(1:params%xps_n_samples), &
                      & overwrite_condition, &
                      &"xps_prediction.dat",&
                      &"core_electron_be xps")
              end if

        end if

        !##############################################################!
        !###---   (Partial) Pair distribution functions and XRD   ---###!
        !##############################################################!

        ! We use these to calculate the (partial) structure factors, which
        ! can be used for X-Ray scattering and (in the future)
        ! neutron scattering.
        !
        ! > We use the formalism which was detailed by
        !   Gutierrez and Johansson, Physical Review B, Volume 65, 104202 (2002)
        !   such that there is consistency between the rdfs and the scattering we calculate.
        !
        ! > Furthermore, calculating the pair distribution function and
        !   the structure factor/XRD becomes much faster, as there is just
        !   a sum over species rather than a double sum over all the
        !   atomic species.
        !
        !   (The ASE implementation of XRD intensity is problematic and
        !   uses the sinc function implemented by DSP
        !   (sin(pi*x)/(pi*x)) which is not what is in the
        !   literature).
        !
        ! ***--- Steps for calculation ---***
        ! 1) We calculate the partial pair distribution functions g_ab
        ! 2) Partial static structure factors S_ab are then calculated from this by Fourier transform
        !    > This calculation can includes a window function ( sin(pi*rij/r_cut) / (pi*rij/r_cut) )
        !      such that termination effects of large sinusoids which come from the cutoff are removed.
        ! 3) If X-Ray Diffraction (xrd) is specified, then the intensity is
        !    calculated from these partial structure factors by the
        !    inclusion of the X-Ray form factors.
        !
        ! ***--- Definitions ---***
        ! There are many definitions of these various functions
        ! in the literature, however, we shall use similar ones
        ! to those in the paper of Gutierrez
        !
        ! Definiton of the structure factors given by Ashcroft and Langreth
        ! N. W. Ashcroft and D. C. Langreth. Phys. Rev., 156(3):685-692 (1967)
        !
        ! PDFs can be smoothed by using a kernel density estimate by a gaussian function when kde_sigma > 0.d0
        !
        ! R(r)    = Radial Distribution Function     === A histogram of atomic distances divided by N, goes as r^2
        ! g(r)    = Pair Distribution Function (PCF) === Scales R(r) by 1/(4 pi r^2) such that it lays flat, converges to 1.
        !         = (N_{r_l < r < r_h} / N) / ( 4 pi r^2 dr  * ( N / V ) )
        !         = n_{r_l < r < r_h} / ( dV * rho_0 )
        !           > n_{r_l < r < r_h} is the average number of particles between r_l and r_h
        !           > rho_0 is the density
        !           > dV is the differential volume between shells
        !
        ! g_ab(r) = Partial Pair Distribution Func.  === Same as above but only for particles a and b
        !         =  (N^{b}_{r_l < r < r_h} / N_b ) / ( 4 pi r^2 dr ) * ( N_a / V )
        !         With kde
        !         =  ( sum_i sum_i/=j exp( - (r - r_ij)^2 / sigma^2 / 2 ) / N_b ) / ( 4 pi r^2 dr ) * ( V / N_a ) &
        !              & * ( (r_max - r_min) / sigma / sqrt(2pi) )
        !
        !           The full pair distribution function is given by the sum (say for the binary system, with a, b atoms)
        !             g(r) = (N_a/N) * g_aa(r) + 2(N_a/N * N_b/N)g_ab + (N_b/N)g_bb
        !
        ! S_ab(q) = delta_ab + 4 pi rho (ca cb)^1/2 int_0^r_cut dr r^2 [ g_ab(r) - 1 ] sin(qr)/(qr) * sin( pi r / R )/ (pi r /R)
        !
        !
        ! XRD(q) = 1/N ( d cross_section/ d Omega )
        !
        ! Total scattering function F^X(q)
        ! F^x(q) = [ XRD(q) - \sum_n c_i f_i(q)^2 ] / [ \sum_n c_i f_i(q) ]^2

        ! First get the number of species in actuality
        n_species_actual = 0
        do i = 1, n_sites
           if (species(i) > n_species_actual) n_species_actual = n_species_actual + 1
        end do

        ! Now find the unique species ids
        if (allocated( species_types_actual)) deallocate( species_types_actual )
        allocate(species_types_actual(1:n_species_actual))

        n_species_actual = 0
        do i = 1, n_sites
           if (species(i) > n_species_actual)then
              n_species_actual = n_species_actual + 1
              species_types_actual(n_species_actual) = params%species_types( species(i) )
           end if
        end do


        if (params%do_pair_distribution)then
           call cpu_time(time_pdf(1))

           call calculate_pair_distribution( params, x_pair_distribution&
                &, y_pair_distribution, y_pair_distribution_temp,&
                & pair_distribution_partial, pair_distribution_partial_temp, &
                & n_species_actual, species_types_actual, n_atoms_of_species, n_sites, a_box,&
                & b_box, c_box, indices, md_istep, mc_istep, i_beg, i_end,&
                & j_beg, j_end, ierr , rjs, xyz, neighbors_list,&
                & n_neigh, neighbor_species, species, rank, params%exp_forces, &
                & pair_distribution_der,&
                & pair_distribution_partial_der,&
#ifdef _MPIF90
                & pair_distribution_partial_temp_der, this_energies_pdf, this_forces_pdf, this_virial_pdf)
#else
           & pair_distribution_partial_temp_der, energies_pdf, forces_pdf, virial_pdf)
#endif

           call cpu_time(time_pdf(2))
           time_pdf(3) = time_pdf(3) + time_pdf(2) - time_pdf(1)
           !           if (rank == 0) print *, rank, " TIME_PDF = ", time_pdf(3)



        end if

        ! Now calculate the structure factors
        if (params%do_structure_factor )then
           call cpu_time(time_sf(1))
           call calculate_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
                & y_structure_factor, y_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp,&
                & x_pair_distribution, y_pair_distribution, &
                & pair_distribution_partial, n_species_actual, species_types_actual, n_atoms_of_species,&
                & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
                & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
                & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
#ifdef _MPIF90
                & pair_distribution_partial_der, this_energies_sf, this_forces_sf, this_virial_sf, &
                & params%structure_factor_matrix_forces)
#else
           & pair_distribution_partial_der, energies_sf, forces_sf, virial_sf,params%structure_factor_matrix_forces)
#endif


           call cpu_time(time_sf(2))
           time_sf(3) = time_sf(3) + time_sf(2) - time_sf(1)



        end if

        if ( params%do_xrd )then
           call cpu_time(time_xrd(1))
           call calculate_xrd( params, x_xrd, x_xrd_temp,&
                & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp,&
                & n_species_actual, species_types_actual, n_atoms_of_species,&
                & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
                & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
                & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
#ifdef _MPIF90
                & pair_distribution_partial_der, this_energies_xrd,&
                & this_forces_xrd, this_virial_xrd, .false., params&
                &%structure_factor_matrix_forces)
#else
           & pair_distribution_partial_der, energies_xrd, forces_xrd,&
                & virial_xrd, .false., params&
                &%structure_factor_matrix_forces )
#endif


           call cpu_time(time_xrd(2))
           time_xrd(3) = time_xrd(3) + time_xrd(2) - time_xrd(1)

           !           if (rank == 0) print *, rank, " TIME_XRD = ", time_xrd(3)

        end if


        if ( params%do_nd )then
           call cpu_time(time_nd(1))
           call calculate_xrd( params, x_nd, x_nd_temp,&
                & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp,&
                & n_species_actual, species_types_actual, n_atoms_of_species,&
                & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
                & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
                & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
#ifdef _MPIF90
                & pair_distribution_partial_der, this_energies_nd,&
                & this_forces_nd, this_virial_nd, .true., params&
                &%structure_factor_matrix_forces)
#else
           & pair_distribution_partial_der, energies_nd, forces_nd,&
                & virial_nd, .true., params&
                &%structure_factor_matrix_forces )
#endif


           call cpu_time(time_nd(2))
           time_nd(3) = time_nd(3) + time_nd(2) - time_nd(1)

           !           if (rank == 0) print *, rank, " TIME_XRD = ", time_xrd(3)

        end if


        !################################################################!
        !###---   Compute similarity of experimental predictions   ---###!
        !################################################################!

        if ( params%do_exp )then
           do i = 1, params%n_exp
              ! First normalize the spectrum if it matches some type of experimental data

              ! Allocate the prediction data
              ! if (.not. allocated( params%exp_data(i)%y_pred ) )then
              !    allocate( params%exp_data(i)%y_pred( 1:size(params%exp_data(i)%y, 1) ))
              ! end if

              if ( trim(params%exp_data(i)%label) == 'pair_distribution' )then
                 params%exp_data(i)%y_pred = y_pair_distribution
              elseif ( trim(params%exp_data(i)%label) == 'structure_factor' )then
                 params%exp_data(i)%y_pred = y_structure_factor
              elseif ( trim(params%exp_data(i)%label) == 'xrd' )then
                 params%exp_data(i)%y_pred = y_xrd
              elseif ( trim(params%exp_data(i)%label) == 'nd' )then
                 params%exp_data(i)%y_pred = y_nd
              end if

              if ( params%exp_data(i)%compute_similarity .and. allocated(params%exp_data(i)%y) )then
                 call get_data_similarity(params%exp_data(i)%y, params&
                      &%exp_data(i)%y_pred, params&
                      &%exp_data(i)%similarity, params&
                      &%exp_similarity_type)
              end if

              ! deallocate(params%exp_data(i)%x)
              ! deallocate(params%exp_data(i)%y)
              ! deallocate(params%exp_data(i)%y_pred)

           end do
        end if





        !##############################################!
        !###---   Finalize experimental arrays   ---###!
        !##############################################!

        if (params%do_pair_distribution)then
           call finalize_pair_distribution( params, x_pair_distribution&
                &, y_pair_distribution, y_pair_distribution_temp,&
                & pair_distribution_partial, pair_distribution_partial_temp, params&
                &%exp_forces, pair_distribution_der,&
                & pair_distribution_partial_der,&
                & pair_distribution_partial_temp_der, n_atoms_of_species, rank)

        end if

        ! Now calculate the structure factors
        if (params%do_structure_factor )then

           call finalize_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
                & y_structure_factor, y_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp,&
                & x_pair_distribution, y_pair_distribution, &
                & pair_distribution_partial, sinc_factor_matrix)

        end if

        if ( params%do_xrd )then
           call finalize_xrd( params, x_xrd, x_xrd_temp,&
                & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp)
        end if

        if ( params%do_nd )then
           call finalize_xrd( params, x_nd, x_nd_temp,&
                & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
                & structure_factor_partial, structure_factor_partial_temp)
        end if

        deallocate(species_types_actual)





        if( params%do_prediction )then
           !       Loop through distance_2b descriptors
           call cpu_time(time1)

           do i = 1, n_distance_2b
              call cpu_time(time_2b(1))
              this_energies = 0.d0
              if( params%do_forces )then
                 this_forces = 0.d0
                 this_virial = 0.d0
              end if
              call get_2b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), distance_2b_hypers(i)%alphas, &
                   distance_2b_hypers(i)%cutoff, &
                   distance_2b_hypers(i)%rcut, 0.5d0, distance_2b_hypers(i)%delta, &
                   distance_2b_hypers(i)%sigma, 0.d0, distance_2b_hypers(i)%Qs(:,1), &
                   n_neigh(i_beg:i_end), params%do_forces, params%do_timing, &
                   species(i_beg:i_end), neighbor_species(j_beg:j_end), &
                   distance_2b_hypers(i)%species1, distance_2b_hypers(i)%species2, &
                   params%species_types, this_energies(i_beg:i_end), this_forces(1:3, i_beg:i_end), &
                   this_virial )
              energies_2b = energies_2b + this_energies
              if( params%do_forces )then
                 forces_2b = forces_2b + this_forces
                 virial_2b = virial_2b + this_virial
              end if
              call cpu_time(time_2b(2))
              time_2b(3) = time_2b(3) + time_2b(2) - time_2b(1)
           end do





           !       Loop through core_pot descriptors
           do i = 1, n_core_pot
              call cpu_time(time_core_pot(1))
              this_energies = 0.d0
              if( params%do_forces )then
                 this_forces = 0.d0
                 this_virial = 0.d0
              end if
              call get_core_pot_energy_and_forces(rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), &
                   core_pot_hypers(i)%x, core_pot_hypers(i)%V, &
                   core_pot_hypers(i)%yp1, core_pot_hypers(i)%ypn, &
                   core_pot_hypers(i)%dVdx2, n_neigh(i_beg:i_end), params%do_forces, &
                   params%do_timing, species(i_beg:i_end), neighbor_species(j_beg:j_end), &
                   core_pot_hypers(i)%species1, core_pot_hypers(i)%species2, &
                   params%species_types, this_energies(i_beg:i_end), this_forces(1:3, i_beg:i_end), &
                   this_virial )
              energies_core_pot = energies_core_pot + this_energies
              if( params%do_forces )then
                 forces_core_pot = forces_core_pot + this_forces
                 virial_core_pot = virial_core_pot + this_virial
              end if
              call cpu_time(time_core_pot(2))
              time_core_pot(3) = time_core_pot(3) + time_core_pot(2) - time_core_pot(1)
           end do




           !       Loop through angle_3b descriptors
           do i = 1, n_angle_3b
              call cpu_time(time_3b(1))
              this_energies = 0.d0
              if( params%do_forces )then
                 this_forces = 0.d0
                 this_virial = 0.d0
              end if
              call get_3b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), angle_3b_hypers(i)%alphas, &
                   angle_3b_hypers(i)%cutoff, &
                   angle_3b_hypers(i)%rcut, 0.5d0, angle_3b_hypers(i)%delta, &
                   angle_3b_hypers(i)%sigma, 0.d0, angle_3b_hypers(i)%Qs, n_neigh(i_beg:i_end), &
                   neighbors_list(j_beg:j_end), &
                   params%do_forces, params%do_timing, angle_3b_hypers(i)%kernel_type, &
                   species(i_beg:i_end), neighbor_species(j_beg:j_end), angle_3b_hypers(i)%species_center, &
                   angle_3b_hypers(i)%species1, angle_3b_hypers(i)%species2, params%species_types, &
                   this_energies(i_beg:i_end), this_forces, this_virial)
              energies_3b = energies_3b + this_energies
              if( params%do_forces )then
                 forces_3b = forces_3b + this_forces
                 virial_3b = virial_3b + this_virial
              end if
              call cpu_time(time_3b(2))
              time_3b(3) = time_3b(3) + time_3b(2) - time_3b(1)
           end do


           call cpu_time(time2)
           time_gap = time_gap + time2 - time1

           !       Communicate all energies and forces here for all
           !       terms
#ifdef _MPIF90
           call cpu_time(time_mpi_ef(1))
           counter2 = 0
           if( n_soap_turbo > 0 )then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_vdw) )then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_lp) )then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_pdf) .and. params%valid_pdf )then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_sf) .and. params%valid_sf)then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_xrd) .and. params%valid_xrd)then
              counter2 = counter2 + 1
           end if
           if( allocated(this_energies_nd) .and. params%valid_nd)then
              counter2 = counter2 + 1
           end if
           if( n_distance_2b > 0 )then
              counter2 = counter2 + 1
           end if
           if( n_core_pot > 0 )then
              counter2 = counter2 + 1
           end if
           if( n_angle_3b > 0 )then
              counter2 = counter2 + 1
           end if


           !       It would probably be faster to use pointers for this
           allocate( all_energies(1:n_sites, 1:counter2) )
           allocate( all_this_energies(1:n_sites, 1:counter2) )
           if( params%do_forces )then
              allocate( all_forces(1:3, 1:n_sites, 1:counter2) )
              allocate( all_this_forces(1:3, 1:n_sites, 1:counter2) )
              allocate( all_virial(1:3, 1:3, 1:counter2) )
              allocate( all_this_virial(1:3, 1:3, 1:counter2) )
           end if
           counter2 = 0
           if( n_soap_turbo > 0 )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = energies_soap(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = forces_soap(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = virial_soap(1:3, 1:3)
              end if
           end if
           if( allocated(this_energies_vdw) )then
              counter2 = counter2 + 1
              !         Note the vdw things have "this" in front
              all_energies(1:n_sites, counter2) = this_energies_vdw(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_vdw(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_vdw(1:3, 1:3)
              end if
           end if

           if( allocated(this_energies_lp) )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = this_energies_lp(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_lp(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_lp(1:3, 1:3)
              end if
           end if

           if( allocated(this_energies_pdf) .and. params%valid_pdf )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = this_energies_pdf(1:n_sites)
              if( params%do_forces .and. params%exp_forces)then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_pdf(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_pdf(1:3, 1:3)
              end if
           end if

           if( allocated(this_energies_sf) .and. params%valid_sf)then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = this_energies_sf(1:n_sites)
              if( params%do_forces .and. params%exp_forces)then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_sf(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_sf(1:3, 1:3)
              end if
           end if

           if( allocated(this_energies_xrd) .and. params%valid_xrd )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = this_energies_xrd(1:n_sites)
              if( params%do_forces .and. params%exp_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_xrd(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_xrd(1:3, 1:3)
              end if
           end if

           if( allocated(this_energies_nd) .and. params%valid_nd )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = this_energies_nd(1:n_sites)
              if( params%do_forces .and. params%exp_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = this_forces_nd(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = this_virial_nd(1:3, 1:3)
              end if
           end if


           if( n_distance_2b > 0 )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = energies_2b(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = forces_2b(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = virial_2b(1:3, 1:3)
              end if
           end if
           if( n_core_pot > 0 )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = energies_core_pot(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = forces_core_pot(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = virial_core_pot(1:3, 1:3)
              end if
           end if
           if( n_angle_3b > 0 )then
              counter2 = counter2 + 1
              all_energies(1:n_sites, counter2) = energies_3b(1:n_sites)
              if( params%do_forces )then
                 all_forces(1:3, 1:n_sites, counter2) = forces_3b(1:3, 1:n_sites)
                 all_virial(1:3, 1:3, counter2) = virial_3b(1:3, 1:3)
              end if
           end if

           !       Here we communicate
           call mpi_reduce(all_energies, all_this_energies, n_sites&
                &*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                & MPI_COMM_WORLD, ierr)
           if( params%do_forces )then
              call mpi_reduce(all_forces, all_this_forces, 3*n_sites&
                   &*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                   & MPI_COMM_WORLD, ierr)
              call mpi_reduce(all_virial, all_this_virial, 9*counter2&
                   &, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                   & MPI_COMM_WORLD, ierr)
           end if

           !       Here we give proper names to the quantities - again, pointers would probably be faster
           counter2 = 0
           if( n_soap_turbo > 0 )then
              counter2 = counter2 + 1
              energies_soap(1:n_sites) = all_this_energies(1:n_sites, counter2)
              if( params%do_forces )then
                 forces_soap(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_soap(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
              end if
           end if
           if( allocated(this_energies_vdw) )then
              counter2 = counter2 + 1
              !         Note the vdw things DO NOT have "this" in front anymore
              energies_vdw(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_vdw)
              if( params%do_forces )then
                 forces_vdw(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_vdw(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_vdw)
              end if
           end if
           if( allocated(this_energies_lp) )then
              counter2 = counter2 + 1
              energies_lp(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_lp)
              if( params%do_forces )then
                 forces_lp(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_lp(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_lp)
              end if
           end if
           if( allocated(this_energies_pdf) .and. params%valid_pdf)then
              counter2 = counter2 + 1
              energies_pdf(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_pdf)
              if( params%do_forces .and. params%exp_forces)then
                 forces_pdf(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_pdf(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_pdf)
              end if
           end if
           if( allocated(this_energies_sf) .and. params%valid_sf)then
              counter2 = counter2 + 1
              energies_sf(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_sf)
              if( params%do_forces .and. params%exp_forces)then
                 forces_sf(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_sf(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_sf)
              end if
           end if
           if( allocated(this_energies_xrd) .and. params%valid_xrd)then
              counter2 = counter2 + 1
              energies_xrd(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_xrd)
              if( params%do_forces  .and. params%exp_forces)then
                 forces_xrd(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_xrd(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_xrd)
              end if
           end if
           if( allocated(this_energies_nd) .and. params%valid_nd)then
              counter2 = counter2 + 1
              energies_nd(1:n_sites) = all_this_energies(1:n_sites, counter2)
              deallocate(this_energies_nd)
              if( params%do_forces  .and. params%exp_forces)then
                 forces_nd(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_nd(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
                 deallocate(this_forces_nd)
              end if
           end if

           if( n_distance_2b > 0 )then
              counter2 = counter2 + 1
              energies_2b(1:n_sites) = all_this_energies(1:n_sites, counter2)
              if( params%do_forces)then
                 forces_2b(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_2b(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
              end if
           end if
           if( n_core_pot > 0 )then
              counter2 = counter2 + 1
              energies_core_pot(1:n_sites) = all_this_energies(1:n_sites, counter2)
              if( params%do_forces )then
                 forces_core_pot(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_core_pot(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
              end if
           end if
           if( n_angle_3b > 0 )then
              counter2 = counter2 + 1
              energies_3b(1:n_sites) = all_this_energies(1:n_sites, counter2)
              if( params%do_forces )then
                 forces_3b(1:3, 1:n_sites) = all_this_forces(1:3, 1:n_sites, counter2)
                 virial_3b(1:3, 1:3) = all_this_virial(1:3, 1:3, counter2)
              end if
           end if

           !       Clean up
           deallocate( all_energies, all_this_energies )
           if( params%do_forces )then
              deallocate( all_forces, all_this_forces, all_virial, all_this_virial )
           end if

           call cpu_time(time_mpi_ef(2))
           time_mpi_ef(3) = time_mpi_ef(3) + time_mpi_ef(2) - time_mpi_ef(1)
#endif


           !       Add up all the energy terms
           energies = energies + energies_soap + energies_2b +&
                & energies_3b + energies_core_pot + energies_vdw !+energies_lp

           if ( valid_xps )                                          energies_exp = energies_exp + energies_lp
           if ( params%valid_pdf .and. params%do_pair_distribution ) energies_exp = energies_exp + energies_pdf
           if ( params%valid_sf .and. params%do_structure_factor )   energies_exp = energies_exp + energies_sf
           if ( params%valid_xrd .and. params%do_xrd )               energies_exp = energies_exp + energies_xrd
           if ( params%valid_nd .and. params%do_nd )               energies_exp = energies_exp + energies_nd

           if (params%exp_energies) energies = energies + energies_exp

           energy_prev = energy
           instant_pressure_prev = instant_pressure
           energy = sum(energies)
           energy_exp = sum( energies_exp )

        end if


        if( .not. params%do_md .and. .not. params%do_mc) then
#ifdef _MPIF90
           IF( rank == 0 )then
#endif
              write(*,*)'                                       |'
              write(*,'(A,1X,F22.8,1X,A)')' SOAP energy:', sum(energies_soap), 'eV |'
              write(*,'(A,1X,F24.8,1X,A)')' 2b energy:', sum(energies_2b), 'eV |'
              write(*,'(A,1X,F24.8,1X,A)')' 3b energy:', sum(energies_3b), 'eV |'
              write(*,'(A,1X,F18.8,1X,A)')' core_pot energy:', sum(energies_core_pot), 'eV |'
              write(*,'(A,1X,F23.8,1X,A)')' vdw energy:', sum(energies_vdw), 'eV |'
              write(*,'(A,1X,F22.8,1X,A)')' Exp. energy:', sum(energies_exp), 'eV |'
              if (valid_xps) write(*,'(A,1X,F23.8,1X,A)')' xps energy:', sum(energies_lp), 'eV |'
              if ( params%valid_pdf .and. params%do_pair_distribution )&
                   & write(*,'(A,1X,F23.8,1X,A)')' pdf energy:',&
                   & sum(energies_pdf), 'eV |'
              if ( params%valid_sf .and. params%do_structure_factor )&
                   & write(*,'(A,1X,F24.8,1X,A)')' sf energy:',&
                   & sum(energies_sf), 'eV |'
              if ( params%valid_xrd .and. params%do_xrd )&
                   & write(*,'(A,1X,F23.8,1X,A)')' xrd energy:',&
                   & sum(energies_xrd), 'eV |'
              if ( params%valid_nd .and. params%do_nd )&
                   & write(*,'(A,1X,F23.8,1X,A)')' nd energy:',&
                   & sum(energies_nd), 'eV |'


              if (.not. params%do_mc .or. (params%do_mc .and.  mc_istep <= 1 ))then
                 write(*,'(A,1X,F21.8,1X,A)')' Total energy:', sum(energies), 'eV |'
              else
                 write(*,'(A,1X,F21.8,1X,A)')' Total energy:', sum(images(i_trial_image)%energies), 'eV |'
              end if

              if ( .not. params%do_mc)then
                 write(*,*)'                                       |'
                 write(*,*)'Energy & forces in "trajectory_out.xyz"|'
                 write(*,*)'                                       |'
                 write(*,*)'.......................................|'
              else if ( mc_istep == 0 )then
                 write(*,*)'                                       |'
                 write(*,*)' MC configs in "mc_current.xyz" and    |'
                 write(*,*)'               "mc_trial.xyz"          |'
                 write(*,*)'               "mc_all.xyz"            |'
                 write(*,*)'.......................................|'
              end if
#ifdef _MPIF90
           END IF
#endif
        end if

        if( params%do_forces )then
           forces = forces_soap + forces_2b + forces_3b + forces_core_pot + forces_vdw
           virial = virial_soap + virial_2b + virial_3b + virial_core_pot + virial_vdw

           if (params%exp_forces .and. valid_xps)        forces = forces + forces_lp
           if (params%exp_forces .and. valid_xps)        virial = virial + virial_lp

           if (params%exp_forces .and. params%valid_pdf) forces = forces + forces_pdf
           if (params%exp_forces .and. params%valid_pdf) virial = virial + virial_pdf

           if (params%exp_forces .and. params%valid_sf ) forces = forces + forces_sf
           if (params%exp_forces .and. params%valid_sf ) virial = virial + virial_sf

           if (params%exp_forces .and. params%valid_xrd) forces = forces + forces_xrd
           if (params%exp_forces .and. params%valid_xrd) virial = virial + virial_xrd

           if (params%exp_forces .and. params%valid_nd) forces = forces + forces_nd
           if (params%exp_forces .and. params%valid_nd) virial = virial + virial_nd

           if ( params%print_vdw_forces )then
              open(unit=90, file="forces_vdw", status="unknown")
              do i = 1, n_sites
                 write(90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                      forces_vdw(1,i), forces_vdw(2,i), forces_vdw(3,i)
              end do
              close(90)

           end if

        end if
        ! For debugging the virial implementation
        if( rank == 0 .and. .false. )then
           write(*,*) "pressure_soap: ", virial_soap / 3.d0 / v_uc
           write(*,*) "pressure_vdw: ", virial_vdw / 3.d0 / v_uc
           write(*,*) "pressure_lp: ", virial_lp / 3.d0 / v_uc
           write(*,*) "pressure_2b: ", virial_2b / 3.d0 / v_uc
           write(*,*) "pressure_3b: ", virial_3b / 3.d0 / v_uc
           write(*,*) "pressure_core_pot: ", virial_core_pot / 3.d0 / v_uc
        end if



        if( params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc)then
#ifdef _MPIF90
           IF( rank == 0 )then
#endif
              !       Write energy and forces if we're just doing static predictions
              !       The masses should be divided by 103.6426965268d0 to have amu units, but
              !       since masses is not allocated for single point calculations, it would
              !       likely lead to a segfault
              call wrap_pbc(positions(1:3,1:n_sites), a_box&
                   &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                   & c_box/dfloat(indices(3)))
              call get_xyz_energy_string(energies_soap, energies_2b,&
                   & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                   &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                   & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                   & params%do_structure_factor, params%do_xrd, params%do_nd, string)

              call write_extxyz( n_sites, -n_xyz, md_time, time_step,&
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
        IF( rank == 0 )then
#endif
           !     Do nothing
           write(*,*)'                                       |'
           write(*,*)'You didn''t ask me to do anything!      |'
           write(*,*)'                                       |'
           write(*,*)'.......................................|'
#ifdef _MPIF90
        END IF
#endif
     end if
     !**************************************************************************






     !**************************************************************************
     !   Do MD stuff here
#ifdef _MPIF90
     IF( rank == 0 )THEN
#endif
        if( params%do_md .and. md_istep > -1)then
           call cpu_time(time_md(1))
           !     Define the time_step and md_time prior to possible scaling (see variable_time_step below)
           if( md_istep > 0 )then
              md_time = md_time + time_step
           else
              md_time = 0.d0
              time_step = params%md_step
           end if
           !     We wrap the positions and remoce CM velocity
           call wrap_pbc(positions(1:3,1:n_sites), a_box&
                &/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box&
                &/dfloat(indices(3)))
           call remove_cm_vel(velocities(1:3,1:n_sites),&
                & masses(1:n_sites))

           !     First we check if this is a variable time step simulation
           if( params%variable_time_step )then
              call variable_time_step(md_istep == 0, velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
                   params%target_pos_step, params%tau_dt, params%md_step, time_step)
           end if


	   !! ------- option for radiation cascade simulation with electronic stopping

	  if ( params%electronic_stopping ) then
		call electron_stopping_velocity_dependent (md_istep, n_species, params%eel_cut, params%eel_freq_out, &
					velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
					params%masses_types, time_step, md_time, nrows, allelstopdata, cum_EEL, 'forces')		
	  end if

	   !! -----------------------------------	******** until here for electronic stopping


	   !! ------- option for electronic stopping based on eph model

	  if ( params%nonadiabatic_processes ) then
		call ephlsc%eph_LangevinForces (velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), &
					masses(1:n_sites), params%masses_types, md_istep, time_step, md_time, &
					positions(1:3, 1:n_sites), n_species, ephbeta, ephfdm)	
	  end if

	  !! -----------------------------------	******** until here for electronic stopping basd on eph model


	  !! ------- option for doing simulation with adaptive time step

	  if ( params%adaptive_time ) then
		if (MOD(md_istep, params%adapt_tstep_interval) == 0) then
			call variable_time_step_adaptive (md_istep == 0, velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), &
						masses(1:n_sites), params%adapt_tmin, params%adapt_tmax, params%adapt_xmax, &
						params%adapt_emax, params%md_step, time_step)
		end if
	  end if

	  !! ----------------------------------	******** until here for adaptive time
           
           !     This takes care of NVE
           !     Velocity Verlet takes positions for t, positions_prev for t-dt, and velocities for t-dt and returns everything
           !     dt later. forces are taken at t, and forces_prev at t-dt. forces is left unchanged by the routine, and
           !     forces_prev is returned as equal to forces (both arrays contain the same information on return)
           if( params%optimize == "vv")then
              call velocity_verlet(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                   forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), time_step, time_step_prev, &
                   md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                   fix_atom(1:3, 1:n_sites))
           else if( params%optimize == "gd" )then
              call gradient_descent(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                   forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), &
                   params%max_opt_step, md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), &
                   c_box/dfloat(indices(3)), fix_atom(1:3, 1:n_sites), energy)
           else if( (params%optimize == "gd-box" .or. params%optimize == "gd-box-ortho") .and. gd_box_do_pos)then
              !       We propagate the positions
              call gradient_descent(positions(1:3, 1:n_sites),&
                   & positions_prev(1:3, 1:n_sites), velocities(1:3,&
                   & 1:n_sites), forces(1:3, 1:n_sites),&
                   & forces_prev(1:3, 1:n_sites), masses(1:n_sites),&
                   & params%max_opt_step, gd_istep == 0, a_box&
                   &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                   & c_box/dfloat(indices(3)), fix_atom(1:3,&
                   & 1:n_sites), energy)
              if( gd_istep > 1 .and. abs(energy-energy_prev) < params&
                   &%e_tol*dfloat(n_sites) .and. maxval(forces) <&
                   & params%f_tol )then
                 !         If the position optimization is converged
                 !         (energy only) we set the code to do the
                 !         box relaxation (below)
                 gd_box_do_pos = .false.
                 gd_istep = 0
              else
                 gd_istep = gd_istep + 1
              end if
           else
              !       If nothing happens we still update these variables
              positions_prev(1:3, 1:n_sites) = positions(1:3, 1:n_sites)
              forces_prev(1:3, 1:n_sites) = forces(1:3, 1:n_sites)
           end if
           
	!! ------- option for radiation cascade simulation with electronic stopping

	  if ( params%electronic_stopping ) then
		call electron_stopping_velocity_dependent (md_istep, n_species, params%eel_cut, params%eel_freq_out, &
					velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
					params%masses_types, time_step, md_time, nrows, allelstopdata, cum_EEL, 'energy')
	  end if

	!! -----------------------------------		******** until here for electronic stopping


	!! ------- option for electronic stopping based on eph model

	  if ( params%nonadiabatic_processes ) then
		call ephlsc%eph_LangevinEnergyDissipation (md_istep, md_time, velocities(1:3, 1:n_sites), &
				positions(1:3, 1:n_sites), time_step, ephfdm)	
	  end if

	!! -----------------------------------		******** until here for electronic stopping basd on eph model

           !     Compute kinetic energy from current velocities. Because Velocity Verlet
           !     works with the velocities at t-dt (except for the first time step) we
           !     have to compute the velocities after call Verlet
           E_kinetic = 0.d0
           do i = 1, n_sites
              E_kinetic = E_kinetic + 0.5d0 * masses(i) * dot_product(velocities(1:3, i), velocities(1:3, i))
           end do
           instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic

           !     Instant pressure in bar
           instant_pressure = (kB*dfloat(n_sites-1)*instant_temp&
                &+(virial(1,1) + virial(2,2) + virial(3,3))/3.d0)&
                &/v_uc*eVperA3tobar
           instant_pressure_tensor(1:3, 1:3) = virial(1:3,1:3)/v_uc&
                &*eVperA3tobar
           do i = 1, 3
              instant_pressure_tensor(i, i) =&
                   & instant_pressure_tensor(i, i) + (kB&
                   &*dfloat(n_sites-1)*instant_temp)/v_uc*eVperA3tobar
           end do

           !     Here we write thermodynamic information -> THIS NEEDS CLEAN UP AND IMPROVEMENT
           if( md_istep == 0 .and. .not. params%do_nested_sampling )then
              open(unit=10, file="thermo.log", status="unknown")
              write(10,"(A,A)") "#     Step             Time      Temperature                E_kin                     E_pot", &
                                "             Pressure"
           else if( md_istep == 0 .and. i_nested == 1 )then
              open(unit=10, file="thermo.log", status="unknown")
              write(10,"(A,A)") "#     Step             Time      Temperature                E_kin                     E_pot", &
                                "             Pressure"
           else
              open(unit=10, file="thermo.log", status="old", position="append")
           end if
           if( .not. params%do_mc .and. (md_istep == 0 .or. md_istep == params%md_nsteps &
                .or. modulo(md_istep, params%write_thermo) == 0) )then
              !       Organize this better so that the user can have more freedom about what gets printed to thermo.log
              !       There should also be a header preceded by # specifying what gets printed
              if (params%do_exp) then
                 write(10, "(I10, 1X, F16.4, 1X, F16.4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8)", advance="no") &
                      md_istep, md_time, instant_temp, E_kinetic, sum(energies), sum(energies_exp), instant_pressure
              else
                 write(10, "(I10, 1X, F16.4, 1X, F16.4, 1X, F20.8, 1X, F20.8, 1X, F20.8)", advance="no") &
                      md_istep, md_time, instant_temp, E_kinetic, sum(energies), instant_pressure
              end if

              if( params%write_lv )then
                 write(10, "(1X, 9F20.8)", advance="no") a_box(1:3)/dfloat(indices(1)), &
                      b_box(1:3)/dfloat(indices(2)), &
                      c_box(1:3)/dfloat(indices(3))
              end if
              !       Further printouts should go here
              !       <<HERE>>
              !
              !       This is to make the pointer advance
              write(10, *)
           end if
           close(10)
           !
           !     Check if we have converged a relaxation calculation
           !     Check if we have converged a relaxation calculation
           if( params%do_md .and. params%optimize == "gd" .and. md_istep > 0 .and. &
                abs(energy-energy_prev) < params%e_tol*dfloat(n_sites) .and. &
                maxval(forces) < params%f_tol .and. rank == 0 )then
              exit_loop = .true.
              if (params%do_mc) exit_loop=.false.
              !     THIS CONDITION ON INSTANT PRESSURE WILL NEED TO BE FINE TUNED, TO ACCOUNT FOR ARBITRARY TARGET PRESSURES
              !     BUT ALSO TO ACCOMMODATE NON-TRICLINIC TARGET BOX SHAPES, WHERE IT MIGHT NOT BE POSSIBLE TO CONVERGE THE
              !     TOTAL PRESSURE BELOW A CERTAIN MINIMUM (DUE TO THE BOX SHAPE CONSTRAINTS)
           else if( params%do_md .and. (params%optimize == "gd-box"&
                & .or. params%optimize == "gd-box-ortho") .and.&
                & gd_istep > 1 .and. abs(energy-energy_prev) < params&
                &%e_tol*dfloat(n_sites) .and. abs(instant_pressure -&
                & instant_pressure_prev) < params%p_tol .and.&
                & maxval(abs(forces)) < params%f_tol .and. rank == 0&
                & )then
              exit_loop = .true.
              if (params%do_mc ) exit_loop=.false.
           end if

           !     We write out the trajectory file. We write positions_prev which is the one for which we have computed
           !     the properties. positions_prev and velocities are synchronous
           if( (md_istep == 0 .and. .not. params%do_nested_sampling) .or. &
                (md_istep == params%md_nsteps .and. .not. params%do_nested_sampling) &
                .or. (modulo(md_istep, params%write_xyz) == 0  .and. .not. params%do_nested_sampling) .or. &
                exit_loop )then
              call wrap_pbc(positions_prev(1:3,1:n_sites), a_box&
                   &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                   & c_box/dfloat(indices(3)))
              call get_xyz_energy_string(energies_soap, energies_2b,&
                   & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                   &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                   & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                   & params%do_structure_factor, params%do_xrd, params%do_nd, string)

              call write_extxyz( n_sites, md_istep, md_time, time_step,&
                   & instant_temp, instant_pressure, a_box&
                   &/dfloat(indices(1)), b_box/dfloat(indices(2)),&
                   & c_box/dfloat(indices(3)), virial, xyz_species,&
                   & positions_prev(1:3, 1:n_sites), velocities,&
                   & forces, energies(1:n_sites), masses, params&
                   &%write_property, params %write_array_property,&
                   & params %write_local_properties,&
                   & local_property_labels, local_properties,&
                   & fix_atom(1:3, 1:n_sites), "trajectory_out.xyz", string, &
                   & md_istep == 0 )
           else if( md_istep == params%md_nsteps .and. params%do_nested_sampling )then
              write(cjunk,'(I8)') i_image
              write(filename,'(A,A,A)') "walkers/", trim(adjustl(cjunk)), ".xyz"
              call wrap_pbc(positions_prev(1:3,1:n_sites), &
                   a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
              call get_xyz_energy_string(energies_soap, energies_2b,&
                   & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                   &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                   & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                   & params%do_structure_factor, params%do_xrd, params%do_nd, string)

              call write_extxyz( n_sites, md_istep, md_time, time_step, instant_temp, instant_pressure, &
                   a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                   virial, xyz_species, &
                   positions_prev(1:3, 1:n_sites), velocities, &
                   forces, energies(1:n_sites), masses, &
                   params%write_property, params%write_array_property&
                   &,params%write_local_properties,&
                   & local_property_labels, local_properties,&
                   & fix_atom(1:3, 1:n_sites), filename, string, .true. )

           end if
           !
           !     If there are pressure/box rescaling operations they happen here
           if( params%scale_box )then
              call box_scaling(positions(1:3, 1:n_sites), a_box(1:3), b_box(1:3), c_box(1:3), &
                   indices, md_istep, params%md_nsteps, params%box_scaling_factor)
           else if( params%barostat == "berendsen" )then
              lv(1:3, 1) = a_box(1:3)
              lv(1:3, 2) = b_box(1:3)
              lv(1:3, 3) = c_box(1:3)
              call berendsen_barostat(lv(1:3,1:3), &
                   params%p_beg + (params%p_end-params%p_beg)*dfloat(md_istep+1)/float(params%md_nsteps), &
                   instant_pressure_tensor, params%barostat_sym, params%tau_p, params%gamma_p, time_step)
              a_box(1:3) = lv(1:3, 1)
              b_box(1:3) = lv(1:3, 2)
              c_box(1:3) = lv(1:3, 3)
              call berendsen_barostat(positions(1:3, 1:n_sites), &
                   params%p_beg + (params%p_end-params%p_beg)*dfloat(md_istep+1)/float(params%md_nsteps), &
                   instant_pressure_tensor, params%barostat_sym, params%tau_p, params%gamma_p, time_step)
           else if( (params%optimize == "gd-box" .or. params%optimize == "gd-box-ortho") &
                .and. .not. gd_box_do_pos )then
              if( gd_istep > 1 .and. ( ( abs(energy-energy_prev) < params%e_tol*dfloat(n_sites) &
                   .and. abs(instant_pressure - instant_pressure_prev) < params%p_tol ) &
                   .or. restart_box_optim) )then
                 gd_box_do_pos = .true.
                 gd_istep = 0
              else
                 !         We rewind positions and forces because they were already updated above
                 positions(1:3, 1:n_sites) = positions_prev(1:3, 1:n_sites)
                 forces(1:3, 1:n_sites) = forces_prev(1:3, 1:n_sites)
                 !
                 a_box = a_box/dfloat(indices(1))
                 b_box = b_box/dfloat(indices(2))
                 c_box = c_box/dfloat(indices(3))
                 call gradient_descent_box(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), &
                      velocities(1:3, 1:n_sites), &
                      forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), &
                      params%max_opt_step_eps, gd_istep == 0, a_box, b_box, c_box, energy, &
                      [virial(1,1), virial(2,2), virial(3,3), virial(2,3), virial(1,3), virial(1,2)], &
                      params%optimize, restart_box_optim )
                 a_box = a_box*dfloat(indices(1))
                 b_box = b_box*dfloat(indices(2))
                 c_box = c_box*dfloat(indices(3))
                 gd_istep = gd_istep + 1
              end if
           end if
           !     If there are thermostating operations they happen here
           if( params%thermostat == "berendsen" )then
              call get_target_temp(  params%t_beg,  params%t_end,&
                   & md_istep,  params%md_nsteps,  params%n_t_hold, &
                   & params%t_hold, target_temp )
              call berendsen_thermostat(velocities(1:3, 1:n_sites), &
                   target_temp, &
                   instant_temp, params%tau_t, time_step)
           else if( params%thermostat == "bussi" )then
              call get_target_temp(  params%t_beg,  params%t_end,&
                   & md_istep,  params%md_nsteps,  params%n_t_hold, &
                   & params%t_hold, target_temp )
              velocities(1:3, 1:n_sites) = velocities(1:3, 1:n_sites) * dsqrt(resamplekin(E_kinetic, &
                   target_temp, &
                   3*n_sites-3,params%tau_t, time_step) / E_kinetic)
           end if
           !     Check what's the maximum atomic displacement since last neighbors build
           positions_diff = positions_diff + positions(1:3, 1:n_sites) - positions_prev(1:3, 1:n_sites)
           rebuild_neighbors_list = .false.
           !--------
           ! CHECK THIS OUT and fix it at some point
           ! Here we set the neighbors list rebuild to always true if the supercell and the primitive unit cell are not
           ! the same. This is because of how atoms get wrapped around the PBC during MD (they get wrapped around the
           ! primitive unit cell) making the neighbors lists obsolete. This is an issue with wrapping, not with the neighbors
           ! lists. A possible solution would be to wrap around the supercell, instead of the unit cell, to maintain the
           ! internal consistency of the positions(:,:) array, and then do a wrapping around the primitive unit cell for
           ! printing the XYZ coordinates only (i.e., keeping the other wrapping convention internally for positions(:,:))
           if( any(indices > 1) )then
              rebuild_neighbors_list = .true.
           end if
           !--------
           do i = 1, n_sites
              if( positions_diff(1, i)**2 + positions_diff(2, i)**2 + positions_diff(3, i)**2 > params%neighbors_buffer/2.d0 )then
                 rebuild_neighbors_list = .true.
                 positions_diff = 0.d0
                 exit
              end if
           end do
           !       We make sure the atoms in the supercell have the same positions and velocities as in the unit cell
           j = 0
           do i2 = 1, indices(1)
              do j2 = 1, indices(2)
                 do k2 = 1, indices(3)
                    do i = 1, n_sites
                       j = j + 1
                       if( j > n_sites )then
                          positions(1:3, j) = positions(1:3, i) + dfloat(i2-1)/dfloat(indices(1))*a_box &
                               + dfloat(j2-1)/dfloat(indices(2))*b_box &
                               + dfloat(k2-1)/dfloat(indices(3))*c_box
                          velocities(1:3, j) = velocities(1:3, i)
                       end if
                    end do
                 end do
              end do
           end do
           call cpu_time(time_md(2))
           time_md(3) = time_md(3) + time_md(2) - time_md(1)
        end if
#ifdef _MPIF90
     END IF
     call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     !   Make sure all ranks have correct positions and velocities
#ifdef _MPIF90
     if( params%do_md )then
        call cpu_time(time_mpi_positions(1))
        n_pos = size(positions,2)
        call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(velocities, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call cpu_time(time_mpi_positions(2))
        time_mpi_positions(3) = time_mpi_positions(3) + time_mpi_positions(2) - time_mpi_positions(1)
     end if
#endif

     !**************************************************************************
     !   Nested sampling
     !   PUT THIS INTO A MODULE!!!!!!!!!!!!!!

     !   This runs at the beginning to read in the initial images
     if( params%do_nested_sampling .and. n_xyz > i_image .and. .not. params%do_md )then
        i_image = i_image + 1
        if( .not. allocated( images ) )then
           allocate( images(1:i_image) )
        else
           allocate( images_temp(1:i_image) )
           images_temp(1:i_image-1) = images(1:i_image-1)
           deallocate( images )
           allocate( images(1:i_image) )
           images = images_temp
           deallocate(images_temp)
        end if
        !     Save initial pool of structures
        velocities = 0.d0
        call from_properties_to_image(images(i_image), positions, velocities, masses, &
             forces, a_box, b_box, c_box, energy, energies, energy_exp,  E_kinetic, &
             species, species_supercell, n_sites, indices, fix_atom, &
             xyz_species, xyz_species_supercell, local_properties)
     end if

     !   This handles the nested sampling iterations after all images have
     !   been read and their energies computed
     if( params%do_nested_sampling .and. .not. repeat_xyz )then
        if( i_nested == 0 )then
           md_istep = -1
           params%write_xyz = params%md_nsteps
           params%do_md = .true.
           if( rank == 0 )then
              write(*,*)'                                       |'
              write(*,*)'Running nested sampling algorithm with |'
              write(*,'(1X,I6,A)') n_xyz, ' walkers.                        |'
              write(*,*)'                                       |'
              write(*,*)'Target pressure in nested sampling:    |'
              write(*,'(A,ES15.7,A)') ' P = ', params%p_nested, ' bar.               |'
              write(*,*)'                                       |'
              write(*,*)'[P = 0 means total energy, rather than |'
              write(*,*)'total enthalphy, simulation]           |'
           end if
        end if
        !     At the end of the MD/MC moves we add the image to the pool if its energy has decreased
        if( md_istep == params%md_nsteps )then
           md_istep = -1
           velocities = 0.d0
           !       Unit cell volume
           v_uc = dot_product( cross_product(a_box, b_box), c_box ) / (dfloat(indices(1)*indices(2)*indices(3)))
           !       We check enthalpy, not internal energy (they are the same for P = 0)
           if( energy + E_kinetic + params%p_nested/eVperA3tobar*v_uc < e_max )then
              call from_properties_to_image(images(i_image), positions, velocities, masses, &
                   forces, a_box, b_box, c_box, energy, energies, energy_exp,  E_kinetic, &
                   species, species_supercell, n_sites, indices, fix_atom, &
                   xyz_species, xyz_species_supercell, local_properties)
           end if
        end if
        !     This selects the highest energy image from the pool
        if( md_istep == -1 .and. i_nested < params%n_nested )then
           i_nested = i_nested + 1
           rebuild_neighbors_list = .true.
           i_max = 0
           e_max = -1.d100
           do i = 1, n_xyz
              v_uc = dot_product( cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box ) / &
                   (dfloat(images(i)%indices(1)*images(i)%indices(2)*images(i)%indices(3)))
              !         We check enthalpy, not potential energy (they are the same for P = 0)
              if( images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*v_uc > e_max )then
                 e_max = images(i)%energy + images(i)%e_kin + params%p_nested/eVperA3tobar*v_uc
                 i_max = i
              end if
           end do
           i_image = i_max
           deallocate( positions, velocities, masses, forces, species, &
                species_supercell, fix_atom, xyz_species, xyz_species_supercell )
           !       Make a copy of a randonmly chosen image which is not i_image
           if( n_xyz == 1 )then
              i = i_image
           else
              i = i_image
              do while( i == i_image )
                 i = mod(irand(), n_xyz) + 1
              end do
           end if
           if( rank == 0 )then
              counter = 1
              !          write(*,*)
              write(*,*)'                                       |'
              write(*,'(A,I8,A,I8,A)') "Nested sampling iter.:", i_nested, "/", params%n_nested, " |"
              write(*,'(A,I8,A)') " - Highest enthalpy walker:    ", i_image, " |"
              write(*,'(A,I8,A)') " - Walker selected for cloning:", i, " |"
              write(*,'(A,F15.7,A)') " - Max. enthalpy: ", e_max, " eV |"
           end if
           call from_image_to_properties(images(i), positions, velocities, masses, &
                forces, a_box, b_box, c_box, energy, energies, energy_exp,  E_kinetic, &
                species, species_supercell, n_sites, indices, fix_atom, &
                xyz_species, xyz_species_supercell, local_properties)
           v_uc = dot_product( cross_product(images(i)%a_box, images(i)%b_box), images(i)%c_box ) / &
                (dfloat(images(i)%indices(1)*images(i)%indices(2)*images(i)%indices(3)))
           !       This only gets triggered if we are doing box rescaling, i.e., if the target nested sampling pressure (*not* the
           !       actual pressure for the atomic configuration) is > 0
!!!!!!!!!!!!!!!!!!!!!!!!!! Temporary hack
           if( params%scale_box_nested )then
              params%scale_box = .true.
              call random_number(rand_scale)
!!!!!!!!!!!!!!! The size of the scaling should also decrease as we reach convergence (otherwise all trial moves will be rejected)
!!!!!!!!!!!!!!! Finally, there should be a limit for the acceptable aspect ratio of the simulation box
              rand_scale = 2.d0*(rand_scale - 0.5d0) * params%nested_max_strain
              params%box_scaling_factor = reshape([1.d0+rand_scale(1), rand_scale(6)/2.d0, rand_scale(5)/2.d0, &
                   rand_scale(6)/2.d0, 1.d0+rand_scale(2), rand_scale(4)/2.d0, &
                   rand_scale(5)/2.d0, rand_scale(4)/2.d0, 1.d0+rand_scale(3)], [3,3])
              ! Make the transformation volume-preserving
              call volume_preserving_strain_transformation(a_box, b_box, c_box, params%box_scaling_factor)
              ! Volume scaling
              call get_ns_unbiased_volume_proposal(1.d0-params%nested_max_volume_change, &
                   1.d0+params%nested_max_volume_change, n_sites, rand)
              params%box_scaling_factor = params%box_scaling_factor * (rand)**(1.d0/3.d0)
              ! Each MPI process has a different set of random numbers so we need to broadcast
#ifdef _MPIF90
              call mpi_bcast(params%box_scaling_factor, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
           end if
           !       This is the so-called total enthalpy Hamiltonian Montecarlo approach (with physical masses)
           !       We do not need to broadcast the velocities here since they get broadcasted later on; otherwise
           !       we would have to do it since each MPI rank may see a different random number
           call random_number( velocities )
           call remove_cm_vel(velocities(1:3,1:n_sites), masses(1:n_sites))
           e_kin = 0.d0
           do i = 1, n_sites
              e_kin = e_kin + 0.5d0 * masses(i) * dot_product(velocities(1:3, i), velocities(1:3, i))
           end do
           call random_number( rand )
           !        rand = rand * 4.d0/3.d0 - 1.d0/3.d0
           !        velocities = velocities / sqrt(e_kin) * sqrt(e_max - energy - params%p_nested/eVperA3tobar*v_uc + &
           !                                                     1.5d0*real(n_sites-1)*kB*params%t_extra*max(0.d0, rand))
           velocities = velocities / sqrt(e_kin) * sqrt(rand*(e_max - energy - params%p_nested/eVperA3tobar*v_uc))
        else if( i_nested == params%n_nested )then
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
     IF( rank == 0 )THEN
#endif

        if(params%do_mc)then
           if (mc_istep == params%mc_nsteps) then
              exit_loop = .true.
           else
              exit_loop = .false.
           end if


           if ( .not. exit_loop .and. ( &
                (md_istep == -1) .or. &
                ( params%do_md .and. ( &
                   (md_istep == params%md_nsteps)  .or. &
                   ( (abs(energy-energy_prev) < params%e_tol*dfloat(n_sites)) .and. (maxval(forces) < params%f_tol) ) &
                   ) ) ) ) then
                 !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
                 !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
                 !       >> First generate a random number in the range of the number of


              call cpu_time(time_mc(1))


              !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
              !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
              !       >> First generate a random number in the range of the number of

              if (mc_istep > 0)then
                 !       Evaluate the conditions for acceptance
                 !       > We have the mc conditions in mc.f90
                 !       > We care about comparing e_store to the energy of the new configuration based on the mc_movw

                 ! Reset the parameters for md / relaxation
                 if (params%do_md)then
                    md_istep = -1
                    params%do_md = .false.
                    do_mc_relax = .false.
                    ! Assume that the number of steps has already been set.
                 end if


                 if (.not. params%mc_hamiltonian) E_kinetic = 0.d0

                 call from_properties_to_image(images(i_trial_image), positions, velocities, masses, &
                      forces, a_box, b_box, c_box,  energy, energies, energy_exp, E_kinetic, &
                      species, species_supercell, n_sites, indices, fix_atom, &
                      xyz_species, xyz_species_supercell, local_properties)

                 if (params%verb > 50) write(*,*)'.......................................|'
                 if (params%verb > 50) write(*,'(A,1X,I0)')   ' MC Iteration:', mc_istep
                 if (params%verb > 50) write(*,'(A,1X,A)')    '    Move type:', mc_move
                 if (params%verb > 50) write(*,'(A,1X,F22.8)')'   &
                      & Etot_prev:', images(i_current_image)&
                      &%energy + images(i_current_image)%e_kin
                 if (params%verb > 50) write(*,'(A,1X,F22.8)')'   &
                      & Etot_new :', images(i_trial_image)%energy &
                      &+ images(i_trial_image)%e_kin

                 v_uc = dot_product( cross_product(a_box, b_box), c_box ) / (dfloat(indices(1)*indices(2)*indices(3)))

                 if (params%accessible_volume)then
                    call get_accessible_volume(v_uc, v_a_uc, species, params%radii)
                    if (params%verb > 50) write(*,'(A,F12.6,A,F12.6&
                         &,1X,A)') ' V_acc new: ', v_a_uc, ' A^3&
                         & V_acc old ', v_a_uc_prev, 'A^3 |'
                 else
                    v_a_uc = v_uc
                 end if


                 call get_mc_acceptance(mc_move, p_accept, &
                      energy + E_kinetic, &
                      images(i_current_image)%energy + images(i_current_image)%e_kin, &
                      params%t_beg, mc_id, mc_mu_id, &
                      params%mc_mu, n_mc_species, v_uc, v_uc_prev,&
                      & v_a_uc, v_a_uc_prev, params&
                      &%masses_types, params%p_beg)

!
!                 call get_mc_acceptance(mc_move, p_accept, &
!                      energy + E_kinetic, &
!                      images(i_current_image)%energy + images(i_current_image)%e_kin, &
!                      params%t_beg, &
!                      params%mc_mu(mc_mu_id), n_mc_species(mc_mu_id), v_uc, v_uc_prev,&
!                      & v_a_uc, v_a_uc_prev, params&
!                      &%masses_types(mc_id(mc_mu_id)), params%p_beg)


                 call random_number(ranf)

                 if (mc_move == "insertion") n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) +1
                 if (mc_move == "removal"  ) n_mc_species(mc_mu_id) = n_mc_species(mc_mu_id) -1

                 !    ACCEPT OR REJECT
                 if (params%verb > 50) write(*, '(A,1X,A,1X,A,L4,1X&
                      &,A,ES12.6,1X,A,1X,ES12.6)') 'Is ',&
                      & trim(mc_move), 'accepted?', p_accept >&
                      & ranf, ' p_accept =', p_accept, ' ranf = ',&
                      & ranf

                 if ( mc_istep == 1 )then
                    open(unit=200, file="mc.log", status="unknown")
                    if (energy_exp > 0.d0)then
                       write(200, '(A)') '# mc_istep  mc_move &
                            & accepted  E_trial              E_current             E_exp_trial&
                            &          E_exp_current  N_tot_trial &
                            & N_mc_species_trial'
                    else
                       write(200, '(A)') '# mc_istep  mc_move &
                            & accepted  E_trial              E_current &
                            &          N_tot_trial  N_mc_species_trial'
                    end if

                 end if
                 if ( mc_istep > 1  )then
                    open(unit=200, file="mc.log", status="old", position="append")
                 end if

                 ! collect the strings for the species etc
                 temp_string = ""
                 temp_string2 = ""

                 do i = 1, params%n_mc_mu
                    temp_string=""
                    write(temp_string, "(A,1X,I8)")  trim(params%mc_species(i)), n_mc_species(i)
                    temp_string2 = trim(temp_string2) // " " // trim(temp_string)
                 end do

                 if (energy_exp > 0.d0 )then

                    write(200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                         mc_istep, trim(adjustl(mc_move)), p_accept > ranf, energy + E_kinetic, &
                         images(i_current_image)%energy +&
                         & images(i_current_image)%e_kin, energy_exp,&
                         & images(i_current_image)%energy_exp,&
                         & images(i_trial_image)%n_sites,&
                         & trim(temp_string2)
                 else
                    write(200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                         mc_istep, trim(adjustl(mc_move)), p_accept > ranf, energy + E_kinetic, &
                         images(i_current_image)%energy + images(i_current_image)%e_kin, &
                         images(i_trial_image)%n_sites,  trim(temp_string2)

                 end if


                 if (mc_istep >= 1 ) close(200)


                 if (p_accept > ranf)then
                    !             Accept
                    ! Set variables
                    n_sites_prev = n_sites
                    v_uc_prev = v_uc
                    v_a_uc_prev = v_a_uc
                    virial_prev = virial
                    !   Assigning the default image with the accepted one
                    images(i_current_image) = images(i_trial_image)

                    if (params%n_mc_mu > 0)then
                       n_mc_species_prev = n_mc_species
                    end if


                 end if
                 instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic
                 instant_pressure = (kB*dfloat(n_sites-1)*instant_temp&
                      &+(virial(1,1) + virial(2,2) + virial(3,3))/3.d0)&
                      &/v_uc*eVperA3tobar


                 if ((params%mc_write_xyz .or. mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                      modulo(mc_istep, params%write_xyz) == 0))then
                    if (params%verb > 50) write(*,'(1X,A)')'&
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
                         & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                         &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                         & params%valid_pdf, params%valid_sf,&
                         & params%valid_xrd, params%valid_nd,&
                         & params%do_pair_distribution, params&
                         &%do_structure_factor, params%do_xrd,&
                         & params%do_nd, string)

                    call write_extxyz( images(i_current_image)%n_sites, 0, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                         images(i_current_image)%a_box/dfloat(indices(1)), &
                         images(i_current_image)%b_box/dfloat(indices(2)), &
                         images(i_current_image)%c_box/dfloat(indices(3)), &
                         virial_prev, images(i_current_image)%xyz_species, &
                         images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites),&
                         images(i_current_image)%velocities, &
                         images(i_current_image)%forces, &
                         images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                         images(i_current_image)%masses,  &
                         params%write_property, params&
                         &%write_array_property, params&
                         &%write_local_properties,&
                         & local_property_labels,&
                         & images(i_current_image)%local_properties&
                         &,images(i_current_image)%fix_atom,&
                         & "mc_current.xyz", string, .true. )

                    call write_extxyz( images(i_current_image)%n_sites, 1, 1.0d0, 0.d0, instant_temp, instant_pressure, &
                         images(i_current_image)%a_box/dfloat(indices(1)), &
                         images(i_current_image)%b_box/dfloat(indices(2)), &
                         images(i_current_image)%c_box/dfloat(indices(3)), &
                         virial_prev, images(i_current_image)%xyz_species, &
                         images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites),&
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
                         & "mc_all.xyz", string, .false. )

                 end if

                 !          Add acceptance to the log file else dont
                 call cpu_time(time_mc(2))
                 time_mc(3) = time_mc(3) + time_mc(2) - time_mc(1)


              else ! if (mc_istep == 0)
                 temp_md_nsteps = params%md_nsteps
                 if (params%verb > 50) write(*,*) '                                       |'
                 if (params%verb > 50) write(*,*) 'Starting MC, using parameters:         |'
                 if (params%verb > 50) write(*,*) '                                       |'
                 if (params%verb > 50) write(*,'(1X,A,1X,I8,1X,A)')  &
                      &  'mc_nsteps     = ', params%mc_nsteps, '     &
                      &        |'
                 if (params%verb > 50) write(*,'(1X,A,1X,I8,1X,A)')  &
                      &  'n_mc_types    = ', params%n_mc_types, '    &
                      &         |'
                 if (params%verb > 50) write(*,'(1X,A)') 'mc_types:                              |'
                 do i = 1, params%n_mc_types
                    if (params%verb > 50) write(*,'(1X,A,1X,A,1X,A)')&
                         & '     ', params%mc_types(i), '|'
                 end do
                 if (params%verb > 50) write(*,'(1X,A)') 'mc_accept_ratio:                       |'
                 do i = 1, params%n_mc_types
                    if (params%verb > 50) write(*,'(1X,A,1X,F12.8,1X&
                         &,A)') '   ', params%mc_acceptance(i), '    &
                         &                  |'
                 end do
                 if (params%verb > 50) write(*,'(1X,A,1X,I8,1X,A)')  &
                      &  'n_mc_swaps    = ', params%n_mc_swaps, '    &
                      &         |'
                 if (params%verb > 50) write(*,'(1X,A)') 'mc_swaps:  &
                      &                            |'
                 do i = 1, 2*params%n_mc_swaps
                    if (params%verb > 50) write(*,'(1X,A,1X,A,1X,A)')&
                         & '   ', params%mc_swaps(i), '              &
                         &        |'
                 end do
                 if (params%verb > 50) write(*,'(1X,A,1X,F17.8,1X&
                      &,A)') 'mc_move_max   = ', params%mc_move_max, &
                      & 'A   |'


                 do i = 1, params%n_mc_mu
                    write(*,'(1X,A,1X,F17.8,1X,A)') 'mc_mu         = ', params%mc_mu(1),        'eV  |'
                    write(*,'(1X,A,1X,A,1X,A)')     'mc_species    = ', trim(params%mc_species(i)),   '                    |'
                 end do

                 if (params%verb > 50) write(*,'(1X,A,1X,F17.8,1X&
                      &,A)') 'mc_min_dist   = ', params%mc_min_dist, &
                      & 'A   |'
                 if (params%verb > 50) write(*,'(1X,A,1X,F17.8,1X&
                      &,A)') 'mc_lnvol_max  = ', params%mc_lnvol_max,&
                      & '    |'
                 if (params%verb > 50) write(*,'(1X,A,1X,L8,1X,A)')  &
                      &  'mc_write_xyz  = ', params%mc_write_xyz, '  &
                      &           |'
                 if (params%verb > 50) write(*,'(1X,A,1X,L8,1X,A)')  &
                      &  'mc_relax      = ', params%mc_relax,     '  &
                      &           |'
                 if (params%verb > 50) write(*,'(1X,A,1X,I8,1X,A)')  &
                      &  'mc_nrelax     = ', params%mc_nrelax,    '  &
                      &           |'
                 if (params%verb > 50) write(*,'(1X,A,1X,A,1X,A)')   &
                      &  'mc_relax_opt  = ', params%mc_relax_opt, '  &
                      &   |'
                 if (params%verb > 50) write(*,'(1X,A,1X,A,1X,A)')   &
                      &  'mc_hybrid_opt = ', params%mc_hybrid_opt,'  &
                      &   |'
                 if (params%verb > 50) write(*,'(1X,A,1X,L8,1X,A)')  &
                      &  'mc_optimize_exp = ', params%mc_optimize_exp&
                      &,'  |'
                 if (params%verb > 50) write(*,'(1X,A,1X,L8,1X,A)')  &
                      &  'mc_hamiltonian = ', params%mc_hamiltonian,'&
                      &  |'
                 if (params%verb > 50) write(*,'(1X,A,1X,L8,1X,A)')  &
                      &   'mc_reverse    = ', params%mc_reverse,'    &
                      & |'
                 if (params%verb > 50) write(*,'(1X,A,1X,F12.6,1X&
                      &,A)') 'mc_reverse_lambda = ', params&
                      &%mc_reverse_lambda,'|'

                 if (params%verb > 50) write(*,*) '                                       |'
                 ! t_beg must

                 if( .not. allocated( images ) .and. .not. params%do_nested_sampling )then
                    allocate( images(1:2) )
                 else if (.not. allocated(images) .and. params%do_nested_sampling )then
                    allocate(images(1:2*i_image))
                 end if


                 if( .not. allocated(mc_id) .and. params%n_mc_mu > 0) then
                    allocate(mc_id(1:params%n_mc_mu))
                    allocate(n_mc_species(1:params%n_mc_mu))
                    allocate(n_mc_species_prev(1:params%n_mc_mu))

                    mc_id = 1
                    n_mc_species = 0

                    !    get the mc species types

                    do j = 1, params%n_mc_mu
                       do i = 1, n_species
                          if (params%species_types(i) == params%mc_species(j) )then
                             mc_id(j)=i
                          end if
                       end do
                    end do
                 end if

                 !       Now use the image construct to store this as the image to compare to
                 call from_properties_to_image(images(i_current_image), positions, velocities, masses, &
                      forces, a_box, b_box, c_box,  energy, energies, energy_exp, E_kinetic, &
                      species, species_supercell, n_sites, indices, fix_atom, &
                      xyz_species, xyz_species_supercell, local_properties)


                 instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic
                 instant_pressure = (kB*dfloat(n_sites-1)*instant_temp&
                      &+(virial(1,1) + virial(2,2) + virial(3,3))/3.d0)&
                      &/v_uc*eVperA3tobar

                 if ((mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                      modulo(mc_istep, params%write_xyz) == 0))then
                    if (params%verb > 50) write(*,'(1X,A)')' Writing mc_current.xyz and mc_all.xyz '
                    call wrap_pbc(images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites), &
                         images(i_current_image)%a_box/dfloat(indices(1)), &
                         images(i_current_image)%b_box/dfloat(indices(2)),&
                         images(i_current_image)%c_box/dfloat(indices(3)))
                    call get_xyz_energy_string(energies_soap, energies_2b,&
                         & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                         &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                         & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                         & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                    call write_extxyz( images(i_current_image)%n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                         images(i_current_image)%a_box/dfloat(indices(1)), &
                         images(i_current_image)%b_box/dfloat(indices(2)), &
                         images(i_current_image)%c_box/dfloat(indices(3)), &
                         virial_prev, images(i_current_image)%xyz_species, &
                         images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites),&
                         images(i_current_image)%velocities, &
                         images(i_current_image)%forces, &
                         images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                         images(i_current_image)%masses,  &
                         params%write_property, params&
                         &%write_array_property, params&
                         &%write_local_properties,&
                         & local_property_labels, images(i_current_image)%local_properties&
                         &, images(i_current_image)%fix_atom,&
                         & "mc_current.xyz", string, .true. )

                    call write_extxyz( images(i_current_image)%n_sites, 1, 1.0d0, 0.0d0,  instant_temp, instant_pressure, &
                         images(i_current_image)%a_box/dfloat(indices(1)), &
                         images(i_current_image)%b_box/dfloat(indices(2)), &
                         images(i_current_image)%c_box/dfloat(indices(3)), &
                         virial_prev, images(i_current_image)%xyz_species, &
                         images(i_current_image)%positions(1:3, 1:images(i_current_image)%n_sites),&
                         images(i_current_image)%velocities, &
                         images(i_current_image)%forces, &
                         images(i_current_image)%energies(1:images(i_current_image)%n_sites), &
                         images(i_current_image)%masses, &
                         params%write_property, params&
                         &%write_array_property, params&
                         &%write_local_properties,&
                         & local_property_labels, images(i_current_image)%local_properties&
                         &, images(i_current_image)%fix_atom,&
                         & "mc_all.xyz", string, .true. )


                    v_uc_prev = dot_product( cross_product(a_box, b_box), c_box ) / (dfloat(indices(1)*indices(2)*indices(3)))
                    if (params%accessible_volume)then
                       call get_accessible_volume(v_uc_prev, v_a_uc_prev, species, params%radii)
                    else
                       v_a_uc_prev = v_uc_prev
                    end if
                 end if


              end if

              !  Now start the mc logic: first, use the stored images properties
              call from_image_to_properties(images(i_current_image), positions, velocities, masses, &
                   forces, a_box, b_box, c_box, energy, energies, energy_exp,  E_kinetic, &
                   species, species_supercell, n_sites, indices, fix_atom, &
                   xyz_species, xyz_species_supercell, local_properties)


                 call perform_mc_step(&
                      & positions, species, xyz_species, masses, fix_atom,&
                      & velocities, positions_prev, positions_diff, disp, d_disp, params%n_local_properties,&
                      & params%mc_acceptance, params%mc_mu_acceptance, local_properties, &
                      images(i_current_image)%local_properties, energies,&
                      & forces, forces_prev, n_sites, params%n_mc_mu, mc_mu_id, n_mc_species,&
                      & mc_move, params %mc_species,&
                      & params%mc_move_max, params%mc_min_dist, params%mc_max_insertion_trials, &
                      params%mc_lnvol_max, params%mc_types, params%masses_types, species_idx,&
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
                      &%mc_relax_after, do_mc_relax, params%verb, &
                      params%mc_n_planes, params%mc_planes, params%mc_max_dist_to_planes, params%mc_planes_restrict_to_polyhedron)


              rebuild_neighbors_list = .true.
              ! end if

              ! NOTE: the species_supercell and xyz_species_supercell are
              ! not commensurate with the new image as these have not been
              ! calculated. If reading from an outputted xyz file, then it
              ! should be okay but really the new atoms should be added to
              ! the supercell in the usual way, but for convenience, one has
              ! not done that.


              if(params%mc_relax .and. do_mc_relax)then
                 ! Set the parameters for relaxatrino
                 md_istep = -1
                 params%do_md = .true.
                 params%optimize = params%mc_relax_opt
                 params%md_nsteps = params%mc_nrelax
                 ! Note, that this may override md steps if the same is chosen! More testing needed
              end if
              ! If doing md, don't relax
              if( mc_move == 'md')then
                 ! Set the parameters for relaxatrino
                 md_istep = -1
                 params%do_md = .true.
                 params%optimize = params%mc_hybrid_opt
                 params%md_nsteps = temp_md_nsteps
                 ! Note, that this may override md steps if the same is chosen! More testing needed
              end if


              if ((params%mc_write_xyz .or. mc_istep == 0 .or. mc_istep == params%mc_nsteps .or. &
                   modulo(mc_istep, params%write_xyz) == 0))then

                 call wrap_pbc(positions(1:3,1:n_sites), &
                      a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
                 call get_xyz_energy_string(energies_soap, energies_2b,&
                      & energies_3b, energies_core_pot, energies_vdw, energies_exp&
                      &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
                      & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                      & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                 call write_extxyz( n_sites, 0, 1.0d0, 0.0d0, instant_temp, instant_pressure, &
                      a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                      virial, xyz_species, &
                      positions(1:3, 1:n_sites), velocities, &
                      forces, energies(1:n_sites), masses, &
                      params%write_property, params&
                      &%write_array_property, params&
                      &%write_local_properties,&
                      & local_property_labels, local_properties&
                      &,fix_atom, mc_file, string,  .true. )
              end if
              ! As we have moved/added/removed, we must check the supercell and  broadcast the results

              call read_xyz(mc_file, .true., params%all_atoms, params%do_timing, &
                   n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                   positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                   xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                   n_sites, .true., fix_atom, params%t_beg, &
                   params%write_array_property(6), .true. )

           else
              if( mc_move == 'md')then
                 if( params%print_progress .and. md_istep == 0 )then
                    write(*,*)'                                       |'
                    write(*,*)'Progress:                              |'
                    write(*,*)'                                       |'
                    write(*,'(1X,A)',advance='no')'[                                    ] |'
                    update_bar = params%md_nsteps/36
                    if( update_bar < 1 )then
                       update_bar = 1
                    end if
                    counter = 1
                 else if( md_istep == params%md_nsteps-1 .or. &
                      (abs(energy-energy_prev) < params%e_tol*dfloat(n_sites) .and. &
                      maxval(forces) < params%f_tol) .and. md_istep > 0 )then
                    write(*,*)
                 else if( params%print_progress .and. counter == update_bar .and. md_istep < params%md_nsteps-1 )then
                    do j = 1, 36+3
                       write(*,"(A)", advance="no") creturn
                    end do
                    write (*,"(1X,A)",advance="no") "["
                    do i = 1, 36*(md_istep+1)/params%md_nsteps
                       write (*,"(A)",advance="no") "."
                    end do
                    do i = 36*(md_istep+1)/params%md_nsteps+1, 36
                       write (*,"(A)",advance="no") " "
                    end do
                    write (*,"(A)",advance="no") "] |"
                    counter = 1
                 else
                    counter = counter + 1
                 end if


                 if (params%mc_hamiltonian)then
                    if (params%verb > 50) write(*,'(1X,A,1X,F20.8,1X&
                         &,A,1X,I8,1X,A,1X,I8)')"Hybrid md step: H =&
                         & T + V = ", energy + E_kinetic, ",&
                         & iteration ", md_istep, "/", params&
                         &%md_nsteps
                 else
                    if (params%verb > 50) write(*,'(1X,A,1X,F20.8,1X&
                         &,A,1X,I8,1X,A,1X,I8)')"Hybrid md step:&
                         & energy = ", energy , ", iteration ",&
                         & md_istep, "/", params%md_nsteps
                 end if

                 if (params%verb > 50) write(*,'(A,1X,F22.8,1X,A)')' SOAP energy:', sum(energies_soap), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F24.8,1X,A)')' 2b energy:', sum(energies_2b), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F24.8,1X,A)')' 3b energy:', sum(energies_3b), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F18.8,1X,A)')'&
                      & core_pot energy:', sum(energies_core_pot),&
                      & 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F23.8,1X,A)')'&
                      & vdw energy:', sum(energies_vdw), 'eV |'
                 if (params%verb > 50 .and. valid_xps) write(*,'(A,1X,F23.8,1X,A)')' xps energy:', sum(energies_lp), 'eV |'

                 if ( params%valid_pdf .and. params%do_pair_distribution .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' pdf energy:',&
                      & sum(energies_pdf), 'eV |'
                 if ( params%valid_sf .and. params%do_structure_factor .and. params%verb > 50)&
                      & write(*,'(A,1X,F24.8,1X,A)')' sf energy:',&
                      & sum(energies_sf), 'eV |'
                 if ( params%valid_xrd .and. params%do_xrd .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' xrd energy:',&
                      & sum(energies_xrd), 'eV |'
                 if ( params%valid_nd .and. params%do_nd .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' nd energy:',&
                      & sum(energies_nd), 'eV |'



              else
                 if( params%print_progress .and. mc_istep == 0 )then
                    write(*,*)'                                       |'
                    write(*,*)'Progress:                              |'
                    write(*,*)'                                       |'
                    write(*,'(1X,A)',advance='no')'[                                    ] |'
                    update_bar = params%mc_nsteps/36
                    if( update_bar < 1 )then
                       update_bar = 1
                    end if
                    counter = 1
                 else if( mc_istep == params%mc_nsteps-1 .and. mc_istep > 0 )then
                    write(*,*)
                 else if( params%print_progress .and. counter == update_bar .and. mc_istep < params%mc_nsteps-1 )then
                    do j = 1, 36+3
                       write(*,"(A)", advance="no") creturn
                    end do
                    write (*,"(1X,A)",advance="no") "["
                    do i = 1, 36*(mc_istep+1)/params%mc_nsteps
                       write (*,"(A)",advance="no") "."
                    end do
                    do i = 36*(mc_istep+1)/params%mc_nsteps+1, 36
                       write (*,"(A)",advance="no") " "
                    end do
                    write (*,"(A)",advance="no") "] |"
                    counter = 1
                 else
                    counter = counter + 1
                 end if


                 if (params%verb > 50 .and. do_mc_relax) write(*,'(1X,A,1X,F20.8,1X,A&
                      &,1X,I8,1X,A,1X,I8)')"MC Relax md step: energy &
                      &= ", energy, ", iteration ", md_istep, "/",&
                      & params%mc_nrelax
                 if (params%verb > 50) write(*,'(A,1X,F22.8,1X,A)')' SOAP energy:', sum(energies_soap), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F24.8,1X,A)')' 2b energy:', sum(energies_2b), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F24.8,1X,A)')' 3b energy:', sum(energies_3b), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F18.8,1X,A)')' core_pot energy:', sum(energies_core_pot), 'eV |'
                 if (params%verb > 50) write(*,'(A,1X,F23.8,1X,A)')' vdw energy:', sum(energies_vdw), 'eV |'
                 if (params%verb > 50 .and. valid_xps) write(*,'(A,1X,F23.8,1X,A)')' xps energy:', sum(energies_lp), 'eV |'

                 if ( params%valid_pdf .and. params%do_pair_distribution .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' pdf energy:',&
                      & sum(energies_pdf), 'eV |'
                 if ( params%valid_sf .and. params%do_structure_factor .and. params%verb > 50)&
                      & write(*,'(A,1X,F24.8,1X,A)')' sf energy:',&
                      & sum(energies_sf), 'eV |'
                 if ( params%valid_xrd .and. params%do_xrd .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' xrd energy:',&
                      & sum(energies_xrd), 'eV |'
                 if ( params%valid_nd .and. params%do_nd .and. params%verb > 50)&
                      & write(*,'(A,1X,F23.8,1X,A)')' nd energy:',&
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
     IF( params%do_mc .and. md_istep == -1 .and. rank == 0 )THEN
        n_pos = size(positions,2)
        n_sp = size(xyz_species,1)
        n_sp_sc = size(xyz_species_supercell,1)
     END IF
     call cpu_time(time_mpi(1))
     call mpi_bcast(n_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(n_sp_sc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(params%do_md, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(md_istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call cpu_time(time_mpi(2))
     time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
     IF( rank /= 0 )THEN !.and. (mc_move == "insertion" .or. mc_move == "removal")
        if(allocated(positions))deallocate(positions)
        allocate( positions(1:3, n_pos) )
        if( params%do_md .or. params%do_nested_sampling  .or. params%do_mc)then
           if(allocated(velocities))deallocate(velocities)
           allocate( velocities(1:3, n_pos) )
           if(allocated(masses))deallocate(masses)
           allocate( masses(1:n_sp) )
           if(allocated(fix_atom))deallocate(fix_atom)
           allocate( fix_atom(1:3, 1:n_sp) )

           ! if(allocated(forces_prev))deallocate(forces_prev)
           ! allocate( forces_prev(1:3, 1:n_sites) )
           ! if(allocated(positions_prev))deallocate(positions_prev)
           ! allocate( positions_prev(1:3, 1:n_sites) )
           ! if(allocated(positions_diff))deallocate(positions_diff)
           ! allocate( positions_diff(1:3, 1:n_sites) )
           ! positions_diff = 0.d0

        end if
        if(allocated(xyz_species))deallocate(xyz_species)
        allocate( xyz_species(1:n_sp) )
        if(allocated(species))deallocate(species)
        allocate( species(1:n_sp) )
        if(allocated(xyz_species_supercell))deallocate(xyz_species_supercell)
        allocate( xyz_species_supercell(1:n_sp_sc) )
        if(allocated(species_supercell))deallocate(species_supercell)
        allocate( species_supercell(1:n_sp_sc) )
     END IF
     call cpu_time(time_mpi_positions(1))
     call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     if( params%do_md .or. params%do_nested_sampling .or. params%do_mc )then
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
     call cpu_time(time_mpi_positions(2))
     time_mpi_positions(3) = time_mpi_positions(3) + time_mpi_positions(2) - time_mpi_positions(1)
#endif
     !   Now that all ranks know the size of n_sites, we allocate do_list
     if( .not. params%do_md .or. (params%do_md .and. md_istep == 0) .or. &
          ( params%do_mc ) )then
        if( allocated(do_list))deallocate(do_list)
        allocate( do_list(1:n_sites) )
        do_list = .true.
     end if
     !
     call cpu_time(time1)
#ifdef _MPIF90
     !   Parallel neighbors list build
     call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif


     if( rebuild_neighbors_list )then
        deallocate( rjs, xyz, thetas, phis, neighbor_species )
        deallocate( neighbors_list, n_neigh )
#ifdef _MPIF90
        deallocate( n_neigh_local )
#endif
     end if
     if( ( params%do_nested_sampling .and. .not. params%do_mc) .and. &
          (params%do_md .and. (md_istep == params%md_nsteps .or. exit_loop)))then
        deallocate( positions, xyz_species, xyz_species_supercell, species, species_supercell, do_list )
        if( allocated(velocities) )deallocate( velocities )
     end if
     if(params%do_mc .and. params%do_md)then
        if (params%do_mc .and. (mc_istep == params%mc_nsteps .or. exit_loop))then
           deallocate( positions, xyz_species, xyz_species_supercell, species, species_supercell, do_list )
           if( allocated(velocities) )deallocate( velocities )
        end if
     end if

     if( (params%do_md .and. .not. params%do_mc) .and. &
          (md_istep == params%md_nsteps .or. exit_loop) .and. rank == 0 )then
        deallocate( positions_prev, forces_prev )
     end if
     if( params%do_mc .and. (mc_istep == params%mc_nsteps .or. exit_loop) .and. rank == 0 )then
        if( allocated(forces_prev)) deallocate(forces_prev)
        if( allocated(positions_prev)) deallocate(positions_prev)
     end if

     if( params%exp_forces .and. (md_istep == params%md_nsteps .or.&
          & mc_istep == params%mc_nsteps .or. exit_loop) )then
        do i = 1, params%n_exp
           if( allocated(params%exp_data(i)%x)) deallocate(params%exp_data(i)%x)
           if( allocated(params%exp_data(i)%y)) deallocate(params%exp_data(i)%y)
           if( allocated(params%exp_data(i)%y_pred)) deallocate(params%exp_data(i)%y_pred)
        end do
     end if


     if (.not. params%do_mc )n_sites_prev = n_sites
     n_atom_pairs_by_rank_prev = n_atom_pairs_by_rank(rank+1)

#ifdef _MPIF90
     call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     if( exit_loop )exit
     ! End of loop through structures in the xyz file or MD steps
  end do



  if( params%do_md .or. params%do_prediction .or. params%do_mc)then
     call cpu_time(time2)
#ifdef _MPIF90
     IF( rank == 0 )then
#endif
        if( params%do_md .and. .not. params%do_nested_sampling )then
           !      write(*,'(A)')'] |'
           !      write(*,*)
           write(*,*)'                                       |'
           !      write(*,'(I8,A,F13.3,A)') params%md_nsteps, ' MD steps:', time2-time3, ' seconds |'
           write(*,'(I8,A,F13.3,A)') md_istep, ' MD steps:', time2-time3, ' seconds |'
        end if
        if( params%do_mc )then
           !      write(*,'(A)')'] |'
           write(*,*)
           write(*,*)'                                       |'
           !      write(*,'(I8,A,F13.3,A)') params%md_nsteps, ' MD steps:', time2-time3, ' seconds |'
           write(*,'(I8,A,F13.3,A)') mc_istep, ' MC steps:', time2-time3, ' seconds |'
        end if

        write(*,*)'                                       |'
        write(*,'(A,F13.3,A)') ' *     Read input:', time_read_input(3), ' seconds |'
        write(*,'(A,F13.3,A)') ' * Read XYZ files:', time_read_xyz(3), ' seconds |'
        write(*,'(A,F13.3,A)') ' * Neighbor lists:', time_neigh, ' seconds |'
        write(*,'(A,F13.3,A)') ' *  GAP desc/pred:', time_gap, ' seconds |'
        write(*,'(A,F13.3,A)') '     - soap_turbo:', time_soap(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -         2b:', time_2b(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -         3b:', time_3b(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -   core_pot:', time_core_pot(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -        vdw:', time_vdw(3), ' seconds |'
        if (valid_xps .or. params%do_pair_distribution .or. params&
             &%do_structure_factor .or. params%do_xrd .or. params%do_nd) write(*,'(A&
             &,F13.3,A)')      ' *  Exp. pred.   :', time_pdf(3) + time_sf(3) + time_xrd(3) + time_nd(3), ' seconds&
             & |'
        if( valid_xps ) write(*,'(A,F13.3,A)') '     -        xps:',&
             & time_xps(3), ' seconds |'
        if( params%do_pair_distribution ) write(*,'(A,F13.3,A)') '     -        pdf:', time_pdf(3), ' seconds |'
        if( params%do_structure_factor  ) write(*,'(A,F13.3,A)') '     -         sf:', time_sf(3), ' seconds |'
        if( params%do_xrd  )              write(*,'(A,F13.3,A)') '     -        xrd:', time_xrd(3), ' seconds |'
        if( params%do_nd  )               write(*,'(A,F13.3,A)') '     -         nd:', time_nd(3), ' seconds |'

        if( params%do_md )then
           write(*,'(A,F13.3,A)') ' *  MD algorithms:', time_md(3), ' seconds |'
        end if
        if( params%do_mc )then
           write(*,'(A,F13.3,A)') ' *  MC algorithms:', time_mc(3), ' seconds |'
        end if

#ifdef _MPIF90
        write(*,'(A,F13.3,A)') ' *  MPI comms.   :', time_mpi(3) + time_mpi_positions(3) + time_mpi_ef(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -  pos & vel:', time_mpi_positions(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     - E & F brc.:', time_mpi_ef(3), ' seconds |'
        write(*,'(A,F13.3,A)') '     -  MPI misc.:', time_mpi(3), ' seconds |'
        write(*,'(A,F13.3,A)') ' *  Miscellaneous:', time2-time3 - time_neigh - time_gap - time_read_input(3) &
             - time_read_xyz(3) - time_mpi(3) - time_mpi_positions(3)&
             & - time_mpi_ef(3) - time_md(3) - time_xps(3) -&
             & time_pdf(3) - time_sf(3) - time_xrd(3) - time_nd(3), ' seconds |'
#else
        write(*,'(A,F13.3,A)') ' *  Miscellaneous:', time2-time3 - time_neigh - time_gap - time_read_input(3) &
             - time_read_xyz(3) - time_md(3)- time_xps(3) -&
             & time_pdf(3) - time_sf(3) - time_xrd(3)  - time_nd(3), ' seconds |'
#endif
        write(*,*)'                                       |'
        write(*,'(A,F13.3,A)') ' *     Total time:', time2-time3, ' seconds |'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
#ifdef _MPIF90
     END IF
#endif
  end if

  if ( allocated( fix_atom ))    deallocate( fix_atom )
  if ( allocated( positions ))    deallocate( positions )
  if ( allocated( velocities ))    deallocate( velocities )
  if ( allocated( positions_diff ))    deallocate( positions_diff )

  if ( allocated( energies ))          deallocate( energies )
  if ( allocated( energies_soap) )     deallocate( energies_soap )
  if ( allocated( energies_2b) )       deallocate( energies_2b )
  if ( allocated( energies_3b) )       deallocate( energies_3b )
  if ( allocated( energies_core_pot) ) deallocate( energies_core_pot )
  if ( allocated( energies_vdw) )      deallocate( energies_vdw )
  if ( allocated( energies_exp) )      deallocate( energies_exp )
  if ( allocated( energies_lp) )       deallocate( energies_lp )
  if ( allocated( energies_pdf ) )     deallocate( energies_pdf)
  if ( allocated( energies_sf ) )      deallocate( energies_sf)
  if ( allocated( energies_xrd ) )     deallocate( energies_xrd)
  if ( allocated( energies_nd ) )      deallocate( energies_nd)

  if ( allocated( this_energies ))          deallocate( this_energies )
  if ( allocated( this_energies_vdw) )      deallocate( this_energies_vdw )
  if ( allocated( this_energies_lp) )       deallocate( this_energies_lp )
  if ( allocated( this_energies_pdf ) )     deallocate( this_energies_pdf)
  if ( allocated( this_energies_sf ) )      deallocate( this_energies_sf)
  if ( allocated( this_energies_xrd ) )     deallocate( this_energies_xrd)
  if ( allocated( this_energies_nd ) )      deallocate( this_energies_nd)



  if ( allocated( forces ))          deallocate( forces )
  if ( allocated( forces_soap) )     deallocate( forces_soap )
  if ( allocated( forces_2b) )       deallocate( forces_2b )
  if ( allocated( forces_3b) )       deallocate( forces_3b )
  if ( allocated( forces_core_pot) ) deallocate( forces_core_pot )
  if ( allocated( forces_vdw) )      deallocate( forces_vdw )
  if ( allocated( forces_lp) )       deallocate( forces_lp )
  if ( allocated( forces_pdf ) )     deallocate( forces_pdf)
  if ( allocated( forces_sf ) )      deallocate( forces_sf)
  if ( allocated( forces_xrd ) )     deallocate( forces_xrd)
  if ( allocated( forces_nd ) )      deallocate( forces_nd)

  if ( allocated( this_forces ))          deallocate( this_forces )
  if ( allocated( this_forces_vdw) )      deallocate( this_forces_vdw )
  if ( allocated( this_forces_lp) )       deallocate( this_forces_lp )
  if ( allocated( this_forces_pdf ) )     deallocate( this_forces_pdf)
  if ( allocated( this_forces_sf ) )      deallocate( this_forces_sf)
  if ( allocated( this_forces_xrd ) )     deallocate( this_forces_xrd)
  if ( allocated( this_forces_nd ) )      deallocate( this_forces_nd)


  if( allocated(local_properties) ) deallocate(local_properties)
  if( allocated(local_properties_cart_der) ) deallocate(local_properties_cart_der)
  if( allocated(this_local_properties) ) deallocate(this_local_properties)
  if( allocated(this_local_properties_cart_der) ) deallocate(this_local_properties_cart_der)


  if( allocated(soap_turbo_hypers) )deallocate(soap_turbo_hypers)
  if( allocated(distance_2b_hypers) )deallocate(distance_2b_hypers)
  if( allocated(angle_3b_hypers) )deallocate(angle_3b_hypers)
  if( allocated(core_pot_hypers) )deallocate(core_pot_hypers)

  deallocate( n_atom_pairs_by_rank )
  if( allocated( n_local_properties_mpi ) ) deallocate( n_local_properties_mpi )
  if( allocated( local_properties_n_sparse_mpi_soap_turbo ) ) deallocate( local_properties_n_sparse_mpi_soap_turbo )
  if( allocated( local_properties_dim_mpi_soap_turbo ) ) deallocate( local_properties_dim_mpi_soap_turbo )
  if( allocated( has_local_properties_mpi ) ) deallocate( has_local_properties_mpi )

  if (allocated(local_property_labels) ) deallocate( local_property_labels)
  if (allocated(local_property_indexes) ) deallocate( local_property_indexes)
  if( allocated(do_list))deallocate(do_list)
  if (allocated( params%write_local_properties )) deallocate(params%write_local_properties)

#ifdef _MPIF90
  IF( rank == 0 )then
#endif
     write(*,*)'                                       |'
     write(*,*)'End of execution                       |'
     write(*,*)'_______________________________________/'
#ifdef _MPIF90
  END IF
#endif



#ifdef _MPIF90
  call mpi_finalize(ierr)
#endif


end program turbogap
