! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap.f90, is copyright (c) 2019-2021, Miguel A. Caro
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
  use gap_interface
  use types
  use vdw
  use soap_turbo_functions
  use F_B_C
  use iso_c_binding
#ifdef _MPIF90
  use mpi
  use mpi_helper
#endif
  use bussi
  use xyz_module
  use F_B_C

!$  use omp_lib


  
  implicit none


!**************************************************************************
! Variable definitions
!
  real*8, allocatable :: rjs(:), thetas(:), phis(:), xyz(:,:), sph_temp(:), sph_temp3(:,:)
  real*8, allocatable :: positions(:,:), positions_prev(:,:), soap(:,:), soap_cart_der(:,:,:), &
                         positions_diff(:,:), forces_prev(:,:)
  real*8 :: rcut_max, a_box(1:3), b_box(1:3), c_box(1:3), max_displacement, energy, energy_prev
  real*8 :: virial(1:3, 1:3), this_virial(1:3, 1:3), virial_soap(1:3, 1:3), virial_2b(1:3, 1:3), &
            virial_3b(1:3,1:3), virial_core_pot(1:3, 1:3), virial_vdw(1:3, 1:3), &
            this_virial_vdw(1:3, 1:3), v_uc, eVperA3tobar = 1602176.6208d0
  real*8, allocatable :: energies(:), forces(:,:), energies_soap(:), forces_soap(:,:), this_energies(:), &
                         this_forces(:,:), &
                         energies_2b(:), forces_2b(:,:), energies_3b(:), forces_3b(:,:), &
                         energies_core_pot(:), forces_core_pot(:,:), &
                         velocities(:,:), masses_types(:), masses(:), hirshfeld_v(:), &
                         hirshfeld_v_cart_der(:,:)
  real*8, allocatable, target :: this_hirshfeld_v(:), this_hirshfeld_v_cart_der(:,:)
  real*8, pointer :: this_hirshfeld_v_pt(:), this_hirshfeld_v_cart_der_pt(:,:)
  real*8, allocatable :: all_energies(:,:), all_forces(:,:,:), all_virial(:,:,:)
  real*8, allocatable :: all_this_energies(:,:), all_this_forces(:,:,:), all_this_virial(:,:,:)
  real*8 :: instant_temp, kB = 8.6173303d-5, E_kinetic, time1, time2, time3, time_neigh, &
            time_gap, time_soap(1:3), time_2b(1:3), time_3b(1:3), time_read_input(1:3), time_read_xyz(1:3), &
            time_mpi(1:3) = 0.d0, time_core_pot(1:3), time_vdw(1:3), instant_pressure, lv(1:3,1:3), &
            time_mpi_positions(1:3) = 0.d0, time_mpi_ef(1:3) = 0.d0, time_md(3) = 0.d0, &
            instant_pressure_tensor(1:3, 1:3), time_step, md_time, &
            solo_time_soap=0.d0, soap_time_soap(3)=0.d0, time_get_soap=0.d0
  integer, allocatable :: displs(:), displs2(:), counts(:), counts2(:)
  integer :: update_bar, n_sparse
  logical, allocatable :: do_list(:), has_vdw_mpi(:), fix_atom(:,:)
  logical :: rebuild_neighbors_list = .true., exit_loop = .true.
  character*1 :: creturn = achar(13)

! Clean up these variables after code refactoring !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, allocatable :: n_neigh(:), neighbors_list(:), alpha_max(:), species(:), species_supercell(:), &
                          neighbor_species(:), sph_temp_int(:), der_neighbors(:), der_neighbors_list(:), &
                          i_beg_list(:), i_end_list(:), j_beg_list(:), j_end_list(:)
  integer :: n_sites, i, j, k, i2, j2, n_soap, k2, k3, l, n_sites_this, ierr, rank, ntasks, dim, n_sp, &
             n_pos, n_sp_sc, this_i_beg, this_i_end, this_j_beg, this_j_end, this_n_sites_mpi, n_sites_prev = 0, &
             n_atom_pairs_by_rank_prev
  integer :: l_max, n_atom_pairs, n_max, ijunk, central_species = 0, n_atom_pairs_total
  integer :: iostatus, counter = 0, counter2
  integer :: which_atom = 0, n_species = 1, n_xyz, indices(1:3)
  integer :: radial_enhancement = 0
  integer :: md_istep

  logical :: repeat_xyz = .true.

  character*1024 :: filename, cjunk, file_compress_soap, file_alphas, file_soap, file_2b, file_alphas_2b, &
                    file_3b, file_alphas_3b, file_gap = "none"
  character*64 :: keyword
  character*16 :: lattice_string(1:9)
  character*8 :: i_char
  character*8, allocatable :: species_types(:), xyz_species(:), xyz_species_supercell(:)
  character*1 :: keyword_first

! This is the mode in which we run TurboGAP
  character*16 :: mode = "none"

! Here we store the input parameters
  type(input_parameters) :: params

! These are the containers for the hyperparameters of descriptors and GAPs
  integer :: n_soap_turbo = 0, n_distance_2b = 0, n_angle_3b = 0, n_core_pot = 0
  type(soap_turbo), allocatable :: soap_turbo_hypers(:)
  type(distance_2b), allocatable :: distance_2b_hypers(:)
  type(angle_3b), allocatable :: angle_3b_hypers(:)
  type(core_pot), allocatable :: core_pot_hypers(:)

!vdw crap
  real*8, allocatable :: v_neigh_vdw(:), energies_vdw(:), forces_vdw(:,:), this_energies_vdw(:), this_forces_vdw(:,:)


! MPI stuff
  real*8, allocatable :: temp_1d(:), temp_1d_bis(:), temp_2d(:,:)
  integer, allocatable :: temp_1d_int(:), n_atom_pairs_by_rank(:), displ(:)
  integer :: i_beg, i_end, n_sites_mpi, j_beg, j_end, size_soap_turbo&
       &, size_distance_2b, size_angle_3b, n_omp, omp_task, omp_n_sites
  integer, allocatable :: n_species_mpi(:),&
       & n_sparse_mpi_soap_turbo(:), dim_mpi(:),&
       & n_sparse_mpi_distance_2b(:), n_sparse_mpi_angle_3b(:),&
       & n_mpi_core_pot(:), vdw_n_sparse_mpi_soap_turbo(:),&
       & n_neigh_local(:)
  logical, allocatable :: compress_soap_mpi(:)
  integer, allocatable :: i_beg_omp(:), i_end_omp(:), j_beg_omp(:), j_end_omp(:)
  type(c_ptr) :: cublas_handle, gpu_stream
!**************************************************************************

  call create_cublas_handle(cublas_handle, gpu_stream)
! integer :: n_ii,i_ii, j_jj,k_ii
! real*8 :: CC(1:1024,1:1024), BB(1:1024,1:1024), AA(1:1024,1:1024), dut1, dut2


! call random_number(AA)
! call random_number(BB)

!**************************************************************************
! Start recording the time
  call cpu_time(time1)
  time3 = time1
!**************************************************************************






!**************************************************************************
! MPI stuff
#ifdef _MPIF90
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  time1=MPI_Wtime()
  time3 = time1
!  allocate( displs(1:ntasks) )
!  allocate( displs2(1:ntasks) )
!  allocate( counts(1:ntasks) )
!  allocate( counts2(1:ntasks) )
  allocate( displ(1:ntasks) )
#else
  rank = 0
  ntasks = 1
#endif
  allocate( n_atom_pairs_by_rank(1:ntasks) )
!**************************************************************************



  call gpu_set_device(rank) 
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
  write(*,*)'                     Welcome to TurboGAP v0.1                    |'
  write(*,*)'                    Written and maintained by                    |'
  write(*,*)'                                                                 |'
  write(*,*)'                         Miguel A. Caro                          |'
  write(*,*)'                       mcaroba@gmail.com                         |'
  write(*,*)'                      miguel.caro@aalto.fi                       |'
  write(*,*)'                                                                 |'
  write(*,*)'       Department of Electrical Engineering and Automation       |'
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
  write(*,*)'*) Miguel A. Caro (Aalto University)                             |'
  write(*,*)'*) Patricia Hernández-León (Aalto University)                    |'
  write(*,*)'*) Suresh Kondati Natarajan (now @ QuantumWise, formerly @ Aalto)|'
  write(*,*)'*) Albert P. Bartók-Pártay (Warwick University)                  |'
  write(*,*)'*) Eelis V. Mielonen (now @ EPFL, formerly @ Aalto)              |'
  write(*,*)'*) Heikki Muhli (Aalto University)                               |'
  write(*,*)'*) Mikhail Kuklin (Aalto University)                             |'
  write(*,*)'*) Gábor Csányi (University of Cambridge)                        |'
  write(*,*)'*) Jan Kloppenburg (Aalto University)                            |'
  write(*,*)'                                                                 |'
  write(*,*)'.................................................................|'
  write(*,*)'                                                                 |'
  write(*,*)'                     Last updated: Jul. 2021                     |'
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
  !call cpu_time(time_read_input(1))
  time_read_input(1)=MPI_Wtime()
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
!   We increase rcut_max by the neighbors buffer
    rcut_max = rcut_max + params%neighbors_buffer
#ifdef _MPIF90
    END IF
#endif
!   THIS CHUNK HERE DISTRIBUTES THE INPUT DATA AMONG ALL THE PROCESSES
!   Broadcast number of descriptors to other processes
#ifdef _MPIF90
    !call cpu_time(time_mpi(1))
    time_mpi(1)=MPI_Wtime()
    call mpi_bcast(n_soap_turbo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_distance_2b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_angle_3b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_core_pot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!   Broadcast the maximum cutoff distance
    call mpi_bcast(rcut_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!   Processes other than 0 need to allocate the data structures on their own
    !call cpu_time(time_mpi(2))
    time_mpi(2)=MPI_Wtime()
    time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
    allocate( n_species_mpi(1:n_soap_turbo) )
    allocate( n_sparse_mpi_soap_turbo(1:n_soap_turbo) )
    allocate( dim_mpi(1:n_soap_turbo) )
    allocate( vdw_n_sparse_mpi_soap_turbo(1:n_soap_turbo) )
    allocate( has_vdw_mpi(1:n_soap_turbo) )
    allocate( compress_soap_mpi(1:n_soap_turbo) )
    allocate( n_sparse_mpi_distance_2b(1:n_distance_2b) )
    allocate( n_sparse_mpi_angle_3b(1:n_angle_3b) )
    allocate( n_mpi_core_pot(1:n_core_pot) )
    IF( rank == 0 )THEN
      n_species_mpi = soap_turbo_hypers(1:n_soap_turbo)%n_species
      n_sparse_mpi_soap_turbo = soap_turbo_hypers(1:n_soap_turbo)%n_sparse
      dim_mpi = soap_turbo_hypers(1:n_soap_turbo)%dim
      vdw_n_sparse_mpi_soap_turbo = soap_turbo_hypers(1:n_soap_turbo)%vdw_n_sparse
      has_vdw_mpi = soap_turbo_hypers(1:n_soap_turbo)%has_vdw
      compress_soap_mpi = soap_turbo_hypers(1:n_soap_turbo)%compress_soap
      n_sparse_mpi_distance_2b = distance_2b_hypers(1:n_distance_2b)%n_sparse
      n_sparse_mpi_angle_3b = angle_3b_hypers(1:n_angle_3b)%n_sparse
      n_mpi_core_pot = core_pot_hypers(1:n_core_pot)%n
    END IF
    !call cpu_time(time_mpi(1))
    time_mpi(1)=MPI_Wtime()
    call mpi_bcast(n_species_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_sparse_mpi_soap_turbo, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(dim_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vdw_n_sparse_mpi_soap_turbo, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(has_vdw_mpi, n_soap_turbo, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(compress_soap_mpi, n_soap_turbo, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_sparse_mpi_distance_2b, n_distance_2b, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_sparse_mpi_angle_3b, n_angle_3b, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_mpi_core_pot, n_core_pot, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call cpu_time(time_mpi(2))
    time_mpi(2)=MPI_Wtime()
    time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
    IF( rank /= 0 )THEN
      call allocate_soap_turbo_hypers(n_soap_turbo, n_species_mpi, n_sparse_mpi_soap_turbo, dim_mpi, &
                                      vdw_n_sparse_mpi_soap_turbo, has_vdw_mpi, compress_soap_mpi, &
                                      soap_turbo_hypers)
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
    !call cpu_time(time_mpi(1))
    time_mpi(1)=MPI_Wtime()
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
        call mpi_bcast(soap_turbo_hypers(i)%compress_soap_indices(1:dim), dim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
      call mpi_bcast(soap_turbo_hypers(i)%has_vdw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if( soap_turbo_hypers(i)%has_vdw )then
        n_sparse = soap_turbo_hypers(i)%vdw_n_sparse
        call mpi_bcast(soap_turbo_hypers(i)%vdw_delta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%vdw_zeta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%vdw_V0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%vdw_alphas(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(soap_turbo_hypers(i)%vdw_Qs(1:dim, 1:n_sparse), n_sparse*dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
    end do
    do i = 1, n_distance_2b
      n_sparse = distance_2b_hypers(i)%n_sparse
      call mpi_bcast(distance_2b_hypers(i)%alphas(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%cutoff(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%Qs(1:n_sparse, 1:1), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%delta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%sigma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%rcut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%buffer, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%species1, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(distance_2b_hypers(i)%species2, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    end do
    do i = 1, n_angle_3b
      n_sparse = angle_3b_hypers(i)%n_sparse
      call mpi_bcast(angle_3b_hypers(i)%alphas(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(angle_3b_hypers(i)%cutoff(1:n_sparse), n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(angle_3b_hypers(i)%Qs(1:n_sparse, 1:3), 3*n_sparse, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(angle_3b_hypers(i)%delta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(angle_3b_hypers(i)%sigma(1:3), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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
    !call cpu_time(time_mpi(2))
    time_mpi(2)=MPI_Wtime()
    time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
!   Clean up
    deallocate( n_species_mpi, n_sparse_mpi_soap_turbo, dim_mpi, compress_soap_mpi, n_sparse_mpi_distance_2b, &
                n_sparse_mpi_angle_3b, n_mpi_core_pot )
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
  !call cpu_time(time_read_input(2))
  time_read_input(2)=MPI_Wtime()
  time_read_input(3) = time_read_input(3) + time_read_input(2) - time_read_input(1)
!**************************************************************************






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
  md_istep = -1
  n_xyz = 0
  do while( repeat_xyz .or. ( params%do_md .and. md_istep < params%md_nsteps) )

    exit_loop = .false.

    if( params%do_md )then
      md_istep = md_istep + 1
    else
      n_xyz = n_xyz + 1
    end if

!   Update progress bar
    if( params%print_progress .and. counter == update_bar )then
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
    if( params%do_md .and. md_istep == 0 )then
      !call cpu_time(time_read_xyz(1))
      time_read_xyz(1)=MPI_Wtime()
#ifdef _MPIF90
      IF( rank == 0 )THEN
#endif
      call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                    n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                    positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                    xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                    n_sites, .false., fix_atom )
!     Only rank 0 handles these variables
!      allocate( positions_prev(1:3, 1:size(positions,2)) )
!      allocate( positions_diff(1:3, 1:size(positions,2)) )
      allocate( forces_prev(1:3, 1:n_sites) )
      allocate( positions_prev(1:3, 1:n_sites) )
      allocate( positions_diff(1:3, 1:n_sites) )
      positions_diff = 0.d0
      rebuild_neighbors_list = .true.
#ifdef _MPIF90
      END IF
#endif
      !call cpu_time(time_read_xyz(2))
      time_read_xyz(2)=MPI_Wtime()
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
      !call cpu_time(time_read_xyz(1))
      time_read_xyz(1)=MPI_Wtime()
#ifdef _MPIF90
      IF( rank == 0 )THEN
#endif
      call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                    n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                    positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                    xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                    n_sites, .false., fix_atom )
#ifdef _MPIF90
      END IF
#endif
      !call cpu_time(time_read_xyz(2))
      time_read_xyz(2)=MPI_Wtime()
      time_read_xyz(3) = time_read_xyz(3) + time_read_xyz(2) - time_read_xyz(1)
#ifdef _MPIF90
      !call cpu_time(time_mpi(1))
      time_mpi(1)=MPI_Wtime()
      call mpi_bcast(repeat_xyz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      !call cpu_time(time_mpi(2))
      time_mpi(2)=MPI_Wtime()
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
    END IF
    !call cpu_time(time_mpi(1))
    time_mpi(1)=MPI_Wtime()
    call mpi_bcast(n_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_sp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(n_sp_sc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call cpu_time(time_mpi(2))
    time_mpi(2)=MPI_Wtime()
    time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
    IF( rank /= 0 .and. (.not. params%do_md .or. md_istep == 0) )THEN
    allocate( positions(1:3, n_pos) )
    if( params%do_md )then
      allocate( velocities(1:3, n_pos) )
!      allocate( masses(n_pos) )
      allocate( masses(1:n_sp) )
      allocate( fix_atom(1:3, 1:n_sp) )
    end if
    allocate( xyz_species(1:n_sp) )
    allocate( species(1:n_sp) )
    allocate( xyz_species_supercell(1:n_sp_sc) )
    allocate( species_supercell(1:n_sp_sc) )
    END IF
    !call cpu_time(time_mpi_positions(1))
    time_mpi_positions(1)=MPI_Wtime()
    call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if( params%do_md )then
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
    !call cpu_time(time_mpi_positions(2))
    time_mpi_positions(2)=MPI_Wtime()
    time_mpi_positions(3) = time_mpi_positions(3) + time_mpi_positions(2) - time_mpi_positions(1)
#endif
!   Now that all ranks know the size of n_sites, we allocate do_list
    if( .not. params%do_md .or. (params%do_md .and. md_istep == 0) )then
      allocate( do_list(1:n_sites) )
      do_list = .true.
    end if
!
    !call cpu_time(time1)
    time1=MPI_Wtime()
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
    if( rebuild_neighbors_list )then
      call read_xyz(params%atoms_file, .true., params%all_atoms, params%do_timing, &
                    n_species, params%species_types, repeat_xyz, rcut_max, params%which_atom, &
                    positions, params%do_md, velocities, params%masses_types, masses, xyz_species, &
                    xyz_species_supercell, species, species_supercell, indices, a_box, b_box, c_box, &
                    n_sites, .true., fix_atom )
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

    ! ### --- OPENMP Parallelization --- ###    
    
    ! Now we split the neighbours that each omp_task will run
    
    ! > We know the total number of atom_pairs that each rank has,
    ! > so lets get the number of neighbours for the number of atoms
    ! > in the task
    
    ! We need to loop over the neighbours with the k index
    ! We want the total number of neighbours and the split the neighbour indices accordingly
    ! Set a pseudo number of omp tasks
    ! > n_omp: The number of omp tasks

    n_omp = 1
    omp_task = 0

!$omp parallel      
!$ n_omp = omp_get_num_threads()
!$omp end parallel     

    allocate(i_beg_omp(1:n_omp))
    allocate(i_end_omp(1:n_omp))
    allocate(j_beg_omp(1:n_omp))
    allocate(j_end_omp(1:n_omp))

    omp_n_sites = i_end - i_beg
    
!$omp parallel private(omp_task, k)
!$ omp_task = omp_get_thread_num()
    
    print *, "Allocated omp index arrays"
    if( omp_task < mod( omp_n_sites, n_omp ) )then
       i_beg_omp(omp_task+1) = i_beg + omp_task*(omp_n_sites / n_omp + 1)
    else
       i_beg_omp(omp_task+1) = i_beg + mod(omp_n_sites, n_omp)*(omp_n_sites / n_omp + 1) + (omp_task - mod( omp_n_sites, n_omp))*(omp_n_sites / n_omp)
    end if

    if( omp_task < mod( omp_n_sites, n_omp ) )then
       i_end_omp(omp_task+1) = i_beg + (omp_task+1)*(omp_n_sites / n_omp + 1)
    else
       i_end_omp(omp_task+1) = i_beg_omp(omp_task+1) + omp_n_sites/n_omp - 1
    end if

    print *,"omp_task = ", omp_task,  " i_beg_omp = ", i_beg_omp(omp_task+1)
    print *,"omp_task = ", omp_task,  " i_end_omp = ", i_end_omp(omp_task+1)

    
! --- Allocating the neighbors indices for reach openmp task 
    
    j_beg_omp = 1
    j_end_omp = n_atom_pairs_by_rank(rank+1)

    k = 0
    do i = i_beg, i_end
       if( i == i_beg_omp(omp_task+1) ) then 
          j_beg_omp(omp_task+1) = k+1
       end if

       k = k + n_neigh(i)

       if( i == i_end_omp(omp_task+1)) then 
          j_end_omp(omp_task+1) = k                       
       end if

    end do


    print *, "omp_task = ", omp_task, " j_beg_omp = ", j_beg_omp(omp_task+1)
    print *, "omp_task = ", omp_task, " j_end_omp = ", j_end_omp(omp_task+1)
    print *, "omp_task = ", omp_task, " k         = ", k
    print *, "omp_task = ", omp_task, " n_at_per_r= ", n_atom_pairs_by_rank(rank+1)
!$omp end parallel     
    
!   Compute the volume of the "primitive" unit cell
    v_uc = dot_product( cross_product(a_box, b_box), c_box ) / (dfloat(indices(1)*indices(2)*indices(3)))
    !call cpu_time(time2)
    time2=MPI_Wtime()
    time_neigh = time_neigh + time2 - time1
!**************************************************************************








!**************************************************************************
!   If we are doing prediction, we run this chunk of code
    if( params%do_prediction .or. params%write_soap .or. params%write_derivatives )then
      !call cpu_time(time1)
      time1=MPI_Wtime()

!     We only need to reallocate the arrays if the number of sites changes
      if( n_sites /= n_sites_prev )then
        if( allocated(energies) )deallocate( energies, energies_soap, energies_2b, energies_3b, energies_core_pot, &
                                             this_energies, energies_vdw, this_forces )
        allocate( energies(1:n_sites) )
        allocate( this_energies(1:n_sites) )
        allocate( energies_soap(1:n_sites) )
        allocate( energies_2b(1:n_sites) )
        allocate( energies_3b(1:n_sites) )
        allocate( energies_core_pot(1:n_sites) )
        allocate( energies_vdw(1:n_sites) )
!       This needs to be allocated even if no force prediction is needed:
        allocate( this_forces(1:3, 1:n_sites) )
      end if
      energies = 0.d0
      energies_soap = 0.d0
      energies_2b = 0.d0
      energies_3b = 0.d0
      energies_core_pot = 0.d0
      energies_vdw = 0.d0
      if( any( soap_turbo_hypers(:)%has_vdw ) )then
        if( n_sites /= n_sites_prev )then
          if( allocated(hirshfeld_v) )then
            nullify( this_hirshfeld_v_pt )
            deallocate( this_hirshfeld_v, hirshfeld_v )
            if( params%do_forces )then
              nullify( this_hirshfeld_v_cart_der_pt )
              deallocate( this_hirshfeld_v_cart_der, hirshfeld_v_cart_der )
            end if
          end if
          allocate( hirshfeld_v(1:n_sites) )
          allocate( this_hirshfeld_v(1:n_sites) )
!         I don't remember why this needs a pointer <----------------------------------------- CHECK
          this_hirshfeld_v_pt => this_hirshfeld_v
        end if
        hirshfeld_v = 0.d0
        if( params%do_forces )then
          if( n_atom_pairs_by_rank(rank+1) /= n_atom_pairs_by_rank_prev )then
            if( allocated(hirshfeld_v_cart_der) )deallocate( hirshfeld_v_cart_der, this_hirshfeld_v_cart_der )
            allocate( hirshfeld_v_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1)) )
            allocate( this_hirshfeld_v_cart_der(1:3, 1:n_atom_pairs_by_rank(rank+1)) )
          end if
          hirshfeld_v_cart_der = 0.d0
        end if
      end if
      if( params%do_forces )then
        if( n_sites /= n_sites_prev )then
          if( allocated(forces) )deallocate( forces, forces_soap, forces_2b, forces_3b, forces_core_pot, forces_vdw )
          allocate( forces(1:3, 1:n_sites) )
          allocate( forces_soap(1:3, 1:n_sites) )
          allocate( forces_2b(1:3, 1:n_sites) )
          allocate( forces_3b(1:3, 1:n_sites) )
          allocate( forces_core_pot(1:3, 1:n_sites) )
          allocate( forces_vdw(1:3,1:n_sites) )
        end if
        forces = 0.d0
        forces_soap = 0.d0
        forces_2b = 0.d0
        forces_3b = 0.d0
        forces_core_pot = 0.d0
        forces_vdw = 0.d0
        virial = 0.d0
        virial_soap = 0.d0
        virial_2b = 0.d0
        virial_3b = 0.d0
        virial_core_pot = 0.d0
        virial_vdw = 0.d0
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
      !call cpu_time(time_mpi_ef(1))
      time_mpi_ef(1)=MPI_Wtime()
      call mpi_reduce(energies, this_energies, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      !call cpu_time(time_mpi_ef(2))
      time_mpi_ef(2)=MPI_Wtime()
      time_mpi_ef(3) = time_mpi_ef(3) + time_mpi_ef(2) - time_mpi_ef(1)
      energies = this_energies
#endif

!     Loop through soap_turbo descriptors - we always call this routine, even if we don't want to do prediction
          !call cpu_time(soap_time_soap(1))
          soap_time_soap(1)=MPI_Wtime()
      do i = 1, n_soap_turbo
        !call cpu_time(time_soap(1))
        time_soap(1)=MPI_Wtime()
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
          if( soap_turbo_hypers(i)%has_vdw )then
            this_hirshfeld_v = 0.d0
            if( params%do_forces )then
              this_hirshfeld_v_cart_der = 0.d0
!             I don't remember why this needs a pointer <----------------------------------------- CHECK
              this_hirshfeld_v_cart_der_pt => this_hirshfeld_v_cart_der(1:3, this_j_beg:this_j_end)
            end if
          end if
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
                            soap_turbo_hypers(i)%compress_soap, soap_turbo_hypers(i)%compress_soap_indices, &
                            soap_turbo_hypers(i)%delta, soap_turbo_hypers(i)%zeta, soap_turbo_hypers(i)%central_species, &
                            xyz_species(this_i_beg:this_i_end), xyz_species_supercell, soap_turbo_hypers(i)%alphas, &
                            soap_turbo_hypers(i)%Qs, params%all_atoms, params%which_atom, indices, soap, soap_cart_der, &
                            der_neighbors, der_neighbors_list, &
                            soap_turbo_hypers(i)%has_vdw, soap_turbo_hypers(i)%vdw_Qs, soap_turbo_hypers(i)%vdw_alphas, &
                            soap_turbo_hypers(i)%vdw_zeta, soap_turbo_hypers(i)%vdw_delta, soap_turbo_hypers(i)%vdw_V0, &
                            this_energies, this_forces, this_hirshfeld_v_pt, this_hirshfeld_v_cart_der_pt, &
                            this_virial, solo_time_soap, time_get_soap, cublas_handle, gpu_stream)

          energies_soap = energies_soap + this_energies
          if( soap_turbo_hypers(i)%has_vdw )then
            hirshfeld_v = hirshfeld_v + this_hirshfeld_v
            if( params%do_forces )then
              hirshfeld_v_cart_der = hirshfeld_v_cart_der + this_hirshfeld_v_cart_der
            end if
          end if
          if( params%do_forces )then
            forces_soap = forces_soap + this_forces
            virial_soap = virial_soap + this_virial
            
            ! write(*,*) 
            ! write(*,*) this_virial(1,:)
            ! write(*,*) this_virial(2,:)
            ! write(*,*) this_virial(3,:) 
            ! write(*,*)
          end if
        end do
        !call cpu_time(soap_time_soap(2))
        soap_time_soap(2)=MPI_Wtime()
        deallocate( i_beg_list, i_end_list, j_beg_list, j_end_list )
        soap_time_soap(3)=soap_time_soap(3)+soap_time_soap(2)-soap_time_soap(1)


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

        !call cpu_time(time_soap(2))
        time_soap(2)=MPI_Wtime()
        time_soap(3) = time_soap(3) + time_soap(2) - time_soap(1)
      end do


#ifdef _MPIF90
      if( any( soap_turbo_hypers(:)%has_vdw ) )then
        !call cpu_time(time_mpi(1))
        time_mpi(1)=MPI_Wtime()
        call mpi_reduce(hirshfeld_v, this_hirshfeld_v, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!        if( params%do_forces )then
!         I'm not sure if this is necessary at all... CHECK
!          call mpi_reduce(hirshfeld_v_cart_der, this_hirshfeld_v_cart_der, 3*n_atom_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, &
!                          0, MPI_COMM_WORLD, ierr)
!          hirshfeld_v_cart_der = this_hirshfeld_v_cart_der
!        end if
        hirshfeld_v = this_hirshfeld_v
        call mpi_bcast(hirshfeld_v, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        !call cpu_time(time_mpi(2))
        time_mpi(2)=MPI_Wtime()
        time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
      end if
#endif



!     Compute vdW energies and forces
      if( any( soap_turbo_hypers(:)%has_vdw ) .and. params%do_prediction .and. params%vdw_type == "ts" )then
        !call cpu_time(time_vdw(1))
        time_vdw(1)=MPI_Wtime()
#ifdef _MPIF90
        allocate( this_energies_vdw(1:n_sites) )
        this_energies_vdw = 0.d0
        if( params%do_forces )then
          allocate( this_forces_vdw(1:3,1:n_sites) )
          this_forces_vdw = 0.d0
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
            v_neigh_vdw(k) = hirshfeld_v(j2)
          end do
        end do
        call get_ts_energy_and_forces( hirshfeld_v(i_beg:i_end), hirshfeld_v_cart_der(1:3, j_beg:j_end), &
                                       n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                                       neighbor_species(j_beg:j_end), &
                                       params%vdw_rcut, params%vdw_buffer, &
                                       params%vdw_rcut_inner, params%vdw_buffer_inner, &
                                       rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                                       params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
                                       params%vdw_alpha0_ref, params%do_forces, &
#ifdef _MPIF90
                                       this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw )
#else
                                       energies_vdw(i_beg:i_end), forces_vdw, virial_vdw )
#endif
        !call cpu_time(time_vdw(2))
        time_vdw(2)=MPI_Wtime()
        time_vdw(3) = time_vdw(2) - time_vdw(1)

        deallocate(v_neigh_vdw)
      end if





      if( params%do_prediction )then
!       Loop through distance_2b descriptors
        do i = 1, n_distance_2b
          !call cpu_time(time_2b(1))
          time_2b(1)=MPI_Wtime()
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
          !call cpu_time(time_2b(2))
          time_2b(2)=MPI_Wtime()
          time_2b(3) = time_2b(3) + time_2b(2) - time_2b(1)
        end do





!       Loop through core_pot descriptors
        do i = 1, n_core_pot
          !call cpu_time(time_core_pot(1))
          time_core_pot(1)=MPI_Wtime()
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
          !call cpu_time(time_core_pot(2))
          time_core_pot(2)=MPI_Wtime()
          time_core_pot(3) = time_core_pot(3) + time_core_pot(2) - time_core_pot(1)
        end do




!       Loop through angle_3b descriptors
        do i = 1, n_angle_3b
          !call cpu_time(time_3b(1))
          time_3b(1)=MPI_Wtime()
          this_energies = 0.d0
          if( params%do_forces )then
            this_forces = 0.d0
            this_virial = 0.d0
         end if

         ! Create custom indices to check that splitting into the specific indices works for openmp

         !
         !$omp parallel reduction(+: this_energies, this_forces, this_virial)
         !$ omp_task = omp_get_thread_num()
         print *, "doing omp parallel, with omp_task ", omp_task  
         call get_3b_energy_and_forces(rjs(j_beg_omp(omp_task&
              &+1):j_end_omp(omp_task+1)), xyz(1:3,j_beg_omp(omp_task&
              &+1):j_end_omp(omp_task+1)), angle_3b_hypers(i)%alphas,&
              & angle_3b_hypers(i)%cutoff, angle_3b_hypers(i)%rcut,&
              & 0.5d0, angle_3b_hypers(i)%delta, angle_3b_hypers(i)&
              &%sigma, 0.d0, angle_3b_hypers(i)%Qs,&
              & n_neigh(i_beg_omp(omp_task+1):i_end_omp(omp_task+1)),&
              & neighbors_list(j_beg_omp(omp_task&
              &+1):j_end_omp(omp_task+1)), params%do_forces, params&
              &%do_timing, angle_3b_hypers(i)%kernel_type,&
              & species(i_beg_omp(omp_task+1):i_end_omp(omp_task+1)),&
              & neighbor_species(j_beg_omp(omp_task&
              &+1):j_end_omp(omp_task+1)), angle_3b_hypers(i)&
              &%species_center, angle_3b_hypers(i)%species1,&
              & angle_3b_hypers(i)%species2, params%species_types,&
              & this_energies(i_beg_omp(omp_task&
              &+1):i_end_omp(omp_task+1)), this_forces, this_virial)
       !$omp end parallel
         
            ! call get_3b_energy_and_forces(rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), angle_3b_hypers(i)%alphas, &
            !                             angle_3b_hypers(i)%cutoff, &
            !                             angle_3b_hypers(i)%rcut, 0.5d0, angle_3b_hypers(i)%delta, &
            !                             angle_3b_hypers(i)%sigma, 0.d0, angle_3b_hypers(i)%Qs, n_neigh(i_beg:i_end), &
            !                             neighbors_list(j_beg:j_end), &
            !                             params%do_forces, params%do_timing, angle_3b_hypers(i)%kernel_type, &
            !                             species(i_beg:i_end), neighbor_species(j_beg:j_end), angle_3b_hypers(i)%species_center, &
            !                             angle_3b_hypers(i)%species1, angle_3b_hypers(i)%species2, params%species_types, &
            !                             this_energies(i_beg:i_end), this_forces, this_virial)
          
          energies_3b = energies_3b + this_energies
          if( params%do_forces )then
            forces_3b = forces_3b + this_forces
            virial_3b = virial_3b + this_virial
          end if
          !call cpu_time(time_3b(2))
          time_3b(2)=MPI_Wtime()
          time_3b(3) = time_3b(3) + time_3b(2) - time_3b(1)
        end do


        !call cpu_time(time2)
        time2=MPI_Wtime()
        time_gap = time_gap + time2 - time1



!       Communicate all energies and forces here for all terms
#ifdef _MPIF90
        !call cpu_time(time_mpi_ef(1))
        time_mpi_ef(1)=MPI_Wtime()
        counter2 = 0
        if( n_soap_turbo > 0 )then
          counter2 = counter2 + 1
        end if
        if( allocated(this_energies_vdw) )then
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
        call mpi_reduce(all_energies, all_this_energies, n_sites*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if( params%do_forces )then
          call mpi_reduce(all_forces, all_this_forces, 3*n_sites*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call mpi_reduce(all_virial, all_this_virial, 9*counter2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
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
        if( n_distance_2b > 0 )then
          counter2 = counter2 + 1
          energies_2b(1:n_sites) = all_this_energies(1:n_sites, counter2)
          if( params%do_forces )then
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

        !call cpu_time(time_mpi_ef(2))
        time_mpi_ef(2)=MPI_Wtime()
        time_mpi_ef(3) = time_mpi_ef(3) + time_mpi_ef(2) - time_mpi_ef(1)
#endif




!       Add up all the energy terms
        energies = energies + energies_soap + energies_2b + energies_3b + energies_core_pot + energies_vdw
        energy_prev = energy
        energy = sum(energies)
      end if

      if( .not. params%do_md )then
#ifdef _MPIF90
        IF( rank == 0 )then
#endif
        write(*,*)'                                       |'
        write(*,'(A,1X,F22.8,1X,A)')' SOAP energy:', sum(energies_soap), 'eV |'
        write(*,'(A,1X,F24.8,1X,A)')' 2b energy:', sum(energies_2b), 'eV |'
        write(*,'(A,1X,F24.8,1X,A)')' 3b energy:', sum(energies_3b), 'eV |'
        write(*,'(A,1X,F18.8,1X,A)')' core_pot energy:', sum(energies_core_pot), 'eV |'
        write(*,'(A,1X,F23.8,1X,A)')' vdw energy:', sum(energies_vdw), 'eV |'
        write(*,'(A,1X,F21.8,1X,A)')' Total energy:', sum(energies), 'eV |'
        write(*,*)'                                       |'
        write(*,*)'Energy & forces in "trajectory_out.xyz"|'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
#ifdef _MPIF90
        END IF
#endif
      end if


      if( params%do_forces )then
        forces = forces_soap + forces_2b + forces_3b + forces_core_pot + forces_vdw
        virial = virial_soap + virial_2b + virial_3b + virial_core_pot + virial_vdw
      end if


! For debugging the virial implementation
if( rank == 0 .and. .false. )then
write(*,*) "pressure_soap: ", virial_soap / 3.d0 / v_uc
write(*,*) "pressure_vdw: ", virial_vdw / 3.d0 / v_uc
write(*,*) "pressure_2b: ", virial_2b / 3.d0 / v_uc
write(*,*) "pressure_3b: ", virial_3b / 3.d0 / v_uc
write(*,*) "pressure_core_pot: ", virial_core_pot / 3.d0 / v_uc
end if



      if( params%do_prediction .and. .not. params%do_md )then
#ifdef _MPIF90
        IF( rank == 0 )then
#endif
!       Write energy and forces if we're just doing static predictions
!       The masses should be divided by 103.6426965268d0 to have amu units, but
!       since masses is not allocated for single point calculations, it would
!       likely lead to a segfault
        call wrap_pbc(positions(1:3,1:n_sites), a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
        call write_extxyz( n_sites, -n_xyz, time_step, instant_temp, instant_pressure, &
                           a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                           virial, xyz_species, &
                           positions(1:3, 1:n_sites), velocities, &
                           forces, energies(1:n_sites), masses, hirshfeld_v, &
                           params%write_property, params%write_array_property, fix_atom )
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
    if( params%do_md )then
      !call cpu_time(time_md(1))
      time_md(1)=MPI_Wtime()
!     Define the time_step and md_time prior to possible scaling (see variable_time_step below)
      if( md_istep > 0 )then
        md_time = md_time + time_step
      else
        md_time = 0.d0
        time_step = params%md_step
      end if
!     We wrap the positions and remoce CM velocity
      call wrap_pbc(positions(1:3,1:n_sites), a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
      call remove_cm_vel(velocities(1:3,1:n_sites), masses(1:n_sites))

!     First we check if this is a variable time step simulation
      if( params%variable_time_step )then
        call variable_time_step(md_istep == 0, velocities(1:3, 1:n_sites), forces(1:3, 1:n_sites), masses(1:n_sites), &
                                params%target_pos_step, params%tau_dt, params%md_step, time_step)
      end if
!     This takes care of NVE
!     Velocity Verlet takes positions for t, positions_prev for t-dt, and velocities for t-dt and returns everything
!     dt later. forces are taken at t, and forces_prev at t-dt. forces is left unchanged by the routine, and
!     forces_prev is returned as equal to forces (both arrays contain the same information on return)
      if( params%optimize == "vv" )then
        call velocity_verlet(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                             forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), time_step, &
                             md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                             fix_atom(1:3, 1:n_sites))
      else if( params%optimize == "gd" )then
        call gradient_descent(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), velocities(1:3, 1:n_sites), &
                              forces(1:3, 1:n_sites), forces_prev(1:3, 1:n_sites), masses(1:n_sites), params%gamma0, &
                              params%max_opt_step, md_istep == 0, a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), &
                              c_box/dfloat(indices(3)), fix_atom(1:3, 1:n_sites), energy)
      end if
!     Compute kinetic energy from current velocities. Because Velocity Verlet
!     works with the velocities at t-dt (except for the first time step) we
!     have to compute the velocities after call Verlet
      E_kinetic = 0.d0
      do i = 1, n_sites
        E_kinetic = E_kinetic + 0.5d0 * masses(i) * dot_product(velocities(1:3, i), velocities(1:3, i))
      end do
      instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic

!     Instant pressure in bar
      instant_pressure = (kB*dfloat(n_sites-1)*instant_temp+(virial(1,1) + virial(2,2) + virial(3,3))/3.d0)/v_uc*eVperA3tobar
      instant_pressure_tensor(1:3, 1:3) = virial(1:3,1:3)/v_uc*eVperA3tobar
      do i = 1, 3
        instant_pressure_tensor(i, i) = instant_pressure_tensor(i, i) + (kB*dfloat(n_sites-1)*instant_temp)/v_uc*eVperA3tobar
      end do

!     Here we write thermodynamic information -> THIS NEEDS CLEAN UP AND IMPROVEMENT
      if( md_istep == 0 )then
        open(unit=10, file="thermo.log", status="unknown")
      else
        open(unit=10, file="thermo.log", status="old", position="append")
      end if
      if( md_istep == 0 .or. md_istep == params%md_nsteps .or. modulo(md_istep, params%write_thermo) == 0 )then
!       Organize this better so that the user can have more freedom about what gets printed to thermo.log
!       There should also be a header preceded by # specifying what gets printed
        write(10, "(I10, 1X, F16.4, 1X, F16.4, 1X, F20.8, 1X, F20.8, 1X, F20.8)", advance="no") &
             md_istep, md_time, instant_temp, E_kinetic, sum(energies), instant_pressure
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
      if( params%do_md .and. params%optimize == "gd" .and. md_istep > 0 .and. &
          abs(energy-energy_prev) < params%e_tol .and. rank == 0 )then
        exit_loop = .true.
      end if

!     We write out the trajectory file. We write positions_prev which is then one for which we have computed
!     the properties. positions_prev and velocities are synchronous
      if( md_istep == 0 .or. md_istep == params%md_nsteps .or. modulo(md_istep, params%write_xyz) == 0 .or. &
          exit_loop )then
        call wrap_pbc(positions_prev(1:3,1:n_sites), a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
        call write_extxyz( n_sites, md_istep, time_step, instant_temp, instant_pressure, &
                           a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), c_box/dfloat(indices(3)), &
                           virial, xyz_species, &
                           positions_prev(1:3, 1:n_sites), velocities, &
                           forces, energies(1:n_sites), masses/103.6426965268d0, hirshfeld_v, &
                           params%write_property, params%write_array_property, fix_atom(1:3, 1:n_sites) )
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
      end if
!     If there are thermostating operations they happen here
      if( params%thermostat == "berendsen" )then
        call berendsen_thermostat(velocities(1:3, 1:n_sites), &
                                  params%t_beg + (params%t_end-params%t_beg)*dfloat(md_istep+1)/float(params%md_nsteps), &
                                  instant_temp, params%tau_t, time_step)
      else if( params%thermostat == "bussi" )then
        velocities(1:3, 1:n_sites) = velocities(1:3, 1:n_sites) * dsqrt(resamplekin(E_kinetic, &
                                     params%t_beg + (params%t_end-params%t_beg)*dfloat(md_istep+1)/float(params%md_nsteps), &
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
      !call cpu_time(time_md(2))
      time_md(2)=MPI_Wtime()
      time_md(3) = time_md(3) + time_md(2) - time_md(1)
    end if
#ifdef _MPIF90
    END IF
    call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
!   Make sure all ranks have correct positions and velocities
#ifdef _MPIF90
    if( params%do_md )then
      !call cpu_time(time_mpi_positions(1))
      time_mpi_positions(1)=MPI_Wtime()
      n_pos = size(positions,2)
      call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(velocities, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      !call cpu_time(time_mpi_positions(2))
      time_mpi_positions(2)=MPI_Wtime()
      time_mpi_positions(3) = time_mpi_positions(3) + time_mpi_positions(2) - time_mpi_positions(1)
    end if
#endif


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


    if( rebuild_neighbors_list )then
      deallocate( rjs, xyz, thetas, phis, neighbor_species )
      deallocate( neighbors_list, n_neigh )
#ifdef _MPIF90
      deallocate( n_neigh_local )
#endif
    end if
    if( .not. params%do_md .or. (params%do_md .and. md_istep == params%md_nsteps) )then
      deallocate( positions, xyz_species, xyz_species_supercell, species, species_supercell, do_list )
    end if
    if( params%do_md .and. md_istep == params%md_nsteps .and. rank == 0 )then
      deallocate( positions_prev, forces_prev )
    end if

    n_sites_prev = n_sites
    n_atom_pairs_by_rank_prev = n_atom_pairs_by_rank(rank+1)

#ifdef _MPIF90
    call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
    if( exit_loop )exit
! End of loop through structures in the xyz file or MD steps
  end do


  if( params%do_md .or. params%do_prediction )then
    !call cpu_time(time2)
    time2=MPI_Wtime()
#ifdef _MPIF90
    IF( rank == 0 )then
#endif
    if( params%do_md )then
!      write(*,'(A)')'] |'
      write(*,*)
      write(*,*)'                                       |'
!      write(*,'(I8,A,F13.3,A)') params%md_nsteps, ' MD steps:', time2-time3, ' seconds |'
      write(*,'(I8,A,F13.3,A)') md_istep, ' MD steps:', time2-time3, ' seconds |'
    end if
    write(*,*)'                                       |'
    write(*,'(A,F13.3,A)') ' *     Read input:', time_read_input(3), ' seconds |'
    write(*,'(A,F13.3,A)') ' * Read XYZ files:', time_read_xyz(3), ' seconds |'
    write(*,'(A,F13.3,A)') ' * Neighbor lists:', time_neigh, ' seconds |'
    write(*,'(A,F13.3,A)') ' *  GAP desc/pred:', time_gap, ' seconds |'
    write(*,'(A,F13.3,A)') '     - soap_turbo:', time_soap(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     - lolo__soap:', soap_time_soap(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     - get___soap:', time_get_soap, ' seconds |'
    write(*,'(A,F13.3,A)') '     - lin__turbo:', solo_time_soap, ' seconds |'
    write(*,'(A,F13.3,A)') '     -         2b:', time_2b(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     -         3b:', time_3b(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     -   core_pot:', time_core_pot(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     -        vdw:', time_vdw(3), ' seconds |'
    if( params%do_md )then
      write(*,'(A,F13.3,A)') ' *  MD algorithms:', time_md(3), ' seconds |'
    end if
#ifdef _MPIF90
    write(*,'(A,F13.3,A)') ' *  MPI comms.   :', time_mpi(3) + time_mpi_positions(3) + time_mpi_ef(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     -  pos & vel:', time_mpi_positions(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     - E & F brc.:', time_mpi_ef(3), ' seconds |'
    write(*,'(A,F13.3,A)') '     -  MPI misc.:', time_mpi(3), ' seconds |'
    write(*,'(A,F13.3,A)') ' *  Miscellaneous:', time2-time3 - time_neigh - time_gap - time_read_input(3) &
      - time_read_xyz(3) - time_mpi(3) - time_mpi_positions(3) - time_mpi_ef(3) - time_md(3), ' seconds |'
#else
    write(*,'(A,F13.3,A)') ' *  Miscellaneous:', time2-time3 - time_neigh - time_gap - time_read_input(3) &
      - time_read_xyz(3) - time_md(3), ' seconds |'
#endif
    write(*,*)'                                       |'
    write(*,'(A,F13.3,A)') ' *     Total time:', time2-time3, ' seconds |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
#ifdef _MPIF90
    END IF
#endif
  end if




  if( allocated(soap_turbo_hypers) )deallocate(soap_turbo_hypers)
  if( allocated(distance_2b_hypers) )deallocate(distance_2b_hypers)
  if( allocated(angle_3b_hypers) )deallocate(angle_3b_hypers)
  if( allocated(core_pot_hypers) )deallocate(core_pot_hypers)

  !call destroy_cublas_handle(cublas_handle)
  call gpu_device_reset()

! write(*,*) "Starting dummy kernel"

! dut1= MPI_Wtime()
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
! dut2= MPI_Wtime()
! write(*,*) "Ending dummy region"
! write(*,*) "Time spent in dummy region", dut2-dut1
   

#ifdef _MPIF90
  IF( rank == 0 )then
#endif
  write(*,*)'                                       |'
  write(*,*)'End of execution                       |'
  write(*,*)'_______________________________________/'
#ifdef _MPIF90
  END IF
#endif


write(*,*) "    - lin__turbo:", solo_time_soap, rank


#ifdef _MPIF90
  call mpi_finalize(ierr)
#endif

!call destroy_cublas_handle(cublas_handle, gpu_stream)

end program turbogap
