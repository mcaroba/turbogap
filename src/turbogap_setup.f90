! Copyright (c) 2020-2026 by Albert Bartok and Miguel Caro
!
! Start-up phase of the TurboGAP driver: read the input file, read and
! broadcast the GAP hyperparameters, and initialise the electronic-stopping
! models.
!
! Lifted verbatim from turbogap.f90, where it ran inline between the "Read
! input file and other files" banner and the printout section. The computation
! is unchanged. What the move buys is 26 fewer variables in the driver's shared
! scope: the per-rank MPI staging arrays for broadcasting hypers
! (n_sparse_mpi_*, dim_mpi, compress_*), the parsing scratch (keyword,
! file_gap, the local_property_labels_temp pair) and the loop counters were
! only ever touched here.
!
! Seven of those -- i, j, ierr, iostatus, cjunk, n_sp, n_lp_count -- are shared
! scratch elsewhere in the driver. They are safe to localise because each is
! written before it is next read after this block; that was checked
! individually, since no test would catch a leaked value.

module turbogap_setup

   use kinds

   use timing
   use types
   use read_files
   use local_prop
   use adaptive_time
   use electronic_stopping
   use eph_fdm
   use eph_beta
   use eph_electronic_stopping
#ifdef _MPIF90
   use mpi
   use mpi_helper
#endif

   implicit none

   private
   public :: read_input_and_gap_files

contains

   subroutine read_input_and_gap_files(mode, rank, ntasks, params, &
                                       soap_turbo_hypers, distance_2b_hypers, angle_3b_hypers, core_pot_hypers, &
                                       n_soap_turbo, n_distance_2b, n_angle_3b, n_core_pot, n_species, rcut_max, &
                                       valid_xps, xps_idx, vdw_lp_index, core_be_lp_index, &
                                       local_property_labels, local_property_indexes, n_local_properties_mpi, &
                                       has_local_properties_mpi, local_properties_n_sparse_mpi_soap_turbo, &
                                       local_properties_dim_mpi_soap_turbo, nrows, allelstopdata, &
                                       ephbeta, ephfdm, ephlsc, time_read_input, time_mpi)

      implicit none

!   Input
      character*16, intent(in) :: mode
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks

!   Everything below is filled in by this routine
      type(input_parameters), intent(inout) :: params
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)
      type(distance_2b), allocatable, intent(inout) :: distance_2b_hypers(:)
      type(angle_3b), allocatable, intent(inout) :: angle_3b_hypers(:)
      type(core_pot), allocatable, intent(inout) :: core_pot_hypers(:)
      integer, intent(inout) :: n_soap_turbo
      integer, intent(inout) :: n_distance_2b
      integer, intent(inout) :: n_angle_3b
      integer, intent(inout) :: n_core_pot
      integer, intent(inout) :: n_species
      real(dp), intent(inout) :: rcut_max
      logical, intent(inout) :: valid_xps
      integer, intent(inout) :: xps_idx
      integer, intent(inout) :: vdw_lp_index
      integer, intent(inout) :: core_be_lp_index
      character*1024, allocatable, intent(inout) :: local_property_labels(:)
      integer, allocatable, intent(inout) :: local_property_indexes(:)
      integer, allocatable, intent(inout) :: n_local_properties_mpi(:)
      logical, allocatable, intent(inout) :: has_local_properties_mpi(:)
      integer, allocatable, intent(inout) :: local_properties_n_sparse_mpi_soap_turbo(:)
      integer, allocatable, intent(inout) :: local_properties_dim_mpi_soap_turbo(:)
      integer, intent(inout) :: nrows
      real(dp), allocatable, intent(inout) :: allelstopdata(:)
      type(EPH_Beta_class), intent(inout) :: ephbeta
      type(EPH_FDM_class), intent(inout) :: ephfdm
      type(EPH_LangevinSpatialCorrelation_class), intent(inout) :: ephlsc
      real(dp), intent(inout) :: time_read_input(1:3)
      real(dp), intent(inout) :: time_mpi(1:3)

!   Local. All of these were variables of the main program that nothing
!   outside this block referenced.
      character*1024 :: file_gap = "none"
      character*1024 :: cjunk
      character*1024, allocatable :: local_property_labels_temp(:)
      character*1024, allocatable :: local_property_labels_temp2(:)
      character*64 :: keyword
      character*1 :: keyword_first
      integer, allocatable :: n_species_mpi(:)
      integer, allocatable :: n_sparse_mpi_soap_turbo(:)
      integer, allocatable :: dim_mpi(:)
      integer, allocatable :: n_sparse_mpi_distance_2b(:)
      integer, allocatable :: n_sparse_mpi_angle_3b(:)
      integer, allocatable :: n_mpi_core_pot(:)
      integer, allocatable :: compress_P_nonzero_mpi(:)
      logical, allocatable :: compress_soap_mpi(:)
      logical :: valid_vdw = .false.
      integer :: n_sparse
      integer :: dim
      integer :: cPnz
      integer :: n_nonzero
      integer :: n_local_properties_tot = 0
      integer :: i
      integer :: j
      integer :: ierr
      integer :: iostatus
      integer :: n_sp
      integer :: n_lp_count
      time_read_input(3) = 0.d0
      call get_time(time_read_input(1))
      open (unit=10, file='input', status='old', iostat=iostatus)
      ! Check for existence of input file
#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         write (*, *) '                                       |'
         write (*, *) 'Checking input file...                 |'
#ifdef _MPIF90
      END IF
#endif
      if (iostatus /= 0) then
         close (10)
#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            write (*, *) '                                       |'
            write (*, *) 'ERROR: input file could not be found   |  <-- ERROR'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
            write (*, *) '                                       |'
            write (*, *) 'End of execution                       |'
            write (*, *) '_______________________________________/'
#ifdef _MPIF90
         END IF
         call mpi_finalize(ierr)
#endif
         stop
      end if
      !
      ! First, we look for n_species, which determines how we allocate the species-specific arrays
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'n_species') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, cjunk, n_species
            if (n_species < 1) then
#ifdef _MPIF90
               IF (rank == 0) THEN
#endif
                  write (*, *) '                                       |'
                  write (*, *) 'ERROR: n_species must be > 0           |  <-- ERROR'
                  write (*, *) '                                       |'
                  write (*, *) '.......................................|'
#ifdef _MPIF90
               END IF
               call mpi_finalize(ierr)
#endif
               stop
            end if
         end if
      end do
      ! Let's look for those and other options in the input file
      rewind (10)
      call read_input_file(n_species, mode, params, rank)

! Make randomized initial velocities (and any other use of random_number)
! reproducible when the input asks for it, so runs can be compared.
      call init_random_seed(params%random_seed_value, rank)

!! If electronic stopping is required to be done, then read the stopping data file for once
!! Reading and storing the elctronic stopping data
      if (params%electronic_stopping) then
         if (params%estop_filename == 'NULL') then
            write (*, *) "ERROR: No stopping data file is provided."
            stop
         else
            call read_electronic_stopping_file(n_species, params%species_types, params%estop_filename, nrows, allelstopdata)
         end if
      end if
!! -----------------------------                                ---- untill here for reading stopping data

      ! TEMPORARY ERROR, FIX THE UNDERLYING ISSUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _MPIF90
      IF (rank == 0 .and. ntasks > 1 .and. (params%write_soap .or. params%write_derivatives)) THEN
         write (*, *) '                                       |'
         write (*, *) 'ERROR: writing of SOAP and/or SOAP     |  <-- ERROR'
         write (*, *) 'derivatives cannot currently be done in|'
         write (*, *) 'parallel. Try using the serial code or |'
         write (*, *) 'running the MPI code on a single CPU.  |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
         write (*, *) '                                       |'
         write (*, *) 'End of execution                       |'
         write (*, *) '_______________________________________/'
      END IF
#endif
      !
      ! Second, we look for pot_file, which contains the GAP difinitions
      rewind (10)
      iostatus = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'pot_file') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, cjunk, file_gap
            exit
         end if
      end do
      close (10)
      ! Now, read file_gap and register each GAP, including the hypers for its descriptor
      if (file_gap /= "none") then
#ifdef _MPIF90
         IF (rank == 0) then
#endif
            call read_gap_hypers(file_gap, &
                                 n_soap_turbo, soap_turbo_hypers, &
                                 n_distance_2b, distance_2b_hypers, &
                                 n_angle_3b, angle_3b_hypers, &
                                 n_core_pot, core_pot_hypers, &
                                 rcut_max, params%do_prediction, &
                                 params)
            !   Check if vdw_rcut is bigger
            if (params%vdw_rcut > rcut_max) then
               rcut_max = params%vdw_rcut
            end if
            if (params%xrd_rcut > rcut_max) then
               rcut_max = params%xrd_rcut
            end if
            if (params%nd_rcut > rcut_max) then
               rcut_max = params%nd_rcut
            end if
            if (params%pair_distribution_rcut > rcut_max) then
               rcut_max = params%pair_distribution_rcut
            end if

            !   We increase rcut_max by the neighbors buffer
            rcut_max = rcut_max + params%neighbors_buffer

            ! Check that the local properties we want to compute are valid,
            ! so things are commensurate between input file and .gap model

            ! Need to set the number of local properties, and get an array of labels and sizes for broadcasting

            write (*, *) '                                       |'
            write (*, *) '.......................................|'
            write (*, *) '                                       |'

            call get_irreducible_local_properties(params, n_local_properties_tot, n_soap_turbo, soap_turbo_hypers, &
                           local_property_labels, local_property_labels_temp, local_property_labels_temp2, local_property_indexes, &
                                                  valid_vdw, vdw_lp_index, core_be_lp_index, valid_xps, xps_idx)

            if (params%n_local_properties > 0) then
               write (*, *) '                                       |'
               write (*, *) ' Irreducible local properties:         |'
               do i = 1, params%n_local_properties
                  write (*, '(A41)') trim(local_property_labels(i))//' |'
               end do

               allocate (params%write_local_properties(1:params%n_local_properties))
               params%write_local_properties = .true.
            end if

#ifdef _MPIF90
         END IF
#endif
         !   THIS CHUNK HERE DISTRIBUTES THE INPUT DATA AMONG ALL THE PROCESSES
         !   Broadcast number of descriptors to other processes
#ifdef _MPIF90
         call get_time(time_mpi(1))
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
         call get_time(time_mpi(2))
         time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
         allocate (n_species_mpi(1:n_soap_turbo))
         allocate (n_sparse_mpi_soap_turbo(1:n_soap_turbo))
         allocate (dim_mpi(1:n_soap_turbo))
         allocate (n_local_properties_mpi(1:n_soap_turbo))
         if (n_local_properties_tot > 0) then
            allocate (local_properties_n_sparse_mpi_soap_turbo(1:n_local_properties_tot))
            allocate (local_properties_dim_mpi_soap_turbo(1:n_local_properties_tot))
            if (rank /= 0) allocate (local_property_indexes(1:n_local_properties_tot))
            if (rank /= 0) allocate (params%write_local_properties(1:params%n_local_properties))
         end if
         allocate (has_local_properties_mpi(1:n_soap_turbo))
         allocate (compress_soap_mpi(1:n_soap_turbo))
         allocate (n_sparse_mpi_distance_2b(1:n_distance_2b))
         allocate (n_sparse_mpi_angle_3b(1:n_angle_3b))
         allocate (n_mpi_core_pot(1:n_core_pot))
         allocate (compress_P_nonzero_mpi(1:n_soap_turbo))
         IF (rank == 0) THEN
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
            if (any(soap_turbo_hypers(:)%has_local_properties)) then
               n_lp_count = 0
               do i = 1, n_soap_turbo
                  if (n_local_properties_mpi(i) > 0) then
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
         call get_time(time_mpi(1))
         call mpi_bcast(n_species_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(n_sparse_mpi_soap_turbo, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(dim_mpi, n_soap_turbo, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (n_local_properties_tot > 0) then
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
         call get_time(time_mpi(2))
         time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)

         IF (rank /= 0) THEN
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
         call get_time(time_mpi(1))
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
            if (soap_turbo_hypers(i)%compress_soap) then
               cPnz = soap_turbo_hypers(i)%compress_P_nonzero
               call mpi_bcast(soap_turbo_hypers(i)%compress_P_el(1:cPnz), cPnz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
               call mpi_bcast(soap_turbo_hypers(i)%compress_P_i(1:cPnz), cPnz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
               call mpi_bcast(soap_turbo_hypers(i)%compress_P_j(1:cPnz), cPnz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            end if
            call mpi_bcast(soap_turbo_hypers(i)%has_local_properties, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(soap_turbo_hypers(i)%has_core_electron_be, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            if (valid_xps .or. params%do_xps) call mpi_bcast(core_be_lp_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            if (valid_vdw) call mpi_bcast(vdw_lp_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call mpi_bcast(soap_turbo_hypers(i)%has_vdw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(soap_turbo_hypers(i)%n_local_properties, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            if (soap_turbo_hypers(i)%has_local_properties) then
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
                       &%local_property_models(j)%alphas(1:n_sparse),&
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
         call get_time(time_mpi(2))
         time_mpi(3) = time_mpi(3) + time_mpi(2) - time_mpi(1)
         !   Clean up
         deallocate (n_species_mpi, n_sparse_mpi_soap_turbo, dim_mpi, compress_soap_mpi, n_sparse_mpi_distance_2b, &
                    n_sparse_mpi_angle_3b, n_mpi_core_pot, compress_P_nonzero_mpi, n_local_properties_mpi, has_local_properties_mpi)
         if (allocated(local_properties_dim_mpi_soap_turbo)) deallocate (&
              & local_properties_dim_mpi_soap_turbo,&
              & local_properties_n_sparse_mpi_soap_turbo)
#endif
      else
#ifdef _MPIF90
         IF (rank == 0) THEN
#endif
            write (*, *) '                                       |'
            write (*, *) 'ERROR: you must provide a "pot_file"   |  <-- ERROR'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
#ifdef _MPIF90
         END IF
         call mpi_finalize(ierr)
#endif
         stop
      end if
      call get_time(time_read_input(2))
      time_read_input(3) = time_read_input(3) + time_read_input(2) - time_read_input(1)
      !**************************************************************************

!! If electronic stopping based on eph model is to be calculated, these data structures are required to be
!! initialized first.
      if (params%nonadiabatic_processes) then
         if (params%eph_Tinfile /= "NULL") then
            call ephfdm%EPH_FDM_input_file(params%eph_Tinfile, params%eph_md_last_step)
         end if
         if (params%eph_Tinfile == "NULL") then
            call ephfdm%EPH_FDM_input_params(params%eph_md_last_step, params%eph_gsx, &
                                           params%eph_gsy, params%eph_gsz, params%in_x0, params%in_x1, params%in_y0, params%in_y1, &
                                             params%in_z0, params%in_z1, params%eph_Ti_e, params%eph_C_e, params%eph_rho_e, &
                                             params%eph_kappa_e, params%eph_fdm_steps)
         end if
         call ephbeta%beta_parameters(params%eph_betafile, n_species)
         call ephlsc%eph_InitialValues(params%eph_friction_option, params%eph_random_option, params%eph_fdm_option, &
                                       params%eph_Toutfile, params%eph_freq_Tout, params%eph_freq_mesh_Tout, params%model_eph, &
                                       params%eph_E_prev_time, params%eph_md_prev_time)
      end if
!! -------------------------                                ---- untill here for initializing eph model elec. stopping

   end subroutine read_input_and_gap_files

end module turbogap_setup
