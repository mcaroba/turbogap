! HND XXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   Copyright (c) 2026, Miguel A. Caro and Tigany Zarrouk
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXX

! How many OpenMP threads a rank takes.
!
! The host loops that carry !$omp directives -- the neighbour build and the
! per-step pair geometry in src/neighbors.f90 -- were unreachable until the
! build started passing -fopenmp: with one rank per GPU a single core built the
! whole list while the device idled. On a 125k-atom cell the neighbour bucket
! is 3.14 s on one thread and 0.47 s on twelve.
!
! The hazard runs the other way. OpenMP's own default is one thread per core,
! so every rank on a node would take every core: four ranks on twelve cores
! would run thirty-six threads over twelve, and the build would be slower than
! serial. So unless the user has said what they want, the cores are divided by
! the ranks that share the node.
!
! OMP_NUM_THREADS, when set, is left alone. It is the user saying what they
! want, and a batch system that sets it is saying it on their behalf.

module threads

#ifdef _MPIF90
   use mpi
#endif
!$ use omp_lib

   implicit none

   private
   public :: threads_init
   public :: threads_in_use

contains

!  Called once, after MPI is up.
   subroutine threads_init()

      implicit none

      character(len=32) :: requested
      integer :: status
      integer :: length

!$    call get_environment_variable("OMP_NUM_THREADS", requested, length, status)
!$    if (status == 0 .and. length > 0) return
!$    call omp_set_num_threads(max(1, omp_get_num_procs()/ranks_sharing_node()))

   end subroutine threads_init

!  One without OpenMP, so callers can print it unconditionally.
   integer function threads_in_use()

      implicit none

      threads_in_use = 1
!$    threads_in_use = omp_get_max_threads()

   end function threads_in_use

!  The ranks that can see this node's cores. MPI_COMM_TYPE_SHARED answers it
!  exactly; the same question is asked of memory in gpu_context_cpu.f90, and
!  for the same reason -- dividing by the total rank count is right only when
!  the job runs on one node.
   integer function ranks_sharing_node()

      implicit none

#ifdef _MPIF90
      integer :: shared_comm
      integer :: ierr
      integer :: n
#endif

      ranks_sharing_node = 1

#ifdef _MPIF90
      call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                               MPI_INFO_NULL, shared_comm, ierr)
      if (ierr == MPI_SUCCESS) then
         call mpi_comm_size(shared_comm, n, ierr)
         if (ierr == MPI_SUCCESS .and. n > 0) ranks_sharing_node = n
         call mpi_comm_free(shared_comm, ierr)
      end if
#endif

   end function ranks_sharing_node

end module threads
