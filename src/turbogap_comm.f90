! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_comm.f90, is copyright (c) 2026, Tigany Zarrouk
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

!  Transport: the only place the driver's collectives touch MPI. Root is always
!  rank 0. Array routines take the element count explicitly, as MPI does, so a
!  converted call keeps its buffer, count and type visibly unchanged. The
!  serial stubs do what a one-rank MPI run does.
module turbogap_comm

   use kinds, only: dp
#ifdef _MPIF90
   use mpi
#endif

   implicit none

   private
   public :: comm_t
   public :: comm_init
   public :: comm_finalize
   public :: comm_bcast
   public :: comm_sum_to_root
   public :: comm_sum_all
   public :: comm_allgather
   public :: comm_with_mpi

#ifdef _MPIF90
   logical, parameter :: comm_with_mpi = .true.
#else
   logical, parameter :: comm_with_mpi = .false.
#endif

   type comm_t
      integer :: handle = 0
      integer :: rank = 0
      integer :: size = 1
      logical :: is_root = .true.
   end type comm_t

   interface comm_bcast
      module procedure bcast_l0
      module procedure bcast_i0
      module procedure bcast_r0
      module procedure bcast_i1
      module procedure bcast_r1
      module procedure bcast_r2
      module procedure bcast_l2
      module procedure bcast_c1
   end interface comm_bcast

!  mpi_reduce to rank 0; recv is undefined on the other ranks.
   interface comm_sum_to_root
      module procedure sum_to_root_i1
      module procedure sum_to_root_r1
      module procedure sum_to_root_r2
      module procedure sum_to_root_r3
   end interface comm_sum_to_root

!  mpi_allreduce in place.
   interface comm_sum_all
      module procedure sum_all_r0
      module procedure sum_all_r2
      module procedure sum_all_r3
   end interface comm_sum_all

   interface comm_allgather
      module procedure allgather_i0
   end interface comm_allgather

contains

   subroutine comm_init(comm)
      type(comm_t), intent(out) :: comm
      integer :: ierr
#ifdef _MPIF90
      call mpi_init(ierr)
      comm%handle = MPI_COMM_WORLD
      call mpi_comm_size(comm%handle, comm%size, ierr)
      call mpi_comm_rank(comm%handle, comm%rank, ierr)
#else
      ierr = 0
#endif
      comm%is_root = (comm%rank == 0)
   end subroutine comm_init

   subroutine comm_finalize(comm)
      type(comm_t), intent(in) :: comm
      integer :: ierr
#ifdef _MPIF90
      call mpi_finalize(ierr)
#else
      ierr = comm%rank
#endif
   end subroutine comm_finalize

   subroutine bcast_l0(comm, x)
      type(comm_t), intent(in) :: comm
      logical, intent(inout) :: x
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, 1, MPI_LOGICAL, 0, comm%handle, ierr)
#endif
   end subroutine bcast_l0

   subroutine bcast_i0(comm, x)
      type(comm_t), intent(in) :: comm
      integer, intent(inout) :: x
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, 1, MPI_INTEGER, 0, comm%handle, ierr)
#endif
   end subroutine bcast_i0

   subroutine bcast_r0(comm, x)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: x
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, 1, MPI_DOUBLE_PRECISION, 0, comm%handle, ierr)
#endif
   end subroutine bcast_r0

   subroutine bcast_i1(comm, x, n)
      type(comm_t), intent(in) :: comm
      integer, contiguous, intent(inout) :: x(:)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, n, MPI_INTEGER, 0, comm%handle, ierr)
#endif
   end subroutine bcast_i1

   subroutine bcast_r1(comm, x, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(inout) :: x(:)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, n, MPI_DOUBLE_PRECISION, 0, comm%handle, ierr)
#endif
   end subroutine bcast_r1

   subroutine bcast_r2(comm, x, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(inout) :: x(:, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, n, MPI_DOUBLE_PRECISION, 0, comm%handle, ierr)
#endif
   end subroutine bcast_r2

   subroutine bcast_l2(comm, x, n)
      type(comm_t), intent(in) :: comm
      logical, contiguous, intent(inout) :: x(:, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, n, MPI_LOGICAL, 0, comm%handle, ierr)
#endif
   end subroutine bcast_l2

!  n counts characters, not strings.
   subroutine bcast_c1(comm, x, n)
      type(comm_t), intent(in) :: comm
      character(len=*), contiguous, intent(inout) :: x(:)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_bcast(x, n, MPI_CHARACTER, 0, comm%handle, ierr)
#endif
   end subroutine bcast_c1

   subroutine sum_to_root_i1(comm, send, recv, n)
      type(comm_t), intent(in) :: comm
      integer, contiguous, intent(in) :: send(:)
      integer, contiguous, intent(inout) :: recv(:)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_reduce(send, recv, n, MPI_INTEGER, MPI_SUM, 0, comm%handle, ierr)
#else
      call copy_i(send, recv, n)
#endif
   end subroutine sum_to_root_i1

   subroutine sum_to_root_r1(comm, send, recv, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:)
      real(dp), contiguous, intent(inout) :: recv(:)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_reduce(send, recv, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm%handle, ierr)
#else
      call copy_r(send, recv, n)
#endif
   end subroutine sum_to_root_r1

   subroutine sum_to_root_r2(comm, send, recv, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:, :)
      real(dp), contiguous, intent(inout) :: recv(:, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_reduce(send, recv, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm%handle, ierr)
#else
      call copy_r(send, recv, n)
#endif
   end subroutine sum_to_root_r2

   subroutine sum_to_root_r3(comm, send, recv, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:, :, :)
      real(dp), contiguous, intent(inout) :: recv(:, :, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_reduce(send, recv, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm%handle, ierr)
#else
      call copy_r(send, recv, n)
#endif
   end subroutine sum_to_root_r3

   subroutine sum_all_r0(comm, x)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: x
      integer :: ierr
#ifdef _MPIF90
      call mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm%handle, ierr)
#endif
   end subroutine sum_all_r0

   subroutine sum_all_r2(comm, x, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(inout) :: x(:, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_allreduce(MPI_IN_PLACE, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm%handle, ierr)
#endif
   end subroutine sum_all_r2

   subroutine sum_all_r3(comm, x, n)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(inout) :: x(:, :, :)
      integer, intent(in) :: n
      integer :: ierr
#ifdef _MPIF90
      call mpi_allreduce(MPI_IN_PLACE, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm%handle, ierr)
#endif
   end subroutine sum_all_r3

   subroutine allgather_i0(comm, x, all)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: x
      integer, contiguous, intent(inout) :: all(:)
      integer :: ierr
#ifdef _MPIF90
      call mpi_allgather(x, 1, MPI_INTEGER, all, 1, MPI_INTEGER, comm%handle, ierr)
#else
      all(1) = x
#endif
   end subroutine allgather_i0

!  Element-sequence copies for the serial stubs, whatever the rank.
   subroutine copy_i(send, recv, n)
      integer, intent(in) :: send(*)
      integer, intent(inout) :: recv(*)
      integer, intent(in) :: n
      recv(1:n) = send(1:n)
   end subroutine copy_i

   subroutine copy_r(send, recv, n)
      real(dp), intent(in) :: send(*)
      real(dp), intent(inout) :: recv(*)
      integer, intent(in) :: n
      recv(1:n) = send(1:n)
   end subroutine copy_r

end module turbogap_comm
