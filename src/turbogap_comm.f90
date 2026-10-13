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
   public :: comm_split
   public :: comm_free
   public :: comm_cart_create
   public :: comm_cart_shift
   public :: comm_cart_coords
   public :: comm_cart_rank
   public :: comm_sendrecv
   public :: comm_scatterv
   public :: comm_gatherv
   public :: comm_gather

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

   interface comm_sendrecv
      module procedure sendrecv_i0
      module procedure sendrecv_r1
   end interface comm_sendrecv

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

!  A child communicator per colour, ranked by key. A negative colour leaves
!  this rank out: it gets a communicator of size 0 that nothing may use.
   subroutine comm_split(parent, color, key, child)
      type(comm_t), intent(in) :: parent
      integer, intent(in) :: color
      integer, intent(in) :: key
      type(comm_t), intent(out) :: child
      integer :: mpi_color
      integer :: ierr
#ifdef _MPIF90
      mpi_color = color
      if (color < 0) mpi_color = MPI_UNDEFINED
      call mpi_comm_split(parent%handle, mpi_color, key, child%handle, ierr)
      if (child%handle == MPI_COMM_NULL) then
         child%rank = -1
         child%size = 0
      else
         call mpi_comm_size(child%handle, child%size, ierr)
         call mpi_comm_rank(child%handle, child%rank, ierr)
      end if
#else
      mpi_color = color
      ierr = key
      child = parent
      if (mpi_color < 0) then
         child%rank = -1
         child%size = 0
      end if
#endif
      child%is_root = (child%rank == 0)
   end subroutine comm_split

   subroutine comm_free(comm)
      type(comm_t), intent(inout) :: comm
      integer :: ierr
#ifdef _MPIF90
      if (comm%size > 0 .and. comm%handle /= MPI_COMM_WORLD) call mpi_comm_free(comm%handle, ierr)
#else
      ierr = 0
#endif
      comm%size = 0
      comm%rank = -1
      comm%is_root = .false.
   end subroutine comm_free

!  A periodic 3-D grid over parent's ranks, in parent's rank order.
   subroutine comm_cart_create(parent, dims, cart)
      type(comm_t), intent(in) :: parent
      integer, intent(in) :: dims(1:3)
      type(comm_t), intent(out) :: cart
      logical :: periods(1:3)
      integer :: ierr
      periods = .true.
#ifdef _MPIF90
      call mpi_cart_create(parent%handle, 3, dims, periods, .false., cart%handle, ierr)
      call mpi_comm_size(cart%handle, cart%size, ierr)
      call mpi_comm_rank(cart%handle, cart%rank, ierr)
#else
      ierr = dims(1)
      cart = parent
#endif
      cart%is_root = (cart%rank == 0)
   end subroutine comm_cart_create

!  axis counts from 1, as Fortran does; MPI counts from 0.
   subroutine comm_cart_shift(cart, axis, lo, hi)
      type(comm_t), intent(in) :: cart
      integer, intent(in) :: axis
      integer, intent(out) :: lo
      integer, intent(out) :: hi
      integer :: ierr
#ifdef _MPIF90
      call mpi_cart_shift(cart%handle, axis - 1, 1, lo, hi, ierr)
#else
      ierr = axis
      lo = cart%rank
      hi = cart%rank
#endif
   end subroutine comm_cart_shift

!  Grid coordinates count from 0.
   subroutine comm_cart_coords(cart, rank, coords)
      type(comm_t), intent(in) :: cart
      integer, intent(in) :: rank
      integer, intent(out) :: coords(1:3)
      integer :: ierr
#ifdef _MPIF90
      call mpi_cart_coords(cart%handle, rank, 3, coords, ierr)
#else
      ierr = rank
      coords = 0
#endif
   end subroutine comm_cart_coords

   subroutine comm_cart_rank(cart, coords, rank)
      type(comm_t), intent(in) :: cart
      integer, intent(in) :: coords(1:3)
      integer, intent(out) :: rank
      integer :: ierr
#ifdef _MPIF90
      call mpi_cart_rank(cart%handle, coords, rank, ierr)
#else
      ierr = coords(1)
      rank = 0
#endif
   end subroutine comm_cart_rank

   subroutine sendrecv_i0(comm, send, dest, recv, source)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: send
      integer, intent(in) :: dest
      integer, intent(out) :: recv
      integer, intent(in) :: source
      integer :: ierr
#ifdef _MPIF90
      call mpi_sendrecv(send, 1, MPI_INTEGER, dest, 0, recv, 1, MPI_INTEGER, source, 0, &
                        comm%handle, MPI_STATUS_IGNORE, ierr)
#else
      ierr = dest + source
      recv = send
#endif
   end subroutine sendrecv_i0

   subroutine sendrecv_r1(comm, send, n_send, dest, recv, n_recv, source)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:)
      integer, intent(in) :: n_send
      integer, intent(in) :: dest
      real(dp), contiguous, intent(inout) :: recv(:)
      integer, intent(in) :: n_recv
      integer, intent(in) :: source
      integer :: ierr
#ifdef _MPIF90
      call mpi_sendrecv(send, n_send, MPI_DOUBLE_PRECISION, dest, 1, recv, n_recv, MPI_DOUBLE_PRECISION, source, 1, &
                        comm%handle, MPI_STATUS_IGNORE, ierr)
#else
      ierr = dest + source + n_recv
      call copy_r(send, recv, n_send)
#endif
   end subroutine sendrecv_r1

!  Root sends counts(i) elements from offset displs(i) to rank i-1.
   subroutine comm_scatterv(comm, send, counts, displs, recv, n_recv)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:)
      integer, contiguous, intent(in) :: counts(:)
      integer, contiguous, intent(in) :: displs(:)
      real(dp), contiguous, intent(inout) :: recv(:)
      integer, intent(in) :: n_recv
      integer :: ierr
#ifdef _MPIF90
      call mpi_scatterv(send, counts, displs, MPI_DOUBLE_PRECISION, recv, n_recv, MPI_DOUBLE_PRECISION, 0, &
                        comm%handle, ierr)
#else
      ierr = displs(1)
      call copy_r(send, recv, counts(1))
#endif
   end subroutine comm_scatterv

!  The reverse: rank i-1 sends n_send elements, root stores them at displs(i).
   subroutine comm_gatherv(comm, send, n_send, recv, counts, displs)
      type(comm_t), intent(in) :: comm
      real(dp), contiguous, intent(in) :: send(:)
      integer, intent(in) :: n_send
      real(dp), contiguous, intent(inout) :: recv(:)
      integer, contiguous, intent(in) :: counts(:)
      integer, contiguous, intent(in) :: displs(:)
      integer :: ierr
#ifdef _MPIF90
      call mpi_gatherv(send, n_send, MPI_DOUBLE_PRECISION, recv, counts, displs, MPI_DOUBLE_PRECISION, 0, &
                       comm%handle, ierr)
#else
      ierr = displs(1) + counts(1)
      call copy_r(send, recv, n_send)
#endif
   end subroutine comm_gatherv

   subroutine comm_gather(comm, x, all)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: x
      integer, contiguous, intent(inout) :: all(:)
      integer :: ierr
#ifdef _MPIF90
      call mpi_gather(x, 1, MPI_INTEGER, all, 1, MPI_INTEGER, 0, comm%handle, ierr)
#else
      ierr = 0
      all(1) = x
#endif
   end subroutine comm_gather

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
