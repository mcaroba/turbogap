
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mpi_utils.f90, is copyright (c) 2019-2025, Miguel A. Caro and
! HND X   Tigany Zarrouk
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
module mpi_utils
#ifdef _MPIF90
   use mpi
#endif
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   interface synchronize_array
      module procedure synchronize_array_int_1
      module procedure synchronize_array_dp_1
      module procedure synchronize_array_dp_2
      module procedure synchronize_array_char_1
      module procedure synchronize_array_logical_1
      module procedure synchronize_array_logical_2
   end interface synchronize_array

contains

   ! subroutine compare_arrays_from_rank( array, rank )
   !   real( dp ), intent(in) :: array(:,:)
   !   integer, intent(in) :: rank
   !   integer :: n_1, n_2
   !   real( dp ), allocatable :: sent_array(:, :)
   !

   subroutine synchronize_array_char_1(array, rank)
      character*8, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, 8*n1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_char_1

   subroutine synchronize_array_logical_1(array, rank)
      logical, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_logical_1

   subroutine synchronize_array_logical_2(array, rank)
      logical, allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1, n2

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
            n2 = size(array, 2)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_bcast(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1 .or. size(array, 2) /= n2) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1, n2))
            end if
         end if

         call MPI_bcast(array, n1*n2, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_logical_2

   subroutine synchronize_array_int_1(array, rank)
      integer, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_int_1

   subroutine synchronize_array_dp_1(array, rank)
      real(dp), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_dp_1

   subroutine synchronize_array_dp_2(array, rank)
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1, n2

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
            n2 = size(array, 2)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_bcast(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1 .or. size(array, 2) /= n2) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1, n2))
            end if
         end if

         call MPI_bcast(array, n1*n2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      end if
#endif

   end subroutine synchronize_array_dp_2
end module mpi_utils
