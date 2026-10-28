! Helpers for managing the ever-growing sparse set of the metadynamics bias

module grow_array
   use kinds
   implicit none

!  Smallest capacity ever allocated, so that the first few appends do not each
!  trigger a reallocation.
   integer, parameter, private :: min_capacity = 16

contains

!  Make sure arr(:) holds at least "needed" elements.
   subroutine ensure_capacity_1d(arr, needed)
      real(dp), allocatable, intent(inout) :: arr(:)
      integer, intent(in) :: needed
      real(dp), allocatable :: tmp(:)

      if (.not. allocated(arr)) then
         allocate (arr(max(needed, min_capacity)))
         return
      end if
      if (needed <= size(arr)) return

      allocate (tmp(max(needed, 2*size(arr))))
      tmp(1:size(arr)) = arr
      call move_alloc(tmp, arr)
   end subroutine ensure_capacity_1d

!  Integer version of ensure_capacity_1d.
   subroutine ensure_capacity_1d_int(arr, needed)
      integer, allocatable, intent(inout) :: arr(:)
      integer, intent(in) :: needed
      integer, allocatable :: tmp(:)

      if (.not. allocated(arr)) then
         allocate (arr(max(needed, min_capacity)))
         return
      end if
      if (needed <= size(arr)) return

      allocate (tmp(max(needed, 2*size(arr))))
      tmp(1:size(arr)) = arr
      call move_alloc(tmp, arr)
   end subroutine ensure_capacity_1d_int

   subroutine ensure_capacity_2d(arr, dim, needed)
      real(dp), allocatable, intent(inout) :: arr(:, :)
      integer, intent(in) :: dim, needed
      real(dp), allocatable :: tmp(:, :)
      integer :: old_cap

      if (.not. allocated(arr)) then
         allocate (arr(dim, max(needed, min_capacity)))
         return
      end if
      if (size(arr, 1) /= dim) error stop "ensure_capacity_2d: vector length mismatch"

      old_cap = size(arr, 2)
      if (needed <= old_cap) return

      allocate (tmp(dim, max(needed, 2*old_cap)))
      tmp(:, 1:old_cap) = arr          ! contiguous block copy
      call move_alloc(tmp, arr)
   end subroutine ensure_capacity_2d


   subroutine ensure_capacity_square(arr, needed)
      real(dp), allocatable, intent(inout) :: arr(:, :)
      integer, intent(in) :: needed
      real(dp), allocatable :: tmp(:, :)
      integer :: old_cap

      if (.not. allocated(arr)) then
         allocate (arr(max(needed, min_capacity), max(needed, min_capacity)))
         return
      end if

      old_cap = size(arr, 1)
      if (needed <= old_cap) return

      allocate (tmp(max(needed, 2*old_cap), max(needed, 2*old_cap)))
      tmp(1:old_cap, 1:old_cap) = arr
      call move_alloc(tmp, arr)
   end subroutine ensure_capacity_square

!  Append one column vector v to arr(:, 1:n)
   subroutine push_vector(arr, n, v)
      real(dp), allocatable, intent(inout) :: arr(:, :)
      integer, intent(inout) :: n
      real(dp), intent(in) :: v(:)

      call ensure_capacity_2d(arr, size(v), n + 1)
      n = n + 1
      arr(:, n) = v
   end subroutine push_vector

!  Append several column vectors
   subroutine push_vectors(arr, n, vs)
      real(dp), allocatable, intent(inout) :: arr(:, :)
      integer, intent(inout) :: n
      real(dp), intent(in) :: vs(:, :)
      integer :: k

      k = size(vs, 2)
      call ensure_capacity_2d(arr, size(vs, 1), n + k)
      arr(:, n + 1:n + k) = vs
      n = n + k
   end subroutine push_vectors

!  Append the elements of x to arr(1:n).
   subroutine push_real_elements(arr, n, x)
      real(dp), allocatable, intent(inout) :: arr(:)
      integer, intent(inout) :: n
      real(dp), intent(in) :: x(:)

      call ensure_capacity_1d(arr, n + size(x))
      arr(n + 1:n + size(x)) = x
      n = n + size(x)
   end subroutine push_real_elements

!  Append the elements of x to arr(1:n).
   subroutine push_integer_elements(arr, n, x)
      integer, allocatable, intent(inout) :: arr(:)
      integer, intent(inout) :: n
      integer, intent(in) :: x(:)

      call ensure_capacity_1d_int(arr, n + size(x))
      arr(n + 1:n + size(x)) = x
      n = n + size(x)
   end subroutine push_integer_elements

end module grow_array

