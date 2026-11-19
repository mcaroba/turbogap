! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, misc.f90, is copyright (c) 2021, Miguel A. Caro
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

module misc

   use kinds

   !use psb_base_mod
   !use psb_prec_mod
   !use psb_krylov_mod
   !use psb_util_mod

!  The CSR layout sparse_mul builds on every call, built once for a pattern that
!  is multiplied many times with different values or vectors: the k-th stored
!  element of the CSR form is A(perm(k)), in sparse_mul's order.
   type :: csr_index_t
      integer :: n = 0
      integer :: dim = 0
      integer, allocatable :: perm(:)
      integer, allocatable :: j_csr(:)
      integer, allocatable :: i_csr(:)
   end type csr_index_t

contains

   subroutine integrate(method, x, y, a, b, integral)

      implicit none

      real(dp), intent(in) :: y(:), x(:), a, b
      character(len=*), intent(in) :: method
      real(dp), intent(out) :: integral

      real(dp), allocatable :: yc(:), xc(:)
      real(dp) :: x_start, x_end, f, x_temp, y_temp
      integer :: n, i, j

      if (size(x) /= size(y)) then
         write (*, *) "ERROR in integrate subroutine: x and y should have the same dimension!"
         stop
      end if

      if (b == a) then
         integral = 0.d0
         return
      else if (b < a) then
         x_start = b
         x_end = a
         f = -1.d0
      else
         x_start = a
         x_end = b
         f = 1.d0
      end if

      n = 0
      do i = 1, size(x)
         if (x(i) >= x_start .and. x(i) <= x_end) then
            n = n + 1
         end if
      end do

      if (n < 2) then
         write (*, *) "ERROR in integrate subroutine: not enough samples within integration domain!"
         stop
      end if

      allocate (xc(1:n))
      allocate (yc(1:n))

      n = 0
      do i = 1, size(x)
         if (x(i) >= x_start .and. x(i) <= x_end) then
            n = n + 1
            xc(n) = x(i)
            yc(n) = y(i)
         end if
      end do
!   Order the arrays
      do i = 1, n
         j = i - 1 + minloc(xc(i:n), 1)
         x_temp = xc(j)
         y_temp = yc(j)
         xc(j) = xc(i)
         yc(j) = yc(i)
         xc(i) = x_temp
         yc(i) = y_temp
      end do

      integral = 0.d0
      if (method == "trapezoidal") then
!     Implements trapezoidal rule
!     If the domain walls are outside of the integration domain we extrapolate
         integral = yc(1)*(xc(1) - x_start) - 0.5d0*(yc(2) - yc(1))*(xc(1) - x_start)**2/(xc(2) - xc(1))
         integral = integral + yc(n)*(x_end - xc(n)) + 0.5d0*(yc(n) - yc(n - 1))*(x_end - xc(n))**2/(xc(n) - xc(n - 1))
         do i = 1, n - 1
            integral = integral + 0.5d0*(yc(i) + yc(i + 1))*(xc(i + 1) - xc(i))
         end do
      else if (method == "simpson") then
!     Implements Simpson's rule
         write (*, *) "ERROR in integrate subroutine: method", method, "not implemented!"
         stop
      else
         write (*, *) "ERROR in integrate subroutine: method", method, "not implemented!"
         stop
      end if

      deallocate (xc, yc)

!   Fix the sign
      integral = f*integral

   end subroutine

   subroutine sparse_mul(A, v, dim, iA, jA, res)

      implicit none

      real(dp), intent(in) :: A(:), v(:)
      integer, intent(in) :: dim
      real(dp), intent(inout) :: res(:)
      integer(i64), intent(in) :: iA(:), jA(:)
      !real(dp), allocatable :: resd(:), ress(:)
      real(dp), allocatable :: A_CSR(:), this_v(:)
      integer, allocatable :: j_CSR(:), i_CSR(:), ipointer(:)
      real(dp) :: t1, t2, td, ts, tcsr, s
      !integer, allocatable :: iA(:), jA(:)
      integer :: i, j, k, n, dim1, dim2, max_el, m

      !dim1 = size(A_dense,1)
      !dim2 = size(A_dense,2)
      dim1 = dim
      dim2 = dim

      n = size(iA, 1)

      allocate (ipointer(1:dim1))

      allocate (A_CSR(1:n))
      allocate (j_CSR(1:n))
      allocate (i_CSR(1:dim1 + 1))

      !call cpu_time(t1)
! Make CSR version of A (this is very inefficient; in an actual production code it would be optimized)
      !max_el = 0
      !k = 1
      !do i = 1, dim1
      !  i_CSR(i) = k
      !  m = 0
      !  do j = 1, dim2
      !    if( A_dense(i, j) /= 0.d0 )then
      !      A_CSR(k) = A_dense(i, j)
      !      j_CSR(k) = j
      !      k = k + 1
      !      m = m + 1
      !    end if
      !  end do
      !  if( m > max_el )then
      !    max_el = m
      !  end if
      !end do
      !i_CSR(dim1+1) = k
      !call cpu_time(t2)
      !write(*,*) t2-t1, "seconds to make A_CSR from dense matrix"

      ipointer = 0
! Count how many elements in each row
      do k = 1, n
         ipointer(iA(k)) = ipointer(iA(k)) + 1
      end do
! Set a pointer to where in the A(1:n) array the elements of i start
      k = ipointer(1)
      ipointer(1) = 1
      do i = 1, dim1 - 1
         j = ipointer(i + 1)
         ipointer(i + 1) = ipointer(i) + k
         k = j
      end do
      do i = 1, dim1
         i_CSR(i) = ipointer(i)
      end do
      i_CSR(dim1 + 1) = n + 1

! Start to populate A_CSR, i_CSR and j_CSR
      ipointer = i_CSR(1:dim1)
      do k = 1, n
         i = iA(k)
         m = ipointer(i)
         A_CSR(m) = A(k)
         j_CSR(m) = jA(k)
         ipointer(i) = m + 1
      end do

! Dense A calculation
!  call cpu_time(t1)
!  res = matmul(A_dense, v)
!  call cpu_time(t2)
!  write(*, *) "dense matmul in", t2-t1, "seconds"

      !resd = res

! Sparse A calculation with randomly ordered elements in A
!  call cpu_time(t1)
!  res = 0.d0
!  do k = 1, n
!    i = iA(k)
!    j = jA(k)
!    res(i) = res(i) + A(k)*v(j)
!  end do
!  call cpu_time(t2)
!  write(*, *) "naive sparse matmul in", t2-t1, "seconds"

      !ress = res

! Sparse A calculation with A expressed using CSR convention (contiguous memory access)
      !call cpu_time(t1)
      res = 0.d0
      do i = 1, dim1
         s = 0.d0
         do k = i_CSR(i), i_CSR(i + 1) - 1
            j = j_CSR(k)
            s = s + A_CSR(k)*v(j)
         end do
         res(i) = s
      end do

      !write(*,*) "Timing: "
      !write(*,*) ts, "(naive sparse)"
      !write(*,*) td, "(dense)"
      !write(*,*) tcsr, "(CSR)"
      !write(*,*) "matmul diff", maxval(abs(res-resd))
      !open(unit=10, file="dat.dat", status="unknown")
      !do i = 1, dim2
      !  write(*,*) resd(i), ress(i), res(i)
      !end do

      deallocate (ipointer, A_CSR, j_CSR, i_CSR)

   end subroutine

!  The CSR index of the pattern (iA, jA), and with diag of that pattern followed
!  by the diagonal (p, p), p = 1..dim -- exactly as sparse_mul places them.
   subroutine csr_index_build(iA, jA, dim, idx, diag)

      implicit none

      integer(i64), intent(in) :: iA(:)
      integer(i64), intent(in) :: jA(:)
      integer, intent(in) :: dim
      type(csr_index_t), intent(inout) :: idx
      logical, intent(in), optional :: diag
      integer, allocatable :: ipointer(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: m
      integer :: n0
      logical :: with_diag

      with_diag = .false.
      if (present(diag)) with_diag = diag
      n0 = size(iA, 1)
      idx%dim = dim
      idx%n = n0
      if (with_diag) idx%n = n0 + dim
      if (allocated(idx%perm)) deallocate (idx%perm, idx%j_csr, idx%i_csr)
      allocate (idx%perm(1:idx%n))
      allocate (idx%j_csr(1:idx%n))
      allocate (idx%i_csr(1:dim + 1))
      allocate (ipointer(1:dim))

      ipointer = 0
      do k = 1, idx%n
         ipointer(row(k)) = ipointer(row(k)) + 1
      end do
      k = ipointer(1)
      ipointer(1) = 1
      do i = 1, dim - 1
         j = ipointer(i + 1)
         ipointer(i + 1) = ipointer(i) + k
         k = j
      end do
      do i = 1, dim
         idx%i_csr(i) = ipointer(i)
      end do
      idx%i_csr(dim + 1) = idx%n + 1
      ipointer = idx%i_csr(1:dim)
      do k = 1, idx%n
         i = row(k)
         m = ipointer(i)
         idx%perm(m) = k
         idx%j_csr(m) = col(k)
         ipointer(i) = m + 1
      end do
      deallocate (ipointer)

   contains

      integer function row(k)
         integer, intent(in) :: k
         if (k <= n0) then
            row = int(iA(k))
         else
            row = k - n0
         end if
      end function row

      integer function col(k)
         integer, intent(in) :: k
         if (k <= n0) then
            col = int(jA(k))
         else
            col = k - n0
         end if
      end function col

   end subroutine csr_index_build

!  sparse_mul with the CSR index built beforehand: the same products, summed in
!  the same order.
   subroutine sparse_mul_csr(A, v, idx, res)

      implicit none

      real(dp), intent(in) :: A(:)
      real(dp), intent(in) :: v(:)
      type(csr_index_t), intent(in) :: idx
      real(dp), intent(inout) :: res(:)
      real(dp) :: s
      integer :: i
      integer :: k

      res = 0.d0
      do i = 1, idx%dim
         s = 0.d0
         do k = idx%i_csr(i), idx%i_csr(i + 1) - 1
            s = s + A(idx%perm(k))*v(idx%j_csr(k))
         end do
         res(i) = s
      end do

   end subroutine sparse_mul_csr

!  sparse_mul_csr on the three columns of v, reading A and the index once.
   subroutine sparse_mul3_csr(A, v, idx, res)

      implicit none

      real(dp), intent(in) :: A(:)
      real(dp), intent(in) :: v(:, :)
      type(csr_index_t), intent(in) :: idx
      real(dp), intent(inout) :: res(:, :)
      real(dp) :: a_k
      real(dp) :: s1
      real(dp) :: s2
      real(dp) :: s3
      integer :: i
      integer :: j
      integer :: k

      res(:, 1:3) = 0.d0
      do i = 1, idx%dim
         s1 = 0.d0
         s2 = 0.d0
         s3 = 0.d0
         do k = idx%i_csr(i), idx%i_csr(i + 1) - 1
            a_k = A(idx%perm(k))
            j = idx%j_csr(k)
            s1 = s1 + a_k*v(j, 1)
            s2 = s2 + a_k*v(j, 2)
            s3 = s3 + a_k*v(j, 3)
         end do
         res(i, 1) = s1
         res(i, 2) = s2
         res(i, 3) = s3
      end do

   end subroutine sparse_mul3_csr

!  power_iteration with the CSR index built beforehand.
   subroutine power_iteration_csr(val, idx, n_iter, b)

      implicit none

      real(dp), intent(in) :: val(:)
      type(csr_index_t), intent(in) :: idx
      integer, intent(in) :: n_iter
      real(dp), intent(inout) :: b(:)
      real(dp), allocatable :: b_k(:)
      real(dp), allocatable :: b_k1(:)
      real(dp) :: b_k1_norm
      integer :: k

      allocate (b_k(1:size(b, 1)))
      allocate (b_k1(1:size(b, 1)))
      b_k = 0.5d0
      do k = 1, n_iter
         call sparse_mul_csr(val, b_k, idx, b_k1)
         b_k1_norm = sqrt(dot_product(b_k1, b_k1))
         b_k = b_k1/b_k1_norm
      end do
      b = b_k
      deallocate (b_k, b_k1)

   end subroutine power_iteration_csr

   subroutine power_iteration(val, ia, ja, dim, n_iter, b) !, myidx, nnz, n_iter, b)

      implicit none

      !real(psb_dpk_),
      real(dp), intent(in) :: val(:)
      !integer(psb_lpk_),
      integer(i64), intent(in) :: ia(:), ja(:) !, myidx(:)
      integer, intent(in) :: dim !, nnz
      integer, intent(in) :: n_iter
      real(dp), intent(inout) :: b(:)
      integer :: k, N
      real(dp), allocatable :: b_k(:), b_k1(:)
      real(dp) :: b_k1_norm, time1, time2

      !type(psb_ctxt_type) :: icontxt
      !integer(psb_ipk_) ::  iam, np, ip, jp, idummy, nr, info_psb
      !type(psb_desc_type) :: desc_a
      !type(psb_dspmat_type) :: A_sp

      N = size(b, 1)
      allocate (b_k(1:N))
      allocate (b_k1(1:N))

      b_k = 0.5d0

      do k = 1, n_iter

         call sparse_mul(val, b_k, dim, ia, ja, b_k1)
         b_k1_norm = sqrt(dot_product(b_k1, b_k1))
         b_k = b_k1/b_k1_norm
      end do

      b = b_k
      deallocate (b_k)

   end subroutine

   subroutine lanczos_algorithm(val, ia, ja, dim, n_iter, l_min, l_max) !, myidx, nnz, n_iter, b)

      implicit none

      !real(psb_dpk_),
      real(dp), intent(in) :: val(:)
      !integer(psb_lpk_),
      integer(i64), intent(in) :: ia(:), ja(:) !, myidx(:)
      integer, intent(in) :: dim !, nnz
      integer, intent(in) :: n_iter
      !real(dp), intent(inout) :: b(:)
      real(dp), intent(inout) :: l_min, l_max
      integer :: k, N, info
      real(dp), allocatable :: v(:, :), w(:, :), alpha(:), beta(:), T(:, :), work_arr(:)
      real(dp) :: b_k1_norm, time1, time2

      allocate (w(1:dim, 1:n_iter))
      allocate (v(1:dim, 1:n_iter))
      allocate (alpha(1:n_iter))
      allocate (beta(2:n_iter))
      allocate (T(1:n_iter, 1:n_iter))
      allocate (work_arr(1:4*n_iter))

      call random_number(v(:, 1))
      v(:, 1) = v(:, 1) - 0.5d0
      v(:, 1) = v(:, 1)/norm2(v(:, 1))

      call sparse_mul(val, v(:, 1), dim, ia, ja, w(:, 1))
      alpha(1) = dot_product(v(:, 1), w(:, 1))
      w(:, 1) = w(:, 1) - alpha(1)*v(:, 1)

      do k = 2, n_iter
         beta(k) = norm2(w(:, k - 1))
         if (beta(k) .ne. 0.d0) then
            v(:, k) = w(:, k - 1)/beta(k)
         else
            call random_number(v(:, k))
            v(:, k) = v(:, k) - 0.5d0
            v(:, k) = v(:, k)/norm2(v(:, k))
         end if
         call sparse_mul(val, v(:, k), dim, ia, ja, w(:, k))
         alpha(k) = dot_product(v(:, k), w(:, k))
         w(:, k) = w(:, k) - alpha(k)*v(:, k) - beta(k)*v(:, k - 1)
      end do

      T(1, 1) = alpha(1)
      T(1, 2) = beta(2)

      do k = 2, n_iter - 1
         T(k, k) = alpha(k)
         T(k, k - 1) = beta(k)
         T(k, k + 1) = beta(k + 1)
      end do

      T(n_iter, n_iter) = alpha(n_iter)
      T(n_iter, n_iter - 1) = beta(n_iter)

      call dpteqr('N', n_iter, alpha, beta, 0.d0, 1, work_arr, info)

      l_min = minval(alpha)
      l_max = maxval(alpha)

      deallocate (w, v, alpha, beta, T, work_arr)

   end subroutine

   function screened_one_over_x(val) result(res)

      implicit none

      real(dp), intent(in) :: val(:)
      real(dp), dimension(1:size(val)) :: res
      real(dp), parameter :: pi = 3.141592653589793d0
      integer :: i

      do i = 1, size(val)
         if (dabs(val(i)) > 1.d0) then
            res(i) = 1.d0/val(i)
         else
            res(i) = dsin(pi/2.d0*val(i))
         end if
      end do

   end function

   function sgn(val) result(res)

      implicit none

      real(dp), intent(in) :: val(:)
      real(dp), dimension(1:size(val)) :: res
      integer :: i

      do i = 1, size(val)
         if (val(i) >= 0.d0) then
            res(i) = 1.d0
         else
            res(i) = -1.d0
         end if
      end do

   end function

   function smooth_ratio(val1, val2, min_val, max_val) result(res)

      implicit none

      real(dp), intent(in) :: val1(:), val2(:), min_val, max_val
      real(dp), dimension(1:size(val1)) :: res
      real(dp) :: r
      integer :: i

      do i = 1, size(val1)
         if (abs(val2(i)) < 1.d-10) then
            res(i) = 0.d0
         else
            r = val1(i)/val2(i)
            if (r < min_val) then
               res(i) = min_val
            else if (r > max_val) then
               res(i) = max_val
            else
               res(i) = r
            end if
         end if
      end do

   end function

   function clip(val, min_val, max_val) result(res)

      implicit none

      real(dp), intent(in) :: val(:), min_val, max_val
      real(dp), dimension(1:size(val)) :: res
      integer :: i

      do i = 1, size(val)
         if (val(i) < min_val) then
            res(i) = min_val
         else if (val(i) > max_val) then
            res(i) = max_val
         else
            res(i) = val(i)
         end if
      end do

   end function

   function poly_cut(x, x_min, x_max) result(res)

      implicit none

      real(dp), intent(in) :: x, x_min, x_max
      real(dp) :: res

      if (x < x_min) then
         res = 0.d0
      else if (x >= x_max) then
         res = 1.d0
      else
         res = 0.5d0 + 1.5d0*(x - (x_max + x_min)/2.d0)/(x_max - x_min) - 2.d0*(x - (x_max + x_min)/2.d0)**3/(x_max - x_min)**3
      end if

   end function

   function std(val) result(res)

      implicit none

      real(dp), intent(in) :: val(:)
      real(dp) :: res, ave
      integer :: i

      ave = sum(val)/size(val)
      res = 0.d0
      do i = 1, size(val)
         res = res + (val(i) - ave)**2
      end do
      res = dsqrt(res/size(val))

   end function

end module
