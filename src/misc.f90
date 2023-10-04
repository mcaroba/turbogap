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

  !use psb_base_mod
  !use psb_prec_mod
  !use psb_krylov_mod
  !use psb_util_mod

  contains

  subroutine integrate(method, x, y, a, b, integral)

    implicit none

!   Input variables
    real*8, intent(in) :: y(:), x(:), a, b
    character(len=*), intent(in) :: method
!   Output variables
    real*8, intent(out) :: integral

!   Internal variables
    real*8, allocatable :: yc(:), xc(:)
    real*8 :: x_start, x_end, f, x_temp, y_temp
    integer :: n, i, j

    if( size(x) /= size(y) )then
      write(*,*) "ERROR in integrate subroutine: x and y should have the same dimension!"
      stop
    end if

    if( b == a )then
      integral = 0.d0
      return
    else if( b < a )then
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
      if( x(i) >= x_start .and. x(i) <= x_end )then
        n = n + 1
      end if
    end do

    if( n < 2 )then
      write(*,*) "ERROR in integrate subroutine: not enough samples within integration domain!"
      stop
    end if

    allocate( xc(1:n) )
    allocate( yc(1:n) )

    n = 0
    do i = 1, size(x)
      if( x(i) >= x_start .and. x(i) <= x_end )then
        n = n + 1
        xc(n) = x(i)
        yc(n) = y(i)
      end if
    end do
!   Order the arrays
    do i = 1, n
      j = i - 1 + minloc(xc(i:n),1)
      x_temp = xc(j)
      y_temp = yc(j)
      xc(j) = xc(i)
      yc(j) = yc(i)
      xc(i) = x_temp
      yc(i) = y_temp
    end do

    integral = 0.d0
    if( method == "trapezoidal" )then
!     Implements trapezoidal rule
!     If the domain walls are outside of the integration domain we extrapolate
      integral = yc(1) * (xc(1)-x_start) - 0.5d0 * (yc(2)-yc(1)) * (xc(1)-x_start)**2 / (xc(2)-xc(1))
      integral = integral + yc(n) * (x_end-xc(n)) + 0.5d0 * (yc(n)-yc(n-1)) * (x_end - xc(n))**2 / (xc(n)-xc(n-1))
      do i = 1, n-1
        integral = integral + 0.5d0 * (yc(i)+yc(i+1)) * (xc(i+1)-xc(i))
      end do
    else if( method == "simpson" )then
!     Implements Simpson's rule
      write(*,*) "ERROR in integrate subroutine: method", method, "not implemented!"
      stop
    else
      write(*,*) "ERROR in integrate subroutine: method", method, "not implemented!"
      stop
    end if

    deallocate( xc, yc )

!   Fix the sign
    integral = f * integral

  end subroutine
  
  subroutine sparse_mul(A,v,dim,iA,jA,res)

  implicit none

  real*8, intent(in) :: A(:), v(:)
  integer, intent(in) :: dim
  real*8, intent(inout) :: res(:)
  integer*8, intent(in) :: iA(:), jA(:)
  real*8, allocatable :: resd(:), ress(:)
  real*8, allocatable :: A_CSR(:), this_v(:)
  integer, allocatable :: j_CSR(:), i_CSR(:), ipointer(:)
  real*8 :: t1, t2, td, ts, tcsr, s
  !integer, allocatable :: iA(:), jA(:)
  integer :: i, j, k, n, dim1, dim2, max_el, m

  !dim1 = size(A_dense,1)
  !dim2 = size(A_dense,2)
  dim1 = dim
  dim2 = dim

  n = size(iA,1)

  !allocate( res(1:dim2) )
  allocate( resd(1:dim2) )
  allocate( ress(1:dim2) )

  allocate( ipointer(1:dim1) )

  allocate( A_CSR(1:n) )
  allocate( j_CSR(1:n) )
  allocate( i_CSR(1:dim1+1) )

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


  !call cpu_time(t1)
  ipointer = 0
! Count how many elements in each row
  do k = 1, n
    ipointer(iA(k)) = ipointer(iA(k)) + 1
  end do
! Set a pointer to where in the A(1:n) array the elements of i start
  k = ipointer(1)
  ipointer(1) = 1
  do i = 1, dim1-1
    j = ipointer(i+1)
    ipointer(i+1) = ipointer(i) + k
    k = j
  end do
  do i = 1, dim1
    i_CSR(i) = ipointer(i)
  end do
  i_CSR(dim1+1) = n + 1

! Start to populate A_CSR, i_CSR and j_CSR
  ipointer = i_CSR(1:dim1)
  do k = 1, n
    i = iA(k)
    m = ipointer(i)
    A_CSR(m) = A(k)
    j_CSR(m) = jA(k)
    ipointer(i) = m + 1
  end do
  !call cpu_time(t2)
  !write(*,*) t2-t1, "seconds to make A_CSR from sparse matrix"








! Dense A calculation
!  call cpu_time(t1)
!  res = matmul(A_dense, v)
!  call cpu_time(t2)
!  write(*, *) "dense matmul in", t2-t1, "seconds"

  resd = res

  !td=t2-t1

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

  ress = res

  !ts=t2-t1

! Sparse A calculation with A expressed using CSR convention (contiguous memory access)
  !call cpu_time(t1)
  res = 0.d0
  do i = 1, dim1
    s = 0.d0
    do k = i_CSR(i), i_CSR(i+1)-1
      j = j_CSR(k)
      s = s + A_CSR(k)*v(j)
    end do
    res(i) = s
  end do
  !call cpu_time(t2)
  !write(*, *) "CSR sparse matmul in", t2-t1, "seconds"

  !tcsr=t2-t1

  !write(*,*) "Timing: "
  !write(*,*) ts, "(naive sparse)"
  !write(*,*) td, "(dense)"
  !write(*,*) tcsr, "(CSR)"
  !write(*,*) "matmul diff", maxval(abs(res-resd))
  !open(unit=10, file="dat.dat", status="unknown")
  !do i = 1, dim2
  !  write(*,*) resd(i), ress(i), res(i)
  !end do

  deallocate( resd, ress, ipointer, A_CSR, j_CSR, i_CSR )


  end subroutine



  subroutine power_iteration(val, ia, ja, dim, n_iter, b) !, myidx, nnz, n_iter, b)

    implicit none

!   Input variables
    !real(psb_dpk_), 
    real*8, intent(in) :: val(:)
    !integer(psb_lpk_),
    integer*8, intent(in) :: ia(:), ja(:) !, myidx(:)
    integer, intent(in) :: dim !, nnz
    integer, intent(in) :: n_iter
!   Output variables
    real*8, intent(inout) :: b(:)
!   Internal variables
    integer :: k, N
    real*8, allocatable :: b_k(:), b_k1(:)
    real*8 :: b_k1_norm, time1, time2
    
    !type(psb_ctxt_type) :: icontxt
    !integer(psb_ipk_) ::  iam, np, ip, jp, idummy, nr, info_psb
    !type(psb_desc_type) :: desc_a
    !type(psb_dspmat_type) :: A_sp

    N = size(b,1)
    allocate( b_k(1:N) )
    allocate( b_k1(1:N) )

!    call random_number( b_k )
    b_k = 0.5d0

      !call psb_init(icontxt)
      !call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
      !call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
      !call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
      !call psb_cdasb(desc_a, info_psb)
      !call psb_spasb(A_sp, desc_a, info_psb)

    do k = 1, n_iter

      !call psb_spmm(1.d0, A_sp, b_k, 0.d0, b_k1, desc_a, info_psb, 'N')
      call sparse_mul(val, b_k, dim, ia, ja, b_k1)
      !call dgemm('N', 'N', size(A,1), 1, size(A,2), 1.d0, A, size(A,1), b_k, size(A,2), 0.d0, b_k1, size(A,1))
      b_k1_norm = sqrt(dot_product(b_k1,b_k1))
      b_k = b_k1 / b_k1_norm
    end do

    b = b_k
    deallocate( b_k )

  end subroutine

end module
