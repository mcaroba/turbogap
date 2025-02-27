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
  use types
  use iso_c_binding
  use F_B_C
  contains



  subroutine gpu_sparse_mul(A,v,dim,iA,jA,res, gpu_stream_blas)

  implicit none

  real*8, intent(in), target :: A(:), v(:)
  integer, intent(in) :: dim
  real*8, intent(inout), target :: res(:)
  integer*8, intent(in) :: iA(:), jA(:)
  !real*8, allocatable :: resd(:), ress(:)
  real*8, allocatable, target :: A_CSR_temp(:), this_v_temp(:)
  real*8, allocatable, target :: A_CSR(:), this_v(:)
  integer, allocatable, target :: j_CSR(:), i_CSR(:), ipointer(:)
  integer, allocatable, target :: j_CSR_temp(:), i_CSR_temp(:), ipointer_temp(:)
  real*8 :: t1, t2, td, ts, tcsr, s
  !integer, allocatable :: iA(:), jA(:)
  integer :: i, j, k, n, dim1, dim2, max_el, m
  
  type(c_ptr) :: gpu_stream, gpu_stream_blas
  integer(c_size_t) :: st_A_CSR_d
  integer(c_size_t) :: st_j_CSR_d
  integer(c_size_t) :: st_i_CSR_d
  integer(c_size_t) :: st_v_d
  integer(c_size_t) :: st_ipointer_d
  integer(c_size_t) :: st_res_d

  type(c_ptr) :: A_CSR_d
  type(c_ptr) :: j_CSR_d
  type(c_ptr) :: i_CSR_d
  type(c_ptr) :: ipointer_d
  type(c_ptr) :: res_d
  type(c_ptr) :: v_d
  type(c_ptr) :: cusparse_handle 
  !dim1 = size(A_dense,1)
  !dim2 = size(A_dense,2)
  dim1 = dim
  dim2 = dim


  ! Creating the stream here used by cusparse explicitly. 
    call create_cusparse_handle(cusparse_handle, gpu_stream)

  n = size(iA,1)

  allocate( ipointer(1:dim1) )
  allocate( A_CSR(1:n) )
  allocate( j_CSR(1:n) )
  allocate( i_CSR(1:dim1+1) )

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
 

! Zero indexing, does it just work like this?
  i_CSR = i_CSR - 1
  j_CSR = j_CSR - 1
!    deallocate( i_CSR, j_CSR, A_CSR ) 
!    allocate( i_CSR(1:5) )
!    allocate( j_CSR(1:9) )
!    allocate( A_CSR(1:9) )


!    dim1 = 4
!    dim2 = 4 
!    n = 9 
!    i_CSR = (/ 0, 3, 4, 7, 9 /)
!    j_CSR = (/ 0, 2, 3, 1, 0, 2, 3, 1, 3 /)
    !i_CSR = i_CSR + 1
    !j_CSR = j_CSR + 1
!    A_CSR = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, &
!                                  6.d0, 7.d0, 8.d0, 9.d0 /)
!    v = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)
!    res = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
!    hY_result[]     = { 19.0f, 8.0f, 51.0f, 52.0f };
   
  st_A_CSR_d     = n * c_double  
  st_j_CSR_d     = n * c_int  
  st_i_CSR_d     = ( dim1 + 1) * c_int  
  st_v_d         = dim1 * c_double  
  st_ipointer_d  = dim1 * c_int  
  st_res_d       = dim1 * c_double  
 


    call gpu_device_sync()
 
!  print *,"    A_CSR    " , size( A_CSR   ),  A_CSR   , st_A_CSR_d
!  print *,"    j_CSR    " , size( j_CSR   ),  j_CSR   , st_j_CSR_d
!  print *,"    i_CSR    " , size( i_CSR   ),  i_CSR   , st_i_CSR_d
!  print *,"    v        " , size( v       ),  v       , st_v_d
!  print *,"    ipointer " , size( ipointer),  ipointer, st_ipointer_d
!  print *,"    res      " , size( res     ),  res     , st_res_d    


  call gpu_malloc_all_blocking( A_CSR_d, st_A_CSR_d)
  call gpu_malloc_all_blocking( j_CSR_d, st_j_CSR_d)
  call gpu_malloc_all_blocking( i_CSR_d, st_i_CSR_d)
  call gpu_malloc_all_blocking( v_d,     st_v_d)
  call gpu_malloc_all_blocking( res_d, st_res_d)
  
!  print *,"    A_CSR    " , size( A_CSR   ),  A_CSR   , st_A_CSR_d
!  print *,"    j_CSR    " , size( j_CSR   ),  j_CSR   , st_j_CSR_d
!  print *,"    i_CSR    " , size( i_CSR   ),  i_CSR   , st_i_CSR_d
!  print *,"    v        " , size( v       ),  v       , st_v_d
!  print *,"    ipointer " , size( ipointer),  ipointer, st_ipointer_d
!  print *,"    res      " , size( res     ),  res     , st_res_d    




! j
! j  allocate( ipointer_temp(1:dim1) )
! j  allocate( A_CSR_temp(1:n) )
! j  allocate( j_CSR_temp(1:n) )
! j  allocate( i_CSR_temp(1:dim1+1) )
! j  allocate( this_v_temp(dim1))

    call gpu_device_sync()
  call cpy_htod_blocking( c_loc( A_CSR ), A_CSR_d, st_A_CSR_d)
  call cpy_htod_blocking( c_loc( j_CSR ), j_CSR_d, st_j_CSR_d)
  call cpy_htod_blocking( c_loc( i_CSR ), i_CSR_d, st_i_CSR_d)
  call cpy_htod_blocking( c_loc( v ),     v_d,     st_v_d)
!  call cpy_htod( c_loc( res ), res_d, st_res_d, gpu_stream)
  
! j
! j  call cpy_dtoh(  A_CSR_d,    c_loc( A_CSR_temp ),   st_A_CSR_d, gpu_stream)
! j  call cpy_dtoh(  j_CSR_d,    c_loc( j_CSR_temp ),   st_j_CSR_d, gpu_stream)
! j  call cpy_dtoh(  i_CSR_d,    c_loc( i_CSR_temp ),   st_i_CSR_d, gpu_stream)
! j  call cpy_dtoh(  v_d,        c_loc( this_v_temp ),       st_v_d, gpu_stream)
! j
! j 
! j  print *,"    A_CSR    temp copied back " , size( A_CSR   ),  A_CSR_temp
! j  print *,"    j_CSR    temp copied back " , size( j_CSR   ),  j_CSR_temp
! j  print *,"    i_CSR    temp copied back " , size( i_CSR   ),  i_CSR_temp
! j  print *,"    v        temp copied back " , size( v       ),  this_v_temp 
! j  print *,"    ipointer temp copied back " , size( ipointer),  ipointer_temp



    call gpu_device_sync()
  call gpu_memset_blocking(res_d, 0,  st_res_d)

  
  

    call gpu_device_sync()
!  write(*,*) "Before matmul"
  call gpu_sparse_matrix_mul_kernel(dim1, dim2, n, &
          A_CSR_d, &
          j_CSR_d, &
          i_CSR_d, &
          ipointer_d, &
          res_d, &
          v_d, gpu_stream, cusparse_handle )
!  write(*,*) "Before cpy_dtoh"
    call gpu_device_sync()
  ! call gpu_stream_sync(gpu_stream)
  ! call gpu_device_sync()
  call cpy_dtoh_blocking( res_d, c_loc( res ), st_res_d)
  !call cpy_dtoh( c_loc( res ), res_d, st_res_d, gpu_stream)
!  write(*,*) "cpy_dtoh done"

  call gpu_free( A_CSR_d)
  call gpu_free( j_CSR_d)
  call gpu_free( i_CSR_d)
  call gpu_free( v_d)
  call gpu_free( res_d)

   deallocate( ipointer, A_CSR, j_CSR, i_CSR )

    call destroy_cusparse_handle(cusparse_handle, gpu_stream)
  end subroutine



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
  !real*8, allocatable :: resd(:), ress(:)
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
  !allocate( resd(1:dim2) )
  !allocate( ress(1:dim2) )

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

  !resd = res

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

  !ress = res

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

  deallocate( ipointer, A_CSR, j_CSR, i_CSR )

  end subroutine



  subroutine power_iteration(val, ia, ja, dim, n_iter, b, do_gpu) !, myidx, nnz, n_iter, b, do_gpu)

    implicit none

!   Input variables
    !real(psb_dpk_), 
    real*8, intent(in) :: val(:)
    !integer(psb_lpk_),
    integer*8, intent(in) :: ia(:), ja(:) !, myidx(:)
    integer, intent(in) :: dim !, nnz
    integer, intent(in) :: n_iter
    logical, intent(in) :: do_gpu
!   Output variables
    real*8, intent(inout) :: b(:)
!   Internal variables
    integer :: k, N
    real*8, allocatable :: b_k(:), b_k1(:)
    real*8 :: b_k1_norm, time1, time2
    
    type(c_ptr) :: gpu_stream

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
      if ( do_gpu ) then
        call gpu_sparse_mul(val, b_k, dim, ia, ja, b_k1, gpu_stream)
      else
        call sparse_mul(val, b_k, dim, ia, ja, b_k1)
      end if
      !call dgemm('N', 'N', size(A,1), 1, size(A,2), 1.d0, A, size(A,1), b_k, size(A,2), 0.d0, b_k1, size(A,1))
      b_k1_norm = sqrt(dot_product(b_k1,b_k1))
      b_k = b_k1 / b_k1_norm
    end do

    b = b_k
    deallocate( b_k )

  end subroutine


  subroutine lanczos_algorithm(val, ia, ja, dim, n_iter, l_min, l_max, do_gpu) !, myidx, nnz, n_iter, b, do_gpu)

    implicit none

!   Input variables
    !real(psb_dpk_), 
    real*8, intent(in) :: val(:)
    !integer(psb_lpk_),
    integer*8, intent(in) :: ia(:), ja(:) !, myidx(:)
    integer, intent(in) :: dim !, nnz
    integer, intent(in) :: n_iter
    logical, intent(in) :: do_gpu
!   Output variables
    !real*8, intent(inout) :: b(:)
    real*8, intent(inout) :: l_min, l_max
!   Internal variables
    integer :: k, N, info
    real*8, allocatable :: v(:,:), w(:,:), alpha(:), beta(:), T(:,:), work_arr(:)
    real*8 :: b_k1_norm, time1, time2

    type(c_ptr) :: gpu_stream


    allocate( w(1:dim,1:n_iter) )
    allocate( v(1:dim,1:n_iter) )
    allocate( alpha(1:n_iter) )
    allocate( beta(2:n_iter) )
    allocate( T(1:n_iter,1:n_iter) )
    allocate( work_arr(1:4*n_iter) )

    call random_number(v(:,1))
    v(:,1) = v(:,1)-0.5d0
    v(:,1) = v(:,1)/norm2(v(:,1))

    if ( do_gpu ) then
      call gpu_sparse_mul(val, v(:,1), dim, ia, ja, w(:,1), gpu_stream)
    else
      call sparse_mul(val, v(:,1), dim, ia, ja, w(:,1))
    end if
    alpha(1) = dot_product(v(:,1),w(:,1))
    w(:,1) = w(:,1) - alpha(1) * v(:,1)

    do k = 2, n_iter     
      beta(k) = norm2(w(:,k-1))
      if ( beta(k) .ne. 0.d0 ) then
        v(:,k) = w(:,k-1)/beta(k)
      else
        call random_number(v(:,k))
        v(:,k) = v(:,k)-0.5d0
        v(:,k) = v(:,k)/norm2(v(:,k))
      end if
      if ( do_gpu ) then
        call gpu_sparse_mul(val, v(:,k), dim, ia, ja, w(:,k), gpu_stream)
      else
        call sparse_mul(val, v(:,k), dim, ia, ja, w(:,k))
      end if
      alpha(k) = dot_product(v(:,k),w(:,k))
      w(:,k) = w(:,k) - alpha(k) * v(:,k) - beta(k) * v(:,k-1)
    end do

    T(1,1) = alpha(1)
    T(1,2) = beta(2)

    do k = 2, n_iter-1
      T(k,k) = alpha(k)
      T(k,k-1) = beta(k)
      T(k,k+1) = beta(k+1)
    end do

    T(n_iter,n_iter) = alpha(n_iter)
    T(n_iter,n_iter-1) = beta(n_iter)    

    call dpteqr( 'N', n_iter, alpha, beta, 0.d0, 1, work_arr ,info )

    !write(*,*) "dpteqr info", info
    !write(*,*) "Lanczos eigs", alpha-1.d0

    l_min = minval(alpha)
    l_max = maxval(alpha)

    deallocate( w, v, alpha, beta, T, work_arr)


  end subroutine




  function screened_one_over_x(val) result(res)

    implicit none

    real*8, intent(in) :: val(:)
    real*8, dimension(1:size(val)) :: res
    real*8, parameter :: pi = 3.141592653589793d0
    integer :: i

    do i = 1, size(val)
      if( dabs(val(i)) > 1.d0 )then
        res(i) = 1.d0/val(i)
      else
        res(i) = dsin(pi/2.d0*val(i))
      end if
    end do

  end function




  function sgn(val) result(res)

    implicit none

    real*8, intent(in) :: val(:)
    real*8, dimension(1:size(val)) :: res
    integer :: i

    do i = 1, size(val)
      if( val(i) >= 0.d0 )then
        res(i) = 1.d0
      else
        res(i) = -1.d0
      end if
    end do

  end function




  function smooth_ratio(val1, val2, min_val, max_val) result(res)

    implicit none

    real*8, intent(in) :: val1(:), val2(:), min_val, max_val
    real*8, dimension(1:size(val1)) :: res
    real*8 :: r
    integer :: i

    do i = 1, size(val1)
      if( abs(val2(i)) < 1.d-10 )then
        res(i) = 0.d0
      else
        r = val1(i)/val2(i)
        if( r < min_val )then
          res(i) = min_val
        else if( r > max_val )then
       	  res(i) = max_val
        else
          res(i) = r
        end if
      end if
    end do

  end function




  function clip(val, min_val, max_val) result(res)

    implicit none

    real*8, intent(in) :: val(:), min_val, max_val
    real*8, dimension(1:size(val)) :: res
    integer :: i

    do i = 1, size(val)
      if( val(i) < min_val )then
        res(i) = min_val
      else if( val(i) > max_val )then
        res(i) = max_val
      else
        res(i) = val(i)
      end if
    end do

  end function




  function poly_cut(x, x_min, x_max) result(res)

    implicit none

    real*8, intent(in) :: x, x_min, x_max
    real*8 :: res

    if( x < x_min )then
      res = 0.d0
    else if( x >= x_max )then
      res = 1.d0
    else
      res = 0.5d0 + 1.5d0*(x-(x_max+x_min)/2.d0)/(x_max-x_min) - 2.d0*(x-(x_max+x_min)/2.d0)**3/(x_max-x_min)**3
    end if

  end function





  function std(val) result(res)

    implicit none

    real*8, intent(in) :: val(:)
    real*8 :: res, ave
    integer :: i

    ave = sum(val) / size(val)
    res = 0.d0
    do i = 1, size(val)
      res = res + (val(i) - ave)**2
    end do
    res = dsqrt(res/size(val))

  end function


end module
