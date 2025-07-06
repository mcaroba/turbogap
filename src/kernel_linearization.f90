! Note that this module only makes sense for zeta integer

module kernel_linearization

  contains

!**************************************************************************
!
! This function checks whether using a linearized kernel is faster than the
! usual kernel.
!
  function check_kernel_vs_linear(zeta, n_sites, n_sparse, n_soap, n_linear, &
                                  do_timing) result(do_linear)

    implicit none

    integer, intent(in) :: zeta, n_sparse, n_soap, n_sites, n_linear
    real*8 :: alphas(1:n_sparse), alphas_linear(1:n_linear), E(1:n_sites)
    real*8 :: Qs(1:n_soap, 1:n_sparse), desc(1:n_soap, 1:n_sites), &
              kernels(1:n_sites, 1:n_sparse), kernels_copy(1:n_sites, 1:n_sparse), &
              desc_copy(1:n_sites, 1:n_soap), desc_linear(1:n_sites, 1:n_linear)
    real*8 :: t1, t2, t_kernel, t_linear
    logical :: do_linear, do_timing

    call random_number(alphas)
    call random_number(alphas_linear)
    call random_number(Qs)
    call random_number(desc)
    desc_linear = 0.d0
    kernels = 0.d0
    do_linear = .false.

    call cpu_time(t1)
    call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, desc, n_soap, &
                Qs, n_soap, 0.d0, kernels, n_sites)
    kernels_copy = kernels**zeta
    E = matmul(kernels_copy, alphas)
    call cpu_time(t2)
    t_kernel = t2 - t1

    call cpu_time(t1)
    desc_copy = transpose(desc)
    call get_linearized_desc(zeta, desc_copy, desc_linear)
    E = matmul(desc_linear, alphas_linear)
    call cpu_time(t2)
    t_linear = t2 - t1

    if( do_timing )then
      write(*,*) "Time energies (usual kernel): ", t_kernel
      write(*,*) "Time energies (linearized kernel): ", t_linear
    end if

!   Add forces check, since it might change the decision chain
    if( t_linear < t_kernel )then
      do_linear = .true.
    else
!     Check here if forces would be faster
      do_linear = check_forces_kernel_vs_linear(zeta, n_sites, n_sparse, n_soap, &
                                                n_linear, do_timing)
    end if

  end function
!**************************************************************************

!**************************************************************************
!
! This function checks whether using a linearized kernel is faster than the
! usual kernel when calculating forces.
!
  function check_forces_kernel_vs_linear(zeta, n_sites, n_sparse, n_soap, &
                                         n_linear, do_timing) result(do_linear)

    implicit none

    integer, intent(in) :: zeta, n_sparse, n_soap, n_sites, n_linear
    integer :: n_neigh(1:n_sites)
    real*8 :: alphas(1:n_sparse), this_Qss(1:n_soap), this_force(1:3), n_neigh0(1:n_sites)
    real*8 :: Qs(1:n_soap, 1:n_sparse), Qs_copy(1:n_soap, 1:n_sparse), Qss(1:n_sites, 1:n_soap), &
              kernels(1:n_sites, 1:n_sparse), kernels_der(1:n_sites, 1:n_sparse), &
              desc(1:n_sites, 1:n_soap), desc_linear(1:n_sites, 1:n_linear), &
              Ms_linear(1:n_soap,1:n_linear), kernels_linear(1:n_soap, 1:n_sites)
    real*8, allocatable :: desc_der(:,:,:), F(:,:)
    integer :: nn = 100, n_atom_pairs, l, i, j, k
    real*8 :: t1, t2, t_kernel, t_linear
    logical :: do_linear, do_timing

    do_linear = .false.

!   Initialize arrays for kernel calculation of forces
    call random_number(alphas)
    call random_number(Qs)
    call random_number(kernels)
    Qss = 0.d0
    Qs_copy = Qs

    call random_number(n_neigh0)
    n_atom_pairs = 0
    do i = 1, n_sites
      n_neigh(i) = 1 + floor(nn*n_neigh0(i))
      n_atom_pairs = n_atom_pairs + n_neigh(i)
    end do
    allocate( desc_der(1:3, 1:n_soap, 1:n_atom_pairs) )
    allocate( F(1:3, 1:n_atom_pairs) )
    call random_number(desc_der)
    F = 0.d0

!   Performance for the usual kernel implementation of forces
    t_kernel = 0.d0
    call cpu_time(t1)
    kernels_der = kernels**(zeta - 1)
    if( n_sites < n_soap )then
      do i = 1, n_sites
        kernels_der(i,:) = kernels_der(i,:)*alphas(:)
      end do
    else
      do i = 1, n_soap
        Qs_copy(i,:) = Qs(i,:)*alphas(:)
      end do
    end if
    if( n_sites > 0 )then
      call dgemm("n", "t", n_sites, n_soap, n_sparse, 1.d0, kernels_der, n_sites, &
                Qs, n_soap, 0.d0, Qss, n_sites)
    end if

    if( do_timing )then
      call cpu_time(t2)
      t_kernel = t2 - t1 
      write(*,*) "Time forces (usual kernel part 1): ", t_kernel
      call cpu_time(t1)
    end if

    l = 0
    do i = 1, n_sites
      this_Qss = Qss(i, 1:n_soap)
      do j = 1, n_neigh(i)
        l = l + 1
        do k = 1, 3
          this_force(k) = dot_product(this_Qss, desc_der(k, 1:n_soap, l))
          F(k, l) = F(k, l) + this_force(k)
        end do
      end do
    end do
    call cpu_time(t2)
    t_kernel = t_kernel + (t2 - t1)

    if( do_timing )then
      write(*,*) "Time forces (usual kernel part 2): ", t2 - t1
      write(*,*) "Time forces (usual kernel): ", t_kernel
    end if

!   Initialize arrays for linearized calculation of forces
    call random_number(Ms_linear)
    call random_number(desc)
    desc_linear = 0.d0

!   Time the performance of linearized implementation of forces
    call cpu_time(t1)
    call get_linearized_desc(zeta, desc, desc_linear)
    call dgemm('n', 't', n_soap, n_sites, n_linear, 1.d0, Ms_linear, n_soap, &
               desc_linear, n_sites, 0.d0, kernels_linear, n_soap)
    l = 1
    do i = 1, n_sites
      this_Qss = kernels_linear(1:n_soap, i)
      do j = 1, n_neigh(i)
        do k = 1, 3
          this_force(k) = dot_product(this_Qss, desc_der(k, 1:n_soap, l))
          F(k, l) = F(k, l) + this_force(k)
        end do
        l = l + 1
      end do
    end do
    call cpu_time(t2)
    t_linear = t2 - t1

    if( do_timing )then
      write(*,*) "Time forces (linearized kernel): ", t_linear
    end if

    if( t_linear < t_kernel )then
      do_linear = .true.
    end if

    deallocate( F, desc_der )

  end function
!**************************************************************************


!**************************************************************************
!
! This subroutine gives the linearized version of a vector v**zeta,
! for zeta integer. That is, the original vector v = (q1, ..., qN) 
! is mapped to a linear model based on its multinomial expansion, 
! returning
!               v_lin = (q1_lin, ..., qN'_lin) 
! where
!             N' = (zeta + N - 1) ··· N / (zeta)!
! and
!     qi_lin = (zeta)! · prod_{k=1}^N qk^alpha_k / (alpha_k)!
! for i --> (alpha_1, ..., alpha_k) such that sum_k alpha_k = zeta
!
! It is possible to pass a matrix of vectors V, to get V_lin, with
! V**T = (v1 v2 ... ) and V_lin**T = (v1_lin v2_lin ... )
!
! NOTE: current implementation only works for zeta == 2, 3 !!!!
!       also, it does NOT include the multinomial coeff.
!
  subroutine get_linearized_desc(zeta, V, V_lin)

    implicit none

    integer, intent(in) :: zeta
    real*8, intent(in) :: V(:,:)
    integer :: n_soap, m, i, j, k, l
    real*8, intent(out) :: V_lin(:,:)

    m = size(V, 1) ! this dim. does not change
    n_soap = size(V, 2)

    if( zeta == 1 )then
      V_lin = V
    else if( zeta == 2 )then
      V_lin(1:m, 1:n_soap) = V**zeta
      k = n_soap + 1
      do i = 1, n_soap
        do j = i + 1, n_soap
          V_lin(1:m, k) = V(1:m, i) * V(1:m, j)
          k = k + 1
        end do
      end do
    else if( zeta == 3 )then
      V_lin(1:m, 1:n_soap) = V**zeta
      k = n_soap + 1
      do i = 1, n_soap
        do j = i + 1, n_soap
          V_lin(1:m, k) = V(1:m, i)**2 * V(1:m, j)
          V_lin(1:m, k + 1) = V(1:m, i) * V(1:m, j)**2
          k = k + 2
        end do
      end do
      do i = 1, n_soap
        do j = i + 1, n_soap
          do l = j + 1, n_soap
            V_lin(1:m, k) = V(1:m, i) * V(1:m, j) * V(1:m, l)
            k = k + 1
          end do
        end do
      end do
    end if

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine gives the multinomial coefficients of a SOAP kernel
! It assumes the that M has size (1:n_sparse, 1:n_linear)
!
  subroutine get_multinomial_coeffs(zeta, n_soap, M)

    implicit none

    integer, intent(in) :: zeta, n_soap
    real*8, intent(out) :: M(:,:)
    integer :: n_sparse, n_linear

    n_sparse = size(M, 1)
    n_linear = size(M, 2)

    M = 1.d0
    if( zeta == 2 )then
      M(1:n_sparse, n_soap + 1:n_linear) = 2.d0
    else if( zeta == 3 )then
      M(1:n_sparse, n_soap + 1:n_soap**2) = 3.d0
      M(1:n_sparse, n_soap**2 + 1:n_linear) = 6.d0
    end if

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This function computes the linearized sparse contribution for energies,
! given by
!                 alphas_linear = Qs_linear · alphas
! For forces, use get_linearized_sparse_forces() 
!
  function get_linearized_sparse(zeta, alphas, Qs) result(alphas_linear)

    implicit none

    integer, intent(in) :: zeta
    real*8, intent(in) :: alphas(:), Qs(:,:)
    integer :: n_soap, n_sparse, n_linear, i
    real*8, allocatable :: Qs_linear(:,:), multinomial_coeffs(:,:)
    real*8, allocatable :: alphas_linear(:)

    n_soap = size(Qs, 2) ! here we assume Qs(1:n_sparse, 1:n_soap), I could check it using alphas size
    n_sparse = size(alphas)
!   Size of the linearized model
    n_linear = 1
    do i = 1, zeta
      n_linear = n_linear*(n_soap + i - 1)/i
    end do

!   Allocate arrays
    allocate( alphas_linear(1:n_linear) )
    allocate( Qs_linear(1:n_sparse, 1:n_linear) )
    allocate( multinomial_coeffs(1:n_sparse, 1:n_linear) )
    alphas_linear = 0.d0
    Qs_linear = 0.d0

!   Get alphas_linear
    call get_linearized_desc(zeta, Qs, Qs_linear)
!   add multiplicity using multinomial coefficients
    call get_multinomial_coeffs(zeta, n_soap, multinomial_coeffs)
    alphas_linear = matmul(alphas, multinomial_coeffs * Qs_linear)

    deallocate( Qs_linear, multinomial_coeffs )

  end function
!**************************************************************************

!**************************************************************************
!
! This function computes the linearized sparse contribution for forces, 
! given in this case by the matrix
!                 Ms_linear = (Qs * diag(alphas)) · (Qs_linear)**T
! Note that Qs_linear is computed here using zeta' = zeta - 1
!
  function get_linearized_sparse_forces(zeta, alphas, Qs) result(Ms_linear)

    implicit none

    integer, intent(in) :: zeta
    real*8, intent(in) :: alphas(:), Qs(:,:)
    integer :: zeta0, n_soap, n_sparse, n_linear, i
    real*8, allocatable :: Qss(:,:), Qs_linear(:,:), multinomial_coeffs(:,:)
    real*8, allocatable :: Ms_linear(:,:)

    n_soap = size(Qs, 2) ! here we assume Qs(1:n_sparse, 1:n_soap)
    n_sparse = size(alphas)
!   Size of the linearized model
    zeta0 = zeta - 1
    n_linear = 1
    do i = 1, zeta0
      n_linear = n_linear*(n_soap + i - 1)/i
    end do

!   Initialize arrays
    allocate( Ms_linear(1:n_soap, 1:n_linear) )
    allocate( Qs_linear(1:n_sparse, 1:n_linear) )
    allocate( Qss(1:n_sparse, 1:n_soap) )
    allocate( multinomial_coeffs(1:n_sparse, 1:n_linear) )
    Ms_linear = 0.d0
    Qs_linear = 0.d0
   
!   Get the linearized contribution
    call get_linearized_desc(zeta0, Qs, Qs_linear)
    call get_multinomial_coeffs(zeta0, n_soap, multinomial_coeffs)
    Qs_linear = Qs_linear * multinomial_coeffs
!   hadamard product with a vector 
    do i=1, n_soap
      Qss(1:n_sparse, i) = Qs(1:n_sparse, i) * alphas
    end do
!   (1:n_soap, 1:n_linear) = (1:n_sparse, 1:n_soap)^T * (1:n_sparse, 1:n_linear)
    call dgemm('t','n', n_soap, n_linear, n_sparse, 1.d0, Qss, n_sparse, Qs_linear, &
                n_sparse, 0.d0, Ms_linear, n_soap)
    deallocate( Qs_linear, Qss )

  end function
!**************************************************************************


end module
