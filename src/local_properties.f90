! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, local_properties.f90, is copyright (c) 2019-2023, Miguel A.
! HND X   Caro, Heikki Muhli and Tigany Zarrouk
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

module local_prop

  use F_B_C
  use iso_c_binding

  contains

    subroutine gpu_local_property_predict( n_sparse, soap, soap_d, &
         & Qs_d, alphas_d, e0, delta, zeta0, local_properties, local_properties_d, &
         & do_derivatives, soap_der_d, local_properties_cart_der,&
         & local_properties_cart_der_d, n_pairs, l_index_d,&
         & cublas_handle, gpu_stream  )

    implicit none

    !   Input variables
    integer(c_int), intent(in) :: n_sparse
    real(c_double), intent(in), target :: soap(:,:), delta, e0, zeta0
    type(c_ptr), intent(in) :: soap_d, soap_der_d
    logical, intent(in) :: do_derivatives 
    type(c_ptr), intent(inout) :: cublas_handle, gpu_stream, alphas_d, Qs_d
    real(c_double), intent(out), target:: local_properties(:), local_properties_cart_der(:,:)
    
    real(c_double) ::  zeta, cdelta_ene, mzetam, cdelta_force
    logical :: is_zeta_int = .false.
    integer(c_int) :: n_sites, n_soap, i, j, k, l, j2, zeta_int, n_sites0, k1, k2

    integer(c_int), intent(in) :: n_pairs 
    real(c_double), allocatable, target :: kernels(:,:), &
                            Qss(:,:), Qs_copy(:,:), this_Qss(:), &
                           kernels_copy(:,:), this_force_h(:,:)
    integer(c_size_t) :: st_alphas, st_Qs, st_kernels, st_local_properties, st_soap, st_local_properties_cart_der

    integer(c_int) :: size_kernels, size_soap, size_Qs, size_alphas, size_local_properties, size_local_properties_cart_der, maxnn
    type(c_ptr) :: kernels_copy_d, kernels_d
    type(c_ptr), intent(inout) :: local_properties_d, local_properties_cart_der_d
    type(c_ptr) :: kernels_der_d, Qss_d, Qs_copy_d !, this_Qss_d
    type(c_ptr), intent(inout) :: l_index_d
    integer :: n1local_properties_cart_der, n2local_properties_cart_der

    
    cdelta_ene=delta*delta
    if( dabs(zeta0-dfloat(int(zeta0))) < 1.d-5 )then
      is_zeta_int = .true.
      zeta_int = int(zeta0)
      zeta = dfloat(zeta_int)
    else
      zeta = zeta0
    end if

    ! n_sparse = size(alphas)
    n_soap = size(soap, 1)
    n_sites = size(soap, 2)
!    n_sites0 = size(forces, 2)

    allocate( kernels(1:n_sites, 1:n_sparse) )
    kernels = 0.d0
    allocate( kernels_copy(1:n_sites, 1:n_sparse) )

    size_kernels = n_sites * n_sparse
    size_soap    = n_soap  * n_sites
    size_Qs      = n_soap  * n_sparse
    ! size_alphas=n_sparse

    size_local_properties=n_sites

    st_kernels=size_kernels*(sizeof(kernels(1,1)))
    st_Qs=size_Qs*(sizeof(e0))
    ! st_alphas=size_alphas*(sizeof(alphas(1)))
    st_local_properties=size_local_properties*(sizeof(local_properties(1)))

    call gpu_malloc_all(kernels_d, st_kernels, gpu_stream)
    call gpu_malloc_all(kernels_copy_d, st_kernels, gpu_stream)
    !call gpu_malloc_all(local_properties_d, st_local_properties, gpu_stream)

    call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
    call gpu_kernels_pow( kernels_d, kernels_copy_d,  zeta, size_kernels, gpu_stream)
    call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, local_properties_d, n_sites, n_sparse)

    call gpu_axpe(local_properties_d, cdelta_ene, e0, size_local_properties, gpu_stream)

    call cpy_dtoh(local_properties_d, c_loc(local_properties), st_local_properties, gpu_stream)

    ! Now we do the derivatives
    if(do_derivatives) then

       
       
      call gpu_malloc_all(kernels_der_d, st_kernels, gpu_stream)
      st_soap = size_soap * sizeof(local_properties(1))
      call gpu_malloc_all(Qss_d, st_soap, gpu_stream)
      call gpu_malloc_all(Qs_copy_d, st_Qs, gpu_stream)
      call cpy_dtod(Qs_d, Qs_copy_d, st_Qs, gpu_stream) 
      
      mzetam=zeta-1
      call gpu_kernels_pow( kernels_d, kernels_der_d,  mzetam, size_kernels, gpu_stream)

      if(n_sites<n_soap) then
        call gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse, gpu_stream)
      else
        call gpu_matvect(Qs_copy_d, alphas_d, n_soap, n_sparse, gpu_stream)
      endif

      
      cdelta_force=-zeta*delta**2 
      call gpu_blas_mmul_n_t(cublas_handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
           n_soap, n_sites, cdelta_force)


      local_properties_cart_der = 0.d0

      n1local_properties_cart_der = size(local_properties_cart_der, 1)
      n2local_properties_cart_der = size(local_properties_cart_der, 2)
      size_local_properties_cart_der=n1local_properties_cart_der*n2local_properties_cart_der
      
      st_local_properties_cart_der = size_local_properties_cart_der*sizeof(local_properties_cart_der(1,1))
      !call gpu_malloc_all(local_properties_cart_der_d, st_local_properties_cart_der, gpu_stream)

      call gpu_local_property_derivatives(n_sites, &
                                              Qss_d, n_soap, l_index_d, &
                                              soap_der_d, &
                                              local_properties_cart_der_d, &
                                              n_pairs, gpu_stream)

      
      call cpy_dtoh(local_properties_cart_der_d, c_loc(local_properties_cart_der), st_local_properties_cart_der, gpu_stream)

      call gpu_free_async(kernels_der_d, gpu_stream)
      call gpu_free_async(Qss_d, gpu_stream)
      call gpu_free_async(Qs_copy_d, gpu_stream)
   end if



   call gpu_free_async(kernels_d, gpu_stream)
   call gpu_free(kernels_copy_d)


   deallocate( kernels, kernels_copy )

  end subroutine gpu_local_property_predict

    subroutine gpu_local_property_predict_single( n_sparse, soap, soap_d, &
         & Qs_d, alphas_d, e0, delta, zeta0, local_properties,&
         & cublas_handle, gpu_stream )

    implicit none

    !   Input variables
    integer(c_int), intent(in) :: n_sparse
    real(c_double), intent(in), target :: soap(:,:), delta, e0, zeta0
    type(c_ptr), intent(in) :: soap_d
    
    type(c_ptr), intent(inout) :: cublas_handle, gpu_stream, alphas_d, Qs_d
    real(c_double), intent(out), target:: local_properties(:)
    
    real(c_double) ::  zeta, cdelta_ene
    logical :: is_zeta_int = .false.
    integer(c_int) :: n_sites, n_soap, i, j, k, l, j2, zeta_int, n_sites0, k1, k2

    real(c_double), allocatable, target :: kernels(:,:), &
                            Qss(:,:), Qs_copy(:,:), this_Qss(:), &
                           kernels_copy(:,:), this_force_h(:,:)
    integer(c_size_t) :: st_alphas, st_Qs, st_kernels, st_local_properties, st_soap

    integer(c_int) :: size_kernels, size_soap, size_Qs, size_alphas, size_local_properties, maxnn
    type(c_ptr) :: kernels_copy_d, kernels_d, local_properties_d
    type(c_ptr) :: kernels_der_d, Qss_d, Qs_copy_d !, this_Qss_d
        
    cdelta_ene=delta*delta
    if( dabs(zeta0-dfloat(int(zeta0))) < 1.d-5 )then
      is_zeta_int = .true.
      zeta_int = int(zeta0)
      zeta = dfloat(zeta_int)
    else
      zeta = zeta0
    end if

    ! n_sparse = size(alphas)
    n_soap = size(soap, 1)
    n_sites = size(soap, 2)
!    n_sites0 = size(forces, 2)

    allocate( kernels(1:n_sites, 1:n_sparse) )
    kernels = 0.d0
    allocate( kernels_copy(1:n_sites, 1:n_sparse) )

    size_kernels = n_sites * n_sparse
    size_soap    = n_soap  * n_sites
    size_Qs      = n_soap  * n_sparse
    ! size_alphas=n_sparse

    size_local_properties=n_sites

    st_kernels=size_kernels*(sizeof(kernels(1,1)))
    st_Qs=size_Qs*(sizeof(e0))
    ! st_alphas=size_alphas*(sizeof(alphas(1)))
    st_local_properties=size_local_properties*(sizeof(local_properties(1)))

    call gpu_malloc_all(kernels_d, st_kernels, gpu_stream)
    call gpu_malloc_all(kernels_copy_d, st_kernels, gpu_stream)
    call gpu_malloc_all(local_properties_d, st_local_properties, gpu_stream)


    call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
    call gpu_kernels_pow( kernels_d, kernels_copy_d,  zeta, size_kernels, gpu_stream)
    call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, local_properties_d, n_sites, n_sparse)

    call gpu_axpe(local_properties_d, cdelta_ene, e0, size_local_properties, gpu_stream)

    call cpy_dtoh(local_properties_d, c_loc(local_properties), st_local_properties, gpu_stream)

    call gpu_free_async(kernels_d, gpu_stream)
    call gpu_free_async(kernels_copy_d, gpu_stream)
    call gpu_free(local_properties_d)

    deallocate( kernels, kernels_copy )

  end subroutine gpu_local_property_predict_single
!**************************************************************************

    

!     subroutine gpu_local_property_predict( n_sparse, soap, Qs_d, alphas_d, e0, delta, zeta0, local_properties, &
!          do_derivatives, soap_der, n_neigh, local_properties_der, soap_d, &
!          soap_der_d, cublas_handle, gpu_stream )

!     implicit none

!     !   Input variables
!     integer(c_int), intent(in) :: n_sparse
!     real(c_double), intent(in), target :: soap(:,:), soap_der(:,:,:), delta, e0, zeta0
!     logical, intent(in) :: do_derivatives
    
!     type(c_ptr), intent(inout) :: cublas_handle, gpu_stream, k2_i_site_d, alphas_d, Qs_d
!     real(c_double), intent(out):: local_properties(:)

!     real(c_double) ::  zeta, cdelta_ene, cdelta_force, mzetam
!     logical :: is_zeta_int = .false.
!     integer(c_int) :: n_sites, n_soap, i, j, k, l, j2, zeta_int, n_sites0, k1, k2

    
!     real(c_double), allocatable, target :: kernels(:,:), kernels_der(:,:), &
!                             Qss(:,:), Qs_copy(:,:), this_Qss(:), &
!                            kernels_copy(:,:), this_force_h(:,:)
!     integer(c_size_t) :: st_alphas, st_Qs, st_kernels, st_local_properties, st_soap

!     integer(c_int) :: size_kernels, size_soap, size_Qs, size_alphas, size_local_properties, maxnn
!     type(c_ptr) :: kernels_copy_d, kernels_d, local_properties_d
!     type(c_ptr) :: kernels_der_d, Qss_d, Qs_copy_d !, this_Qss_d


    
!     ! real*8, intent(in) :: soap(:,:), Qs(:,:), alphas(:), V0, delta, zeta, soap_cart_der(:,:,:)
! !     integer, intent(in) :: n_neigh(:)
! !     logical, intent(in) :: do_derivatives
! ! !   Output variables
! !     real*8, intent(out) :: V(:), V_der(:,:)
! ! !   Internal variables
! !     real*8, allocatable :: K(:,:), K_der(:,:), Qss(:,:), Qs_copy(:,:)

! !     integer :: i, j, i2, cart

        
!     cdelta_ene=delta*delta
!     if( dabs(zeta0-dfloat(int(zeta0))) < 1.d-5 )then
!       is_zeta_int = .true.
!       zeta_int = int(zeta0)
!       zeta = dfloat(zeta_int)
!     else
!       zeta = zeta0
!     end if

!     ! n_sparse = size(alphas)
!     n_soap = size(soap, 1)
!     n_sites = size(soap, 2)
! !    n_sites0 = size(forces, 2)

!     allocate( kernels(1:n_sites, 1:n_sparse) )
!     kernels = 0.d0
!     allocate( kernels_copy(1:n_sites, 1:n_sparse) )

!     size_kernels=n_sites*n_sparse
!     size_soap=n_soap*n_sites
!     size_Qs = n_soap*n_sparse
!     ! size_alphas=n_sparse

!     size_local_properties=n_sites

!     st_kernels=size_kernels*(sizeof(kernels(1,1)))
!     st_Qs=size_Qs*(sizeof(e0))
!     ! st_alphas=size_alphas*(sizeof(alphas(1)))
!     st_local_properties=size_local_properties*(sizeof(local_properties(1)))

!     call gpu_malloc_all(kernels_d,st_kernels, gpu_stream)
!     call gpu_malloc_all(kernels_copy_d,st_kernels, gpu_stream)

!     call gpu_malloc_all(local_properties_d,st_local_properties, gpu_stream)


!     call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
!     call gpu_kernels_pow( kernels_d, kernels_copy_d,  zeta, size_kernels, gpu_stream)
!     call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, local_properties_d, n_sites, n_sparse)

!     call gpu_axpe(local_properties_d, cdelta_ene, e0, size_local_properties, gpu_stream)

!     call cpy_dtoh(local_properties_d, c_loc(local_properties),st_local_properties, gpu_stream)

!     call gpu_free_async(kernels_d, gpu_stream)
!     call gpu_free_async(kernels_copy_d, gpu_stream)
!     call gpu_free_async(local_properties_d, gpu_stream)

    
!     ! call cpy_dtoh(energies_d,c_loc(tmp_energies),st_energies,gpu_stream)
!     ! energies=tmp_energies
    
    
! !     n_sparse = size(alphas)
! !     n_soap = size(soap, 1)
! !     n_sites = size(soap, 2)

! !     allocate( K(1:n_sites, 1:n_sparse) )
! !     if( do_derivatives )then
! !       n_pairs = size(soap_cart_der, 3)
! !       allocate( K_der(1:n_sites, 1:n_sparse) )
! !       allocate( Qss(1:n_sites, 1:n_soap) )
! !       allocate( Qs_copy(1:n_soap, 1:n_sparse) )
! !     end if


! !     if( n_sites > 0 )then
! !       call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
! !                 K, n_sites)
! !     end if

! !     zeta_int = nint(zeta)
! !     if( dabs( dfloat(zeta_int) - zeta ) < 1.d-10 )then
! !       if( do_derivatives )then
! !         K_der = zeta_int*K**(zeta_int-1)
! !         if( n_sites < n_soap)then
! !           do i = 1, n_sites
! !             K_der(i,:) = K_der(i,:) * alphas(:)
! !             Qs_copy = Qs
! !           end do
! !         else
! !           do i = 1, n_soap
! !             Qs_copy(i,:) = Qs(i,:) * alphas(:)
! !           end do
! !         end if
! !       end if
! !       K = K**zeta_int
! !     else
! !       if( do_derivatives )then
! !         K_der = zeta*K**(zeta-1.d0)
! !         if( n_sites < n_soap)then
! !           do i = 1, n_sites
! !             K_der(i,:) = K_der(i,:) * alphas(:)
! !             Qs_copy = Qs
! !           end do
! !         else
! !           do i = 1, n_soap
! !             Qs_copy(i,:) = Qs(i,:) * alphas(:)
! !           end do
! !         end if
! !       end if
! !       K = K**zeta
! !     end if

! ! !    V = delta**2 * matmul( K, alphas ) + V0
! !     if( n_sites > 0 )then
! !       call dgemm( "n", "n", n_sites, 1, n_sparse, delta**2, K, n_sites, alphas, n_sparse, 0.d0, V, n_sites)
! !     end if
! !     V = V + V0

! ! !   Make sure all V are >= 0
! !     do i = 1, size(V)
! !       if( V(i) < 0.d0 )then
! !         V(i) = 0.d0
! !       end if
! !     end do

! !     if( do_derivatives)then
! !       if( n_sites > 0 )then
! !         call dgemm("n", "t", n_sites, n_soap, n_sparse, delta**2, K_der, n_sites, &
! !                    Qs_copy, n_soap, 0.d0, Qss, n_sites)
! !       end if
! !       j = 1
! !       do i = 1, n_sites
! !         do i2 = 1, n_neigh(i)
! !           if( V(i) == 0.d0 )then
! !             V_der(1:3, j) = 0.d0
! !           else
! !              do cart = 1, 3
! !                 V_der(cart, j) = dot_product( Qss(i,:), soap_cart_der(cart, :, j) )
! !             end do
! !           end if
! !           j = j + 1
! !         end do
! !       end do
! !       deallocate(Qs_copy)
! !       deallocate(Qss)
! !       deallocate(K_der)
! !     end if
! !   deallocate(K)

! end subroutine gpu_local_property_predict
! !**************************************************************************

  

  subroutine local_property_predict( soap, Qs, alphas, V0, delta, zeta, V, &
                                do_derivatives, soap_cart_der, n_neigh, V_der )

    implicit none

!   Input variables
    real*8, intent(in) :: soap(:,:), Qs(:,:), alphas(:), V0, delta, zeta, soap_cart_der(:,:,:)
    integer, intent(in) :: n_neigh(:)
    logical, intent(in) :: do_derivatives
!   Output variables
    real*8, intent(out) :: V(:), V_der(:,:)
!   Internal variables
    real*8, allocatable :: K(:,:), K_der(:,:), Qss(:,:), Qs_copy(:,:)
    integer :: n_sites, n_soap, n_sparse, zeta_int, n_pairs
    integer :: i, j, i2, cart

    n_sparse = size(alphas)
    n_soap = size(soap, 1)
    n_sites = size(soap, 2)

    allocate( K(1:n_sites, 1:n_sparse) )
    if( do_derivatives )then
      n_pairs = size(soap_cart_der, 3)
      allocate( K_der(1:n_sites, 1:n_sparse) )
      allocate( Qss(1:n_sites, 1:n_soap) )
      allocate( Qs_copy(1:n_soap, 1:n_sparse) )
    end if


    if( n_sites > 0 )then
      call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                K, n_sites)
    end if

    zeta_int = nint(zeta)
    if( dabs( dfloat(zeta_int) - zeta ) < 1.d-10 )then
      if( do_derivatives )then
        K_der = zeta_int*K**(zeta_int-1)
        if( n_sites < n_soap)then
          do i = 1, n_sites
            K_der(i,:) = K_der(i,:) * alphas(:)
            Qs_copy = Qs
          end do
        else
          do i = 1, n_soap
            Qs_copy(i,:) = Qs(i,:) * alphas(:)
          end do
        end if
      end if
      K = K**zeta_int
    else
      if( do_derivatives )then
        K_der = zeta*K**(zeta-1.d0)
        if( n_sites < n_soap)then
          do i = 1, n_sites
            K_der(i,:) = K_der(i,:) * alphas(:)
            Qs_copy = Qs
          end do
        else
          do i = 1, n_soap
            Qs_copy(i,:) = Qs(i,:) * alphas(:)
          end do
        end if
      end if
      K = K**zeta
    end if

!    V = delta**2 * matmul( K, alphas ) + V0
    if( n_sites > 0 )then
      call dgemm( "n", "n", n_sites, 1, n_sparse, delta**2, K, n_sites, alphas, n_sparse, 0.d0, V, n_sites)
    end if
    V = V + V0

!   Make sure all V are >= 0
    do i = 1, size(V)
      if( V(i) < 0.d0 )then
        V(i) = 0.d0
      end if
    end do

    if( do_derivatives)then
      if( n_sites > 0 )then
        call dgemm("n", "t", n_sites, n_soap, n_sparse, delta**2, K_der, n_sites, &
                   Qs_copy, n_soap, 0.d0, Qss, n_sites)
      end if
      j = 1
      do i = 1, n_sites
        do i2 = 1, n_neigh(i)
          if( V(i) == 0.d0 )then
            V_der(1:3, j) = 0.d0
          else
             do cart = 1, 3
                V_der(cart, j) = dot_product( Qss(i,:), soap_cart_der(cart, :, j) )
            end do
          end if
          j = j + 1
        end do
      end do
      deallocate(Qs_copy)
      deallocate(Qss)
      deallocate(K_der)
    end if
  deallocate(K)

  end subroutine
!**************************************************************************

end module
