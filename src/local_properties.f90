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

   use kinds

#ifdef _GPU
   use F_B_C
   use iso_c_binding
#endif
contains

#ifdef _GPU
   subroutine gpu_local_property_predict(n_sparse, soap, soap_d, &
        & Qs_d, alphas_d, e0, delta, zeta0, local_properties, local_properties_d, &
        & do_derivatives, soap_der_d, local_properties_cart_der,&
        & local_properties_cart_der_d, n_pairs, l_index_d,&
        & cublas_handle, gpu_stream)

      implicit none

      integer(c_int), intent(in) :: n_sparse
      real(c_double), intent(in), target :: soap(:, :)
      real(c_double), intent(in), target :: delta
      real(c_double), intent(in), target :: e0
      real(c_double), intent(in), target :: zeta0
      type(c_ptr), intent(in) :: soap_d
      type(c_ptr), intent(in) :: soap_der_d
      logical, intent(in) :: do_derivatives
      type(c_ptr), intent(inout) :: cublas_handle
      type(c_ptr), intent(inout) :: gpu_stream
      type(c_ptr), intent(inout) :: alphas_d
      type(c_ptr), intent(inout) :: Qs_d
      real(c_double), intent(out), target:: local_properties(:)
      real(c_double), intent(out), target:: local_properties_cart_der(:, :)

      real(c_double) :: zeta
      real(c_double) :: cdelta_ene
      real(c_double) :: mzetam
      real(c_double) :: cdelta_force
      logical :: is_zeta_int = .false.
      integer(c_int) :: n_sites
      integer(c_int) :: n_soap
      integer(c_int) :: i
      integer(c_int) :: j
      integer(c_int) :: k
      integer(c_int) :: l
      integer(c_int) :: j2
      integer(c_int) :: zeta_int
      integer(c_int) :: n_sites0
      integer(c_int) :: k1
      integer(c_int) :: k2

      integer(c_int), intent(in) :: n_pairs
      real(c_double), allocatable, target :: kernels(:, :)
      real(c_double), allocatable, target :: Qss(:, :)
      real(c_double), allocatable, target :: Qs_copy(:, :)
      real(c_double), allocatable, target :: this_Qss(:)
      real(c_double), allocatable, target :: kernels_copy(:, :)
      real(c_double), allocatable, target :: this_force_h(:, :)
      integer(c_size_t) :: st_alphas
      integer(c_size_t) :: st_Qs
      integer(c_size_t) :: st_kernels
      integer(c_size_t) :: st_local_properties
      integer(c_size_t) :: st_soap
      integer(c_size_t) :: st_local_properties_cart_der

      integer(c_int) :: size_kernels
      integer(c_int) :: size_soap
      integer(c_int) :: size_Qs
      integer(c_int) :: size_alphas
      integer(c_int) :: size_local_properties
      integer(c_int) :: size_local_properties_cart_der
      integer(c_int) :: maxnn
      type(c_ptr) :: kernels_copy_d
      type(c_ptr) :: kernels_d
      type(c_ptr), intent(inout) :: local_properties_d
      type(c_ptr), intent(inout) :: local_properties_cart_der_d
      type(c_ptr) :: kernels_der_d
      type(c_ptr) :: Qss_d
      type(c_ptr) :: Qs_copy_d !
      type(c_ptr) :: this_Qss_d
      type(c_ptr), intent(inout) :: l_index_d
      integer :: n1local_properties_cart_der
      integer :: n2local_properties_cart_der

      cdelta_ene = delta*delta
      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
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

      allocate (kernels(1:n_sites, 1:n_sparse))
      kernels = 0.d0
      allocate (kernels_copy(1:n_sites, 1:n_sparse))

      size_kernels = n_sites*n_sparse
      size_soap = n_soap*n_sites
      size_Qs = n_soap*n_sparse
      ! size_alphas=n_sparse

      size_local_properties = n_sites

      st_kernels = size_kernels*(sizeof(kernels(1, 1)))
      st_Qs = size_Qs*(sizeof(e0))
      ! st_alphas=size_alphas*(sizeof(alphas(1)))
      st_local_properties = size_local_properties*(sizeof(local_properties(1)))

      call gpu_malloc_async(kernels_d, st_kernels, gpu_stream)
      call gpu_malloc_async(kernels_copy_d, st_kernels, gpu_stream)

      call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
      call gpu_kernels_pow(kernels_d, kernels_copy_d, zeta, size_kernels, gpu_stream)
      call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, local_properties_d, n_sites, n_sparse)

      call gpu_axpe(local_properties_d, cdelta_ene, e0, size_local_properties, gpu_stream)

      call cpy_dtoh(local_properties_d, c_loc(local_properties), st_local_properties, gpu_stream)

      ! Now we do the derivatives
      if (do_derivatives) then

         call gpu_malloc_async(kernels_der_d, st_kernels, gpu_stream)
         st_soap = size_soap*sizeof(local_properties(1))
         call gpu_malloc_async(Qss_d, st_soap, gpu_stream)
         call gpu_malloc_async(Qs_copy_d, st_Qs, gpu_stream)
         call cpy_dtod(Qs_d, Qs_copy_d, st_Qs, gpu_stream)

         mzetam = zeta - 1
         call gpu_kernels_pow(kernels_d, kernels_der_d, mzetam, size_kernels, gpu_stream)

         if (n_sites < n_soap) then
            call gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse, gpu_stream)
         else
            call gpu_matvect(Qs_copy_d, alphas_d, n_soap, n_sparse, gpu_stream)
         end if

         cdelta_force = -zeta*delta**2
         call gpu_blas_mmul_n_t(cublas_handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                                n_soap, n_sites, cdelta_force)

         local_properties_cart_der = 0.d0

         n1local_properties_cart_der = size(local_properties_cart_der, 1)
         n2local_properties_cart_der = size(local_properties_cart_der, 2)
         size_local_properties_cart_der = n1local_properties_cart_der*n2local_properties_cart_der

         st_local_properties_cart_der = size_local_properties_cart_der*sizeof(local_properties_cart_der(1, 1))

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
      call gpu_free_async(kernels_copy_d, gpu_stream)

      deallocate (kernels, kernels_copy)

   end subroutine gpu_local_property_predict

#endif
   ! subroutine get_local_property_details( n_soap_turbo, soap_turbo_hypers, n_local_properties_tot, local_property_labels, write_local_properties )
   !   implicit none
   !   integer, intent(in) :: n_soap_turbo
   !   integer, intent(out) :: n_local_properties_tot
   !   character*1024, allocatable, intent(inout) ::  local_property_labels(:)
   !   character*1024, allocatable :: local_property_labels_temp(:), local_property_labels_temp2(:)
   !   integer :: i, j, k, i2

   !  i2 = 1 ! using this as a counter for the labels
   !  do j = 1, n_soap_turbo
   !     if( soap_turbo_hypers(j)%has_local_properties )then
   !        ! This property has the labels of the quantities to
   !        ! compute. We must specify the number of local properties, for the sake of coding simplicity

   !        if(.not. allocated(local_property_labels))then
   !           allocate(local_property_labels(1:n_local_properties_tot))
   !           do i = 1, n_local_properties_tot
   !              local_property_labels(i) = soap_turbo_hypers(j)%local_property_models(i)%label
   !              write(*,*)' Local property found                  |'
   !              write(*,'(A,1X,I8,1X,A,1X,A)')' Descriptor ', j,&
   !                   & trim(soap_turbo_hypers(j)&
   !                   &%local_property_models(i)%label),  ' |'
   !           end do
   !        else
   !           allocate( local_property_labels_temp( 1:n_local_properties_tot - soap_turbo_hypers(j)%n_local_properties ))
   !           local_property_labels_temp = local_property_labels
   !           deallocate(local_property_labels)
   !           allocate(local_property_labels(1:n_local_properties_tot))

   !           do i = 1, nprop
   !              local_property_labels(i + n_local_properties_tot -&
   !                   & nprop) = soap_turbo_hypers(j)&
   !                   &%local_property_models(i)%label
   !              write(*,*)' Local property found                  |'
   !              write(*,'(A,1X,I8,1X,A,1X,A)')' Descriptor ', j,&
   !                   & trim(soap_turbo_hypers(j)&
   !                   &%local_property_models(i)%label),  ' |'

   !  ! Now we create an irreducible list of the labels
   !  if (n_local_properties_tot > 0)then
   !     allocate( local_property_labels_temp( 1:1 ))
   !     local_property_labels_temp(1) = local_property_labels(1)
   !     i2 = 1
   !     if (n_local_properties_tot > 1)then
   !        do i = 2, n_local_properties_tot
   !           label_in_list = .false.
   !           ! Iterate through irreducible list to see if there is a mismatch
   !           do j = 1, size( local_property_labels_temp, 1 )
   !              if (trim( local_property_labels_temp(j) ) == trim( local_property_labels(i) )) label_in_list = .true.
   !           end do
   !           if (.not. label_in_list) then
   !              i2 = i2 + 1
   !              allocate(local_property_labels_temp2(1:i2))
   !              local_property_labels_temp2(1:i2-1) = local_property_labels_temp(1:i2-1)
   !              local_property_labels_temp2(i2)     = local_property_labels(i)
   !              deallocate(local_property_labels_temp)
   !              allocate(local_property_labels_temp(1:i2))
   !              local_property_labels_temp(1:i2) = local_property_labels_temp2(1:i2)
   !              deallocate(local_property_labels_temp2)
   !           end if
   !        end do
   !     end if

   !     params%n_local_properties = i2

   !     ! Now we can have an array which has a soap turbo index as an input and it can give us the corresponding label
   !     allocate(local_property_indexes(1:n_local_properties_tot))
   !     i2 = 1
   !     do i = 1, params%n_local_properties
   !        do j = 1, n_local_properties_tot
   !           if ( trim(local_property_labels(j)) == trim( local_property_labels_temp(i) ) )then

   !              local_property_indexes(j) = i

   !                 ! Check if there is experimental data for one to do xps fitting
   !                 do i2 = 1, params%n_exp
   !                    if(( trim(params%exp_data(i2)%label) == "xps" .and.  &
   !                         .not. ( trim(params%exp_data(i2)%file_data) == "none" )))then
   !                       valid_xps = .true.
   !                       do k = 1, soap_turbo_hypers(j)%n_local_properties
   !                          if (trim(soap_turbo_hypers(j)%local_property_models(k)%label) == "xps")then
   !                             soap_turbo_hypers(j)%local_property_models(k)%do_derivatives = .false.
   !                             if( params%exp_forces .and. params%do_derivatives)then
   !                                soap_turbo_hypers(j)%local_property_models(k)%do_derivatives = .true.
   !                             end if
   !                          end if
   !                       end do
   !                    end if
   !                 end do
   !              end if

   !     print *, "n_local_properties ", params%n_local_properties
   !     print *, "n_local_properties_tot ", n_local_properties_tot
   !     ! print *, "local_property_labels ", local_property_labels
   !     ! print *, "local_property_labels_temp (irreducible) ", local_property_labels_temp

   subroutine local_property_predict(soap, Qs, alphas, V0, delta, zeta, V, &
                                     do_derivatives, soap_cart_der, n_neigh, V_der, &
                                     zero_trunc, label)

      implicit none

      real(dp), intent(in) :: soap(:, :)
      real(dp), intent(in) :: Qs(:, :)
      real(dp), intent(in) :: alphas(:)
      real(dp), intent(in) :: V0
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: zeta
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      integer, intent(in) :: n_neigh(:)
      logical, intent(in) :: do_derivatives
      real(dp), intent(out) :: V(:)
!     Clamp negative predictions at zero. Right for a Hirshfeld volume, where a
!     negative value is meaningless; wrong for an atomic charge, where it is
!     half the atoms. Absent means clamp, which is what the standalone caller
!     below has always done.
      logical, intent(in), optional :: zero_trunc
!     Named only for the warning below, so that "a negative value was floored"
!     says which property it was.
      character(len=*), intent(in), optional :: label
      logical :: truncate
      integer :: n_floored
      logical, save :: warned = .false.
      real(dp), intent(out) :: V_der(:, :)
      real(dp), allocatable :: K(:, :)
      real(dp), allocatable :: K_der(:, :)
      real(dp), allocatable :: Qss(:, :)
      real(dp), allocatable :: Qs_copy(:, :)
      integer :: n_sites
      integer :: n_soap
      integer :: n_sparse
      integer :: zeta_int
      integer :: n_pairs
      integer :: i
      integer :: j
      integer :: i2
      integer :: cart

      n_sparse = size(alphas)
      n_soap = size(soap, 1)
      n_sites = size(soap, 2)

      allocate (K(1:n_sites, 1:n_sparse))
      if (do_derivatives) then
         n_pairs = size(soap_cart_der, 3)
         allocate (K_der(1:n_sites, 1:n_sparse))
         allocate (Qss(1:n_sites, 1:n_soap))
         allocate (Qs_copy(1:n_soap, 1:n_sparse))
      end if

      if (n_sites > 0) then
         call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                    K, n_sites)
      end if

      zeta_int = nint(zeta)
      if (dabs(dfloat(zeta_int) - zeta) < 1.d-10) then
         if (do_derivatives) then
            K_der = zeta_int*K**(zeta_int - 1)
            if (n_sites < n_soap) then
               do i = 1, n_sites
                  K_der(i, :) = K_der(i, :)*alphas(:)
                  Qs_copy = Qs
               end do
            else
               do i = 1, n_soap
                  Qs_copy(i, :) = Qs(i, :)*alphas(:)
               end do
            end if
         end if
         K = K**zeta_int
      else
         if (do_derivatives) then
            K_der = zeta*K**(zeta - 1.d0)
            if (n_sites < n_soap) then
               do i = 1, n_sites
                  K_der(i, :) = K_der(i, :)*alphas(:)
                  Qs_copy = Qs
               end do
            else
               do i = 1, n_soap
                  Qs_copy(i, :) = Qs(i, :)*alphas(:)
               end do
            end if
         end if
         K = K**zeta
      end if

!    V = delta**2 * matmul( K, alphas ) + V0
      if (n_sites > 0) then
         call dgemm("n", "n", n_sites, 1, n_sparse, delta**2, K, n_sites, alphas, n_sparse, 0.d0, V, n_sites)
      end if
      V = V + V0

      truncate = .true.
      if (present(zero_trunc)) truncate = zero_trunc

      if (truncate) then
         n_floored = 0
         do i = 1, size(V)
            if (V(i) < 0.d0) then
               V(i) = 0.d0
               n_floored = n_floored + 1
            end if
         end do
!        A Hirshfeld volume or a binding energy cannot be negative, so this
!        floor should never do anything. When it does, the prediction is wrong
!        and the floored value is not a repair -- it is the wrong answer,
!        rounded up. Say so once rather than let it pass.
         if (n_floored > 0 .and. .not. warned) then
            warned = .true.
            if (present(label)) then
               write (*, *) "WARNING: ", n_floored, " negative values of "//trim(label)// &
                  " were floored at zero."
            else
               write (*, *) "WARNING: ", n_floored, " negative local property values were floored at zero."
            end if
            write (*, *) "         That quantity cannot be negative, so the model is"
            write (*, *) "         predicting something impossible; the floored values are"
            write (*, *) "         not a repair. Reported once per run."
         end if
      end if

      if (do_derivatives) then
         if (n_sites > 0) then
            call dgemm("n", "t", n_sites, n_soap, n_sparse, delta**2, K_der, n_sites, &
                       Qs_copy, n_soap, 0.d0, Qss, n_sites)
         end if
         j = 1
         do i = 1, n_sites
            do i2 = 1, n_neigh(i)
!              A clamped site has no gradient, because its value no longer
!              depends on the descriptor. Without the clamp a zero is an
!              ordinary value and its gradient is the ordinary one.
               if (truncate .and. V(i) == 0.d0) then
                  V_der(1:3, j) = 0.d0
               else
                  do cart = 1, 3
                     V_der(cart, j) = dot_product(Qss(i, :), soap_cart_der(cart, :, j))
                  end do
               end if
               j = j + 1
            end do
         end do
         deallocate (Qs_copy)
         deallocate (Qss)
         deallocate (K_der)
      end if
      deallocate (K)

   end subroutine

end module
