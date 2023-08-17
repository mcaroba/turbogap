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


  contains


    ! !**************************************************************************
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

    !        n_local_properties_tot = n_local_properties_tot + soap_turbo_hypers(j)%n_local_properties

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

    !           nprop = soap_turbo_hypers(j)%n_local_properties
    !           do i = 1, n_local_properties_tot - nprop
    !              local_property_labels(i) = local_property_labels_temp(i)
    !           end do

    !           deallocate(local_property_labels_temp)

    !           do i = 1, nprop
    !              local_property_labels(i + n_local_properties_tot -&
    !                   & nprop) = soap_turbo_hypers(j)&
    !                   &%local_property_models(i)%label
    !              write(*,*)' Local property found                  |'
    !              write(*,'(A,1X,I8,1X,A,1X,A)')' Descriptor ', j,&
    !                   & trim(soap_turbo_hypers(j)&
    !                   &%local_property_models(i)%label),  ' |'

    !           end do
    !        end if
    !     end if
    !  end do

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


    !              if ( trim(local_property_labels(j)) == "hirshfeld_v" )then
    !                 vdw_lp_index = i
    !                 valid_vdw = .true.
    !              end if


    !              if ( trim(local_property_labels(j)) == "core_electron_be" )then
    !                 core_be_lp_index = i

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


    !           end if
    !        end do
    !     end do

    !     print *, "n_local_properties ", params%n_local_properties
    !     print *, "n_local_properties_tot ", n_local_properties_tot
    !     ! print *, "local_property_labels ", local_property_labels
    !     ! print *, "local_property_labels_temp (irreducible) ", local_property_labels_temp

    !     deallocate(local_property_labels)
    !     allocate(local_property_labels(1:size(local_property_labels_temp,1)))
    !     local_property_labels = local_property_labels_temp
    !     deallocate(local_property_labels_temp)

    !     do i = 1, params%n_local_properties
    !        write(*,*)' Irreducible local properties                   |'
    !        write(*,'(1X,A)') trim( local_property_labels(i) )
    !     end do

    !     allocate( params%write_local_properties(1:params%n_local_properties) )
    !     params%write_local_properties = .true.
    !  end if


    ! end subroutine get_local_property_details


  subroutine local_property_predict(soap, Qs, alphas, V0, delta, zeta, V, &
                                    do_derivatives, soap_cart_der, n_neigh, V_der, &
                                    do_zero_floor)

    implicit none

!   Input variables
    real*8, intent(in) :: soap(:,:), Qs(:,:), alphas(:), V0, delta, zeta, soap_cart_der(:,:,:)
    integer, intent(in) :: n_neigh(:)
    logical, intent(in) :: do_derivatives
    logical, intent(in), optional :: do_zero_floor
!   Output variables
    real*8, intent(out) :: V(:), V_der(:,:)
!   Internal variables
    real*8, allocatable :: K(:,:), K_der(:,:), Qss(:,:), Qs_copy(:,:)
    integer :: n_sites, n_soap, n_sparse, zeta_int, n_pairs
    integer :: i, j, i2, cart

    ! The default behaviour should be _not_ to truncate predictions to zero!
    if (.not. present(do_zero_floor)) then
        do_zero_floor = .false.
    end if

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

    !   Make sure all V are >= 0, if requested
    if (do_zero_floor)
      do i = 1, size(V)
        if( V(i) < 0.d0 )then
          V(i) = 0.d0
        end if
      end do
    end if

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
