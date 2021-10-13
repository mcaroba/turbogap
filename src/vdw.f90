! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, vdw.f90, is copyright (c) 2019-2021, Miguel A. Caro and Heikki
! HND X   Muhli
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

module vdw

  use misc

  contains

!**************************************************************************
  subroutine hirshfeld_predict( soap, Qs, alphas, V0, delta, zeta, V, &
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
    logical :: is_zeta_int = .false.

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

    call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
              K, n_sites)

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
    call dgemm( "n", "n", n_sites, 1, n_sparse, delta**2, K, n_sites, alphas, n_sparse, 0.d0, V, n_sites)
    V = V + V0

!   Make sure all V are >= 0
    do i = 1, size(V)
      if( V(i) < 0.d0 )then
        V(i) = 0.d0
      end if
    end do

    if( do_derivatives)then
      call dgemm("n", "t", n_sites, n_soap, n_sparse, delta**2, K_der, n_sites, &
                 Qs_copy, n_soap, 0.d0, Qss, n_sites)
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





!**************************************************************************
  subroutine get_ts_energy_and_forces( hirshfeld_v, hirshfeld_v_cart_der, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, hirshfeld_v_neigh, &
                                       sR, d, c6_ref, r0_ref, alpha0_ref, do_forces, &
                                       energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v(:), hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), hirshfeld_v_neigh(:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_forces
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), neighbor_c6_ij(:), r0_ii(:), r0_ij(:), &
                           exp_damp(:), f_damp(:), c6_ij_free(:), neighbor_alpha0(:), &
                           pref_force1(:), pref_force2(:), r6(:), r6_der(:)
    real*8 :: time1, time2, c6_ii, c6_jj, r0_i, r0_j, alpha0_i, alpha0_j, rbuf, this_force(1:3)
    integer, allocatable:: i_buffer(:)
    integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0
    integer :: i, j, i2, j2, k, n_in_buffer, k1, k2
    logical, allocatable :: is_in_buffer(:)
    logical :: do_timing = .false.

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!   We precompute the C6 coefficients of all the neighbors
    allocate( neighbor_c6_ii(1:n_pairs) )
    allocate( neighbor_c6_ij(1:n_pairs) )
    allocate( r0_ij(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( exp_damp(1:n_pairs) )
    allocate( f_damp(1:n_pairs) )
    allocate( c6_ij_free(1:n_pairs) )
    allocate( r6(1:n_pairs) )
    allocate( is_in_buffer(1:n_pairs) )
    is_in_buffer = .false.
    if( do_forces )then
      allocate( pref_force1(1:n_sites) )
      allocate( pref_force2(1:n_sites) )
      allocate( r6_der(1:n_pairs) )
    end if

    if( do_timing) then
      call cpu_time(time1)
    end if

!   Check which atoms are in the buffer region
    do k = 1, n_pairs
      j = neighbor_species(k)
      neighbor_c6_ii(k) = c6_ref(j)
      r0_ii(k) = r0_ref(j)
      neighbor_alpha0(k) = alpha0_ref(j)
    end do
    n_in_buffer = 0
    if( buffer > 0.d0 .or. buffer_inner > 0.d0 )then
      do k = 1, n_pairs
        if( (rjs(k) > rcut-buffer .and. rjs(k) < rcut) .or. &
            (rjs(k) < rcut_inner+buffer_inner .and. rjs(k) > rcut_inner) )then
          n_in_buffer = n_in_buffer + 1
          is_in_buffer(k) = .true.
        end if
      end do
    end if
    allocate( i_buffer(1:n_in_buffer) )
    if( buffer > 0.d0 .or. buffer_inner > 0.d0 )then
      i = 0
      do k = 1, n_pairs
        if( is_in_buffer(k) )then
          i = i + 1
          i_buffer(i) = k
        end if
      end do
    end if

!   Precompute r6 = 1/r^6 * cutoff(r) and its derivative
    r6 = 1.d0 / rjs**6
!   We do the forces first, so that we have not modified r6 yet
    if( do_forces )then
      r6_der = -6.d0 / rjs * r6
      do i = 1, n_in_buffer
        k = i_buffer(i)
        if( rjs(k) <= rcut_inner + buffer_inner )then
          rbuf = -(rjs(k)-rcut_inner-buffer_inner)/buffer_inner
          r6_der(k) = r6_der(k) + 6.d0 * r6(k) / (-buffer_inner) * (-rbuf + rbuf**2)
        else
          rbuf = (rjs(k)-rcut+buffer)/buffer
          r6_der(k) = r6_der(k) + 6.d0 * r6(k) / buffer * (-rbuf + rbuf**2)
        end if
      end do
    end if
!   Now we do these terms
    do i = 1, n_in_buffer
      k = i_buffer(i)
      if( rjs(k) <= rcut_inner + buffer_inner )then
        rbuf = (rjs(k)-rcut_inner-buffer_inner)/(-buffer_inner)
      else
        rbuf = (rjs(k)-rcut+buffer)/buffer
      end if
      r6(k) = r6(k) * ( 1.d0 - 3.d0*rbuf**2 + 2.d0*rbuf**3 )
    end do

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 1:", time2-time1
      call cpu_time(time1)
    end if

!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 2:", time2-time1
      call cpu_time(time1)
    end if

!   Computing pair parameters
    k = 0
    do i = 1, n_sites
      k = k + 1
      c6_ii = neighbor_c6_ii(k)
      r0_i = r0_ii(k)
      alpha0_i = neighbor_alpha0(k)
!     We don't really need these ones but whatever
      neighbor_c6_ij(k) = c6_ii
      r0_ij(k) = r0_i
      do j = 2, n_neigh(i)
        k = k + 1
        c6_jj = neighbor_c6_ii(k)
        alpha0_j = neighbor_alpha0(k)
        if( c6_ii == 0.d0 .or. c6_jj == 0.d0 )then
          neighbor_c6_ij(k) = 0.d0
        else
          neighbor_c6_ij(k) = (2.d0 * c6_ii * c6_jj) / (alpha0_j/alpha0_i * c6_ii + alpha0_i/alpha0_j * c6_jj)
        end if
        r0_j = r0_ii(k)
        r0_ij(k) = r0_i + r0_j
      end do
    end do

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 3:", time2-time1
      call cpu_time(time1)
    end if

!   Compute the TS local energies
    f_damp = 0.d0
    k = 0
    do i = 1, n_sites
      k = k + 1
      do j = 2, n_neigh(i)
        k = k + 1
        if( rjs(k) > rcut_inner .and. rjs(k) < rcut )then
          exp_damp(k) = exp( -d*(rjs(k)/(sR*r0_ij(k)) - 1.d0) )
          f_damp(k) = 1.d0/( 1.d0 + exp_damp(k) )
          energies(i) = energies(i) + neighbor_c6_ij(k) * r6(k) * f_damp(k)
        end if         
      end do
    end do
    energies = -0.5d0 * energies

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing energies:", time2-time1
      call cpu_time(time1)
    end if

!   Compute TS forces
    if( do_forces )then
      forces0 = 0.d0
      virial = 0.d0
!     First, we compute the forces acting on all the "SOAP-neighbors" of atom i
!     (including i itself) due to the gradients of i's Hirshfeld volume wrt the
!     positions of its neighbors 
      pref_force1 = 0.d0
      pref_force2 = 0.d0
      k = 0
      do i = 1, n_sites
        k = k + 1
        i2 = neighbor_species(k)
        do j = 2, n_neigh(i)
          k = k + 1
          if( rjs(k) > rcut_inner .and. rjs(k) < rcut )then
!           For the dependence of C6_ij on the V_A
            pref_force1(i) = pref_force1(i) + neighbor_c6_ij(k) * f_damp(k) * r6(k)
!           For the dependence of f_damp on the V_A
            pref_force2(i) = pref_force2(i) - neighbor_c6_ij(k) * f_damp(k)**2 * rjs(k) * r6(k) * &
                             exp_damp(k) * d / sR / r0_ij(k)**2
          end if
        end do
!       Make sure no NaNs arise for zero Hirshfeld volumes
        if( hirshfeld_v(i) == 0.d0 )then
          pref_force1(i) = 0.d0
          pref_force2(i) = 0.d0
        else
          pref_force1(i) = pref_force1(i) / hirshfeld_v(i)
          pref_force2(i) = pref_force2(i) * r0_ref(i2) / 3.d0 / hirshfeld_v(i)**(2.d0/3.d0)
        end if
      end do
      k = 0
      do i = 1, n_sites
        i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
        do j = 1, n_neigh(i)
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
! NOT SURE IF ALL THESE IFS ARE VERY EFFICIENT
          if( rjs(k) > rcut_inner .and. rjs(k) < rcut )then
!             SOAP neighbors
! MAY NEED TO CHANGE THIS TO ACCOUNT FOR MACHINE PRECISION
            if( .not. all(hirshfeld_v_cart_der(1:3, k) == 0.d0) )then
              this_force(1:3) = hirshfeld_v_cart_der(1:3, k) * ( pref_force1(i) + pref_force2(i) )
              forces0(1:3, j2) = forces0(1:3, j2) + this_force(1:3)
!             Sign is plus because this force is acting on j2. Factor of one is because this is
!             derived from a local energy
!              virial = virial + dot_product(this_force(1:3), xyz(1:3,k))
              do k1 = 1, 3
                do k2 =1, 3
                  virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                end do
              end do
            end if
            if( r0_ij(k) == 0.d0 )then
              this_force(1:3) = 0.d0
            else
              this_force(1:3) = neighbor_c6_ij(k) / rjs(k) * xyz(1:3,k) &
                               * f_damp(k) * ( -r6_der(k) - r6(k) * f_damp(k) * exp_damp(k) &
                                               * d / sR / r0_ij(k) )
            end if
!           There is no net force acting on i2 from its periodic replicas...
            if( j2 /= i2 )then
              forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
            end if
!           ... but the periodic replicas DO contribute to the virial
!           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
!           derived from a pair energy
!            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
            do k1 = 1, 3
              do k2 =1, 3
                virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
              end do
            end do
          end if
        end do
      end do
    end if

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing forces:", time2-time1
    end if
 
!   Clean up
    deallocate( neighbor_c6_ii, neighbor_c6_ij, r0_ij, exp_damp, f_damp, c6_ij_free, r6, is_in_buffer, i_buffer )
    if( do_forces )then
      deallocate( pref_force1, pref_force2, r6_der )
    end if

  end subroutine
!**************************************************************************






!**************************************************************************
  subroutine get_mbd_energy_and_forces( hirshfeld_v, hirshfeld_v_cart_der, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, hirshfeld_v_neigh, &
                                       sR, d, c6_ref, r0_ref, alpha0_ref, do_forces, &
                                       energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v(:), hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), hirshfeld_v_neigh(:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_forces
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), neighbor_c6_ij(:), r0_ii(:), r0_ij(:), &
                           exp_damp(:), f_damp(:), c6_ij_free(:), neighbor_alpha0(:), &
                           pref_force1(:), pref_force2(:), r6(:), r6_der(:), &
                           T_func(:,:), T_func_der(:,:,:,:), h_func(:,:,:,:,:), g_func(:,:,:), &
                           h_func_der(:,:,:,:), g_func_der(:,:,:,:), omegas(:), omega_i(:), &
                           alpha_i(:,:), sigma_i(:,:), sigma_ij(:), T_SR(:,:,:), B_mat(:,:,:), &
                           rjs_H(:), xyz_H(:,:), A_mat(:,:,:), work_arr(:), A_i(:,:), alpha_SCS(:,:), &
                           A_LR(:,:,:), T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), I_mat(:,:), AT(:,:,:), &
                           logIAT(:,:,:), VR(:,:), logMapo(:,:), VRinv(:,:), WR(:), WI(:), VL(:,:), &
                           integrand(:)
    real*8 :: time1, time2, c6_ii, c6_jj, r0_i, r0_j, alpha0_i, alpha0_j, rbuf, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD
    integer, allocatable:: i_buffer(:), ipiv(:)
    integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0, info
    integer :: i, i2, j, j2, k, k2, a, n_in_buffer, c1, c2, c3, lwork
    logical, allocatable :: is_in_buffer(:)
    logical :: do_timing = .false.

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

!   We precompute the C6 coefficients of all the neighbors
    allocate( neighbor_c6_ii(1:n_pairs) )
    allocate( neighbor_c6_ij(1:n_pairs) )
    allocate( r0_ij(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( exp_damp(1:n_pairs) )
    allocate( f_damp(1:n_pairs) )
    allocate( c6_ij_free(1:n_pairs) )
    allocate( r6(1:n_pairs) )
    allocate( is_in_buffer(1:n_pairs) )
    allocate( T_func(1:3*n_sites,1:3*n_sites) )
    allocate( h_func(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( g_func(1:n_sites,1:n_sites,1:11) )
    allocate( omegas(1:11) )
    allocate( omega_i(1:n_pairs) )
    allocate( sigma_i(1:n_sites,1:11) )
    allocate( alpha_i(1:n_sites,1:11) )
    allocate( sigma_ij(1:11) )
    allocate( T_SR(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( B_mat(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( xyz_H(1:3,1:n_pairs) )
    allocate( rjs_H(1:n_pairs) )
    allocate( A_mat(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( work_arr(1:12*n_sites) )
    allocate( ipiv(1:3*n_sites) )
    allocate( A_i(1:3,1:3) )
    allocate( alpha_SCS(1:n_sites,1:11) )
    allocate( A_LR(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( T_LR(1:3*n_sites,1:3*n_sites) )
    allocate( f_damp_SCS(1:n_pairs) )
    allocate( r0_ii_SCS(1:n_pairs) )
    allocate( I_mat(1:3*n_sites,1:3*n_sites) )
    allocate( AT(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( logIAT(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( VR(1:3*n_sites,1:3*n_sites) )
    allocate( logMapo(1:3*n_sites,1:3*n_sites) )
    allocate( VRinv(1:3*n_sites,1:3*n_sites) )
    allocate( WR(1:3*n_sites) )
    allocate( WI(1:3*n_sites) )
    allocate( VL(1:3*n_sites,1:3*n_sites) )
    allocate( integrand(1:11) )
    is_in_buffer = .false.
    if( do_forces )then
      allocate( pref_force1(1:n_sites) )
      allocate( pref_force2(1:n_sites) )
      allocate( r6_der(1:n_pairs) )
      allocate( T_func_der(1:3*n_sites,1:3*n_sites,1:n_sites,1:3) )
      allocate( h_func_der(1:n_sites,1:n_sites,1:n_sites,1:3) )
      allocate( g_func_der(1:n_sites,1:n_sites,1:n_sites,1:3) )
    end if

    if( do_timing) then
      call cpu_time(time1)
    end if
    
!   Frequencies used for integration:
    omega = 0.d0
    do i = 1, 11
      omegas(i) = omega
      omega = omega + 0.4d0 
    end do

!   Check which atoms are in the buffer region
    do k = 1, n_pairs
      j = neighbor_species(k)
      neighbor_c6_ii(k) = c6_ref(j) / (Hartree*Bohr**6)
      r0_ii(k) = r0_ref(j) / Bohr
      neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
      xyz_H(:,k) = xyz(:,k)/Bohr
      rjs_H(k) = rjs(k)/Bohr
    end do
    n_in_buffer = 0
    if( buffer > 0.d0 .or. buffer_inner > 0.d0 )then
      do k = 1, n_pairs
        if( (rjs(k) > rcut-buffer .and. rjs(k) < rcut) .or. &
            (rjs(k) < rcut_inner+buffer_inner .and. rjs(k) > rcut_inner) )then
          n_in_buffer = n_in_buffer + 1
          is_in_buffer(k) = .true.
        end if
      end do
    end if
    allocate( i_buffer(1:n_in_buffer) )
    if( buffer > 0.d0 .or. buffer_inner > 0.d0 )then
      i = 0
      do k = 1, n_pairs
        if( is_in_buffer(k) )then
          i = i + 1
          i_buffer(i) = k
        end if
      end do
    end if

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 1:", time2-time1
      call cpu_time(time1)
    end if

!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    omega_i = (4.d0 * neighbor_c6_ii)/(3*neighbor_alpha0**2) 
    
    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 2:", time2-time1
      call cpu_time(time1)
    end if

    do i = 1, n_sites
      k = (i-1)*n_sites+1
      do j = 1, 11
        alpha_i(i,j) = neighbor_alpha0(k)/(1.d0 + omegas(j)/omega_i(k))**2
        sigma_i(i,j) = (sqrt(2.d0/pi) * alpha_i(i,j)/3.d0)**(1.d0/3.d0)
      end do
    end do

    write(*,*) sigma_i(1:3,1:11)

!   Computing dipole interaction tensor
!   Requires the complete supercell to get correct dimensions for T_func! 3*N_at x 3*N_at, where N_at are atoms in supercell
    T_func = 0.d0
    f_damp = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k) ! NOTE: mod not necessary because we are using only single C60 for now
        if( rjs(k) < rcut )then
          f_damp(k) = 1.d0/( 1.d0 + exp( -d*( rjs(k)/(sR*(r0_ii(n_sites*i+1) + r0_ii(k))) - 1.d0 ) ) )
          do c1 = 1, 3
            T_func(3*(i-1)+c1,3*(j-1)+c1) = (3*xyz_H(c1,k) * xyz_H(c1,k) - rjs_H(k)**2)/rjs_H(k)**5
            T_func(3*(j-1)+c1,3*(i-1)+c1) = T_func(3*(i-1)+c1,3*(j-1)+c1)
            do c2 = c1+1, 3
              T_func(3*(i-1)+c1,3*(j-1)+c2) = (3*xyz_H(c1,k) * xyz_H(c2,k) - rjs_H(k)**2)/rjs_H(k)**5
              T_func(3*(j-1)+c1,3*(i-1)+c2) = T_func(3*(i-1)+c1,3*(j-1)+c2)
              T_func(3*(i-1)+c2,3*(j-1)+c1) = T_func(3*(i-1)+c1,3*(j-1)+c2)
              T_func(3*(j-1)+c2,3*(i-1)+c1) = T_func(3*(i-1)+c1,3*(j-1)+c2)
            end do
          end do
        end if
      end do
    end do

    write(*,*) "f_damp:", f_damp(1:6)
    write(*,*) "Size of neighbor list:", size(neighbors_list)
    write(*,*) "Size:", size(T_func, 1), size(T_func, 2)
    write(*,*) "T_ij:", T_func(1,1:9)
    write(*,*) "T_ij:", T_func(2,1:9)
    write(*,*) "T_ij:", T_func(3,1:9)
    write(*,*) "T_ij:", T_func(4,1:9)
    write(*,*) "T_ij:", T_func(5,1:9)
    write(*,*) "T_ij:", T_func(6,1:9)
    write(*,*) "T_ij:", T_func(7,1:9)
    write(*,*) "T_ij:", T_func(8,1:9)
    write(*,*) "T_ij:", T_func(9,1:9)

    g_func = 0.d0
    h_func = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        if( rjs(k) < rcut )then
          sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
          g_func(i,j,:) = erf(rjs_H(k)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k)/sigma_ij) * exp(-rjs_H(k)**2.d0/sigma_ij**2)
          g_func(j,i,:) = g_func(i,j,:)
          do c1 = 1, 3
            do c2 = 1, 3
              h_func(i,j,c1,c2,:) = 4.d0/sqrt(pi) * (rjs_H(k)/sigma_ij)**3 * &
                                    xyz_H(c1,k)*xyz_H(c2,k)/rjs_H(k)**5 * exp(-rjs_H(k)**2/sigma_ij**2)
              h_func(i,j,c2,c1,:) = h_func(i,j,c1,c2,:)
              h_func(j,i,c1,c2,:) = h_func(i,j,c1,c2,:)
              h_func(j,i,c2,c1,:) = h_func(i,j,c1,c2,:)
            end do
          end do 
        end if
      end do
    end do

    write(*,*) c6_ref(1)/(Hartree*Bohr**6), alpha0_ref(1)/Bohr**3
    write(*,*) "g_func:", g_func(1,1:6,1)
    write(*,*) "g_func:", g_func(2,1:6,1)
    write(*,*) "g_func:", g_func(3,1:6,1)
    write(*,*) "g_func:", g_func(4,1:6,1)
    write(*,*) "g_func:", g_func(5,1:6,1)
    write(*,*) "g_func:", g_func(6,1:6,1)
    write(*,*) "h_func:", h_func(1,1:6,1,1,1)
    write(*,*) "h_func:", h_func(2,1:6,1,1,1)
    write(*,*) "h_func:", h_func(3,1:6,1,1,1)
    write(*,*) "h_func:", h_func(4,1:6,1,1,1)
    write(*,*) "h_func:", h_func(5,1:6,1,1,1)
    write(*,*) "h_func:", h_func(6,1:6,1,1,1)
    write(*,*) "h_func max:", maxval(h_func)

    T_SR = 0.d0
    B_mat = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = i+1, n_sites
        k = k+1
        j = neighbors_list(k)
        do c1 = 1, 3
          B_mat(3*(i-1)+c1,3*(i-1)+c1,:) = 1.d0/alpha_i(i,:)
          do c2 = 1, 3
            T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = (1-f_damp(k)) * (-T_func(3*(i-1)+c1,3*(j-1)+c2) * &
                                          g_func(i,j,:) + h_func(i,j,c1,c2,:))
            T_SR(3*(j-1)+c1,3*(i-1)+c2,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
            T_SR(3*(i-1)+c2,3*(j-1)+c1,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
            T_SR(3*(j-1)+c2,3*(i-1)+c1,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
          end do
        end do
      end do
    end do
    B_mat = B_mat + T_SR

    write(*,*) "T_SR:", T_SR(1,1:9,1)
    write(*,*) "T_SR:", T_SR(2,1:9,1)
    write(*,*) "T_SR:", T_SR(3,1:9,1)
    write(*,*) "T_SR:", T_SR(4,1:9,1)
    write(*,*) "T_SR:", T_SR(5,1:9,1)
    write(*,*) "T_SR:", T_SR(6,1:9,1)
    write(*,*) "T_SR:", T_SR(7,1:9,1)
    write(*,*) "T_SR:", T_SR(8,1:9,1)
    write(*,*) "T_SR:", T_SR(9,1:9,1)
    write(*,*) "B_mat:", B_mat(1,1:9,1)
    write(*,*) "B_mat:", B_mat(2,1:9,1)
    write(*,*) "B_mat:", B_mat(3,1:9,1)
    write(*,*) "B_mat:", B_mat(4,1:9,1)
    write(*,*) "B_mat:", B_mat(5,1:9,1)
    write(*,*) "B_mat:", B_mat(6,1:9,1)
    write(*,*) "B_mat:", B_mat(7,1:9,1)
    write(*,*) "B_mat:", B_mat(8,1:9,1)
    write(*,*) "B_mat:", B_mat(9,1:9,1)

    A_mat = B_mat

    do i = 1, 11
      call dgetrf(3*n_sites, 3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, info)
      call dgetri(3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
    end do

    write(*,*) "A_mat:", A_mat(1,1:6,1)
    write(*,*) "A_mat:", A_mat(2,1:6,1)
    write(*,*) "A_mat:", A_mat(3,1:6,1)
    write(*,*) "A_mat:", A_mat(4,1:6,1)
    write(*,*) "A_mat:", A_mat(5,1:6,1)
    write(*,*) "A_mat:", A_mat(6,1:6,1)

    do k = 1, 11
      do i = 1, n_sites
        A_i = 0.d0
        do j = 1, n_sites
          A_i = A_i + A_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,k)
        end do
        alpha_SCS(i,k) = A_i(1,1)+A_i(2,2)+A_i(3,3)
      end do
    end do

    A_LR = 0.d0

    do k = 1, 11
      do i = 1, n_sites
        do c1 = 1, 3
          A_LR(3*(i-1)+c1,3*(i-1)+c1,k) = alpha_SCS(i,k)
        end do
      end do
    end do

    write(*,*) "alpha_SCS:", alpha_SCS(:,1)
    write(*,*) "A_LR:", A_LR(1,1:6,1)
    write(*,*) "A_LR:", A_LR(2,1:6,1)
    write(*,*) "A_LR:", A_LR(3,1:6,1)
    write(*,*) "A_LR:", A_LR(4,1:6,1)
    write(*,*) "A_LR:", A_LR(5,1:6,1)
    write(*,*) "A_LR:", A_LR(6,1:6,1)

    T_LR = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j = i+1, n_sites
        k = k+1
        T_LR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = f_damp(k) * T_func(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
        T_LR(3*(j-1)+1:3*(j-1)+3,3*(i-1)+1:3*(i-1)+3) = T_LR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
      end do
    end do

    write(*,*) "T_LR:", T_LR(1,1:6)
    write(*,*) "T_LR:", T_LR(2,1:6)
    write(*,*) "T_LR:", T_LR(3,1:6)
    write(*,*) "T_LR:", T_LR(4,1:6)
    write(*,*) "T_LR:", T_LR(5,1:6)
    write(*,*) "T_LR:", T_LR(6,1:6)

    I_mat = 0.d0
    do i = 1, 3*n_sites
      I_mat(i,i) = 1.d0
    end do

    AT = 0.d0
    do k = 1, 11
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_mat(:,:,k), 3*n_sites, T_LR, 3*n_sites, &
                 0.d0, AT(:,:,k), 3*n_sites) 
    end do

    write(*,*) "AT:", AT(1,1:6,1)
    write(*,*) "AT:", AT(2,1:6,1)
    write(*,*) "AT:", AT(3,1:6,1)
    write(*,*) "AT:", AT(4,1:6,1)
    write(*,*) "AT:", AT(5,1:6,1)
    write(*,*) "AT:", AT(6,1:6,1)

    do k = 1,11
      WR = 0.d0
      VL = 0.d0
      VR = 0.d0
      call dgeev('n', 'v', 3*n_sites, I_mat-AT(:,:,k), 3*n_sites, WR, WI, VL, 3*n_sites, VR, 3*n_sites, &
                 work_arr, 12*n_sites, info) 
      logMapo = 0.d0
      write(*,*) "WR:"
      do i = 1, 3*n_sites
        write(*,*) WR(i)
        logMapo(i,i) = log(WR(i))
      end do
      VRinv = VR
      VL = 0.d0
      call dgetrf(3*n_sites, 3*n_sites, VRinv, 3*n_sites, ipiv, info)
      call dgetri(3*n_sites, VRinv, 3*n_sites, ipiv, work_arr, 12*n_sites, info)
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, VR, 3*n_sites, logMapo, 3*n_sites, &
                 0.d0, VL, 3*n_sites)
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, VL, 3*n_sites, VRinv, 3*n_sites, &
                 0.d0, logIAT(:,:,k), 3*n_sites)
    end do

    write(*,*) "logIAT:", logIAT(1,1:6,1)
    write(*,*) "logIAT:", logIAT(2,1:6,1)
    write(*,*) "logIAT:", logIAT(3,1:6,1)
    write(*,*) "logIAT:", logIAT(4,1:6,1)
    write(*,*) "logIAT:", logIAT(5,1:6,1)
    write(*,*) "logIAT:", logIAT(6,1:6,1)

    integrand = 0.d0
    do k = 1,11
      do i = 1,3*n_sites
        integrand(k) = integrand(k) + logIAT(i,i,k)
      end do 
    end do

!    NOTE: Include misc.mod before commenting this out
!    call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
!    E_MBD = integral/(2.d0*pi)

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: creating pair quantities 3:", time2-time1
      call cpu_time(time1)
    end if

!   Compute the TS local energies

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing energies:", time2-time1
      call cpu_time(time1)
    end if

!   Compute MBD forces
    if( do_forces )then
      forces0 = 0.d0
      virial = 0.d0
    end if

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing forces:", time2-time1
    end if
 
!   Clean up
    deallocate( neighbor_c6_ii, neighbor_c6_ij, r0_ij, exp_damp, f_damp, c6_ij_free, r6, is_in_buffer, i_buffer, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, A_LR, r0_ii_SCS, f_damp_SCS, I_mat, AT, logIAT, logMapo, VR, VRinv, WR, WI, VL, integrand )
    if( do_forces )then
      deallocate( pref_force1, pref_force2, r6_der, T_func_der, h_func_der, g_func_der )
    end if

  end subroutine
!**************************************************************************





!**************************************************************************
  subroutine get_vdw_ref_params(element, C6, R0, alpha0, rank)

    implicit none

    character*8, intent(in) :: element
    integer, intent(in) :: rank
    real*8, intent(out) :: C6, R0, alpha0
    real*8 :: Hartree = 27.211386024367243d0, Bohr = 0.5291772105638411d0

    C6 = 0.d0
    R0 =  0.d0
    alpha0 = 0.d0

!   These should be in the correct units (enegy in eV, distances in Angstrom)
!   Variables to help convert between units are provided in this subroutine

    if( element == "H" )then
!     This is the value provided by VASP, for which they give "private comm."
!     as reference in the TS implementation paper:
      R0 = 1.64d0
!     These values are given by Chu and Dalgarno (J Chem Phys 121, 4083 [2004])
      alpha0 = 4.5d0 * Bohr**3
      C6 = 6.5d0 * Hartree * Bohr**6
    else if( element == "C" )then
!     This is the value provided by VASP, for which they give "private comm."
!     as reference in the TS implementation paper:
      R0 = 1.900d0
!     This is the one given by Grimme (J Comput Chem 27, 1787 [2006]):
!      R0 = 1.452d0
!     These values are given by Chu and Dalgarno (J Chem Phys 121, 4083 [2004])
      alpha0 = 12.d0 * Bohr**3
      C6 = 46.6d0 * Hartree * Bohr**6
    else
      if( rank == 0 )then
        write(*,*)'                                       |'
        write(*,*)'WARNING: default vdW reference parame- |'
        write(*,*)'ters not available for element:        |'
        write(*,*)'                                       |'
        write(*,'(1X,A8,A)') element, '                               |'
        write(*,*)'                                       |'
        write(*,*)'You can safely disregard this warning  |'
        write(*,*)'if your potential does not use van der |'
        write(*,*)'Waals corrections. Otherwise, or if you|'
        write(*,*)'want to overrride the defaults, you    |'
        write(*,*)'need to provide your own definitions of|'
        write(*,*)'vdw_c6_ref, vdw_r0_ref and             |'
        write(*,*)'vdw_alpha0_ref in the input file.      |'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
      end if
    end if

  end subroutine
!**************************************************************************


end module
