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
    real*8, intent(in) :: hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:), hirshfeld_v(:), hirshfeld_v_neigh(:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_forces
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), &
                           f_damp(:), neighbor_alpha0(:), &
                           pref_force1(:), pref_force2(:), &
                           T_func(:), h_func(:,:), g_func(:,:), &
                           omegas(:), omega_i(:), &
                           alpha_i(:,:), sigma_i(:,:), sigma_ij(:), T_SR(:,:,:), B_mat(:,:,:), &
                           rjs_H(:), xyz_H(:,:), A_mat(:,:,:), work_arr(:), A_i(:,:), alpha_SCS(:,:), &
                           A_LR(:,:,:), T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), I_mat(:,:), AT(:,:,:), &
                           logIAT(:,:,:), VR(:,:), logMapo(:,:), VRinv(:,:), WR(:), WI(:), VL(:,:), &
                           integrand(:), dT(:,:), dT_SR(:,:,:,:,:), f_damp_der(:,:), &
                           g_func_der(:,:,:), h_func_der(:,:,:), dA_mat(:,:,:,:,:), &
                           dalpha(:,:,:,:), dA_LR(:,:,:,:,:), dT_LR(:,:,:,:), f_damp_der_SCS(:,:,:), &
                           invIAT(:,:,:), G_mat(:,:,:,:,:), force_integrand(:,:,:), forces_MBD(:,:), &
                           coeff_h_der(:), terms(:), dT_SR_A_mat(:,:), dT_SR_v(:,:,:,:,:), &
                           dB_mat(:,:,:,:,:), dB_mat_v(:,:,:,:,:), &
                           dv_i(:), dv_j(:), &
                           coeff_der(:,:,:,:,:), coeff_fdamp(:,:,:,:,:), dg(:), &
                           dh(:), hirshfeld_v_cart_der_H(:,:), AT_n(:,:,:,:), A_LR_k(:,:,:,:), AT_k(:,:,:,:), &
                           integrand_k(:), E_MBD_k(:), alpha_k(:,:), sigma_k(:,:), s_i(:), s_j(:), &
                           T_SR_i(:,:,:), B_mat_i(:,:,:), alpha_SCS_i(:,:), A_mat_i(:,:,:), A_LR_i(:,:,:), &
                           T_LR_i(:,:), AT_i(:,:,:), series(:,:,:), integrand_2(:), I_mat_n(:,:), &
                           logIAT_n(:,:,:), WR_n(:), WI_n(:), VL_n(:,:), VR_n(:,:), VRinv_n(:,:), &
                           logMapo_n(:,:), dB_mat_n(:,:,:,:), dA_mat_n(:,:,:,:), dBA_n(:,:), &
                           dA_LR_n(:,:,:,:), dalpha_n(:,:,:), f_damp_der_ij_n(:), f_damp_der_SCS_ij_n(:), &
                           dT_LR_n(:,:,:), force_integrand_n(:,:), forces_MBD_k(:,:), invIAT_n(:,:,:), &
                           G_mat_n(:,:,:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, dist, f_damp_SCS_ij
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, n_order, n_tot
    integer :: i, i2, j, j2, k, k2, k3, a, a2, c1, c2, c3, lwork, b, p, q, kf
    logical :: do_timing = .false., do_hirshfeld_gradients = .true., nonlocal = .false.

!   Change these to be input variables (NOTE THAT THEY ARE IN ANGSTROMS!):
    write(*,*) "rcut", rcut
!    rcut_vdw = 4.d0
    n_order = 100

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!   Number of neighbors check for finite difference:
!    do i = 1, n_sites
!      write(*,*) "n_neigh:", i, n_neigh(i)
!    end do
!    write(*,*) "Size of neighbors_list:", size(neighbors_list)
!    k = 0
!    do i = 1, n_sites
!      write(*,*) neighbors_list(k+1:k+n_neigh(i))
!      k = k+n_neigh(i)
!      write(*,*) i, hirshfeld_v(i)
!    end do

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

!   Allocate all the necessary stuff
    allocate( neighbor_c6_ii(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( f_damp(1:n_pairs) )
    allocate( T_func(1:9*n_pairs) )
    allocate( h_func(1:9*n_pairs,1:11) )
    allocate( g_func(1:n_pairs,1:11) )
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
    allocate( dT(1:9*n_pairs,1:3) )
    allocate( dT_SR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( f_damp_der(1:n_pairs,1:3) )
    allocate( g_func_der(1:n_pairs,1:3,1:11) )
    allocate( h_func_der(1:9*n_pairs,1:3,1:11) )
    allocate( dA_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dalpha(1:n_sites,1:n_sites,1:3,1:11) )
    allocate( dA_LR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dT_LR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3) )
    allocate( f_damp_der_SCS(1:n_pairs,1:n_sites,1:3) )
    allocate( invIAT(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( G_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( force_integrand(1:n_sites,1:3,1:11) )
    allocate( forces_MBD(1:n_sites,1:3) )
    allocate( coeff_h_der(1:11) )
    allocate( terms(1:11) )
    allocate( dT_SR_A_mat(1:3*n_sites,1:3*n_sites) )
    allocate( dT_SR_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dB_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dB_mat_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dv_i(1:3) )
    allocate( dv_j(1:3) )
    allocate( coeff_der(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( coeff_fdamp(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( dg(1:11) )
    allocate( dh(1:11) )
    allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
    allocate( A_LR_k(1:3*n_sites,1:3*n_sites,1:11,1:n_sites) )
    allocate( AT_k(1:3*n_sites,1:3*n_sites,1:11,1:n_sites) )
    allocate( integrand_k(1:11) )
    allocate( E_MBD_k(1:n_sites) )
    allocate( alpha_k(1:n_pairs,1:11) )
    allocate( sigma_k(1:n_pairs,1:11) )
    allocate( s_i(1:11) )
    allocate( s_j(1:11) )
    allocate( integrand_2(1:11) )
    allocate( f_damp_der_ij_n(1:3) )
    allocate( f_damp_der_SCS_ij_n(1:3) )
    allocate( force_integrand_n(1:3,1:11) )
    allocate( forces_MBD_k(1:n_sites,1:3) )

    if( do_timing) then
      call cpu_time(time1)
    end if
    
!   Frequencies used for integration:
    omega = 0.d0
    do i = 1, 11
      omegas(i) = omega
      omega = omega + 0.4d0 
    end do

    do k = 1, n_pairs
      j = neighbor_species(k)
      neighbor_c6_ii(k) = c6_ref(j) / (Hartree*Bohr**6)
      r0_ii(k) = r0_ref(j) / Bohr
      neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
      xyz_H(:,k) = xyz(:,k)/Bohr
      rjs_H(k) = rjs(k)/Bohr
      hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,k)*Bohr
    end do

!    do i = 1, n_neigh(i)
!      write(*,*) "i, neighbors_list(i), dv_i", i, neighbors_list(i), rjs(i), hirshfeld_v_cart_der_H(:,i)
!    end do

!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    omega_i = (4.d0 * neighbor_c6_ii)/(3.d0*neighbor_alpha0**2)

    k2 = 1
    do i = 1, n_sites
      do k = 1, 11
        alpha_i(i,k) = neighbor_alpha0(k2)/(1.d0 + omegas(k)**2/omega_i(k2)**2)
        sigma_i(i,k) = (sqrt(2.d0/pi) * alpha_i(i,k)/3.d0)**(1.d0/3.d0)
      end do
      k2 = k2+n_neigh(i)
    end do

    do k = 1, 11
      alpha_k(:,k) = neighbor_alpha0/(1.d0 + omegas(k)**2/omega_i**2)
      sigma_k(:,k) = (sqrt(2.d0/pi) * alpha_k(:,k)/3.d0)**(1.d0/3.d0)
    end do

!   Computing dipole interaction tensor
!   Requires the complete supercell to get correct dimensions for T_func! 3*N_at x 3*N_at, where N_at are atoms in supercell
    T_func = 0.d0
    f_damp = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      r_vdw_i = r0_ii(k)
      do j2 = 2, n_neigh(i)
        k = k+1
        r_vdw_j = r0_ii(k)
        if( rjs(k) < rcut )then
          f_damp(k) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2 + 1
              if (c1 == c2) then
                T_func(k2) = (3*xyz_H(c1,k) * xyz_H(c1,k) - rjs_H(k)**2)/rjs_H(k)**5
              else
                T_func(k2) = (3*xyz_H(c1,k) * xyz_H(c2,k))/rjs_H(k)**5
              end if
            end do
          end do
        end if
      end do
    end do

    g_func = 0.d0
    h_func = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      s_i = sigma_k(k,:)
      do j2 = 2, n_neigh(i)
        k = k+1
        s_j = sigma_k(k,:)
        j = neighbors_list(k)
        if( rjs(k) < rcut )then
          sigma_ij = sqrt(s_i**2 + s_j**2)
          g_func(k,:) = erf(rjs_H(k)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k)/sigma_ij) * exp(-rjs_H(k)**2.d0/sigma_ij**2)
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2+1
              h_func(k2,:) = 4.d0/sqrt(pi) * (rjs_H(k)/sigma_ij)**3 * &
                                    xyz_H(c1,k)*xyz_H(c2,k)/rjs_H(k)**5 * exp(-rjs_H(k)**2/sigma_ij**2)
            end do
          end do 
        end if
      end do
    end do

    T_SR = 0.d0
    B_mat = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do c1 = 1, 3
        B_mat(3*(i-1)+c1,3*(i-1)+c1,:) = 1.d0/alpha_k(k,:)
      end do
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        if( rjs(k) < rcut )then
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2+1
              T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(k2) * &
                                              g_func(k,:) + h_func(k2,:))
            end do
          end do
        end if
      end do
    end do
    B_mat = B_mat + T_SR

    A_mat = B_mat

    do i = 1, 11
      call dgetrf(3*n_sites, 3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, info)
      call dgetri(3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
    end do

    do k = 1, 11
      do i = 1, n_sites
        A_i = 0.d0
        do j = 1, n_sites
          A_i = A_i + A_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,k)
        end do
        alpha_SCS(i,k) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
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

    r0_ii_SCS = 0.d0

    do k = 1, n_pairs
      j = neighbors_list(k)
      r0_ii_SCS(k) = r0_ii(k) * (alpha_SCS(j,1)/neighbor_alpha0(k))**(1.d0/3.d0)
    end do

    f_damp_SCS = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      r_vdw_i = r0_ii_SCS(k)
      do j2 = 2, n_neigh(i)
        k=k+1
        if (rjs(k) < rcut) then
          r_vdw_j = r0_ii_SCS(k)
          j = neighbors_list(k)
          f_damp_SCS(k) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
        end if
      end do
    end do

    T_LR = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        if (rjs(k) < rcut) then
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2+1
              T_LR(3*(i-1)+c1,3*(j-1)+c2) = f_damp_SCS(k) * T_func(k2)
              T_LR(3*(j-1)+c1,3*(i-1)+c2) = T_LR(3*(i-1)+c1,3*(j-1)+c2)
            end do
          end do
        end if
      end do
    end do

    I_mat = 0.d0
    do i = 1, 3*n_sites
      I_mat(i,i) = 1.d0
    end do

    AT = 0.d0
    do k = 1, 11
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR(:,:,k), 3*n_sites, T_LR, 3*n_sites, &
                0.d0, AT(:,:,k), 3*n_sites)
    end do

    dT = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        if (rjs(k) < rcut) then
          k2 = 9*(k-1)
          do c1 = 1,3
            do c2 = 1,3
              k2 = k2+1
              do c3 = 1,3
                dT(k2,c3) = (-15.d0 * xyz_H(c1,k) * xyz_H(c2,k) * xyz_H(c3,k))/rjs_H(k)**7
                if (c1 == c2) then
                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c3,k)
                end if
                if (c2 == c3) then
                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c1,k)
                end if
                if (c1 == c3) then
                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c2,k)
                end if
              end do
            end do
          end do
        end if
      end do
    end do

    f_damp_der = 0.d0
    g_func_der = 0.d0
    h_func_der = 0.d0
    dT_SR = 0.d0

    do a = 1, n_sites
      k = 0
      do i = 1, n_sites
        k = k+1
        r_vdw_i = r0_ii(k)
        do j2 = 2, n_neigh(i)
          k = k+1
          r_vdw_j = r0_ii(k)
          j = neighbors_list(k)
          if (rjs(k) < rcut) then
            if (a == i .or. a == j) then
              sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
              do c3 = 1, 3
                f_damp_der(k,c3) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k)**2 * &
                                     exp( -d*(rjs_H(k)/(sR*(r_vdw_i + &
                                     r_vdw_j)) - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
                g_func_der(k,c3,:) = 4.d0/sqrt(pi) * rjs_H(k)/sigma_ij**3 * xyz_H(c3,k) * &
                                       exp(-rjs_H(k)**2/sigma_ij**2)
!                if (a == 1 .and. c3 == 1 .and. i == 1 .and. j == 33) then
!                  write(*,*) "a, c3, i, j, f_damp_der, g_func_der", a, c3, i, j, f_damp_der(k,c3), g_func_der(k,c3,1)
!                end if
!                if (a == 1 .and. c3 == 1 .and. i == 33 .and. j == 1) then
!                  write(*,*) "a, c3, i, j, f_damp_der, g_func_der", a, c3, i, j, f_damp_der(k,c3), g_func_der(k,c3,1)
!                end if
                k2 = 9*(k-1)
                do c1 = 1, 3
                  do c2 = 1, 3
                    k2 = k2+1
                    coeff_h_der = 4.d0/sqrt(pi) * 1.d0/sigma_ij**3 * exp(-rjs_H(k)**2/sigma_ij**2)
                    terms = 0.d0
                    if (c1 == c3) then
                      terms = terms + xyz_H(c2,k)/rjs_H(k)**2
                    end if
                    if (c2 == c3) then
                      terms = terms + xyz_H(c1,k)/rjs_H(k)**2
                    end if
                    terms = terms + -2.d0*(xyz_H(c1,k)*xyz_H(c2,k)*xyz_H(c3,k))/rjs_H(k)**2 * &
                            (1.d0/rjs_H(k)**2 + 1.d0/sigma_ij**2)
                    h_func_der(k2,c3,:) = coeff_h_der * terms
                    dT_SR(3*(i-1)+c1,3*(j-1)+c2,a,c3,:) = f_damp_der(k,c3) * T_func(k2) * &
                                                          g_func(k,:) - (1.d0 - f_damp(k)) * dT(k2,c3) * &
                                                          g_func(k,:) - g_func_der(k,c3,:) * (1.d0 - f_damp(k)) * &
                                                          T_func(k2) - f_damp_der(k,c3) * h_func(k2,:) + &
                                                          h_func_der(k2,c3,:) * (1.d0 - f_damp(k))
                  end do
                end do
              end do
              if (a == i) then
                dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:) = -dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:)
              end if
            end if
          end if
        end do
      end do
    end do

!    write(*,*) "dT_SR:"
!    write(*,*) dT_SR(3*(33-1)+1,1,1,1,1)
!    write(*,*) dT_SR(1,3*(33-1)+1,1,1,1)

    dB_mat = dT_SR
    if (do_hirshfeld_gradients) then
      dB_mat_v = 0.d0
      dT_SR_v = 0.d0
      coeff_fdamp = 0.d0
      coeff_der = 0.d0
      k = 0
      do i = 1, n_sites
        k = k+1
        r_vdw_i = r0_ii(k)
        do j2 = 2, n_neigh(i)
          k = k+1
          r_vdw_j = r0_ii(k)
          j = neighbors_list(k)
          if (rjs(k) < rcut) then
            sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
            S_vdW_ij = sR*(r_vdw_i + r_vdw_j)
            exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
            k2 = 9*(k-1)
            do c1 = 1, 3
              do c2 = 1, 3
                k2 = k2+1
                coeff_fdamp(i,j,c1,c2,:) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                           (-T_func(k2) * g_func(k,:) + h_func(k2,:))
                dg = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * (rjs_H(k)/sigma_ij)**2 * &
                     (-rjs_H(k)/(3.d0 * sigma_ij**3))
                dh = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * xyz_H(c1,k)*xyz_H(c2,k) / &
                     (sigma_ij**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij)**2)
                coeff_der(i,j,c1,c2,:) = (1.d0-f_damp(k)) * (-T_func(k2)*dg+dh)
              end do
            end do
          end if
        end do
      end do

      k2 = 0
      k3 = 0
      do i = 1, n_sites
        r_vdw_i = r0_ii(k3+1)
        do a2 = 1, n_neigh(i)
          k2 = k2+1
          a = neighbors_list(k2)
          do c3 = 1, 3
            do c1 = 1, 3
              dB_mat_v(3*(i-1)+c1,3*(i-1)+c1,a,c3,:) = -1.d0/(hirshfeld_v(i)*alpha_i(i,:)) * &
                            hirshfeld_v_cart_der_H(c3,k2)
            end do
            do k = 1, 11
              do j2 = 2, n_neigh(i)
                if (rjs(k3+j2) < rcut) then
                  j = neighbors_list(k3+j2)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                        coeff_der(i,j,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                        coeff_fdamp(i,j,c1,c2,k) * r_vdw_i/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
                    end do
                  end do
                end if      
              end do
            end do
          end do
        end do
        k3 = k3 + n_neigh(i)
      end do
      k2 = 0
      k3 = 0
      do j = 1, n_sites
        r_vdw_j = r0_ii(k3+1)
        do a2 = 1, n_neigh(j)
          k2 = k2+1
          a = neighbors_list(k2)
          do c3 = 1, 3
            do k = 1, 11
              do i2 = 2, n_neigh(j)
                if (rjs(k3+i2) < rcut) then
                  i = neighbors_list(k3+i2)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                        coeff_der(i,j,c1,c2,k) * (sigma_i(j,k)**2/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2))
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                        coeff_fdamp(i,j,c1,c2,k) * r_vdw_j/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2)
                    end do
                  end do
                end if
              end do
            end do
          end do
        end do
        k3 = k3 + n_neigh(j)
      end do

      dB_mat_v = dB_mat_v + dT_SR_v

      dB_mat = dB_mat + dB_mat_v
    end if

!    write(*,*) "dT_SR_v:"
!    write(*,*) dT_SR_v(3*(33-1)+1,1,1,1,1)
!    write(*,*) dT_SR_v(1,3*(33-1)+1,1,1,1)

!    write(*,*) "dB_mat:"
!    write(*,*) dB_mat(3*(33-1)+1,1,1,1,1)
!    write(*,*) dB_mat(1,3*(33-1)+1,1,1,1)
!    do p = 1, 3*n_sites
!      write(*,*) dB_mat(p,:,1,1,1)
!    end do

    if (nonlocal) then
      logIAT = 0.d0
      do k = 1,11
        WR = 0.d0
        WI = 0.d0
        VL = 0.d0
        VR = 0.d0
        call dgeev('n', 'v', 3*n_sites, I_mat-AT(:,:,k), 3*n_sites, WR, WI, VL, 3*n_sites, VR, 3*n_sites, &
                   work_arr, 12*n_sites, info)
        logMapo = 0.d0
        do i = 1, 3*n_sites
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

      integrand = 0.d0
      do k = 1,11
        do i = 1,3*n_sites
          integrand(k) = integrand(k) + logIAT(i,i,k)
        end do 
      end do

!     Local energy test:
!      E_MBD_k = 0.d0
!      do i = 1, n_sites
!        integrand = 0.d0
!        do k = 1,11
!          do i2 = 1,3*n_sites
!            integrand(k) = integrand(k) + logIAT(i2,i2,k)
!          end do
!          integrand(k) = alpha_SCS(i,k)/sum(alpha_SCS(:,k)) * integrand(k)
!        end do
!        integral = 0.d0
!        call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
!        E_MBD_k(i,1) = integral/(2.d0*pi)
!        E_MBD_k(i,1) = E_MBD_k(i,1) * 27.211386245988
!        write(*,*) "E_MBD:", i, E_MBD_k(i,1)
!      end do
!      write(*,*) "Total MBD energy:", sum(E_MBD_k(:,1))


      call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)

      E_MBD = integral/(2.d0*pi)

      ! Conversion to eV:
      E_MBD = E_MBD * 27.211386245988
      write(*,*) "E_MBD:", E_MBD

!      write(*,*) "Local energy attempt:"
!      do i = 1, n_sites
!        write(*,*) "E_MBD local for atom:", i, alpha_SCS(i,1)/sum(alpha_SCS(:,1)) * E_MBD
!      end do

      dA_mat = 0.d0
      do k = 1, 11
        do a = 1, n_sites
          do c3 = 1, 3
            dT_SR_A_mat = 0.d0
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, dB_mat(:,:,a,c3,k), 3*n_sites, A_mat(:,:,k), 3*n_sites, &
                       0.d0, dT_SR_A_mat, 3*n_sites)
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, -A_mat(:,:,k), 3*n_sites, dT_SR_A_mat, 3*n_sites, &
                       0.d0, dA_mat(:,:,a,c3,k), 3*n_sites)
          end do
        end do
      end do

!      write(*,*) "dA_mat:"
!      do c1 = 1, 3
!        write(*,*) dA_mat(c1,3*(29-1)+1:3*(29-1)+3,1,1,1)
!      end do

      dalpha = 0.d0
      do k = 1, 11
        do i = 1, n_sites
          do a = 1, n_sites
            do c3 = 1, 3
              A_i = 0.d0
              do j = 1, n_sites
                A_i = A_i + dA_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,c3,k)
              end do
              do j2 = 1, 3
                dalpha(i,a,c3,k) = dalpha(i,a,c3,k) + A_i(j2,j2)
              end do
              dalpha(i,a,c3,k) = 1.d0/3.d0 * dalpha(i,a,c3,k)
            end do
          end do
        end do
      end do

!      write(*,*) "dalpha:", dalpha(:,1,1,1)

      dA_LR = 0.d0

      do k = 1, 11
        do a = 1, n_sites
          do c3 = 1, 3
            do i = 1, n_sites
              do c1 = 1, 3
                dA_LR(3*(i-1)+c1,3*(i-1)+c1,a,c3,k) = dalpha(i,a,c3,k)
              end do
            end do
          end do
        end do
      end do

!      do a = 1, 3
!        write(*,*) a
!        do i2 = 1, 6
!          write(*,*) "dA_LR:", dA_LR(i2,1:6,a,1,1)
!        end do
!      end do

      f_damp_der_SCS = 0.d0
      f_damp_der = 0.d0 ! This is cleared so we can recalculate it with SCS values
      dT_LR = 0.d0
      do a = 1, n_sites
        k = 0
        do i = 1, n_sites
          k = k+1
          r_vdw_i = r0_ii_SCS(k)
          do j2 = 2, n_neigh(i)
            k = k+1
            if (rjs(k) < rcut) then
              j = neighbors_list(k)
              r_vdw_j = r0_ii_SCS(k)
              R_vdW_SCS_ij = r_vdw_i + r_vdw_j
              S_vdW_ij = sR*R_vdW_SCS_ij
              do c3 = 1, 3
                dS_vdW_ij = sR/3.d0 * ( r0_ii_SCS(n_sites*(i-1)+1)/alpha_SCS(i,1) * dalpha(i,a,c3,1) + &
                                        r0_ii_SCS(k)/alpha_SCS(j,1) * dalpha(j,a,c3,1) )
                f_damp_der_SCS(k,a,c3) = -(d*rjs_H(k))/S_vdW_ij**2 * f_damp_SCS(k)**2 * &
                                           exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                k2 = 9*(k-1)
                do c1 = 1, 3
                  do c2 = 1, 3
                    k2 = k2+1
                    dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) + &
                                      T_func(k2) * f_damp_der_SCS(k,a,c3) 
                  end do
                end do
              end do
              if (a == i .or. a == j) then
                do c3 = 1, 3
                f_damp_der(k,c3) = d/S_vdW_ij * f_damp_SCS(k)**2 * &
                                     exp( -d*(rjs_H(k)/S_vdW_ij - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
                  k2 = 9*(k-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k2 = k2+1
                      if (a == i) then
                        dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) - &
                                           T_func(k2) * f_damp_der(k,c3) - &
                                           dT(k2,c3) * f_damp_SCS(k)
                      else if (a == j) then
                        dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) + &
                                           T_func(k2) * f_damp_der(k,c3) + &
                                           dT(k2,c3) * f_damp_SCS(k)
                      end if
                    end do
                  end do
                end do
              end if
            end if
          end do
        end do
      end do

!      write(*,*) "dT_LR:"
!      do c1 = 1, 3
!        write(*,*) dT_LR(c1,3*(29-1)+1:3*(29-1)+3,1,1)
!      end do

      invIAT = 0.d0
      G_mat = 0.d0
      force_integrand = 0.d0
 
      invIAT = 0.d0
      do k = 1, 11
        invIAT(:,:,k) = I_mat - AT(:,:,k)
        call dgetrf(3*n_sites, 3*n_sites, invIAT(:,:,k), 3*n_sites, ipiv, info)
        call dgetri(3*n_sites, invIAT(:,:,k), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
        do a = 1, n_sites
          do c3 = 1, 3
            ! call dgemm for G_mat
            VL = 0.d0
            VR = 0.d0
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR(:,:,k), 3*n_sites, dT_LR(:,:,a,c3), & 
                       3*n_sites, 0.d0, VL, 3*n_sites)
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, dA_LR(:,:,a,c3,k), 3*n_sites, T_LR, &
                       3*n_sites, 0.d0, VR, 3*n_sites)
            G_mat(:,:,a,c3,k) = VL + VR
            VL = 0.d0
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, invIAT(:,:,k), 3*n_sites, G_mat(:,:,a,c3,k), &
                       3*n_sites, 0.d0, VL, 3*n_sites)
            ! Take trace of invIAT * G_mat to get force integrand
            do i = 1, 3*n_sites
              force_integrand(a,c3,k) = force_integrand(a,c3,k) + VL(i,i)
            end do
          end do
        end do
      end do

      forces_MBD = 0.d0

      do a = 1, n_sites
        do c3 = 1, 3
          call integrate("trapezoidal", omegas, force_integrand(a,c3,:), 0.d0, 10.d0, integral)
          forces_MBD(a,c3) = 1.d0/(2.d0*pi) * integral
        end do
      end do

      ! Conversion to eV/A:
      forces_MBD = 51.42208619083232 * forces_MBD
  
      write(*,*) "MBD forces:"
      do a = 1, n_sites
        write(*,*) forces_MBD(a,:)
      end do

    else
      E_MBD_k = 0.d0
      forces_MBD_k = 0.d0
      n_tot = 0
      kf = 0 ! Counter for force calculation
      do i = 1, n_sites
        allocate( T_SR_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( B_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( A_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( alpha_SCS_i(1:n_neigh(i),1:11) )
        allocate( A_LR_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( T_LR_i(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( AT_n(1:3*n_neigh(i),1:3,1:11,1:n_order) )
        allocate( AT_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( I_mat_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( series(1:3*n_neigh(i),1:3,1:11) )
        allocate( WR_n(1:3*n_neigh(i)) )
        allocate( WI_n(1:3*n_neigh(i)) )
        allocate( VR_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( VL_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( VRinv_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( logMapo_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( logIAT_n(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( dB_mat_n(1:3*n_neigh(i),1:3*n_neigh(i),1:3,1:11) )
        allocate( dA_mat_n(1:3*n_neigh(i),1:3*n_neigh(i),1:3,1:11) )
        allocate( dBA_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( dalpha_n(1:n_neigh(i),1:3,1:11) )
        allocate( dA_LR_n(1:3*n_neigh(i),1:3*n_neigh(i),1:3,1:11) )
        allocate( dT_LR_n(1:3*n_neigh(i),1:3*n_neigh(i),1:3) )
        allocate( invIAT_n(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( G_mat_n(1:3*n_neigh(i),1:3*n_neigh(i),1:3,1:11) )
        T_SR_i = 0.d0
        B_mat_i = 0.d0
        A_mat_i = 0.d0
        alpha_SCS_i = 0.d0
        A_LR_i = 0.d0
        T_LR_i = 0.d0
        k = 0
        do i2 = 1, n_sites
          if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
            p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
!            write(*,*) "i2, p, neighbors_list(n_tot+p)", i2, p, neighbors_list(n_tot+p)
            k = k+1
            do c1 = 1, 3
              B_mat_i(3*(p-1)+c1,3*(p-1)+c1,:) = 1.d0/alpha_k(n_tot+p,:)
            end do
            do j2 = 2, n_neigh(i2)
              k = k+1
              j = neighbors_list(k)
!              write(*,*) "i2, j, rjs(k)", i2, j, rjs(k)
              if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
                 if ( rjs(k) < rcut) then
                  q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
                  k2 = 9*(k-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k2 = k2+1
                      T_SR_i(3*(p-1)+c1,3*(q-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(k2) * &
                                                        g_func(k,:) + h_func(k2,:))
                    end do
                  end do
                end if
              end if
            end do
          else
            k = k+n_neigh(i2)
          end if
        end do
        B_mat_i = B_mat_i + T_SR_i
        A_mat_i = B_mat_i

!        write(*,*) "B_mat_i:"
!        do p = 1, 6
!          write(*,*) B_mat_i(p,1:6,1)
!        end do


        do k3 = 1, 11
          call dgetrf(3*n_neigh(i), 3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, info)
          call dgetri(3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, work_arr, 12*n_sites, info)
        end do

        do k3 = 1, 11
          do p = 1, n_neigh(i)
            A_i = 0.d0
            do q = 1, n_neigh(i)
              A_i = A_i + A_mat_i(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,k3)
            end do
            alpha_SCS_i(p,k3) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
          end do
        end do

!        write(*,*) "alpha_SCS_i:", alpha_SCS_i(:,1)

        do k3 = 1, 11
          do p = 1, n_neigh(i)
            do c1 = 1, 3
              A_LR_i(3*(p-1)+c1,3*(p-1)+c1,k3) = alpha_SCS_i(p,k3)
            end do
          end do
        end do

        k = 0
        do i2 = 1, n_sites
          if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
            p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
            k = k+1
            r_vdw_i = r0_ii(k) * (alpha_SCS_i(p,1)/neighbor_alpha0(k))**(1.d0/3.d0)
            do j2 = 2, n_neigh(i2)
              k = k+1
              j = neighbors_list(k)
              if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
                if (rjs(k) < rcut) then
                  q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
                  r_vdw_j = r0_ii(k) * (alpha_SCS_i(q,1)/neighbor_alpha0(k))**(1.d0/3.d0)
                  f_damp_SCS_ij = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                  k2 = 9*(k-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k2 = k2+1
                      T_LR_i(3*(p-1)+c1,3*(q-1)+c2) = f_damp_SCS_ij * T_func(k2)
                      T_LR_i(3*(q-1)+c1,3*(p-1)+c2) = T_LR_i(3*(p-1)+c1,3*(q-1)+c2)
                    end do
                  end do
                end if
              end if
            end do
          else
            k = k + n_neigh(i2)
          end if
        end do

!        write(*,*) "T_LR_i:"
!        do p = 1, 6
!          write(*,*) T_LR_i(p,1:6)
!        end do

        AT_i = 0.d0
        do k = 1, 11
          call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, A_LR_i(:,:,k), 3*n_neigh(i), &
                     T_LR_i, 3*n_neigh(i), 0.d0, AT_i(:,:,k), 3*n_neigh(i))
        end do

!        write(*,*) "AT_i:"
!        do p = 1, 6
!          write(*,*) AT_i(p,1:6,1)
!        end do

        I_mat_n = 0.d0
        do i2 = 1, 3*n_neigh(i)
          I_mat_n(i2,i2) = 1.d0
        end do

        logIAT_n = 0.d0
        do k = 1,11
          WR_n = 0.d0
          WI_n = 0.d0
          VL_n = 0.d0
          VR_n = 0.d0
          VRinv_n = 0.d0
          call dgeev('n', 'v', 3*n_neigh(i), I_mat_n-AT_i(:,:,k), 3*n_neigh(i), WR_n, WI_n, VL_n, &
                     3*n_neigh(i), VR_n, 3*n_neigh(i), work_arr, 12*n_sites, info)
          logMapo_n = 0.d0
          do i2 = 1, 3*n_neigh(i)
            logMapo_n(i2,i2) = log(WR_n(i2))
          end do
          VRinv_n = VR_n
          VL_n = 0.d0
          call dgetrf(3*n_neigh(i), 3*n_neigh(i), VRinv_n, 3*n_neigh(i), ipiv, info)
          call dgetri(3*n_neigh(i), VRinv_n, 3*n_neigh(i), ipiv, work_arr, 12*n_sites, info)
          call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, VR_n, 3*n_neigh(i), logMapo_n, &
                     3*n_neigh(i), 0.d0, VL_n, 3*n_neigh(i))
          call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, VL_n, 3*n_neigh(i), VRinv_n, &
                     3*n_neigh(i), 0.d0, logIAT_n(:,:,k), 3*n_neigh(i))
        end do
!       Local energy test:
        integrand = 0.d0
        integrand_k = 0.d0
!        write(*,*) "alpha_SCS_i:", i, alpha_SCS_i(1,1)
        do k = 1,11
          do i2 = 1,3*n_neigh(i)
            integrand(k) = integrand(k) + logIAT_n(i2,i2,k)
          end do
          integrand_k(k) = alpha_SCS_i(1,k)/sum(alpha_SCS_i(:,k)) * integrand(k)
        end do
!        write(*,*) "Integrand:", integrand
        integral = 0.d0
        call integrate("trapezoidal", omegas, integrand_k, 0.d0, 10.d0, integral)
        E_MBD_k(i) = integral/(2.d0*pi)
        E_MBD_k(i) = E_MBD_k(i) * 27.211386245988
        integral = 0.d0
        call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
        E_MBD = integral/(2.d0*pi)
        E_MBD = E_MBD * 27.211386245988
        write(*,*) "E_MBD_k:", i, E_MBD_k(i)
        write(*,*) "E_MBD_tot for k:", i, E_MBD
!        write(*,*) "Total MBD energy:", sum(E_MBD_k(:,1))



!        series = 0.d0
!        AT_n = 0.d0
!        integrand = 0.d0
!        do k = 1, 11
!          call dgemm('n', 'n', 3*n_neigh(i), 3, 3*n_neigh(i), 1.d0, A_LR_i(:,:,k), 3*n_neigh(i), &
!                      T_LR_i(:,1:3), 3*n_neigh(i), 0.d0, AT_n(:,:,k,1), 3*n_neigh(i))
!          series(:,:,k) = -1.d0/2 * AT_n(:,:,k,1)
!          do c1 = 1, 3
!            integrand(k) = integrand(k) + alpha_SCS(1,k) * &
!              dot_product(T_LR_i(c1,:), series(:,c1,k))
!          end do
!        end do
!        integrand = integrand / (2.d0*pi)
!        call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, E_MBD_k(i,1))

!        do k2 = 2, n_order
!          integrand = 0.d0
!          do k = 1, 11
!            call dgemm('n', 'n', 3*n_neigh(i), 3, 3*n_neigh(i), 1.d0, AT_i(:,:,k), 3*n_neigh(i), &
!                       AT_n(:,:,k,k2-1), 3*n_neigh(i), 0.d0, AT_n(:,:,k,k2), 3*n_neigh(i))
!            series(:,:,k) = series(:,:,k) - 1.d0/(k2+1) * AT_n(:,:,k,k2)
!            do c1 = 1, 3
!              integrand(k) = integrand(k) + alpha_SCS(1,k) * &
!                dot_product(T_LR_i(c1,:), series(:,c1,k))
!            end do
!          end do
!          integrand = integrand / (2.d0*pi)
!          call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, E_MBD_k(i,k2))
!        end do

!        write(*,*) "Local energy for site", i, E_MBD_k(i,n_order)*27.211386245988

!        write(*,*) "dB_mat calc:"
        dB_mat_n = 0.d0
        k = 0
        do i2 = 1, n_sites
          if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
            k = k+1
            p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
!            write(*,*) p, i2
            do c1 = 1, 3
              do c3 = 1, 3
                dB_mat_n(3*(p-1)+c1,3*(p-1)+c1,c3,:) = dB_mat(3*(i2-1)+c1,3*(i2-1)+c1,i,c3,:)
              end do
            end do
            do j2 = 2, n_neigh(i2)
              k = k+1
              j = neighbors_list(k)
!              write(*,*) "i2, j, rjs(k)", i2, j, rjs(k)
              if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
                if (rjs(k) < rcut) then
!                  write(*,*) "i2, j:", i2, j
                  q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
!                  if (p == 1 .and. q == 2) then
!                    write(*,*) "p, q, i2, j:", p, q, i2, j, dB_mat(3*(i2-1)+1,3*(j-1)+1,i,1,1)
!                  end if
!                  if (p == 2 .and. q == 1) then
!                    write(*,*) "p, q, i2, j:", p, q, i2, j, dB_mat(3*(i2-1)+1,3*(j-1)+1,i,1,1)
!                  end if
                  do c1 = 1, 3
                    do c2 = 1, 3
                      do c3 = 1, 3
                        dB_mat_n(3*(p-1)+c1,3*(q-1)+c2,c3,:) = dB_mat(3*(i2-1)+c1,3*(j-1)+c2,i,c3,:)
                      end do
                    end do
                  end do
                end if
              end if
            end do
          else
            k = k + n_neigh(i2)
          end if
        end do

!        write(*,*) "dB_mat_n:"
!        do p = 1, 6
!          write(*,*) dB_mat_n(p,1:6,1,1)
!        end do


!        write(*,*) "dB_mat_n:"
!        do p = 1, 3*n_neigh(i)
!          write(*,*) dB_mat_n(p,:,1,1)
!        end do

!        write(*,*) "dB_mat_n:"
!        do p = 1, 6
!          write(*,*) dB_mat_n(p,1:6,1,1)
!        end do

!        write(*,*) "dB_mat_n:"
!        write(*,*) dB_mat_n(3*(1-1)+1,1:6,1,1)
!        write(*,*) dB_mat_n(3*(1-1)+2,1:6,1,1)
!        write(*,*) dB_mat_n(3*(1-1)+3,1:6,1,1)
!        write(*,*) dB_mat_n(3*(42-1)+1,3*(42-1)+1,1,1)
!        write(*,*) dB_mat_n(3*(42-1)+2,3*(42-1)+2,1,1)
!        write(*,*) dB_mat_n(3*(42-1)+3,3*(42-1)+3,1,1)

        dA_mat_n = 0.d0
        do k = 1, 11
          do c3 = 1, 3
            dBA_n = 0.d0
            call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, dB_mat_n(:,:,c3,k), 3*n_neigh(i), &
                       A_mat_i(:,:,k), 3*n_neigh(i), 0.d0, dBA_n, 3*n_neigh(i))
            call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, -A_mat_i(:,:,k), 3*n_neigh(i), dBA_n, &
                       3*n_neigh(i), 0.d0, dA_mat_n(:,:,c3,k), 3*n_neigh(i))
          end do
        end do

!        write(*,*) "dA_mat_n:"
!        write(*,*) dA_mat_n(3*(1-1)+1,1:6,1,1)
!        write(*,*) dA_mat_n(3*(1-1)+2,1:6,1,1)
!        write(*,*) dA_mat_n(3*(1-1)+3,1:6,1,1)

        dalpha_n = 0.d0
        do k = 1, 11
          do p = 1, n_neigh(i)
            do c3 = 1, 3
              A_i = 0.d0
              do q = 1, n_neigh(i)
                A_i = A_i + dA_mat_n(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,c3,k)
              end do
              do c1 = 1, 3
                dalpha_n(p,c3,k) = dalpha_n(p,c3,k) + A_i(c1,c1)
              end do
              dalpha_n(p,c3,k) = 1.d0/3.d0 * dalpha_n(p,c3,k)
            end do
          end do
        end do

!        write(*,*) "dalpha_n:", dalpha_n(:,1,1)
!        write(*,*) "dalpha_n(42):", dalpha_n(42,1,1)

        dA_LR_n = 0.d0

        do k = 1, 11
          do c3 = 1, 3
            do p = 1, n_neigh(i)
              do c1 = 1, 3
                dA_LR_n(3*(p-1)+c1,3*(p-1)+c1,c3,k) = dalpha_n(p,c3,k)
              end do
            end do
          end do
        end do

!        do i2 = 1, 6
!          write(*,*) "dA_LR_n:", dA_LR_n(i2,1:6,1,1)
!        end do

        dT_LR_n = 0.d0
!        f_damp_der_n = 0.d0
!        f_damp_der_SCS_n = 0.d0
        k = 0
!        write(*,*) "dT_LR_n calc:"
        do i2 = 1, n_sites
          if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
            p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
            k = k+1
            r_vdw_i = r0_ii(k) * (alpha_SCS_i(p,1)/neighbor_alpha0(k))**(1.d0/3.d0)
            do j2 = 2, n_neigh(i2)
              k = k+1
              j = neighbors_list(k)
!              write(*,*) "i2, j, rjs(k)", i2, j, rjs(k)
              if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
                if (rjs(k) < rcut) then
                  q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
                  r_vdw_j = r0_ii(k) * (alpha_SCS_i(q,1)/neighbor_alpha0(k))**(1.d0/3.d0)
                  R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                  S_vdW_ij = sR*R_vdW_SCS_ij
                  f_damp_SCS_ij = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                  do c3 = 1, 3
                    dS_vdW_ij = sR/3.d0 * ( r_vdw_i/alpha_SCS_i(p,1) * dalpha_n(p,c3,1) + &
                                            r_vdw_j/alpha_SCS_i(q,1) * dalpha_n(q,c3,1) )
                    f_damp_der_SCS_ij_n(c3) = -(d*rjs_H(k))/S_vdW_ij**2 * f_damp_SCS_ij**2 * &
                                               exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                    k2 = 9*(k-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k2 = k2+1
                        dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) = dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) + &
                                      T_func(k2) * f_damp_der_SCS_ij_n(c3)
                      end do
                    end do
                  end do
                  if (i == i2 .or. i == j) then
                    do c3 = 1, 3
                    f_damp_der_ij_n(c3) = d/S_vdW_ij * f_damp_SCS_ij**2 * &
                                         exp( -d*(rjs_H(k)/S_vdW_ij - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
                      k2 = 9*(k-1)
                      do c1 = 1, 3
                        do c2 = 1, 3
                          k2 = k2+1
                          if (i == i2) then
                            dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) = dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) - &
                                               T_func(k2) * f_damp_der_ij_n(c3) - &
                                               dT(k2,c3) * f_damp_SCS_ij
                          else if (i == j) then
                            dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) = dT_LR_n(3*(p-1)+c1,3*(q-1)+c2,c3) + &
                                               T_func(k2) * f_damp_der_ij_n(c3) + &
                                               dT(k2,c3) * f_damp_SCS_ij
                          end if
                        end do
                      end do
                    end do
                  end if                  
                end if
              end if
            end do
          else
            k = k + n_neigh(i2)
          end if
        end do

!        write(*,*) "dT_LR_n:"
!        do p = 1, 6
!          write(*,*) dT_LR_n(p,1:6,1)
!        end do

!        write(*,*) "dT_LR_n:"
!        write(*,*) dT_LR_n(3*(1-1)+1,1:6,1)
!        write(*,*) dT_LR_n(3*(1-1)+2,1:6,1)
!        write(*,*) dT_LR_n(3*(1-1)+3,1:6,1)

        invIAT_n = 0.d0
        G_mat_n = 0.d0
        force_integrand_n = 0.d0

        do k = 1, 11
          invIAT_n(:,:,k) = I_mat_n - AT_i(:,:,k)
          call dgetrf(3*n_neigh(i), 3*n_neigh(i), invIAT_n(:,:,k), 3*n_neigh(i), ipiv, info)
          call dgetri(3*n_neigh(i), invIAT_n(:,:,k), 3*n_neigh(i), ipiv, work_arr, 12*n_sites, info)
          do c3 = 1, 3
            ! call dgemm for G_mat
            VL_n = 0.d0
            VR_n = 0.d0
            call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, A_LR_i(:,:,k), 3*n_neigh(i), &
                       dT_LR_n(:,:,c3), 3*n_neigh(i), 0.d0, VL_n, 3*n_neigh(i))
            call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, dA_LR_n(:,:,c3,k), 3*n_neigh(i), &
                       T_LR_i, 3*n_neigh(i), 0.d0, VR_n, 3*n_neigh(i))
            G_mat_n(:,:,c3,k) = VL_n + VR_n
            VL_n = 0.d0
            call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, invIAT_n(:,:,k), 3*n_neigh(i), &
                       G_mat_n(:,:,c3,k), 3*n_neigh(i), 0.d0, VL_n, 3*n_neigh(i))
            ! Take trace of invIAT * G_mat to get force integrand
            do i2 = 1, 3*n_neigh(i)
              force_integrand_n(c3,k) = force_integrand_n(c3,k) + VL_n(i2,i2)
            end do
          end do
        end do

!        forces_MBD = 0.d0

        integral = 0.d0
        do c3 = 1, 3
          call integrate("trapezoidal", omegas, force_integrand_n(c3,:), 0.d0, 10.d0, integral)
          forces_MBD_k(i,c3) = 1.d0/(2.d0*pi) * integral
        end do
        forces_MBD_k(i,:) = forces_MBD_k(i,:) * 51.42208619083232
        write(*,*) "force_k:", i, forces_MBD_k(i,:)

        deallocate( T_SR_i, B_mat_i, A_mat_i, alpha_SCS_i, A_LR_i, T_LR_i, AT_n, AT_i, series, I_mat_n, &
                    WR_n, WI_n, VR_n, VL_n, VRinv_n, logMapo_n, logIAT_n, dB_mat_n, dA_mat_n, dBA_n, dalpha_n, &
                    dA_LR_n, dT_LR_n, invIAT_n, G_mat_n )
        n_tot = n_tot + n_neigh(i)
      end do
!      write(*,*) "n_order | Total MBD energy:"
!      do k2 = 1, n_order
      write(*,*) "Total MBD energy:", sum(E_MBD_k)
!      end do
    end if

!      A_LR_k = 0.d0
!      do i = 1, n_sites
!        A_LR_k(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3,:,i) = A_LR(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3,:)
!      end do
!      AT_n = 0.d0
!      AT_n(:,:,:,1) = AT
!      logIAT = 0.d0
!      do k = 1, 11
!        logIAT(:,:,k) = -I_mat
!      end do
!      do k2 = 1, n_order-1
!        integrand_k = 0.d0
!        do k = 1, 11
!          call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, AT(:,:,k), 3*n_sites, AT_n(:,:,k,k2), &
!                     3*n_sites, 0.d0, AT_n(:,:,k,k2+1), 3*n_sites)
!          logIAT(:,:,k) = logIAT(:,:,k) - 1.d0/(k2+1) * AT_n(:,:,k,k2)
!          do i =1, n_sites
!            VL = 0.d0
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, T_LR, 3*n_sites, logIAT(:,:,k), &
!                       3*n_sites, 0.d0, VL, 3*n_sites)
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR_k(:,:,k,i), 3*n_sites, VL, &
!                       3*n_sites, 0.d0, AT_k(:,:,k,i), 3*n_sites)
!            do i2 = 1, 3*n_sites
!              integrand_k(k,i) = integrand_k(k,i) + AT_k(i2,i2,k,i)
!            end do
!            do c1 = 1, 3
!              integrand_k(k,i) = integrand_k(k,i) + alpha_SCS(i,k) * &
!                dot_product(T_LR(3*(i-1)+c1,:), logIAT(:,3*(i-1)+c1,k))
!            end do
!          end do
!          do i = 1, 3*n_sites
!            integrand(k) = integrand(k) + logIAT(i,i,k)
!          end do 
!        end do
!        write(*,*) "integrand_k", integrand_k(:,1)
!        E_MBD_k = 0.d0
!        do i = 1, n_sites
!          call integrate("trapezoidal", omegas, integrand_k(:,i), 0.d0, 10.d0, E_MBD_k(i))
!        end do
!        E_MBD_k = E_MBD_k/(2.d0*pi)
!        write(*,*) "Local energies for order", k2+1, E_MBD_k*27.211386245988
!        write(*,*) "Sum of local energies", sum(E_MBD_k)*27.211386245988
!        call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
!        write(*,*) "n_order, E_MBD:", k2+1, integral/(2.d0*pi) * 27.211386245988
!      end do
!      integrand = 0.d0
!      do k = 1,11
!        do i = 1,3*n_sites
!          integrand(k) = integrand(k) + logIAT(i,i,k)
!        end do
!      end do

!      call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
!    end if

!   Force calculation starts here:
!    dT = 0.d0
!    k = 0
!    do i = 1, n_sites
!      k = k+1
!      do j2 = 2, n_neigh(i)
!        k = k+1
!        if (rjs(k) < rcut_vdw) then
!          k2 = 9*(k-1)
!          do c1 = 1,3
!            do c2 = 1,3
!              k2 = k2+1
!              do c3 = 1,3
!                dT(k2,c3) = (-15.d0 * xyz_H(c1,k) * xyz_H(c2,k) * xyz_H(c3,k))/rjs_H(k)**7
!                if (c1 == c2) then
!                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c3,k)
!                end if
!                if (c2 == c3) then
!                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c1,k)
!                end if
!                if (c1 == c3) then
!                  dT(k2,c3) = dT(k2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c2,k)
!                end if
!              end do
!            end do
!          end do
!        end if
!      end do
!    end do
!
!    f_damp_der = 0.d0
!    g_func_der = 0.d0
!    h_func_der = 0.d0
!    dT_SR = 0.d0
!
!    do a = 1, n_sites
!      k = 0
!      do i = 1, n_sites
!        k = k+1
!        do j2 = 2, n_neigh(i)
!          k = k+1
!          j = neighbors_list(k)
!          if (rjs(k) < rcut_vdw) then
!            if (a == i .or. a == j) then
!              sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
!              do c3 = 1, 3
!                f_damp_der(k,c3) = d/(sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))) * f_damp(k)**2 * &
!                                     exp( -d*(rjs_H(k)/(sR*(r0_ii(n_sites*(i-1)+1) + &
!                                     r0_ii(k))) - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
!                g_func_der(k,c3,:) = 4.d0/sqrt(pi) * rjs_H(k)/sigma_ij**3 * xyz_H(c3,k) * &
!                                       exp(-rjs_H(k)**2/sigma_ij**2) 
!                k2 = 9*(k-1)
!                do c1 = 1, 3
!                  do c2 = 1, 3
!                    k2 = k2+1
!                    coeff_h_der = 4.d0/sqrt(pi) * 1.d0/sigma_ij**3 * exp(-rjs_H(k)**2/sigma_ij**2)
!                    terms = 0.d0
!                    if (c1 == c3) then
!                      terms = terms + xyz_H(c2,k)/rjs_H(k)**2
!                    end if
!                    if (c2 == c3) then
!                      terms = terms + xyz_H(c1,k)/rjs_H(k)**2
!                    end if
!                    terms = terms + -2.d0*(xyz_H(c1,k)*xyz_H(c2,k)*xyz_H(c3,k))/rjs_H(k)**2 * &
!                            (1.d0/rjs_H(k)**2 + 1.d0/sigma_ij**2)
!                    h_func_der(k2,c3,:) = coeff_h_der * terms
!                    dT_SR(3*(i-1)+c1,3*(j-1)+c2,a,c3,:) = f_damp_der(k,c3) * T_func(k2) * &
!                                                          g_func(k,:) - (1.d0 - f_damp(k)) * dT(k2,c3) * &
!                                                          g_func(k,:) - g_func_der(k,c3,:) * (1.d0 - f_damp(k)) * &
!                                                          T_func(k2) - f_damp_der(k,c3) * h_func(k2,:) + &
!                                                          h_func_der(k2,c3,:) * (1.d0 - f_damp(k)) 
!                  end do
!                end do
!              end do
!              if (a == i) then
!                dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:) = -dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:)
!              end if
!            end if
!          end if
!        end do
!      end do
!    end do
!
!    dB_mat = dT_SR
!    if (do_hirshfeld_gradients) then
!      dB_mat_v = 0.d0
!      dT_SR_v = 0.d0
!      coeff_fdamp = 0.d0
!      coeff_der = 0.d0
!      k = 0
!      do i = 1, n_sites
!        k = k+1
!        do j2 = 2, n_neigh(i)
!          k = k+1
!          j = neighbors_list(k)
!          if (rjs(k) < rcut_vdw) then
!            sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
!            S_vdW_ij = sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))
!            exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
!            k2 = 9*(k-1)
!            do c1 = 1, 3
!              do c2 = 1, 3
!                k2 = k2+1
!                coeff_fdamp(i,j,c1,c2,:) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
!                                           (-T_func(k2) * g_func(k,:) + h_func(k2,:))
!                dg = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * (rjs_H(k)/sigma_ij)**2 * &
!                     (-rjs_H(k)/(3.d0 * sigma_ij**3))
!                dh = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * xyz_H(c1,k)*xyz_H(c2,k) / &
!                     (sigma_ij**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij)**2)
!                coeff_der(i,j,c1,c2,:) = (1.d0-f_damp(k)) * (-T_func(k2)*dg+dh)
!              end do
!            end do
!          end if
!        end do
!      end do
!
!      k2 = 0
!      k3 = 0
!      do i = 1, n_sites
!        do a2 = 1, n_neigh(i)
!          k2 = k2+1
!          a = neighbors_list(k2)
!          do c3 = 1, 3
!            do c1 = 1, 3
!              dB_mat_v(3*(i-1)+c1,3*(i-1)+c1,a,c3,:) = -1.d0/(hirshfeld_v(i)*alpha_i(i,:)) * &
!                            hirshfeld_v_cart_der_H(c3,k2)
!            end do
!            do k = 1, 11
!              do j2 = 2, n_neigh(i)
!                if (rjs(k3+j2) < rcut_vdw) then
!                  j = neighbors_list(k3+j2)
!                  do c1 = 1, 3
!                    do c2 = 1, 3
!                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                        coeff_der(i,j,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
!                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                        coeff_fdamp(i,j,c1,c2,k) * r0_ii(n_sites*(i-1)+1)/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
!                    end do
!                  end do
!                end if      
!              end do
!            end do
!          end do
!        end do
!        k3 = k3 + n_neigh(i)
!      end do
!      k2 = 0
!      k3 = 0
!      do j = 1, n_sites
!        do a2 = 1, n_neigh(j)
!          k2 = k2+1
!          a = neighbors_list(k2)
!          do c3 = 1, 3
!            do k = 1, 11
!              do i2 = 2, n_neigh(j)
!                if (rjs(k3+i2) < rcut_vdw) then
!                  i = neighbors_list(k3+i2)
!                  do c1 = 1, 3
!                    do c2 = 1, 3
!                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                        coeff_der(i,j,c1,c2,k) * (sigma_i(j,k)**2/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2))
!                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                        coeff_fdamp(i,j,c1,c2,k) * r0_ii(n_sites*(j-1)+1)/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2)
!                    end do
!                  end do
!                end if
!              end do
!            end do
!          end do
!        end do
!        k3 = k3 + n_neigh(j)
!      end do
!
!      dB_mat_v = dB_mat_v + dT_SR_v
!
!      dB_mat = dB_mat + dB_mat_v
!    end if


    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing forces:", time2-time1
    end if
 
!   Clean up
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, A_LR, r0_ii_SCS, f_damp_SCS, I_mat, AT, logIAT, logMapo, VR, VRinv, WR, WI, VL, integrand, &
                dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, dalpha, dA_LR, dT_LR, f_damp_der_SCS, &
                invIAT, G_mat, force_integrand, forces_MBD, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                dB_mat, dB_mat_v, dv_i, dv_j, &
                coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H, A_LR_k, alpha_k, sigma_k, s_i, s_j, &
                integrand_2, E_MBD_k, f_damp_der_ij_n, f_damp_der_SCS_ij_n, force_integrand_n, forces_MBD_k, integrand_k )

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
