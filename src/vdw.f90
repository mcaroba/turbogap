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
                          alpha0_ref(:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_forces
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
!   TEMPORARILY ADDED hirshfeld_v AND hirshfeld_v_neigh AS INOUT VARIABLES FOR DEBUGGIN!!! (originally in vars)
    real*8, intent(inout) :: energies(:), forces0(:,:), hirshfeld_v(:), hirshfeld_v_neigh(:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), neighbor_c6_ij(:), r0_ii(:), r0_ij(:), &
                           exp_damp(:), f_damp(:), c6_ij_free(:), neighbor_alpha0(:), &
                           pref_force1(:), pref_force2(:), r6(:), r6_der(:), &
                           T_func(:,:), h_func(:,:,:,:,:), g_func(:,:,:), &
                           omegas(:), omega_i(:), &
                           alpha_i(:,:), sigma_i(:,:), sigma_ij(:), T_SR(:,:,:), B_mat(:,:,:), &
                           rjs_H(:), xyz_H(:,:), A_mat(:,:,:), work_arr(:), A_i(:,:), alpha_SCS(:,:), &
                           A_LR(:,:,:), T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), I_mat(:,:), AT(:,:,:), &
                           logIAT(:,:,:), VR(:,:), logMapo(:,:), VRinv(:,:), WR(:), WI(:), VL(:,:), &
                           integrand(:), v_arr(:), dT(:,:,:), dT_SR(:,:,:,:,:), f_damp_der(:,:,:), &
                           g_func_der(:,:,:,:), h_func_der(:,:,:,:,:,:), dA_mat(:,:,:,:,:), &
                           dalpha(:,:,:,:), dA_LR(:,:,:,:,:), dT_LR(:,:,:,:), f_damp_der_SCS(:,:,:,:), &
                           invIAT(:,:,:), G_mat(:,:,:,:,:), force_integrand(:,:,:), forces_MBD(:,:), &
                           coeff_h_der(:), terms(:), dT_SR_A_mat(:,:), dT_SR_v(:,:,:,:,:), &
                           f_damp_der_v(:,:,:,:), g_func_der_v(:,:,:,:,:), h_func_der_v(:,:,:,:,:,:,:), &
                           dB_mat(:,:,:,:,:), dB_mat_v(:,:,:,:,:), dA_mat_v(:,:,:,:,:), dalpha_v(:,:,:,:), &
                           dalpha_v_r(:,:,:,:), dR_vdW(:), dv_i(:), dv_j(:), dsigma(:,:), &
                           g_func_der_v_coeff(:), coeff_der(:,:,:,:,:), coeff_fdamp(:,:,:,:,:), dg(:), &
                           dh(:), hirshfeld_v_cart_der_H(:,:)
    real*8 :: time1, time2, c6_ii, c6_jj, r0_i, r0_j, alpha0_i, alpha0_j, rbuf, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              f_damp_der_v_coeff
    integer, allocatable:: i_buffer(:), ipiv(:)
    integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0, info
    integer :: i, i2, j, j2, k, k2, k3, a, n_in_buffer, c1, c2, c3, lwork, b, a2
    logical, allocatable :: is_in_buffer(:)
    logical :: do_timing = .false., do_hirshfeld_gradients = .true.

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!    write(*,*) "hirshfeld_v:", hirshfeld_v
!    write(*,*) "hirshfeld_v_neigh:", hirshfeld_v_neigh

!   Debugging for dimer (remove later):
!   plus:
!    hirshfeld_v = (/ 1.0154463183569753, 1.0154463183569753 /)
!    hirshfeld_v_neigh = (/ 1.0154463183569753, 1.0154463183569753, &
!                           1.0154463183569753, 1.0154463183569753 /)
!   minus:
!    hirshfeld_v = (/ 1.0168126496003040, 1.0168126496003040 /)
!    hirshfeld_v_neigh = (/ 1.0168126496003040, 1.0168126496003040, &
!                           1.0168126496003040, 1.0168126496003040 /)


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
    allocate( v_arr(1:n_sites) )
    allocate( dT(1:3*n_sites,1:3*n_sites,1:3) )
    allocate( dT_SR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( f_damp_der(1:n_sites,1:n_sites,1:3) )
    allocate( g_func_der(1:n_sites,1:n_sites,1:3,1:11) )
    allocate( h_func_der(1:n_sites,1:n_sites,1:3,1:3,1:3,1:11) )
    allocate( dA_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dalpha(1:n_sites,1:n_sites,1:3,1:11) )
    allocate( dA_LR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dT_LR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3) )
    allocate( f_damp_der_SCS(1:n_sites,1:n_sites,1:n_sites,1:3) )
    allocate( invIAT(1:3*n_sites,1:3*n_sites,1:11) )
    allocate( G_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( force_integrand(1:n_sites,1:3,1:11) )
    allocate( forces_MBD(1:n_sites,1:3) )
    allocate( coeff_h_der(1:11) )
    allocate( terms(1:11) )
    allocate( dT_SR_A_mat(1:3*n_sites,1:3*n_sites) )
    allocate( dT_SR_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( f_damp_der_v(1:n_sites,1:n_sites,1:n_sites,1:3) )
    allocate( g_func_der_v(1:n_sites,1:n_sites,1:n_sites,1:3,1:11) )
    allocate( h_func_der_v(1:n_sites,1:n_sites,1:n_sites,1:3,1:3,1:3,1:11) )
    allocate( dB_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dB_mat_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dalpha_v(1:n_sites,1:n_sites,1:3,1:11) )
    allocate( dA_mat_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:11) )
    allocate( dR_vdW(1:3) )
    allocate( dv_i(1:3) )
    allocate( dv_j(1:3) )
    allocate( dsigma(1:3,1:11) )
    allocate( g_func_der_v_coeff(1:11) )
    allocate( coeff_der(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( coeff_fdamp(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( dg(1:11) )
    allocate( dh(1:11) )
    allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
!    allocate( dalpha_v_r(1:n_sites,1:n_sites,1:3,1:11) )
    is_in_buffer = .false.
    if( do_forces )then
      allocate( pref_force1(1:n_sites) )
      allocate( pref_force2(1:n_sites) )
      allocate( r6_der(1:n_pairs) )
    end if

!   Hirshfeld derivative test:
!    k = 0
!    write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh
!    do i = 1, n_neigh(1)
!      k = k+1
!      j = neighbors_list(k)
!      write(*,*) i, j, hirshfeld_v_neigh(k)
!    end do
!    k = 0
!    write(*,*) "hirshfeld_v_cart_der"
!    do i = 1, n_sites
!      do a2 = 1, n_neigh(i)
!        k = k+1
!        a = neighbors_list(k)
!        write(*,*) i, a, hirshfeld_v_cart_der(1:3,k)
!      end do
!    end do

!   Temporary fixed array for Hirshfeld volumes:
    v_arr = (/ 0.85308339, &
               0.85305028, &
               0.85307521, &
               0.85306433, &
               0.85309976, &
               0.85300157, &
               0.85300707, &
               0.85305648, &
               0.85298019, &
               0.85303452, &
               0.85295826, &
               0.85306464, &
               0.85306557, &
               0.85305240, &
               0.85299413, &
               0.85309070, &
               0.85299418, &
               0.85310156, &
               0.85296718, &
               0.85312167, &
               0.85317982, &
               0.85305540, &
               0.85310142, &
               0.85317319, &
               0.85309425, &
               0.85311298, &
               0.85310860, &
               0.85305652, &
               0.85300463, &
               0.85308386, &
               0.85306393, &
               0.85300010, &
               0.85312190, &
               0.85306154, &
               0.85315732, &
               0.85299361, &
               0.85298426, &
               0.85317383, &
               0.85318570, &
               0.85305124, &
               0.85308609, &
               0.85306547, &
               0.85303523, &
               0.85306989, &
               0.85312138, &
               0.85306431, &
               0.85304930, &
               0.85315359, &
               0.85309796, &
               0.85300060, &
               0.85306619, &
               0.85310520, &
               0.85300061, &
               0.85306354, &
               0.85312192, &
               0.85296821, &
               0.85298307, &
               0.85312609, &
               0.85310662, &
               0.85295523 /)
!    do k = 1, n_sites
!      write(*,*) "v_arr:", v_arr(k)
!    end do

!    write(*,*) "v_arr:", v_arr(1)
!   THIS IS FOR TESTING:
!    v_arr(1) = v_arr(1) - 0.01d0
!   REMOVE AFTERWARDS
!    write(*,*) "v_arr:", v_arr(1)

    if( do_timing) then
      call cpu_time(time1)
    end if
    
!   Frequencies used for integration:
    omega = 0.d0
    do i = 1, 11
      omegas(i) = omega
      omega = omega + 0.4d0 
    end do
!    write(*,*) "omegas", omegas

!   Check which atoms are in the buffer region
    do k = 1, n_pairs
      j = neighbor_species(k)
      neighbor_c6_ii(k) = c6_ref(j) / (Hartree*Bohr**6)
      r0_ii(k) = r0_ref(j) / Bohr
      neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
      xyz_H(:,k) = xyz(:,k)/Bohr
      rjs_H(k) = rjs(k)/Bohr
      hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,k)*Bohr
    end do
!    write(*,*) "r0_ii:", r0_ii(3090), r0_ii(3564)

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

!    if( do_timing) then
!      call cpu_time(time2)
!      write(*,*) "vdw: creating pair quantities 1:", time2-time1
!      call cpu_time(time1)
!    end if

!   UNCOMMENT THIS
!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    omega_i = (4.d0 * neighbor_c6_ii)/(3*neighbor_alpha0**2) 
!   UNCOMMENT END

!    write(*,*) "i, R_vdW:", 1, r0_ii(1)
!    write(*,*) "i, R_vdW:", 2, r0_ii(4)
!    write(*,*) "i, a, dR_vdW:", 1, 1, 1.d0/3.d0 * r0_ii(1)/hirshfeld_v(1) * hirshfeld_v_cart_der(1:3,1)
!    write(*,*) "i, a, dR_vdW:", 2, 1, 1.d0/3.d0 * r0_ii(4)/hirshfeld_v(2) * hirshfeld_v_cart_der(1:3,4)

!   Temporary volume assignment:
!    do k = 1, n_pairs
!      j = neighbors_list(k)
!      neighbor_c6_ii(k) = neighbor_c6_ii(k) * v_arr(j)**2
!      r0_ii(k) = r0_ii(k) * v_arr(j)**(1.d0/3.d0)
!      neighbor_alpha0(k) = neighbor_alpha0(k) * v_arr(j)
!    end do
!    write(*,*) "r0_ii:", r0_ii(3090), r0_ii(3564)
!   Temporary assignment ends here.
    omega_i = (4.d0 * neighbor_c6_ii)/(3.d0*neighbor_alpha0**2)

    
!    if( do_timing) then
!      call cpu_time(time2)
!      write(*,*) "vdw: creating pair quantities 2:", time2-time1
!      call cpu_time(time1)
!    end if

    do i = 1, n_sites
      k = (i-1)*n_sites+1
      do j = 1, 11
        alpha_i(i,j) = neighbor_alpha0(k)/(1.d0 + omegas(j)**2/omega_i(k)**2)
        sigma_i(i,j) = (sqrt(2.d0/pi) * alpha_i(i,j)/3.d0)**(1.d0/3.d0)
      end do
    end do

!    write(*,*) "sigma_21:", sqrt(sigma_i(2,1)**2 + sigma_i(1,1)**2)

!    write(*,*) "alpha_i:", alpha_i(:,1)
!    write(*,*) "sigma_i:", sigma_i(:,1)

!   Computing dipole interaction tensor
!   Requires the complete supercell to get correct dimensions for T_func! 3*N_at x 3*N_at, where N_at are atoms in supercell
    T_func = 0.d0
    f_damp = 0.d0
    k = 0
!    write(*,*) "f_damp calc:"
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k) ! NOTE: mod not necessary because we are using only single C60 for now
!        if (i == j) then
!          write(*,*) i, j
!        end if
        if( rjs(k) < rcut )then
          f_damp(k) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))) - 1.d0 ) ) )
!          if ( i == 2 .and. j == 1) then
!            write(*,*) "i, j, f_damp:", i, j, f_damp(k)
!          end if
!          if (i < 7 .and. j < 7) then
!            write(*,*) "f_damp"
!            write(*,*) i, j
!            write(*,*) f_damp(k)
!          end if
!          write(*,*) i, j, f_damp(k)
          do c1 = 1, 3
            do c2 = 1, 3
              if (c1 == c2) then
                T_func(3*(i-1)+c1,3*(j-1)+c1) = (3*xyz_H(c1,k) * xyz_H(c1,k) - rjs_H(k)**2)/rjs_H(k)**5
              else
                T_func(3*(i-1)+c1,3*(j-1)+c2) = (3*xyz_H(c1,k) * xyz_H(c2,k))/rjs_H(k)**5
!                if (3*(i-1)+c1 == 3*(j-1)+c2) then
!                  write(*,*) i,j,c1,c2
!                end if
              end if
!            T_func(3*(i-1)+(c1+1),3*(j-1)+(c1+1)) = (3*xyz_H(c1+1,k) * xyz_H(c1+1,k) - rjs_H(k)**2)/rjs_H(k)**5
!            T_func(3*(j-1)+(c1+1),3*(i-1)+(c1+1)) = T_func(3*(i-1)+(c1+1),3*(j-1)+(c1+1))
!            do c2 = c1, 2
!              T_func(3*(i-1)+(c1+1),3*(j-1)+(c2+1)) = (3*xyz_H(c1+1,k) * xyz_H(c2+1,k))/rjs_H(k)**5
!              T_func(3*(j-1)+c1,3*(i-1)+c2) = T_func(3*(i-1)+c1,3*(j-1)+c2)
!              T_func(3*(i-1)+(c2+1),3*(j-1)+(c1+1)) = T_func(3*(i-1)+(c1+1),3*(j-1)+(c2+1))
!              T_func(3*(j-1)+c2,3*(i-1)+c1) = T_func(3*(i-1)+c1,3*(j-1)+c2)
            end do
          end do
        end if
      end do
    end do

!    write(*,*) "f_damp:", f_damp
!    write(*,*) "Size of neighbor list:", size(neighbors_list)
!    write(*,*) "Size:", size(T_func, 1), size(T_func, 2)
!    write(*,*) "T_ij:", T_func(1,1:9)
!    write(*,*) "T_ij:", T_func(2,1:9)
!    write(*,*) "T_ij:", T_func(3,1:9)
!    write(*,*) "T_ij:", T_func(4,1:9)
!    write(*,*) "T_ij:", T_func(5,1:9)
!    write(*,*) "T_ij:", T_func(6,1:9)
!    write(*,*) "T_ij:", T_func(7,1:9)
!    write(*,*) "T_ij:", T_func(8,1:9)
!    write(*,*) "T_ij:", T_func(9,1:9)
!    write(*,*) "T_ij end"
!    write(*,*) "T_ij:", T_func(172,172:180)
!    write(*,*) "T_ij:", T_func(173,172:180)
!    write(*,*) "T_ij:", T_func(174,172:180)
!    write(*,*) "T_ij:", T_func(175,172:180)
!    write(*,*) "T_ij:", T_func(176,172:180)
!    write(*,*) "T_ij:", T_func(177,172:180)
!    write(*,*) "T_ij:", T_func(178,172:180)
!    write(*,*) "T_ij:", T_func(179,172:180)
!    write(*,*) "T_ij:", T_func(180,172:180)
!    write(*,*) "T_ij full:", T_func


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
!          g_func(j,i,:) = g_func(i,j,:)
!          write(*,*) i, j, g_func(i,j,1)
          do c1 = 1, 3
            do c2 = 1, 3
              h_func(i,j,c1,c2,:) = 4.d0/sqrt(pi) * (rjs_H(k)/sigma_ij)**3 * &
                                    xyz_H(c1,k)*xyz_H(c2,k)/rjs_H(k)**5 * exp(-rjs_H(k)**2/sigma_ij**2)
!              h_func(i,j,c2,c1,:) = h_func(i,j,c1,c2,:)
!              h_func(j,i,c1,c2,:) = h_func(i,j,c1,c2,:)
!              h_func(j,i,c2,c1,:) = h_func(i,j,c1,c2,:)
            end do
          end do 
!          write(*,*) i,j,h_func(i,j,1,3,1)
        end if
      end do
    end do

!    write(*,*) "g_ij:", g_func(2,1,1)
!    write(*,*) "h_ij:", h_func(2,1,1,1,1)
!    write(*,*) "-T_ij * g(...) + h(...)", -T_func(4,1) * g_func(2,1,1) + h_func(2,1,1,1,1)
!    write(*,*) "T_func:"
!    do i = 1, 6
!      write(*,*) T_func(i,1:6)
!    end do

!    write(*,*) c6_ref(1)/(Hartree*Bohr**6), alpha0_ref(1)/Bohr**3
!    write(*,*) "g_func:", g_func(1,1:6,1)
!    write(*,*) "g_func:", g_func(2,1:6,1)
!    write(*,*) "g_func:", g_func(3,1:6,1)
!    write(*,*) "g_func:", g_func(4,1:6,1)
!    write(*,*) "g_func:", g_func(5,1:6,1)
!    write(*,*) "g_func:", g_func(6,1:6,1)
!    write(*,*) "h_func:", h_func(1,1:6,1,1,1)
!    write(*,*) "h_func:", h_func(2,1:6,1,1,1)
!    write(*,*) "h_func:", h_func(3,1:6,1,1,1)
!    write(*,*) "h_func:", h_func(4,1:6,1,1,1)
!    write(*,*) "h_func:", h_func(5,1:6,1,1,1)
!    write(*,*) "h_func:", h_func(6,1:6,1,1,1)
!    write(*,*) "h_func max:", maxval(h_func)

!    write(*,*) "g_func:"
!    do i = 1,6
!      write(*,*) g_func(i,1:6,1)
!    end do
!    write(*,*) "h_func:"
!    do i = 1,6
!      write(*,*) h_func(i,1:6,1,1,1)
!    end do


    T_SR = 0.d0
    B_mat = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do c1 = 1, 3
        B_mat(3*(i-1)+c1,3*(i-1)+c1,:) = 1.d0/alpha_i(i,:)
      end do
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        do c1 = 1, 3
          do c2 = 1, 3
            T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(3*(i-1)+c1,3*(j-1)+c2) * &
                                          g_func(i,j,:) + h_func(i,j,c1,c2,:))
!            if (3*(i-1)+c1 == 178 .and. 3*(j-1)+c2 == 154) then
!              write(*,*) "T_SR value:", T_SR(3*(i-1)+c1,3*(j-1)+c2,1)
!              write(*,*) f_damp(k)
!              write(*,*) T_func(3*(i-1)+c1,3*(j-1)+c2)
!              write(*,*) g_func(i,j,1)
!              write(*,*) h_func(i,j,c1,c2,1)
!              write(*,*) "k = ", k
!            end if
!            if (3*(i-1)+c1 == 154 .and. 3*(j-1)+c2 == 178) then
!              write(*,*) "T_SR value:", T_SR(3*(i-1)+c1,3*(j-1)+c2,1)
!              write(*,*) f_damp(k)
!              write(*,*) T_func(3*(i-1)+c1,3*(j-1)+c2)
!              write(*,*) g_func(i,j,1)
!              write(*,*) h_func(i,j,c1,c2,1)
!              write(*,*) "k = ", k
!            end if
 !           T_SR(3*(j-1)+c1,3*(i-1)+c2,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
 !           T_SR(3*(i-1)+c2,3*(j-1)+c1,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
 !           T_SR(3*(j-1)+c2,3*(i-1)+c1,:) = T_SR(3*(i-1)+c1,3*(j-1)+c2,:)
          end do
        end do
      end do
    end do
    B_mat = B_mat + T_SR

!    write(*,*) "B_mat:"
!    do i = 1, 6
!      write(*,*) B_mat(i,1:6,1)
!    end do

!    write(*,*) "T_SR:", T_SR(1,1:9,1)
!    write(*,*) "T_SR:", T_SR(2,1:9,1)
!    write(*,*) "T_SR:", T_SR(3,1:9,1)
!    write(*,*) "T_SR:", T_SR(4,1:9,1)
!    write(*,*) "T_SR:", T_SR(5,1:9,1)
!    write(*,*) "T_SR:", T_SR(6,1:9,1)
!    write(*,*) "T_SR:", T_SR(7,1:9,1)
!    write(*,*) "T_SR:", T_SR(8,1:9,1)
!    write(*,*) "T_SR:", T_SR(9,1:9,1)
!    write(*,*) "B_mat:", B_mat(1,1:9,1)
!    write(*,*) "B_mat:", B_mat(2,1:9,1)
!    write(*,*) "B_mat:", B_mat(3,1:9,1)
!    write(*,*) "B_mat:", B_mat(4,1:9,1)
!    write(*,*) "B_mat:", B_mat(5,1:9,1)
!    write(*,*) "B_mat:", B_mat(6,1:9,1)
!    write(*,*) "B_mat:", B_mat(7,1:9,1)
!    write(*,*) "B_mat:", B_mat(8,1:9,1)
!    write(*,*) "B_mat:", B_mat(9,1:9,1)
!    write(*,*) "B_mat end:", B_mat(172,172:180,1)
!    write(*,*) "B_mat end:", B_mat(173,172:180,1)
!    write(*,*) "B_mat end:", B_mat(174,172:180,1)
!    write(*,*) "B_mat end:", B_mat(175,172:180,1)
!    write(*,*) "B_mat end:", B_mat(176,172:180,1)
!    write(*,*) "B_mat end:", B_mat(177,172:180,1)
!    write(*,*) "B_mat end:", B_mat(178,172:180,1)
!    write(*,*) "B_mat end:", B_mat(179,172:180,1)
!    write(*,*) "B_mat end:", B_mat(180,172:180,1)

!    write(*,*) "T_SR full:"
!    do i = 1, 3*n_sites
!      write(*,*) T_SR(i,:,1)
!    end do

!    write(*,*) "B_mat full:"
!    do k = 1,3*n_sites
!      write(*,*) B_mat(k,:,1)
!    end do

!    write(*,*) "B_mat:"
!    do i = 1, 6
!      write(*,*) B_mat(i,1:6,1)
!    end do

    A_mat = B_mat

    do i = 1, 11
      call dgetrf(3*n_sites, 3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, info)
      call dgetri(3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
    end do

!    VL = 0.d0
!    call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_mat(:,:,1), 3*n_sites, B_mat(:,:,1), 3*n_sites, &
!                 0.d0, VL, 3*n_sites)

!    write(*,*) "A_mat:"
!    do i = 1, 3*n_sites
!      write(*,*) A_mat(i,1:3*n_sites,1)
!    end do

    do k = 1, 11
      do i = 1, n_sites
        A_i = 0.d0
        do j = 1, n_sites
          A_i = A_i + A_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,k)
        end do
        alpha_SCS(i,k) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
      end do
    end do

!    write(*,*) "alpha_SCS:"
!    write(*,*) alpha_SCS(:,1)

    A_LR = 0.d0

    do k = 1, 11
      do i = 1, n_sites
        do c1 = 1, 3
          A_LR(3*(i-1)+c1,3*(i-1)+c1,k) = alpha_SCS(i,k)
        end do
      end do
    end do

!    write(*,*) "alpha_SCS:", alpha_SCS(:,1)
!    write(*,*) "A_LR:", A_LR(1,1:6,1)
!    write(*,*) "A_LR:", A_LR(2,1:6,1)
!    write(*,*) "A_LR:", A_LR(3,1:6,1)
!    write(*,*) "A_LR:", A_LR(4,1:6,1)
!    write(*,*) "A_LR:", A_LR(5,1:6,1)
!    write(*,*) "A_LR:", A_LR(6,1:6,1)

    do k = 1, n_pairs
      j = neighbors_list(k)
      r0_ii_SCS(k) = r0_ii(k) * (alpha_SCS(j,1)/neighbor_alpha0(k))**(1.d0/3.d0)
    end do

    k = 0
    do i = 1, n_sites
      k=k+1
      do j2 = 2, n_neigh(i)
        k=k+1
        j = neighbors_list(k)
        f_damp_SCS(k) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/(sR*(r0_ii_SCS(n_sites*(i-1)+1) + r0_ii_SCS(k))) - 1.d0 ) ) )
      end do
    end do

    T_LR = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        T_LR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = f_damp_SCS(k) * T_func(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
        T_LR(3*(j-1)+1:3*(j-1)+3,3*(i-1)+1:3*(i-1)+3) = T_LR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
      end do
    end do

!    do i = 1, n_sites
!      write(*,*) "T_LR diag:", T_LR(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3)
!    end do
!    write(*,*) "T_LR:", T_LR(1,1:6)
!    write(*,*) "T_LR:", T_LR(2,1:6)
!    write(*,*) "T_LR:", T_LR(3,1:6)
!    write(*,*) "T_LR:", T_LR(4,1:6)
!    write(*,*) "T_LR:", T_LR(5,1:6)
!    write(*,*) "T_LR:", T_LR(6,1:6)

    I_mat = 0.d0
    do i = 1, 3*n_sites
      I_mat(i,i) = 1.d0
    end do

    AT = 0.d0
    do k = 1, 11
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR(:,:,k), 3*n_sites, T_LR, 3*n_sites, &
                0.d0, AT(:,:,k), 3*n_sites)
!      AT(:,:,k) = matmul(A_mat(:,:,k),T_LR) 
    end do

!    write(*,*) "AT:", AT(1,1:6,1)
!    write(*,*) "AT:", AT(2,1:6,1)
!    write(*,*) "AT:", AT(3,1:6,1)
!    write(*,*) "AT:", AT(4,1:6,1)
!    write(*,*) "AT:", AT(5,1:6,1)
!    write(*,*) "AT:", AT(6,1:6,1)

    do k = 1,11
      WR = 0.d0
      VL = 0.d0
      VR = 0.d0
      call dgeev('n', 'v', 3*n_sites, I_mat-AT(:,:,k), 3*n_sites, WR, WI, VL, 3*n_sites, VR, 3*n_sites, &
                 work_arr, 12*n_sites, info)
!      write(*,*) "WI:", WI 
      logMapo = 0.d0
!      write(*,*) "WR:"
      do i = 1, 3*n_sites
!        write(*,*) WR(i)
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

!    write(*,*) "logIAT:", logIAT(1,1:6,1)
!    write(*,*) "logIAT:", logIAT(2,1:6,1)
!    write(*,*) "logIAT:", logIAT(3,1:6,1)
!    write(*,*) "logIAT:", logIAT(4,1:6,1)
!    write(*,*) "logIAT:", logIAT(5,1:6,1)
!    write(*,*) "logIAT:", logIAT(6,1:6,1)
!    write(*,*) "logIAT full:"
!    do i = 1, 3*n_sites
!      write(*,*) logIAT(i,:,1)
!    end do

    integrand = 0.d0
    do k = 1,11
      do i = 1,3*n_sites
        integrand(k) = integrand(k) + logIAT(i,i,k)
      end do 
    end do

!    write(*,*) "Integrand:", integrand

!    NOTE: Include misc.mod before commenting this out
    call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
!    E_MBD = integral/(2.d0*pi)

!    integral = 0.d0
!    do k = 1,10
!      integral = integral + 0.5d0 * (integrand(k)+integrand(k+1)) * (omegas(k+1)-omegas(k))
!    end do
    E_MBD = integral/(2.d0*pi)

    ! Conversion to eV:
    E_MBD = E_MBD * 27.211386245988
    write(*,*) "E_MBD:", E_MBD

!    if( do_timing) then
!      call cpu_time(time2)
!      write(*,*) "vdw: creating pair quantities 3:", time2-time1
!      call cpu_time(time1)
!    end if


!    if( do_timing) then
!      call cpu_time(time2)
!      write(*,*) "vdw: computing energies:", time2-time1
!      call cpu_time(time1)
!    end if

!   Force calculation starts here:
    dT = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        do c1 = 1,3
          do c2 = 1,3
            do c3 = 1,3
              dT(3*(i-1)+c1,3*(j-1)+c2,c3) = (-15.d0 * xyz_H(c1,k) * xyz_H(c2,k) * xyz_H(c3,k))/rjs_H(k)**7
              if (c1 == c2) then
                dT(3*(i-1)+c1,3*(j-1)+c2,c3) = dT(3*(i-1)+c1,3*(j-1)+c2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c3,k)
              end if
              if (c2 == c3) then
                dT(3*(i-1)+c1,3*(j-1)+c2,c3) = dT(3*(i-1)+c1,3*(j-1)+c2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c1,k)
              end if
              if (c1 == c3) then
                dT(3*(i-1)+c1,3*(j-1)+c2,c3) = dT(3*(i-1)+c1,3*(j-1)+c2,c3) + 3.d0/rjs_H(k)**5 * xyz_H(c2,k)
              end if
            end do
          end do
        end do
      end do
    end do

!    write(*,*) "dT:", dT(1,1:6,1)
!    write(*,*) "dT:", dT(2,1:6,1)
!    write(*,*) "dT:", dT(3,1:6,1)
!    write(*,*) "dT:", dT(4,1:6,1)
!    write(*,*) "dT:", dT(5,1:6,1)
!    write(*,*) "dT:", dT(6,1:6,1)

    f_damp_der = 0.d0
    g_func_der = 0.d0
    h_func_der = 0.d0
    dT_SR = 0.d0

    do a = 1, n_sites
      k = 0
      do i = 1, n_sites
        k = k+1
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          if (a == i .or. a == j) then
!            write(*,*) i, j
            sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
            do c3 = 1, 3
              f_damp_der(i,j,c3) = d/(sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))) * f_damp(k)**2 * &
                                   exp( -d*(rjs_H(k)/(sR*(r0_ii(n_sites*(i-1)+1) + &
                                   r0_ii(k))) - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
!              if ( a == 1 .and. i == 2 .and. j == 1 .and. c3 == 1) then
!                write(*,*) "i, j, f_damp_der:", i, j, f_damp_der(i,j,c3)
!              end if
!                              f_damp(k) =               exp( -d*(rjs_H(k)/(sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))) - 1.d0 ) )
              g_func_der(i,j,c3,:) = 4.d0/sqrt(pi) * rjs_H(k)/sigma_ij**3 * xyz_H(c3,k) * &
                                     exp(-rjs_H(k)**2/sigma_ij**2) 
              do c1 = 1, 3
                do c2 = 1,3
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
                  h_func_der(i,j,c1,c2,c3,:) = coeff_h_der * terms
                  dT_SR(3*(i-1)+c1,3*(j-1)+c2,a,c3,:) = f_damp_der(i,j,c3) * T_func(3*(i-1)+c1,3*(j-1)+c2) * &
                                                        g_func(i,j,:) - (1.d0 - f_damp(k)) * dT(3*(i-1)+c1,3*(j-1)+c2,c3) * &
                                                        g_func(i,j,:) - g_func_der(i,j,c3,:) * (1.d0 - f_damp(k)) * &
                                                        T_func(3*(i-1)+c1,3*(j-1)+c2) - f_damp_der(i,j,c3) * h_func(i,j,c1,c2,:) + &
                                                        h_func_der(i,j,c1,c2,c3,:) * (1.d0 - f_damp(k)) 
!                  if (a == 1 .and. c3 == 1 .and. c1 == 2 .and. c2 == 3) then
!                    write(*,*) i-1, j-1
!                    write(*,*) "1st term", f_damp_der(i,j,c3) * T_func(3*(i-1)+c1,3*(j-1)+c2) * g_func(i,j,:)
!                    write(*,*) "2nd term", -(1.d0 - f_damp(k)) * dT(3*(i-1)+c1,3*(j-1)+c2,c3) * g_func(i,j,:)
!                    write(*,*) "3rd term", -g_func_der(i,j,c3,:) * (1.d0 - f_damp(k)) * T_func(3*(i-1)+c1,3*(j-1)+c2)
!                    write(*,*) "4th term", -f_damp_der(i,j,c3) * h_func(i,j,c1,c2,:)
!                    write(*,*) "5th term", h_func_der(i,j,c1,c2,c3,:) * (1.d0 + f_damp(k))
!                    write(*,*) "T_ij", T_func(3*(i-1)+c1,3*(j-1)+c2)
!                    write(*,*) "T_ij_der", dT(3*(i-1)+c1,3*(j-1)+c2,c3)
!                    write(*,*) "damp", f_damp(k)
!                    write(*,*) "damp_der", f_damp_der(i,j,c3)
!                    write(*,*) "R_vdW_ij", r0_ii(n_sites*(i-1)+1) + r0_ii(k)
!                    write(*,*) "r_ij", rjs_H(k)
!                    write(*,*) "r_vec", xyz_H(:,k)
!                    write(*,*) "d, sR", d, sR
!                    write(*,*) "g_ij", g_func(i,j,:)
!                    write(*,*) "g_ij_der", g_func_der(i,j,c3,:)
!                    write(*,*) "h_ij", h_func(i,j,c1,c2,:)
!                    write(*,*) "h_ij_der", h_func_der(i,j,c1,c2,c3,:)
!                  end if
                end do
              end do
            end do
            if (a == i) then
              dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:) = -dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,:,:)
            end if
          end if
        end do
      end do
    end do

!    dA_mat = 0.d0
!    do k = 1, 11
!      do a = 1, n_sites
!        do c3 = 1, 3
!          dT_SR_A_mat = 0.d0
!          call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, A_mat(:,:,k), 3*n_sites, &
!                     0.d0, dT_SR_A_mat, 3*n_sites)
!          call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, -A_mat(:,:,k), 3*n_sites, dT_SR_A_mat, 3*n_sites, &
!                     0.d0, dA_mat(:,:,a,c3,k), 3*n_sites)
!        end do
!      end do
!    end do

!    write(*,*) "dA_mat:", dA_mat(1,1:6,1,1,1)
!    write(*,*) "dA_mat:", dA_mat(2,1:6,1,1,1)
!    write(*,*) "dA_mat:", dA_mat(3,1:6,1,1,1)
!    write(*,*) "dA_mat:", dA_mat(4,1:6,1,1,1)
!    write(*,*) "dA_mat:", dA_mat(5,1:6,1,1,1)
!    write(*,*) "dA_mat:", dA_mat(6,1:6,1,1,1)

!    dalpha = 0.d0
!    do k = 1, 11
!      do i = 1, n_sites
!        do a = 1, n_sites
!          do c3 = 1, 3
!            A_i = 0.d0
!            do j = 1, n_sites
!              A_i = A_i + dA_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,c3,k)
!            end do
!            do j2 = 1, 3
!              dalpha(i,a,c3,k) = dalpha(i,a,c3,k) + A_i(j2,j2)
!            end do
!            dalpha(i,a,c3,k) = 1.d0/3.d0 * dalpha(i,a,c3,k)
!          end do
!        end do
!      end do
!    end do

!    write(*,*) "d, sR", d, sR
!    write(*,*) "r0_ii:", r0_ii

!    write(*,*) "T_func:"
!    do i = 1, 6
!      write(*,*) T_func(i,1:6)
!    end do

    dB_mat = dT_SR
    if (do_hirshfeld_gradients) then
      dB_mat_v = 0.d0
      dT_SR_v = 0.d0
      coeff_fdamp = 0.d0
      coeff_der = 0.d0
      k = 0
      do i = 1, n_sites
        k = k+1
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
          S_vdW_ij = sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))
          exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
          ! This should includ T_ij * g_func + h_func
!          if ( i == 2 .and. j == 1) then
!            write(*,*) "i, j, df_damp:", i, j, -(d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
!              ( r0_ii(4)/hirshfeld_v(2) * hirshfeld_v_cart_der_H(1,4) + r0_ii(1)/hirshfeld_v(1) * hirshfeld_v_cart_der_H(1,1) )
!            write(*,*) "i, j, dg:", 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij(1)**2) * (rjs_H(k)/sigma_ij(1))**2 * &
!              (-rjs_H(k)/(3.d0 * sigma_ij(1)**3)) * &
!              (sigma_i(2,1)**2/hirshfeld_v(2) * hirshfeld_v_cart_der_H(1,4) + sigma_i(1,1)**2/hirshfeld_v(1) * &
!              hirshfeld_v_cart_der_H(1,1))
!            write(*,*) "i, j, dh:", 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij(1)**2) * xyz_H(1,k)*xyz_H(1,k) / &
!              (sigma_ij(1)**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij(1))**2) * &
!              (sigma_i(2,1)**2/hirshfeld_v(2) * hirshfeld_v_cart_der_H(1,4) + sigma_i(1,1)**2/hirshfeld_v(1) * &
!              hirshfeld_v_cart_der_H(1,1))
!          end if
          do c1 = 1, 3
            do c2 = 1, 3
              coeff_fdamp(i,j,c1,c2,:) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                 (-T_func(3*(i-1)+c1,3*(j-1)+c2) * g_func(i,j,:) + h_func(i,j,c1,c2,:))
!              if (c1 == c2) then
!                coeff_h_der = 2*xyz_H(c1,k)*xyz_H(c2,k)/sigma_ij**2 - 1.d0
!              else
!                coeff_h_der = 2*xyz_H(c1,k)*xyz_H(c2,k)/sigma_ij**2                
!              end if
!              coeff_der(i,j,c1,c2,:) = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * 1.d0/(3*sigma_ij**7) * &
!                                 (1.d0-f_damp(k)) * 1.d0/3.d0 * rjs_H(k)**2/sigma_ij**2 * coeff_h_der
              dg = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * (rjs_H(k)/sigma_ij)**2 * &
                   (-rjs_H(k)/(3.d0 * sigma_ij**3))
              dh = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * xyz_H(c1,k)*xyz_H(c2,k) / &
                   (sigma_ij**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij)**2)
!              if ( i == 2 .and. j == 1 .and. c1 == 1 .and. c2 == 1 ) then
!                write(*,*) "dg:", dg
!                write(*,*) "dh:", dh
!                write(*,*) "dT_SR_v:", (1.d0-f_damp(k)) * (-T_func(3*(i-1)+c1,3*(j-1)+c2)*dg+dh) * &
!                (sigma_i(2,1)**2/hirshfeld_v(2) * hirshfeld_v_cart_der_H(1,4) + sigma_i(1,1)**2/hirshfeld_v(1) * &
!                hirshfeld_v_cart_der_H(1,1)) + coeff_fdamp(i,j,c1,c2,:) * &
!                ( r0_ii(4)/hirshfeld_v(2) * hirshfeld_v_cart_der_H(1,4) + r0_ii(1)/hirshfeld_v(1) * hirshfeld_v_cart_der_H(1,1) )
!              end if
              coeff_der(i,j,c1,c2,:) = (1.d0-f_damp(k)) * (-T_func(3*(i-1)+c1,3*(j-1)+c2)*dg+dh)
            end do
          end do
!          coeff_der(i,j,:) = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * 1.d0/(3*sigma_ij**5) * &
!                             (1.d0-f_damp(k)) * (-1.d0 + 2.d0/3.d0 * rjs_H(k)**2/sigma_ij**2)
        end do
      end do

!      write(*,*) "coeff_fdamp:", coeff_fdamp
!      write(*,*) "coeff_der:", coeff_der(:,:,1)

      k2 = 0
      k3 = 0
      do i = 1, n_sites
        do a2 = 1, n_neigh(i)
          k2 = k2+1
          a = neighbors_list(k2)
!          write(*,*) i, a, hirshfeld_v_cart_der_H(1,k2)
          do c3 = 1, 3
            do c1 = 1, 3
              dB_mat_v(3*(i-1)+c1,3*(i-1)+c1,a,c3,:) = -1.d0/(hirshfeld_v(i)*alpha_i(i,:)) * &
                            hirshfeld_v_cart_der_H(c3,k2)
            end do
            do k = 1, 11
              do j2 = 2, n_neigh(i)
                j = neighbors_list(k3+j2)
!                write(*,*) "i, j:", i, j
!                if (k == 1) then
!                  write(*,*) i, j, coeff_der(i,j,1) * r0_ii(n_sites*(i-1)+1) * hirshfeld_v_cart_der(c3,k2)
!                end if
                do c1 = 1, 3
!                  dT_SR_v(3*(i-1)+c1,3*(j-1)+c1,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c1,a,c3,k) + &
!                    coeff_der(i,j,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der(c3,k2))
                  do c2 = 1, 3
                    dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                      coeff_der(i,j,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
                    dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                      coeff_fdamp(i,j,c1,c2,k) * r0_ii(n_sites*(i-1)+1)/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
                  end do
                end do      
              end do
            end do
          end do
        end do
        k3 = k3 + n_neigh(i)
      end do
      k2 = 0
      k3 = 0
      do j = 1, n_sites
        do a2 = 1, n_neigh(j)
          k2 = k2+1
          a = neighbors_list(k2)
          do c3 = 1, 3
            do k = 1, 11
              do i2 = 2, n_neigh(j)
                i = neighbors_list(k3+i2)
!                write(*,*) "j, i:", j, i
!                if (k == 1) then
!                  write(*,*) i, j, coeff_der(i,j,1) * r0_ii(n_sites*(j-1)+1) * hirshfeld_v_cart_der(c3,k2)
!                end if
                do c1 = 1, 3
!                  dT_SR_v(3*(i-1)+c1,3*(j-1)+c1,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c1,a,c3,k) + &
!                    coeff_der(i,j,k) * (sigma_i(j,k)**2/hirshfeld_v(j) * hirshfeld_v_cart_der(c3,k2))
                  do c2 = 1, 3
                    dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                      coeff_der(i,j,c1,c2,k) * (sigma_i(j,k)**2/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2))
                    dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                      coeff_fdamp(i,j,c1,c2,k) * r0_ii(n_sites*(j-1)+1)/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2)
                  end do
                end do
              end do
            end do
          end do
        end do
        k3 = k3 + n_neigh(j)
      end do

      dB_mat_v = dB_mat_v + dT_SR_v

!      k = 0
!      do i = 1, n_sites
!        k = k+1
!        do j2 = 2, n_neigh(i)
!          k = k+1
!          j = neighbors_list(k)
!          sigma_ij = sqrt(sigma_i(i,:)**2 + sigma_i(j,:)**2)
!          S_vdW_ij = sR*(r0_ii(n_sites*(i-1)+1) + r0_ii(k))
!          exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
!          f_damp_der_v_coeff = -d*sR*rjs_H(k)/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2
!          g_func_der_v_coeff = -4.d0/(3*sqrt(pi)) * rjs_H(k)**3/sigma_ij**5 * exp(-rjs_H(k)**2/sigma_ij**2)
!          k2 = 0
!          dv_i = 0.d0
!          dv_j = 0.d0
!          do a = 1, n_sites
!            k2 = k2+1
!            if (neighbors_list(n_sites*(i-1)+k2) == a) then
!              dv_i = hirshfeld_v_cart_der(:,n_sites*(i-1)+k2)
!            end if
!            if (neighbors_list(n_sites*(j-1)+k2) == a) then
!              dv_j = hirshfeld_v_cart_der(:,n_sites*(j-1)+k2)
!            end if
!            dR_vdW(1:3) = r0_ii(n_sites*(i-1)+1)/hirshfeld_v_neigh(n_sites*(i-1)+1) * dv_i + &                
!                          r0_ii(k)/hirshfeld_v_neigh(k) * dv_j
!            f_damp_der_v(i,j,a,:) = f_damp_der_v_coeff * dR_vdW
!            dsigma = 0.d0
!            do c3 = 1, 3
!              dsigma(c3,:) = sigma_i(i,:)**2/hirshfeld_v_neigh(n_sites*(i-1)+1) * dv_i(c3) + &
!                             sigma_i(j,:)**2/hirshfeld_v_neigh(k) * dv_j(c3)
!              g_func_der_v(i,j,a,c3,:) = g_func_der_v_coeff * dsigma(c3,:)
!              do c1 = 1, 3
!                if (a == i) then
!                  dB_mat_v(3*(i-1)+c1,3*(i-1)+c1,a,c3,:) = -1.d0/(hirshfeld_v_neigh(n_sites*(i-1)+1)*alpha_i(i,:)) * &
!                                                           dv_i(c3)
!                end if
!                do c2 = 1, 3
!                  coeff_h_der = 4.d0/sqrt(pi) * xyz_H(c1,k) * xyz_H(c2,k)/(sigma_ij**5*rjs_H(k)**2) * &
!                              exp(-rjs_H(k)**2/sigma_ij**2) * (-1.d0 + 2.d0/3.d0 * rjs_H(k)**2 / sigma_ij**2)
!                  h_func_der_v(i,j,a,c1,c2,c3,:) = coeff_h_der * dsigma(c3,:)
!                  dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,:) = -f_damp_der_v(i,j,a,c3) * (-T_func(3*(i-1)+c1,3*(j-1)+c2) * &
!                                                       g_func(i,j,:) + h_func(i,j,c1,c2,:)) + (1.d0 - f_damp(k)) * &
!                                                       (-T_func(3*(i-1)+c1,3*(j-1)+c2) * g_func_der_v(i,j,a,c3,:) + &
!                                                       h_func_der_v(i,j,a,c1,c2,c3,:))
!                end do
!              end do
!            end do
!          end do
!        end do
!      end do

!      dB_mat_v = dB_mat_v + dT_SR_v
!      dA_mat_v = 0.d0
!      do k = 1, 11
!        do a = 1, n_sites
!          do c3 = 1, 3
!            dT_SR_A_mat = 0.d0
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, dB_mat_v(:,:,a,c3,k), 3*n_sites, A_mat(:,:,k), 3*n_sites, &
!                       0.d0, dT_SR_A_mat, 3*n_sites)
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, -A_mat(:,:,k), 3*n_sites, dT_SR_A_mat, 3*n_sites, &
!                       0.d0, dA_mat_v(:,:,a,c3,k), 3*n_sites)
!          end do
!        end do
!      end do
!      dalpha_v = 0.d0
!      do k = 1, 11
!        do i = 1, n_sites
!          do a = 1, n_sites
!            do c3 = 1, 3
!              A_i = 0.d0
!              do j = 1, n_sites
!                A_i = A_i + dA_mat_v(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,a,c3,k)
!              end do
!              do j2 = 1, 3
!                dalpha_v(i,a,k,c3) = dalpha_v(i,a,k,c3) + A_i(j2,j2)
!              end do
!            end do
!          end do
!        end do
!      end do
!      dalpha_v = dalpha_v/3.d0
!      dalpha_v = 0.d0
!      do i = 1, n_sites
!        do a = 1, n_sites
!          do j = 1, n_sites
!            do c1 = 1, 3
!              do i2 = 1, 3*n_sites
!                do j2 = 1, 3*n_sites
!                  dalpha_v(i,a,:) = dalpha_v(i,a,:) - A_mat(3*(i-1)+c1,i2,:) * dB_mat_v(i2,j2,a,:) * A_mat(j2,3*(j-1)+c1,:) 
!                end do
!              end do
!            end do
!          end do
!        end do
!      end do
!     k = 0
!      dalpha_v_r = 0.d0
!      do b = 1, n_sites
!        do j2 = 1, n_neigh(b)
!          k = k+1
!          a = neighbors_list(k)
!          do i = 1, n_sites
!            do c3 = 1, 3
!              dalpha_v_r(i,a,c3,:) = dalpha_v_r(i,a,c3,:) + dalpha_v(i,b,:) * hirshfeld_v_cart_der(c3,k)
!              write(*,*) "b, a, i, c3, dalpha_v, hirshfeld_v_cart_der"
!              write(*,*) b, a, i, c3, dalpha_v(i,b,1), hirshfeld_v_cart_der(c3,k)
!            end do
!          end do
!        end do
!      end do
!      write(*,*) "k, size of hirshfeld der:", k, size(hirshfeld_v_cart_der(1,:))
!      dalpha = dalpha + dalpha_v
      dB_mat = dB_mat + dB_mat_v
    end if

!    write(*,*) "dB_mat:"
!    do i = 1, 6
!      write(*,*) dB_mat(i,1:6,1,1,1)
!    end do

!    write(*,*) "dT_SR:"
!    do i = 1, 6
!      write(*,*) dT_SR(i,1:6,1,1,1)
!    end do

!    write(*,*) "dT_SR_v:"
!    do i = 1, 6
!      write(*,*) dT_SR_v(i,1:6,1,1,1)
!    end do

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

!    do i = 55, 60
!      write(*,*) "f_damp_der_v:", f_damp_der_v(i,55:60,60)
!    end do
!    do i = 55, 60
!      write(*,*) "g_func_der_v:", g_func_der_v(i,55:60,60,1)
!    end do
!    do i = 55, 60
!      write(*,*) "h_func_der_v:", h_func_der_v(i,55:60,60,1,1,1)
!    end do
!    do i = 175, 180
!      write(*,*) "dT_SR_v:", dT_SR_v(i,175:180,60,1)
!    end do
!    write(*,*) "dB_mat_v:"
!    do i = 1, 3*n_sites
!      write(*,*) dB_mat_v(i,1:3*n_sites,1,1)
!    end do
!    write(*,*) "dA_mat_v:"
!    do i = 1, 3*n_sites
!      write(*,*) dA_mat_v(i,1:3*n_sites,1,1)
!    end do
!    write(*,*) "dalpha_v:", dalpha_v(:,1,1)
!    write(*,*) "dalpha_v_r:", dalpha_v_r(:,1,1,1)
!    write(*,*) "dalpha:", dalpha(:,1,1,1)

!    do a = 1, n_sites
!      write(*,*) "dalpha", dalpha(1,a,:,1)
!    end do

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

!    do i = 1, 9
!      write(*,*) "dA_LR", dA_LR(i,1:9,1,1,1)
!    end do

    f_damp_der_SCS = 0.d0
    f_damp_der = 0.d0 ! This is cleared so we can recalculate it with SCS values
    dT_LR = 0.d0
    do a = 1, n_sites
      k = 0
      do i = 1, n_sites
        k = k+1
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          R_vdW_SCS_ij = r0_ii_SCS(n_sites*(i-1)+1) + r0_ii_SCS(k)
          S_vdW_ij = sR*R_vdW_SCS_ij
          do c3 = 1, 3
            dS_vdW_ij = sR/3.d0 * ( r0_ii_SCS(n_sites*(i-1)+1)/alpha_SCS(i,1) * dalpha(i,a,c3,1) + &
                                    r0_ii_SCS(k)/alpha_SCS(j,1) * dalpha(j,a,c3,1) )
            f_damp_der_SCS(i,j,a,c3) = -(d*rjs_H(k))/S_vdW_ij**2 * f_damp_SCS(k)**2 * &
                                       exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0)) * dS_vdW_ij
!            if (a == 1 .and. c3 == 1 .and. i == 1 .and. j == 2) then
!              write(*,*) i-1, j-1
!              write(*,*) "f_damp_der_SCS", f_damp_der_SCS(i,j,a,1)
!              write(*,*) "dS_vdW_ij", dS_vdW_ij
!              write(*,*) "R_vdW_SCS_ij", R_vdW_SCS_ij
!              write(*,*) "S_vdW_ij", S_vdW_ij
!            end if
!            if (a == 1 .and. c3 == 1 .and. i == 2 .and. j == 1) then
!              write(*,*) i-1, j-1
!              write(*,*) "f_damp_der_SCS", f_damp_der_SCS(i,j,a,1)
!              write(*,*) "dS_vdW_ij", dS_vdW_ij
!              write(*,*) "R_vdW_SCS_ij", R_vdW_SCS_ij
!              write(*,*) "S_vdW_ij", S_vdW_ij
!            end if
            do c1 = 1, 3
              do c2 = 1, 3
                dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) + &
                                  T_func(3*(i-1)+c1,3*(j-1)+c2) * f_damp_der_SCS(i,j,a,c3) 
              end do
            end do
          end do
          if (a == i .or. a == j) then
            do c3 = 1, 3
            f_damp_der(i,j,c3) = d/S_vdW_ij * f_damp_SCS(k)**2 * &
                                 exp( -d*(rjs_H(k)/S_vdW_ij - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
              do c1 = 1, 3
                do c2 = 1, 3
                  if (a == i) then
                    dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) - &
                                       T_func(3*(i-1)+c1,3*(j-1)+c2) * f_damp_der(i,j,c3) - &
                                       dT(3*(i-1)+c1,3*(j-1)+c2,c3) * f_damp_SCS(k)
                  else if (a == j) then
                    dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) = dT_LR(3*(i-1)+c1,3*(j-1)+c2,a,c3) + &
                                       T_func(3*(i-1)+c1,3*(j-1)+c2) * f_damp_der(i,j,c3) + &
                                       dT(3*(i-1)+c1,3*(j-1)+c2,c3) * f_damp_SCS(k)
                  end if
                end do
              end do
            end do
          end if
        end do
      end do
    end do

!    do i = 1, n_sites
!      write(*,*) "dT_LR", dT_LR(i,1:6,1,1)
!    end do

!    do i = 1, n_sites
!      write(*,*) "dA_LR", dA_LR(i,1:6,1,2,1)
!    end do

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

!    do i = 1, 9
!      write(*,*) "G_mat:", G_mat(i,1:6,1,1,1)
!    end do

!    do a = 1, n_sites
!      write(*,*) "Force integrand:", force_integrand(a,:,1)
!    end do

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
!   Compute MBD forces
!    if( do_forces )then
!      forces0 = 0.d0
!      virial = 0.d0
!    end if

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "vdw: computing forces:", time2-time1
    end if
 
!   Clean up
    deallocate( neighbor_c6_ii, neighbor_c6_ij, r0_ij, exp_damp, f_damp, c6_ij_free, r6, is_in_buffer, i_buffer, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, A_LR, r0_ii_SCS, f_damp_SCS, I_mat, AT, logIAT, logMapo, VR, VRinv, WR, WI, VL, integrand, &
                v_arr, dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, dalpha, dA_LR, dT_LR, f_damp_der_SCS, &
                invIAT, G_mat, force_integrand, forces_MBD, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                f_damp_der_v, g_func_der_v, h_func_der_v, dB_mat, dB_mat_v, dalpha_v, dA_mat_v, dR_vdW, dv_i, dv_j, &
                dsigma, g_func_der_v_coeff, coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H )
    if( do_forces )then
      deallocate( pref_force1, pref_force2, r6_der )
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
