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
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod

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
  subroutine get_scs_polarizabilities( hirshfeld_v, hirshfeld_v_cart_der, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, hirshfeld_v_neigh, &
                                       sR, d, c6_ref, r0_ref, alpha0_ref, do_derivatives, alpha_SCS0, dalpha_full, &
                                       energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:), hirshfeld_v(:), hirshfeld_v_neigh(:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_derivatives
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:), alpha_SCS0(:,:), dalpha_full(:,:,:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), &
                           f_damp(:), neighbor_alpha0(:), &
                           T_func(:), h_func(:,:), g_func(:,:), &
                           omegas(:), omega_i(:), &
                           alpha_i(:,:), sigma_i(:,:), sigma_ij(:), T_SR(:,:,:), B_mat(:,:,:), &
                           rjs_H(:), xyz_H(:,:), A_mat(:,:,:), work_arr(:), A_i(:,:), &
                           dT(:,:), dT_SR(:,:,:,:,:), f_damp_der(:,:), &
                           g_func_der(:,:,:), h_func_der(:,:,:), dA_mat(:,:,:,:,:), &
                           dalpha(:,:,:), &
                           coeff_h_der(:), terms(:), dT_SR_A_mat(:,:), dT_SR_v(:,:,:,:,:), &
                           dB_mat(:,:,:,:,:), dB_mat_v(:,:,:,:,:), &
                           coeff_der(:,:,:,:,:), coeff_fdamp(:,:,:,:,:), dg(:), &
                           dh(:), hirshfeld_v_cart_der_H(:,:), &
                           alpha_k(:,:), sigma_k(:,:), s_i(:), s_j(:), &
                           T_SR_i(:,:,:), B_mat_i(:,:,:), alpha_SCS_i(:,:), A_mat_i(:,:,:), alpha_SCS(:,:), &
                           dB_mat_n(:,:,:,:,:), dA_mat_n(:,:,:,:,:), dBA_n(:,:), I_mat(:,:), &
                           dalpha_n(:,:,:), hirshfeld_v_neigh_H(:), hirshfeld_v_cart_der2(:,:), &
                           rjs2(:), xyz2(:,:), hirshfeld_v_neigh2(:), rjs_central2(:), rjs_central(:), &
                           xyz_central2(:,:), xyz_central(:,:), dalpha2(:,:,:), I_aT(:,:,:), a_vec(:,:,:), &
                           regularization(:,:), I_aT_a_vec(:,:), I_aT2(:,:), BTB_reg(:,:,:), B_reg(:,:,:), &
                           BTB_reg_copy(:,:), a_SCS(:,:,:), da_vec(:,:,:,:,:), dBTB(:,:), vect1(:,:), &
                           vect2(:,:), vect3(:,:), vect4(:,:), vect_temp(:,:), da_SCS(:,:,:,:,:), AdB_n(:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, dist, f_damp_SCS_ij, t1, t2, &
              reg_param
    integer, allocatable :: ipiv(:), n_sneigh(:), sneighbors_list(:), p_to_i(:), i_to_p(:), &
                            neighbors_list2(:), local_neighbors(:), neighbor_species2(:), &
                            neighbor_species_H(:), neighbors_list3(:), n_sneigh_vder(:), &
                            sneighbors_list_vder(:), n_sneigh_list2(:), &
                            n_sneigh_list(:), n_ssites_list(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, n_order, n_tot, n_spairs, &
               n_ssites, n_spairs_vder, n_freq
    integer :: i, i2, j, j2, k, k2, k3, a, a2, c1, c2, c3, lwork, b, p, q, n_count
    logical :: do_timing = .true., do_hirshfeld_gradients = .true., &
               total_energy = .true.
    logical, allocatable :: i0_buffer2(:), ij_buffer2(:), i0_buffer(:), ij_buffer(:)
    type(psb_ctxt_type) :: icontxt
    integer(psb_ipk_) ::  iam, np, ip, jp, idummy, nr, nnz, info_psb
    type(psb_desc_type) :: desc_a
    type(psb_dspmat_type) :: A_sp
    real*8, allocatable :: x_vec(:,:), b_vec(:,:)
    integer(psb_lpk_), allocatable :: ia(:), ja(:), myidx(:)
    real(psb_dpk_), allocatable :: val(:), val_xv(:,:), val_bv(:,:)
    type(psb_dprec_type)  :: prec
    character(len=20) :: ptype

!   IMPORTANT NOTE ABOUT THE DERIVATIVES:
!   If rcut < rcut_soap, the derivatives in the new implementation omit the terms that fall outside of rcut.
!   This means that the finite difference and the analytical derivative do not match in the present version
!   in this special case. The old implementation gives the correct derivatives even in this case. I will
!   probably fix this at some point but it will require using the full n_neigh(i) again, instead of 
!   n_ssites, to construct the sneighbors_list.

    write(*,*) "rcut", rcut

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

!   This should allow to only take a subset of atoms for parallelization:

!   Number of frequencies
    n_freq = size(alpha_SCS0, 2)

!   This implementation assumes that rcut is the largest cutoff, that is, the neigbhbors_list contains only the atoms within the vdW cutoff.
!   The implementation matches with the implementation above if rcut is the largest cutoff or the cutoff is so small that the only neighbor
!   the atoms see are themselves (n_neigh(i) = 1 for all i).

    if( do_timing) then
      call cpu_time(time1)
    end if

!   Allocate all the necessary stuff
    allocate( neighbor_c6_ii(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( f_damp(1:n_pairs) )
    allocate( T_func(1:9*n_pairs) )
    allocate( h_func(1:9*n_pairs,1:n_freq) )
    allocate( g_func(1:n_pairs,1:n_freq) )
    allocate( omegas(1:n_freq) )
    allocate( omega_i(1:n_pairs) )
    allocate( sigma_i(1:n_sites,1:n_freq) )
    allocate( alpha_i(1:n_sites,1:n_freq) )
    allocate( sigma_ij(1:n_freq) )
    allocate( T_SR(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( B_mat(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( xyz_H(1:3,1:n_pairs) )
    allocate( rjs_H(1:n_pairs) )
    allocate( A_mat(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( work_arr(1:12*n_sites) )
    allocate( ipiv(1:3*n_sites) )
    allocate( A_i(1:3,1:3) )
    allocate( alpha_SCS(1:n_sites,1:n_freq) )
    if ( do_derivatives ) then
      allocate( dT(1:9*n_pairs,1:3) )
      allocate( dT_SR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( f_damp_der(1:n_pairs,1:3) )
      allocate( g_func_der(1:n_pairs,1:3,1:n_freq) )
      allocate( h_func_der(1:9*n_pairs,1:3,1:n_freq) )
      allocate( dA_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( coeff_h_der(1:n_freq) )
      allocate( terms(1:n_freq) )
      allocate( dT_SR_A_mat(1:3*n_sites,1:3*n_sites) )
      allocate( dT_SR_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( dB_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( dB_mat_v(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( coeff_der(1:n_sites,1:n_sites,1:3,1:3,1:n_freq) )
      allocate( coeff_fdamp(1:n_sites,1:n_sites,1:3,1:3,1:n_freq) )
      allocate( dg(1:n_freq) )
      allocate( dh(1:n_freq) )
      allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
    end if
    allocate( alpha_k(1:n_pairs,1:n_freq) )
    allocate( sigma_k(1:n_pairs,1:n_freq) )
    allocate( s_i(1:n_freq) )
    allocate( s_j(1:n_freq) )

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "Timing for allocation:", time2-time1
      call cpu_time(time1)
    end if
  
!   Frequencies used for integration:
    omega = 0.d0
    do i = 1, n_freq
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
      if ( do_derivatives ) then
        hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,k)*Bohr
      end if
    end do

!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    omega_i = (4.d0 * neighbor_c6_ii)/(3.d0*neighbor_alpha0**2)

    k2 = 1
    do i = 1, n_sites
      do k = 1, n_freq
        alpha_i(i,k) = neighbor_alpha0(k2)/(1.d0 + omegas(k)**2/omega_i(k2)**2)
        sigma_i(i,k) = (sqrt(2.d0/pi) * alpha_i(i,k)/3.d0)**(1.d0/3.d0)
      end do
      k2 = k2+n_neigh(i)
    end do

    do k = 1, n_freq
      alpha_k(:,k) = neighbor_alpha0/(1.d0 + omegas(k)**2/omega_i**2)
      sigma_k(:,k) = (sqrt(2.d0/pi) * alpha_k(:,k)/3.d0)**(1.d0/3.d0)
    end do

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "Timing for initialization:", time2-time1
      call cpu_time(time1)
    end if

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

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "Timing for T_func, g_func, h_func:", time2-time1
      call cpu_time(time1)
    end if

    allocate( ia(1:9*n_sites*n_sites) )
    allocate( ja(1:9*n_sites*n_sites) )
    allocate( val(1:9*n_sites*n_sites) )

    T_SR = 0.d0
    B_mat = 0.d0
    k = 0
    nnz = 0
    do i = 1, n_sites
      k = k+1
      do c1 = 1, 3
        B_mat(3*(i-1)+c1,3*(i-1)+c1,:) = 1.d0
        nnz = nnz+1
        ia(nnz) = 3*(i-1)+c1
        ja(nnz) = 3*(i-1)+c1
        val(nnz) = B_mat(3*(i-1)+c1,3*(i-1)+c1,1)
      end do
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        j = modulo(j-1,n_sites)+1
        if( rjs(k) < rcut )then
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2+1
              nnz = nnz+1
              do k3 = 1, n_freq
                T_SR(3*(i-1)+c1,3*(j-1)+c2,k3) = alpha_i(i,k3) * (1.d0-f_damp(k)) * (-T_func(k2) * &
                                                g_func(k,k3) + h_func(k2,k3))
              end do
              ia(nnz) = 3*(i-1)+c1
              ja(nnz) = 3*(j-1)+c2
              val(nnz) = T_SR(3*(i-1)+c1,3*(j-1)+c2,1)
            end do
          end do
        end if
      end do
    end do

    B_mat = B_mat + T_SR

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "Timing for B_mat:", time2-time1
      call cpu_time(time1)
    end if

    if ( do_derivatives ) then
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

      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Timing for dT:", time2-time1
        call cpu_time(time1)
      end if

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

      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Timing for dT_SR (no hirshfeld):", time2-time1
        call cpu_time(time1)
      end if

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
              do k = 1, n_freq
                do j2 = 2, n_neigh(i)
                  if (rjs(k3+j2) < rcut) then
                    j = neighbors_list(k3+j2)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                          (coeff_der(i,j,c1,c2,k) * sigma_i(i,k)**2/hirshfeld_v(i) + &
                          coeff_fdamp(i,j,c1,c2,k) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,k2)
                        dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) = dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) + &
                          (coeff_der(j,i,c1,c2,k) * sigma_i(i,k)**2/hirshfeld_v(i) + &
                          coeff_fdamp(j,i,c1,c2,k) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,k2)
! TEST:
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                          coeff_der(i,j,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
!                          dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                          coeff_fdamp(i,j,c1,c2,k) * r_vdw_i/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
! TEST:
!                        dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) = dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) + &
!                          coeff_der(j,i,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
!                        dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) = &
!                          dT_SR_v(3*(j-1)+c1,3*(i-1)+c2,a,c3,k) + &
!                          coeff_fdamp(j,i,c1,c2,k) * r_vdw_i/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
! END TEST.
                      end do
                    end do
                  end if      
                end do
              end do
            end do
          end do
          k3 = k3 + n_neigh(i)
        end do

!        k2 = 0
!        k3 = 0
!        do j = 1, n_sites
!          r_vdw_j = r0_ii(k3+1)
!          do a2 = 1, n_neigh(j)
!            k2 = k2+1
!            a = neighbors_list(k2)
!            do c3 = 1, 3
!              do k = 1, n_freq
!                do i2 = 2, n_neigh(j)
!                  if (rjs(k3+i2) < rcut) then
!                    i = neighbors_list(k3+i2)
!                    do c1 = 1, 3
!                      do c2 = 1, 3
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                          coeff_der(i,j,c1,c2,k) * (sigma_i(j,k)**2/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2))
!                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
!                          dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
!                          coeff_fdamp(i,j,c1,c2,k) * r_vdw_j/hirshfeld_v(j) * hirshfeld_v_cart_der_H(c3,k2)
!                      end do
!                    end do
!                  end if
!                end do
!              end do!
!            end do
!          end do
!          k3 = k3 + n_neigh(j)
!        end do

        k2 = 0
        do i = 1, n_sites
          do a2 = 1, n_neigh(i)
            k2 = k2+1
            a = neighbors_list(k2)
            do k = 1, n_freq
              do c3 = 1, 3
                dT_SR_v(3*(i-1)+1:3*(i-1)+3,:,a,c3,k) = dT_SR_v(3*(i-1)+1:3*(i-1)+3,:,a,c3,k) + &
                  1.d0/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2) * &
                  T_SR(3*(i-1)+1:3*(i-1)+3,:,k) / alpha_i(i,k)
              end do
            end do
          end do
        end do

        dT_SR = dT_SR + dT_SR_v
      end if

      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Timing for dT_SR (hirshfeld):", time2-time1
        call cpu_time(time1)
      end if

      do i = 1, n_sites
        do k = 1, n_freq
          dT_SR(3*(i-1)+1:3*(i-1)+3,:,:,:,k) = alpha_i(i,k) * dT_SR(3*(i-1)+1:3*(i-1)+3,:,:,:,k) 
        end do
      end do

      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Timing for dT_SR (alpha multiplication):", time2-time1
        call cpu_time(time1)
      end if

    end if ! do_derivatives

    !if ( iterative ) then

      !call psb_init(icontxt)
      !call psb_cdall(icontxt, desc_a, info_psb, nl=3*n_sites)
      !write(*,*) "cdall", info_psb
      !call psb_spall(A_sp, desc_a, info_psb, nnz)
      !write(*,*) "spall", info_psb
      !call psb_geall(x_vec, desc_a, info_psb, 3, 1)
      !write(*,*) "geall x", info_psb, size(x_vec,1), size(x_vec,2)
      !call psb_geall(b_vec, desc_a, info_psb, 3, 1)
      !write(*,*) "geall b", info_psb, size(b_vec,1), size(b_vec,2)
      !write(*,*) "size of b", size(b_vec,1)
      !call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
      !write(*,*) "spins", info_psb
      !allocate( val_xv(1:3*n_sites,1:3) )
      !allocate( val_bv(1:3*n_sites,1:3) )
      !allocate( myidx(1:3*n_sites) )
      !val_xv = 0.d0
      !val_bv = 0.d0
      !do p = 1, n_sites
      !  do c1 = 1, 3
      !    myidx(3*(p-1)+c1) = 3*(p-1)+c1
      !    val_xv(3*(p-1)+c1,c1) = 1.d0
      !    val_bv(3*(p-1)+c1,c1) = 1.d0
      !  end do
      !end do
      !call psb_geins(3*n_sites, myidx, val_xv, x_vec, desc_a, info_psb)
      !write(*,*) "x_vec", info_psb
      !call psb_geins(3*n_sites, myidx, val_bv, b_vec, desc_a, info_psb)
      !write(*,*) "b_vec", info_psb
      !call psb_cdasb(desc_a, info_psb)
      !write(*,*) "cdasb", info_psb
      !call psb_spasb(A_sp, desc_a, info_psb)
      !write(*,*) "spasb", info_psb
      !call psb_geasb(x_vec, desc_a, info_psb)
      !write(*,*) "geasb x", info_psb
      !call psb_geasb(b_vec, desc_a, info_psb)
      !write(*,*) "geasb b", info_psb
      !ptype="DIAG"
      !call prec%init(icontxt, ptype, info_psb)
      !call prec%build(A_sp, desc_a, info_psb)
      !write(*,*) "prec build", info_psb
      !call psb_krylov("BICGSTAB", A_sp, prec, b_vec, x_vec, 1.d-6, desc_a, info)
      !deallocate( val_xv, val_bv )
      !call psb_exit(icontxt)
      deallocate( ia, ja, val )

      !n_iter = 100
      allocate( I_mat(1:3*n_sites,1:3*n_sites) )
      allocate( a_vec(1:3*n_sites,1:3,1:n_freq) )
      I_mat = 0.d0
      do p = 1, n_sites
        do c1 = 1, 3
          I_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0
        end do
      end do
      a_vec = 0.d0
      do k = 1, n_freq
        do p = 1, n_sites
          do c1 = 1, 3
            a_vec(3*(p-1)+c1,c1,k) = alpha_i(p,k)
          end do
        end do
      end do

      allocate( a_SCS(1:3*n_sites,1:3,1:n_freq) )
      allocate( BTB_reg_copy(1:3*n_sites,1:3*n_sites) )

      reg_param = 0.01d0
      allocate( BTB_reg(1:3*n_sites,1:3*n_sites,1:n_freq) )
      allocate( B_reg(1:3*n_sites,1:3*n_sites,1:n_freq) )
      do k = 1, n_freq
        B_reg(:,:,k) = B_mat(:,:,k) + reg_param * I_mat
        call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_reg(:,:,k), 3*n_sites, &
                    a_vec(:,:,k), 3*n_sites, 0.d0, a_SCS(:,:,k), 3*n_sites)
        call dgemm('t', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, B_mat(:,:,k), 3*n_sites, &
                    B_mat(:,:,k), 3*n_sites, 0.d0, BTB_reg(:,:,k), 3*n_sites)
        BTB_reg(:,:,k) = BTB_reg(:,:,k) + reg_param * I_mat
        BTB_reg_copy = BTB_reg(:,:,k)
        call dsysv('U', 3*n_sites, 3, BTB_reg_copy, 3*n_sites, ipiv, &
                    a_SCS(:,:,k), 3*n_sites, work_arr, 12*n_sites, info)
      end do

      alpha_SCS0 = 0.d0
      do k = 1, n_freq
        do p = 1, n_sites
          do c1 = 1, 3
            alpha_SCS0(p,k) = alpha_SCS0(p,k) + a_SCS(3*(p-1)+c1,c1,k)
          end do
        end do
      end do
      alpha_SCS0 = alpha_SCS0/3.d0

      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Timing for alpha_SCS:", time2-time1
        call cpu_time(time1)
      end if

      if ( do_derivatives ) then

        allocate( da_vec(1:3*n_sites,1:3,1*n_sites,1:3,1:n_freq) )
        da_vec = 0.d0
        do k = 1, n_freq
          k2 = 0
          do i = 1, n_sites
            do a2 = 1, n_neigh(i)
              k2 = k2+1
              a = neighbors_list(k2)
              do c3 = 1, 3
                do c1 = 1, 3
                  da_vec(3*(i-1)+c1,c1,a,c3,k) = 1.d0/hirshfeld_v(i) * alpha_i(i,k) * &
                    hirshfeld_v_cart_der_H(c3,k2)
                end do
              end do
            end do
          end do
        end do

        if( do_timing) then
          call cpu_time(time2)
          write(*,*) "Timing for da_vec:", time2-time1
          call cpu_time(time1)
        end if

        allocate( dBTB(1:3*n_sites,1:3*n_sites) )
        allocate( vect1(1:3*n_sites,1:3) )
        allocate( vect2(1:3*n_sites,1:3) )
        allocate( vect3(1:3*n_sites,1:3) )
        allocate( vect4(1:3*n_sites,1:3) )
        allocate( vect_temp(1:3*n_sites,1:3) )
        allocate( da_SCS(1:3*n_sites,1:3,1:n_sites,1:3,1:n_freq) )

        if( do_timing) then
          call cpu_time(time2)
          write(*,*) "Timing for alpha_SCS_grad allocation:", time2-time1
          call cpu_time(time1)
        end if

      ! What happens here: regularized equation
      ! (B^T * B + reg_param * I) a_SCS = (B^T + reg_param*I) a_vec
      ! is differentiated:
      ! (B^T * B + reg_param * I) a_SCS' = (B')^T a_vec + (B^T + reg_param*I) a_vec' - ((B')^T * B) a_SCS - (B^T * B') a_SCS
      !                                    |---vect1--|   |---------vect2----------|   |-----vect3------|   |----vect4-----|
      ! a_SCS' is solved with dsysv.

!        write(*,*) "dT_SR"
!        write(*,*) dT_SR(1,:,1,1,1)
!        write(*,*) "every third"
!        write(*,*) dT_SR(1,1::3,1,1,1)
!        write(*,*) "a_vec"
!        write(*,*) a_vec(:,1,1)
!        write(*,*) "every third"
!        write(*,*) a_vec(1::3,1,1)

        do k = 1, n_freq
          do a = 1, n_sites
            do c3 = 1, 3
!TEST
              do i = 1, 3*n_sites
                do c1 = 1, 3
                  vect1(i,c1) = dot_product(dT_SR(c1::3,i,a,c3,k), a_vec(c1::3,c1,k))
                  vect2(i,c1) = dot_product(B_reg(c1::3,i,k), da_vec(c1::3,c1,a,c3,k))
                end do
              end do
!              write(*,*) "vect1"
!              do i = 1, 6
!                write(*,*) vect1(i,:)
!              end do
!              call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, &
!                a_vec(:,:,k), 3*n_sites, 0.d0, vect1, 3*n_sites)
!              call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_reg(:,:,k), 3*n_sites, &
!                da_vec(:,:,a,c3,k), 3*n_sites, 0.d0, vect2, 3*n_sites)
!              write(*,*) "vect1 after"
!              do i = 1, 6
!                write(*,*) vect1(i,:)
!              end do
!TEST ENDS
              call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_mat(:,:,k), 3*n_sites, &
                a_SCS(:,:,k), 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, &
                vect_temp, 3*n_sites, 0.d0, vect3, 3*n_sites)
              call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, &
                a_SCS(:,:,k), 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_mat(:,:,k), 3*n_sites, &
                vect_temp, 3*n_sites, 0.d0, vect4, 3*n_sites)
              da_SCS(:,:,a,c3,k) = vect1 + vect2 - vect3 - vect4
              BTB_reg_copy = BTB_reg(:,:,k)
              call dsysv('U', 3*n_sites, 3, BTB_reg_copy, 3*n_sites, ipiv, &
                da_SCS(:,:,a,c3,k), 3*n_sites, work_arr, 12*n_sites, info)
            end do
          end do
        end do

        if( do_timing) then
          call cpu_time(time2)
          write(*,*) "Timing for da_SCS:", time2-time1
          call cpu_time(time1)
        end if

        do k = 1, n_freq
          do i = 1, n_sites
            do a = 1, n_sites
              do c3 = 1, 3
                dalpha_full(i,a,c3,k) = 1.d0/3.d0 * (da_SCS(3*(i-1)+1,1,a,c3,k) + &
                  da_SCS(3*(i-1)+2,2,a,c3,k) + da_SCS(3*(i-1)+3,3,a,c3,k))
              end do
            end do
          end do
        end do

        deallocate( da_vec, dBTB, vect1, vect2, vect3, vect4, vect_temp, da_SCS )

        if( do_timing) then
          call cpu_time(time2)
          write(*,*) "Timing for alpha_SCS_grad:", time2-time1
          call cpu_time(time1)
        end if


      end if

      deallocate( BTB_reg, B_reg, BTB_reg_copy, I_mat, a_vec, a_SCS )

      write(*,*) "alpha_SCS:" !,  alpha_SCS0(1,1)
      do p = 1, n_sites
        write(*,*) p, alpha_SCS0(p,1)
      end do


      if ( do_derivatives ) then
        write(*,*) "dalpha_SCS w.r.t. atom 1 in x direction:"
        do p = 1, n_sites
          write(*,*) p, dalpha_full(p,1,1,1)
        end do
      end if

!   Clean up
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, alpha_k, sigma_k, s_i, s_j, )
    if ( do_derivatives ) then
      deallocate( dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                  dB_mat, dB_mat_v, coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H )
    
    end if

  end subroutine

!**************************************************************************





!**************************************************************************

!NOTE: T_func and dT remain unchanged. Should they be directly passed to this subroutine?

  subroutine get_mbd_energies_and_forces( alpha_SCS0, alpha_SCS_grad, n_neigh, neighbors_list, neighbor_species, &
                                          rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, &
                                          sR, d, c6_ref, r0_ref, alpha0_ref, do_derivatives, &
                                          energies, forces0, virial)

    implicit none
    real*8, intent(inout) :: alpha_SCS0(:,:), alpha_SCS_grad(:,:,:,:), energies(:), forces0(:,:)
    real*8, intent(in) :: rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:)
    real*8, intent(out) :: virial(1:3, 1:3)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_derivatives
    real*8, allocatable :: T_LR(:,:), r0_ii_SCS(:), r0_ii(:), neighbor_alpha0(:), &
                           xyz_H(:,:), rjs_H(:), f_damp_SCS(:), T_func(:), AT(:,:,:), AT_n(:,:,:,:), &
                           energy_series(:,:,:), integrand(:,:), omegas(:), f_damp_der(:,:), &
                           f_damp_der_SCS(:,:,:), dT_LR(:,:,:,:), dT(:,:), &
                           G_mat(:,:,:,:,:), force_integrand(:,:,:), force_series(:,:,:), &
                           VL(:,:)
    integer :: n_order, n_freq, n_sites, n_pairs, n_species, n_sites0
    integer :: k, k2, i, j, j2, c1, c2, c3, a
    real*8 :: Bohr, Hartree, pi, r_vdw_i, r_vdw_j, E_MBD, integral, omega, R_vdW_SCS_ij, S_vdW_ij, &
              dS_vdW_ij, time1, time2
    logical :: series_average = .true., do_timing = .true.

    write(*,*) "rcut (MBD)", rcut
!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

    n_order = 3

    n_freq = size(alpha_SCS0, 2)
    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

    call cpu_time(time1)

    allocate( omegas(1:n_freq) )
    allocate( xyz_H(1:3,1:n_pairs) )
    allocate( rjs_H(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( T_func(1:9*n_pairs) )
    allocate( T_LR(1:3*n_sites,1:3*n_sites) )
    allocate( r0_ii_SCS(1:n_pairs) )
    allocate( f_damp_SCS(1:n_pairs) )
    allocate( AT(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( AT_n(1:3*n_sites,1:3*n_sites, 1:n_order-1, 1:n_freq) )
    allocate( energy_series(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( integrand(1:n_sites,1:n_freq) )
    if ( do_derivatives ) then
      allocate( dT(1:9*n_pairs,1:3) )
      allocate( f_damp_der(1:n_pairs,1:3) )
      allocate( f_damp_der_SCS(1:n_pairs,1:n_sites,1:3) )
      allocate( dT_LR(1:3*n_sites,1:3*n_sites,1:n_sites,1:3) )
      allocate( G_mat(1:3*n_sites,1:3*n_sites,1:n_sites,1:3,1:n_freq) )
      allocate( force_integrand(1:n_sites,1:3,1:n_freq) ) 
      allocate( force_series(1:3*n_sites,1:3*n_sites,1:n_freq) )
      allocate( VL(1:3*n_sites,1:3*n_sites) )
    end if

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for allocation:", time2-time1
      call cpu_time(time1)
    end if

    omega = 0.d0
    do i = 1, n_freq
      omegas(i) = omega
      omega = omega + 0.4d0
    end do

    do k = 1, n_pairs
      j = neighbor_species(k)
      r0_ii(k) = r0_ref(j) / Bohr
      neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
      xyz_H(:,k) = xyz(:,k)/Bohr
      rjs_H(k) = rjs(k)/Bohr
    end do

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for initialization:", time2-time1
      call cpu_time(time1)
    end if

    T_func = 0.d0
    k = 0
    do i = 1, n_sites
      k = k+1
      do j2 = 2, n_neigh(i)
        k = k+1
        if( rjs(k) < rcut )then
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

    do k = 1, n_pairs
      j = neighbors_list(k)
      r0_ii_SCS(k) = r0_ii(k) * (alpha_SCS0(j,1)/neighbor_alpha0(k))**(1.d0/3.d0)
    end do

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for T_func, r0_ii_SCS:", time2-time1
      call cpu_time(time1)
    end if

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

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for f_damp_SCS:", time2-time1
      call cpu_time(time1)
    end if

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
            end do
          end do
        end if
      end do
    end do

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for T_LR:", time2-time1
      call cpu_time(time1)
    end if

    if ( do_derivatives ) then

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
    
      if( do_timing ) then
        call cpu_time(time2)
        write(*,*) "Timing for dT:", time2-time1
        call cpu_time(time1)
      end if

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
                dS_vdW_ij = sR/3.d0 * ( r_vdw_i/alpha_SCS0(i,1) * alpha_SCS_grad(i,a,c3,1) + &
                                        r_vdw_j/alpha_SCS0(j,1) * alpha_SCS_grad(j,a,c3,1) )
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

      if( do_timing ) then
        call cpu_time(time2)
        write(*,*) "Timing for dT_LR:", time2-time1
        call cpu_time(time1)
      end if

    end if

    AT = 0.d0
    do k = 1, n_freq
      do i = 1, n_sites
        AT(3*(i-1)+1:3*(i-1)+3,:,k) = alpha_SCS0(i,k) * T_LR(3*(i-1)+1:3*(i-1)+3,:)
      end do
    end do

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for AT:", time2-time1
      call cpu_time(time1)
    end if

    AT_n = 0.d0
    integrand = 0.d0
    energy_series = 0.d0

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for initialization of AT_n, integrand, energy_series:", time2-time1
      call cpu_time(time1)
    end if

    if ( n_order > 1 ) then
      do k = 1, n_freq
        do k2 = 1, n_order-1
          ! Precalculate the full AT_n for forces:
          if ( k2 == 1 ) then
            AT_n(:,:,k2,k) = AT(:,:,k)
          else
            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, AT(:,:,k), 3*n_sites, &
                       AT_n(:,:,k2-1,k), 3*n_sites, 0.d0, AT_n(:,:,k2,k), 3*n_sites)
          end if
          if (series_average .and. k2 == n_order-1) then
            energy_series(:,:,k) = energy_series(:,:,k) - 1.d0/(k2+1)*AT_n(:,:,k2,k)/2.d0
          else
            energy_series(:,:,k) = energy_series(:,:,k) - 1.d0/(k2+1)*AT_n(:,:,k2,k)
          end if
        end do
        do i = 1, n_sites
          do c1 = 1, 3
            integrand(i,k) = integrand(i,k) + alpha_SCS0(i,k) * dot_product(T_LR(3*(i-1)+c1,:), &
                           energy_series(:,3*(i-1)+c1,k))
          end do
        end do
      end do
    else
      write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
      write(*,*) "and forces will be zero."
    end if

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for energy integrand:", time2-time1
      call cpu_time(time1)
    end if

    do i = 1, n_sites
      integral = 0.d0
      call integrate("trapezoidal", omegas, integrand(i,:), omegas(1), omegas(n_freq), integral)
      E_MBD = integral/(2.d0*pi)
      energies(i) = E_MBD * 27.211386245988
      write(*,*) i, energies(i)
    end do
    write(*,*) "MBD energy:", sum(energies)

    if( do_timing ) then
      call cpu_time(time2)
      write(*,*) "Timing for MBD energy integration:", time2-time1
      call cpu_time(time1)
    end if

    if ( do_derivatives ) then

      G_mat = 0.d0
      force_integrand = 0.d0
      force_series = 0.d0

      if( do_timing ) then
        call cpu_time(time2)
        write(*,*) "Timing for initialization of G_mat, force_integrand, force_series:", time2-time1
        call cpu_time(time1)
      end if

      if ( n_order > 1 ) then
        do k = 1, n_freq
          if( do_timing ) then
            call cpu_time(time2)
            write(*,*) "force_series start:", time2-time1
            call cpu_time(time1)
          end if
          do k2 = 1, n_order-1
            if (series_average .and. k2 == n_order-1) then
              force_series(:,:,k) = force_series(:,:,k) + AT_n(:,:,k2,k)/2.d0
            else
              force_series(:,:,k) = force_series(:,:,k) + AT_n(:,:,k2,k)
            end if
          end do
          if( do_timing ) then
            call cpu_time(time2)
            write(*,*) "Timing for force_series freq:", k, time2-time1
            call cpu_time(time1)
          end if
          do a = 1, n_sites
            do c3 = 1, 3
              ! call dgemm for G_mat
              do i = 1, n_sites
                G_mat(3*(i-1)+1:3*(i-1)+3,:,a,c3,k) = G_mat(3*(i-1)+1:3*(i-1)+3,:,a,c3,k) + &
                  alpha_SCS0(i,k) * dT_LR(3*(i-1)+1:3*(i-1)+3,:,a,c3) + &
                  alpha_SCS_grad(i,a,c3,k) * T_LR(3*(i-1)+1:3*(i-1)+3,:)
              end do
              !call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, force_series(:,:,k), 3*n_sites, &
              !           G_mat(:,:,a,c3,k), 3*n_sites, 0.d0, VL, 3*n_sites)
              do i = 1, 3*n_sites
                !force_integrand(a,c3,k) = force_integrand(a,c3,k) + VL(i,i)
                force_integrand(a,c3,k) = force_integrand(a,c3,k) + &
                  dot_product(G_mat(i,:,a,c3,k), force_series(:,i,k))
              end do
            end do
          end do
          if( do_timing ) then
            call cpu_time(time2)
            write(*,*) "Timing for force_integrand freq:", k, time2-time1
            call cpu_time(time1)
          end if
        end do
      end if

      if( do_timing ) then
        call cpu_time(time2)
        write(*,*) "Timing for force_integrand:", time2-time1
        call cpu_time(time1)
      end if

      do a = 1, n_sites
        do c3 = 1, 3
          integral = 0.d0
          call integrate("trapezoidal", omegas, force_integrand(a,c3,:), omegas(1), omegas(n_freq), integral)
          forces0(c3,a) = 1.d0/(2.d0*pi) * integral
        end do
      end do
      forces0 = forces0 * 51.42208619083232

      if( do_timing ) then
        call cpu_time(time2)
        write(*,*) "Timing for MBD force integration:", time2-time1
        call cpu_time(time1)
      end if

      write(*,*) "MBD forces:"
      do i = 1, n_sites
        write(*,*) i, forces0(:,i)
      end do

    end if

    deallocate( omegas, xyz_H, rjs_H, r0_ii, neighbor_alpha0, T_func, T_LR, r0_ii_SCS, &
                f_damp_SCS, AT, AT_n, energy_series, integrand )
    if ( do_derivatives ) then
      deallocate( dT, dT_LR, f_damp_der, f_damp_der_SCS, force_series, force_integrand, G_mat, VL )
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
