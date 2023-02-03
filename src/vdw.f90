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
  !use psb_base_mod
  !use psb_prec_mod
  !use psb_krylov_mod
  !use psb_util_mod

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
                                       sR, d, c6_ref, r0_ref, alpha0_ref, c6_scs, r0_scs, alpha0_scs, do_forces, &
                                       energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v(:), hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), hirshfeld_v_neigh(:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:), c6_scs(:), r0_scs(:), alpha0_scs(:)
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
    logical :: do_timing = .false., read_hirshfeld = .false.

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
    if ( read_hirshfeld ) then
      neighbor_c6_ii = c6_scs
      r0_ii = r0_scs
      neighbor_alpha0 = alpha0_scs
    else
      neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!     This is slow, could replace by Taylor expansion maybe
      r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
      neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    end if

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
    energies = 0.d0
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

    write(*,*) "TS energy", sum(energies)

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
  subroutine get_scs_polarizabilities( hirshfeld_v, hirshfeld_v_cart_der, hirshfeld_v_cart_der_ji, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, rcut_mbd, rcut_2b, r_buffer, rcut_inner, buffer_inner, rjs, xyz, &
                                       hirshfeld_v_neigh, sR, d, c6_ref, r0_ref, alpha0_ref, do_derivatives, &
                                       do_hirshfeld_gradients, polynomial_expansion, n_freq, n_order, &
                                       central_pol, central_omega, dalpha_full, &
                                       c6_scs, r0_scs, alpha0_scs, energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v_cart_der(:,:), rcut, r_buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), rcut_mbd, rcut_2b, hirshfeld_v_cart_der_ji(:,:), &
                          alpha0_ref(:) !, hirshfeld_v(:), hirshfeld_v_neigh(:) !NOTE: uncomment this in final implementation
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:), n_freq, n_order
    logical, intent(in) :: do_derivatives, do_hirshfeld_gradients, polynomial_expansion
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:), central_pol(:), central_omega(:), dalpha_full(:,:), hirshfeld_v(:), &
                             hirshfeld_v_neigh(:), c6_scs(:), r0_scs(:), alpha0_scs(:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), f_damp(:), neighbor_alpha0(:), T_func(:), h_func(:), g_func(:), &
                           omegas(:), omega_i(:), B_mat(:,:), rjs_H(:), xyz_H(:,:), work_arr(:), dT(:), f_damp_der(:), &
                           g_func_der(:), h_func_der(:), coeff_der(:), coeff_fdamp(:), hirshfeld_v_cart_der_H(:,:), &
                           a_SCS(:,:), da_SCS(:,:), b_der(:,:), alpha_SCS_full(:,:,:), dB_mat(:,:), alpha_test(:,:), &
                           neighbor_sigma(:), hirshfeld_v_sub_der(:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, f_damp_SCS_ij, t1, t2, &
              sigma_ij, coeff_h_der, dg, dh, s_i, s_j, terms, omega_ref, xyz_i(1:3), xyz_j(1:3)
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, om, n_tot
    integer :: i, i0, i1, i2, i3, j, j1, j2, j3, k, k2, k3, k4, a, a2, c1, c2, c3, lwork, b, p, q, r, n_count, n_max, k_i, k_j
    logical :: read_hirshfeld = .false.
               
!    LOCAL TEST stuff:
    integer, allocatable :: sub_neighbors_list(:), n_sub_neigh(:), p_list(:)
    integer :: n_sub_sites, n_sub_pairs, s

!   MBD stuff:
    real*8, allocatable :: T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), AT(:,:,:), AT_n(:,:,:,:), energy_series(:,:), &
                           integrand(:), AT_n_f(:,:,:,:), f_damp_der_SCS(:), dT_LR(:,:), G_mat(:,:,:), force_series(:,:), &
                           total_energy_series(:,:), total_integrand(:), rjs_0(:), &
                           T_mbd(:), r0_ii_mbd(:), neighbor_alpha0_mbd(:), omegas_mbd(:), rjs_0_mbd(:), xyz_0_mbd(:,:), &
                           xyz_mbd(:,:), rjs_mbd(:), d_der(:,:), dT_mbd(:), f_damp_der_mbd(:), a_mbd(:), da_mbd(:), &
                           a_iso(:,:), o_p(:), da_iso(:,:,:),  o_mbd(:), sub_2b_list(:), &
                           xyz_2b(:,:), rjs_2b(:), r0_ii_2b(:), neighbor_alpha0_2b(:), f_damp_SCS_2b(:), &
                           a_2b(:), r0_ii_SCS_2b(:), C6_2b(:), da_2b(:), T_SR(:), T_SR_mult(:), d_arr_i(:), d_arr_o(:), &
                           d_mult_i(:), d_mult_o(:), dT_SR_mult(:,:), d_dmult_i(:,:), d_dmult_o(:,:), do_mbd(:), &
                           hirshfeld_sub_neigh(:), o_2b(:), do_2b(:)
    real*8 :: a_mbd_i, a_mbd_j, da_i, da_j, pol1, E_TS, f_damp_der_2b, dr_vdw_i, &
              dr_vdw_j, forces_TS, dC6_2b, mult1_i, mult1_j, mult2, dmult1_i(1:3), dmult1_j(1:3), dmult2(1:3), hv_p_der, &
              hv_q_der, do_pref
    integer :: n_mbd_sites, n_mbd_pairs, n_2b_sites
    integer, allocatable :: n_mbd_neigh(:), mbd_neighbors_list(:), p_mbd(:)


    !PSBLAS stuff:
    !type(psb_ctxt_type) :: icontxt
    !integer(psb_ipk_) ::  iam, np, ip, jp, idummy, nr, nnz, info_psb
    !type(psb_desc_type) :: desc_a
    !type(psb_dspmat_type) :: A_sp
    !type(psb_d_vect_type) :: x_vec, b_vec
    !integer(psb_lpk_), allocatable :: ia(:), ja(:), myidx(:)
    !real(psb_dpk_), allocatable :: val(:), val_xv(:,:), b_i(:,:), d_vec(:,:)
    !type(psb_dprec_type) :: prec
    !character(len=20) :: ptype
    real*8 :: polyfit(1:15)
    integer :: n_degree
    real*8, allocatable :: B_pol(:,:), B_mult(:,:), b_i(:,:), d_vec(:,:), val_xv(:,:)

!   IMPORTANT NOTE ABOUT THE DERIVATIVES:
!   If rcut < rcut_soap, the derivatives in the new implementation omit the terms that fall outside of rcut.
!   This means that the finite difference and the analytical derivative do not match in the present version
!   in this special case. The old implementation gives the correct derivatives even in this case. I will
!   probably fix this at some point but it will require using the full n_neigh(i) again, instead of 
!   n_ssites, to construct the sneighbors_list.

    !write(*,*) "rcut (SCS)", rcut

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

    !open(unit=89, file="hirshfeld.dat", status="new")
    !do i = 1, n_sites
    !  write(89,*) hirshfeld_v(i)
    !end do
    !close(89)


    ! HACK FOR HIRSHFELD DERIVATIVES
    allocate( hirshfeld_v_cart_der_H(1:3, n_pairs) )

    hirshfeld_v_cart_der_H = 0.d0

    !write(*,*) "dv"
    !k = 0
    !do i = 1, n_sites
    !  do j2 = 1, n_neigh(i)
    !    k = k+1
    !    j = neighbors_list(k)
    !    n_tot = sum(n_neigh(1:j)) - n_neigh(j)
    !    a = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(j)),i,1)
    !    hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,n_tot+a)
    !  end do
    !end do

    !do i = 1, n_neigh(1)
    !  write(*,*) hirshfeld_v_cart_der_H(1:3,i), hirshfeld_v_cart_der_ji(1:3,i)
    !end do
    if ( do_derivatives .and. do_hirshfeld_gradients ) then
    hirshfeld_v_cart_der_H = hirshfeld_v_cart_der_ji
    !write(*,*) "hirshfeld_der"
    !write(*,*) hirshfeld_v_cart_der_H(1:n_neigh(1))
    end if

    ! HACK END


    if ( read_hirshfeld ) then

      open(unit=89, file="hirshfeld.dat", status="old")
      do i = 1, n_sites
        read(89,*) hirshfeld_v(i)
      end do
      close(89)

      hirshfeld_v_neigh = 0.d0
      k = 0
      do i = 1, n_sites
        do j = 1, n_neigh(i)
          k = k+1
          i2 = modulo(neighbors_list(k)-1, n_sites) + 1 !neighbors_list(k)
          hirshfeld_v_neigh(k) = hirshfeld_v(i2)
        end do
      end do

      !write(*,*) "hirshfeld_v", hirshfeld_v
      !write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh 

    end if

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

    write(*,*) "rcut", rcut/Bohr

!   Number of frequencies
    !n_freq = size(alpha_SCS0, 2)
      
    write(*,*) "Starting local calculation"

    E_MBD = 0.d0

    !alpha_SCS0 = 0.d0
    !allocate( alpha_SCS_full(1:3*n_sites,1:3,1:2) )
    if ( do_derivatives ) then
      dalpha_full = 0.d0
    end if

    !alpha_SCS_full = 0.d0
    
    !allocate( central_pol(1:n_sites) )
    !central_pol = 0.d0
    !allocate( central_omega(1:n_sites) )
    !central_omega = 0.d0
      
    !r_buffer = 0.5d0

    call cpu_time(time1)

    do i = 1, n_sites      

      n_tot = sum(n_neigh(1:i))-n_neigh(i)
      i0 = modulo(neighbors_list(n_tot+1)-1, n_sites0) + 1
      s = neighbor_species(n_tot+1)
      omega_ref = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
      n_sub_sites = 0
      n_sub_pairs = 0
      k_i = 0
      do i3 = 1, n_neigh(i)
        k_i = k_i + 1
        if (rjs(n_tot+k_i) < rcut ) then
          n_sub_sites = n_sub_sites + 1
          n_sub_pairs = n_sub_pairs + 1
          xyz_i = xyz(:,n_tot+k_i)/Bohr
          k_j = 0
          do j3 = 1, n_neigh(i)
            k_j = k_j + 1
            if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
              if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                  n_sub_pairs = n_sub_pairs + 1
                end if
              end if
            end if
          end do
        end if
      end do
        
      allocate( sub_neighbors_list(1:n_sub_pairs) )
      allocate( n_sub_neigh(1:n_sub_sites) )
      allocate( p_list(1:n_sub_pairs) )
      allocate( xyz_H(1:3,1:n_sub_pairs) )
      allocate( rjs_H(1:n_sub_pairs) )
      allocate( r0_ii(1:n_sub_pairs) )
      allocate( neighbor_alpha0(1:n_sub_pairs) )
      allocate( neighbor_sigma(1:n_sub_pairs) )
      allocate( omegas(1:n_sub_pairs) )
      allocate( T_func(1:9*n_sub_pairs) )
      allocate( B_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
      allocate( f_damp(1:n_sub_pairs) )
      allocate( g_func(1:n_sub_pairs) )
      allocate( h_func(1:9*n_sub_pairs) )
      allocate( a_SCS(1:3*n_sub_sites,1:3) )
      allocate( ipiv(1:3*n_sub_sites) )
      allocate( work_arr(1:12*n_sub_sites) )
      allocate( rjs_0(1:n_sub_pairs) )
      !allocate( a_iso(1:n_sub_sites,1:2) )
      allocate( o_p(1:n_sub_sites) )
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        allocate( hirshfeld_v_sub_der(1:3,1:n_sub_sites) )
        hirshfeld_v_sub_der = 0.d0
      end if
        
      !allocate( ia(1:9*n_sub_pairs) )
      !allocate( ja(1:9*n_sub_pairs) )
      !allocate( val(1:9*n_sub_pairs) )
      allocate( b_i(1:3*n_sub_sites,1:3) )
      allocate( d_vec(1:3*n_sub_sites,1:3) )
        
      !a_iso = 0.d0

      do om = 1, 2

        xyz_H = 0.d0
        rjs_H = 0.d0
        r0_ii = 0.d0
        neighbor_alpha0 = 0.d0
        neighbor_sigma = 0.d0

        f_damp = 0.d0
        g_func = 0.d0
        h_func = 0.d0

        !nnz = 0
        k2 = 0
        T_func = 0.d0
        B_mat = 0.d0

        a_SCS = 0.d0
        b_i = 0.d0

        d_vec = 0.d0
        do p = 1, n_sub_sites
          do c1 = 1, 3
            d_vec(3*(p-1)+c1,c1) = 1.d0
          end do
        end do

        p_list = 0
        n_sub_neigh = 0

        k_i = 0
        p = 0

        do i3 = 1, n_neigh(i)
          k_i = k_i+1
          i2 = neighbors_list(n_tot+k_i)
          if ( rjs(n_tot+k_i) < rcut ) then
            p = p+1 
            k2 = k2+1
            rjs_0(k2) = rjs(n_tot+k_i)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = i2
            if ( do_derivatives .and. do_hirshfeld_gradients ) then
              hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
            end if
            neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
            xyz_H(:,k2) = xyz(:,n_tot+k_i)/Bohr
            xyz_i = xyz_H(:,k2)
            rjs_H(k2) = rjs(n_tot+k_i)/Bohr
            r_vdw_i = r0_ii(k2)
            s_i = neighbor_sigma(k2)
            if (p == 1) then
              do c1 = 1, 3
                b_i(c1,c1) = 1.d0/neighbor_alpha0(k2)
              end do
            end if
            do c1 = 1, 3
              B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/neighbor_alpha0(k2)
              !nnz = nnz+1
              !ia(nnz) = 3*(p-1)+c1
              !ja(nnz) = 3*(p-1)+c1
              !val(nnz) = 1.d0/neighbor_alpha0(k2)
            end do
            mult1_i = 1.d0
            if ( rjs(n_tot+k_i) .ge. rcut-r_buffer ) then
              mult1_i = mult1_i * &
                              (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**2 &
                               + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**3 )
            end if
            k_j = 0
            q = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
                if ( rjs(n_tot+k_j) < rcut ) then
                  q = q+1
                end if
                if (i3 .ne. j3) then
                  xyz_j = xyz(:,n_tot+k_j)/Bohr
                  if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                    n_sub_neigh(p) = n_sub_neigh(p) + 1
                    j = neighbors_list(n_tot+k_j)
                    k2 = k2+1
                    rjs_0(k2) = rjs(n_tot+k_j)
                    if ( rjs(n_tot+k_j) < rcut ) then
                      p_list(k2) = q
                    end if    
                    s = neighbor_species(n_tot+k_j)
                    sub_neighbors_list(k2) = j                        
                    r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                    omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
                    if ( om == 2 ) then
                      neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                    else
                      neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)) / &
                        (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
                    end if
                    neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
                    xyz_H(:,k2) = xyz_j-xyz_i
                    rjs_H(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                    r_vdw_j = r0_ii(k2)
                    f_damp(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                    s_j = neighbor_sigma(k2)
                    sigma_ij = sqrt(s_i**2 + s_j**2)
                    g_func(k2) = erf(rjs_H(k2)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij) * exp(-rjs_H(k2)**2.d0/sigma_ij**2)
                    k3 = 9*(k2-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k3 = k3 + 1
                        if (c1 == c2) then
                          T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c1,k2) - rjs_H(k2)**2)/rjs_H(k2)**5
                        else
                          T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c2,k2))/rjs_H(k2)**5
                        end if
                        h_func(k3) = 4.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij)**3 * &
                                        xyz_H(c1,k2)*xyz_H(c2,k2)/rjs_H(k2)**5 * exp(-rjs_H(k2)**2/sigma_ij**2)
                        if ( rjs(n_tot+k_j) < rcut ) then
                          mult1_j = 1.d0
                          mult2 = 1.d0
                          if ( rjs(n_tot+k_j) .ge. rcut-r_buffer ) then
                            mult1_j = mult1_j * &
                                           (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**2 &
                                                 + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**3 )
                          end if
                          if ( rjs_H(k2) .ge. (rcut-r_buffer)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                            mult2 = mult2 * &
                                           (1.d0 - 3.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 &
                                                 + 2.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**3 )
                          end if
                          B_mat(3*(p-1)+c1,3*(q-1)+c2) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                        g_func(k2) + h_func(k3)) * mult1_i * mult1_j * mult2
                          if ( p == 1 ) then
                            b_i(3*(q-1)+c2,c1) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                  g_func(k2) + h_func(k3)) * mult1_i * mult1_j * mult2
                          end if
                          !nnz = nnz+1
                          !ia(nnz) = 3*(p-1)+c1
                          !ja(nnz) = 3*(q-1)+c2
                          !val(nnz) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                          !           g_func(k2) + h_func(k3)) * mult1_j * mult2
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                    g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                    mult1_i * ( 1.d0 - mult1_j ) * mult2
                        else ! rjs(n_tot+k_j) < rcut
                          mult1_j = 1.d0
                          mult2 = 1.d0
                          if ( rjs(n_tot+k_j) .ge. 2.d0*rcut-r_buffer ) then
                            mult1_j = mult1_j * &
                                           (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer))**2 &
                                                 + 2.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer))**3 )
                          end if
                          if ( rjs_H(k2) .ge. (rcut-r_buffer)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                            mult2 = mult2 * &
                                          (1.d0 - 3.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 &
                                                + 2.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**3 )
                          end if
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                      mult1_j * mult2
                        end if
                      end do
                    end do
                  end if
                end if
              end if
            end do
          end if
        end do
        
        !call cpu_time(time1)

        if ( polynomial_expansion ) then
        
          if (om == 2) then
            polyfit = (/ 3.237385145550585d+02, -4.241125470183307d+04, 3.008572712845031d+06, &
                        -1.309430416378132d+08, 3.756106046665028d+09, -7.433108326602138d+10, &
                         1.045457946646248d+12, -1.064278184563909d+13, 7.908875986032283d+13, &
                        -4.286093180295281d+14, 1.673744363742281d+15, -4.582925299269683d+15, &
                         8.343821283398570d+15, -9.066835011532402d+15, 4.447833222479864d+15 /)

          else
            polyfit = (/ 3.464569029392560e+01, -5.541785287730104e+02, 5.429135883990769e+03, &
                        -3.643266484189429e+04, 1.774056111172961e+05, -6.476342669175469e+05, &
                         1.805107059718787e+06, -3.873561223705132e+06, 6.400647427099578e+06, &
                        -8.079084665579021e+06, 7.651382671329762e+06, -5.264071084777111e+06, &
                         2.484267601324048e+06, -7.192678344511647e+05, 9.633392189452502e+04 /)
          end if
         
          n_degree = size(polyfit)-1

          !allocate( myidx(1:3*n_sub_sites) )
          allocate( val_xv(1:3*n_sub_sites,1:3) )
          val_xv = 0.d0
          !k2 = 0
          !do c1 = 1, 3
          !  do p = 1, n_sub_sites
          !    k2 = k2+1
          !    myidx(k2) = k2
          !  end do
          !end do

          !call psb_init(icontxt)
          a_SCS = polyfit(2)*b_i
          !call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
          !call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
          !call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
          !call psb_cdasb(desc_a, info_psb)
          !call psb_spasb(A_sp, desc_a, info_psb)
          do k2 = 3, n_degree+1
            !call psb_spmm(1.d0, A_sp, b_i, 0.d0, val_xv, desc_a, info_psb, 'T')
            call dgemm( "n", "n", 3*n_sub_sites, 3, 3*n_sub_sites, 1.d0, B_mat, 3*n_sub_sites, b_i, &
                        3*n_sub_sites, 0.d0, val_xv, 3*n_sub_sites)
            a_SCS = a_SCS + polyfit(k2) * val_xv
            b_i = val_xv 
          end do
          if ( om == 1 ) then
            pol1 = polyfit(1)
            do c1 = 1, 3
              pol1 = pol1 + dot_product(a_SCS(:,c1),d_vec(:,c1))/3.d0
            end do
            !write(*,*) "pol1", pol1
          else
            central_pol(i) = polyfit(1)
            do c1 = 1, 3
              central_pol(i) = central_pol(i) + dot_product(a_SCS(:,c1),d_vec(:,c1))/3.d0
            end do
            central_omega(i) = 2.d0*omega_ref/sqrt(central_pol(i)/pol1-1.d0)
            write(*,*) "Polynomial central polarizability", i0, central_pol(i), central_omega(i)
          end if
      
          deallocate( val_xv )
          !deallocate( myidx )
        
        else
        
          a_SCS = d_vec

          call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, a_SCS, 3*n_sub_sites, work_arr, &
                      12*n_sub_sites, info )
                    

          if ( om == 1 ) then
            pol1 = 0.d0
            do c1 = 1, 3
              pol1 = pol1 + a_SCS(c1,c1) 
            end do
            pol1 = pol1/3.d0
            central_pol(i) = 0.d0
            do c1 = 1, 3
              central_pol(i) = central_pol(i) + a_SCS(c1,c1)
            end do
            central_pol(i) = central_pol(i)/3.d0
          else
            central_pol(i) = 0.d0
            do c1 = 1, 3
              central_pol(i) = central_pol(i) + a_SCS(c1,c1)
            end do
            central_pol(i) = central_pol(i)/3.d0
            write(*,*) i0, central_pol(i)
            central_omega(i) = 2.d0*omega_ref/sqrt(central_pol(i)/pol1-1.d0)
          end if
        
        end if

        !call cpu_time(time2)

        !write(*,*) "Timing for solving the polarizability of central atom:", time2-time1

      end do
        
      deallocate( sub_neighbors_list, n_sub_neigh, p_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, &
                  omegas, T_func, B_mat, f_damp, g_func, h_func, a_SCS, ipiv, work_arr, rjs_0, o_p )
      !deallocate( a_iso )
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        deallocate( hirshfeld_v_sub_der )
      end if
        
      deallocate( b_i, d_vec )
      !deallocate( ia, ja, val )        

    end do 

    call cpu_time(time2)

    write(*,*) "Central polarizabilities timing", time2-time1

      
!******************************************************************************************
! This is where you would break this into get_scs and get_mbd for two separate subroutines:
! Store central_pol and central_omega and pass them to get_mbd
!******************************************************************************************

! TEST !!!!!!!!!!!!!!!!
central_pol = 10.d0
central_omega = 0.5d0


    call cpu_time(time1)

    do i = 1, n_sites      

      n_tot = sum(n_neigh(1:i))-n_neigh(i)
      i0 = modulo(neighbors_list(n_tot+1)-1, n_sites0) + 1
      s = neighbor_species(n_tot+1)
      omega_ref = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
      n_sub_sites = 0
      n_sub_pairs = 0
      k_i = 0
      do i3 = 1, n_neigh(i)
        k_i = k_i + 1
        if (rjs(n_tot+k_i) < rcut ) then
          n_sub_sites = n_sub_sites + 1
          n_sub_pairs = n_sub_pairs + 1
          xyz_i = xyz(:,n_tot+k_i)/Bohr
          k_j = 0
          do j3 = 1, n_neigh(i)
            k_j = k_j + 1
            if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
              if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                  n_sub_pairs = n_sub_pairs + 1
                end if
              end if
            end if
          end do
        end if
      end do
        
      allocate( sub_neighbors_list(1:n_sub_pairs) )
      allocate( n_sub_neigh(1:n_sub_sites) )
      allocate( p_list(1:n_sub_pairs) )
      allocate( xyz_H(1:3,1:n_sub_pairs) )
      allocate( rjs_H(1:n_sub_pairs) )
      allocate( r0_ii(1:n_sub_pairs) )
      allocate( neighbor_alpha0(1:n_sub_pairs) )
      allocate( neighbor_sigma(1:n_sub_pairs) )
      allocate( omegas(1:n_sub_pairs) )
      allocate( T_func(1:9*n_sub_pairs) )
      allocate( B_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
      allocate( f_damp(1:n_sub_pairs) )
      allocate( g_func(1:n_sub_pairs) )
      allocate( h_func(1:9*n_sub_pairs) )
      allocate( a_SCS(1:3*n_sub_sites,1:3) )
      allocate( ipiv(1:3*n_sub_sites) )
      allocate( work_arr(1:12*n_sub_sites) )
      allocate( rjs_0(1:n_sub_pairs) )
      allocate( a_iso(1:n_sub_sites,1:2) )
      allocate( o_p(1:n_sub_sites) )
      allocate( T_SR(1:9*n_sub_pairs) )
      allocate( T_SR_mult(1:n_sub_pairs) )
      allocate( d_arr_i(1:9*n_sub_pairs) )
      allocate( d_arr_o(1:9*n_sub_pairs) )
      allocate( d_mult_i(1:n_sub_pairs) )
      allocate( d_mult_o(1:n_sub_pairs) )
      allocate( dT_SR_mult(1:n_sub_pairs,1:3) )
      allocate( d_dmult_i(1:n_sub_pairs,1:3) )
      allocate( d_dmult_o(1:n_sub_pairs,1:3) )
      allocate( hirshfeld_sub_neigh(1:n_sub_pairs) )
      if ( do_derivatives ) then
        allocate( da_iso(1:n_sub_sites,1:3,1:2) )
        da_iso = 0.d0
      end if
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        allocate( hirshfeld_v_sub_der(1:3,1:n_sub_sites) )
        hirshfeld_v_sub_der = 0.d0
      end if
        
      !allocate( ia(1:9*n_sub_pairs) )
      !allocate( ja(1:9*n_sub_pairs) )
      !allocate( val(1:9*n_sub_pairs) )
      allocate( b_i(1:3*n_sub_sites,1:3) )
      allocate( d_vec(1:3*n_sub_sites,1:3) )
      if ( polynomial_expansion ) then
        allocate( B_pol(1:3*n_sub_sites,1:3*n_sub_sites) )
      end if

      a_iso = 0.d0
      T_SR = 0.d0
      T_SR_mult = 0.d0
      d_arr_i = 0.d0
      d_arr_o = 0.d0
      d_mult_i = 0.d0
      d_mult_o = 0.d0
      dT_SR_mult = 0.d0
      d_dmult_i = 0.d0
      d_dmult_o = 0.d0

      do om = 1, 2

        xyz_H = 0.d0
        rjs_H = 0.d0
        r0_ii = 0.d0
        neighbor_alpha0 = 0.d0
        neighbor_sigma = 0.d0
        hirshfeld_sub_neigh = 0.d0

        f_damp = 0.d0
        g_func = 0.d0
        h_func = 0.d0

        !nnz = 0
        k2 = 0
        T_func = 0.d0
        B_mat = 0.d0

        a_SCS = 0.d0
        b_i = 0.d0

        d_vec = 0.d0
        do p = 1, n_sub_sites
          do c1 = 1, 3
            d_vec(3*(p-1)+c1,c1) = 1.d0
          end do
        end do

        p_list = 0
        n_sub_neigh = 0
        k_i = 0
        p = 0

        do i3 = 1, n_neigh(i)
          k_i = k_i+1
          i2 = neighbors_list(n_tot+k_i)
          if ( rjs(n_tot+k_i) < rcut ) then
            p = p+1 
            k2 = k2+1
            rjs_0(k2) = rjs(n_tot+k_i)
            hirshfeld_sub_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = i2
            if ( do_derivatives .and. do_hirshfeld_gradients ) then
              hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
            end if
            neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
            xyz_H(:,k2) = xyz(:,n_tot+k_i)/Bohr
            xyz_i = xyz_H(:,k2)
            rjs_H(k2) = rjs(n_tot+k_i)/Bohr
            r_vdw_i = r0_ii(k2)
            s_i = neighbor_sigma(k2)
            if (p == 1) then
              do c1 = 1, 3
                b_i(c1,c1) = 1.d0/neighbor_alpha0(k2)
              end do
            end if
            k3 = 9*(k2-1)
            do c1 = 1, 3
              do c2 = 1, 3
                k3 = k3+1
                if ( c1 == c2 ) then
                B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/neighbor_alpha0(k2)
                !nnz = nnz+1
                !ia(nnz) = 3*(p-1)+c1
                !ja(nnz) = 3*(p-1)+c1
                !val(nnz) = 1.d0/neighbor_alpha0(k2)
                T_SR(k3) = 1.d0/neighbor_alpha0(k2)
                end if
              end do
            end do
            mult1_i = 1.d0
            dmult1_i = 0.d0
            if ( rjs(n_tot+k_i) .ge. rcut-r_buffer ) then
              mult1_i = mult1_i * &
                              (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**2 &
                               + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**3 )
              dmult1_i = (- 6.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer)) &
                              + 6.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**2 ) &
                              * ( -xyz(:,n_tot+k_i)/rjs(n_tot+k_i)/(r_buffer/Bohr))
            end if
            k_j = 0
            q = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
                if ( rjs(n_tot+k_j) < rcut ) then
                  q = q+1
                end if
                if (i3 .ne. j3) then
                  xyz_j = xyz(:,n_tot+k_j)/Bohr
                  if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                    n_sub_neigh(p) = n_sub_neigh(p) + 1
                    j = neighbors_list(n_tot+k_j)
                    k2 = k2+1
                    rjs_0(k2) = rjs(n_tot+k_j)
                    hirshfeld_sub_neigh(k2) = hirshfeld_v_neigh(n_tot+k_j)
                    if ( rjs(n_tot+k_j) < rcut ) then
                      p_list(k2) = q
                    end if    
                    s = neighbor_species(n_tot+k_j)
                    sub_neighbors_list(k2) = j                        
                    r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                    omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
                    if ( om == 2 ) then
                      neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                    else
                      neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)) / &
                        (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
                    end if
                    neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
                    xyz_H(:,k2) = xyz_j-xyz_i
                    rjs_H(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                    r_vdw_j = r0_ii(k2)
                    f_damp(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                    s_j = neighbor_sigma(k2)
                    sigma_ij = sqrt(s_i**2 + s_j**2)
                    g_func(k2) = erf(rjs_H(k2)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij) * exp(-rjs_H(k2)**2.d0/sigma_ij**2)
                    k3 = 9*(k2-1)
                    if ( rjs(n_tot+k_j) < rcut ) then
                      mult1_j = 1.d0
                      dmult1_j = 0.d0
                      mult2 = 1.d0
                      dmult2 = 0.d0
                      if ( rjs(n_tot+k_j) .ge. rcut-r_buffer ) then
                        mult1_j = mult1_j * &
                                       (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**2 &
                                             + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**3 )
                        dmult1_j = (- 6.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer)) &
                                + 6.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**2 ) &
                                * ( -xyz(:,n_tot+k_j)/rjs(n_tot+k_j)/(r_buffer/Bohr))
                      end if
                      if ( rjs_H(k2) .ge. (rcut-r_buffer)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                        mult2 = mult2 * &
                                       (1.d0 - 3.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 &
                                             + 2.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**3 )
                        if ( i2 == i ) then
                          dmult2 = (- 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr)) &
                            + 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 ) &
                            * ( -xyz_H(:,k2)/rjs_H(k2)/(r_buffer/Bohr))
                        end if
                        if ( j == i ) then
                          dmult2 = (- 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr)) &
                            + 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 ) &
                            * ( xyz_H(:,k2)/rjs_H(k2)/(r_buffer/Bohr))
                        end if
                      end if
                      T_SR_mult(k2) = mult1_i * mult1_j * mult2
                      dT_SR_mult(k2,:) = dmult1_i * mult1_j * mult2 + mult1_i * dmult1_j * mult2 + &
                                         mult1_i * mult1_j * dmult2
                      d_mult_i(k2) = mult1_i * ( 1.d0 - mult1_j ) * mult2
                      d_dmult_i(k2,:) = dmult1_i * (1.d0 - mult1_j ) * mult2 + mult1_i * (-dmult1_j) * mult2 &
                                      + mult1_i * ( 1.d0 - mult1_j ) * dmult2
                    else ! rjs(n_tot+k_j) < rcut
                      mult1_j = 1.d0
                      dmult1_j = 0.d0
                      mult2 = 1.d0
                      dmult2 = 0.d0
                      if ( rjs(n_tot+k_j) .ge. 2.d0*rcut-r_buffer ) then
                        mult1_j = mult1_j * &
                                       (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer))**2 &
                                             + 2.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer))**3 )
                        dmult1_j = (- 6.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer)) &
                                + 6.d0 * ((rjs(n_tot+k_j)-2.d0*rcut+r_buffer)/(r_buffer))**2 ) &
                                * ( -xyz(:,n_tot+k_j)/rjs(n_tot+k_j)/(r_buffer/Bohr))
                      end if
                      if ( rjs_H(k2) .ge. (rcut-r_buffer)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                        mult2 = mult2 * &
                                      (1.d0 - 3.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 &
                                            + 2.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**3 )
                        if ( i2 == i ) then
                          dmult2 = (- 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr)) &
                            + 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 ) &
                            * ( -xyz_H(:,k2)/rjs_H(k2)/(r_buffer/Bohr))
                        end if
                        if ( j == i ) then
                          dmult2 = (- 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr)) &
                            + 6.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 ) &
                            * ( xyz_H(:,k2)/rjs_H(k2)/(r_buffer/Bohr))
                        end if
                      end if
                      d_mult_o(k2) = mult1_j * mult2
                      d_dmult_o(k2,:) = dmult1_j * mult2 + mult1_j * dmult2
                    end if
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k3 = k3 + 1
                        if (c1 == c2) then
                          T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c1,k2) - rjs_H(k2)**2)/rjs_H(k2)**5
                        else
                          T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c2,k2))/rjs_H(k2)**5
                        end if
                        h_func(k3) = 4.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij)**3 * &
                                        xyz_H(c1,k2)*xyz_H(c2,k2)/rjs_H(k2)**5 * exp(-rjs_H(k2)**2/sigma_ij**2)
                        if ( rjs(n_tot+k_j) < rcut ) then
                          B_mat(3*(p-1)+c1,3*(q-1)+c2) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                        g_func(k2) + h_func(k3)) * T_SR_mult(k2)
                          T_SR(k3) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                        g_func(k2) + h_func(k3))
                          if ( p == 1 ) then
                            b_i(3*(q-1)+c2,c1) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                  g_func(k2) + h_func(k3)) * T_SR_mult(k2)
                          end if
                          !nnz = nnz+1
                          !ia(nnz) = 3*(p-1)+c1
                          !ja(nnz) = 3*(q-1)+c2
                          !val(nnz) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                          !           g_func(k2) + h_func(k3)) * T_SR_mult(k2)
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                    g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                    d_mult_i(k2)
                          d_arr_i(k3) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                    g_func(k2) + h_func(k3)) * neighbor_alpha0(k2)
                        else
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                      d_mult_o(k2)
                          d_arr_o(k3) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * neighbor_alpha0(k2)
                        end if
                      end do
                    end do
                  end if
                end if
              end if
            end do
          end if
        end do

        !if ( i == 1 .and. om == 2 ) then
        !  write(*,*) "B_mat"
        !  do p = 1, 3*n_sub_sites
        !    write(*,*) B_mat(p,:)
        !  end do
        !end if
        
        if ( polynomial_expansion ) then

          if (om == 2) then
            polyfit = (/ 3.237385145550585d+02, -4.241125470183307d+04, 3.008572712845031d+06, &
                        -1.309430416378132d+08, 3.756106046665028d+09, -7.433108326602138d+10, &
                         1.045457946646248d+12, -1.064278184563909d+13, 7.908875986032283d+13, &
                        -4.286093180295281d+14, 1.673744363742281d+15, -4.582925299269683d+15, &
                         8.343821283398570d+15, -9.066835011532402d+15, 4.447833222479864d+15 /)

          else
            polyfit = (/ 3.464569029392560e+01, -5.541785287730104e+02, 5.429135883990769e+03, &
                        -3.643266484189429e+04, 1.774056111172961e+05, -6.476342669175469e+05, &
                         1.805107059718787e+06, -3.873561223705132e+06, 6.400647427099578e+06, &
                        -8.079084665579021e+06, 7.651382671329762e+06, -5.264071084777111e+06, &
                         2.484267601324048e+06, -7.192678344511647e+05, 9.633392189452502e+04 /)
          end if
         
          n_degree = size(polyfit)-1

          !allocate( myidx(1:3*n_sub_sites) )
          allocate( val_xv(1:3*n_sub_sites,1:3*n_sub_sites) )
          !allocate( B_pol(1:3*n_sub_sites,1:3*n_sub_sites) )
          allocate( B_mult(1:3*n_sub_sites,1:3*n_sub_sites) )
          B_mult = B_mat
          val_xv = 0.d0
          !k2 = 0
          !do c1 = 1, 3
          !  do p = 1, n_sub_sites
          !    k2 = k2+1
          !    myidx(k2) = k2
          !  end do
          !end do

          !call psb_init(icontxt)
          B_pol = polyfit(2)*B_mat
          do p = 1, 3*n_sub_sites
            B_pol(p,p) = B_pol(p,p) + polyfit(1)
          end do
          !call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
          !call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
          !call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
          !call psb_cdasb(desc_a, info_psb)
          !call psb_spasb(A_sp, desc_a, info_psb)

          !call cpu_time(time1)

          do k2 = 3, n_degree+1
            ! ATTENTION: For some reason psb_spmm slows down after the first iteration and dgemm perfoms faster.
            !call psb_spmm(1.d0, A_sp, B_mult, 0.d0, val_xv, desc_a, info_psb, 'T')
            call dgemm( "n", "n", 3*n_sub_sites, 3*n_sub_sites, 3*n_sub_sites, 1.d0, B_mult, 3*n_sub_sites, B_mat, &
                        3*n_sub_sites, 0.d0, val_xv, 3*n_sub_sites)
            B_pol = B_pol + polyfit(k2) * val_xv
            B_mult = val_xv 
          end do

          !call cpu_time(time2)

          !write(*,*) "Timing for matrix multiplications", time2-time1

          a_SCS = 0.d0
          do p = 1, n_sub_sites
            do c1 = 1, 3
              do c2 = 1, 3
                a_SCS(3*(p-1)+c1,c2) = a_SCS(3*(p-1)+c1,c2) + dot_product(B_pol(3*(p-1)+c1,:),d_vec(:,c2))
              end do
            end do
          end do
          do p = 1, n_sub_sites
            do c1 = 1, 3
              a_iso(p,om) = a_iso(p,om) + a_SCS(3*(p-1)+c1,c1)
            end do
          end do
          a_iso(:,om) = a_iso(:,om)/3.d0

          if ( om == 2 ) then
            do p = 1, n_sub_sites
              o_p(p) = 2.d0*omega_ref/sqrt(a_iso(p,2)/a_iso(p,1)-1.d0)
            end do
          end if

          deallocate( val_xv, B_mult )
          !deallocate( myidx )

        else

          ! This is the exact solution:
          a_SCS = d_vec

          call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, a_SCS, 3*n_sub_sites, work_arr, &
                      12*n_sub_sites, info )

          do p = 1, n_sub_sites
            do c1 = 1, 3
              a_iso(p,om) = a_iso(p,om) + a_SCS(3*(p-1)+c1,c1)
            end do
          end do
          a_iso(:,om) = a_iso(:,om)/3.d0
        
          if ( om == 2 ) then
            do p = 1, n_sub_sites
              o_p(p) = 2.d0*omega_ref/sqrt(a_iso(p,2)/a_iso(p,1)-1.d0)
            end do
          end if

          !do c1 = 1, 3
          !  do c2 = 1, 3
          !    alpha_SCS_full(3*(i-1)+c1,c2,om) = a_SCS(c1,c2)
          !  end do
          !end do
          
          !alpha_SCS0(i,om) = 0.d0
          !do c1 = 1, 3
          !  alpha_SCS0(i,om) = alpha_SCS0(i,om) + a_SCS(c1,c1)
          !end do
          !alpha_SCS0(i,om) = alpha_SCS0(i,om)/3.d0
          ! Exact solution ends
        
        end if

        if ( om == 2 ) then

          ! MBD for local polarizabilities:
          ! At least for now: rcut <= rcut_mbd <= rcut_2b
          !rcut_mbd = 5.d0
          !rcut_2b = 9.d0
          n_mbd_sites = 0
          n_mbd_pairs = 0
        
          k_i = 0
          do i3 = 1, n_neigh(i)
            k_i = k_i + 1
            if (rjs(n_tot+k_i) < rcut_mbd ) then
              n_mbd_sites = n_mbd_sites + 1
              n_mbd_pairs = n_mbd_pairs + 1
              xyz_i = xyz(:,n_tot+k_i)/Bohr
              k_j = 0
              do j3 = 1, n_neigh(i)
                k_j = k_j + 1
                if ( rjs(n_tot+k_j) < rcut_mbd ) then
                  if (i3 .ne. j3) then
                    xyz_j = xyz(:,n_tot+k_j)/Bohr
                    if ( sqrt(sum((xyz_j-xyz_i)**2)) < (rcut_mbd)/Bohr ) then
                      n_mbd_pairs = n_mbd_pairs + 1
                    end if
                  end if
                end if
              end do
            end if
          end do
        
          k_i = 0
          n_2b_sites = 0
          do j2 = 1, n_neigh(i)
            k_i = k_i + 1
            if (rjs(n_tot+k_i) < rcut_2b .and. rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer) then
              n_2b_sites = n_2b_sites + 1
            end if
          end do

          !n_order = 4
          !n_freq = 12

          allocate( T_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( r0_ii_SCS(1:n_mbd_pairs) )
          allocate( f_damp_SCS(1:n_mbd_pairs) )
          allocate( AT(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_freq) )
          allocate( AT_n(1:3*n_mbd_sites,1:3,n_order-1,1:n_freq) )
          allocate( energy_series(1:3*n_mbd_sites,1:3) )
          allocate( omegas_mbd(1:n_freq) )
          allocate( integrand(1:n_freq) )
          allocate( n_mbd_neigh(1:n_mbd_sites) )
          allocate( mbd_neighbors_list(1:n_mbd_pairs) )
          allocate( p_mbd(1:n_mbd_pairs) )
          allocate( r0_ii_mbd(1:n_mbd_pairs) )
          allocate( neighbor_alpha0_mbd(1:n_mbd_pairs) )
          allocate( xyz_mbd(1:3,n_mbd_pairs) )
          allocate( rjs_mbd(n_mbd_pairs) )
          allocate( T_mbd(1:9*n_mbd_pairs) )
          allocate( a_mbd(1:n_mbd_pairs) )
          allocate( o_mbd(1:n_mbd_pairs) )
          allocate( rjs_0_mbd(1:n_mbd_pairs) )
          allocate( xyz_0_mbd(1:3,n_mbd_pairs) )
          if ( do_derivatives ) then
            allocate( da_mbd(1:n_mbd_pairs) )
            allocate( AT_n_f(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_order-1,1:n_freq) )
            allocate( dT_mbd(1:9*n_mbd_pairs) )
            allocate( f_damp_der_mbd(1:n_mbd_pairs) )
            allocate( f_damp_der_SCS(1:n_mbd_pairs) )
            allocate( dT_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
            allocate( G_mat(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_freq) )
            allocate( force_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
            allocate( total_energy_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
            allocate( total_integrand(1:n_freq) )
            allocate( do_mbd(1:n_mbd_pairs) )
          end if
        
          allocate( sub_2b_list(1:n_2b_sites) )
          allocate( xyz_2b(1:3,1:n_2b_sites) )
          allocate( rjs_2b(1:n_2b_sites) )
          allocate( r0_ii_2b(1:n_2b_sites) )
          allocate( neighbor_alpha0_2b(1:n_2b_sites) )
          allocate( a_2b(1:n_2b_sites) )
          allocate( o_2b(1:n_2b_sites) )
          allocate( r0_ii_SCS_2b(1:n_2b_sites) )
          allocate( f_damp_SCS_2b(1:n_2b_sites) )
          allocate( C6_2b(1:n_2b_sites) )
          if ( do_derivatives ) then
            allocate( da_2b(1:n_2b_sites) )
            allocate( do_2b(1:n_2b_sites) )
          end if
        
          a_mbd = 0.d0
        
          omegas_mbd = 0.d0
          omega = 4.d0 ! Max frequency for integration
          ! Using non-uniform distribution because the polarizabilities change more rapidly with low frequency values:
          do i2 = 1, n_freq
            omegas_mbd(i2) = (i2-1)**2 * omega/(n_freq-1)**2
          end do

          T_LR = 0.d0
          AT = 0.d0
          AT_n = 0.d0
          r0_ii_SCS = 0.d0
          f_damp_SCS = 0.d0
        
          n_mbd_neigh = 0
          mbd_neighbors_list = 0
          p_mbd = 0
          r0_ii_mbd = 0.d0
          neighbor_alpha0_mbd = 0.d0
          xyz_mbd = 0.d0
          rjs_mbd = 0.d0
          T_mbd = 0.d0

          k2 = 0
          k_i = 0
          p = 0
          do i3 = 1, n_neigh(i)
            k_i = k_i+1
            i2 = neighbors_list(n_tot+k_i)
            i1 = modulo(i2-1, n_sites0) + 1
            if ( rjs(n_tot+k_i) < rcut_mbd ) then
              p = p+1
              k2 = k2+1
              s = neighbor_species(n_tot+k_i)
              mbd_neighbors_list(k2) = i2
              n_mbd_neigh(p) = n_mbd_neigh(p) + 1
              p_mbd(k2) = p
              r0_ii_mbd(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
              xyz_mbd(:,k2) = xyz(:,n_tot+k_i)/Bohr
              xyz_i = xyz_mbd(:,k2)
              rjs_mbd(k2) = rjs(n_tot+k_i)/Bohr
              rjs_0_mbd(k2) = rjs(n_tot+k_i)/Bohr
              xyz_0_mbd(:,k2) = xyz(:,n_tot+k_i)/Bohr
              neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
              if ( rjs(n_tot+k_i) < rcut-r_buffer ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                a_mbd(k2) = a_iso(r,2)
                o_mbd(k2) = o_p(r)
              else if ( rjs(n_tot+k_i) < rcut .and. rjs(n_tot+k_i) .ge. rcut-r_buffer ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                a_mbd(k2) = a_iso(r,2) * &
                            (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                  + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3) + &
                            central_pol(i1) * & 
                                  ( + 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                  - 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3)
                o_mbd(k2) = o_p(r) * &
                            (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                  + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3) + &
                            central_omega(i1) * & 
                                  ( + 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                  - 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3)
              else if ( rjs(n_tot+k_i) .ge. rcut .and. rjs(n_tot+k_i) < rcut_mbd-r_buffer) then
                a_mbd(k2) = central_pol(i1)
                o_mbd(k2) = central_omega(i1)
              else if ( rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer .and. rjs(n_tot+k_i) < rcut_mbd ) then
                a_mbd(k2) = central_pol(i1) * (1.d0 & 
                                    - 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                    + 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
                o_mbd(k2) = central_omega(i1) * (1.d0 & 
                                    - 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                    + 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
              end if
              r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_mbd(k2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
              r_vdw_i = r0_ii_SCS(k2)
              k_j = 0
              q = 0
              do j3 = 1, n_neigh(i)
                k_j = k_j+1
                if ( rjs(n_tot+k_j) < rcut_mbd ) then
                  q = q+1
                  if (i3 .ne. j3) then
                    xyz_j = xyz(:,n_tot+k_j)/Bohr
                    if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut_mbd/Bohr ) then
                      n_mbd_neigh(p) = n_mbd_neigh(p) + 1
                      j = neighbors_list(n_tot+k_j)
                      j1 = modulo(j-1, n_sites0) + 1
                      k2 = k2+1    
                      s = neighbor_species(n_tot+k_j)
                      mbd_neighbors_list(k2) = j
                      p_mbd(k2) = q                       
                      r0_ii_mbd(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                      xyz_mbd(:,k2) = xyz_j-xyz_i
                      rjs_mbd(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                      rjs_0_mbd(k2) = rjs(n_tot+k_j)/Bohr
                      xyz_0_mbd(:,k2) = xyz(:,n_tot+k_j)/Bohr
                      neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                      if ( rjs(n_tot+k_j) < rcut-r_buffer ) then
                        r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                        a_mbd(k2) = a_iso(r,2)
                        o_mbd(k2) = o_p(r)
                      else if ( rjs(n_tot+k_j) < rcut .and. rjs(n_tot+k_j) .ge. rcut-r_buffer ) then
                        r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                        a_mbd(k2) = a_iso(r,2) * &
                              (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                    + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3) + &
                        central_pol(j1) * & 
                                    ( + 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                      - 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3)
                        o_mbd(k2) = o_p(r) * &
                              (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                    + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3) + &
                        central_omega(j1) * & 
                                    ( + 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                    - 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3)
                      else if ( rjs(n_tot+k_j) .ge. rcut .and. rjs(n_tot+k_j) < rcut_mbd-r_buffer) then
                        a_mbd(k2) = central_pol(j1)
                        o_mbd(k2) = central_omega(j1)
                      else
                        a_mbd(k2) = central_pol(j1) * (1.d0 & 
                                            - 3.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                            + 2.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**3)
                        o_mbd(k2) = central_omega(j1) * (1.d0 & 
                                            - 3.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                            + 2.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**3)
                      end if
                      r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_mbd(k2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
                      r_vdw_j = r0_ii_SCS(k2)
                      f_damp_SCS(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_mbd(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                      k3 = 9*(k2-1)
                      do c1 = 1, 3
                        do c2 = 1, 3
                          k3 = k3 + 1
                          if (c1 == c2) then
                            T_mbd(k3) = (3*xyz_mbd(c1,k2) * xyz_mbd(c1,k2) - rjs_mbd(k2)**2)/rjs_mbd(k2)**5
                          else
                            T_mbd(k3) = (3*xyz_mbd(c1,k2) * xyz_mbd(c2,k2))/rjs_mbd(k2)**5
                          end if
                          T_LR(3*(p-1)+c1,3*(q-1)+c2) = f_damp_SCS(k2) * T_mbd(k3)
                        end do
                      end do
                    end if
                  end if
                end if
              end do
            end if
          end do
        
          E_TS = 0.d0
          k2 = 0
          k_i = 0
          a_2b = 0.d0
          o_2b = 0.d0
          rjs_2b = 0.d0
          xyz_2b = 0.d0
          sub_2b_list = 0
          r0_ii_2b = 0.d0
          neighbor_alpha0_2b = 0.d0
          C6_2b = 0.d0
          r0_ii_SCS_2b = 0.d0
          f_damp_SCS_2b = 0.d0
          s = neighbor_species(n_tot+1)
          r_vdw_i = r0_ref(s) / Bohr * (a_iso(1,2)/(alpha0_ref(s)/Bohr**3))**(1.d0/3.d0)
          do i3 = 1, n_neigh(i)
            k_i = k_i+1
            i2 = neighbors_list(n_tot+k_i)
            i1 = modulo(i2-1, n_sites0) + 1
            if ( rjs(n_tot+k_i) < rcut_2b .and. rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer ) then
              k2 = k2+1
              s = neighbor_species(n_tot+k_i)
              sub_2b_list(k2) = neighbors_list(n_tot+k_i)
              r0_ii_2b(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
              xyz_2b(:,k2) = xyz(:,n_tot+k_i)/Bohr
              xyz_i = xyz_2b(:,k2)
              rjs_2b(k2) = rjs(n_tot+k_i)/Bohr
              neighbor_alpha0_2b(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_i)
              if ( rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer .and. rjs(n_tot+k_i) < rcut_mbd ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                a_2b(k2) = central_pol(i1) * &
                                  ( + 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                  - 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
                o_2b(k2) = central_omega(i1) * &
                                  ( + 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                  - 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
              else if ( rjs(n_tot+k_i) .ge. rcut_mbd .and. rjs(n_tot+k_i) < rcut_2b-r_buffer) then
                a_2b(k2) = central_pol(i1)
                o_2b(k2) = central_omega(i1)
              else if ( rjs(n_tot+k_i) .ge. rcut_2b-r_buffer .and. rjs(n_tot+k_i) < rcut_2b) then
                a_2b(k2) = central_pol(i1) * (1.d0 & 
                           - 3.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**2 &
                           + 2.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**3)
                o_2b(k2) = central_omega(i1) * (1.d0 &
                           - 3.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**2 &
                           + 2.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**3)
              end if
              C6_2b(k2) = 3.d0/2.d0 * a_iso(1,2) * a_2b(k2) * (central_omega(i0) * o_2b(k2)) / &
                      (central_omega(i0) + o_2b(k2))
              r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (a_2b(k2)/neighbor_alpha0_2b(k2))**(1.d0/3.d0)
              r_vdw_j = r0_ii_SCS_2b(k2)
              f_damp_SCS_2b(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
              E_TS = E_TS - C6_2b(k2)/rjs_2b(k2)**6 * f_damp_SCS_2b(k2)
            end if
          end do
          E_TS = 1.d0/2.d0 * E_TS

          if ( i == 1 .and. om == 2 ) then
            write(*,*) "E_TS", E_TS
          end if

          do i2 = 1, n_freq
            k3 = 0
            do p = 1, n_mbd_sites
              AT(3*(p-1)+1:3*(p-1)+3,:,i2) = a_mbd(k3+1)/(1.d0+(omegas_mbd(i2)/o_mbd(k3+1))**2) &
                                             * T_LR(3*(p-1)+1:3*(p-1)+3,:)
              k3 = k3 + n_mbd_neigh(p)
            end do

            if ( n_order > 1 ) then
              do k2 = 1, n_order-1
                ! Precalculate the full AT_n for forces:
                if ( k2 == 1 ) then
                  AT_n(:,:,k2,i2) = AT(:,1:3,i2)
                  if ( do_derivatives ) then
                    AT_n_f(:,:,k2,i2) = AT(:,:,i2)
                  end if
                else
                  call dgemm('n', 'n', 3*n_mbd_sites, 3, 3*n_mbd_sites, 1.d0, AT(:,:,i2), &
                             3*n_mbd_sites, AT_n(:,:,k2-1,i2), 3*n_mbd_sites, 0.d0, AT_n(:,:,k2,i2), &
                             3*n_mbd_sites)
                  if ( do_derivatives ) then
                    call dgemm('n', 'n', 3*n_mbd_sites, 3*n_mbd_sites, 3*n_mbd_sites, 1.d0, AT(:,:,i2), &
                               3*n_mbd_sites, AT_n_f(:,:,k2-1,i2), 3*n_mbd_sites, 0.d0, AT_n_f(:,:,k2,i2), &
                               3*n_mbd_sites)
                  end if
                end if
              end do
            else
              write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
              write(*,*) "and forces will be zero."
            end if
          end do

          !if ( i == 1 .and. om == 2 ) then
          !  write(*,*) "AT"
          !  do p = 1, 3*n_sub_sites
          !    write(*,*) AT(p,:,1)
          !  end do
          !end if

          integrand = 0.d0
          do i2 = 1, n_freq
            energy_series = 0.d0
            do k2 = 1, n_order-1
              energy_series = energy_series - 1.d0/(k2+1) * AT_n(:,1:3,k2,i2) 
            end do
            do c1 = 1, 3
              integrand(i2) = integrand(i2) + a_mbd(1)/(1.d0 + (omegas_mbd(i2)/o_mbd(1))**2) & 
                              * dot_product(T_LR(c1,:),energy_series(:,c1))
            end do
          end do

          integral = 0.d0
          call integrate("trapezoidal", omegas_mbd, integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
          integral = integral/(2.d0*pi)
          energies(i) = (integral + E_TS) * Hartree
          write(*,*) "MBD energy", i, energies(i)
          E_MBD = E_MBD + energies(i)          
        
        end if ! om loop
        
        ! Derivatives:
        
        if (do_derivatives) then
        
          allocate( da_SCS(1:3*n_sub_sites,1:3) )
          allocate( dT(1:9*n_sub_pairs) )
          allocate( dB_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
          dB_mat = 0.d0
          allocate( b_der(1:3*n_sub_sites,1:3) )
          allocate( f_damp_der(1:n_sub_pairs) )
          allocate( g_func_der(1:n_sub_pairs) )
          allocate( h_func_der(1:9*n_sub_pairs) )
          allocate( d_der(1:3*n_sub_sites,1:3) )        

          do c3 = 1, 3
          
            dT = 0.d0
            k2 = 0
            do p = 1, n_sub_sites 
              k2 = k2+1
              do j3 = 2, n_sub_neigh(p)
                k2 = k2+1
                k3 = 9*(k2-1)
                do c1 = 1,3
                  do c2 = 1,3
                    k3 = k3+1
                    dT(k3) = (-15.d0 * xyz_H(c1,k2) * xyz_H(c2,k2) * xyz_H(c3,k2))/rjs_H(k2)**7
                    if (c1 == c2) then
                      dT(k3) = dT(k3) + 3.d0/rjs_H(k2)**5 * xyz_H(c3,k2)
                    end if
                    if (c2 == c3) then
                      dT(k3) = dT(k3) + 3.d0/rjs_H(k2)**5 * xyz_H(c1,k2)
                    end if
                    if (c1 == c3) then
                      dT(k3) = dT(k3) + 3.d0/rjs_H(k2)**5 * xyz_H(c2,k2)
                    end if
                  end do
                end do
              end do
            end do
                
            ! Derivatives are always w.r.t. the central atom:
            a = i

            dB_mat = 0.d0
            f_damp_der = 0.d0
            g_func_der = 0.d0
            h_func_der = 0.d0
            b_der = 0.d0
            d_der = 0.d0
          
            k3 = 0
            do p = 1, n_sub_sites
              k3 = k3+1
              r_vdw_i = r0_ii(k3)
              s_i = neighbor_sigma(k3)
              i2 = sub_neighbors_list(k3)
              do j2 = 2, n_sub_neigh(p)
                k3 = k3+1
                r_vdw_j = r0_ii(k3)
                s_j = neighbor_sigma(k3)
                j = sub_neighbors_list(k3)
                if (a == i2 .or. a == j) then
                  if ( rjs_0(k3) < rcut ) then
                    q = p_list(k3)
                  end if
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  f_damp_der(k3) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k3)**2 * &
                                     exp( -d*(rjs_H(k3)/(sR*(r_vdw_i + &
                                     r_vdw_j)) - 1.d0) ) * xyz_H(c3,k3)/rjs_H(k3)
                  g_func_der(k3) = 4.d0/sqrt(pi) * rjs_H(k3)/sigma_ij**3 * xyz_H(c3,k3) * &
                                       exp(-rjs_H(k3)**2/sigma_ij**2)
                  k4 = 9*(k3-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k4 = k4+1
                      coeff_h_der = 4.d0/sqrt(pi) * 1.d0/sigma_ij**3 * exp(-rjs_H(k3)**2/sigma_ij**2)
                      terms = 0.d0
                      if (c1 == c3) then
                        terms = terms + xyz_H(c2,k3)/rjs_H(k3)**2
                      end if
                      if (c2 == c3) then
                        terms = terms + xyz_H(c1,k3)/rjs_H(k3)**2
                      end if
                      terms = terms -2.d0*(xyz_H(c1,k3)*xyz_H(c2,k3)*xyz_H(c3,k3))/rjs_H(k3)**2 * &
                              (1.d0/rjs_H(k3)**2 + 1.d0/sigma_ij**2)
                      h_func_der(k4) = coeff_h_der * terms
                      if ( rjs_0(k3) < rcut ) then
                        if (a == i2) then
                          b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            T_SR_mult(k3) * a_SCS(3*(q-1)+c2,:)
                          dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) - (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * T_SR_mult(k3)                   
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) + (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * d_mult_i(k3) * &
                                                            neighbor_alpha0(k3)
                        else
                          b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            T_SR_mult(k3) * a_SCS(3*(q-1)+c2,:)
                          dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * T_SR_mult(k3)
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - (f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * d_mult_i(k3) * &
                                                            neighbor_alpha0(k3)
                        end if      
                      else
                        if ( a == i2 ) then
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) + ((f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            neighbor_alpha0(k3)) * d_mult_o(k3)
                        else
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - ((f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            neighbor_alpha0(k3)) * d_mult_o(k3)

                        end if
                      end if
                    end do
                  end do
                end if
              end do
            end do

            if (do_hirshfeld_gradients) then

              allocate( coeff_der(1:9*n_sub_pairs) )
              allocate( coeff_fdamp(1:9*n_sub_pairs) )
            
              coeff_fdamp = 0.d0
              coeff_der = 0.d0

              k3 = 0
              do p = 1, n_sub_sites
                k3 = k3+1
                r_vdw_i = r0_ii(k3)
                s_i = neighbor_sigma(k3)
                i2 = sub_neighbors_list(k3)
                do j2 = 2, n_sub_neigh(p)
                  k3 = k3+1
                  r_vdw_j = r0_ii(k3)
                  s_j = neighbor_sigma(k3)
                  j = sub_neighbors_list(k3)
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  S_vdW_ij = sR*(r_vdw_i + r_vdw_j)
                  exp_term = exp(-d*(rjs_H(k3)/S_vdW_ij - 1.d0))
                  k4 = 9*(k3-1)
                  if ( rjs_0(k3) < rcut ) then
                    q = p_list(k3)
                  end if
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k4 = k4+1
                      coeff_fdamp(k4) = (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                               (-T_func(k4) * g_func(k3) + h_func(k4))
                      dg = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * (rjs_H(k3)/sigma_ij)**2 * &
                           (-rjs_H(k3)/(3.d0 * sigma_ij**3))
                      dh = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * xyz_H(c1,k3)*xyz_H(c2,k3) / &
                           (sigma_ij**5 * rjs_H(k3)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k3)/sigma_ij)**2)
                      coeff_der(k4) = (1.d0-f_damp(k3)) * (-T_func(k4)*dg+dh)
                    end do
                  end do
                end do
              end do

              k3 = 0
              do p = 1, n_sub_sites
                i2 = sub_neighbors_list(k3+1)
                hv_p_der = hirshfeld_v_sub_der(c3,p)*Bohr
                r_vdw_i = r0_ii(k3+1)
                s_i = neighbor_sigma(k3+1)
                do c1 = 1, 3
                  b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_sub_neigh(k3+1)) * &
                                        hirshfeld_v_sub_der(c3,p)*Bohr * a_SCS(3*(p-1)+c1,:) 
                  dB_mat(3*(p-1)+c1,3*(p-1)+c1) = dB_mat(3*(p-1)+c1,3*(p-1)+c1) - 1.d0/(neighbor_alpha0(k3+1) &
                                        * hirshfeld_sub_neigh(k3+1)) * &
                                        hirshfeld_v_sub_der(c3,p)*Bohr
                end do
                do j2 = 2, n_sub_neigh(p)
                  j = sub_neighbors_list(k3+j2)
                  r_vdw_j = r0_ii(k3+j2)
                  s_j = neighbor_sigma(k3+j2)
                  if ( rjs_0(k3+j2) < rcut ) then
                    q = p_list(k3+j2)
                    hv_q_der = hirshfeld_v_sub_der(c3,q)*Bohr
                  else
                    hv_q_der = 0.d0
                  end if
                  k4 = 9*(k3+j2-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k4 = k4+1
                      if ( rjs_0(k3+j2) < rcut ) then                      
                        b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * &
                          a_SCS(3*(q-1)+c2,:) * T_SR_mult(k3+j2)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * T_SR_mult(k3+j2)
                        d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * &
                          neighbor_alpha0(k3+j2) * d_mult_i(k3+j2)
                      else
                        d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * &
                          neighbor_alpha0(k3+j2) * d_mult_o(k3+j2)
                      end if

                    end do
                  end do
                end do
                k3 = k3 + n_sub_neigh(p)
              end do         

              deallocate( coeff_der, coeff_fdamp )
              
            end if
            
            k3 = 0
            do p = 1, n_sub_sites
              k3 = k3+1
              do j2 = 2, n_sub_neigh(p)
                k3 = k3+1
                if ( rjs_0(k3) < rcut ) then
                  q = p_list(k3)
                end if
                k4 = 9*(k3-1)
                do c1 = 1, 3
                  do c2 = 1, 3
                    k4 = k4+1
                    if ( rjs_0(k3) < rcut ) then
                      dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + &
                          dT_SR_mult(k3,c3) * T_SR(k4)
                      b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + dT_SR_mult(k3,c3) * T_SR(k4) * &
                                            a_SCS(3*(q-1)+c2,:) 
                      d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - d_dmult_i(k3,c3) * d_arr_i(k4)
                    else
                      d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - d_dmult_o(k3,c3) * d_arr_o(k4)
                    end if
                  end do
                end do
              end do
            end do

            !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
            !  write(*,*) "dB_mat"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) dB_mat(p,:)
            !  end do
            !end if

            b_der = -b_der+d_der

            !call cpu_time(time1)

            if ( polynomial_expansion ) then
            
              call dgemm( "n", "n", 3*n_sub_sites, 3, 3*n_sub_sites, 1.d0, B_pol, 3*n_sub_sites, b_der, &
                          3*n_sub_sites, 0.d0, da_SCS, 3*n_sub_sites)
            
            else

              call dsytrs( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, b_der, 3*n_sub_sites, work_arr, &
                           12*n_sub_sites, info )

              da_SCS = b_der
            
            end if

            !call cpu_time(time2)

            !write(*,*) "Timing for solving force component", time2-time1
            
            do p = 1, n_sub_sites
              do c1 = 1, 3
                da_iso(p,c3,om) = da_iso(p,c3,om) + da_SCS(3*(p-1)+c1,c1)
              end do
            end do
            da_iso(:,c3,om) = da_iso(:,c3,om)/3.d0
            
            if ( om == 2 ) then
            
              f_damp_der_SCS = 0.d0
              f_damp_der_mbd = 0.d0 ! This is cleared so we can recalculate it with SCS values
              dT_mbd = 0.d0
              dT_LR = 0.d0
              da_mbd = 0.d0
              do_mbd = 0.d0
              k2 = 0
              do p = 1, n_mbd_sites
                k2 = k2+1
                r_vdw_i = r0_ii_SCS(k2)
                i2 = mbd_neighbors_list(k2)
                i1 = modulo(i2-1, n_sites0) + 1
                if ( rjs_0_mbd(k2) < (rcut-r_buffer)/Bohr ) then
                  r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                  da_mbd(k2) = da_iso(r,c3,2)
                  do_pref = -omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                               ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                  do_mbd(k2) = do_pref
                else if ( rjs_0_mbd(k2) .ge. (rcut-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut/Bohr ) then
                  r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                  do_pref = -omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                               ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                  da_mbd(k2) = da_iso(r,c3,2) * &
                              (1.d0 - 3.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2 &
                                    + 2.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**3) + &
                                    a_iso(r,2) * &
                              (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                               + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr)) + &
                               central_pol(i1) * &
                               ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                 - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                  do_mbd(k2) = do_pref * & 
                              (1.d0 - 3.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2 &
                                    + 2.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**3) + &
                                    o_p(r) * &
                              (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                               + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr)) + &
                               central_omega(i1) * &
                               ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                 - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                else if ( rjs_0_mbd(k2) .ge. rcut/Bohr .and. rjs_0_mbd(k2) < (rcut_mbd-r_buffer)/Bohr ) then
                  da_mbd(k2) = 0.d0
                  do_mbd(k2) = 0.d0
                else if ( rjs_0_mbd(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut_mbd/Bohr ) then
                  da_mbd(k2) = central_pol(i1) &
                       * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                       + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                       * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                  do_mbd(k2) = central_omega(i1) &
                       * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                       + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                       * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                end if
                da_i = da_mbd(k2)
                a_mbd_i = a_mbd(k2)
                do j3 = 2, n_mbd_neigh(p)
                  k2 = k2+1
                  j = mbd_neighbors_list(k2)
                  j1 = modulo(j-1, n_sites0) + 1
                  q = p_mbd(k2)
                  if ( rjs_0_mbd(k2) < (rcut-r_buffer)/Bohr ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                    da_mbd(k2) = da_iso(r,c3,2)
                    do_pref = -omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                               ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                  else if ( rjs_0_mbd(k2) .ge. (rcut-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut/Bohr ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                    do_pref = -omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                               ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                    da_mbd(k2) = da_iso(r,c3,2) * &
                                (1.d0 - 3.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2 &
                                      + 2.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**3) + &
                                      a_iso(r,2) * &
                                (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                 + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr)) + &
                                 central_pol(j1) * &
                                  ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                   - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                   * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr)) 
                    do_mbd(k2) = do_pref * &
                                (1.d0 - 3.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2 &
                                      + 2.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**3) + &
                                      o_p(r) * &
                                (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                 + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr)) + &
                                 central_omega(j1) * &
                                  ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                   - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                   * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                  else if ( rjs_0_mbd(k2) .ge. rcut/Bohr .and. rjs_0_mbd(k2) < (rcut_mbd-r_buffer)/Bohr ) then
                    da_mbd(k2) = 0.d0
                    do_mbd(k2) = 0.d0
                  else if ( rjs_0_mbd(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut_mbd/Bohr ) then
                    da_mbd(k2) = central_pol(j1) &
                         * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                         + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                         * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                    do_mbd(k2) = central_omega(j1) &
                         * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                         + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                         * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                  end if
                  da_j = da_mbd(k2)
                  a_mbd_j = a_mbd(k2)
                  r_vdw_j = r0_ii_SCS(k2)
                  R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                  S_vdW_ij = sR*R_vdW_SCS_ij
                  dS_vdW_ij = sR/3.d0 * ( r_vdw_i/a_mbd_i * da_i + &
                                          r_vdw_j/a_mbd_j * da_j )
                  f_damp_der_SCS(k2) = -(d*rjs_mbd(k2))/S_vdW_ij**2 * f_damp_SCS(k2)**2 * &
                                         exp(-d*(rjs_mbd(k2)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                  k3 = 9*(k2-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k3 = k3+1
                      dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                        T_mbd(k3) * f_damp_der_SCS(k2)
                    end do
                  end do
                  if (i == i2 .or. i == j) then
                    f_damp_der_mbd(k2) = d/S_vdW_ij * f_damp_SCS(k2)**2 * &
                                     exp( -d*(rjs_mbd(k2)/S_vdW_ij - 1.d0) ) * xyz_mbd(c3,k2)/rjs_mbd(k2)
                    k3 = 9*(k2-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k3 = k3+1
                        dT_mbd(k3) = (-15.d0 * xyz_mbd(c1,k2) * xyz_mbd(c2,k2) * xyz_mbd(c3,k2))/rjs_mbd(k2)**7
                        if (c1 == c2) then
                          dT_mbd(k3) = dT_mbd(k3) + 3.d0/rjs_mbd(k2)**5 * xyz_mbd(c3,k2)
                        end if
                        if (c2 == c3) then
                          dT_mbd(k3) = dT_mbd(k3) + 3.d0/rjs_mbd(k2)**5 * xyz_mbd(c1,k2)
                        end if
                        if (c1 == c3) then
                          dT_mbd(k3) = dT_mbd(k3) + 3.d0/rjs_mbd(k2)**5 * xyz_mbd(c2,k2)
                        end if
                        if (i == i2) then
                          dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) - &
                                             T_mbd(k3) * f_damp_der_mbd(k2) - &
                                             dT_mbd(k3) * f_damp_SCS(k2)
                        else if (i == j) then
                          dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                             T_mbd(k3) * f_damp_der_mbd(k2) + &
                                             dT_mbd(k3) * f_damp_SCS(k2)
                        end if
                      end do
                    end do
                  end if
                end do
              end do
            
              forces_TS = 0.d0
              k2 = 0
              da_2b = 0.d0
              do_2b = 0.d0
              s = neighbor_species(n_tot+1)
              r_vdw_i = r0_ref(s) / Bohr * (a_iso(1,2)/(alpha0_ref(s)/Bohr**3))**(1.d0/3.d0)
              dr_vdw_i = r_vdw_i / (3.d0 * a_iso(1,2)) * da_iso(1,c3,2)
              do_pref = -omega_ref * (a_iso(1,1) * da_iso(1,c3,2) - a_iso(1,2) * da_iso(1,c3,1)) / &
                               ( a_iso(1,1)**2 * (a_iso(1,2)/a_iso(1,1) - 1.d0)**(3.d0/2.d0) )
              do p = 1, n_2b_sites
                k2 = k2+1
                i2 = sub_2b_list(k2)
                i1 = modulo(i2-1, n_sites0) + 1
                if ( rjs_2b(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_2b(k2) < rcut_mbd/Bohr ) then
                  da_2b(k2) = central_pol(i1) * &
                                ( + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer)) &
                                  - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer))**2) &
                                  * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
                  do_2b(k2) = central_omega(i1) * &
                                ( + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer)) &
                                  - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer))**2) &
                                  * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
                else if ( rjs_2b(k2) .ge. rcut_mbd/Bohr .and. rjs_2b(k2) < (rcut_2b-r_buffer)/Bohr ) then
                  da_2b(k2) = 0.d0
                  do_2b(k2) = 0.d0
                else if ( rjs_2b(k2) .ge. (rcut_2b-r_buffer)/Bohr .and. rjs_2b(k2) < rcut_2b/Bohr ) then
                  da_2b(k2) = central_pol(i1) * &
                         ( - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer)) &
                           + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer))**2) &
                         * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
                  do_2b(k2) = central_omega(i1) * &
                         ( - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer)) &
                           + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer))**2) &
                         * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
                end if
                r_vdw_j = r0_ii_SCS_2b(k2)
                dr_vdw_j = r_vdw_j / (3.d0 * a_2b(k2)) * da_2b(k2)
                dC6_2b = 3.d0/2.d0 * (central_omega(i0)*o_2b(k2) &
                            / (central_omega(i0)+o_2b(k2)) &
                            * (da_iso(1,c3,2)*a_2b(k2) + a_iso(1,2)*da_2b(k2)) &
                            + a_iso(1,2) * a_2b(k2) / (central_omega(i0)+o_2b(k2))**2 &
                            * (do_pref * o_2b(k2)**2 + central_omega(i0)**2 * do_2b(k2)))
                f_damp_der_2b = d * f_damp_SCS_2b(k2)**2 * &
                                        exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (-xyz_2b(c3,k2)/rjs_2b(k2)) &
                                        - rjs_2b(k2)/0.97d0 * 1.d0/(r_vdw_i+r_vdw_j)**2 &
                                        * (dr_vdw_i + dr_vdw_j))
                forces_TS = forces_TS + ( dC6_2b * f_damp_SCS_2b(k2) / rjs_2b(k2)**6 &
                                        + 6.d0/rjs_2b(k2)**8 * xyz_2b(c3,k2) * C6_2b(k2) * f_damp_SCS_2b(k2) &
                                        + C6_2b(k2)/rjs_2b(k2)**6 * f_damp_der_2b )
              end do
              forces_TS = 1.d0/2.d0 * forces_TS

              if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
                write(*,*) "forces_TS", forces_TS
              end if

              G_mat = 0.d0

              do j = 1, n_freq
                k3 = 0
                do p = 1, n_mbd_sites
                  i2 = mbd_neighbors_list(k3+1)
                  G_mat(3*(p-1)+1:3*(p-1)+3,:,j) = G_mat(3*(p-1)+1:3*(p-1)+3,:,j) + &
                    a_mbd(k3+1)/(1.d0 + (omegas_mbd(j)/o_mbd(k3+1))**2) * &
                    dT_LR(3*(p-1)+1:3*(p-1)+3,:) + &
                    da_mbd(k3+1)/(1.d0 + (omegas_mbd(j)/o_mbd(k3+1))**2) * &
                    T_LR(3*(p-1)+1:3*(p-1)+3,:) + &
                    a_mbd(k3+1) * (2.d0 * omegas_mbd(j)**2 * o_mbd(k3+1)) * &
                    do_mbd(k3+1) / ( o_mbd(k3+1)**2 + omegas_mbd(j)**2 )**2 * &
                    T_LR(3*(p-1)+1:3*(p-1)+3,:)
                    !if ( p == 29 .and. i == 1 .and. c3 == 1 .and. om == 2 .and. j == 1 ) then
                    !  write(*,*) "G_mat", G_mat(3*(p-1)+3,3*(21-1)+1,1)
                    !  write(*,*) "a_mbd", a_mbd(k3+1)
                    !  write(*,*) "o_mbd", o_mbd(k3+1)
                    !  write(*,*) "dT_LR", dT_LR(3*(p-1)+3,3*(21-1)+1)
                    !  write(*,*) "da_mbd", da_mbd(k3+1)
                    !  write(*,*) "T_LR", T_LR(3*(p-1)+3,3*(21-1)+1)
                    !  write(*,*) "do_mbd", do_mbd(k3+1)
                    !  write(*,*) "k3+1", k3+1
                    !  write(*,*) "rjs_0_mbd", rjs_0_mbd(k3+1)*Bohr
                    !end if
                  k3 = k3+n_mbd_neigh(p)
                end do
              end do


              !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
              !  write(*,*) "G_mat"
              !  do p = 1, 3*n_sub_sites
              !    write(*,*) G_mat(p,:,1)
              !  end do
              !end if

              if ( n_order > 1 ) then

                integrand = 0.d0
                total_integrand = 0.d0
                do j = 1, n_freq

                  force_series = 0.d0
                  total_energy_series = 0.d0
                  do k2 = 1, n_order-1
                    force_series = force_series + AT_n_f(:,:,k2,j) 
                    total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2,j) 
                  end do
                  k3 = 0
                  do p = 1, n_mbd_sites
                    i2 = mbd_neighbors_list(k3+1)
                    do c1 = 1, 3
                      integrand(j) = integrand(j) + & !1.d0/(1.d0 + (omegas_mbd(j)/0.5d0)**2) * &
                      dot_product(G_mat(3*(p-1)+c1,:,j),force_series(:,3*(p-1)+c1))
                      total_integrand(j) = total_integrand(j) + a_mbd(k3+1) / &
                            (1.d0 + (omegas_mbd(j)/o_mbd(k3+1))**2) &
                            * dot_product(T_LR(3*(p-1)+c1,:), &
                            total_energy_series(:,3*(p-1)+c1))
                    end do
                    k3 = k3 + n_mbd_neigh(p)
                  end do
                end do

                integral = 0.d0
                call integrate("trapezoidal", omegas_mbd, integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
                forces0(c3,i) = forces0(c3,i) + (1.d0/(2.d0*pi) * integral + forces_TS) * Hartree/Bohr

                write(*,*) "integral", integral
                write(*,*) "MBD force", i, c3, forces0(c3,i)
                integral = 0.d0
                if (c3 == 1 ) then
                  call integrate("trapezoidal", omegas_mbd, total_integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
                  write(*,*) "MBD total energy of sphere", i, (integral / (2.d0*pi) + E_TS) * Hartree
                end if

              end if
            
            end if ! om loop(?)

          end do ! c3 loop
        
          deallocate( da_SCS, dT, b_der, g_func_der, h_func_der, d_der, f_damp_der )
          deallocate( dB_mat )

        end if ! do_derivatives
        
        if ( om == 2 ) then
          deallocate( T_LR, r0_ii_SCS, f_damp_SCS, AT, AT_n, energy_series, omegas_mbd, integrand, n_mbd_neigh, &
                      mbd_neighbors_list, p_mbd, r0_ii_mbd, neighbor_alpha0_mbd, xyz_mbd, rjs_mbd, T_mbd, a_mbd, &
                      rjs_0_mbd, xyz_0_mbd, o_mbd, sub_2b_list, xyz_2b, rjs_2b, r0_ii_2b, neighbor_alpha0_2b, &
                      a_2b, o_2b, r0_ii_SCS_2b, f_damp_SCS_2b, C6_2b )
        end if
                    
        if ( do_derivatives .and. om == 2 ) then
          deallocate( da_mbd, AT_n_f, dT_mbd, f_damp_der_mbd, f_damp_der_SCS, dT_LR, G_mat, force_series, &
                      total_energy_series, total_integrand, da_2b, do_2b, do_mbd )
        end if

      end do ! om

      deallocate( n_sub_neigh, sub_neighbors_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, omegas, &
                  T_func, b_i, d_vec, g_func, h_func, f_damp, a_SCS, ipiv, p_list, work_arr, rjs_0, &
                  a_iso, o_p, T_SR, T_SR_mult, d_arr_i, d_arr_o, d_mult_i, d_mult_o, dT_SR_mult, d_dmult_i, &
                  d_dmult_o, hirshfeld_sub_neigh )
      !deallocate( ia, ja, val )

      deallocate( B_mat )
      if ( polynomial_expansion ) then
        deallocate( B_pol )
      end if

      if ( do_derivatives ) then
        deallocate( da_iso )
      end if
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        deallocate( hirshfeld_v_sub_der )
      end if
           
    end do

    call cpu_time(time2)

    write(*,*) "Local energies and forces timing", time2-time1

    write(*,*) "E_MBD", E_MBD
      
    !deallocate( central_pol, central_omega )
    deallocate( hirshfeld_v_cart_der_H )
    !deallocate( alpha_SCS_full )

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
