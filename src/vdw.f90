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
  subroutine get_scs_polarizabilities( hirshfeld_v, hirshfeld_v_cart_der, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, hirshfeld_v_neigh, &
                                       sR, d, c6_ref, r0_ref, alpha0_ref, do_derivatives, alpha_SCS0, dalpha_full, &
                                       c6_scs, r0_scs, alpha0_scs, energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: hirshfeld_v_cart_der(:,:), rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:) !, hirshfeld_v(:), hirshfeld_v_neigh(:) !NOTE: uncomment this in final implementation
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_derivatives
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:), alpha_SCS0(:,:), dalpha_full(:,:), hirshfeld_v(:), &
                             hirshfeld_v_neigh(:), c6_scs(:), r0_scs(:), alpha0_scs(:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), &
                           f_damp(:), neighbor_alpha0(:), &
                           T_func(:), h_func(:), g_func(:), &
                           omegas(:), omega_i(:), &
                           alpha_i(:), sigma_i(:), B_mat(:,:), &
                           rjs_H(:), xyz_H(:,:), work_arr(:), &
                           dT(:), & !dT_SR(:,:)
                           f_damp_der(:), &
                           g_func_der(:), h_func_der(:), &
                           dalpha(:,:,:), &
                           !dT_SR_v(:,:), &
                           !coeff_der(:,:,:,:), coeff_fdamp(:,:,:,:), &
                           coeff_der(:), coeff_fdamp(:), & ! Assuming that coeff_der and coeff_fdamp are symmetric
                           hirshfeld_v_cart_der_H(:,:), I_mat(:,:), &
                           a_vec(:,:), &
                           BTB_reg(:,:), B_reg(:,:), &
                           a_SCS(:,:), da_vec(:,:), vect1(:,:), &
                           vect2(:,:), vect3(:,:), vect4(:,:), vect_temp(:,:), da_SCS(:,:), &
                           c6_nsites(:), b_der(:,:), alpha_SCS_full(:,:,:), dB_mat(:,:), alpha_test(:,:), &
                           test_vector(:,:), neighbor_sigma(:), hirshfeld_v_sub_der(:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, dist, f_damp_SCS_ij, t1, t2, &
              reg_param, sigma_ij, coeff_h_der, dg, dh, s_i, s_j, terms, omega_ref, xyz_i(1:3), xyz_j(1:3), xyz_a(1:3), &
              da_p(1:3,1:3), terms3(1:3,1:3)
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, n_order, n_freq, om, n_tot
    integer :: i, i2, i3, j, j2, j3, k, k2, k3, k4, a, a2, c1, c2, c3, lwork, b, p, q, r, n_count, n_max, k_i, k_j, &
               k_a
    logical :: do_timing = .false., do_hirshfeld_gradients = .true., &
               total_energy = .true., regularization = .false., read_hirshfeld = .false., &
               psblas = .false., polynomial_expansion = .true.
               
!    LOCAL TEST stuff:
    integer, allocatable :: p_to_i(:), i_to_p(:), sub_neighbors_list(:), n_sub_neigh(:), p_list(:)
    logical, allocatable :: in_cutoff(:), in_force_cutoff(:)
    integer :: n_sub_sites, n_sub_pairs, n_tot2, s, n_degree
    logical :: local = .true.

!   MBD stuff:
    real*8, allocatable :: T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), AT(:,:,:), AT_n(:,:,:,:), energy_series(:,:), &
                           integrand(:), AT_n_f(:,:,:,:), f_damp_der_SCS(:), dT_LR(:,:), G_mat(:,:,:), force_series(:,:), &
                           VL(:,:), total_energy_series(:,:), alpha_grad(:,:), total_integrand(:), B_inv(:,:), rjs_0(:), &
                           T_mbd(:), r0_ii_mbd(:), neighbor_alpha0_mbd(:), omegas_mbd(:), rjs_0_mbd(:), xyz_0_mbd(:,:), &
                           xyz_mbd(:,:), rjs_mbd(:), d_der(:,:), dT_mbd(:), f_damp_der_mbd(:), a_mbd(:), da_mbd(:), &
                           a_iso(:,:), o_p(:), central_pol(:), da_iso(:,:,:), central_omega(:), o_mbd(:), sub_2b_list(:), &
                           xyz_2b(:,:), rjs_2b(:), r0_ii_2b(:), neighbor_alpha0_2b(:), f_damp_SCS_2b(:), &
                           a_2b(:), r0_ii_SCS_2b(:), C6_2b(:), da_2b(:), T_SR(:), T_SR_mult(:), d_arr_i(:), d_arr_o(:), &
                           d_mult_i(:), d_mult_o(:), dT_SR_mult(:,:), d_dmult_i(:,:), d_dmult_o(:,:), do_mbd(:), &
                           hirshfeld_sub_neigh(:)
    real*8 :: rcut_mbd, xyz_l(1:3), a_mbd_i, a_mbd_j, da_i, da_j, pol1, E_TS, rcut_2b, r_buffer, f_damp_der_2b, dr_vdw_i, &
              dr_vdw_j, forces_TS, dC6_2b, mult1_i, mult1_j, mult2, dmult1_i(1:3), dmult1_j(1:3), dmult2(1:3), hv_p_der, &
              hv_q_der, do_pref
    logical :: derivative_expansion = .false.
    integer :: n_mbd_sites, n_mbd_pairs, l, l_cent, n_2b_sites
    integer, allocatable :: n_mbd_neigh(:), mbd_neighbors_list(:), p_mbd(:)

!   DERIVATIVE TEST stuff:
    !real*8 :: polyfit(1:7)
    !real*8 :: polyfit(1:12)
    real*8 :: polyfit(1:15)
    real*8, allocatable :: eigval(:) 

!    PSBLAS stuff:
    type(psb_ctxt_type) :: icontxt
    integer(psb_ipk_) ::  iam, np, ip, jp, idummy, nr, nnz, info_psb
    type(psb_desc_type) :: desc_a
    type(psb_dspmat_type) :: A_sp
    !real*8, allocatable :: x_vec(:,:), b_vec(:,:)
    type(psb_d_vect_type) :: x_vec, b_vec
    integer(psb_lpk_), allocatable :: ia(:), ja(:), myidx(:)
    real(psb_dpk_), allocatable :: val(:), val_xv(:,:), val_bv(:,:), b_i(:,:), d_vec(:,:), A_i(:,:)
    type(psb_dprec_type) :: prec
    character(len=20) :: ptype


    !write(*,*) "neighbors_list"
    !write(*,*) neighbors_list

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
    
    n_degree = size(polyfit)-1
    !n_degree = 3

    !open(unit=89, file="hirshfeld.dat", status="new")
    !do i = 1, n_sites
    !  write(89,*) hirshfeld_v(i)
    !end do
    !close(89)


    ! HACK FOR HIRSHFELD DERIVATIVES
    allocate( hirshfeld_v_cart_der_H(1:3, n_pairs) )

    hirshfeld_v_cart_der_H = 0.d0

    write(*,*) "dv"
    k = 0
    do i = 1, n_sites
      do j2 = 1, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        n_tot = sum(n_neigh(1:j)) - n_neigh(j)
        a = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(j)),i,1) 
        hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,n_tot+a)
        !write(*,*) "i, j, der, der_H", i, j, hirshfeld_v_cart_der(:,k), hirshfeld_v_cart_der_H(:,k)
      end do
    end do

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
    !omega_ref = (4.d0 * c6_ref(1)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(1)/Bohr**3)**2)
    pi = acos(-1.d0)

    write(*,*) "rcut", rcut/Bohr

!   This should allow to only take a subset of atoms for parallelization:

!   Number of frequencies
    n_freq = size(alpha_SCS0, 2)
    

!   LOCAL TEST:
    if ( local ) then
      
      !open(unit=69, file="alpha_SCS_local_inv.dat", status="new")
      !open(unit=79, file="alpha_SCS_local_pol.dat", status="new")
      !write(*,*) "hirshfeld_v"
      !do i = 1, n_sites
      !  write(*,*) hirshfeld_v(i), ", &"
      !end do
      write(*,*) "Starting local calculation"

      E_MBD = 0.d0

      alpha_SCS0 = 0.d0
      !allocate( sigma_i(1:n_sites) )
      !allocate( alpha_i(1:n_sites) )
      allocate( alpha_SCS_full(1:3*n_sites,1:3,1:2) )
      !write(*,*) "test 1"
      if ( do_derivatives ) then
        dalpha_full = 0.d0
        !allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
        !do k = 1, n_pairs
        !  hirshfeld_v_cart_der_H(:,k) = hirshfeld_v_cart_der(:,k)*Bohr
        !end do
      end if
      !write(*,*) "test 2"
      !write(*,*) "hirshfeld der"
      !do p = 1, n_neigh(1)
      !  write(*,*) neighbors_list(p), hirshfeld_v_cart_der_H(1,p)
      !end do
      !write(*,*) "the other"
      !do p = sum(n_neigh(1:2))-n_neigh(2)+1, sum(n_neigh(1:2))
      !  write(*,*) neighbors_list(p), hirshfeld_v_cart_der_H(1,p)
      !end do

      alpha_SCS_full = 0.d0
      
      !k2 = 1
      !do i = 1, n_sites
      !  alpha_i(i) = alpha0_ref(1) / Bohr**3 * hirshfeld_v(i)
      !  sigma_i(i) = (sqrt(2.d0/pi) * alpha_i(i)/3.d0)**(1.d0/3.d0)
      !  k2 = k2+n_neigh(i)
      !end do
      
      n_max = maxval(neighbors_list)
    
      allocate( central_pol(1:n_sites) )
      central_pol = 0.d0
      allocate( central_omega(1:n_sites) )
      central_omega = 0.d0
      allocate( in_cutoff(1:n_max) )
      allocate( p_to_i(1:n_max) )
      allocate( i_to_p(1:n_max) )
      allocate( A_i(1:3*n_max,1:3*n_sites) )
      !write(*,*) "test 3"
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        allocate( in_force_cutoff(1:n_max) )
      end if
      !write(*,*) "test 4"
      A_i = 0.d0 ! This is a temporary solution to store a_SCS's for each atom
      r_buffer = 0.d0

      do i = 1, n_sites      

        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        s = neighbor_species(n_tot+1)
        omega_ref = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
        n_sub_sites = 0
        n_sub_pairs = 0
        k_i = 0
        do i3 = 1, n_neigh(i)
          !k = k+1
          k_i = k_i + 1
          !j = p_to_i(p)
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
        !allocate( a_vec(1:3*n_sub_sites,1:3) )
        allocate( a_SCS(1:3*n_sub_sites,1:3) )
        allocate( ipiv(1:3*n_sub_sites) )
        allocate( work_arr(1:12*n_sub_sites) )
        allocate( rjs_0(1:n_sub_pairs) )
        allocate( a_iso(1:n_sub_sites,1:2) )
        allocate( o_p(1:n_sub_sites) )
        if ( do_derivatives .and. do_hirshfeld_gradients ) then
          allocate( hirshfeld_v_sub_der(1:3,1:n_sub_sites) )
          hirshfeld_v_sub_der = 0.d0
        end if
        
        allocate( ia(1:9*n_sub_pairs) )
        allocate( ja(1:9*n_sub_pairs) )
        allocate( val(1:9*n_sub_pairs) )
        allocate( b_i(1:3*n_sub_sites,1:3) )
        allocate( d_vec(1:3*n_sub_sites,1:3) )
        
        a_iso = 0.d0

        !write(*,*) "allocation successful"

        do om = 1, 2

        xyz_H = 0.d0
        rjs_H = 0.d0
        r0_ii = 0.d0
        neighbor_alpha0 = 0.d0
        neighbor_sigma = 0.d0

        f_damp = 0.d0
        g_func = 0.d0
        h_func = 0.d0

        nnz = 0
        k2 = 0
        T_func = 0.d0
        B_mat = 0.d0

        a_SCS = 0.d0
        b_i = 0.d0
        !a_iso = 0.d0

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
          !if (i == 1 .and. i3 == 8) then
          !write(*,*) "k2", k2
          !end if
          k_i = k_i+1
          i2 = neighbors_list(n_tot+k_i)
          if ( rjs(n_tot+k_i) < rcut ) then
            p = p+1
            !p = i_to_p(i2)
            !write(*,*) "p, i_to_p", p, i_to_p(i2) 
            k2 = k2+1
            rjs_0(k2) = rjs(n_tot+k_i)
            !n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
            !s = neighbor_species(n_tot2+1)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = neighbors_list(n_tot+k_i)
            if ( do_derivatives .and. do_hirshfeld_gradients ) then
              hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            !r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot2+1)**(1.d0/3.d0)
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              !neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot2+1)
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              !neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot2+1)) / &
              !  (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
            end if
            neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
            !xyz_H(:,k2) = xyz(:,n_tot2+1)/Bohr
            !rjs_H(k2) = rjs(n_tot2+1)/Bohr
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
                nnz = nnz+1
                ia(nnz) = 3*(p-1)+c1
                ja(nnz) = 3*(p-1)+c1
                val(nnz) = 1.d0/neighbor_alpha0(k2)
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
                ! TWO-BODY TEST
                !write(*,*) "q, i_to_p", q, i_to_p(neighbors_list(n_tot+k_j))
                if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                !j = neighbors_list(n_tot+k_j)
                !if ( rjs(n_tot+k_j) < 2*rcut ) then
                  !if (i == 1 .and. i3 == 7) then
                  !write(*,*) "i, j", i2, j
                  !end if
                  n_sub_neigh(p) = n_sub_neigh(p) + 1
                  j = neighbors_list(n_tot+k_j)
                  !q = i_to_p(j)
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
                  !rjs_H(k2) = rjs(n_tot2+j3)/Bohr
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
                        !if ( rjs(n_tot+k_j) < rcut-r_buffer ) then)
                        if ( rjs(n_tot+k_j) .ge. rcut-r_buffer ) then
                          mult1_j = mult1_j * &
                                         (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**2 &
                                               + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/(r_buffer))**3 )
                        end if
                        if ( rjs_H(k2) .ge. (rcut-r_buffer)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                          mult2 = mult2 * &
                                         (1.d0 - 3.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**2 &
                                               + 2.d0 * ((rjs_H(k2)-(rcut-r_buffer)/Bohr)/(r_buffer/Bohr))**3 )
                        !else if ( rjs_H(k2) > (rcut-r_buffer)/Bohr ) then
                        !  mult2 = 0.d0
                        end if
                        !write(*,*) "mult1, mult2", mult1, mult2
                        B_mat(3*(p-1)+c1,3*(q-1)+c2) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * mult1_i * mult1_j * mult2
                        if ( p == 1 ) then
                          b_i(3*(q-1)+c2,c1) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                g_func(k2) + h_func(k3)) * mult1_i * mult1_j * mult2
                        end if
                        nnz = nnz+1
                        ia(nnz) = 3*(p-1)+c1
                        ja(nnz) = 3*(q-1)+c2
                        val(nnz) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                   g_func(k2) + h_func(k3)) * mult1_j * mult2
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
                        !else if ( rjs_H(k2) > (rcut-r_buffer)/Bohr ) then
                        !  mult2 = 0.d0
                        end if
                        !write(*,*) "d mult1, mult2", mult1, mult2
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
        
        a_SCS = d_vec

        call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, a_SCS, 3*n_sub_sites, work_arr, &
                    12*n_sub_sites, info )

        !if ( i == 1 ) then
        !  write(*,*) "a_SCS"
        !  k2 = 0
        !  do p = 1, n_sub_sites
        !    write(*,*) sub_neighbors_list(k2+1), 1.d0/3.d0 * ( a_SCS(3*(p-1)+1,1) + a_SCS(3*(p-1)+2,2) + &
        !                                                       a_SCS(3*(p-1)+3,3) )
        !    k2 = k2+n_sub_neigh(p)
        !  end do
        !end if

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
          write(*,*) "central polarizability, om", i, om, central_pol(i)
        else
          central_pol(i) = 0.d0
          do c1 = 1, 3
            central_pol(i) = central_pol(i) + a_SCS(c1,c1)
          end do
          central_pol(i) = central_pol(i)/3.d0
          write(*,*) "central polarizability, om", i, om, central_pol(i)
          central_omega(i) = 2.d0*omega_ref/sqrt(central_pol(i)/pol1-1.d0)
          !write(*,*) "central omega", i, central_omega(i)
        end if

        end do
        
        deallocate( sub_neighbors_list, n_sub_neigh, p_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, &
                    omegas, T_func, B_mat, f_damp, g_func, h_func, a_SCS, ipiv, work_arr, rjs_0, a_iso, o_p )
        if ( do_derivatives .and. do_hirshfeld_gradients ) then
          deallocate( hirshfeld_v_sub_der )
        end if
        
        deallocate( ia, ja, val, b_i, d_vec )        
        
       ! end do

      end do  

      !k = 0
      do i = 1, n_sites      

        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !write(*,*) i, "/", n_sites
        !p_to_i = 0
        !i_to_p = 0
        !in_cutoff = .false.
        !k = k+1
        s = neighbor_species(n_tot+1)
        omega_ref = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
        !p = 1
        !p_to_i(p) = i
        !i_to_p(i) = p
        !in_cutoff(i) = .true.
        !do j2 = 2, n_neigh(i)
        !  k = k+1
        !  j = neighbors_list(k)
        !  if (rjs(k) < 2*rcut) then
        !    in_cutoff(j) = .true.
        !    p = p+1
        !    p_to_i(p) = j
        !    i_to_p(j) = p
        !  end if
        !end do
        !write(*,*) "in cutoff"
        !write(*,*) in_cutoff
        !write(*,*) "p_to_i"
        !write(*,*) p_to_i
        !write(*,*) "i_to_p"
        !write(*,*) i_to_p
        !n_sub_sites = p
        !write(*,*) "n_sub_sites", n_sub_sites
        n_sub_sites = 0
        n_sub_pairs = 0
        !n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !allocate( n_sub_neigh(1:n_sub_sites) )
        !n_sub_neigh = 0
        k_i = 0
        do i3 = 1, n_neigh(i)
          !k = k+1
          k_i = k_i + 1
          !j = p_to_i(p)
          if (rjs(n_tot+k_i) < rcut ) then
            n_sub_sites = n_sub_sites + 1
            n_sub_pairs = n_sub_pairs + 1
            !n_sub_neigh(p) = n_sub_neigh(p) + 1
            !n_tot2 = sum(n_neigh(1:j))-n_neigh(j)
            !write(*,*) "n_tot2", n_tot2
            xyz_i = xyz(:,n_tot+k_i)/Bohr
            k_j = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j + 1
              if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
                if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                  !if (i == 1 .and. p == 7) then
                  !write(*,*) "i, j", j, neighbors_list(n_tot2+j3)
                  !end if
                  !n_sub_neigh(p) = n_sub_neigh(p) + 1
                  n_sub_pairs = n_sub_pairs + 1
                end if
                end if
              end if
            end do
          end if
        end do

        !write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_sub_pairs", n_sub_pairs
     
        !write(*,*) "p_to_i", p_to_i
        !write(*,*) "i_to_p", i_to_p
        !write(*,*) "n_sub_pairs", n_sub_pairs
        !write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_sub_neigh", n_sub_neigh
        
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
        !allocate( a_vec(1:3*n_sub_sites,1:3) )
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
        
        allocate( ia(1:9*n_sub_pairs) )
        allocate( ja(1:9*n_sub_pairs) )
        allocate( val(1:9*n_sub_pairs) )
        allocate( b_i(1:3*n_sub_sites,1:3) )
        allocate( d_vec(1:3*n_sub_sites,1:3) )
        
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

        !write(*,*) "allocation successful"
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

        nnz = 0
        k2 = 0
        T_func = 0.d0
        B_mat = 0.d0

        a_SCS = 0.d0
        b_i = 0.d0
        !a_iso = 0.d0

        d_vec = 0.d0
        do p = 1, n_sub_sites
          do c1 = 1, 3
            d_vec(3*(p-1)+c1,c1) = 1.d0
          end do
        end do

        !do p = 1, n_sub_sites
        !  do c1 = 1, 3
        !    !B_mat(j2,j2) = 1.d0
        !    !B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/alpha_i(p_to_i(p))
        !    nnz = nnz+1
        !    ia(nnz) = 3*(p-1)+c1
        !    ja(nnz) = 3*(p-1)+c1
        !    !val(nnz) = 1.d0
        !    val(nnz) = 1.d0/alpha_i(p_to_i(p))
        !  end do
        !end do


        !n_tot = sum(n_neigh(1:i))-n_neigh(i)

        p_list = 0
        n_sub_neigh = 0
        !do j2 = 1, n_neigh(i)
        !  i2 = modulo(neighbors_list(n_tot+j2)-1,n_sites)+1

        !NEIGHBORS_LIST TEST
        !do p = 1, n_sub_sites
        !  i2 = p_to_i(p)
        !write(*,*) "n_neigh", n_neigh(i)
        k_i = 0
        p = 0

        do i3 = 1, n_neigh(i)
          !if (i == 1 .and. i3 == 8) then
          !write(*,*) "k2", k2
          !end if
          k_i = k_i+1
          i2 = neighbors_list(n_tot+k_i)
          if ( rjs(n_tot+k_i) < rcut ) then
            p = p+1
            !p = i_to_p(i2)
            !write(*,*) "p, i_to_p", p, i_to_p(i2) 
            k2 = k2+1
            rjs_0(k2) = rjs(n_tot+k_i)
            hirshfeld_sub_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
            !n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
            !s = neighbor_species(n_tot2+1)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = neighbors_list(n_tot+k_i)
            if ( do_derivatives .and. do_hirshfeld_gradients ) then
              hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            !r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot2+1)**(1.d0/3.d0)
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              !neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot2+1)
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              !neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot2+1)) / &
              !  (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (2.d0 * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
            end if
            neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
            !xyz_H(:,k2) = xyz(:,n_tot2+1)/Bohr
            !rjs_H(k2) = rjs(n_tot2+1)/Bohr
            xyz_H(:,k2) = xyz(:,n_tot+k_i)/Bohr
            xyz_i = xyz_H(:,k2)
            rjs_H(k2) = rjs(n_tot+k_i)/Bohr
            r_vdw_i = r0_ii(k2)
            s_i = neighbor_sigma(k2)
            !T_SR(k2) = 1.d0/neighbor_alpha0(k2)
            if (p == 1) then
              do c1 = 1, 3
                b_i(c1,c1) = 1.d0/neighbor_alpha0(k2)
              end do
            end if
            k3 = 9*(k2-1)
            do c1 = 1, 3
              !if ( rjs(n_tot+k_i) < rcut ) then
              do c2 = 1, 3
                k3 = k3+1
                if ( c1 == c2 ) then
                B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/neighbor_alpha0(k2)
                nnz = nnz+1
                ia(nnz) = 3*(p-1)+c1
                ja(nnz) = 3*(p-1)+c1
                val(nnz) = 1.d0/neighbor_alpha0(k2)
              !else
              !  B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/(neighbor_alpha0(k2))
              !  nnz = nnz+1
              !  ia(nnz) = 3*(p-1)+c1
              !  ja(nnz) = 3*(p-1)+c1
              !  val(nnz) = 1.d0/(neighbor_alpha0(k2))
              !end if
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
              !if ( i2 == i ) then
                dmult1_i = (- 6.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer)) &
                            + 6.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/(r_buffer))**2 ) &
                           * ( -xyz(:,n_tot+k_i)/rjs(n_tot+k_i)/(r_buffer/Bohr))
              !end if
            end if
            k_j = 0
            q = 0
            !do j3 = 2, n_neigh(i2)
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              if ( rjs(n_tot+k_j) < 2*rcut ) then !2*rcut
                if ( rjs(n_tot+k_j) < rcut ) then
                  q = q+1
                end if
                ! TWO-BODY TEST
                !write(*,*) "q, i_to_p", q, i_to_p(neighbors_list(n_tot+k_j))
                if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                !j = neighbors_list(n_tot+k_j)
                !if ( rjs(n_tot+k_j) < 2*rcut ) then
                  !if (i == 1 .and. i3 == 7) then
                  !write(*,*) "i, j", i2, j
                  !end if
                  n_sub_neigh(p) = n_sub_neigh(p) + 1
                  j = neighbors_list(n_tot+k_j)
                  !q = i_to_p(j)
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
                  !rjs_H(k2) = rjs(n_tot2+j3)/Bohr
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
                    !if ( rjs(n_tot+k_j) < rcut-r_buffer ) then)
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
                    !d_mult_i(k2) = 0.d0
                    !d_dmult_i(k2,:) = 0.d0
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
                        nnz = nnz+1
                        ia(nnz) = 3*(p-1)+c1
                        ja(nnz) = 3*(q-1)+c2
                        val(nnz) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                   g_func(k2) + h_func(k3)) * T_SR_mult(k2)
                        d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                  g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                  d_mult_i(k2)
                        d_arr_i(k3) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                  g_func(k2) + h_func(k3)) * neighbor_alpha0(k2)
                      else
                        !write(*,*) "d mult1, mult2", mult1, mult2
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

       if ( i == 1 .and. om == 1 ) then
         write(*,*) "d_vec"
         do p = 1, 3*n_sub_sites
           write(*,*) d_vec(p,:)
         end do
       end if

       if ( .false. ) then

         if (om == 1) then
polyfit = (/ 3.237385145550585d+02, -4.241125470183307d+04, 3.008572712845031d+06, -1.309430416378132d+08, &
3.756106046665028d+09, -7.433108326602138d+10, 1.045457946646248d+12, &
-1.064278184563909d+13, 7.908875986032283d+13, -4.286093180295281d+14, &
1.673744363742281d+15, -4.582925299269683d+15, 8.343821283398570d+15, &
-9.066835011532402d+15, 4.447833222479864d+15 /)

         else
polyfit = (/ 3.464569029392560e+01, -5.541785287730104e+02, 5.429135883990769e+03, -3.643266484189429e+04, &
1.774056111172961e+05, -6.476342669175469e+05, 1.805107059718787e+06, &
-3.873561223705132e+06, 6.400647427099578e+06, -8.079084665579021e+06, &
7.651382671329762e+06, -5.264071084777111e+06, 2.484267601324048e+06, &
-7.192678344511647e+05, 9.633392189452502e+04 /)
         end if

!call cpu_time(time1)
        allocate( myidx(1:3*n_sub_sites) )
        allocate( val_xv(1:3*n_sub_sites,1:3) )
        allocate( val_bv(1:3*n_sub_sites,1:3) )
        val_xv = 0.d0
        !val_bv = 0.d0
        k2 = 0
        do c1 = 1, 3
          do p = 1, n_sub_sites
            k2 = k2+1
            myidx(k2) = k2
        !    i2 = p_to_i(p)
        !    val_xv(3*(p-1)+c1,c1) = 1.d0
        !    val_bv(3*(p-1)+c1,c1) = 1.d0
          end do
        end do
        
        !write(*,*) "i according to p = 1"
        !write(*,*) p_to_i(1)

        !write(*,*) "val_xv"
        !do p = 1, 6
        !  write(*,*) val_xv(p,:)
        !end do


        call cpu_time(time1)
        call psb_init(icontxt)
        a_SCS = polyfit(2)*b_i
        !write(*,*) "polyfit(n_degree+1)"
        !write(*,*) polyfit(n_degree+1)
        !write(*,*) "a_SCS"
        !write(*,*) a_SCS
        call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
          !write(*,*) "cdall", info_psb
        call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
          !write(*,*) "spall", info_psb
          !call psb_geall(x_vec,desc_a,info_psb)
          !write(*,*) "geall x", info_psb
          !call psb_geall(b_vec,desc_a,info_psb)
          !write(*,*) "geall b", info_psb
        call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
          !write(*,*) "spins", info_psb
          !call psb_geins(3*n_sites, myidx, val_xv(:,c1), x_vec, desc_a, info_psb)
          !write(*,*) "x_vec", info_psb
          !call psb_geins(3*n_sites, myidx, val_bv(:,c1), b_vec, desc_a, info_psb)
          !write(*,*) "b_vec", info_psb
        call psb_cdasb(desc_a, info_psb)
          !write(*,*) "cdasb", info_psb
        call psb_spasb(A_sp, desc_a, info_psb)
        call cpu_time(time1)
        do k2 = 3, n_degree+1
          call psb_spmm(1.d0, A_sp, b_i, 0.d0, val_xv, desc_a, info_psb, 'T')
          !write(*,*) "psb_spmm", val_xv
          a_SCS = a_SCS + polyfit(k2) * val_xv
          b_i = val_xv 
          !write(*,*) "k2, polyfit", k2, polyfit(k2)
          !write(*,*) "a_SCS"
          !write(*,*) a_SCS
          !val_bv = val_xv
        end do
        call cpu_time(time2)
        !write(*,*) "psb_spmm timing", time2-time1

        !allocate( d_vec(1:3,1:3*n_sub_sites) )
        !d_vec = 0.d0
        !do p = 1, n_sub_sites
        !  do c1 = 1, 3
        !    d_vec(c1,3*(p-1)+c1) = 1.d0
        !  end do
        !end do
        do c1 = 1, 3
          do c2 = 1, 3
            alpha_SCS_full(3*(i-1)+c2,c1,om) = dot_product(d_vec(c1,:),a_SCS(:,c2))
            if ( c1 == c2) then
              alpha_SCS_full(3*(i-1)+c2,c1,om) = alpha_SCS_full(3*(i-1)+c2,c1,om) + polyfit(1)
            end if
          end do
        end do
      end if

        ! This is the exact solution:
        a_SCS = d_vec

        !if ( i == 17 ) then
        !  write(*,*) "B_mat"
        !  do p = 1, 3*n_sub_sites
        !    write(*,*) B_mat(p,:)
        !  end do
        !  write(*,*) "d_vec"
        !  do p = 1, 3*n_sub_sites
        !    write(*,*) d_vec(p,:)
        !  end do
        !end if

        call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, a_SCS, 3*n_sub_sites, work_arr, &
                    12*n_sub_sites, info )

        !if ( i == 1 ) then
        !  write(*,*) "a_SCS"
        !  k2 = 0
        !  do p = 1, n_sub_sites
        !    write(*,*) sub_neighbors_list(k2+1), 1.d0/3.d0 * ( a_SCS(3*(p-1)+1,1) + a_SCS(3*(p-1)+2,2) + &
        !                                                       a_SCS(3*(p-1)+3,3) )
        !    k2 = k2+n_sub_neigh(p)
        !  end do
        !end if

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
        
          if ( i == 1 ) then
          write(*,*) "a_iso"
          k2 = 0
          do p = 1, n_sub_sites
            write(*,*) sub_neighbors_list(k2+1), a_iso(p,om)
            k2 = k2+n_sub_neigh(p)
          end do
          end if
          !write(*,*) i, a_iso(1,2), neighbor_alpha0(1)
          !if ( a_iso(1,2) > neighbor_alpha0(1) ) then
          !  write(*,*) "#####################################################################################"
          !end if
          !write(*,*) "omega SCS", o_p(1)
          !write(*,*) "omega average", sum(o_p)/n_sub_sites

        do c1 = 1, 3
          do c2 = 1, 3
            alpha_SCS_full(3*(i-1)+c1,c2,om) = a_SCS(c1,c2)
          end do
        end do

        alpha_SCS0(i,om) = 0.d0
        do c1 = 1, 3
          alpha_SCS0(i,om) = alpha_SCS0(i,om) + a_SCS(c1,c1)
        end do
        alpha_SCS0(i,om) = alpha_SCS0(i,om)/3.d0
        ! Exact solution ends
        
        if ( om == 2 ) then

        ! MBD for local polarizabilities:
        ! At least for now: rcut <= rcut_mbd <= rcut_2b
        rcut_mbd = 4.5d0
        rcut_2b = 4.5d0
        !r_buffer = 0.5d0
        n_mbd_sites = 0
        n_mbd_pairs = 0
        
        k_i = 0
        do i3 = 1, n_neigh(i)
          !k = k+1
          k_i = k_i + 1
          !j = p_to_i(p)
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
          !k = k+1
          k_i = k_i + 1
          if (rjs(n_tot+k_i) < rcut_2b .and. rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer) then
            n_2b_sites = n_2b_sites + 1
          end if
        end do

        n_order = 4
        n_freq = 12

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
          allocate( VL(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( total_energy_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( total_integrand(1:n_freq) )
          allocate( do_mbd(1:n_mbd_pairs) )
          allocate( alpha_grad(1:3,1:n_mbd_sites) )
          alpha_grad = 0.d0
        end if
        !if ( do_derivatives ) then
        !  allocate( AT_n_f(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_order-1) )
        !  allocate( dT_mbd(1:9*n_mbd_pairs) )
        !  allocate( f_damp_der_mbd(1:n_mbd_pairs) )
        !  allocate( f_damp_der_SCS(1:n_mbd_pairs) )
        !  allocate( dT_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        !  allocate( G_mat(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        !  allocate( force_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        !  allocate( VL(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        !  allocate( total_energy_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        !  allocate( total_integrand(1:n_freq) )
        !  allocate( alpha_grad(1:3,1:n_mbd_sites) )
        !  alpha_grad = 0.d0
        !end if
        
        allocate( sub_2b_list(1:n_2b_sites) )
        allocate( xyz_2b(1:3,1:n_2b_sites) )
        allocate( rjs_2b(1:n_2b_sites) )
        allocate( r0_ii_2b(1:n_2b_sites) )
        allocate( neighbor_alpha0_2b(1:n_2b_sites) )
        allocate( a_2b(1:n_2b_sites) )
        allocate( r0_ii_SCS_2b(1:n_2b_sites) )
        allocate( f_damp_SCS_2b(1:n_2b_sites) )
        allocate( C6_2b(1:n_2b_sites) )
        if ( do_derivatives ) then
          allocate( da_2b(1:n_2b_sites) )
        end if
        
        a_mbd = 0.d0
        
        omegas_mbd = 0.d0
        omega = 0.d0
        do i2 = 1, n_freq
          omegas_mbd(i2) = omega
          omega = omega + 0.2d0
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

        write(*,*) "o_mbd"

        k2 = 0
        k_i = 0
        p = 0
        do i3 = 1, n_neigh(i)
          k_i = k_i+1
          i2 = neighbors_list(n_tot+k_i)
          if ( rjs(n_tot+k_i) < rcut_mbd ) then
            p = p+1
            k2 = k2+1
            s = neighbor_species(n_tot+k_i)
            mbd_neighbors_list(k2) = neighbors_list(n_tot+k_i)
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
              if ( i == 1 ) then
                !write(*,*) a_iso(r,1)
                write(*,*) o_mbd(k2)
              end if
            else if ( rjs(n_tot+k_i) < rcut .and. rjs(n_tot+k_i) .ge. rcut-r_buffer ) then
              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
              a_mbd(k2) = a_iso(r,2) * &
                          (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3) + &
                          central_pol(i2) * & !neighbor_alpha0_mbd(k2) * &
                                ( + 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                - 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3)
              o_mbd(k2) = o_p(r) * &
                          (1.d0 - 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                + 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3) + &
                          central_omega(i2) * & !neighbor_alpha0_mbd(k2) * &
                                ( + 3.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**2 &
                                - 2.d0 * ((rjs(n_tot+k_i)-rcut+r_buffer)/r_buffer)**3)
              !write(*,*) "a_mbd", a_mbd(k2)
            else if ( rjs(n_tot+k_i) .ge. rcut .and. rjs(n_tot+k_i) < rcut_mbd-r_buffer) then
              a_mbd(k2) = central_pol(i2)
              o_mbd(k2) = central_omega(i2)
            else if ( rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer .and. rjs(n_tot+k_i) < rcut_mbd ) then
              a_mbd(k2) = central_pol(i2) * (1.d0 & !neighbor_alpha0_mbd(k2) * (1.d0 &
                                  - 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                  + 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
              o_mbd(k2) = central_omega(i2) * (1.d0 & !neighbor_alpha0_mbd(k2) * (1.d0 &
                                  - 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                  + 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
              !write(*,*) "a_mbd larger", a_mbd(k2)
            end if
            !o_mbd(k2) = central_omega(i2)
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
                    central_pol(j) * & !neighbor_alpha0_mbd(k2) * &
                                ( + 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                - 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3)
                    o_mbd(k2) = o_p(r) * &
                          (1.d0 - 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                + 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3) + &
                    central_omega(j) * & !neighbor_alpha0_mbd(k2) * &
                                ( + 3.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**2 &
                                - 2.d0 * ((rjs(n_tot+k_j)-rcut+r_buffer)/r_buffer)**3)
                  else if ( rjs(n_tot+k_j) .ge. rcut .and. rjs(n_tot+k_j) < rcut_mbd-r_buffer) then
                    a_mbd(k2) = central_pol(j)
                    o_mbd(k2) = central_omega(j)
                  else
                    a_mbd(k2) = central_pol(j) * (1.d0 & !neighbor_alpha0_mbd(k2) * (1.d0 &
                                        - 3.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                        + 2.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**3)
                    o_mbd(k2) = central_omega(j) * (1.d0 & !neighbor_alpha0_mbd(k2) * (1.d0 &
                                        - 3.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                        + 2.d0 * ((rjs(n_tot+k_j)-rcut_mbd+r_buffer)/(r_buffer))**3)
                   !write(*,*) "a_mbd larger", a_mbd(k2)
                  end if
                  !o_mbd(k2) = central_omega(j)
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
          if ( rjs(n_tot+k_i) < rcut_2b .and. rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer ) then
            k2 = k2+1
            s = neighbor_species(n_tot+k_i)
            sub_2b_list(k2) = neighbors_list(n_tot+k_i)
            r0_ii_2b(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            xyz_2b(:,k2) = xyz(:,n_tot+k_i)/Bohr
            xyz_i = xyz_2b(:,k2)
            rjs_2b(k2) = rjs(n_tot+k_i)/Bohr
            neighbor_alpha0_2b(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            if ( rjs(n_tot+k_i) .ge. rcut_mbd-r_buffer .and. rjs(n_tot+k_i) < rcut_mbd ) then
              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
              a_2b(k2) = central_pol(i2) * &
                                ( + 3.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**2 &
                                - 2.d0 * ((rjs(n_tot+k_i)-rcut_mbd+r_buffer)/(r_buffer))**3)
              !write(*,*) "a_mbd", a_mbd(k2)
            else if ( rjs(n_tot+k_i) .ge. rcut_mbd .and. rjs(n_tot+k_i) < rcut_2b-r_buffer) then
              a_2b(k2) = central_pol(i2)
            else if ( rjs(n_tot+k_i) .ge. rcut_2b-r_buffer .and. rjs(n_tot+k_i) < rcut_2b) then
              a_2b(k2) = central_pol(i2) * (1.d0 & 
                         - 3.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**2 &
                         + 2.d0 * ((rjs(n_tot+k_i)-rcut_2b+r_buffer)/(r_buffer))**3)
              !write(*,*) "a_mbd larger", a_mbd(k2)
            end if
            C6_2b(k2) = 3.d0/2.d0 * a_iso(1,2) * a_2b(k2) * (central_omega(i) * central_omega(i2)) / &
                    (central_omega(i) + central_omega(i2))
            r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (a_2b(k2)/neighbor_alpha0_2b(k2))**(1.d0/3.d0)
            r_vdw_j = r0_ii_SCS_2b(k2)
            f_damp_SCS_2b(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
            E_TS = E_TS - C6_2b(k2)/rjs_2b(k2)**6 * f_damp_SCS_2b(k2)
          end if
        end do
        E_TS = 1.d0/2.d0 * E_TS

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
            !energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,1:3,k2)
            !if ( do_derivatives ) then
            !  force_series = force_series + AT_n_f(:,:,k2)
            !  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2)
            !end if
          end do
        else
          write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
          write(*,*) "and forces will be zero."
        end if
        end do

          integrand = 0.d0
          do i2 = 1, n_freq
            energy_series = 0.d0
            do k2 = 1, n_order-1
              energy_series = energy_series - 1.d0/(k2+1) * AT_n(:,1:3,k2,i2) !* &
                              !(1.d0/(1.d0 + (omegas_mbd(i2)/0.5d0)**2))**k2 !alpha_SCS0(i,3))**2))**k2
            end do
            do c1 = 1, 3
              integrand(i2) = integrand(i2) + a_mbd(1)/(1.d0 + (omegas_mbd(i2)/o_mbd(1))**2) & 
                              * dot_product(T_LR(c1,:),energy_series(:,c1))
            end do
          end do

          integral = 0.d0
          !write(*,*) "omegas_mbd", omegas_mbd
          !write(*,*) "o_mbd", o_mbd(1)
          call integrate("trapezoidal", omegas_mbd, integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
          integral = integral/(2.d0*pi)
          energies(i) = (integral + E_TS) * 27.211386245988
          write(*,*) "MBD energy", i, energies(i)
          E_MBD = E_MBD + energies(i)          
        !else
        !  write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
        !  write(*,*) "and forces will be zero."
        !end if
        
        end if ! om loop

        
        ! Derivatives:
        
        if (do_derivatives) then
        
        allocate( da_SCS(1:3*n_sub_sites,1:3) )
        !allocate( da_iso(1:n_sub_sites,1:2) )
        !da_iso = 0.d0
        allocate( dT(1:9*n_sub_pairs) )
        allocate( dB_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
        dB_mat = 0.d0
        allocate( b_der(1:3*n_sub_sites,1:3) )
        allocate( f_damp_der(1:n_sub_pairs) )
        allocate( g_func_der(1:n_sub_pairs) )
        allocate( h_func_der(1:9*n_sub_pairs) )
        allocate( d_der(1:3*n_sub_sites,1:3) )        


        if ( do_timing ) then
        call cpu_time(time2)
        write(*,*) "gradient allocation", time2-time1
        end if

        do c3 = 1, 3  ! REMOVE NEIGHBOR-NEIGHBOR LIST DEPENDENCIES INSIDE THIS LOOP!!!!!
          if ( do_timing ) then
          call cpu_time(time1)
          end if
          dT = 0.d0
          k2 = 0
          do p = 1, n_sub_sites
            !i2 = p_to_i(p)
            !if ( in_cutoff(i2) ) then 
            k2 = k2+1
            do j3 = 2, n_sub_neigh(p)
              k2 = k2+1
              !if ( rjs_0(k2) < rcut ) then
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
              !end if
            end do
            !end if
          end do
        !write(*,*) "k2", k2
        !write(*,*) "n_sub_pairs", n_sub_pairs
        
          if ( do_timing ) then
          call cpu_time(time2)
          write(*,*) "dT timing", time2-time1
          end if        

!TEST!!!!!
          !k_a = sum(n_neigh(1:i))-n_neigh(i)
          !do r = 1, n_neigh(i)
            !k_a = k_a+1
            !if ( rjs(k_a) < 2*rcut ) then
            !a = neighbors_list(k_a)
            a = i
            if ( do_timing ) then
            call cpu_time(time1)
            end if
!TEST!!!!!
            !dB_mat = 0.d0
            f_damp_der = 0.d0
            g_func_der = 0.d0
            h_func_der = 0.d0
            b_der = 0.d0
            d_der = 0.d0
            
            if ( do_timing ) then
            call cpu_time(time2)
            write(*,*) "initialization of gradients", time2-time1
            call cpu_time(time1)
            end if
          
            k3 = 0
            do p = 1, n_sub_sites
              k3 = k3+1
              r_vdw_i = r0_ii(k3)
              s_i = neighbor_sigma(k3)
              i2 = sub_neighbors_list(k3)
              do j2 = 2, n_sub_neigh(p)
                k3 = k3+1
                !if ( rjs_0(k3) < rcut ) then
                r_vdw_j = r0_ii(k3)
                s_j = neighbor_sigma(k3)
                j = sub_neighbors_list(k3)
                j3 = modulo(sub_neighbors_list(k3)-1, n_sites) + 1
!TEST!!!!!!!!!!!
                if (a == i2 .or. a == j) then
!                if (i == i2 .or. i == j) then
!TEST!!!!!!!!!!!!
                  if ( rjs_0(k3) < rcut ) then
                    q = p_list(k3)
                  end if
                  !sigma_ij = sqrt(sigma_i(i2)**2 + sigma_i(j)**2)
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  f_damp_der(k3) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k3)**2 * &
                                     exp( -d*(rjs_H(k3)/(sR*(r_vdw_i + &
                                     r_vdw_j)) - 1.d0) ) * xyz_H(c3,k3)/rjs_H(k3)
                  g_func_der(k3) = 4.d0/sqrt(pi) * rjs_H(k3)/sigma_ij**3 * xyz_H(c3,k3) * &
                                       exp(-rjs_H(k3)**2/sigma_ij**2)
                  !if ( i == 2 .and. c3 == 1 .and. a == 1 .and. p == 40 .and. q == 41 ) then
                  !  write(*,*) "f_damp_der", f_damp_der(k3)
                  !end if
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
!TEST!!!!!!!!!!!!!!!!
                      if ( rjs_0(k3) < rcut ) then
                      if (a == i2) then
!                      if ( i == i2 ) then
!TEST!!!!!!!!!!!!!!!
                        b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                          T_SR_mult(k3) * a_SCS(3*(q-1)+c2,:)
                                                          !alpha_SCS_full(3*(j3-1)+c2,:,1)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) - (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3))) * T_SR_mult(k3)
                        !dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) - (f_damp_der(k3) * T_func(k4) * &
                        !                                  g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                  g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                  T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                  h_func_der(k4) * (1.d0 - f_damp(k3)))
                        
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
                                                          !alpha_SCS_full(3*(j3-1)+c2,:,1)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3))) * T_SR_mult(k3)
                        !dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + (f_damp_der(k3) * T_func(k4) * &
                        !                                  g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                  g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                  T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                  h_func_der(k4) * (1.d0 - f_damp(k3)))
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
                          !write(*,*) "d_der i2", d_der(3*(p-1)+c1,c2)
                        else
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - ((f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            neighbor_alpha0(k3)) * d_mult_o(k3)

                          !write(*,*) "d_der j"
                        end if
                      end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. i2 == 2 .and. j == 1 .and. c1 == 1 .and. c2 == 1) then
                      !write(*,*) "T_func", T_func(k4)
                      !write(*,*) "g_func", g_func(k3)
                      !write(*,*) "h_func", h_func(k4)
                      !write(*,*) "f_damp", f_damp(k3)
                      !write(*,*) "f_damp_der", f_damp_der(k3)
                      !write(*,*) "dT", dT(k4)
                      !write(*,*) "g_func_der", g_func_der(k3)
                      !write(*,*) "h_func_der", h_func_der(k4)
                      !end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. p == 1 .and. q == 3 ) then
                      !  write(*,*) "f_damp_der", f_damp_der(k3)
                      !  write(*,*) "g_func_der", g_func_der(k3)
                      !  write(*,*) "h_func_der", h_func_der(k4)
                      !  write(*,*) "dT", dT(k4)
                      !  write(*,*) "T_func, g_func, h_func, f_damp", T_func(k4), g_func(k3), h_func(k4), f_damp(k3)
                      !end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. p == 3 .and. q == 1 ) then
                      !  write(*,*) "f_damp_der", f_damp_der(k3)
                      !  write(*,*) "g_func_der", g_func_der(k3)
                      !  write(*,*) "h_func_der", h_func_der(k4)
                      !  write(*,*) "dT", dT(k4)
                      !  write(*,*) "T_func, g_func, h_func, f_damp", T_func(k4), g_func(k3), h_func(k4), f_damp(k3)
                      !end if
                      !write(*,*) "a, i2, j", a, i2, j
                      !write(*,*) "b_der", (f_damp_der(k3) * T_func(k4) * &
                      !                                    g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                      !                                    g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                      !                                    T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                      !                                    h_func_der(k4) * (1.d0 - f_damp(k3)))
                    end do
                  end do
                end if
                !end if
              end do
            end do
            
            !if ( a == 1 .and. c3 == 1 .and. i == 1) then
            !  open(unit=69, file="dB_mat.dat", status="new")
            !  write(*,*) "dB_mat"
            !  write(*,*) a, c3, i
            !  write(*,*) "before"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) b_der(p,:)
            !  end do
            !  close(69)
            !  write(*,*) "after"
            !  b_der = matmul(dB_mat,alpha_SCS_full(:,:,1))
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) b_der(p,:)
            !  end do
            !end if

            if ( do_timing ) then
            call cpu_time(time2)
            write(*,*) "gradients without hirshfeld timing", time2-time1
            end if            

            !if ( i == 1 .and. c3 == 1 .and. om == 1) then
            !write(*,*) "dB_mat"
            !p = 40
            !write(*,*) p, p_to_i(p)
            !do p = 1, 3*n_sub_sites
            !  write(*,*) dB_mat(p,:)
            !end do
            !end if

            if (do_hirshfeld_gradients) then
           
              if ( do_timing ) then
              call cpu_time(time1)
              end if

              !allocate( coeff_der(1:n_sub_sites,1:n_sub_sites,1:3,1:3) )
              !allocate( coeff_fdamp(1:n_sub_sites,1:n_sub_sites,1:3,1:3) )
              allocate( coeff_der(1:9*n_sub_pairs) )
              allocate( coeff_fdamp(1:9*n_sub_pairs) )
              
              !in_force_cutoff = .false.
              !k3 = sum(n_neigh(1:a))-n_neigh(a)
              !do j2 = 1, n_neigh(a)
              !  k3 = k3+1
              !  j = neighbors_list(k3)
              !  if ( rjs(k3) < 2 * rcut ) then
              !    in_force_cutoff(j) = .true.
              !  end if
              !end do
              !xyz_a = xyz(:,k_a)/Bohr
            
              coeff_fdamp = 0.d0
              coeff_der = 0.d0
              
              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "allocation and force cutoff array timing", time2-time1
              call cpu_time(time1)
              end if
              
              k3 = 0
              do p = 1, n_sub_sites
                k3 = k3+1
                r_vdw_i = r0_ii(k3)
                s_i = neighbor_sigma(k3)
                i2 = sub_neighbors_list(k3)
                do j2 = 2, n_sub_neigh(p)
                  k3 = k3+1
                  !if ( rjs_0(k3) < rcut ) then
                  r_vdw_j = r0_ii(k3)
                  s_j = neighbor_sigma(k3)
                  j = sub_neighbors_list(k3)
                  !sigma_ij = sqrt(sigma_i(i2)**2 + sigma_i(j)**2)
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
                      !coeff_fdamp(p,q,c1,c2) = (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                      !                         (-T_func(k4) * g_func(k3) + h_func(k4))
                      coeff_fdamp(k4) = (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                               (-T_func(k4) * g_func(k3) + h_func(k4))
                      dg = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * (rjs_H(k3)/sigma_ij)**2 * &
                           (-rjs_H(k3)/(3.d0 * sigma_ij**3))
                      dh = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * xyz_H(c1,k3)*xyz_H(c2,k3) / &
                           (sigma_ij**5 * rjs_H(k3)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k3)/sigma_ij)**2)
                      !coeff_der(p,q,c1,c2) = (1.d0-f_damp(k3)) * (-T_func(k4)*dg+dh)
                      coeff_der(k4) = (1.d0-f_damp(k3)) * (-T_func(k4)*dg+dh)
                      !if ( i == 2 .and. p == 40 .and. q == 41 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                      !       .and. c2 == 2 ) then
                      !  write(*,*) "coeff_fdamp first", (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2
                      !end if
                    end do
                  end do
                  !end if
                end do
              end do

              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "calculation of coeff_fdamp and coeff_der timing", time2-time1
              call cpu_time(time1)
              end if

              k3 = 0
              do p = 1, n_sub_sites
                !write(*,*) k3, "/", n_sub_pairs
                i2 = sub_neighbors_list(k3+1)
                !write(*,*) "I2", i2
                !n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
                !xyz_i = xyz_H(:,k3+1)
                !if ( in_cutoff(a) .and. in_force_cutoff(i2) ) then
                !if ( sqrt(sum((xyz_a-xyz_i)**2)) < rcut/Bohr ) then ! This should be SOAP cut-off (at the moment rcut_scs = rcut_soap)!
                !if ( dabs(hirshfeld_v_sub_der(c3,p)) > 1.d-10 ) then
                  hv_p_der = hirshfeld_v_sub_der(c3,p)*Bohr
                  r_vdw_i = r0_ii(k3+1)
                  s_i = neighbor_sigma(k3+1)
                  do c1 = 1, 3
                    !dB_mat(3*(p-1)+c1,3*(p-1)+c1) = - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v(i2)) * &
                    !                         hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr
                    b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_sub_neigh(k3+1)) * &
                                          hirshfeld_v_sub_der(c3,p)*Bohr * a_SCS(3*(p-1)+c1,:) !alpha_SCS_full(3*(i2-1)+c1,:,1)
                                          !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                    dB_mat(3*(p-1)+c1,3*(p-1)+c1) = dB_mat(3*(p-1)+c1,3*(p-1)+c1) - 1.d0/(neighbor_alpha0(k3+1) &
                                          * hirshfeld_sub_neigh(k3+1)) * &
                                          hirshfeld_v_sub_der(c3,p)*Bohr
                  end do
                  !j3 = findloc(neighbors_list(n_tot2+1:n_tot2+n_neigh(i2)),a,1)
                  do j2 = 2, n_sub_neigh(p)
                    j = sub_neighbors_list(k3+j2)
                    r_vdw_j = r0_ii(k3+j2)
                    s_j = neighbor_sigma(k3+j2)
                    if ( rjs_0(k3+j2) < rcut ) then
                      q = p_list(k3+j2)
                      hv_q_der = hirshfeld_v_sub_der(c3,q)*Bohr
                    !write(*,*) "k3+j2", k3+j2, "/", n_sub_pairs
                    !write(*,*) "j, q", j, q
                    else
                      hv_q_der = 0.d0
                    end if
                    !if ( i2 == j ) then
                    !  do c1 = 1, 3
                        !dB_mat(3*(p-1)+c1,3*(p-1)+c1) = - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v(i2)) * &
                        !                         hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr
                    !    b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v_neigh(k3+1)) * &
                    !                          hirshfeld_v_sub_der(c3,p)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                                              !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                    !  end do
                    !end if
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
                        if ( p == 2 .and. q == 1 .and. c1 == 1 .and. c2 == 1 .and. om == 2 .and. i == 1 .and. c3 == 1 ) then
                          write(*,*) "p, q", p, q
                          write(*,*) "coeff_der", coeff_der(k4)
                          write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                          write(*,*) "T_SR_mult", T_SR_mult(k3+j2)
                          write(*,*) "s_i, s_j", s_i, s_j
                          write(*,*) "r_vdw_i, r_vdw_j", r_vdw_i, r_vdw_j
                          write(*,*) "hv_p_der, hv_q_der", hv_p_der, hv_q_der
                          write(*,*) "hirshfeld_i, hirshfeld_j", hirshfeld_sub_neigh(k3+1), hirshfeld_sub_neigh(k3+j2)
                        end if
                        if ( p == 1 .and. q == 2 .and. c1 == 1 .and. c2 == 1 .and. om == 2 .and. i == 1 .and. c3 == 1 ) then
                          write(*,*) "p, q", p, q
                          write(*,*) "coeff_der", coeff_der(k4)
                          write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                          write(*,*) "T_SR_mult", T_SR_mult(k3+j2)
                          write(*,*) "s_i, s_j", s_i, s_j
                          write(*,*) "r_vdw_i, r_vdw_j", r_vdw_i, r_vdw_j
                          write(*,*) "hv_p_der, hv_q_der", hv_p_der, hv_q_der
                          write(*,*) "hirshfeld_i, hirshfeld_j", hirshfeld_sub_neigh(k3+1), hirshfeld_sub_neigh(k3+j2)
                        end if
                        d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * &
                          neighbor_alpha0(k3+j2) * d_mult_i(k3+j2)
!                        b_der(3*(q-1)+c1,:) = b_der(3*(q-1)+c1,:) + &
!                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
!                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
!                          hirshfeld_v_sub_der(c3,p)*Bohr) * &
!                          a_SCS(3*(p-1)+c2,:) * T_SR_mult(k3+j2)
!                          !alpha_SCS_full(3*(i2-1)+c2,:,1)
!                        dB_mat(3*(q-1)+c1,3*(p-1)+c2) = dB_mat(3*(q-1)+c1,3*(p-1)+c2) + &
!                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
!                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
!                          hirshfeld_v_sub_der(c3,p)*Bohr) * T_SR_mult(k3+j2)
!
!                       d_der(3*(q-1)+c1,c2) = d_der(3*(q-1)+c1,c2) - &
!                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
!                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
!                          hirshfeld_v_sub_der(c3,p)*Bohr) * &
!                          neighbor_alpha0(k3+j2) * d_mult_i(k3+j2)                        

                      !  if ( i == 1 .and. c3 == 1 .and. l == 1 .and. p == 50 .and. q == 48 .and. c1 == 3 .and. c2 == 3 ) then
                      !  write(*,*) "coeff_der", coeff_der(k4)
                      !  write(*,*) "s_i", s_i
                      !  write(*,*) "r_vdw_i", r_vdw_i
                      !  write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh(k3+1)
                      !  write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                      !  write(*,*) "hirshfeld_v_sub_der", hirshfeld_v_sub_der(c3,p)*Bohr
                      !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                      !  end if

                      !  if ( i == 1 .and. c3 == 1 .and. l == 1 .and. p == 48 .and. q == 50 .and. c1 == 3 .and. c2 == 3 ) then
                      !  write(*,*) "coeff_der", coeff_der(k4)
                      !  write(*,*) "s_i", s_i
                      !  write(*,*) "r_vdw_i", r_vdw_i
                      !  write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh(k3+1)
                      !  write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                      !  write(*,*) "hirshfeld_v_sub_der", hirshfeld_v_sub_der(c3,p)*Bohr
                      !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                      !  end if


                          !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr) * &
                          !alpha_SCS_full(3*(i2-1)+c2,:,1)
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

                        !b_der(3*(q-1)+c1,:) = b_der(3*(q-1)+c1,:) + &
                        !  ((coeff_der(q,p,c1,c2) * sigma_i(i2)**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(q,p,c1,c2) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)) * &
                        !  alpha_SCS_full(3*(i2-1)+c2,:)
                          
                        !dB_mat(3*(q-1)+c1,3*(p-1)+c2) = dB_mat(3*(q-1)+c1,3*(p-1)+c2) + &
                        !  ((coeff_der(k4) * s_i**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(k4) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr) 
                        !if ( i == 2 .and. p == 40 .and. q == 41 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                        !     .and. c2 == 2) then
                        !  write(*,*) "p, q, i2, j", p, q, i2, j
                        !  write(*,*) "coeff_fdamp", coeff_fdamp(p,q,c1,c2) !/(-T_func(k4) * g_func(k3) + h_func(k4)) !* &
                        !                            !r_vdw_i/hirshfeld_v(i2) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                        !  write(*,*) "hirshfeld_v_cart_der term", r_vdw_i/hirshfeld_v(i2) * &
                        !                  hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !end if
                        !if ( i == 2 .and. p == 41 .and. q == 40 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                        !     .and. c2 == 2 ) then
                        !  write(*,*) "p, q", p, q
                        !  write(*,*) "coeff_fdamp", coeff_fdamp(q,p,c1,c2) !/(-T_func(k4) * g_func(k3) + h_func(k4)) * &
                        !                            !r_vdw_i/hirshfeld_v(i2) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !  write(*,*) "dB_mat", dB_mat(3*(q-1)+c1,3*(p-1)+c2)
                        !  write(*,*) "hirshfeld_v_cart_der term", r_vdw_i/hirshfeld_v(i2) * &
                        !                 hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !end if
                      end do
                    end do
                  end do
                !end if
                k3 = k3 + n_sub_neigh(p)
              end do
              
              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "hirshfeld gradient timing", time2-time1
              end if            

              deallocate( coeff_der, coeff_fdamp )
              
            end if
            
            if ( do_timing ) then
            call cpu_time(time1)
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


            if ( i == 1 .and. c3 == 1 .and. om == 1 ) then

            write(*,*) "neighbor_alpha0", neighbor_alpha0

            write(*,*) "a_SCS"
            do p = 1, 3*n_sub_sites
              write(*,*) a_SCS(p,:)
            end do

            write(*,*) "d_der"

            do p = 1, 3*n_sub_sites
              write(*,*) d_der(p,:)
            end do

            end if

            b_der = -b_der+d_der

            call dsytrs( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, b_der, 3*n_sub_sites, work_arr, &
                         12*n_sub_sites, info )
            !call dgetrs( 'N', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, b_der, 3*n_sub_sites, info )

            da_SCS = b_der

            ! MBD forces
            if ( i == 1 .and. c3 == 1 .and. om == 1 ) then
              write(*,*) "da_SCS"
              do p = 1, 3*n_sub_sites
                write(*,*) da_SCS(p,:)
              end do
            end if
            

            if ( i == 1 .and. c3 == 1) then
              write(*,*) "om", om
            end if
            
            do p = 1, n_sub_sites
              do c1 = 1, 3
                da_iso(p,c3,om) = da_iso(p,c3,om) + da_SCS(3*(p-1)+c1,c1)
              end do
            end do
            da_iso(:,c3,om) = da_iso(:,c3,om)/3.d0
        

            if ( i == 1 .and. c3 == 1 ) then
              write(*,*) "da_SCS"
              k2 = 0
              do p = 1, n_sub_sites
                write(*,*) sub_neighbors_list(k2+1), 1.d0/3.d0 * ( da_SCS(3*(p-1)+1,1) + da_SCS(3*(p-1)+2,2)+ &
                                                                   da_SCS(3*(p-1)+3,3) )
                k2 = k2+n_sub_neigh(p) 
              end do
            end if
            
            if ( om == 2 ) then
            
            write(*,*) "do_mbd"
            
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
              if ( rjs_0_mbd(k2) < (rcut-r_buffer)/Bohr ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                da_mbd(k2) = da_iso(r,c3,2)
                do_pref = -omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                             ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                do_mbd(k2) = do_pref
                if ( i == 1 .and. c3 == 1 ) then
                  !write(*,*) da_iso(r,c3,1)
                  write(*,*) do_mbd(k2)
                end if
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
                             central_pol(i2) * &
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
                             central_omega(i2) * &
                             ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                               - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                             * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
              else if ( rjs_0_mbd(k2) .ge. rcut/Bohr .and. rjs_0_mbd(k2) < (rcut_mbd-r_buffer)/Bohr ) then
                da_mbd(k2) = 0.d0
                do_mbd(k2) = 0.d0
              else if ( rjs_0_mbd(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut_mbd/Bohr ) then
                da_mbd(k2) = central_pol(i2) &
                     * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                     + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                     * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                do_mbd(k2) = central_omega(i2) &
                     * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                     + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                     * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
              end if
              da_i = da_mbd(k2)
              a_mbd_i = a_mbd(k2)
              do j3 = 2, n_mbd_neigh(p)
                k2 = k2+1
                j = mbd_neighbors_list(k2)
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
                               central_pol(j) * &
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
                               central_omega(j) * &
                                ( + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer) &
                                 - 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut+r_buffer)/r_buffer)**2) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                else if ( rjs_0_mbd(k2) .ge. rcut/Bohr .and. rjs_0_mbd(k2) < (rcut_mbd-r_buffer)/Bohr ) then
                  da_mbd(k2) = 0.d0
                  do_mbd(k2) = 0.d0
                else if ( rjs_0_mbd(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_0_mbd(k2) < rcut_mbd/Bohr ) then
                  da_mbd(k2) = central_pol(j) &
                       * (- 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer) &
                       + 6.d0 * ((rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buffer)/r_buffer)**2) &
                       * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buffer/Bohr))
                  do_mbd(k2) = central_omega(j) &
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
            s = neighbor_species(n_tot+1)
            r_vdw_i = r0_ref(s) / Bohr * (a_iso(1,2)/(alpha0_ref(s)/Bohr**3))**(1.d0/3.d0)
            dr_vdw_i = r_vdw_i / (3.d0 * a_iso(1,2)) * da_iso(1,c3,2)
            do p = 1, n_2b_sites
              k2 = k2+1
              i2 = sub_2b_list(k2)
              if ( rjs_2b(k2) .ge. (rcut_mbd-r_buffer)/Bohr .and. rjs_2b(k2) < rcut_mbd/Bohr ) then
                da_2b(k2) = a_2b(k2) * &
                              ( + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer)) &
                                - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_mbd+r_buffer)/(r_buffer))**2) &
                                * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
              else if ( rjs_2b(k2) .ge. rcut_mbd/Bohr .and. rjs_2b(k2) < (rcut_2b-r_buffer)/Bohr ) then
                da_2b(k2) = 0.d0
              else if ( rjs_2b(k2) .ge. (rcut_2b-r_buffer)/Bohr .and. rjs_2b(k2) < rcut_2b/Bohr ) then
                da_2b(k2) = a_2b(k2) * &
                       ( - 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer)) &
                         + 6.d0 * ((rjs_2b(k2)*Bohr-rcut_2b+r_buffer)/(r_buffer))**2) &
                       * ( -xyz_2b(c3,k2)/rjs_2b(k2)/(r_buffer/Bohr))
              end if
              r_vdw_j = r0_ii_SCS_2b(k2)
              dr_vdw_j = r_vdw_j / (3.d0 * a_2b(k2)) * da_2b(k2)
              !f_damp_SCS_2b(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
              dC6_2b = 3.d0/2.d0*central_omega(i)*central_omega(i2) &
                          / (central_omega(i)+central_omega(i2)) &
                          * (da_iso(1,c3,2)*a_2b(k2) + a_iso(1,2)*da_2b(k2))
              f_damp_der_2b = d/0.97d0 * f_damp_SCS_2b(k2)**2 * &
                                      exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                      * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (-xyz_2b(c3,k2)/rjs_2b(k2)) &
                                      - rjs_2b(k2)/0.97d0 * 1.d0/(r_vdw_i+r_vdw_j)**2 &
                                      * (dr_vdw_i + dr_vdw_j))
              forces_TS = forces_TS + ( dC6_2b * f_damp_SCS_2b(k2) / rjs_2b(k2)**6 &
                                      + 6.d0/rjs_2b(k2)**8 * xyz_2b(c3,k2) * C6_2b(k2) * f_damp_SCS_2b(k2) &
                                      + C6_2b(k2)/rjs_2b(k2)**6 * f_damp_der_2b )
            end do
            forces_TS = 1.d0/2.d0 * forces_TS


            !if ( i == 1 .and. c3 == 1 ) then
            !  write(*,*) "dT_LR"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) dT_LR(p,:)
            !  end do
            !end if

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
                a_mbd(k3+1) * (2.d0 * omegas_mbd(j) * o_mbd(k3+1)) / &
                ( o_mbd(k3+1)**2 + omegas_mbd(j)**2 ) * &
                T_LR(3*(p-1)+1:3*(p-1)+3,:)
              k3 = k3+n_mbd_neigh(p)
            end do
            end do


            if ( n_order > 1 ) then


              integrand = 0.d0
              total_integrand = 0.d0
              do j = 1, n_freq

                force_series = 0.d0
                total_energy_series = 0.d0
                do k2 = 1, n_order-1
                  force_series = force_series + AT_n_f(:,:,k2,j) !* &
                                 !(1.d0/(1.d0 + (omegas_mbd(j)/0.5d0)**2))**k2 !/alpha_SCS0(i,3))**2))**k2
                  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2,j) !* &
                                 !(1.d0/(1.d0 + (omegas_mbd(j)/0.5d0)**2))**k2 !/alpha_SCS0(i,3))**2))**k2
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
              forces0(c3,i) = forces0(c3,i) + (1.d0/(2.d0*pi) * integral + forces_TS) * 51.42208619083232
              !forces0(c3,i) = 1.d0/(2.d0*pi) * integral * 51.42208619083232
              !forces0(c3,i) = forces0(c3,i) * 51.42208619083232

              write(*,*) "MBD force", i, c3, forces0(c3,i)
              !write(*,*) i, c3, forces0(c3,i)
              integral = 0.d0
              if (c3 == 1 ) then
              call integrate("trapezoidal", omegas_mbd, total_integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
              write(*,*) "MBD total energy of sphere", i, (integral / (2.d0*pi) + E_TS) * 27.211386245988
              !write(*,*) "Total TS energy of sphere", i, integral/(2.d0*pi) * 27.211386245988
              end if

            end if
            
          end if ! om loop(?)

        end do ! c3 loop
        
        deallocate( da_SCS, dT, dB_mat, b_der, g_func_der, h_func_der, d_der, f_damp_der )

        end if ! do_derivatives
        
        if ( om == 2 ) then
        deallocate( T_LR, r0_ii_SCS, f_damp_SCS, AT, AT_n, energy_series, omegas_mbd, integrand, n_mbd_neigh, &
                    mbd_neighbors_list, p_mbd, r0_ii_mbd, neighbor_alpha0_mbd, xyz_mbd, rjs_mbd, T_mbd, a_mbd, &
                    rjs_0_mbd, xyz_0_mbd, o_mbd, sub_2b_list, xyz_2b, rjs_2b, r0_ii_2b, neighbor_alpha0_2b, &
                    a_2b, r0_ii_SCS_2b, f_damp_SCS_2b, C6_2b )
        end if
                    
        if ( do_derivatives .and. om == 2 ) then
          deallocate( da_mbd, AT_n_f, dT_mbd, f_damp_der_mbd, f_damp_der_SCS, dT_LR, G_mat, force_series, VL, &
                      total_energy_series, total_integrand, alpha_grad, da_2b, do_mbd )
        end if
        
        !end if

        !deallocate( val_xv, val_bv, myidx )

!call cpu_time(time2)
!write(*,*) "polynomial timing", time2-time1


        end do ! om

        ! store omega_p for MBD calculation
        !alpha_SCS0(i,3) = 2.d0*omega_ref/sqrt(alpha_SCS0(i,1)/alpha_SCS0(i,2)-1.d0)
        !write(*,*) "omega_p", alpha_SCS0(i,3)

        !write(*,*) "omega_p", alpha_SCS0(i,3)

        deallocate( n_sub_neigh, sub_neighbors_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, omegas, &
                    T_func, b_i, d_vec, g_func, h_func, f_damp, a_SCS, ipiv, ia, ja, val, p_list, work_arr, rjs_0, &
                    a_iso, o_p, T_SR, T_SR_mult, d_arr_i, d_arr_o, d_mult_i, d_mult_o, dT_SR_mult, d_dmult_i, &
                    d_dmult_o, hirshfeld_sub_neigh )

        deallocate( B_mat)

        if ( do_derivatives ) then
          deallocate( da_iso )
        end if
        if ( do_derivatives .and. do_hirshfeld_gradients ) then
          deallocate( hirshfeld_v_sub_der )
        end if
           
      end do

      write(*,*) "E_MBD", E_MBD

      !write(*,*) "alpha_SCS_full"
      !open(unit=89, file="alpha_SCS_full.dat", status="new")
      !do p = 1, 3*n_sites
      !  write(*,*) alpha_SCS_full(p,:,1)
      !end do
      !close(89)

      !write(*,*) "A_i"
      !do p = 1, 3*n_sites
      !  write(*,*) A_i(p,:)
      !end do
      
      !close(69)
      !close(79)
      !write(*,*) "Writing polarizabilities"
      !open(unit=69, file="alpha_SCS_polynomial.dat", status="new")
      !do i = 1, n_sites
      !  write(69,*) alpha_SCS0(i,1)
      !end do
      !close(69)

      ! Next we do the force calculation separately after having calculated A_i and alpha_SCS0 for all atoms first.
      ! We should have alpha_SCS0 and A_i stored now.
      ! We only need dB/dr for each atom within the MBD cut-off of atom i.
      
      !write(*,*) "alpha_SCS_full"
      !alpha_SCS_full = 0.d0
      !do i = 1, n_sites
      !  do c1 = 1, 3
      !    alpha_SCS_full(3*(i-1)+c1,c1,1) = alpha_SCS0(i,1)
      !  end do
      !end do

      !if ( polynomial_expansion ) then
      if ( .false. ) then

      !write(*,*) "Calculating gradients"

polyfit = (/ 3.237385145550585d+02, -4.241125470183307d+04, 3.008572712845031d+06, -1.309430416378132d+08, &
3.756106046665028d+09, -7.433108326602138d+10, 1.045457946646248d+12, &
-1.064278184563909d+13, 7.908875986032283d+13, -4.286093180295281d+14, &
1.673744363742281d+15, -4.582925299269683d+15, 8.343821283398570d+15, &
-9.066835011532402d+15, 4.447833222479864d+15 /)

         !  polyfit = (/ 1.488424328064546e+02, -1.013335812271306e+04, 4.172679550725457e+05, &
         !                        -1.157445216098117e+07, 2.278478136982700e+08, &
         ! -3.263921412009479e+09, 3.428207402483869e+10, -2.620296686535526e+11, 1.421359676944887e+12, &
         ! -5.193951907119399e+12, 1.148018668128465e+13, -1.160719499531639e+13 /)

      !polyfit = (/ 2.619652273544946e+02, -2.541165321210297e+04, 1.238702036672193e+06, &
      ! -3.484275967677591e+07, 6.121187495909173e+08, &
      ! -7.039412728611286e+09, 5.432799199845714e+10, -2.831213305568388e+11, 9.822909546821646e+11, &
      ! -2.172784949806289e+12, 2.770968854860906e+12, -1.549957105267603e+12 /)


      k = 0
      do i = 1, n_sites
      
        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !write(*,*) i, "/", n_sites
        !p_to_i = 0
        !i_to_p = 0
        !in_cutoff = .false.
        !k = k+1
        s = neighbor_species(n_tot+1)
        omega_ref = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
        !p = 1
        !p_to_i(p) = i
        !i_to_p(i) = p
        !in_cutoff(i) = .true.
        !do j2 = 2, n_neigh(i)
        !  k = k+1
        !  j = neighbors_list(k)
        !  if (rjs(k) < 2*rcut) then
        !    in_cutoff(j) = .true.
        !    p = p+1
        !    p_to_i(p) = j
        !    i_to_p(j) = p
        !  end if
        !end do
        !write(*,*) "in cutoff"
        !write(*,*) in_cutoff
        !write(*,*) "p_to_i"
        !write(*,*) p_to_i
        !write(*,*) "i_to_p"
        !write(*,*) i_to_p
        !n_sub_sites = p
        !write(*,*) "n_sub_sites", n_sub_sites
        
        call cpu_time(time1)

        rcut_mbd = 9.d0
        n_mbd_sites = 0
        n_mbd_pairs = 0
        
        k_i = 0
        do i3 = 1, n_neigh(i)
          k = k+1
          k_i = k_i + 1
          !j = p_to_i(p)
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
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut_mbd/Bohr ) then
                  n_mbd_pairs = n_mbd_pairs + 1
                end if
                end if
              end if
            end do
          end if
        end do

        n_order = 6
        n_freq = 12

        allocate( T_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
        allocate( r0_ii_SCS(1:n_mbd_pairs) )
        allocate( f_damp_SCS(1:n_mbd_pairs) )
        allocate( AT(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_freq) )
        allocate( AT_n(1:3*n_mbd_sites,1:3,n_order-1,1:n_freq) )
        allocate( energy_series(1:3*n_mbd_sites,1:3) )
        allocate( omegas(1:n_freq) )
        allocate( integrand(1:n_freq) )
        allocate( n_mbd_neigh(1:n_mbd_sites) )
        allocate( mbd_neighbors_list(1:n_mbd_pairs) )
        allocate( p_mbd(1:n_mbd_pairs) )
        allocate( r0_ii_mbd(1:n_mbd_pairs) )
        allocate( neighbor_alpha0_mbd(1:n_mbd_pairs) )
        allocate( xyz_mbd(1:3,n_mbd_pairs) )
        allocate( rjs_mbd(n_mbd_pairs) )
        allocate( T_mbd(1:9*n_mbd_pairs) )
        if ( do_derivatives ) then
          allocate( AT_n_f(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_order-1,1:n_freq) )
          allocate( dT_mbd(1:9*n_mbd_pairs) )
          allocate( f_damp_der_mbd(1:n_mbd_pairs) )
          allocate( f_damp_der_SCS(1:n_mbd_pairs) )
          allocate( dT_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( G_mat(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_freq) )
          allocate( force_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( VL(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( total_energy_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          allocate( total_integrand(1:n_freq) )
          allocate( alpha_grad(1:3,1:n_mbd_sites) )
          alpha_grad = 0.d0
        end if
        
        omegas = 0.d0
        omega = 0.d0
        do i2 = 1, n_freq
          omegas(i2) = omega
          omega = omega + 0.4d0
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
          if ( rjs(n_tot+k_i) < rcut_mbd ) then
            p = p+1
            k2 = k2+1
            s = neighbor_species(n_tot+k_i)
            mbd_neighbors_list(k2) = neighbors_list(n_tot+k_i)
            n_mbd_neigh(p) = n_mbd_neigh(p) + 1
            p_mbd(k2) = p
            r0_ii_mbd(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            r0_ii_SCS(k2) = r0_ii_mbd(k2) * (alpha_SCS0(i2,1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
            xyz_mbd(:,k2) = xyz(:,n_tot+k_i)/Bohr
            xyz_i = xyz_mbd(:,k2)
            rjs_mbd(k2) = rjs(n_tot+k_i)/Bohr
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
                  k2 = k2+1    
                  s = neighbor_species(n_tot+k_j)
                  mbd_neighbors_list(k2) = j
                  p_mbd(k2) = q                       
                  r0_ii_mbd(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                  neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                  r0_ii_SCS(k2) = r0_ii_mbd(k2) * (alpha_SCS0(j,1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
                  xyz_mbd(:,k2) = xyz_j-xyz_i
                  rjs_mbd(k2) = sqrt(sum((xyz_j-xyz_i)**2))
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

        k3 = 0
        do p = 1, n_mbd_sites
          i2 = mbd_neighbors_list(k3+1)
          AT(3*(p-1)+1:3*(p-1)+3,:,1) = alpha_SCS0(i2,1) * T_LR(3*(p-1)+1:3*(p-1)+3,:)
          k3 = k3 + n_mbd_neigh(p)
        end do

        if ( n_order > 1 ) then
          do k2 = 1, n_order-1
            ! Precalculate the full AT_n for forces:
            if ( k2 == 1 ) then
              AT_n(:,:,k2,1) = AT(:,1:3,1)
              if ( do_derivatives ) then
                AT_n_f(:,:,k2,1) = AT(:,:,1)
              end if
            else
              call dgemm('n', 'n', 3*n_mbd_sites, 3, 3*n_mbd_sites, 1.d0, AT(:,:,1), &
                         3*n_mbd_sites, AT_n(:,:,k2-1,1), 3*n_mbd_sites, 0.d0, AT_n(:,:,k2,1), &
                         3*n_mbd_sites)
              if ( do_derivatives ) then
                call dgemm('n', 'n', 3*n_mbd_sites, 3*n_mbd_sites, 3*n_mbd_sites, 1.d0, AT(:,:,1), &
                           3*n_mbd_sites, AT_n_f(:,:,k2-1,1), 3*n_mbd_sites, 0.d0, AT_n_f(:,:,k2,1), &
                           3*n_mbd_sites)
              end if
            end if
            !energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,1:3,k2)
            !if ( do_derivatives ) then
            !  force_series = force_series + AT_n_f(:,:,k2)
            !  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2)
            !end if
          end do
          integrand = 0.d0
          do om = 1, n_freq
            energy_series = 0.d0
            do k2 = 1, n_order-1
              energy_series = energy_series - 1.d0/(k2+1) * AT_n(:,1:3,k2,1) * &
                              (1.d0/(1.d0 + (omegas(om)/0.5d0)**2))**k2 !alpha_SCS0(i,3))**2))**k2
            end do
            do c1 = 1, 3
              integrand(om) = integrand(om) + alpha_SCS0(i,1)/(1.d0 + (omegas(om)/0.5d0)**2) & !/alpha_SCS0(i,3))**2) &
                              * dot_product(T_LR(c1,:),energy_series(:,c1))
            end do
          end do

          integral = 0.d0
          call integrate("trapezoidal", omegas, integrand, omegas(1), omegas(n_freq), integral)
          integral = integral/(2.d0*pi)
          energies(i) = integral * 27.211386245988
          !write(*,*) "MBD energy", i, energies(i)

          !if ( do_derivatives ) then
          !  do p = 1, n_sub_sites
          !    i2 = p_to_i(p)
          !    do c1 = 1, 3
          !      total_integrand(i,om) = total_integrand(i,om) + alpha_SCS0(i,1)/(1.d0 + &
          !           (omegas(om)/alpha_SCS0(i2,3))**2) * dot_product(T_LR(3*(p-1)+c1,:), &
          !           total_energy_series(:,3*(p-1)+c1))
          !      end do
          !  end do
          !end if

        else
          write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
          write(*,*) "and forces will be zero."
        end if

        call cpu_time(time2)
        !write(*,*) "MBD energy timing", time2-time1
        
    if ( do_derivatives ) then
        
      call cpu_time(time1)

      do l = 1, n_neigh(i)
      !do l = 1, 1

      if ( rjs(n_tot+l) < rcut ) then

        xyz_l = xyz(1:3,n_tot+l)/Bohr

        n_sub_sites = 0
        n_sub_pairs = 0
        !n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !allocate( n_sub_neigh(1:n_sub_sites) )
        !n_sub_neigh = 0
        k_i = 0
        do i3 = 1, n_neigh(i)
          k = k+1
          k_i = k_i + 1
          xyz_i = xyz(1:3,n_tot+k_i)/Bohr
          !j = p_to_i(p)
          if (sqrt(sum((xyz_i-xyz_l)**2)) < rcut/Bohr ) then
            n_sub_sites = n_sub_sites + 1
            n_sub_pairs = n_sub_pairs + 1
            !n_sub_neigh(p) = n_sub_neigh(p) + 1
            !n_tot2 = sum(n_neigh(1:j))-n_neigh(j)
            !write(*,*) "n_tot2", n_tot2
            !xyz_i = xyz(:,n_tot+k_i)/Bohr
            k_j = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j + 1
              xyz_j = xyz(1:3,n_tot+k_j)/Bohr
              if ( sqrt(sum((xyz_j-xyz_l)**2)) < rcut/Bohr ) then
                if (i3 .ne. j3) then
                !xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                  !if (i == 1 .and. p == 7) then
                  !write(*,*) "i, j", j, neighbors_list(n_tot2+j3)
                  !end if
                  !n_sub_neigh(p) = n_sub_neigh(p) + 1
                  n_sub_pairs = n_sub_pairs + 1
                end if
                end if
              end if
            end do
          end if
        end do

        !write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_sub_pairs", n_sub_pairs
     
        !write(*,*) "p_to_i", p_to_i
        !write(*,*) "i_to_p", i_to_p
        !write(*,*) "n_sub_pairs", n_sub_pairs
        !write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_sub_neigh", n_sub_neigh
        
        allocate( sub_neighbors_list(1:n_sub_pairs) )
        allocate( p_list(1:n_sub_pairs) )
        allocate( n_sub_neigh(1:n_sub_sites) )
        allocate( xyz_H(1:3,1:n_sub_pairs) )
        allocate( rjs_H(1:n_sub_pairs) )
        allocate( r0_ii(1:n_sub_pairs) )
        allocate( neighbor_alpha0(1:n_sub_pairs) )
        allocate( neighbor_sigma(1:n_sub_pairs) )
        !allocate( omegas(1:n_sub_pairs) )
        allocate( T_func(1:9*n_sub_pairs) )
        allocate( B_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
        allocate( f_damp(1:n_sub_pairs) )
        allocate( g_func(1:n_sub_pairs) )
        allocate( h_func(1:9*n_sub_pairs) )
        allocate( hirshfeld_v_sub_der(1:3,1:n_sub_sites) )
        !allocate( a_vec(1:3*n_sub_sites,1:3) )
        !allocate( a_SCS(1:3*n_sub_sites,1:3) )
        allocate( ipiv(1:3*n_sub_sites) )
        allocate( work_arr(1:12*n_sub_sites) )        

        allocate( ia(1:9*n_sub_pairs) )
        allocate( ja(1:9*n_sub_pairs) )
        allocate( val(1:9*n_sub_pairs) )
        allocate( d_vec(1:3*n_sub_sites,1:3) )

        allocate( rjs_0(1:n_sub_pairs) )
        !allocate( b_i(1:3*n_sub_sites,1:3) )
        
        !write(*,*) "allocation successful"

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
        nnz = 0
        ia = 0
        ja = 0
        val = 0.d0

        ipiv = 0
        work_arr = 0.d0

        !a_SCS = 0.d0
        !b_i = 0.d0
        !do p = 1, n_sub_sites
        !  do c1 = 1, 3
        !    !B_mat(j2,j2) = 1.d0
        !    !B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/alpha_i(p_to_i(p))
        !    nnz = nnz+1
        !    ia(nnz) = 3*(p-1)+c1
        !    ja(nnz) = 3*(p-1)+c1
        !    !val(nnz) = 1.d0
        !    val(nnz) = 1.d0/alpha_i(p_to_i(p))
        !  end do
        !end do

        !n_tot = sum(n_neigh(1:i))-n_neigh(i)

        d_vec = 0.d0
        do p = 1, n_sub_sites
          do c1 = 1, 3
            d_vec(3*(p-1)+c1,c1) = 1.d0
          end do
        end do


        hirshfeld_v_sub_der = 0.d0
        n_sub_neigh = 0
        sub_neighbors_list = 0
        p_list = 0
        !do j2 = 1, n_neigh(i)
        !  i2 = modulo(neighbors_list(n_tot+j2)-1,n_sites)+1

        !NEIGHBORS_LIST TEST
        !do p = 1, n_sub_sites
        !  i2 = p_to_i(p)
        !write(*,*) "n_neigh", n_neigh(i)
        k_i = 0
        p = 0
        do i3 = 1, n_neigh(i)
          !if (i == 1 .and. i3 == 8) then
          !write(*,*) "k2", k2
          !end if
          k_i = k_i+1
          xyz_i = xyz(1:3,n_tot+k_i)/Bohr
          i2 = neighbors_list(n_tot+k_i)
          if ( sqrt(sum((xyz_i-xyz_l)**2)) < rcut/Bohr ) then
            p = p+1
            if ( dabs(sum(xyz_i-xyz_l)) < 1.d-10 ) then
              l_cent = p
            end if
            !p = i_to_p(i2)
            !write(*,*) "p, i_to_p", p, i_to_p(i2) 
            k2 = k2+1
            !n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
            !s = neighbor_species(n_tot2+1)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = neighbors_list(n_tot+k_i)
            hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            !r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot2+1)**(1.d0/3.d0)
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            !omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
            !xyz_H(:,k2) = xyz(:,n_tot2+1)/Bohr
            !rjs_H(k2) = rjs(n_tot2+1)/Bohr
            xyz_H(:,k2) = xyz(:,n_tot+k_i)/Bohr
            !xyz_i = xyz_H(:,k2)
            rjs_H(k2) = rjs(n_tot+k_i)/Bohr
            rjs_0(k2) = rjs(n_tot+k_i)
            r_vdw_i = r0_ii(k2)
            s_i = neighbor_sigma(k2)
            do c1 = 1, 3
              B_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0/neighbor_alpha0(k2)
              nnz = nnz+1
              ia(nnz) = 3*(p-1)+c1
              ja(nnz) = 3*(p-1)+c1
              val(nnz) = 1.d0/neighbor_alpha0(k2)
            end do
            k_j = 0
            q = 0
            !do j3 = 2, n_neigh(i2)
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              xyz_j = xyz(1:3,n_tot+k_j)/Bohr
              if ( sqrt(sum((xyz_j-xyz_l)**2)) < rcut/Bohr ) then
                !if ( rjs(n_tot+k_j) < rcut) then
                q = q+1
                !end if
                !write(*,*) "q, i_to_p", q, i_to_p(neighbors_list(n_tot+k_j))
                if (i3 .ne. j3) then
                !xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) < rcut/Bohr ) then
                !j = neighbors_list(n_tot+k_j)
                !if ( rjs(n_tot+k_j) < 2*rcut ) then
                  !if (i == 1 .and. i3 == 7) then
                  !write(*,*) "i, j", i2, j
                  !end if
                  n_sub_neigh(p) = n_sub_neigh(p) + 1
                  j = neighbors_list(n_tot+k_j)
                  !q = i_to_p(j)
                  k2 = k2+1    
                  s = neighbor_species(n_tot+k_j)
                  sub_neighbors_list(k2) = j
                  !if ( rjs(n_tot+k_j) < rcut ) then
                  p_list(k2) = q
                  !end if                        
                  r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                  !omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
                  neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                  neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
                  xyz_H(:,k2) = xyz_j-xyz_i
                  !rjs_H(k2) = rjs(n_tot2+j3)/Bohr
                  rjs_H(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                  rjs_0(k2) = rjs(n_tot+k_j)
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
                      !if ( rjs(n_tot+k_j) < rcut ) then
                        !if ( rjs_H(k2) < rcut/Bohr ) then
                          B_mat(3*(p-1)+c1,3*(q-1)+c2) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3))
                          nnz = nnz+1
                          ia(nnz) = 3*(p-1)+c1
                          ja(nnz) = 3*(q-1)+c2
                          val(nnz) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                     g_func(k2) + h_func(k3))
                        !end if
                      !else ! rjs(n_tot+k_j) < rcut
                      !  d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                      !                            g_func(k2) + h_func(k3)) * neighbor_alpha0(k2)
                      !end if
                    end do
                  end do
                end if
              end if
              end if
            end do
          end if
        end do

        !call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, d_vec, 3*n_sub_sites, work_arr, &
        !            12*n_sub_sites, info )

        !call dsytrf( 'U', 3*n_sub_sites, B_mat, 3*n_sub_sites, ipiv, work_arr, 12*n_sub_sites, info)
        call dgetrf( 3*n_sub_sites, 3*n_sub_sites, B_mat, 3*n_sub_sites, ipiv, info )

        call dgetrs( 'N', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, d_vec, 3*n_sub_sites, info )         

        if ( do_timing ) then
        call cpu_time(time2)
        write(*,*) "pre-calculation stuff timing", time2-time1
        call cpu_time(time1)
        end if

        allocate( dT(1:9*n_sub_pairs) )
        allocate( dB_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
        dB_mat = 0.d0
        allocate( b_der(1:3*n_sub_sites,1:3) )
        allocate( f_damp_der(1:n_sub_pairs) )
        allocate( g_func_der(1:n_sub_pairs) )
        allocate( h_func_der(1:9*n_sub_pairs) )
        allocate( d_der(1:3*n_sub_sites,1:3) )        


        if ( do_timing ) then
        call cpu_time(time2)
        write(*,*) "gradient allocation", time2-time1
        end if

        do c3 = 1, 3  ! REMOVE NEIGHBOR-NEIGHBOR LIST DEPENDENCIES INSIDE THIS LOOP!!!!!
          if ( do_timing ) then
          call cpu_time(time1)
          end if
          dT = 0.d0
          k2 = 0
          do p = 1, n_sub_sites
            !i2 = p_to_i(p)
            !if ( in_cutoff(i2) ) then 
            k2 = k2+1
            do j3 = 2, n_sub_neigh(p)
              k2 = k2+1
              !if ( rjs_0(k2) < rcut ) then
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
              !end if
            end do
            !end if
          end do
        !write(*,*) "k2", k2
        !write(*,*) "n_sub_pairs", n_sub_pairs
        
          if ( do_timing ) then
          call cpu_time(time2)
          write(*,*) "dT timing", time2-time1
          end if        

!TEST!!!!!
          !k_a = sum(n_neigh(1:i))-n_neigh(i)
          !do r = 1, n_neigh(i)
            !k_a = k_a+1
            !if ( rjs(k_a) < 2*rcut ) then
            !a = neighbors_list(k_a)
            a = i
            if ( do_timing ) then
            call cpu_time(time1)
            end if
!TEST!!!!!
            !dB_mat = 0.d0
            f_damp_der = 0.d0
            g_func_der = 0.d0
            h_func_der = 0.d0
            b_der = 0.d0
            d_der = 0.d0
            
            if ( do_timing ) then
            call cpu_time(time2)
            write(*,*) "initialization of gradients", time2-time1
            call cpu_time(time1)
            end if
          
            k3 = 0
            do p = 1, n_sub_sites
              k3 = k3+1
              r_vdw_i = r0_ii(k3)
              s_i = neighbor_sigma(k3)
              i2 = sub_neighbors_list(k3)
              do j2 = 2, n_sub_neigh(p)
                k3 = k3+1
                !if ( rjs_0(k3) < rcut ) then
                r_vdw_j = r0_ii(k3)
                s_j = neighbor_sigma(k3)
                j = sub_neighbors_list(k3)
                j3 = modulo(sub_neighbors_list(k3)-1, n_sites) + 1
!TEST!!!!!!!!!!!
                if (a == i2 .or. a == j) then
!                if (i == i2 .or. i == j) then
!TEST!!!!!!!!!!!!
                  !if ( rjs_0(k3) < rcut ) then
                  q = p_list(k3)
                  !end if
                  !sigma_ij = sqrt(sigma_i(i2)**2 + sigma_i(j)**2)
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  f_damp_der(k3) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k3)**2 * &
                                     exp( -d*(rjs_H(k3)/(sR*(r_vdw_i + &
                                     r_vdw_j)) - 1.d0) ) * xyz_H(c3,k3)/rjs_H(k3)
                  g_func_der(k3) = 4.d0/sqrt(pi) * rjs_H(k3)/sigma_ij**3 * xyz_H(c3,k3) * &
                                       exp(-rjs_H(k3)**2/sigma_ij**2)
                  !if ( i == 2 .and. c3 == 1 .and. a == 1 .and. p == 40 .and. q == 41 ) then
                  !  write(*,*) "f_damp_der", f_damp_der(k3)
                  !end if
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
!TEST!!!!!!!!!!!!!!!!
                      !if ( rjs_0(k3) < rcut ) then
                      if (a == i2) then
!                      if ( i == i2 ) then
!TEST!!!!!!!!!!!!!!!
                        b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                          alpha_SCS_full(3*(j3-1)+c2,:,1)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) - (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3)))
                        !dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) - (f_damp_der(k3) * T_func(k4) * &
                        !                                  g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                  g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                  T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                  h_func_der(k4) * (1.d0 - f_damp(k3)))
                      else
                        b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                          alpha_SCS_full(3*(j3-1)+c2,:,1)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + (f_damp_der(k3) * T_func(k4) * &
                                                          g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                          g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                          T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                          h_func_der(k4) * (1.d0 - f_damp(k3)))
                        !dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + (f_damp_der(k3) * T_func(k4) * &
                        !                                  g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                  g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                  T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                  h_func_der(k4) * (1.d0 - f_damp(k3)))                                
                      end if
                      !else
                        !if ( a == i2 ) then
                        !  d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - (f_damp_der(k3) * T_func(k4) * &
                        !                                    g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                    g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                    T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                    h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                        !                                    neighbor_alpha0(k3)
                          !write(*,*) "d_der i2", d_der(3*(p-1)+c1,c2)
                        !else
                        !  d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) + (f_damp_der(k3) * T_func(k4) * &
                        !                                    g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                        !                                    g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                        !                                    T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                        !                                    h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                        !                                    neighbor_alpha0(k3)
                          !write(*,*) "d_der j"
                        !end if
                      !end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. i2 == 2 .and. j == 1 .and. c1 == 1 .and. c2 == 1) then
                      !write(*,*) "T_func", T_func(k4)
                      !write(*,*) "g_func", g_func(k3)
                      !write(*,*) "h_func", h_func(k4)
                      !write(*,*) "f_damp", f_damp(k3)
                      !write(*,*) "f_damp_der", f_damp_der(k3)
                      !write(*,*) "dT", dT(k4)
                      !write(*,*) "g_func_der", g_func_der(k3)
                      !write(*,*) "h_func_der", h_func_der(k4)
                      !end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. p == 1 .and. q == 3 ) then
                      !  write(*,*) "f_damp_der", f_damp_der(k3)
                      !  write(*,*) "g_func_der", g_func_der(k3)
                      !  write(*,*) "h_func_der", h_func_der(k4)
                      !  write(*,*) "dT", dT(k4)
                      !  write(*,*) "T_func, g_func, h_func, f_damp", T_func(k4), g_func(k3), h_func(k4), f_damp(k3)
                      !end if
                      !if ( i == 1 .and. a == 1 .and. c3 == 1 .and. p == 3 .and. q == 1 ) then
                      !  write(*,*) "f_damp_der", f_damp_der(k3)
                      !  write(*,*) "g_func_der", g_func_der(k3)
                      !  write(*,*) "h_func_der", h_func_der(k4)
                      !  write(*,*) "dT", dT(k4)
                      !  write(*,*) "T_func, g_func, h_func, f_damp", T_func(k4), g_func(k3), h_func(k4), f_damp(k3)
                      !end if
                      !write(*,*) "a, i2, j", a, i2, j
                      !write(*,*) "b_der", (f_damp_der(k3) * T_func(k4) * &
                      !                                    g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                      !                                    g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                      !                                    T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                      !                                    h_func_der(k4) * (1.d0 - f_damp(k3)))
                    end do
                  end do
                end if
                !end if
              end do
            end do

            !if ( a == 1 .and. c3 == 1 .and. i == 1) then
            !  open(unit=69, file="dB_mat.dat", status="new")
            !  write(*,*) "dB_mat"
            !  write(*,*) a, c3, i
            !  write(*,*) "before"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) b_der(p,:)
            !  end do
            !  close(69)
            !  write(*,*) "after"
            !  b_der = matmul(dB_mat,alpha_SCS_full(:,:,1))
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) b_der(p,:)
            !  end do
            !end if

            if ( do_timing ) then
            call cpu_time(time2)
            write(*,*) "gradients without hirshfeld timing", time2-time1
            end if            

            !if ( i == 2 .and. a == 1 .and. c3 == 1) then
            !write(*,*) "dB_mat"
            !p = 40
            !write(*,*) p, p_to_i(p)
            !do c1 = 1, 3
            !  write(*,*) dB_mat(3*(p-1)+c1,3*(41-1)+1:3*(41-1)+3)
            !end do
            !end if

            if (do_hirshfeld_gradients) then
           
              if ( do_timing ) then
              call cpu_time(time1)
              end if

              !allocate( coeff_der(1:n_sub_sites,1:n_sub_sites,1:3,1:3) )
              !allocate( coeff_fdamp(1:n_sub_sites,1:n_sub_sites,1:3,1:3) )
              allocate( coeff_der(1:9*n_sub_pairs) )
              allocate( coeff_fdamp(1:9*n_sub_pairs) )
              
              !in_force_cutoff = .false.
              !k3 = sum(n_neigh(1:a))-n_neigh(a)
              !do j2 = 1, n_neigh(a)
              !  k3 = k3+1
              !  j = neighbors_list(k3)
              !  if ( rjs(k3) < 2 * rcut ) then
              !    in_force_cutoff(j) = .true.
              !  end if
              !end do
              !xyz_a = xyz(:,k_a)/Bohr
            
              coeff_fdamp = 0.d0
              coeff_der = 0.d0
              
              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "allocation and force cutoff array timing", time2-time1
              call cpu_time(time1)
              end if
              
              k3 = 0
              do p = 1, n_sub_sites
                k3 = k3+1
                r_vdw_i = r0_ii(k3)
                s_i = neighbor_sigma(k3)
                i2 = sub_neighbors_list(k3)
                do j2 = 2, n_sub_neigh(p)
                  k3 = k3+1
                  !if ( rjs_0(k3) < rcut ) then
                  r_vdw_j = r0_ii(k3)
                  s_j = neighbor_sigma(k3)
                  j = sub_neighbors_list(k3)
                  !sigma_ij = sqrt(sigma_i(i2)**2 + sigma_i(j)**2)
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  S_vdW_ij = sR*(r_vdw_i + r_vdw_j)
                  exp_term = exp(-d*(rjs_H(k3)/S_vdW_ij - 1.d0))
                  k4 = 9*(k3-1)
                  q = p_list(k3)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k4 = k4+1
                      !coeff_fdamp(p,q,c1,c2) = (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                      !                         (-T_func(k4) * g_func(k3) + h_func(k4))
                      coeff_fdamp(k4) = (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                               (-T_func(k4) * g_func(k3) + h_func(k4))
                      dg = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * (rjs_H(k3)/sigma_ij)**2 * &
                           (-rjs_H(k3)/(3.d0 * sigma_ij**3))
                      dh = 4.d0/sqrt(pi) * exp(-rjs_H(k3)**2/sigma_ij**2) * xyz_H(c1,k3)*xyz_H(c2,k3) / &
                           (sigma_ij**5 * rjs_H(k3)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k3)/sigma_ij)**2)
                      !coeff_der(p,q,c1,c2) = (1.d0-f_damp(k3)) * (-T_func(k4)*dg+dh)
                      coeff_der(k4) = (1.d0-f_damp(k3)) * (-T_func(k4)*dg+dh)
                      !if ( i == 2 .and. p == 40 .and. q == 41 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                      !       .and. c2 == 2 ) then
                      !  write(*,*) "coeff_fdamp first", (d*sR*rjs_H(k3))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2
                      !end if
                    end do
                  end do
                  !end if
                end do
              end do

              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "calculation of coeff_fdamp and coeff_der timing", time2-time1
              call cpu_time(time1)
              end if

              k3 = 0
              do p = 1, n_sub_sites
                !write(*,*) k3, "/", n_sub_pairs
                i2 = sub_neighbors_list(k3+1)
                !write(*,*) "I2", i2
                !n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
                !xyz_i = xyz_H(:,k3+1)
                !if ( in_cutoff(a) .and. in_force_cutoff(i2) ) then
                !if ( sqrt(sum((xyz_a-xyz_i)**2)) < rcut/Bohr ) then ! This should be SOAP cut-off (at the moment rcut_scs = rcut_soap)!
                if ( dabs(hirshfeld_v_sub_der(c3,p)) > 1.d-10 ) then
                  r_vdw_i = r0_ii(k3+1)
                  s_i = neighbor_sigma(k3+1)
                  do c1 = 1, 3
                    !dB_mat(3*(p-1)+c1,3*(p-1)+c1) = - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v(i2)) * &
                    !                         hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr
                    b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v_neigh(k3+1)) * &
                                          hirshfeld_v_sub_der(c3,p)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                                          !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                    dB_mat(3*(p-1)+c1,3*(p-1)+c1) = dB_mat(3*(p-1)+c1,3*(p-1)+c1) - 1.d0/(neighbor_alpha0(k3+1) &
                                          * hirshfeld_v_neigh(k3+1)) * &
                                          hirshfeld_v_sub_der(c3,p)*Bohr
                  end do
                  !j3 = findloc(neighbors_list(n_tot2+1:n_tot2+n_neigh(i2)),a,1)
                  do j2 = 2, n_sub_neigh(p)
                    j = sub_neighbors_list(k3+j2)
                    !if ( rjs_0(k3+j2) < rcut ) then
                    q = p_list(k3+j2)
                    !write(*,*) "k3+j2", k3+j2, "/", n_sub_pairs
                    !write(*,*) "j, q", j, q
                    !end if
                    !if ( i2 == j ) then
                    !  do c1 = 1, 3
                        !dB_mat(3*(p-1)+c1,3*(p-1)+c1) = - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v(i2)) * &
                        !                         hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr
                    !    b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_v_neigh(k3+1)) * &
                    !                          hirshfeld_v_sub_der(c3,p)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                                              !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr * alpha_SCS_full(3*(i2-1)+c1,:,1)
                    !  end do
                    !end if
                    k4 = 9*(k3+j2-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k4 = k4+1
                        !if ( rjs_0(k3+j2) < rcut ) then                      
                        b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
                          hirshfeld_v_sub_der(c3,p)*Bohr) * &
                          alpha_SCS_full(3*(j-1)+c2,:,1)
                        dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
                          hirshfeld_v_sub_der(c3,p)*Bohr)  

                          !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr) * &
                          !alpha_SCS_full(3*(j-1)+c2,:,1)



                        !b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + &
                        !  ((coeff_der(p,q,c1,c2) * sigma_i(i2)**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(p,q,c1,c2) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)) * &
                        !  alpha_SCS_full(3*(j-1)+c2,:)
                        
                        !dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + &
                        !  ((coeff_der(k4) * s_i**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(k4) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr)
                        b_der(3*(q-1)+c1,:) = b_der(3*(q-1)+c1,:) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
                          hirshfeld_v_sub_der(c3,p)*Bohr) * &
                          alpha_SCS_full(3*(i2-1)+c2,:,1)
                        dB_mat(3*(q-1)+c1,3*(p-1)+c2) = dB_mat(3*(q-1)+c1,3*(p-1)+c2) + &
                          ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
                          hirshfeld_v_sub_der(c3,p)*Bohr)



                      !  if ( i == 1 .and. c3 == 1 .and. l == 1 .and. p == 50 .and. q == 48 .and. c1 == 3 .and. c2 == 3 ) then
                      !  write(*,*) "coeff_der", coeff_der(k4)
                      !  write(*,*) "s_i", s_i
                      !  write(*,*) "r_vdw_i", r_vdw_i
                      !  write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh(k3+1)
                      !  write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                      !  write(*,*) "hirshfeld_v_sub_der", hirshfeld_v_sub_der(c3,p)*Bohr
                      !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                      !  end if

                      !  if ( i == 1 .and. c3 == 1 .and. l == 1 .and. p == 48 .and. q == 50 .and. c1 == 3 .and. c2 == 3 ) then
                      !  write(*,*) "coeff_der", coeff_der(k4)
                      !  write(*,*) "s_i", s_i
                      !  write(*,*) "r_vdw_i", r_vdw_i
                      !  write(*,*) "hirshfeld_v_neigh", hirshfeld_v_neigh(k3+1)
                      !  write(*,*) "coeff_fdamp", coeff_fdamp(k4)
                      !  write(*,*) "hirshfeld_v_sub_der", hirshfeld_v_sub_der(c3,p)*Bohr
                      !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                      !  end if


                          !hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr) * &
                          !alpha_SCS_full(3*(i2-1)+c2,:,1)
          !              else
          !              d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) + &
          !                ((coeff_der(k4) * s_i**2/hirshfeld_v_neigh(k3+1) + &
          !                coeff_fdamp(k4) * r_vdw_i/hirshfeld_v_neigh(k3+1)) * &
          !                hirshfeld_v_sub_der(c3,p)*Bohr) * &
          !                neighbor_alpha0(k3+j2)
                        !end if

                        !b_der(3*(q-1)+c1,:) = b_der(3*(q-1)+c1,:) + &
                        !  ((coeff_der(q,p,c1,c2) * sigma_i(i2)**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(q,p,c1,c2) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)) * &
                        !  alpha_SCS_full(3*(i2-1)+c2,:)
                          
                        !dB_mat(3*(q-1)+c1,3*(p-1)+c2) = dB_mat(3*(q-1)+c1,3*(p-1)+c2) + &
                        !  ((coeff_der(k4) * s_i**2/hirshfeld_v(i2) + &
                        !  coeff_fdamp(k4) * r_vdw_i/hirshfeld_v(i2)) * hirshfeld_v_cart_der(c3,n_tot2+j3)*Bohr) 
                        !if ( i == 2 .and. p == 40 .and. q == 41 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                        !     .and. c2 == 2) then
                        !  write(*,*) "p, q, i2, j", p, q, i2, j
                        !  write(*,*) "coeff_fdamp", coeff_fdamp(p,q,c1,c2) !/(-T_func(k4) * g_func(k3) + h_func(k4)) !* &
                        !                            !r_vdw_i/hirshfeld_v(i2) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !  write(*,*) "dB_mat", dB_mat(3*(p-1)+c1,3*(q-1)+c2)
                        !  write(*,*) "hirshfeld_v_cart_der term", r_vdw_i/hirshfeld_v(i2) * &
                        !                  hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !end if
                        !if ( i == 2 .and. p == 41 .and. q == 40 .and. a == 1 .and. c3 == 1 .and. c1 == 3 &
                        !     .and. c2 == 2 ) then
                        !  write(*,*) "p, q", p, q
                        !  write(*,*) "coeff_fdamp", coeff_fdamp(q,p,c1,c2) !/(-T_func(k4) * g_func(k3) + h_func(k4)) * &
                        !                            !r_vdw_i/hirshfeld_v(i2) * hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !  write(*,*) "dB_mat", dB_mat(3*(q-1)+c1,3*(p-1)+c2)
                        !  write(*,*) "hirshfeld_v_cart_der term", r_vdw_i/hirshfeld_v(i2) * &
                        !                 hirshfeld_v_cart_der_H(c3,n_tot2+j3)
                        !end if
                      end do
                    end do
                  end do
                end if
                k3 = k3 + n_sub_neigh(p)
              end do
              
              if ( do_timing ) then
              call cpu_time(time2)
              write(*,*) "hirshfeld gradient timing", time2-time1
              end if            

              deallocate( coeff_der, coeff_fdamp )
              
            end if

            if ( do_timing ) then
            call cpu_time(time1)
            end if

            !if ( i == 1 .and. c3 == 1 .and. l == 1 ) then
            !  write(*,*) "dB_mat"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) dB_mat(p,:)
            !  end do
            !end if


            !if ( i == 1 .and. c3 == 1 .and. l == 1 ) then

            !write(*,*) "d_der", d_der

            !write(*,*) "b_der"
            !do p = 1, 3*n_sub_sites
            !  write(*,*) b_der(p,:)
            !end do

            !end if

            b_der = b_der+d_der

            !call dsytrs( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, b_der, 3*n_sub_sites, work_arr, &
            !             12*n_sub_sites, info )
            call dgetrs( 'N', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, b_der, 3*n_sub_sites, info )


            ! MBD forces

            k2 = 0
            do p = 1, n_sub_sites
              i2 = sub_neighbors_list(k2+1)
              !write(*,*) "i2", i2
              if ( p == l_cent ) then
                q = findloc(mbd_neighbors_list(1:n_mbd_neigh(1)),i2,1)
                do c1 = 1, 3
                  alpha_grad(c3,q) = alpha_grad(c3,q) - 1.d0/3.d0 * b_der(3*(p-1)+c1,c1)
                end do
                if ( i == 2 .and. c3 == 1 ) then
                !  write(*,*) "l_cent", l_cent
                !  write(*,*) "alpha_grad", i2, alpha_grad(c3,q)
                end if
              end if
              k2 = k2+n_sub_neigh(p)
            end do

        end do ! c3 loop

        deallocate( dT, b_der, g_func_der, h_func_der, hirshfeld_v_sub_der, d_der, f_damp_der )
        deallocate( dB_mat )

        deallocate( n_sub_neigh, sub_neighbors_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, T_func, &
                    g_func, h_func, f_damp, p_list )

        deallocate( B_mat, ipiv, work_arr, ia, ja, val, rjs_0, d_vec )
        
        end if ! if (rjs(n_tot+l) < rcut)
        
        end do ! l loop

        call cpu_time(time2)
        !write(*,*) "polarizability gradient timing", time2-time1
        call cpu_time(time1)

        do c3 = 1, 3

            f_damp_der_SCS = 0.d0
            f_damp_der_mbd = 0.d0 ! This is cleared so we can recalculate it with SCS values
            dT_mbd = 0.d0
            dT_LR = 0.d0
            k2 = 0
            do p = 1, n_mbd_sites
              k2 = k2+1
              r_vdw_i = r0_ii_SCS(k2)
              i2 = mbd_neighbors_list(k2)
              do j3 = 2, n_mbd_neigh(p)
                k2 = k2+1
                j = mbd_neighbors_list(k2)
                q = p_mbd(k2)
                r_vdw_j = r0_ii_SCS(k2)
                R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                S_vdW_ij = sR*R_vdW_SCS_ij
                dS_vdW_ij = sR/3.d0 * ( r_vdw_i/alpha_SCS0(i2,1) * alpha_grad(c3,p) + &
                                        r_vdw_j/alpha_SCS0(j,1) * alpha_grad(c3,q) )
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

            !if ( i == 1 .and. c3 == 1 ) then
            !  write(*,*) "dT_LR"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) dT_LR(p,:)
            !  end do
            !end if

            G_mat = 0.d0

            k3 = 0
            do p = 1, n_mbd_sites
              i2 = mbd_neighbors_list(k3+1)
              G_mat(3*(p-1)+1:3*(p-1)+3,:,1) = G_mat(3*(p-1)+1:3*(p-1)+3,:,1) + &
                alpha_SCS0(i2,1) * dT_LR(3*(p-1)+1:3*(p-1)+3,:) + &
                alpha_grad(c3,p) * T_LR(3*(p-1)+1:3*(p-1)+3,:)
              k3 = k3+n_mbd_neigh(p)
            end do


            if ( n_order > 1 ) then


              integrand = 0.d0
              total_integrand = 0.d0
              do om = 1, n_freq

                force_series = 0.d0
                total_energy_series = 0.d0
                do k2 = 1, n_order-1
                  force_series = force_series + AT_n_f(:,:,k2,1) * &
                                 (1.d0/(1.d0 + (omegas(om)/0.5d0)**2))**k2 !/alpha_SCS0(i,3))**2))**k2
                  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2,1) * &
                                 (1.d0/(1.d0 + (omegas(om)/0.5d0)**2))**k2 !/alpha_SCS0(i,3))**2))**k2
                end do
                k3 = 0
                do p = 1, n_mbd_sites
                  i2 = mbd_neighbors_list(k3+1)
                  do c1 = 1, 3
                    integrand(om) = integrand(om) + 1.d0/(1.d0 + (omegas(om)/0.5d0)**2) * &
                    dot_product(G_mat(3*(p-1)+c1,:,1),force_series(:,3*(p-1)+c1))
                    total_integrand(om) = total_integrand(om) + alpha_SCS0(i2,1) / &
                          (1.d0 + (omegas(om)/0.5d0)**2) & !/alpha_SCS0(i,3))**2) &
                          * dot_product(T_LR(3*(p-1)+c1,:), &
                          total_energy_series(:,3*(p-1)+c1))
                  end do
                  k3 = k3 + n_mbd_neigh(p)
                end do
              end do

              integral = 0.d0
              call integrate("trapezoidal", omegas, integrand, omegas(1), omegas(n_freq), integral)
              forces0(c3,i) = forces0(c3,i) + 1.d0/(2.d0*pi) * integral
              forces0(c3,i) = forces0(c3,i) * 51.42208619083232

             ! write(*,*) "MBD force", i, c3, forces0(c3,i)
              integral = 0.d0
              call integrate("trapezoidal", omegas, total_integrand, omegas(1), omegas(n_freq), integral)
             ! write(*,*) "MBD total energy of sphere", i, integral / (2.d0*pi) * 27.211386245988

            end if
           
            
        end do ! c3
            
        call cpu_time(time2)
        !write(*,*) "MBD forces timing", time2-time1
            
        end if ! do_derivatives
          

        if (derivative_expansion) then
          deallocate( B_inv )
        end if
        ! MBD stuff
        deallocate( T_LR, r0_ii_SCS, f_damp_SCS, AT, AT_n, energy_series, omegas, integrand, n_mbd_neigh, &
                    mbd_neighbors_list, p_mbd, r0_ii_mbd, neighbor_alpha0_mbd, xyz_mbd, rjs_mbd, T_mbd )    
        if ( do_derivatives ) then
          deallocate( AT_n_f, f_damp_der_mbd, f_damp_der_SCS, dT_LR, G_mat, force_series, VL, total_energy_series, &
                      total_integrand, alpha_grad, dT_mbd ) 
        end if

               
      end do

      !write(*,*) "dalpha w.r.t. atom 1 in x direction"
      !do p = 1, n_sites
      !  write(*,*) p, dalpha_full(p,1,1,1)
      !end do

      end if ! do_derivatives
      
      deallocate( central_pol, central_omega )
      deallocate( alpha_SCS_full, in_cutoff, p_to_i, i_to_p, A_i, hirshfeld_v_cart_der_H )
      !deallocate( alpha_i, sigma_i, hirshfeld_v_cart_der_H )
      if ( do_derivatives .and. do_hirshfeld_gradients ) then
        deallocate( in_force_cutoff )
      end if
      
    else ! local

!   This implementation assumes that rcut is the largest cutoff, that is, the neigbhbors_list contains only the atoms within the vdW cutoff.
!   The implementation matches with the implementation above if rcut is the largest cutoff or the cutoff is so small that the only neighbor
!   the atoms see are themselves (n_neigh(i) = 1 for all i).

    if( do_timing ) then
      call cpu_time(time1)
    end if

!   Allocate all the necessary stuff
    allocate( neighbor_c6_ii(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( f_damp(1:n_pairs) )
    allocate( T_func(1:9*n_pairs) )
    allocate( h_func(1:9*n_pairs) )
    allocate( g_func(1:n_pairs) )
    allocate( omegas(1:n_freq) )
    allocate( omega_i(1:n_pairs) )
    allocate( sigma_i(1:n_sites) )
    allocate( alpha_i(1:n_sites) )
!    allocate( T_SR(1:3*n_sites,1:3*n_sites) )
    allocate( B_mat(1:3*n_sites,1:3*n_sites) )
    allocate( xyz_H(1:3,1:n_pairs) )
    allocate( rjs_H(1:n_pairs) )
    allocate( work_arr(1:12*n_sites) )
    allocate( ipiv(1:3*n_sites) )
    allocate( I_mat(1:3*n_sites,1:3*n_sites) )
    allocate( a_vec(1:3*n_sites,1:3) )
    allocate( a_SCS(1:3*n_sites,1:3) )
!    allocate( BTB_reg_copy(1:3*n_sites,1:3*n_sites) )
    allocate( BTB_reg(1:3*n_sites,1:3*n_sites) )
    allocate( B_reg(1:3*n_sites,1:3*n_sites) )
    allocate( c6_nsites(1:n_sites) )
    if ( do_derivatives ) then
      allocate( dT(1:9*n_pairs) )
      !allocate( dT_SR(1:3*n_sites,1:3*n_sites) )
      allocate( f_damp_der(1:n_pairs) )
      allocate( g_func_der(1:n_pairs) )
      allocate( h_func_der(1:9*n_pairs) )
      !allocate( dT_SR_v(1:3*n_sites,1:3*n_sites) )
      !allocate( coeff_der(1:n_sites,1:n_sites,1:3,1:3) )   ! THIS IS CHANGED RECENTLY, MIGHT CHANGE THE VALUES!!!!!!!!!!!!
      !allocate( coeff_fdamp(1:n_sites,1:n_sites,1:3,1:3) ) ! THIS TOO!!!
      allocate( coeff_der(1:9*n_pairs) )
      allocate( coeff_fdamp(1:9*n_pairs) )
      allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
      allocate( da_vec(1:3*n_sites,1:3) )
      allocate( vect1(1:3*n_sites,1:3) )
      allocate( vect2(1:3*n_sites,1:3) )
      allocate( vect3(1:3*n_sites,1:3) )
      allocate( vect4(1:3*n_sites,1:3) )
      allocate( vect_temp(1:3*n_sites,1:3) )
      allocate( da_SCS(1:3*n_sites,1:3) )
      allocate( b_der(1:3*n_sites,1:3) )
      dalpha_full = 0.d0
    end if

!   Frequencies used for integration:
    !omega = 0.d0
    !do i = 1, n_freq
    !  omegas(i) = omega
    !  omega = omega + 0.4d0 
    !end do

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
    
    omegas(1) = 0.d0
    omegas(2) = 2.d0 * omega_i(1)

    write(*,*) "omega_i", omega_i(1)

!   Identity matrix
    I_mat = 0.d0
    do p = 1, n_sites
      do c1 = 1, 3
        I_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0
      end do
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

    if( do_timing) then
      call cpu_time(time2)
      write(*,*) "Timing for initial stuff:", time2-time1
    end if
    
    allocate( ia(1:9*n_sites*n_sites) )
    allocate( ja(1:9*n_sites*n_sites) )
    allocate( val(1:9*n_sites*n_sites) )

    do om = 1, 2

      write(*,*) "Doing frequency", om

      if( do_timing) then
        call cpu_time(time1)
      end if

      k2 = 1
      do i = 1, n_sites
        alpha_i(i) = neighbor_alpha0(k2)/(1.d0 + omegas(om)**2/omega_i(k2)**2)
        sigma_i(i) = (sqrt(2.d0/pi) * alpha_i(i)/3.d0)**(1.d0/3.d0)
        k2 = k2+n_neigh(i)
      end do

      g_func = 0.d0
      h_func = 0.d0
      k = 0
      do i = 1, n_sites
        k = k+1
        s_i = sigma_i(i)
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          s_j = sigma_i(j)
          if( rjs(k) < rcut )then
            sigma_ij = sqrt(s_i**2 + s_j**2)
            g_func(k) = erf(rjs_H(k)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k)/sigma_ij) * exp(-rjs_H(k)**2.d0/sigma_ij**2)
            k2 = 9*(k-1)
            do c1 = 1, 3
              do c2 = 1, 3
                k2 = k2+1
                h_func(k2) = 4.d0/sqrt(pi) * (rjs_H(k)/sigma_ij)**3 * &
                                      xyz_H(c1,k)*xyz_H(c2,k)/rjs_H(k)**5 * exp(-rjs_H(k)**2/sigma_ij**2)
              end do
            end do 
          end if
        end do
      end do

!      T_SR = 0.d0
      B_mat = 0.d0
      k = 0
      nnz = 0
      do i = 1, n_sites
        k = k+1
        do c1 = 1, 3
          !B_mat(3*(i-1)+c1,3*(i-1)+c1) = 1.d0
          B_mat(3*(i-1)+c1,3*(i-1)+c1) = 1.d0/alpha_i(i)
          nnz = nnz+1
          ia(nnz) = 3*(i-1)+c1
          ja(nnz) = 3*(i-1)+c1
          val(nnz) = B_mat(3*(i-1)+c1,3*(i-1)+c1)
        end do
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          !write(*,*) "j", j
          !j = modulo(j-1,n_sites)+1
          !write(*,*) "j2", j
          if( rjs(k) < rcut )then
            k2 = 9*(k-1)
            do c1 = 1, 3
              do c2 = 1, 3
                k2 = k2+1
                nnz = nnz+1
                !B_mat(3*(i-1)+c1,3*(j-1)+c2) = alpha_i(i) * (1.d0-f_damp(k)) * (-T_func(k2) * &
                !                                  g_func(k) + h_func(k2))
                B_mat(3*(i-1)+c1,3*(j-1)+c2) = (1.d0-f_damp(k)) * (-T_func(k2) * &
                                                  g_func(k) + h_func(k2))
                ia(nnz) = 3*(i-1)+c1
                ja(nnz) = 3*(j-1)+c2
                val(nnz) = B_mat(3*(i-1)+c1,3*(j-1)+c2)
              end do
            end do
          end if
        end do
      end do

        !if (om == 1) then
        !allocate( eigval(1:3*n_sites) )
        !call dsyev('n', 'u', 3*n_sites, B_mat, 3*n_sites, eigval, work_arr, 12*n_sites, info)
        !write(*,*) "i, min eig", i, minval(eigval)
        !open(unit=79, file="eigs_aC.dat", status="new")
        !write(*,*) "B_mat eigenvalues"
        !write(79,*) eigval
        !close(79)
        !deallocate( eigval)
        !end if

          !write(*,*) "B_mat", B_mat(1,:)
          write(*,*) "Solving full polarizabilities"
          !allocate( alpha_test(1:3*n_sites,1:3) )
          a_SCS = 0.d0
          !write(*,*) "n_sub_sites", n_sub_sites
          do p = 1, n_sites
            do c1 = 1, 3
              a_SCS(3*(p-1)+c1,c1) = 1.d0
            end do
          end do
          call dsysv('U', 3*n_sites, 3, B_mat, 3*n_sites, ipiv, a_SCS, 3*n_sites, &
                      work_arr, 12*n_sites, info)
          do p = 1, n_sites
            alpha_SCS0(p,om) = 1.d0/3.d0 * (a_SCS(3*(p-1)+1,1) &
                       + a_SCS(3*(p-1)+2,2) + a_SCS(3*(p-1)+3,3))
            write(*,*) p, alpha_SCS0(p,om)
          end do
          !deallocate( alpha_test )
          !write(*,*) "B_mat after", B_mat(1,:)
      

      !allocate( eigval(1:3*n_sites) )

      !if ( om == 1 ) then
      !  write(*,*) hirshfeld_v
      !  call dsyev('v', 'u', 3*n_sites, B_mat, 3*n_sites, eigval, work_arr, 12*n_sites, info)
      !  open(unit=79, file="eigs.dat", status="new")
      !  open(unit=89, file="eigvec.dat",status="new")
      !  write(*,*) "B_mat"
      !  write(79,*) eigval
      !  do p = 1, 3*n_sites
      !    write(89,*) B_mat(p,:)
      !  end do
      !  close(79)
      !  close(89)
      !end if

      !deallocate( eigval )

      !write(*,*) "Matrix norm"
      !do p = 1, 3*n_sites
      !  write(*,*) sum(dabs(-B_mat(p,:)+I_mat(p,:)))
      !end do

      !write(*,*) "Writing B_mat"
      !open(unit=79, file="B_mat.dat", status="new")
      !do p = 1, 3*n_sites
      !  write(79,*) B_mat(p,:)
      !end do
      !close(79)
      !write(*,*) "Writing done"

!      B_mat = B_mat + T_SR

      !EIGENVALUE TEST
      !if ( .false. ) then

      if ( psblas ) then

              !polyfit = (/ -0.011111d0,    0.10555365d0, -0.42261792d0,  0.93352476d0, -1.27295025d0,  1.20591827d0, &
              !   -1.00960849d0,  0.96270989d0, -0.99307669d0,  1.00191393d0, -1.00020537d0,  0.99999289d0 /)

       ! polyfit = (/ -8.52126705e+11,  1.54164485e+12, -1.22546032e+12,  5.62886097e+11, &
       !-1.65311654e+11,  3.24458270e+10, -4.32224085e+09,  3.89166719e+08, &
       !-2.31700842e+07,  8.74333211e+05, -1.94556850e+04,  2.24501478e+02 /)

       !polyfit = (/ -11465.40348552d0, 6509.62160556d0, -1165.56417327d0, 70.55762083d0 /)

      allocate( val_xv(1:3*n_sites,1:3) )
      allocate( val_bv(1:3*n_sites,1:3) )
      allocate( myidx(1:3*n_sites) )
      val_xv = 0.d0
      val_bv = 0.d0
      k2 = 0
      do c1 = 1, 3
        do p = 1, n_sites
          k2 = k2+1
          myidx(k2) = k2
          val_xv(3*(p-1)+c1,c1) = 1.d0
          val_bv(3*(p-1)+c1,c1) = 1.d0
        end do
      end do

      call cpu_time(time1)
      call psb_init(icontxt)
      a_SCS = polyfit(n_degree+1)*val_bv
      call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
        !write(*,*) "cdall", info_psb
      call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
        !write(*,*) "spall", info_psb
        !call psb_geall(x_vec,desc_a,info_psb)
        !write(*,*) "geall x", info_psb
        !call psb_geall(b_vec,desc_a,info_psb)
        !write(*,*) "geall b", info_psb
      call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
        !write(*,*) "spins", info_psb
        !call psb_geins(3*n_sites, myidx, val_xv(:,c1), x_vec, desc_a, info_psb)
        !write(*,*) "x_vec", info_psb
        !call psb_geins(3*n_sites, myidx, val_bv(:,c1), b_vec, desc_a, info_psb)
        !write(*,*) "b_vec", info_psb
      call psb_cdasb(desc_a, info_psb)
        !write(*,*) "cdasb", info_psb
      call psb_spasb(A_sp, desc_a, info_psb)
      call cpu_time(time1)
      do k2 = 1, n_degree
        call psb_spmm(1.d0, A_sp, val_bv, 0.d0, val_xv, desc_a, info_psb)
        write(*,*) "psb_spmm", info_psb
        a_SCS = a_SCS + polyfit(n_degree+1-k2) * val_xv 
        val_bv = val_xv
        !write(*,*) "spasb", info_psb
        !call psb_geasb(x_vec, desc_a, info_psb)
        !write(*,*) "geasb x", info_psb
        !call psb_geasb(b_vec, desc_a, info_psb)
        !write(*,*) "geasb b", info_psb
        !ptype="DIAG"
        ! NOTE: Everything works fine until preconditioner has to be set. Then the compilation fails.
        !call prec%init(icontxt, ptype, info_psb)
        !write(*,*) "prec init", info_psb
        !call prec%build(A_sp, desc_a, info_psb)
        !write(*,*) "prec build", info_psb
        !call psb_krylov("BICGSTAB", A_sp, prec, b_vec, x_vec, 0.000001d0, desc_a, info_psb)
        !write(*,*) "krylov", info_psb
        !val_xv(:,c1) = x_vec%get_vect()
      end do
      call cpu_time(time2)
      write(*,*) "psb_spmm timing", time2-time1
      !call psb_exit(icontxt)
      !write(*,*) "val_xv"
      !a_SCS = val_xv
      !alpha_SCS0(:,om) = 0.d0
      do p = 1, n_sites
        alpha_SCS0(p,om) = 1.d0/3.d0 * (a_SCS(3*(p-1)+1,1) + a_SCS(3*(p-1)+2,2) + a_SCS(3*(p-1)+3,3))
      end do
      deallocate( val_xv, val_bv, myidx )
      !call psb_exit(icontxt)
      !deallocate( ia, ja, val )
      call cpu_time(time2)
      write(*,*) "Timing for PSBLAS", time2-time1
      write(*,*) alpha_SCS0(:,om)
      else ! psblas

      !n_iter = 100
      
      !a_vec = 0.d0
      !do p = 1, n_sites
      !  do c1 = 1, 3
      !    !a_vec(3*(p-1)+c1,c1) = alpha_i(p)
      !    a_vec(3*(p-1)+c1,c1) = 1.d0
      ! end do
      !end do

      !if ( om == 1 ) then
      !write(*,*) "a_vec"
      !do p = 1, 3*n_sites
      !  write(*,*) a_vec(p,:)
      !end do
      !end if

      !if( do_timing) then
      !  call cpu_time(time2)
      !  write(*,*) "Energies: timing for everything else:", time2-time1
      !  call cpu_time(time1)
      !end if     

      !if ( regularization ) then
      !  reg_param = 0.01d0
      !  B_reg = B_mat + reg_param * I_mat
      !  call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_reg, 3*n_sites, &
      !              a_vec, 3*n_sites, 0.d0, a_SCS, 3*n_sites)
      !  call dgemm('t', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, B_mat, 3*n_sites, &
      !              B_mat, 3*n_sites, 0.d0, BTB_reg, 3*n_sites)
      !  BTB_reg = BTB_reg + reg_param * I_mat
!     !   BTB_reg_copy = BTB_reg
      !  call dsysv('U', 3*n_sites, 3, BTB_reg, 3*n_sites, ipiv, &
      !              a_SCS, 3*n_sites, work_arr, 12*n_sites, info)
      !else
      !  BTB_reg = B_mat
      !  a_SCS = a_vec
      !  call dgesv(3*n_sites, 3, BTB_reg, 3*n_sites, ipiv, &
      !              a_SCS, 3*n_sites, info)
      !end if

!      if ( om == 1 ) then
!      write(*,*) "a_SCS"
!      do i = 1, 3*n_sites
!        write(*,*) a_SCS(i,:)
!      end do

!      write(*,*) "B_mat"
!      do i = 1, 6
!        write(*,*) B_mat(i,1:6)
!      end do
!      end if

      !alpha_SCS0(:,om) = 0.d0
      !do p = 1, n_sites
      !  do c1 = 1, 3
      !    alpha_SCS0(p,om) = alpha_SCS0(p,om) + a_SCS(3*(p-1)+c1,c1)
      !  end do
      !end do
      !alpha_SCS0(:,om) = alpha_SCS0(:,om)/3.d0
      
      if( do_timing) then
        call cpu_time(time2)
        write(*,*) "Energies: timing for solving alpha_SCS:", time2-time1
      end if
      
      end if ! psblas

      if ( do_derivatives ) then

        if( do_timing) then
          call cpu_time(time1)
        end if

        do c3 = 1, 3
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
                    dT(k2) = (-15.d0 * xyz_H(c1,k) * xyz_H(c2,k) * xyz_H(c3,k))/rjs_H(k)**7
                    if (c1 == c2) then
                      dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c3,k)
                    end if
                    if (c2 == c3) then
                      dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c1,k)
                    end if
                    if (c1 == c3) then
                      dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c2,k)
                    end if
                  end do
                end do
              end if
            end do
          end do

          do a = 1, n_sites

            write(*,*) "frequency, cartesian, atom:", om, c3, a

            f_damp_der = 0.d0
            g_func_der = 0.d0
            h_func_der = 0.d0
            !dT_SR = 0.d0
            b_der = 0.d0
          
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
                    sigma_ij = sqrt(sigma_i(i)**2 + sigma_i(j)**2)
                    f_damp_der(k) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k)**2 * &
                                       exp( -d*(rjs_H(k)/(sR*(r_vdw_i + &
                                       r_vdw_j)) - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
                    g_func_der(k) = 4.d0/sqrt(pi) * rjs_H(k)/sigma_ij**3 * xyz_H(c3,k) * &
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
                        h_func_der(k2) = coeff_h_der * terms
                        !dT_SR(3*(i-1)+c1,3*(j-1)+c2) = f_damp_der(k) * T_func(k2) * &
                        !                                    g_func(k) - (1.d0 - f_damp(k)) * dT(k2) * &
                        !                                    g_func(k) - g_func_der(k) * (1.d0 - f_damp(k)) * &
                        !                                    T_func(k2) - f_damp_der(k) * h_func(k2) + &
                        !                                    h_func_der(k2) * (1.d0 - f_damp(k))
                        if (a == i) then
                          !b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) + alpha_i(i) * (f_damp_der(k) * T_func(k2) * &
                          !                                  g_func(k) - (1.d0 - f_damp(k)) * dT(k2) * &
                          !                                  g_func(k) - g_func_der(k) * (1.d0 - f_damp(k)) * &
                          !                                  T_func(k2) - f_damp_der(k) * h_func(k2) + &
                          !                                  h_func_der(k2) * (1.d0 - f_damp(k))) * a_SCS(3*(j-1)+c2,:)
                          b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) + (f_damp_der(k) * T_func(k2) * &
                                                            g_func(k) - (1.d0 - f_damp(k)) * dT(k2) * &
                                                            g_func(k) - g_func_der(k) * (1.d0 - f_damp(k)) * &
                                                            T_func(k2) - f_damp_der(k) * h_func(k2) + &
                                                            h_func_der(k2) * (1.d0 - f_damp(k))) * a_SCS(3*(j-1)+c2,:)
                        else
                          b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) - (f_damp_der(k) * T_func(k2) * &
                                                            g_func(k) - (1.d0 - f_damp(k)) * dT(k2) * &
                                                            g_func(k) - g_func_der(k) * (1.d0 - f_damp(k)) * &
                                                            T_func(k2) - f_damp_der(k) * h_func(k2) + &
                                                            h_func_der(k2) * (1.d0 - f_damp(k))) * a_SCS(3*(j-1)+c2,:)
                        end if
                      end do
                    end do
                    !if (a == i) then
                    !   dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = -dT_SR(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
                    !end if
                  end if
                end if
              end do
            end do

            if (do_hirshfeld_gradients) then
              !dT_SR_v = 0.d0
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
                    sigma_ij = sqrt(sigma_i(i)**2 + sigma_i(j)**2)
                    S_vdW_ij = sR*(r_vdw_i + r_vdw_j)
                    exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
                    k2 = 9*(k-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k2 = k2+1
                        !coeff_fdamp(i,j,c1,c2) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                        !                         (-T_func(k2) * g_func(k) + h_func(k2))
                        coeff_fdamp(k2) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                                 (-T_func(k2) * g_func(k) + h_func(k2))
                        dg = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * (rjs_H(k)/sigma_ij)**2 * &
                             (-rjs_H(k)/(3.d0 * sigma_ij**3))
                        dh = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * xyz_H(c1,k)*xyz_H(c2,k) / &
                             (sigma_ij**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij)**2)
                        !coeff_der(i,j,c1,c2) = (1.d0-f_damp(k)) * (-T_func(k2)*dg+dh)
                        coeff_der(k2) = (1.d0-f_damp(k)) * (-T_func(k2)*dg+dh)
                      end do
                    end do
                  end if
                end do
              end do

            !k2 = 0
            !k3 = 0
            
            !do i = 1, n_sites
            !  r_vdw_i = r0_ii(k3+1)
            !  do a2 = 1, n_neigh(i)
              k3 = 0
              do i = 1, n_sites
                n_tot = sum(n_neigh(1:i))-n_neigh(i)
                if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == a) ) then
                  p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),a,1)
            !    k2 = k2+1
            !    a = neighbors_list(k2)
                  do j2 = 1, n_neigh(i)
                    if (rjs(k3+j2) < rcut) then
                      j = neighbors_list(k3+j2)
                      if ( i == j ) then
                        do c1 = 1, 3
                          !do c2 = 1, 3
                          b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) + 1.d0/(alpha_i(i) * hirshfeld_v(i)) * &
                                                   hirshfeld_v_cart_der_H(c3,n_tot+p) * a_SCS(3*(i-1)+c1,:)
                          !end do 
                        end do
                      end if
                      k4 = 9*(k3+j2-1)
                      do c1 = 1, 3
                        do c2 = 1, 3
                          k4 = k4+1
                          !dT_SR_v(3*(i-1)+c1,3*(j-1)+c2) = dT_SR_v(3*(i-1)+c1,3*(j-1)+c2) + &
                          !  (coeff_der(i,j,c1,c2) * sigma_i(i)**2/hirshfeld_v(i) + &
                          !  coeff_fdamp(i,j,c1,c2) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,n_tot+p) ! n_tot+p = k2
                          !dT_SR_v(3*(j-1)+c1,3*(i-1)+c2) = dT_SR_v(3*(j-1)+c1,3*(i-1)+c2) + &
                          !  (coeff_der(j,i,c1,c2) * sigma_i(i)**2/hirshfeld_v(i) + &
                          !  coeff_fdamp(j,i,c1,c2) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,n_tot+p) ! n_tot+p = k2
                          b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) - &
                            ((coeff_der(k4) * sigma_i(i)**2/hirshfeld_v(i) + &
                            coeff_fdamp(k4) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,n_tot+p)) * &
                            a_SCS(3*(j-1)+c2,:)
                          ! THE FOLLOWING ONLY WORKS ON THE ASSUMPTION THAT COEFF_DER AND COEFF_FDAMP ARE SYMMETRIC!
                          ! THAT IS, B HAS TO BE SYMMETRIC!
                          b_der(3*(j-1)+c1,:) = b_der(3*(j-1)+c1,:) - &
                            ((coeff_der(k4) * sigma_i(i)**2/hirshfeld_v(i) + &
                            coeff_fdamp(k4) * r_vdw_i/hirshfeld_v(i)) * hirshfeld_v_cart_der_H(c3,n_tot+p)) * &
                            a_SCS(3*(i-1)+c2,:)
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
                end if
                k3 = k3 + n_neigh(i)
              end do
            !  k3 = k3 + n_neigh(i)

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

            !k2 = 0
            !do i = 1, n_sites
            !  do a2 = 1, n_neigh(i)
            !    k2 = k2+1
            !    a = neighbors_list(k2)
       !       do i = 1, n_sites
       !         n_tot = sum(n_neigh(1:i))-n_neigh(i)
       !         if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == a) ) then
       !           p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),a,1)
! NOTE: There is probably a better way to do this:
                  !dT_SR_v(3*(i-1)+1:3*(i-1)+3,:) = dT_SR_v(3*(i-1)+1:3*(i-1)+3,:) + &
                  !  1.d0/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,n_tot+p) * &
                  !  (B_mat(3*(i-1)+1:3*(i-1)+3,:)-I_mat(3*(i-1)+1:3*(i-1)+3,:)) / alpha_i(i)
       !           do j = 1, n_sites
       !             do c1 = 1, 3
       !               do c2 = 1, 3
       !                 b_der(3*(i-1)+c1,:) = b_der(3*(i-1)+c1,:) - alpha_i(i) * (1.d0/hirshfeld_v(i) * &
       !                   hirshfeld_v_cart_der_H(c3,n_tot+p) * (B_mat(3*(i-1)+c1,3*(j-1)+c2) &
       !                   -I_mat(3*(i-1)+c1,3*(j-1)+c2)) / alpha_i(i)) * a_SCS(3*(j-1)+c2,:)
       !               end do
       !             end do
       !           end do
       !         end if
       !       end do

              !dT_SR = dT_SR + dT_SR_v
            end if

            !do i = 1, n_sites
            !  dT_SR(3*(i-1)+1:3*(i-1)+3,:) = alpha_i(i) * dT_SR(3*(i-1)+1:3*(i-1)+3,:) 
            !end do

            da_vec = 0.d0
            !k2 = 0
        !    do i = 1, n_sites
        !      n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !      if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == a) ) then
        !        p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),a,1)
        !        do c1 = 1, 3
        !          da_vec(3*(i-1)+c1,c1) = 1.d0/hirshfeld_v(i) * alpha_i(i) * &
        !            hirshfeld_v_cart_der_H(c3,n_tot+p)
        !        end do
        !      end if
        !    end do
            
      ! What happens here: regularized equation
      ! (B^T * B + reg_param * I) a_SCS = (B^T + reg_param*I) a_vec
      ! is differentiated:
      ! (B^T * B + reg_param * I) a_SCS' = (B')^T a_vec + (B^T + reg_param*I) a_vec' - ((B')^T * B) a_SCS - (B^T * B') a_SCS
      !                                    |---vect1--|   |---------vect2----------|   |-----vect3------|   |----vect4-----|
      ! a_SCS' is solved with dsysv.

            if( do_timing) then
              call cpu_time(time2)
              write(*,*) "Gradients: timing for everything else:", time2-time1
              call cpu_time(time1)
            end if
            
      if ( psblas ) then
      
      !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR, 3*n_sites, &
      !           a_SCS, 3*n_sites, 0.d0, vect_temp, 3*n_sites)
      !da_SCS = da_vec - vect_temp
      !b_der = b_der + da_vec


      !write(*,*) "da_SCS"
      !do p = 1, 6
      !  write(*,*) da_SCS(p,:)
      !end do
      !write(*,*) "b_der"
      !do p = 1, 6
      !  write(*,*) b_der(p,:)
      !end do

      !allocate( val_xv(1:3*n_sites,1:3) )
      !allocate( val_bv(1:3*n_sites,1:3) )
      !allocate( myidx(1:3*n_sites) )
      !val_xv = 0.d0
      !val_bv = 0.d0
      !k2 = 0
      !do p = 1, 3*n_sites
      !    k2 = k2+1
      !    myidx(k2) = k2
      !    !val_xv(p,:) = da_SCS(p,:)
      !    !val_bv(p,:) = da_SCS(p,:)
      !    val_xv(p,:) = b_der(p,:)
      !    val_bv(p,:) = b_der(p,:)
      !end do

      call cpu_time(time1)
      !call psb_init(icontxt)
      !do c1 = 1, 3
      !  call psb_cdall(icontxt, desc_a, info_psb, vl=myidx)
      !  !write(*,*) "cdall", info_psb
      !  call psb_spall(A_sp, desc_a, info_psb, nnz=nnz)
      !  !write(*,*) "spall", info_psb
      !  call psb_geall(x_vec,desc_a,info_psb)
      !  !write(*,*) "geall x", info_psb
      !  call psb_geall(b_vec,desc_a,info_psb)
      !  !write(*,*) "geall b", info_psb
      !  call psb_spins(nnz, ia(1:nnz), ja(1:nnz), val(1:nnz), A_sp, desc_a, info_psb)
      !  !write(*,*) "spins", info_psb
      !  call psb_geins(3*n_sites, myidx, val_xv(:,c1), x_vec, desc_a, info_psb)
      !  !write(*,*) "x_vec", info_psb
      !  call psb_geins(3*n_sites, myidx, val_bv(:,c1), b_vec, desc_a, info_psb)
      !  !write(*,*) "b_vec", info_psb
      !  call psb_cdasb(desc_a, info_psb)
      !  !write(*,*) "cdasb", info_psb
      !  call psb_spasb(A_sp, desc_a, info_psb)
      !  !write(*,*) "spasb", info_psb
      !  call psb_geasb(x_vec, desc_a, info_psb)
      !  !write(*,*) "geasb x", info_psb
      !  call psb_geasb(b_vec, desc_a, info_psb)
      !  !write(*,*) "geasb b", info_psb
      !  ptype="DIAG"
      !  ! NOTE: Everything works fine until preconditioner has to be set. Then the compilation fails.
      !  call prec%init(icontxt, ptype, info_psb)
      !  !write(*,*) "prec init", info_psb
      !  call prec%build(A_sp, desc_a, info_psb)
      !  !write(*,*) "prec build", info_psb
      !  call psb_krylov("BICGSTAB", A_sp, prec, b_vec, x_vec, 0.000001d0, desc_a, info_psb)
      !  !write(*,*) "krylov", info_psb
      !  val_xv(:,c1) = x_vec%get_vect()
      !end do
      !call psb_exit(icontxt)
      !write(*,*) "val_xv"
      !dalpha_full(:,a,c3,om) = 0.d0
      !do p = 1, n_sites
      !  dalpha_full(p,a,c3,om) = 1.d0/3.d0 * (val_xv(3*(p-1)+1,1) + val_xv(3*(p-1)+2,2) + val_xv(3*(p-1)+3,3))
      !end do
      !deallocate( val_xv, val_bv, myidx )
      !call psb_exit(icontxt)
      !deallocate( ia, ja, val )
      !call cpu_time(time2)
      !write(*,*) "Timing for PSBLAS", time2-time1
      
      else ! psblas

            if ( regularization ) then
              do i = 1, 3*n_sites
                do c1 = 1, 3
                  vect1 = 0.d0
                  !vect1(i,c1) = dot_product(dT_SR(c1::3,i), a_vec(c1::3,c1))
                  vect2(i,c1) = dot_product(B_reg(c1::3,i), da_vec(c1::3,c1))
                end do
              end do
              call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_mat, 3*n_sites, &
                a_SCS, 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              !call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR, 3*n_sites, &
              !  vect_temp, 3*n_sites, 0.d0, vect3, 3*n_sites)
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR, 3*n_sites, &
              !  a_SCS, 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_mat, 3*n_sites, &
                vect_temp, 3*n_sites, 0.d0, vect4, 3*n_sites)
              da_SCS = vect1 + vect2 - vect3 - vect4
              call dsytrs('U', 3*n_sites, 3, BTB_reg, 3*n_sites, ipiv, &
                da_SCS, 3*n_sites, info)
            else
              !polyfit = (/ -0.01252063d0,  0.16690202d0, -0.90556109d0,  2.52934965d0, -3.713878d0, 2.3867889d0, &
              !              0.03266611d0, -0.19889042d0, -0.97601187d0,  1.19685207d0, -1.01171277d0,  0.99889649d0 /)
              !polyfit = (/ -0.011111d0,    0.10555365d0, -0.42261792d0,  0.93352476d0, -1.27295025d0,  1.20591827d0, &
              !   -1.00960849d0,  0.96270989d0, -0.99307669d0,  1.00191393d0, -1.00020537d0,  0.99999289d0 /)
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR, 3*n_sites, &
              !  a_SCS, 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              !da_SCS = da_vec - vect_temp
              !call dgetrs('n', 3*n_sites, 3, BTB_reg, 3*n_sites, ipiv, &
              !  da_SCS, 3*n_sites, info)
              !b_der = da_vec + b_der
              ! DERIVATIVE TEST HACK
              !call dgetrs('n', 3*n_sites, 3, BTB_reg, 3*n_sites, ipiv, &
              !  b_der, 3*n_sites, info)
              !da_SCS = polyfit(n_degree) * b_der
              !vect_temp = b_der
              !do p = 1, n_degree
              !  call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_mat-I_mat, 3*n_sites, &
              !    vect_temp, 3*n_sites, 0.d0, vect1, 3*n_sites)
              !    da_SCS = da_SCS + polyfit(n_degree+1-p)*vect1
              !    vect_temp = vect1
              !end do
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, -B_mat+I_mat, 3*n_sites, &
              !  b_der, 3*n_sites, 0.d0, vect1, 3*n_sites)
              !da_SCS = da_SCS - 1.13904254d0 * vect1
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, -B_mat+I_mat, 3*n_sites, &
              !  vect1, 3*n_sites, 0.d0, vect2, 3*n_sites)
              !da_SCS = da_SCS + 0.92858092d0 * vect2
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, -B_mat+I_mat, 3*n_sites, &
              !  vect2, 3*n_sites, 0.d0, vect3, 3*n_sites)
              !da_SCS = da_SCS - 0.37501278d0 * vect3
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, -B_mat+I_mat, 3*n_sites, &
              !  vect3, 3*n_sites, 0.d0, vect4, 3*n_sites)
              !da_SCS = da_SCS + 0.05495458d0 * vect4
              !call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, -B_mat+I_mat, 3*n_sites, &
              !  vect4, 3*n_sites, 0.d0, vect_temp, 3*n_sites)
              !da_SCS = da_SCS + vect_temp/120.d0
              !write(*,*) "b_der"
              !do p = 1, 3*n_sites
              !  write(*,*) b_der(p,:)
              !end do
              !write(*,*) "b_der", b_der
              call dsytrs('U', 3*n_sites, 3, B_mat, 3*n_sites, ipiv, &
                b_der, 3*n_sites, info)
              !write(*,*) "b_der", b_der
              !b_der = da_SCS
              !write(*,*) "just for fun"
              ! END
            end if

            !n_tot = sum(n_neigh(1:a))-n_neigh(a)
            !do i2 = 1, n_neigh(a)
            !do i = 1, n_sites
                ! DERIVATIVE TEST: THIS NEEDS TO BE CHANGED INTO FORM dalpha(k_a,c3) TO WORK
                !dalpha_full(i,a,c3,om) = 1.d0/3.d0 * (b_der(3*(i-1)+1,1) + &
                ! b_der(3*(i-1)+2,2) + b_der(3*(i-1)+3,3))
            !end do

            if( do_timing) then
              call cpu_time(time2)
              write(*,*) "Gradients: timing for solving alpha_SCS_grad:", time2-time1
              call cpu_time(time1)
            end if
            
            end if ! psblas

          end do ! a loop
          
        end do ! c3 loop

      end if ! do_derivatives

      !end if ! EIGENVALUE TEST if (.false.)

    end do ! om loop

    do i = 1, n_sites
      alpha_SCS0(i,3) = 2.d0*omega_i(1)/sqrt(alpha_SCS0(i,1)/alpha_SCS0(i,2)-1.d0)
    end do

            ! Cut-off test: This is a hack to test the validity of the approximation;
            ! Implement this in the loops above
            !do i = 1, n_sites
            !  n_tot = sum(n_neigh(1:i))-n_neigh(i)
            !  do a = 1, n_sites
            !    if (.not.( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == a) )) then   
            !      dalpha_full(i,a,:,:) = 0.d0
            !    end if
            !  end do
            !end do

    !EIGENVALUE TEST
    !if ( .false. ) then

    end if ! local

    !write(*,*) "alpha_SCS:" !,  alpha_SCS0(1,1)
    !do p = 1, n_sites
    !  write(*,*) p, alpha_SCS0(p,:)
    !end do

    if ( .false. ) then
      write(*,*) "atom, alpha_SCS, C6_SCS:"
      do p = 1, n_sites
        call integrate("trapezoidal", omegas, alpha_SCS0(p,:)**2, omegas(1), omegas(n_freq), c6_nsites(p))
        c6_nsites(p) = 3.d0/pi*c6_nsites(p)
        write(*,*) p, alpha_SCS0(p,1) * Bohr**3, c6_nsites(p)
      end do

      k = 0
      do i = 1, n_sites
        do j = 1,  n_neigh(i)
          k = k+1
          i2 = neighbors_list(k)
          c6_scs(k) = c6_nsites(i2) * (Hartree*Bohr**6)
          r0_scs(k) = (alpha_SCS0(i2,1)/neighbor_alpha0(k))**(1.d0/3.d0) * r0_ii(k) * Bohr
          alpha0_scs(k) = alpha_SCS0(i2,1) * Bohr**3
        end do
      end do
    end if

    !if ( do_derivatives ) then
      !write(*,*) "dalpha_SCS w.r.t. atom 1 in x direction:"
    !  write(*,*) "dalpha_SCS of atom 1 w.r.t. its neighbors:"
      !do i2 = 1, n_neigh(1)
      !  write(*,*) neighbors_list(i2), dalpha_full(i2,1)
      !end do
      !k = 0
    !  do i2 = 1, n_sites
      !write(*,*) "gradients of atom x w.r.t. atom 1"
        !do a2 = 1, n_neigh(i2)
        !  k = k+1
        !  if ( neighbors_list(k) == 1 .or. neighbors_list(k) == 2) then
    !        write(*,*) i2, dalpha_full(i2,1)
        !  end if
        !end do
    !  end do
    !end if

    write(*,*) "Full MBD energy", sum(energies)

    !deallocate( ia, ja, val )

    !end if ! EIGENVALUE TEST if ( .false. )

!   Clean up

    if ( .not. local ) then
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, B_mat, xyz_H, rjs_H, &
                BTB_reg, B_reg, I_mat, a_vec, a_SCS, c6_nsites )
    if ( do_derivatives ) then
      deallocate( dT, f_damp_der, g_func_der, h_func_der, &
                  coeff_der, coeff_fdamp, hirshfeld_v_cart_der_H, da_vec, &
                  vect1, vect2, vect3, vect4, vect_temp, da_SCS, b_der ) !dT_SR, dT_SR_v
    
    end if
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
                           xyz_H(:,:), rjs_H(:), f_damp_SCS(:), T_func(:), AT(:,:), AT_n(:,:,:), &
                           AT_n_f(:,:,:), energy_series(:,:), integrand(:,:), omegas(:), f_damp_der(:), &
                           f_damp_der_SCS(:), dT_LR(:,:), dT(:), &
                           G_mat(:,:), force_integrand(:,:,:), force_series(:,:), &
                           VL(:,:), total_energy_series(:,:), total_integrand(:,:), &
                           xyz_2b(:,:), rjs_2b(:), r0_ii_2b(:), neighbor_alpha0_2b(:), &
                           T_func_2b(:), t_vec(:,:), r0_ii_SCS_2b(:), f_damp_SCS_2b(:), g_vec(:,:)
    integer :: n_order, n_freq, n_sites, n_pairs, n_species, n_sites0, n_sub_sites, n_sub_pairs, &
               n_2b_sites, n_2b_pairs, n_tot_sites
    integer :: k, k2, k3, i, i2, j, j2, j3, c1, c2, c3, a, a2, om, p, q, r, n_tot, n_tot2, s
    real*8 :: Bohr, Hartree, pi, r_vdw_i, r_vdw_j, E_MBD, integral, omega, R_vdW_SCS_ij, S_vdW_ij, &
              dS_vdW_ij, time1, time2, rcut_mbd, C6_2b, E_TS
    logical :: series_average = .false., do_timing = .false., local = .true.
    logical, allocatable :: in_cutoff(:), in_2b_cutoff(:) 
    integer, allocatable :: p_to_i(:), i_to_p(:), sub_neighbors_list(:), n_sub_neigh(:), sub_2b_list(:), &
                            r_to_i(:), i_to_r(:)


    rcut_mbd = 4.5d0
    write(*,*) "Full cutoff", rcut
    write(*,*) "MBD cutoff", rcut_mbd
!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

    n_order = 2

    !n_freq = size(alpha_SCS0, 2)
    n_freq = 12
    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

    call cpu_time(time1)

    allocate( omegas(1:n_freq) )
    allocate( integrand(1:n_sites,1:n_freq) )
    if ( do_derivatives ) then
      allocate( force_integrand(1:n_sites,1:3,1:n_freq) )
    end if
    
    omega = 0.d0
    do i = 1, n_freq
      omegas(i) = omega
      omega = omega + 0.4d0
    end do

    integrand = 0.d0
    if ( do_derivatives ) then
      force_integrand = 0.d0
    end if

    if ( local ) then

      if ( do_derivatives ) then
        allocate( total_integrand(1:n_sites,1:n_freq) )
        total_integrand = 0.d0
      end if
    
      allocate( in_cutoff(1:n_sites) )
      allocate( in_2b_cutoff(1:n_sites) )
      allocate( p_to_i(1:n_sites) )
      allocate( i_to_p(1:n_sites) )
      allocate( r_to_i(1:n_sites) )
      allocate( i_to_r(1:n_sites) )

      k = 0
      do i = 1, n_sites
      
        p_to_i = 0
        i_to_p = 0
        r_to_i = 0
        i_to_r = 0
        in_cutoff = .false.
        in_2b_cutoff = .false.
        k = k+1
        p = 1
        r = 0
        p_to_i(p) = i
        i_to_p(i) = p
        in_cutoff(i) = .true.
        do j2 = 2, n_neigh(i)
          k = k+1
          j = neighbors_list(k)
          if (rjs(k) < rcut_mbd) then
            in_cutoff(j) = .true.
            p = p+1
            p_to_i(p) = j
            i_to_p(j) = p
          end if
          if (rjs(k) .ge. rcut_mbd .and. rjs(k) < rcut) then
            in_2b_cutoff(j) = .true.
            r = r+1
            r_to_i(r) = j
            i_to_r(j) = r
          end if
        end do
        n_sub_sites = p
        n_2b_sites = r
        !write(*,*) "in_cutoff", in_cutoff
        n_sub_pairs = 0
        !n_2b_pairs = 0
        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        allocate( n_sub_neigh(1:n_sub_sites) )
        n_sub_neigh = 0
        p = 0
        !r = 0
        do j2 = 1, n_neigh(i)
          j = neighbors_list(n_tot+j2)
          if ( in_cutoff(j) ) then
            p = p + 1
            n_sub_pairs = n_sub_pairs + 1
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            n_tot2 = sum(n_neigh(1:j))-n_neigh(j)
            do j3 = 2, n_neigh(j)
              if ( rjs(n_tot2+j3) < rcut_mbd ) then
                if ( in_cutoff(neighbors_list(n_tot2+j3)) ) then
                  n_sub_neigh(p) = n_sub_neigh(p) + 1
                  n_sub_pairs = n_sub_pairs + 1
                end if
              end if
            end do
          end if
          !if ( in_2b_cutoff(j) ) then
          !  r = r + 1
          !  n_2b_pairs = n_2b_pairs + 1
          !end if
        end do

        write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_2b_pairs", n_2b_pairs
        write(*,*) "n_2b_sites", n_2b_sites        

        !write(*,*) "p_to_i", p_to_i
        !write(*,*) "i_to_p", i_to_p
        !write(*,*) "n_sub_pairs", n_sub_pairs
        !write(*,*) "n_sub_sites", n_sub_sites
        !write(*,*) "n_sub_neigh", n_sub_neigh
        
        allocate( sub_neighbors_list(1:n_sub_pairs) )
        allocate( xyz_H(1:3,1:n_sub_pairs) )
        allocate( rjs_H(1:n_sub_pairs) )
        allocate( r0_ii(1:n_sub_pairs) )
        allocate( neighbor_alpha0(1:n_sub_pairs) )
        allocate( T_func(1:9*n_sub_pairs) )
        allocate( T_LR(1:3*n_sub_sites,1:3*n_sub_sites) )
        allocate( r0_ii_SCS(1:n_sub_pairs) )
        allocate( f_damp_SCS(1:n_sub_pairs) )
        allocate( AT(1:3*n_sub_sites,1:3*n_sub_sites) )
        !allocate( AT_n(1:3*n_sub_sites,1:3*n_sub_sites,1:n_order-1) )
        allocate( AT_n(1:3*n_sub_sites,1:3,1:n_order-1) )
        allocate( energy_series(1:3*n_sub_sites,1:3) )
        if ( do_derivatives ) then
          allocate( AT_n_f(1:3*n_sub_sites,1:3*n_sub_sites,1:n_order-1) )
          allocate( dT(1:9*n_sub_pairs) )
          allocate( f_damp_der(1:n_sub_pairs) )
          allocate( f_damp_der_SCS(1:n_sub_pairs) )
          allocate( dT_LR(1:3*n_sub_sites,1:3*n_sub_sites) )
          allocate( G_mat(1:3*n_sub_sites,1:3*n_sub_sites) )
          !allocate( force_series(1:3*n_sub_sites,1:3*n_sub_sites) )
          allocate( force_series(1:3*n_sub_sites,1:3*n_sub_sites) )
          allocate( VL(1:3*n_sub_sites,1:3*n_sub_sites) )
          allocate( total_energy_series(1:3*n_sub_sites,1:3*n_sub_sites) )
        end if
        
        n_tot_sites = n_sub_sites+n_2b_sites

        allocate( sub_2b_list(1:n_2b_sites) )
        allocate( xyz_2b(1:3,1:n_2b_sites) )
        allocate( rjs_2b(1:n_2b_sites) )
        allocate( r0_ii_2b(1:n_2b_sites) )
        allocate( neighbor_alpha0_2b(1:n_2b_sites) )
        allocate( T_func_2b(1:9*n_2b_sites) )
        allocate( t_vec(1:3,1:3*n_2b_sites) )
        allocate( r0_ii_SCS_2b(1:n_2b_sites) )
        allocate( f_damp_SCS_2b(1:n_2b_sites) )
        allocate( g_vec(1:3,1:3*n_2b_sites) )


        k2 = 0
        T_func = 0.d0
        T_LR = 0.d0
        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        !n_sub_neigh = 0
        do j2 = 1, n_neigh(i)
          i2 = neighbors_list(n_tot+j2)
          if ( in_cutoff(i2) ) then 
            k2 = k2+1
            s = neighbor_species(n_tot+j2)
            sub_neighbors_list(k2) = i2
            !n_sub_neigh(j2) = n_sub_neigh(i2) + 1
            r0_ii(k2) = r0_ref(s) / Bohr
            neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3
            xyz_H(:,k2) = xyz(:,n_tot+j2)/Bohr
            rjs_H(k2) = rjs(n_tot+j2)/Bohr
            r0_ii_SCS(k2) = r0_ii(k2) * (alpha_SCS0(i2,1)/neighbor_alpha0(k2))**(1.d0/3.d0)
            r_vdw_i = r0_ii_SCS(k2)
            n_tot2 = sum(n_neigh(1:i2))-n_neigh(i2)
            do j3 = 2, n_neigh(i2)
              if ( rjs(n_tot2+j3) < rcut_mbd ) then
                if ( in_cutoff(neighbors_list(n_tot2+j3)) ) then
                  !n_sub_neigh(j2) = n_sub_neigh(j2) + 1
                  k2 = k2+1    
                  s = neighbor_species(n_tot2+j3)
                  j = neighbors_list(n_tot2+j3)
                  sub_neighbors_list(k2) = j                        
                  r0_ii(k2) = r0_ref(s) / Bohr
                  neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3
                  xyz_H(:,k2) = xyz(:,n_tot2+j3)/Bohr
                  rjs_H(k2) = rjs(n_tot2+j3)/Bohr
                  r0_ii_SCS(k2) = r0_ii(k2) * (alpha_SCS0(j,1)/neighbor_alpha0(k2))**(1.d0/3.d0)
                  r_vdw_j = r0_ii_SCS(k2)
                  f_damp_SCS(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                  k3 = 9*(k2-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k3 = k3 + 1
                      if (c1 == c2) then
                        T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c1,k2) - rjs_H(k2)**2)/rjs_H(k2)**5
                      else
                        T_func(k3) = (3*xyz_H(c1,k2) * xyz_H(c2,k2))/rjs_H(k2)**5
                      end if
                      p = i_to_p(i2)
                      q = i_to_p(j)
                      !write(*,*) "f_damp_SCS", f_damp_SCS(k2)
                      !write(*,*) "T_func", T_func(k3)
                      T_LR(3*(p-1)+c1,3*(q-1)+c2) = f_damp_SCS(k2) * T_func(k3)
                    end do
                  end do
                end if
              end if
            end do
          end if
        end do

        E_TS = 0.d0
        k2 = 0
        T_func_2b = 0.d0
        t_vec = 0.d0
        !t_vec(1:3,1:3*n_sub_sites) = T_LR(1:3,1:3*n_sub_sites)
        n_tot = sum(n_neigh(1:i))-n_neigh(i)
        s = neighbor_species(n_tot+1)
        r_vdw_i = r0_ref(s) / Bohr * (alpha_SCS0(i,1)/(alpha0_ref(s)/Bohr**3))**(1.d0/3.d0)
        do j2 = 1, n_neigh(i)
          i2 = neighbors_list(n_tot+j2)
          if ( in_2b_cutoff(i2) ) then
            k2 = k2+1
            s = neighbor_species(n_tot+j2)
            sub_2b_list(k2) = i2
            r0_ii_2b(k2) = r0_ref(s) / Bohr
            neighbor_alpha0_2b(k2) = alpha0_ref(s) / Bohr**3
            xyz_2b(:,k2) = xyz(:,n_tot+j2)/Bohr
            rjs_2b(k2) = rjs(n_tot+j2)/Bohr
            r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (alpha_SCS0(i2,1)/neighbor_alpha0_2b(k2))**(1.d0/3.d0)
            r_vdw_j = r0_ii_SCS_2b(k2)
            ! Here sR=0.97 for TS-SCS!
            f_damp_SCS_2b(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
            C6_2b = 3.d0/2.d0 * alpha_SCS0(i,1) * alpha_SCS0(i2,1) * (alpha_SCS0(i,3) * alpha_SCS0(i2,3)) / &
                    (alpha_SCS0(i,3) + alpha_SCS0(i2,3))
            E_TS = E_TS - C6_2b/rjs_2b(k2)**6 * f_damp_SCS_2b(k2)
            !MBD 2b approach below
            !k3 = 9*(k2-1)
            !do c1 = 1, 3
            !  do c2 = 1, 3
            !    k3 = k3 + 1
            !    if (c1 == c2) then
            !      T_func_2b(k3) = (3*xyz_2b(c1,k2) * xyz_2b(c1,k2) - rjs_2b(k2)**2)/rjs_2b(k2)**5
            !    else
            !      T_func_2b(k3) = (3*xyz_2b(c1,k2) * xyz_2b(c2,k2))/rjs_2b(k2)**5
            !    end if
            !    q = i_to_r(i2)
            !    !q = i_to_r(j)
            !    t_vec(c1,3*(q-1)+c2) = f_damp_SCS_2b(k2) * T_func_2b(k3)
            !  end do
            !end do
          end if
        end do
        E_TS = 1.d0/2.d0 * E_TS

        !write(*,*) "t_vec"
        !do p = 1, 3*n_2b_sites
        !  write(*,*) p, t_vec(:,p)
        !end do



        !write(*,*) "T_LR"
        !do p = 1, n_sub_sites
        !  write(*,*) T_LR(p,:)
        !end do
            
        do om = 1, n_freq
    
          AT = 0.d0
          do p = 1, n_sub_sites
            i2 = p_to_i(p)
            AT(3*(p-1)+1:3*(p-1)+3,:) = alpha_SCS0(i2,1)/(1.d0 + (omegas(om)/alpha_SCS0(i2,3))**2) * T_LR(3*(p-1)+1:3*(p-1)+3,:)
          end do
          g_vec = 0.d0
          do p = 1, n_2b_sites
            i2 = r_to_i(p)
            g_vec(:,3*(p-1)+1:3*(p-1)+3) = alpha_SCS0(i2,1)/(1.d0 + (omegas(om)/alpha_SCS0(i2,3))**2) * t_vec(:,3*(p-1)+1:3*(p-1)+3)
          end do

          !if ( i == 1 .and. om == 1 ) then
          !  write(*,*) "3*n_sub_sites", 3*n_sub_sites
          !  write(*,*) "AT"
          !  do p = 1, 3*n_sub_sites
          !    write(*,*) AT(p,:)
          !  end do
          !end if
 
          AT_n = 0.d0
          energy_series = 0.d0
          if ( do_derivatives ) then
            AT_n_f = 0.d0
            force_series = 0.d0
            total_energy_series = 0.d0
          end if

          if ( n_order > 1 ) then
            do k2 = 1, n_order-1
              ! Precalculate the full AT_n for forces:
              if ( k2 == 1 ) then
                AT_n(:,:,k2) = AT(:,1:3)
                if ( do_derivatives ) then
                  AT_n_f(:,:,k2) = AT
                end if
              else
                call dgemm('n', 'n', 3*n_sub_sites, 3, 3*n_sub_sites, 1.d0, AT, &
                           3*n_sub_sites, AT_n(:,:,k2-1), 3*n_sub_sites, 0.d0, AT_n(:,:,k2), &
                           3*n_sub_sites)
                if ( do_derivatives ) then
                  call dgemm('n', 'n', 3*n_sub_sites, 3*n_sub_sites, 3*n_sub_sites, 1.d0, AT, &
                             3*n_sub_sites, AT_n_f(:,:,k2-1), 3*n_sub_sites, 0.d0, AT_n_f(:,:,k2), &
                             3*n_sub_sites)
                end if
              end if
              if (series_average .and. k2 == n_order-1) then
                energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,:,k2)/2.d0
                if ( do_derivatives ) then
                  force_series = force_series + AT_n_f(:,:,k2)/2.d0
                  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2)/2.d0
                end if
              else
                energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,1:3,k2)
                if ( do_derivatives ) then
                  force_series = force_series + AT_n_f(:,:,k2)
                  total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2)
                end if
              end if
            end do
            do c1 = 1, 3
              integrand(i,om) = integrand(i,om) + alpha_SCS0(i,1)/(1.d0 + (omegas(om)/alpha_SCS0(i2,3))**2) &
                                * dot_product(T_LR(c1,:),energy_series(:,c1))
              !write(*,*) "T_LR", T_LR(c1,:)
              !write(*,*) "energy_series", energy_series(:,c1)
            end do
            if ( do_derivatives ) then
              do p = 1, n_sub_sites
                i2 = p_to_i(p)
                do c1 = 1, 3
                  total_integrand(i,om) = total_integrand(i,om) + alpha_SCS0(i,1)/(1.d0 + &
                       (omegas(om)/alpha_SCS0(i2,3))**2) * dot_product(T_LR(3*(p-1)+c1,:), &
                       total_energy_series(:,3*(p-1)+c1))
                  end do
              end do
            end if
            ! Two-body stuff
            !do c1 = 1, 3
            !  integrand(i,om) = integrand(i,om) - alpha_SCS0(i,1)/(1.d0 + (omegas(om)/alpha_SCS0(i,3))**2) * &
            !                    dot_product(t_vec(c1,:),g_vec(c1,:))/2.d0
            !end do


          else
            write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
            write(*,*) "and forces will be zero."
          end if

          
          
          if ( do_derivatives ) then
      
          do c3 = 1, 3
          
            dT = 0.d0
            k2 = 0
            do p = 1, n_sub_sites
              k2 = k2+1
              i2 = sub_neighbors_list(k2)
              if ( n_sub_neigh(p) > 1 ) then
                do j3 = 2, n_sub_neigh(p)
                  k2 = k2+1
                  j = sub_neighbors_list(k2)    
                  k3 = 9*(k2-1)
                  do c1 = 1,3
                    do c2 = 1,3
                      k3 = k3+1
                      !write(*,*) "k3", k3
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
              end if
            end do
            
            !do a2 = 1, n_neigh(i)

              f_damp_der_SCS = 0.d0
              f_damp_der = 0.d0 ! This is cleared so we can recalculate it with SCS values
              dT_LR = 0.d0
              !a = neighbors_list(n_tot+a2)
              k2 = 0
              do p = 1, n_sub_sites
                k2 = k2+1
                r_vdw_i = r0_ii_SCS(k2)
                i2 = sub_neighbors_list(k2)
                if ( n_sub_neigh(p) > 1 ) then
                  do j3 = 2, n_sub_neigh(p)
                    k2 = k2+1
                    j = sub_neighbors_list(k2)
                    r_vdw_j = r0_ii_SCS(k2)
                    R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                    S_vdW_ij = sR*R_vdW_SCS_ij
                    dS_vdW_ij = sR/3.d0 * ( r_vdw_i/alpha_SCS0(i2,1) * alpha_SCS_grad(i2,i,c3,1) + &
                                              r_vdw_j/alpha_SCS0(j,1) * alpha_SCS_grad(j,i,c3,1) )
                    f_damp_der_SCS(k2) = -(d*rjs_H(k2))/S_vdW_ij**2 * f_damp_SCS(k2)**2 * &
                                           exp(-d*(rjs_H(k2)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                    ! DERIVATIVE TEST HACK:
                    !f_damp_der_SCS(k2) = 0.d0
                    ! END
                    k3 = 9*(k2-1)
                    !p = i_to_p(i2)
                    q = i_to_p(j)
                    !write(*,*) "p, j2", p, j2
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k3 = k3+1
                        dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                          T_func(k3) * f_damp_der_SCS(k2) 
                      end do
                    end do
                    if (i == i2 .or. i == j) then
                      f_damp_der(k2) = d/S_vdW_ij * f_damp_SCS(k2)**2 * &
                                      exp( -d*(rjs_H(k2)/S_vdW_ij - 1.d0) ) * xyz_H(c3,k2)/rjs_H(k2)
                      k3 = 9*(k2-1)
                      do c1 = 1, 3
                        do c2 = 1, 3
                          k3 = k3+1
                          if (i == i2) then
                            dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) - &
                                               T_func(k3) * f_damp_der(k2) - &
                                               dT(k3) * f_damp_SCS(k2)
                          else if (i == j) then
                            dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                               T_func(k3) * f_damp_der(k2) + &
                                               dT(k3) * f_damp_SCS(k2)
                          end if
                        end do
                      end do
                    end if
                  end do
                end if
              end do

              G_mat = 0.d0
!              force_series = 0.d0

              if ( n_order > 1 ) then

                do p = 1, n_sub_sites
                  i2 = p_to_i(p)
                  G_mat(3*(p-1)+1:3*(p-1)+3,:) = G_mat(3*(p-1)+1:3*(p-1)+3,:) + &
                    alpha_SCS0(i2,1)/(1.d0 + (omegas(om)/alpha_SCS0(i2,3))**2) * &
                    dT_LR(3*(p-1)+1:3*(p-1)+3,:) + &
                    alpha_SCS_grad(i2,i,c3,1)/(1.d0 + (omegas(om)/alpha_SCS0(i2,3))**2) * &
                    T_LR(3*(p-1)+1:3*(p-1)+3,:)
                end do
                !if ( i == 1 .and. om == 1 .and. c3 == 1 ) then
                !  write(*,*) "G_mat"
                !  do p = 1, 3*n_sub_sites
                !    write(*,*) G_mat(p,:)
                !  end do
                !end if
                do i2 = 1, 3
                !do i2 = 1, 3*n_sub_sites
                  force_integrand(i,c3,om) = force_integrand(i,c3,om) + &
                    dot_product(G_mat(i2,:), force_series(:,i2))
                end do
                
              end if
              
            !end do ! a loop
          
          end do ! c3 loop
        
        end if ! do derivatives

        end do ! om loop

        !write(*,*) "integrand", integrand(i,:)
        !write(*,*) "omegas", omegas
        integral = 0.d0
        call integrate("trapezoidal", omegas, integrand(i,:), omegas(1), omegas(n_freq), integral)
        E_MBD = integral/(2.d0*pi)
        energies(i) = (E_MBD+E_TS) * 27.211386245988
        write(*,*) "Local energy", i, energies(i)
        write(*,*) "MBD:", E_MBD * 27.211386245988
        write(*,*) "TS:", E_TS * 27.211386245988        

        if ( do_derivatives ) then
          do c3 = 1, 3
            integral = 0.d0
            call integrate("trapezoidal", omegas, force_integrand(i,c3,:), omegas(1), omegas(n_freq), integral)
            forces0(c3,i) = forces0(c3,i) + 1.d0/(2.d0*pi) * integral         
          end do
          forces0(:,i) = forces0(:,i) * 51.42208619083232
          write(*,*) "Force", i, forces0(:,i)
          integral = 0.d0
          call integrate("trapezoidal", omegas, total_integrand(i,:), omegas(1), omegas(n_freq), integral)
          E_MBD = integral/(2.d0*pi) * 27.211386245988
          write(*,*) "Total energy", i, E_MBD
        end if


        deallocate( n_sub_neigh, sub_neighbors_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, T_func, T_LR, &
                    r0_ii_SCS, f_damp_SCS, AT, AT_n, energy_series )
        deallocate( sub_2b_list, xyz_2b, rjs_2b, r0_ii_2b, neighbor_alpha0_2b, T_func_2b, t_vec, r0_ii_SCS_2b, &
                    f_damp_SCS_2b, g_vec )
        if ( do_derivatives ) then
          deallocate( AT_n_f, dT, f_damp_der, f_damp_der_SCS, dT_LR, G_mat, force_series, VL, total_energy_series )
        end if


      end do
      
      write(*,*) "MBD sum of local energies:", sum(energies)

      if ( do_derivatives ) then
        write(*,*) "MBD forces:"
        do i = 1, n_sites
          write(*,*) i, forces0(:,i)
        end do
      end if

      deallocate( in_cutoff, in_2b_cutoff, p_to_i, i_to_p, r_to_i, i_to_r )
      if ( do_derivatives ) then
        deallocate( total_integrand )
      end if

    else

      allocate( xyz_H(1:3,1:n_pairs) )
      allocate( rjs_H(1:n_pairs) )
      allocate( r0_ii(1:n_pairs) )
      allocate( neighbor_alpha0(1:n_pairs) )
      allocate( T_func(1:9*n_pairs) )
      allocate( T_LR(1:3*n_sites,1:3*n_sites) )
      allocate( r0_ii_SCS(1:n_pairs) )
      allocate( f_damp_SCS(1:n_pairs) )
      allocate( AT(1:3*n_sites,1:3*n_sites) )
      allocate( AT_n(1:3*n_sites,1:3*n_sites, 1:n_order-1) )
      allocate( energy_series(1:3*n_sites,1:3*n_sites) )
      if ( do_derivatives ) then
        allocate( dT(1:9*n_pairs) )
        allocate( f_damp_der(1:n_pairs) )
        allocate( f_damp_der_SCS(1:n_pairs) )
        allocate( dT_LR(1:3*n_sites,1:3*n_sites) )
        allocate( G_mat(1:3*n_sites,1:3*n_sites) )
        allocate( force_series(1:3*n_sites,1:3*n_sites) )
        allocate( VL(1:3*n_sites,1:3*n_sites) )
      end if

      do k = 1, n_pairs
        j = neighbor_species(k)
        r0_ii(k) = r0_ref(j) / Bohr
        neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
        xyz_H(:,k) = xyz(:,k)/Bohr
        rjs_H(k) = rjs(k)/Bohr
      end do

      T_func = 0.d0
      k = 0
      do i = 1, n_sites
        k = k+1
        do j2 = 2, n_neigh(i)
          k = k+1
          if( rjs(k) < rcut_mbd )then
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

      f_damp_SCS = 0.d0
      k = 0
      do i = 1, n_sites
        k = k+1
        r_vdw_i = r0_ii_SCS(k)
        do j2 = 2, n_neigh(i)
          k=k+1
          if (rjs(k) < rcut_mbd) then
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
          if (rjs(k) < rcut_mbd) then
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
    
      integrand = 0.d0
      if ( do_derivatives ) then
        force_integrand = 0.d0
      end if

      do om = 1, n_freq
    
        AT = 0.d0
        do i = 1, n_sites
          AT(3*(i-1)+1:3*(i-1)+3,:) = alpha_SCS0(i,om) * T_LR(3*(i-1)+1:3*(i-1)+3,:)
        end do

        AT_n = 0.d0
        energy_series = 0.d0
        if ( do_derivatives ) then
          force_series = 0.d0
        end if

        if ( n_order > 1 ) then
          do k2 = 1, n_order-1
            ! Precalculate the full AT_n for forces:
            if ( k2 == 1 ) then
              AT_n(:,:,k2) = AT
            else
              call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, AT, 3*n_sites, &
                         AT_n(:,:,k2-1), 3*n_sites, 0.d0, AT_n(:,:,k2), 3*n_sites)
            end if
            if (series_average .and. k2 == n_order-1) then
              energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,:,k2)/2.d0
              if ( do_derivatives ) then
                force_series = force_series + AT_n(:,:,k2)/2.d0
              end if
            else
              energy_series = energy_series - 1.d0/(k2+1)*AT_n(:,:,k2)
              if ( do_derivatives ) then
                force_series = force_series + AT_n(:,:,k2)
              end if
            end if
          end do
          do i = 1, n_sites
            do c1 = 1, 3
              integrand(i,om) = integrand(i,om) + alpha_SCS0(i,om) * dot_product(T_LR(3*(i-1)+c1,:), &
                             energy_series(:,3*(i-1)+c1))
            end do
          end do
        else
          write(*,*) "WARNING: Series expansion requires that vdw_mbd_order > 1 or the resulting energies"
          write(*,*) "and forces will be zero."
        end if

        if ( do_derivatives ) then
      
          do c3 = 1, 3

            dT = 0.d0
            k = 0
            do i = 1, n_sites
              k = k+1
              do j2 = 2, n_neigh(i)
                k = k+1
                if (rjs(k) < rcut_mbd) then
                  k2 = 9*(k-1)
                  do c1 = 1,3
                    do c2 = 1,3
                      k2 = k2+1
                      dT(k2) = (-15.d0 * xyz_H(c1,k) * xyz_H(c2,k) * xyz_H(c3,k))/rjs_H(k)**7
                      if (c1 == c2) then
                        dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c3,k)
                      end if
                      if (c2 == c3) then
                        dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c1,k)
                      end if
                      if (c1 == c3) then
                        dT(k2) = dT(k2) + 3.d0/rjs_H(k)**5 * xyz_H(c2,k)
                      end if
                    end do
                  end do
                end if
              end do
            end do

            do a = 1, n_sites

              f_damp_der_SCS = 0.d0
              f_damp_der = 0.d0 ! This is cleared so we can recalculate it with SCS values
              dT_LR = 0.d0

              k = 0
              do i = 1, n_sites
                k = k+1
                r_vdw_i = r0_ii_SCS(k)
                do j2 = 2, n_neigh(i)
                  k = k+1
                  if (rjs(k) < rcut_mbd) then
                    j = neighbors_list(k)
                    r_vdw_j = r0_ii_SCS(k)
                    R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                    S_vdW_ij = sR*R_vdW_SCS_ij
                    dS_vdW_ij = sR/3.d0 * ( r_vdw_i/alpha_SCS0(i,1) * alpha_SCS_grad(i,a,c3,1) + &
                                            r_vdw_j/alpha_SCS0(j,1) * alpha_SCS_grad(j,a,c3,1) )
                    f_damp_der_SCS(k) = -(d*rjs_H(k))/S_vdW_ij**2 * f_damp_SCS(k)**2 * &
                                         exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                    k2 = 9*(k-1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        k2 = k2+1
                        dT_LR(3*(i-1)+c1,3*(j-1)+c2) = dT_LR(3*(i-1)+c1,3*(j-1)+c2) + &
                                          T_func(k2) * f_damp_der_SCS(k) 
                      end do
                    end do
                    if (a == i .or. a == j) then
                      f_damp_der(k) = d/S_vdW_ij * f_damp_SCS(k)**2 * &
                                      exp( -d*(rjs_H(k)/S_vdW_ij - 1.d0) ) * xyz_H(c3,k)/rjs_H(k)
                      k2 = 9*(k-1)
                      do c1 = 1, 3
                        do c2 = 1, 3
                          k2 = k2+1
                          if (a == i) then
                            dT_LR(3*(i-1)+c1,3*(j-1)+c2) = dT_LR(3*(i-1)+c1,3*(j-1)+c2) - &
                                               T_func(k2) * f_damp_der(k) - &
                                               dT(k2) * f_damp_SCS(k)
                          else if (a == j) then
                            dT_LR(3*(i-1)+c1,3*(j-1)+c2) = dT_LR(3*(i-1)+c1,3*(j-1)+c2) + &
                                               T_func(k2) * f_damp_der(k) + &
                                               dT(k2) * f_damp_SCS(k)
                          end if
                        end do
                      end do
                    end if
                  end if
                end do
              end do

              G_mat = 0.d0
!              force_series = 0.d0

              if ( n_order > 1 ) then
!              do k2 = 1, n_order-1
!                if (series_average .and. k2 == n_order-1) then
!                  force_series = force_series + AT_n(:,:,k2)/2.d0
!                else
!                  force_series = force_series + AT_n(:,:,k2)
!                end if
!              end do
                ! call dgemm for G_mat
                do i = 1, n_sites
                  G_mat(3*(i-1)+1:3*(i-1)+3,:) = G_mat(3*(i-1)+1:3*(i-1)+3,:) + &
                    alpha_SCS0(i,om) * dT_LR(3*(i-1)+1:3*(i-1)+3,:) + &
                    alpha_SCS_grad(i,a,c3,om) * T_LR(3*(i-1)+1:3*(i-1)+3,:)
                end do
                !call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, force_series(:,:,k), 3*n_sites, &
                !           G_mat(:,:,a,c3,k), 3*n_sites, 0.d0, VL, 3*n_sites)
                do i = 1, 3*n_sites
                  !force_integrand(a,c3,k) = force_integrand(a,c3,k) + VL(i,i)
                  force_integrand(a,c3,om) = force_integrand(a,c3,om) + &
                    dot_product(G_mat(i,:), force_series(:,i))
                end do
              end if
            
            end do ! a loop
          
          end do ! c3 loop
        
        end if ! do derivatives
      
      end do ! om loop
      
      do i = 1, n_sites
        integral = 0.d0
        call integrate("trapezoidal", omegas, integrand(i,:), omegas(1), omegas(n_freq), integral)
        E_MBD = integral/(2.d0*pi)
        energies(i) = E_MBD * 27.211386245988
        write(*,*) i, energies(i)
      end do
      write(*,*) "MBD energy full system:", sum(energies)

      if ( do_derivatives ) then
        do a = 1, n_sites
          do c3 = 1, 3
            integral = 0.d0
            call integrate("trapezoidal", omegas, force_integrand(a,c3,:), omegas(1), omegas(n_freq), integral)
            forces0(c3,a) = 1.d0/(2.d0*pi) * integral
          end do
        end do
        forces0 = forces0 * 51.42208619083232

        write(*,*) "MBD forces:"
        do i = 1, n_sites
          write(*,*) i, forces0(:,i)
        end do

      end if

      deallocate( xyz_H, rjs_H, r0_ii, neighbor_alpha0, T_func, T_LR, r0_ii_SCS, &
                  f_damp_SCS, AT, AT_n, energy_series )
      if ( do_derivatives ) then
        deallocate( dT, dT_LR, f_damp_der, f_damp_der_SCS, force_series, G_mat, VL )
      end if
      
    end if ! local
    
    deallocate( omegas, integrand )
    if ( do_derivatives ) then
      deallocate( force_integrand )
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
