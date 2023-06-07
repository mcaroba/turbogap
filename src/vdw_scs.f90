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
  use nonneg_leastsq
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

    !write(*,*) "TS energy", sum(energies)

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
  subroutine get_scs_polarizabilities( n_neigh, neighbors_list, neighbor_species, &
                                       rcut, r_buffer, rjs, xyz, &
                                       hirshfeld_v_neigh, sR, d, c6_ref, r0_ref, alpha0_ref, &
                                       polynomial_expansion, vdw_omega_ref, &
                                       central_pol, central_omega, forces0 )

    implicit none

!   Input variables
    real*8, intent(in) :: rcut, r_buffer, hirshfeld_v_neigh(:), &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:), vdw_omega_ref !, hirshfeld_v(:), hirshfeld_v_neigh(:) !NOTE: uncomment this in final implementation
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: polynomial_expansion
!   Output variables
!   In-Out variables
    real*8, intent(inout) :: central_pol(:), central_omega(:), forces0(:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), f_damp(:), neighbor_alpha0(:), T_func(:), h_func(:), g_func(:), &
                           omegas(:), B_mat(:,:), rjs_H(:), xyz_H(:,:), work_arr(:), &
                           a_SCS(:,:), inner_damp(:), &
                           neighbor_sigma(:) !, rjs_0(:)
    real*8 :: time1, time2, Bohr, Hartree, &
              omega, pi, &
              r_vdw_i, r_vdw_j, t1, t2, &
              sigma_ij, s_i, s_j, omega_ref, xyz_i(1:3), xyz_j(1:3), mult1_i, mult1_j, mult2, pol1, r_buf, rb
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, om, n_tot
    integer :: i, i0, i1, i2, i3, j, j1, j2, j3, k, k2, k3, k4, a, a2, c1, c2, c3, lwork, b, p, q, r, k_i, k_j
               
!    LOCAL TEST stuff:
    integer, allocatable :: sub_neighbors_list(:), n_sub_neigh(:), p_list(:)
    integer :: n_sub_sites, n_sub_pairs, s

!   MBD stuff:
!    real*8, allocatable :: T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), AT(:,:,:), AT_n(:,:,:,:), energy_series(:,:), &
!                           integrand(:), AT_n_f(:,:,:,:), f_damp_der_SCS(:), dT_LR(:,:), G_mat(:,:,:), force_series(:,:), &
!                           total_energy_series(:,:), total_integrand(:), rjs_0(:), &
!                           T_mbd(:), r0_ii_mbd(:), neighbor_alpha0_mbd(:), omegas_mbd(:), rjs_0_mbd(:), xyz_0_mbd(:,:), &
!                           xyz_mbd(:,:), rjs_mbd(:), d_der(:,:), dT_mbd(:), f_damp_der_mbd(:), a_mbd(:), da_mbd(:), &
!                           a_iso(:,:), o_p(:), da_iso(:,:,:),  o_mbd(:), sub_2b_list(:), &
!                           xyz_2b(:,:), rjs_2b(:), r0_ii_2b(:), neighbor_alpha0_2b(:), f_damp_SCS_2b(:), &
!                           a_2b(:), r0_ii_SCS_2b(:), C6_2b(:), da_2b(:), T_SR(:), T_SR_mult(:), d_arr_i(:), d_arr_o(:), &
!                           d_mult_i(:), d_mult_o(:), dT_SR_mult(:,:), d_dmult_i(:,:), d_dmult_o(:,:), do_mbd(:), &
!                           hirshfeld_sub_neigh(:), o_2b(:), do_2b(:)
!    real*8 :: a_mbd_i, a_mbd_j, da_i, da_j, pol1, E_TS, f_damp_der_2b, dr_vdw_i, &
!              dr_vdw_j, forces_TS, dC6_2b, mult1_i, mult1_j, mult2, dmult1_i(1:3), dmult1_j(1:3), dmult2(1:3), hv_p_der, &
!              hv_q_der, do_pref
!    integer :: n_mbd_sites, n_mbd_pairs, n_2b_sites
!    integer, allocatable :: n_mbd_neigh(:), mbd_neighbors_list(:), p_mbd(:)


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
    !allocate( hirshfeld_v_cart_der_H(1:3, n_pairs) )

    !hirshfeld_v_cart_der_H = 0.d0

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
    !if ( do_derivatives .and. do_hirshfeld_gradients ) then
    !hirshfeld_v_cart_der_H = hirshfeld_v_cart_der_ji
    !write(*,*) "hirshfeld_der"
    !write(*,*) hirshfeld_v_cart_der_H(1:n_neigh(1))
    !end if

    ! HACK END

    !write(*,*) "hirshfeld"
    !do p = 1, n_neigh(1)
    !  write(*,*) hirshfeld_v_neigh(p)
    !end do
    !write(*,*) "der"
    !do p = 1, n_neigh(1)
    !  write(*,*) hirshfeld_v_cart_der_H(1,p)
    !end do

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411d0
    Hartree = 27.211386024367243d0
    pi = acos(-1.d0)

    !write(*,*) "rcut", rcut/Bohr

!   Number of frequencies
    !n_freq = size(alpha_SCS0, 2)
      
    !write(*,*) "Starting local calculation"

    !E_MBD = 0.d0

    !alpha_SCS0 = 0.d0
    !allocate( alpha_SCS_full(1:3*n_sites,1:3,1:2) )
    !if ( do_derivatives ) then
    !  dalpha_full = 0.d0
    !end if

    !alpha_SCS_full = 0.d0
    
    !allocate( central_pol(1:n_sites) )
    !central_pol = 0.d0
    !allocate( central_omega(1:n_sites) )
    !central_omega = 0.d0
    
    if ( rcut < r_buffer ) then
      r_buf = rcut
    else
      r_buf = r_buffer
    end if
      
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
        if (rjs(n_tot+k_i) .le. rcut ) then
          n_sub_sites = n_sub_sites + 1
          n_sub_pairs = n_sub_pairs + 1
          xyz_i = xyz(:,n_tot+k_i)/Bohr
          k_j = 0
          do j3 = 1, n_neigh(i)
            k_j = k_j + 1
            if ( rjs(n_tot+k_j) .le. 2*rcut ) then !2*rcut
              if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut/Bohr ) then
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
      allocate( inner_damp(1:n_sub_pairs) )
      !allocate( a_iso(1:n_sub_sites,1:2) )
      !allocate( o_p(1:n_sub_sites) )
      !if ( do_derivatives .and. do_hirshfeld_gradients ) then
      !  allocate( hirshfeld_v_sub_der(1:3,1:n_sub_sites) )
      !  hirshfeld_v_sub_der = 0.d0
      !end if
        
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

        inner_damp = 1.d0
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
          if ( rjs(n_tot+k_i) .le. rcut ) then
            p = p+1 
            k2 = k2+1
            !rjs_0(k2) = rjs(n_tot+k_i)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = i2
            !if ( do_derivatives .and. do_hirshfeld_gradients ) then
            !  hirshfeld_v_sub_der(1:3,p) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            !end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (vdw_omega_ref * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
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
            if ( rjs(n_tot+k_i) > rcut-r_buf ) then
              rb = (rjs(n_tot+k_i)-rcut+r_buf)/(r_buf)
              mult1_i = mult1_i * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
            end if
            k_j = 0
            q = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              if ( rjs(n_tot+k_j) .le. 2*rcut ) then !2*rcut
                if ( rjs(n_tot+k_j) .le. rcut ) then
                  q = q+1
                end if
                if (i3 .ne. j3) then
                  xyz_j = xyz(:,n_tot+k_j)/Bohr
                  if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut/Bohr ) then
                    n_sub_neigh(p) = n_sub_neigh(p) + 1
                    j = neighbors_list(n_tot+k_j)
                    k2 = k2+1
                    !rjs_0(k2) = rjs(n_tot+k_j)
                    if ( rjs(n_tot+k_j) .le. rcut ) then
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
                        (1.d0 + (vdw_omega_ref * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
                    end if
                    neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
                    xyz_H(:,k2) = xyz_j-xyz_i
                    rjs_H(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                    r_vdw_j = r0_ii(k2)
                    f_damp(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                    s_j = neighbor_sigma(k2)
                    sigma_ij = sqrt(s_i**2 + s_j**2)
                    g_func(k2) = erf(rjs_H(k2)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij) * exp(-rjs_H(k2)**2.d0/sigma_ij**2)
                    if ( rjs_H(k2)*Bohr < 2.d0 .and. rjs(n_tot+k_j) > rcut-r_buf ) then
                      rb = (rjs_H(k2)*Bohr)/2.d0
                      inner_damp(k2) = 10.d0*rb**3-15.d0*rb**4+6.d0*rb**5
                    end if
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
                        if ( rjs(n_tot+k_j) .le. rcut ) then
                          mult1_j = 1.d0
                          mult2 = 1.d0
                          if ( rjs(n_tot+k_j) > rcut-r_buf ) then
                            rb = (rjs(n_tot+k_j)-rcut+r_buf)/(r_buf)
                            mult1_j = mult1_j * &
                                           (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                          end if
                          if ( rjs_H(k2) > (rcut-r_buf)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                            rb = (rjs_H(k2)-(rcut-r_buf)/Bohr)/(r_buf/Bohr)
                            mult2 = mult2 * &
                                           (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
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
                                                    mult1_i * ( 1.d0 - mult1_j ) * mult2 * inner_damp(k2)
                        else ! rjs(n_tot+k_j) < rcut
                          mult1_j = 1.d0
                          mult2 = 1.d0
                          if ( rjs(n_tot+k_j) > 2.d0*rcut-r_buf ) then
                            rb = (rjs(n_tot+k_j)-2.d0*rcut+r_buf)/(r_buf)
                            mult1_j = mult1_j * &
                                           (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                          end if
                          if ( rjs_H(k2) > (rcut-r_buf)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                            rb = (rjs_H(k2)-(rcut-r_buf)/Bohr)/(r_buf/Bohr)
                            mult2 = mult2 * &
                                          (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                          end if
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                      mult1_j * mult2 * inner_damp(k2)
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
            if (central_pol(i) .le. pol1) then
              write(*,*) "ERROR: frequency dependency approximation fails for"
              write(*,*) "one of the central polarizabilities."
              write(*,*) "You can try using a larger value for"
              write(*,*) "vdw_omega_ref or larger vdw_scs_rcut."
              return
            else
              central_omega(i) = vdw_omega_ref*omega_ref/sqrt(central_pol(i)/pol1-1.d0)
            end if  
            !write(*,*) "Polynomial central polarizability", i0, central_pol(i), central_omega(i)
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
            if ( pol1 .le. 0.d0 ) then
              write(*,*) "ERROR: At least one of the SCS polarizabilities"
              write(*,*) "received a non-positive value. This is probably due to"
              write(*,*) "very small interatomic distances where the local SCS"
              write(*,*) "cycle for MBD fails. We recommend switching off"
              write(*,*) "vdW corrections, using Tkatchenko-Scheffler or"
              write(*,*) "trying to increase the value of vdw_scs_rcut."
              return
            end if
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
            if ( central_pol(i) .le. 0.d0 ) then
              write(*,*) "ERROR: At least one of the SCS polarizabilities"
              write(*,*) "received a non-positive value. This is probably due to"
              write(*,*) "very small interatomic distances where the SCS"
              write(*,*) "cycle for MBD fails. We recommend switching off"
              write(*,*) "vdW corrections, using Tkatchenko-Scheffler or"
              write(*,*) "trying to increase the value of vdw_scs_rcut."
              return
            end if
            central_pol(i) = central_pol(i)/3.d0
            !write(*,*) i0, central_pol(i)
            if (central_pol(i) .le. pol1) then
              write(*,*) "ERROR: frequency dependency approximation fails for"
              write(*,*) "one of the central polarizabilities."
              write(*,*) "You can try using a larger value for"
              write(*,*) "vdw_omega_ref."
              return
            else
              central_omega(i) = vdw_omega_ref*omega_ref/sqrt(central_pol(i)/pol1-1.d0)
            end if
          end if
        
        end if

        !call cpu_time(time2)

        !write(*,*) "Timing for solving the polarizability of central atom:", time2-time1

      end do
        
      deallocate( sub_neighbors_list, n_sub_neigh, p_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, &
                  omegas, T_func, B_mat, f_damp, g_func, h_func, a_SCS, ipiv, work_arr, inner_damp )
      !deallocate( a_iso )
      !if ( do_derivatives .and. do_hirshfeld_gradients ) then
      !  deallocate( hirshfeld_v_sub_der )
      !end if
        
      deallocate( b_i, d_vec )
      !deallocate( ia, ja, val )        

    end do 

    call cpu_time(time2)

    !write(*,*) "Central polarizabilities timing", time2-time1



  end subroutine
!**************************************************************************







      
!******************************************************************************************
! This is where you would break this into get_scs and get_mbd for two separate subroutines:
! Store central_pol and central_omega and pass them to get_mbd
!******************************************************************************************

! TEST !!!!!!!!!!!!!!!!
!central_pol = 10.d0
!central_omega = 0.5d0

!**************************************************************************
  subroutine get_mbd_energies_and_forces( hirshfeld_v_cart_der_ji, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, rcut_loc, rcut_mbd, rcut_2b, r_buffer, rjs, xyz, &
                                       hirshfeld_v_neigh, sR, d, c6_ref, r0_ref, alpha0_ref, do_derivatives, &
                                       do_hirshfeld_gradients, polynomial_expansion, do_nnls, n_freq, n_order, &
                                       vdw_omega_ref, central_pol, central_omega, &
                                       energies, forces0, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: rcut, r_buffer, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), rcut_loc, rcut_mbd, rcut_2b, &
                          hirshfeld_v_cart_der_ji(:,:), &
                          alpha0_ref(:), vdw_omega_ref !, hirshfeld_v(:), hirshfeld_v_neigh(:) !NOTE: uncomment this in final implementation
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:), n_freq, n_order
    logical, intent(in) :: do_derivatives, do_hirshfeld_gradients, polynomial_expansion, do_nnls
!   Output variables
    real*8, intent(out) :: virial(1:3, 1:3)
!   In-Out variables
    real*8, intent(inout) :: energies(:), forces0(:,:), central_pol(:), central_omega(:), &
                             hirshfeld_v_neigh(:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), f_damp(:), neighbor_alpha0(:), T_func(:), h_func(:), g_func(:), &
                           omegas(:), B_mat(:,:), rjs_H(:), xyz_H(:,:), work_arr(:), dT(:), f_damp_der(:), &
                           g_func_der(:), h_func_der(:), coeff_der(:), coeff_fdamp(:), hirshfeld_v_cart_der_H(:,:), &
                           a_SCS(:,:), da_SCS(:,:), b_der(:,:), dB_mat(:,:), &
                           neighbor_sigma(:), hirshfeld_v_sub_der(:,:), inner_damp(:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              r_vdw_i, r_vdw_j, t1, t2, r_buf_scs, r_buf_mbd, r_buf_2b, r_buf_loc, &
              sigma_ij, coeff_h_der, dg, dh, s_i, s_j, terms, omega_ref, xyz_i(1:3), xyz_j(1:3), rjs_i, rjs_j
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, om, n_tot
    integer :: i, i0, i1, i2, i3, j, j1, j2, j3, k, k2, k3, k4, a, a2, c1, c2, c3, lwork, b, p, q, r, k_i, k_j

!    LOCAL TEST stuff:
    integer, allocatable :: sub_neighbors_list(:), n_sub_neigh(:), p_list(:)
    integer :: n_sub_sites, n_sub_pairs, s

!   MBD stuff:
    real*8, allocatable :: T_LR(:,:), r0_ii_SCS(:), f_damp_SCS(:), AT(:,:,:), AT_n(:,:,:,:), energy_series(:,:), &
                           integrand(:), AT_n_f(:,:,:,:), f_damp_der_SCS(:), dT_LR(:,:), G_mat(:,:,:), force_series(:,:), &
                           total_energy_series(:,:), total_integrand(:), rjs_0(:), &
                           T_mbd(:), r0_ii_mbd(:), neighbor_alpha0_mbd(:), omegas_mbd(:), rjs_0_mbd(:), xyz_0_mbd(:,:), &
                           xyz_mbd(:,:), rjs_mbd(:), d_der(:,:), dT_mbd(:), f_damp_der_mbd(:), a_mbd(:), da_mbd(:), &
                           a_iso(:,:), o_p(:), da_iso(:,:,:),  o_mbd(:), sub_2b_list(:), hirshfeld_2b_neigh(:), &
                           xyz_2b(:,:), rjs_2b(:), r0_ii_2b(:), neighbor_alpha0_2b(:), f_damp_SCS_2b(:), &
                           a_2b(:), r0_ii_SCS_2b(:), C6_2b(:), da_2b(:), T_SR(:), T_SR_mult(:), d_arr_i(:), d_arr_o(:), &
                           d_mult_i(:), d_mult_o(:), dT_SR_mult(:,:), d_dmult_i(:,:), d_dmult_o(:,:), do_mbd(:), &
                           hirshfeld_sub_neigh(:), o_2b(:), do_2b(:), hirshfeld_v_mbd_der(:,:), hirshfeld_mbd_neigh(:), &
                           hirshfeld_v_2b_der(:,:), dr0_ii_SCS(:), dr0_ii_SCS_2b(:), V_int(:,:), T_LR_mult_0i(:), &
                           T_LR_mult_0j(:), T_LR_mult_ij(:), dT_LR_mult_0i(:), dT_LR_mult_0j(:), dT_LR_mult_ij(:), &
                           r6_mult(:), dr6_mult(:), T_LR_mult_0ij(:), T_LR_mult_0ji(:), dT_LR_mult_0ij(:), &
                           dT_LR_mult_0ji(:), xyz_2b_tot(:,:), rjs_2b_tot(:), r0_ii_2b_tot(:), neighbor_alpha0_2b_tot(:), &
                           f_damp_SCS_2b_tot(:), hirshfeld_2b_tot_neigh(:), a_2b_tot(:), o_2b_tot(:), &
                           r0_ii_SCS_2b_tot(:), c6_2b_tot(:), r6_mult_2b_tot(:), r6_mult_0i(:), r6_mult_0j(:), &
                           dr6_mult_0i(:), dr6_mult_0j(:)
    real*8 :: a_mbd_i, a_mbd_j, da_i, da_j, pol1, E_TS, f_damp_der_2b, dr_vdw_i, &
              dr_vdw_j, forces_TS, dC6_2b, mult1_i, mult1_j, mult2, dmult1_i(1:3), dmult1_j(1:3), dmult2(1:3), hv_p_der, &
              hv_q_der, do_pref, rb, inner_damp_der, rjs_0_i, rcut_forces, o_i, o_j, do_i, do_j, T_LR_mult_i, T_LR_mult_j, &
              dT_LR_mult_i, dT_LR_mult_j, dT_LR_mult_0ij_2, dT_LR_mult_0ji_2, a_i, a_j, E_TS_tot, r6_der, ac2, ac3, ac4
              
    integer :: n_mbd_sites, n_mbd_pairs, n_2b_sites, n_2b_tot_sites, n_2b_tot_pairs
    integer, allocatable :: n_mbd_neigh(:), mbd_neighbors_list(:), p_mbd(:), sub_2b_tot_list(:), n_2b_tot_neigh(:), &
                            p_2b_tot(:)
    real*8 :: polyfit(1:15)
    integer :: n_degree
    real*8, allocatable :: B_pol(:,:), B_mult(:,:), b_i(:,:), d_vec(:,:), val_xv(:,:)
!   NNLS stuff;
    
    real*8, allocatable :: A_nnls(:,:), b_nnls(:), coeff_nnls(:), work_nnls(:), omegas_nnls(:), integrand_nnls(:), &
                           total_integrand_nnls(:), work_integrand(:), trace_nnls(:), force_series0(:,:), full_integrand(:)
    integer, allocatable :: ind_nnls(:)
    real*8 :: res_nnls, E_tot, denom
    integer :: mode_nnls
    logical :: do_total_energy = .true. ! Finite difference testing purposes

!central_pol = 10.d0
!central_omega = 0.5d0

!hirshfeld_v_neigh = 1.d0

ac4 = 8.d0*(2.d0)**(1.d0/2.d0) - 8.d0
ac3 = -2.d0*ac4 -2.d0
ac2 = 3.d0 + ac4

!central_omega(1:size(central_omega)/3) = 0.5d0
!central_omega(size(central_omega)/3:2*size(central_omega)/3) = 0.8d0
!central_omega(2*size(central_omega)/3:size(central_omega)) = 1.4d0

! TESTING: Change neighbor_alpha0_mbd... and neighbor_alpha0_2b... back to central_pol in a_2b and a_mbd parts
! Test G_mat with finite difference; it should be equal to the analytical solution now

    if ( rcut_loc > rcut ) then
      write(*,*) "ERROR: vdw_loc_rcut should be at most vdw_scs_rcut!"
      return
    end if
    if ( rcut_2b < rcut_mbd ) then
      write(*,*) "ERROR: vdw_2b_rcut should be larger than vdw_mbd_rcut!"
      return
    end if
    if ( n_freq < 2 ) then
      write(*,*) "ERROR: not enough frequency points for integration."
      write(*,*) "Try vdw_mbd_nfreq > 1."
      return
    end if
    if ( n_order < 2 ) then
      write(*,*) "WARNING: vdw_norder < 2 results in zero energies"
      write(*,*) "and forces. Try larger value (vdw_norder = n"
      write(*,*) "results in interactions up to n-th order)."
    end if

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

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

    !write(*,*) "hirshfeld"
    !do p = 1, n_neigh(1)
    !  write(*,*) hirshfeld_v_neigh(p)
    !end do
    !write(*,*) "der"
    !do p = 1, n_neigh(1)
    !  write(*,*) hirshfeld_v_cart_der_H(1,p)
    !end do
    Bohr = 0.5291772105638411d0
    Hartree = 27.211386024367243d0
    pi = acos(-1.d0)

    !write(*,*) "rcut", rcut/Bohr

!   Number of frequencies
    !n_freq = size(alpha_SCS0, 2)

    !write(*,*) "Starting local calculation"

    E_MBD = 0.d0

    !alpha_SCS0 = 0.d0
    !allocate( alpha_SCS_full(1:3*n_sites,1:3,1:2) )

    !alpha_SCS_full = 0.d0

    !allocate( central_pol(1:n_sites) )
    !central_pol = 0.d0
    !allocate( central_omega(1:n_sites) )
    !central_omega = 0.d0

    !r_buffer = 0.5d0
    
    if ( rcut < r_buffer ) then
      r_buf_scs = rcut
    else
      r_buf_scs = r_buffer
    end if

    if ( rcut_loc < r_buffer ) then
      r_buf_loc = rcut_loc
    else
      r_buf_loc = r_buffer
    end if
    
    if ( rcut_mbd < r_buffer ) then
      r_buf_mbd = rcut_mbd
    else
      r_buf_mbd = r_buffer
    end if

    r_buf_mbd = 0.d0
    
    if ( rcut_2b-rcut_mbd < r_buffer ) then
      r_buf_2b = rcut_2b-rcut_mbd
    else
      r_buf_2b = r_buffer
    end if

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
        if (rjs(n_tot+k_i) .le. rcut ) then
          n_sub_sites = n_sub_sites + 1
          n_sub_pairs = n_sub_pairs + 1
          xyz_i = xyz(:,n_tot+k_i)/Bohr
          k_j = 0
          do j3 = 1, n_neigh(i)
            k_j = k_j + 1
            if ( rjs(n_tot+k_j) .le. 2*rcut ) then !2*rcut
              if (i3 .ne. j3) then
                xyz_j = xyz(:,n_tot+k_j)/Bohr
                if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut/Bohr ) then
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
      allocate( inner_damp(1:n_sub_pairs) )
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
        allocate( hirshfeld_v_sub_der(1:3,1:n_sub_pairs) )
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
      !T_SR = 0.d0
      !T_SR_mult = 0.d0
      !d_arr_i = 0.d0
      !d_arr_o = 0.d0
      !d_mult_i = 0.d0
      !d_mult_o = 0.d0
      !dT_SR_mult = 0.d0
      !d_dmult_i = 0.d0
      !d_dmult_o = 0.d0

      do om = 1, 2

        T_SR = 0.d0
        T_SR_mult = 0.d0
        d_arr_i = 0.d0
        d_arr_o = 0.d0
        d_mult_i = 0.d0
        d_mult_o = 0.d0
        dT_SR_mult = 0.d0
        d_dmult_i = 0.d0
        d_dmult_o = 0.d0

        xyz_H = 0.d0
        rjs_H = 0.d0
        r0_ii = 0.d0
        neighbor_alpha0 = 0.d0
        neighbor_sigma = 0.d0
        hirshfeld_sub_neigh = 0.d0

        f_damp = 0.d0
        g_func = 0.d0
        h_func = 0.d0
        inner_damp = 1.d0

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
          if ( rjs(n_tot+k_i) .le. rcut ) then
            p = p+1 
            k2 = k2+1
            rjs_0(k2) = rjs(n_tot+k_i)
            hirshfeld_sub_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
            s = neighbor_species(n_tot+k_i)
            sub_neighbors_list(k2) = i2
            if ( do_derivatives .and. do_hirshfeld_gradients ) then
              hirshfeld_v_sub_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
            end if
            n_sub_neigh(p) = n_sub_neigh(p) + 1
            p_list(k2) = p
            r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
            omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
            if ( om == 2 ) then
              neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)
            else
              neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_i)) / &
                (1.d0 + (vdw_omega_ref * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
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
            if ( rjs(n_tot+k_i) > rcut-r_buf_scs ) then
              rb = (rjs(n_tot+k_i)-rcut+r_buf_scs)/(r_buf_scs)
              mult1_i = mult1_i * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
              dmult1_i = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                              * ( -xyz(:,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_scs/Bohr))
            end if
            k_j = 0
            q = 0
            do j3 = 1, n_neigh(i)
              k_j = k_j+1
              if ( rjs(n_tot+k_j) .le. 2*rcut ) then !2*rcut
                if ( rjs(n_tot+k_j) .le. rcut ) then
                  q = q+1
                end if
                if (i3 .ne. j3) then
                  xyz_j = xyz(:,n_tot+k_j)/Bohr
                  if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut/Bohr ) then
                    n_sub_neigh(p) = n_sub_neigh(p) + 1
                    j = neighbors_list(n_tot+k_j)
                    k2 = k2+1
                    rjs_0(k2) = rjs(n_tot+k_j)
                    hirshfeld_sub_neigh(k2) = hirshfeld_v_neigh(n_tot+k_j)
                    if ( rjs(n_tot+k_j) .le. rcut ) then
                      p_list(k2) = q
                    end if
                    if ( do_derivatives .and. do_hirshfeld_gradients ) then
                      hirshfeld_v_sub_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_j)
                    end if    
                    s = neighbor_species(n_tot+k_j)
                    sub_neighbors_list(k2) = j                        
                    r0_ii(k2) = r0_ref(s) / Bohr * hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                    omegas(k2) = (4.d0 * c6_ref(s)/(Hartree*Bohr**6)) / (3.d0*(alpha0_ref(s)/Bohr**3)**2)
                    if ( om == 2 ) then
                      neighbor_alpha0(k2) = alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)
                    else
                      neighbor_alpha0(k2) = (alpha0_ref(s) / Bohr**3 * hirshfeld_v_neigh(n_tot+k_j)) / &
                        (1.d0 + (vdw_omega_ref * omega_ref/omegas(k2))**2) ! omega = 2*omega_p
                    end if
                    neighbor_sigma(k2) = (sqrt(2.d0/pi) * neighbor_alpha0(k2)/3.d0)**(1.d0/3.d0)
                    xyz_H(:,k2) = xyz_j-xyz_i
                    rjs_H(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                    r_vdw_j = r0_ii(k2)
                    f_damp(k2) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k2)/(sR*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                    s_j = neighbor_sigma(k2)
                    sigma_ij = sqrt(s_i**2 + s_j**2)
                    g_func(k2) = erf(rjs_H(k2)/sigma_ij) - 2.d0/sqrt(pi) * (rjs_H(k2)/sigma_ij) * exp(-rjs_H(k2)**2.d0/sigma_ij**2)
                    if ( rjs_H(k2)*Bohr < 2.d0 .and. rjs(n_tot+k_j) > rcut-r_buf_scs ) then
                      rb = (rjs_H(k2)*Bohr)/2.d0
                      inner_damp(k2) = 10.d0*rb**3-15.d0*rb**4+6.d0*rb**5
                    end if
                    k3 = 9*(k2-1)
                    if ( rjs(n_tot+k_j) .le. rcut ) then
                      mult1_j = 1.d0
                      dmult1_j = 0.d0
                      mult2 = 1.d0
                      dmult2 = 0.d0
                      if ( rjs(n_tot+k_j) > rcut-r_buf_scs ) then
                        rb = (rjs(n_tot+k_j)-rcut+r_buf_scs)/(r_buf_scs)
                        mult1_j = mult1_j * &
                                       (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                        dmult1_j = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                                * ( -xyz(:,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_scs/Bohr))
                      end if
                      if ( rjs_H(k2) > (rcut-r_buf_scs)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                        rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                        mult2 = mult2 * &
                                       (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                        if ( i2 == i ) then
                          rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                          dmult2 = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                            * ( -xyz_H(:,k2)/rjs_H(k2)/(r_buf_scs/Bohr))
                        end if
                        if ( j == i ) then
                          rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                          dmult2 = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                            * ( xyz_H(:,k2)/rjs_H(k2)/(r_buf_scs/Bohr))
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
                      if ( rjs(n_tot+k_j) > 2.d0*rcut-r_buf_scs ) then
                        rb = (rjs(n_tot+k_j)-2.d0*rcut+r_buf_scs)/(r_buf_scs)
                        mult1_j = mult1_j * &
                                       (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                        dmult1_j = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                                * ( -xyz(:,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_scs/Bohr))
                      end if
                      if ( rjs_H(k2) > (rcut-r_buf_scs)/Bohr ) then !.and. rjs_H(k2) < (rcut-r_buffer)/Bohr ) then
                        rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                        mult2 = mult2 * &
                                      (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 )
                        if ( i2 == i ) then
                          rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                          dmult2 = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                            * ( -xyz_H(:,k2)/rjs_H(k2)/(r_buf_scs/Bohr))
                        end if
                        if ( j == i ) then
                          rb = (rjs_H(k2)-(rcut-r_buf_scs)/Bohr)/(r_buf_scs/Bohr)
                          dmult2 = (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4 ) &
                            * ( xyz_H(:,k2)/rjs_H(k2)/(r_buf_scs/Bohr))
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
                        if ( rjs(n_tot+k_j) .le. rcut ) then
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
                                                    d_mult_i(k2) * inner_damp(k2)
                          d_arr_i(k3) = (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                    g_func(k2) + h_func(k3)) * neighbor_alpha0(k2)
                        else
                          d_vec(3*(p-1)+c1,c2) = d_vec(3*(p-1)+c1,c2) - (1.d0-f_damp(k2)) * (-T_func(k3) * &
                                                      g_func(k2) + h_func(k3)) * neighbor_alpha0(k2) * &
                                                      d_mult_o(k2) * inner_damp(k2)
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
        !  write(*,*) "d_vec"
        !  write(*,*) "B_mat"
        !  do p = 1, 3*n_sub_sites
        !    write(*,*) B_mat(p,:)
        !    write(*,*) d_vec(p,:)
        !  end do
        !  write(*,*) "d_mult_o", d_mult_o
        !  write(*,*) "d_dmult_o", d_dmult_o
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
            if ( a_iso(p,om) .le. 0.d0 ) then
              write(*,*) "ERROR: At least one of the SCS polarizabilities"
              write(*,*) "received a non-positive value. This is probably due to"
              write(*,*) "very small interatomic distances where the SCS"
              write(*,*) "cycle for MBD fails. We recommend switching off"
              write(*,*) "vdW corrections, using Tkatchenko-Scheffler or"
              write(*,*) "trying to increase the value of vdw_scs_rcut."
              return
            end if
          end do
          a_iso(:,om) = a_iso(:,om)/3.d0

          if ( om == 2 ) then
            k2 = 0
            do p = 1, n_sub_sites
              if (a_iso(p,2) .le. a_iso(p,1)) then
                write(*,*) "WARNING: frequency dependency approximation fails for"
                write(*,*) "one of the local atoms. Using central characteristic"
                write(*,*) "frequency for now. You can try using a larger value for"
                write(*,*) "vdw_omega_ref or vdw_scs_rcut. If you see this warning,"
                write(*,*) "your input structure probably has at least one very small"
                write(*,*) "intertomic distance where the SCS cycle in MBD"
                write(*,*) "calculation fails. For high energy calculations"
                write(*,*) "we recommend switching off vdW corrections." 
                i2 = p_list(k2)
                i1 = modulo(i2-1, n_sites0) + 1
                o_p(p) = central_omega(i1)
              else
                o_p(p) = vdw_omega_ref*omega_ref/sqrt(a_iso(p,2)/a_iso(p,1)-1.d0)
              end if
              k2 = k2 + n_sub_neigh(p)
            end do
          end if

          deallocate( val_xv, B_mult )
          !deallocate( myidx )

        else


          !if ( i == 1 .and. om == 2 ) then
          !  write(*,*) "B_mat"
          !  do p = 1, n_sub_sites
          !    do c1 = 1, 3
          !      write(*,*) B_mat(3*(p-1)+c1,:)
          !    end do
          !  end do
          !end if
          ! This is the exact solution:
          a_SCS = d_vec

          call dsysv( 'U', 3*n_sub_sites, 3, B_mat, 3*n_sub_sites, ipiv, a_SCS, 3*n_sub_sites, work_arr, &
                      12*n_sub_sites, info )

          do p = 1, n_sub_sites
            do c1 = 1, 3
              a_iso(p,om) = a_iso(p,om) + a_SCS(3*(p-1)+c1,c1)
            end do
            if ( a_iso(p,om) .le. 0.d0 ) then
              write(*,*) "ERROR: At least one of the SCS polarizabilities"
              write(*,*) "received a non-positive value. This is probably due to"
              write(*,*) "very small interatomic distances where the SCS"
              write(*,*) "cycle for MBD fails. We recommend switching off"
              write(*,*) "vdW corrections, using Tkatchenko-Scheffler or"
              write(*,*) "trying to increase the value of vdw_scs_rcut."
              write(*,*) "p, om, a_iso", p, om, a_iso(p,om)
              return
            end if
          end do
          a_iso(:,om) = a_iso(:,om)/3.d0
        
          if ( om == 2 ) then
            k2 = 0
            do p = 1, n_sub_sites
              if (a_iso(p,2) .le. a_iso(p,1)) then
                write(*,*) "WARNING: frequency dependency approximation fails for"
                write(*,*) "one of the local atoms. Using central characteristic"
                write(*,*) "frequency for now. You can try using a larger value for"
                write(*,*) "vdw_omega_ref or vdw_scs_rcut. If you see this warning,"
                write(*,*) "your input structure probably has at least one very small"
                write(*,*) "intertomic distance where the SCS cycle in MBD"
                write(*,*) "calculation fails. For high energy calculations"                          
                write(*,*) "we recommend switching off vdW corrections."
                i2 = p_list(k2)
                i1 = modulo(i2-1, n_sites0) + 1
                o_p(p) = central_omega(i1)
              else
                o_p(p) = vdw_omega_ref*omega_ref/sqrt(a_iso(p,2)/a_iso(p,1)-1.d0)
              end if
              k2 = k2 + n_sub_neigh(p)
            end do
          end if

          if (i == 1 .and. om == 2) then
            write(*,*) "a_iso"
            do p = 1, n_sub_sites
              write(*,*) a_iso(p,2)
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

if ( .false. ) then

          ! MBD for local polarizabilities:
          ! At least for now: rcut <= rcut_mbd <= rcut_2b
          !rcut_mbd = 5.d0
          !rcut_2b = 9.d0
          n_mbd_sites = 0
          n_mbd_pairs = 0
        
          k_i = 0
          do i3 = 1, n_neigh(i)
            k_i = k_i + 1
            if (rjs(n_tot+k_i) .le. rcut_mbd ) then
              n_mbd_sites = n_mbd_sites + 1
              n_mbd_pairs = n_mbd_pairs + 1
              xyz_i = xyz(:,n_tot+k_i)/Bohr
              k_j = 0
              do j3 = 1, n_neigh(i)
                k_j = k_j + 1
                if ( rjs(n_tot+k_j) .le. rcut_mbd ) then
                  if (i3 .ne. j3) then
                    xyz_j = xyz(:,n_tot+k_j)/Bohr
                    if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. (rcut_mbd)/Bohr ) then
                      n_mbd_pairs = n_mbd_pairs + 1
                    end if
                  end if
                end if
              end do
            end if
          end do
        
end if

          k_i = 0
          n_2b_sites = 0
          do j2 = 1, n_neigh(i)
            k_i = k_i + 1
            if (rjs(n_tot+k_i) .le. rcut_2b .and. rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd) then
              n_2b_sites = n_2b_sites + 1
            end if
          end do

          !n_order = 4
          !n_freq = 12

if ( .false. ) then

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
          allocate( hirshfeld_mbd_neigh(1:n_mbd_pairs) )
          allocate( T_LR_mult_0i(1:n_mbd_pairs) )
          allocate( T_LR_mult_0j(1:n_mbd_pairs) )
          allocate( T_LR_mult_ij(1:n_mbd_pairs) )
          allocate( T_LR_mult_0ij(1:n_mbd_pairs) )
          allocate( T_LR_mult_0ji(1:n_mbd_pairs) )
          hirshfeld_mbd_neigh = 0.d0
          if ( do_derivatives ) then
            allocate( da_mbd(1:n_mbd_pairs) )
            allocate( AT_n_f(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_order-1,1:n_freq) )
            allocate( dT_mbd(1:9*n_mbd_pairs) )
            allocate( f_damp_der_mbd(1:n_mbd_pairs) )
            allocate( f_damp_der_SCS(1:n_mbd_pairs) )
            allocate( dT_LR(1:3*n_mbd_sites,1:3*n_mbd_sites) )
            allocate( force_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
            if ( do_total_energy ) then
              allocate( total_energy_series(1:3*n_mbd_sites,1:3*n_mbd_sites) )
              allocate( total_integrand(1:n_freq) )
            end if
            allocate( do_mbd(1:n_mbd_pairs) )
            allocate( dr0_ii_SCS(1:n_mbd_pairs) )
            allocate( dT_LR_mult_0i(1:n_mbd_pairs) )
            allocate( dT_LR_mult_0j(1:n_mbd_pairs) )
            allocate( dT_LR_mult_ij(1:n_mbd_pairs) )
            allocate( dT_LR_mult_0ij(1:n_mbd_pairs) )
            allocate( dT_LR_mult_0ji(1:n_mbd_pairs) )
          end if
          if ( do_derivatives .and. do_hirshfeld_gradients ) then
            allocate( hirshfeld_v_mbd_der(1:3,1:n_mbd_pairs) )
            hirshfeld_v_mbd_der = 0.d0
          end if

end if

        
          allocate( sub_2b_list(1:n_2b_sites) )
          allocate( xyz_2b(1:3,1:n_2b_sites) )
          allocate( rjs_2b(1:n_2b_sites) )
          allocate( r0_ii_2b(1:n_2b_sites) )
          allocate( neighbor_alpha0_2b(1:n_2b_sites) )
          allocate( hirshfeld_2b_neigh(1:n_2b_sites) )
          allocate( a_2b(1:n_2b_sites) )
          allocate( o_2b(1:n_2b_sites) )
          allocate( r0_ii_SCS_2b(1:n_2b_sites) )
          allocate( f_damp_SCS_2b(1:n_2b_sites) )
          allocate( C6_2b(1:n_2b_sites) )
          allocate( r6_mult(1:n_2b_sites) )

if ( .false. ) then

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
          
          T_LR_mult_0i = 1.d0
          T_LR_mult_0j = 1.d0
          T_LR_mult_ij = 1.d0
          T_LR_mult_0ij = 1.d0
          T_LR_mult_0ji = 1.d0        

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
            if ( rjs(n_tot+k_i) .le. rcut_mbd ) then
              p = p+1
              k2 = k2+1
              s = neighbor_species(n_tot+k_i)
              mbd_neighbors_list(k2) = i2
              n_mbd_neigh(p) = n_mbd_neigh(p) + 1
              p_mbd(k2) = p
              if ( do_derivatives .and. do_hirshfeld_gradients ) then
                hirshfeld_v_mbd_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
              end if
              r0_ii_mbd(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
              xyz_mbd(:,k2) = xyz(:,n_tot+k_i)/Bohr
              xyz_i = xyz_mbd(:,k2)
              rjs_mbd(k2) = rjs(n_tot+k_i)/Bohr
              rjs_0_mbd(k2) = rjs(n_tot+k_i)/Bohr
              xyz_0_mbd(:,k2) = xyz(:,n_tot+k_i)/Bohr
              hirshfeld_mbd_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
              neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 ! * hirshfeld_v_neigh(n_tot+k_i)
              if ( rjs(n_tot+k_i) .le. rcut_loc-r_buf_loc ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                a_mbd(k2) = a_iso(r,2)
                o_mbd(k2) = o_p(r)
                r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
              else if ( rjs(n_tot+k_i) .le. rcut_loc .and. rjs(n_tot+k_i) > rcut_loc-r_buf_loc ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                rb = (rjs(n_tot+k_i)-rcut_loc+r_buf_loc)/r_buf_loc
                a_mbd(k2) = a_iso(r,2) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 ) + &
                            central_pol(i1) * &
                            !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2) * &
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                o_mbd(k2) = o_p(r) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                            central_omega(i1) * & 
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                r0_ii_mbd(k2) * (central_pol(i1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                                 !r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0) * &
                                ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
              else if ( rjs(n_tot+k_i) > rcut_loc .and. rjs(n_tot+k_i) .le. rcut_mbd) then
                a_mbd(k2) =  central_pol(i1)
                             !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2)
                o_mbd(k2) = central_omega(i1)
                r0_ii_SCS(k2) = r0_ii_mbd(k2) * (central_pol(i1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) 
                                             !* (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0)
              end if
              if ( rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd .and. &
                   rjs(n_tot+k_i) .le. rcut_mbd ) then
                rb = (rjs(n_tot+k_i)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                T_LR_mult_0i(k2) = (1.d0 & 
                               - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5)**(1.d0/n_order)
                !if ( rb < 0.5d0 ) then
                !T_LR_mult_0i(k2) = (1.d0 - ( ac2*rb**2 + ac3*rb**3 + ac4*rb**4 ) )
                !else
                !  T_LR_mult_0i(k2) = (1.d0 -3.d0*rb**2 +2.d0*rb**3)**(1.d0/2.d0)
                !end if
              else if ( rjs(n_tot+k_i) > rcut_mbd ) then
                T_LR_mult_0i(k2) = 0.d0
              end if
              T_LR_mult_i = T_LR_mult_0i(k2)
              !r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_mbd(k2)/neighbor_alpha0_mbd(k2)/hirshfeld_mbd_neigh(k2))**(1.d0/3.d0)
              r_vdw_i = r0_ii_SCS(k2)
              k_j = 0
              q = 0
              do j3 = 1, n_neigh(i)
                k_j = k_j+1
                if ( rjs(n_tot+k_j) .le. rcut_mbd ) then
                  q = q+1
                  if (i3 .ne. j3) then
                    xyz_j = xyz(:,n_tot+k_j)/Bohr
                    if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. (rcut_mbd)/Bohr ) then
                      n_mbd_neigh(p) = n_mbd_neigh(p) + 1
                      j = neighbors_list(n_tot+k_j)
                      j1 = modulo(j-1, n_sites0) + 1
                      k2 = k2+1
                      s = neighbor_species(n_tot+k_j)
                      mbd_neighbors_list(k2) = j
                      p_mbd(k2) = q
                      if ( do_derivatives .and. do_hirshfeld_gradients ) then
                        hirshfeld_v_mbd_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_j)
                      end if                     
                      r0_ii_mbd(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_j)**(1.d0/3.d0)
                      xyz_mbd(:,k2) = xyz_j-xyz_i
                      rjs_mbd(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                      rjs_0_mbd(k2) = rjs(n_tot+k_j)/Bohr
                      xyz_0_mbd(:,k2) = xyz(:,n_tot+k_j)/Bohr
                      hirshfeld_mbd_neigh(k2) = hirshfeld_v_neigh(n_tot+k_j)
                      neighbor_alpha0_mbd(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_j)
                      if ( rjs(n_tot+k_j) .le. rcut_loc-r_buf_loc ) then
                        r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                        a_mbd(k2) = a_iso(r,2)
                        o_mbd(k2) = o_p(r)
                        r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
                      else if ( rjs(n_tot+k_j) .le. rcut_loc .and. rjs(n_tot+k_j) > rcut-r_buf_loc ) then
                        r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                        rb = (rjs(n_tot+k_j)-rcut_loc+r_buf_loc)/r_buf_loc
                        a_mbd(k2) = a_iso(r,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                        !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2) * &
                        central_pol(j1) * & 
                                    ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                        o_mbd(k2) = o_p(r) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                        central_omega(j1) * & 
                                    ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                        r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                r0_ii_mbd(k2) * (central_pol(j1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                                 !r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0) * &
                                ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                      else if ( rjs(n_tot+k_j) > rcut_loc .and. rjs(n_tot+k_j) .le. rcut_mbd) then
                        a_mbd(k2) = central_pol(j1) 
                                    !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2)
                        o_mbd(k2) = central_omega(j1)
                        r0_ii_SCS(k2) = r0_ii_mbd(k2)  * (central_pol(j1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0)
                                              ! * (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0)
                      end if
                      if ( rjs(n_tot+k_j) > rcut_mbd-r_buf_mbd &
                           .and. rjs(n_tot+k_j) .le. rcut_mbd ) then
                        rb = (rjs(n_tot+k_j)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                        T_LR_mult_0j(k2) = (1.d0 & 
                               - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5)**(1.d0/n_order)
                        !if ( rb < 0.5d0 ) then
                        !T_LR_mult_0j(k2) = (1.d0 - ( ac2*rb**2 + ac3*rb**3 + ac4*rb**4 ) )
                        !else
                        !  T_LR_mult_0i(k2) = (1.d0 -3.d0*rb**2 +2.d0*rb**3)**(1.d0/2.d0)
                        !end if

                      else if ( rjs(n_tot+k_j) > rcut_mbd ) then
                        T_LR_mult_0j(k2) = 0.d0
                      end if
                      T_LR_mult_j = T_LR_mult_0j(k2)
                      if ( p .ne. 1 .and. q .ne. 1 .and. rjs_mbd(k2)*Bohr > rcut_mbd-r_buf_mbd &
                          .and. rjs_mbd(k2)*Bohr .le. rcut_mbd ) then
                        rb = (rjs_mbd(k2)*Bohr-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                        T_LR_mult_ij(k2) = (1.d0 & 
                               - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5)**(1.d0/n_order)
                        !if ( rb < 0.5d0 ) then
                        !T_LR_mult_ij(k2) = (1.d0 - ( ac2*rb**2 + ac3*rb**3 + ac4*rb**4 ) )
                        !else
                        !  T_LR_mult_0i(k2) = (1.d0 -3.d0*rb**2 +2.d0*rb**3)**(1.d0/2.d0)
                        !end if

                      else if ( rjs_mbd(k2)*Bohr > rcut_mbd ) then
                        T_LR_mult_ij(k2) = 0.d0
                      end if
                      !if ( rjs_mbd(k2)*Bohr+rjs(n_tot+k_i) > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                      !     rjs_mbd(k2)*Bohr+rjs(n_tot+k_i) .le. 1.0d0*rcut_mbd ) then
                      !  rb = (rjs_mbd(k2)*Bohr+rjs(n_tot+k_i)-1.0d0*rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                      !  T_LR_mult_0ij(k2) = (1.d0 & 
                      !         - 10.d0 * rb**3 &
                      !         + 15.d0 * rb**4 &
                      !         - 6.d0 * rb**5) 
                      !else if ( rjs_mbd(k2)*Bohr+rjs(n_tot+k_i) > 1.0d0*rcut_mbd ) then
                      !  T_LR_mult_0ij(k2) = 0.d0
                      !end if
                      !if ( rjs_mbd(k2)*Bohr+rjs(n_tot+k_j) > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                      !     rjs_mbd(k2)*Bohr+rjs(n_tot+k_j) .le. 1.0d0*rcut_mbd ) then
                      !  rb = (rjs_mbd(k2)*Bohr+rjs(n_tot+k_j)-1.0d0*rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                      !  T_LR_mult_0ji(k2) = (1.d0 & 
                      !         - 10.d0 * rb**3 &
                      !         + 15.d0 * rb**4 &
                      !         - 6.d0 * rb**5)
                      !else if ( rjs_mbd(k2)*Bohr+rjs(n_tot+k_j) > 1.0d0*rcut_mbd ) then
                      !  T_LR_mult_0ji(k2) = 0.d0
                      !end if
                      !!r0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_mbd(k2)/neighbor_alpha0_mbd(k2)/hirshfeld_mbd_neigh(k2))**(1.d0/3.d0)
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
                          T_LR(3*(p-1)+c1,3*(q-1)+c2) = f_damp_SCS(k2) * T_mbd(k3) * T_LR_mult_i * T_LR_mult_j &
                                                        * T_LR_mult_ij(k2)
                          !write(*,*) "rcut_mbd, r_buf_mbd", rcut_mbd, r_buf_mbd
                        end do
                      end do
                    end if
                  end if
                end if
              end do
            end if
          end do

end if


          !write(*,*) "a_mbd", a_mbd
          !write(*,*) "o_mbd", o_mbd
        
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
          hirshfeld_2b_neigh = 0.d0
          C6_2b = 0.d0
          r0_ii_SCS_2b = 0.d0
          f_damp_SCS_2b = 0.d0
          r6_mult = 1.d0
          s = neighbor_species(n_tot+1)
          r_vdw_i = r0_ref(s) / Bohr * (a_iso(1,2)/(alpha0_ref(s)/Bohr**3))**(1.d0/3.d0)
          do i3 = 1, n_neigh(i)
            k_i = k_i+1
            i2 = neighbors_list(n_tot+k_i)
            i1 = modulo(i2-1, n_sites0) + 1
            if ( rjs(n_tot+k_i) .le. rcut_2b .and. rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd ) then
              k2 = k2+1
              s = neighbor_species(n_tot+k_i)
              sub_2b_list(k2) = neighbors_list(n_tot+k_i)
              !if ( do_derivatives .and. do_hirshfeld_gradients ) then
              !  hirshfeld_v_2b_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
              !end if
              r0_ii_2b(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
              xyz_2b(:,k2) = xyz(:,n_tot+k_i)/Bohr
              xyz_i = xyz_2b(:,k2)
              rjs_2b(k2) = rjs(n_tot+k_i)/Bohr
              neighbor_alpha0_2b(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_i)
              hirshfeld_2b_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
              if ( rjs(n_tot+k_i) .le. rcut_loc-r_buf_loc ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                a_2b(k2) = a_iso(r,2)
                o_2b(k2) = o_p(r)
                r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (a_iso(r,2)/neighbor_alpha0_2b(k2))**(1.d0/3.d0)
              else if ( rjs(n_tot+k_i) .le. rcut_loc .and. rjs(n_tot+k_i) > rcut_loc-r_buf_loc ) then
                r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                rb = (rjs(n_tot+k_i)-rcut_loc+r_buf_loc)/r_buf_loc
                a_2b(k2) = a_iso(r,2) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 ) + &
                            central_pol(i1) * &
                            !neighbor_alpha0_2b(k2) * hirshfeld_v_neigh(n_tot+k_i) * &
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                o_2b(k2) = o_p(r) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                            central_omega(i1) * &
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (a_iso(r,2)/neighbor_alpha0_2b(k2))**(1.d0/3.d0) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                r0_ii_2b(k2) * (central_pol(i1)/neighbor_alpha0_2b(k2))**(1.d0/3.d0) * &
                                !r0_ii_2b(k2) * (hirshfeld_2b_neigh(k2))**(1.d0/3.d0) * &
                                ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
              else if ( rjs(n_tot+k_i) > rcut_loc ) then
                a_2b(k2) = central_pol(i1) 
                            !neighbor_alpha0_2b(k2) * hirshfeld_v_neigh(n_tot+k_i)
                o_2b(k2) = central_omega(i1)
                r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (hirshfeld_2b_neigh(k2))**(1.d0/3.d0)
              end if
              if ( rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd .and. rjs(n_tot+k_i) .le. rcut_mbd ) then
                rb = (rjs(n_tot+k_i)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                r6_mult(k2) = ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
              else if ( rjs(n_tot+k_i) > rcut_2b-r_buf_2b .and. rjs(n_tot+k_i) .le. rcut_2b) then
                rb = (rjs(n_tot+k_i)-rcut_2b+r_buf_2b)/(r_buf_2b)
                r6_mult(k2) = (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5)
              end if
              C6_2b(k2) = 3.d0/2.d0 * a_iso(1,2) * a_2b(k2) * (o_p(1) * o_2b(k2)) / &
                      (o_p(1) + o_2b(k2))
              !r0_ii_SCS_2b(k2) = r0_ii_2b(k2) * (a_2b(k2)/neighbor_alpha0_2b(k2))**(1.d0/3.d0)
              r_vdw_j = r0_ii_SCS_2b(k2)
              f_damp_SCS_2b(k2) = 1.d0/( 1.d0 + exp( -20.d0*( rjs_2b(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
              E_TS = E_TS - C6_2b(k2)/rjs_2b(k2)**6 * f_damp_SCS_2b(k2) * r6_mult(k2)
            end if
          end do
          E_TS = 1.d0/2.d0 * E_TS
          !write(*,*) "E_TS", E_TS
          energies(i) = E_TS * Hartree

if ( .false. ) then

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

          !write(*,*) "n_mbd_sites", n_mbd_sites

          !if ( i == 1 .and. om == 2 ) then
          !  write(*,*) "AT"
          !  do p = 1, 3*n_mbd_sites
          !    write(*,*) AT(p,:,1)
          !  end do
          !end if
          
          !TEST!!!!!!!!!!!!!!!!!
          !allocate( V_int(1:3*n_mbd_sites,1:3*n_mbd_sites) )
          !V_int = 0.d0
          !k3 = 0
          !do p = 1, n_mbd_sites
          !  k3 = k3+1
          !  i2 = mbd_neighbors_list(k3)
          !  i1 = modulo(i2-1, n_sites0) + 1
          !  o_i = central_omega(i1)
          !  a_mbd_i = central_pol(i1)
          !  do c1 = 1, 3
          !    V_int(3*(p-1)+c1,3*(p-1)+c1) = o_i**2
          !  end do
          !  do j3 = 2, n_mbd_neigh(p)
          !    k3 = k3+1
          !    q = p_mbd(k3)
          !    j = mbd_neighbors_list(k3)
          !    j1 = modulo(j-1, n_sites0) + 1
          !    V_int(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3) = o_i * central_omega(j1) * &
          !          sqrt(a_mbd_i * central_pol(j1)) * T_LR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3)
          !  end do
          !end do
          !write(*,*) "omega_sum", sum(central_omega)
          !if ( i == 1 .and. om == 2 ) then
          !  open(unit=89, file="V_int.dat", status="new")
          !  write(*,*) "V_int"
          !  do p = 1, 3*n_mbd_sites
          !    write(89,*) V_int(p,:)
          !  end do
          !  close(89)
          !end if
          !
          !deallocate( V_int )
          !TEST!!!!!!!!!!!!!!!!!

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

          if ( do_nnls ) then

          allocate( A_nnls(1:n_freq,1:n_order+1) )
          allocate( b_nnls(1:n_freq) )
          allocate( coeff_nnls(1:n_order+1) )
          allocate( work_nnls(1:n_order+1) )
          allocate( ind_nnls(1:n_order+1) )
          allocate( omegas_nnls(1:201) )
          allocate( integrand_nnls(1:201) )

          A_nnls = 0.d0
          b_nnls = 0.d0
          coeff_nnls = 0.d0
          res_nnls = 0.d0
          work_nnls = 0.d0
          ind_nnls = 0.d0

          do i2 = 1, n_freq
            b_nnls(i2) = abs(integrand(i2))
            do j2 = 1, n_order+1
              if ( j2 == 1 ) then
                A_nnls(i2,j2) = 1.d0
              else
                A_nnls(i2,j2) = -abs(integrand(i2))*omegas_mbd(i2)**(2*(j2-1))
              end if
            end do
          end do

          call nnls(A_nnls, n_freq, n_order+1, b_nnls, coeff_nnls, res_nnls, work_nnls, ind_nnls, mode_nnls)

          integrand_nnls = 1.d0
          omegas_nnls = 0.d0
          do j2 = 2, n_order+1
            integrand_nnls(1) = integrand_nnls(1) + coeff_nnls(j2)*omegas_nnls(1)**(2.d0*(j2-1))
          end do
          integrand_nnls(1) = sign(coeff_nnls(1),integrand(1))/integrand_nnls(1)
          do i2 = 2, 201
            omegas_nnls(i2) = omegas_nnls(i2-1)+0.02d0
            do j2 = 2, n_order+1
              integrand_nnls(i2) = integrand_nnls(i2) + coeff_nnls(j2)*omegas_nnls(i2)**(2.d0*(j2-1)) 
            end do
            integrand_nnls(i2) = sign(coeff_nnls(1),integrand(1))/integrand_nnls(i2)
          end do

          integral = 0.d0
          call cpu_time(time1)
          call integrate("trapezoidal", omegas_nnls, integrand_nnls, omegas_nnls(1), omegas_nnls(size(omegas_nnls)), integral)
          call cpu_time(time2)
          write(*,*) "Integration time", time2-time1
          integral = integral/(2.d0*pi)
         
          deallocate( A_nnls, b_nnls, coeff_nnls, work_nnls, ind_nnls, omegas_nnls, integrand_nnls )

          else

          integral = 0.d0
          call integrate("trapezoidal", omegas_mbd, integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)
          integral = integral/(2.d0*pi)

          end if ! do_nnls

          energies(i) = (integral + E_TS) * Hartree
          write(*,*) "MBD energy", i, energies(i)
          E_MBD = E_MBD + energies(i)          
        
        end if ! om loop
        
end if
        
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

            !dB_mat = 0.d0
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
                  if ( rjs_0(k3) .le. rcut ) then
                    q = p_list(k3)
                  end if
                  sigma_ij = sqrt(s_i**2 + s_j**2)
                  f_damp_der(k3) = d/(sR*(r_vdw_i + r_vdw_j)) * f_damp(k3)**2 * &
                                     exp( -d*(rjs_H(k3)/(sR*(r_vdw_i + &
                                     r_vdw_j)) - 1.d0) ) * xyz_H(c3,k3)/rjs_H(k3)
                  g_func_der(k3) = 4.d0/sqrt(pi) * rjs_H(k3)/sigma_ij**3 * xyz_H(c3,k3) * &
                                       exp(-rjs_H(k3)**2/sigma_ij**2)
                  inner_damp_der = 0.d0
                  if ( rjs_H(k3)*Bohr < 2.d0 .and. rjs_0(k3) > rcut-r_buf_scs ) then ! This is probably never true but it's still here just in case.
                    rb = (rjs_H(k3)*Bohr)/2.d0                                       ! This is also material dependent! Change it at some point. Used to avoid singularities at short distances.
                    if ( a == i2 ) then
                      inner_damp_der = (30.d0*rb**2-60.d0*rb**3+30.d0*rb**4) * &
                                             (-xyz_H(c3,k3)/rjs_H(k3))*(Bohr/2.d0)
                    else
                      inner_damp_der = (30.d0*rb**2-60.d0*rb**3+30.d0*rb**4) * &
                                             (xyz_H(c3,k3)/rjs_H(k3))*(Bohr/2.d0)
                    end if
                  end if                  
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
                      if ( rjs_0(k3) .le. rcut ) then
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
                                                            neighbor_alpha0(k3) * inner_damp(k3) - d_arr_i(k4) * &
                                                            d_mult_i(k3) * inner_damp_der
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
                                                            neighbor_alpha0(k3) * inner_damp(k3) - d_arr_i(k4) * &
                                                            d_mult_i(k3) * inner_damp_der
                        end if      
                      else
                        if ( a == i2 ) then
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) + ((f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            neighbor_alpha0(k3)) * d_mult_o(k3) * inner_damp(k3) - &
                                                            d_arr_o(k4) * d_mult_o(k3) * &
                                                            inner_damp_der
                        else
                          d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - ((f_damp_der(k3) * T_func(k4) * &
                                                            g_func(k3) - (1.d0 - f_damp(k3)) * dT(k4) * &
                                                            g_func(k3) - g_func_der(k3) * (1.d0 - f_damp(k3)) * &
                                                            T_func(k4) - f_damp_der(k3) * h_func(k4) + &
                                                            h_func_der(k4) * (1.d0 - f_damp(k3))) * &
                                                            neighbor_alpha0(k3)) * d_mult_o(k3) * inner_damp(k3) - &
                                                            d_arr_o(k4) * d_mult_o(k3) * &
                                                            inner_damp_der
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
                  if ( rjs_0(k3) .le. rcut ) then
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
                hv_p_der = hirshfeld_v_sub_der(c3,k3+1)*Bohr
                r_vdw_i = r0_ii(k3+1)
                s_i = neighbor_sigma(k3+1)
                do c1 = 1, 3
                  b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) - 1.d0/(neighbor_alpha0(k3+1) * hirshfeld_sub_neigh(k3+1)) * &
                                        hirshfeld_v_sub_der(c3,k3+1)*Bohr * a_SCS(3*(p-1)+c1,:) 
                  dB_mat(3*(p-1)+c1,3*(p-1)+c1) = dB_mat(3*(p-1)+c1,3*(p-1)+c1) - 1.d0/(neighbor_alpha0(k3+1) &
                                        * hirshfeld_sub_neigh(k3+1)) * &
                                        hirshfeld_v_sub_der(c3,k3+1)*Bohr
                end do
                do j2 = 2, n_sub_neigh(p)
                  j = sub_neighbors_list(k3+j2)
                  r_vdw_j = r0_ii(k3+j2)
                  s_j = neighbor_sigma(k3+j2)
                  hv_q_der = hirshfeld_v_sub_der(c3,k3+j2)*Bohr
                  if ( rjs_0(k3+j2) .le. rcut ) then
                    q = p_list(k3+j2)
                  !else ! This needs to be changed for potential case rcut < rcut_soap
                  !  hv_q_der = 0.d0
                  end if
                  k4 = 9*(k3+j2-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k4 = k4+1
                      if ( rjs_0(k3+j2) .le. rcut ) then                      
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
                          neighbor_alpha0(k3+j2) * d_mult_i(k3+j2) * inner_damp(k3+j2) - &
                          d_arr_i(k4)/hirshfeld_sub_neigh(k3+j2) * hv_q_der * d_mult_i(k3+j2) * &
                          inner_damp(k3+j2)
                      else
                        d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - &
                          ((coeff_der(k4) * s_i**2/hirshfeld_sub_neigh(k3+1) + &
                          coeff_fdamp(k4) * r_vdw_i/hirshfeld_sub_neigh(k3+1)) * &
                          hv_p_der + &
                          (coeff_der(k4) * s_j**2/hirshfeld_sub_neigh(k3+j2) + &
                          coeff_fdamp(k4) * r_vdw_j/hirshfeld_sub_neigh(k3+j2)) * &
                          hv_q_der) * &
                          neighbor_alpha0(k3+j2) * d_mult_o(k3+j2) * inner_damp(k3+j2) - &
                          d_arr_o(k4)/hirshfeld_sub_neigh(k3+j2) * hv_q_der * d_mult_o(k3+j2) * &
                          inner_damp(k3+j2)
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
                if ( rjs_0(k3) .le. rcut ) then
                  q = p_list(k3)
                end if
                k4 = 9*(k3-1)
                do c1 = 1, 3
                  do c2 = 1, 3
                    k4 = k4+1
                    if ( rjs_0(k3) .le. rcut ) then
                      dB_mat(3*(p-1)+c1,3*(q-1)+c2) = dB_mat(3*(p-1)+c1,3*(q-1)+c2) + &
                          dT_SR_mult(k3,c3) * T_SR(k4)
                      b_der(3*(p-1)+c1,:) = b_der(3*(p-1)+c1,:) + dT_SR_mult(k3,c3) * T_SR(k4) * &
                                            a_SCS(3*(q-1)+c2,:) 
                      d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - d_dmult_i(k3,c3) * d_arr_i(k4) * &
                                                  inner_damp(k3)
                    else
                      d_der(3*(p-1)+c1,c2) = d_der(3*(p-1)+c1,c2) - d_dmult_o(k3,c3) * d_arr_o(k4) * &
                                                  inner_damp(k3)
                    end if
                  end do
                end do
              end do
            end do

            !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
             ! write(*,*) "d_der"
             ! do c1 = 1, 3
             !   write(*,*) d_der(3*(55-1)+c1,:)
             ! end do
             ! write(*,*) "dB_mat"
             ! do p = 1, 3*n_sub_sites
             !   write(*,*) dB_mat(p,:)
             ! end do
            !  write(*,*) "n_sub_sites", n_sub_sites
            !  write(*,*) "a_SCS"
            !  do p = 1, 3*n_sub_sites
            !    write(*,*) a_SCS(p,:)
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

            if ( om == 2 .and. c3 == 3 .and. i == 1 ) then            
            write(*,*) "da_iso", i
            do p = 1, n_sub_sites
              write(*,*) da_iso(p,:,2)
            end do
            end if

            !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
            !  write(*,*) "da_iso"
            !  !write(*,*) "da_SCS"
            !  do p = 1, n_sub_sites
            !    write(*,*) da_iso(p,1,2)
            !    !write(*,*) da_SCS(p,:)
            !  end do
            !end if



            if ( om == 2 ) then

if ( .false. ) then
            
              f_damp_der_SCS = 0.d0
              f_damp_der_mbd = 0.d0 ! This is cleared so we can recalculate it with SCS values
              dT_mbd = 0.d0
              dT_LR = 0.d0
              dT_LR_mult_0i = 0.d0
              dT_LR_mult_0j = 0.d0
              dT_LR_mult_ij = 0.d0
              da_mbd = 0.d0
              do_mbd = 0.d0
              dr0_ii_SCS = 0.d0
              k2 = 0
              do p = 1, n_mbd_sites
                k2 = k2+1
                r_vdw_i = r0_ii_SCS(k2)
                i2 = mbd_neighbors_list(k2)
                i1 = modulo(i2-1, n_sites0) + 1
                rjs_i = rjs_0_mbd(k2)
                xyz_i = xyz_0_mbd(:,k2)
                if ( rjs_0_mbd(k2) .le. (rcut_loc-r_buf_loc)/Bohr ) then
                  r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                  da_mbd(k2) = da_iso(r,c3,2)
                  if ( a_iso(r,2) > a_iso(r,1) ) then
                    do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                                 ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                  else
                    write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                    do_pref = 0.d0
                  end if
                  do_mbd(k2) = do_pref
                  dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2)
                else if ( rjs_0_mbd(k2) > (rcut_loc-r_buf_loc)/Bohr .and. rjs_0_mbd(k2) .le. rcut_loc/Bohr ) then
                  r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                  if ( a_iso(r,2) > a_iso(r,1) ) then
                    do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                                 ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                  else
                    write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                    do_pref = 0.d0
                  end if
                  rb = (rjs_0_mbd(k2)*Bohr-rcut_loc+r_buf_loc)/r_buf_loc
                  da_mbd(k2) = da_iso(r,c3,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                    a_iso(r,2) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                              central_pol(i1) * & 
                               !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2) * &
                               ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                  dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                   r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                                r0_ii_mbd(k2) &
                                * (central_pol(i1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * &
                                !* (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0) * &
                              (+30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                  if ( do_hirshfeld_gradients ) then
                    da_mbd(k2) = da_mbd(k2) + neighbor_alpha0_mbd(k2) * &
                                    hirshfeld_v_mbd_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                    dr0_ii_SCS(k2) = dr0_ii_SCS(k2) + &
                               r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                               hirshfeld_v_mbd_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)       
                  end if
                  do_mbd(k2) = do_pref * & 
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                    o_p(r) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                               central_omega(i1) * &
                               ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                else if ( rjs_0_mbd(k2) > rcut_loc/Bohr .and. rjs_0_mbd(k2) .le. rcut_mbd/Bohr ) then
                  if ( do_hirshfeld_gradients ) then
                    da_mbd(k2) = neighbor_alpha0_mbd(k2) * &
                                    hirshfeld_v_mbd_der(c3,k2)*Bohr
                    dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                               hirshfeld_v_mbd_der(c3,k2)*Bohr
                  else
                    da_mbd(k2) = 0.d0
                    dr0_ii_SCS(k2) = 0.d0
                  end if
                  do_mbd(k2) = 0.d0
                end if
                if ( rjs_0_mbd(k2) > (rcut_mbd-r_buf_mbd)/Bohr .and. &
                     rjs_0_mbd(k2) .le. rcut_mbd/Bohr ) then
                  rb = (rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buf_mbd)/r_buf_mbd
                  dT_LR_mult_0i(k2) = (1.d0/n_order) * ( 1.d0 - 10.d0*rb**3 &
                               +15.d0 * rb**4 &
                               - 6.d0 * rb**5 )**(1.d0/n_order-1.d0) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                              * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_mbd/Bohr))
                end if
                T_LR_mult_i = T_LR_mult_0i(k2)
                dT_LR_mult_i = dT_LR_mult_0i(k2)
                dr_vdw_i = dr0_ii_SCS(k2)
                a_mbd_i = a_mbd(k2)
                do j3 = 2, n_mbd_neigh(p)
                  k2 = k2+1
                  j = mbd_neighbors_list(k2)
                  j1 = modulo(j-1, n_sites0) + 1
                  q = p_mbd(k2)
                  rjs_j = rjs_0_mbd(k2)
                  xyz_j = xyz_0_mbd(:,k2)
                  if ( rjs_0_mbd(k2) .le. (rcut_loc-r_buf_loc)/Bohr ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                    da_mbd(k2) = da_iso(r,c3,2)
                    if ( a_iso(r,2) > a_iso(r,1) ) then
                      do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                                 ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                    else
                      write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                      do_pref = 0.d0
                    end if
                    dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2)
                  else if ( rjs_0_mbd(k2) > (rcut_loc-r_buf_loc)/Bohr .and. rjs_0_mbd(k2) .le. rcut_loc/Bohr ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j,1)
                    if ( a_iso(r,2) > a_iso(r,1) ) then
                      do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) - a_iso(r,2) * da_iso(r,c3,1)) / &
                                 ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                    else
                      write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                      do_pref = 0.d0
                    end if
                    rb = (rjs_0_mbd(k2)*Bohr-rcut_loc+r_buf_loc)/r_buf_loc
                    da_mbd(k2) = da_iso(r,c3,2) * &
                                (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                      a_iso(r,2) * &
                                (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                                 central_pol(j1) * & 
                                  !neighbor_alpha0_mbd(k2) * hirshfeld_mbd_neigh(k2) * &
                                  ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                   * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                    dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                   r0_ii_mbd(k2) * (a_iso(r,2)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) *  &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                                r0_ii_mbd(k2) &
                                * (central_pol(j1)/neighbor_alpha0_mbd(k2))**(1.d0/3.d0) * & 
                                !* (hirshfeld_mbd_neigh(k2))**(1.d0/3.d0) * &
                              (+30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                               * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                    if ( do_hirshfeld_gradients ) then
                      da_mbd(k2) = da_mbd(k2) + neighbor_alpha0_mbd(k2) * &
                                    hirshfeld_v_mbd_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                      dr0_ii_SCS(k2) = dr0_ii_SCS(k2) + &
                               r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                               hirshfeld_v_mbd_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5) 
                    end if 
                    do_mbd(k2) = do_pref * &
                                (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                      o_p(r) * &
                                (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                                 * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr)) + &
                                 central_omega(j1) * &
                                  ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                   * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_loc/Bohr))
                  else if ( rjs_0_mbd(k2) > rcut_loc/Bohr .and. rjs_0_mbd(k2) .le. rcut_mbd/Bohr ) then
                    if ( do_hirshfeld_gradients ) then
                      da_mbd(k2) = neighbor_alpha0_mbd(k2) * &
                                    hirshfeld_v_mbd_der(c3,k2)*Bohr
                      dr0_ii_SCS(k2) = r0_ii_mbd(k2) * (hirshfeld_mbd_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                               hirshfeld_v_mbd_der(c3,k2)*Bohr
                    else
                      da_mbd(k2) = 0.d0
                      dr0_ii_SCS(k2) = 0.d0
                    end if
                    do_mbd(k2) = 0.d0
                  end if
                  if ( rjs_0_mbd(k2) > (rcut_mbd-r_buf_mbd)/Bohr .and. &
                       rjs_0_mbd(k2) .le. rcut_mbd/Bohr ) then
                    rb = (rjs_0_mbd(k2)*Bohr-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                    dT_LR_mult_0j(k2) = (1.d0/n_order) * ( 1.d0 - 10.d0*rb**3 &
                               +15.d0 * rb**4 &
                               - 6.d0 * rb**5 )**(1.d0/n_order-1.d0) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                              * ( -xyz_0_mbd(c3,k2)/rjs_0_mbd(k2)/(r_buf_mbd/Bohr))
                  end if
                  T_LR_mult_j = T_LR_mult_0j(k2)
                  dT_LR_mult_j = dT_LR_mult_0j(k2)
                  !if ( rjs_mbd(k2)*Bohr+rjs_i > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                  !     rjs_mbd(k2)*Bohr+rjs_i .le. 1.0d0*rcut_mbd ) then
                  !  rb = (rjs_mbd(k2)*Bohr+rjs_i-1.0d0*rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                  !  dT_LR_mult_0ij(k2) = 1.d0/3.d0 * (1.d0 - 10.d0 * rb**3 &  
                  !             + 15.d0 * rb**4 &
                  !             - 6.d0 * rb**5)**(-2.d0/3.d0) * &
                  !            (- 30.d0 * rb**2 &
                  !             + 60.d0 * rb**3 &
                  !             - 30.d0 * rb**4 ) &
                  !             * ( -xyz_i(c3)/rjs_i/(r_buf_mbd/Bohr))
                  !end if
                  !if ( rjs_mbd(k2)*Bohr+rjs_j > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                  !         rjs_mbd(k2)*Bohr+rjs_j .le. 1.0d0*rcut_mbd ) then
                  !  rb = (rjs_mbd(k2)*Bohr+rjs_j-1.0d0*rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                  !  dT_LR_mult_0ji(k2) = (- 30.d0 * rb**2 &
                  !             + 60.d0 * rb**3 &
                  !             - 30.d0 * rb**4) &
                  !             * ( -xyz_j(c3)/rjs_j/(r_buf_mbd/Bohr))
                  !end if
                  dr_vdw_j = dr0_ii_SCS(k2)
                  a_mbd_j = a_mbd(k2)
                  r_vdw_j = r0_ii_SCS(k2)
                  R_vdW_SCS_ij = r_vdw_i + r_vdw_j
                  S_vdW_ij = sR*R_vdW_SCS_ij
                  dS_vdW_ij = sR * ( dr_vdw_i + dr_vdw_j )
                  f_damp_der_SCS(k2) = -(d*rjs_mbd(k2))/S_vdW_ij**2 * f_damp_SCS(k2)**2 * &
                                         exp(-d*(rjs_mbd(k2)/S_vdW_ij - 1.d0)) * dS_vdW_ij
                  k3 = 9*(k2-1)
                  do c1 = 1, 3
                    do c2 = 1, 3
                      k3 = k3+1
                      dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                        T_mbd(k3) * f_damp_der_SCS(k2) * T_LR_mult_i * &
                                        T_LR_mult_j * T_LR_mult_ij(k2) + &
                                        f_damp_SCS(k2) * T_mbd(k3) * T_LR_mult_ij(k2) * &
                                        T_LR_mult_j * dT_LR_mult_i + &
                                        f_damp_SCS(k2) * T_mbd(k3) * T_LR_mult_ij(k2) * &
                                        T_LR_mult_i * dT_LR_mult_j
                    end do
                  end do
                  if (i == i2 .or. i == j) then
                    f_damp_der_mbd(k2) = d/S_vdW_ij * f_damp_SCS(k2)**2 * &
                                     exp( -d*(rjs_mbd(k2)/S_vdW_ij - 1.d0) ) * xyz_mbd(c3,k2)/rjs_mbd(k2)
                    if ( p .ne. 1 .and. q .ne. 1 .and. rjs_mbd(k2) > (rcut_mbd-r_buf_mbd)/Bohr &
                        .and. rjs_mbd(k2) .le. rcut_mbd/Bohr ) then
                      rb = (rjs_mbd(k2)*Bohr-rcut_mbd+r_buf_mbd)/r_buf_mbd
                      if ( i == i2 ) then
                        dT_LR_mult_ij(k2) = (1.d0/n_order) * ( 1.d0 - 10.d0*rb**3 &
                               +15.d0 * rb**4 &
                               - 6.d0 * rb**5 )**(1.d0/n_order-1.d0) * &
                             (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                              * ( -xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                      else if ( i == j ) then
                        dT_LR_mult_ij(k2) = (1.d0/n_order) * ( 1.d0 - 10.d0*rb**3 &
                               +15.d0 * rb**4 &
                               - 6.d0 * rb**5 )**(1.d0/n_order-1.d0) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                              * ( xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                      end if
                    end if
                    !if ( rjs_mbd(k2)*Bohr+rjs_i > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                    !     rjs_mbd(k2)*Bohr+rjs_i .le. 1.0d0*rcut_mbd ) then
                    !  rb = (rjs_mbd(k2)*Bohr+rjs_i-1.0d0*rcut_mbd+r_buf_mbd)/r_buf_mbd
                    !  if ( i == i2 ) then
                    !    dT_LR_mult_0ij_2 = (-30.d0 * rb**2 &
                    !          + 60.d0 * rb**3 &
                    !          - 30.d0 * rb**4) &
                    !          * ( -xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                    !    dT_LR_mult_0ij(k2) = dT_LR_mult_0ij(k2) + dT_LR_mult_0ij_2
                    !  else if ( i == j ) then
                    !    dT_LR_mult_0ij_2 = (-30.d0 * rb**2 &
                    !          + 60.d0 * rb**3 &
                    !          - 30.d0 * rb**4) &
                    !          * ( xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                    !    dT_LR_mult_0ji(k2) = dT_LR_mult_0ji(k2) + dT_LR_mult_0ji_2
                    !  end if
                    !end if
                    !if ( rjs_mbd(k2)*Bohr+rjs_j > 1.0d0*rcut_mbd-r_buf_mbd .and. &
                    !     rjs_mbd(k2)*Bohr+rjs_j .le. 1.0d0*rcut_mbd ) then
                    !  rb = (rjs_mbd(k2)*Bohr+rjs_j-1.0d0*rcut_mbd+r_buf_mbd)/r_buf_mbd
                    !  if ( i == i2 ) then
                    !    dT_LR_mult_0ij_2 = (-30.d0 * rb**2 &
                    !          + 60.d0 * rb**3 &
                    !          - 30.d0 * rb**4) &
                    !          * ( -xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                    !  else if ( i == j ) then
                    !    dT_LR_mult_0ji_2 = (-30.d0 * rb**2 &
                    !          + 60.d0 * rb**3 &
                    !          - 30.d0 * rb**4) &
                    !          * ( xyz_mbd(c3,k2)/rjs_mbd(k2)/(r_buf_mbd/Bohr))
                    !  end if
                    !end if
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
                                             (T_mbd(k3) * f_damp_der_mbd(k2) + &
                                             dT_mbd(k3) * f_damp_SCS(k2)) * &
                                             (T_LR_mult_i * T_LR_mult_j * T_LR_mult_ij(k2))
                        else if (i == j) then
                          dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                             (T_mbd(k3) * f_damp_der_mbd(k2) + &
                                             dT_mbd(k3) * f_damp_SCS(k2)) * &
                                             (T_LR_mult_i * T_LR_mult_j * T_LR_mult_ij(k2))
                        end if
                        dT_LR(3*(p-1)+c1,3*(q-1)+c2) = dT_LR(3*(p-1)+c1,3*(q-1)+c2) + &
                                             f_damp_SCS(k2) * T_mbd(k3) * T_LR_mult_i * &
                                             T_LR_mult_j * dT_LR_mult_ij(k2)
                      end do
                    end do
                  end if
                end do
              end do

              !if ( i == 1 .and. c3 == 1 ) then
              !  write(*,*) "dT_LR"
              !  do p = 1, 3*n_mbd_sites
              !    write(*,*) dT_LR(p,:)
              !  end do
              !end if
            
              !if ( i == 1 .and. c3 == 1 ) then
              !write(*,*) "r0_ii_SCS"
              !do k3 = 1, n_mbd_pairs
              !  write(*,*) k3, r0_ii_SCS(k3)
              !end do
              !write(*,*) "dr0_ii_SCS"
              !do k3 = 1, n_mbd_pairs
              !  write(*,*) k3, dr0_ii_SCS(k3)
              !end do
              !end if

              !write(*,*) "da_mbd", da_mbd
              !write(*,*) "do_mbd", do_mbd
              
end if

              n_tot = sum(n_neigh(1:i))-n_neigh(i)
              n_2b_tot_sites = 0
              n_2b_tot_pairs = 0
              k_i = 0
              do i3 = 1, n_neigh(i)
                k_i = k_i + 1
                if (rjs(n_tot+k_i) .le. rcut_2b ) then
                  n_2b_tot_sites = n_2b_tot_sites + 1
                  n_2b_tot_pairs = n_2b_tot_pairs + 1
                  xyz_i = xyz(:,n_tot+k_i)/Bohr
                  k_j = 0
                  do j3 = 1, n_neigh(i)
                    k_j = k_j + 1
                    if ( rjs(n_tot+k_i) .le. rcut_mbd-r_buf_mbd) then
                      !if ( rjs(n_tot+k_j) .le. rcut_2b .and. rjs(n_tot+k_j) > rcut_mbd-r_buf_mbd ) then
                        if (i3 .ne. j3) then
                          xyz_j = xyz(:,n_tot+k_j)/Bohr
                          if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut_2b/Bohr ) then
                            n_2b_tot_pairs = n_2b_tot_pairs + 1
                          end if
                        end if
                      !end if
                    else
                      !if ( rjs(n_tot+k_j) .le. rcut_2b ) then
                        if (i3 .ne. j3) then
                          xyz_j = xyz(:,n_tot+k_j)/Bohr
                          if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut_2b/Bohr ) then
                            n_2b_tot_pairs = n_2b_tot_pairs + 1
                          end if
                        end if
                      !end if
                    end if
                  end do
                end if
              end do

              allocate( sub_2b_tot_list(1:n_2b_tot_pairs) )
              allocate( n_2b_tot_neigh(1:n_2b_tot_sites) )
              allocate( p_2b_tot(1:n_2b_tot_pairs) )
              allocate( xyz_2b_tot(1:3,1:n_2b_tot_pairs) )
              allocate( rjs_2b_tot(1:n_2b_tot_pairs) )
              allocate( r0_ii_2b_tot(1:n_2b_tot_pairs) )
              allocate( neighbor_alpha0_2b_tot(1:n_2b_tot_pairs) )
              allocate( f_damp_SCS_2b_tot(1:n_2b_tot_pairs) )
              allocate( hirshfeld_2b_tot_neigh(1:n_2b_tot_pairs) )
              allocate( hirshfeld_v_2b_der(1:3,n_2b_tot_pairs) )
              allocate( a_2b_tot(1:n_2b_tot_pairs) )
              allocate( o_2b_tot(1:n_2b_tot_pairs) )
              allocate( r0_ii_SCS_2b_tot(1:n_2b_tot_pairs) )
              allocate( c6_2b_tot(1:n_2b_tot_pairs) )
              allocate( r6_mult_2b_tot(1:n_2b_tot_pairs) )
              allocate( dr6_mult(1:n_2b_tot_pairs) )
              allocate( da_2b(1:n_2b_tot_pairs) )
              allocate( do_2b(1:n_2b_tot_pairs) )
              allocate( dr0_ii_SCS_2b(1:n_2b_tot_pairs) )
              allocate( r6_mult_0i(1:n_2b_tot_pairs) )
              allocate( r6_mult_0j(1:n_2b_tot_pairs) )
              allocate( dr6_mult_0i(1:n_2b_tot_pairs) )
              allocate( dr6_mult_0j(1:n_2b_tot_pairs) )

              E_TS_tot = 0.d0
              forces_TS = 0.d0
              
              sub_2b_tot_list = 0
              n_2b_tot_neigh = 0
              p_2b_tot = 0
              xyz_2b_tot = 0.d0
              rjs_2b_tot = 0.d0
              r0_ii_2b_tot = 0.d0
              neighbor_alpha0_2b_tot = 0.d0
              f_damp_SCS_2b_tot = 0.d0
              hirshfeld_2b_tot_neigh = 0.d0
              hirshfeld_v_2b_der = 0.d0
              a_2b_tot = 0.d0
              o_2b_tot = 0.d0
              r0_ii_SCS_2b_tot = 0.d0
              c6_2b_tot = 0.d0
              r6_mult_2b_tot = 1.d0
              dr6_mult = 0.d0
              da_2b = 0.d0
              do_2b = 0.d0
              dr0_ii_SCS_2b = 0.d0
              r6_mult_0i = 1.d0
              r6_mult_0j = 1.d0
              dr6_mult_0i = 0.d0
              dr6_mult_0j = 0.d0

        

              k2 = 0
              k_i = 0     
              do i3 = 1, n_neigh(i)
                k_i = k_i+1
                i2 = neighbors_list(n_tot+k_i)
                i1 = modulo(i2-1, n_sites0) + 1
                if ( rjs(n_tot+k_i) .le. rcut_2b ) then
                  k2 = k2+1
                  s = neighbor_species(n_tot+k_i)
                  sub_2b_tot_list(k2) = neighbors_list(n_tot+k_i)
                  hirshfeld_v_2b_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_i)
                  r0_ii_2b_tot(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
                  xyz_2b_tot(:,k2) = xyz(:,n_tot+k_i)/Bohr
                  xyz_i = xyz_2b_tot(:,k2)
                  rjs_2b_tot(k2) = rjs(n_tot+k_i)/Bohr
                  neighbor_alpha0_2b_tot(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_i)
                  hirshfeld_2b_tot_neigh(k2) = hirshfeld_v_neigh(n_tot+k_i)
                  if ( rjs(n_tot+k_i) .le. rcut_loc-r_buf_loc ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                    a_2b_tot(k2) = a_iso(r,2)
                    o_2b_tot(k2) = o_p(r)
                    r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0)
                    da_2b(k2) = da_iso(r,c3,2)
                    if ( a_iso(r,2) > a_iso(r,1) ) then
                      do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                  - a_iso(r,2) * da_iso(r,c3,1)) / &
                                  ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                    else
                      write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                      do_pref = 0.d0
                    end if
                    do_2b(k2) = do_pref
                    dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2)
                  else if ( rjs(n_tot+k_i) .le. rcut_loc .and. rjs(n_tot+k_i) > rcut_loc-r_buf_loc ) then
                    r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),i2,1)
                    rb = (rjs(n_tot+k_i)-rcut_loc+r_buf_loc)/r_buf_loc
                    a_2b_tot(k2) = a_iso(r,2) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5 ) + &
                            central_pol(i1) * &
                            !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                    o_2b_tot(k2) = o_p(r) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                            central_omega(i1) * &
                                  ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                    r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                            (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                r0_ii_2b_tot(k2) * (central_pol(i1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                !r0_ii_2b_tot(k2) * (hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                                ( + 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                    if ( a_iso(r,2) > a_iso(r,1) ) then
                      do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                 - a_iso(r,2) * da_iso(r,c3,1)) / &
                                 ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                    else
                      write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                      do_pref = 0.d0
                    end if
                    da_2b(k2) = da_iso(r,c3,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                    a_iso(r,2) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr)) + &
                              central_pol(i1) * &
                              !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                               ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr))
                    dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2) * &
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                   r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr)) + &
                                r0_ii_2b_tot(k2) &
                                * (central_pol(i1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                !* (hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                              (+30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                               * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr))
                    if ( do_hirshfeld_gradients ) then
                      da_2b(k2) = da_2b(k2) + neighbor_alpha0_2b_tot(k2) * &
                                    hirshfeld_v_2b_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)
                      dr0_ii_SCS_2b(k2) = dr0_ii_SCS_2b(k2) + &
                               r0_ii_2b_tot(k2) * (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                               hirshfeld_v_2b_der(c3,k2)*Bohr * &
                              (+ 10.d0 *rb**3 &
                               - 15.d0 * rb**4 &
                               + 6.d0 * rb**5)       
                    end if
                    do_2b(k2) = do_pref * & 
                              (1.d0 - 10.d0 * rb**3 &
                               + 15.d0 * rb**4 &
                               - 6.d0 * rb**5) + &
                                    o_p(r) * &
                              (-30.d0 * rb**2 &
                              + 60.d0 * rb**3 &
                              - 30.d0 * rb**4) &
                               * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr)) + &
                               central_omega(i1) * &
                               ( +30.d0 * rb**2 &
                              - 60.d0 * rb**3 &
                              + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_loc/Bohr))
                  else if ( rjs(n_tot+k_i) > rcut_loc ) then
                    a_2b_tot(k2) = &
                          central_pol(i1)
                          !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2)
                    o_2b_tot(k2) = central_omega(i1)
                    r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * &
                            (central_pol(i1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) 
                            !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0)
                    if ( do_hirshfeld_gradients ) then
                      da_2b(k2) = neighbor_alpha0_2b_tot(k2) * hirshfeld_v_2b_der(c3,k2)*Bohr
                      dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * &
                              (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * hirshfeld_v_2b_der(c3,k2)*Bohr
                    else
                      da_2b(k2) = 0.d0
                      dr0_ii_SCS_2b(k2) = 0.d0
                    end if
                    do_2b(k2) = 0.d0
                  end if
                  !if ( rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd .and. rjs(n_tot+k_i) .le. rcut_mbd ) then
                  !  rb = (rjs(n_tot+k_i)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                  !  r6_mult_2b_tot(k2) = ( + 10.d0 *rb**3 &
                  !                - 15.d0 * rb**4 &
                  !                + 6.d0 * rb**5)
                  !else if ( rjs(n_tot+k_i) > rcut_2b-r_buf_2b .and. rjs(n_tot+k_i) .le. rcut_2b) then
                  !  rb = (rjs(n_tot+k_i)-rcut_2b+r_buf_2b)/(r_buf_2b)
                  !  r6_mult_2b_tot(k2) = (1.d0 - 10.d0 * rb**3 &
                  !                 + 15.d0 * rb**4 &
                  !                 - 6.d0 * rb**5)
                  !end if
                  r_vdw_i = r0_ii_SCS_2b_tot(k2)
                  a_i = a_2b_tot(k2)
                  o_i = o_2b_tot(k2)
                  da_i = da_2b(k2)
                  do_i = do_2b(k2)
                  dr_vdw_i = dr0_ii_SCS_2b(k2)
                  k_j = 0
                  do j3 = 1, n_neigh(i)
                    k_j = k_j + 1
                    j2 = neighbors_list(n_tot+k_j)
                    j1 = modulo(j2-1, n_sites0) + 1
                    if (.false.) then
                    !if ( rjs(n_tot+k_i) .le. rcut_mbd .and. rjs(n_tot+k_j) .le. rcut_2b &
                    !     .and. rjs(n_tot+k_j) > rcut_mbd-r_buf_mbd ) then
                        if (i3 .ne. j3) then
                          xyz_j = xyz(:,n_tot+k_j)/Bohr
                          if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut_2b/Bohr ) then
                            k2 = k2+1
                            s = neighbor_species(n_tot+k_j)
                            sub_2b_tot_list(k2) = neighbors_list(n_tot+k_j)
                            hirshfeld_v_2b_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_j)
                            r0_ii_2b_tot(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
                            xyz_2b_tot(:,k2) = xyz_j-xyz_i
                            rjs_2b_tot(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                            neighbor_alpha0_2b_tot(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_i)
                            hirshfeld_2b_tot_neigh(k2) = hirshfeld_v_neigh(n_tot+k_j)
                            if ( rjs(n_tot+k_j) .le. rcut_loc-r_buf_loc ) then
                              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j2,1)
                              a_2b_tot(k2) = a_iso(r,2)
                              o_2b_tot(k2) = o_p(r)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0)
                              da_2b(k2) = da_iso(r,c3,2)
                             if ( a_iso(r,2) > a_iso(r,1) ) then
                               do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                  - a_iso(r,2) * da_iso(r,c3,1)) / &
                                  ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                             else
                               write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                               do_pref = 0.d0
                             end if
                             do_2b(k2) = do_pref
                             dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2)
                            else if ( rjs(n_tot+k_j) .le. rcut_loc .and. rjs(n_tot+k_j) > rcut_loc-r_buf_loc ) then
                              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j2,1)
                              rb = (rjs(n_tot+k_j)-rcut_loc+r_buf_loc)/r_buf_loc
                              a_2b_tot(k2) = a_iso(r,2) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                  + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5 ) + &
                                   central_pol(j1) * &
                                   !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                                     ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              o_2b_tot(k2) = o_p(r) * &
                                  (1.d0 - 10.d0 * rb**3 &
                                    + 15.d0 * rb**4 &
                                    - 6.d0 * rb**5) + &
                                  central_omega(j1) * &
                                     ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                  + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5) + &
                                   r0_ii_2b_tot(k2) * &
                                   (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                    !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                                   ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              if ( a_iso(r,2) > a_iso(r,1) ) then
                                do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                     - a_iso(r,2) * da_iso(r,c3,1)) / &
                                     ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                              else
                                write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                                do_pref = 0.d0
                              end if
                              da_2b(k2) = da_iso(r,c3,2) * &
                                  (1.d0 - 10.d0 * rb**3 &
                                   + 15.d0 * rb**4 &
                                   - 6.d0 * rb**5) + &
                                    a_iso(r,2) * &
                                 (-30.d0 * rb**2 &
                                 + 60.d0 * rb**3 &
                                 - 30.d0 * rb**4) &
                                  * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                 central_pol(j1) * & 
                                 !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                                 ( +30.d0 * rb**2 &
                                - 60.d0 * rb**3 &
                                + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                              dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                    + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5) + &
                                   r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                   (-30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                   r0_ii_2b_tot(k2) * &
                                   (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                   !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                                  (+30.d0 * rb**2 &
                                  - 60.d0 * rb**3 &
                                  + 30.d0 * rb**4) &
                                  * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                              if ( do_hirshfeld_gradients ) then
                                da_2b(k2) = da_2b(k2) + neighbor_alpha0_2b_tot(k2) * &
                                    hirshfeld_v_2b_der(c3,k2)*Bohr * &
                                    (+ 10.d0 *rb**3 &
                                    - 15.d0 * rb**4 &
                                    + 6.d0 * rb**5)
                                dr0_ii_SCS_2b(k2) = dr0_ii_SCS_2b(k2) + &
                                   r0_ii_2b_tot(k2) * (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                                   hirshfeld_v_2b_der(c3,k2)*Bohr * &
                                   (+ 10.d0 *rb**3 &
                                   - 15.d0 * rb**4 &
                                   + 6.d0 * rb**5)       
                              end if
                              do_2b(k2) = do_pref * & 
                                  (1.d0 - 10.d0 * rb**3 &
                                   + 15.d0 * rb**4 &
                                   - 6.d0 * rb**5) + &
                                    o_p(r) * &
                                  (-30.d0 * rb**2 &
                                  + 60.d0 * rb**3 &
                                  - 30.d0 * rb**4) &
                                   * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                  central_omega(j1) * &
                                  ( +30.d0 * rb**2 &
                                - 60.d0 * rb**3 &
                                 + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                            else if ( rjs(n_tot+k_j) > rcut_loc ) then
                              a_2b_tot(k2) = central_pol(j1) 
                                    !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2)
                              o_2b_tot(k2) = central_omega(j1)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0)
                                           !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0)
                              if ( do_hirshfeld_gradients ) then
                                da_2b(k2) = neighbor_alpha0_2b_tot(k2) * hirshfeld_v_2b_der(c3,k2)*Bohr
                                dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * &
                                     (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * hirshfeld_v_2b_der(c3,k2)*Bohr
                              else
                                da_2b(k2) = 0.d0
                                dr0_ii_SCS_2b(k2) = 0.d0
                              end if
                              do_2b(k2) = 0.d0
                            end if
                            if ( rjs_2b_tot(k2)*Bohr > rcut_2b-r_buf_2b .and. rjs_2b_tot(k2)*Bohr .le. rcut_2b) then
                              rb = (rjs_2b_tot(k2)*Bohr-rcut_2b+r_buf_2b)/(r_buf_2b)
                              r6_mult_2b_tot(k2) = (1.d0 - 10.d0 * rb**3 &
                                      + 15.d0 * rb**4 &
                                      - 6.d0 * rb**5)
                              if ( i2 == i ) then 
                                dr6_mult(k2) = (- 30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * ( -xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)/(r_buf_2b/Bohr))
                              else if ( j2 == i ) then
                                dr6_mult(k2) = (- 30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * (xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)/(r_buf_2b/Bohr))
                              end if
                            end if
                            if ( rjs(n_tot+k_j) .le. rcut_mbd .and. rjs(n_tot+k_j) > rcut_mbd-r_buf_mbd ) then
                              rb = (rjs(n_tot+k_j)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                              r6_mult_0j(k2) = (10.d0 * rb**3 &
                                      - 15.d0 * rb**4 &
                                      + 6.d0 * rb**5)
                              dr6_mult_0j(k2) = (30.d0 * rb**2 &
                                      - 60.d0 * rb**3 &
                                      + 30.d0 * rb**4) * &
                                      (-xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_mbd/Bohr))                                      
                            else if ( rjs(n_tot+k_j) .le. rcut_2b .and. rjs(n_tot+k_j) > rcut_2b-r_buf_2b ) then
                              rb = (rjs(n_tot+k_j)-rcut_2b+r_buf_2b)/(r_buf_2b)
                              r6_mult_0j(k2) = (1.d0 - 10.d0 * rb**3 &
                                      + 15.d0 * rb**4 &
                                      - 6.d0 * rb**5)
                              !dr6_mult_0j(k2) = (-30.d0 * rb**2 &
                              !        + 60.d0 * rb**3 &
                              !        - 30.d0 * rb**4) * &
                              !        (-xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_2b/Bohr))
                            end if
                            if ( rjs(n_tot+k_i) .le. rcut_mbd .and. rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd ) then
                              rb = (rjs(n_tot+k_i)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                              r6_mult_0i(k2) = (10.d0 * rb**3 &
                                      - 15.d0 * rb**4 &
                                      + 6.d0 * rb**5)
                              dr6_mult_0i(k2) = (30.d0 * rb**2 &
                                      - 60.d0 * rb**3 &
                                      + 30.d0 * rb**4) * &
                                      (-xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_mbd/Bohr))
                              !if ( i2 == i ) then
                              !  write(*,*) "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
                              !  dr6_mult_0i(k2) = dr6_mult_0i(k2) + (30.d0 * rb**2 &
                              !        - 60.d0 * rb**3 &
                              !        + 30.d0 * rb**4) * &
                              !        (xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_mbd/Bohr))
                              !end if
                            end if
                            a_j = a_2b_tot(k2)
                            r_vdw_j = r0_ii_SCS_2b_tot(k2)
                            o_j = o_2b_tot(k2)
                            dr_vdw_j = dr0_ii_SCS_2b(k2)
                            da_j = da_2b(k2)
                            do_j = do_2b(k2)
                            c6_2b_tot(k2) = 3.d0/2.d0 * a_i * a_j * (o_i * o_j) / &
                               (o_i + o_j)
                            f_damp_SCS_2b_tot(k2) = 1.d0/( 1.d0 + &
                                exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                            E_TS_tot = E_TS_tot - c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) & ! NOTE: f_damp is squared here only to match MBD. Delete it later.
                                       * r6_mult_2b_tot(k2) * r6_mult_0i(k2) * r6_mult_0j(k2)
                            dC6_2b = 3.d0/2.d0 * (o_i*o_j &
                                / (o_i+o_j) &
                                * (da_i*a_j + a_i*da_j) &
                                + a_i * a_j / (o_i+o_j)**2 &
                                * (do_i * o_j**2 + o_i**2 * do_j))
                            f_damp_der_2b = 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (- rjs_2b_tot(k2)/0.97d0 * 1.d0/(r_vdw_i+r_vdw_j)**2 &
                                        * (dr_vdw_i + dr_vdw_j))
                            if ( i2 == i ) then        
                              f_damp_der_2b = f_damp_der_2b + 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (-xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)))
                            else if ( j2 == i ) then
                              f_damp_der_2b = f_damp_der_2b + 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)))
                            end if
                            r6_der = 0.d0
                            if ( i2 == i ) then
                              r6_der = 6.d0/rjs_2b_tot(k2)**8 * xyz_2b_tot(c3,k2)
                            else if ( j2 == i ) then
                              r6_der = -6.d0/rjs_2b_tot(k2)**8 * xyz_2b_tot(c3,k2)
                            end if
                            forces_TS = forces_TS + dC6_2b * f_damp_SCS_2b_tot(k2) / rjs_2b_tot(k2)**6 * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_der_2b * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * dr6_mult(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + r6_der * c6_2b_tot(k2) * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * dr6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * dr6_mult_0j(k2)
                            !if ( i == 1 .and. c3 == 1 .and. abs(dr6_mult_0j(k2)) > 0.d0 ) then
                            !  write(*,*) "forces_TS", dC6_2b * f_damp_SCS_2b_tot(k2) / rjs_2b_tot(k2)**6 * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_der_2b * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * dr6_mult(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + r6_der * c6_2b_tot(k2) * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * dr6_mult_0j(k2), "upper", k2
                            !  write(*,*) "buffer term", c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * &
                            !      f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * dr6_mult_0j(k2)
                            !  write(*,*) "buffer", dr6_mult_0j(k2)
                            !end if
                          end if
                        end if
                      !end if
                    !else if ( rjs(n_tot+k_i) .le. rcut_2b .and. rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd .and. &
                    !          rjs(n_tot+k_j) .le. rcut_2b ) then
                    else if ( rjs(n_tot+k_i) .le. rcut_2b ) then
                        if (i3 .ne. j3) then
                          xyz_j = xyz(:,n_tot+k_j)/Bohr
                          if ( sqrt(sum((xyz_j-xyz_i)**2)) .le. rcut_2b/Bohr ) then
                            k2 = k2+1
                            s = neighbor_species(n_tot+k_j)
                            sub_2b_tot_list(k2) = neighbors_list(n_tot+k_j)
                            hirshfeld_v_2b_der(1:3,k2) = hirshfeld_v_cart_der_H(1:3,n_tot+k_j)
                            r0_ii_2b_tot(k2) = r0_ref(s) / Bohr !* hirshfeld_v_neigh(n_tot+k_i)**(1.d0/3.d0)
                            xyz_2b_tot(:,k2) = xyz_j-xyz_i
                            rjs_2b_tot(k2) = sqrt(sum((xyz_j-xyz_i)**2))
                            neighbor_alpha0_2b_tot(k2) = alpha0_ref(s) / Bohr**3 !* hirshfeld_v_neigh(n_tot+k_i)
                            hirshfeld_2b_tot_neigh(k2) = hirshfeld_v_neigh(n_tot+k_j)
                            if ( rjs(n_tot+k_j) .le. rcut_loc-r_buf_loc ) then
                              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j2,1)
                              a_2b_tot(k2) = a_iso(r,2)
                              o_2b_tot(k2) = o_p(r)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0)
                              da_2b(k2) = da_iso(r,c3,2)
                             if ( a_iso(r,2) > a_iso(r,1) ) then
                               do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                  - a_iso(r,2) * da_iso(r,c3,1)) / &
                                  ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                             else
                               write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                               do_pref = 0.d0
                             end if
                             do_2b(k2) = do_pref
                             dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2)
                            else if ( rjs(n_tot+k_j) .le. rcut_loc .and. rjs(n_tot+k_j) > rcut_loc-r_buf_loc ) then
                              r = findloc(sub_neighbors_list(1:n_sub_neigh(1)),j2,1)
                              rb = (rjs(n_tot+k_j)-rcut_loc+r_buf_loc)/r_buf_loc
                              a_2b_tot(k2) = a_iso(r,2) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                  + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5 ) + &
                                  central_pol(j1) * &
                                   !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                                     ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              o_2b_tot(k2) = o_p(r) * &
                                  (1.d0 - 10.d0 * rb**3 &
                                    + 15.d0 * rb**4 &
                                    - 6.d0 * rb**5) + &
                                  central_omega(j1) * &
                                     ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                  + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5) + &
                                   r0_ii_2b_tot(k2) * (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                    !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                                   ( + 10.d0 *rb**3 &
                                  - 15.d0 * rb**4 &
                                  + 6.d0 * rb**5)
                              if ( a_iso(r,2) > a_iso(r,1) ) then
                                do_pref = -0.5d0*vdw_omega_ref*omega_ref * (a_iso(r,1) * da_iso(r,c3,2) &
                                     - a_iso(r,2) * da_iso(r,c3,1)) / &
                                     ( a_iso(r,1)**2 * (a_iso(r,2)/a_iso(r,1) - 1.d0)**(3.d0/2.d0) )
                              else
                                write(*,*) "WARNING: frequency dependency failure. Use larger vdw_omega_ref."
                                do_pref = 0.d0
                              end if
                              da_2b(k2) = da_iso(r,c3,2) * &
                                  (1.d0 - 10.d0 * rb**3 &
                                   + 15.d0 * rb**4 &
                                   - 6.d0 * rb**5) + &
                                    a_iso(r,2) * &
                                 (-30.d0 * rb**2 &
                                 + 60.d0 * rb**3 &
                                 - 30.d0 * rb**4) &
                                  * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                 central_pol(j1) * & 
                                 !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2) * &
                                 ( +30.d0 * rb**2 &
                                - 60.d0 * rb**3 &
                                + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                              dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) / &
                                   (3.d0 * a_iso(r,2)) * da_iso(r,c3,2) * &
                                   (1.d0 - 10.d0 * rb**3 &
                                    + 15.d0 * rb**4 &
                                  - 6.d0 * rb**5) + &
                                   r0_ii_2b_tot(k2) * (a_iso(r,2)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * &
                                   (-30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                   r0_ii_2b_tot(k2) * &
                                   (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0) * & 
                                   !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0) * &
                                  (+30.d0 * rb**2 &
                                  - 60.d0 * rb**3 &
                                  + 30.d0 * rb**4) &
                                  * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                              if ( do_hirshfeld_gradients ) then
                                da_2b(k2) = da_2b(k2) + neighbor_alpha0_2b_tot(k2) * &
                                    hirshfeld_v_2b_der(c3,k2)*Bohr * &
                                    (+ 10.d0 *rb**3 &
                                    - 15.d0 * rb**4 &
                                    + 6.d0 * rb**5)
                                dr0_ii_SCS_2b(k2) = dr0_ii_SCS_2b(k2) + &
                                   r0_ii_2b_tot(k2) * (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * &
                                   hirshfeld_v_2b_der(c3,k2)*Bohr * &
                                   (+ 10.d0 *rb**3 &
                                   - 15.d0 * rb**4 &
                                   + 6.d0 * rb**5)       
                              end if
                              do_2b(k2) = do_pref * & 
                                  (1.d0 - 10.d0 * rb**3 &
                                   + 15.d0 * rb**4 &
                                   - 6.d0 * rb**5) + &
                                    o_p(r) * &
                                  (-30.d0 * rb**2 &
                                  + 60.d0 * rb**3 &
                                  - 30.d0 * rb**4) &
                                   * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr)) + &
                                  central_omega(j1) * &
                                  ( +30.d0 * rb**2 &
                                - 60.d0 * rb**3 &
                                 + 30.d0 * rb**4) &
                                 * ( -xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_loc/Bohr))
                            else if ( rjs(n_tot+k_j) > rcut_loc ) then
                              a_2b_tot(k2) = central_pol(j1)
                                    !neighbor_alpha0_2b_tot(k2) * hirshfeld_2b_tot_neigh(k2)
                              o_2b_tot(k2) = central_omega(j1)
                              r0_ii_SCS_2b_tot(k2) = r0_ii_2b_tot(k2) * &
                                   (central_pol(j1)/neighbor_alpha0_2b_tot(k2))**(1.d0/3.d0)
                                   !(hirshfeld_2b_tot_neigh(k2))**(1.d0/3.d0)
                              if ( do_hirshfeld_gradients ) then
                                da_2b(k2) = neighbor_alpha0_2b_tot(k2) * hirshfeld_v_2b_der(c3,k2)*Bohr
                                dr0_ii_SCS_2b(k2) = r0_ii_2b_tot(k2) * &
                                     (hirshfeld_2b_tot_neigh(k2))**(-2.d0/3.d0)/3.d0 * hirshfeld_v_2b_der(c3,k2)*Bohr
                              else
                                da_2b(k2) = 0.d0
                                dr0_ii_SCS_2b(k2) = 0.d0
                              end if
                              do_2b(k2) = 0.d0
                            end if
                            if ( rjs_2b_tot(k2)*Bohr > rcut_2b-r_buf_2b .and. rjs_2b_tot(k2)*Bohr .le. rcut_2b) then
                              rb = (rjs_2b_tot(k2)*Bohr-rcut_2b+r_buf_2b)/(r_buf_2b)
                              r6_mult_2b_tot(k2) = (1.d0 - 10.d0 * rb**3 &
                                      + 15.d0 * rb**4 &
                                      - 6.d0 * rb**5)
                              if ( i2 == i ) then 
                                dr6_mult(k2) = (- 30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * ( -xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)/(r_buf_2b/Bohr))
                              else if ( j2 == i ) then
                                dr6_mult(k2) = (- 30.d0 * rb**2 &
                                   + 60.d0 * rb**3 &
                                   - 30.d0 * rb**4) &
                                   * (xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)/(r_buf_2b/Bohr))
                              end if
                            end if
                            !if ( rjs(n_tot+k_i) .le. rcut_mbd .and. rjs(n_tot+k_i) > rcut_mbd-r_buf_mbd ) then
                            !  rb = (rjs(n_tot+k_i)-rcut_mbd+r_buf_mbd)/(r_buf_mbd)
                            !  r6_mult_0i(k2) = (10.d0 * rb**3 &
                            !          - 15.d0 * rb**4 &
                            !          + 6.d0 * rb**5)
                            !  dr6_mult_0i(k2) = (30.d0 * rb**2 &
                            !          - 60.d0 * rb**3 &
                            !          + 30.d0 * rb**4) * &
                            !          (-xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_mbd/Bohr))
                            !else if ( rjs(n_tot+k_i) .le. rcut_2b .and. rjs(n_tot+k_i) > rcut_2b-r_buf_2b ) then
                            !  rb = (rjs(n_tot+k_i)-rcut_2b+r_buf_2b)/(r_buf_2b)
                            !  r6_mult_0i(k2) = (1.d0 - 10.d0 * rb**3 &
                            !          + 15.d0 * rb**4 &
                            !          - 6.d0 * rb**5)
                              !dr6_mult_0i(k2) = (-30.d0 * rb**2 &
                              !        + 60.d0 * rb**3 &
                              !        - 30.d0 * rb**4) * &
                              !        (-xyz(c3,n_tot+k_i)/rjs(n_tot+k_i)/(r_buf_2b/Bohr))
                            !end if
                            !if ( rjs(n_tot+k_j) .le. rcut_2b .and. rjs(n_tot+k_j) > rcut_2b-r_buf_2b ) then
                            !  rb = (rjs(n_tot+k_j)-rcut_2b+r_buf_2b)/(r_buf_2b)
                            !  r6_mult_0j(k2) = (1.d0 - 10.d0 * rb**3 &
                            !          + 15.d0 * rb**4 &
                            !          - 6.d0 * rb**5)
                              !dr6_mult_0j(k2) = (-30.d0 * rb**2 &
                              !        + 60.d0 * rb**3 &
                              !        - 30.d0 * rb**4) * &
                              !        (-xyz(c3,n_tot+k_j)/rjs(n_tot+k_j)/(r_buf_2b/Bohr))
                            !end if
                            a_j = a_2b_tot(k2)
                            r_vdw_j = r0_ii_SCS_2b_tot(k2)
                            o_j = o_2b_tot(k2)
                            dr_vdw_j = dr0_ii_SCS_2b(k2)
                            da_j = da_2b(k2)
                            do_j = do_2b(k2)
                            c6_2b_tot(k2) = 3.d0/2.d0 * a_i * a_j * (o_i * o_j) / &
                               (o_i + o_j)
                            f_damp_SCS_2b_tot(k2) = 1.d0/( 1.d0 + &
                                exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) )
                            E_TS_tot = E_TS_tot - c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) & ! NOTE: f_damp is squared here only to match MBD. Delete it later.
                                       * r6_mult_2b_tot(k2) * r6_mult_0i(k2) * r6_mult_0j(k2)
                            dC6_2b = 3.d0/2.d0 * (o_i*o_j &
                                / (o_i+o_j) &
                                * (da_i*a_j + a_i*da_j) &
                                + a_i * a_j / (o_i+o_j)**2 &
                                * (do_i * o_j**2 + o_i**2 * do_j))
                            f_damp_der_2b = 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (- rjs_2b_tot(k2)/0.97d0 * 1.d0/(r_vdw_i+r_vdw_j)**2 &
                                        * (dr_vdw_i + dr_vdw_j))
                            if ( i2 == i ) then        
                              f_damp_der_2b = f_damp_der_2b + 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (-xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)))
                            else if ( j2 == i ) then
                              f_damp_der_2b = f_damp_der_2b + 20.d0 * f_damp_SCS_2b_tot(k2)**2 * &
                                        exp( -20.d0*( rjs_2b_tot(k2)/(0.97d0*(r_vdw_i + r_vdw_j)) - 1.d0 ) ) &
                                        * (1.d0/(0.97d0 * (r_vdw_i+r_vdw_j)) * (xyz_2b_tot(c3,k2)/rjs_2b_tot(k2)))
                            end if
                            r6_der = 0.d0
                            if ( i2 == i ) then
                              r6_der = 6.d0/rjs_2b_tot(k2)**8 * xyz_2b_tot(c3,k2)                    
                            else if ( j2 == i ) then
                              r6_der = -6.d0/rjs_2b_tot(k2)**8 * xyz_2b_tot(c3,k2)
                            end if
                            forces_TS = forces_TS + dC6_2b * f_damp_SCS_2b_tot(k2) / rjs_2b_tot(k2)**6 * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_der_2b * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * dr6_mult(k2) & 
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + r6_der * c6_2b_tot(k2) * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * dr6_mult_0i(k2) * r6_mult_0j(k2) &
                                  + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                                  * r6_mult_0i(k2) * dr6_mult_0j(k2)
                            !if ( i == 1 .and. c3 == 1 .and. abs(dr6_mult_0j(k2)) > 0.d0 ) then
                            !  write(*,*) "forces_TS", dC6_2b * f_damp_SCS_2b_tot(k2) / rjs_2b_tot(k2)**6 * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_der_2b * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * dr6_mult(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + r6_der * c6_2b_tot(k2) * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * r6_mult_0j(k2) &
                            !      + c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * dr6_mult_0j(k2), "upper", k2
                            !  write(*,*) "buffer term", c6_2b_tot(k2)/rjs_2b_tot(k2)**6 * &
                            !      f_damp_SCS_2b_tot(k2) * r6_mult_2b_tot(k2) &
                            !      * dr6_mult_0j(k2)
                            !  write(*,*) "buffer", dr6_mult_0j(k2)
                            !end if
                          end if
                        end if
                      !end if
                    end if
                  end do
                end if
              end do
              E_TS_tot = 1.d0/2.d0 * E_TS_tot
              forces_TS = 1.d0/2.d0 * forces_TS
              forces0(c3,i) = forces0(c3,i) + forces_TS * Hartree/Bohr
              
              !write(*,*) "E_TS_tot", E_TS_tot*Hartree
              !write(*,*) "forces_TS", forces_TS * Hartree/Bohr

              deallocate( sub_2b_tot_list, n_2b_tot_neigh, p_2b_tot, xyz_2b_tot, rjs_2b_tot, r0_ii_2b_tot, &
                          neighbor_alpha0_2b_tot, f_damp_SCS_2b_tot, hirshfeld_2b_tot_neigh, hirshfeld_v_2b_der, &
                          a_2b_tot, o_2b_tot, r0_ii_SCS_2b_tot, c6_2b_tot, r6_mult_2b_tot, dr6_mult, da_2b, &
                          do_2b, dr0_ii_SCS_2b, r6_mult_0i, r6_mult_0j, dr6_mult_0i, dr6_mult_0j )
              !if ( i == 1 .and. c3 == 1 ) then
              !write(*,*) "r6_mult, dr6_mult"
              !do k3 = 1, n_2b_sites
              !  write(*,*) k3, r6_mult(k3), dr6_mult(k3)
              !end do
              !end if


              !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
              !  write(*,*) "forces_TS", forces_TS
              !end if

              !if ( i == 1 .and. c3 == 1 .and. om == 2 ) then
              !  write(*,*) "o_mbd"
              !  k3 = 0
              !  do p = 1, n_mbd_sites
              !    write(*,*) o_mbd(k3+1)
              !    k3 = k3+n_mbd_neigh(p)
              !  end do
              !end if
      
if ( .false. ) then
        
              allocate( G_mat(1:3*n_mbd_sites,1:3*n_mbd_sites,1:n_freq) )

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
                  k3 = k3+n_mbd_neigh(p)
                end do
              end do

              !if ( i == 1 .and. c3 == 1 ) then
              !  write(*,*) "G_mat"
              !  open(unit=89, file="G_mat.dat", status="new")
              !  do p = 1, 3*n_mbd_sites
              !    write(89,*) G_mat(p,:,1)
              !  end do
              !  close(89)
              !end if

              integrand = 0.d0
              if ( c3 == 1 .and. do_total_energy ) then
                total_integrand = 0.d0
              end if
              do j = 1, n_freq
              
                force_series = 0.d0
                if ( c3 == 1 .and. do_total_energy ) then
                  total_energy_series = 0.d0
                end if
                do k2 = 1, n_order-1
                  force_series = force_series + AT_n_f(:,:,k2,j) 
                  if ( c3 == 1 .and. do_total_energy ) then
                    total_energy_series = total_energy_series - 1.d0/(k2+1)*AT_n_f(:,:,k2,j)
                  end if 
                end do
                k3 = 0
                do p = 1, n_mbd_sites
                  i2 = mbd_neighbors_list(k3+1)
                  do c1 = 1, 3
                    integrand(j) = integrand(j) + & !1.d0/(1.d0 + (omegas_mbd(j)/0.5d0)**2) * &
                    dot_product(G_mat(3*(p-1)+c1,:,j),force_series(:,3*(p-1)+c1))
                    if ( c3 == 1 .and. do_total_energy ) then
                      total_integrand(j) = total_integrand(j) + a_mbd(k3+1) / &
                            (1.d0 + (omegas_mbd(j)/o_mbd(k3+1))**2) &
                            * dot_product(T_LR(3*(p-1)+c1,:), &
                            total_energy_series(:,3*(p-1)+c1))
                    end if
                  end do
                  k3 = k3 + n_mbd_neigh(p)
                end do
              end do
              
              deallocate ( G_mat )

              if ( n_order > 1 ) then

                if ( do_nnls ) then
               
                allocate( A_nnls(1:n_freq,1:2*n_order+1) )
                allocate( b_nnls(1:n_freq) )
                allocate( coeff_nnls(1:2*n_order+1) )
                allocate( work_nnls(1:2*n_order+1) )
                allocate( ind_nnls(1:2*n_order+1) )
                allocate( omegas_nnls(1:21) )
                allocate( integrand_nnls(1:21) )
                allocate( work_integrand(1:size(integrand)) )

                A_nnls = 0.d0
                b_nnls = 0.d0
                coeff_nnls = 0.d0
                res_nnls = 0.d0
                work_nnls = 0.d0
                ind_nnls = 0.d0
                
                if (integrand(1) < 0.d0 ) then
                  work_integrand = -integrand
                else
                  work_integrand = integrand
                end if

                do i2 = 1, n_freq
                  b_nnls(i2) = work_integrand(i2)
                  A_nnls(i2,1) = 1.d0
                  do j2 = 1, n_order
                    A_nnls(i2,j2+1) = (-1.d0 * omegas_mbd(i2))**j2
                    A_nnls(i2,n_order+j2+1) = -work_integrand(i2) * omegas_mbd(i2)**(2*j2)
                  end do
                end do
                
                call nnls(A_nnls, n_freq, n_order+3, b_nnls, coeff_nnls, res_nnls, work_nnls, ind_nnls, mode_nnls)

                if ( integrand(1) < 0.d0 ) then
                  integrand_nnls = -coeff_nnls(1)
                else
                  integrand_nnls = coeff_nnls(1)
                end if
                omegas_nnls = 0.d0
                do j2 = 1, n_order
                  if ( integrand(1) < 0.d0 ) then
                    integrand_nnls(1) = integrand_nnls(1) - coeff_nnls(j2+1)*(-1.d0*omegas_nnls(1))**j2
                  else
                    integrand_nnls(1) = integrand_nnls(1) + coeff_nnls(j2+1)*(-1.d0*omegas_nnls(1))**j2
                  end if
                end do
                denom = 1.d0
                do j2 = 1, n_order
                  denom = denom + coeff_nnls(n_order+j2+1) * omegas_nnls(1)**(2**j2)
                end do
                integrand_nnls(1) = integrand_nnls(1)/denom
                do i2 = 2, 21
                  omegas_nnls(i2) = omegas_nnls(i2-1)+0.2d0
                  do j2 = 1, n_order
                    if ( integrand(1) < 0.d0 ) then
                      integrand_nnls(i2) = integrand_nnls(i2) - coeff_nnls(j2+1)*(-1.d0*omegas_nnls(i2))**j2
                    else
                      integrand_nnls(i2) = integrand_nnls(i2) + coeff_nnls(j2+1)*(-1.d0*omegas_nnls(i2))**j2
                    end if
                  end do
                  denom = 1.d0
                  do j2 = 1, n_order
                    denom = denom + coeff_nnls(n_order+j2+1) * omegas_nnls(i2)**(2**j2)
                  end do
                  integrand_nnls(i2) = integrand_nnls(i2)/denom
                end do

                integral = 0.d0
                call integrate("trapezoidal", omegas_nnls, integrand_nnls, omegas_nnls(1), omegas_nnls(size(omegas_nnls)), integral)

                deallocate( A_nnls, b_nnls, coeff_nnls, work_nnls, ind_nnls, omegas_nnls, integrand_nnls, work_integrand )
          
                else

                integral = 0.d0
                call integrate("trapezoidal", omegas_mbd, integrand, omegas_mbd(1), omegas_mbd(n_freq), integral)

                if ( c3 == 1 .and. do_total_energy ) then
                  E_tot = 0.d0
                  call integrate("trapezoidal", omegas_mbd, total_integrand, omegas_mbd(1), omegas_mbd(n_freq), E_tot)
                  E_tot = (E_tot/(2.d0*pi)+ E_TS_tot) * Hartree
                end if

                end if

                forces0(c3,i) = forces0(c3,i) + (1.d0/(2.d0*pi) * integral + forces_TS) * Hartree/Bohr
                
                if ( c3 == 1 .and. do_total_energy ) then
                  write(*,*) "Total energy of sphere", i, E_tot
                end if

                !write(*,*) "MBD force", i, c3, 1.d0/(2.d0*pi) * integral * Hartree/Bohr
                write(*,*) "Total force", i, c3, forces0(c3,i)

              end if
            
            end if ! om loop(?)
            
end if

          end do ! c3 loop
        
          deallocate( da_SCS, dT, b_der, g_func_der, h_func_der, d_der, f_damp_der )
          deallocate( dB_mat )

        end if ! do_derivatives
      
      
        if ( om == 2 ) then
          deallocate( sub_2b_list, xyz_2b, rjs_2b, r0_ii_2b, neighbor_alpha0_2b, &
                      hirshfeld_2b_neigh, a_2b, o_2b, r0_ii_SCS_2b, f_damp_SCS_2b, C6_2b, r6_mult )
        end if

      
if ( .false. ) then  
        if ( om == 2 ) then
          deallocate( T_LR, r0_ii_SCS, f_damp_SCS, AT, AT_n, energy_series, omegas_mbd, integrand, n_mbd_neigh, &
                      mbd_neighbors_list, p_mbd, r0_ii_mbd, neighbor_alpha0_mbd, xyz_mbd, rjs_mbd, T_mbd, a_mbd, &
                      rjs_0_mbd, xyz_0_mbd, o_mbd, hirshfeld_mbd_neigh, &
                      T_LR_mult_0i, T_LR_mult_0j, T_LR_mult_ij, T_LR_mult_0ij, T_LR_mult_0ji )
        end if
                    
        if ( do_derivatives .and. om == 2 ) then
          deallocate( da_mbd, AT_n_f, dT_mbd, f_damp_der_mbd, f_damp_der_SCS, dT_LR, force_series, &
                      do_mbd, dr0_ii_SCS, dT_LR_mult_0i, dT_LR_mult_0j,&
                      dT_LR_mult_ij, dT_LR_mult_0ij, dT_LR_mult_0ji )
          
          if ( do_total_energy ) then
            deallocate( total_energy_series, total_integrand )
          end if
          if ( do_hirshfeld_gradients ) then
            deallocate( hirshfeld_v_mbd_der )
          end if
        end if
end if

      end do ! om

      deallocate( n_sub_neigh, sub_neighbors_list, xyz_H, rjs_H, r0_ii, neighbor_alpha0, neighbor_sigma, omegas, &
                  T_func, b_i, d_vec, g_func, h_func, f_damp, a_SCS, ipiv, p_list, work_arr, rjs_0, &
                  a_iso, o_p, T_SR, T_SR_mult, d_arr_i, d_arr_o, d_mult_i, d_mult_o, dT_SR_mult, d_dmult_i, &
                  d_dmult_o, hirshfeld_sub_neigh, inner_damp )
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

    !write(*,*) "Local energies and forces timing", time2-time1

    !write(*,*) "E_MBD", E_MBD
      
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
