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
  use psb_d_base_prec_mod
  use psb_d_ilu_fact_mod
  use psb_d_ainv_fact_mod
  use psb_d_invk_fact_mod
  use psb_d_invt_fact_mod

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
  subroutine get_scs_polarizabilities( hirshfeld_v, hirshfeld_v_cart_der, &
                                       n_neigh, neighbors_list, neighbor_species, &
                                       rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, hirshfeld_v_neigh, &
                                       sR, d, c6_ref, r0_ref, alpha0_ref, do_forces, alpha_SCS0, &
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
    real*8, intent(inout) :: energies(:), forces0(:,:), alpha_SCS0(:,:)
!   Internal variables
    real*8, allocatable :: neighbor_c6_ii(:), r0_ii(:), &
                           f_damp(:), neighbor_alpha0(:), &
                           T_func(:), h_func(:,:), g_func(:,:), &
                           omegas(:), omega_i(:), &
                           alpha_i(:,:), sigma_i(:,:), sigma_ij(:), T_SR(:,:,:), B_mat(:,:,:), &
                           rjs_H(:), xyz_H(:,:), A_mat(:,:,:), work_arr(:), A_i(:,:), &
                           dT(:,:), dT_SR(:,:,:,:,:), f_damp_der(:,:), &
                           g_func_der(:,:,:), h_func_der(:,:,:), dA_mat(:,:,:,:,:), &
                           dalpha(:,:,:), dalpha_full(:,:,:,:), &
                           coeff_h_der(:), terms(:), dT_SR_A_mat(:,:), dT_SR_v(:,:,:,:,:), &
                           dB_mat(:,:,:,:,:), dB_mat_v(:,:,:,:,:), &
                           coeff_der(:,:,:,:,:), coeff_fdamp(:,:,:,:,:), dg(:), &
                           dh(:), hirshfeld_v_cart_der_H(:,:), &
                           alpha_k(:,:), sigma_k(:,:), s_i(:), s_j(:), &
                           T_SR_i(:,:,:), B_mat_i(:,:,:), alpha_SCS_i(:,:), A_mat_i(:,:,:), alpha_SCS(:,:), &
                           dB_mat_n(:,:,:,:,:), dA_mat_n(:,:,:,:,:), dBA_n(:,:), I_mat(:,:), &
                           dalpha_n(:,:,:), AdB_n(:,:), hirshfeld_v_neigh_H(:), hirshfeld_v_cart_der2(:,:), &
                           rjs2(:), xyz2(:,:), hirshfeld_v_neigh2(:), rjs_central2(:), rjs_central(:), &
                           xyz_central2(:,:), xyz_central(:,:), dalpha2(:,:,:), I_aT(:,:,:), a_vec(:,:,:), &
                           regularization(:,:), I_aT_a_vec(:,:), I_aT2(:,:), BTB_reg(:,:,:), B_reg(:,:,:), &
                           BTB_reg_copy(:,:), a_SCS(:,:,:), da_vec(:,:,:,:,:), dBTB(:,:), vect1(:,:), &
                           vect2(:,:), vect3(:,:), da_SCS(:,:,:,:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, dist, f_damp_SCS_ij, t1, t2, rcut_H, buffer_H, rbuf, fcut, dfcut, &
              reg_param
    integer, allocatable :: ipiv(:), n_sneigh(:), sneighbors_list(:), p_to_i(:), i_to_p(:), &
                            neighbors_list2(:), local_neighbors(:), neighbor_species2(:), &
                            neighbor_species_H(:), neighbors_list3(:), n_sneigh_vder(:), &
                            sneighbors_list_vder(:), n_sneigh_list2(:), &
                            n_sneigh_list(:), n_ssites_list(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, n_order, n_tot, n_spairs, n_beg, n_end, &
               n_ssites, n_spairs_vder, n_freq, n_iter
    integer :: i, i2, j, j2, k, k2, k3, a, a2, c1, c2, c3, lwork, b, p, q, n_count
    logical :: do_timing = .false., do_hirshfeld_gradients = .true., nonlocal = .false., &
               series_expansion = .true., total_energy = .true., series_average = .true., &
               new_implementation = .false., do_derivatives = .true., iterative = .true.
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

!   Change these to be input variables (NOTE THAT THEY ARE IN ANGSTROMS!):
    write(*,*) "rcut", rcut
!TEST2 : CHANGE EACH OCCURENCE OF rcut_vdw BACK TO rcut
!    rcut_vdw = rcut-7.d0
!    write(*,*) "rcut_vdw", rcut_vdw
    ! n_order has to be at least 2
!    n_order = 6
!    n_order = 3

!    write(*,*) "neighbors_list", neighbors_list

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

!   This should allow to only take a subset of atoms for parallelization:
    n_beg = 1
    n_end = 1

!   Number of frequencies
    n_freq = size(alpha_SCS0, 2)

    if ( new_implementation ) then
    allocate( alpha_SCS(n_beg:n_end,1:n_freq) )
    if ( do_derivatives ) then
      allocate( dalpha2(1:n_pairs,1:3,1:n_freq) )
    end if
    allocate( n_ssites_list(n_beg:n_end) )
    allocate( n_sneigh_list2(1:n_pairs) )
    alpha_SCS = 0
    dalpha2 = 0
    n_ssites_list = 0
    n_sneigh_list2 = 0
    n_count = 0

!   Frequencies used for integration:
    allocate( omegas(1:n_freq) )
    omega = 0.d0
    do i = 1, n_freq
      omegas(i) = omega
      omega = omega + 0.4d0
    end do
!   Some related quantities:
    allocate( s_i(1:n_freq) )
    allocate( s_j(1:n_freq) )
    allocate( sigma_ij(1:n_freq) )
!    allocate( alpha_i(n_beg:n_end,1:11) )
!    allocate( sigma_i(n_beg:n_end,1:11) )
    allocate( coeff_h_der(1:n_freq) )
    allocate( terms(1:n_freq) )
    allocate( dg(1:n_freq) )
    allocate( dh(1:n_freq) )
    allocate( A_i(1:3,1:3) )


    write(*,*) "max of neighbors list:", maxval(neighbors_list)
!   Let's try to build a sub neighbor list for one atomic environment:
    call cpu_time(time1)
    allocate( n_sneigh(1:maxval(neighbors_list)) )
    allocate( neighbors_list2(1:n_pairs) )
    allocate( hirshfeld_v_cart_der2(1:3,1:n_pairs) )
    allocate( rjs2(1:n_pairs) )
    allocate( xyz2(1:3,1:n_pairs) )
    allocate( hirshfeld_v_neigh2(1:n_pairs) )
    allocate( neighbor_species2(1:n_pairs) )
    allocate( p_to_i(1:maxval(neighbors_list)) )
    allocate( i_to_p(1:maxval(neighbors_list)) )
    allocate( neighbors_list3(1:n_pairs) )
    allocate( n_sneigh_vder(1:maxval(neighbors_list)) )
    allocate( ij_buffer2(1:n_pairs) )
    n_tot = sum(n_neigh(1:n_beg))-n_neigh(n_beg)  ! Can you just do sum(n_neigh(1:n_beg-1)) or what happens if n_beg = 1?
!    write(*,*) "n_tot init:", n_tot

    do i = n_beg, n_end
      allocate( local_neighbors(1:n_neigh(i)) )
      allocate( i0_buffer2(1:n_neigh(i)) )
      allocate( rjs_central2(1:n_neigh(i)) )
      allocate( xyz_central2(1:3,1:n_neigh(i)) )
      local_neighbors = neighbors_list(n_tot+1:n_tot+n_neigh(i))
      n_ssites = 0
      n_sneigh = 0
      n_sneigh_vder = 0
      neighbors_list2 = 0
      neighbors_list3 = 0
      hirshfeld_v_cart_der2 = 0
      rjs2 = 0
      xyz2 = 0
      hirshfeld_v_neigh2 = 0
      neighbor_species2 = 0
      i0_buffer2 = .false.
      ij_buffer2 = .false.
      rjs_central2 = 0.d0
      xyz_central2 = 0.d0
      p_to_i = 0
      i_to_p = 0
      k = 0
      k2 = 0
      k3 = 0
      n_spairs_vder = 0
      write(*,*) "initialization ok"
      do p = 1, n_neigh(i)
!TEST2
        write(*,*) "new p", p
!        if ( rjs(n_tot+p) < rcut ) then
        if ( rjs(n_tot+p) < rcut) then
          n_ssites = n_ssites + 1
          i2 = local_neighbors(p)
          write(*,*) "i2", i2
          k2 = k2+1
          k3 = k3+1
          p_to_i(n_ssites) = i2
          i_to_p(i2) = n_ssites
          k = sum(n_neigh(1:i2)) - n_neigh(i2) + 1
          neighbors_list2(k2) = neighbors_list(k)
!          hirshfeld_v_cart_der2(1:3,k2) = hirshfeld_v_cart_der(1:3,k)
          rjs2(k2) = rjs(k)
          xyz2(1:3,k2) = xyz(1:3,k)
          hirshfeld_v_neigh2(k2) = hirshfeld_v_neigh(k)
          neighbor_species2(k2) = neighbor_species(k)
          n_sneigh(n_ssites) = n_sneigh(n_ssites) + 1
          n_sneigh_vder(n_ssites) = n_sneigh_vder(n_ssites) + 1
          neighbors_list3(k3) = neighbors_list(k)
          hirshfeld_v_cart_der2(1:3,k3) = hirshfeld_v_cart_der(1:3,k)
          n_spairs_vder = n_spairs_vder + 1
          ! Go through the neighbors of i2
!buffer test
          if ( rjs(n_tot+p) > rcut - buffer ) then
!          if ( rjs(n_tot+p) > rcut/2.d0 ) then
!          if ( rjs(n_tot+p) > 0.d0 ) then
            i0_buffer2(n_ssites) = .true.
            rjs_central2(n_ssites) = rjs(n_tot+p)
            xyz_central2(1:3,n_ssites) = xyz(1:3,n_tot+p)
!            write(*,*) "rjs in buffer:", rjs(n_tot+p)
          end if
!TEST2
          do j2 = 2, n_neigh(i2)
            k = k+1
            j = neighbors_list(k)
!            write(*,*) "j", j
!            write(*,*) "local neighbors:", local_neighbors
            ! Is j a neighbor of i as well?
            if ( any(local_neighbors == j) ) then
!              write(*,*) "hello?"
              q = findloc(local_neighbors,j,1)
              ! Is it within cutoff of central atom?
!              write(*,*) "check if within central cutoff"
!              if ( rjs(n_tot+q) < rcut ) then
              if ( rjs(n_tot+q) < rcut) then
                k3 = k3+1
                n_spairs_vder = n_spairs_vder + 1
                n_sneigh_vder(n_ssites) = n_sneigh_vder(n_ssites) + 1
                neighbors_list3(k3) = neighbors_list(k)
                hirshfeld_v_cart_der2(1:3,k3) = hirshfeld_v_cart_der(1:3,k)
                ! Is it also within cutoff of i2?
!TEST2
!                write(*,*) "check for neighbor-neighbor"
!                if ( rjs(k) < rcut ) then
                  n_sneigh(n_ssites) = n_sneigh(n_ssites) + 1
                  k2 = k2+1
                  neighbors_list2(k2) = j
                  rjs2(k2) = rjs(k)
                  xyz2(1:3,k2) = xyz(1:3,k)
                  hirshfeld_v_neigh2(k2) = hirshfeld_v_neigh(k)
                  neighbor_species2(k2) = neighbor_species(k)
!TEST2
!                  write(*,*) "check if buffer"
!                  if ( rjs(k) > rcut - buffer ) then
!                  if ( rjs(k) > rcut/2.d0 ) then
!                  if ( rjs(k) > 0.d0 ) then
!                    ij_buffer2(k2) = .true.
!                  end if
!                end if
              end if
            end if
          end do
        end if
      end do

      write(*,*) "n_sneigh", n_sneigh

      call cpu_time(time2)

      write(*,*) "neighbor list stuff 1:", time2-time1
      call cpu_time(time1)

      n_spairs = sum(n_sneigh)
      allocate( sneighbors_list(1:n_spairs) )
      allocate( hirshfeld_v_cart_der_H(1:3,1:n_spairs_vder) )
      allocate( rjs_H(1:n_spairs) )
      allocate( xyz_H(1:3,1:n_spairs) )
      allocate( hirshfeld_v_neigh_H(1:n_spairs) )
      allocate( neighbor_species_H(1:n_spairs) )
      allocate( sneighbors_list_vder(1:n_spairs_vder) )
      allocate( i0_buffer(1:n_ssites) )
      allocate( ij_buffer(1:n_spairs) )
      allocate( rjs_central(1:n_ssites) )
      allocate( xyz_central(1:3,1:n_ssites) )

      sneighbors_list = 0
      hirshfeld_v_cart_der_H = 0
      rjs_H = 0
      xyz_H = 0
      hirshfeld_v_neigh_H = 0
      neighbor_species_H = 0
      sneighbors_list_vder = 0
      i0_buffer = .false.
      ij_buffer = .false.
      rjs_central = 0.d0
      xyz_central = 0.d0

      do k = 1, n_ssites
        i0_buffer(k) = i0_buffer2(k)
        rjs_central(k) = rjs_central2(k)
        xyz_central(1:3,k) = xyz_central2(1:3,k)
      end do

!      write(*,*) "i0_buffer", i0_buffer

      do k = 1, n_spairs
        if (neighbors_list2(k) > 0) then
          sneighbors_list(k) = neighbors_list2(k)
!          hirshfeld_v_cart_der_H(1:3,k) = hirshfeld_v_cart_der2(1:3,k)
          rjs_H(k) = rjs2(k)
          xyz_H(1:3,k) = xyz2(1:3,k)
          hirshfeld_v_neigh_H(k) = hirshfeld_v_neigh2(k)
          neighbor_species_H(k) = neighbor_species2(k)
          ij_buffer(k) = ij_buffer2(k)
        end if
      end do

!      write(*,*) "sum of n_sneigh_vder", sum(n_sneigh_vder)
!      write(*,*) "n_spairs_vder", n_spairs_vder
!      write(*,*) "n_spairs", n_spairs

      do k = 1, n_spairs_vder
        if (neighbors_list3(k) > 0) then
          sneighbors_list_vder(k) = neighbors_list3(k)
          hirshfeld_v_cart_der_H(1:3,k) = hirshfeld_v_cart_der2(1:3,k)
        end if
      end do
!      k2 = 0
!      do p = 1, n_ssites
!        write(*,*) "atom", p_to_i(p)
!        do a2 = 1, n_sneigh_vder(p)
!          k2 = k2+1
!          write(*,*) sneighbors_list_vder(k2)
!        end do
!      end do
!      write(*,*) "what"
!      k2 = 0
!      do p = 1, n_ssites
!        write(*,*) "atom", p_to_i(p)
!        do j2 = 1, n_sneigh(p)
!          k2 = k2+1
!          write(*,*) sneighbors_list(k2)
!        end do
!      end do

!      write(*,*) "sneighbors_list"
!      k = 0
!      do p = 1, n_ssites
!        write(*,*) "p", p
!        write(*,*) "p to i", p_to_i(p)
!        do q = 1, n_sneigh(p)
!          k = k+1
!          write(*,*) sneighbors_list(k)
!        end do
!      end do
!      write(*,*) "k", k
!      write(*,*) size(sneighbors_list)
!      write(*,*) "n_neigh(i)", n_neigh(i)
!      write(*,*) "n_ssites", n_ssites

      call cpu_time(time2)
      write(*,*) "neighbor list stuff 2:", time2-time1
      call cpu_time(time1)

!TEST2
      ! rcut_H = rcut/Bohr
      rcut_H = rcut/Bohr
!TEST2
      buffer_H = buffer/Bohr
!      buffer_H = rcut_H
!      buffer_H = rcut_H/2.d0
      rjs_H = rjs_H/Bohr
      xyz_H = xyz_H/Bohr
      hirshfeld_v_cart_der_H = hirshfeld_v_cart_der_H*Bohr
      rjs_central = rjs_central/Bohr
      xyz_central = xyz_central/Bohr

      allocate( neighbor_c6_ii(1:n_spairs) )
      allocate( r0_ii(1:n_spairs) )
      allocate( neighbor_alpha0(1:n_spairs) )
      allocate( omega_i(1:n_spairs) )

      do k = 1, n_spairs
        j = neighbor_species_H(k)
        neighbor_c6_ii(k) = c6_ref(j) / (Hartree*Bohr**6)
        r0_ii(k) = r0_ref(j) / Bohr
        neighbor_alpha0(k) = alpha0_ref(j) / Bohr**3
      end do

!     Precompute some other pair quantities
      neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh_H**2
!     This is slow, could replace by Taylor expansion maybe
      r0_ii = r0_ii * hirshfeld_v_neigh_H**(1.d0/3.d0)
      neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh_H
      omega_i = (4.d0 * neighbor_c6_ii)/(3.d0*neighbor_alpha0**2)

!      write(*,*) "omega_i:", omega_i

      allocate( alpha_k(1:n_spairs,1:n_freq) )
      allocate( sigma_k(1:n_spairs,1:n_freq) )

!      k2 = 1
!      do i2 = n_beg, n_end
!        do k = 1, 11
!          alpha_i(i2,k) = neighbor_alpha0(k2)/(1.d0 + omegas(k)**2/omega_i(k2)**2)
!          sigma_i(i2,k) = (sqrt(2.d0/pi) * alpha_i(i2,k)/3.d0)**(1.d0/3.d0)
!        end do
!        k2 = k2+n_neigh(i)
!      end do

      do k = 1, n_freq
        alpha_k(:,k) = neighbor_alpha0/(1.d0 + omegas(k)**2/omega_i**2)
        sigma_k(:,k) = (sqrt(2.d0/pi) * alpha_k(:,k)/3.d0)**(1.d0/3.d0)
      end do

      call cpu_time(time2)
      write(*,*) "pair quantities:", time2-time1
      call cpu_time(time1)

!      write(*,*) "T_func"
      
      allocate( f_damp(1:n_spairs) )
      allocate( T_func(1:9*n_spairs) )
      allocate( h_func(1:9*n_spairs,1:n_freq) )
      allocate( g_func(1:n_spairs,1:n_freq) )
      
      T_func = 0.d0
      f_damp = 0.d0
      k = 0
      do p = 1, n_ssites
        k = k+1
        r_vdw_i = r0_ii(k)
!        if (n_sneigh(p) > 1) then
        do q = 2, n_sneigh(p)
          k = k+1
          r_vdw_j = r0_ii(k)
!TEST2
!          if( rjs_H(k) < rcut_H )then
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
!TEST2
!          end if
        end do
!        end if
      end do

      g_func = 0.d0
      h_func = 0.d0
      k = 0
      do p = 1, n_ssites
        k = k+1
        s_i = sigma_k(k,:)
!        if (n_sneigh(p) > 1) then
        do q = 2, n_sneigh(p)
          k = k+1
          s_j = sigma_k(k,:)
!          j = sneighbors_list(k)
!TEST2
!          if( rjs_H(k) < rcut_H )then
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
!TEST2
!          end if
        end do
!        end if
      end do
      
!      write(*,*) "T_SR"

      call cpu_time(time2)
      write(*,*) "more pair quantities:", time2-time1
      call cpu_time(time1)

      allocate( T_SR(1:3*n_ssites,1:3*n_ssites,1:n_freq) )
      allocate( B_mat(1:3*n_ssites,1:3*n_ssites,1:n_freq) )

      T_SR = 0.d0
      B_mat = 0.d0
      k = 0
      do p = 1, n_ssites
        k = k+1
        do c1 = 1, 3
          B_mat(3*(p-1)+c1,3*(p-1)+c1,:) = 1.d0/alpha_k(k,:)
        end do
!        if (n_sneigh(p) > 1) then
        do j2 = 2, n_sneigh(p)
          k = k+1
          j = sneighbors_list(k)
          q = i_to_p(j)
!TEST2
!          if( rjs_H(k) < rcut_H )then
            k2 = 9*(k-1)
            do c1 = 1, 3
              do c2 = 1, 3
                k2 = k2+1
                T_SR(3*(p-1)+c1,3*(q-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(k2) * &
                                                g_func(k,:) + h_func(k2,:))
              end do
            end do
!TEST2
            if ( ij_buffer(k) ) then
              rbuf = (rjs_H(k)-rcut_H+buffer_H)/buffer_H
              T_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:) = T_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:) * &
                  (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
            end if
!TEST2
!          end if
        end do
!        end if
      end do

      do p = 1, n_ssites 
        if ( i0_buffer(p) ) then
          rbuf = (rjs_central(p)-rcut_H+buffer_H)/buffer_H
          T_SR(3*(p-1)+1:3*(p-1)+3,:,:) = T_SR(3*(p-1)+1:3*(p-1)+3,:,:) * (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
          T_SR(:,3*(p-1)+1:3*(p-1)+3,:) = T_SR(:,3*(p-1)+1:3*(p-1)+3,:) * (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
        end if
      end do

!      write(*,*) "T_SR"
!      do p = 1, 3*n_ssites
!        write(*,*) T_SR(p,:,1)
!      end do

      B_mat = B_mat + T_SR
      
!      write(*,*) "dT"
      call cpu_time(time2)
      write(*,*) "B_mat:", time2-time1
      call cpu_time(time1)

      if ( do_derivatives ) then

      allocate( dT(1:9*n_spairs,1:3) )

      dT = 0.d0
      k = 0
      do p = 1, n_ssites
        k = k+1
        if (n_sneigh(p) > 1) then
        do j2 = 2, n_sneigh(p)
          k = k+1
!TEST2
!          if (rjs_H(k) < rcut_H) then
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
!TEST2
!          end if
        end do
        end if
      end do

      call cpu_time(time2)
      write(*,*) "dT:", time2-time1
      call cpu_time(time1)

!      write(*,*) "dT_SR"
      
      allocate( dT_SR(1:3*n_ssites,1:3*n_ssites,1:n_ssites,1:3,1:n_freq) )
      allocate( f_damp_der(1:n_spairs,1:3) )
      allocate( g_func_der(1:n_spairs,1:3,1:n_freq) )
      allocate( h_func_der(1:9*n_spairs,1:3,1:n_freq) )

      f_damp_der = 0.d0
      g_func_der = 0.d0
      h_func_der = 0.d0
      dT_SR = 0.d0

      do a = 1, n_ssites
        k = 0
        k3 = 0
        do p = 1, n_ssites
          k = k+1
!          if (n_sneigh(p) > 1) then
          r_vdw_i = r0_ii(k)
          do j2 = 2, n_sneigh(p)
            k = k+1
            r_vdw_j = r0_ii(k)
            j = sneighbors_list(k)
            q = i_to_p(j)
!TEST2
!            if (rjs_H(k) < rcut_H) then
              if (a == p .or. a == q) then
                sigma_ij = sqrt(sigma_k(k3+1,:)**2 + sigma_k(k,:)**2)
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
                      dT_SR(3*(p-1)+c1,3*(q-1)+c2,a,c3,:) = f_damp_der(k,c3) * T_func(k2) * &
                                                            g_func(k,:) - (1.d0 - f_damp(k)) * dT(k2,c3) * &
                                                            g_func(k,:) - g_func_der(k,c3,:) * (1.d0 - f_damp(k)) * &
                                                            T_func(k2) - f_damp_der(k,c3) * h_func(k2,:) + &
                                                            h_func_der(k2,c3,:) * (1.d0 - f_damp(k))
                    end do
                  end do
                end do
                if (a == p) then
                  dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,:,:) = -dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,:,:)
                end if
              end if
!TEST2
!            end if
          end do
!          end if
          k3 = k3+n_sneigh(p)
        end do
      end do

      call cpu_time(time2)
      write(*,*) "dT_SR:", time2-time1
      call cpu_time(time1)
      
!      write(*,*) "dB_mat"

      allocate( dB_mat(1:3*n_ssites,1:3*n_ssites,1:n_ssites,1:3,1:n_freq) )
      allocate( dT_SR_v(1:3*n_ssites,1:3*n_ssites,1:n_ssites,1:3,1:n_freq) )
      allocate( dB_mat_v(1:3*n_ssites,1:3*n_ssites,1:n_ssites,1:3,1:n_freq) )
      allocate( coeff_der(1:n_ssites,1:n_ssites,1:3,1:3,1:n_freq) )
      allocate( coeff_fdamp(1:n_ssites,1:n_ssites,1:3,1:3,1:n_freq) )

      call cpu_time(time2)
      write(*,*) "allocation:", time2-time1
      call cpu_time(time1)

      dB_mat = 0.d0
      if (do_hirshfeld_gradients) then
        dB_mat_v = 0.d0
        dT_SR_v = 0.d0
        coeff_fdamp = 0.d0
        coeff_der = 0.d0
        k = 0
        k3 = 0
        do p = 1, n_ssites
          k = k+1
!          if (n_sneigh(p) > 1) then
          r_vdw_i = r0_ii(k)
          do j2 = 2, n_sneigh(p)
            k = k+1
            r_vdw_j = r0_ii(k)
            j = sneighbors_list(k)
            q = i_to_p(j)
            !if (rjs_H(k) < rcut_H) then
              sigma_ij = sqrt(sigma_k(k3+1,:)**2 + sigma_k(k,:)**2)
              S_vdW_ij = sR*(r_vdw_i + r_vdw_j)
              exp_term = exp(-d*(rjs_H(k)/S_vdW_ij - 1.d0))
              k2 = 9*(k-1)
              do c1 = 1, 3
                do c2 = 1, 3
                  k2 = k2+1
                  coeff_fdamp(p,q,c1,c2,:) = (d*sR*rjs_H(k))/(3*S_vdW_ij**2) * exp_term/(1.d0+exp_term)**2 * &
                                             (-T_func(k2) * g_func(k,:) + h_func(k2,:))
                  dg = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * (rjs_H(k)/sigma_ij)**2 * &
                       (-rjs_H(k)/(3.d0 * sigma_ij**3))
                  dh = 4.d0/sqrt(pi) * exp(-rjs_H(k)**2/sigma_ij**2) * xyz_H(c1,k)*xyz_H(c2,k) / &
                       (sigma_ij**5 * rjs_H(k)**2) * (-1.d0 + 2.d0/3.d0 * (rjs_H(k)/sigma_ij)**2)
                  coeff_der(p,q,c1,c2,:) = (1.d0-f_damp(k)) * (-T_func(k2)*dg+dh)
                end do
              end do
            !end if
          end do
!          end if
          k3 = k3 + n_sneigh(p)
        end do

        call cpu_time(time2)
        write(*,*) "dB_mat_v 1:", time2-time1
        call cpu_time(time1)

!        write(*,*) "coeff_fdamp", coeff_fdamp(1,3,1,1,1)
!        write(*,*) "coeff_der", coeff_der(1,3,1,1,1)

        k2 = 0
        k3 = 0
!        write(*,*) "n_ssites", n_ssites
        do p = 1, n_ssites
          r_vdw_i = r0_ii(k3+1)
          i2 = p_to_i(p)
          do a2 = 1, n_sneigh_vder(p)  !NOTE: THIS NEEDS TO GO FROM 1 TO N_SSITES!!!!!!!!! hirshfeld_v_cart_der_H should be longer, "symmetry" of dT_SR_v requires it
            k2 = k2+1
            a = i_to_p(sneighbors_list_vder(k2))
!            write(*,*) "a", a            
!            write(*,*) "k2, sneighbor", k2, sneighbors_list(k2)
!            write(*,*) "p, a:", p, a
            do c3 = 1, 3
              do c1 = 1, 3
                dB_mat_v(3*(p-1)+c1,3*(p-1)+c1,a,c3,:) = -1.d0/(hirshfeld_v_neigh_H(k3+1)*alpha_k(k3+1,:)) * &
                              hirshfeld_v_cart_der_H(c3,k2)
              end do
!              if (n_sneigh(p) > 1) then
              do k = 1, n_freq
                do j2 = 2, n_sneigh(p)
                  !if (rjs_H(k3+j2) < rcut_H) then
                    j = sneighbors_list(k3+j2)
                    q = i_to_p(j)
!                    write(*,*) "q", q
                    do c1 = 1, 3
                      do c2 = 1, 3
                        dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) = dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) + &
                          coeff_der(p,q,c1,c2,k) * (sigma_k(k3+1,k)**2/hirshfeld_v_neigh_H(k3+1) * &
                          hirshfeld_v_cart_der_H(c3,k2))
                        dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) = &
                          dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) + &
                          coeff_fdamp(p,q,c1,c2,k) * r_vdw_i/hirshfeld_v_neigh_H(k3+1) * hirshfeld_v_cart_der_H(c3,k2)
!                        if ( a == 2 .and. c3 == 1 .and. k == 1 .and. p == 1 .and. q == 3 .and. c1 == 1 .and. c2 == 1) then
!                          write(*,*) "r_vdw_i", r_vdw_i
!                          write(*,*) "hirshfeld_v_neigh_H", hirshfeld_v_neigh_H(k3+1)
!                          write(*,*) "sigma_k", sigma_k(k3+1,k)
!                          write(*,*) "hirshfeld_v_cart_der_H", hirshfeld_v_cart_der_H(c3,k2)
!                          write(*,*) "dT_SR_v", dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k)
!                        end if 
                      end do
                    end do
                  !end if      
                end do
              end do
!              end if
            end do
          end do
          k3 = k3 + n_sneigh(p)
        end do
        k2 = 0
        k3 = 0
        do q = 1, n_ssites
          r_vdw_j = r0_ii(k3+1)
          do a2 = 1, n_sneigh_vder(q)
            k2 = k2+1
!            if (n_sneigh(p) > 1) then
            a = i_to_p(sneighbors_list_vder(k2))
!            write(*,*) "k2, sneighbor", k2, sneighbors_list(k2)
!            write(*,*) "q, a:", q, a
            do c3 = 1, 3
              do k = 1, n_freq
                do j2 = 2, n_sneigh(q)
                  !if (rjs_H(k3+j2) < rcut_H) then
                    i2 = sneighbors_list(k3+j2)
                    p = i_to_p(i2)
!                    write(*,*) "p", p
                    do c1 = 1, 3
                      do c2 = 1, 3
                        dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) = dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) + &
                          coeff_der(p,q,c1,c2,k) * (sigma_k(k3+1,k)**2/hirshfeld_v_neigh_H(k3+1) * &
                          hirshfeld_v_cart_der_H(c3,k2))
                        dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) = &
                          dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k) + &
                          coeff_fdamp(p,q,c1,c2,k) * r_vdw_j/hirshfeld_v_neigh_H(k3+1) * hirshfeld_v_cart_der_H(c3,k2)
!                        if ( a == 2 .and. c3 == 1 .and. k == 1 .and. p == 1 .and. q == 3 .and. c1 == 1 .and. c2 == 1) then
!                          write(*,*) "r_vdw_j", r_vdw_j
!                          write(*,*) "hirshfeld_v_neigh_H", hirshfeld_v_neigh_H(k3+1)
!                          write(*,*) "sigma_k", sigma_k(k3+1,k)
!                          write(*,*) "hirshfeld_v_cart_der_H", hirshfeld_v_cart_der_H(c3,k2)
!                          write(*,*) "dT_SR_v", dT_SR_v(3*(p-1)+c1,3*(q-1)+c2,a,c3,k)
!                        end if
                      end do
                    end do
                  !end if
                end do
              end do
            end do
!            end if
          end do
          k3 = k3 + n_sneigh(q)
        end do

!        write(*,*) "dT_SR_v:", dT_SR_v(1,3*(3-1)+1,2,1,1)

        dT_SR = dT_SR + dT_SR_v

        dB_mat = dB_mat_v
        call cpu_time(time2)
        write(*,*) "dB_mat_v 2:", time2-time1
        call cpu_time(time1)
      end if

      ! Buffer region:
!      allocate( dfB(1:3*n_ssites,1:3*n_ssites,1:n_ssites,1:3,1:11) )
      k = 0
!      k2 = 0
!      k3 = 0
      do p = 1, n_ssites
        k = k+1
        if ( i0_buffer(p) ) then
          rbuf = (rjs_central(p)-rcut_H+buffer_H)/buffer_H
          dT_SR(3*(p-1)+1:3*(p-1)+3,:,:,:,:) = dT_SR(3*(p-1)+1:3*(p-1)+3,:,:,:,:) * &
                 (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3) 
          dT_SR(:,3*(p-1)+1:3*(p-1)+3,:,:,:) = dT_SR(:,3*(p-1)+1:3*(p-1)+3,:,:,:) * &
                 (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
        end if
!TEST1
        do j2 = 2, n_sneigh(p)
          k = k+1
          j = sneighbors_list(k)
          q = i_to_p(j)
!TEST2
          if( rjs_H(k) < rcut_H )then
            if ( ij_buffer(k) ) then
              rbuf = (rjs_H(k)-rcut_H+buffer_H)/buffer_H
              dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:,:,:) = dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:,:,:) * &
                  (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
            end if
          end if
        end do
      end do

      do a = 1, n_ssites
        do p = 1, n_ssites
          if ( i0_buffer(p) ) then
            rbuf = (rjs_central(p)-rcut_H+buffer_H)/buffer_H
            fcut = 1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3
            do c3 = 1, 3
              dfcut = -6.d0 * (rbuf-rbuf**2) * 1.d0/buffer_H * xyz_central(c3,p)/rjs_central(p)
              if ( a == 1 ) then
                dT_SR(3*(p-1)+1:3*(p-1)+3,:,a,c3,:) = dT_SR(3*(p-1)+1:3*(p-1)+3,:,a,c3,:) - &
                    dfcut/fcut * T_SR(3*(p-1)+1:3*(p-1)+3,:,:)
                dT_SR(:,3*(p-1)+1:3*(p-1)+3,a,c3,:) = dT_SR(:,3*(p-1)+1:3*(p-1)+3,a,c3,:) - &
                    dfcut/fcut * T_SR(:,3*(p-1)+1:3*(p-1)+3,:)
              end if
              if ( a == p ) then
                dT_SR(3*(p-1)+1:3*(p-1)+3,:,a,c3,:) = dT_SR(3*(p-1)+1:3*(p-1)+3,:,a,c3,:) + &
                    dfcut/fcut * T_SR(3*(p-1)+1:3*(p-1)+3,:,:)
                dT_SR(:,3*(p-1)+1:3*(p-1)+3,a,c3,:) = dT_SR(:,3*(p-1)+1:3*(p-1)+3,a,c3,:) + &
                    dfcut/fcut * T_SR(:,3*(p-1)+1:3*(p-1)+3,:)
              end if
            end do
          end if
        end do
      end do

      do a = 1, n_ssites
        k = 0
        do p = 1, n_ssites
          k = k+1
          do j2 = 2, n_sneigh(p)
            k = k+1
            j = sneighbors_list(k)
            q = i_to_p(j)
            if( rjs_H(k) < rcut_H )then
              if ( ij_buffer(k) ) then
                rbuf = (rjs_H(k)-rcut_H+buffer_H)/buffer_H
                fcut = 1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3
                do c3 = 1, 3
                  dfcut = -6.d0 * (rbuf-rbuf**2) * 1.d0/buffer_H * xyz_H(c3,k)/rjs_H(k)
                  if ( a == p ) then
                    dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,c3,:) = &
                        dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,c3,:) - &
                        dfcut/fcut * T_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:)
                  end if
                  if ( a == q ) then
                    dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,c3,:) = &
                        dT_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a,c3,:) + &
                        dfcut/fcut * T_SR(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,:)
                  end if
                end do
              end if
            end if
          end do
        end do
      end do

      call cpu_time(time2)
      write(*,*) "cutoff function derivatives:", time2-time1
      call cpu_time(time1)

      dB_mat = dB_mat + dT_SR

      call cpu_time(time2)
      write(*,*) "add dB_mat+dT_SR:", time2-time1
      call cpu_time(time1)

      ! end if for do_derivatives
      end if

      allocate( A_mat(1:3*n_ssites,1:3*n_ssites,1:n_freq) )
      A_mat = B_mat

      allocate( ipiv(1:3*n_ssites) )
      allocate( work_arr(1:12*n_ssites) )
      do k3 = 1, n_freq
        call dgetrf(3*n_ssites, 3*n_ssites, A_mat(:,:,k3), 3*n_ssites, ipiv, info)
        call dgetri(3*n_ssites, A_mat(:,:,k3), 3*n_ssites, ipiv, work_arr, 12*n_ssites, info)
      end do

      do k3 = 1, n_freq
!        do p = 1, n_neigh(i)
        A_i = 0.d0
        do q = 1, n_ssites
          A_i = A_i + A_mat(1:3,3*(q-1)+1:3*(q-1)+3,k3)
        end do
        alpha_SCS(i,k3) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
!        end do
      end do

! TEST FOR SURROUNDING POLARIZABILITIES #########################
      write(*,*) "rjs_central", rjs_central
      do k3 = 1, n_freq
        !A_i = 0.d0
        do p = 1, n_ssites
          A_i = 0.d0
          do q = 1, n_ssites
            A_i = A_i + A_mat(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,k3)
          end do
          if (k3 == 1) then
            write(*,*) "alpha SCS for atom:", p_to_i(p), rjs_central(p)*Bohr,  1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
          end if
        end do
      end do
! ###############################################################

      call cpu_time(time2)
      write(*,*) "alpha SCS (timing):", time2-time1
      call cpu_time(time1)

      if ( do_derivatives ) then

      allocate( dA_mat(1:3,1:3*n_ssites,1:n_ssites,1:3,1:n_freq) )
      allocate( AdB_n(1:3,1:3*n_ssites) )
      dA_mat = 0.d0
      do a2 = 1, n_ssites
        do k = 1, n_freq
          do c3 = 1, 3
            AdB_n = 0.d0
            call dgemm('n', 'n', 3, 3*n_ssites, 3*n_ssites, 1.d0, -A_mat(1:3,:,k), 3, &
                       dB_mat(:,:,a2,c3,k), 3*n_ssites, 0.d0, AdB_n, 3)
            call dgemm('n', 'n', 3, 3*n_ssites, 3*n_ssites, 1.d0, AdB_n, 3, A_mat(:,:,k), &
                       3*n_ssites, 0.d0, dA_mat(:,:,a2,c3,k), 3)
          end do
        end do
      end do

      call cpu_time(time2)
      write(*,*) "dA_mat", time2-time1
      call cpu_time(time1)

!      write(*,*) "dB_mat for atom", p_to_i(3)
!      do p = 1, 3*n_ssites
!        write(*,*) dB_mat(p,:,2,1,1)
!      end do
!      write(*,*) dB_mat(1,3*(3-1)+1:3*(3-1)+3,2,1,1)

!      write(*,*) "A_mat"
!      do p = 1, 3*n_ssites
!        write(*,*) A_mat(p,:,1)
!      end do
!      write(*,*) A_mat(1,:,1)

!      write(*,*) "dA_mat"
!      do p = 1, 3
!        write(*,*) dA_mat(p,:,2,1,1)
!      end do

!     Remember that the a2 here are just the permuted indices: p_to_i should be part of the output array
      do a2 = 1, n_ssites
        do k = 1, n_freq
          do c3 = 1, 3
            A_i = 0.d0
            do q = 1, n_ssites
              A_i = A_i + dA_mat(1:3,3*(q-1)+1:3*(q-1)+3,a2,c3,k)
            end do
            do c1 = 1, 3
              dalpha2(n_count+a2,c3,k) = dalpha2(n_count+a2,c3,k) + A_i(c1,c1)
            end do
            dalpha2(n_count+a2,c3,k) = 1.d0/3.d0 * dalpha2(n_count+a2,c3,k)
          end do
        end do
      end do

      end if

      n_ssites_list(i) = n_ssites
      n_sneigh_list2(n_count+1:n_count+n_ssites) = n_sneigh(1:n_ssites)
!TEST1
      write(*,*) "atom:", i
      write(*,*) "alpha_SCS:", alpha_SCS(i,1)

      if ( do_derivatives ) then
      write(*,*) "dalpha:"
      do p = 1, n_ssites
        write(*,*) p_to_i(p), dalpha2(n_count+p,1,1)
      end do
      end if
!      write(*,*) p_to_i(2), dalpha(n_tot+2,1,1)

!      write(*,*) "Central atom:", i
!      k = 0
!      do p = 1, n_neigh(i)
!        if (n_sneigh(p) > 0) then
!          write(*,*) sneighbors_list(k+1:k+n_sneigh(p))
!          k = k + n_sneigh(p)
!        end if
!      end do
      n_count = n_count + n_ssites

      deallocate( local_neighbors, sneighbors_list, rjs_H, xyz_H, hirshfeld_v_neigh_H, &
                  neighbor_species_H, neighbor_c6_ii, r0_ii, neighbor_alpha0, omega_i, alpha_k, sigma_k, f_damp, &
                  g_func, h_func, T_func, T_SR, B_mat, A_mat, ipiv, work_arr, &
                  sneighbors_list_vder, hirshfeld_v_cart_der_H, i0_buffer2, i0_buffer, ij_buffer, rjs_central2, &
                  rjs_central, xyz_central2, xyz_central )
      if ( do_derivatives ) then
        deallocate( dT, dT_SR, f_damp_der, g_func_der, h_func_der, dB_mat, &
                  dT_SR_v, dB_mat_v, coeff_der, coeff_fdamp, dA_mat, AdB_n ) 
      end if

      n_tot = n_tot+n_neigh(i)
    end do

    ! These should be output variables (with alpha_SCS and n_ssites_list)
    if ( do_derivatives ) then
      allocate( dalpha(1:n_count,1:3,1:n_freq) )
      dalpha(1:n_count,:,:) = dalpha2(1:n_count,:,:)
    end if
    allocate( n_sneigh_list(1:n_count) )
    n_sneigh_list = n_sneigh_list2(1:n_count)


    deallocate( omegas, s_i, s_j, sigma_ij, coeff_h_der, terms, n_sneigh, p_to_i, i_to_p, &
                neighbors_list2, hirshfeld_v_cart_der2, rjs2, xyz2, hirshfeld_v_neigh2, neighbor_species2, &
                dg, dh, A_i, neighbors_list3, n_sneigh_vder, ij_buffer2, &
                n_sneigh_list2 )
    if ( do_derivatives ) then
      deallocate( dalpha2 )
    end if
    call cpu_time(time2)
    write(*,*) "the rest (deallocating and writing)", time2-time1
!TEST1

!    write(*,*) "sub neighbor list timing:", time2-time1
    else

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OLD IMPLEMENTATION STARTS HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(*,*) "Start of the old implementation; values are precalculated and matrices"
    write(*,*) "are larger but contain the same number of non-zero elements."
    write(*,*) "Eventually a switch will be implemented to choose which one"
    write(*,*) "to use, or the old one will be scrapped."

!   This implementation assumes that rcut is the largest cutoff, that is, the neigbhbors_list contains only the atoms within the vdW cutoff.
!   The implementation matches with the implementation above if rcut is the largest cutoff or the cutoff is so small that the only neighbor
!   the atoms see are themselves (n_neigh(i) = 1 for all i).

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
!   T_SR SHOULD BE ALLOCATED FOR EACH ATOM SEPARATELY USING THE NUMBER OF NEIGHBORS INSIDE THE RCUT_SCS
!   THE SAME APPLIES TO A LOT OF THESE VARIABLES HERE. N_SITES CAN BE AS SMALL AS 1, DEPENDING ON HOW
!   PARALLELISM IS HANDLED
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
    allocate( dalpha_full(1:n_sites,1:n_sites,1:3,1:n_freq) )
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
    allocate( AdB_n(1:3*n_sites,1:3*n_sites) )
    end if
    allocate( alpha_k(1:n_pairs,1:n_freq) )
    allocate( sigma_k(1:n_pairs,1:n_freq) )
    allocate( s_i(1:n_freq) )
    allocate( s_j(1:n_freq) )

    if( do_timing) then
      call cpu_time(time1)
    end if

!    write(*,*) "hirshfeld_v:"
!    do p = 1, n_sites
!      write(*,*) p, hirshfeld_v(p)
!    end do
    
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

    rcut_H = rcut/Bohr
    buffer_H = buffer/Bohr

!   Precompute some other pair quantities
    neighbor_c6_ii = neighbor_c6_ii * hirshfeld_v_neigh**2
!   This is slow, could replace by Taylor expansion maybe
    r0_ii = r0_ii * hirshfeld_v_neigh**(1.d0/3.d0)
    neighbor_alpha0 = neighbor_alpha0 * hirshfeld_v_neigh
    omega_i = (4.d0 * neighbor_c6_ii)/(3.d0*neighbor_alpha0**2)

!    write(*,*) "omega_i:", omega_i

    k2 = 1
    do i = 1, n_sites
      do k = 1, n_freq
        alpha_i(i,k) = neighbor_alpha0(k2)/(1.d0 + omegas(k)**2/omega_i(k2)**2)
        sigma_i(i,k) = (sqrt(2.d0/pi) * alpha_i(i,k)/3.d0)**(1.d0/3.d0)
      end do
      k2 = k2+n_neigh(i)
    end do

!    write(*,*) "alpha_i", alpha_i(:,1)

    do k = 1, n_freq
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
! FDAMP TEST:
!          f_damp(k) = 1.d0/( 1.d0 + exp( -d*( rjs_H(k)/((0.5d0*rcut)) - 1.d0 ) ) )
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
          !write(*,*) "sigma_ij:", sigma_ij
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
!REGULARIZATION TEST:
!        B_mat(3*(i-1)+c1,3*(i-1)+c1,:) = alpha_k(k,:)
        nnz = nnz+1
        ia(nnz) = 3*(i-1)+c1
        ja(nnz) = 3*(i-1)+c1
        val(nnz) = B_mat(3*(i-1)+c1,3*(i-1)+c1,1)
        !write(*,*) "i, nnz", i, nnz
      end do
      do j2 = 2, n_neigh(i)
        k = k+1
        j = neighbors_list(k)
        j = modulo(j-1,n_sites)+1
        if( rjs(k) < rcut )then
!          write(*,*) i, j
!          if (i == 8 .and. j == 1) then
!            write(*,*) "i, j", i, j
!            write(*,*) "f_damp", f_damp(k)
!            write(*,*) "T_func", T_func(9*(k-1)+1:9*(k-1)+9)
!            write(*,*) "g_func", g_func(k,1)
!            write(*,*) "h_func", h_func(9*(k-1)+1:9*(k-1)+9,1)
!          end if
!          if (i == 1 .and. j == 8) then
!            write(*,*) "i, j", i, j
!            write(*,*) "f_damp", f_damp(k)
!            write(*,*) "T_func", T_func(9*(k-1)+1:9*(k-1)+9)
!            write(*,*) "g_func", g_func(k,1)
!            write(*,*) "h_func", h_func(9*(k-1)+1:9*(k-1)+9,1)
!          end if
          k2 = 9*(k-1)
          do c1 = 1, 3
            do c2 = 1, 3
              k2 = k2+1
              !if ( iterative ) then
                nnz = nnz+1
                !write(*,*) "j, nnz", j, nnz
                !if ( rjs(k) > rcut-buffer ) then
                !  rbuf = (rjs_H(k)-rcut_H+buffer_H)/buffer_H
                !  T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = alpha_i(i,:) * (1.d0-f_damp(k)) * (-T_func(k2) * &
                !                                  g_func(k,:) + h_func(k2,:)) * (1.d0 - 3.d0 * rbuf**2 + 2.d0 * rbuf**3)
                !else
                do k3 = 1, n_freq
                  T_SR(3*(i-1)+c1,3*(j-1)+c2,k3) = alpha_i(i,k3) * (1.d0-f_damp(k)) * (-T_func(k2) * &
                                                  g_func(k,k3) + h_func(k2,k3))
                end do
!              T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(k2) * &
!                                              g_func(k,:) + h_func(k2,:))
                !end if
                ia(nnz) = 3*(i-1)+c1
                ja(nnz) = 3*(j-1)+c2
                val(nnz) = T_SR(3*(i-1)+c1,3*(j-1)+c2,1)
              !else
              !  T_SR(3*(i-1)+c1,3*(j-1)+c2,:) = (1.d0-f_damp(k)) * (-T_func(k2) * &
              !                                  g_func(k,:) + h_func(k2,:))
              !end if
            end do
          end do
        end if
      end do
    end do

    write(*,*) "alpha_i", size(alpha_i,1)
    write(*,*) "hirshfeld_v", size(hirshfeld_v,1)
    write(*,*) "T_SR:"
    do i = 1, 6
      write(*,*) T_SR(i,1:6,1)
    end do

    !write(*,*) "nnz", nnz

    B_mat = B_mat + T_SR

!REGULARIZATION TEST
!    do k2 = 1, n_freq
!    B_mat(:,:,k2) = matmul(B_mat(:,:,k2),T_SR(:,:,k2))
!    end do

!    write(*,*) "B_mat"
!    do p = 1, 3*n_sites
!      write(*,*) B_mat(p,:,1)
!    end do


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

!      write(*,*) "coeff_fdamp", coeff_fdamp(1,3,1,1,1)
!      write(*,*) "coeff_der", coeff_der(1,3,1,1,1)

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
                        coeff_der(i,j,c1,c2,k) * (sigma_i(i,k)**2/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2))
                      dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) = &
                        dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k) + &
                        coeff_fdamp(i,j,c1,c2,k) * r_vdw_i/hirshfeld_v(i) * hirshfeld_v_cart_der_H(c3,k2)
!                        if ( a == 12 .and. c3 == 1 .and. k == 1 .and. i == 1 .and. j == 3 .and. c1 == 1 .and. c2 == 1) then
!                          write(*,*) "r_vdw_i", r_vdw_i
!                          write(*,*) "hirshfeld_v_neigh_H", hirshfeld_v(i)
!                          write(*,*) "sigma_k", sigma_i(i,k)
!                          write(*,*) "hirshfeld_v_cart_der_H", hirshfeld_v_cart_der_H(c3,k2)
!                          write(*,*) "dT_SR_v", dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k)
!                        end if
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
            do k = 1, n_freq
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
!                        if ( a == 12 .and. c3 == 1 .and. k == 1 .and. i == 1 .and. j == 3 .and. c1 == 1 .and. c2 == 1) then
!                          write(*,*) "r_vdw_j", r_vdw_j
!                          write(*,*) "hirshfeld_v_neigh_H", hirshfeld_v(j)
!                          write(*,*) "sigma_k", sigma_i(j,k)
!                          write(*,*) "hirshfeld_v_cart_der_H", hirshfeld_v_cart_der_H(c3,k2)
!                          write(*,*) "dT_SR_v", dT_SR_v(3*(i-1)+c1,3*(j-1)+c2,a,c3,k)
!                        end if
                    end do
                  end do
                end if
              end do
            end do
          end do
        end do
        k3 = k3 + n_neigh(j)
      end do

    !if ( iterative ) then
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
              if ( i == 1 .and. a == 1 .and. k == 1 .and. c3 == 1) then
                write(*,*) "hirshfeld_v", hirshfeld_v(i)
                write(*,*) "dv/dr", hirshfeld_v_cart_der_H(c3,k2)
              end if
            end do
          end do
        end do
      end do
    !end if

    dT_SR = dT_SR + dT_SR_v
    end if

    !if ( iterative ) then
      do i = 1, n_sites
        do k = 1, n_freq
          dT_SR(3*(i-1)+1:3*(i-1)+3,:,:,:,k) = alpha_i(i,k) * dT_SR(3*(i-1)+1:3*(i-1)+3,:,:,:,k) 
        end do
      end do
    !end if

    write(*,*) "dT_SR:"
    do i = 1, 6
      write(*,*) dT_SR(i,1:6,1,1,1)
    end do

    end if ! do_derivatives

    
!    A_mat = B_mat

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


      n_iter = 100
      allocate( I_mat(1:3*n_sites,1:3*n_sites) )
      !allocate( I_aT(1:3*n_sites,1:3*n_sites,1:n_freq) )
      allocate( a_vec(1:3*n_sites,1:3,1:n_freq) )
      !allocate( regularization(1:3*n_sites,1:3*n_sites) )
      I_mat = 0.d0
!      regularization = 0.d0
      do p = 1, n_sites
        do c1 = 1, 3
          I_mat(3*(p-1)+c1,3*(p-1)+c1) = 1.d0
!          regularization(3*(p-1)+c1,3*(p-1)+c1) = 0.0d0
        end do
      end do
      !I_aT = 0.d0
      a_vec = 0.d0
      do k = 1, n_freq
        !I_aT(:,:,k) = I_mat + T_SR(:,:,k)
        do p = 1, n_sites
          do c1 = 1, 3
            a_vec(3*(p-1)+c1,c1,k) = alpha_i(p,k)
          end do
        end do
      end do
      !write(*,*) "a_TS:"
      !do p = 1, 3*n_sites
      !  write(*,*) a_next(p,:,1)
      !end do

      allocate( a_SCS(1:3*n_sites,1:3,1:n_freq) )
      allocate( BTB_reg_copy(1:3*n_sites,1:3*n_sites) )
      !allocate( I_aT2(1:3*n_sites,1:3*n_sites) )

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
!        if (k == 1) then
!          write(*,*) "I_aT2 before"
!          do p = 1, 6
!            write(*,*) I_aT2(p,1:6)
!          end do
!        end if
        call dsysv('U', 3*n_sites, 3, BTB_reg_copy, 3*n_sites, ipiv, &
                    a_SCS(:,:,k), 3*n_sites, work_arr, 12*n_sites, info)
!        if (k == 1) then
!          write(*,*) "I_aT2 after"
!          do p = 1, 6
!            write(*,*) I_aT2(p,1:6)
!          end do
!        end if  
      end do

      !deallocate( BTB_reg, B_reg, BTB_reg_copy )
      !deallocate( I_aT_a_vec, I_aT2 )

      alpha_SCS0 = 0.d0
      do k = 1, n_freq
        do p = 1, n_sites
          do c1 = 1, 3
            alpha_SCS0(p,k) = alpha_SCS0(p,k) + a_SCS(3*(p-1)+c1,c1,k)
          end do
        end do
      end do
      alpha_SCS0 = alpha_SCS0/3.d0

      write(*,*) "a_SCS"
      do i = 1, n_sites
        write(*,*) a_SCS(i,:,1)
      end do

!      deallocate( BTB_reg, B_reg, BTB_reg_copy, I_mat, a_vec, a_SCS )
    !else
!      deallocate( ia, ja, val)
!      allocate( a_vec(1:3*n_sites,1:3,1:n_freq) )
!
!      a_vec = 0.d0
!      do k = 1, n_freq
!        do p = 1, n_sites
!          do c1 = 1, 3
!            a_vec(3*(p-1)+c1,c1,k) = 1.d0
!          end do
!        end do
!      end do
!      do k = 1, n_freq
!        call dsysv('U', 3*n_sites, 3, B_mat(:,:,k), 3*n_sites, ipiv, &
!                    a_vec(:,:,k), 3*n_sites, work_arr, 12*n_sites, info)
!      end do
!      alpha_SCS0 = 0.d0
!      do k = 1, n_freq
!        do p = 1, n_sites
!          do c1 = 1, 3
!            alpha_SCS0(p,k) = alpha_SCS0(p,k) + a_vec(3*(p-1)+c1,c1,k)
!          end do
!        end do
!      end do
!      alpha_SCS0 = alpha_SCS0/3.d0
      !write(*,*) "alpha_SCS0"
      !do p = 1, n_sites
      !  write(*,*) alpha_SCS0(p,1)
      !end do

!      deallocate( a_vec )

!      A_mat = B_mat
!      do k3 = 1, n_freq
!        call dgetrf(3*n_sites, 3*n_sites, A_mat(:,:,k3), 3*n_sites, ipiv, info)
!        call dgetri(3*n_sites, A_mat(:,:,k3), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
!      end do

!      alpha_SCS0 = 0.d0
!      do k3 = 1, n_freq
!        do p = 1, n_sites
!          A_i = 0.d0
!          do q = 1, n_sites
!            A_i = A_i + A_mat(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,k3)
!          end do
!          alpha_SCS0(p,k3) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
!        end do
!      end do
!    end if

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

      !write(*,*) "da_vec"
      !do i = 1, 3*n_sites
      !  write(*,*) da_vec(i,:,1,1,1)
      !end do

      allocate( dBTB(1:3*n_sites,1:3*n_sites) )
      allocate( vect1(1:3*n_sites,1:3) )
      allocate( vect2(1:3*n_sites,1:3) )
      allocate( vect3(1:3*n_sites,1:3) )
      allocate( da_SCS(1:3*n_sites,1:3,1:n_sites,1:3,1:n_freq) )

      ! What happens here: regularized equation
      ! (B^T * B + reg_param * I) a_SCS = (B^T + reg_param*I) a_vec
      ! is differentiated:
      ! (B^T * B + reg_param * I) a_SCS' = (B')^T a_vec + (B^T + reg_param*I) a_vec' - ((B')^T * B + B^T * B') a_SCS
      !                                    |---vect1--|   |---------vect2----------|   |--------vect3--------------|
      ! a_SCS' is solved with dsysv.

      write(*,*) "Doing: freq, atom"
      do k = 1, n_freq
        do a = 1, n_sites
          write(*,*) k, a
          do c3 = 1, 3
            call dgemm('t', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, &
              B_mat(:,:,k), 3*n_sites, 0.d0, dBTB, 3*n_sites)
            call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dT_SR(:,:,a,c3,k), 3*n_sites, &
              a_vec(:,:,k), 3*n_sites, 0.d0, vect1, 3*n_sites)
            call dgemm('t', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, B_reg(:,:,k), 3*n_sites, &
              da_vec(:,:,a,c3,k), 3*n_sites, 0.d0, vect2, 3*n_sites)
            call dgemm('n', 'n', 3*n_sites, 3, 3*n_sites, 1.d0, dBTB+transpose(dBTB), 3*n_sites, &
              a_SCS(:,:,k), 3*n_sites, 0.d0, vect3, 3*n_sites)
            da_SCS(:,:,a,c3,k) = vect1 + vect2 - vect3
            BTB_reg_copy = BTB_reg(:,:,k)
            call dsysv('U', 3*n_sites, 3, BTB_reg_copy, 3*n_sites, ipiv, &
              da_SCS(:,:,a,c3,k), 3*n_sites, work_arr, 12*n_sites, info)
          end do
        end do
      end do

      write(*,*) "da_SCS"
      do p = 1, n_sites
        write(*,*) da_SCS(p,:,1,1,1)
      end do

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



      end if

      deallocate( BTB_reg, B_reg, BTB_reg_copy, I_mat, a_vec, a_SCS, da_vec, dBTB, vect1, vect2, &
                  vect3, da_SCS )

!      dA_mat = 0.d0
!      do a2 = 1, n_sites
!        do k = 1, n_freq
!          do c3 = 1, 3
!            AdB_n = 0.d0
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, -A_mat(:,:,k), 3*n_sites, &
!                       dB_mat(:,:,a2,c3,k), 3*n_sites, 0.d0, AdB_n, 3*n_sites)
!            call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, AdB_n, 3*n_sites, A_mat(:,:,k), &
!                       3*n_sites, 0.d0, dA_mat(:,:,a2,c3,k), 3*n_sites)
!          end do
!        end do
!      end do
!
!      dalpha_full = 0.d0
!      do a2 = 1, n_sites
!        do k = 1, n_freq
!          do c3 = 1, 3
!            do p = 1, n_sites
!              A_i = 0.d0
!              do q = 1, n_sites
!                A_i = A_i + dA_mat(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,a2,c3,k)
!              end do
!              do c1 = 1, 3
!                dalpha_full(p,a2,c3,k) = dalpha_full(p,a2,c3,k) + A_i(c1,c1)
!              end do
!              dalpha_full(p,a2,c3,k) = 1.d0/3.d0 * dalpha_full(p,a2,c3,k)
!            end do
!          end do
!        end do
!      end do
!
!      end if
! TEST1
      write(*,*) "alpha_SCS:" !,  alpha_SCS0(1,1)
      do p = 1, n_sites
        write(*,*) p, alpha_SCS0(p,1)
      end do


      if ( do_derivatives ) then
      write(*,*) "dalpha_SCS w.r.t. atom 1 in x direction:"
      do p = 1, n_sites
        write(*,*) p, dalpha_full(p,1,1,1)
      end do
!      write(*,*) neighbors_list(n_tot+9), dalpha_n(n_tot+9,1,1)
      end if

    ! TODO: This subroutine should output alpha_SCS_i and dalpha_n: they can be read using the neighbors_list
    ! The alpha_SCS_i should probably be reduced to the size of n_sites (only the SCS polarizability of each
    ! central atom).

!   Clean up
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, alpha_k, sigma_k, s_i, s_j, )
    if ( do_derivatives ) then
    deallocate( dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                dB_mat, dB_mat_v, dalpha_full, coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H, AdB_n )
    end if

  if ( .false. ) then
    n_tot = 0
    allocate( alpha_SCS_i(1:n_pairs,1:11) )
    allocate( dalpha_n(1:n_pairs,1:3,1:11) )
    alpha_SCS_i = 0.d0
    dalpha_n = 0.d0
! TEST1
    do i = 1, n_sites
!    do i = 1, n_sites
!      write(*,*) "allocate"
      allocate( T_SR_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
      allocate( B_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
      allocate( A_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
!      allocate( alpha_SCS_i(1:n_neigh(i),1:11) )
      allocate( dB_mat_n(1:3*n_neigh(i),1:3*n_neigh(i),1:n_neigh(i),1:3,1:11) )
      allocate( dA_mat_n(1:3,1:3*n_neigh(i),1:n_neigh(i),1:3,1:11) )
      allocate( AdB_n(1:3,1:3*n_neigh(i)) )
!      allocate( dalpha_n(1:n_neigh(i),1:3,1:11) )
      T_SR_i = 0.d0
      B_mat_i = 0.d0
      A_mat_i = 0.d0
      alpha_SCS_i = 0.d0

!      write(*,*) "rjs"

!      do p = 1, n_neigh(i)
!        write(*,*) rjs(n_tot+p)
!      end do


      k = 0
      ! Changing the order and logic in some of the loops here might make a difference. Have to test.
      ! The problem is having to go through the entire neighbors list for each atom (to get interactions
      ! between neighbors correctly, the interactions with central atom are straightforward).
      do i2 = 1, n_sites
        if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
          p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
!          write(*,*) "i2, p, neighbors_list(n_tot+p)", i2, p, neighbors_list(n_tot+p)
!          if ( rjs(n_tot+p) < rcut ) then
!            write(*,*) p, alpha_k(n_tot+p,:)
          k = k+1
          do c1 = 1, 3
            B_mat_i(3*(p-1)+c1,3*(p-1)+c1,:) = 1.d0/alpha_k(n_tot+p,:)
          end do
          do j2 = 2, n_neigh(i2)
            k = k+1
            j = neighbors_list(k)
!              write(*,*) "i2, j, rjs(k)", i2, j, rjs(k)
            if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
              q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
              if ( rjs(n_tot+q) < rcut) then
!                write(*,*) "q, rjs(k):", q, rjs(k)
                if ( rjs(k) < rcut) then
!                  write(*,*) "q:", q
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
            end if
          end do
!          else
!            k = k+n_neigh(i2)
!          end if
        else
          k = k+n_neigh(i2)
        end if
      end do
      B_mat_i = B_mat_i + T_SR_i
      A_mat_i = B_mat_i

!      write(*,*) "B_mat_i:"
!      do p = 1, 6
!        write(*,*) B_mat_i(p,1:6,1)
!      end do


      do k3 = 1, 11
        call dgetrf(3*n_neigh(i), 3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, info)
!        write(*,*) "dgetrf:", k3, info
        call dgetri(3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, work_arr, 12*n_sites, info)
!        write(*,*) "dgetri:", k3, info
      end do

!      write(*,*) "A_mat_i:"
!      do p = 1, 6
!        write(*,*) A_mat_i(p,1:6,1)
!      end do
        
      ! Have to check, but probably only the central one is needed (p = 1)
      do k3 = 1, 11
        do p = 1, n_neigh(i)
          A_i = 0.d0
          do q = 1, n_neigh(i)
            A_i = A_i + A_mat_i(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,k3)
          end do
          alpha_SCS_i(n_tot+p,k3) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
        end do
      end do

!      write(*,*) "dB_mat_n"
      dB_mat_n = 0.d0
      do a2 = 1, n_neigh(i)
        a = neighbors_list(n_tot+a2)
!        write(*,*) "a2, a:", a2, a
        k = 0
        do i2 = 1, n_sites
          if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == i2) ) then
            k = k+1
            p = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),i2,1)
!            write(*,*) p, i2
            do c1 = 1, 3
              do c3 = 1, 3
                dB_mat_n(3*(p-1)+c1,3*(p-1)+c1,a2,c3,:) = dB_mat(3*(i2-1)+c1,3*(i2-1)+c1,a,c3,:)
              end do
            end do
            do j2 = 2, n_neigh(i2)
              k = k+1
              j = neighbors_list(k)
!              write(*,*) "i2, j, rjs(k)", i2, j, rjs(k)
              if ( any(neighbors_list(n_tot+1:n_tot+n_neigh(i)) == j) ) then
                q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
                if (rjs(n_tot+q) < rcut) then
                  if (rjs(k) < rcut) then
!                  write(*,*) "i2, j:", i2, j
                    q = findloc(neighbors_list(n_tot+1:n_tot+n_neigh(i)),j,1)
                    do c1 = 1, 3
                      do c2 = 1, 3
                        do c3 = 1, 3
                          dB_mat_n(3*(p-1)+c1,3*(q-1)+c2,a2,c3,:) = dB_mat(3*(i2-1)+c1,3*(j-1)+c2,a,c3,:)
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
      end do

!      write(*,*) "dA_mat_n"
      dA_mat_n = 0.d0
      do a2 = 1, n_neigh(i)
        do k = 1, 11
          do c3 = 1, 3
            AdB_n = 0.d0
            call dgemm('n', 'n', 3, 3*n_neigh(i), 3*n_neigh(i), 1.d0, -A_mat_i(1:3,:,k), 3, &
                       dB_mat_n(:,:,a2,c3,k), 3*n_neigh(i), 0.d0, AdB_n, 3)
            call dgemm('n', 'n', 3, 3*n_neigh(i), 3*n_neigh(i), 1.d0, AdB_n, 3, A_mat_i(:,:,k), &
                       3*n_neigh(i), 0.d0, dA_mat_n(:,:,a2,c3,k), 3)
          end do
        end do
      end do

!      write(*,*) "dB_mat for atom:", neighbors_list(n_tot+12)
!      do p = 1, 3*n_neigh(i)
!        write(*,*) dB_mat_n(p,:,9,1,1)
!      end do
!      write(*,*) dB_mat_n(1,3*(12-1)+1:3*(12-1)+3,9,1,1)

!      write(*,*) "A_mat"
!      do p = 1, 3*n_neigh(i)
!        write(*,*) A_mat_i(p,:,1)
!      end do
!      write(*,*) A_mat_i(1,:,1)

!      write(*,*) "dA_mat"
!      do p = 1, 3
!        write(*,*) dA_mat_n(p,:,9,1,1)
!      end do

!      write(*,*) "dalpha_n"
!      dalpha_n = 0.d0
      do a2 = 1, n_neigh(i)
        do k = 1, 11
          do c3 = 1, 3
            A_i = 0.d0
            do q = 1, n_neigh(i)
              A_i = A_i + dA_mat_n(1:3,3*(q-1)+1:3*(q-1)+3,a2,c3,k)
            end do
            do c1 = 1, 3
              dalpha_n(n_tot+a2,c3,k) = dalpha_n(n_tot+a2,c3,k) + A_i(c1,c1)
            end do
            dalpha_n(n_tot+a2,c3,k) = 1.d0/3.d0 * dalpha_n(n_tot+a2,c3,k)
          end do
        end do
      end do
! TEST1
      write(*,*) "atom:", i
      ! SCS polarizability of the central atom
      write(*,*) "alpha_SCS:", alpha_SCS_i(n_tot+1,1)
      write(*,*) "dalpha_SCS:"
      do p = 1, n_neigh(i)
        write(*,*) neighbors_list(n_tot+p), dalpha_n(n_tot+p,1,1)
      end do
!      write(*,*) neighbors_list(n_tot+9), dalpha_n(n_tot+9,1,1)



      n_tot = n_tot + n_neigh(i)
      
!      write(*,*) "deallocate"
      ! MEMORY CORRUPTION HERE: double free or corruption
      deallocate( T_SR_i, B_mat_i, A_mat_i, dB_mat_n, dA_mat_n, AdB_n )
!      write(*,*) "end"
    end do
    ! TODO: This subroutine should output alpha_SCS_i and dalpha_n: they can be read using the neighbors_list
    ! The alpha_SCS_i should probably be reduced to the size of n_sites (only the SCS polarizability of each
    ! central atom).
    deallocate( alpha_SCS_i, dalpha_n )

!   Clean up
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, &
                dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                dB_mat, dB_mat_v, &
                coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H, alpha_k, sigma_k, s_i, s_j )
    
    end if
    
    end if

  if ( new_implementation ) then
    deallocate( alpha_SCS, n_ssites_list, n_sneigh_list )
    if ( do_derivatives ) then
      deallocate( dalpha )
    end if
  end if

  end subroutine

!**************************************************************************

!NOTE: T_func and dT remain unchanged. Should they be directly passed to this subroutine?

  subroutine get_mbd( alpha_SCS0, n_neigh, neighbors_list, neighbor_species, &
                      rcut, buffer, rcut_inner, buffer_inner, rjs, xyz, &
                      sR, d, c6_ref, r0_ref, alpha0_ref, do_forces, &
                      energies, forces0, virial)

    implicit none
    real*8, intent(inout) :: alpha_SCS0(:,:), energies(:), forces0(:,:)
    real*8, intent(in) :: rcut, buffer, rcut_inner, buffer_inner, &
                          rjs(:), xyz(:,:), sR, d, c6_ref(:), r0_ref(:), &
                          alpha0_ref(:)
    real*8, intent(out) :: virial(1:3, 1:3)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    logical, intent(in) :: do_forces
    real*8, allocatable :: A_LR(:,:,:), T_LR(:,:), r0_ii_SCS(:), r0_ii(:), neighbor_alpha0(:), &
                           xyz_H(:,:), rjs_H(:), f_damp_SCS(:), T_func(:), AT(:,:,:), AT_n(:,:,:,:), &
                           energy_series(:,:,:), integrand(:,:), omegas(:)
    integer :: n_order, n_freq, n_sites, n_pairs, n_species, n_sites0
    integer :: k, k2, i, j, j2, c1, c2
    real*8 :: Bohr, Hartree, pi, r_vdw_i, r_vdw_j, E_MBD, integral, omega
    logical :: series_average = .true.

!   Hartree units (calculations done in Hartree units for simplicity)
    Bohr = 0.5291772105638411
    Hartree = 27.211386024367243
    pi = acos(-1.d0)

    n_order = 8

    n_freq = size(alpha_SCS0, 2)
    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_species = size(c6_ref)
    n_sites0 = size(forces0, 2)

    allocate( omegas(1:n_freq) )
    allocate( xyz_H(1:3,1:n_pairs) )
    allocate( rjs_H(1:n_pairs) )
    allocate( r0_ii(1:n_pairs) )
    allocate( neighbor_alpha0(1:n_pairs) )
    allocate( T_func(1:9*n_pairs) )
    allocate( A_LR(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( T_LR(1:3*n_sites,1:3*n_sites) )
    allocate( r0_ii_SCS(1:n_pairs) )
    allocate( f_damp_SCS(1:n_pairs) )
    allocate( AT(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( AT_n(1:3*n_sites,1:3*n_sites, 1:n_order-1, 1:n_freq) )
    allocate( energy_series(1:3*n_sites,1:3*n_sites,1:n_freq) )
    allocate( integrand(1:n_sites,1:n_freq) )

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

    A_LR = 0.d0

    do k = 1, n_freq
      do i = 1, n_sites
        do c1 = 1, 3
          A_LR(3*(i-1)+c1,3*(i-1)+c1,k) = alpha_SCS0(i,k)
        end do
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
            end do
          end do
        end if
      end do
    end do

    AT = 0.d0
    do k = 1, n_freq
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR(:,:,k), 3*n_sites, &
                 T_LR, 3*n_sites, 0.d0, AT(:,:,k), 3*n_sites)
    end do

    AT_n = 0.d0
    integrand = 0.d0
    energy_series = -0.5d0 * AT(:,:,:)
    do k = 1, n_freq
      AT_n(:,:,1,k) = AT(:,:,k)
      do k2 = 1, n_order-2
        ! Precalculate the full AT_n for forces:
        call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, AT(:,:,k), 3*n_sites, &
                   AT_n(:,:,k2,k), 3*n_sites, 0.d0, AT_n(:,:,k2+1,k), 3*n_sites)
        if (series_average .and. k2 == n_order-2) then
          energy_series(:,:,k) = energy_series(:,:,k) - 1.d0/(k2+2)*AT_n(:,:,k2+1,k)/2.d0
        else
          energy_series(:,:,k) = energy_series(:,:,k) - 1.d0/(k2+2)*AT_n(:,:,k2+1,k)
        end if
      end do
      do i = 1, n_sites
        do c1 = 1, 3
          integrand(i,k) = integrand(i,k) + alpha_SCS0(i,k) * dot_product(T_LR(3*(i-1)+c1,:), &
                         energy_series(:,3*(i-1)+c1,k))
        end do
      end do
    end do
    do i = 1, n_sites
      write(*,*) "integrand", i, integrand(i,:)
    end do

    do i = 1, n_sites
      integral = 0.d0
      call integrate("trapezoidal", omegas, integrand(i,:), 0.d0, 10.d0, integral)
      E_MBD = integral/(2.d0*pi)
      energies(i) = E_MBD * 27.211386245988
      write(*,*) "local energy", i, energies(i)
    end do
    write(*,*) "MBD energy:", sum(energies)

    deallocate( omegas, xyz_H, rjs_H, r0_ii, neighbor_alpha0, T_func, A_LR, T_LR, r0_ii_SCS, &
                f_damp_SCS, AT, AT_n, energy_series, integrand )

  end subroutine

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
                           coeff_der(:,:,:,:,:), coeff_fdamp(:,:,:,:,:), dg(:), &
                           dh(:), hirshfeld_v_cart_der_H(:,:), &
                           integrand_k(:), E_MBD_k(:), alpha_k(:,:), sigma_k(:,:), s_i(:), s_j(:), &
                           T_SR_i(:,:,:), B_mat_i(:,:,:), alpha_SCS_i(:,:), A_mat_i(:,:,:), A_LR_i(:,:,:), &
                           T_LR_i(:,:), AT_i(:,:,:), energy_series(:,:,:), I_mat_n(:,:), &
                           logIAT_n(:,:,:), WR_n(:), WI_n(:), VL_n(:,:), VR_n(:,:), VRinv_n(:,:), &
                           logMapo_n(:,:), dB_mat_n(:,:,:,:), dA_mat_n(:,:,:,:), dBA_n(:,:), &
                           dA_LR_n(:,:,:,:), dalpha_n(:,:,:), f_damp_der_ij_n(:), f_damp_der_SCS_ij_n(:), &
                           dT_LR_n(:,:,:), force_integrand_n(:,:), forces_MBD_k(:,:), invIAT_n(:,:,:), &
                           G_mat_n(:,:,:,:), AT_n(:,:,:,:), force_series(:,:,:), energy_term(:,:)
    real*8 :: time1, time2, this_force(1:3), Bohr, Hartree, &
              omega, pi, integral, E_MBD, R_vdW_ij, R_vdW_SCS_ij, S_vdW_ij, dS_vdW_ij, exp_term, &
              rcut_vdw, r_vdw_i, r_vdw_j, dist, f_damp_SCS_ij, t1, t2
    integer, allocatable :: ipiv(:)
    integer :: n_sites, n_pairs, n_species, n_sites0, info, n_order, n_tot
    integer :: i, i2, j, j2, k, k2, k3, a, a2, c1, c2, c3, lwork, b, p, q
    logical :: do_timing = .false., do_hirshfeld_gradients = .false., nonlocal = .false., &
               series_expansion = .true., total_energy = .true., series_average = .true.

! WE SHOULD ENFORCE RCUT_MBD, RCUT_SCS <= RCUT - BUFFER
! 4.5, 4.5, 0.5, 0.5; make sure 
real*8 :: rcut_scs, rcut_mbd, buffer_scs, buffer_mbd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TO DO: SPLIT THIS SUBROUTINE INTO TWO (OR, MAYBE HAVE A KEYWORD TO SELECT IN WHICH "MODE" THE SUBROUTINE IS
!        BEING CALLED, THAT WAY WE CAN USE SAVE STATEMENTS, E.G., FOR THE PAIR QUANTITIES):
! 1) THE FIRST ONE IS USED TO COMPUTE THE SELF-CONSISTENT POLARIZABILITIES. THESE ARE REDUCED/BROADCASTED BY
!    RANK 0 AND PASSED ON TO THE SECOND SUBROUTINE SO THAT ALL THE ATOMS ARE CHARACTERIZED BY THE MORE ACCURATE
!    POLARIZABILIES RATHER THAN TS POLARIZABILITIES. THIS PART SHOULD BE REASONABLY CHEAP (COMPARABLE TO
!    BROADCASTING THE HIRSHFELD VOLUMES IN THE TS IMPLEMENTATION)
! 2) THE SECOND ONE DOES THE SCREENED DIPOLE-DIPOLE COUPLING CALCULATION USING THE SCS POLARIZABILITIES
!
! NOTES:
! *) THE SCS POLARIZABILITIES HAVE THEIR OWN CUTOFF RADIUS AND CUTOFF FUNCTION, WHICH IS ENFORCED THROUGH THE
!    SHORT-RANGE DIPOLE INTERACTION TENSOR CUTOFF
! *) T_LR SHOULD HAVE A CUTOFF FOR FULL COUPLING (ALL ATOMS INTERACT) AND ANOTHER FOR PAIR INTERACTIONS (ONLY
!    CENTRAL ATOM AND OTHER ONE INTERACT; THERE ARE THE T_1J = T_J1 TERMS). THIS WILL MAKE T_LR SPARSE.
!    HOPEFULLY, ACCURATE MODELS CAN BE OBTAINED WITH FULL COUPLING CUTOFFS MUCH SMALLER THAN OVERALL CUTOFFS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !There are some flags and other parameters that should eventually be made function arguments:
    !do_hirshfeld_gradients = .true. -> If this is true, the Hirshfeld gradients are included in the forces
    !nonlocal = .false. -> If this is true, the full MBD energy and forces are calculated for the input .xyz 
    !                      (instead of the local version)
    !series_expansion = .false. -> If this is true, series expansion is used for the local stuff, if false, 
    !                              the local energies are not that well defined (just 1/N_neigh of the 
    !                              entire cutoff sphere)
    !total_energy = .true. -> This is used to calculate the total energies for the cutoff spheres; mostly 
    !                         used for doing finite difference for the forces
    !series_average = .true. -> If this is true, the series expansion is averaged for n_order and
    !                           n_order-1, which should give better results because the series expansion
    !                           seems to be a Cauchy sequence of n_order
    !
    !n_order -> This gives the order of the series expansion for the local part and has to be at least 2


!   Change these to be input variables (NOTE THAT THEY ARE IN ANGSTROMS!):
!    write(*,*) "rcut", rcut
!    rcut_vdw = 4.d0
    ! n_order has to be at least 2
!    n_order = 6
    n_order = 8

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
!   T_SR SHOULD BE ALLOCATED FOR EACH ATOM SEPARATELY USING THE NUMBER OF NEIGHBORS INSIDE THE RCUT_SCS
!   THE SAME APPLIES TO A LOT OF THESE VARIABLES HERE. N_SITES CAN BE AS SMALL AS 1, DEPENDING ON HOW
!   PARALLELISM IS HANDLED
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
    allocate( coeff_der(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( coeff_fdamp(1:n_sites,1:n_sites,1:3,1:3,1:11) )
    allocate( dg(1:11) )
    allocate( dh(1:11) )
    allocate( hirshfeld_v_cart_der_H(1:3,1:n_pairs) )
    allocate( integrand_k(1:11) )
    allocate( E_MBD_k(1:n_sites) )
    allocate( alpha_k(1:n_pairs,1:11) )
    allocate( sigma_k(1:n_pairs,1:11) )
    allocate( s_i(1:11) )
    allocate( s_j(1:11) )
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

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "Computing basic pair quantities:", t2-t1, "seconds"

!   Computing dipole interaction tensor
!   Requires the complete supercell to get correct dimensions for T_func! 3*N_at x 3*N_at, where N_at are atoms in supercell
call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "Computing dipole-dipole tensor:", t2-t1, "seconds"

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "g and h functions, whatever that is:", t2-t1, "seconds"

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "computing T_SR, A_mat, B_mat:", t2-t1, "seconds"

call cpu_time(t1)
    do i = 1, 11
      call dgetrf(3*n_sites, 3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, info)
      call dgetri(3*n_sites, A_mat(:,:,i), 3*n_sites, ipiv, work_arr, 12*n_sites, info)
    end do
call cpu_time(t2)
!write(*,*) "matrix multiplications with LAPACK:", t2-t1, "seconds"

call cpu_time(t1)
    do k = 1, 11
      do i = 1, n_sites
        A_i = 0.d0
        do j = 1, n_sites
          A_i = A_i + A_mat(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3,k)
        end do
        alpha_SCS(i,k) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
      end do
    end do
call cpu_time(t2)

!write(*,*) "alpha_SCS old method"
!do i = 1, n_sites
!  write(*,*) alpha_SCS(i,1)
!end do


!write(*,*) "A_i:", t2-t1, "seconds"

!write(*,*) "alpha_TS(1) =", neighbor_alpha0(1)
!write(*,*) "alpha_scs:", alpha_SCS(:,1)

call cpu_time(t1)
    A_LR = 0.d0

    do k = 1, 11
      do i = 1, n_sites
        do c1 = 1, 3
          A_LR(3*(i-1)+c1,3*(i-1)+c1,k) = alpha_SCS(i,k)
        end do
      end do
    end do
call cpu_time(t2)
!write(*,*) "A_LR:", t2-t1, "seconds"


call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "damping function:", t2-t1, "seconds"

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "T_LR", t2-t1, "seconds"

    write(*,*) "T_LR:"
    do p = 1, 6
      write(*,*) T_LR(p,1:6)
    end do



    I_mat = 0.d0
    do i = 1, 3*n_sites
      I_mat(i,i) = 1.d0
    end do

call cpu_time(t1)
    AT = 0.d0
    do k = 1, 11
      call dgemm('n', 'n', 3*n_sites, 3*n_sites, 3*n_sites, 1.d0, A_LR(:,:,k), 3*n_sites, T_LR, 3*n_sites, &
                0.d0, AT(:,:,k), 3*n_sites)
    end do
call cpu_time(t2)
!write(*,*) "A_LR x T_LR:", t2-t1, "seconds"

    write(*,*) "AT"
    do p = 1, 6
      write(*,*) AT(p,1:6,1)
    end do


call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "ugly loop:", t2-t1, "seconds"

    f_damp_der = 0.d0
    g_func_der = 0.d0
    h_func_der = 0.d0
    dT_SR = 0.d0

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "very ugly loop:", t2-t1, "seconds"

!    write(*,*) "dT_SR:"
!    write(*,*) dT_SR(3*(33-1)+1,1,1,1,1)
!    write(*,*) dT_SR(1,3*(33-1)+1,1,1,1)

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "very ugly gradients loop:", t2-t1, "seconds"

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
      do i = 1, n_sites
call cpu_time(t1)
        allocate( T_SR_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( B_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( A_mat_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( alpha_SCS_i(1:n_neigh(i),1:11) )
        allocate( A_LR_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( T_LR_i(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( AT_i(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( I_mat_n(1:3*n_neigh(i),1:3*n_neigh(i)) )
        allocate( energy_series(1:3*n_neigh(i),1:3,1:11) )
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
        allocate( AT_n(1:3*n_neigh(i),1:3*n_neigh(i), 1:n_order-1, 1:11) )
        allocate( force_series(1:3*n_neigh(i),1:3*n_neigh(i),1:11) )
        allocate( energy_term(1:3*n_neigh(i),1:3) )
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
call cpu_time(t2)
!write(*,*) "loop within loop for atom", i, "done in:", t2-t1, "seconds"

!        write(*,*) "B_mat_i:"
!        do p = 1, 6
!          write(*,*) B_mat_i(p,1:6,1)
!        end do


call cpu_time(t1)
        do k3 = 1, 11
          call dgetrf(3*n_neigh(i), 3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, info)
          call dgetri(3*n_neigh(i), A_mat_i(:,:,k3), 3*n_neigh(i), ipiv, work_arr, 12*n_sites, info)
        end do
call cpu_time(t2)
!write(*,*) "LAPACK for atom", i, "done in:", t2-t1, "seconds"

call cpu_time(t1)
        do k3 = 1, 11
          do p = 1, n_neigh(i)
            A_i = 0.d0
            do q = 1, n_neigh(i)
              A_i = A_i + A_mat_i(3*(p-1)+1:3*(p-1)+3,3*(q-1)+1:3*(q-1)+3,k3)
            end do
            alpha_SCS_i(p,k3) = 1.d0/3.d0 * (A_i(1,1)+A_i(2,2)+A_i(3,3))
          end do
        end do

!write(*,*) "alpha_scs_i", i, alpha_SCS_i(1,1)

!        write(*,*) "alpha_SCS_i:", alpha_SCS_i(:,1)

        do k3 = 1, 11
          do p = 1, n_neigh(i)
            do c1 = 1, 3
              A_LR_i(3*(p-1)+c1,3*(p-1)+c1,k3) = alpha_SCS_i(p,k3)
            end do
          end do
        end do
call cpu_time(t2)
!write(*,*) "Polarizability matrices for atom", i, "done in:", t2-t1, "seconds"

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "dipole-dipole tensor for atom", i, "done in:", t2-t1, "seconds"

!        write(*,*) "T_LR_i:"
!        do p = 1, 6
!          write(*,*) T_LR_i(p,1:6)
!        end do

call cpu_time(t1)
        AT_i = 0.d0
        do k = 1, 11
          call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, A_LR_i(:,:,k), 3*n_neigh(i), &
                     T_LR_i, 3*n_neigh(i), 0.d0, AT_i(:,:,k), 3*n_neigh(i))
        end do
call cpu_time(t2)
!write(*,*) "A_LR x T_LR for atom", i, "done in:", t2-t1, "seconds"

!        write(*,*) "AT_i:"
!        do p = 1, 6
!          write(*,*) AT_i(p,1:6,1)
!        end do

        I_mat_n = 0.d0
        do i2 = 1, 3*n_neigh(i)
          I_mat_n(i2,i2) = 1.d0
        end do

        if (series_expansion) then
call cpu_time(t1)
!open(unit=25, file="at.dat", status="unknown")
!do k = 1, 3*n_neigh(i)
!do k2 = 1, 3*n_neigh(i)
!write(25,*) AT_i(k2,k,1)
!end do
!end do
!close(25)
          AT_n = 0.d0
          integrand = 0.d0
          energy_series = -0.5d0 * AT_i(:,1:3,:)
          do k = 1, 11
            AT_n(:,:,1,k) = AT_i(:,:,k)
            do k2 = 1, n_order-2
              ! Precalculate the full AT_n for forces:
              call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, AT_i(:,:,k), 3*n_neigh(i), &
                         AT_n(:,:,k2,k), 3*n_neigh(i), 0.d0, AT_n(:,:,k2+1,k), 3*n_neigh(i))
              ! Use only slice for local energies:
              energy_term = 0.d0
              call dgemm('n', 'n', 3*n_neigh(i), 3, 3*n_neigh(i), 1.d0, AT_n(:,:,k2,k), 3*n_neigh(i), &
                         AT_i(:,1:3,k), 3*n_neigh(i), 0.d0, energy_term, 3*n_neigh(i))
              if (series_average .and. k2 == n_order-2) then
                energy_series(:,:,k) = (2.d0 * energy_series(:,:,k) - 1.d0/(k2+2)*energy_term)/2.d0 
              else
                energy_series(:,:,k) = energy_series(:,:,k) - 1.d0/(k2+2)*energy_term
              end if
            end do
            do c1 = 1, 3
              integrand(k) = integrand(k) + alpha_SCS_i(1,k) * dot_product(T_LR_i(c1,:), &
                             energy_series(:,c1,k))
            end do
!write(*,*) k, integrand(k)
          end do
          write(*,*) "integrand:", i, integrand
          integral = 0.d0
          call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
          E_MBD_k(i) = integral/(2.d0*pi)
          E_MBD_k(i) = E_MBD_k(i) * 27.211386245988
call cpu_time(t2)
!write(*,*) "MBD energy for atom", i, "done in:", t2-t1, "seconds"
          write(*,*) "E_MBD_k:", i, E_MBD_k(i)
          ! Calculate total MBD energy inside the cutoff sphere (mostly for checking finite difference)
          if (total_energy) then
            force_series = 0.d0
            integrand = 0.d0
            do k = 1, 11
              do k2 = 2, n_order-1
                force_series(:,:,k) = force_series(:,:,k) - 1.d0/k2 * AT_n(:,:,k,k2)
              end do
              VL_n = 0.d0
              call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, AT_i(:,:,k), 3*n_neigh(i), &
                         AT_n(:,:,k,n_order-1), 3*n_neigh(i), 0.d0, VL_n, 3*n_neigh(i))
              if ( series_average ) then
                force_series(:,:,k) = (2.d0*force_series(:,:,k) -1.d0/n_order * VL_n)/2.d0
              else
                force_series(:,:,k) = force_series(:,:,k) -1.d0/n_order * VL_n
              end if
              do i2 = 1, 3*n_neigh(i)
                integrand(k) = integrand(k) + force_series(i2,i2,k)
              end do
            end do
            integral = 0.d0
            call integrate("trapezoidal", omegas, integrand, 0.d0, 10.d0, integral)
            E_MBD = integral/(2.d0*pi)
            E_MBD = E_MBD * 27.211386245988
            write(*,*) "E_MBD_k_tot:", i, E_MBD
          end if
        else
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
!         Local energy test:
          integrand = 0.d0
          integrand_k = 0.d0
!          write(*,*) "alpha_SCS_i:", i, alpha_SCS_i(1,1)
          do k = 1,11
            do i2 = 1,3*n_neigh(i)
              integrand(k) = integrand(k) + logIAT_n(i2,i2,k)
            end do
            integrand_k(k) = 1.d0/n_neigh(i) * integrand(k)
          end do
!          write(*,*) "Integrand:", integrand
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
!          write(*,*) "Total MBD energy:", sum(E_MBD_k(:,1))
        end if


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
call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "dB_mat stuff:", t2-t1, "seconds"

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

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "LAPACK with dA_mat stuff:", t2-t1, "seconds"

!        write(*,*) "dA_mat_n:"
!        write(*,*) dA_mat_n(3*(1-1)+1,1:6,1,1)
!        write(*,*) dA_mat_n(3*(1-1)+2,1:6,1,1)
!        write(*,*) dA_mat_n(3*(1-1)+3,1:6,1,1)

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "dalpha:", t2-t1, "seconds"

        write(*,*) "dalpha_n:", dalpha_n(1,1,1)
!        write(*,*) "dalpha_n(42):", dalpha_n(42,1,1)

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "dA_LR:", t2-t1, "seconds"

!        do i2 = 1, 6
!          write(*,*) "dA_LR_n:", dA_LR_n(i2,1:6,1,1)
!        end do

call cpu_time(t1)
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
call cpu_time(t2)
!write(*,*) "dT_LR:", t2-t1, "seconds"

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
        force_series = 0.d0

        if (series_expansion) then
call cpu_time(t1)
          do k = 1, 11
            do k2 = 1, n_order-1
              if (series_average .and. k2 == n_order-1) then
                force_series(:,:,k) = (2.d0 * force_series(:,:,k) + AT_n(:,:,k,k2))/2.d0
              else
                force_series(:,:,k) = force_series(:,:,k) + AT_n(:,:,k,k2)
              end if
            end do
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
              call dgemm('n', 'n', 3*n_neigh(i), 3*n_neigh(i), 3*n_neigh(i), 1.d0, force_series(:,:,k), 3*n_neigh(i), &
                         G_mat_n(:,:,c3,k), 3*n_neigh(i), 0.d0, VL_n, 3*n_neigh(i))
              do i2 = 1, 3*n_neigh(i)
                force_integrand_n(c3,k) = force_integrand_n(c3,k) + VL_n(i2,i2)
              end do
            end do
          end do
call cpu_time(t2)
!write(*,*) "LAPACK for force:", t2-t1, "seconds"
        else
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
        end if

call cpu_time(t1)
        integral = 0.d0
        do c3 = 1, 3
          call integrate("trapezoidal", omegas, force_integrand_n(c3,:), 0.d0, 10.d0, integral)
          forces_MBD_k(i,c3) = 1.d0/(2.d0*pi) * integral
        end do
        forces_MBD_k(i,:) = forces_MBD_k(i,:) * 51.42208619083232
        write(*,*) "force_k:", i, forces_MBD_k(i,:)
call cpu_time(t2)
!write(*,*) "force integration:", t2-t1, "seconds"

        deallocate( T_SR_i, B_mat_i, A_mat_i, alpha_SCS_i, A_LR_i, T_LR_i, AT_i, I_mat_n, &
                    WR_n, WI_n, VR_n, VL_n, VRinv_n, logMapo_n, logIAT_n, dB_mat_n, dA_mat_n, dBA_n, dalpha_n, &
                    dA_LR_n, dT_LR_n, invIAT_n, G_mat_n, energy_series, AT_n, force_series, energy_term )
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


!    if( do_timing) then
!      call cpu_time(time2)
!      write(*,*) "vdw: computing forces:", time2-time1
!    end if
 
!   Clean up
    deallocate( neighbor_c6_ii, f_damp, T_func, &
                h_func, g_func, omegas, omega_i, alpha_i, sigma_i, sigma_ij, T_SR, B_mat, xyz_H, rjs_H, A_mat, A_i, &
                alpha_SCS, A_LR, r0_ii_SCS, f_damp_SCS, I_mat, AT, logIAT, logMapo, VR, VRinv, WR, WI, VL, integrand, &
                dT, dT_SR, f_damp_der, g_func_der, h_func_der, dA_mat, dalpha, dA_LR, dT_LR, f_damp_der_SCS, &
                invIAT, G_mat, force_integrand, forces_MBD, coeff_h_der, terms, dT_SR_A_mat, dT_SR_v, &
                dB_mat, dB_mat_v, &
                coeff_der, coeff_fdamp, dg, dh, hirshfeld_v_cart_der_H, alpha_k, sigma_k, s_i, s_j, &
                E_MBD_k, f_damp_der_ij_n, f_damp_der_SCS_ij_n, force_integrand_n, forces_MBD_k, integrand_k )

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
