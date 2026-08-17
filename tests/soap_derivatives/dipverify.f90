! End-to-end validation of the dipole gradient.
!
! A GAP dipole model predicts mu_i = dE_i/dr_i for a fictitious local energy
!
!     E_i = delta^2 sum_s alpha_s (q_i . q_s)^zeta,
!
! so the dipole is w_i . dq/dr_i with w_i = dE_i/dq. Differentiating that with
! respect to a neighbour's position gives two terms, and both have to be right:
!
!   A   the weight vector moves:  sum_de (d2E/dq_d dq_e)(dq_d/dr_ia)(dq_e/dr_jb)
!       = zeta(zeta-1) delta^2 sum_s alpha_s k_s^(zeta-2)
!         (q_s . dq/dr_ia)(q_s . dq/dr_jb)
!
!   B   the descriptor curves:    sum_d w_d d2q_d/dr_ia dr_jb
!
! Term B is what get_soap_central_hessian supplies, and hverify checks it
! against finite differences of soap_cart_der at fixed w. That test says nothing
! about term A, because it holds w still. This one differences the dipole
! itself, so it is the first test in which both terms have to be simultaneously
! correct, with the right relative weight and sign.
!
! Term A is written so it stays linear in the number of pairs: the naive form
! needs (q_s . dq/dr_jb) for every sparse point and every pair, which is
! n_sparse times the cost of a force contraction. Collapsing the sparse sum
! into three descriptor-space vectors per site first,
!
!     V_a = sum_s beta_s (q_s . dq/dr_ia) q_s,      A_ab = V_a . dq/dr_jb,
!
! moves the n_sparse work to once per site and leaves 9*n_soap per pair.
!
program dipverify
   use gharness
   use soap_turbo_desc
   use soap_turbo_radial, only: legacy_filter_seed
   implicit none

   integer, parameter :: n_sparse = 12
   integer :: n_soap, i, j, k, s, a, b, k2, kself, ih, ntested
   real(dp) :: delta, zeta, u, h, worst, scale
   real(dp), allocatable :: soap(:, :), der(:, :, :), hess(:, :, :, :)
   real(dp), allocatable :: Qsp(:, :), alph(:), kern(:, :), wts(:, :, :)
   real(dp), allocatable :: Vv(:, :, :), dmu(:, :, :), mu(:, :)
   real(dp), allocatable :: sp(:, :), dsp(:, :, :), sm(:, :), dsm(:, :, :), mup(:, :), mum(:, :)
   logical, parameter :: comp = .true.

   call setup_system(40, 20260817, 8.0d0)
   call build_neighbours()
   legacy_filter_seed = .false.
   n_soap = n_soap_compressed
   delta = 1.3d0
   zeta = 4.d0

   allocate (soap(n_soap, n_sites), der(3, n_soap, n_atom_pairs))
   allocate (sp(n_soap, n_sites), dsp(3, n_soap, n_atom_pairs))
   allocate (sm(n_soap, n_sites), dsm(3, n_soap, n_atom_pairs))
   allocate (hess(3, 3, 1, n_atom_pairs), wts(n_soap, 1, n_sites))
   allocate (Qsp(n_soap, n_sparse), alph(n_sparse), kern(n_sites, n_sparse))
   allocate (Vv(n_soap, 3, n_sites), dmu(3, 3, n_atom_pairs))
   allocate (mu(3, n_sites), mup(3, n_sites), mum(3, n_sites))

! A stand-in sparse set. Normalised, like real sparse points.
   do s = 1, n_sparse
      do k = 1, n_soap
         u = dble(mod(104729*(k + 37*s), 7919))/7919.d0 - 0.5d0
         Qsp(k, s) = u
      end do
      Qsp(:, s) = Qsp(:, s)/dsqrt(dot_product(Qsp(:, s), Qsp(:, s)))
      alph(s) = 0.3d0 + dble(mod(7*s, 5))*0.11d0
   end do

   write (*, '(A,I0,A,I0,A,I0,A,I0,A,F4.1)') "system: ", n_atoms, " atoms, ", n_atom_pairs, &
      " pairs, n_soap = ", n_soap, ", n_sparse = ", n_sparse, ",  zeta = ", zeta
   write (*, *)

   call model(.true.)

   scale = maxval(abs(dmu))
   write (*, '(A,ES11.3)') "max |d mu / d r| = ", scale

! Automatic, but it catches pair-index bookkeeping: a rigid shift cannot change
! a dipole, so the blocks of one site sum to zero.
   worst = 0.d0
   k2 = 0
   do i = 1, n_sites
      kself = k2 + 1
      do j = 1, n_neigh(i)
         k2 = k2 + 1
      end do
      do b = 1, 3
         do a = 1, 3
            worst = max(worst, abs(sum(dmu(a, b, kself:kself + n_neigh(i) - 1))))
         end do
      end do
   end do
   write (*, '(A,ES11.3)') "translational sum rule (automatic), abs = ", worst
   write (*, *)

   write (*, '(A)') "finite differences of the dipole in a neighbour's position:"
   do ih = 1, 7
      h = 1.d-2*10.d0**(-0.5d0*(ih - 1))
      call fdcheck(h, worst, ntested)
      write (*, '(A,ES9.2,A,I0,A,ES11.3)') "   h = ", h, "  over ", ntested, &
         " blocks:  max rel err ", worst
   end do

contains

! Descriptor, dipole and, if asked, the analytic dipole gradient.
   subroutine model(want_grad)
      logical, intent(in) :: want_grad
      real(dp) :: bet, gsa(3), kk
      integer :: ii, jj, aa, bb, ss, kb

      soap_hessian_enabled = want_grad
      soap = 0.d0; der = 0.d0
      call run(soap, der)

!   kernels and the weight vector w_i = dE_i/dq
      kern = matmul(transpose(soap), Qsp)
      do ii = 1, n_sites
         wts(:, 1, ii) = 0.d0
         do ss = 1, n_sparse
            wts(:, 1, ii) = wts(:, 1, ii) + (zeta*delta**2*alph(ss)*kern(ii, ss)**(nint(zeta) - 1))*Qsp(:, ss)
         end do
      end do

!   mu_i = w_i . dq/dr_i, the self pair
      kb = 0
      do ii = 1, n_sites
         kb = kb + 1
         do aa = 1, 3
            mu(aa, ii) = dot_product(wts(:, 1, ii), der(aa, :, kb))
         end do
         kb = kb + n_neigh(ii) - 1
      end do
      if (.not. want_grad) return

      call get_soap_central_hessian(n_sites, n_neigh, n_atom_pairs, n_species, mask, rjs, thetas, phis, &
                                    l_max, rcut_hard, atom_sigma_t, atom_sigma_t_scaling, &
                                    comp, compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, &
                                    soap, 1, wts, hess)
      soap_hessian_enabled = .false.

!   V_a, three descriptor-space vectors per site, so that term A costs
!   9*n_soap per pair instead of n_sparse*n_soap
      kb = 0
      do ii = 1, n_sites
         kb = kb + 1
         Vv(:, :, ii) = 0.d0
         do ss = 1, n_sparse
            kk = kern(ii, ss)
            bet = zeta*(zeta - 1.d0)*delta**2*alph(ss)*kk**(nint(zeta) - 2)
            do aa = 1, 3
               gsa(aa) = dot_product(Qsp(:, ss), der(aa, :, kb))
            end do
            do aa = 1, 3
               Vv(:, aa, ii) = Vv(:, aa, ii) + (bet*gsa(aa))*Qsp(:, ss)
            end do
         end do
         kb = kb + n_neigh(ii) - 1
      end do

      kb = 0
      do ii = 1, n_sites
         do jj = 1, n_neigh(ii)
            kb = kb + 1
            do bb = 1, 3
               do aa = 1, 3
                  dmu(aa, bb, kb) = dot_product(Vv(:, aa, ii), der(bb, :, kb)) + hess(aa, bb, 1, kb)
               end do
            end do
         end do
      end do
   end subroutine model

   subroutine run(s, d)
      real(dp), intent(out) :: s(:, :), d(:, :, :)
      s = 0.d0; d = 0.d0
      call get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                    thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                    atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, &
                    radial_enhancement, central_weight, basis, scaling_mode, .false., .true., comp, &
                    compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, s, d)
   end subroutine run

   logical function stable(ai)
      integer, intent(in) :: ai
      integer :: jj
      real(dp) :: dd(3), rr
      stable = .true.
      do jj = 1, n_atoms
         if (jj == ai) cycle
         dd = pos(:, jj) - pos(:, ai)
         rr = dsqrt(sum(dd*dd))
         if (abs(rr - rcut_max) < 0.05d0) stable = .false.
         if (rr < 0.5d0) stable = .false.
      end do
   end function stable

! Displace the NEIGHBOUR of a pair and difference mu of the site.
   subroutine fdcheck(hh, worst, ndone)
      real(dp), intent(in) :: hh
      real(dp), intent(out) :: worst
      integer, intent(out) :: ndone
      integer :: ii, jj, aa, bb, kk, kb, aj
      real(dp) :: save_pos(3), fd

      worst = 0.d0
      ndone = 0
      kb = 0
      do ii = 1, n_sites
         kself = kb + 1
         if (mod(ii, 9) /= 1) then
            kb = kb + n_neigh(ii)
            cycle
         end if
         do jj = 1, n_neigh(ii)
            kk = kself + jj - 1
            if (jj > 1 .and. mod(jj, 5) /= 2) cycle
            aj = nl_atom(kk)
            if (.not. stable(aj)) cycle
            do bb = 1, 3
               save_pos = pos(:, aj)
               pos(bb, aj) = save_pos(bb) + hh
               call build_neighbours()
               call model(.false.)
               mup = mu
               pos(:, aj) = save_pos
               pos(bb, aj) = save_pos(bb) - hh
               call build_neighbours()
               call model(.false.)
               mum = mu
               pos(:, aj) = save_pos
               call build_neighbours()
               do aa = 1, 3
                  fd = (mup(aa, ii) - mum(aa, ii))/(2.d0*hh)
                  worst = max(worst, abs(fd - dmu(aa, bb, kk))/scale)
               end do
            end do
            ndone = ndone + 1
         end do
         kb = kb + n_neigh(ii)
      end do
   end subroutine fdcheck

end program dipverify
