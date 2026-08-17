! Validation of the central-atom SOAP Hessian.
!
! get_soap_central_hessian returns
!
!     hess(a,b,v,k2) = sum_d vecs(d,v,i) d^2 q_d / d r_ia d r_jb
!
! for pair k2 = (site i, neighbour j). soap_cart_der(b,:,k2) is d q / d r_jb,
! so this is the derivative of something the first-derivative path already
! returns, with respect to the central atom's position. The test is therefore a
! central difference in r_i of v . soap_cart_der(b,:,k2). It exercises the
! whole chain at once: the radial second derivatives, the pole-free angular
! Hessian, the S/Z/ZG contraction and the normalisation algebra.
!
! Two checks, and one non-check worth being explicit about:
!
!   1. FD agreement, scanned over h. A single h proves little -- a wrong
!      derivative can sit within truncation error at one step size. What
!      identifies a correct derivative is the shape: h^2 decay until round-off
!      takes over, and a clear minimum. The block set is held fixed across the
!      scan (the neighbour-list guard uses a fixed margin, not one proportional
!      to h) so the numbers are comparable to each other.
!
!   2. The reconstruction self-test inside the routine. The Hessian rebuilds
!      the angular coefficient from scratch via Q_l^m, so the same contraction
!      also yields the FIRST derivative, which is already known independently
!      from soap_rad_der/azi/pol. Agreement there tests Omega, dOmega, F, F',
!      the S and Z build, compression and the species mask with no finite
!      differences involved, which is what makes it useful for localising a
!      failure to the second-order terms.
!
!   NOT a check: symmetry of the self block, and the translational sum rule.
!   Both are automatic. The self block is built as minus the sum over
!   neighbours, and every term in that sum -- including the D.ZG cross term,
!   which is asymmetric pair by pair -- becomes symmetric once summed, because
!   sum_j D^j = G. They are printed below as bookkeeping checks only. Do not
!   read them as evidence that the second derivatives are right; an early
!   version of this code had the D.ZG term double-counted and passed both at
!   1e-16 while failing the finite differences by 3 per cent.
!
program hverify
   use gharness
   use soap_turbo_desc
   use soap_turbo_radial, only: legacy_filter_seed
   implicit none

   integer, parameter :: nvec = 2
   integer :: ic
   logical :: comp
! Shared by every internal procedure below through host association: Fortran
! does not allow a contained procedure to contain further procedures, so these
! live here rather than inside one_case.
   real(dp), allocatable :: soap(:, :), der(:, :, :), hess(:, :, :, :), vecs(:, :, :)
   real(dp), allocatable :: sp(:, :), dsp(:, :, :), sm(:, :), dsm(:, :, :)
   integer :: n_soap, kbase
   real(dp) :: scale

   call setup_system(40, 20260817, 8.0d0)
   call build_neighbours()
   legacy_filter_seed = .false.

   write (*, '(A,I0,A,I0)') "system: ", n_atoms, " atoms, ", n_atom_pairs, " pairs"
   write (*, *)

   do ic = 1, 2
      comp = (ic == 1)
      call one_case(comp)
   end do

contains

   subroutine one_case(comp)
      logical, intent(in) :: comp
      integer :: i, j, k, k2, k3, a, b, v, ih, ntested
      real(dp) :: h, worst, u, worstsym, worstsum, acc(3, 3)

      if (comp) then
         n_soap = n_soap_compressed
      else
         n_soap = n_soap_uncompressed
      end if
      write (*, '(A,L1,A,I0)') "=== compress_soap = ", comp, "   n_soap = ", n_soap

      allocate (soap(n_soap, n_sites), der(3, n_soap, n_atom_pairs))
      allocate (sp(n_soap, n_sites), dsp(3, n_soap, n_atom_pairs))
      allocate (sm(n_soap, n_sites), dsm(3, n_soap, n_atom_pairs))
      allocate (hess(3, 3, nvec, n_atom_pairs), vecs(n_soap, nvec, n_sites))

!   Fixed weight vectors. They stand in for dE_i/dq, and must not move when the
!   atoms do -- their own motion is a separate term in a dipole gradient.
      do i = 1, n_sites
         do v = 1, nvec
            do k = 1, n_soap
               u = dble(mod(7919*(k + 131*v + 17*i), 2003))/2003.d0 - 0.5d0
               vecs(k, v, i) = u
            end do
         end do
      end do

      soap_hessian_selftest = .true.
      call run_hess(soap, der, hess)
      soap_hessian_selftest = .false.

      scale = maxval(abs(hess))
      write (*, '(A,ES11.3)') "  max |hess|                                    ", scale
      write (*, '(A,ES11.3)') "  reconstruction self-test (no FD), rel err     ", soap_hessian_selftest_err

      worstsym = 0.d0
      worstsum = 0.d0
      k2 = 0
      do i = 1, n_sites
         k3 = k2 + 1
         acc = 0.d0
         do j = 1, n_neigh(i)
            k2 = k2 + 1
            acc(:, :) = acc(:, :) + hess(:, :, 1, k2)
         end do
         worstsum = max(worstsum, maxval(abs(acc)))
         do v = 1, nvec
            do a = 1, 3
               do b = a + 1, 3
                  worstsym = max(worstsym, abs(hess(a, b, v, k3) - hess(b, a, v, k3)))
               end do
            end do
         end do
      end do
      write (*, '(A,ES11.3)') "  self-block symmetry (automatic), abs          ", worstsym
      write (*, '(A,ES11.3)') "  translational sum rule (automatic), abs       ", worstsum
      write (*, '(A)') "  finite differences in the central atom position:"
      do ih = 1, 7
         h = 1.d-2*10.d0**(-0.5d0*(ih - 1))
         call fdcheck(h, worst, ntested)
         write (*, '(A,ES9.2,A,I0,A,ES11.3)') "     h = ", h, "  over ", ntested, &
            " blocks:  max rel err ", worst
      end do
      write (*, *)

      deallocate (soap, der, sp, dsp, sm, dsm, hess, vecs)
   end subroutine one_case

   subroutine run_hess(s, d, hh)
      real(dp), intent(out) :: s(:, :), d(:, :, :), hh(:, :, :, :)
      soap_hessian_enabled = .true.
      s = 0.d0; d = 0.d0
      call get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                    thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                    atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, &
                    radial_enhancement, central_weight, basis, scaling_mode, .false., .true., comp, &
                    compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, s, d)
      call get_soap_central_hessian(n_sites, n_neigh, n_atom_pairs, n_species, mask, rjs, thetas, phis, &
                                    l_max, rcut_hard, atom_sigma_t, atom_sigma_t_scaling, &
                                    comp, compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, &
                                    s, nvec, vecs, hh)
      soap_hessian_enabled = .false.
   end subroutine run_hess

   subroutine run_plain(s, d)
      real(dp), intent(out) :: s(:, :), d(:, :, :)
      s = 0.d0; d = 0.d0
      call get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                    thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                    atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, &
                    radial_enhancement, central_weight, basis, scaling_mode, .false., .true., comp, &
                    compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, s, d)
   end subroutine run_plain

!   Moving atom ai must not change anybody's neighbour list, or the pair
!   indices shift under the difference. The margin is fixed, not proportional
!   to h, so that every step size tests the same set of blocks.
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

   subroutine fdcheck(hh, worst, ndone)
      real(dp), intent(in) :: hh
      real(dp), intent(out) :: worst
      integer, intent(out) :: ndone
      integer :: ii, jj, aa, bb, vv, kk, kself
      real(dp) :: gp, gm, fd, an, save_pos(3)

      worst = 0.d0
      ndone = 0
      kbase = 0
      do ii = 1, n_sites
         kself = kbase + 1
!       one site in seven, so the scan over h stays cheap
         if (mod(ii, 7) /= 1 .or. .not. stable(ii)) then
            kbase = kbase + n_neigh(ii)
            cycle
         end if

         do aa = 1, 3
            save_pos = pos(:, ii)
            pos(aa, ii) = save_pos(aa) + hh
            call build_neighbours()
            call run_plain(sp, dsp)
            pos(:, ii) = save_pos
            pos(aa, ii) = save_pos(aa) - hh
            call build_neighbours()
            call run_plain(sm, dsm)
            pos(:, ii) = save_pos
            call build_neighbours()

!         the self pair, which is the block a dipole gradient needs, plus a
!         sample of the neighbour blocks
            do jj = 1, n_neigh(ii)
               kk = kself + jj - 1
               if (jj > 1 .and. mod(jj, 3) /= 2) cycle
               do vv = 1, nvec
                  do bb = 1, 3
                     gp = dot_product(vecs(:, vv, ii), dsp(bb, :, kk))
                     gm = dot_product(vecs(:, vv, ii), dsm(bb, :, kk))
                     fd = (gp - gm)/(2.d0*hh)
                     an = hess(aa, bb, vv, kk)
                     worst = max(worst, abs(fd - an)/scale)
                  end do
               end do
               ndone = ndone + 1
            end do
         end do
         kbase = kbase + n_neigh(ii)
      end do
   end subroutine fdcheck

end program hverify
