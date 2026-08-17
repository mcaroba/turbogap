! First derivatives: soap_cart_der against central differences of soap.
!
! soap_cart_der(:,:,k2) for pair k2 = (site i, neighbour j) is d q_i / d r_j,
! so displacing atom j and differencing q_i tests it directly. The self pair
! j = 1 is included: soap_turbo builds it as minus the sum over the others, and
! it is the one a dipole model reads, so it is worth testing explicitly rather
! than trusting the identity.
!
! Both compression settings are covered. The uncompressed path cannot be
! reached end to end with any fitted potential in turbogap_tests -- every GAP
! there was fitted with compression on -- so driving get_soap directly is the
! only way to exercise it.
!
program dverify
   use gharness
   use soap_turbo_desc
   use soap_turbo_radial, only: legacy_filter_seed
   implicit none

   integer :: ic, ih, ntested
   logical :: comp, legacy
   real(dp) :: h, worst
   real(dp), allocatable :: soap(:, :), der(:, :, :), sp(:, :), sm(:, :), dd(:, :, :)
   integer :: n_soap
   real(dp) :: scale

   call setup_system(40, 20260817, 8.0d0)
   call build_neighbours()

   write (*, '(A,I0,A,I0)') "system: ", n_atoms, " atoms, ", n_atom_pairs, " pairs"
   write (*, *)

! Both filter seeds. The legacy seed drops a surface term at rcut_hard, so its
! derivatives are inconsistent with its own values at the 1e-6 level; that is a
! known, documented property of the pre-2026 convention and this test is what
! makes it visible rather than a surprise.
   do ic = 1, 4
      comp = (mod(ic, 2) == 1)
      legacy = (ic <= 2)
      call one_case(comp, legacy)
   end do

contains

   subroutine one_case(comp, legacy)
      logical, intent(in) :: comp, legacy

      legacy_filter_seed = legacy
      if (comp) then
         n_soap = n_soap_compressed
      else
         n_soap = n_soap_uncompressed
      end if
      write (*, '(A,L1,A,L1,A,I0)') "=== compress_soap = ", comp, &
         "   legacy_filter_seed = ", legacy, "   n_soap = ", n_soap

      allocate (soap(n_soap, n_sites), der(3, n_soap, n_atom_pairs))
      allocate (sp(n_soap, n_sites), sm(n_soap, n_sites), dd(3, n_soap, n_atom_pairs))
      call run(comp, soap, der)
      scale = maxval(abs(der))
      write (*, '(A,ES11.3)') "  max |soap_cart_der| ", scale
      do ih = 1, 6
         h = 1.d-2*10.d0**(-0.5d0*(ih - 1))
         call fdcheck(comp, h, worst, ntested)
         write (*, '(A,ES9.2,A,I0,A,ES11.3)') "     h = ", h, "  over ", ntested, &
            " pairs:  max rel err ", worst
      end do
      write (*, *)
      deallocate (soap, der, sp, sm, dd)
   end subroutine one_case

   subroutine run(comp, s, d)
      logical, intent(in) :: comp
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
      real(dp) :: v(3), rr
      stable = .true.
      do jj = 1, n_atoms
         if (jj == ai) cycle
         v = pos(:, jj) - pos(:, ai)
         rr = dsqrt(sum(v*v))
         if (abs(rr - rcut_max) < 0.05d0) stable = .false.
         if (rr < 0.5d0) stable = .false.
      end do
   end function stable

   subroutine fdcheck(comp, hh, worst, ndone)
      logical, intent(in) :: comp
      real(dp), intent(in) :: hh
      real(dp), intent(out) :: worst
      integer, intent(out) :: ndone
      integer :: ii, jj, bb, kk, kself, kb, aj, nn
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
               call run(comp, sp, dd)
               pos(:, aj) = save_pos
               pos(bb, aj) = save_pos(bb) - hh
               call build_neighbours()
               call run(comp, sm, dd)
               pos(:, aj) = save_pos
               call build_neighbours()
               do nn = 1, n_soap
                  fd = (sp(nn, ii) - sm(nn, ii))/(2.d0*hh)
                  worst = max(worst, abs(fd - der(bb, nn, kk))/scale)
               end do
            end do
            ndone = ndone + 1
         end do
         kb = kb + n_neigh(ii)
      end do
   end subroutine fdcheck

end program dverify
