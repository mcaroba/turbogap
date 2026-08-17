! The complete MAD IR force, end to end.
!
!     f_jb = - dL/dr_jb = - sum_a lambda_a  dmu_a/dr_jb
!
! Every link in the chain is exercised at once, and against the real production
! routines: this program links libturbogap.a rather than recompiling anything,
! so get_soap_dipole_weights and accumulate_dmu_dr are the versions that ship.
!
!   get_soap                       -> q and dq/dr
!   get_soap_dipole_weights        -> w = dE/dq, and V (the kernel-curvature term)
!   get_soap_central_hessian       -> d2q/dr_i dr_j contracted against w
!   accumulate_dmu_dr              -> dmu_a/dr_jb, scattered per atom
!   mad_ir_evaluate                -> lambda = dE/dmu(newest)
!   mad_ir_forces                  -> f
!
! The reference is a central difference of the loss itself with respect to an
! atom's position, with the older frames of the ensemble held fixed -- which is
! exactly the approximation the scheme makes, so the test measures the thing
! the code claims to compute rather than an idealisation of it.
!
! A short window is used on purpose. Resolution is irrelevant to a gradient
! test and a long one would only make it slow; what matters is that the newest
! frame enters C(0) and every C(tau), which it does for any n_lag > 1.
!
program madverify
!  gharness declares dp itself (soap_turbo is a submodule and must not depend on
!  the parent's kinds module), so importing kinds as well would make dp
!  ambiguous. They are the same real64 either way.
   use gharness
   use soap_turbo_desc
   use soap_turbo_radial, only: legacy_filter_seed
   use types, only: soap_turbo
   use gap, only: soap_backend_begin, get_soap_dipole_weights, accumulate_dmu_dr
   use mad_ir, only: mad_ir_type, mad_ir_init, mad_ir_push, mad_ir_evaluate, mad_ir_forces
   implicit none

   integer, parameter :: n_sparse = 10
   integer, parameter :: n_lag = 60, n_window = 120, n_freq = 80
   type(soap_turbo), target :: hyp
   type(mad_ir_type) :: st
   integer :: n_soap, i, k, s, a, b, ih, ntested, aj
   real(dp) :: delta, zeta, u, h, worst, scale, loss, lambda(3), dt, escale
   real(dp), allocatable :: soap(:, :), der(:, :, :), hess(:, :, :, :)
   real(dp), allocatable :: w(:, :), wv(:, :, :), V(:, :, :), dmu_dr(:, :, :)
   real(dp), allocatable :: nu(:), Iexp(:), wgt(:), forces(:, :), mu_hist(:, :)
   real(dp) :: mu_tot(3)
   logical, parameter :: comp = .true.

   call setup_system(40, 20260817, 9.5d0)
   call build_neighbours()
   legacy_filter_seed = .false.
   n_soap = n_soap_compressed
   delta = 1.3d0
   zeta = 4.d0
   dt = 2.d0
   escale = 3.7d0        ! exp_energy_scales for this observable

   allocate (soap(n_soap, n_sites), der(3, n_soap, n_atom_pairs))
   allocate (hess(3, 3, 1, n_atom_pairs), wv(n_soap, 1, n_sites))
   allocate (w(n_soap, n_sites), V(n_soap, 3, n_sites))
   allocate (dmu_dr(3, 3, n_atoms), forces(3, n_atoms))
   allocate (nu(n_freq), Iexp(n_freq), wgt(n_freq), mu_hist(3, n_window))

!  A stand-in sparse set, bound the way the driver binds a real one.
   allocate (hyp%alphas(n_sparse), hyp%Qs(n_soap, n_sparse))
   do s = 1, n_sparse
      do k = 1, n_soap
         u = dble(mod(104729*(k + 37*s), 7919))/7919.d0 - 0.5d0
         hyp%Qs(k, s) = u
      end do
      hyp%Qs(:, s) = hyp%Qs(:, s)/dsqrt(dot_product(hyp%Qs(:, s), hyp%Qs(:, s)))
      hyp%alphas(s) = 0.3d0 + dble(mod(7*s, 5))*0.11d0
   end do
   call soap_backend_begin(hyp)

   do k = 1, n_freq
      nu(k) = 4000.d0*dble(k)/dble(n_freq)
      wgt(k) = 1.d0
   end do

!  Fill the ensemble: n_window-1 frames of arbitrary history, then the frame
!  belonging to the actual configuration.
   do k = 1, n_window - 1
      mu_hist(1, k) = 0.9d0*dcos(0.31d0*dble(k)) + 0.2d0*dcos(0.07d0*dble(k))
      mu_hist(2, k) = 0.5d0*dcos(0.53d0*dble(k) + 0.4d0)
      mu_hist(3, k) = 0.7d0*dcos(0.17d0*dble(k) + 1.1d0)
   end do

   write (*, '(A,I0,A,I0,A,I0)') "system: ", n_atoms, " atoms, ", n_atom_pairs, &
      " pairs, n_soap = ", n_soap
   write (*, '(A,I0,A,I0,A,F6.2,A)') "ensemble: n_window = ", n_window, ", n_lag = ", n_lag, &
      ", dt = ", dt, " fs"

!  A target built from a spectrum this configuration actually produces, so the
!  fitted scale is O(1) rather than absorbing ten orders of magnitude.
   call init_state()
   call dipole(mu_tot)
   st%mu_hist(:, st%head) = mu_tot
   call mad_ir_evaluate(st, escale, loss, lambda)
   Iexp = 1.2d0*st%I_calc + 0.1d0*maxval(abs(st%I_calc))
   call init_state()

   call analytic(forces, loss, lambda)
   scale = maxval(abs(forces))
   write (*, '(A,ES11.3,A,ES11.3)') "loss = ", loss, "   max |f_MAD| = ", scale
   write (*, '(A,3ES12.4)') "lambda = ", lambda
   write (*, *)

   write (*, '(A)') "finite differences of the loss in an atom position:"
   do ih = 1, 7
      h = 1.d-2*10.d0**(-0.5d0*(ih - 1))
      call fdcheck(h, worst, ntested)
      write (*, '(A,ES9.2,A,I0,A,ES11.3)') "   h = ", h, "  over ", ntested, &
         " atoms:  max rel err ", worst
   end do

contains

   subroutine init_state()
      integer :: kk
      call mad_ir_init(st, dt, n_lag, n_window, nu, Iexp, wgt, .true., 2.d0, "hann")
      do kk = 1, n_window - 1
         call mad_ir_push(st, mu_hist(:, kk))
      end do
      call mad_ir_push(st, [0.d0, 0.d0, 0.d0])
   end subroutine init_state

   subroutine run_soap(want_hess)
      logical, intent(in) :: want_hess
      soap_hessian_enabled = want_hess
      soap = 0.d0; der = 0.d0
      call get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                    thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                    atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, &
                    radial_enhancement, central_weight, basis, scaling_mode, .false., .true., comp, &
                    compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, soap, der)
   end subroutine run_soap

!  Total dipole of the current configuration.
   subroutine dipole(mt)
      real(dp), intent(out) :: mt(3)
      integer :: ii, kb, aa
      call run_soap(.false.)
      call get_soap_dipole_weights(soap, der, delta, zeta, n_neigh, w, V)
      mt = 0.d0
      kb = 0
      do ii = 1, n_sites
         kb = kb + 1
         do aa = 1, 3
            mt(aa) = mt(aa) + dot_product(w(:, ii), der(aa, :, kb))
         end do
         kb = kb + n_neigh(ii) - 1
      end do
      soap_hessian_enabled = .false.
   end subroutine dipole

   subroutine analytic(f, l_out, lam)
      real(dp), intent(out) :: f(:, :), l_out, lam(3)
      integer :: ii
      call run_soap(.true.)
      call get_soap_dipole_weights(soap, der, delta, zeta, n_neigh, w, V)
      mu_tot = 0.d0
      k = 0
      do ii = 1, n_sites
         k = k + 1
         do a = 1, 3
            mu_tot(a) = mu_tot(a) + dot_product(w(:, ii), der(a, :, k))
         end do
         k = k + n_neigh(ii) - 1
      end do
      wv(:, 1, :) = w
      call get_soap_central_hessian(n_sites, n_neigh, n_atom_pairs, n_species, mask, rjs, thetas, phis, &
                                    l_max, rcut_hard, atom_sigma_t, atom_sigma_t_scaling, &
                                    comp, compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, &
                                    soap, 1, wv, hess)
      soap_hessian_enabled = .false.
      dmu_dr = 0.d0
      call accumulate_dmu_dr(V, der, hess, n_neigh, nl_atom, n_atoms, dmu_dr)

      st%mu_hist(:, st%head) = mu_tot
      call mad_ir_evaluate(st, escale, l_out, lam)
      f = 0.d0
      call mad_ir_forces(lam, dmu_dr, f)
   end subroutine analytic

!  Loss at the current positions, older frames untouched.
   real(dp) function loss_now()
      real(dp) :: mt(3), lam(3), l
      call dipole(mt)
      st%mu_hist(:, st%head) = mt
      call mad_ir_evaluate(st, escale, l, lam)
      loss_now = l
   end function loss_now

   logical function stable(ai)
      integer, intent(in) :: ai
      integer :: jj
      real(dp) :: v3(3), rr
      stable = .true.
      do jj = 1, n_atoms
         if (jj == ai) cycle
         v3 = pos(:, jj) - pos(:, ai)
         rr = dsqrt(sum(v3*v3))
         if (abs(rr - rcut_max) < 0.05d0) stable = .false.
         if (rr < 0.5d0) stable = .false.
      end do
   end function stable

   subroutine fdcheck(hh, worst, ndone)
      real(dp), intent(in) :: hh
      real(dp), intent(out) :: worst
      integer, intent(out) :: ndone
      integer :: bb
      real(dp) :: save_pos(3), lp, lm, fd

      worst = 0.d0
      ndone = 0
      do aj = 1, n_atoms
!        every atom whose displacement cannot reorder a neighbour list
         if (.not. stable(aj)) cycle
         do bb = 1, 3
            save_pos = pos(:, aj)
            pos(bb, aj) = save_pos(bb) + hh
            call build_neighbours()
            lp = loss_now()
            pos(:, aj) = save_pos
            pos(bb, aj) = save_pos(bb) - hh
            call build_neighbours()
            lm = loss_now()
            pos(:, aj) = save_pos
            call build_neighbours()
            fd = -(lp - lm)/(2.d0*hh)
            worst = max(worst, abs(fd - forces(bb, aj))/scale)
         end do
         ndone = ndone + 1
      end do
   end subroutine fdcheck

end program madverify
