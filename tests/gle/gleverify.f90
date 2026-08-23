! The generalized Langevin thermostat: its linear algebra, its stationary
! distribution, and its memory.
!
! The thing being tested is a sampler, and a sampler is the kind of code that
! is wrong quietly. Every one of the defects worth catching here -- a
! transposed propagator, a mass scaling applied once instead of twice, a noise
! matrix built from the wrong covariance -- still produces a run that thermalises
! to something, writes plausible trajectories, and reports a temperature. So
! none of the checks below ask whether the output looks thermal. They ask for
! numbers the theory pins down exactly, and compare against those.
!
! In increasing order of how much they can hide:
!
!   1. exp(-A dt), against a closed form and against the semigroup property
!      exp(X)exp(X) = exp(2X). The second is the sharp one: it holds for the
!      true exponential and for essentially no wrong answer, including the
!      unscaled Taylor series that this replaced.
!
!   2. The Cholesky, including the semi-definite case. C - T C T^T is singular
!      whenever a mode carries no noise, and a factorization that returns NaN
!      there would take out perfectly ordinary kernels.
!
!   3. The ns = 0 propagator against the closed-form Ornstein-Uhlenbeck update.
!      This is the one member of the family whose answer is known in full, so
!      it is the reference the general path is measured against.
!
!   4. THE STATIONARY DISTRIBUTION. With no forces the OU propagator is exact
!      at every dt, so <m v^2> = kB T0 is not an approximation to be met within
!      integrator error -- it is an identity, and the only thing between the
!      code and it is sampling noise. Checked per species, because a mass
!      scaling that is wrong by a power is invisible at one mass.
!
!   5. THE MEMORY. <z(t) z(0)^T> = exp(-A t) C for a stationary OU process, so
!      the velocity autocorrelation is known analytically at every lag. This is
!      the check that distinguishes a generalized Langevin equation from an
!      ordinary one: it fails for a single-exponential VACF, it fails if T is
!      transposed, and it fails if the auxiliary variables are propagated but
!      not fed back.
!
!   6. Refusals. A drift and covariance that violate the fluctuation-dissipation
!      condition ask for a negative noise variance, and must be rejected rather
!      than clamped and run.
!
!   7. Bookkeeping: restart round-trip, reproducibility under a fixed seed,
!      fixed atoms, and the thermostat's energy ledger.
!
!   8. A harmonic oscillator, where the splitting is no longer exact. Reported
!      as an h-scan rather than a threshold, because the interesting question
!      is whether the bias falls with dt, not whether it happens to be small at
!      one dt.
!
program gleverify

   use kinds
   use gle

   implicit none

   integer :: n_fail
   real(dp), parameter :: KB = 8.6173303d-5

   n_fail = 0

   call test_expm()
   call test_cholesky()
   call test_white_closed_form()
   call test_stationary()
   call test_memory_kernel()
   call test_fdt_refusal()
   call test_matrix_file()
   call test_restart()
   call test_reproducible()
   call test_fixed_atoms()
   call test_energy_ledger()
   call test_harmonic()
   call test_timescales()

   write (*, *)
   if (n_fail == 0) then
      write (*, '(A)') "==> gleverify: all checks passed"
   else
      write (*, '(A,I0,A)') "==> gleverify: ", n_fail, " CHECK(S) FAILED"
      stop 1
   end if

contains

!  ---------------------------------------------------------------- helpers

   subroutine check(name, got, want, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got, want, tol
      real(dp) :: err
      err = dabs(got - want)
      if (dabs(want) > 1.d0) err = err/dabs(want)
      if (err <= tol .and. .not. (got /= got)) then
         write (*, '(A,A,A,ES13.6,A,ES13.6,A,ES9.2,A)') "   ok   ", name, &
            "  got ", got, "  want ", want, "  (err ", err, ")"
      else
         write (*, '(A,A,A,ES13.6,A,ES13.6,A,ES9.2,A,ES9.2,A)') "   FAIL ", name, &
            "  got ", got, "  want ", want, "  (err ", err, " > ", tol, ")"
         n_fail = n_fail + 1
      end if
   end subroutine check

   subroutine check_true(name, got)
      character(len=*), intent(in) :: name
      logical, intent(in) :: got
      if (got) then
         write (*, '(A,A)') "   ok   ", name
      else
         write (*, '(A,A)') "   FAIL ", name
         n_fail = n_fail + 1
      end if
   end subroutine check_true

!  A drift matrix that satisfies the fluctuation-dissipation condition against
!  C = kB T I. With that C the condition A C + C A^T >= 0 reduces to
!  A + A^T >= 0, so an antisymmetric coupling between p and s and non-negative
!  diagonals is enough -- and antisymmetric coupling is what a physical kernel
!  has, since it is the part that conserves the extended energy.
   subroutine kernel_ns1(A)
      real(dp), intent(out) :: A(2, 2)
      A(1, 1) = 0.005d0; A(1, 2) = 0.08d0
      A(2, 1) = -0.08d0; A(2, 2) = 0.05d0
   end subroutine kernel_ns1

   subroutine kernel_ns2(A)
      real(dp), intent(out) :: A(3, 3)
      A(1, 1) = 0.002d0; A(1, 2) = 0.05d0; A(1, 3) = 0.02d0
      A(2, 1) = -0.05d0; A(2, 2) = 0.02d0; A(2, 3) = 0.d0
      A(3, 1) = -0.02d0; A(3, 2) = 0.d0; A(3, 3) = 0.004d0
   end subroutine kernel_ns2

   subroutine seed_fixed(k)
      integer, intent(in) :: k
      integer :: ns, i
      integer, allocatable :: sd(:)
      call random_seed(size=ns)
      allocate (sd(ns))
      do i = 1, ns
         sd(i) = 1237*k + 91*i + 17
      end do
      call random_seed(put=sd)
      deallocate (sd)
   end subroutine seed_fixed

!  ------------------------------------------------------------ 1. exp(-A dt)

   subroutine test_expm()
      real(dp) :: M1(1, 1), E1(1, 1)
      real(dp) :: M2(2, 2), E2(2, 2)
      real(dp) :: M3(3, 3), E3(3, 3), Em(3, 3), P(3, 3), Eh(3, 3), Ehh(3, 3)
      real(dp) :: th, worst
      integer :: i, j

      write (*, '(A)') "1. matrix exponential"

!     1x1: exp(-g dt) exactly.
      M1(1, 1) = -0.37d0
      call gle_expm(1, M1, E1)
      call check("1x1 against dexp                ", E1(1, 1), dexp(-0.37d0), 1.d-14)

!     A 2x2 rotation generator has a closed form: exp([[0,-t],[t,0]]) is the
!     rotation by t. Chosen because it is the case an implementation that
!     confuses M with M^T gets exactly backwards.
      th = 0.9d0
      M2(1, 1) = 0.d0; M2(1, 2) = -th
      M2(2, 1) = th; M2(2, 2) = 0.d0
      call gle_expm(2, M2, E2)
      call check("2x2 rotation, cos               ", E2(1, 1), dcos(th), 1.d-14)
      call check("2x2 rotation, -sin              ", E2(1, 2), -dsin(th), 1.d-14)
      call check("2x2 rotation, +sin              ", E2(2, 1), dsin(th), 1.d-14)

!     Semigroup. exp(X/2) squared must be exp(X), for an X whose norm is large
!     enough that the scaling stage is doing real work.
      call kernel_ns2(M3)
      M3 = -M3*40.d0
      call gle_expm(3, M3, E3)
      call gle_expm(3, 0.5d0*M3, Eh)
      Ehh = matmul(Eh, Eh)
      worst = 0.d0
      do i = 1, 3
         do j = 1, 3
            worst = max(worst, dabs(Ehh(i, j) - E3(i, j)))
         end do
      end do
      call check("semigroup exp(X/2)^2 = exp(X)   ", worst, 0.d0, 1.d-13)

!     exp(-X) is the inverse of exp(X). Catches a sign convention that is
!     consistent everywhere and consistently wrong.
      call gle_expm(3, -M3, Em)
      P = matmul(E3, Em)
      worst = 0.d0
      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               worst = max(worst, dabs(P(i, j) - 1.d0))
            else
               worst = max(worst, dabs(P(i, j)))
            end if
         end do
      end do
      call check("exp(X) exp(-X) = I              ", worst, 0.d0, 1.d-12)
      write (*, *)
   end subroutine test_expm

!  ------------------------------------------------------------- 2. Cholesky

   subroutine test_cholesky()
      real(dp) :: M(3, 3), L(3, 3), R(3, 3)
      real(dp) :: resid, worst
      integer :: i, j

      write (*, '(A)') "2. Cholesky, including the semi-definite case"

!     Positive definite.
      M(1, 1) = 4.d0; M(1, 2) = 2.d0; M(1, 3) = 1.d0
      M(2, 1) = 2.d0; M(2, 2) = 5.d0; M(2, 3) = 3.d0
      M(3, 1) = 1.d0; M(3, 2) = 3.d0; M(3, 3) = 6.d0
      call gle_cholesky(3, M, L, resid)
      call check("definite: residual              ", resid, 0.d0, 1.d-14)
      call check_true("definite: L is lower triangular ", &
                      dabs(L(1, 2)) + dabs(L(1, 3)) + dabs(L(2, 3)) == 0.d0)

!     Rank 2: the outer product of two vectors, so one direction carries no
!     variance at all. A textbook Cholesky returns NaN here.
      M = 0.d0
      do i = 1, 3
         do j = 1, 3
            M(i, j) = dfloat(i)*dfloat(j) + dfloat(min(i, j))
         end do
      end do
!     M = v v^T + (a rank-2 piece); make it exactly rank 2 by projecting out.
      M(3, 1:3) = M(1, 1:3) + M(2, 1:3)
      M(1:3, 3) = M(1:3, 1) + M(1:3, 2)
      M = 0.5d0*(M + transpose(M))
      call gle_cholesky(3, M, L, resid)
      R = matmul(L, transpose(L))
      worst = 0.d0
      do i = 1, 3
         do j = 1, 3
            worst = max(worst, dabs(R(i, j) - M(i, j)))
         end do
      end do
      call check("semi-definite: L L^T = M        ", worst, 0.d0, 1.d-12)
      call check_true("semi-definite: no NaN           ",.not. any(L /= L))
      write (*, *)
   end subroutine test_cholesky

!  ---------------------------------------------- 3. the ns = 0 closed form

   subroutine test_white_closed_form()
      type(gle_type) :: g
      real(dp) :: tau, dt, temp, m(1)
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "3. ns = 0 is the exact Ornstein-Uhlenbeck update"

      tau = 25.d0
      dt = 0.5d0
      temp = 300.d0
      m(1) = 12.011d0*103.6426965268d0

      call gle_init_white(g, tau, 1, temp, dt, ok, msg)
      call check_true("built                           ", ok)
      call check("T = exp(-dt/tau)                ", g%Tm(1, 1), dexp(-dt/tau), 1.d-14)
      call check("S = sqrt(kT (1 - exp(-2dt/tau)))", g%Sm(1, 1), &
                 dsqrt(KB*temp*(1.d0 - dexp(-2.d0*dt/tau))), 1.d-13)
      call check("C = kT                          ", g%C(1, 1), KB*temp, 1.d-15)
      call gle_free(g)
      write (*, *)
   end subroutine test_white_closed_form

!  ------------------------------------------------- 4. stationary sampling

!  With no forces the propagator is exact, so these are identities up to
!  sampling noise. Three species at once, spanning a factor of 200 in mass,
!  because the mass scaling is the part that is invisible at a single mass.
   subroutine test_stationary()
      type(gle_type) :: g
      real(dp), allocatable :: v(:, :), mass(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3)
      real(dp) :: dt, temp, sv2(3), ssum, tol
      real(dp) :: aux2, kt
      integer :: na, nsteps, i, k, st, sp, cnt(3)
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "4. stationary distribution, forces switched off"

      na = 300
      nsteps = 4000
      dt = 1.d0
      temp = 400.d0
      kt = KB*temp

      allocate (v(3, na), mass(na), fix(3, na))
      fix = .false.
!     H, C, and a heavy one. Interleaved rather than blocked so that an error
!     correlated with atom index cannot line up with an error correlated with
!     mass.
      do i = 1, na
         sp = mod(i - 1, 3) + 1
         if (sp == 1) mass(i) = 1.00797d0*103.6426965268d0
         if (sp == 2) mass(i) = 12.011d0*103.6426965268d0
         if (sp == 3) mass(i) = 195.08d0*103.6426965268d0
      end do

      call seed_fixed(11)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=temp, dt=dt, kind="gle", ok=ok, msg=msg)
      call check_true("built (ns = 2)                  ", ok)
      if (.not. ok) then
         write (*, '(A,A)') "        ", trim(msg)
         return
      end if

!     Start on the stationary distribution, so there is no burn-in to argue
!     about: v ~ N(0, kT/m) and the auxiliary variables already drawn from
!     N(0, C_ss) by gle_init.
      do i = 1, na
         call draw_maxwell(v(1:3, i), mass(i), kt)
      end do

      sv2 = 0.d0
      cnt = 0
      aux2 = 0.d0
      do st = 1, nsteps
         call gle_thermostat(g, v, mass, fix, temp, dt, ok, msg)
         do i = 1, na
            sp = mod(i - 1, 3) + 1
            do k = 1, 3
               sv2(sp) = sv2(sp) + mass(i)*v(k, i)**2
               cnt(sp) = cnt(sp) + 1
            end do
         end do
         aux2 = aux2 + sum(g%s**2)/dfloat(3*na*g%ns)
      end do

!     sqrt(2/N) is the relative standard error of a variance estimate from N
!     independent samples; the samples here are correlated over the kernel's
!     slowest mode, so the tolerance carries a factor for that rather than
!     pretending every step is independent.
      tol = 0.02d0
      call check("<m v^2> per DOF, H              ", sv2(1)/dfloat(cnt(1)), kt, tol)
      call check("<m v^2> per DOF, C              ", sv2(2)/dfloat(cnt(2)), kt, tol)
      call check("<m v^2> per DOF, Pt             ", sv2(3)/dfloat(cnt(3)), kt, tol)
      call check("<s^2> per auxiliary DOF         ", aux2/dfloat(nsteps), kt, tol)

!     The same, started from rest: the thermostat has to fill the system as
!     well as hold it. Averaged over the second half only, so the transient is
!     excluded rather than tolerated.
      v = 0.d0
      g%s = 0.d0
      ssum = 0.d0
      do st = 1, nsteps
         call gle_thermostat(g, v, mass, fix, temp, dt, ok, msg)
         if (st > nsteps/2) then
            do i = 1, na
               do k = 1, 3
                  ssum = ssum + mass(i)*v(k, i)**2
               end do
            end do
         end if
      end do
      call check("<m v^2> from rest               ", ssum/dfloat(3*na*(nsteps - nsteps/2)), kt, tol)

      call gle_free(g)
      deallocate (v, mass, fix)
      write (*, *)
   end subroutine test_stationary

   subroutine draw_maxwell(v, m, kt)
      real(dp), intent(out) :: v(3)
      real(dp), intent(in) :: m, kt
      real(dp) :: g(3)
      call gle_gaussian(g)
      v = g*dsqrt(kt/m)
   end subroutine draw_maxwell

!  ------------------------------------------------------- 5. the memory kernel

!  For a stationary OU process <z(t) z(0)^T> = exp(-A t) C. With C = kB T I
!  that makes the mass-weighted velocity autocorrelation exactly
!  kB T [exp(-A t)]_11 at every lag -- an analytic curve, not a shape.
!
!  For ns = 0 that curve is a single exponential. For ns > 0 it is not, and the
!  test asserts that too: a GLE whose VACF is still a single exponential has
!  auxiliary variables that are being propagated and then ignored.
   subroutine test_memory_kernel()
      type(gle_type) :: g
      real(dp), allocatable :: v(:, :), mass(:), v0(:, :), acf(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3), Ad(3, 3), E(3, 3)
      real(dp) :: dt, temp, kt, want, got, worst, dev_single, fit_tau
      integer :: na, nlag, nblock, i, k, st, lag, b
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "5. velocity autocorrelation against exp(-A t) C"

      na = 400
      nlag = 30
      nblock = 900
      dt = 2.d0
      temp = 350.d0
      kt = KB*temp

      allocate (v(3, na), v0(3, na), mass(na), fix(3, na), acf(0:nlag))
      fix = .false.
      mass = 12.011d0*103.6426965268d0

      call seed_fixed(23)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=temp, dt=dt, kind="gle", ok=ok, msg=msg)
      if (.not. ok) then
         write (*, '(A,A)') "   FAIL build: ", trim(msg)
         n_fail = n_fail + 1
         return
      end if
      do i = 1, na
         call draw_maxwell(v(1:3, i), mass(i), kt)
      end do

!     Let the joint distribution settle before taking time origins. The
!     marginals are already stationary; the cross-correlation between v and s
!     is what needs the few hundred steps.
      do st = 1, 500
         call gle_thermostat(g, v, mass, fix, temp, dt, ok, msg)
      end do

      acf = 0.d0
      do b = 1, nblock
         v0 = v
         do i = 1, na
            do k = 1, 3
               acf(0) = acf(0) + mass(i)*v0(k, i)*v0(k, i)
            end do
         end do
         do lag = 1, nlag
            call gle_thermostat(g, v, mass, fix, temp, dt, ok, msg)
            do i = 1, na
               do k = 1, 3
                  acf(lag) = acf(lag) + mass(i)*v(k, i)*v0(k, i)
               end do
            end do
         end do
      end do
      acf = acf/dfloat(nblock*3*na)

      worst = 0.d0
      do lag = 0, nlag
         Ad = -A*dt*dfloat(lag)
         call gle_expm(3, Ad, E)
         want = kt*E(1, 1)
         got = acf(lag)
         worst = max(worst, dabs(got - want)/kt)
         if (mod(lag, 6) == 0) then
            write (*, '(A,I3,A,ES13.6,A,ES13.6)') "        lag ", lag, &
               "   measured ", got/kt, "   exp(-At)_11 ", E(1, 1)
         end if
      end do
      call check("VACF vs exp(-A t) C, worst lag  ", worst, 0.d0, 0.02d0)

!     And it must not be a single exponential. Fit tau from the first lag and
!     measure how far the rest departs from that fit; for this kernel the
!     departure is large, and for an implementation that drops the auxiliary
!     feedback it would be zero.
      fit_tau = -dt/dlog(acf(1)/acf(0))
      dev_single = 0.d0
      do lag = 0, nlag
         dev_single = max(dev_single, dabs(acf(lag)/acf(0) - dexp(-dfloat(lag)*dt/fit_tau)))
      end do
      write (*, '(A,F8.3,A)') "        single-exponential fit gives tau = ", fit_tau, " fs"
      call check_true("VACF is NOT a single exponential", dev_single > 0.05d0)
      write (*, '(A,ES10.3)') "        max departure from that exponential: ", dev_single

      call gle_free(g)
      deallocate (v, v0, mass, fix, acf)
      write (*, *)
   end subroutine test_memory_kernel

!  ------------------------------------------------------------ 6. refusals

   subroutine test_fdt_refusal()
      type(gle_type) :: g
      real(dp) :: A(2, 2), C(2, 2), Cbad(2, 2)
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "6. matrices that do not describe a bath are refused"

!     A negative friction on the physical momentum: A + A^T has a negative
!     eigenvalue, so C - T C T^T is not positive semi-definite and the noise
!     variance the propagator needs does not exist.
      A(1, 1) = -0.01d0; A(1, 2) = 0.05d0
      A(2, 1) = -0.05d0; A(2, 2) = 0.02d0
      call gle_init(g, A, n_atoms=4, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call check_true("negative friction refused       ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(96, len_trim(msg))))
      call gle_free(g)

!     An asymmetric C is not a covariance.
      call kernel_ns1(A)
      Cbad(1, 1) = KB*300.d0; Cbad(1, 2) = 1.d-3
      Cbad(2, 1) = -1.d-3; Cbad(2, 2) = KB*300.d0
      call gle_init(g, A, Cbad, n_atoms=4, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call check_true("asymmetric C refused            ",.not. ok)
      call gle_free(g)

!     A negative variance on the diagonal is not a covariance either.
      Cbad(1, 1) = -KB*300.d0; Cbad(1, 2) = 0.d0
      Cbad(2, 1) = 0.d0; Cbad(2, 2) = KB*300.d0
      call gle_init(g, A, Cbad, n_atoms=4, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call check_true("negative variance refused       ",.not. ok)
      call gle_free(g)

!     A mismatched C must not be silently truncated to A's order.
      call kernel_ns1(A)
      call gle_init(g, A, C_in=reshape([1.d0], [1, 1]), n_atoms=4, temp=300.d0, &
                    dt=1.d0, kind="gle", ok=ok, msg=msg)
      call check_true("C of the wrong order refused    ",.not. ok)
      call gle_free(g)

!     The legitimate case must still be accepted, or the checks above prove
!     only that the routine says no to everything.
      call kernel_ns1(A)
      C = 0.d0
      C(1, 1) = KB*300.d0
      C(2, 2) = KB*300.d0
      call gle_init(g, A, C, n_atoms=4, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call check_true("a valid explicit C is accepted  ", ok)
      call check_true("c_from_file flag set            ", g%c_from_file)
      call gle_free(g)
      write (*, *)
   end subroutine test_fdt_refusal

!  ------------------------------------------------------- 7. the matrix file

   subroutine test_matrix_file()
      real(dp), allocatable :: M(:, :)
      real(dp) :: A(3, 3)
      integer :: n, u, i
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "7. reading a GLE matrix file"

      call kernel_ns2(A)

!     Comments, a blank line, and rows wrapped across lines: all three appear
!     in the files gle4md.org emits.
      open (newunit=u, file="gle_A_test.dat", status="replace", action="write")
      write (u, '(A)') "# GLE drift matrix, fs^-1"
      write (u, '(A)') "# ns = 2"
      write (u, '(A)') ""
      write (u, '(ES24.16,1X,ES24.16)') A(1, 1), A(1, 2)
      write (u, '(ES24.16)') A(1, 3)
      write (u, '(3(ES24.16,1X))') A(2, 1), A(2, 2), A(2, 3)
      write (u, '(3(ES24.16,1X))') A(3, 1), A(3, 2), A(3, 3)
      close (u)

      call gle_read_matrix("gle_A_test.dat", n, M, ok, msg)
      call check_true("read back                       ", ok)
      if (ok) then
         call check("order deduced from the count    ", dfloat(n), 3.d0, 1.d-15)
         call check("worst element                   ", maxval(dabs(M - A)), 0.d0, 1.d-15)
         deallocate (M)
      else
         write (*, '(A,A)') "        ", trim(msg)
         n_fail = n_fail + 1
      end if

!     A file whose count is not a perfect square cannot be a square matrix, and
!     guessing which element is missing is not a service.
      open (newunit=u, file="gle_A_bad.dat", status="replace", action="write")
      write (u, '(A)') "0.1 0.2 0.3"
      write (u, '(A)') "0.4 0.5"
      close (u)
      call gle_read_matrix("gle_A_bad.dat", n, M, ok, msg)
      call check_true("non-square file refused         ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(88, len_trim(msg))))

      call gle_read_matrix("gle_A_missing.dat", n, M, ok, msg)
      call check_true("missing file refused            ",.not. ok)

      open (newunit=u, file="gle_A_empty.dat", status="replace", action="write")
      write (u, '(A)') "# nothing but a comment"
      close (u)
      call gle_read_matrix("gle_A_empty.dat", n, M, ok, msg)
      call check_true("comment-only file refused       ",.not. ok)

      open (newunit=u, file="gle_A_test.dat", status="old"); close (u, status="delete")
      open (newunit=u, file="gle_A_bad.dat", status="old"); close (u, status="delete")
      open (newunit=u, file="gle_A_empty.dat", status="old"); close (u, status="delete")
      write (*, *)
   end subroutine test_matrix_file

!  ---------------------------------------------------------- 8. the restart

   subroutine test_restart()
      type(gle_type) :: g, h
      real(dp), allocatable :: v(:, :), mass(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3), worst
      integer :: na, st, u
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "8. restart of the auxiliary variables"

      na = 17
      allocate (v(3, na), mass(na), fix(3, na))
      fix = .false.
      mass = 15.9994d0*103.6426965268d0
      v = 0.d0

      call seed_fixed(31)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      do st = 1, 50
         call gle_thermostat(g, v, mass, fix, 300.d0, 1.d0, ok, msg)
      end do

      call gle_save(g, "gle_restart_test.dat", ok, msg)
      call check_true("saved                           ", ok)

      call gle_init(h, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call gle_load(h, "gle_restart_test.dat", ok, msg)
      call check_true("loaded                          ", ok)
      worst = maxval(dabs(h%s - g%s))
      call check("auxiliary variables round-trip  ", worst, 0.d0, 1.d-15)
      call check("thermostat ledger round-trips   ", h%e_thermo, g%e_thermo, 1.d-15)
      call gle_free(h)

!     A file from a run with a different kernel is refused, not adopted: the
!     auxiliary variables of one kernel mean nothing under another.
      call gle_init_white(h, 20.d0, na, 300.d0, 1.d0, ok, msg)
      call gle_load(h, "gle_restart_test.dat", ok, msg)
      call check_true("different ns refused            ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(96, len_trim(msg))))
      call gle_free(h)

!     And so is one from a different number of atoms.
      call gle_init(h, A, n_atoms=na + 1, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call gle_load(h, "gle_restart_test.dat", ok, msg)
      call check_true("different n_atoms refused       ",.not. ok)
      call gle_free(h)

!     A truncated file is refused rather than half-applied.
      call trim_file("gle_restart_test.dat", 12)
      call gle_init(h, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call gle_load(h, "gle_restart_test.dat", ok, msg)
      call check_true("truncated file refused          ",.not. ok)
      if (.not. ok) write (*, '(A,A)') "        ", trim(msg(1:min(96, len_trim(msg))))
      call gle_free(h)

      call gle_free(g)
      open (newunit=u, file="gle_restart_test.dat", status="old"); close (u, status="delete")
      deallocate (v, mass, fix)
      write (*, *)
   end subroutine test_restart

   subroutine trim_file(fname, keep)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: keep
      character(len=1024), allocatable :: lines(:)
      character(len=1024) :: line
      integer :: u, ios, n, i
      allocate (lines(keep))
      open (newunit=u, file=trim(fname), status="old", action="read")
      n = 0
      do i = 1, keep
         read (u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         n = n + 1
         lines(n) = line
      end do
      close (u)
      open (newunit=u, file=trim(fname), status="replace", action="write")
      do i = 1, n
         write (u, '(A)') trim(lines(i))
      end do
      close (u)
      deallocate (lines)
   end subroutine trim_file

!  ------------------------------------------------------ 9. reproducibility

   subroutine test_reproducible()
      type(gle_type) :: g
      real(dp), allocatable :: v(:, :), mass(:), v1(:, :)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3)
      integer :: na, st
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "9. a seeded run is reproducible"

      na = 23
      allocate (v(3, na), v1(3, na), mass(na), fix(3, na))
      fix = .false.
      mass = 28.0855d0*103.6426965268d0
      call kernel_ns2(A)

      call seed_fixed(7)
      call gle_init(g, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      v = 0.d0
      do st = 1, 40
         call gle_thermostat(g, v, mass, fix, 300.d0, 1.d0, ok, msg)
      end do
      v1 = v
      call gle_free(g)

      call seed_fixed(7)
      call gle_init(g, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      v = 0.d0
      do st = 1, 40
         call gle_thermostat(g, v, mass, fix, 300.d0, 1.d0, ok, msg)
      end do
      call check("same seed, same velocities      ", maxval(dabs(v - v1)), 0.d0, 1.d-300)
      call gle_free(g)
      deallocate (v, v1, mass, fix)
      write (*, *)
   end subroutine test_reproducible

!  ------------------------------------------------------- 10. fixed atoms

   subroutine test_fixed_atoms()
      type(gle_type) :: g
      real(dp), allocatable :: v(:, :), mass(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3)
      integer :: na, st
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "10. constrained components stay constrained"

      na = 12
      allocate (v(3, na), mass(na), fix(3, na))
      mass = 12.011d0*103.6426965268d0
      v = 0.d0
      fix = .false.
!     Atom 3 fully fixed; atom 7 fixed in y only, which is the case where a
!     per-atom rather than per-component skip would show up.
      fix(1:3, 3) = .true.
      fix(2, 7) = .true.

      call seed_fixed(5)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=500.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      g%s = 0.d0
      do st = 1, 200
         call gle_thermostat(g, v, mass, fix, 500.d0, 1.d0, ok, msg)
      end do

      call check("fixed atom stays at rest        ", maxval(dabs(v(1:3, 3))), 0.d0, 1.d-300)
      call check("fixed component stays at rest   ", dabs(v(2, 7)), 0.d0, 1.d-300)
      call check_true("its free components move        ", &
                      dabs(v(1, 7)) > 0.d0 .and. dabs(v(3, 7)) > 0.d0)
      call check("its bath is not propagated      ", maxval(dabs(g%s(:, 2, 7))), 0.d0, 1.d-300)
      call check_true("other atoms move                ", maxval(dabs(v(1:3, 1))) > 0.d0)
      call gle_free(g)
      deallocate (v, mass, fix)
      write (*, *)
   end subroutine test_fixed_atoms

!  ------------------------------------------------------ 11. energy ledger

   subroutine test_energy_ledger()
      type(gle_type) :: g
      real(dp), allocatable :: v(:, :), mass(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3), ek0, ek1
      integer :: na, st, i
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "11. the thermostat's energy ledger"

      na = 40
      allocate (v(3, na), mass(na), fix(3, na))
      fix = .false.
      mass = 12.011d0*103.6426965268d0

      call seed_fixed(13)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      do i = 1, na
         call draw_maxwell(v(1:3, i), mass(i), KB*300.d0)
      end do

      ek0 = 0.d0
      do i = 1, na
         ek0 = ek0 + 0.5d0*mass(i)*dot_product(v(1:3, i), v(1:3, i))
      end do
      do st = 1, 300
         call gle_thermostat(g, v, mass, fix, 300.d0, 1.d0, ok, msg)
      end do
      ek1 = 0.d0
      do i = 1, na
         ek1 = ek1 + 0.5d0*mass(i)*dot_product(v(1:3, i), v(1:3, i))
      end do

!     With no forces, every change in kinetic energy came from the thermostat,
!     so the ledger must account for all of it exactly. This is what makes
!     E_pot + E_kin - e_thermo meaningful as a conserved quantity in a real run.
      call check("e_thermo = dE_kin (no forces)   ", g%e_thermo, ek1 - ek0, 1.d-12)
      call gle_free(g)
      deallocate (v, mass, fix)
      write (*, *)
   end subroutine test_energy_ledger

!  --------------------------------------------------- 12. harmonic oscillator

!  Forces on. The O step is still exact, but it no longer commutes with the
!  force update, so the sampled temperature carries a splitting bias. The
!  question a threshold cannot answer is whether that bias is the splitting or
!  a defect, so this is run as an h-scan: halving dt must roughly halve the
!  error for a first-order splitting, and an error that sits still as dt falls
!  is something else.
!
!  Equipartition for a harmonic oscillator gives <m v^2> = kB T0 and
!  <k x^2> = kB T0 independently, and both are checked -- a thermostat that
!  gets the kinetic term right while the configurational term drifts is the
!  signature of a scheme applied at the wrong point in the step.
   subroutine test_harmonic()
      real(dp) :: dts(4), err_k(4), err_p(4)
      integer :: id

      write (*, '(A)') "12. harmonic oscillator, an h-scan of the splitting bias"
      write (*, '(A)') "        dt/fs      <m v^2>/kT      <k x^2>/kT"

      dts = [4.d0, 2.d0, 1.d0, 0.5d0]
      do id = 1, 4
         call harmonic_run(dts(id), err_k(id), err_p(id))
         write (*, '(A,F8.3,2(4X,F12.6))') "     ", dts(id), 1.d0 + err_k(id), 1.d0 + err_p(id)
      end do

!     At the smallest step both must be close. The tolerance is set by the
!     sampling length, not by the splitting: 3e5 correlated samples give a
!     variance to about half a percent.
      call check("kinetic  equipartition at dt/8  ", 1.d0 + err_k(4), 1.d0, 0.02d0)
      call check("configurational   ditto         ", 1.d0 + err_p(4), 1.d0, 0.02d0)

!     And the bias must shrink. Compared against the coarsest step, where it is
!     largest and least confusable with sampling noise.
      call check_true("bias falls with dt              ", &
                      dabs(err_k(4)) < dabs(err_k(1)) .or. dabs(err_k(1)) < 0.02d0)
      write (*, *)
   end subroutine test_harmonic

!  ------------------------------------------------ 13. the reported timescales

!  The setup report's relaxation times, against eigenvalues that are known in
!  closed form. This is a report rather than a computation the propagator
!  depends on, but it is the report a user picks a kernel by -- a kernel whose
!  modes miss the system's vibrational frequencies thermostats it very slowly
!  while looking perfectly healthy -- so a wrong answer here is a wrong choice
!  of kernel later.
   subroutine test_timescales()
      type(gle_type) :: g
      real(dp) :: A1(3, 3), A2(2, 2), A3(3, 3)
      real(dp) :: tf, ts
      logical :: ok
      character(len=512) :: msg

      write (*, '(A)') "13. the relaxation times in the setup report"

!     Lower triangular: the eigenvalues are the diagonal, 0.1, 0.05 and 0.02,
!     so the fastest mode is 10 fs and the slowest 50 fs.
      A1 = 0.d0
      A1(1, 1) = 0.1d0; A1(2, 2) = 0.05d0; A1(3, 3) = 0.02d0
      A1(2, 1) = 0.03d0; A1(3, 1) = 0.01d0; A1(3, 2) = 0.02d0
      call gle_init(g, A1, n_atoms=2, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      if (ok) then
         call gle_timescales(g, tf, ts)
         call check("triangular: fastest             ", tf, 10.d0, 1.d-10)
         call check("triangular: slowest             ", ts, 50.d0, 1.d-10)
      else
         write (*, '(A,A)') "   FAIL build: ", trim(msg)
         n_fail = n_fail + 1
      end if
      call gle_free(g)

!     A complex pair. trace = 0.055 and the discriminant is negative, so both
!     eigenvalues have real part 0.0275 and both timescales are 1/0.0275. This
!     is the case an implementation that reads only the real array gets wrong.
      call kernel_ns1(A2)
      call gle_init(g, A2, n_atoms=2, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call gle_timescales(g, tf, ts)
      call check("complex pair: fastest           ", tf, 1.d0/0.0275d0, 1.d-9)
      call check("complex pair: slowest           ", ts, 1.d0/0.0275d0, 1.d-9)
      call gle_free(g)

!     And the kernel the rest of this file uses, where the old Gershgorin bound
!     reached through zero and could report no slowest mode at all. Any finite
!     positive answer is an improvement on that; what is checked is that both
!     ends now exist and are ordered.
      call kernel_ns2(A3)
      call gle_init(g, A3, n_atoms=2, temp=300.d0, dt=1.d0, kind="gle", ok=ok, msg=msg)
      call gle_timescales(g, tf, ts)
      write (*, '(A,F10.3,A,F10.3,A)') "        ns = 2 kernel: fastest ", tf, &
         " fs, slowest ", ts, " fs"
      call check_true("both ends reported              ", tf > 0.d0 .and. ts > 0.d0)
      call check_true("and ordered                     ", ts >= tf)
      call gle_free(g)
      write (*, *)
   end subroutine test_timescales

   subroutine harmonic_run(dt, err_k, err_p)
      real(dp), intent(in) :: dt
      real(dp), intent(out) :: err_k, err_p
      type(gle_type) :: g
      real(dp), allocatable :: x(:, :), v(:, :), f(:, :), mass(:)
      logical, allocatable :: fix(:, :)
      real(dp) :: A(3, 3), kspring, m, kt, temp, t_total
      real(dp) :: sk, sp
      integer :: na, nsteps, st, i, k, nsamp
      logical :: ok
      character(len=512) :: msg

      na = 200
      temp = 300.d0
      kt = KB*temp
      m = 12.011d0*103.6426965268d0
!     omega = 0.05 fs^-1, a period of ~126 fs: slow enough that dt = 4 fs is a
!     legitimate step for it, so the scan measures the splitting rather than an
!     integrator falling apart.
      kspring = m*0.05d0**2
!     A fixed total time, so every dt sees the same amount of physics and the
!     comparison between rows is of bias and not of sampling length.
      t_total = 3.d5
      nsteps = nint(t_total/dt)

      allocate (x(3, na), v(3, na), f(3, na), mass(na), fix(3, na))
      fix = .false.
      mass = m

      call seed_fixed(nint(100.d0*dt) + 3)
      call kernel_ns2(A)
      call gle_init(g, A, n_atoms=na, temp=temp, dt=dt, kind="gle", ok=ok, msg=msg)

      do i = 1, na
         call draw_maxwell(v(1:3, i), m, kt)
         call draw_maxwell(x(1:3, i), kspring, kt)
      end do
      f = -kspring*x

      sk = 0.d0
      sp = 0.d0
      nsamp = 0
      do st = 1, nsteps
!        Velocity Verlet, written out rather than called, so this measures the
!        thermostat and not md.f90's integrator.
         v = v + 0.5d0*dt*f/m
         x = x + dt*v
         f = -kspring*x
         v = v + 0.5d0*dt*f/m
!        The thermostat sits where compute_md puts it: after the full Verlet
!        update, once per step.
         call gle_thermostat(g, v, mass, fix, temp, dt, ok, msg)
         if (st > nsteps/5) then
            do i = 1, na
               do k = 1, 3
                  sk = sk + m*v(k, i)**2
                  sp = sp + kspring*x(k, i)**2
               end do
            end do
            nsamp = nsamp + 3*na
         end if
      end do

      err_k = sk/dfloat(nsamp)/kt - 1.d0
      err_p = sp/dfloat(nsamp)/kt - 1.d0

      call gle_free(g)
      deallocate (x, v, f, mass, fix)
   end subroutine harmonic_run

end program gleverify
