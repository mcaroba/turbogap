! Generalized Langevin dynamics by Markovian embedding.
!
! THE PROBLEM. A Langevin thermostat couples each momentum to a bath through
! an instantaneous friction and white noise. The generalized Langevin equation
! replaces the friction with a memory kernel and the white noise with coloured
! noise correlated on the same timescale,
!
!   m rdot'(t) = -dV/dr - int_{-inf}^{t} K(t - t') rdot(t') dt' + eta(t)      (1)
!   <eta(t) eta(t')> = kB T K(t - t')                                          (2)
!
! (2) is the fluctuation-dissipation theorem, and it is what makes (1) sample
! the canonical distribution rather than merely damp. The trouble is the
! integral: evaluated literally it needs the whole velocity history at every
! step, so cost and storage grow with the length of the run.
!
! THE FIX. Couple the physical momentum to a handful of extra variables that
! obey ordinary, memoryless dynamics, and let the memory emerge from having
! integrated them out. For one Cartesian degree of freedom, with p the
! mass-scaled momentum and s a vector of ns auxiliary momenta,
!
!   d ( p )   = ( -dV/dr ) dt  -  A ( p ) dt  +  B dW                          (3)
!     ( s )     (    0   )         p ( s )
!
! where A_p is (ns+1) x (ns+1). Nothing in (3) looks back. Eliminating s from
! (3) reproduces (1) with a kernel that is a sum of ns decaying exponentials
! (complex pairs give damped oscillations), so any kernel that can be fitted by
! such a sum is available at the cost of ns extra numbers per degree of freedom.
! ns = 4..8 covers several decades of frequency; the fitted matrices published
! at gle4md.org are of exactly this shape and can be used directly.
!
! WHAT IS PROPAGATED. Between force evaluations (3) is a linear
! Ornstein-Uhlenbeck process, and an OU process has an exact finite-time
! propagator -- there is no integrator error in the thermostat itself, at any
! dt:
!
!   z <- T z + S xi,   T = exp(-A_p dt),   S S^T = C_p - T C_p T^T,  xi ~ N(0,1)
!
! with z = (p, s) and C_p the stationary covariance the process is asked to
! have. For canonical sampling at temperature T0, C_p = kB T0 I. Supplying
! some other C_p is how the quantum thermostats work: they give each normal
! mode the energy of a quantum oscillator rather than kB T0.
!
! B never appears. The pair (A_p, C_p) determines it through
! B B^T = A_p C_p + C_p A_p^T, and going through S = chol(C_p - T C_p T^T)
! instead is both cheaper and better behaved -- it is exact for any dt, whereas
! integrating B dW over a step is not.
!
! MASS SCALING. z(0) is sqrt(m) v, not m v and not v. In that variable the
! stationary covariance is kB T0 for every atom regardless of mass, so one
! (A_p, C_p) pair thermostats every species correctly and the matrices from a
! published fit mean the same thing here as they did there. The auxiliary
! variables are carried in the same units.
!
! WHERE IT ACTS. gle_thermostat is called once per step, at the point in
! compute_md where berendsen and bussi are applied -- after the velocity-Verlet
! update, when the velocities are synchronous with the positions of the
! previous step. That makes the splitting first order in dt (an "OBA" scheme),
! which is what the other thermostats at that call site are as well. The O part
! is exact whatever dt is, so a run with no forces samples kB T0 exactly; with
! forces there is the usual O(dt) splitting bias in the sampled kinetic energy,
! and tests/gle measures it rather than assuming it away.
!
! ns = 0 IS ORDINARY LANGEVIN. With no auxiliary variables A_p is the 1x1
! matrix [gamma] and the propagator collapses to
!
!   v <- e^{-gamma dt} v + sqrt( kB T0 (1 - e^{-2 gamma dt}) / m ) xi
!
! which is the exact Langevin/OU update. `thermostat = langevin` builds that
! case internally from tau_t (gamma = 1/tau_t) and takes no matrix files. It
! is not a degenerate afterthought: it is the one member of the family with a
! closed-form answer for everything, which is what makes it the reference the
! tests measure the general path against.
!
! MPI. compute_md runs the whole MD step on rank 0 and broadcasts positions and
! velocities afterwards, so the auxiliary variables and the random draws live
! on rank 0 alone. There is nothing to reduce and nothing to keep in step, and
! a run on n ranks reproduces the same seeded trajectory as a run on one.
!
! RANDOM NUMBERS. Drawn from the intrinsic generator, so they belong to the
! stream random_seed fixes and a seeded run is reproducible. md.f90's
! gaussian_deviates does the same thing for the same reason; this module keeps
! its own copy rather than using it so that gle.f90 uses no other module in
! this tree, which is what lets tests/gle build the module under test on its
! own. The only external it needs is LAPACK's dgeev, and that only to report
! the kernel's relaxation times at setup -- nothing the propagator does depends
! on it.
module gle

   use kinds

   implicit none

!  eV/K. The same value the rest of the code uses.
   real(dp), parameter :: GLE_KB = 8.6173303d-5

!  Terms in the Taylor series for the scaled matrix exponential. With the
!  argument scaled to infinity-norm <= 1/2 the truncation error after 18 terms
!  is below the double-precision epsilon, so the squaring stage sees a matrix
!  that is exact to round-off.
   integer, parameter :: GLE_EXP_TERMS = 18

!  Guard on the scaling loop, so a matrix carrying a NaN cannot spin forever.
   integer, parameter :: GLE_EXP_MAX_SCALE = 64

   type :: gle_type
      logical :: active = .false.
!     ns auxiliary momenta per Cartesian degree of freedom; nd = ns + 1 is the
!     dimension of everything below.
      integer :: ns = 0
      integer :: nd = 1
      integer :: n_atoms = 0
!     Which flavour built this: "gle" from matrix files, "langevin" from tau_t.
      character(len=32) :: kind = "gle"
!     Was C supplied by the user? If it was, it is a statement about the bath
!     at one temperature and must not be rescaled when t_beg /= t_end; if it
!     was not, C = kB T0 I and follows the ramp.
      logical :: c_from_file = .false.
!     What T and S were last built for. Rebuilt when either moves, which is
!     O(nd^3) on a matrix of order ten and therefore free.
      real(dp) :: dt_built = -1.d0
      real(dp) :: temp_built = -1.d0
!     Running sum of the kinetic energy the thermostat has put in (positive) or
!     taken out (negative), in eV. E_pot + E_kin - e_thermo is the conserved
!     quantity of the extended system and is what a GLE run should be watched
!     with; plain E_pot + E_kin is not conserved and is not meant to be.
      real(dp) :: e_thermo = 0.d0
!     Largest |L L^T - M| seen when factoring C - T C T^T, relative to max|C|.
!     Kept because it is the one number that says whether the supplied A and C
!     are compatible; see gle_rebuild.
      real(dp) :: chol_residual = 0.d0
      real(dp), allocatable :: A(:, :)       ! (nd,nd) drift, fs^-1
      real(dp), allocatable :: C(:, :)       ! (nd,nd) stationary covariance, eV
      real(dp), allocatable :: Tm(:, :)      ! (nd,nd) exp(-A dt)
      real(dp), allocatable :: Sm(:, :)      ! (nd,nd) S S^T = C - T C T^T
      real(dp), allocatable :: s(:, :, :)    ! (ns,3,n_atoms) aux momenta, mass-scaled
   end type gle_type

   type(gle_type), save :: gle_state

contains

!
! Standard normal deviates by the polar Box-Muller transform. See the module
! header for why this is a copy of md.f90's rather than a use of it.
   subroutine gle_gaussian(g)

      implicit none

      real(dp), intent(out) :: g(:)
      real(dp) :: u(1:2), r, f
      integer :: n, i

      n = size(g)
      i = 1
      do while (i <= n)
         r = 2.d0
         do while (r >= 1.d0 .or. r == 0.d0)
            call random_number(u)
            u = 2.d0*u - 1.d0
            r = u(1)*u(1) + u(2)*u(2)
         end do
         f = dsqrt(-2.d0*dlog(r)/r)
         g(i) = u(1)*f
         if (i + 1 <= n) g(i + 1) = u(2)*f
         i = i + 2
      end do

   end subroutine gle_gaussian

!
! exp(M) for a small dense M, by scaling and squaring around a Taylor series.
!
! The scaling is not optional. exp of a matrix with a large norm computed
! directly from its Taylor series loses every digit to cancellation between
! terms that are individually enormous -- and a GLE drift matrix spanning
! several decades of frequency is exactly that kind of matrix, since the whole
! point of it is that the fastest and slowest modes are far apart. Halving
! until the norm is below 1/2 and then squaring back up keeps every
! intermediate O(1).
   subroutine gle_expm(n, M, E)

      implicit none

      integer, intent(in) :: n
      real(dp), intent(in) :: M(n, n)
      real(dp), intent(out) :: E(n, n)
      real(dp) :: Ms(n, n), Tk(n, n), nrm
      integer :: i, k, j

      nrm = 0.d0
      do i = 1, n
         nrm = max(nrm, sum(dabs(M(i, 1:n))))
      end do

      j = 0
      do while (nrm > 0.5d0 .and. j < GLE_EXP_MAX_SCALE)
         nrm = 0.5d0*nrm
         j = j + 1
      end do

      Ms = M/(2.d0**j)

      E = 0.d0
      Tk = 0.d0
      do i = 1, n
         E(i, i) = 1.d0
         Tk(i, i) = 1.d0
      end do
      do k = 1, GLE_EXP_TERMS
         Tk = matmul(Tk, Ms)/dfloat(k)
         E = E + Tk
      end do

      do k = 1, j
         E = matmul(E, E)
      end do

   end subroutine gle_expm

!
! Cholesky factor of a symmetric positive-semi-definite M, L L^T = M, with L
! lower triangular.
!
! Semi-definite, not definite: C - T C T^T is singular whenever the process has
! a mode the noise does not reach, and it tends to singular as dt -> 0. A
! textbook Cholesky takes a square root of a pivot that is zero to round-off
! and returns a NaN, or fails outright, on cases that are perfectly well posed.
! A pivot at or below the tolerance means that direction carries no noise, so
! the honest factor has a zero column there -- and any L with L L^T = M will do,
! since only the product enters the propagator.
!
! resid comes back as max |L L^T - M|. It is not a diagnostic for the
! factorization, which is exact when it succeeds; it is a test of the CALLER's
! matrices. C - T C T^T is positive semi-definite for every dt if and only if
! A C + C A^T is, which is the fluctuation-dissipation condition (2). A user
! who writes down an A and a C that do not satisfy it has asked for a process
! with negative noise variance, and the only place that shows up is here, as a
! clamped negative pivot and a residual that does not go away.
   subroutine gle_cholesky(n, M, L, resid)

      implicit none

      integer, intent(in) :: n
      real(dp), intent(in) :: M(n, n)
      real(dp), intent(out) :: L(n, n)
      real(dp), intent(out) :: resid
      real(dp) :: acc, tol, dmax
      integer :: i, j, k

      dmax = 0.d0
      do i = 1, n
         dmax = max(dmax, dabs(M(i, i)))
      end do
!     Scale-free floor. An absolute one would be meaningless: C is in eV and
!     its entries are ~1e-2 at room temperature.
      tol = dmax*1.d-13

      L = 0.d0
      do i = 1, n
         acc = M(i, i)
         do k = 1, i - 1
            acc = acc - L(i, k)*L(i, k)
         end do
         if (acc <= tol) then
!           No noise in this direction. Leave column i of L at zero.
            cycle
         end if
         L(i, i) = dsqrt(acc)
         do j = i + 1, n
            acc = M(j, i)
            do k = 1, i - 1
               acc = acc - L(j, k)*L(i, k)
            end do
            L(j, i) = acc/L(i, i)
         end do
      end do

      resid = 0.d0
      do i = 1, n
         do j = 1, n
            acc = 0.d0
            do k = 1, n
               acc = acc + L(i, k)*L(j, k)
            end do
            resid = max(resid, dabs(acc - M(i, j)))
         end do
      end do
      if (dmax > 0.d0) resid = resid/dmax

   end subroutine gle_cholesky

!
! Read a square matrix from a text file.
!
! The format is the one the fitted matrices at gle4md.org come in: free-format
! numbers, row-major, whitespace or newline separated, with `#` comment lines
! carrying the metadata. The order is not stated in the file and is not asked
! for on the input line -- it is deduced from the count, because a matrix file
! whose declared size disagrees with its contents is a failure mode with no
! upside, and the count settles it with no room for disagreement.
!
! Units are NOT converted and NOT guessed. A is in fs^-1 and C is in eV. The
! gle4md generator emits matrices in a unit system chosen at download time, and
! a wrong guess here would produce a thermostat that runs at the wrong
! temperature or the wrong rate while looking entirely healthy -- so the
! keyword documentation states the units and this routine takes the file at its
! word.
   subroutine gle_read_matrix(fname, n, M, ok, msg)

      implicit none

      character(len=*), intent(in) :: fname
      integer, intent(out) :: n
      real(dp), allocatable, intent(out) :: M(:, :)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: buf(:)
      real(dp) :: x
      integer :: u, ios, cnt, cap, i, j, k
      character(len=1024) :: line
      character(len=1024) :: rest
      logical :: ex

      ok = .false.
      msg = ""
      n = 0

      inquire (file=trim(fname), exist=ex)
      if (.not. ex) then
         msg = "GLE matrix file not found: "//trim(fname)
         return
      end if

      cap = 256
      allocate (buf(cap))
      cnt = 0

      open (newunit=u, file=trim(fname), status="old", action="read", iostat=ios)
      if (ios /= 0) then
         msg = "could not open GLE matrix file: "//trim(fname)
         deallocate (buf)
         return
      end if

      do
         read (u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == "#" .or. line(1:1) == "!") cycle
!        Free-format read of however many numbers the line holds. Reading them
!        one at a time keeps the row length out of the format entirely, so a
!        file wrapped at 80 columns and one with a row per line both work.
         rest = line
         do
            read (rest, *, iostat=ios) x
            if (ios /= 0) exit
            cnt = cnt + 1
            if (cnt > cap) then
               call gle_grow(buf, cap)
            end if
            buf(cnt) = x
!           Drop the token just consumed and go round again.
            call gle_drop_token(rest)
            if (len_trim(rest) == 0) exit
         end do
      end do
      close (u)

      if (cnt == 0) then
         msg = "GLE matrix file holds no numbers: "//trim(fname)
         deallocate (buf)
         return
      end if

      n = nint(dsqrt(dfloat(cnt)))
      if (n*n /= cnt) then
         write (msg, '(A,I0,A,A)') "GLE matrix file is not square: ", cnt, &
            " numbers in ", trim(fname)
         deallocate (buf)
         n = 0
         return
      end if

      allocate (M(n, n))
      k = 0
      do i = 1, n
         do j = 1, n
            k = k + 1
            M(i, j) = buf(k)
         end do
      end do
      deallocate (buf)
      ok = .true.

   end subroutine gle_read_matrix

!  Double the scratch buffer. Split out only to keep gle_read_matrix readable.
   subroutine gle_grow(buf, cap)
      implicit none
      real(dp), allocatable, intent(inout) :: buf(:)
      integer, intent(inout) :: cap
      real(dp), allocatable :: tmp(:)
      allocate (tmp(2*cap))
      tmp(1:cap) = buf(1:cap)
      call move_alloc(tmp, buf)
      cap = 2*cap
   end subroutine gle_grow

!  Remove the leading whitespace-delimited token from s, in place.
   subroutine gle_drop_token(s)
      implicit none
      character(len=*), intent(inout) :: s
      integer :: i
      s = adjustl(s)
      i = index(trim(s), " ")
      if (i == 0) then
         s = ""
      else
         s = adjustl(s(i:))
      end if
   end subroutine gle_drop_token

!
! Build the propagator for the current dt and target temperature.
!
! Called from gle_thermostat whenever either has moved, which on a run with a
! fixed step and a fixed temperature is once. C follows the temperature ramp
! only when this module built it; a C that came from a file describes a bath at
! one temperature and rescaling it would silently turn a quantum thermostat
! into something with no name.
   subroutine gle_rebuild(this, dt, temp, ok, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: temp
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: W(:, :), M(:, :)
      integer :: i, nd

      ok = .true.
      msg = ""
      nd = this%nd

      if (.not. this%c_from_file) then
         this%C = 0.d0
         do i = 1, nd
            this%C(i, i) = GLE_KB*temp
         end do
      end if

      allocate (W(nd, nd), M(nd, nd))

!     T = exp(-A dt)
      W = -this%A*dt
      call gle_expm(nd, W, this%Tm)

!     M = C - T C T^T, symmetrized. The symmetrization is not cosmetic: the
!     factorization below reads only the lower triangle of M via M(j,i), so an
!     asymmetry of order round-off in the product would be taken as real and
!     would show up in the residual as if the user's matrices were at fault.
      M = this%C - matmul(this%Tm, matmul(this%C, transpose(this%Tm)))
      M = 0.5d0*(M + transpose(M))

      call gle_cholesky(nd, M, this%Sm, this%chol_residual)

      if (this%chol_residual > 1.d-8) then
         ok = .false.
         write (msg, '(A,ES10.3,A)') &
            "GLE: C - T C T^T is not positive semi-definite (residual ", &
            this%chol_residual, "). The supplied A and C do not satisfy the "// &
            "fluctuation-dissipation condition A C + C A^T >= 0."
      end if

      deallocate (W, M)

      this%dt_built = dt
      this%temp_built = temp

   end subroutine gle_rebuild

!
! Set up from an explicit drift matrix, and optionally an explicit covariance.
!
! A_in is (ns+1) x (ns+1) with the physical momentum first; C_in, if present,
! must match it. n_atoms sizes the auxiliary array, which is (ns,3,n_atoms) and
! therefore 8*ns*3*n_atoms bytes -- 1.7 MB for ns = 8 at 7000 atoms, which is
! why the auxiliary variables are stored and not recomputed.
   subroutine gle_init(this, A_in, C_in, n_atoms, temp, dt, kind, ok, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      real(dp), intent(in) :: A_in(:, :)
      real(dp), intent(in), optional :: C_in(:, :)
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: temp
      real(dp), intent(in) :: dt
      character(len=*), intent(in) :: kind
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: nd, i

      ok = .false.
      msg = ""

      nd = size(A_in, 1)
      if (size(A_in, 2) /= nd) then
         msg = "GLE: drift matrix A is not square."
         return
      end if
      if (nd < 1) then
         msg = "GLE: drift matrix A is empty."
         return
      end if
      if (present(C_in)) then
         if (size(C_in, 1) /= nd .or. size(C_in, 2) /= nd) then
            write (msg, '(A,I0,A,I0,A)') "GLE: C is ", size(C_in, 1), "x", &
               size(C_in, 2), " but A sets the order; they must agree."
            return
         end if
      end if
      if (n_atoms < 1) then
         msg = "GLE: n_atoms must be positive."
         return
      end if
      if (dt <= 0.d0) then
         msg = "GLE: the time step must be positive."
         return
      end if
      if (.not. present(C_in) .and. temp <= 0.d0) then
         msg = "GLE: the target temperature must be positive."
         return
      end if

      call gle_free(this)

      this%ns = nd - 1
      this%nd = nd
      this%n_atoms = n_atoms
      this%kind = kind
      this%c_from_file = present(C_in)

      allocate (this%A(nd, nd), this%C(nd, nd), this%Tm(nd, nd), this%Sm(nd, nd))
      this%A = A_in
      if (present(C_in)) then
         this%C = C_in
!        A covariance is a covariance: symmetric, with non-negative variances.
!        Both are cheap to check and both are things a hand-edited file gets
!        wrong in ways that do not otherwise surface until the residual test in
!        gle_rebuild, where the message points at the wrong thing.
         do i = 1, nd
            if (this%C(i, i) < 0.d0) then
               write (msg, '(A,I0,A)') "GLE: C has a negative variance on the diagonal (element ", &
                  i, "); it is not a covariance matrix."
               call gle_free(this)
               return
            end if
         end do
         if (maxval(dabs(this%C - transpose(this%C))) > 1.d-10*max(1.d0, maxval(dabs(this%C)))) then
            msg = "GLE: C is not symmetric; it is not a covariance matrix."
            call gle_free(this)
            return
         end if
      else
         this%C = 0.d0
      end if

      if (this%ns > 0) then
         allocate (this%s(this%ns, 3, n_atoms))
         this%s = 0.d0
      end if

      call gle_rebuild(this, dt, temp, ok, msg)
      if (.not. ok) then
         call gle_free(this)
         return
      end if

      call gle_draw_aux(this)

      this%e_thermo = 0.d0
      this%active = .true.
      ok = .true.

   end subroutine gle_init

!
! The ns = 0 case, built from a friction rather than a file: ordinary Langevin
! with gamma = 1/tau_t. See the module header for why this shares every line
! below gle_init with the general path instead of being its own integrator.
   subroutine gle_init_white(this, tau_t, n_atoms, temp, dt, ok, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      real(dp), intent(in) :: tau_t
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: temp
      real(dp), intent(in) :: dt
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp) :: A(1, 1)

      ok = .false.
      msg = ""
      if (tau_t <= 0.d0) then
         msg = "GLE: tau_t must be positive for thermostat = langevin."
         return
      end if

      A(1, 1) = 1.d0/tau_t
      call gle_init(this, A, n_atoms=n_atoms, temp=temp, dt=dt, kind="langevin", &
                    ok=ok, msg=msg)

   end subroutine gle_init_white

!
! Draw the auxiliary variables from their stationary marginal.
!
! Starting them at zero is a transient the length of the slowest mode, during
! which the thermostat is not the one that was asked for; drawing them from
! N(0, C_ss) starts the bath already equilibrated. The mass-scaled convention
! is what makes this possible without knowing the masses.
   subroutine gle_draw_aux(this)

      implicit none

      type(gle_type), intent(inout) :: this
      real(dp), allocatable :: Css(:, :), L(:, :), g(:)
      real(dp) :: resid
      integer :: i, k, a, ns

      ns = this%ns
      if (ns < 1) return

      allocate (Css(ns, ns), L(ns, ns), g(ns))
      Css = this%C(2:this%nd, 2:this%nd)
      call gle_cholesky(ns, Css, L, resid)

      do a = 1, this%n_atoms
         do k = 1, 3
            call gle_gaussian(g)
            do i = 1, ns
               this%s(i, k, a) = dot_product(L(i, 1:i), g(1:i))
            end do
         end do
      end do

      deallocate (Css, L, g)

   end subroutine gle_draw_aux

!
! One thermostat step: z <- T z + S xi, for every Cartesian degree of freedom.
!
! target_temp and dt are passed in rather than stored because both can move
! during a run -- a temperature ramp changes the first and variable_time_step
! the second -- and a propagator built for a dt the run is no longer taking is
! wrong in a way nothing downstream can detect.
!
! Fixed components are left alone. velocity_verlet holds them at zero, and a
! thermostat that heated them would be quietly undoing the constraint; their
! auxiliary variables are not propagated either, so a component that is later
! released does not carry a bath that has been running without it.
   subroutine gle_thermostat(this, velocities, masses, fix_atom, target_temp, dt, ok, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      real(dp), intent(inout) :: velocities(:, :)
      real(dp), intent(in) :: masses(:)
      logical, intent(in) :: fix_atom(:, :)
      real(dp), intent(in) :: target_temp
      real(dp), intent(in) :: dt
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: z(:), zn(:), g(:)
      real(dp) :: sq, ek_before, ek_after
      integer :: a, k, i, j, n, nd

      ok = .true.
      msg = ""
      if (.not. this%active) return

      n = min(size(masses), this%n_atoms)
      nd = this%nd

!     Rebuild only when something the propagator depends on has actually moved.
!     The comparison is exact on purpose: these are copies of the same values,
!     not the results of separate arithmetic, so a tolerance would only make it
!     possible to miss a small deliberate change.
      if (dt /= this%dt_built .or. &
          (.not. this%c_from_file .and. target_temp /= this%temp_built)) then
         call gle_rebuild(this, dt, target_temp, ok, msg)
         if (.not. ok) return
      end if

      ek_before = 0.d0
      do a = 1, n
         ek_before = ek_before + 0.5d0*masses(a)*dot_product(velocities(1:3, a), velocities(1:3, a))
      end do

      allocate (z(nd), zn(nd), g(nd))

      do a = 1, n
         sq = dsqrt(masses(a))
         do k = 1, 3
            if (fix_atom(k, a)) cycle
            z(1) = sq*velocities(k, a)
            do i = 2, nd
               z(i) = this%s(i - 1, k, a)
            end do

            call gle_gaussian(g)
!           zn = T z + S g. S is lower triangular, so the noise sum runs only
!           to i.
            do i = 1, nd
               zn(i) = 0.d0
               do j = 1, nd
                  zn(i) = zn(i) + this%Tm(i, j)*z(j)
               end do
               do j = 1, i
                  zn(i) = zn(i) + this%Sm(i, j)*g(j)
               end do
            end do

            velocities(k, a) = zn(1)/sq
            do i = 2, nd
               this%s(i - 1, k, a) = zn(i)
            end do
         end do
      end do

      deallocate (z, zn, g)

      ek_after = 0.d0
      do a = 1, n
         ek_after = ek_after + 0.5d0*masses(a)*dot_product(velocities(1:3, a), velocities(1:3, a))
      end do
      this%e_thermo = this%e_thermo + (ek_after - ek_before)

   end subroutine gle_thermostat

!
! Release everything. Safe to call on a state that was never initialised, which
! is what lets the error paths in gle_init use it as an unwind.
   subroutine gle_free(this)

      implicit none

      type(gle_type), intent(inout) :: this

      if (allocated(this%A)) deallocate (this%A)
      if (allocated(this%C)) deallocate (this%C)
      if (allocated(this%Tm)) deallocate (this%Tm)
      if (allocated(this%Sm)) deallocate (this%Sm)
      if (allocated(this%s)) deallocate (this%s)
      this%active = .false.
      this%ns = 0
      this%nd = 1
      this%n_atoms = 0
      this%dt_built = -1.d0
      this%temp_built = -1.d0
      this%e_thermo = 0.d0
      this%chol_residual = 0.d0
      this%c_from_file = .false.

   end subroutine gle_free

!
! Write the auxiliary variables so a restart continues the same bath.
!
! Without this a restarted run begins with s drawn afresh, which is not wrong
! -- the draw is from the stationary distribution -- but it discards the
! correlation between the bath and the atoms, and for a memory kernel whose
! slowest mode is longer than the restart interval that correlation is the
! whole of what the thermostat was doing.
   subroutine gle_save(this, fname, ok, msg)

      implicit none

      type(gle_type), intent(in) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: u, ios, a, k

      ok = .false.
      msg = ""
      if (.not. this%active) then
         msg = "GLE: nothing to save."
         return
      end if

      open (newunit=u, file=trim(fname), status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         msg = "GLE: could not write "//trim(fname)
         return
      end if

      write (u, '(A)') "# TurboGAP GLE thermostat state"
      write (u, '(A)') "# ns  n_atoms  kind"
      write (u, '(I0,1X,I0,1X,A)') this%ns, this%n_atoms, trim(this%kind)
      write (u, '(A)') "# accumulated thermostat energy (eV)"
      write (u, '(ES24.16)') this%e_thermo
      write (u, '(A)') "# auxiliary momenta, mass-scaled: atom, direction, s(1:ns)"
      do a = 1, this%n_atoms
         do k = 1, 3
            if (this%ns > 0) then
               write (u, '(I0,1X,I0,1X,*(ES24.16,1X))') a, k, this%s(1:this%ns, k, a)
            else
               write (u, '(I0,1X,I0)') a, k
            end if
         end do
      end do
      close (u)

      ok = .true.

   end subroutine gle_save

!
! Read auxiliary variables back.
!
! A file that does not describe this run is REFUSED, not adopted. ns is the
! shape of the memory kernel and n_atoms the shape of the system; either one
! differing means the file belongs to a different simulation, and continuing
! with it would run a bath fitted for one kernel against the propagator of
! another. Refusal is not fatal -- the caller falls back to a fresh draw from
! the stationary distribution, which is a legitimate start -- so the message
! matters more than the status.
   subroutine gle_load(this, fname, ok, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      character(len=*), intent(in) :: fname
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
      integer :: u, ios, a, k, ns_f, na_f, aa, kk, i
      character(len=32) :: kind_f
      character(len=1024) :: line
      real(dp), allocatable :: tmp(:)
      logical :: ex

      ok = .false.
      msg = ""

      inquire (file=trim(fname), exist=ex)
      if (.not. ex) then
         msg = "GLE: no restart file "//trim(fname)
         return
      end if

      open (newunit=u, file=trim(fname), status="old", action="read", iostat=ios)
      if (ios /= 0) then
         msg = "GLE: could not read "//trim(fname)
         return
      end if

      call gle_next_data_line(u, line, ios)
      if (ios /= 0) then
         msg = "GLE: restart file is empty: "//trim(fname)
         close (u)
         return
      end if
      kind_f = ""
      read (line, *, iostat=ios) ns_f, na_f, kind_f
      if (ios /= 0) read (line, *, iostat=ios) ns_f, na_f
      if (ios /= 0) then
         msg = "GLE: restart header is not 'ns n_atoms [kind]': "//trim(fname)
         close (u)
         return
      end if

      if (ns_f /= this%ns .or. na_f /= this%n_atoms) then
         write (msg, '(A,I0,A,I0,A,I0,A,I0,A)') &
            "GLE: restart file describes a different run (ns = ", ns_f, &
            ", n_atoms = ", na_f, "; this run has ns = ", this%ns, &
            ", n_atoms = ", this%n_atoms, "). Starting a fresh bath."
         close (u)
         return
      end if

      call gle_next_data_line(u, line, ios)
      if (ios == 0) read (line, *, iostat=ios) this%e_thermo

      allocate (tmp(max(this%ns, 1)))
      do a = 1, this%n_atoms
         do k = 1, 3
            call gle_next_data_line(u, line, ios)
            if (ios /= 0) then
               write (msg, '(A,I0,A)') "GLE: restart file ends after ", &
                  3*(a - 1) + k - 1, " degrees of freedom. Starting a fresh bath."
               close (u)
               deallocate (tmp)
               return
            end if
            if (this%ns > 0) then
               read (line, *, iostat=ios) aa, kk, tmp(1:this%ns)
            else
               read (line, *, iostat=ios) aa, kk
            end if
            if (ios /= 0 .or. aa /= a .or. kk /= k) then
               write (msg, '(A,I0,A,I0,A)') &
                  "GLE: restart file is out of order at atom ", a, ", direction ", k, &
                  ". Starting a fresh bath."
               close (u)
               deallocate (tmp)
               return
            end if
            if (this%ns > 0) then
               do i = 1, this%ns
                  this%s(i, k, a) = tmp(i)
               end do
            end if
         end do
      end do
      close (u)
      deallocate (tmp)

      ok = .true.
      msg = "GLE: bath resumed from "//trim(fname)

   end subroutine gle_load

!  Next line that is neither blank nor a comment.
   subroutine gle_next_data_line(u, line, ios)
      implicit none
      integer, intent(in) :: u
      character(len=*), intent(out) :: line
      integer, intent(out) :: ios
      do
         read (u, '(A)', iostat=ios) line
         if (ios /= 0) return
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == "#" .or. line(1:1) == "!") cycle
         return
      end do
   end subroutine gle_next_data_line

!
! Everything a driver needs to start a GLE run: read the matrices, build the
! propagator, and resume the bath if there is one to resume.
!
! Called lazily on the first thermostat step rather than from turbogap_setup,
! for the same reason mad_ir's setup is: n_atoms is not known until the system
! is read, and the report belongs next to the thing that consumes it.
!
! `resumed` is separate from `ok` on purpose. Failing to read a restart file is
! not a failure of the run -- a fresh bath drawn from the stationary
! distribution is a legitimate start -- so it comes back as resumed = .false.
! with a message, while ok = .false. is reserved for the matrices themselves
! being unusable, which is fatal.
   subroutine gle_setup(this, kind, a_file, c_file, restart_file, do_restart, &
                        tau_t, n_atoms, temp, dt, ok, resumed, msg)

      implicit none

      type(gle_type), intent(inout) :: this
      character(len=*), intent(in) :: kind
      character(len=*), intent(in) :: a_file
      character(len=*), intent(in) :: c_file
      character(len=*), intent(in) :: restart_file
      logical, intent(in) :: do_restart
      real(dp), intent(in) :: tau_t
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: temp
      real(dp), intent(in) :: dt
      logical, intent(out) :: ok
      logical, intent(out) :: resumed
      character(len=*), intent(out) :: msg
      real(dp), allocatable :: A(:, :), C(:, :)
      integer :: na, nc
      logical :: lok
      character(len=512) :: lmsg

      ok = .false.
      resumed = .false.
      msg = ""

      if (trim(kind) == "langevin") then
         call gle_init_white(this, tau_t, n_atoms, temp, dt, ok, msg)
      else
         call gle_read_matrix(a_file, na, A, lok, lmsg)
         if (.not. lok) then
            msg = lmsg
            return
         end if
         if (len_trim(c_file) > 0 .and. trim(c_file) /= "none") then
            call gle_read_matrix(c_file, nc, C, lok, lmsg)
            if (.not. lok) then
               msg = lmsg
               deallocate (A)
               return
            end if
            call gle_init(this, A, C, n_atoms, temp, dt, "gle", ok, msg)
            deallocate (C)
         else
            call gle_init(this, A, n_atoms=n_atoms, temp=temp, dt=dt, kind="gle", &
                          ok=ok, msg=msg)
         end if
         deallocate (A)
      end if
      if (.not. ok) return

      if (do_restart .and. len_trim(restart_file) > 0 .and. trim(restart_file) /= "none") then
         call gle_load(this, restart_file, resumed, lmsg)
         msg = lmsg
      end if

   end subroutine gle_setup

!
! The setup report.
!
! A GLE is specified by a matrix, and the one thing a matrix does not tell its
! user is whether it covers the timescales of the system it is about to
! thermostat. The relaxation-time bounds are printed for that reason, next to
! the step the run is actually taking: a kernel whose slowest mode is shorter
! than dt is being sampled by an integrator that cannot see it, and a kernel
! whose fastest mode is much longer than the run is a thermostat that will not
! act within it. Neither is an error and neither is refused -- both are
! legitimate things to ask for -- but both are things to have been told.
   subroutine gle_report(this, dt, temp, resumed, resume_msg)

      implicit none

      type(gle_type), intent(in) :: this
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: temp
      logical, intent(in) :: resumed
      character(len=*), intent(in) :: resume_msg
      real(dp) :: tau_fast, tau_slow

      call gle_timescales(this, tau_fast, tau_slow)

      write (*, '(A)') '                                             |'
      write (*, '(A)') ' Generalized Langevin thermostat:            |'
      write (*, '(A,A12,A)') '  *) kernel:            ', trim(this%kind), '         |'
      write (*, '(A,I12,A)') '  *) auxiliary DOF:     ', this%ns, '         |'
      write (*, '(A,F12.4,A)') '  *) time step:         ', dt, ' fs      |'
      if (this%c_from_file) then
         write (*, '(A)') '  *) covariance:          from file          |'
      else
         write (*, '(A,F12.2,A)') '  *) target T:          ', temp, ' K       |'
      end if
      if (tau_fast > 0.d0) then
         write (*, '(A,F12.4,A)') '  *) fastest mode:      ', tau_fast, ' fs      |'
      end if
      if (tau_slow > 0.d0) then
         write (*, '(A,F12.4,A)') '  *) slowest mode:      ', tau_slow, ' fs      |'
      end if
      if (resumed) then
         write (*, '(A)') '  *) bath resumed from restart               |'
      else if (len_trim(resume_msg) > 0) then
         write (*, '(A)') '  *) fresh bath                              |'
      end if
      write (*, '(A)') ' ............................................|'
      if (tau_fast > 0.d0 .and. tau_fast < dt) then
         write (*, *) '                                       |'
         write (*, *) 'WARNING: the kernel has a mode faster   |  <-- WARNING'
         write (*, *) 'than the time step. The propagator is   |'
         write (*, *) 'still exact, but the forces cannot see  |'
         write (*, *) 'what it does between steps.             |'
      end if
!     The REASON a bath was not resumed, printed outside the box because it is
!     a sentence and the box is 39 columns. "fresh bath" alone is not enough:
!     a missing file and a file belonging to another run both land here, and
!     only one of those is something the user meant to happen. A refusal that
!     does not say what it refused is indistinguishable from no file at all.
      if (.not. resumed .and. len_trim(resume_msg) > 0) then
         write (*, '(1X,A)') trim(resume_msg)
      end if

   end subroutine gle_report

!
! The timescales the supplied kernel actually covers.
!
! A GLE is specified by a matrix, and the one thing a matrix does not tell its
! user is whether it spans the frequencies of the system it is about to
! thermostat. That is not a cosmetic question. Measured here on 897 atoms of
! amorphous carbon, a kernel whose modes sit at 33 and 167 fs cooled the system
! several times more slowly than a plain Langevin thermostat of nominally
! weaker friction, because carbon's optical phonons are near 25 fs and the
! kernel's spectral weight is nowhere near them. Nothing was wrong; the kernel
! simply did not reach the modes carrying the energy. So the eigenvalues are
! the report that matters, and they are worth an O(nd^3) solve at setup.
!
! Re(lambda) is a rate: 1/Re(lambda) is the relaxation time of that mode. A
! complex pair is an oscillating kernel component whose envelope still decays
! at 1/Re(lambda), so the real parts are what is reported either way.
!
! DGEEV rather than Gershgorin, which this used to use. Gershgorin bounds the
! eigenvalues by discs around the diagonal, and for any kernel with strong
! p-to-s coupling -- which is every kernel worth having, since that coupling IS
! the memory -- the discs reach through zero and the bound degenerates to "the
! slowest mode is somewhere between infinity and nothing". It reported nothing
! useful for exactly the matrices a user needs it for.
   subroutine gle_timescales(this, tau_fast, tau_slow)

      implicit none

      type(gle_type), intent(in) :: this
      real(dp), intent(out) :: tau_fast, tau_slow
      real(dp), allocatable :: Ac(:, :), wr(:), wi(:), work(:), vdummy(:, :)
      real(dp) :: lo, hi
      integer :: nd, lwork, info, i

      tau_fast = -1.d0
      tau_slow = -1.d0
      if (.not. this%active) return

      nd = this%nd
      lwork = 8*nd + 64
      allocate (Ac(nd, nd), wr(nd), wi(nd), work(lwork), vdummy(1, 1))
!     DGEEV overwrites its argument.
      Ac = this%A

      call dgeev('N', 'N', nd, Ac, nd, wr, wi, vdummy, 1, vdummy, 1, work, lwork, info)

      if (info /= 0) then
!        No eigenvalues means no report, not a wrong report. The thermostat
!        itself does not depend on this, so a failure here is not fatal.
         deallocate (Ac, wr, wi, work, vdummy)
         return
      end if

      lo = huge(1.d0)
      hi = 0.d0
      do i = 1, nd
!        A non-positive real part is a mode that does not decay. It cannot
!        occur for a kernel that passed the positive-semi-definite test in
!        gle_rebuild, and if it somehow does it has no timescale to quote.
         if (wr(i) <= 0.d0) cycle
         lo = min(lo, wr(i))
         hi = max(hi, wr(i))
      end do

      if (hi > 0.d0) tau_fast = 1.d0/hi
      if (lo < huge(1.d0)) tau_slow = 1.d0/lo

      deallocate (Ac, wr, wi, work, vdummy)

   end subroutine gle_timescales

end module gle
