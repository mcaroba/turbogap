! The MAD IR spectrum, its loss, and the sensitivity of that loss to the
! newest dipole.
!
! Four things, in increasing order of how much they can hide:
!
!   1. Sizing. n_lag follows from the requested resolution and n_window from
!      n_lag; a sampling interval too coarse for the requested wavenumber must
!      be refused rather than aliased.
!
!   2. Nyquist in practice. A dipole oscillating above the Nyquist limit must
!      show up folded back to 2*nu_nyq - nu0, not near nu0. This is the check
!      that catches a factor-of-two or a 2*pi in the transform phase, which a
!      peak-position test at a single well-behaved frequency will not.
!
!   3. Peak placement. A pure cosine at nu0 must produce its maximum at nu0.
!
!   4. lambda = dL/dmu(newest), against central differences of the loss with
!      respect to the newest stored dipole. This is the quantity the MAD force
!      is built from, and it is the one that has to be right to the last
!      digit -- everything upstream of it is a spectrum nobody differentiates.
!      Both with and without the fitted overall scale, because the scale sits
!      at its own optimum and the envelope-theorem argument that lets dL/dI
!      ignore it is exactly the kind of thing that is wrong silently.
!
program irverify
   use kinds
   use mad_ir
   implicit none

   type(mad_ir_type) :: st
   integer :: n_lag, n_window, n_freq, k, t, ih, c
   real(dp) :: dt, nu_res, nu_max, loss, lambda(3), h, fd, worst, sc, escale
   real(dp) :: lp, lm, nu0, peak, two_pi, amp(3)
   real(dp), allocatable :: nu(:), Iexp(:), wgt(:), mu(:, :)
   logical :: ok
   character(len=256) :: msg

   two_pi = 2.d0*dacos(-1.d0)
!  exp_energy_scales for this observable. The energy is
!  1/2 * escale * sum (I_pred - I_exp)^2, the same form every MAD observable
!  uses, and lambda is its derivative -- so this tests the gradient of the
!  spectral difference, not of something proportional to it.
   escale = 3.7d0

! ------------------------------------------------------------------ 1. sizing
   write (*, '(A)') "1. window sizing"
   dt = 2.d0          ! fs between stored frames
   nu_res = 8.d0      ! cm^-1 resolution wanted
   nu_max = 4000.d0   ! highest wavenumber wanted
   call mad_ir_size_window(dt, nu_res, nu_max, 2, n_lag, n_window, ok, msg)
   write (*, '(A,L1)') "   accepted (dt=2 fs, res=8, max=4000):    ", ok
   write (*, '(A,I0,A,I0)') "   n_lag = ", n_lag, "   n_window = ", n_window
   write (*, '(A,F9.2,A)') "   implied resolution = ", CM_PER_INV_FS/(dfloat(n_lag)*dt), " cm^-1"
   write (*, '(A,F9.1,A)') "   Nyquist            = ", CM_PER_INV_FS/(2.d0*dt), " cm^-1"
!  same resolution but sampling too coarsely for the wanted range
   call mad_ir_size_window(6.d0, nu_res, nu_max, 2, n_lag, n_window, ok, msg)
   write (*, '(A,L1)') "   accepted (dt=6 fs, max=4000):           ", ok
   if (.not. ok) write (*, '(A,A)') "   refused with: ", trim(msg)
   write (*, *)

! ------------------------------------------------------ 2/3. spectral content
!  A grid fine enough to locate a peak, over a range a real fit would use.
   dt = 2.d0
   nu_res = 20.d0
   call mad_ir_size_window(dt, nu_res, 4000.d0, 2, n_lag, n_window, ok, msg)
   n_freq = 400
   allocate (nu(n_freq), Iexp(n_freq), wgt(n_freq))
   do k = 1, n_freq
      nu(k) = 4000.d0*dfloat(k)/dfloat(n_freq)
   end do
   Iexp = 0.d0
   wgt = 1.d0

   allocate (mu(3, n_window))

   write (*, '(A)') "2. a dipole above Nyquist must alias, not vanish"
!  choose the drive so that the fold-back lands at 1500 cm^-1, well inside the
!  grid: a test whose expected answer is off the end of the grid proves nothing
   nu0 = 2.d0*(CM_PER_INV_FS/(2.d0*dt)) - 1500.d0
   call fill_cosine(nu0)
   call load_and_eval()
   peak = nu(maxloc(st%I_calc, 1))
   write (*, '(A,F8.1,A,F8.1,A)') "   driven at ", nu0, " cm^-1 (Nyquist ", &
      CM_PER_INV_FS/(2.d0*dt), ")"
   write (*, '(A,F8.1,A,F8.1)') "   peak found at ", peak, &
      "   expected alias at ", 2.d0*CM_PER_INV_FS/(2.d0*dt) - nu0
   write (*, *)

   write (*, '(A)') "3. peak placement for a resolvable dipole"
   do k = 1, 3
      nu0 = 600.d0*dfloat(k)
      call fill_cosine(nu0)
      call load_and_eval()
      peak = nu(maxloc(st%I_calc, 1))
      write (*, '(A,F8.1,A,F8.1,A,F6.1,A)') "   driven at ", nu0, &
         " cm^-1   peak at ", peak, "   (grid spacing ", nu(2) - nu(1), ")"
   end do
   write (*, *)

! ------------------------------------------------------------- 4. sensitivity
!  The target is a spectrum this code actually produced, from a different
!  trajectory and then distorted. Inventing an I_exp of order unity instead
!  would leave the computed spectrum ten orders of magnitude away, the fitted
!  scale absorbing all of it, and the test measuring almost nothing about the
!  loss. Here the fitted scale comes out O(1), which is the regime a real fit
!  runs in.
   call fill_cosine(1500.d0)
   call load_and_eval_scale(.false.)
   do k = 1, n_freq
      Iexp(k) = 1.15d0*st%I_calc(k) + 0.05d0*maxval(st%I_calc)*dexp(-((nu(k) - 2900.d0)/150.d0)**2)
      wgt(k) = 1.d0
   end do
   call fill_mixed()

   write (*, '(A)') "4. lambda = dE/dmu(newest) vs central differences"
   do c = 1, 2
      call load_and_eval_scale(c == 1)
      write (*, '(A,L1,A,ES11.3,A,ES11.3)') "   match_scale = ", (c == 1), &
         "   energy = ", loss, "   fitted scale = ", st%scale
      do ih = 1, 9
         h = 1.d-1*10.d0**(-0.5d0*(ih - 1))
         worst = 0.d0
         sc = maxval(abs(lambda))
         do t = 1, 3
            st%mu_hist(t, st%head) = mu(t, n_window) + h
            call mad_ir_evaluate(st, escale, lp, lambda)
            st%mu_hist(t, st%head) = mu(t, n_window) - h
            call mad_ir_evaluate(st, escale, lm, lambda)
            st%mu_hist(t, st%head) = mu(t, n_window)
            call mad_ir_evaluate(st, escale, loss, lambda)
            fd = (lp - lm)/(2.d0*h)
            worst = max(worst, abs(fd - lambda(t))/sc)
         end do
         write (*, '(A,ES9.2,A,ES11.3)') "      h = ", h, "   max rel err ", worst
      end do
      write (*, *)
   end do

contains

   subroutine fill_cosine(nu_drive)
      real(dp), intent(in) :: nu_drive
      integer :: tt
!     frame tt is at time (tt-1)*dt; the newest frame is tt = n_window
      do tt = 1, n_window
         mu(1, tt) = dcos(two_pi*nu_drive*dfloat(tt - 1)*dt/CM_PER_INV_FS)
         mu(2, tt) = 0.d0
         mu(3, tt) = 0.d0
      end do
   end subroutine fill_cosine

   subroutine fill_mixed()
      integer :: tt
      do tt = 1, n_window
         amp(1) = dcos(two_pi*1500.d0*dfloat(tt - 1)*dt/CM_PER_INV_FS)
         amp(2) = 0.6d0*dcos(two_pi*2900.d0*dfloat(tt - 1)*dt/CM_PER_INV_FS + 0.7d0)
         amp(3) = 0.3d0*dcos(two_pi*800.d0*dfloat(tt - 1)*dt/CM_PER_INV_FS + 1.9d0)
         mu(1, tt) = amp(1) + 0.2d0*amp(3)
         mu(2, tt) = amp(2) - 0.1d0*amp(1)
         mu(3, tt) = amp(3) + 0.3d0*amp(2)
      end do
   end subroutine fill_mixed

   subroutine load_and_eval()
      call load_and_eval_scale(.true.)
   end subroutine load_and_eval

   subroutine load_and_eval_scale(ms)
      logical, intent(in) :: ms
      integer :: tt
      call mad_ir_init(st, dt, n_lag, n_window, nu, Iexp, wgt, ms, 2.d0, "hann")
      do tt = 1, n_window
         call mad_ir_push(st, mu(:, tt))
      end do
      call mad_ir_evaluate(st, escale, loss, lambda)
   end subroutine load_and_eval_scale

end program irverify
