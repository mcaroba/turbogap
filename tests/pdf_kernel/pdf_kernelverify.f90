! The pdf kernel's recurrence against the Gaussian in quad precision.
!
! pair_kde_row with step builds a pair's row from one exp at the sample c
! nearest r and two multiplications a sample, so it evaluates the Gaussian at
! x(c) + k step; without step it takes an exp at every stored x(l). For each
! grid and width, every pair distance on a fine scan from below the grid to past
! its end: the same window from both; each row within RTOL of the quad Gaussian
! at the points it evaluates, of the row's peak and relative wherever above 1e-6
! of it; and the stored grid within DX_TOL of x(1) + k step. The rows differ
! from each other by more than RTOL (~7e-13 relative at 24 A, sigma 0.05):
! that is the stored x(l)'s rounding seen through the Gaussian's slope, not the
! arithmetic. A grid not evenly spaced must not get a step at all.
program pdf_kernelverify
   use kinds, only: dp
   use exp_utils, only: pdf_kernel_t, make_pdf_kernel, pair_kde_row, uniform_step, PDF_KDE_CUT
   implicit none
   integer, parameter :: qp = selected_real_kind(30)
   real(dp), parameter :: RTOL = 1.d-13
   real(dp), parameter :: DX_TOL = 1.d-13
   integer, parameter :: N_R = 20000
   real(dp), parameter :: r_min(2) = [0.1d0, 0.1d0]
   real(dp), parameter :: r_max(2) = [12.0d0, 24.4d0]
   integer, parameter :: n_samples(2) = [301, 613]
   real(dp), parameter :: sigmas(3) = [0.05d0, 0.1d0, 0.2d0]
   type(pdf_kernel_t) :: kern
   real(dp), allocatable :: direct(:)
   real(dp), allocatable :: recur(:)
   real(dp), allocatable :: bent(:)
   real(qp) :: exact
   real(dp) :: r
   real(dp) :: peak
   real(dp) :: dx
   real(dp) :: worst_abs(2)
   real(dp) :: worst_rel(2)
   integer :: lo_d
   integer :: hi_d
   integer :: lo_r
   integer :: hi_r
   integer :: c
   integer :: g
   integer :: s
   integer :: i
   integer :: l
   integer :: n_fail

   n_fail = 0
   do g = 1, 2
      do s = 1, 3
         call make_pdf_kernel(r_min(g), r_max(g), n_samples(g), sigmas(s), kern)
         if (kern%step /= kern%dr) then
            write (*, '(A,I0,A,F5.2,A)') 'FAIL: grid ', g, ', sigma ', sigmas(s), ' was not found evenly spaced'
            n_fail = n_fail + 1
            cycle
         end if
         dx = 0.d0
         do l = 1, n_samples(g)
            dx = max(dx, abs(kern%x(1) + dfloat(l - 1)*kern%step - kern%x(l)))
         end do
         allocate (direct(1:n_samples(g)), recur(1:n_samples(g)))
         worst_abs = 0.d0
         worst_rel = 0.d0
         do i = 0, N_R
            r = r_min(g) - sigmas(s) + (r_max(g) - r_min(g) + 8.d0*sigmas(s))*dfloat(i)/dfloat(N_R)
            direct = 0.d0
            recur = 0.d0
            call pair_kde_row(r, kern%n_samples, kern%x, kern%r_min, kern%dr, kern%sigma, PDF_KDE_CUT, &
                              lo_d, hi_d, direct)
            call pair_kde_row(r, kern%n_samples, kern%x, kern%r_min, kern%dr, kern%sigma, PDF_KDE_CUT, &
                              lo_r, hi_r, recur, kern%step)
            if (lo_d /= lo_r .or. hi_d /= hi_r) then
               write (*, '(A,F12.6)') 'FAIL: the window moved at r = ', r
               n_fail = n_fail + 1
               exit
            end if
            if (hi_d < lo_d) cycle
            peak = maxval(direct(lo_d:hi_d))
            if (peak <= 0.d0) cycle
!           pair_kde_row's c.
            c = min(hi_d, max(lo_d, lo_d + nint((r - kern%x(lo_d))/kern%step)))
            do l = lo_d, hi_d
               exact = exp(-((real(kern%x(l), qp) - real(r, qp))/real(kern%sigma, qp))**2/2)
               call tally(1, direct(l), exact)
               exact = real(kern%x(c), qp) + real(l - c, qp)*real(kern%step, qp)
               exact = exp(-((exact - real(r, qp))/real(kern%sigma, qp))**2/2)
               call tally(2, recur(l), exact)
            end do
         end do
         write (*, '(A,I4,A,F5.2,A,2ES9.2,A,2ES9.2,A,ES9.2,A)') '  samples ', n_samples(g), ', sigma ', sigmas(s), &
            ': exp, recurrence from quad', worst_abs, ' of the peak,', worst_rel, ' relative; grid', dx, ' A'
         if (maxval(worst_abs) > RTOL .or. maxval(worst_rel) > RTOL) then
            write (*, '(A)') 'FAIL: a row is further from the quad Gaussian than RTOL'
            n_fail = n_fail + 1
         end if
         if (dx > DX_TOL) then
            write (*, '(A)') 'FAIL: the stored grid is further than DX_TOL from evenly spaced'
            n_fail = n_fail + 1
         end if
         deallocate (direct, recur)
      end do
   end do

   call make_pdf_kernel(r_min(1), r_max(1), n_samples(1), 0.1d0, kern)
   bent = kern%x
   bent(5) = bent(5) + 1.d-6
   if (uniform_step(bent, kern%dr) /= 0.d0) then
      write (*, '(A)') 'FAIL: a grid with one sample moved by 1e-6 A was taken as evenly spaced'
      n_fail = n_fail + 1
   end if

   if (n_fail > 0) then
      write (*, '(I0,A)') n_fail, ' failure(s)'
      stop 1
   end if
   write (*, '(A)') 'PASS'

contains

   subroutine tally(k, val, exact)
      integer, intent(in) :: k
      real(dp), intent(in) :: val
      real(qp), intent(in) :: exact

      worst_abs(k) = max(worst_abs(k), real(abs(val - exact), dp)/peak)
      if (val > 1.d-6*peak) worst_rel(k) = max(worst_rel(k), real(abs(val - exact)/exact, dp))
   end subroutine tally

end program pdf_kernelverify
