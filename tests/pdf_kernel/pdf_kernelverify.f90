! The pdf's sample grid and Gaussian kernel against quad precision.
!
! For a 12 A and a 24.4 A pdf grid at three kernel widths: the stored samples
! within DX_TOL of evenly spaced, as r_min + (l - 1/2) dr computed in quad;
! and pair_kde_row, over a fine scan of pair distances from below the grid to
! past its end, within RTOL of the Gaussian in quad precision at the samples
! it uses, of the row's peak and relative wherever above 1e-6 of it. A grid
! built in single precision is ~1e-6 A off even spacing, which at 3 sigma moves
! a pair's kernel by ~1e-5 of itself.
program pdf_kernelverify
   use kinds, only: dp
   use exp_utils, only: pdf_kernel_t, make_pdf_kernel, pair_kde_row, PDF_KDE_CUT
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
   real(dp), allocatable :: row(:)
   real(qp) :: exact
   real(qp) :: ideal
   real(dp) :: r
   real(dp) :: peak
   real(dp) :: dx
   real(dp) :: worst_abs
   real(dp) :: worst_rel
   integer :: lo
   integer :: hi
   integer :: g
   integer :: s
   integer :: i
   integer :: l
   integer :: n_fail

   n_fail = 0
   do g = 1, 2
      do s = 1, 3
         call make_pdf_kernel(r_min(g), r_max(g), n_samples(g), sigmas(s), kern)
         dx = 0.d0
         do l = 1, n_samples(g)
            ideal = real(r_min(g), qp) + (real(l, qp) - 0.5_qp)*(real(r_max(g), qp) - real(r_min(g), qp)) &
                    /real(n_samples(g), qp)
            dx = max(dx, real(abs(real(kern%x(l), qp) - ideal), dp))
         end do
         allocate (row(1:n_samples(g)))
         worst_abs = 0.d0
         worst_rel = 0.d0
         do i = 0, N_R
            r = r_min(g) - sigmas(s) + (r_max(g) - r_min(g) + 8.d0*sigmas(s))*dfloat(i)/dfloat(N_R)
            row = 0.d0
            call pair_kde_row(r, kern%n_samples, kern%x, kern%r_min, kern%dr, kern%sigma, PDF_KDE_CUT, lo, hi, row)
            if (hi < lo) cycle
            peak = maxval(row(lo:hi))
            if (peak <= 0.d0) cycle
            do l = lo, hi
               exact = exp(-((real(kern%x(l), qp) - real(r, qp))/real(kern%sigma, qp))**2/2)
               worst_abs = max(worst_abs, real(abs(row(l) - exact), dp)/peak)
               if (row(l) > 1.d-6*peak) worst_rel = max(worst_rel, real(abs(row(l) - exact)/exact, dp))
            end do
         end do
         write (*, '(A,I4,A,F5.2,A,ES9.2,A,ES9.2,A,ES9.2,A)') '  samples ', n_samples(g), ', sigma ', sigmas(s), &
            ': row from quad', worst_abs, ' of the peak,', worst_rel, ' relative; grid', dx, ' A'
         if (worst_abs > RTOL .or. worst_rel > RTOL) then
            write (*, '(A)') 'FAIL: a row is further from the quad Gaussian than RTOL'
            n_fail = n_fail + 1
         end if
         if (dx > DX_TOL) then
            write (*, '(A)') 'FAIL: the stored grid is further than DX_TOL from evenly spaced'
            n_fail = n_fail + 1
         end if
         deallocate (row)
      end do
   end do

   if (n_fail > 0) then
      write (*, '(I0,A)') n_fail, ' failure(s)'
      stop 1
   end if
   write (*, '(A)') 'PASS'
end program pdf_kernelverify
