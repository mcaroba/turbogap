! Drive ir_fft.f90 from a dipole file, so that the Fortran pipeline can be
! compared line by line with the Python it was translated from
! (TNEP/spectroscopy.py). Also finite-differences the MAD adjoint.
!
! Usage:
!
!   irfftverify spectrum <dipoles.dat> <dt_fs> <window> <acf_ratio> \
!               <max_freq_cm> <smooth_k> <smooth_kind> <qc> <T> <dc_cut> \
!               <out.dat>
!
!       <dipoles.dat> is T lines of "mux muy muz", chronological. The output
!       has one line per bin: nu, intensity(normalised), power(normalised),
!       intensity_raw, power_raw, and a header carrying the sizing.
!
!   irfftverify acf <dipoles.dat> <out.dat>
!
!       C(tau) for every lag, so that the FFT correlation can be checked
!       against a direct sum without the rest of the pipeline in the way.
!
!   irfftverify grad <dipoles.dat> <dt_fs> <exp.dat> <h> [qc] [smooth_kind]
!
!       h-scan of the MAD gradient: compares lambda from ir_fft_loss with a
!       central difference of the loss in each Cartesian component of the
!       NEWEST dipole. An h-scan and not a single h, because a gradient that
!       is wrong by a constant factor passes a single-h check whenever the
!       tolerance is loose enough, and does not pass a scan.
!
!   irfftverify fft <n>
!
!       Round-trip and DFT-identity checks on the radix-2 transform itself.
!
program irfftverify

   use kinds
   use ir_fft

   implicit none

   character(len=1024) :: mode, arg, fname, outname
   integer :: nargs

   nargs = command_argument_count()
   if (nargs < 1) then
      write (*, *) "usage: irfftverify spectrum|acf|grad|fft ..."
      stop 1
   end if
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ("spectrum")
      call do_spectrum()
   case ("acf")
      call do_acf()
   case ("grad")
      call do_grad()
   case ("fft")
      call get_command_argument(2, arg)
      call do_fft(arg)
   case default
      write (*, *) "unknown mode ", trim(mode)
      stop 1
   end select

contains

!  Read a T x 3 text file of dipoles.
   subroutine read_dipoles(fname, mu, n)

      implicit none

      character(len=*), intent(in) :: fname
      real(dp), allocatable, intent(out) :: mu(:, :)
      integer, intent(out) :: n
      integer :: u, ios, i
      real(dp) :: a, b, c

      open (newunit=u, file=trim(fname), status="old", action="read")
      n = 0
      do
         read (u, *, iostat=ios) a, b, c
         if (ios /= 0) exit
         n = n + 1
      end do
      rewind (u)
      allocate (mu(1:3, 1:n))
      do i = 1, n
         read (u, *) mu(1, i), mu(2, i), mu(3, i)
      end do
      close (u)

   end subroutine read_dipoles

!  Read a two-column nu, I file.
   subroutine read_exp(fname, nu, I_e, n)

      implicit none

      character(len=*), intent(in) :: fname
      real(dp), allocatable, intent(out) :: nu(:), I_e(:)
      integer, intent(out) :: n
      integer :: u, ios, i
      real(dp) :: a, b

      open (newunit=u, file=trim(fname), status="old", action="read")
      n = 0
      do
         read (u, *, iostat=ios) a, b
         if (ios /= 0) exit
         n = n + 1
      end do
      rewind (u)
      allocate (nu(1:n), I_e(1:n))
      do i = 1, n
         read (u, *) nu(i), I_e(i)
      end do
      close (u)

   end subroutine read_exp

   subroutine do_spectrum()

      implicit none

      type(ir_fft_config_type) :: cfg
      type(ir_fft_result_type) :: res
      real(dp), allocatable :: mu(:, :)
      integer :: n, u, k
      logical :: ok
      character(len=512) :: msg

      call get_command_argument(2, fname)
      call read_dipoles(fname, mu, n)
      call get_command_argument(3, arg); read (arg, *) cfg%dt_fs
      call get_command_argument(4, cfg%window)
      call get_command_argument(5, arg); read (arg, *) cfg%acf_ratio
      call get_command_argument(6, arg); read (arg, *) cfg%max_freq_cm
      call get_command_argument(7, arg); read (arg, *) cfg%smooth_k
      call get_command_argument(8, cfg%smooth_kind)
      call get_command_argument(9, cfg%quantum_correction)
      call get_command_argument(10, arg); read (arg, *) cfg%temperature
      call get_command_argument(11, arg); read (arg, *) cfg%power_dc_cutoff_cm
      call get_command_argument(12, outname)
      cfg%normalise = .true.
      cfg%subtract_mean = .true.

      call ir_fft_spectrum(mu, n, cfg, res, ok, msg)
      if (.not. ok) then
         write (*, *) "ERROR: ", trim(msg)
         stop 1
      end if

      open (newunit=u, file=trim(outname), status="replace", action="write")
      write (u, '(A,I0,A,I0,A,I0)') "# n_frames = ", res%n_frames, &
         "   n_lag = ", res%n_lag, "   n_freq = ", res%n_freq
      write (u, '(A,F14.6,A,F14.6,A,F14.6)') "# d_nu = ", res%d_nu, &
         "   resolution = ", res%resolution, "   nyquist = ", res%nyquist
      write (u, '(A,ES22.14,A,ES22.14)') "# peak = ", res%peak, &
         "   peak_power = ", res%peak_power
      write (u, '(A)') "# nu   intensity   power   intensity_raw   power_raw"
      do k = 1, res%n_freq
         write (u, '(5ES24.15)') res%freq(k), res%intensity(k), res%power(k), &
            res%intensity_raw(k), res%power_raw(k)
      end do
      close (u)

      write (*, '(A,I0,A,I0,A,F10.4)') "n_lag = ", res%n_lag, "  n_freq = ", &
         res%n_freq, "  d_nu = ", res%d_nu

   end subroutine do_spectrum

   subroutine do_acf()

      implicit none

      real(dp), allocatable :: mu(:, :), acf(:), direct(:)
      real(dp) :: mu_mean(1:3), acc, worst, den
      integer :: n, u, t, a, d
      integer :: n_direct

      call get_command_argument(2, fname)
      call get_command_argument(3, outname)
      call read_dipoles(fname, mu, n)

      allocate (acf(0:n - 1))
      call ir_fft_autocorrelation(mu, n, .true., acf, mu_mean)

!     The same thing by the definition, so that the FFT is checked against
!     what it is supposed to be computing rather than only against the Python
!     (which uses an FFT too, so an error in the padding would agree with
!     itself). Capped at 400 lags: this is O(n * n_direct).
      n_direct = min(n - 1, 400)
      allocate (direct(0:n_direct))
      worst = 0.d0
      do t = 0, n_direct
         acc = 0.d0
         do a = 1, n - t
            do d = 1, 3
               acc = acc + (mu(d, a) - mu_mean(d))*(mu(d, a + t) - mu_mean(d))
            end do
         end do
         direct(t) = acc/dfloat(n)
         den = max(dabs(direct(t)), dabs(acf(0)))
         if (den > 0.d0) worst = max(worst, dabs(direct(t) - acf(t))/den)
      end do

      open (newunit=u, file=trim(outname), status="replace", action="write")
      write (u, '(A)') "# tau   acf_fft   acf_direct"
      do t = 0, n_direct
         write (u, '(I8,2ES24.15)') t, acf(t), direct(t)
      end do
      do t = n_direct + 1, n - 1
         write (u, '(I8,ES24.15)') t, acf(t)
      end do
      close (u)

      write (*, '(A,ES12.4)') "worst relative FFT-vs-direct difference: ", worst
      if (worst > 1.d-11) then
         write (*, *) "FAIL: the FFT correlation does not match the definition"
         stop 1
      end if
      write (*, *) "PASS: FFT correlation matches the direct sum"

   end subroutine do_acf

   subroutine do_grad()

      implicit none

      type(ir_fft_config_type) :: cfg
      real(dp), allocatable :: mu(:, :), nu(:), I_e(:), wgt(:), I_fit(:)
      real(dp) :: h, energy, lambda(1:3), scale, offset, dissim, dissim_ref
      real(dp) :: ep, em, fd(1:3), rel, worst
      real(dp) :: save_mu
      integer :: n, n_exp, i, c, j
      logical :: ok
      character(len=512) :: msg
      character(len=1024) :: qc, sk, sm

      call get_command_argument(2, fname)
      call get_command_argument(3, arg); read (arg, *) cfg%dt_fs
      call get_command_argument(4, outname)
      call get_command_argument(5, arg); read (arg, *) h
      qc = "harmonic"
      sk = "gaussian"
      sm = "mean"
      if (command_argument_count() >= 6) call get_command_argument(6, qc)
      if (command_argument_count() >= 7) call get_command_argument(7, sk)
      if (command_argument_count() >= 8) call get_command_argument(8, sm)

      call read_dipoles(fname, mu, n)
      call read_exp(outname, nu, I_e, n_exp)
      allocate (wgt(1:n_exp), I_fit(1:n_exp))
      wgt = 1.d0

      cfg%window = "hann"
      cfg%acf_ratio = 0.2d0
      cfg%max_freq_cm = 4000.d0
      cfg%smooth_k = 5
      cfg%smooth_kind = trim(sk)
      cfg%quantum_correction = trim(qc)
      cfg%temperature = 300.d0
      cfg%subtract_mean = (trim(sm) == "mean")
      cfg%normalise = .false.

      call ir_fft_loss(mu, n, cfg, nu, I_e, wgt, n_exp, .true., .false., 1.d0, &
                       energy, lambda, I_fit, scale, offset, dissim, dissim_ref, &
                       ok, msg)
      if (.not. ok) then
         write (*, *) "ERROR: ", trim(msg)
         stop 1
      end if

      write (*, '(A,ES16.8,A,ES14.6,A,ES14.6)') "energy = ", energy, &
         "   scale = ", scale, "   dissim = ", dissim
      write (*, '(A,3ES16.8)') "lambda   = ", lambda(1:3)
      write (*, '(A)') ""
      write (*, '(A)') "     h        c      analytic          finite diff        rel.err"

      do j = 0, 7
         worst = 0.d0
         do c = 1, 3
            save_mu = mu(c, n)
            mu(c, n) = save_mu + h
            call ir_fft_loss(mu, n, cfg, nu, I_e, wgt, n_exp, .true., .false., 1.d0, &
                             ep, fd, I_fit, scale, offset, dissim, dissim_ref, ok, msg)
            mu(c, n) = save_mu - h
            call ir_fft_loss(mu, n, cfg, nu, I_e, wgt, n_exp, .true., .false., 1.d0, &
                             em, fd, I_fit, scale, offset, dissim, dissim_ref, ok, msg)
            mu(c, n) = save_mu
            fd(c) = (ep - em)/(2.d0*h)
            rel = dabs(fd(c) - lambda(c))/max(dabs(lambda(c)), 1.d-30)
            worst = max(worst, rel)
            write (*, '(ES10.2,I8,3ES18.8)') h, c, lambda(c), fd(c), rel
         end do
         write (*, '(A,ES10.2,A,ES12.4)') "  -> h = ", h, "   worst rel.err = ", worst
         h = h*0.5d0
      end do

      i = 0

   end subroutine do_grad

!  Check the transform on its own: a forward followed by an inverse must be
!  the identity up to 1/n, and the forward must agree with the DFT written out
!  as a sum.
   subroutine do_fft(arg_n)

      implicit none

      character(len=*), intent(in) :: arg_n
      complex(dp), allocatable :: z(:), z0(:), zd(:)
      real(dp) :: two_pi, th, worst1, worst2, s
      integer :: n, i, k, j
      integer :: seed_size
      integer, allocatable :: seed(:)
      real(dp) :: r1, r2

      read (arg_n, *) n
      n = ir_fft_next_pow2(n)
      allocate (z(0:n - 1), z0(0:n - 1), zd(0:n - 1))

      call random_seed(size=seed_size)
      allocate (seed(seed_size))
      seed = 20260824
      call random_seed(put=seed)
      do i = 0, n - 1
         call random_number(r1)
         call random_number(r2)
         z0(i) = dcmplx(2.d0*r1 - 1.d0, 2.d0*r2 - 1.d0)
      end do

!     round trip
      z = z0
      call ir_fft_transform(z, n, -1)
      call ir_fft_transform(z, n, +1)
      z = z/dfloat(n)
      worst1 = 0.d0
      do i = 0, n - 1
         worst1 = max(worst1, abs(z(i) - z0(i)))
      end do

!     against the definition, on the first min(n,64) bins
      two_pi = 2.d0*dacos(-1.d0)
      z = z0
      call ir_fft_transform(z, n, -1)
      worst2 = 0.d0
      s = 0.d0
      do k = 0, min(n - 1, 63)
         zd(k) = dcmplx(0.d0, 0.d0)
         do j = 0, n - 1
            th = -two_pi*dfloat(j)*dfloat(k)/dfloat(n)
            zd(k) = zd(k) + z0(j)*dcmplx(dcos(th), dsin(th))
         end do
         worst2 = max(worst2, abs(zd(k) - z(k)))
         s = max(s, abs(zd(k)))
      end do

      write (*, '(A,I0)') "n = ", n
      write (*, '(A,ES12.4)') "round-trip max |err|          : ", worst1
      write (*, '(A,ES12.4)') "forward vs definition max |err|: ", worst2/s
      if (worst1 > 1.d-12 .or. worst2/s > 1.d-12) then
         write (*, *) "FAIL"
         stop 1
      end if
      write (*, *) "PASS"

   end subroutine do_fft

end program irfftverify
