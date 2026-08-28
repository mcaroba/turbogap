
! Tigany Zarrouk 26/08/2026
!
! This module is meant to implement the use of extended variables to allow for
! the matching of IR spectra primarily and maybe the inclusion of other
! temporal/frequency-related spectra
!
! The general idea is that we couple the observable directly to a set of fictitious variables.
!
! If we start simply from the case of one set of fictitious variables {x_k}
!
! We can make a modified Hamiltonian
!
! The initial Hamiltonian is
!
! H = T + V
!
! where V = V( {r^a_k} ) and T = T( {r^a_k} ), where r^a_k is an atomic component.
!
! In the initial MAD formulation we have
!
! H = T + V + V'
!
! Where V' is the experimental potential which is given by
!
! V' = V'({r^a_k}; \nu) = 1/2 \sum_\nu ( I_pred(\nu) - I_exp(\nu) )^2
!
! Where \nu is the dependent variable of the observable (which, is frequency in the case of IR spectra).
!
! In the case of an observable that depends on the atomic positions, then we can
! construct a concrete conservative potential which can minimize the loss.
!
! L = V + V'
!
! However, in the case of IR spectra, we have an observable that depends on both
! the atomic coordinates and an ensemble of structures, as the IR spectrum is an
! absorption spectrum which samples the characteristic frequencies of the modes
! present in the experimentally measured sample.
!
! IR: I_pred = I_pred( \nu, { r^a_k } )
!     I_pred = prefactor * FT( autocorrelation_function( m(0), m(t)  ))
!
! Where the prefactor depends on the quantum correction and m(t) is the total dipole predicted at a given time
! In the case of the models developed so far,
!
!
! In the lagrangian we have the normal lagrangian as
! L = T - V
!
! Doing something similar to the fictitious variables, we have
!
! T_ext_i = 1/2 eff_mass_k d_t(X_k_i)^2
! V_ext_i = 1/2 w_k^2 X_k_i^2
! L_ext =  \sum_{i,k} [  P_k_i^2/(2 eff_mass_k) - 1/2 eff_mass_k w_k^2 X_k_i^2 + g_k X_k_i . M_i(q_soap) ]
!
! So this L_ext, which is to be added to the original lagrangian, is a set of harmonic oscillators which are /driven/ by the driving force g_k X_k_i m_i(q_soap)
!
! where mu_k is the fictitious mass associated with the extended variable X_k_i, and where we have k enumerating the oscillators.
!
! Define
!
! R_k^2 = X_k_i^2 + (d( X_k_i )/dt)^2 / w_k^2
!
! which where s_k -> x_k  and d( s_k )/dt -> x_k_dot (the momentum)
!
! So R_k is the amplitude of these oscillators. This encodes the information that we need to read out the spectrum later.
!
! The total Hamiltonian is
!
! H = T + V + \sum_k [  P_k^2/(2 eff_mass_k) + 1/2 eff_mass_k w_k^2 X_k^2 + g_k X_k . M(q_soap) ] + V'
!
! From the Hamiltonian relations we have
!
! For the standard physical coordinates
! dr_i/dt =  v
! dv_i/dt =  F/m + \sum_k g_k grad( M(q_soap) ) . X_k / m
!
! dX_k/dt = P_k / eff_mass_k + dV'/dP_k
! dP_k/dt = -eff_mass_k w_k^2 X_k  + g_k M(q_soap) - dV'/dX_k
!
!
! dV'/dX_k = dV'/dr^a_l * X_k/R_k
! dV'/dP_k = dV'/dr^a_l * P_k/(eff_mass_k^2 w_k^2 R_k)
!
! The Hamiltonian is conserved over time.
!
! During dynamics, the extended variables will couple to the atoms which will change their velocities
!
! And now the bias potential has the explicit form
!
! V' = 1/2 ( R_k - R_k_exp )^2
!
! Where R_k_exp is tuned to the current scale that turbogap estimates the dipole autocorrelation function.
!

! To see why this works, think of Carr-Parrinello molecular dynamics, where a
! set of variables is introduced to propagate the wavefunctions by the auxiliary
! variables.
!
! This has a lagrangian of a similar form, just the "forcing" term is the
! Lagrangian multiplier which restricts the fact that the total wavefunction
! should be normalized. See https://www.acmm.nl/ensing/thesis/node10.html

module ir_auxiliary_dynamics
   use kinds, only: dp
   ! use the calculation of the dipole gradients (the central atom jacobian of
   ! soap), as long as it is turned on
   use read_utils, only: check_file_exists
   use mad_ir, only: mad_ir_dmu_dr
   use error, only: turbogap_abort
   use md, only: randomize_velocities, remove_cm_vel
   implicit none

   private
   public :: calc_amplitude_bias_forces

   !  Wavenumber in cm^-1 of a frequency of 1/fs: 1/(c) with c in cm/fs.
   real(dp), parameter :: CM_PER_INV_FS = 33356.40952d0
   real(dp), parameter :: AMU = 103.6426965268  ! amu -> eV fs^2 / A^2
   real(dp), parameter :: DEFAULT_EFF_MASS = 100.0_dp*AMU ! amu
   real(dp), parameter :: kB = 8.6173303d-5
   real(dp), parameter :: PI = dacos(-1.0_dp)
   real(dp), parameter :: R_MAX = 3.0_dp

   type ir_auxiliary

      ! This is the variable which sets the size of the extended variables
      integer :: n_freq = -1

      ! Extended variables
      ! size is 3 * n_modes
      ! Can think of X and P being like oscillators that react to an either
      ! electric or a magnetic field.
      real(dp), allocatable :: X(:, :)
      real(dp), allocatable :: X_prev(:, :)
      real(dp), allocatable :: P(:, :)
      real(dp), allocatable :: P_pred(:, :)
      real(dp), allocatable :: E_aux(:)

      ! Magnitude of the extended variables gives you R, which is like an energy
      ! density of sorts
      real(dp), allocatable :: R(:)

      real(dp), allocatable :: R_target(:)
      real(dp), allocatable :: normalisation

      ! Effective masses
      ! These should be larger than the masses of the atoms
      ! Therefore ~100amu, such that the auxiliary variables can still track the
      real(dp), allocatable :: eff_mass(:)

      ! Frequencies
      real(dp), allocatable :: omega(:)

      ! The couplings which are the size of the frequencies
      real(dp), allocatable :: g_k(:)

      ! The couplings can be calibrated by a linear response which just needs a
      ! calculation of the bare FT( dipole_autocorrelation )
      ! The bandwidth ( resolution ), is the damping term
      ! damping_k === resolution ~= 5-10cm^-1
      !
      ! C(omega) = \int < M(0) . M(t) > exp( - i omega t)
      !
      ! Susceptibility is
      !
      ! \chi_k ( omega ) = g_k / eff_mass_k / ( omega_k^2 - omega^2  - i omega damping_k )
      !
      ! < R_k^2 > = 1/2\pi \int d_omega |\chi_k(omega)|^2 ( 1 + omega^2/omega_k^2 ) C(omega)
      !
      ! Assuming the resolution is quite sharp, then \chi is essentially a dirac delta giving
      !
      ! < R_k^2 > ~=  g_k^2 / (eff_mass_k^2  omega_k^2 damping_k ) C(omega_k)
      !
      ! So the couplings can be determined by
      !
      ! g_k ~= eff_mass_k omega_k sqrt( damping_k < R_k^2, target > / C(omega_k))
      !

      ! Just setting even damping right now
      ! It should be set to the resolution of experiment, but this is a rough value.
      real(dp), allocatable :: damping_k(:)
      real(dp) :: damping_default = 10.0_dp*CM_PER_INV_FS

      ! Forces on the extended variables
      real(dp), allocatable :: f_aux(:, :)
      real(dp), allocatable :: f_aux_prev(:, :)

      ! Thermo parameters
      real(dp) :: t_beg_aux
      real(dp) :: instant_temp_aux
      real(dp) :: E_kinetic_aux

   end type ir_auxiliary

contains

   pure elemental function gaussian(x, mean, sigma) result(res)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: mean
      real(dp), intent(in) :: sigma
      real(dp) :: res
      res = 1/dsqrt(2.0_dp*PI)*exp(-(x - sigma)**2/2.0_dp/sigma**2)
   end function gaussian

   subroutine remove_center_of_mass_momentum(vel, M)
      ! Here vel is the actually the momentum, this was copied from md.f90
      real(dp), intent(inout) :: vel(:, :)
      real(dp), intent(in) :: M(:)
      real(dp) :: cm_pos(1:3)
      real(dp) :: cm_vel(1:3)
      real(dp) :: total_mass
      integer :: Np
      integer :: i

      Np = size(vel, 2)

      cm_vel = 0.d0
      total_mass = 0.d0
      do i = 1, Np
         ! removing the mass from the line as it's already a momentum
         ! cm_vel(1:3) = cm_vel(1:3) + M(i)*vel(1:3, i)
         cm_vel(1:3) = cm_vel(1:3) + vel(1:3, i)
         total_mass = total_mass + M(i)
      end do
      cm_vel = cm_vel/total_mass
      do i = 1, Np
         vel(1:3, i) = vel(1:3, i) - M(i)*cm_vel(1:3)
      end do

   end subroutine remove_center_of_mass_momentum

   subroutine get_amplitude(this)
      type(ir_auxiliary), intent(inout) :: this
      real(dp) :: X_mag_sq
      real(dp) :: P_mag_sq
      integer :: k
      do k = 1, this%n_freq
         X_mag_sq = this%X(1, k)**2 + this%X(2, k)**2 + this%X(3, k)**2
         P_mag_sq = this%P(1, k)**2 + this%P(2, k)**2 + this%P(3, k)**2
         this%R(k) = dsqrt(X_mag_sq + P_mag_sq/(this%eff_mass(k)**2*this%omega(k)**2))
      end do
   end subroutine get_amplitude

   subroutine ir_auxiliary_restart_write(this, ir_auxiliary_restart_file_name)
      type(ir_auxiliary), intent(inout) :: this
      character*1024, intent(in) :: ir_auxiliary_restart_file_name
      integer :: unit_number
      integer :: k

      ! Format is
      ! n_freq
      ! X11 X12 Z13 P11 P12 P13 eff_mass omega R_exp g_k
      open (newunit=unit_number, file=ir_auxiliary_restart_file, status="replace")
      write (unit_number, '(I10)') this%n_freq
      do k = 1, this%n_freq
         write (unit_number, '(10(F20.8,1X))') &
            this%X(1, k), this%X(2, k), this%X(3, k), &
            this%P(1, k), this%P(2, k), this%P(3, k), &
            this%eff_mass(k), &
            this%omega(k), &
            this%R_target(k), &
            this%g_k(k)
      end do
      close (unit_number)

   end subroutine ir_auxiliary_restart_write

   subroutine ir_auxiliary_restart_read(this, ir_auxiliary_restart_file)
      type(ir_auxiliary), intent(inout) :: this
      character*1024, intent(in) :: ir_auxiliary_restart_file
      integer :: unit_number
      integer :: iostatus
      integer :: n_lines
      integer :: n_freq
      integer :: k

      if (.not. trim(ir_auxiliary_restart_file) == "none") then
         ! We read the restart file and allocate the IR things
         call check_file_exists(ir_auxiliary_restart_file)
         ! As the file must exist, we can read it in
         ! First check the number of lines

         open (newunit=unit_number, file=ir_auxiliary_restart_file, status="old")
         iostatus = 0
         n_lines = -1
         do while (iostatus == 0)
            read (unit_number, *, iostat=iostatus)
            n_lines = n_lines + 1
         end do
         close (unit_number)

         ! Format is
         ! n_freq
         ! X11 X12 Z13 P11 P12 P13 eff_mass omega R_exp g_k

         if (n_lines < 2) then
            write (*, *) 'ERROR: ir_auxiliary_variable: '
            write (*, *) 'The ir auxiliary restart file is invalid! There is no data here'
            call turbogap_abort()
         end if

         open (newunit=unit_number, file=ir_auxiliary_restart_file, status="old")
         iostatus = 0
         read (unit_number, *, iostat=iostatus) n_freq

         if (iostatus /= 0) then
            write (*, *) 'ERROR: ir_auxiliary_variable: There is something wrong with the'
            write (*, *) '       ir_auxiliary_restart_file supplied, we cannot continue!'
            close (unit_number)
            call turbogap_abort()
         end if

         if (n_freq < 1) then
            write (*, *) 'ERROR: ir_auxiliary_variable: There are no frequencies'
            write (*, *) '       or the file is garbled! Make sure the start'
            write (*, *) '       of the file has n_freq equal to something greater than zero.'
            close (unit_number)
            call turbogap_abort()
         end if

         if (n_freq > n_lines - 1) then
            write (*, *) 'ERROR: ir_auxiliary_variable: The number of frequencies in the'
            write (*, *) '       ir_auxiliary_restart_file is more than the number of lines'
            write (*, *) '       we cannot continue with this!'
            close (unit_number)
            call turbogap_abort()
         end if

         if (n_freq < n_lines - 1) then
            write (*, *) 'WARNING: ir_auxiliary_variable: The number of frequencies in the'
            write (*, *) '         ir_auxiliary_restart_file is less than the number of lines'
            write (*, *) '         this is a little weird!'
            write (*, *) '         but we will assume you know what you are doing!'
         end if

         ! n_freq is okay now, therefore allocate and try to read the rest
         this%n_freq = n_freq

         allocate (this%X(3, n_freq))
         allocate (this%X_prev(3, n_freq))
         allocate (this%P(3, n_freq))
         allocate (this%P_pred(3, n_freq))
         allocate (this%f_aux(3, n_freq))
         allocate (this%f_aux_prev(3, n_freq))
         allocate (this%E_aux(n_freq))
         allocate (this%R(n_freq))
         allocate (this%R_target(n_freq))
         allocate (this%omega(n_freq))
         allocate (this%g_k(n_freq))
         allocate (this%damping_k(n_freq))

         do k = 1, n_freq

            if (iostatus /= 0) then
               write (*, *) 'ERROR: ir_auxiliary_variable: There is something wrong with the'
               write (*, *) '       ir_auxiliary_restart_file while reading, we cannot continue!'
               call turbogap_abort()
            end if
            ! Check the number of variables in the line that it matches what we expect
            ! Not doing that yet just trying out
            ! X11 X12 Z13 P11 P12 P13 eff_mass omega R_exp g_k
            read (unit_number, *, iostat=iostatus) &
               this%X(1, k), this%X(2, k), this%X(3, k), &
               this%P(1, k), this%P(2, k), this%P(3, k), &
               this%eff_mass(i), &
               this%omega(k), &
               this%R_target(k), &
               this%g_k(k)
         end do

         close (unit_number)
      end if
   end subroutine ir_auxiliary_restart_read

   ! this routine should be done after an initial ACF calculation
   subroutine ir_auxiliary_dynamics_init(this, n_freq, freqs, ir_exp, dipole_acf, t_beg_aux, eff_mass, damping_k)
      type(ir_auxiliary), intent(inout) :: this
      integer, intent(in) :: n_freq
      ! The initial auxiliary temperature, which is usually the physical temperature
      real(dp), intent(in) :: t_beg_aux
      ! The freqs are the input freqs of the spectrum we want to match
      real(dp), intent(in) :: freqs(n_freq)
      ! Experimental intensities
      real(dp), intent(in) :: ir_exp(n_freq)
      ! autocorrelation function if doing some run to determine the correct parameters
      ! The dipole_acf which is input is the bare, un-normalized one
      real(dp), intent(in) :: dipole_acf(n_freq)

      ! Input effective mass for all the oscillators
      real(dp), intent(in) :: eff_mass
      ! Input damping for all the oscillators
      real(dp), intent(in) :: damping_k
      real(dp) :: max_ir_exp

      real(dp) :: rk_sq
      real(dp) :: inv_eff_mass_omega_sq

      integer :: i, k

      if (n_sites < 1) then
         write (*, *) 'ERROR: ir_auxiliary_dynamics: n_sites is less than one!! Exiting.'
         call turbogap_abort()
      end if

      if (n_freq < 1) then
         write (*, *) 'ERROR: ir_auxiliary_dynamics: n_freq is less than zero!! Exiting.'
         call turbogap_abort()
      end if

      if (t_beg_aux < 0.0_dp) then
         write (*, *) 'ERROR: ir_auxiliary_dynamics: t_beg_aux is less than zero!! Exiting.'
         call turbogap_abort()
      end if

      this%n_freq = n_freq
      this%t_beg_aux = t_beg_aux

      ! Allocate the arrays
      allocate (this%X(3, this%n_freq))
      allocate (this%X_prev(3, this%n_freq))
      allocate (this%P(3, this%n_freq))
      allocate (this%P_pred(3, this%n_freq))

      allocate (this%f_aux(3, this%n_freq))
      allocate (this%f_aux_prev(3, this%n_freq))

      allocate (this%E_aux(this%n_freq))

      ! The magnitude/amplitude of the oscillators
      allocate (this%R(this%n_freq))
      ! The magnitude of the experimental observable at omega_k
      allocate (this%R_target(this%n_freq))

      allocate (this%omega(this%n_freq))
      this%omega = freqs

      ! Set later
      allocate (this%g_k(this%n_freq))

      allocate (this%damping_k(this%n_freq))

      if (damping_k > 1e-12) then
         this%damping_k = damping_k
      else
         this%damping_k = this%damping_default
      end if

      ! initialise the effective masses to the default
      if (eff_mass > 1e-12) then
         allocate (this%eff_mass(this%n_freq), source=eff_mass*AMU)
      else
         allocate (this%eff_mass(this%n_freq), source=DEFAULT_EFF_MASS)
      end if

      ! This takes in velocities, usually and not the momentum, so we will just multiply by the effective masses after
      ! Make the distribution Maxwellian
      call randomize_velocities( &
         this%P, &
         this%n_freq, &
         this%E_kinetic_aux, &
         this%eff_mass, &
         this%instant_temp_aux, &
         this%t_beg_aux, &
         "maxwell")

      do k = 1, this%n_freq
         this%P(1:3, k) = this%P(1:3, k)*eff_mass(k)
      end do

      ! First randomize X and P to be uniform then get distributions: Maxwellian
      call random_number(this%X)
      call random_number(this%P)
      do k = 1, this%n_freq
         this%X(1:3, k) = dsqrt(kB*this%t_beg_aux/(this%eff_mass(k)*this%omega(k)**2))*gaussian(this%X(1:3, k), 0.0_dp, 1.0_dp)
         this%P(1:3, k) = dsqrt(kB*this%t_beg_aux*this%eff_mass(k)*this%omega(k)**2)*gaussian(this%P(1:3, k), 0.0_dp, 1.0_dp)
      end do

      ! call remove_cm_vel(this%P(1:3, 1:n_sites), masses(1:n_sites))

      call remove_center_of_mass_momentum(this%P(1:3, 1:n_sites), masses(1:n_sites))

      ! Get the value of R here, from the magnitude
      call get_amplitude(this)

      ! Set the normalization factor to be relative to the maximum intensity
      ! From here we assume a classical prefactor of omega^2, but we should change this.
      max_ir_exp = -1e6
      do i = 1, this%n_freq
         max_ir_exp = max(max_ir_exp, ir_exp(k)/this%omega(k)**2)
      end do

      this%normalisation = R_MAX**2/max_ir_exp

      ! Now get the < R_k^2_target >
      do k = 1, this%n_freq
         rk_sq = this%normalisation*ir_exp(k)/this%omega(k)**2
         this%R_target(k) = dsqrt(rk_sq)

         ! g_k ~= eff_mass_k omega_k sqrt( damping_k < R_k^2, target > / C(omega_k))
         ! This tunes the coupling to something relative to the scale of the
         ! code in comparison to the experimental spectrum
         this%g_k(k) = this%eff_mass(k)*this%omega(k) &
                       *dsqrt(this%damping_k(k)*rk_sq/dipole_acf(k))
      end do

   end subroutine ir_auxiliary_dynamics_init

   subroutine ir_auxiliary_dynamics_free(this)
      type(ir_auxiliary) :: this

      if (allocated(this%X)) deallocate (this%X)
      if (allocated(this%X_prev)) deallocate (this%X_prev)
      if (allocated(this%P)) deallocate (this%P)
      if (allocated(this%P_pred)) deallocate (this%P_pred)

      if (allocated(this%E_aux)) deallocate (this%E_aux)
      if (allocated(this%f_aux)) deallocate (this%f_aux)
      if (allocated(this%f_aux_prev)) deallocate (this%f_aux_prev)

      if (allocated(this%R)) deallocate (this%R)
      if (allocated(this%R_target)) deallocate (this%R_target)

      if (allocated(this%omega)) deallocate (this%omega)
      if (allocated(this%g_k)) deallocate (this%g_k)
      if (allocated(this%damping_k)) deallocate (this%damping_k)

   end subroutine ir_auxiliary_dynamics_free

   ! The dipole jacobian here is exactly the d mu_a / d r_jb that comes from accumulate_dmu_dr
   subroutine calculate_ir_auxiliary_physical_force_bias( &
      n_sites, n_freq, g_k, x_k, dipole_jacobian, f_phys)
      integer, intent(in)  :: n_sites
      integer, intent(in)  :: n_freq
      real(dp), intent(in) :: g_k(n_freq)
      real(dp), intent(in) :: x_k(3, n_freq)

      ! From accumulate dmu_dr, the jacobian:
      !   dmu_dr(a, b, j) = d mu_a / d r_jb,     mu = sum_i mu_i
      !
      ! The dipole jacobian tensor is : (dipole_component, coordinate_component, atom_index)
      ! J_{alpha, beta} = d m_{alpha} / d q_{beta}
      real(dp), intent(in)    :: dipole_jacobian(3, 3, n_sites)

      ! Physical forces to be updated, these will be mpi reduced
      real(dp), intent(inout) :: f_phys(3, n_sites)

      integer  :: i, k
      real(dp) :: E_eff(3)

      E_eff = 0.0_dp

      !$omp parallel do reduction(+:E_eff) private(k) schedule(static)
      do k = 1, n_freq
         E_eff(1) = E_eff(1) + g_k(k)*x_k(1, k)
         E_eff(2) = E_eff(2) + g_k(k)*x_k(2, k)
         E_eff(3) = E_eff(3) + g_k(k)*x_k(3, k)
      end do
      !$omp end parallel do

      !$omp parallel do private(i) schedule(static)
      do i = 1, n_sites

         ! Unrolled matrix-vector multiplication: F_i = J_i^T * E_eff
         ! We transpose by iterating over the first index (alpha) of the Jacobian
         !
         ! here we have 3x3 blocks for each site i
         ! (1,1) (1,2) (1,3)
         ! (2,1) (2,2) (2,3)
         ! (3,1) (3,2) (3,3)
         !
         ! Transpose is
         ! (1,1) (2,1) (3,1)
         ! (1,2) (2,2) (3,2)
         ! (1,3) (2,3) (3,3)
         !
         ! First mult is
         ! (1,1) * E_eff(1) (2,1) * E_eff(2) (3,1) * E_eff(3)

         f_phys(1, i) = f_phys(1, i) + &
                        dipole_jacobian(1, 1, i)*E_eff(1) + &
                        dipole_jacobian(2, 1, i)*E_eff(2) + &
                        dipole_jacobian(3, 1, i)*E_eff(3)

         f_phys(2, i) = f_phys(2, i) + &
                        dipole_jacobian(1, 2, i)*E_eff(1) + &
                        dipole_jacobian(2, 2, i)*E_eff(2) + &
                        dipole_jacobian(3, 2, i)*E_eff(3)

         f_phys(3, i) = f_phys(3, i) + &
                        dipole_jacobian(1, 3, i)*E_eff(1) + &
                        dipole_jacobian(2, 3, i)*E_eff(2) + &
                        dipole_jacobian(3, 3, i)*E_eff(3)

      end do
      !$omp end parallel do

   end subroutine calculate_ir_auxiliary_physical_force_bias

   ! Calculates the bias forces on the auxiliary spatial (X) and momentum (P) coordinates.
   subroutine calculate_ir_auxiliary_extended_variable_force_bias( &
      n_freq, x_k, p_k, g_k, dipole, eff_mass, omega, r_target, k_bias, E_aux, f_x, f_p)
      ! Number of modes
      integer, intent(in)  :: n_freq
      ! Extended variables
      real(dp), intent(in) :: x_k(3, n_freq)
      real(dp), intent(in) :: p_k(3, n_freq)
      real(dp), intent(in) :: g_k(n_freq)
      real(dp), intent(in) :: E_aux(n_freq)
      real(dp), intent(in) :: dipole(3)
      ! Effective masses
      real(dp), intent(in) :: eff_mass(n_freq)
      ! Frequencies
      real(dp), intent(in) :: omega(n_freq)
      ! Experimental observable
      real(dp), intent(in) :: r_target(n_freq)
      ! The biases ( exp_energy_scale * weights in other MAD stuff)
      real(dp), intent(in) :: k_bias(n_freq)

      real(dp), intent(out) :: f_x(3, n_freq)
      real(dp), intent(out) :: f_p(3, n_freq)

      integer  :: k
      real(dp) :: eff_mass_omega_sq_inv, r_k_sq, r_k, prefactor
      real(dp), parameter :: EPS = 1.0e-12_dp

      !$omp parallel do private(k, eff_mass_omega_sq_inv, r_k_sq, r_k, prefactor)
      do k = 1, n_freq
         eff_mass_omega_sq_inv = 1.0_dp/(eff_mass(k)**2*omega(k)**2)

         r_k_sq = (x_k(1, k)**2 + x_k(2, k)**2 + x_k(3, k)**2) + &
                  (p_k(1, k)**2 + p_k(2, k)**2 + p_k(3, k)**2)*eff_mass_omega_sq_inv

         r_k = sqrt(r_k_sq)

         E_aux(k) = R_k*0.5_dp*eff_mass(k)*omega(k)**2
         if (r_k > EPS) then
            ! prefactor: -K * (R - R_target) / R

            prefactor = -k_bias(k)*(r_k - r_target(k))/r_k

            ! - d V' / d X_k
            f_x(1:3, k) = prefactor*x_k(1:3, k)

            ! gradient of bias     harmonic restoring force   coupling to dipole
            ! - d V' / d P_k   - eff_mass_k omega_k**2 X_k  + g_k M
            f_p(1:3, k) = prefactor*eff_mass_omega_sq_inv*p_k(1:3, k)
            ! Now add in the other term for the momenta force
            ! - mu omega**2 X
            f_p(:, k) = f_p(:, k) - eff_mass(k)*omega(k)**2*x_k(:, k) + g_k(k)*dipole

         else
            ! If amplitude is near zero
            f_x(1:3, k) = 0.0_dp
            f_p(1:3, k) = 0.0_dp
         end if
      end do
      !$omp end parallel do

   end subroutine calculate_ir_auxiliary_extended_variable_force_bias

   subroutine velocity_verlet_auxiliary(X, P, f_X, f_P, eff_mass, dt, first_step)
      real(dp), intent(inout) :: X(:, :)
      real(dp), intent(inout) :: P(:, :)
      real(dp), intent(in)    :: f_X(:, :)
      real(dp), intent(in)    :: f_P(:, :)
      real(dp), intent(in)    :: eff_mass(:)
      real(dp), intent(in)    :: dt
      logical, intent(in)    :: first_step

      integer :: n_freq, i, j

      n_sites = size(X, 2)

      ! 1. Momentum half-step (skipped on the very first initialization step if un-synced)
      if (.not. first_step) then
      do i = 1, n_sites
         do j = 1, 3
            P(j, i) = P(j, i) + 0.5_dp*f_P(j, i)*dt
         end do
      end do
      end if

      ! 2. Full position step
      ! X updates using velocity (P / mass) plus any position-dependent bias forces (f_X)
      do i = 1, n_freq
         do j = 1, 3
            X(j, i) = X(j, i) + dt*((P(j, i)/eff_mass(i)) + f_X(j, i))
         end do
      end do

   end subroutine velocity_verlet_auxiliary

   subroutine propagate_ir_auxiliary(step, dt, this, &
                                     n_sites, forces_atoms_aux, k_bias, dipole, dipole_jacobian)
      integer, intent(in) :: step
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: dipole(3)
      type(ir_auxiliary), intent(inout) :: this
      real(dp), intent(inout) :: forces_atoms_aux(:, :)
      real(dp), intent(inout) :: dipole_jacobian(:, :, :)
      real(dp), intent(in) :: k_bias(:)
      real(dp) :: E_eff(3)
      integer :: i, k, iter
      real(dp) :: R_k_sq
      real(dp) :: R_k
      real(dp) :: prefactor

      do k = 1, this%n_freq

         ! The Hamiltonian formulaton has the momenta velocity updates at the
         ! same time however we want to use the MD input with staggered functions

         if (first_step) then
            this%P_pred(:, k) = this%P(:, k)

            R_k_sq = sum(this%X(:, k)**2) &
                     + sum(this%P_pred(:, k)**2)/(this%eff_mass(k)**2*this%omega_k(k)**2)
            R_k = sqrt(max(R_k_sq, 1.0e-12_dp))

            this%E_aux(k) = R_k*0.5_dp*this%eff_mass(k)*this%omega(k)**2

            prefactor = -k_bias(k)*(R_k - this%R_target(k))/R_k

            this%f_aux(:, k) = -this%eff_mass(k)*this%omega(k)**2*this%X(:, k) &
                               + this%g_k(k)*dipole(:) &
                               + prefactor*this%X(:, k)
         else
            ! For all other steps, P_k is at time (t - dt).
            this%P_pred(:, k) = this%P(:, k) + (this%f_aux_prev(:, k))*dt

            ! Iterate to converge V(t) and F(t) simultaneously
            do iter = 1, 3
               ! Calculate amplitude envelope using the current guess for V(t)
               R_k_sq = sum(this%X(:, k)**2) &
                        + sum(this%P_pred(:, k)**2)/(this%eff_mass(k)**2*this%omega_k(k)**2)
               R_k = sqrt(max(R_k_sq, 1.0e-12_dp))

               this%E_aux(k) = R_k*0.5_dp*this%eff_mass(k)*this%omega(k)**2
               ! Calculate the trial force F(t)
               prefactor = -k_bias(k)*(R_k - this%R_target(k))/R_k
               this%f_aux(:, k) = -this%eff_mass(k)*this%omega(k)**2*this%X(:, k) &
                                  + this%g_k(k)*dipole(:) &
                                  + prefactor*this%X(:, k)

               ! Refine V(t) using the standard Trapezoidal/Verlet momentum rule
               this%P_pred(:, k) = this%P(:, k) &
                                   + 0.5_dp*(this%f_aux(:, k) + this%f_aux_prev(:, k))*dt
            end do
            ! By the end of this loop, f_aux is fully converged and mathematically
            ! consistent with the velocity your Verlet routine will generate.
         end if

      end do

      ! Compute the global effective field from the auxiliary positions E = sum(g_k * X_k)
      E_eff = 0.0_dp
      do k = 1, this%n_freq
         E_eff(:) = E_eff(:) + g_k(k)*X_k(:, k)
      end do

      ! Apply the local Jacobian to project the global field onto local atomic forces
      do i = 1, n_sites
         ! F_i = F_base,i + J_i^T * E_eff
         forces_atoms_aux(1, i) = forces_atoms_aux(1, i) &
                                  + dot_product(dipole_jacobian(1:3, 1, i), E_eff)
         forces_atoms_aux(2, i) = forces_atoms_aux(2, i) &
                                  + dot_product(dipole_jacobian(1:3, 2, i), E_eff)
         forces_atoms_aux(3, i) = forces_atoms_aux(3, i) &
                                  + dot_product(dipole_jacobian(1:3, 3, i), E_eff)
      end do

   end subroutine propagate_ir_auxiliary

! ==============================================================================
end module ir_auxiliary_dynamics
