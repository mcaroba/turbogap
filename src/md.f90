! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, md.f90, is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module md


  use neighbors


  contains


!**************************************************************************
! Verlet is two subroutines
!
! Regular Verlet
  subroutine verlet(x_in, x_in_prev, F, m, dt, x_out)

    implicit none

    real*8, intent(in) :: x_in(1:3), x_in_prev(1:3)
    real*8, intent(out) :: x_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = 2.d0*x_in(1:3) - x_in_prev(1:3) + F(1:3)/m*dt**2

  end subroutine
  subroutine regular_verlet(positions, positions_prev, velocities, forces, masses, dt, first_step, &
                            a_box, b_box, c_box)

    implicit none

    real*8, intent(in) :: forces(:,:), masses(:), dt, a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(inout) :: positions(:,:), positions_prev(:,:)
    real*8, intent(inout) :: velocities(:,:)
    logical, intent(in) :: first_step
    integer :: natoms, i, i_shift(1:3)
    real*8 :: pos(1:3), d

    natoms = size(positions,2)

    if( first_step )then
      do i = 1, natoms
        positions_prev(1:3, i) = positions(1:3, i)
        positions(1:3, i) = positions(1:3, i) + velocities(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2
        velocities(1:3, i) = (positions(1:3, i) - positions_prev(1:3, i)) / dt
      end do
    else
      do i = 1, natoms
!       Make sure that we preserve the minimum image convention for the positions:
        call get_distance(positions_prev(1:3, i), positions(1:3, i), a_box, b_box, c_box, &
                          [.true., .true., .true.], pos(1:3), d, i_shift(1:3))
        positions(1:3, i) = positions_prev(1:3, i) + pos(1:3)
        call verlet(positions(1:3,i), positions_prev(1:3, i), forces(1:3,i), masses(i), dt, pos(1:3))
        positions_prev(1:3, i) = positions(1:3, i)
        positions(1:3, i) = pos(1:3)
!       This gives the velocity at t+dt/2, positions and velocities are not sinchronous in this
!       implementation
        velocities(1:3, i) = (positions(1:3, i) - positions_prev(1:3, i)) / dt
      end do
    end if

  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine velocity_verlet(positions, positions_prev, velocities, &
                             forces, forces_prev, masses, dt, dt_prev, &
                             first_step, a_box, b_box, c_box, fix_atom)

    implicit none

!   Input variables
    real*8, intent(inout) :: positions(:,:), positions_prev(:,:), velocities(:,:), &
                             forces_prev(:,:), dt_prev
    real*8, intent(in) :: forces(:,:), masses(:), dt, a_box(1:3), b_box(1:3), &
                          c_box(1:3)
    logical, intent(in) :: first_step, fix_atom(:,:)
!   Internal variables
    integer :: n_sites, i, j

    n_sites = size(masses)

!   After this whole routine, velocities and positions_prev are synchronous. positions
!   is dt ahead of velocities

!   velocities are given at t-dt (except for the first step, when they're given for t); compute for t
    if( .not. first_step )then
      do i = 1, n_sites
!        velocities(1:3, i) = velocities(1:3, i) + 0.5d0 * (forces(1:3, i) + forces_prev(1:3, i))/masses(i) * dt
        do j = 1, 3
          if( .not. fix_atom(j, i) )then
            velocities(j, i) = velocities(j, i) + 0.5d0 * (forces(j, i) + forces_prev(j, i))/masses(i) * dt_prev	! dt
          else
            velocities(j, i) = 0.d0
          end if
        end do
      end do
    end if
!   positions are given at t; compute for t+dt
    positions_prev = positions
    forces_prev = forces
    do i = 1, n_sites
!     positions(1:3, i) = positions(1:3, i) + velocities(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2
      do j = 1, 3
        if( .not. fix_atom(j, i) )then
          positions(j, i) = positions(j, i) + velocities(j, i)*dt + 0.5d0*forces(j, i)/masses(i)*dt**2
        end if
      end do
    end do
    
    dt_prev = dt	!! minimum modification for variable time-step situations 

  end subroutine
!**************************************************************************





!**************************************************************************
! Berendsen's velocity rescaling thermostat
!
  subroutine berendsen_thermostat(vel, T0, T, tau, dt)

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: T0, T, tau, dt
    real*8 :: f
    integer :: Np, i

    Np = size(vel, 2)

    f = dsqrt(1.d0 + dt/tau * (T0/T - 1.d0))

    if( T > 0.d0 )then
      do i = 1, Np
        vel(1:3, i) = vel(1:3, i) * f
      end do
    end if

  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine remove_cm_vel(vel, M)

!   I should adapt this code to mixed boundary conditions, where
!   the CM velocity can be removed per Cartesian dimension independently

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: M(:)
    real*8 :: cm_pos(1:3), cm_vel(1:3), total_mass
    integer :: Np, i

    Np = size(vel, 2)

    cm_vel = 0.d0
    total_mass = 0.d0
    do i = 1, Np
      cm_vel(1:3) = cm_vel(1:3) + M(i)*vel(1:3,i)
      total_mass = total_mass + M(i)
    end do
    cm_vel = cm_vel / total_mass
    do i = 1, Np
      vel(1:3,i) = vel(1:3,i) - cm_vel(1:3)
    end do

  end subroutine
!**************************************************************************





!**************************************************************************
  subroutine wrap_pbc(positions, a_box, b_box, c_box)

    implicit none

    real*8, intent(inout) :: positions(:,:)
    real*8, intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8 :: dist(1:3), d, mid(1:3)
    integer :: Np, i, i_shift(1:3)

    Np = size(positions,2)

!    mid(1:3) = 0.5d0 * (a_box(1:3) + b_box(1:3) + c_box(1:3))

!    do i = 1, Np
!      call get_distance( mid(1:3), positions(1:3, i), a_box(1:3), b_box(1:3), c_box(1:3), &
!                         [.true., .true., .true.], dist, d, i_shift(1:3))
!      positions(1:3, i) = mid(1:3) + dist(1:3) &
!                          - i_shift(1)*a_box(1:3) - i_shift(2)*b_box(1:3) - i_shift(3)*c_box(1:3)
!    end do

    do i = 1, Np
      call get_distance( [0.d0, 0.d0, 0.d0], positions(1:3, i), a_box(1:3), b_box(1:3), c_box(1:3), &
                         [.true., .true., .true.], dist, d, i_shift(1:3))
      positions(1:3, i) = positions(1:3, i) - i_shift(1)*a_box(1:3) - i_shift(2)*b_box(1:3) - i_shift(3)*c_box(1:3)
    end do

  end subroutine
!**************************************************************************






!**************************************************************************
  subroutine berendsen_barostat(positions, P0, P, sym, tau, gamma, dt)
!   Berendsen barostat that takes the bulk moduli ratio to that of water, gamma,
!   and takes P in bar
!   gamma = B / B_water; i.e., if the materials is very "hard", like diamond,
!   gamma >> 1; if the material is a gas, gamma << 1; if you have a liquid, then
!   gamma = 1 is probably a good choice.
    implicit none

    real*8, intent(inout) :: positions(:, :)
    real*8, intent(in) :: P0, P(1:3,1:3), tau, dt, gamma
    character(*), intent(in) :: sym
    real*8 :: P_iso
    integer :: i, n

    P_iso = (P(1,1) + P(2,2) + P(3,3))/3.d0

    if( sym(1:3) == "iso" )then
      positions = positions * (1.d0 + dt/tau * 4.5d-5/gamma * (P_iso - P0))**(1.d0/3.d0)
    else if( sym(1:4) == "diag" )then
      n = size(positions, 2)
      do i = 1, 3
        positions(i, 1:n) = positions(i, 1:n) * (1.d0 + dt/tau * 4.5d-5/gamma * (P(i,i) - P0))**(1.d0/3.d0)
      end do
    else
      write(*,*) "ERROR: I don't understand the specified barostat_sym keyword"
      stop
    end if

  end subroutine
!**************************************************************************








!**************************************************************************
  subroutine box_scaling(positions, a_box, b_box, c_box, indices, i_step, n_steps, gamma)

    implicit none

    real*8, intent(inout) :: positions(:, :), a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(in) :: gamma(3,3)
    integer, intent(in) :: n_steps, indices(:)
    real*8 :: f(3,3), identity(3,3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3,3])
    real*8, save :: a0(1:3), b0(1:3), c0(1:3)
    integer :: i_step
    integer, save :: indices0(1:3)

    if( i_step == 0 )then
      a0(1:3) = a_box(1:3)
      b0(1:3) = b_box(1:3)
      c0(1:3) = c_box(1:3)
      indices0(1:3) = indices(1:3)
    end if

!    positions = positions * ( 1.d0 + (gamma-1.d0) / dfloat(n_steps) )
!    f = 1.d0 + (gamma-1.d0) * dfloat(i_step+1) / dfloat(n_steps)
!    a_box = a0 * f / dfloat(indices0(1)) * dfloat(indices(1))
!    b_box = b0 * f / dfloat(indices0(2)) * dfloat(indices(2))
!    c_box = c0 * f / dfloat(indices0(3)) * dfloat(indices(3))
    positions = positions + matmul(gamma-identity, positions) / dfloat(n_steps)
    f = identity + (gamma-identity) * dfloat(i_step+1) / dfloat(n_steps)
    a_box = matmul(f, a0) / dfloat(indices0(1)) * dfloat(indices(1))
    b_box = matmul(f, b0) / dfloat(indices0(2)) * dfloat(indices(2))
    c_box = matmul(f, c0) / dfloat(indices0(3)) * dfloat(indices(3))

  end subroutine 
!**************************************************************************






!**************************************************************************
! Custom variable time step algorithm
!
  subroutine variable_time_step(init, vel, forces, masses, target_pos_step, tau_dt, dt0, dt)

    implicit none

    real*8, intent(inout) :: dt
    real*8, intent(in) :: vel(:,:), target_pos_step, dt0, tau_dt, forces(:,:), masses(:)
    logical, intent(in) :: init
    real*8, allocatable :: d(:)
    real*8 :: new_dt, d_max, dt_prev
    integer :: Np, i, i_max
    logical :: optimize_time_step, too_large

    Np = size(vel, 2)

    if( allocated( d ) )then
      if( size(d) /= Np )then
        deallocate( d )
        allocate( d(1:Np) )
      end if
    else
      allocate( d(1:Np) )
    end if


    do i = 1, Np
      d(i) = sqrt( dot_product( vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2, &
                                vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2 ) )
    end do
    d_max = maxval(d)
    if( d_max > target_pos_step )then
      too_large = .true.
    else
      too_large = .false.
    end if

!   new_dt estimates the optimal time step for this snapshot (within 1% accuracy)
    new_dt = dt
    optimize_time_step = .true.
    do while( optimize_time_step )
      do i = 1, Np
        d(i) = sqrt( dot_product( vel(1:3, i)*new_dt + 0.5d0*forces(1:3, i)/masses(i)*new_dt**2, &
                                  vel(1:3, i)*new_dt + 0.5d0*forces(1:3, i)/masses(i)*new_dt**2 ) )
      end do
      d_max = maxval(d)
!      i_max = maxloc(d)
      if( d_max > target_pos_step .and. too_large )then
        new_dt = new_dt * 0.99d0
      else if( d_max < target_pos_step .and. .not. too_large )then
        new_dt = new_dt * 1.01d0
      else
        optimize_time_step = .false.
      end if
    end do

    if( init )then
!     If we're at the first step (init) we use new_dt as estimated above
      dt = min(new_dt, dt0)
    else
!     Otherwise we use a Berendsen approach
      new_dt = dt * tau_dt / (tau_dt + dt - new_dt)
      dt = min(new_dt, dt0)
    end if


  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine gradient_descent(positions, positions_prev, velocities, &
                              forces, forces_prev, masses, max_opt_step, &
                              first_step, a_box, b_box, c_box, fix_atom, energy)

    implicit none

!   Input variables
    real*8, intent(inout) :: positions(:,:), positions_prev(:,:), velocities(:,:), &
                             forces_prev(:,:), forces(:,:)
    real*8, intent(in) :: masses(:), a_box(1:3), b_box(1:3), &
                          c_box(1:3), max_opt_step, energy
    logical, intent(in) :: fix_atom(:,:), first_step
!   Internal variables
    real*8 :: gamma, max_force, this_force, pos(1:3), d
    real*8, save :: gamma_prev, energy0, m_prev, gamma_back0
    real*8, allocatable, save :: positions0(:,:), forces0(:,:)
    integer :: n_sites, i, j, i_shift(1:3)
    logical, save :: backtracking, initialized = .false.

    n_sites = size(masses)

!   Here we always set the velocities to zero
    velocities = 0.d0

    if( first_step )then
      backtracking = .true.
      if( .not. allocated(positions0) )allocate( positions0(1:3, 1:size(positions,2)) )
      if( .not. allocated(forces0) )allocate( forces0(1:3, 1:size(positions,2)) )
      positions0 = positions
      forces0 = forces
      energy0 = energy
!     The first step is (over)estimated from user provided values
      max_force = 0.d0
      do i = 1, n_sites
        this_force = sqrt( dot_product(forces(1:3,i), forces(1:3,i)) )
        if( this_force > max_force )then
          max_force = this_force
        end if
      end do
      if( max_force == 0.d0 )then
        gamma = 0.d0
      else if( initialized )then
        gamma = gamma_back0
      else
!        gamma = max(gamma0, max_opt_step/max_force)
        gamma = max_opt_step/max_force
      end if
    else if( backtracking )then
!     After the first step, we perform backtracking line search until fullfilling the
!     Armijo-Goldstein condition
      if( energy <= energy0 - gamma_prev*0.5d0*m_prev )then
        backtracking = .false.
        initialized = .true.
        gamma_back0 = gamma_prev
      else
!       If the condition is not fulfilled, we restore the original positions and decrease
!       the step by half
        gamma = gamma_prev * 0.5d0
        positions = positions0
        forces = forces0
      end if
    end if

    if( .not. first_step .and. .not. backtracking )then
!     Make sure we use the same image convention for positions and positions_prev
      do i = 1, n_sites
        call get_distance(positions_prev(1:3, i), positions(1:3, i), a_box, b_box, c_box, &
                          [.true., .true., .true.], pos(1:3), d, i_shift(1:3))
        positions_prev(1:3, i) = positions(1:3, i) - pos(1:3)
      end do
!     Barzilai–Borwein method for finding gamma
      gamma = sum( (positions-positions_prev) * (forces-forces_prev)) / sum( (forces-forces_prev)**2 )
      gamma = abs( gamma )
    end if

    positions_prev = positions
    forces_prev = forces

    do i = 1, n_sites
      do j = 1, 3
        if( .not. fix_atom(j, i) )then
          positions(j, i) = positions_prev(j, i) + gamma*forces_prev(j, i)
        end if
      end do
    end do

    gamma_prev = gamma
    m_prev = sum( forces**2 )

  end subroutine
!**************************************************************************






!**************************************************************************
  subroutine volume_preserving_strain_transformation(a_box, b_box, c_box, gamma)

    implicit none

    real*8, intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(inout) :: gamma(3,3)
    real*8 :: identity(3,3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3,3]), &
              a(1:3), b(1:3), c(1:3), vol, vol_p, vol_ratio

    a = a_box + matmul(gamma-identity, a_box)
    b = b_box + matmul(gamma-identity, b_box)
    c = c_box + matmul(gamma-identity, c_box)

    vol = dot_product(cross_product(a_box, b_box), c_box)
    vol_p = dot_product(cross_product(a, b), c)

    vol_ratio = (vol / vol_p)**(1.d0/3.d0)

    gamma = gamma * vol_ratio

  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine get_ns_unbiased_volume_proposal(V1, V2, n_sites, V)

!   V1 is the minimal volume, V2 is the maximal volume, V is the volume proposal
!   The likelihood of a volume proposal should scale as V^n_sites. Therefore
!   larger volumes are favored. We use a homogeneous random distribution within
!   (0,1) and scale it so that the integrated probability density matches that
!   of V^n_sites within [V1, V2], before making the assignment between the two
    !   distributions

    implicit none
    real*8, intent(in) :: V1, V2
    integer, intent(in) :: n_sites
!   Output variables
    real*8, intent(out) :: V
!   Internal variables
    real*8 :: rand, log_V_V2, v_ratio

    call random_number( rand )

    v_ratio = (V1/V2)**(n_sites+1)

    log_V_V2 = dlog(rand + (1.d0-rand)*v_ratio) / dfloat(n_sites + 1)

    V = V2 * dexp(log_V_V2)

  end subroutine






  subroutine gradient_descent_box(positions, positions_prev, velocities, &
                                  forces, forces_prev, masses, max_opt_step_eps, &
                                  first_step, a_box, b_box, c_box, energy, &
                                  virial, optim_mode, restart )
!                                  virial, optim_mode, n_restart, restart )
    implicit none
!   Input variables
    real*8, intent(inout) :: positions(:,:), positions_prev(:,:), velocities(:,:), &
                             forces_prev(:,:), a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, intent(in) :: forces(:,:), masses(:), max_opt_step_eps, energy, virial(1:6)
!    integer, intent(inout) :: n_restart
    integer :: n_restart
    logical, intent(in) :: first_step
    character*16, intent(in) :: optim_mode
!   Output variables
    logical :: restart
!   Internal variables
    real*8 :: max_force, this_force, pos(1:3), d, gamma_eps, t_eps(1:3, 1:3)
    real*8, allocatable, save :: frac_pos(:,:), frac_pos_prev(:,:)
    real*8, save :: energy0, m_prev, a_box0(1:3), b_box0(1:3), c_box0(1:3), &
                    eps(1:6), eps_prev(1:6), gamma_eps_prev, &
                    m_eps_prev, virial_prev(1:6), virial0(1:6), this_virial(1:6), &
                    gamma_back0
    real*8, allocatable, save :: positions0(:,:)
    integer :: n_sites, i, j, i_shift(1:3)
    integer, save :: i_restart
    logical, save :: backtracking, initialized = .false.

    n_sites = size(masses)

!   Here we always set the velocities to zero
    velocities = 0.d0

    this_virial = virial

!   HARDCODED FOR NOW
    n_restart = 10
    if( n_restart < 2 )then
      n_restart = 2
    end if

    if( first_step )then
      i_restart = 0
      backtracking = .true.
      if( .not. allocated(frac_pos) )allocate( frac_pos(1:3, 1:n_sites) )
      if( .not. allocated(frac_pos_prev) )allocate( frac_pos_prev(1:3, 1:n_sites) )
      if( .not. allocated(positions0) )allocate( positions0(1:3, 1:size(positions,2)) )
      positions0 = positions
      virial0 = this_virial
      energy0 = energy
      a_box0 = a_box
      b_box0 = b_box
      c_box0 = c_box
      eps = 0.d0
!     The first step is (over)estimated from user-provided values
      max_force = 0.d0
      do i = 1, 6
        this_force = sqrt( this_virial(i)**2 )
        if( this_force > max_force )then
          max_force = this_force
        end if
      end do
      if( max_force == 0.d0 )then
        gamma_eps = 0.d0
      else if( initialized )then
        gamma_eps = gamma_back0
      else
        gamma_eps = max_opt_step_eps / max_force
      end if
    else if( backtracking )then
!     After the first step, we perform backtracking line search until fullfilling the
!     Armijo-Goldstein condition
      if( energy <= energy0 - gamma_eps_prev*0.5d0*m_eps_prev )then
        backtracking = .false.
        initialized = .true.
        gamma_back0 = gamma_eps_prev
      else
!       If the condition is not fulfilled, we restore the original positions and decrease
!       the step by half
        gamma_eps = gamma_eps_prev * 0.5d0
        gamma_back0 = gamma_eps
        a_box = a_box0
        b_box = b_box0
        c_box = c_box0
        eps = 0.d0
        positions = positions0
        this_virial = virial0
      end if
    end if

!   Transform positions to fractional coordinate system
    call get_fractional_coordinates(positions, a_box, b_box, c_box, frac_pos)

    if( .not. first_step .and. .not. backtracking )then
!     Make sure we use the same image convention for positions and positions_prev
      do i = 1, n_sites
        call get_distance(frac_pos_prev(1:3, i), frac_pos(1:3, i), [1.d0, 0.d0, 0.d0], [0.d0, 1.d0, 0.d0], &
                          [0.d0, 0.d0, 1.d0], [.true., .true., .true.], pos(1:3), d, i_shift(1:3))
        frac_pos_prev(1:3, i) = frac_pos(1:3, i) - pos(1:3)
      end do
!     Barzilai–Borwein method for finding gamma
      gamma_eps = sum( (eps(:)-eps_prev(:)) * (this_virial(:)-virial_prev(:)) ) / &
                  sum( (this_virial(:)-virial_prev(:))**2 )
      gamma_eps = abs( gamma_eps )
    end if

    virial_prev = this_virial
    eps_prev = eps
    frac_pos_prev = frac_pos
    do i = 1, n_sites
      positions(1:3, i) = frac_pos(1, i) * a_box(1:3) + frac_pos(2, i) * b_box(1:3) + &
                          frac_pos(3, i) * c_box(1:3)
    end do

    positions_prev = positions
    forces_prev = forces

    eps(1:6) = eps_prev(1:6) + gamma_eps * virial_prev(1:6)
    gamma_eps_prev = gamma_eps
    m_eps_prev = sum( virial_prev(1:6)**2 )

    if( optim_mode == "gd-box-ortho" )then
      eps(4) = 0.d0
      eps(5) = 0.d0
      eps(6) = 0.d0
    end if

    t_eps(1:3, 1) = [1.d0 + eps(1), eps(6)/2.d0, eps(5)/2.d0]
    t_eps(1:3, 2) = [eps(6)/2.d0, 1.d0 + eps(2), eps(4)/2.d0]
    t_eps(1:3, 3) = [eps(5)/2.d0, eps(4)/2.d0, 1.d0 + eps(3)]
    a_box = matmul( t_eps, a_box0 )
    b_box = matmul( t_eps, b_box0 )
    c_box = matmul( t_eps, c_box0 )

    do i = 1, n_sites
      positions(1:3, i) = frac_pos(1, i) * a_box(1:3) + frac_pos(2, i) * b_box(1:3) + &
                          frac_pos(3, i) * c_box(1:3)
    end do

    i_restart = i_restart + 1
    if( i_restart >= n_restart .and. energy < energy0 )then
      restart = .true.
    else
      restart = .false.
    end if
  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine get_atomic_mass( element, mass, is_in_database )

    implicit none

!   Input variables
    character*8, intent(in) :: element
!   Output variables
    real*8, intent(out) :: mass
    logical, intent(out) :: is_in_database
!   Internal variables
    real*8 :: masses(1:96)
    character*8 :: elements(1:96)
    integer :: i

    elements = [" H", "He", "Li", "Be", " B", " C", " N", " O", &
                " F", "Ne", "Na", "Mg", "Al", "Si", " P", " S", &
                "Cl", " K", "Ar", "Ca", "Sc", "Ti", " V", "Cr", &
                "Mn", "Fe", "Ni", "Co", "Cu", "Zn", "Ga", "Ge", &
                "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", &
                "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", &
                "In", "Sn", "Sb", " I", "Te", "Xe", "Cs", "Ba", &
                "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", &
                "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", &
                "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
                "Th", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", &
                "Ac", "Pa", "Th", "Np", " U", "Pu", "Am", "Cm"]

    masses = [  1.00797,   4.00260,  6.941,    9.01218, 10.81,    12.011,  14.0067,  15.9994, &
               18.998403, 20.179,   22.98977, 24.305,   26.98154, 28.0855, 30.97376, 32.06, &
               35.453,    39.0983,  39.948,   40.08,    44.9559,  47.90,   50.9415,  51.996, &
               54.9380,   55.847,   58.70,    58.9332,  63.546,   65.38,   69.72,    72.59, &
               74.9216,   78.96,    79.904,   83.80,    85.4678,  87.62,   88.9059,  91.22, &
               92.9064,   95.94,    98.,     101.07,   102.9055, 106.4,   107.868,  112.41, &
              114.82,    118.69,   121.75,   126.9045, 127.60,   131.30,  132.9054, 137.33, &
              138.9055,  140.12,   140.9077, 144.24,   145.,     150.4,   151.96,   157.25, &
              158.9254,  162.50,   164.9304, 167.26,   168.9342, 173.04,  174.967,  178.49, &
              180.9479,  183.85,   186.207,  190.2,    192.22,   195.09,  196.9665,  200.59, &
              204.37,    207.2,    208.9804, 209.,     210.,     222.,    223.,      226.0254, &
              227.0278,  231.0359, 232.0381, 237.0482, 238.029,  244.,    243.,      247.]

    is_in_database = .false.

    mass = 0.d0
    do i = 1, size(elements, 1)
      if( trim(adjustl(elements(i))) == trim(adjustl(element)) )then
        mass = masses(i)
        is_in_database = .true.
        exit
      end if
    end do

  end subroutine
!**************************************************************************

end module
