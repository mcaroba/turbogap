! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, md.f90, is copyright (c) 2019-2021, Miguel A. Caro
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
                             forces, forces_prev, masses, dt, &
                             first_step, a_box, b_box, c_box)

    implicit none

!   Input variables
    real*8, intent(inout) :: positions(:,:), positions_prev(:,:), velocities(:,:), &
                             forces_prev(:,:)
    real*8, intent(in) :: forces(:,:), masses(:), dt, a_box(1:3), b_box(1:3), &
                          c_box(1:3)
    logical, intent(in) :: first_step
!   Internal variables
    integer :: n_sites, i

    n_sites = size(masses)

!   After this whole routine, velocities and positions_prev are synchronous. positions
!   is dt ahead of velocities

!   velocities are given at t-dt (except for the first step, when they're given for t); compute for t
    if( .not. first_step )then
      do i = 1, n_sites
        velocities(1:3, i) = velocities(1:3, i) + 0.5d0 * (forces(1:3, i) + forces_prev(1:3, i))/masses(i) * dt
      end do
    end if
!   positions are given at t; compute for t+dt
    positions_prev = positions
    forces_prev = forces
    do i = 1, n_sites
      positions(1:3, i) = positions(1:3, i) + velocities(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2
    end do

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


end module
