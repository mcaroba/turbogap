! ------- option for doing simulation with adaptive time step			
!********* by Uttiyoarnab Saha

!**************************************************************************
! Adaptive Time Step for Upgrading from the existing variable time step 
! algorithm in TurboGap. The fix is applied after every specified number of
! time step(s). It checks whether the displacement (and the energy transfer)
! between atoms in a time step is below a limit specified by the user and 
! modifies the time step, i.e. dt, if these are not so.

module adaptive_time

contains

subroutine variable_time_step_adaptive (init, vel, forces, masses, tmin, tmax, xmax, emax, dt0, dt)

    implicit none

    real*8, intent(inout) :: dt
    real*8, intent(in) :: vel(:,:), dt0, forces(:,:), masses(:), tmin, tmax, xmax, emax
    logical, intent(in) :: init
    real*8, allocatable :: d(:)
    real*8 :: dtmin, vsq, fsq, dte, dtf, dtv
    integer :: Np, i

    Np = size(vel, 2)

    if( allocated( d ) )then
      if( size(d) /= Np )then
        deallocate( d )
        allocate( d(1:Np) )
      end if
    else
      allocate( d(1:Np) )
    end if

	dtmin = 1.0E+20
	
	do i = 1, Np
		vsq = dot_product (vel(1:3, i), vel(1:3, i))
		fsq = dot_product (forces(1:3, i), forces(1:3, i))
		if (vsq > 0.0) dtv = xmax / sqrt(vsq)
		if (fsq > 0.0) dtf = sqrt(2.0 * xmax * masses(i)/ (sqrt(fsq)))
		dt = min(dtv, dtf)
		
		if ((emax > 0.0) .and. (vsq*fsq > 0.0)) then
			dte = emax / sqrt(vsq*fsq)
			dt = min(dt, dte)
		end if
		
		! finally check the allowed maximum displacement per time-step
		! if the allowance in displacement is too large, then the criteria for 
		! maximum energy transfer per time-step will determine the value of time-step, dt 
		
		d(i) = sqrt( dot_product( vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2, &
					vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2 ) )
		
		if (d(i) > xmax) dt = dt * xmax / d(i)
		
		dtmin = min(dtmin, dt)
		
		if (tmin .ne. 0) dt = max(dt, tmin)
		if (tmax .ne. 0) dt = min(dt, tmax)
	end do

!	If we're at the first step (init) we use new dt as estimated above
!	and the initial time-step, dt0
   
    if ( init ) dt = min(dt, dt0)

end subroutine variable_time_step_adaptive

!**************************************************************************

end module adaptive_time

!******** until here for adaptive time

