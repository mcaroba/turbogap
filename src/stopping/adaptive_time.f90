! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, adaptive_time.f90, is copyright (c) 2023, Uttiyoarnab Saha
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

!**************************************************************************

! ------- option for doing simulation with adaptive time step			
! Adaptive Time Step for Upgrading from the existing variable time step 
! algorithm in TurboGap. The is applied after every specified number of
! time step(s). It checks whether the displacement (and the energy transfer)
! between atoms in a time step is below a limit specified by the user and 
! modifies the time step, i.e. dt, if these are not so.
!********* by Uttiyoarnab Saha

!**************************************************************************

module adaptive_time

contains

subroutine variable_time_step_adaptive (init, vel, forces, masses, &
			tmin, tmax, xmax, emax, dt0, dt, adapt_time_for_atoms)

	implicit none

	real*8, intent(inout) :: dt
	real*8, intent(in) :: vel(:,:), dt0, forces(:,:), masses(:), tmin, tmax, xmax, emax
	logical, intent(in) :: init
	integer, intent(in) :: adapt_time_for_atoms(:)
	real*8 :: dist, dtmin, vsq, fsq, dte, dtf, dtv
	integer :: Np, ki, i

	!! checking the input values of calculation parameters
	
	if (xmax <= 0.0) then
		write(*,*) "ERROR: value of xmax must be greater than 0."
		stop
	end if
	if (emax <= 0.0) then
		write(*,*) "WARNING: value of emax is given less or equal to 0."
		write(*,*) "Assuming default value for emax = 10.0 eV."
	end if
	if (tmax <= 0.0) then
		write(*,*) "ERROR: value of tmax must be greater than 0."
		stop
	end if
	if (tmin <= 0.0) then
		write(*,*) "ERROR: value of tmin must be greater than 0."
		stop
	end if
	if (tmax <= tmin) then
		write(*,*) "ERROR: tmax must be greater than tmin."
		stop
	end if
	
	
	Np = size(adapt_time_for_atoms)		!! number of atoms in specified group
	
	!! this will finally keep the smallest time-step required
	
	dtmin = 1.0E+20
	
	!! Initializing time steps to choose the minimum from a proper set
	!! of values of times and avoid getting 0.
	
	if ( init ) then
		dtv = dt0; dtf = dt0; dte = dt0
	else
		dtv = dt; dtf = dt; dte = dt
	end if

	do ki = 1, Np
		i = adapt_time_for_atoms(ki)
		!! velocity is needed for both xmax and emax criteria
		vsq = dot_product (vel(1:3, i), vel(1:3, i))
		!! force is needed for emax criteria
		fsq = dot_product (forces(1:3, i), forces(1:3, i))
		
		!! obtain the smaller dt between xmax/v and sqrt(2xmax/a)
		
		!! time from the velocity, x = vt
		if (vsq > 0.0) dtv = xmax / sqrt(vsq)
		
		!! time from the acceleration, x = 1/2at²
		if (fsq > 0.0) dtf = sqrt(2.0 * xmax * masses(i)/ (sqrt(fsq)))
		
		dt = min(dtv, dtf)
		
		!! apply the criterion emax
		!! time from the energy, t = e/(fv)
		if ((emax > 0.0) .and. (vsq*fsq > 0.0)) then
			dte = emax / sqrt(vsq*fsq)
			dt = min(dt, dte)
		end if
		
		!! finally check the allowed maximum displacement per time-step
		!! if the allowance in displacement is too large, then the criteria for 
		!! maximum energy transfer per time-step will determine the value of time-step, dt 
		
		!! calculate the distance that an atom moves with the present v, f, dt
		dist = sqrt( dot_product( vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2, &
					vel(1:3, i)*dt + 0.5d0*forces(1:3, i)/masses(i)*dt**2 ) )
		
		!! apply the criterion xmax, if necessary
		!! rescale dt by the xmax when the above calculated distance is larger than xmax  
		if (dist > xmax) dt = dt * xmax / dist

		dtmin = min(dtmin, dt)
	end do
	
	!! time-step is always within the maximum and minimum specified limits
	
	if (tmin .ne. 0) dt = max(dtmin, tmin)
	if (tmax .ne. 0) dt = min(dtmin, tmax)
	
	!! warnings for specified time limits, based on the warnings user can change the min. and max. t-limits
	!! no need to print these warnings as these consumes time
	
	!if (dtmin < tmin) then
	!	write(*,*) "WARNING: given xmax or emax criterion demands even lower value of tmin,"
	!	write(*,*) "(", dtmin, ")", " but doing with tmin ...."
	!end if
	!
	!if (dtmin > tmax) then
	!	write(*,*) "WARNING: given xmax or emax criterion demands even larger value of tmax,"
	!	write(*,*) "(", dtmin, ")", " but doing with tmax ...."
	!end if

!!	If we're at the first step (init) we use new dt as estimated above
!!	and the initial time-step, dt0
   
    if ( init ) dt = min(dt, dt0)

end subroutine variable_time_step_adaptive

end module adaptive_time

!******** until here for adaptive time

