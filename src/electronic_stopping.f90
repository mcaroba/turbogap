! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, electronic_stopping.f90, is copyright (c) 2023, Uttiyoarnab Saha
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

! ------- option for radiation cascade simulation with velocity-dependent electronic stopping				
! The electron stopping should be calculated in every time step if the simulation 
! is done by invoking this option of electron stopping. The balance forces on the atoms 
! after electronic energy loss is returned on the timestep. The total electronic
! energy loss on the time step is also written to output file for calculation of total
! cumulative energy loss through electronic stopping during post processing or analyses.
! This calculation is done here based on the tabular data (given through an input data file)
! of energy versus electronic energy loss for each species of atoms in the system.
! Such data can be provided as the input by using the SRIM-2013 software (for example).
! The data must be in a particular (simple) format: 1st line - information;
! 2nd line - number of rows of data points; 3rd line - energy unit (eV) and symbols of
! elements (must be) in order as in input file; 4th line onwards - energy and electronic
! stopping values (eV/Ang) like E_1 dE/dx(element 1) dE/dx(element 2) .. (.. 3) ..
!								E_2  "   (.. 1)			" (.. 2)	   " (.. 3) ..
!								..	..		..			..				..		..
!********* by Uttiyoarnab Saha

!**************************************************************************

module electronic_stopping

contains

subroutine electron_stopping_velocity_dependent (md_istep, natomtypes, Ecut, eel_freq_out, &
			vel, forces, masses, type_mass, dt, md_time, nrows, allelstopdata, to_calculate)

	implicit none
	
	real*8, intent(in) :: vel(:,:), dt, md_time 
	real*8, intent(inout) :: forces(:,:)
	integer, intent(in) :: md_istep, nrows, eel_freq_out
	character*6, intent(in) :: to_calculate
	integer :: Np, i, j, itype
	real*8 :: vsq, energy, Se, Se_lo, Se_hi, E_lo, E_hi, factor, vabs, SeLoss
	
	!! Data from the user-given input script file
	
	integer, intent(in) :: natomtypes
	real*8, intent(in) :: Ecut,  masses(:), type_mass(:)
	
	!! To receive stopping data from file to array containers
	
	integer :: ncols, ndata, irow		!icol 
	real*8, allocatable, intent(in) :: allelstopdata(:)
	real*8, allocatable :: En_elstopfile(:), elstop(:,:)
	character*1024 :: infoline
	
	!! To write the electronic energy loss data to a file after evaluation at each time step
	
	if( md_istep == 0 .or. md_istep == -1 )then
		open (unit = 100, file = "ElectronicEnergyLoss.txt", status = "unknown")
		write(100,*) 'Time (fs)  Electronic energy loss (eV)'
		write(100,1) md_time, 0.0
	else
		if ((to_calculate == 'energy') .and. &
				(MOD(md_istep, eel_freq_out) == 0)) &
					open (unit = 100, file = "ElectronicEnergyLoss.txt", &
											status = "old", position = "append")
	end if


	!! Do all calculations for MD time greater than 0
	
	if ( .not. (md_istep == 0 .or. md_istep == -1) ) then

		ncols = natomtypes + 1
		ndata = nrows*ncols
		allocate (En_elstopfile(nrows), elstop(nrows, ncols-1))
		
		irow = 1
		
		do i = 1, ndata, ncols
			En_elstopfile(irow) = allelstopdata(i)
			elstop(irow,:) = allelstopdata(i+1:i+natomtypes)
			irow = irow + 1	
		end do
		
		Np = size(vel, 2)
		
		SeLoss = 0.0
		
		do i = 1, Np
			Se = 0.0
			do j = 1, natomtypes
				if (masses(i) == type_mass(j)) then
					itype = j
					exit
				end if
			end do
			
			vsq = dot_product (vel(1:3, i), vel(1:3, i))
			energy = 0.5 * masses(i) * vsq
		
			if (energy < Ecut) continue
			if (energy < En_elstopfile(1)) continue
			if (energy > En_elstopfile(nrows)) then
				write (*,*) "ERROR: Kinetic energy of atom is higher than electron stopping data"
				stop
			end if
			!! Here the condition for matching the region can be given
			!! .... like a particular group atoms only .....
			!! .... !
			
			!! Find position of atom K.E in the data file and then corresponding electronic stopping
			!! to apply the friction to the current forces
			
			do j = 1, nrows
				if (energy == En_elstopfile(j)) then 
					Se = elstop(j, itype)
					exit
				end if
				if (En_elstopfile(j) < energy .and. energy < En_elstopfile(j+1)) then
					Se_lo = elstop(j, itype)
					Se_hi = elstop(j+1, itype)
					E_lo = En_elstopfile(j)
					E_hi = En_elstopfile(j+1)
					Se = Se_lo + (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo)
					exit
				end if
			end do
			
			vabs = sqrt(vsq)
			factor = -Se / vabs
			
			if (to_calculate == 'forces') then
				!! The current forces get modified (reduced) 
				forces(1,i) = forces(1,i) + vel(1,i) * factor
				forces(2,i) = forces(2,i) + vel(2,i) * factor
				forces(3,i) = forces(3,i) + vel(3,i) * factor
			end if
			
			if (to_calculate == 'energy') then
				!! roughly, E = E + (dE/dx) * (dx) = (Se) * (vabs*dt)
				SeLoss = SeLoss + (Se * vabs * dt)
			end if
		end do
		
		if ((to_calculate == 'energy') .and. &
							(MOD(md_istep, eel_freq_out) == 0)) then
			write(100, 1) md_time, SeLoss
			close(unit = 100)
		end if
1 FORMAT (E13.6, E14.6)

	end if 	!! when md_time > 0.0

end subroutine electron_stopping_velocity_dependent

!**************************************************************************

end module electronic_stopping

!******** until here for electronic stopping
