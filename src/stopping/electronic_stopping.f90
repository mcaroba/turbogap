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

use mpi

type electronic_stopping_scalar_class
	integer :: nrows_esdata, ncols_esdata, elstop_length, md_last_step
	real*8, allocatable :: En_elstopfile(:), elstop(:)
	real*8 :: cum_EEL, md_prev_time
	contains
	procedure :: read_electronic_stopping_file, electron_stopping_velocity_dependent
	procedure :: electronic_stopping_scalar_broadcast
	procedure :: electronic_stopping_InitialValues
end type electronic_stopping_scalar_class

contains

!! initial values
subroutine electronic_stopping_InitialValues(this, E_eel_cumulative_prev, md_last_step, md_prev_time)
	implicit none
	class (electronic_stopping_scalar_class) :: this
	integer, intent(in) :: md_last_step
	real*8, intent(in) :: E_eel_cumulative_prev, md_prev_time
	this%cum_EEL = E_eel_cumulative_prev + 0.0d0
	this%md_last_step = md_last_step
	this%md_prev_time = md_prev_time
end subroutine electronic_stopping_InitialValues

!! read and store all electronic stopping power data
!! also give error messages if the data in the file is not in proper format
subroutine read_electronic_stopping_file (this, rank, ierr, n_species, species_types, estopfilename)

	implicit none
	class (electronic_stopping_scalar_class) :: this
	
	character*1024, intent(in) :: estopfilename
	integer, intent(in) :: rank, n_species
	integer :: ierr
	character*8, intent(in) :: species_types(n_species)
	real*8, allocatable :: allelstopdata(:)
	
	character*6, allocatable :: infoline(:)
	integer :: i, j, irow, ndata_esdata
	
	open (unit = 1000, file = estopfilename)
	! first line gives information
	! second line gives number of energy-stopping data points, i.e no. of rows of data
	read(1000,*)
	read(1000,*) this%nrows_esdata
	if (this%nrows_esdata <= 0) then
		if (rank == 0) write(*,*) "ERROR: Number of data rows in stopping file is 0 or less."
		call mpi_finalize(ierr)
		stop
	end if
	this%ncols_esdata = n_species + 1
	allocate(infoline(this%ncols_esdata))
	! third line gives energy units, names of elements in order of the atom species types in input file
	read(1000,*) (infoline(i), i = 1, this%ncols_esdata)
	do i = 2, this%ncols_esdata
		if (trim(infoline(i)) /= trim(species_types(i-1))) then
			if (rank == 0) write(*,*) "ERROR: Stopping powers for Elements are not given in order."
			call mpi_finalize(ierr)
			stop
		end if
	end do
	ndata_esdata = this%nrows_esdata*this%ncols_esdata
	allocate (allelstopdata(ndata_esdata))
	read(1000,*) (allelstopdata(i), i = 1, ndata_esdata)

	close(unit = 1000)

	this%elstop_length = this%nrows_esdata*(this%ncols_esdata-1)
	allocate (this%En_elstopfile(this%nrows_esdata), this%elstop(this%elstop_length))

	!! energy values
	irow = 1
	do i = 1, ndata_esdata, this%ncols_esdata
		this%En_elstopfile(irow) = allelstopdata(i)
		irow = irow + 1
	end do
	!! ES values
	do i = 1, n_species
		irow = 1
		do j = i+1, ndata_esdata, this%ncols_esdata
			this%elstop((i-1)*this%nrows_esdata + irow) = allelstopdata(j)
			irow = irow + 1
		end do
	end do

end subroutine read_electronic_stopping_file


!! broadcast required quantities
subroutine electronic_stopping_scalar_broadcast(this, ierr)

	implicit none
	class (electronic_stopping_scalar_class) :: this
	integer :: ierr

	call mpi_bcast (this%nrows_esdata, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call mpi_bcast (this%ncols_esdata, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call mpi_bcast (this%En_elstopfile, this%nrows_esdata, MPI_DOUBLE_PRECISION, 0, &
					MPI_COMM_WORLD, ierr)
	call mpi_bcast (this%elstop, this%elstop_length, MPI_DOUBLE_PRECISION, 0, &
					MPI_COMM_WORLD, ierr)
	call mpi_bcast (this%cum_EEL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

end subroutine electronic_stopping_scalar_broadcast


!! find forces after electronic stopping and electronic energy loss
subroutine electron_stopping_velocity_dependent (this, md_istep, num_md_steps, natomtypes, Ecut, eel_freq_out, &
			vel, forces, masses, type_mass, dt, md_time, eel_for_atoms, Np, to_calculate, rank, ntasks, ierr)

	implicit none
	class (electronic_stopping_scalar_class) :: this
	real*8, intent(in) :: vel(:,:), dt, md_time 
	real*8, intent(inout) :: forces(:,:)
	integer, intent(in) :: md_istep, num_md_steps, eel_freq_out, eel_for_atoms(:), Np, rank, ntasks
	character*6, intent(in) :: to_calculate
	integer :: Np_all, ki, i, j, itype, atom_type_start
	real*8 :: vsq, energy, Se, Se_lo, Se_hi, E_lo, E_hi, factor = 0.0d0, vabs, SeLoss, &
	E_kinetic_atoms, instant_temp_atoms, SeLoss_temp, E_kinetic_atoms_temp
	real*8, parameter :: boltzconst = 8.61733326E-05
	integer :: istart, istop, ierr

	!! Data from the user-given input script file

	integer, intent(in) :: natomtypes
	real*8, intent(in) :: Ecut, masses(:), type_mass(:)

	Np_all = size(vel,2)

	!! To write the electronic energy loss data to a file after evaluation at each time step

IF (rank == 0) THEN
	if (to_calculate == 'energy') then
		if( md_istep == 0 .or. md_istep == -1 )then
			E_kinetic_atoms = 0.0d0
			do ki = 1, Np
				i = eel_for_atoms(ki)
				E_kinetic_atoms = E_kinetic_atoms + &
								0.5d0 * masses(i) * dot_product(vel(1:3, i), vel(1:3, i))
			end do
			instant_temp_atoms = 2.d0/3.d0/dfloat(Np-1)/boltzconst*E_kinetic_atoms

			if ( this%md_last_step == 0 ) then
				open (unit = 100, file = "ElectronicEnergyLoss.txt", status = "unknown")
				write(100,*) 'Time (fs) / EEL (eV) / Cum.EEL (eV) / Kin.E_a (eV) / T_a (K)'
				write(100,1) md_time, 0.0, 0.0, E_kinetic_atoms, instant_temp_atoms
			end if
		else
			if ((MOD(md_istep, eel_freq_out) == 0) .or. md_istep == num_md_steps) &
				open (unit = 100, file = "ElectronicEnergyLoss.txt", &
												status = "old", position = "append")
		end if
	end if
END IF

	!! Do all calculations for MD time greater than 0

	if ( .not. (md_istep == 0 .or. md_istep == -1) ) then

		!! Finding the forces after friction

		if (to_calculate == 'forces') then
IF (rank /= 0) forces = 0.0d0
			call iterationRangesToProcesses (1, Np, ntasks, rank, istart, istop)

			do ki = istart, istop
				i = eel_for_atoms(ki)
				Se = 0.0
				do j = 1, natomtypes
					if (masses(i) == type_mass(j)) then
						itype = j
						atom_type_start = (itype-1)*this%nrows_esdata
						exit
					end if
				end do

				vsq = dot_product (vel(1:3, i), vel(1:3, i))
				energy = 0.5d0 * masses(i) * vsq

				if (energy < Ecut) cycle
				if (energy < this%En_elstopfile(1)) cycle
				if (energy > this%En_elstopfile(this%nrows_esdata)) then
					if (rank == 0) write (*,*) trim("ERROR: Kinetic energy "), energy, &
								trim("eV of atom is higher than electron stopping data")
					call mpi_finalize(ierr)
					stop
				end if

				!! Find position of atom K.E in the data file and then corresponding
				!! electronic stopping to apply the friction to the current forces

				do j = 1, this%nrows_esdata
					if (energy == this%En_elstopfile(j)) then 
						Se = this%elstop(atom_type_start + j)
						exit
					end if
					if (this%En_elstopfile(j) < energy .and. energy < this%En_elstopfile(j+1)) then
						Se_lo = this%elstop(atom_type_start + j)
						Se_hi = this%elstop(atom_type_start + j+1)
						E_lo = this%En_elstopfile(j)
						E_hi = this%En_elstopfile(j+1)
						Se = Se_lo + (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo)
						exit
					end if
				end do

				vabs = sqrt(vsq)
				if (vabs /= 0.0d0) factor = -Se / vabs

				!! The current forces get modified (reduced) 
				forces(1,i) = forces(1,i) + vel(1,i) * factor
				forces(2,i) = forces(2,i) + vel(2,i) * factor
				forces(3,i) = forces(3,i) + vel(3,i) * factor
			end do
			call mpi_allreduce (MPI_IN_PLACE, forces, 3*Np_all, MPI_DOUBLE_PRECISION, &
						MPI_SUM, MPI_COMM_WORLD, ierr)

		end if
		
		!! Finding the electronic energy loss

		if (to_calculate == 'energy') then

			call iterationRangesToProcesses (1, Np, ntasks, rank, istart, istop)

			SeLoss = 0.0d0
			E_kinetic_atoms = 0.0d0	
			do ki = istart, istop
				i = eel_for_atoms(ki)
				Se = 0.0
				do j = 1, natomtypes
					if (masses(i) == type_mass(j)) then
						itype = j
						atom_type_start = (itype-1)*this%nrows_esdata
						exit
					end if
				end do

				vsq = dot_product (vel(1:3, i), vel(1:3, i))
				energy = 0.5d0 * masses(i) * vsq
				E_kinetic_atoms = E_kinetic_atoms + energy

				if (energy < Ecut) cycle
				if (energy < this%En_elstopfile(1)) cycle
				if (energy > this%En_elstopfile(this%nrows_esdata)) then
					if (rank == 0) write (*,*) "ERROR: Kinetic energy ", energy, &
								"eV of atom is higher than electron stopping data"
					call mpi_finalize(ierr)
					stop
				end if

				!! Find position of atom K.E in the data file and then corresponding electronic stopping
				!! to apply the friction to the current forces

				do j = 1, this%nrows_esdata
					if (energy == this%En_elstopfile(j)) then 
						Se = this%elstop(atom_type_start + j)
						exit
					end if
					if (this%En_elstopfile(j) < energy .and. energy < this%En_elstopfile(j+1)) then
						Se_lo = this%elstop(atom_type_start + j)
						Se_hi = this%elstop(atom_type_start + j+1)
						E_lo = this%En_elstopfile(j)
						E_hi = this%En_elstopfile(j+1)
						Se = Se_lo + (Se_hi - Se_lo) / (E_hi - E_lo) * (energy - E_lo)
						exit
					end if
				end do

				vabs = sqrt(vsq)
				!! roughly, E = E + (dE/dx) * (dx) = (Se) * (vabs*dt)
				SeLoss = SeLoss + (Se * vabs * dt)
			end do

			call mpi_reduce (E_kinetic_atoms, E_kinetic_atoms_temp, 1, MPI_DOUBLE_PRECISION, &
							MPI_SUM, 0,MPI_COMM_WORLD, ierr)
			E_kinetic_atoms = E_kinetic_atoms_temp

			instant_temp_atoms = 2.d0/3.d0/dfloat(Np-1)/boltzconst*E_kinetic_atoms

			call mpi_reduce (SeLoss, SeLoss_temp, 1, MPI_DOUBLE_PRECISION, &
							MPI_SUM, 0,MPI_COMM_WORLD, ierr)
			SeLoss = SeLoss_temp

IF (rank == 0) THEN
			this%cum_EEL = this%cum_EEL + SeLoss

			if (MOD(md_istep, eel_freq_out) == 0 .or. md_istep == num_md_steps) then
				write(100, 1) md_time + this%md_prev_time, SeLoss, this%cum_EEL, &
							E_kinetic_atoms, instant_temp_atoms
				close(unit = 100)
			end if
END IF
		end if

1 FORMAT (E13.6, 4E14.6)

	end if 	!! when md_time > 0.0

end subroutine electron_stopping_velocity_dependent

!! delegate tasks to processes 
subroutine iterationRangesToProcesses (nstart, nstop, nprocs, myid, istart, istop)

	implicit none
	integer, intent(in) :: nstart, nstop, nprocs, myid
	integer, intent(out) :: istart, istop
	integer :: value1, value2
	!! this is for ideal equal division
	value1 = (nstop - nstart + 1)/nprocs
	!! this is for the excess to be distributed 
	value2 = mod( (nstop - nstart + 1), nprocs )

	istart = myid*value1 + nstart + min(myid, value2)
	istop = istart + value1 - 1
	!! ranks lower than the excess value gets one additional task each
	if (value2 > myid) istop = istop + 1

end subroutine iterationRangesToProcesses

end module electronic_stopping
