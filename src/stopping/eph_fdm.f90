! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, eph_fdm.f90, is copyright (c) 2023, Uttiyoarnab Saha
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

! ------- option for radiation cascade simulation with electronic stopping 
! based on electron-phonon coupling model.
! This module provides the required data from the FDM input file or the input
! parameters such as the number of grids, their dimensions, T_e, C_e, K_e, etc. 
! that are provided by the user. This data will be used for the electronic
! system in the calculation of electronic stopping according to the eph model.				
!********* by Uttiyoarnab Saha

!**************************************************************************

module eph_fdm
use mpi
type EPH_FDM_class
	character*128 :: FDM_infile
	integer, allocatable :: x_mesh(:), y_mesh(:), z_mesh(:)
	integer :: nx, ny, nz, steps, ntotal, md_last_step
	real*8 :: dx, dy, dz, dV
	real*8 :: x0, x1, y0, y1, z0, z1
	real*8, allocatable :: file_C_e_T(:), file_kappa_e_T(:), file_T_for_kappae(:), &
	file_T_for_Ce(:), C_e_T_deriv(:), K_e_T_deriv(:)
	real*8, allocatable :: T_e(:,:,:), S_e(:,:,:), rho_e(:,:,:), &
	C_e(:,:,:), kappa_e(:,:,:), Q_ei(:,:,:) 
	integer, allocatable :: flag(:,:,:), T_dynamic_flag(:,:,:)
	logical :: T_dependent_parameters = .false., Source_term = .false.
	integer :: Num_C_e, Num_K_e
	
	contains
	procedure :: EPH_FDM_input_params, EPH_FDM_input_file, setMeshIndices
	procedure :: edgesOf3DGrids, feedback_ei_energy, heatDiffusionSolve
	procedure :: whichGrid, Collect_Te, saveOutputToFile, beforeNextFeedback
	procedure :: splineDerivatives, spline_int
	procedure :: EPH_FDM_broadcastQuantities
	
end type EPH_FDM_class

contains

!! The calculation is based on the parameters provided in the input file ....
subroutine EPH_FDM_input_params (this, rank, ierr, md_last_step, &
			in_nx, in_ny, in_nz, in_x0, in_x1, in_y0, in_y1, in_z0, &
			in_z1,in_T_e, in_C_e, in_rho_e, in_kappa_e, in_steps)
	implicit none
	class (EPH_FDM_class) :: this
	integer, intent(in) :: in_nx, in_ny, in_nz, in_steps
	integer, intent(in) :: rank, md_last_step
	real*8, intent(in) :: in_x0, in_x1, in_y0, in_y1, in_z0, in_z1, in_T_e, &
	in_C_e,in_rho_e,in_kappa_e
	integer :: i,j,k,loc, ierr

	this%nx = in_nx; this%ny = in_ny; this%nz = in_nz
	this%x0 = in_x0; this%x1 = in_x1; this%y0 = in_y0; this%y1 = in_y1
	this%z0 = in_z0; this%z1 = in_z1
	
	this%ntotal = in_nx*in_ny*in_nz
	this%steps = in_steps
	this%md_last_step = md_last_step

	allocate(this%x_mesh(this%ntotal), this%y_mesh(this%ntotal), this%z_mesh(this%ntotal))

	call this%edgesOf3DGrids(rank, ierr, in_nx,in_ny,in_nz,in_x0,in_x1,in_y0,in_y1,in_z0,in_z1)
	call this%setMeshIndices()

	allocate(this%T_e(in_nx,in_ny,in_nz), this%S_e(in_nx,in_ny,in_nz), &
	this%rho_e(in_nx,in_ny,in_nz), this%C_e(in_nx,in_ny,in_nz), &
	this%kappa_e(in_nx,in_ny,in_nz), this%flag(in_nx,in_ny,in_nz), &
	this%T_dynamic_flag(in_nx,in_ny,in_nz), this%Q_ei(in_nx,in_ny,in_nz))

	this%S_e = 0.0d0
	this%Q_ei = 0.0d0
	this%T_e = in_T_e
	this%rho_e = in_rho_e
	this%C_e = in_C_e 
	this%kappa_e = in_kappa_e/1000.0d0		!! 1000.0 is for ps <---> fs
	this%flag = 0
	this%T_dynamic_flag = 0

end subroutine EPH_FDM_input_params


!! The calculation is based on the FDM parameters provided through a mesh in parameter file ....
subroutine EPH_FDM_input_file (this, rank, ierr, FDM_infile, md_last_step)
	implicit none
	class (EPH_FDM_class) :: this
	character*128, intent(in) :: FDM_infile
	integer, intent(in) :: rank, md_last_step
	integer :: i, ierr
	real*8 :: T_e_val, S_e_val, rho_e_val, C_e_val, kappa_e_val
	integer :: T_dynamic_flag_val, flag_val
	real*8, parameter :: bignum = 1.1e30
	
	this%FDM_infile = FDM_infile
	this%md_last_step = md_last_step

	open (unit = 10, file = FDM_infile)
	!! First 3 lines are comments		
	read(10,*)
	read(10,*)
	read(10,*)

	!! Grid size is given in next line
	read(10,*) this%nx, this%ny, this%nz, this%steps
	read(10,*) this%x0,this%x1
	read(10,*) this%y0,this%y1
	read(10,*) this%z0,this%z1
	read(10,*)		!! Column headers (i j k T_e S_e rho_e C_e K_e flag T_dyn_flag)

	call this%edgesOf3DGrids(rank, ierr, this%nx,this%ny,this%nz,this%x0,this%x1,this%y0,this%y1, &
	this%z0,this%z1)

	this%ntotal = this%nx*this%ny*this%nz

	!! read grid values
	allocate(this%T_e(this%nx,this%ny,this%nz),this%S_e(this%nx,this%ny,this%nz), &
	this%rho_e(this%nx,this%ny,this%nz), this%C_e(this%nx,this%ny,this%nz), &
	this%kappa_e(this%nx,this%ny,this%nz),this%flag(this%nx,this%ny,this%nz), &
	this%T_dynamic_flag(this%nx,this%ny,this%nz))

	allocate(this%Q_ei(this%nx,this%ny,this%nz)) 	!! for electron-ion energy exchanges  
	this%Q_ei = 0.0d0

	allocate(this%x_mesh(this%ntotal), this%y_mesh(this%ntotal), this%z_mesh(this%ntotal))

	do i = 1, this%ntotal
		read(10,*) this%x_mesh(i), this%y_mesh(i), this%z_mesh(i), &
		T_e_val, S_e_val, rho_e_val, C_e_val, kappa_e_val, &
		flag_val, T_dynamic_flag_val
		
		!! Always give grid indices from 1 to n, not from 0 to n-1
		
		!! bad mesh index
		if ((this%x_mesh(i) < 1).or.(this%x_mesh(i) > this%nx).or.(this%y_mesh(i) < 1).or. &
		(this%y_mesh(i) > this%ny).or.(this%z_mesh(i) < 1).or.(this%z_mesh(i) > this%nz)) then
			if (rank == 0) then
				write(*,*) "ERROR: Index of FDM grid is invalid."
				write(*,*) i, this%x_mesh(i), this%y_mesh(i), this%z_mesh(i)
			end if
			call mpi_finalize(ierr)
			stop
		end if
		!! unphysical parameters
		if (T_e_val<0 .or. rho_e_val<0 .or. C_e_val<0 .or. &
		kappa_e_val<0 .or. flag_val<0 .or. T_dynamic_flag_val<0) then
			if (rank == 0) write(*,*) "ERROR: Negative electronic parameters found."
			call mpi_finalize(ierr)
			stop
		end if

		this%T_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = T_e_val
		this%S_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = S_e_val
		this%rho_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = rho_e_val
		this%C_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = C_e_val 
		this%kappa_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = kappa_e_val/1000.0d0	 !! 1000.0 is for ps <---> fs
		this%flag(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = flag_val
		this%T_dynamic_flag(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = T_dynamic_flag_val
		
		if (T_dynamic_flag_val == 1) this%T_dependent_parameters = .true.
		if (flag_val == 1) this%Source_term = .true.
		
	end do
	close(unit = 10)
	
	!! read T-dependent C_e and K_e from data file
	
	if ( this%T_dependent_parameters ) then
		open (unit = 10, file = 'Te-dependent_e-parameters.txt')
		!! First 3 lines are comments		
		read(10,*)
		read(10,*)
		read(10,*)
		read(10,*) this%Num_C_e
		allocate (this%file_T_for_Ce(this%Num_C_e), this%file_C_e_T(this%Num_C_e))
		do i = 1, this%Num_C_e
			read(10,*) this%file_T_for_Ce(i), this%file_C_e_T(i) 
		end do
		read(10,*) this%Num_K_e
		allocate (this%file_T_for_kappae(this%Num_K_e), this%file_kappa_e_T(this%Num_K_e))
		do i = 1, this%Num_K_e
			read(10,*) this%file_T_for_kappae(i), this%file_kappa_e_T(i)
		end do
		close(unit = 10)

		this%file_kappa_e_T = this%file_kappa_e_T/1000.0d0	 	!! 1000.0 is for ps <---> fs
		
		!! If there is only one or less data point for both C_e(T_e) and K_e(T_e)
		!! and also, if any of these have 0 entries within T-dependent option,
		!! then stop and do with T-independent option.
		
		if (this%Num_C_e <= 1 .and. this%Num_K_e <= 1) then
			if (rank == 0) then
				write(*,*) 'Not enough data for C_e(T_e) and K_e(T_e).'
				write(*,*) 'Try doing with all T_dyn_flag set to 0.'
			end if
			call mpi_finalize(ierr)
			stop
		end if
		if (this%Num_C_e == 0 .or. this%Num_K_e == 0) then
			if (rank == 0) then
				write(*,*) 'Not a single data found for C_e(T_e) or K_e(T_e).'
				write(*,*) 'Try doing with all T_dyn_flag set to 0.'
			end if
			call mpi_finalize(ierr)
			stop
		end if
		
		!! 2nd derivatives for Ce and Ke
		
		allocate(this%C_e_T_deriv(this%Num_C_e), this%K_e_T_deriv(this%Num_K_e))
		this%C_e_T_deriv = 0.0d0
		this%K_e_T_deriv = 0.0d0
		call this%splineDerivatives (this%file_T_for_Ce,this%file_C_e_T,this%Num_C_e,bignum, &
		bignum,this%C_e_T_deriv)
		call this%splineDerivatives (this%file_T_for_kappae,this%file_kappa_e_T,this%Num_K_e, &
		bignum,bignum,this%K_e_T_deriv)
	end if

end subroutine EPH_FDM_input_file


!! broadcast required quantities
subroutine EPH_FDM_broadcastQuantities (this, eph_random_option, eph_fdm_option, ierr)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: eph_random_option, eph_fdm_option, ierr

	if (eph_random_option == 1 .or. eph_fdm_option == 1) then
		call mpi_bcast (this%ntotal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%x0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%x1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%y0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%y1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%z0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%z1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%dy, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%dz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%dV, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%x_mesh, this%ntotal, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%y_mesh, this%ntotal, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%z_mesh, this%ntotal, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		call mpi_bcast (this%T_e, this%ntotal, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	end if
	if (eph_fdm_option == 1) then
		call mpi_bcast (this%Q_ei, this%ntotal, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	end if

end subroutine EPH_FDM_broadcastQuantities

!! set the i,j,k of the coarse mesh
subroutine setMeshIndices(this)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: i, j, k, loc
	do k = 1, this%nz
		do j = 1, this%ny
			do i = 1, this%nx
				loc = i + (j-1)*this%nx + (k-1)*this%nx*this%ny
				this%x_mesh(loc) = i
				this%y_mesh(loc) = j
				this%z_mesh(loc) = k
			end do
		end do
	end do
end subroutine setMeshIndices


!! set the dimensions of the grids in box
subroutine edgesOf3DGrids (this,rank,ierr,nx,ny,nz,x0,x1,y0,y1,z0,z1)
	implicit none
	integer, intent(in) :: rank, nx, ny, nz
	integer :: ierr
	real*8, intent(in) :: x0,x1,y0,y1,z0,z1
	class (EPH_FDM_class) :: this
	if (x0 >= x1 .or. y0 >= y1 .or. z0 >= z1) then
		if (rank == 0) write(*,*) "ERROR: Mesh boundaries are not correct in input or in file"
		call mpi_finalize(ierr)
		stop
	end if
	this%dx = (x1 - x0)/nx
	this%dy = (y1 - y0)/ny
	this%dz = (z1 - z0)/nz
	this%dV = this%dx*this%dy*this%dz
end subroutine edgesOf3DGrids


!! determine the mesh location depending on the position of the atom
integer function whichGrid (this, x, y, z)	result (indx)
	implicit none
	class (EPH_FDM_class) :: this
	real*8, intent(in) :: x, y, z
	integer :: lx, px, ly, py, lz, pz

	lx = floor((x-this%x0) / this%dx)
	px = floor(real(lx / this%nx))
	lx = lx - px * this%nx

	ly = floor((y-this%y0) / this%dy)
	py = floor(real(ly / this%ny))
	ly = ly - py * this%ny

	lz = floor((z-this%z0) / this%dz);
	pz = floor(real(lz / this%nz));
	lz = lz - pz * this%nz

	indx = 1 + lx + ly*this%nx + lz*this%nx*this%ny
end function whichGrid


!! add energy into a cell of mesh
subroutine feedback_ei_energy (this, x, y, z, En, dt)
	implicit none
	class (EPH_FDM_class) :: this
	real*8, intent(in) :: x, y, z, En, dt
	real*8 :: converter, xtrue, ytrue, ztrue
	integer :: indx
	
	xtrue = x; ytrue = y; ztrue = z
	
	!! if points are outside mesh limits, shift them by periodic conditions
	if (x < this%x0) xtrue = this%x1 - abs(x - this%x0)
	if (x > this%x1) xtrue = this%x0 + abs(x - this%x1)
	if (y < this%y0) ytrue = this%y1 - abs(y - this%y0)
	if (y > this%y1) ytrue = this%y0 + abs(y - this%y1)
	if (z < this%z0) ztrue = this%z1 - abs(z - this%z0)
	if (z > this%z1) ztrue = this%z0 + abs(z - this%z1)
	
	indx = this%whichGrid (xtrue, ytrue, ztrue)
	converter = dt * this%dV 
	!! convert to energy per unit time per unit volume
	this%Q_ei(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx)) = &
	this%Q_ei(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx)) + En/converter
end subroutine feedback_ei_energy


!! get temperature of a mesh cell for electrons
real function Collect_Te (this, x, y, z)	result(T_e_indx)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: indx
	real*8, intent(in) :: x, y, z

	indx = this%whichGrid (x, y, z)
	T_e_indx = this%T_e(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx))
end function Collect_Te


!! solve PDE to find the electronic mesh temperatures  
subroutine heatDiffusionSolve(this, dt)
	implicit none
	class (EPH_FDM_class) :: this
	real*8, intent(in) :: dt
	integer :: i,j,k, new_steps, n, indx
	real*8 :: inner_dt, e_sp_heat_min, e_rho_min, e_kappa_max, grad_sq_T, &
	grad_kappa_grad_T, stability, invdx2, invdy2, invdz2, &
	sum_invd, factor, multiply_factor
	integer :: xback, xfront, yback, yfront, zback, zfront

	new_steps = 1
	inner_dt = dt / this%steps
	new_steps = dt / inner_dt
	invdx2 = 1.0d0/this%dx/this%dx
	invdy2 = 1.0d0/this%dy/this%dy
	invdz2 = 1.0d0/this%dz/this%dz
	sum_invd = invdx2 + invdy2 + invdz2
	
	e_sp_heat_min = minval(this%C_e)
	e_rho_min = minval(this%rho_e)
	e_kappa_max = maxval(this%kappa_e)
	
	!! stability criteria (Press, Teukolsky, Vetterling, Flannery, 
	!! Numerical Recipes: The Art of Scientific Computing)	
	
	factor = inner_dt*e_kappa_max*sum_invd / (e_sp_heat_min*e_rho_min)

	if (factor > 0.4d0 .or. factor < 0.0d0) then	!! also works with 0.5
		inner_dt = 0.4d0*inner_dt / factor		!! also works with 0.5
		new_steps = max(int(dt/inner_dt),1)
		inner_dt = dt/new_steps
	end if
	
	!! when parameters are T-dependent and there is no additional source term

	if ( this%T_dependent_parameters .and. (.not. this%Source_term) ) then
		do n = 1, new_steps
			do k = 1, this%nz
				do j = 1, this%ny
					do i = 1, this%nx
						xback = i - 1
						yback = j - 1
						zback = k - 1
						xfront = i + 1
						yfront = j + 1
						zfront = k + 1
						!! check if two sides of x,y,z are within specified mesh limits
						!! otherwise, shift by periodic conditions
						if (xback < 1) xback = this%nx
						if (yback < 1) yback = this%ny
						if (zback < 1) zback = this%nz
						if (xfront > this%nx) xfront = 1
						if (yfront > this%ny) yfront = 1
						if (zfront > this%nz) zfront = 1
						
						grad_sq_T = &
						(this%T_e(k,j,xback) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,j,xfront))*invdx2 + &
						(this%T_e(k,yback,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,yfront,i))*invdy2 + &
						(this%T_e(zback,j,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(zfront,j,i))*invdz2

						!! interpolate for kappa_e(T) which are not on boundaries
						
						if (this%T_e(k,j,xfront) <= this%file_T_for_kappae(1)) then 
							this%kappa_e(k,j,xfront) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,xfront) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,xfront) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,xfront), this%kappa_e(k,j,xfront))
						end if
						
						if (this%T_e(k,j,xback) <= this%file_T_for_kappae(1)) then 
							this%kappa_e(k,j,xback) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,xback) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,xback) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,xback), this%kappa_e(k,j,xback))
						end if
						
						if (this%T_e(k,yfront,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,yfront,i) = this%file_kappa_e_T(1)
						else if (this%T_e(k,yfront,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(k,yfront,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,yfront,i), this%kappa_e(k,yfront,i))
						end if
						
						if (this%T_e(k,yback,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,yback,i) = this%file_kappa_e_T(1)
						else if (this%T_e(k,yback,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(k,yback,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,yback,i), this%kappa_e(k,yback,i))
						end if
						
						if (this%T_e(zfront,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(zfront,j,i) = this%file_kappa_e_T(1)
						else if (this%T_e(zfront,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(zfront,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(zfront,j,i), this%kappa_e(zfront,j,i))
						end if
						
						if (this%T_e(zback,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(zback,j,i) = this%file_kappa_e_T(1)
						else if (this%T_e(zback,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(zback,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(zback,j,i), this%kappa_e(zback,j,i))
						end if

						if (this%T_e(k,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,j,i) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,i), this%kappa_e(k,j,i))
						end if

						!! interpolate for C_e(T) which are not on boundaries
						
						if (this%T_e(k,j,i) <= this%file_T_for_Ce(1)) then
							this%C_e(k,j,i) = this%file_C_e_T(1)
						else if (this%T_e(k,j,i) >= this%file_T_for_Ce(this%Num_C_e)) then
							this%C_e(k,j,i) =  this%file_C_e_T(this%Num_C_e)
						else
							call this%spline_int (this%file_T_for_Ce,this%file_C_e_T, &
						this%C_e_T_deriv, this%Num_C_e, this%T_e(k,j,i), this%C_e(k,j,i))
						end if
						
						multiply_factor = inner_dt/(this%C_e(k,j,i) * this%rho_e(k,j,i))
						
						grad_kappa_grad_T = &
						(this%kappa_e(k,j,xfront) - this%kappa_e(k,j,xback)) * &
						(this%T_e(k,j,xfront) - this%T_e(k,j,xback))*invdx2/4.0 + &
						(this%kappa_e(k,yfront,i) - this%kappa_e(k,yback,i)) * &
						(this%T_e(k,yfront,i) - this%T_e(k,yback,i))*invdy2/4.0 + &
						(this%kappa_e(zfront,j,i) - this%kappa_e(zback,j,i)) * &
						(this%T_e(zfront,j,i) - this%T_e(zback,j,i))*invdz2/4.0
	
						this%T_e(k,j,i) = this%T_e(k,j,i) + multiply_factor * &
								( this%Q_ei(k,j,i) + (this%kappa_e(k,j,i) * grad_sq_T) + &
									grad_kappa_grad_T )
	
						if (this%T_e(k,j,i) < 0.0d0) this%T_e(k,j,i) = 0.0d0
					end do
				end do
			end do
		end do
	end if

	!! when parameters are T-dependent and there is additional source term
					
	if ( this%T_dependent_parameters .and. this%Source_term ) then
		do n = 1, new_steps
			!T_e_prev = this%T_e
			do k = 1, this%nz
				do j = 1, this%ny
					do i = 1, this%nx
						xback = i - 1
						yback = j - 1
						zback = k - 1
						xfront = i + 1
						yfront = j + 1
						zfront = k + 1
						!! check if two sides of x,y,z are within specified mesh limits
						!! otherwise, shift by periodic conditions
						if (xback < 1) xback = this%nx
						if (yback < 1) yback = this%ny
						if (zback < 1) zback = this%nz
						if (xfront > this%nx) xfront = 1
						if (yfront > this%ny) yfront = 1
						if (zfront > this%nz) zfront = 1
						
						grad_sq_T = &
						(this%T_e(k,j,xback) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,j,xfront))*invdx2 + &
						(this%T_e(k,yback,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,yfront,i))*invdy2 + &
						(this%T_e(zback,j,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(zfront,j,i))*invdz2

						!! interpolate for kappa_e(T) which are not on boundaries
						
						if (this%T_e(k,j,xfront) <= this%file_T_for_kappae(1)) then 
							this%kappa_e(k,j,xfront) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,xfront) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,xfront) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,xfront), this%kappa_e(k,j,xfront))
						end if
						
						if (this%T_e(k,j,xback) <= this%file_T_for_kappae(1)) then 
							this%kappa_e(k,j,xback) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,xback) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,xback) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,xback), this%kappa_e(k,j,xback))
						end if
						
						if (this%T_e(k,yfront,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,yfront,i) = this%file_kappa_e_T(1)
						else if (this%T_e(k,yfront,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(k,yfront,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,yfront,i), this%kappa_e(k,yfront,i))
						end if
						
						if (this%T_e(k,yback,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,yback,i) = this%file_kappa_e_T(1)
						else if (this%T_e(k,yback,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(k,yback,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,yback,i), this%kappa_e(k,yback,i))
						end if
						
						if (this%T_e(zfront,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(zfront,j,i) = this%file_kappa_e_T(1)
						else if (this%T_e(zfront,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(zfront,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(zfront,j,i), this%kappa_e(zfront,j,i))
						end if
						
						if (this%T_e(zback,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(zback,j,i) = this%file_kappa_e_T(1)
						else if (this%T_e(zback,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then
							this%kappa_e(zback,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(zback,j,i), this%kappa_e(zback,j,i))
						end if

						if (this%T_e(k,j,i) <= this%file_T_for_kappae(1)) then
							this%kappa_e(k,j,i) =  this%file_kappa_e_T(1)
						else if (this%T_e(k,j,i) >= this%file_T_for_kappae(this%Num_K_e)) then 
							this%kappa_e(k,j,i) =  this%file_kappa_e_T(this%Num_K_e)
						else
							call this%spline_int (this%file_T_for_kappae,this%file_kappa_e_T, &
						this%K_e_T_deriv, this%Num_K_e, this%T_e(k,j,i), this%kappa_e(k,j,i))
						end if

						!! interpolate for C_e(T) which are not on boundaries
						
						if (this%T_e(k,j,i) <= this%file_T_for_Ce(1)) then
							this%C_e(k,j,i) = this%file_C_e_T(1)
						else if (this%T_e(k,j,i) >= this%file_T_for_Ce(this%Num_C_e)) then
							this%C_e(k,j,i) =  this%file_C_e_T(this%Num_C_e)
						else
							call this%spline_int (this%file_T_for_Ce,this%file_C_e_T, &
						this%C_e_T_deriv, this%Num_C_e, this%T_e(k,j,i), this%C_e(k,j,i))
						end if
						
						multiply_factor = inner_dt/(this%C_e(k,j,i) * this%rho_e(k,j,i))
						
						grad_kappa_grad_T = &
						(this%kappa_e(k,j,xfront) - this%kappa_e(k,j,xback)) * &
						(this%T_e(k,j,xfront) - this%T_e(k,j,xback))*invdx2/4.0 + &
						(this%kappa_e(k,yfront,i) - this%kappa_e(k,yback,i)) * &
						(this%T_e(k,yfront,i) - this%T_e(k,yback,i))*invdy2/4.0 + &
						(this%kappa_e(zfront,j,i) - this%kappa_e(zback,j,i)) * &
						(this%T_e(zfront,j,i) - this%T_e(zback,j,i))*invdz2/4.0
	
						this%T_e(k,j,i) = this%T_e(k,j,i) + multiply_factor * &
								( this%Q_ei(k,j,i) + (this%kappa_e(k,j,i) * grad_sq_T) + &
									grad_kappa_grad_T + this%S_e(k,j,i))
	
						if (this%T_e(k,j,i) < 0.0d0) this%T_e(k,j,i) = 0.0d0
					end do
				end do
			end do
		end do
	end if
	
	!! when parameters are not T-dependent and there is additional source term

	if ( (.not. this%T_dependent_parameters) .and. this%Source_term ) then
		do n = 1, new_steps
			do k = 1, this%nz
				do j = 1, this%ny
					do i = 1, this%nx
						xback = i - 1
						yback = j - 1
						zback = k - 1
						xfront = i + 1
						yfront = j + 1
						zfront = k + 1
						!! check if two sides of x,y,z are within specified mesh limits
						!! otherwise, shift by periodic conditions
						if (xback < 1) xback = this%nx
						if (yback < 1) yback = this%ny
						if (zback < 1) zback = this%nz
						if (xfront > this%nx) xfront = 1
						if (yfront > this%ny) yfront = 1
						if (zfront > this%nz) zfront = 1

						grad_sq_T = &
						(this%T_e(k,j,xback) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,j,xfront))*invdx2 + &
						(this%T_e(k,yback,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,yfront,i))*invdy2 + &
						(this%T_e(zback,j,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(zfront,j,i))*invdz2

						multiply_factor = inner_dt/(this%C_e(k,j,i) * this%rho_e(k,j,i))

						this%T_e(k,j,i) = this%T_e(k,j,i) + multiply_factor * &
							( this%Q_ei(k,j,i) + (this%kappa_e(k,j,i) * grad_sq_T) + this%S_e(k,j,i) )

						if (this%T_e(k,j,i) < 0.0d0) this%T_e(k,j,i) = 0.0d0
					end do
				end do
			end do
		end do
	end if

	!! when parameters are not T-dependent and there is no additional source term
	!! this case is most often used	

	if ( (.not. this%T_dependent_parameters) .and. (.not. this%Source_term)) then
		do n = 1, new_steps
			do k = 1, this%nz
				do j = 1, this%ny
					do i = 1, this%nx
						xback = i - 1
						yback = j - 1
						zback = k - 1
						xfront = i + 1
						yfront = j + 1
						zfront = k + 1
						!! check if two sides of x,y,z are within specified mesh limits
						!! otherwise, shift by periodic conditions
						if (xback < 1) xback = this%nx
						if (yback < 1) yback = this%ny
						if (zback < 1) zback = this%nz
						if (xfront > this%nx) xfront = 1
						if (yfront > this%ny) yfront = 1
						if (zfront > this%nz) zfront = 1
						
						grad_sq_T = &
						(this%T_e(k,j,xback) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,j,xfront))*invdx2 + &
						(this%T_e(k,yback,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(k,yfront,i))*invdy2 + &
						(this%T_e(zback,j,i) - 2.0d0*this%T_e(k,j,i) + this%T_e(zfront,j,i))*invdz2

						multiply_factor = inner_dt/(this%C_e(k,j,i) * this%rho_e(k,j,i))

						this%T_e(k,j,i) = this%T_e(k,j,i) + multiply_factor * &
							( this%Q_ei(k,j,i) + (this%kappa_e(k,j,i) * grad_sq_T) )

						if (this%T_e(k,j,i) < 0.0d0) this%T_e(k,j,i) = 0.0d0
					end do
				end do
			end do
		end do
	end if	
	
end subroutine heatDiffusionSolve


!! Save temperature state to ouput file
subroutine saveOutputToFile (this, outputfile, md_istep, dt)
	implicit none
	class (EPH_FDM_class) :: this
	character*128, intent(in) :: outputfile
	integer, intent(in) :: md_istep
	real*8, intent(in) :: dt
	integer :: i,j,k

	open (unit=300, file = outputfile, status = "old", position = "append")
	write(300,*) md_istep + this%md_last_step

	do k = 1, this%nz
		do j = 1, this%ny
			do i = 1, this%nx
				write (300, 1) i, j, k,this%T_e(k,j,i), this%S_e(k,j,i), this%rho_e(k,j,i), &
				this%C_e(k,j,i), this%kappa_e(k,j,i)*1000.0d0, this%flag(k,j,i), this%T_dynamic_flag(k,j,i)
			end do
		end do
	end do
	close(unit = 300)
1 	format(3I6,5e13.6,2I2)
end subroutine saveOutputToFile


!! prepare the holder for energy transfer for next step 
subroutine beforeNextFeedback(this)
	implicit none
	class (EPH_FDM_class) :: this
	this%Q_ei = 0.0d0
end subroutine beforeNextFeedback


!! Find the interpolated value of y corresponding to a given value of x.
!! Data arrays are xarr(1:n) and yarr(1:n). Second derivative of ya is y2arr(1:n).
!! For a given value of x, y is the cubic-spline interpolated value.
subroutine spline_int(this,xarr,yarr,y2arr,n,x,y)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: n
	real*8 :: x, y, xarr(n), yarr(n), y2arr(n)
	integer :: indx, hiindx, loindx
	real*8 :: xwidth, A, B, C, D
	
	!! Interpolate by method of bisection. Find the index limits within which x lies.
	
	loindx = 1
	hiindx = n
	
	do while ((hiindx - loindx) > 1)
		indx = (hiindx + loindx)/2
		if (xarr(indx) > x) then
			hiindx = indx
		else
			loindx = indx
		end if
	end do
	
	!! Cubic spline polynomial value of y for given x.
	!! y = Ay_j + By_j+1 + Cy2_j + Dy2_j+1
	!! A = (x_j+1 - x) / (x_j+1 - x_j)
	!! B = 1 - A
	!! xwidth = x_j+1 - x_j
	!! C = (A³ - A) * xwidth² / 6
	!! D = (B³ - B) * xwidth² / 6
	
	xwidth = xarr(hiindx) - xarr(loindx)
	if (xwidth == 0.0d0) then 
		write(*,*) 'ERROR: The x-values in function for spline interpolation are not distinct.'
		stop
	end if
	
	A = (xarr(hiindx) - x)/xwidth
	B = (x - xarr(loindx))/xwidth
	C = (A**3 - A)*(xwidth**2)/6
	D = (B**3 - B)*(xwidth**2)/6
	
	y = A*yarr(loindx) + B*yarr(hiindx) + C*y2arr(loindx) + D*y2arr(hiindx)

end subroutine spline_int


!! Cubic spline interpolation 
!! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)
!! Find the second derivatives of the interpolating function at tabulated points (x)
!! Find the array of second derivatives of y(x).
subroutine splineDerivatives(this,x,y,n,yp1,ypn,y2)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: n
	real*8 :: yp1, ypn, x(n), y(n), y2(n)
	integer :: i, k
	real*8 :: p, qn, sig, un, u(n)
	
	!! The lower boundary condition is set either to be “natural”, else to have a specified first derivative.
	!! The first derivative is not known, so it is set high value --> natural.
	if (yp1 > 0.99e30) then
		y2(1) = 0.0d0
		u(1) = 0.0d0
	else
		y2(1) = -0.5d0
		u(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	end if
	
	!! This is the decomposition loop of the tridiagonal
	!! algorithm. y2 and u are used for temporary
	!! storage of the decomposed factors.
	do i = 2, n-1
		sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
		p = sig*y2(i-1) + 2.0d0
		y2(i) = (sig-1.0d0)/p
		u(i) = (6.0d0*(( y(i+1)-y(i) ) / (x(i+1)-x(i))-(y(i)-y(i-1)) / &
		( x(i) - x(i-1) ))/( x(i+1) - x(i-1) )-sig*u(i-1)) / p
	end do
	
	!! The upper boundary condition is set either to be “natural”, else to have a specified first derivative.
	!! The first derivative is not known, so it is set high value --> natural.
	if (ypn > 0.99e30) then 
		qn = 0.0d0
		un = 0.0d0
	else
		qn = 0.5d0
		un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	end if
	y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0d0)
	
	!! This is the backsubstitution loop of the tridiagonal algorithm.
	do k = n-1,1,-1
		y2(k) = y2(k)*y2(k+1) + u(k)
	end do

end subroutine splineDerivatives

end module eph_fdm
