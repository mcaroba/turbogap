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
! that are provided by the user. This data will be used for the calculation 
! of electronic stopping according to eph model.				
!********* by Uttiyoarnab Saha

!**************************************************************************

module eph_fdm

type EPH_FDM_class
	character*128 :: FDM_infile, TDP_infile
	integer, allocatable :: x_mesh(:), y_mesh(:), z_mesh(:)
	integer :: nx, ny, nz, steps, ntotal, fdm_option, md_last_step
	real*8 :: dx, dy, dz, dV
	real*8 :: dT, C_e_T, kappa_e_T, E_e_T
	real*8 :: x0, x1, y0, y1, z0, z1
	real*8, allocatable :: file_C_e_T(:), file_kappa_e_T(:)
	real*8, allocatable :: T_e(:,:,:), S_e(:,:,:), rho_e(:,:,:), &
	C_e(:,:,:), kappa_e(:,:,:), Q_ei(:,:,:) 
	integer, allocatable :: flag(:,:,:), T_dynamic_flag(:,:,:)
	
	contains
	procedure :: EPH_FDM_input_params, EPH_FDM_input_file, set_mesh_indices
	procedure :: edges_of_3D_grids, feedback_ei_energy, heat_diffusion_solve
	procedure :: which_grid, Collect_Te, save_output_to_file, beforeNextFeedback 
	
end type EPH_FDM_class

contains

! The calculation is based on the parameters provided in the input file ....
subroutine EPH_FDM_input_params (this,fdm_option,in_nx,in_ny,in_nz,in_x0,in_x1,in_y0, &
		in_y1, in_z0, in_z1,in_T_e, in_C_e, in_rho_e, in_kappa_e, in_steps)
	implicit none
	class (EPH_FDM_class) :: this
	integer, intent(in) :: in_nx, in_ny, in_nz, in_steps
	integer, intent(in) :: fdm_option
	real*8, intent(in) :: in_x0,in_x1,in_y0,in_y1,in_z0,in_z1,in_T_e,in_C_e,in_rho_e,in_kappa_e
	integer :: i,j,k,loc

	this%md_last_step = 0
	this%fdm_option = fdm_option
	this%nx = in_nx; this%ny = in_ny; this%nz = in_nz
	this%x0 = in_x0; this%x1 = in_x1; this%y0 = in_y0; this%y1 = in_y1
	this%z0 = in_z0; this%z1 = in_z1
	
	this%ntotal = in_nx*in_ny*in_nz
	this%steps = in_steps
	
	allocate(this%x_mesh(this%ntotal), this%y_mesh(this%ntotal), this%z_mesh(this%ntotal))
	
	call this%edges_of_3D_grids(in_nx,in_ny,in_nz,in_x0,in_x1,in_y0,in_y1,in_z0,in_z1)
	
	allocate(this%T_e(in_nx,in_ny,in_nz),this%S_e(in_nx,in_ny,in_nz),this%rho_e(in_nx,in_ny,in_nz), &
	this%C_e(in_nx,in_ny,in_nz), this%kappa_e(in_nx,in_ny,in_nz),this%flag(in_nx,in_ny,in_nz), &
	this%T_dynamic_flag(in_nx,in_ny,in_nz))
	allocate(this%Q_ei(in_nx,in_ny,in_nz)) ! for electron-ion energy exchanges
	
	call this%set_mesh_indices()
	
	this%Q_ei = 0.0
	this%T_e = in_T_e 
	this%rho_e = in_rho_e
	this%C_e = in_C_e 
	this%kappa_e = in_kappa_e	!/1000.0		! 1000.0 is for ps <---> fs
	this%flag = 0
	this%T_dynamic_flag = 0
end subroutine EPH_FDM_input_params


! The calculation is based on the FDM parameters provided through a mesh in parameter file ....
subroutine EPH_FDM_input_file (this, fdm_option, FDM_infile, md_last_step, TDP_infile)
	implicit none
	class (EPH_FDM_class) :: this
	character*128, intent(in) :: FDM_infile
	character*128, optional, intent(in) :: TDP_infile
	integer, intent(in) :: fdm_option, md_last_step
	integer :: i
	! some T-dependent parameters will be used later ....
	!real*8 :: dT, C_e_T, kappa_e_T, E_e_T
	!real*8, allocatable :: file_C_e_T(:), file_kappa_e_T(:)
	real*8 :: T_e_val, S_e_val, rho_e_val, C_e_val, kappa_e_val, flag_val
	integer :: T_dynamic_flag_val, iostatus
	
	this%fdm_option = fdm_option
	this%FDM_infile = FDM_infile
	this%md_last_step = md_last_step
	
	open (unit = 10, file = FDM_infile)
	! First 3 lines are comments		
	read(10,*)
	read(10,*)
	read(10,*)
	      
	! Grid size is given in next line
	read(10,*) this%nx, this%ny, this%nz, this%steps
	read(10,*) this%x0,this%x1
	read(10,*) this%y0,this%y1
	read(10,*) this%z0,this%z1
	read(10,*)	! Column headers (i j k T_e S_e rho_e C_e K_e flag T_dyn_flag)
	call this%edges_of_3D_grids(this%nx,this%ny,this%nz,this%x0,this%x1,this%y0,this%y1, &
	this%z0,this%z1)
	
	this%ntotal = this%nx*this%ny*this%nz
	
	! read grid values
	allocate(this%T_e(this%nx,this%ny,this%nz),this%S_e(this%nx,this%ny,this%nz), &
	this%rho_e(this%nx,this%ny,this%nz), this%C_e(this%nx,this%ny,this%nz), &
	this%kappa_e(this%nx,this%ny,this%nz),this%flag(this%nx,this%ny,this%nz), &
	this%T_dynamic_flag(this%nx,this%ny,this%nz))
	
	allocate(this%Q_ei(this%nx,this%ny,this%nz)) ! for electron-ion energy exchanges  
	
	allocate(this%x_mesh(this%ntotal), this%y_mesh(this%ntotal), this%z_mesh(this%ntotal))

	do i = 1, this%ntotal
		read(10,*) this%x_mesh(i), this%y_mesh(i), this%z_mesh(i), &
		T_e_val, S_e_val, rho_e_val, C_e_val, kappa_e_val, &
		flag_val, T_dynamic_flag_val
		
		! Always give grid indices from 1 to n, not from 0 to n-1
		
		! bad mesh index
		if ((this%x_mesh(i) < 1).or.(this%x_mesh(i) > this%nx).or.(this%y_mesh(i) < 1).or. &
		(this%y_mesh(i) > this%ny).or.(this%z_mesh(i) < 1).or.(this%z_mesh(i) > this%nz)) then
			write(*,*) "ERROR: Index of FDM grid is invalid."
			write(*,*) i, this%x_mesh(i), this%y_mesh(i), this%z_mesh(i)
			stop
		end if
		! unphysical parameters
		if (T_e_val<0 .or. S_e_val<0 .or. rho_e_val<0 .or. C_e_val<0 .or. &
		kappa_e_val<0 .or. flag_val<0 .or. T_dynamic_flag_val<0) then
			write(*,*) "ERROR: Negative electronic parameters found."
			stop
		end if
		
		this%T_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = T_e_val
		this%S_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = S_e_val
		this%rho_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = rho_e_val
		this%C_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = C_e_val 
		this%kappa_e(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = kappa_e_val	!/1000.0	! 1000.0 is for ps <---> fs
		this%flag(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = flag_val
		this%T_dynamic_flag(this%z_mesh(i),this%y_mesh(i),this%x_mesh(i)) = T_dynamic_flag_val
	end do
	close(unit = 10)

end subroutine EPH_FDM_input_file


! set the i,j,k of the coarse mesh
subroutine set_mesh_indices(this)
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
end subroutine set_mesh_indices


! set the dimensions of the grids in box
subroutine edges_of_3D_grids (this,nx,ny,nz,x0,x1,y0,y1,z0,z1)
	implicit none
	integer, intent(in) :: nx, ny, nz
	real*8, intent(in) :: x0,x1,y0,y1,z0,z1
	class (EPH_FDM_class) :: this
	if (x0 >= x1 .or. y0 >= y1 .or. z0 >= z1) then
		write(*,*) "ERROR: Box dimensions are not correct in input or in file"
		stop
	end if
	this%dx = (x1 - x0)/nx
	this%dy = (y1 - y0)/ny
	this%dz = (z1 - z0)/nz
	this%dV = this%dx*this%dy*this%dz
end subroutine edges_of_3D_grids


! determine the mesh location depending on the position of the atom
integer function which_grid (this, x, y, z)	result (indx)
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
end function which_grid


! add energy into a cell of mesh
subroutine feedback_ei_energy (this, x, y, z, En, dt)
	implicit none
	class (EPH_FDM_class) :: this
	real*8, intent(in) :: x, y, z, En, dt
	real*8 :: converter
	integer :: indx
	indx = this%which_grid (x, y, z)
	converter = dt * this%dV 
	! convert energy to power
	this%Q_ei(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx)) = &
	this%Q_ei(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx)) + En/converter
end subroutine feedback_ei_energy


! get temperature of a mesh cell for electrons
real function Collect_Te (this, x, y, z)	result(T_e_indx)
	implicit none
	class (EPH_FDM_class) :: this
	integer :: indx
	real*8, intent(in) :: x, y, z

	indx = this%which_grid (x, y, z)
	T_e_indx = this%T_e(this%z_mesh(indx),this%y_mesh(indx),this%x_mesh(indx))
end function Collect_Te


! solve PDE to find the electronic mesh temperatures  
subroutine heat_diffusion_solve(this,dt)
	implicit none
	class (EPH_FDM_class) :: this
	real*8, intent(in) :: dt
	real*8 :: inner_dt, e_sp_heat_min,e_rho_min,e_kappa_max,E_e
	integer :: i,j,k, new_steps, n, indx
	real*8 :: stability, invdx2, invdy2, invdz2, sum_invd, factor
	real*8, allocatable :: T_e_prev(:,:,:)
	integer :: xback, xfront, yback, yfront, zback, zfront
	
	new_steps = 1
	inner_dt = dt / this%steps
	new_steps = dt / inner_dt
	invdx2 = 1.0/this%dx/this%dx
	invdy2 = 1.0/this%dy/this%dy
	invdz2 = 1.0/this%dz/this%dz
	sum_invd = invdx2 + invdy2 + invdz2
	
	e_sp_heat_min = minval(this%C_e)
	e_rho_min = minval(this%rho_e)
	e_kappa_max = maxval(this%kappa_e)
	
	! stability criteria (Press, Teukolsky, Vetterling, Flannery, 
	! Numerical Recipes: The Art of Scientific Computing)	
	
	factor = inner_dt*e_kappa_max*sum_invd / (e_sp_heat_min*e_rho_min)
	
	!stability = 1.0 - 2.0*factor
	
	!if (stability < 0.0) then
	!	inner_dt = 0.5*e_kappa_max*sum_invd / (e_sp_heat_min*e_rho_min)
	!	new_steps = int(dt/inner_dt) + 1
	!	inner_dt = dt/new_steps
	!end if

	if (factor > 0.4 .or. factor < 0.0) then	! also works with 0.5, to see in details
		inner_dt = 0.4*inner_dt / factor		! also works with 0.5. to see in details
		new_steps = max(int(dt/inner_dt),1)
		inner_dt = dt/new_steps
	end if
	
	allocate(T_e_prev(this%nx,this%ny,this%nz))
	
	do n = 1, new_steps
		T_e_prev = this%T_e
		do k = 1, this%nz
			do j = 1, this%ny
				do i = 1, this%nx
					xback = i - 1
					yback = j - 1
					zback = k - 1
					xfront = i + 1
					yfront = j + 1
					zfront = k + 1
					! check if two sides of x,y,z are within specified mesh limits
					if (xback < 1) xback = this%nx
					if (yback < 1) yback = this%ny
					if (zback < 1) zback = this%nz
					if (xfront > this%nx) xfront = 1
					if (yfront > this%ny) yfront = 1
					if (zfront > this%nz) zfront = 1
					
					this%T_e(k,j,i) = T_e_prev(k,j,i) + &
					inner_dt/(this%C_e(k,j,i) * this%rho_e(k,j,i)) * &
					((this%kappa_e(k,j,i) * &
					((T_e_prev(k,j,xback) - 2.0*T_e_prev(k,j,i) + T_e_prev(k,j,xfront))*invdx2 + &
					(T_e_prev(k,yback,i) - 2.0*T_e_prev(k,j,i) + T_e_prev(k,yfront,i))*invdy2 + &
					(T_e_prev(zback,j,i) - 2.0*T_e_prev(k,j,i) + T_e_prev(zfront,j,i))*invdz2)) + &
					this%Q_ei(k,j,i))
							
					if (this%T_e(k,j,i) < 0.0) this%T_e(k,j,i) = 0.0
				end do
			end do
		end do
	end do
end subroutine heat_diffusion_solve


! Save temperature state to ouput file
subroutine save_output_to_file (this, outputfile, md_istep, dt)
	implicit none
	class (EPH_FDM_class) :: this
	character*128, intent(in) :: outputfile
	integer, intent(in) :: md_istep
	real*8, intent(in) :: dt
	integer :: i,j,k

	if( md_istep == 0 .or. md_istep == -1 )then
		open (unit=300, file = outputfile, status = "unknown")
		write(300,*) this%nx, ' ', this%ny, ' ', this%nz
		write(300,*) 'i j k T_e S_e rho_e C_e K_e flag T_dyn_flag'		! Column headers
		! write(300,*) 'i  j  k  T  Q T_average Q_average Q_total time-step' 
		write(300,*) md_istep + this%md_last_step
		else
			open (unit=300, file = outputfile, status = "old", position = "append")
			write(300,*) md_istep + this%md_last_step
	end if

	do k = 1, this%nz
		do j = 1, this%ny
			do i = 1, this%nx
				!write (300, "(I4,I4,I4,E13.6,E14.6)") i, j, k, this%T_e(k,j,i), -this%Q_ei(k,j,i)
				write (300, *) i, j, k,this%T_e(k,j,i), this%S_e(k,j,i), this%rho_e(k,j,i), &
				this%C_e(k,j,i), this%kappa_e(k,j,i), this%flag(k,j,i), this%T_dynamic_flag(k,j,i)
			end do
		end do
	end do
	!write(300, "(38x,E13.6,2E14.6,E13.6)") sum(this%T_e)/this%ntotal, -sum(this%Q_ei)/this%ntotal, -sum(this%Q_ei), dt
	close(unit = 300)
end subroutine save_output_to_file


! prepare the holder for energy transfer for next step 
subroutine beforeNextFeedback(this)
	implicit none
	class (EPH_FDM_class) :: this
	this%Q_ei = 0.0
end subroutine beforeNextFeedback

end module eph_fdm
