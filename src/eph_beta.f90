! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, eph_beta.f90, is copyright (c) 2023, Uttiyoarnab Saha
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
! This module provides the required electron desities and beta(rho) function
! needed for the calculation of the friction forces according to the eph model.
! The data are extracted from the TD-DFT beta file provided by the user. 				
!********* by Uttiyoarnab Saha

!**************************************************************************
module eph_beta

type EPH_Beta_class
	character*128 :: beta_infile
	character*2, dimension(5) :: line
	character*2, allocatable :: element_name(:)
	integer :: n_elements, n_points_rho, n_points_beta
	real*8 :: dr, drho, r_cutoff, rho_cutoff
	integer, allocatable :: element_number(:)
	real*8, allocatable :: r(:), data_rho(:,:), rho(:), data_beta(:,:), data_alpha(:,:)
	real*8, allocatable :: y2rho(:,:), y2beta(:,:), y2alpha(:,:), y2rho_r_sq(:,:)
	integer :: n_species
	
	! --------------------------------------------------------------------
	! ** As fix user-eph does in LAMMPS, but not as shown in the papers **
	real*8 :: dr_sq, r_cutoff_sq
	real*8, allocatable :: data_rho_rsq(:,:), r_sq(:)
	! --------------------------------------------------------------------
	
	contains
	procedure :: beta_parameters, spline_int, spline

end type EPH_Beta_class

contains

! Get the rho, beta and alpha data from .beta file 
subroutine beta_parameters(this,beta_infile,n_species)
	implicit none
	class (EPH_Beta_class) :: this
	integer, intent(in) :: n_species
	integer :: i, j
	character*128, intent(in) :: beta_infile
	real*8, allocatable :: y2(:), z2(:), w2(:), y2r2(:)
	real*8, parameter :: bignum = 1.1e30
	
	this%n_species = n_species
	this%beta_infile = beta_infile
	
	open (unit = 10, file = beta_infile)
	! First 3 lines are comments		
	read(10,*)
	read(10,*)
	read(10,*)
	read(10,*) (this%line(i), i = 1, (this%n_species+1))
	read(this%line(1),'(I2)') this%n_elements
	if (this%n_elements <= 0) then
		write(*,*) "ERROR: Negative or 0 elements in the density file"
	stop
	end if

	allocate(this%element_name(this%n_elements), this%element_number(this%n_elements))

	! read the number of elements and their names
	do i = 1, this%n_elements
		this%element_name(i) = this%line(i+1)
	end do
	
	! Read r, rho and beta parameters from the beta file
	read(10,*) this%n_points_rho, this%dr, this%n_points_beta, this%drho, &
	this%r_cutoff
	
	this%rho_cutoff = this%drho * (this%n_points_beta - 1)
	
	allocate(this%data_rho(this%n_elements,this%n_points_rho), &
		this%data_beta(this%n_elements,this%n_points_beta), &
			this%data_alpha(this%n_elements,this%n_points_beta))

	! It is assumed that data in the beta file for different elements is
	! according to the corresponding types of the species as specified in 
	! the input file of TurboGAP.
	
	write(*,*) ' -- MESSAGE --'
	write(*,*) 'It is assumed that data in the .beta file for different elements is'
	write(*,*) 'according to the corresponding types of the species as specified in'
	write(*,*) 'the input file, i.e the order of elements in both files are same.'
	write(*,*) 'Please check it.'
	
	this%data_rho = 0
	this%data_beta = 0
	do i = 1, this%n_elements
		read(10,*) this%element_number(i)
		do j = 1, this%n_points_rho
			read(10,*) this%data_rho(i,j)
		end do

		do j = 1, this%n_points_beta
			read(10,*) this%data_beta(i,j)
			this%data_beta(i,j) = this%data_beta(i,j)		!*1000.0
		end do
	end do
	close(unit = 10)

	! find the values of alpha from beta
	do i = 1, this%n_elements
		do j = 1, this%n_points_beta
			this%data_alpha(i,j) = sqrt(this%data_beta(i,j))
		end do
	end do
	
	! create the array with values of r for rho vs. r from data in beta file
	if (.not. allocated(this%r)) allocate(this%r(this%n_points_rho))
	this%r(1) = 0.0
	do i = 2, this%n_points_rho
		this%r(i) = this%r(i-1) + this%dr
	end do
	
	! create the array with values of rho for beta vs. rho from data in beta file
	if (.not. allocated(this%rho)) allocate(this%rho(this%n_points_beta))
	this%rho(1) = 0.0
	do i = 2, this%n_points_beta
		this%rho(i) = this%rho(i-1) + this%drho
	end do
	
	! Have the y"(x) in y2rho and y2beta for cubic spline interpolation 
	! for all types of atoms and use them as and when needed afterwards
	
	allocate(this%y2rho(this%n_elements,this%n_points_rho), y2(this%n_points_rho))
	
	allocate(this%y2beta(this%n_elements,this%n_points_beta), z2(this%n_points_beta), &
				this%y2alpha(this%n_elements,this%n_points_beta), w2(this%n_points_beta))

	do i = 1, this%n_elements
		y2 = 0
		call this%spline (this%r,this%data_rho(i,:),this%n_points_rho,bignum,bignum,y2)
		this%y2rho(i,:) = y2(:)
		z2 = 0
		call this%spline (this%rho,this%data_beta(i,:),this%n_points_beta,bignum,bignum,z2)
		this%y2beta(i,:) = z2(:)
		w2 = 0
		call this%spline (this%rho, this%data_alpha(i,:),this%n_points_beta,bignum,bignum,w2)
		this%y2alpha(i,:) = w2(:)
	end do
	
	
	! ------------------------------------------------------------------------
	! ** As fix user-eph does with the function rho(r²) versus r² **
	
	this%r_cutoff_sq = this%r_cutoff * this%r_cutoff
	this%dr_sq = this%r_cutoff_sq / (this%n_points_rho - 1)
	
	! create the array with values of r²
	if (.not. allocated(this%r_sq)) allocate(this%r_sq(this%n_points_rho))
	this%r_sq(1) = 0.0
	do i = 2, this%n_points_rho
		this%r_sq(i) = this%r_sq(i-1) + this%dr_sq
	end do
	
	if (.not. allocated(this%data_rho_rsq)) &
		allocate(this%data_rho_rsq(this%n_elements,this%n_points_rho))
	
	do i = 1, this%n_elements
		do j = 1, this%n_points_rho
			call this%spline_int (this%r, this%data_rho(i,:), this%y2rho(i,:), &
					this%n_points_rho, sqrt((j-1)*this%dr_sq), this%data_rho_rsq(i,j))
		end do
	end do
	
	
	! keep y"(x) for data_rho_rsq
	
	allocate(this%y2rho_r_sq(this%n_elements,this%n_points_rho), y2r2(this%n_points_rho))
	
	do i = 1, this%n_elements
		y2r2 = 0
		call this%spline (this%r_sq,this%data_rho_rsq(i,:),this%n_points_rho,bignum,bignum,y2r2)
		this%y2rho_r_sq(i,:) = y2r2(:)
	end do
	! ------------------------------------------------------------------------
	
	
	! -------------------------------------------------
	! To check rho-r data that is read
	
	!open(unit = 300, file = "rho-data.txt")
	!
	!write(300,*) 'rho versus r'
	!do i = 1, this%n_points_rho
	!	if (this%r(i) < this%r_cutoff) write(300,*) this%r(i), '', this%data_rho(1,i) 
	!end do
	!write(300,*) '======================================='
	!write(300,*) 'rho_r_sq versus r_sq'
	!do i = 1, this%n_points_rho
	!	if (this%r_sq(i) < this%r_cutoff_sq) write(300,*) this%r_sq(i), '', this%data_rho_rsq(1,i)
	!end do
	!close(unit = 300)

	! -------------------------------------------------
	
	
end subroutine beta_parameters

! Find the interpolated value of y corresponding to a given value of x.
! Data arrays are xarr(1:n) and yarr(1:n). Second derivative of ya is y2arr(1:n).
! For a given value of x, y is the cubic-spline interpolated value.
subroutine spline_int(this,xarr,yarr,y2arr,n,x,y)
	implicit none
	class (EPH_Beta_class) :: this
	integer :: n
	real*8 :: x, y, xarr(n), yarr(n), y2arr(n)
	integer :: indx, hiindx, loindx
	real*8 :: xwidth, A, B, C, D
	
	! Interpolate by method of bisection. Find the index limits within which x lies.
	
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
	
	! Cubic spline polynomial value of y for given x.
	! y = Ay_j + By_j+1 + Cy2_j + Dy2_j+1
	! A = (x_j+1 - x) / (x_j+1 - x_j)
	! B = 1 - A
	! xwidth = x_j+1 - x_j
	! C = (A³ - A) * xwidth² / 6
	! D = (B³ - B) * xwidth² / 6
	
	xwidth = xarr(hiindx) - xarr(loindx)
	if (xwidth == 0.0) then 
		write(*,*) 'ERROR: The x-values in function for spline interpolation are not distinct.'
		stop
	end if
	
	A = (xarr(hiindx) - x)/xwidth
	B = (x - xarr(loindx))/xwidth
	C = (A**3 - A)*(xwidth**2)/6
	D = (B**3 - B)*(xwidth**2)/6
	
	y = A*yarr(loindx) + B*yarr(hiindx) + C*y2arr(loindx) + D*y2arr(hiindx)

end subroutine spline_int


! Cubic spline interpolation 
! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)
! Find the second derivatives of the interpolating function at tabulated points (x)
! Find the array of second derivatives of y(x).
subroutine spline(this,x,y,n,yp1,ypn,y2)
	implicit none
	class (EPH_Beta_class) :: this
	integer :: n
	real*8 :: yp1, ypn, x(n), y(n), y2(n)
	integer :: i, k
	real*8 :: p, qn, sig, un, u(n)
	
	! The lower boundary condition is set either to be “natural”, else to have a specified first derivative.
	! The first derivative is not known, so it is set high value --> natural.
	if (yp1 > 0.99e30) then
		y2(1) = 0.0
		u(1) = 0.0
	else
		y2(1) = -0.5
		u(1) = (3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	end if
	
	! This is the decomposition loop of the tridiagonal
	! algorithm. y2 and u are used for temporary
	! storage of the decomposed factors.
	do i = 2, n-1
		sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
		p = sig*y2(i-1) + 2.0
		y2(i) = (sig-1.0)/p
		u(i) = (6.0*(( y(i+1)-y(i) ) / (x(i+1)-x(i))-(y(i)-y(i-1)) / &
		( x(i) - x(i-1) ))/( x(i+1) - x(i-1) )-sig*u(i-1)) / p
	end do
	
	! The upper boundary condition is set either to be “natural”, else to have a specified first derivative.
	! The first derivative is not known, so it is set high value --> natural.
	if (ypn > 0.99e30) then 
		qn = 0.0
		un = 0.0
	else
		qn = 0.5
		un = (3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	end if
	y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0)
	
	!This is the backsubstitution loop of the tridiagonal algorithm.
	do k = n-1,1,-1
		y2(k) = y2(k)*y2(k+1) + u(k)
	end do

end subroutine spline

end module eph_beta
