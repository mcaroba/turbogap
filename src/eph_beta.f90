! ------- option for radiation cascade simulation with electronic stopping 
! based on electron-phonon coupling model.
! This module provides the required electron desities and beta(rho) function
! needed for the calculation of the friction forces according to the eph model.
! The data are extracted from the TDDFT beta file provided by the user. 				
!********* by Uttiyoarnab Saha

!**************************************************************************
module eph_beta

type EPH_Beta_class
	character*128 :: beta_infile
	character*2, dimension(5) :: line
	character*2, allocatable :: element_name(:)
	integer :: n_elements, n_points_rho, n_points_beta
	real*8 :: dr, dr_sq, drho, r_cutoff, rho_cutoff, r_cutoff_sq
	integer, allocatable :: element_number(:)
	real*8, allocatable :: r(:), data_rho(:,:), rho(:), data_beta(:,:)
	integer :: n_species
	
	contains
	procedure :: beta_parameters

end type EPH_Beta_class

contains

subroutine beta_parameters(this,beta_infile,n_species)
	implicit none
	class (EPH_Beta_class) :: this
	integer, intent(in) :: n_species
	integer :: i, j
	character*128, intent(in) :: beta_infile
	
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
	
	this%r_cutoff_sq = this%r_cutoff * this%r_cutoff
	this%rho_cutoff = this%drho * (this%n_points_beta - 1)
	this%dr_sq = this%r_cutoff_sq / (this%n_points_rho - 1)
	
	allocate(this%data_rho(this%n_elements,this%n_points_rho), &
	this%data_beta(this%n_elements,this%n_points_beta))

	! It is assumed that data in the beta file for different atoms is
	! according to the corresponding types of them as specified in input file.
	this%data_rho = 0
	this%data_beta = 0
	do i = 1, this%n_elements
		read(10,*) this%element_number(i)
		do j = 1, this%n_points_rho
			read(10,*) this%data_rho(i,j)
		end do

		do j = 1, this%n_points_beta
			read(10,*) this%data_beta(i,j)
		end do
	end do
	close(unit = 10)
	
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
	
end subroutine beta_parameters

end module eph_beta
