! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, eph_electronic_stopping.f90, is copyright (c) 2023, Uttiyoarnab Saha
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
! based on electron-phonon coupling model				
! The model for electron-phonon coupling for calculation of electronic stopping
! according to the Langevin dynamics with spatial correlations is implemented here.
! The reference papers are PRL 120, 185501 (2018); PRB 99, 174302 (2019);
! PRB 99, 174301 (2019).
!********* by Uttiyoarnab Saha

!**************************************************************************

module eph_electronic_stopping

use eph_fdm
use eph_beta
USE OMP_LIB

type EPH_LangevinSpatialCorrelation_class
	
	character*128 :: T_outfile
	integer :: T_out_freq, isfriction, israndom, model_eph, T_out_mesh_freq, &
			fdm_option
	real*8, allocatable :: forces_fric(:,:), forces_rnd(:,:)
	real*8 :: E_eph_cumulative, md_prev_time
	
	contains
	procedure :: eph_InitialValues	
	procedure :: eph_LangevinForces, eph_LangevinEnergyDissipation
	
end type EPH_LangevinSpatialCorrelation_class

contains


!! Initialize some required variables
subroutine eph_InitialValues (this, isfriction, israndom, fdm_option, T_outfile, &
			T_out_freq, T_out_mesh_freq, model_eph, E_eph_cumulative_prev, md_prev_time)
	implicit none
	
	class (EPH_LangevinSpatialCorrelation_class) :: this
	character*128, intent(in) :: T_outfile
	integer :: T_out_freq, isfriction, israndom, fdm_option, model_eph, T_out_mesh_freq
	real*8 :: E_eph_cumulative_prev, md_prev_time
	
	this%T_outfile = T_outfile
	this%T_out_freq = T_out_freq
	this%isfriction = isfriction 
	this%israndom = israndom
	this%fdm_option = fdm_option
	this%model_eph = model_eph
	this%T_out_mesh_freq = T_out_mesh_freq
	this%E_eph_cumulative = E_eph_cumulative_prev + 0.0d0	!! cumulative energy transferred
	this%md_prev_time = md_prev_time		!! MD time from previous run
	
end subroutine eph_InitialValues


!! Calculate the friction and random forces according to the Langevin spatial correlations
!! and then correct the forces obtained using the interatomic potential. 
subroutine eph_LangevinForces (this, vel, forces, masses, type_mass, &
			md_istep, dt, md_time, positions, natomtypes, eph_for_atoms, beta, fdm)
	implicit none
	
	class (eph_LangevinSpatialCorrelation_class) :: this
	
	type (EPH_Beta_class), intent(in) :: beta
	type (EPH_FDM_class), intent(in) :: fdm
	
	integer, intent(in) :: natomtypes, md_istep, eph_for_atoms(:)
	real*8, intent(in) :: vel(:,:), positions(:,:), masses(:), &
	type_mass(:), dt, md_time
	
	real*8, intent(inout) :: forces(:,:)
	
	integer :: Np, i, j, ki, kj, itype, jtype, atom_type
	
	real*8 :: xi, yi, zi, xj, yj, zj, r_ij, alpha_I, alpha_J, rho_ij, v_Te, rel_ij(3), &
	rel_v_ij(3), r_ij_sq, correl_factor_eta, multiply_factor, multiply_factor1, &
	multiply_factor2

	real*8, allocatable :: rand_vec(:,:), rho_I(:), w_I(:,:)
	real*8, parameter :: boltzconst = 8.61733326E-05

	!! Do all calculations for MD time greater than 0
	
	if ( .not. (md_istep == 0 .or. md_istep == -1) ) then

	Np = size(eph_for_atoms)	!! Number of atoms in specified group
	!! The size of arrays should be same as the total number of atoms because
	!! the specified group of atoms contain their ids.

	if (.not. allocated(this%forces_fric)) allocate(this%forces_fric(3,size(vel,2)))
	if (.not. allocated(this%forces_rnd)) allocate(this%forces_rnd(3,size(vel,2)))

	this%forces_fric = 0.0d0
	this%forces_rnd = 0.0d0
	
	!! --------------------------------------------
	!! 					 Definitions 
	!! --------------------------------------------
	!! *-*-*-*- -*-*-*-*
	!! Implementation of random forces and friction forces 
	!! according to PRL 120, 185501 (2018)
	!!		----------------------------	
	
	!! ---------------------------------------------------
	!! To find forces sigma_I and eta_I ......
	!! ---------------------------------------------------

	if (this%isfriction == 1 .or. this%israndom == 1) then

		! Find total atomic electronic densities at all sites limited by the 
		! value of r_cutoff given in the beta file.

		allocate(rho_I(size(vel,2)))
		rho_I = 0.0d0

!$OMP PARALLEL
!$OMP DO
		do ki = 1, Np
			i = eph_for_atoms(ki)
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			
			do kj = 1, Np
				j = eph_for_atoms(kj)
				if (i /= j) then
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi, xj,yj,zj)
					if (r_ij < beta%r_cutoff) then
						call getAtomType(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						rho_ij = 0.0d0
						call beta%spline_int (beta%r,beta%data_rho(jtype,:), &
						beta%y2rho(jtype,:),beta%n_points_rho, r_ij, rho_ij)

						rho_I(i) = rho_I(i) + rho_ij
					end if
				end if
			end do
		end do
!$OMP END DO
!$OMP END PARALLEL

		!! -------------------------------------------------
		!! When random forces need to be calculated
		!! -------------------------------------------------
		
		if (this%israndom == 1) then
			!! generate uncorrelated random vectors
			allocate(rand_vec(3,size(vel,2)))
			call randomGaussianArray(size(vel,2), 0.0d0, 1.0d0, rand_vec(1,:))
			call randomGaussianArray(size(vel,2), 0.0d0, 1.0d0, rand_vec(2,:))
			call randomGaussianArray(size(vel,2), 0.0d0, 1.0d0, rand_vec(3,:))
			correl_factor_eta = sqrt(2.0d0*boltzconst*1000.0d0/dt)	!! time units fs ----> ps
!$OMP PARALLEL
!$OMP DO
			do ki = 1, Np
				i = eph_for_atoms(ki)
				call getAtomType(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)

				do kj = 1, Np
					j = eph_for_atoms(kj)
					if (i /= j) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						
						r_ij_sq = r_ij*r_ij

						if (r_ij < beta%r_cutoff) then
							call relativeVector(positions(:,j), positions(:,i), rel_ij)

							call getAtomType(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							!! find rho_J
							rho_ij = 0.0d0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:),beta%y2rho(jtype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							!! find alpha_I	
							alpha_I = 0.0d0
							call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
									beta%y2alpha(itype,:),beta%n_points_beta, rho_I(i), alpha_I)

							multiply_factor1 = dot_product( rel_ij, rand_vec(:,i) )
							
							multiply_factor1 = alpha_I * rho_ij * multiply_factor1 / r_ij_sq / rho_I(i)
							
							!! find rho_I
							rho_ij = 0.0d0
							call beta%spline_int (beta%r,beta%data_rho(itype,:),beta%y2rho(itype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							! find alpha_J
							alpha_J = 0.0d0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
									beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
							
							multiply_factor2 = dot_product( rel_ij, rand_vec(:,j) )
							
							multiply_factor2 = alpha_J * rho_ij * multiply_factor2 / r_ij_sq / rho_I(j)
							
							multiply_factor = multiply_factor1 - multiply_factor2

							this%forces_rnd(1,i) = this%forces_rnd(1,i) + multiply_factor*rel_ij(1)
							this%forces_rnd(2,i) = this%forces_rnd(2,i) + multiply_factor*rel_ij(2)
							this%forces_rnd(3,i) = this%forces_rnd(3,i) + multiply_factor*rel_ij(3)
						end if
					end if
				end do
				v_Te = fdm%Collect_Te(xi,yi,zi)
				this%forces_rnd(1,i) = this%forces_rnd(1,i) * correl_factor_eta * sqrt(v_Te)
				this%forces_rnd(2,i) = this%forces_rnd(2,i) * correl_factor_eta * sqrt(v_Te)
				this%forces_rnd(3,i) = this%forces_rnd(3,i) * correl_factor_eta * sqrt(v_Te)
			end do
!$OMP END DO
!$OMP END PARALLEL
			if ( allocated(rand_vec) ) deallocate(rand_vec)
		end if	!! for condition (this%israndom == 1)
	
		!! -------------------------------------------------
		!! When friction forces need to be calculated
		!! -------------------------------------------------
	
		if (this%isfriction == 1) then

			!! Find auxillary set w_I
			allocate(w_I(3,size(vel,2)))
			w_I = 0.0d0
!$OMP PARALLEL
!$OMP DO
			do ki = 1, Np
				i = eph_for_atoms(ki)
				call getAtomType(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)

				!! alpha_I for the itype atom
				alpha_I = 0.0d0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
					beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)
				do kj = 1, Np
					j = eph_for_atoms(kj)
					if (i /= j) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						r_ij_sq = r_ij*r_ij	
						if (r_ij < beta%r_cutoff) then

							call getAtomType(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type

							!! find rho_J
							rho_ij = 0.0d0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)

							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							call relativeVector(vel(:,i), vel(:,j), rel_v_ij)

							rel_v_ij = rel_v_ij*1000.0d0		! A/fs ----> A/ps
							multiply_factor = alpha_I * rho_ij * dot_product( rel_v_ij, rel_ij )

							multiply_factor = multiply_factor / rho_I(i) / r_ij_sq

							w_I(1,i) = w_I(1,i) + multiply_factor * rel_ij(1)
							w_I(2,i) = w_I(2,i) + multiply_factor * rel_ij(2)
							w_I(3,i) = w_I(3,i) + multiply_factor * rel_ij(3)
						end if
					end if
				end do
			end do
!$OMP END DO
!$OMP END PARALLEL

			!! Find sig_I --> F_fric

!$OMP PARALLEL
!$OMP DO
			do ki = 1, Np
				i = eph_for_atoms(ki)
				call getAtomType(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)

				!! alpha_I for the itype atom
				alpha_I = 0.0d0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
						beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)

				do kj = 1, Np
					j = eph_for_atoms(kj)
					if (i /= j) then

						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						r_ij_sq = r_ij*r_ij

						if (r_ij < beta%r_cutoff) then

							call getAtomType(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type

							!! find rho_J
							rho_ij = 0.0d0
							call beta%spline_int (beta%r, beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)

							call relativeVector(positions(:,j), positions(:,i), rel_ij)

							multiply_factor1 = alpha_I * rho_ij * dot_product( w_I(:,i), rel_ij )

							multiply_factor1 = multiply_factor1 / rho_I(i) / r_ij_sq

							!! alpha_J for the jtype atom
							alpha_J = 0.0d0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
								beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)

							!! find rho_I
							rho_ij = 0.0d0
							call beta%spline_int (beta%r, beta%data_rho(itype,:), &
									beta%y2rho(itype,:), beta%n_points_rho, r_ij, rho_ij)

							multiply_factor2 = alpha_J * rho_ij * dot_product( w_I(:,j), rel_ij )

							multiply_factor2 = multiply_factor2 / rho_I(j) / r_ij_sq
	
							multiply_factor = multiply_factor1 - multiply_factor2

							this%forces_fric(1,i) = this%forces_fric(1,i) - multiply_factor * rel_ij(1) 
							this%forces_fric(2,i) = this%forces_fric(2,i) - multiply_factor * rel_ij(2) 
							this%forces_fric(3,i) = this%forces_fric(3,i) - multiply_factor * rel_ij(3)
						end if
					end if
				end do
			end do
!$OMP END DO
!$OMP END PARALLEL
		end if	!! for condition (this%isfriction == 1)

	end if	!! for condition (this%isfriction == 1 .or. this%israndom == 1)

	!! ----- The current forces get modified -------------------

	!! Different cases of switching ON/OFF the friction and random forces	 
	!! this is the full model
	if (this%isfriction == 1 .and. this%israndom == 1) then
!$OMP PARALLEL
!$OMP DO
		do ki = 1, Np
			i = eph_for_atoms(ki)
			forces(1,i) = forces(1,i) + this%forces_rnd(1,i) + this%forces_fric(1,i)
			forces(2,i) = forces(2,i) + this%forces_rnd(2,i) + this%forces_fric(2,i)
			forces(3,i) = forces(3,i) + this%forces_rnd(3,i) + this%forces_fric(3,i)
		end do
!$OMP END DO
!$OMP END PARALLEL
	!! this is with only friction
	else if (this%isfriction == 1 .and. this%israndom == 0) then
!$OMP PARALLEL
!$OMP DO
		do ki = 1, Np
			i = eph_for_atoms(ki)
			forces(1,i) = forces(1,i) + this%forces_fric(1,i)
			forces(2,i) = forces(2,i) + this%forces_fric(2,i)
			forces(3,i) = forces(3,i) + this%forces_fric(3,i)
		end do
!$OMP END DO
!$OMP END PARALLEL
	!! this is with only random	
	else if (this%isfriction == 0 .and. this%israndom == 1) then			
!$OMP PARALLEL
!$OMP DO
		do ki = 1, Np
			i = eph_for_atoms(ki)
			forces(1,i) = forces(1,i) + this%forces_rnd(1,i)
			forces(2,i) = forces(2,i) + this%forces_rnd(2,i)
			forces(3,i) = forces(3,i) + this%forces_rnd(3,i)
		end do
!$OMP END DO
!$OMP END PARALLEL
	!! this is for not including friction and random, this is definitely not needed
	else if (this%isfriction == 0 .and. this%israndom == 0) then
		do ki = 1, Np
			i = eph_for_atoms(ki)
			do j = 1, 3
				forces(j,i) = forces(j,i)
			end do
		end do
	end if

	end if 	!! when md_time > 0.0

end subroutine eph_LangevinForces


!! Calculate the energies that will be transferred between the atomic and electronic systems
!! due to the forces that have been modified according to the Langevin spatial correlations 
subroutine eph_LangevinEnergyDissipation (this, md_istep, md_time, vel, &
				positions_prev, masses, energies, dt, eph_for_atoms, fdm)

	implicit none

	class (EPH_LangevinSpatialCorrelation_class) :: this

	type (EPH_FDM_class), intent(inout) :: fdm

	integer, intent(in) :: md_istep, eph_for_atoms(:)
	real*8, intent(in) :: dt, md_time, vel(:,:), positions_prev(:,:), masses(:), energies(:)
	integer :: Np, i, j, ki
	real*8 :: xi, yi, zi, Energy_val_fric, Energy_val_rnd, &
	E_val_i_fric, E_val_i_eph, E_val_i_rnd, E_kinetic_atoms, &
	E_pot_atoms, instant_temp_atoms
	real*8, parameter :: boltzconst = 8.61733326E-05

	Np = size(eph_for_atoms)	!! Number of atoms in specified group

	E_kinetic_atoms = 0.0d0
	E_pot_atoms = 0.0d0
	do ki = 1, Np
		i = eph_for_atoms(ki)
		E_kinetic_atoms = E_kinetic_atoms + 0.5d0 * masses(i) * dot_product(vel(1:3, i), vel(1:3, i))
		E_pot_atoms = E_pot_atoms + energies(i)
	end do
	instant_temp_atoms = 2.d0/3.d0/dfloat(Np-1)/boltzconst*E_kinetic_atoms

	if ( md_istep == 0 .or. md_istep == -1 ) then

		!! To write x, y, z, T_e, other quantities mesh map
		if ( fdm%md_last_step == 0 ) then
			open (unit=300, file = this%T_outfile, status = "unknown")
			write(300,*) fdm%nx, ' ', fdm%ny, ' ', fdm%nz
			write(300,*) 'i j k T_e S_e rho_e C_e K_e flag T_dyn_flag'		! Column headers
			close(unit = 300)
		end if

		!! To write the energies and T's at some time steps to a file
		if ( fdm%md_last_step == 0 ) then
			open (unit = 100, file = "eph-EnergySharingData.txt", status = "unknown")
			write(100,*) 'Time (fs) / E_fric(eV) / E_rand(eV) / E_net_cum(eV) / &
							T_e (K) / T_a (K) / Kin.E_a (eV) / Pot.E_a (eV)' 
			write(100,1) md_time, 0.0, 0.0, 0.0, sum(fdm%T_e)/fdm%ntotal, &
							instant_temp_atoms, E_kinetic_atoms, E_pot_atoms 
		end if

	else
		if (MOD(md_istep, this%T_out_freq) == 0) open (unit = 100, file = "eph-EnergySharingData.txt", &
											status = "old", position = "append")
	end if

	!! Do all calculations for MD time greater than 0

	if ( .not. (md_istep == 0 .or. md_istep == -1) ) then

	!! ----- Calculate the energies for exchange ----------

	!! Different cases of switching ON/OFF the friction and random forces

	Energy_val_fric = 0.0d0; Energy_val_rnd = 0.0d0

	!! this is the full model
	if (this%isfriction == 1 .and. this%israndom == 1) then
		do ki = 1, Np
			i = eph_for_atoms(ki)
			E_val_i_fric = 0.0d0; E_val_i_rnd = 0.0d0; E_val_i_eph = 0.0d0
			xi = positions_prev(1,i); yi = positions_prev(2,i); zi = positions_prev(3,i)
			E_val_i_fric = E_val_i_fric - dt*dot_product(this%forces_fric(1:3,i), vel(1:3,i))
			E_val_i_rnd = E_val_i_rnd - dt*dot_product(this%forces_rnd(1:3,i), vel(1:3,i))

			E_val_i_eph = E_val_i_fric + E_val_i_rnd

			if (this%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_eph, dt)

			Energy_val_fric = Energy_val_fric + E_val_i_fric
			Energy_val_rnd = Energy_val_rnd + E_val_i_rnd		
		end do

	!! this is with only friction
	else if (this%isfriction == 1 .and. this%israndom == 0) then
		do ki = 1, Np
			i = eph_for_atoms(ki)
			E_val_i_fric = 0.0d0
			xi = positions_prev(1,i); yi = positions_prev(2,i); zi = positions_prev(3,i)
			E_val_i_fric = E_val_i_fric - dt*dot_product(this%forces_fric(1:3,i), vel(1:3,i))

			if (this%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_fric, dt)
			Energy_val_fric = Energy_val_fric + E_val_i_fric
		end do

	!! this is with only random	
	else if (this%isfriction == 0 .and. this%israndom == 1) then	
		do ki = 1, Np
			i = eph_for_atoms(ki)
			E_val_i_rnd = 0.0d0
			xi = positions_prev(1,i); yi = positions_prev(2,i); zi = positions_prev(3,i)
			E_val_i_rnd = E_val_i_rnd - dt*dot_product(this%forces_rnd(1:3,i), vel(1:3,i))

			if (this%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_rnd, dt)
			Energy_val_rnd = Energy_val_rnd + E_val_i_rnd
		end do

	!! this is for not including friction and random, this is definitely not needed
	else if (this%isfriction == 0 .and. this%israndom == 0) then
		do ki = 1, Np
			i = eph_for_atoms(ki)
			xi = positions_prev(1,i); yi = positions_prev(2,i); zi = positions_prev(3,i)
			if (this%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, 0.0d0, dt)
		end do
	end if

	!! To calculate electronic temperature

	if (this%fdm_option == 1) call fdm%heatDiffusionSolve (dt)

	!! To write the electronic mesh temperatures

	if (this%T_out_mesh_freq /= 0) then
		if (MOD(md_istep, this%T_out_mesh_freq) == 0) then
			call fdm%saveOutputToFile ( this%T_outfile, md_istep, dt )
		end if
	end if

	!! To write the energy transfers within atomic-electonic system and temperatures

	this%E_eph_cumulative = this%E_eph_cumulative + Energy_val_fric + Energy_val_rnd

	if (MOD(md_istep, this%T_out_freq) == 0) then
		write(100,1) md_time + this%md_prev_time, Energy_val_fric, Energy_val_rnd, &
					this%E_eph_cumulative, sum(fdm%T_e)/fdm%ntotal, instant_temp_atoms, &
					E_kinetic_atoms, E_pot_atoms
		close(unit=100)
	end if

1	format(1E13.6,3E14.6,3E13.6,1E14.6)

	call fdm%beforeNextFeedback()

	end if 	!! when md_time > 0.0

end subroutine eph_LangevinEnergyDissipation


!! Find distance between two atoms.
real function get_distance(xi,yi,zi,xj,yj,zj) result (r_ij)
	implicit none
	real*8, intent(in) :: xi,yi,zi,xj,yj,zj
	r_ij = sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
end function get_distance


!! Find the relative vectors
subroutine relativeVector(vec1, vec2, vec12)
	implicit none
	real*8, intent(in) :: vec1(3), vec2(3)
	real*8, intent(out) :: vec12(3)
	vec12(1) = vec1(1) - vec2(1)
	vec12(2) = vec1(2) - vec2(2)
	vec12(3) = vec1(3) - vec2(3)
end subroutine relativeVector


!! Determine which type of atom is atom i
!! It is assumed that data in the beta file for different atoms is
!! according to the corresponding types of them as specified in
!! the input file.
subroutine getAtomType(p,natomtypes,masses,type_mass,atom_type)
	implicit none
	integer, intent(in) :: p, natomtypes
	real*8, intent(in) :: masses(:), type_mass(:)
	integer, intent(out) :: atom_type
	integer :: k
	do k = 1, natomtypes
		if (masses(p) == type_mass(k)) then
			atom_type = k
			exit
		end if
	end do
end subroutine getAtomType


!! Box-Muller transmorfation to get (nearly) an array of Gaussian random variates
!! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)
subroutine randomGaussianArray(Np, mean, standard_deviation, rand_array)
	implicit none
	integer, intent(in) :: Np
	real*8, intent(in) :: mean, standard_deviation
	real*8, intent(out) :: rand_array(Np)
	integer :: i
	real*8 :: one_value, mean_actual, sd_actual
	real*8 :: Pi = 4.0d0*atan(1.0d0)

	!! Uniformly distributed random numbers from 0 to 1
	call random_number(rand_array) 

	!! Box-Muller transmorfation to normally distributed random numbers
	do i = 1, Np-1, 2
		one_value = standard_deviation * &
		sqrt(-2.0d0*log(rand_array(i))) * sin(2.0d0*Pi*rand_array(i+1)) + mean

		rand_array(i+1) = standard_deviation * &
		sqrt(-2.0d0*log(rand_array(i))) * cos(2.0d0*Pi*rand_array(i+1)) + mean

		rand_array(i) = one_value
	end do

	!! Mean and standard deviation must be brought to required values 
	!! This works very good when mean is 0
	mean_actual = sum(rand_array)/Np
	sd_actual = sqrt(sum((rand_array - mean_actual)**2)/Np)
	do i = 1, Np
		rand_array(i) = rand_array(i) - mean_actual + mean
		rand_array(i) = rand_array(i) * standard_deviation/sd_actual
	end do
end subroutine randomGaussianArray

end module eph_electronic_stopping
