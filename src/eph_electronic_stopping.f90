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
! PRB 94, 024305 (2016); PRB 99, 174301 (2019).
!********* by Uttiyoarnab Saha

!**************************************************************************

module eph_electronic_stopping

use eph_fdm
use eph_beta

contains

subroutine eph_Langevin_spatial_correlation(isfriction, israndom, vel, forces, masses, &
			type_mass, md_istep, dt, md_time, positions, natomtypes, &
			T_outfile, T_out_freq, T_out_mesh_freq, model_eph, beta, fdm)
	implicit none
	type (EPH_Beta_class), intent(in) :: beta
	type (EPH_FDM_class), intent(in) :: fdm
	
	character*128, intent(in) :: T_outfile
	integer, intent(in) :: natomtypes, md_istep, T_out_freq, isfriction, israndom, model_eph, &
	T_out_mesh_freq
	real*8, intent(in) :: vel(:,:), positions(:,:), masses(:), type_mass(:), dt, md_time
	real*8, intent(inout) :: forces(:,:)
	real*8, allocatable :: forces_fric(:,:), forces_rnd(:,:)
	integer :: Np, i, j, k, itype, jtype, ktype, atom_type, component
	
	real*8 :: beta_rho, xi, yi, zi, xj, yj, zj, r_ij, &
	alpha_I, alpha_J, &
	sig_I(3),B_IJ(3,3),eta_I(3), W_IJ_I(3,3), W_IJ_J(3,3), W_IJ(3,3), &
	outer_p(3,3), rho_ij, v_Te, rel_ij(3), &
	rel_v_ij(3), correl_factor_eta, &
	Energy_val_fric,Energy_val_rnd, Energy_val_eph, E_val_i_fric,E_val_i_rnd,E_val_i_eph, &
	multiply_factor, r_v1, r_v2, r_ik_sq, st_dev_factor, & 
	multiply_factor1, multiply_factor2, &
	sum_of_forces
	
	
	real*8, allocatable :: rand_vec(:,:), rho_I(:), w_I(:,:)
	real*8, parameter :: boltzconst = 8.61733326E-5 
	
	Np = size(vel,2)
	
	! To write the electronic energy loss data to a file after evaluation at each time step

	if ( md_istep == 0 .or. md_istep == -1 ) then
		open (unit = 100, file = "eph-EnergySharingData.txt", status = "unknown")
		write(100,*) 'Time (fs)  E_fric(eV) E_rand(eV) E_net(eV) T_e (K)' 
		else
			open (unit = 100, file = "eph-EnergySharingData.txt", status = "old", position = "append")
	end if
	
	! generate uncorrelated random vectors
	allocate(rand_vec(3,Np))
	do i = 1, Np
		rand_vec(1,i) = random_gaussian()
		rand_vec(2,i) = random_gaussian()
		rand_vec(3,i) = random_gaussian()
	end do
	
	! The random vectors are component by component follow an independent
	! Gaussian distribution, related to a heat bath at temperature T_e
	! PRL 120, 185501 (2018).  
	
	! the random forces need finally to have certain correlations
	  
	!v_Te = sum(fdm%T_e) / fdm%ntotal
	!correl_factor_eta = sqrt(2.0*boltzconst*v_Te/(dt))
	!call get_distribution_factor(rand_vec(1,:), Np, correl_factor_eta, st_dev_factor)
	!rand_vec(1,:) = rand_vec(1,:)*st_dev_factor
	!call get_distribution_factor(rand_vec(2,:), Np, correl_factor_eta, st_dev_factor)
	!rand_vec(2,:) = rand_vec(2,:)*st_dev_factor
	!call get_distribution_factor(rand_vec(3,:), Np, correl_factor_eta, st_dev_factor)
	!rand_vec(3,:) = rand_vec(3,:)*st_dev_factor

	
	allocate(forces_fric(3,Np), forces_rnd(3,Np))

	forces_fric = 0.0
	forces_rnd = 0.0

	! .......... Model 1 ..........
	! Implementation of random forces and friction forces 
	! according to PRL 120, 185501 (2018)
	
	if (model_eph == 1) then
	
		! Find total atomic electronic densities at all sites limited by the 
		! value of r_cutoff given in the beta file.

		allocate(rho_I(Np))
		rho_I = 0.0
		do i = 1, Np
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			do j = 1, Np
				if (i /= j) then
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi,xj,yj,zj)
					if (r_ij <= beta%r_cutoff) then
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						rho_ij = 0.0
						call beta%spline_int (beta%r,beta%data_rho(jtype,:), &
						beta%y2rho(jtype,:),beta%n_points_rho, r_ij, rho_ij)
						
						rho_I(i) = rho_I(i) + rho_ij
					end if
				end if
			end do
		end do
		
		! To find sigma_I and eta_I ......

		if (isfriction == 1 .or. israndom == 1) then
		
		! -------------------------------------------------
		! When random forces need to be calculated
		! -------------------------------------------------
		
		if (israndom == 1) then

			do i = 1, Np
 
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				do j = 1, Np
					W_IJ_J = 0.0
					W_IJ_I = 0.0
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						
						if (r_ij <= beta%r_cutoff) then
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							
							!outer_p = 0.0
							!call get_outer_product(rel_ij, outer_p)
							
							! find rho_I
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(itype,:),beta%y2rho(itype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							! find alpha_J
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
									beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
	
							!W_IJ_J = - alpha_J * outer_p * rho_ij / r_ij / rho_I(j)
							
							!eta_I = 0.0
							!call get_eta_I(W_IJ_J,rand_vec(:,j),eta_I)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,j) )
							multiply_factor = - alpha_J * rho_ij * multiply_factor / r_ij / rho_I(j)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)	!eta_I(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)	!eta_I(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)	!eta_I(3)
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:),beta%y2rho(jtype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							! find alpha_I	
							alpha_I = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
									beta%y2alpha(itype,:),beta%n_points_beta, rho_I(i), alpha_I)

							!W_IJ_I = alpha_I * outer_p * rho_ij / r_ij / rho_I(i)

							!eta_I = 0.0
							!call get_eta_I(W_IJ_I,rand_vec(:,i),eta_I)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,i) )
							multiply_factor = alpha_I * rho_ij * multiply_factor / r_ij / rho_I(i)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)	!eta_I(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)	!eta_I(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)	!eta_I(3)
						end if
					end if
				end do
				
				v_Te = fdm%Collect_Te(xi,yi,zi)
				correl_factor_eta = sqrt(2.0*boltzconst*v_Te/(dt))
				forces_rnd(1,i) = forces_rnd(1,i) * correl_factor_eta
				forces_rnd(2,i) = forces_rnd(2,i) * correl_factor_eta
				forces_rnd(3,i) = forces_rnd(3,i) * correl_factor_eta
			end do
		end if
	
		! -------------------------------------------------
		! When friction forces need to be calculated
		! -------------------------------------------------
	
		if (isfriction == 1) then
	
			! Find auxillary set w_I
			allocate(w_I(3,Np))
	
			do i = 1, Np
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				! alpha_I for the itype atom
				alpha_I = 0.0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
					beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)
				
				do j = 1, Np
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij <= beta%r_cutoff) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							call relative_vector(vel(:,i), vel(:,j), rel_v_ij)
							
							multiply_factor = alpha_I * rho_ij * dot_product( rel_v_ij, rel_ij )
							
							multiply_factor = multiply_factor  / rho_I(i) / r_ij
							
							w_I(1,i) = w_I(1,i) + multiply_factor * rel_ij(1)
							w_I(2,i) = w_I(2,i) + multiply_factor * rel_ij(2)
							w_I(3,i) = w_I(3,i) + multiply_factor * rel_ij(3)
						end if
					end if
				end do
			end do
	
			! Find sig_I --> F_fric
	
			do i = 1, Np
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				! alpha_I for the itype atom
				alpha_I = 0.0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
						beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)
				
				do j = 1, Np
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij <= beta%r_cutoff) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r, beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							
							multiply_factor1 = alpha_I * rho_ij * dot_product( w_I(:,i), rel_ij )
							
							multiply_factor1 = multiply_factor1 / rho_I(i) / r_ij
							
							
							! alpha_J for the jtype atom
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
								beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
							
							! find rho_I
							rho_ij = 0.0
							call beta%spline_int (beta%r, beta%data_rho(itype,:), &
									beta%y2rho(itype,:), beta%n_points_rho, r_ij, rho_ij)
							
							multiply_factor2 = alpha_J * rho_ij * dot_product( w_I(:,j), rel_ij )
						
							multiply_factor2 = multiply_factor2 / rho_I(j) / r_ij
							
							sig_I(1) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(1)
							sig_I(2) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(2)
							sig_I(3) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(3)
							
							forces_fric(1,i) = forces_fric(1,i) - sig_I(1) 
							forces_fric(2,i) = forces_fric(2,i) - sig_I(2) 
							forces_fric(3,i) = forces_fric(3,i) - sig_I(3)
						
						end if
					end if
				end do
			end do
					
		end if	! for condition (isfriction == 1)
		
		end if	! for condition (isfriction == 1 .or. israndom == 1)
	
	end if	! for condition (model == 1)
	
	
	!=========================================================================
	! ***** The following model is as per USER-FIX EPH in LAMMPS ********
	!=========================================================================
	
	! .......... Model 2 ..........
	! As per USER-EPH FIX in LAMMPS which does with "rho_r_sq"
	
	if (model_eph == 2) then
		
		! Find total atomic electronic densities at all sites limited by the 
		! value of r_cutoff given in the beta file.

		allocate(rho_I(Np))
		rho_I = 0.0
		do i = 1, Np
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			do j = 1, Np
				if (i /= j) then
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi,xj,yj,zj)
					if (r_ij**2 <= beta%r_cutoff_sq) then
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						rho_ij = 0.0
						call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
						beta%y2rho_r_sq(jtype,:),beta%n_points_rho, r_ij**2, rho_ij)
						
						rho_I(i) = rho_I(i) + rho_ij
					end if
				end if
			end do
		end do
		
		! To find sigma_I and eta_I ......

		if (isfriction == 1 .or. israndom == 1) then
		
		! -------------------------------------------------
		! When random forces need to be calculated
		! -------------------------------------------------
		
		if (israndom == 1) then

			do i = 1, Np
 
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				do j = 1, Np
					W_IJ_J = 0.0
					W_IJ_I = 0.0
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							
							!outer_p = 0.0
							!call get_outer_product(rel_ij, outer_p)
							
							! find rho_I_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(itype,:), &
								beta%y2rho_r_sq(itype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							! find alpha_J
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
									beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
	
							!W_IJ_J = - alpha_J * outer_p * rho_ij / r_ij**2 / rho_I(j)
							
							!eta_I = 0.0
							!call get_eta_I(W_IJ_J,rand_vec(:,j),eta_I)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,j) )
							multiply_factor = - alpha_J * rho_ij * multiply_factor / r_ij**2 / rho_I(j)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)	!eta_I(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)	!eta_I(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)	!eta_I(3)
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)

							! find alpha_I	
							alpha_I = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
									beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)

							!W_IJ_I = alpha_I * outer_p * rho_ij / r_ij**2 / rho_I(i)

							!eta_I = 0.0
							!call get_eta_I(W_IJ_I,rand_vec(:,i),eta_I)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,i) )
							multiply_factor = alpha_I * rho_ij * multiply_factor / r_ij**2 / rho_I(i)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)	!eta_I(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)	!eta_I(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)	!eta_I(3)
						end if
					end if
				end do
				
				v_Te = fdm%Collect_Te(xi,yi,zi)
				correl_factor_eta = sqrt(2.0*boltzconst*v_Te/(dt))
				forces_rnd(1,i) = forces_rnd(1,i) * correl_factor_eta
				forces_rnd(2,i) = forces_rnd(2,i) * correl_factor_eta
				forces_rnd(3,i) = forces_rnd(3,i) * correl_factor_eta
			end do
		end if
	
		! -------------------------------------------------
		! When friction forces need to be calculated
		! -------------------------------------------------
	
		if (isfriction == 1) then
	
			! Find auxillary set w_I
			allocate(w_I(3,Np))
	
			do i = 1, Np
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				! alpha_I for the itype atom
				alpha_I = 0.0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
					beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)
				
				do j = 1, Np
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							call relative_vector(vel(:,i), vel(:,j), rel_v_ij)
							
							multiply_factor = alpha_I * rho_ij * dot_product( rel_v_ij, rel_ij )
							
							multiply_factor = multiply_factor  / rho_I(i) / r_ij**2
							
							w_I(1,i) = w_I(1,i) + multiply_factor * rel_ij(1)
							w_I(2,i) = w_I(2,i) + multiply_factor * rel_ij(2)
							w_I(3,i) = w_I(3,i) + multiply_factor * rel_ij(3)
						end if
					end if
				end do
			end do
	
			! Find sig_I --> F_fric
	
			do i = 1, Np
				call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
				itype = atom_type
				xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
				
				! alpha_I for the itype atom
				alpha_I = 0.0
				call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
						beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)
				
				do j = 1, Np
					if (j /= i) then
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							call relative_vector(positions(:,i), positions(:,j), rel_ij)
							
							multiply_factor1 = alpha_I * rho_ij * dot_product( w_I(:,i), rel_ij )
							
							multiply_factor1 = multiply_factor1 / rho_I(i) / r_ij**2
							
							
							! alpha_J for the jtype atom
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
								beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
							
							! find rho_I_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(itype,:), &
									beta%y2rho_r_sq(itype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							multiply_factor2 = alpha_J * rho_ij * dot_product( w_I(:,j), rel_ij )
						
							multiply_factor2 = multiply_factor2 / rho_I(j) / r_ij**2
							
							sig_I(1) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(1)
							sig_I(2) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(2)
							sig_I(3) = ( multiply_factor1 - multiply_factor2 ) * rel_ij(3)
							
							forces_fric(1,i) = forces_fric(1,i) - sig_I(1) 
							forces_fric(2,i) = forces_fric(2,i) - sig_I(2) 
							forces_fric(3,i) = forces_fric(3,i) - sig_I(3)
						
						end if
					end if
				end do
			end do
					
		end if	! for condition (isfriction == 1)
		
		end if	! for condition (isfriction == 1 .or. israndom == 1)
	
	end if	! for condition (model == 2)

	
	! The current forces get modified
	
	! First ensure that sum_I(eta_I) = 0, sum_I(sigma_I) = 0

	!if (isfriction == 1) then
	!	sum_of_forces = sum(forces_fric(1,:))/Np
	!	forces_fric(1,:) = forces_fric(1,:) - sum_of_forces
	!	
	!	sum_of_forces = sum(forces_fric(2,:))/Np
	!	forces_fric(2,:) = forces_fric(2,:) - sum_of_forces
	!	
	!	sum_of_forces = sum(forces_fric(3,:))/Np
	!	forces_fric(3,:) = forces_fric(3,:) - sum_of_forces
	!end if
	!
	!if (israndom == 1) then
	!	sum_of_forces = sum(forces_rnd(1,:))/Np
	!	forces_rnd(1,:) = forces_rnd(1,:) - sum_of_forces
	!	
	!	sum_of_forces = sum(forces_rnd(2,:))/Np
	!	forces_rnd(2,:) = forces_rnd(2,:) - sum_of_forces
	!	
	!	sum_of_forces = sum(forces_rnd(3,:))/Np
	!	forces_rnd(3,:) = forces_rnd(3,:) - sum_of_forces
	!end if
	
	
	! Different cases of switching ON/OFF the friction and random forces
		 
	! this is the full model
	if (isfriction == 1 .and. israndom == 1) then
	
		Energy_val_fric = 0.0; Energy_val_rnd = 0.0; Energy_val_eph = 0.0
	
		do i = 1, Np
			E_val_i_fric = 0.0; E_val_i_rnd = 0.0; E_val_i_eph = 0.0
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i) 
	
			do component = 1, 3
		!! ---- forces are adjusted by friction according to the relative signs of forces and velocities ----
		!! when velocity and force are opposite, friction acts along the force
		!		if ((vel(component,i)>0 .and. forces(component,i)<0).or. &
		!		(vel(component,i)<0 .and. forces(component,i)>0)) then
		!			if (forces(component,i)<0) then
		!				forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
		!			end if
		!			if (forces(component,i)>0) then
		!				forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
		!			end if
		!		end if
		!!! when velocity and force are in same direction, friction acts opposite to the force
		!		if ((vel(component,i)>0 .and. forces(component,i)>0).or. &
		!		(vel(component,i)<0 .and. forces(component,i)<0)) then
		!			if (forces(component,i)<0) then
		!				forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
		!			end if
		!			if (forces(component,i)>0) then 
		!				forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
		!			end if
		!		end if
				
				! ------ up to here -------
				
				forces(component,i) = forces(component,i) + forces_rnd(component,i) + forces_fric(component,i)
	
				E_val_i_fric = E_val_i_fric - dt*forces_fric(component,i)*vel(component,i)
				E_val_i_rnd = E_val_i_rnd - dt*forces_rnd(component,i)*vel(component,i)
			end do
			
			E_val_i_eph = E_val_i_eph + E_val_i_fric + E_val_i_rnd
			
			if (fdm%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_eph, dt)
		
			Energy_val_fric = Energy_val_fric + E_val_i_fric
			Energy_val_rnd = Energy_val_rnd + E_val_i_rnd		
		end do
		Energy_val_eph = Energy_val_fric + Energy_val_rnd
	end if
	
	! this is with only friction
	if (isfriction == 1 .and. israndom == 0) then			
		Energy_val_fric = 0.0; Energy_val_rnd = 0.0; Energy_val_eph = 0.0
		do i = 1, Np
			E_val_i_fric = 0.0
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			
			do component = 1, 3
		!! ---- forces are adjusted by friction according to the relative signs of forces and velocities ----	
		!! when velocity and force are opposite, friction acts along the force
		!		if ((vel(component,i)>0 .and. forces(component,i)<0).or. &
		!		(vel(component,i)<0 .and. forces(component,i)>0)) then
		!			if (forces(component,i)<0) then
		!				forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
		!			end if
		!			if (forces(component,i)>0) then
		!				forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
		!			end if
		!		end if
		!! when velocity and force are in same direction, friction acts opposite to the force
		!		if ((vel(component,i)>0 .and. forces(component,i)>0).or. &
		!		(vel(component,i)<0 .and. forces(component,i)<0)) then
		!			if (forces(component,i)<0) then
		!				forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
		!			end if
		!			if (forces(component,i)>0) then 
		!				forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
		!			end if
		!		end if
				forces(component,i) = forces(component,i) + forces_fric(component,i)
				E_val_i_fric = E_val_i_fric - dt*forces_fric(component,i)*vel(component,i)
			end do
			
			if (fdm%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_fric, dt)
			Energy_val_fric = Energy_val_fric + E_val_i_fric
		end do
		Energy_val_eph = Energy_val_fric
	end if
	
	! this is with only random	
	if (isfriction == 0 .and. israndom == 1) then			
		
		Energy_val_fric = 0.0; Energy_val_rnd = 0.0; Energy_val_eph = 0.0
		
		do i = 1, Np
			E_val_i_rnd = 0.0
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
		
			do component = 1, 3
				forces(component,i) = forces(component,i) + forces_rnd(component,i)
				E_val_i_rnd = E_val_i_rnd - dt*forces_rnd(component,i)*vel(component,i)
			end do
			
			if (fdm%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_rnd, dt)
			Energy_val_rnd = Energy_val_rnd + E_val_i_rnd
		end do
		Energy_val_eph = Energy_val_rnd
	end if
	
	! this is for not including friction and random, this is definitely not needed
	if (isfriction == 0 .and. israndom == 0) then
		Energy_val_fric = 0.0; Energy_val_rnd = 0.0; Energy_val_eph = 0.0
		do i = 1, Np
			do component = 1, 3
				forces(component,i) = forces(component,i)
			end do
		end do
	end if

	! To calculate electronic temperature

	if (fdm%fdm_option == 1) call fdm%heat_diffusion_solve (dt)

	! To write the electronic mesh temperatures

	if (T_out_mesh_freq /= 0) then
		if (MOD(md_istep, T_out_mesh_freq) == 0) then
			call fdm%save_output_to_file (T_outfile,md_istep,dt)
		end if
	end if

	! To write the energy transfers within atomic-electonic system and electronic temperature

	if (MOD(md_istep, T_out_freq) == 0) then
		write(100,1) md_time, Energy_val_fric, Energy_val_rnd, Energy_val_eph, sum(fdm%T_e)/fdm%ntotal 
	end if
1	format(1E13.6,3E14.6,1E13.6)	
	close(unit=100)
	
	call fdm%beforeNextFeedback()

end subroutine eph_Langevin_spatial_correlation


! Find distance between two atoms.
real function get_distance(xi,yi,zi,xj,yj,zj) result (r_ij)
	implicit none
	real*8, intent(in) :: xi,yi,zi,xj,yj,zj
	r_ij = sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
end function get_distance

! Find the relative vectors
subroutine relative_vector(vec1, vec2, vec12)
	implicit none
	real*8, intent(in) :: vec1(3), vec2(3)
	real*8, intent(out) :: vec12(3)
	vec12(1) = vec1(1) - vec2(1)
	vec12(2) = vec1(2) - vec2(2)
	vec12(3) = vec1(3) - vec2(3)
end subroutine relative_vector


! Determine which type of atom is atom i
! It is assumed that data in the beta file for different atoms is
! according to the corresponding types of them as specified in
! the input file.
subroutine get_atom_type(p,natomtypes,masses,type_mass,atom_type)
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
end subroutine get_atom_type


! Find the outer product of the unit vectors which are
! used to calculate the random forces, eta_I.
!subroutine get_outer_product(vector,outer_p)
!	implicit none
!	real*8, intent(in) :: vector(3)
!	real*8, intent(out) :: outer_p(3,3)
!	real*8 :: vec(3,1), vecT(1,3)
!	integer :: l, m
!	
!	vec(1,1) = vector(1); vec(2,1) = vector(2); vec(3,1) = vector(3)
!	vecT(1,1) = vector(1); vecT(1,2) = vector(2); vecT(1,3) = vector(3) 
!	do l = 1, 3
!		do m = 1, 3
!			outer_p(l,m) = vec(l,1) * vecT(1,m)
!		end do
!	end do
!end subroutine get_outer_product


! Find eta_I, which are the random forces on the atoms. It is
! found using the projection of the random vectors associated with
! every atom on the (outer product of) relative unit vectors
! of the atoms.
!subroutine get_eta_I(W_ij,vector,eta_I)
!	implicit none
!	real*8, intent(in) :: W_ij(3,3), vector(3)
!	real*8, intent(out) :: eta_I(3)
!	integer :: i, j
!	
!	eta_I = 0.0
!	do i = 1, 3
!		do j = 1, 3
!			eta_I(i) = eta_I(i) + W_ij(i,j) * vector(j)
!		end do
!	end do
!end subroutine get_eta_I

! Find sig_I, which are the friction forces on the atoms. It is
! found using the projection of B_ij associated with W_ij by fluctuation-
! dissipation theorem on the relative velocities of the atoms.
!subroutine get_sigma_I(B_ij,vector,sig_I)
!	implicit none
!	real*8, intent(in) :: B_ij(3,3), vector(3)
!	real*8, intent(out) :: sig_I(3)
!	integer :: i, j
!	
!	sig_I = 0.0
!	do i = 1, 3
!		do j = 1, 3
!			sig_I(i) = sig_I(i) + B_ij(i,j) * vector(j)
!		end do
!	end do
!end subroutine get_sigma_I


! Box-Muller transmorfation to get (nearly) Gaussian random variate
real function random_gaussian() result(randg)
	implicit none
	real*8 :: n1, n2, R, theta
	real*8, parameter :: twoPi = 2.0*3.14285714
	call random_number(n1)
	call random_number(n2)
	R = sqrt(-2.0*log(n1)/n1)
	theta = twoPi*n2
	randg = R * cos(theta)
end function random_gaussian


! Box-Muller transmorfation to get Gaussian random variate
! Choosing n1 and n2 as ordinate and abscissa of a random point
! inside the unit circle around the origin.
! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)

! When the distribution of each of the components x, y, z are checked 
! individually, they do not show a Gaussian distribution if we use this
! method. Neither there is correct variation of the temperatures of the
! electronic and atomic systems observed with this, specifically, the
! temperatures do not show any significant variation either from 
! the beginning or after some time. So, it is commented out.

!real function random_gaussian() result(randg)
!	implicit none
!	
!	integer :: iset
!	real*8 :: n1, n2, Rsq, factor, randn1
!	save :: iset, randn1
!	data iset/0/
!	!real*8, parameter :: twoPi = 2.0*3.14285714
!	Rsq = 1.0
!	if (iset == 0) then
!		do while ( Rsq >= 1.0 .or. Rsq == 0.0)
!			call random_number(n1)
!			call random_number(n2)
!			n1 = 2.0*n1 - 1.0
!			n2 = 2.0*n2 - 1.0
!			Rsq = n1**2 + n2**2	
!		end do
!		factor = sqrt(-2.0*log(Rsq)/Rsq)
!		randn1 = n1*factor
!		randg = n2*factor
!		iset = 1
!	else
!		randg = randn1
!		iset = 0
!	end if
!end function random_gaussian


! The random_gaussian() function generate distributions of each component
! of the vector with different variance from the required values.
! So, the distributions are re-scaled to have the required variance.  

! There is a correlation factor between 'eta', the random forces
! which involves the T_e at position R_I(x,y,z) of the atom and not on the
! overall temperature of the so-called "electronic heat bath". So, the following
! subroutine seems not not to be a necessity. 

!subroutine get_distribution_factor(vect, Np, st_dev_req, st_dev_factor)
!	implicit none
!	integer, intent(in) :: Np
!	real*8, intent(in) :: vect(Np), st_dev_req
!	real*8, intent(out) :: st_dev_factor
!	real*8 :: st_dev_now, mean
!	
!	mean = sum(vect)/ size(vect)
!	st_dev_now = sqrt( sum( (vect(:) - mean)**2 ) / size(vect) )
!	st_dev_factor = st_dev_req / st_dev_now
!end subroutine get_distribution_factor

end module eph_electronic_stopping

