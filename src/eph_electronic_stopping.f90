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

subroutine eph_Langevin_spatial_correlation (isfriction, israndom, vel, forces, masses, &
			type_mass, md_istep, dt, md_time, positions, neighbors_list, n_neigh, &
			natomtypes, T_outfile, T_out_freq, T_out_mesh_freq, model_eph, beta, fdm)
	implicit none
	type (EPH_Beta_class), intent(in) :: beta
	type (EPH_FDM_class), intent(in) :: fdm
	
	character*128, intent(in) :: T_outfile
	integer, intent(in) :: neighbors_list(:), n_neigh(:), natomtypes, md_istep, &
	T_out_freq, isfriction, israndom, model_eph, T_out_mesh_freq
	real*8, intent(in) :: vel(:,:), positions(:,:), masses(:), type_mass(:), dt, md_time
	real*8, intent(inout) :: forces(:,:)
	real*8, allocatable :: forces_fric(:,:), forces_rnd(:,:)
	integer :: Np, i, j, jneigh, jneigh_start, jneigh_stop, itype, jtype, atom_type, component
	
	real*8 :: xi, yi, zi, xj, yj, zj, r_ij, alpha_I, alpha_J, rho_ij, v_Te, rel_ij(3), &
	rel_v_ij(3), correl_factor_eta, Energy_val_fric,Energy_val_rnd, Energy_val_eph, &
	E_val_i_fric, E_val_i_eph, E_val_i_rnd, multiply_factor, multiply_factor1, &
	multiply_factor2, st_dev_factor
	
	
	real*8, allocatable :: rand_vec(:,:), rho_I(:), w_I(:,:)
	real*8, parameter :: boltzconst = 8.61733326E-05 
	
	Np = size(vel,2)
	
	! To write the electronic energy loss data to a file after evaluation at each time step

	if ( md_istep == 0 .or. md_istep == -1 ) then
		! x, y, z, T_e, other quantities mesh map 
		open (unit=300, file = T_outfile, status = "unknown")
		write(300,*) fdm%nx, ' ', fdm%ny, ' ', fdm%nz
		write(300,*) 'i j k T_e S_e rho_e C_e K_e flag T_dyn_flag'		! Column headers
		close(unit = 300)
		! energy transfers and T_e at some time steps
		open (unit = 100, file = "eph-EnergySharingData.txt", status = "unknown")
		write(100,*) 'Time (fs)  E_fric(eV) E_rand(eV) E_net(eV) T_e (K)' 
		write(100,1) md_time, 0.0, 0.0, 0.0, sum(fdm%T_e)/fdm%ntotal 
		else
			open (unit = 100, file = "eph-EnergySharingData.txt", status = "old", position = "append")
	end if
	
	if ( .not. (md_istep == 0 .or. md_istep == -1) ) then
	
	! generate uncorrelated random vectors
	allocate(rand_vec(3,Np))
	call randomGaussianArray(Np, 0.0d0, 1.0d0, rand_vec(1,:))
	call randomGaussianArray(Np, 0.0d0, 1.0d0, rand_vec(2,:))
	call randomGaussianArray(Np, 0.0d0, 1.0d0, rand_vec(3,:))
	
	correl_factor_eta = sqrt(2.0*boltzconst*1000.0/dt)	! time units fs ----> ps

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
			
			jneigh_start = (i-1)*n_neigh(i) + 1
			jneigh_stop = i*n_neigh(i)
			
			do jneigh = jneigh_start, jneigh_stop
				if (i /= neighbors_list(jneigh)) then
					j = neighbors_list(jneigh)
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi, xj,yj,zj)
					if (r_ij < beta%r_cutoff) then
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						
						if (r_ij < beta%r_cutoff) then
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							
							! find rho_I
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(itype,:),beta%y2rho(itype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							! find alpha_J
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
									beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,j) )
							multiply_factor = - alpha_J * rho_ij * multiply_factor / r_ij / rho_I(j)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:),beta%y2rho(jtype,:), &
											beta%n_points_rho, r_ij, rho_ij)
							
							! find alpha_I	
							alpha_I = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
									beta%y2alpha(itype,:),beta%n_points_beta, rho_I(i), alpha_I)

							multiply_factor = dot_product( rel_ij, rand_vec(:,i) )
							multiply_factor = alpha_I * rho_ij * multiply_factor / r_ij / rho_I(i)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)
						end if
					end if
				end do
				
				v_Te = fdm%Collect_Te(xi,yi,zi)
				forces_rnd(1,i) = forces_rnd(1,i) * correl_factor_eta * sqrt(v_Te)
				forces_rnd(2,i) = forces_rnd(2,i) * correl_factor_eta * sqrt(v_Te)
				forces_rnd(3,i) = forces_rnd(3,i) * correl_factor_eta * sqrt(v_Te)
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij < beta%r_cutoff) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r,beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)
							
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							call relativeVector(vel(:,i), vel(:,j), rel_v_ij)

							rel_v_ij = rel_v_ij*1000.0		! A/fs ----> A/ps
							multiply_factor = alpha_I * rho_ij * dot_product( rel_v_ij, rel_ij )
							
							multiply_factor = multiply_factor / rho_I(i) / r_ij
							
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
					
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij < beta%r_cutoff) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J
							rho_ij = 0.0
							call beta%spline_int (beta%r, beta%data_rho(jtype,:), &
									beta%y2rho(jtype,:), beta%n_points_rho, r_ij, rho_ij)
							
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							
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
							
							multiply_factor = multiply_factor1 - multiply_factor2

							forces_fric(1,i) = forces_fric(1,i) - multiply_factor * rel_ij(1) 
							forces_fric(2,i) = forces_fric(2,i) - multiply_factor * rel_ij(2) 
							forces_fric(3,i) = forces_fric(3,i) - multiply_factor * rel_ij(3)

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
			jneigh_start = (i-1)*n_neigh(i) + 1
			jneigh_stop = i*n_neigh(i)
			do jneigh = jneigh_start, jneigh_stop
				if (i /= neighbors_list(jneigh)) then
					j = neighbors_list(jneigh)
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi, xj,yj,zj)
					if (r_ij**2 < beta%r_cutoff_sq) then
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
						
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							
							! find rho_I_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(itype,:), &
								beta%y2rho_r_sq(itype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							! find alpha_J
							alpha_J = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(jtype,:), &
									beta%y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha_J)
							
							multiply_factor = dot_product( rel_ij, rand_vec(:,j) )
							multiply_factor = - alpha_J * rho_ij * multiply_factor / r_ij**2 / rho_I(j)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)

							! find alpha_I	
							alpha_I = 0.0
							call beta%spline_int (beta%rho, beta%data_alpha(itype,:), &
									beta%y2alpha(itype,:), beta%n_points_beta, rho_I(i), alpha_I)

							multiply_factor = dot_product( rel_ij, rand_vec(:,i) )
							multiply_factor = alpha_I * rho_ij * multiply_factor / r_ij**2 / rho_I(i)
							
							forces_rnd(1,i) = forces_rnd(1,i) + multiply_factor*rel_ij(1)
							forces_rnd(2,i) = forces_rnd(2,i) + multiply_factor*rel_ij(2)
							forces_rnd(3,i) = forces_rnd(3,i) + multiply_factor*rel_ij(3)
						end if
					end if
				end do
				
				v_Te = fdm%Collect_Te(xi,yi,zi)
				forces_rnd(1,i) = forces_rnd(1,i) * correl_factor_eta * sqrt(v_Te)
				forces_rnd(2,i) = forces_rnd(2,i) * correl_factor_eta * sqrt(v_Te)
				forces_rnd(3,i) = forces_rnd(3,i) * correl_factor_eta * sqrt(v_Te)
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							call relativeVector(vel(:,i), vel(:,j), rel_v_ij)
							
							rel_v_ij = rel_v_ij*1000.0		! A/fs ----> A/ps
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
				
				jneigh_start = (i-1)*n_neigh(i) + 1
				jneigh_stop = i*n_neigh(i)
				
				do jneigh = jneigh_start, jneigh_stop
					
					if (i /= neighbors_list(jneigh)) then
						
						j = neighbors_list(jneigh)
					
						xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
						r_ij = get_distance(xi,yi,zi,xj,yj,zj)
							
						if (r_ij**2 <= beta%r_cutoff_sq) then
							
							call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
							jtype = atom_type
							
							! find rho_J_r_sq
							rho_ij = 0.0
							call beta%spline_int (beta%r_sq, beta%data_rho_rsq(jtype,:), &
								beta%y2rho_r_sq(jtype,:), beta%n_points_rho, r_ij**2, rho_ij)
							
							call relativeVector(positions(:,j), positions(:,i), rel_ij)
							
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
							
							multiply_factor = multiply_factor1 - multiply_factor2
							
							forces_fric(1,i) = forces_fric(1,i) - multiply_factor * rel_ij(1) 
							forces_fric(2,i) = forces_fric(2,i) - multiply_factor * rel_ij(2) 
							forces_fric(3,i) = forces_fric(3,i) - multiply_factor * rel_ij(3)
						end if
					end if
				end do
			end do
					
		end if	!! for condition (isfriction == 1)
		
		end if	!! for condition (isfriction == 1 .or. israndom == 1)
	
	end if	!! for condition (model == 2)

	
	!! The current forces get modified
	
	!! Different cases of switching ON/OFF the friction and random forces
		 
	!! this is the full model
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
				
				!! ------ up to here -------
				
				forces(component,i) = forces(component,i) + forces_rnd(component,i) + forces_fric(component,i)
	
				E_val_i_fric = E_val_i_fric - dt*forces_fric(component,i)*vel(component,i)
				E_val_i_rnd = E_val_i_rnd - dt*forces_rnd(component,i)*vel(component,i)
			end do
			
			E_val_i_eph = E_val_i_fric + E_val_i_rnd
			
			if (fdm%fdm_option == 1) call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_eph, dt)
		
			Energy_val_fric = Energy_val_fric + E_val_i_fric
			Energy_val_rnd = Energy_val_rnd + E_val_i_rnd		
		end do
		Energy_val_eph = Energy_val_fric + Energy_val_rnd
	end if
	
	!! this is with only friction
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
	
	!! this is with only random	
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
	
	!! this is for not including friction and random, this is definitely not needed
	if (isfriction == 0 .and. israndom == 0) then
		Energy_val_fric = 0.0; Energy_val_rnd = 0.0; Energy_val_eph = 0.0
		do i = 1, Np
			do component = 1, 3
				forces(component,i) = forces(component,i)
			end do
		end do
	end if

	!! To calculate electronic temperature

	if (fdm%fdm_option == 1) call fdm%heat_diffusion_solve (dt)

	!! To write the electronic mesh temperatures

	if (T_out_mesh_freq /= 0) then
		if (MOD(md_istep, T_out_mesh_freq) == 0) then
			call fdm%save_output_to_file (T_outfile,md_istep,dt)
		end if
	end if

	!! To write the energy transfers within atomic-electonic system and electronic temperature

	if (MOD(md_istep, T_out_freq) == 0) then
		write(100,1) md_time, Energy_val_fric, Energy_val_rnd, Energy_val_eph, sum(fdm%T_e)/fdm%ntotal 
	end if
1	format(1E13.6,3E14.6,1E13.6)	
	close(unit=100)
	
	call fdm%beforeNextFeedback()
	
	end if 	! when md_time > 0.0

end subroutine eph_Langevin_spatial_correlation


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


!! Box-Muller transmorfation to get (nearly) Gaussian random variate
!! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)
real function randomGaussian() result(randg)
	implicit none
	real*8 :: n1, n2, R, theta
	real*8, parameter :: twoPi = 2.0*3.14285714
	call random_number(n1)
	call random_number(n2)
	R = sqrt(-2.0*log(n1)/n1)
	theta = twoPi*n2
	randg = R * cos(theta)
end function randomGaussian


!! Box-Muller transmorfation to get (nearly) an array of Gaussian random variates
!! (Ref. Numerical recipes in F77, vol. 1, W.H. Press et al.)
subroutine randomGaussianArray(Np, mean, standard_deviation, rand_array)
	implicit none
	integer, intent(in) :: Np
	real*8, intent(in) :: mean, standard_deviation
	real*8, intent(out) :: rand_array(Np)
	integer :: i
	real*8 :: one_value, mean_actual, sd_actual
	real*8 :: Pi = 4.0*atan(1.0)
	
	!! Uniformly distributed random numbers from 0 to 1
	call random_number(rand_array) 
	
	!! Box-Muller transmorfation to normally distributed random numbers
	do i = 1, Np-1, 2
		one_value = standard_deviation * &
		sqrt(-2.0*log(rand_array(i))) * sin(2*Pi*rand_array(i+1)) + mean
		
		rand_array(i+1) = standard_deviation * &
		sqrt(-2.0*log(rand_array(i))) * cos(2*Pi*rand_array(i+1)) + mean
		
		rand_array(i) = one_value
	end do
	
	!! Mean and standard deviation must be brought to required values 
	!! This works very good when mean is 0
	!mean_actual = sum(rand_array)/Np
	!sd_actual = sqrt(sum((rand_array - mean_actual)**2)/Np)
	!do i = 1, Np
	!	rand_array(i) = rand_array(i) - mean_actual + mean
	!	rand_array(i) = rand_array(i) * standard_deviation/sd_actual
	!end do
end subroutine randomGaussianArray



!! There is a correlation factor between 'eta', the random forces
!! which involves the T_e at position R_I(x,y,z) of the atom and not on the
!! overall temperature of the so-called "electronic heat bath". So, the following
!! subroutine seems not not to be necessary. 

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

