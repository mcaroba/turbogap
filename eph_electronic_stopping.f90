! ------- option for radiation cascade simulation with electronic stopping 
! based on electron-phonon coupling model				
!********* by Uttiyoarnab Saha

!**************************************************************************
! The model for electron-phonon coupling for calculation of electronic stopping
! according to the Langevin dynamics with spatial correlations is implemented here.
! The reference papers are PRL 120, 185501 (2018); PRB 99, 174302 (2019);
! PRB 94, 024305 (2016); PRB 99, 174301 (2019).

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
	real*8, allocatable :: v_relative(:,:), y2rho(:,:), y2(:), y2beta(:,:), z2(:), &
	forces_fric(:,:), forces_rnd(:,:), y2alpha(:,:), w2(:)
	integer :: Np, i, j, k, itype, jtype, ktype, atom_type, component
	
	real*8 :: beta_rho, xi,yi,zi,xj,yj,zj,xk,yk,zk, r_ij, r_ik, alpha, alpha_i, alpha_k, &
	sig_I(3),B_IJ(3,3),eta_I(3),W_IJ(3,3),sum_outer_p(3,3),outer_p(3,3),rho_ij, &
	rho_ik, rho_ki, v_ij_x, v_ij_y, v_ij_z, v_Te, rel_ij(3), rel_ik(3), rel_v_ij(3), WTv(3), &
	var1, var2, dvar, ezi1, ezi2, st_dev, &
	Energy_val_eph,Energy_val_fric,Energy_val_rnd,E_val_i_fric,E_val_i_rnd,E_val_i_eph, &
	multiply_factor, r_v1, r_v2, r_ik_sq
	
	
	real*8, allocatable :: rand_vec(:,:), rho_I(:), w_i(:,:)
	real*8, parameter :: bignum = 1.1e30
	real*8, parameter :: boltzconst = 8.61733326E-5 
	
	Np = size(vel,2)
	allocate(v_relative(3,Np))
	
	! To write the electronic energy loss data to a file after evaluation at each time step

	if( md_istep == 0 .or. md_istep == -1 )then
		open (unit = 100, file = "eph-EnergySharingData.txt", status = "unknown")
		write(100,*) 'Time (fs)  energy_fric(eV) energy_rand(eV) energy_net(eV) T_e (K)' 
		else
			open (unit = 100, file = "eph-EnergySharingData.txt", status = "old", position = "append")
	end if
	
	! to generate random forces
	allocate(rand_vec(3,Np))
	do i = 1, Np
		xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
		v_Te = fdm%get_T(xi,yi,zi)
		st_dev = sqrt(2*boltzconst*v_Te/(dt))
		rand_vec(1,i) = random_gaussian(0d1,st_dev)
		rand_vec(2,i) = random_gaussian(0d1,st_dev)
		rand_vec(3,i) = random_gaussian(0d1,st_dev)
	end do
	
	! have the y"(x) in y2rho and y2beta for cubic spline interpolation 
	! for all types of atoms
	
	allocate(y2rho(natomtypes,beta%n_points_rho), y2(beta%n_points_rho))
	
	allocate(y2beta(natomtypes,beta%n_points_beta), &
	z2(beta%n_points_beta),	y2alpha(natomtypes,beta%n_points_beta), w2(beta%n_points_beta))

	do i = 1, natomtypes
		y2 = 0
		call spline (beta%r,beta%data_rho(i,:),beta%n_points_rho,bignum,bignum,y2)
		y2rho(i,:) = y2(:)
		z2 = 0
		call spline (beta%rho,beta%data_beta(i,:),beta%n_points_beta,bignum,bignum,z2)
		y2beta(i,:) = z2(:)
		w2 = 0
		call spline (beta%rho,sqrt(beta%data_beta(i,:)),beta%n_points_beta,bignum,bignum,w2)
		y2alpha(i,:) = w2(:)
	end do
	
	allocate(forces_fric(3,Np), forces_rnd(3,Np))
	
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
					call spline_int (beta%r,beta%data_rho(jtype,:), &
					y2rho(jtype,:),beta%n_points_rho,r_ij,rho_ij)
					
					rho_I(i) = rho_I(i) + rho_ij
				end if
			end if
		end do
	end do
	
	
	forces_fric = 0.0
	forces_rnd = 0.0
	
	! friction forces according to (PRB 94, 024305 (2016))
	! Have to decide on this .........
	
	!do i = 1, Np  
	!	call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
	!	itype = atom_type
	!	xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
	!	do j = 1, Np
	!		if (i /= j) then
	!			xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
	!			r_ij = get_distance(xi,yi,zi,xj,yj,zj)
	!			if (r_ij <= beta%r_cutoff) then
	!				call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
	!				jtype = atom_type
	!				rho_ij = 0.0
	!				call spline_int (beta%r,beta%data_rho(jtype,:), &
	!				y2rho(jtype,:),beta%n_points_rho,r_ij,rho_ij)
	!				
	!				! Relative velocities for friction forces
	!				call relative_vector(vel(:,j), vel(:,i), rel_v_ij)
	!				v_relative(1,i) = v_relative(1,i) + rho_ij*rel_v_ij(1)
	!				v_relative(2,i) = v_relative(2,i) + rho_ij*rel_v_ij(2)
	!				v_relative(3,i) = v_relative(3,i) + rho_ij*rel_v_ij(3)
	!			end if
	!		end if
	!	end do
	!	
	!	! Friction forces depending on relative electron densities and velocities
	!	v_relative(1,i) = v_relative(1,i)/rho_I(i)
	!	v_relative(2,i) = v_relative(2,i)/rho_I(i)
	!	v_relative(3,i) = v_relative(3,i)/rho_I(i)
	!	
	!	call spline_int (beta%rho,beta%data_beta(itype,:),y2beta(itype,:),beta%n_points_beta,rho_I(i),beta_rho)
	!	
	!	forces_fric(1,i) = forces_fric(1,i) - beta_rho*v_relative(1,i)
	!	forces_fric(2,i) = forces_fric(2,i) - beta_rho*v_relative(2,i)
	!	forces_fric(3,i) = forces_fric(3,i) - beta_rho*v_relative(3,i)
	!end do
	
	
	
	! .......... Model 1 ..........
	! friction and random forces (PRL 120, 185501 (2018))
	
	allocate(w_i(3,Np))
	
	if (model_eph == 1) then
		
		! friction forces
		! find the w_i vectors for each atom
		
		do i = 1, Np
			
			call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
			itype = atom_type
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			
			! get alpha_i
			alpha_i = 0.0
			call spline_int (beta%rho,sqrt(beta%data_beta(itype,:)), &
					y2alpha(itype,:),beta%n_points_beta,rho_I(i),alpha_i)
			
			do k = 1, Np
				if (i /= k) then
					xk = positions(1,k); yk = positions(2,k); zk = positions(3,k)
					r_ik = get_distance(xi,yi,zi,xk,yk,zk)
		
					if (r_ik <= beta%r_cutoff) then
						call get_atom_type(k,natomtypes,masses,type_mass,atom_type)
						ktype = atom_type
						r_ik_sq = r_ik*r_ik
						
						! first sum
						rho_ki = 0
						call spline_int (beta%r,beta%data_rho(ktype,:), &
							y2rho(ktype,:),beta%n_points_rho,r_ik,rho_ki)

						multiply_factor = alpha_i * rho_ki / (rho_I(i) * r_ik_sq)
			
						call relative_vector(positions(:,k), positions(:,i), rel_ik)
		
						r_v1 = dot_product(rel_ik, vel(:,i))
						var1 = multiply_factor * r_v1
						r_v2 = dot_product(rel_ik, vel(:,k))
						var2 = multiply_factor * r_v2

						dvar = var1 - var2

						w_i(1,i) = w_i(1,i) + dvar * rel_ik(1)
						w_i(2,i) = w_i(2,i) + dvar * rel_ik(2)
						w_i(3,i) = w_i(3,i) + dvar * rel_ik(3)
					end if
				end if
			end do
		end do

		! calculate the forces by using the w_i vector arrays
		! f_i = W_ij w_j

		do i = 1, Np
	
			call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
			itype = atom_type
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			
			! get alpha_i
			alpha_i = 0.0
			call spline_int (beta%rho,sqrt(beta%data_beta(itype,:)), &
							y2alpha(itype,:),beta%n_points_beta,rho_I(i),alpha_i)
		
			do k = 1, Np
				if (i /= k) then
					xk = positions(1,k); yk = positions(2,k); zk = positions(3,k)
					r_ik = get_distance(xi,yi,zi,xk,yk,zk)
	
					if (r_ik <= beta%r_cutoff) then
						call get_atom_type(k,natomtypes,masses,type_mass,atom_type)
						
						ktype = atom_type
						r_ik_sq = r_ik*r_ik
						
						! get alpha_k
						alpha_k = 0.0
						call spline_int (beta%rho,sqrt(beta%data_beta(ktype,:)), &
								y2alpha(ktype,:),beta%n_points_beta,rho_I(k),alpha_k)
			
						rho_ki = 0
						call spline_int (beta%r,beta%data_rho(ktype,:), &
									y2rho(ktype,:),beta%n_points_rho,r_ik,rho_ki)
						
						call relative_vector(positions(:,k), positions(:,i), rel_ik)
	
						r_v1 = dot_product(rel_ik, w_i(:,i))
						var1 = alpha_i * rho_ki * r_v1 / (rho_I(i) * r_ik_sq)
	
						rho_ik = 0
						call spline_int (beta%r,beta%data_rho(itype,:), &
									y2rho(itype,:),beta%n_points_rho,r_ik,rho_ik)
	
						r_v2 = dot_product(rel_ik, w_i(:,k))
	
						var2 = alpha_k * rho_ik * r_v2 / (rho_I(k) * r_ik_sq)
	
						dvar = var1 - var2
	
						! frictional forces are negative
						
						forces_fric(1,i) = forces_fric(1,i) - dvar * rel_ik(1)
						forces_fric(2,i) = forces_fric(2,i) - dvar * rel_ik(2)
						forces_fric(3,i) = forces_fric(3,i) - dvar * rel_ik(3)
					end if
				end if
			end do
		end do

		! random forces
	
		do i = 1, Np
			call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
			itype = atom_type
			
			alpha_i = 0.0
			! get alpha_i
			call spline_int (beta%rho,sqrt(beta%data_beta(itype,:)), &
			y2alpha(itype,:),beta%n_points_beta,rho_I(i),alpha_i)
			
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			
			do k = 1, Np
				
				if (i /= k) then
					xk = positions(1,k); yk = positions(2,k); zk = positions(3,k)  
					r_ik = get_distance(xi,yi,zi,xk,yk,zk)
					r_ik_sq = r_ik*r_ik
					if (r_ik <= beta%r_cutoff) then
						call relative_vector(positions(:,k), positions(:,i), rel_ik)
						
						call get_atom_type(k,natomtypes,masses,type_mass,atom_type)
						ktype = atom_type
						
						alpha_k = 0.0
						! get alpha_k
						call spline_int (beta%rho,sqrt(beta%data_beta(ktype,:)), &
								y2alpha(ktype,:),beta%n_points_beta,rho_I(k),alpha_k)
						
						rho_ki = 0
						
						call spline_int (beta%r,beta%data_rho(ktype,:), &
								y2rho(ktype,:),beta%n_points_rho,r_ik,rho_ki)
						
						ezi1 = dot_product(rel_ik, rand_vec(:,i))
						
						var1 = alpha_i*rho_ki*ezi1/(r_ik_sq * rho_I(i))		! r_ik_sq
						
						rho_ik = 0
						
						call spline_int (beta%r,beta%data_rho(itype,:), &
						y2rho(itype,:),beta%n_points_rho,r_ik,rho_ik)
						
						ezi2 = dot_product(rel_ik, rand_vec(:,k))
						
						var2 = alpha_k*rho_ik*ezi2/(r_ik_sq * rho_I(k))		! r_ik_sq
						
						dvar = var1 - var2
						
						forces_rnd(1,i) = forces_rnd(1,i) + dvar * rel_ik(1)
						forces_rnd(2,i) = forces_rnd(2,i) + dvar * rel_ik(2)
						forces_rnd(3,i) = forces_rnd(3,i) + dvar * rel_ik(3)
					end if
				end if
			end do
		end do
	end if


	! .......... Model 2 ..........
	! Implementation of random forces and friction forces according to (PRB 99, 174302 (2019))
	! [supposed to be similar to that in PRL 120, 185501 (2018)]
	
	if (model_eph == 2) then
	
		do i = 1, Np
			xi = positions(1,i); yi = positions(2,i); zi = positions(3,i)
			do j = 1, Np
				! (I not = J)
				if (j /= i) then
					xj = positions(1,j); yj = positions(2,j); zj = positions(3,j) 
					r_ij = get_distance(xi,yi,zi,xj,yj,zj)
					if (r_ij <= beta%r_cutoff) then
						call get_atom_type(j,natomtypes,masses,type_mass,atom_type)
						jtype = atom_type
						rho_ij = 0
						call spline_int (beta%r,beta%data_rho(jtype,:),y2rho(jtype,:), &
						beta%n_points_rho,r_ij,rho_ij)
						
						call relative_vector(positions(:,j), positions(:,i), rel_ij)
						
						outer_p = 0
						call get_outer_product(rel_ij, outer_p)
						outer_p = outer_p * rho_ij / r_ij**2
						
						alpha = 0.0
						call spline_int (beta%rho, sqrt(beta%data_beta(jtype,:)), &
						y2alpha(jtype,:), beta%n_points_beta, rho_I(j), alpha)
						
						W_IJ = - alpha*outer_p/rho_I(j)
						WTv = 0.0
						call get_sigma_I(transpose(W_IJ), vel(:,j), WTv)	! fluctuation-dissipation
						sig_I = 0.0
						call get_sigma_I(W_IJ, WTv, sig_I)
	
						forces_fric(1,i) = forces_fric(1,i) - sig_I(1) 
						forces_fric(2,i) = forces_fric(2,i) - sig_I(2) 
						forces_fric(3,i) = forces_fric(3,i) - sig_I(3)
						
						eta_I = 0.0
						call get_eta_I(W_IJ,rand_vec(:,j),eta_I)
						forces_rnd(1,i) = forces_rnd(1,i) + eta_I(1)
						forces_rnd(2,i) = forces_rnd(2,i) + eta_I(2)
						forces_rnd(3,i) = forces_rnd(3,i) + eta_I(3)
					end if
				end if
				
				! (I = J)
				! since i and j are same we find everything for i th. atom 
				! in the following (as per the notation in references)
				
				if (j == i) then
					sum_outer_p = 0.0
					do k = 1, Np
						if (k /= i) then
							xk = positions(1,k); yk = positions(2,k); zk = positions(3,k) 
							r_ik = get_distance(xi,yi,zi,xk,yk,zk)
							if (r_ik <= beta%r_cutoff) then
								call get_atom_type(k,natomtypes,masses,type_mass,atom_type)
								ktype = atom_type
								rho_ik = 0
								call spline_int (beta%r,beta%data_rho(ktype,:), &
								y2rho(ktype,:),beta%n_points_rho,r_ik,rho_ik)
								
								call relative_vector(positions(:,k), positions(:,i), rel_ik)
								
								outer_p = 0
								call get_outer_product(rel_ik, outer_p)
								outer_p = outer_p * rho_ik / r_ik**2
								sum_outer_p = sum_outer_p + outer_p
							end if
						end if
					end do
					call get_atom_type(i,natomtypes,masses,type_mass,atom_type)
					itype = atom_type
					alpha = 0.0
					call spline_int (beta%rho,sqrt(beta%data_beta(itype,:)), &
					y2alpha(itype,:),beta%n_points_beta,rho_I(i),alpha)
					
					W_IJ = alpha*sum_outer_p/rho_I(i)
					
					! fluctuation-dissipation
					WTv = 0.0
					call get_sigma_I(transpose(W_IJ), vel(:,i), WTv)	
					sig_I = 0.0
					call get_sigma_I(W_IJ, WTv, sig_I)
				
					forces_fric(1,i) = forces_fric(1,i) - sig_I(1) 
					forces_fric(2,i) = forces_fric(2,i) - sig_I(2) 
					forces_fric(3,i) = forces_fric(3,i) - sig_I(3)
					
					eta_I = 0.0
					call get_eta_I(W_IJ,rand_vec(:,i),eta_I)
					forces_rnd(1,i) = forces_rnd(1,i) + eta_I(1)
					forces_rnd(2,i) = forces_rnd(2,i) + eta_I(2)
					forces_rnd(3,i) = forces_rnd(3,i) + eta_I(3)
				end if
			end do
		end do
	end if


	! The current forces get modified 
	Energy_val_eph = 0.0; Energy_val_fric = 0.0; Energy_val_rnd = 0.0
	do i = 1, Np
		E_val_i_fric = 0.0; E_val_i_rnd = 0.0; E_val_i_eph = 0.0
		xi = positions(1,i); yi = positions(2,i); zi = positions(3,i) 
		
		! Different cases of switching ON/OFF the friction and random forces
		 
		! this is the full model
		if (isfriction == 1 .and. israndom == 1) then		
			
			do component = 1, 3
				
				! ---- forces are adjusted by friction according to the relative signs of forces and velocities ----
				
				! when velocity and force are opposite, friction acts along the force
				!if ((vel(component,i)>0 .and. forces(component,i)<0).or. &
				!(vel(component,i)<0 .and. forces(component,i)>0)) then
				!	if (forces(component,i)<0) then
				!		forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
				!	end if
				!	if (forces(component,i)>0) then
				!		forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
				!	end if
				!end if
				!! when velocity and force are in same direction, friction acts opposite to the force
				!if ((vel(component,i)>0 .and. forces(component,i)>0).or. &
				!(vel(component,i)<0 .and. forces(component,i)<0)) then
				!	if (forces(component,i)<0) then
				!		forces(component,i) = forces(component,i) + abs(forces_fric(component,i))
				!	end if
				!	if (forces(component,i)>0) then 
				!		forces(component,i) = forces(component,i) - abs(forces_fric(component,i))
				!	end if
				!end if
				
				! ------ up to here -------
				
				forces(component,i) = forces(component,i) + &
										forces_rnd(component,i) + &
											forces_fric(component,i)

				E_val_i_fric = E_val_i_fric - dt*forces_fric(component,i)*vel(component,i)
				E_val_i_rnd = E_val_i_rnd - dt*forces_rnd(component,i)*vel(component,i)
			end do

			!E_val_i_fric = E_val_i_fric - dt*sqrt((forces_fric(1,i)**2 + &
			!forces_fric(2,i)**2 + forces_fric(3,i)**2)*(vel(1,i)**2 + vel(2,i)**2 + &
			!vel(3,i)**2))
						
			!E_val_i_rnd = E_val_i_rnd + dt*(forces_rnd(1,i)*vel(1,i) + &
			!forces_rnd(2,i)*vel(2,i) + forces_rnd(3,i)*vel(3,i))
			
			E_val_i_eph = E_val_i_eph + E_val_i_fric + E_val_i_rnd
			
			call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_eph, dt)
		end if
		
		if (isfriction == 1 .and. israndom == 0) then			
			
			do component = 1, 3
				forces(component,i) = forces(component,i) + forces_fric(component,i)
				E_val_i_fric = E_val_i_fric - dt*forces_fric(component,i)*vel(component,i)
			end do
			
			!E_val_i_fric = E_val_i_fric - dt*sqrt((forces_fric(1,i)**2 + &
			!forces_fric(2,i)**2 + forces_fric(3,i)**2)*(vel(1,i)**2 + vel(2,i)**2 + &
			!vel(3,i)**2))
			
			call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_fric, dt)
		end if
		
		if (isfriction == 0 .and. israndom == 1) then			
			
			do component = 1, 3
				forces(component,i) = forces(component,i) + forces_rnd(component,i)
				E_val_i_rnd = E_val_i_rnd - dt*forces_rnd(component,i)*vel(component,i)
			end do

			!E_val_i_rnd = E_val_i_rnd + dt*(forces_rnd(1,i)*vel(1,i) + &
			!forces_rnd(2,i)*vel(2,i) + forces_rnd(3,i)*vel(3,i))
			
			call fdm%feedback_ei_energy(xi, yi, zi, E_val_i_rnd, dt)
		end if
		
		if (isfriction == 0 .and. israndom == 0) then
			do component = 1, 3
				forces(component,i) = forces(component,i)
			end do
		end if

		Energy_val_fric = Energy_val_fric + E_val_i_fric
		Energy_val_rnd = Energy_val_rnd + E_val_i_rnd
		Energy_val_eph = Energy_val_eph + E_val_i_eph
		
	end do

	! To calculate electronic temperature

	if (fdm%fdm_option == 1) call fdm%heat_diffusion_solve (dt)
	
	! To write the electronic mesh temperatures
	
	if (T_out_mesh_freq /= 0) then
		if (MOD(md_istep, T_out_mesh_freq) == 0) then
			call fdm%save_output_to_file (T_outfile,md_istep,dt)
		end if
	end if
	
	call fdm%before_next_feedback()
	
	! To write the energy transfers from atomic system and electronic temperature
	
	if (MOD(md_istep, T_out_freq) == 0) then
		write(100,1) md_time, Energy_val_fric, Energy_val_rnd, Energy_val_eph, sum(fdm%T_e)/fdm%ntotal 
	end if
1	format(1E13.6,3E14.6,1E13.6)	

	close(unit=100)

end subroutine eph_Langevin_spatial_correlation


! Find distance between two atoms.
real function get_distance(xi,yi,zi,xj,yj,zj) result (r_ij)
	implicit none
	real*8, intent(in) :: xi,yi,zi,xj,yj,zj
	r_ij = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5
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
subroutine get_outer_product(vector,outer_p)
	implicit none
	real*8, intent(in) :: vector(3)
	real*8, intent(out) :: outer_p(3,3)
	real*8 :: vec(3,1), vecT(1,3)
	integer :: l, m
	
	vec(1,1) = vector(1); vec(2,1) = vector(2); vec(3,1) = vector(3)
	vecT(1,1) = vector(1); vecT(1,2) = vector(2); vecT(1,3) = vector(3) 
	do l = 1, 3
		do m = 1, 3
			outer_p(l,m) = vec(l,1) * vecT(1,m)
		end do
	end do
end subroutine get_outer_product


! Find eta_I, which are the random forces on the atoms. It is
! found using the projection of the random vectors associated with
! every atom on the (outer product of) relative unit vectors
! of the atoms.
subroutine get_eta_I(W_ij,vector,eta_I)
	implicit none
	real*8, intent(in) :: W_ij(3,3), vector(3)
	real*8, intent(out) :: eta_I(3)
	integer :: i, j
	
	eta_I = 0.0
	do i = 1, 3
		do j = 1, 3
			eta_I(i) = eta_I(i) + W_ij(i,j) * vector(j)
		end do
	end do
end subroutine get_eta_I

! Find sig_I, which are the friction forces on the atoms. It is
! found using the projection of B_ij associated with W_ij by fluctuation-
! dissipation theorem on the relative velocities of the atoms.
subroutine get_sigma_I(B_ij,vector,sig_I)
	implicit none
	real*8, intent(in) :: B_ij(3,3), vector(3)
	real*8, intent(out) :: sig_I(3)
	integer :: i, j
	
	sig_I = 0.0
	do i = 1, 3
		do j = 1, 3
			sig_I(i) = sig_I(i) + B_ij(i,j) * vector(j)
		end do
	end do
end subroutine get_sigma_I


! Box-Muller transmorfation to get Gaussian random variate
real function random_gaussian(mean, sigma) result(randg)
	implicit none
	real*8 :: n1, n2, R, theta, mean, sigma
	real*8, parameter :: twoPi = 2*3.14285714
	call random_number(n1)
	call random_number(n2)
	R = sqrt(-2*log(n1)/n1)
	theta = twoPi*n2
	randg = sigma * R * cos(theta) + mean
end function random_gaussian

! Find the interpolated value of y corresponding to a given value of x.
! Data arrays are xarr(1:n) and yarr(1:n). Second derivative of ya is y2arr(1:n).
! For a given value of x, y is the cubic-spline interpolated value.
subroutine spline_int(xarr,yarr,y2arr,n,x,y)
	implicit none
	integer :: n
	real*8 :: x, y, xarr(n), yarr(n), y2arr(n)
	integer :: indx, hiindx, loindx
	real*8 :: xwidth, A, B, C, D
	
	! Interpolate by method of bisection. Find the index limits within which x lies.
	
	loindx = 1
	hiindx = n
	
	do while (hiindx - loindx > 1)
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
subroutine spline(x,y,n,yp1,ypn,y2)
	implicit none
	integer :: n
	real*8 :: yp1, ypn, x(n), y(n), y2(n)
	integer :: i, k
	real*8 :: p, qn, sig, un, u(n)
	
	!The lower boundary condition is set either to be “natural”, else to have a specified first derivative.
	if (yp1 > 0.99e30) then
		y2(1) = 0.0
		u(1) = 0.0
	else
		y2(1) = -0.5
		u(1) = (3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	end if
	
	!This is the decomposition loop of the tridiagonal
	!algorithm. y2 and u are used for temporary
	!storage of the decomposed factors.
	do i = 2, n-1
		sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
		p = sig*y2(i-1) + 2.0
		y2(i) = (sig-1.0)/p
		u(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/ &
		(x(i) - x(i-1)))/(x(i+1) - x(i-1))-sig*u(i-1))/p
	end do
	
	!The upper boundary condition is set either to be “natural”, else to have a specified first derivative.
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

end module eph_electronic_stopping

