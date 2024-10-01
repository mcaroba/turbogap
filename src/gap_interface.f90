! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, gap_interface.f90, is copyright (c) 2019-2022, Miguel A. Caro
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

module gap_interface

  use soap_turbo_desc
  use gap
  use read_files
  use types
  use local_prop
  use F_B_C
  use iso_c_binding
  use mpi


  contains





!**************************************************************************
  subroutine get_gap_soap(n_sparse, n_total_sites, n_sites0, n_neigh0, neighbors_list0, n_species, species_types, &
                          rjs0, thetas0, phis0, xyz0, alpha_max_d, alpha_max, l_max, n_soap, &
                          rcut_hard_d, rcut_hard, rcut_soft_d, nf_d, global_scaling_d, atom_sigma_r_d, &
                          atom_sigma_r, atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
                          amplitude_scaling_d, radial_enhancement, central_weight_d, central_weight, basis, &
                          scaling_mode, do_timing, do_derivatives, do_forces, do_prediction, &
                          write_soap, write_derivatives, compress_soap, &
                          compress_soap_indices, &
                          delta, zeta, central_species, &
                          xyz_species, xyz_species_supercell, alphas_d, Qs_d, all_atoms, &
                          which_atom, indices, soap, soap_cart_der, der_neighbors, der_neighbors_list, &
                          has_local_properties, n_local_properties, local_property_models, lp_index,&
                          energies0, forces0, local_properties0,&
                          & local_properties_cart_der0,&
                          & local_property_indexes,  virial, solo_time_soap, &
                          time_get_soap, W_d, S_d, multiplicity_array_d,&
                          & st_W_d, st_S_d, st_multiplicity_array_d,&
                          & recompute_basis,cublas_handle , gpu_stream)

    implicit none


    !   Input variables
    type(local_property_soap_turbo), allocatable :: local_property_models(:)

!   real*8, intent(in) :: rjs0(:), thetas0(:), phis0(:), xyz0(:,:), rcut_hard(:), rcut_soft(:), &
    real*8, intent(in) :: rjs0(:), thetas0(:), phis0(:), xyz0(:,:), rcut_hard(:), &
!                         nf(:), global_scaling(:), atom_sigma_r(:), atom_sigma_r_scaling(:), &
                          atom_sigma_r(:), &
!                         atom_sigma_t(:), atom_sigma_t_scaling(:), amplitude_scaling(:), &
!                         central_weight(:), delta, zeta, Qs(:,:), vdw_Qs(:,:), &
                          central_weight(:), delta, zeta
    integer, intent(in) :: n_sites0, n_neigh0(:), neighbors_list0(:), n_species, central_species, &
                           radial_enhancement, compress_soap_indices(:), which_atom, &
                           indices(1:3), alpha_max(:), l_max, n_total_sites, n_sparse,&
                           & n_local_properties,&
                           & local_property_indexes(:), lp_index

    logical, intent(in) :: do_timing, do_derivatives, compress_soap, do_forces, do_prediction, &
                           all_atoms, write_soap, write_derivatives,  has_local_properties

    character*64, intent(in) :: basis
    character*32, intent(in) :: scaling_mode
    character*8, intent(in) :: xyz_species(:), xyz_species_supercell(:), species_types(:)

!   Output variables
    real*8, allocatable, intent(out) :: soap(:,:), soap_cart_der(:,:,:)
    real*8, intent(out) :: virial(1:3, 1:3)
    real*8, intent(inout)  :: solo_time_soap, time_get_soap

!   Inout variables
    real*8, intent(inout) :: energies0(:), forces0(:,:), local_properties0(:,:), local_properties_cart_der0(:,:,:)

!   Internal variables
    real*8, allocatable :: rjs(:), thetas(:), phis(:), energies(:), forces(:,:), soap_temp(:,:), &
                           local_properties(:), local_properties_cart_der(:,:), xyz(:,:)

    integer, allocatable :: in_to_out_site(:), n_neigh(:), neighbors_list(:), species_multiplicity(:), &
                            species(:,:), species0(:,:), species_multiplicity0(:), out_to_in_site(:), &
                            der_neighbors(:), der_neighbors_list(:), in_to_out_pairs(:)
    integer, allocatable :: species_multiplicity_supercell(:)
    integer :: n_sites, i, j, n_atom_pairs, k, k2, i2, j2, n_soap, max_species_multiplicity, i3, n_all_sites, &
               i4, n_sites_supercell, j3
    integer :: ierr, rank
    logical, allocatable :: mask(:,:)
    logical, allocatable ::  mask0(:,:), is_atom_seen(:)
!   CLEAN THIS UP
    real*8 :: time1, time2, rcut_max
    real*8 :: ttt(2)
    type(c_ptr) :: soap_cart_der_d, soap_d, nf_d, rcut_hard_d, rcut_soft_d, global_scaling_d, atom_sigma_r_d, atom_sigma_r_scaling_d
    type(c_ptr) :: atom_sigma_t_d, atom_sigma_t_scaling_d, amplitude_scaling_d, alpha_max_d, central_weight_d, alphas_d,Qs_d
    type(c_ptr) :: cublas_handle, gpu_stream
    type(c_ptr) :: n_neigh_d

    type(c_ptr), intent(inout) :: W_d, S_d, multiplicity_array_d
    integer(c_size_t), intent(inout) :: st_W_d, st_S_d, st_multiplicity_array_d
    logical, intent(inout) ::  recompute_basis
    
    ! For local properties
    type(c_ptr) :: alphas_lp_d, Qs_lp_d, l_index_d, local_properties_d, local_properties_cart_der_d
    integer(c_int) :: n_pairs

    ! Check for lp prediction
    integer :: n_sparse_lp, dim_lp
    integer(c_size_t) :: st_soap, st_soap_cart_der, st_size_nf



    
    !--- TODO: CHANGE THE LOCAL PROPERTIES PARAMS TO DOUBLES TOO ---!

    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! call gpu_set_device(rank) ! Every node has 4 GPUs. Even if there are more than 1 nodes used. This will assing the ranks to GPU in a roundbin fashion
    
    !call create_cublas_handle(cublas_handle, gpu_stream)
    n_sites_supercell = size(xyz_species_supercell)


!   Work out the multiplicity stuff
!   Check if the user wants to represent one or more of the species more than once
    max_species_multiplicity = 0
    do i = 1, n_species
      k = 1
      do j = i+1, n_species
        if( species_types(i) == species_types(j) )then
          k = k + 1
        end if
      end do
      if( k > max_species_multiplicity )then
        max_species_multiplicity = k
      end if
    end do
    n_atom_pairs = size(rjs0)
    call assign_species_multiplicity(max_species_multiplicity, species_types, xyz_species, &
                                     xyz_species_supercell, n_species, all_atoms, which_atom, &
                                     indices, n_neigh0, n_atom_pairs, neighbors_list0, mask0, &
                                     species0, species_multiplicity0, &
                                     species_multiplicity_supercell )

!   We need to sort out which atoms from the general neighbors list we actually keep
!
!   First, we just count the number of sites to keep and the number of atom pairs for allocation purposes
    n_sites = 0
!   n_all_sites are all the sites (central or else) to which this descriptor is not blind (including neighbors)
!   and on which forces might be acting
    n_all_sites = 0
    n_atom_pairs = 0
    k = 0
    allocate( is_atom_seen(1:n_total_sites) )
    is_atom_seen = .false.
    rcut_max=maxval(rcut_hard)
    do i = 1, n_sites0
      if( species0(1, i) == central_species .or. central_species == 0 )then
        n_sites = n_sites + 1
        do j = 1, n_neigh0(i)
          k = k + 1
          j2 = mod(neighbors_list0(k)-1, n_total_sites)+1
          if( rjs0(k) < rcut_max .and. species_multiplicity_supercell(j2) > 0 )then
            n_atom_pairs = n_atom_pairs + 1
            is_atom_seen(j2) = .true.
          end if
        end do
      else
        k = k + n_neigh0(i)
      end if
    end do
    do i = 1, n_total_sites
      if( is_atom_seen(i) )then
        n_all_sites = n_all_sites + 1
      end if
    end do

!   Second, we build all the new neighbors list stuff with the correct subset of information
    allocate( n_neigh(1:n_sites) )
    allocate( species_multiplicity(1:n_sites) )
    allocate( species(1:max_species_multiplicity, 1:n_sites) )
    allocate( neighbors_list(1:n_atom_pairs) )
    allocate( rjs(1:n_atom_pairs) )
    allocate( thetas(1:n_atom_pairs) )
    allocate( phis(1:n_atom_pairs) )
    allocate( xyz(1:3, 1:n_atom_pairs) )
    allocate( mask(1:n_atom_pairs, 1:n_species) )
    allocate( in_to_out_site(1:n_all_sites) )
    allocate( in_to_out_pairs(1:n_atom_pairs) )
!    allocate( out_to_in_site(1:n_sites0) )
    allocate( out_to_in_site(1:n_total_sites) )
    in_to_out_site = 0
    in_to_out_pairs = 0
    out_to_in_site = 0
    k2 = 0
    i4 = 0
    k = 0
    i2 = 0
!   The first n_sites sites are assigned to the central sites
    do i = 1, n_sites0
      k = k + 1
      j = mod(neighbors_list0(k)-1, n_total_sites) + 1
      if( species0(1, i) == central_species .or. central_species == 0 )then
        i2 = i2 + 1
        in_to_out_site(i2) = j
        out_to_in_site(j) = i2
        species_multiplicity(i2) = species_multiplicity0(i)
        species(:, i2) = species0(:, i)
        k = k + n_neigh0(i) - 1
      else
        k = k + n_neigh0(i) - 1
      end if
    end do
    k = 0
    i2 = 0
!   Now we loop over neighbors and make sure that if a neighbor has already been registered, which
!   means its out_to_in_site /= 0, we do not update the out_to_in_site variable further
    do i = 1, n_sites0
      if( species0(1, i) == central_species .or. central_species == 0 )then
        i2 = i2 + 1
        j2 = 0
        do j = 1, n_neigh0(i)
          k = k + 1
!         We fold neighbors in the supercell back to the central unit cell
          j3 = mod(neighbors_list0(k)-1, n_total_sites)+1
          if( rjs0(k) <= rcut_max .and. species_multiplicity_supercell(j3) > 0 )then
            j2 = j2 + 1
            k2 = k2 + 1
            rjs(k2) = rjs0(k)
            phis(k2) = phis0(k)
            thetas(k2) = thetas0(k)
            xyz(1:3, k2) = xyz0(1:3, k)
            in_to_out_pairs(k2) = k
            do i3 = 1, n_species
              mask(k2, i3) = mask0(k, i3)
            end do
            if( j > 1 .and. out_to_in_site(j3) == 0 )then
              i4 = i4 + 1
              in_to_out_site(n_sites + i4) = j3
              out_to_in_site(j3) = n_sites + i4
              neighbors_list(k2) = out_to_in_site(j3)
            else if( j > 1 .and. out_to_in_site(j3) /= 0 )then
              neighbors_list(k2) = out_to_in_site(j3)
            else if( j == 1 )then
              neighbors_list(k2) = i2
            end if
          end if
        end do
        n_neigh(i2) = j2
      else
        k = k + n_neigh0(i)
      end if
    end do

!   Get SOAP vectors and derivatives:
    allocate( soap(1:n_soap, 1:n_sites) )
    soap = 0.d0
    st_soap=n_soap*n_sites*sizeof(soap(1,1))
    call gpu_malloc_all(soap_d, st_soap, gpu_stream)

    
    if( do_derivatives )then
      allocate( soap_cart_der(1:3, 1:n_soap, 1:n_atom_pairs) )
      soap_cart_der = 0.d0
      st_soap_cart_der =   3*n_soap*n_atom_pairs*sizeof(soap_cart_der(1,1,1))
      call gpu_malloc_all(soap_cart_der_d, st_soap_cart_der, gpu_stream) !call gpu_malloc_all_blocking(soap_cart_der_d, st_soap_cart_de      
    end if

    print *, " starting get soap "
    
    if( n_sites > 0 )then
!      call cpu_time(ttt(1))
      call get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                    thetas, phis, alpha_max_d, alpha_max, l_max, rcut_hard_d, rcut_hard, rcut_soft_d, nf_d, global_scaling_d, &
                    atom_sigma_r_d, atom_sigma_r, atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
                    amplitude_scaling_d, radial_enhancement, central_weight_d, central_weight, basis, scaling_mode, do_timing, &
                    do_derivatives, compress_soap, compress_soap_indices, soap,&
                    & soap_cart_der, time_get_soap, soap_d, &
                    & soap_cart_der_d, n_neigh_d, W_d, S_d, multiplicity_array_d,&
                    & st_W_d, st_S_d, st_multiplicity_array_d,&
                    & recompute_basis, &
                    & cublas_handle, gpu_stream)


!      call cpu_time(ttt(2))

    end if
   ! write(*,*) "n_total_sites", n_total_sites, "n_sites", n_sites, "nsites0", n_sites0
   !time_get_soap=time_get_soap+ttt(2)-ttt(1)






    if( do_prediction )then
!     Get energies and forces
      allocate( energies(1:n_sites) )
      energies = 0.d0
      if( do_forces )then
        allocate( forces(1:3, 1:n_all_sites) )
        forces = 0.d0
      end if
      if( n_sites > 0 )then
        call get_soap_energy_and_forces(n_sparse, soap, soap_cart_der, alphas_d, delta, zeta, 0.d0, Qs_d, &
                                        n_neigh, neighbors_list, xyz, do_forces, do_timing, &
                                        energies, forces, virial, solo_time_soap,  &
                                        soap_d,  soap_cart_der_d, n_neigh_d, n_pairs, l_index_d,  cublas_handle, gpu_stream )
        
      end if

      do i = 1, n_sites
        i2 = in_to_out_site(i)
        energies0(i2) = energies(i)
      end do
      if( do_forces )then
        do i = 1, n_all_sites
          i2 = in_to_out_site(i)
          forces0(1:3, i2) = forces(1:3, i)
        end do
      end if
    end if



!     !###########################################!
!     !###---   Local property prediction   ---###!
!     !###########################################!

    
    
    print *, " starting local properties (cpu) "
    !--- TODO: CONVERT THIS FUNCTION INTO GPU KERNEL ---!
    if( has_local_properties )then
       !call cpu_time(time1)

       allocate( local_properties( 1:n_sites ) )
       if( do_derivatives )then
          allocate( local_properties_cart_der(1:3, 1:n_atom_pairs) )
       end if

       ! Allocate the arrays here so we dont have to reallocate 
       ! Allocate gpu memory
       st_size_nf = n_sites * sizeof( local_properties(1) )
       call gpu_malloc_all(local_properties_d, st_size_nf, gpu_stream)
       st_size_nf = n_atom_pairs * 3 * sizeof( local_properties_cart_der(1,1) )
       call gpu_malloc_all(local_properties_cart_der_d, st_size_nf, gpu_stream)
       
       
       ! We need to iterate over the number of local properties       
       do i4 = 1, n_local_properties

          call gpu_local_property_predict( local_property_models(i4)%n_sparse, soap,&
               & soap_d, local_property_models(i4)%Qs_d,&
               & local_property_models(i4)%alphas_d,&
               & local_property_models(i4)%V0,&
               & local_property_models(i4)%delta,&
               & local_property_models(i4)%zeta, local_properties, local_properties_d,&
               & do_derivatives, soap_cart_der_d,&
               & local_properties_cart_der,&
               & local_properties_cart_der_d, n_pairs, l_index_d ,&
               & cublas_handle, gpu_stream )


          do i = 1, n_sites
             i2 = in_to_out_site(i)
             local_properties0(i2, local_property_indexes(lp_index + i4)) = local_properties(i)
!             write(*,'(A,1X,I8,1X,A,1X,I8,1X,A,1X,F20.12)')  "rank ", rank, "lp, index  ", i2, " : ", local_properties0(i2, local_property_indexes(lp_index + i4))  
          end do
          if( do_derivatives )then
             do k = 1, n_atom_pairs
                k2 = in_to_out_pairs(k)
                local_properties_cart_der0(1:3,  k2, local_property_indexes(lp_index + i4)) &
                     & = local_properties_cart_der(1:3, k)
!                write(*,'(A,1X,I8,1X,A,1X,I8,1X,A,1X,F20.12,1X,F20.12,1X,F20.12)')  "rank ", rank,  "lpder, index  ", k2, " : ",&
                     ! & local_properties_cart_der(1, k),&
                     ! & local_properties_cart_der(2, k),&
                     ! & local_properties_cart_der(3, k)
             end do
          end if
       end do

       deallocate( local_properties )
       call gpu_free_async(local_properties_d, gpu_stream)
       if( do_derivatives )then
          deallocate( local_properties_cart_der )
          call gpu_free_async(local_properties_cart_der_d, gpu_stream)
       end if

    end if


    call gpu_free_async(l_index_d,gpu_stream)


!! commenting out to here 
    
!    call gpu_free_async(j2_index_d,gpu_stream) 
    call gpu_free_async(soap_d,gpu_stream)
    call gpu_free(soap_cart_der_d)

    
!   This is slow, but writing to disk is even slower so who cares
    if( write_soap )then
      allocate( soap_temp(1:n_soap, 1:n_sites0) )
      soap_temp = 0.d0
      do i = 1, n_sites
        i2 = in_to_out_site(i)
        soap_temp(1:n_soap, i2) = soap(1:n_soap, i)
      end do
      deallocate( soap )
      allocate( soap(1:n_soap, 1:n_sites0) )
      soap = soap_temp
      deallocate( soap_temp )
    end if
    if( do_derivatives .and. write_derivatives )then
      allocate( der_neighbors(1:n_sites) )
      allocate( der_neighbors_list(1:n_atom_pairs) )
      der_neighbors = n_neigh
      k = 0
      do i = 1, n_sites
        do j = 1, n_neigh(i)
          k = k + 1
          j2 = neighbors_list(k)
          der_neighbors_list(k) = in_to_out_site(j2)
        end do
      end do
    end if



    if( do_prediction )then
      deallocate( energies )
    end if
    if( do_forces )then
      deallocate( forces )
    end if

    deallocate( neighbors_list, n_neigh, rjs, thetas, phis, mask, mask0, is_atom_seen )
    deallocate( species0, species_multiplicity0, species, species_multiplicity, species_multiplicity_supercell )
    deallocate( in_to_out_site, out_to_in_site , in_to_out_pairs )
    if( .not. write_soap )then
      deallocate( soap )
    end if
    if( do_derivatives .and. .not. write_derivatives )then
      deallocate( soap_cart_der )
    end if

  !call destroy_cublas_handle(cublas_handle, gpu_stream)
end subroutine get_gap_soap
!**************************************************************************

    

! !**************************************************************************
!   subroutine get_gap_soap(n_total_sites, n_sites0, n_neigh0, neighbors_list0, n_species, species_types, &
!                           rjs0, thetas0, phis0, xyz0, alpha_max, l_max, n_soap, &
!                           rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
!                           atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, &
!                           amplitude_scaling, radial_enhancement, central_weight, basis, &
!                           scaling_mode, do_timing, do_derivatives, do_forces, do_prediction, &
!                           write_soap, write_derivatives, compress_soap, &
!                           compress_P_nonzero, compress_P_i, compress_P_j, compress_P_el, &
!                           delta, zeta, central_species, &
!                           xyz_species, xyz_species_supercell, alphas, Qs, all_atoms, &
!                           which_atom, indices, soap, soap_cart_der, der_neighbors, der_neighbors_list, &
!                           energies0, forces0, virial,&
!                           & has_local_properties, n_atom_pairs, in_to_out_pairs,&
!                           & n_all_sites, in_to_out_site, n_neigh_out,&
!                           & n_sites_out, this_j_beg, this_j_end,&
!                           & local_properties0,&
!                           & local_properties_cart_der0  )

!     implicit none

! !   Input variables
!     real*8, intent(in) :: rjs0(:), thetas0(:), phis0(:), xyz0(:,:), rcut_hard(:), rcut_soft(:), &
!                           nf(:), global_scaling(:), atom_sigma_r(:), atom_sigma_r_scaling(:), &
!                           atom_sigma_t(:), atom_sigma_t_scaling(:), amplitude_scaling(:), &
!                           central_weight(:), delta, zeta, alphas(:), Qs(:,:), compress_P_el(:)
!     integer, intent(in) :: n_sites0, n_neigh0(:), neighbors_list0(:), n_species, central_species, &
!                            radial_enhancement, which_atom, indices(1:3), alpha_max(:), l_max, &
!                            n_total_sites, compress_P_nonzero,&
!                            & compress_P_i(:), compress_P_j(:)
!     logical, intent(in) :: do_timing, do_derivatives, compress_soap,&
!          & do_forces, do_prediction, all_atoms, write_soap,&
!          & write_derivatives,  has_local_properties
!     character*64, intent(in) :: basis
!     character*32, intent(in) :: scaling_mode
!     character*8, intent(in) :: xyz_species(:), xyz_species_supercell(:), species_types(:)

! !   Output variables
!     real*8, allocatable, intent(out) :: soap(:,:), soap_cart_der(:,:,:)
!     real*8, intent(out) :: virial(1:3, 1:3)

! !   Inout variables
!     real*8, intent(inout) :: energies0(:), forces0(:,:), local_properties0(:), local_properties_cart_der0(:,:)

! !   Internal variables
!     real*8, allocatable :: rjs(:), thetas(:), phis(:), energies(:), forces(:,:), soap_temp(:,:), &
!                            & xyz(:,:)
!     real*8 :: rcut_max
!     integer, allocatable :: neighbors_list(:), species_multiplicity(:), &
!                             species(:,:), species0(:,:), species_multiplicity0(:), out_to_in_site(:), &
!                             der_neighbors(:), der_neighbors_list(:)
!     integer, allocatable, intent(out) :: in_to_out_pairs(:), in_to_out_site(:), n_neigh_out(:)
!     integer, allocatable :: species_multiplicity_supercell(:)
!     integer :: n_sites, i, j, k, k2, i2, j2, n_soap, max_species_multiplicity, i3, &
!          i4, n_sites_supercell, j3, this_j_beg, this_j_end
!     integer, intent(out) :: n_atom_pairs, n_all_sites, n_sites_out
!     logical, allocatable :: mask(:,:), mask0(:,:), is_atom_seen(:)


! !   CLEAN THIS UP
!     real*8 :: time1, time2

!     n_sites_supercell = size(xyz_species_supercell)


! !   Work out the multiplicity stuff
! !   Check if the user wants to represent one or more of the species more than once
!     max_species_multiplicity = 0
!     do i = 1, n_species
!       k = 1
!       do j = i+1, n_species
!         if( species_types(i) == species_types(j) )then
!           k = k + 1
!         end if
!       end do
!       if( k > max_species_multiplicity )then
!         max_species_multiplicity = k
!       end if
!     end do
!     n_atom_pairs = size(rjs0)
!     call assign_species_multiplicity(max_species_multiplicity, species_types, xyz_species, &
!                                      xyz_species_supercell, n_species, all_atoms, which_atom, &
!                                      indices, n_neigh0, n_atom_pairs, neighbors_list0, mask0, &
!                                      species0, species_multiplicity0, &
!                                      species_multiplicity_supercell )

! !   We need to sort out which atoms from the general neighbors list we actually keep
! !
! !   First, we just count the number of sites to keep and the number of atom pairs for allocation purposes
!     rcut_max = maxval(rcut_hard)
!     n_sites_out = 0
! !   n_all_sites are all the sites (central or else) to which this descriptor is not blind (including neighbors)
! !   and on which forces might be acting
!     n_all_sites = 0
!     n_atom_pairs = 0
!     k = 0
!     allocate( is_atom_seen(1:n_total_sites) )
!     is_atom_seen = .false.
!     do i = 1, n_sites0
!       if( species0(1, i) == central_species .or. central_species == 0 )then
!         n_sites_out = n_sites_out + 1
!         do j = 1, n_neigh0(i)
!           k = k + 1
!           j2 = mod(neighbors_list0(k)-1, n_total_sites)+1
!           if( rjs0(k) < rcut_max .and. species_multiplicity_supercell(j2) > 0 )then
!             n_atom_pairs = n_atom_pairs + 1
!             is_atom_seen(j2) = .true.
!           end if
!         end do
!       else
!         k = k + n_neigh0(i)
!       end if
!     end do
!     do i = 1, n_total_sites
!       if( is_atom_seen(i) )then
!         n_all_sites = n_all_sites + 1
!       end if
!     end do

! !   Second, we build all the new neighbors list stuff with the correct subset of information
!     allocate( n_neigh_out(1:n_sites_out) )
!     allocate( species_multiplicity(1:n_sites_out) )
!     allocate( species(1:max_species_multiplicity, 1:n_sites_out) )
!     allocate( neighbors_list(1:n_atom_pairs) )
!     allocate( rjs(1:n_atom_pairs) )
!     allocate( thetas(1:n_atom_pairs) )
!     allocate( phis(1:n_atom_pairs) )
!     allocate( xyz(1:3, 1:n_atom_pairs) )
!     allocate( mask(1:n_atom_pairs, 1:n_species) )
!     allocate( in_to_out_site(1:n_all_sites) )
!     allocate( in_to_out_pairs(1:n_atom_pairs) )
! !    allocate( out_to_in_site(1:n_sites_out0) )
!     allocate( out_to_in_site(1:n_total_sites) )
!     in_to_out_site = 0
!     in_to_out_pairs = 0
!     out_to_in_site = 0
!     k2 = 0
!     i4 = 0
!     k = 0
!     i2 = 0
! !   The first n_sites_out sites are assigned to the central sites
!     do i = 1, n_sites0
!       k = k + 1
!       j = mod(neighbors_list0(k)-1, n_total_sites) + 1
!       if( species0(1, i) == central_species .or. central_species == 0 )then
!         i2 = i2 + 1
!         in_to_out_site(i2) = j
!         out_to_in_site(j) = i2
!         species_multiplicity(i2) = species_multiplicity0(i)
!         species(:, i2) = species0(:, i)
!         k = k + n_neigh0(i) - 1
!       else
!         k = k + n_neigh0(i) - 1
!       end if
!     end do
!     k = 0
!     i2 = 0
! !   Now we loop over neighbors and make sure that if a neighbor has already been registered, which
! !   means its out_to_in_site /= 0, we do not update the out_to_in_site variable further
!     do i = 1, n_sites0
!       if( species0(1, i) == central_species .or. central_species == 0 )then
!         i2 = i2 + 1
!         j2 = 0
!         do j = 1, n_neigh0(i)
!           k = k + 1
! !         We fold neighbors in the supercell back to the central unit cell
!           j3 = mod(neighbors_list0(k)-1, n_total_sites)+1
!           if( rjs0(k) < rcut_max .and. species_multiplicity_supercell(j3) > 0 )then
!             j2 = j2 + 1
!             k2 = k2 + 1
!             rjs(k2) = rjs0(k)
!             phis(k2) = phis0(k)
!             thetas(k2) = thetas0(k)
!             xyz(1:3, k2) = xyz0(1:3, k)
!             in_to_out_pairs(k2) = k
!             do i3 = 1, n_species
!               mask(k2, i3) = mask0(k, i3)
!             end do
!             if( j > 1 .and. out_to_in_site(j3) == 0 )then
!               i4 = i4 + 1
!               in_to_out_site(n_sites_out + i4) = j3
!               out_to_in_site(j3) = n_sites_out + i4
!               neighbors_list(k2) = out_to_in_site(j3)
!             else if( j > 1 .and. out_to_in_site(j3) /= 0 )then
!               neighbors_list(k2) = out_to_in_site(j3)
!             else if( j == 1 )then
!               neighbors_list(k2) = i2
!             end if
!           end if
!         end do
!         n_neigh_out(i2) = j2
!       else
!         k = k + n_neigh0(i)
!       end if
!     end do

! !   Get SOAP vectors and derivatives:
!     allocate( soap(1:n_soap, 1:n_sites_out) )
!     soap = 0.d0
!     if( do_derivatives )then
!       allocate( soap_cart_der(1:3, 1:n_soap, 1:n_atom_pairs) )
!       soap_cart_der = 0.d0
!     end if

!     if( n_sites_out > 0 )then
!       call get_soap(n_sites_out, n_neigh_out, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
!                     thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, &
!                     atom_sigma_r, atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, &
!                     amplitude_scaling, radial_enhancement, central_weight, basis, scaling_mode, do_timing, &
!                     do_derivatives, compress_soap, compress_P_nonzero, compress_P_i, compress_P_j, &
!                     compress_P_el, soap, soap_cart_der)
!     end if

! !     if( has_vdw )then
! ! !call cpu_time(time1)
! !       allocate( hirshfeld_v( 1:n_sites_out ) )
! !       hirshfeld_v = 0.d0
! !       if( do_derivatives )then
! !         allocate( hirshfeld_v_cart_der(1:3, 1:n_atom_pairs) )
! !         hirshfeld_v_cart_der = 0.d0
! !       end if

! !       call hirshfeld_predict( soap, vdw_Qs, vdw_alphas, vdw_V0, vdw_delta, vdw_zeta, &
! !                               hirshfeld_v, do_derivatives, soap_cart_der, n_neigh_out, &
! !                               hirshfeld_v_cart_der )
! !       do i = 1, n_sites_out
! !         i2 = in_to_out_site(i)
! !         hirshfeld_v0(i2) = hirshfeld_v(i)
! !       end do
! !       if( do_derivatives )then
! !         do k = 1, n_atom_pairs
! !           k2 = in_to_out_pairs(k)
! !           hirshfeld_v_cart_der0(1:3, k2) = hirshfeld_v_cart_der(1:3, k)
! !         end do
! !       end if

! !       deallocate( hirshfeld_v )
! !       if( do_derivatives )then
! !         deallocate( hirshfeld_v_cart_der )
! !       end if
! ! !call cpu_time(time2)
! ! !write(*,*) "hirshfeld_v time =", time2-time1, "seconds"
! !     end if



!     if( has_local_properties )then
! !call cpu_time(time1)
!        allocate( local_properties( 1:n_sites_out, 1:n_local_properties ) )

!       do i = 1, n_local_properties
!          if(local_property_models(i)%do_derivatives)then
!             allocate( local_properties_cart_der(1:3, 1:n_sites_out, 1:n_local_properties ) )
!          end if
!       end do

!       do i = 1, n_local_properties
!          local_properties = 0.d0
!          if(local_property_models(i)%do_derivatives)then
!             local_properties_cart_der = 0.d0
!          end if

!          call local_property_predict( soap, &
!               & local_property_models(i)%Qs,&
!               & local_property_models(i)%alphas,&
!               & local_property_models(i)%V0,&
!               & local_property_models(i)%delta,&
!               & local_property_models(i)%zeta,&
!               & local_properties(i,1:n_sites_out),&
!               & local_property_models(i)%do_derivatives,&
!               & soap_cart_der, n_neigh_out, local_properties_cart_der(1:3, 1:n_sites_out, i) )
!          do j = 1, n_sites_out
!             i2 = in_to_out_site(j)
!             local_properties0(i2) = local_properties(j)
!          end do
!          if( local_property_models(i)%do_derivatives )then
!             do k = 1, n_atom_pairs
!                k2 = in_to_out_pairs(k)
!                local_properties_cart_der0(1:3, k2) = local_properties_cart_der(1:3, k)
!             end do
!          end if
!       end do

!       deallocate( local_properties )
!       do i = 1, n_local_properties
!          if(local_property_models(i)%do_derivatives) deallocate( local_properties_cart_der )
!       end do

! !call cpu_time(time2)
! !write(*,*) "hirshfeld_v time =", time2-time1, "seconds"
!    end if



!     if( do_prediction )then
! !     Get energies and forces
!       allocate( energies(1:n_sites_out) )
!       energies = 0.d0
!       if( do_forces )then
!         allocate( forces(1:3, 1:n_all_sites) )
!         forces = 0.d0
!       end if
!       if( n_sites_out > 0 )then
!         call get_soap_energy_and_forces(soap, soap_cart_der, alphas, delta, zeta, 0.d0, Qs, &
!                                         n_neigh_out, neighbors_list, xyz, do_forces, do_timing, &
!                                         energies, forces, virial)
!       end if

!       do i = 1, n_sites_out
!         i2 = in_to_out_site(i)
!         energies0(i2) = energies(i)
!       end do
!       if( do_forces )then
!         do i = 1, n_all_sites
!           i2 = in_to_out_site(i)
!           forces0(1:3, i2) = forces(1:3, i)
!         end do
!       end if
!     end if

! !   This is slow, but writing to disk is even slower so who cares
!     if( write_soap )then
!       allocate( soap_temp(1:n_soap, 1:n_sites0) )
!       soap_temp = 0.d0
!       do i = 1, n_sites_out
!         i2 = in_to_out_site(i)
!         soap_temp(1:n_soap, i2) = soap(1:n_soap, i)
!       end do
!       deallocate( soap )
!       allocate( soap(1:n_soap, 1:n_sites0) )
!       soap = soap_temp
!       deallocate( soap_temp )
!     end if
!     if( do_derivatives .and. write_derivatives )then
!       allocate( der_neighbors(1:n_sites_out) )
!       allocate( der_neighbors_list(1:n_atom_pairs) )
!       der_neighbors = n_neigh_out
!       k = 0
!       do i = 1, n_sites_out
!         do j = 1, n_neigh_out(i)
!           k = k + 1
!           j2 = neighbors_list(k)
!           der_neighbors_list(k) = in_to_out_site(j2)
!         end do
!       end do
!     end if



!     if( do_prediction )then
!       deallocate( energies )
!     end if
!     if( do_forces )then
!       deallocate( forces )
!     end if


!     deallocate( neighbors_list, rjs, thetas, phis, mask, mask0, is_atom_seen )
!     deallocate( species0, species_multiplicity0, species, species_multiplicity, species_multiplicity_supercell )
!     if(.not. has_local_properties )then
!        deallocate( in_to_out_site, out_to_in_site , in_to_out_pairs, n_neigh_out )
!     else
!        deallocate(out_to_in_site)
!     end if

!     if( .not. write_soap .and. .not. has_local_properties )then
!        deallocate( soap )
!     end if
!     if( do_derivatives .and. .not. write_derivatives .and. .not. has_local_properties )then
!       deallocate( soap_cart_der )
!    end if

!   end subroutine
! !**************************************************************************


subroutine get_local_properties( soap, Qs, alphas, V0, delta, zeta, &
                              local_property0, do_derivatives, soap_cart_der, n_neigh_out, &
                              local_property_cart_der0, n_atom_pairs,&
                              & in_to_out_pairs, n_all_sites,&
                              & in_to_out_site, n_sites_out )

  real*8, allocatable, intent(in) :: soap(:,:), soap_cart_der(:,:,:)
  real*8, intent(in) ::  Qs(:,:), alphas(:), zeta, delta, V0
  real*8, intent(inout) :: local_property0(:), local_property_cart_der0(:,:)
  real*8, allocatable :: local_property(:), local_property_cart_der(:,:)
  integer, intent(in) :: n_atom_pairs, in_to_out_pairs(:), n_all_sites,&
       & in_to_out_site(:), n_neigh_out(:), n_sites_out
  logical, intent(in) :: do_derivatives

  allocate( local_property( 1:n_sites_out ) )
  local_property = 0.d0
  if( do_derivatives )then
     allocate( local_property_cart_der(1:3, 1:n_atom_pairs) )
     local_property_cart_der = 0.d0
  end if

  call local_property_predict( soap, Qs, alphas, V0, delta, zeta, &
       local_property, do_derivatives, soap_cart_der, n_neigh_out, &
       local_property_cart_der )
  do i = 1, n_sites_out
     i2 = in_to_out_site(i)
     local_property0(i2) = local_property(i)
  end do
  if( do_derivatives )then
     do k = 1, n_atom_pairs
        k2 = in_to_out_pairs(k)
        local_property_cart_der0(1:3, k2) = local_property_cart_der(1:3, k)
     end do
  end if

  deallocate( local_property )
  if( do_derivatives )then
     deallocate( local_property_cart_der )
  end if
!call cpu_time(time2)
!write(*,*) "hirshfeld_v time =", time2-time1, "seconds"
end subroutine get_local_properties






end module
