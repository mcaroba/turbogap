! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_soap.f90, is copyright (c) 2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
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

!  The soap_turbo descriptors and the GAPs on them: one pass over every
!  descriptor, split into batches that fit the memory budget, accumulating
!  energies, forces, virial, local properties and dipoles, and writing the
!  descriptors and their derivatives out when asked.
module turbogap_soap

   use kinds, only: dp
   use types, only: input_parameters
   use timing, only: times_t, time_start, time_end, get_time
   use neighbors, only: get_number_of_atom_pairs
   use gap, only: soap_backend_begin, soap_backend_end
   use gap_interface, only: get_gap_soap
   use turbogap_comm, only: comm_t
   use turbogap_structure, only: state_t
   use turbogap_domain, only: domain_t, neighbors_t
   use turbogap_setup, only: model_t
   use turbogap_results, only: results_t
   use turbogap_loop, only: loop_t
#ifdef _GPU
   use iso_c_binding
   use F_B_C, only: gpu_malloc_async, cpy_htod, gpu_free_async
   use neighbors, only: get_number_of_atom_pairs_batches
   use gpu_context, only: gpu_stream
#endif

   implicit none

   private
   public :: compute_soap
   public :: soap_free_device

contains

   subroutine compute_soap(res, state, nl, dom, model, params, loop, comm, time)
      type(results_t), target, intent(inout) :: res
      type(state_t), intent(in) :: state
      type(neighbors_t), intent(inout) :: nl
      type(domain_t), intent(in) :: dom
      type(model_t), target, intent(inout) :: model
      type(input_parameters), intent(inout) :: params
      type(loop_t), intent(in) :: loop
      type(comm_t), intent(in) :: comm
      type(times_t), intent(inout) :: time
      real(dp), allocatable :: soap(:, :)
      real(dp), allocatable :: soap_cart_der(:, :, :)
      integer, allocatable :: der_neighbors(:)
      integer, allocatable :: der_neighbors_list(:)
      integer, allocatable :: i_beg_list(:)
      integer, allocatable :: i_end_list(:)
      integer, allocatable :: j_beg_list(:)
      integer, allocatable :: j_end_list(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: i2
      integer :: k2
      integer :: n_lp_count
      integer :: n_sites_this
      integer :: n_soap
      integer :: this_i_beg
      integer :: this_i_end
      integer :: this_j_beg
      integer :: this_j_end
      integer :: this_n_sites_mpi
      character*8 :: i_char
#ifdef _GPU
      integer(c_size_t) :: st_size_nf
      type(c_ptr) :: nf_d
      type(c_ptr) :: rcut_hard_d
      type(c_ptr) :: rcut_soft_d
      type(c_ptr) :: global_scaling_d
      type(c_ptr) :: atom_sigma_r_d
      type(c_ptr) :: atom_sigma_r_scaling_d
      type(c_ptr) :: atom_sigma_t_d
      type(c_ptr) :: atom_sigma_t_scaling_d
      type(c_ptr) :: amplitude_scaling_d
      type(c_ptr) :: alpha_max_d
      type(c_ptr) :: central_weight_d
      integer :: n_sparse
      integer :: dim
      integer :: n_sp
#endif

      !     Loop through soap_turbo descriptors - we always call this routine, even if we don't want to do prediction
      n_lp_count = 0 ! This counts the local properties
      call time_start(time%gap)
      do i = 1, model%n_soap_turbo
         call time_start(time%soap)
         !       Compute number of pairs for this SOAP. SOAP has in general a different cutoff than overall max
         !       cutoff, so the number of pairs may be a lot smaller for the SOAP subset.
         !       This subroutine splits the load optimally so as to not use more memory per MPI process than available.
         !       TurboGAP does not check how much memory is available, it just relies on heuristics and a user provided
         !       max_Gbytes_per_process (default = 1.d0)
#ifdef _GPU
         if (params%n_batches > 0) then
            call get_number_of_atom_pairs_batches(params%n_batches, nl%n_neigh(dom%i_beg:dom%i_end), &
                                                  nl%rjs(dom%j_beg:dom%j_end), model%soap_turbo_hypers(i)%rcut_max, &
                                                  model%soap_turbo_hypers(i)%l_max, &
                                                  model%soap_turbo_hypers(i)%n_max, &
                                                  model%soap_turbo_hypers(i)%dim, &
                                                  model%soap_turbo_hypers(i)%n_species, &
                                                  params%max_Gbytes_per_process, i_beg_list, &
                                                  i_end_list, j_beg_list, j_end_list)
         else
            call get_number_of_atom_pairs(nl%n_neigh(dom%i_beg:dom%i_end), nl%rjs(dom%j_beg:dom%j_end), &
                                          model%soap_turbo_hypers(i)%rcut_max, &
                                          model%soap_turbo_hypers(i)%l_max, &
                                          model%soap_turbo_hypers(i)%n_max, &
                                          model%soap_turbo_hypers(i)%dim, &
                                          model%soap_turbo_hypers(i)%n_species, &
                                          params%max_Gbytes_per_process, i_beg_list, &
                                          i_end_list, j_beg_list, j_end_list)
         end if

         n_sp = model%soap_turbo_hypers(i)%n_species

         st_size_nf = n_sp*sizeof(model%soap_turbo_hypers(i)%nf(1))
         call gpu_malloc_async(nf_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%nf), nf_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(rcut_hard_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%rcut_hard), rcut_hard_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(rcut_soft_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%rcut_soft), rcut_soft_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(global_scaling_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%global_scaling), global_scaling_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(atom_sigma_r_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_r), atom_sigma_r_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_r_scaling), atom_sigma_r_scaling_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(atom_sigma_t_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_t), atom_sigma_t_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%atom_sigma_t_scaling), atom_sigma_t_scaling_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(amplitude_scaling_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%amplitude_scaling), amplitude_scaling_d, st_size_nf, gpu_stream)
         call gpu_malloc_async(central_weight_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%central_weight), central_weight_d, st_size_nf, gpu_stream)
         st_size_nf = n_sp*sizeof(model%soap_turbo_hypers(i)%alpha_max(1))
         call gpu_malloc_async(alpha_max_d, st_size_nf, gpu_stream)
         call cpy_htod(c_loc(model%soap_turbo_hypers(i)%alpha_max), alpha_max_d, st_size_nf, gpu_stream)
         n_sparse = model%soap_turbo_hypers(i)%n_sparse
         dim = model%soap_turbo_hypers(i)%dim
         call soap_backend_begin(model%soap_turbo_hypers(i))

         if (model%soap_turbo_hypers(i)%has_local_properties) then
            ! Allocate gpu memory
            do j = 1, model%soap_turbo_hypers(i)%n_local_properties
               model%soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas = &
                  model%soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                  sizeof(model%soap_turbo_hypers(i)%local_property_models(j)%alphas(1))
               call gpu_malloc_async(model%soap_turbo_hypers(i)%local_property_models(j)%alphas_d, &
                                     model%soap_turbo_hypers(i)%local_property_models(j)%st_size_alphas, gpu_stream)
               call cpy_htod(c_loc(model%soap_turbo_hypers(i)&
                 &%local_property_models(j)%alphas), &
                 & model%soap_turbo_hypers(i)%local_property_models(j)&
                 &%alphas_d, model%soap_turbo_hypers(i)&
                 &%local_property_models(j)%st_size_alphas,&
                 & gpu_stream)

               model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs = &
                  model%soap_turbo_hypers(i)%local_property_models(j)%n_sparse* &
                  model%soap_turbo_hypers(i)%local_property_models(j)%dim* &
                  sizeof(model%soap_turbo_hypers(i)%local_property_models(j)%Qs(1, 1))

               call gpu_malloc_async(model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                                     model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, gpu_stream)
               call cpy_htod(c_loc(model%soap_turbo_hypers(i)%local_property_models(j)%Qs), &
                             model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, &
                             model%soap_turbo_hypers(i)%local_property_models(j)%st_size_Qs, &
                             gpu_stream)

            end do

         end if
#else
         call get_number_of_atom_pairs(nl%n_neigh(dom%i_beg:dom%i_end), nl%rjs(dom%j_beg:dom%j_end), &
                                       model%soap_turbo_hypers(i)%rcut_max, &
                                       model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%n_max, &
                                       model%soap_turbo_hypers(i)%dim, model%soap_turbo_hypers(i)%n_species, &
                                       params%max_Gbytes_per_process, i_beg_list, i_end_list, j_beg_list, j_end_list)
#endif

         do j = 1, size(i_beg_list)
            this_i_beg = dom%i_beg - 1 + i_beg_list(j)
            this_i_end = dom%i_beg - 1 + i_end_list(j)
            this_j_beg = dom%j_beg - 1 + j_beg_list(j)
            this_j_end = dom%j_beg - 1 + j_end_list(j)
            this_n_sites_mpi = this_i_end - this_i_beg + 1
            res%this_energies = 0.d0
            if (params%do_forces) then
               res%this_forces = 0.d0
               res%this_virial = 0.d0
            end if
            if (model%soap_turbo_hypers(i)%is_dipole_model) then
               res%this_local_dipoles = 0.d0
               res%this_energies_dipole = 0.d0
            end if
            if (model%soap_turbo_hypers(i)%has_local_properties) then
               res%this_local_properties = 0.d0
               if (params%do_forces) then
                  res%this_local_properties_cart_der = 0.d0
                  !             I don't remember why this needs a pointer <----------------------------------------- CHECK
                  nullify (res%this_local_properties_cart_der_pt)
                  res%this_local_properties_cart_der_pt =>&
                       & res%this_local_properties_cart_der(1:3,&
                       & this_j_beg:this_j_end, 1:params&
                       &%n_local_properties)
               end if
            end if

#ifdef _GPU
            call get_gap_soap( &
               n_sparse, state%n_sites, this_n_sites_mpi, nl%n_neigh(this_i_beg:this_i_end), &
               nl%neighbors_list(this_j_beg:this_j_end), model%soap_turbo_hypers(i)%n_species, &
               model%soap_turbo_hypers(i)%species_types, nl%rjs(this_j_beg:this_j_end), nl%thetas(this_j_beg:this_j_end), &
               nl%phis(this_j_beg:this_j_end), nl%xyz(1:3, this_j_beg:this_j_end), alpha_max_d, &
               model%soap_turbo_hypers(i)%alpha_max, model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%dim, &
               rcut_hard_d, &
               model%soap_turbo_hypers(i)%rcut_hard, rcut_soft_d, nf_d, global_scaling_d, atom_sigma_r_d, &
               model%soap_turbo_hypers(i)%atom_sigma_r, atom_sigma_r_scaling_d, atom_sigma_t_d, atom_sigma_t_scaling_d, &
               amplitude_scaling_d, model%soap_turbo_hypers(i)%radial_enhancement, central_weight_d, &
               model%soap_turbo_hypers(i)%central_weight, model%soap_turbo_hypers(i)%basis, &
               model%soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
               params%do_prediction, params%write_soap, params%write_derivatives, &
               model%soap_turbo_hypers(i)%compress_soap, model%soap_turbo_hypers(i)%compress_soap_indices, &
               model%soap_turbo_hypers(i)%delta, model%soap_turbo_hypers(i)%zeta, model%soap_turbo_hypers(i)%central_species, &
               state%xyz_species(this_i_beg:this_i_end), state%xyz_species_supercell, params%all_atoms, &
               params%which_atom, state%indices, soap, soap_cart_der, der_neighbors, der_neighbors_list, &
               model%soap_turbo_hypers(i)%has_local_properties, model%soap_turbo_hypers(i)%n_local_properties, &
               model%soap_turbo_hypers(i)%local_property_models, n_lp_count, res%energies_soap, res%forces_soap, &
               res%this_local_properties_pt, res%this_local_properties_cart_der_pt, model%local_property_indexes, &
               res%this_virial, &
               time%soap_lin(3), time%get_soap(3), model%soap_turbo_hypers(i)%W_d, model%soap_turbo_hypers(i)%S_d, &
               model%soap_turbo_hypers(i)%multiplicity_array_d, model%soap_turbo_hypers(i)%st_W_d, &
               model%soap_turbo_hypers(i)%st_S_d, model%soap_turbo_hypers(i)%st_multiplicity_array_d, &
               model%soap_turbo_hypers(i)%recompute_basis, time%local_prop, &
               model%soap_turbo_hypers(i)%is_dipole_model, res%local_dipoles, res%energies_dipole)
#else
            call soap_backend_begin(model%soap_turbo_hypers(i))
            call get_gap_soap(state%n_sites, this_n_sites_mpi, nl%n_neigh(this_i_beg:this_i_end), &
               nl%neighbors_list(this_j_beg:this_j_end), &
                 model%soap_turbo_hypers(i)%n_species, model%soap_turbo_hypers(i)%species_types, &
                 nl%rjs(this_j_beg:this_j_end), nl%thetas(this_j_beg:this_j_end), nl%phis(this_j_beg:this_j_end), &
                 nl%xyz(1:3, this_j_beg:this_j_end), &
                 model%soap_turbo_hypers(i)%alpha_max, &
                 model%soap_turbo_hypers(i)%l_max, model%soap_turbo_hypers(i)%dim, model%soap_turbo_hypers(i)%rcut_hard, &
                 model%soap_turbo_hypers(i)%rcut_soft, model%soap_turbo_hypers(i)%nf, &
                    model%soap_turbo_hypers(i)%global_scaling, &
                 model%soap_turbo_hypers(i)%atom_sigma_r, model%soap_turbo_hypers(i)%atom_sigma_r_scaling, &
                 model%soap_turbo_hypers(i)%atom_sigma_t, model%soap_turbo_hypers(i)%atom_sigma_t_scaling, &
                 model%soap_turbo_hypers(i)%amplitude_scaling, model%soap_turbo_hypers(i)%radial_enhancement, &
                 model%soap_turbo_hypers(i)%central_weight, model%soap_turbo_hypers(i)%basis, &
                 model%soap_turbo_hypers(i)%scaling_mode, params%do_timing, params%do_derivatives, params%do_forces, &
                 params%do_prediction, params%write_soap, params%write_derivatives, &
                 model%soap_turbo_hypers(i)%compress_soap, model%soap_turbo_hypers(i)%compress_P_nonzero, &
                 model%soap_turbo_hypers(i)%compress_P_i, model%soap_turbo_hypers(i)%compress_P_j, &
                 model%soap_turbo_hypers(i)%compress_P_el, &
                 model%soap_turbo_hypers(i)%delta, model%soap_turbo_hypers(i)%zeta, model%soap_turbo_hypers(i)%central_species, &
                 state%xyz_species(this_i_beg:this_i_end), state%xyz_species_supercell, &
                 params%all_atoms, params%which_atom, state%indices, soap, soap_cart_der, &
                 der_neighbors, der_neighbors_list, &
                 & model%soap_turbo_hypers(i)%has_local_properties,&
                 & model%soap_turbo_hypers(i)%n_local_properties,&
                 & model%soap_turbo_hypers(i)%local_property_models,&
                 & res%this_energies, res%this_forces, res%this_local_properties_pt,&
                 & res%this_local_properties_cart_der_pt,&
                 & model%local_property_indexes, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                 & res%this_virial, n_lp_count, model%soap_turbo_hypers(i)%is_dipole_model, &
                 & res%this_local_dipoles, res%this_energies_dipole)

            call soap_backend_end()
#endif

            ! We can have a pointer to specific parts of this_local_properties array to then

!              A dipole descriptor leaves this_energies and this_forces at the
!              zero they were set to above -- get_gap_soap never writes them --
!              so its fictitious energy stays out of energies_soap and its
!              gradient out of forces_soap. It is carried separately.
            res%energies_soap = res%energies_soap + res%this_energies

            if (model%soap_turbo_hypers(i)%is_dipole_model) then
               res%local_dipoles = res%local_dipoles + res%this_local_dipoles
               res%energies_dipole = res%energies_dipole + res%this_energies_dipole
            end if

            if (model%soap_turbo_hypers(i)%has_local_properties) then

               res%local_properties(:, :) = res%local_properties(:, :) + res%this_local_properties(:, :)
               if (any(model%soap_turbo_hypers(i)&
                    &%local_property_models(:)%do_derivatives) &
                    & .and. params%do_derivatives) then
                  res%local_properties_cart_der(:, :, :) =&
                       & res%local_properties_cart_der(:, :, :) +&
                       & res%this_local_properties_cart_der(:, :, :)
               end if

            end if
            if (params%do_forces) then
               res%forces_soap = res%forces_soap + res%this_forces
               res%virial_soap = res%virial_soap + res%this_virial
            end if
         end do
         n_lp_count = n_lp_count + model%soap_turbo_hypers(i)%n_local_properties

#ifdef _GPU
         call gpu_free_async(nf_d, gpu_stream)
         call gpu_free_async(rcut_hard_d, gpu_stream)
         call gpu_free_async(rcut_soft_d, gpu_stream)
         call gpu_free_async(global_scaling_d, gpu_stream)
         call gpu_free_async(atom_sigma_r_d, gpu_stream)
         call gpu_free_async(atom_sigma_r_scaling_d, gpu_stream)
         call gpu_free_async(atom_sigma_t_d, gpu_stream)
         call gpu_free_async(atom_sigma_t_scaling_d, gpu_stream)
         call gpu_free_async(amplitude_scaling_d, gpu_stream)
         call gpu_free_async(alpha_max_d, gpu_stream)
         call gpu_free_async(central_weight_d, gpu_stream)

         if (model%soap_turbo_hypers(i)%has_local_properties) then
            do j = 1, model%soap_turbo_hypers(i)%n_local_properties
               call gpu_free_async(model%soap_turbo_hypers(i)%local_property_models(j)%alphas_d, gpu_stream)
               call gpu_free_async(model%soap_turbo_hypers(i)%local_property_models(j)%Qs_d, gpu_stream)
            end do
         end if

         call soap_backend_end()

         call get_time(time%soap_solo(2))
#endif
         deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)
#ifdef _GPU
         time%soap_solo(3) = time%soap_solo(3) + time%soap_solo(2) - time%soap_solo(1)
#endif

         ! THIS WON'T WORK! THE SOAP AND SOAP DERIVATIVES NEED TO BE COLLECTED FROM ALL RANKS <--------------------- FIX THIS!!!!
         ! AT THE MOMENT I'M MAKING THE CODE PRINT AN ERROR MESSAGE AND STOP EXECUTION IF THE USER TRIES TO WRITE OUT THESE
         ! FILES WITH MORE THAN ONE MPI TASK
         if (comm%rank == 0) then
            !       Write out stuff - THIS SHOULD PROBABLY BE PUT IN A MODULE
            if (model%n_soap_turbo == 1) then
               i_char = ""
            else
               write (i_char, '(I7)') i
               i_char = "_"//adjustl(i_char)
            end if
            !       Write the SOAP vectors - NOT THE OPTIMAL STRATEGY IN TERMS OF DISK SPACE SINCE SOME ATOMS HAVE SOAP = 0
            if (params%write_soap) then
               if (loop%n_xyz == 1 .or. loop%md_istep == 0) then
                  open (unit=10, file="soap"//trim(i_char)//".dat", status="unknown")
               else
                  open (unit=10, file="soap"//trim(i_char)//".dat", status="old", position="append")
               end if
               if (.not. params%do_md .or. &
                   (params%do_md .and. (loop%md_istep == 0 .or. loop%md_istep == params%md_nsteps .or. &
                                        modulo(loop%md_istep, params%write_xyz) == 0))) then
                  n_sites_this = size(soap, 2)
                  n_soap = size(soap, 1)
                  write (10, *) n_sites_this, n_soap
                  do i2 = 1, n_sites_this
                     write (10, '(*(ES24.15))') soap(1:n_soap, i2)
                  end do
               end if
               close (10)
            end if
            if (allocated(soap)) deallocate (soap)

            !       Optionally, write out the derivatives (might take a lot of disk space)
            if ((params%do_derivatives .or. params%do_derivatives_fd) .and. params%write_derivatives) then
               if (loop%n_xyz == 1 .or. loop%md_istep == 0) then
                  open (unit=10, file="soap_der"//trim(i_char)//".dat", status="unknown")
               else
                  open (unit=10, file="soap_der"//trim(i_char)//".dat", status="old", position="append")
               end if
               if (.not. params%do_md .or. &
                   (params%do_md .and. (loop%md_istep == 0 .or. loop%md_istep == params%md_nsteps .or. &
                                        modulo(loop%md_istep, params%write_xyz) == 0))) then
                  !           Note, this n_sites is not the same as the total number of sites, it's just the total number
                  !           of sites that have a derivative, since the first neighbor of each site is itself, the site
                  !           ID can always be retrieved from there. Note also that the sites are not necessarily given in
                  !           order
                  n_sites_this = size(der_neighbors, 1)
                  n_soap = size(soap_cart_der, 2)
                  nl%n_atom_pairs = size(der_neighbors_list, 1)
                  write (10, *) state%n_sites, n_soap, nl%n_atom_pairs
                  k = 1
                  k2 = 0
                  do i2 = 1, n_sites_this
                     write (10, *) der_neighbors_list(k), der_neighbors(i2), der_neighbors_list(k:k + der_neighbors(i2) - 1)
                     k = k + der_neighbors(i)
                     do j = 1, der_neighbors(i)
                        k2 = k2 + 1
                        write (10, '(*(ES24.15))') soap_cart_der(1, 1:n_soap, k2)
                        write (10, '(*(ES24.15))') soap_cart_der(2, 1:n_soap, k2)
                        write (10, '(*(ES24.15))') soap_cart_der(3, 1:n_soap, k2)
                     end do
                  end do
               end if
               close (10)
            end if
            if (params%write_derivatives) then
               deallocate (soap_cart_der, der_neighbors, der_neighbors_list)
            end if
         end if

         call time_end(time%soap)

      end do
      call time_end(time%gap)
   end subroutine compute_soap

!  Release the device copies of the compression basis that outlive a step.
   subroutine soap_free_device(model)
      type(model_t), intent(inout) :: model
#ifdef _GPU
      integer :: i

      do i = 1, model%n_soap_turbo
         if (.not. model%soap_turbo_hypers(i)%recompute_basis) then
            call gpu_free_async(model%soap_turbo_hypers(i)%W_d, gpu_stream)
            call gpu_free_async(model%soap_turbo_hypers(i)%S_d, gpu_stream)
            call gpu_free_async(model%soap_turbo_hypers(i)%multiplicity_array_d, gpu_stream)
         end if
      end do
#endif
   end subroutine soap_free_device

end module turbogap_soap
