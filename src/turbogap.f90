! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap.f90, is copyright (c) 2019-2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X   Uttiyoarnab saha
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

program turbogap

   use kinds

   use timing
   use neighbors
   use soap_turbo_desc
   use gap
   use read_files
   use md
   use adaptive_time                        ! for adaptive time simulation (TurboGAP will use these five modules for radiation cascades)
   use electronic_stopping                ! for electronic stopping correction in radiation cascades
   use eph_fdm                                ! for T - dependent parameters - elec. stop. - eph model
   use eph_beta                                ! for the atomic electronic densities  - elec. stop. - eph model
   use eph_electronic_stopping                ! for electronic stopping based in radiation cascades on the eph model
   use mc
   use gap_interface
   use types
   use vdw
   use electrostatics, only: compute_coulomb_direct, compute_coulomb_dsf, compute_coulomb_lamichhane
   use turbogap_setup
   use turbogap_structure, only: state_t, structure_acquire
   use turbogap_domain, only: domain_t, neighbors_t, domain_sync_state, domain_build, &
                              domain_complete_sites, domain_complete_e0, &
                              domain_complete_contributions, domain_sync_after_md, &
                              domain_sync_after_ipi
   use turbogap_results, only: results_t, results_prepare
   use turbogap_soap, only: compute_soap
   use turbogap_sampling, only: sampling_t, mc_prepare_step, nested_step
   use turbogap_ir, only: ir_run_t, ir_init, ir_step_begin, ir_before_evaluate, ir_push_frame, &
                          ir_after_forces, ir_step_end, ir_finish, ir_report
   use turbogap_loop, only: loop_t, loop_init, loop_continues, loop_begin_step, creturn
   use turbogap_output, only: print_banner, print_options
   use turbogap_exp
   use turbogap_md
   use ipi_driver, only: ipi_driver_open, ipi_driver_exchange, ipi_driver_close
   use gap_backend
   use gpu_context
   use turbogap_vdw
   use turbogap_estat
   use exp_utils
   use exp_interface
   use soap_turbo_functions
   use mad_ir
   use mad_ir_xl
   use ir_auxiliary_dynamics, only: ir_aux_state, ir_aux_active, ir_aux_setup, &
                                    ir_aux_advance, ir_aux_evaluate, ir_aux_forces, &
                                    ir_aux_calibrate, ir_aux_calibrated, ir_aux_save, &
                                    ir_aux_write_spectrum, ir_aux_bank_energy, &
                                    ir_aux_energy_pumped, ir_aux_stability, &
                                    ir_aux_escale_max
   use ir_fft
   use ir_fft_io
   use turbogap_comm, only: comm_t, comm_init, comm_finalize, comm_bcast, comm_sum_to_root, &
                            comm_sum_all, comm_allgather, comm_with_mpi
   use bussi
   use xyz_module
   use keyword_help
#ifdef _GPU
   use F_B_C
   use iso_c_binding
#endif

   implicit none

   ! Variable definitions

   real(dp) :: time1
   real(dp) :: time2
   real(dp) :: time3
!   Every wall-clock bucket lives in one times_t (src/timing.f90), so the
!   extracted modules take a single argument instead of thirteen and the two
!   branches' signatures agree.
   type(times_t) :: time
   logical :: write_condition = .false.
   logical :: overwrite_condition = .false.

   ! Clean up these variables after code refactoring !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer, allocatable :: alpha_max(:)
   integer, allocatable :: i_beg_list(:)
   integer, allocatable :: i_end_list(:)
   integer, allocatable :: j_beg_list(:)
   integer, allocatable :: j_end_list(:)
   integer :: i
   integer :: j
   integer :: ierr
   integer :: rank
   integer :: ntasks
   type(comm_t) :: comm
   type(state_t) :: state
   type(domain_t) :: dom
   type(neighbors_t) :: nl
   type(model_t), target :: model
   type(results_t), target :: res
   logical :: resized
   type(loop_t) :: loop
   type(dynamics_t) :: dyn
   type(ir_run_t) :: ir
   type(sampling_t) :: smp
   integer :: n_pos
   integer :: this_i_beg
   integer :: this_i_end
   integer :: this_j_beg
   integer :: this_j_end

   integer :: l_max
   integer :: n_max
   integer :: central_species = 0
   integer :: iostatus

   type(perform_t) :: perform
   integer :: which_atom = 0
   integer :: n_omp = 1
   integer :: radial_enhancement = 0

   character*1024 :: filename
   character*1024 :: string
   character*1024 :: temp_string
   character*1024 :: temp_string2

   ! This is the mode in which we run TurboGAP
   character*16 :: mode = "none"
   character*16 :: help_topic = ""
   character*32 :: exp_output = "none"

   ! Here we store the input parameters
   type(input_parameters) :: params

   logical :: do_electrostatics = .true.
! Persistent ts+mbd correction state, owned by turbogap_vdw
   type(vdw_state) :: vdw_ws

   character*32 :: implemented_exp_observables(1:5)
#ifdef _GPU
   integer :: omp_task
   type(c_ptr) :: alphas_d
   type(c_ptr) :: qs_d
   integer :: n_pairs_temp
   integer :: n_sites_temp
   character*8, allocatable, target :: species_types_actual(:)
#endif

!  --help answers from the generated keyword reference and exits. It is the
!  first thing the program does because it must work with no input file, no
!  GPU and no MPI: everything below this point assumes at least one of those.
   call get_command_argument(1, mode)
   if (mode == "--help" .or. mode == "-h" .or. mode == "help") then
      call get_command_argument(2, help_topic)
!     Validated against the SAME list the error message prints, which
!     keyword_help.f90 generates from tools/keyword_docs.py. It used to be a
!     hardcoded chain of comparisons beside a generated message, and the two
!     drifted the moment a mode was added: --help ipi was rejected by a message
!     that listed ipi as valid. Slashes on both sides so that a topic cannot
!     match a substring of another.
      if (len_trim(help_topic) > 0 .and. &
          index("/"//trim(keyword_help_topics())//"/", "/"//trim(help_topic)//"/") == 0) then
         write (*, '(A)') 'ERROR: unknown help topic "'//trim(help_topic)// &
            '". turbogap --help ['//trim(keyword_help_topics())//']'
         stop 1
      end if
      call print_keyword_help(help_topic)
      stop
   end if
   mode = "none"

   implemented_exp_observables(1) = "xps"
   implemented_exp_observables(2) = "xrd"
   implemented_exp_observables(3) = "saxs"
   implemented_exp_observables(4) = "pair_distribution"
   implemented_exp_observables(5) = "structure_factor"

!  Bring the device context up. Empty on this branch (src/gpu_context.f90);
!  on the GPU branch the same two names create the streams and cuBLAS handles.
   call time_start(time%create_streams)
   call gpu_context_init(params, rank, n_omp)
   call gap_backend_init()
   call time_end(time%create_streams)

   ! Start recording the time
   call get_time(time1)
   time3 = time1
!  Everything before the first pass of the main loop: the input file, the
!  potential files, and the allocation and broadcast that follow them. Without
!  this bucket the pre-loop cost fell into Miscellaneous, which is why a run
!  whose real work took 1.2 s reported 0.4 s "miscellaneous" and a 31 s run
!  reported the same 0.4 s -- a constant, and therefore obviously a setup cost,
!  but not one the report could name.
   call time_start(time%setup)
   ! Start random seed
   call srand(int(time1*1000))

   call comm_init(comm)
   rank = comm%rank
   ntasks = comm%size
   allocate (dom%n_atom_pairs_by_rank(1:ntasks))

   ! Read the mode. It should be "soap", "predict" or "md"
   call get_command_argument(1, mode)
   if (mode == "" .or. mode == "none") then
      write (*, *) "ERROR: you need to run 'turbogap md', 'turbogap mc', 'turbogap predict'"
      write (*, *) "       or 'turbogap ipi' (forces for an i-PI server; see ipi_address)"
      write (*, *) "       'turbogap --help [predict|md|mc|soap|gap]' lists the keywords"
      stop
      ! THIS SHOULD BE FIXED, IN CASE THE USER JUST WANT TO OUTPUT THE SOAP DESCRIPTORS
      mode = "soap"
   end if

   call print_banner(comm)

   ! Read input file and other files
   call read_input_and_gap_files(mode, rank, ntasks, params, &
                                 model%soap_turbo_hypers, model%distance_2b_hypers, model%angle_3b_hypers, model%core_pot_hypers, &
                                 model%n_soap_turbo, model%n_distance_2b, model%n_angle_3b, model%n_core_pot, model%n_species, &
                                 model%rcut_max, &
                                 model%valid_xps, model%xps_idx, model%vdw_lp_index, model%core_be_lp_index, &
                                 model%valid_estat_charges, model%charge_lp_index, &
                                 model%local_property_labels, model%local_property_indexes, model%n_local_properties_mpi, &
                                 model%has_local_properties_mpi, model%local_properties_n_sparse_mpi_soap_turbo, &
                                 model%local_properties_dim_mpi_soap_turbo, dyn%nrows, dyn%allelstopdata, &
                                 dyn%ephbeta, dyn%ephfdm, dyn%ephlsc, time)

!  The host memory budget, which has to sit exactly here.
!
!  After read_input_and_gap_files, because it reads mem_fraction and writes
!  max_Gbytes_per_process -- placed next to gpu_context_init it would run before
!  the input existed and size the loop from defaults whatever the input said.
!
!  ntasks is passed for the case where MPI cannot say how the ranks are laid
!  out; the routine prefers to ask MPI which ranks share a node, because that is
!  the set that shares the memory it is dividing.
!
!  The GPU branch calls the same name in the same place, where it budgets from
!  the device instead of from the node.
   call gpu_memory_budget_init(params, rank, ntasks)

   call print_options(comm, params, model)

   ! Print progress bar and initialize timers

   model%xps_idx = params%xps_idx
   smp%i_nested = 0
   smp%i_image = 0

   call loop_init(loop, params, comm)

   ! This checks if we need to do the SOAP calculation more than once, if there are several concatenated
   ! structures in the xyz file provided or we're doing molecular dynamics

!   The exp-observable decisions, evaluated once.  Every input is a params
!   field or valid_xps, none of which changes inside the main loop.
!
!   This closes a defect.  The allocation guards asked do_X .and. valid_X, the
!   zeroing guards asked do_X .and. exp_forces .and. valid_X, and the force
!   accumulation asked only exp_forces .and. valid_X -- so a deck supplying an
!   experimental dataset for an observable it had not switched on, with
!   exp_forces set, accumulated forces_X and virial_X that the allocation
!   guard had skipped.  do_X and valid_X are independent: valid_X is set from
!   a label in the experimental data file, do_X is its own input keyword.
!   Same shape as the electrostatics guard and as has_vdw against
!   has_local_properties.
   perform%pdf = params%do_pair_distribution .and. params%valid_pdf
   perform%sf = params%do_structure_factor .and. params%valid_sf
   perform%xrd = params%do_xrd .and. params%valid_xrd
   perform%nd = params%do_nd .and. params%valid_nd

   perform%pdf_forces = perform%pdf .and. params%exp_forces
   perform%sf_forces = perform%sf .and. params%exp_forces
   perform%xrd_forces = perform%xrd .and. params%exp_forces
   perform%nd_forces = perform%nd .and. params%exp_forces
   perform%xps_forces = model%valid_xps .and. params%exp_forces

   call ir_init(ir, params, comm)

   call time_end(time%setup)

!  Connect before the first force call, so that a missing or unstarted i-PI
!  server is reported now rather than after the first GAP evaluation.
   if (mode == "ipi") call ipi_driver_open(params%ipi_address, rank)

   do while (loop_continues(loop, params))
      loop%exit_loop = .false.

      call ir_step_begin(ir, params)

      call loop_begin_step(loop, params, comm)

      call structure_acquire(state, nl%rebuild_neighbors_list, loop, params, model, comm, smp%mc_file, time)
      !   Broadcast the info in the XYZ file: positions, velocities, masses, xyz_species, xyz_species_supercell,
      !   species, species_supercell, indices, a_box, b_box, c_box and n_sites. I should put this into a module!!!!!!!

      call md_prepare_velocities(dyn, state, params, loop, comm)
      call mc_prepare_step(smp, dyn, state, params, loop, comm)
      call domain_sync_state(dom, comm, state, params, time)
      call domain_build(dom, nl, comm, state, params, model, loop, smp%mc_file, time)
      !   Compute the volume of the "primitive" unit cell
      state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                               state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))

      !   If we are doing prediction, we run this chunk of code
      if (params%do_prediction .or. params%write_soap .or. params%write_derivatives) then

         call results_prepare(res, state, params, model, perform, loop, &
                              dom%n_atom_pairs_by_rank(rank + 1), dom%n_atom_pairs_by_rank_prev, resized)
         if (resized) call vdw_read_ts_scaling(res, state, params, loop, comm)
         if (params%do_forces) call ir_before_evaluate(ir, params, state, res, loop, comm)

         if (params%do_prediction) then
            !       Assign the e0 to each atom according to its species
            do i = dom%i_beg, dom%i_end
               do j = 1, model%n_species
                  if (state%xyz_species(i) == params%species_types(j)) then
                     res%energies(i) = params%e0(j)
                  end if
               end do
            end do
         end if
         !     Collect all energies
         call domain_complete_e0(dom, comm, res, state, time)

         call compute_soap(res, state, nl, dom, model, params, loop, comm, time)

         call domain_complete_sites(dom, comm, res, state, params, model, time)

         if (params%do_dipole) then
            res%dipole(1) = sum(res%local_dipoles(1, 1:state%n_sites))
            res%dipole(2) = sum(res%local_dipoles(2, 1:state%n_sites))
            res%dipole(3) = sum(res%local_dipoles(3, 1:state%n_sites))
         end if

!        IR PREDICTION FROM A TRAJECTORY. This frame's total dipole joins the
!        ensemble, with the time its comment line claimed. Nothing is
!        transformed yet: the file's length is not known until the end of it,
!        and the resolution follows from that length.
!
!        Every rank keeps the same buffer. local_dipoles was all-reduced and
!        broadcast just above, so the sums agree bit for bit, and having the
!        ensemble replicated means the final transform needs no communication.
!        Three doubles a frame; a 100 ps trajectory at 1 fs is 2.4 MB.
         call ir_push_frame(ir, res, state)

         !     Compute vdW energies and forces

!        Compute ELECTROSTATIC energies and forces
!
!        Ported from the GPU branch. That branch additionally routes the gsf
!        method through a batched device implementation when params%gpu_batched
!        is set; here gsf always takes the compute_coulomb_lamichhane path,
!        which is what the GPU branch itself falls back to.
!        valid_estat_charges is part of the guard, not an afterthought: without it a
!        deck that asks for electrostatics against a GAP with no atomic_charge local
!        property indexes local_properties with an uninitialised charge_lp_index and
!        segfaults. Same shape as the has_vdw/has_local_properties defect.
!        Writes this rank's partial sums into the this_ arrays, which the
!        reduction below completes.
#ifdef _GPU
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, state%species, nl%neighbor_species, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, n_omp, &
                            res%this_energies_estat, res%this_forces_estat, res%this_virial_estat, time)
#else
         call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                            state%n_sites, nl%n_neigh, nl%neighbors_list, nl%rjs, nl%xyz, &
                            res%local_properties, res%local_properties_cart_der, &
                            dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, &
                            res%this_energies_estat, res%this_forces_estat, res%this_virial_estat, time)
#endif

         call compute_vdw(params, any_has_vdw(model%soap_turbo_hypers), state%n_sites, &
                          nl%n_neigh, nl%neighbors_list, nl%neighbor_species, nl%rjs, nl%xyz, &
                          res%local_properties, res%local_properties_cart_der, model%vdw_lp_index, &
                          dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, dom%n_atom_pairs_by_rank, dom%site_in_rank, &
                          state%indices, rank, ntasks, loop%md_istep, vdw_ws, &
                          res%energies_vdw, res%forces_vdw, res%virial_vdw, res%local_virial_vdw_diag, &
                          res%this_energies_vdw, res%this_forces_vdw, res%this_virial_vdw, &
                          res%this_local_virial_vdw_diag, res%energies_vdw_corr, res%forces_vdw_corr, &
                          res%local_virial_vdw_diag_corr, res%mbd_ts_scaling, res%this_mbd_ts_scaling, &
                          res%update_mbd_ts_scaling, time)

         !--- EXPERIMENTAL SPECTRUM CALCULATION AND FORCES ---!

         ! --- Changing the implementation:
         !     > All experimental prediction should be done here
         !     > do_exp is the variable which says whether calculation should be done
         !     > experimental_forces = .true. will add forces to the calculation

         !###---   Compute Experimental Data Interpolation   ---###!

         if (params%do_exp) then
            do i = 1, params%n_exp
               ! If we want to compute the experimental interpolation, we do it now.

               call get_write_condition(params%do_mc, params%do_md&
                    &, loop%mc_istep, loop%md_istep, params%write_xyz,&
                    & write_condition)

               if (params%exp_data(i)%compute_exp) then
                  if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
                  if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
                  call calculate_exp_interpolation(params%exp_data(i)&
                       &%x, params%exp_data(i)%y, params%exp_data(i)&
                       &%n_samples, params%exp_data(i)%data)

!                 The weights live on the same grid as the experiment, so they
!                 are built here, from the same x, every time it is rebuilt.
                  call build_exp_weights(params%exp_data(i)%x, params%exp_data(i)%w, &
                                         params%exp_data(i)%weights_data, &
                                         params%exp_data(i)%n_weights, &
                                         params%exp_data(i)%data, &
                                         params%exp_data(i)%data_weights, &
                                         params%exp_data(i)%n_data_weights, &
                                         params%exp_data(i)%n_data, &
                                         trim(params%exp_data(i)%file_data_weights))

                  call preprocess_exp_data(params, params%exp_data(i)%x,&
                       & params%exp_data(i)%y, params%exp_data(i)%label,&
                       & state%n_sites, dot_product(cross_product(state%a_box,&
                       & state%b_box), state%c_box)/(dfloat(state%indices(1)*state%indices(2) &
                       &*state%indices(3))), params%exp_data(i)%input, exp_output, .true.)

                  if (params%write_exp .and. .not. params&
                       &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition) then

                     call get_overwrite_condition(params%do_mc,&
                          & params%do_md, loop%mc_istep, loop%md_istep, params&
                          &%write_xyz, overwrite_condition)

                     call write_exp_data(params%exp_data(i)%x, params&
                          &%exp_data(i)%y, overwrite_condition,&
                          & trim(params%exp_data(i)%label)//&
                          & "_exp.dat", params%exp_data(i)%label)
                  end if

               end if

               if (params%exp_data(i)%compute_exp .and. .not. params&
                    &%exp_data(i)%wrote_exp .and. rank == 0 .and. write_condition) then

                  if (params%write_exp) then
                     write (filename, '(A)')&
                          & trim(params%exp_data(i)%label)//"_exp_fit.dat"

                     call get_overwrite_condition(params%do_mc,&
                          & params%do_md, loop%mc_istep, loop%md_istep, params&
                          &%write_xyz, overwrite_condition)

                     call write_exp_data(params%exp_data(i)%x, params%exp_data(i)%y,&
                          & overwrite_condition, trim(params&
                          &%exp_data(i)%label)//"_exp_fit.dat",&
                          & trim(params%exp_data(i)%label)//" : output = "&
                          & //trim(exp_output))

                  end if

               end if

               params%exp_data(i)%wrote_exp = .true.
               params%exp_data(i)%compute_exp = .true.

            end do
         end if

         !###---   XPS Forces and Spectra Prediction   ---###!

         !     Compute core_electron_be energies and forces
         !
         ! Partial sums into the this_ arrays, as for electrostatics.
         call compute_exp_xps(params, state%n_sites, loop%n_xyz, nl%xyz, nl%neighbors_list, nl%n_neigh, &
                              res%local_properties, res%local_properties_cart_der, model%soap_turbo_hypers, &
                              state%a_box, state%b_box, state%c_box, state%indices, dom%i_beg, dom%i_end, dom%j_beg, &
                              dom%j_end, rank, &
                              loop%md_istep, loop%mc_istep, model%valid_xps, model%xps_idx, model%core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              res%this_energies_lp, res%this_forces_lp, res%this_virial_lp, time)

         !###---   (Partial) Pair distribution functions and XRD   ---###!
         !
         ! Partial sums into the this_ arrays; exp_interface allocates them.
#ifdef _GPU
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, &
                                  dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, res%this_energies_pdf, &
                                  res%this_forces_pdf, res%this_virial_pdf, res%this_energies_sf, &
                                  res%this_forces_sf, res%this_virial_sf, res%this_energies_xrd, res%this_forces_xrd, &
                                  res%this_virial_xrd, res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, time, &
                                  i_beg_list, i_end_list, j_beg_list, &
                                  j_end_list, n_omp, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                                  n_sites_temp, n_pairs_temp, write_condition, overwrite_condition, &
                                  temp_string, species_types_actual, state%v_uc)
#else
         call compute_exp_spectra(params, state%n_sites, state%species, state%positions, nl%rjs, nl%xyz, nl%neighbors_list, &
                                  nl%n_neigh, nl%neighbor_species, state%indices, state%a_box, state%b_box, state%c_box, &
                                  dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, rank, ntasks, ierr, loop%md_istep, loop%mc_istep, &
                                  res%this_energies_pdf, res%this_forces_pdf, res%this_virial_pdf, &
                                  res%this_energies_sf, res%this_forces_sf, res%this_virial_sf, &
                                  res%this_energies_xrd, res%this_forces_xrd, res%this_virial_xrd, &
                                  res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, &
                                  time)
#endif

         if (params%do_prediction) then
            !       Two-body, core-potential and three-body contributions, via the
            !       gap_backend seam. The CPU implementation is in
            !       src/gap_backend_cpu.f90; the GPU branch provides the same three
            !       names from src/gap_backend_gpu.f90 and the Makefile picks one.
            call time_start(time%gap)

            call gap_backend_begin(params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                   nl%neighbors_list, dom%i_beg, dom%i_end, dom%j_beg, dom%j_end)

            call add_2b_contribution(model%n_distance_2b, model%distance_2b_hypers, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%energies_2b, res%forces_2b, res%virial_2b, time)

            call add_core_pot_contribution(model%n_core_pot, model%core_pot_hypers, &
                                           params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                           dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                           res%this_virial, &
                                           res%energies_core_pot, res%forces_core_pot, res%virial_core_pot, time)

#ifdef _GPU
            call add_3b_contribution(model%n_angle_3b, model%angle_3b_hypers, nl%neighbors_list, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%forces, res%energies_3b, res%forces_3b, res%virial_3b, time)
#else
            call add_3b_contribution(model%n_angle_3b, model%angle_3b_hypers, nl%neighbors_list, &
                                     params, nl%rjs, nl%xyz, nl%n_neigh, state%species, nl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%energies_3b, res%forces_3b, res%virial_3b, time)
#endif

            call gap_backend_end()

            call time_end(time%gap)
            !       Communicate all energies and forces here for all
            !       terms
            call domain_complete_contributions(dom, comm, res, state, params, model, time)

            !       Add up all the energy terms
            res%energies = res%energies + res%energies_soap + res%energies_2b +&
                 & res%energies_3b + res%energies_core_pot + res%energies_vdw + res%energies_estat !+energies_lp

            if (model%valid_xps) res%energies_exp = res%energies_exp + res%energies_lp
            if (perform%pdf) res%energies_exp = res%energies_exp + res%energies_pdf
            if (perform%sf) res%energies_exp = res%energies_exp + res%energies_sf
            if (perform%xrd) res%energies_exp = res%energies_exp + res%energies_xrd
            if (perform%nd) res%energies_exp = res%energies_exp + res%energies_nd

            if (params%exp_energies) res%energies = res%energies + res%energies_exp

            res%energy_prev = res%energy
            dyn%instant_pressure_prev = dyn%instant_pressure
            res%energy = sum(res%energies)
            res%energy_exp = sum(res%energies_exp)

         end if

         if (.not. params%do_md .and. .not. params%do_mc) then
            if (rank == 0) then
               write (*, *) '                                       |'
               write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
               write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
               write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
               write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
               write (*, '(A,1X,F21.8,1X,A)') ' estat energy:', sum(res%energies_estat), 'eV |'
               write (*, '(A,1X,F22.8,1X,A)') ' Exp. energy:', sum(res%energies_exp), 'eV |'
               if (model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', sum(res%energies_lp), 'eV |'
               if (perform%pdf)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                    & sum(res%energies_pdf), 'eV |'
               if (perform%sf)&
                    & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                    & sum(res%energies_sf), 'eV |'
               if (perform%xrd)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                    & sum(res%energies_xrd), 'eV |'
               if (perform%nd)&
                    & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                    & sum(res%energies_nd), 'eV |'

               if (.not. params%do_mc .or. (params%do_mc .and. loop%mc_istep <= 1)) then
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(res%energies), 'eV |'
               else
                  write (*, '(A,1X,F21.8,1X,A)') ' Total energy:', sum(smp%images(smp%i_trial_image)%energies), 'eV |'
               end if

               if (.not. params%do_mc) then
                  write (*, *) '                                       |'
                  write (*, *) 'Energy & forces in "trajectory_out.xyz"|'
                  write (*, *) '                                       |'
                  write (*, *) '.......................................|'
               else if (loop%mc_istep == 0) then
                  write (*, *) '                                       |'
                  write (*, *) ' MC configs in "mc_current.xyz" and    |'
                  write (*, *) '               "mc_trial.xyz"          |'
                  write (*, *) '               "mc_all.xyz"            |'
                  write (*, *) '.......................................|'
               end if
            end if
         end if

         if (params%do_forces) then
            res%forces = res%forces_soap + res%forces_2b + res%forces_3b + res%forces_core_pot + res%forces_vdw
            res%virial = res%virial_soap + res%virial_2b + res%virial_3b + res%virial_core_pot + res%virial_vdw

!           MAD IR bias. The dipole of this configuration joins the ensemble,
!           the spectrum is compared with the experiment, and the gradient of
!           the mismatch with respect to THIS configuration is added to the
!           forces. Nothing is applied until the ensemble is full, because a
!           partly filled one has a resolution that changes step to step.
!
!           No virial: the bias is a function of the dipole, not of the cell,
!           and a stress from it would be wrong rather than merely missing.
            call ir_after_forces(ir, params, state, res, loop, comm, dyn%md_time, time3, write_condition, time)

            if (model%valid_estat_charges) res%forces = res%forces + res%forces_estat
            if (model%valid_estat_charges) res%virial = res%virial + res%virial_estat

            if (perform%xps_forces) res%forces = res%forces + res%forces_lp
            if (perform%xps_forces) res%virial = res%virial + res%virial_lp

            if (perform%pdf_forces) res%forces = res%forces + res%forces_pdf
            if (perform%pdf_forces) res%virial = res%virial + res%virial_pdf

            if (perform%sf_forces) res%forces = res%forces + res%forces_sf
            if (perform%sf_forces) res%virial = res%virial + res%virial_sf

            if (perform%xrd_forces) res%forces = res%forces + res%forces_xrd
            if (perform%xrd_forces) res%virial = res%virial + res%virial_xrd

            if (perform%nd_forces) res%forces = res%forces + res%forces_nd
            if (perform%nd_forces) res%virial = res%virial + res%virial_nd

            if (rank == 0 .and. params%print_vdw_forces) then
               print *, "> Virial ESTAT "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_estat(i, j)
                  end do
               end do

               print *, "> Virial soap "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_soap(i, j)
                  end do
               end do

               print *, "> Virial 2b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_2b(i, j)
                  end do
               end do

               print *, "> Virial 3b "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_3b(i, j)
                  end do
               end do

               print *, "> Virial core_pot "
               do i = 1, 3
                  do j = 1, 3
                     print *, " i, ", i, " j ", j, " ", res%virial_core_pot(i, j)
                  end do
               end do

               if (perform%xrd_forces) then
                  print *, "> Virial xrd "
                  do i = 1, 3
                     do j = 1, 3
                        print *, " i, ", i, " j ", j, " ", res%virial_xrd(i, j)
                     end do
                  end do
                  temp_string = ""
                  temp_string2 = ""
                  write (temp_string, "(I8)") loop%md_istep
                  write (temp_string2, "(A)") "forces_xrd_"//trim(adjustl(temp_string))
                  open (unit=90, file=temp_string2, status="unknown")
                  do i = 1, state%n_sites
                     write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                        res%forces_xrd(1, i), res%forces_xrd(2, i), res%forces_xrd(3, i)
                  end do
                  close (90)

               end if

            end if

            if (params%print_vdw_forces) then
               open (unit=90, file="forces_vdw", status="unknown")
               do i = 1, state%n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%forces_vdw(1, i), res%forces_vdw(2, i), res%forces_vdw(3, i)
               end do
               close (90)

            end if

            if (rank == 0 .and. params%print_estat_forces) then
               open (unit=90, file="forces_estat", status="unknown")
               do i = 1, state%n_sites
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%forces_estat(1, i), res%forces_estat(2, i), res%forces_estat(3, i)
               end do
               close (90)

               open (unit=90, file="charge_gradients_estat", status="unknown")
               do i = 1, dom%n_atom_pairs_by_rank(rank + 1)
                  write (90, "(F20.8, 1X, F20.8, 1X, F20.8)") &
                     res%local_properties_cart_der(1, i, model%charge_lp_index), &
                     res%local_properties_cart_der(2, i, model%charge_lp_index), &
                     res%local_properties_cart_der(3, i, model%charge_lp_index)
               end do
               close (90)

            end if

         end if
         ! For debugging the virial implementation
         if (rank == 0 .and. .false.) then
            write (*, *) "pressure_soap: ", res%virial_soap/3.d0/state%v_uc
            write (*, *) "pressure_vdw: ", res%virial_vdw/3.d0/state%v_uc
            write (*, *) "pressure_lp: ", res%virial_lp/3.d0/state%v_uc
            write (*, *) "pressure_2b: ", res%virial_2b/3.d0/state%v_uc
            write (*, *) "pressure_3b: ", res%virial_3b/3.d0/state%v_uc
            write (*, *) "pressure_core_pot: ", res%virial_core_pot/3.d0/state%v_uc
         end if
! For debugging the virial implementation
         if (rank == 0 .and. .false.) then
            write (*, *) "pressure_soap: ", res%virial_soap/3.d0/state%v_uc
            write (*, *) "pressure_vdw: ", res%virial_vdw/3.d0/state%v_uc
            do i = 1, 3
               write (*, *) res%virial_vdw(i, :)/state%v_uc
            end do
            write (*, *) "Trace of vdw pressure:", (res%virial_vdw(1, 1) + res%virial_vdw(2, 2) + res%virial_vdw(3, &
                                                                                                                 3))/3.d0/state%v_uc
            write (*, *) "pressure_2b: ", res%virial_2b/3.d0/state%v_uc
            write (*, *) "pressure_3b: ", res%virial_3b/3.d0/state%v_uc
            write (*, *) "pressure_core_pot: ", res%virial_core_pot/3.d0/state%v_uc
            write (*, *) "full vdw forces"
            do i = 1, state%n_sites
               write (*, *) i, res%forces_vdw(1:3, i)
            end do
            write (*, *) "Local virial", res%local_virial_vdw_diag
         end if

         if (params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc) then
            if (rank == 0) then
               !       Write energy and forces if we're just doing static predictions
               !       The masses should be divided by 103.6426965268d0 to have amu units, but
               !       since masses is not allocated for single point calculations, it would
               !       likely lead to a segfault
               call wrap_pbc(state%positions(1:3, 1:state%n_sites), state%a_box&
                    &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
                    & state%c_box/dfloat(state%indices(3)))
               call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                    & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                    &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                    & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                    & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                    & params%do_dipole, res%dipole, res%energies_dipole)

               call write_extxyz(state%n_sites, -loop%n_xyz, dyn%md_time, dyn%time_step,&
                    & dyn%instant_temp, dyn%instant_pressure, state%a_box&
                    &/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)),&
                    & state%c_box/dfloat(state%indices(3)), res%virial, state%xyz_species,&
                    & state%positions(1:3, 1:state%n_sites), state%velocities, res%forces,&
                    & res%energies(1:state%n_sites), state%masses, params&
                    &%write_property, params%write_array_property,&
                    & params%write_local_properties, model%local_property_labels, res%local_properties, &
                    & state%fix_atom, "trajectory_out.xyz", string, .false.,&
                    & params%do_dipole, res%local_dipoles(1:3, 1:state%n_sites))

            end if
         end if
      else
         if (rank == 0) then
            !     Do nothing
            write (*, *) '                                       |'
            write (*, *) 'You didn''t ask me to do anything!      |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if

      !   Do MD stuff here. Moved to src/turbogap_md.f90; the rank guard and the
      !   position broadcast moved with it.
!     In i-PI mode the integrator is i-PI's, so the forces just computed go
!     out over the socket and the next coordinates come back in. Everything
!     compute_md does AROUND the integration -- the skin accounting, the
!     supercell refresh, the broadcast -- happens inside the exchange.
      if (mode == "ipi") then
         call ipi_driver_exchange(rank, state%n_sites, state%positions, state%positions_prev, state%positions_diff, &
                                  state%velocities, state%a_box, state%b_box, state%c_box, state%indices, params%neighbors_buffer, &
                                  res%forces, res%energy, res%virial, loop%exit_loop, nl%rebuild_neighbors_list)
         call domain_sync_after_ipi(dom, comm, state, nl, loop)
      else
         call compute_md(params, rank, ierr, state%n_sites, model%n_species, loop%md_istep, dyn%md_time, dyn%time_step, &
                         state%positions, state%positions_prev, state%positions_diff, state%velocities, res%forces, &
                         state%forces_prev, state%masses, &
                         dyn%masses_types, nl%xyz, state%xyz_species, state%a_box, state%b_box, state%c_box, state%indices, &
                         state%v_uc, res%virial, res%energy, &
                         res%energy_prev, res%energies, res%energies_soap, res%energies_2b, res%energies_3b, &
                         res%energies_core_pot, &
                         res%energies_vdw, res%energies_lp, res%energies_exp, res%energies_pdf, res%energies_sf, res%energies_xrd, &
                         res%energies_nd, res%local_properties, model%local_property_labels, dyn%instant_temp, &
                         dyn%instant_pressure, dyn%instant_pressure_prev, dyn%e_kin, dyn%e_kinetic, dyn%kb, dyn%evpera3tobar, &
                         state%fix_atom, loop%exit_loop, nl%rebuild_neighbors_list, smp%i_image, smp%i_nested, n_pos, dyn%nrows, &
                         filename, string, dyn%allelstopdata, dyn%ephbeta, dyn%ephfdm, dyn%ephlsc, time, &
                         dyn%cum_eel, dyn%gd_istep, &
                         dyn%target_temp, dyn%time_step_prev, res%dipole, res%local_dipoles, res%energies_dipole)
         call domain_sync_after_md(dom, comm, state, nl, params, time)
      end if

      call nested_step(smp, state, res, dyn, nl, params, loop, comm)

      if (rank == 0) then

         if (params%do_mc) then
            if (loop%mc_istep == params%mc_nsteps) then
               loop%exit_loop = .true.
            else
               loop%exit_loop = .false.
            end if

            if (.not. loop%exit_loop .and. ( &
                (loop%md_istep == -1) .or. &
                (params%do_md .and. ( &
                 (loop%md_istep == params%md_nsteps) .or. &
                 ((abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites)) .and. (maxval(abs(res%forces)) < &
                                                                                                  params%f_tol)) &
                 )))) then
               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               call time_start(time%mc)

               !       Now we do a monte-carlo step: we choose what the steps are from the available list and then choose a random number
               !       -- We have the list of move types in params%mc_types and the number params%n_mc_types --
               !       >> First generate a random number in the range of the number of

               if (loop%mc_istep > 0) then
                  !       Evaluate the conditions for acceptance
                  !       > We have the mc conditions in mc.f90
                  !       > We care about comparing e_store to the energy of the new configuration based on the mc_movw

                  ! Reset the parameters for md / relaxation
                  smp%trial_came_from_md = params%do_md
                  if (params%do_md) then
                     loop%md_istep = -1
                     params%do_md = .false.
                     smp%do_mc_relax = .false.
                     ! Assume that the number of steps has already been set.
                  end if

                  if (.not. params%mc_hamiltonian) dyn%E_kinetic = 0.d0

!                 The trial configuration has to be the one `energy` belongs to.
!                 An "md" move, or a relaxation after any other move, leaves
!                 `positions` one integrator step *past* the last force
!                 evaluation: compute_md advances them after the energy was
!                 computed, and stashes the configuration it was computed at in
!                 positions_prev. md.f90 says as much -- "velocities and
!                 positions_prev are synchronous, positions is dt ahead of
!                 velocities" -- and compute_md writes positions_prev to the
!                 trajectory for exactly this reason.
!
!                 Storing `positions` here accepted or rejected x_(n+1) on the
!                 strength of E(x_n), and wrote a frame to mc_all.xyz whose
!                 energy was not the energy of its own coordinates: ~1 eV out on
!                 512 atoms after a 0.5 fs velocity-Verlet burst. Rewind by the
!                 one step, which also pairs the stored positions with the
!                 stored velocities.
                  if (smp%trial_came_from_md) then
                     state%positions(1:3, 1:state%n_sites) = state%positions_prev(1:3, 1:state%n_sites)
                  end if

                  call from_properties_to_image(smp%images(smp%i_trial_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, dyn%E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

                  if (params%verb > 50) write (*, *) '.......................................|'
                  if (params%verb > 50) write (*, '(A,1X,I0)') ' MC Iteration:', loop%mc_istep
                  if (params%verb > 50) write (*, '(A,1X,A)') '    Move type:', smp%mc_move

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_prev:', smp%images(smp%i_current_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_prev:', smp%images(smp%i_current_image)&
                       &%energy + smp%images(smp%i_current_image)%e_kin

                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Ekin_new:', smp%images(smp%i_trial_image)%e_kin
                  if (params%verb > 50) write (*, '(A,1X,F22.8)') '   &
                       & Etot_new :', smp%images(smp%i_trial_image)%energy &
                       &+ smp%images(smp%i_trial_image)%e_kin

                  state%v_uc = dot_product(cross_product(state%a_box, state%b_box), &
                                           state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))

                  if (params%accessible_volume) then
                     call get_accessible_volume(state%v_uc, smp%v_a_uc, state%species, params%radii)
                     if (params%verb > 50) write (*, '(A,F12.6,A,F12.6&
                          &,1X,A)') ' V_acc new: ', smp%v_a_uc, ' A^3&
                          & V_acc old ', smp%v_a_uc_prev, 'A^3 |'
                  else
                     smp%v_a_uc = state%v_uc
                  end if

                  call get_mc_acceptance(smp%mc_move, smp%p_accept, &
                       res%energy + dyn%E_kinetic, &
                       smp%images(smp%i_current_image)%energy + smp%images(smp%i_current_image)%e_kin, &
                       params%t_beg, smp%mc_mu_id, &
                       params%mc_mu, smp%n_mc_species, state%v_uc, smp%v_uc_prev,&
                       & smp%v_a_uc, smp%v_a_uc_prev, params%mc_exchange_mass, &
                       & params%mc_exchange_e0, params%mc_mu_reference, &
                       & params%p_beg, state%n_sites)

!                 call get_mc_acceptance(mc_move, p_accept, &
!                      energy + E_kinetic, &
!                      images(i_current_image)%energy + images(i_current_image)%e_kin, &
!                      params%t_beg, &
!                      params%mc_mu(mc_mu_id), n_mc_species(mc_mu_id), v_uc, v_uc_prev,&
!                      & v_a_uc, v_a_uc_prev, params&
!                      &%masses_types(mc_id(mc_mu_id)), params%p_beg)

                  call random_number(smp%ranf)

                  if (smp%mc_move == "insertion") smp%n_mc_species(smp%mc_mu_id) = smp%n_mc_species(smp%mc_mu_id) + 1
                  if (smp%mc_move == "removal") smp%n_mc_species(smp%mc_mu_id) = smp%n_mc_species(smp%mc_mu_id) - 1

                  !    ACCEPT OR REJECT
                  if (params%verb > 50) write (*, '(A,1X,A,1X,A,L4,1X&
                       &,A,ES12.6,1X,A,1X,ES12.6)') 'Is ',&
                       & trim(smp%mc_move), 'accepted?', smp%p_accept >&
                       & smp%ranf, ' p_accept =', smp%p_accept, ' ranf = ',&
                       & smp%ranf

                  if (loop%mc_istep == 1) then
                     open (unit=200, file="mc.log", status="unknown")
                     if (res%energy_exp > 0.d0) then
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current             E_exp_trial&
                             &          E_exp_current  N_tot_trial &
                             & N_mc_species_trial'
                     else
                        write (200, '(A)') '# mc_istep  mc_move &
                             & accepted  E_trial              E_current &
                             &          N_tot_trial  N_mc_species_trial'
                     end if

                  end if
                  if (loop%mc_istep > 1) then
                     open (unit=200, file="mc.log", status="old", position="append")
                  end if

                  ! collect the strings for the species etc
                  temp_string = ""
                  temp_string2 = ""

                  do i = 1, params%n_mc_mu
                     temp_string = ""
                     write (temp_string, "(A,1X,I8)") trim(params%mc_species(i)), smp%n_mc_species(i)
                     temp_string2 = trim(temp_string2)//" "//trim(temp_string)
                  end do

                  if (res%energy_exp > 0.d0) then

                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                          loop%mc_istep, trim(adjustl(smp%mc_move)), smp%p_accept > smp%ranf, res%energy + dyn%E_kinetic, &
                          smp%images(smp%i_current_image)%energy +&
                          & smp%images(smp%i_current_image)%e_kin, res%energy_exp,&
                          & smp%images(smp%i_current_image)%energy_exp,&
                          & smp%images(smp%i_trial_image)%n_sites,&
                          & trim(temp_string2)
                  else
                     write (200, "(I8, 1X, A10, 1X, L4, 1X, F20.8, 1X, F20.8, 1X, I8, 1X, A)") &
                        loop%mc_istep, trim(adjustl(smp%mc_move)), smp%p_accept > smp%ranf, res%energy + dyn%E_kinetic, &
                        smp%images(smp%i_current_image)%energy + smp%images(smp%i_current_image)%e_kin, &
                        smp%images(smp%i_trial_image)%n_sites, trim(temp_string2)

                  end if

                  if (loop%mc_istep >= 1) close (200)

                  if (smp%p_accept > smp%ranf) then
                     !             Accept
                     ! Set variables
                     loop%n_sites_prev = state%n_sites
                     smp%v_uc_prev = state%v_uc
                     smp%v_a_uc_prev = smp%v_a_uc
                     smp%virial_prev = res%virial
                     !   Assigning the default image with the accepted one
                     smp%images(smp%i_current_image) = smp%images(smp%i_trial_image)

                     if (params%n_mc_mu > 0) then
                        smp%n_mc_species_prev = smp%n_mc_species
                     end if

                  end if
                  if (state%n_sites > 1) then
                     dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
                     dyn%instant_pressure = (dyn%kB*dfloat(state%n_sites - 1)*dyn%instant_temp&
                          &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                          &/state%v_uc*dyn%eVperA3tobar
                  else
                     dyn%instant_temp = 0.0d0
                     dyn%instant_pressure = 0.0d0
                  end if

                  if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') '&
                          & Writing mc_current.xyz and&
                          & mc_all.xyz '
                     call wrap_pbc(smp%images(smp%i_current_image)&
                          &%positions(1:3,&
                          & 1:smp%images(smp%i_current_image)%n_sites),&
                          & smp%images(smp%i_current_image)%a_box&
                          &/dfloat(state%indices(1)),&
                          & smp%images(smp%i_current_image)%b_box&
                          &/dfloat(state%indices(2)),&
                          & smp%images(smp%i_current_image)%c_box&
                          &/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf,&
                          & params%valid_xrd, params%valid_nd,&
                          & params%do_pair_distribution, params&
                          &%do_structure_factor, params%do_xrd,&
                          & params%do_nd, string, params%do_dipole,&
                          & smp%images(smp%i_current_image)%dipole,&
                          & smp%images(smp%i_current_image)%energies_dipole)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 0, 1.0d0, 0.d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 1, 1.0d0, 0.d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels,&
                          & smp%images(smp%i_current_image)%local_properties,&
                          & smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .false., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                  end if

                  !          Add acceptance to the log file else dont
                  call time_end(time%mc)

               else ! if (mc_istep == 0)
                  smp%temp_md_nsteps = params%md_nsteps
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, *) 'Starting MC, using parameters:         |'
                  if (params%verb > 50) write (*, *) '                                       |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nsteps     = ', params%mc_nsteps, '     &
                       &        |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_types    = ', params%n_mc_types, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_types:                              |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '     ', params%mc_types(i), '|'
                  end do
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_accept_ratio:                       |'
                  do i = 1, params%n_mc_types
                     if (params%verb > 50) write (*, '(1X,A,1X,F12.8,1X&
                          &,A)') '   ', params%mc_acceptance(i), '    &
                          &                  |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'n_mc_swaps    = ', params%n_mc_swaps, '    &
                       &         |'
                  if (params%verb > 50) write (*, '(1X,A)') 'mc_swaps:  &
                       &                            |'
                  do i = 1, 2*params%n_mc_swaps
                     if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')&
                          & '   ', params%mc_swaps(i), '              &
                          &        |'
                  end do
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_move_max   = ', params%mc_move_max, &
                       & 'A   |'

                  do i = 1, params%n_mc_mu
                     write (*, '(1X,A,1X,F17.8,1X,A)') 'mc_mu         = ', params%mc_mu(1), 'eV  |'
                     write (*, '(1X,A,1X,A,1X,A)') 'mc_species    = ', trim(params%mc_species(i)), '                    |'
                  end do

                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_min_dist   = ', params%mc_min_dist, &
                       & 'A   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,F17.8,1X&
                       &,A)') 'mc_lnvol_max  = ', params%mc_lnvol_max,&
                       & '    |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_write_xyz  = ', params%mc_write_xyz, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_relax      = ', params%mc_relax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,I8,1X,A)')  &
                       &  'mc_nrelax     = ', params%mc_nrelax, '  &
                       &           |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_relax_opt  = ', params%mc_relax_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,A,1X,A)')   &
                       &  'mc_hybrid_opt = ', params%mc_hybrid_opt, '  &
                       &   |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_optimize_exp = ', params%mc_optimize_exp&
                       &, '  |'
                  if (params%verb > 50) write (*, '(1X,A,1X,L8,1X,A)')  &
                       &  'mc_hamiltonian = ', params%mc_hamiltonian, '&
                       &  |'

                  if (params%verb > 50) write (*, *) '                                       |'
                  ! t_beg must

                  if (.not. allocated(smp%images) .and. .not. params%do_nested_sampling) then
                     allocate (smp%images(1:2))
                  else if (.not. allocated(smp%images) .and. params%do_nested_sampling) then
                     allocate (smp%images(1:2*smp%i_image))
                  end if

                  if (.not. allocated(smp%mc_mol_id)) then
                     allocate (smp%mc_mol_id(1:state%n_sites), smp%mc_mol_mu(1:state%n_sites))
                     smp%mc_mol_id = 0
                     smp%mc_mol_mu = 0
                  end if

                  if (.not. allocated(smp%mc_id) .and. params%n_mc_mu > 0) then
                     allocate (smp%mc_id(1:params%n_mc_mu))
                     allocate (smp%n_mc_species(1:params%n_mc_mu))
                     allocate (smp%n_mc_species_prev(1:params%n_mc_mu))

                     smp%mc_id = 1
                     smp%n_mc_species = 0

                     !    get the mc species types

                     do j = 1, params%n_mc_mu
                        do i = 1, model%n_species
                           if (params%species_types(i) == params%mc_species(j)) then
                              smp%mc_id(j) = i
                           end if
                        end do
                     end do
                  end if

                  !       Now use the image construct to store this as the image to compare to
                  call from_properties_to_image(smp%images(smp%i_current_image), state%positions, state%velocities, state%masses, &
                                                res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                                res%energy_exp, dyn%E_kinetic, &
                                                state%species, state%species_supercell, state%n_sites, state%indices, &
                                                state%fix_atom, &
                                                state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                                res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

                  dyn%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/dyn%kB*dyn%E_kinetic
                  dyn%instant_pressure = (dyn%kB*dfloat(state%n_sites - 1)*dyn%instant_temp&
                       &+ (res%virial(1, 1) + res%virial(2, 2) + res%virial(3, 3))/3.d0)&
                       &/state%v_uc*dyn%eVperA3tobar

                  if ((loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                       modulo(loop%mc_istep, params%write_xyz) == 0)) then
                     if (params%verb > 50) write (*, '(1X,A)') ' Writing mc_current.xyz and mc_all.xyz '
                     call wrap_pbc(smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                                   smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                                   smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                                   smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)))
                     call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                          & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                          &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                          & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                          & params%do_structure_factor, params%do_xrd, params%do_nd, string,&
                          & params%do_dipole, smp%images(smp%i_current_image)%dipole,&
                          & smp%images(smp%i_current_image)%energies_dipole)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 0, 1.0d0, 0.0d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_current.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     call write_extxyz(smp%images(smp%i_current_image)%n_sites, 1, 1.0d0, 0.0d0, dyn%instant_temp, &
                        dyn%instant_pressure, &
                          smp%images(smp%i_current_image)%a_box/dfloat(state%indices(1)), &
                          smp%images(smp%i_current_image)%b_box/dfloat(state%indices(2)), &
                          smp%images(smp%i_current_image)%c_box/dfloat(state%indices(3)), &
                          smp%virial_prev, smp%images(smp%i_current_image)%xyz_species, &
                          smp%images(smp%i_current_image)%positions(1:3, 1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%velocities, &
                          smp%images(smp%i_current_image)%forces, &
                          smp%images(smp%i_current_image)%energies(1:smp%images(smp%i_current_image)%n_sites), &
                          smp%images(smp%i_current_image)%masses, &
                          params%write_property, params&
                          &%write_array_property, params&
                          &%write_local_properties,&
                          & model%local_property_labels, smp%images(smp%i_current_image)%local_properties&
                          &, smp%images(smp%i_current_image)%fix_atom,&
                          & "mc_all.xyz", string, .true., &
                          & params%do_dipole,&
                          & smp%images(smp%i_current_image)%local_dipoles)

                     smp%v_uc_prev = dot_product(cross_product(state%a_box, state%b_box), &
                                                 state%c_box)/(dfloat(state%indices(1)*state%indices(2)*state%indices(3)))
                     if (params%accessible_volume) then
                        call get_accessible_volume(smp%v_uc_prev, smp%v_a_uc_prev, state%species, params%radii)
                     else
                        smp%v_a_uc_prev = smp%v_uc_prev
                     end if
                  end if

               end if

               !  Now start the mc logic: first, use the stored images properties
               call from_image_to_properties(smp%images(smp%i_current_image), state%positions, state%velocities, state%masses, &
                                             res%forces, state%a_box, state%b_box, state%c_box, res%energy, res%energies, &
                                             res%energy_exp, dyn%E_kinetic, &
                                             state%species, state%species_supercell, state%n_sites, state%indices, state%fix_atom, &
                                             state%xyz_species, state%xyz_species_supercell, res%local_properties, &
                                             res%local_dipoles, res%energies_dipole, res%dipole, smp%mc_mol_id, smp%mc_mol_mu)

               call perform_mc_step(&
                    & state%positions, state%species, state%xyz_species, state%masses, state%fix_atom,&
                    & state%velocities, state%positions_prev, state%positions_diff, smp%disp, smp%d_disp, &
                       params%n_local_properties,&
                    & params%mc_acceptance, params%mc_mu_acceptance, res%local_properties, &
                    smp%images(smp%i_current_image)%local_properties, res%energies,&
                    & res%forces, state%forces_prev, state%n_sites, params%n_mc_mu, smp%mc_mu_id, smp%n_mc_species,&
                    & smp%mc_move, params%mc_species,&
                    & params%mc_move_max, params%mc_min_dist, params%mc_max_dist, params%mc_max_insertion_trials, &
                    params%mc_lnvol_max, params%mc_types, params%masses_types, smp%species_idx,&
                    & smp%images(smp%i_current_image)%positions,&
                    & smp%images(smp%i_current_image)%species,&
                    & smp%images(smp%i_current_image)%xyz_species,&
                    & smp%images(smp%i_current_image)%fix_atom,&
                    & smp%images(smp%i_current_image)%masses, state%a_box(1:3), state%b_box(1:3),&
                    & state%c_box(1:3), state%indices, params%do_md, params%mc_relax,&
                    & loop%md_istep, smp%mc_id, dyn%E_kinetic, dyn%instant_temp, params%t_beg,&
                    & params%n_mc_swaps, params%mc_swaps, params%mc_swaps_id, &
                    & params%species_types, params%mc_hamiltonian,&
                    & params%n_mc_relax_after, params&
                    &%mc_relax_after, smp%do_mc_relax, params%verb, &
                    params%mc_n_planes, params%mc_planes, params%mc_max_dist_to_planes, &
                    params%mc_planes_restrict_to_polyhedron, &
                    params%mc_molecules, smp%mc_mol_id, smp%mc_mol_mu, &
                    smp%images(smp%i_current_image)%mc_mol_id, smp%images(smp%i_current_image)%mc_mol_mu, smp%mc_mol_next)

               nl%rebuild_neighbors_list = .true.

               ! NOTE: the species_supercell and xyz_species_supercell are
               ! not commensurate with the new image as these have not been
               ! calculated. If reading from an outputted xyz file, then it
               ! should be okay but really the new atoms should be added to
               ! the supercell in the usual way, but for convenience, one has
               ! not done that.

               if (params%mc_relax .and. smp%do_mc_relax) then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_relax_opt
                  params%md_nsteps = params%mc_nrelax

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, dyn%E_kinetic, state%masses, dyn%instant_temp, &
                                            params%t_beg, &
                                            params%velocity_distribution)

                  if (params%mc_hamiltonian) dyn%E_kinetic_prev = dyn%E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if
               ! If doing md, don't relax
               if (smp%mc_move == 'md') then
                  ! Set the parameters for relaxatrino
                  loop%md_istep = -1
                  params%do_md = .true.
                  params%optimize = params%mc_hybrid_opt
                  params%md_nsteps = smp%temp_md_nsteps

                  if (state%n_sites == 1) then
                     params%do_md = .false.
                  end if

                  call randomize_velocities(state%velocities, state%n_sites, dyn%E_kinetic, state%masses, dyn%instant_temp, &
                                            params%t_beg, &
                                            params%velocity_distribution)
                  if (params%mc_hamiltonian) dyn%E_kinetic_prev = dyn%E_kinetic
                  ! Note, that this may override md steps if the same is chosen! More testing needed
               end if

               if ((params%mc_write_xyz .or. loop%mc_istep == 0 .or. loop%mc_istep == params%mc_nsteps .or. &
                    modulo(loop%mc_istep, params%write_xyz) == 0)) then

                  call wrap_pbc(state%positions(1:3, 1:state%n_sites), &
                                state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                                state%c_box/dfloat(state%indices(3)))
                  call get_xyz_energy_string(res%energies_soap, res%energies_2b,&
                       & res%energies_3b, res%energies_core_pot, res%energies_vdw, res%energies_exp&
                       &, res%energies_lp, res%energies_pdf, res%energies_sf, res%energies_xrd, res%energies_nd,&
                       & params%valid_pdf, params%valid_sf, params%valid_xrd, params%valid_nd, params%do_pair_distribution,&
                       & params%do_structure_factor, params%do_xrd, params%do_nd, string)

                  call write_extxyz(state%n_sites, 0, 1.0d0, 0.0d0, dyn%instant_temp, dyn%instant_pressure, &
                       state%a_box/dfloat(state%indices(1)), state%b_box/dfloat(state%indices(2)), &
                          state%c_box/dfloat(state%indices(3)), &
                       res%virial, state%xyz_species, &
                       state%positions(1:3, 1:state%n_sites), state%velocities, &
                       res%forces, res%energies(1:state%n_sites), state%masses, &
                       params%write_property, params&
                       &%write_array_property, params&
                       &%write_local_properties,&
                       & model%local_property_labels, res%local_properties&
                       &, state%fix_atom, smp%mc_file, string, .true.)
               end if
               ! As we have moved/added/removed, we must check the supercell and  broadcast the results

               call read_xyz(smp%mc_file, .true., params%all_atoms, params%do_timing, &
                             model%n_species, params%species_types, loop%repeat_xyz, model%rcut_max, params%which_atom, &
                             state%positions, params%do_md, state%velocities, params%masses_types, state%masses, &
                             state%xyz_species, &
                             state%xyz_species_supercell, state%species, state%species_supercell, state%indices, state%a_box, &
                             state%b_box, state%c_box, &
                             state%n_sites, .true., state%fix_atom, params%t_beg, &
                             params%write_array_property(6), .true., params%randomize_velocities)

            else
               if (smp%mc_move == 'md') then
                  if (params%print_progress .and. loop%md_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%md_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%md_istep == params%md_nsteps - 1 .or. &
                           (abs(res%energy - res%energy_prev) < params%e_tol*dfloat(state%n_sites) .and. &
                            maxval(abs(res%forces)) < params%f_tol) .and. loop%md_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%md_istep < params%md_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%md_istep + 1)/params%md_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%md_istep + 1)/params%md_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%mc_hamiltonian) then
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step: H =&
                          & T + V = ", res%energy + dyn%E_kinetic, ",&
                          & iteration ", loop%md_istep, "/", params&
                          &%md_nsteps
                  else
                     if (params%verb > 50) write (*, '(1X,A,1X,F20.8,1X&
                          &,A,1X,I8,1X,A,1X,I8)') "Hybrid md step:&
                          & energy = ", res%energy, ", iteration ",&
                          & loop%md_istep, "/", params%md_nsteps
                  end if

                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') '&
                       & core_pot energy:', sum(res%energies_core_pot),&
                       & 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') '&
                       & vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               else
                  if (params%print_progress .and. loop%mc_istep == 0) then
                     write (*, *) '                                       |'
                     write (*, *) 'Progress:                              |'
                     write (*, *) '                                       |'
                     write (*, '(1X,A)', advance='no') '[                                    ] |'
                     loop%update_bar = params%mc_nsteps/36
                     if (loop%update_bar < 1) then
                        loop%update_bar = 1
                     end if
                     loop%counter = 1
                  else if (loop%mc_istep == params%mc_nsteps - 1 .and. loop%mc_istep > 0) then
                     write (*, *)
                  else if (params%print_progress .and. loop%counter == loop%update_bar .and. loop%mc_istep < params%mc_nsteps &
                           - 1) then
                     do j = 1, 36 + 3
                        write (*, "(A)", advance="no") creturn
                     end do
                     write (*, "(1X,A)", advance="no") "["
                     do i = 1, 36*(loop%mc_istep + 1)/params%mc_nsteps
                        write (*, "(A)", advance="no") "."
                     end do
                     do i = 36*(loop%mc_istep + 1)/params%mc_nsteps + 1, 36
                        write (*, "(A)", advance="no") " "
                     end do
                     write (*, "(A)", advance="no") "] |"
                     loop%counter = 1
                  else
                     loop%counter = loop%counter + 1
                  end if

                  if (params%verb > 50 .and. smp%do_mc_relax) write (*, '(1X,A,1X,F20.8,1X,A&
                       &,1X,I8,1X,A,1X,I8)') "MC Relax md step: energy &
                       &= ", res%energy, ", iteration ", loop%md_istep, "/",&
                       & params%mc_nrelax
                  if (params%verb > 50) write (*, '(A,1X,F22.8,1X,A)') ' SOAP energy:', sum(res%energies_soap), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 2b energy:', sum(res%energies_2b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F24.8,1X,A)') ' 3b energy:', sum(res%energies_3b), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F18.8,1X,A)') ' core_pot energy:', sum(res%energies_core_pot), 'eV |'
                  if (params%verb > 50) write (*, '(A,1X,F23.8,1X,A)') ' vdw energy:', sum(res%energies_vdw), 'eV |'
                  if (params%verb > 50 .and. model%valid_xps) write (*, '(A,1X,F23.8,1X,A)') ' xps energy:', &
                     sum(res%energies_lp), 'eV |'

                  if (perform%pdf .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' pdf energy:',&
                       & sum(res%energies_pdf), 'eV |'
                  if (perform%sf .and. params%verb > 50)&
                       & write (*, '(A,1X,F24.8,1X,A)') ' sf energy:',&
                       & sum(res%energies_sf), 'eV |'
                  if (perform%xrd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' xrd energy:',&
                       & sum(res%energies_xrd), 'eV |'
                  if (perform%nd .and. params%verb > 50)&
                       & write (*, '(A,1X,F23.8,1X,A)') ' nd energy:',&
                       & sum(res%energies_nd), 'eV |'

               end if
            end if
         end if

      end if

      ! NOTE!! One tried for far far too long to be smart and implement some
      ! sort of conditional broadcasting: having a logical array named
      ! broadcast, which perform_mc_step would then to set values to
      ! true. Specific indexes referenced specific quantities to be
      ! broadcasted, which allowed for the broadcasting amount to be
      ! dependent on the step, e.g. if it were an insertion step then
      ! positions, masses, n_sites, etc would have to be broadcast, whereas
      ! for a simple move only positions had to be broadcasted. This array
      ! would then subsequently be broadcast to all other ranks, thereby
      ! allowing for the minimum number of allocations and
      ! communication. BUT, for some reason, this led to segfaults
      ! (corrupted unsorted chunks or something of that sort).

      ! This doesn't make sense to be as all ranks have the same broadcast
      ! array (as it is broadcasted before) so it seems like it should work
      ! but it does not! Hence, in the following broadcasting, everything is
      ! transmitted.

      ! This can be optimised, so please do if you are smarter than me

      call time_start(time%mpi)
      call comm_bcast(comm, params%do_md)
      call comm_bcast(comm, loop%md_istep)
      call time_end(time%mpi)
      call domain_sync_state(dom, comm, state, params, time)
      !   Now that all ranks know the size of n_sites, we allocate do_list
      if (.not. params%do_md .or. (params%do_md .and. loop%md_istep == 0) .or. &
          (params%do_mc)) then
         if (allocated(dom%do_list)) deallocate (dom%do_list)
         allocate (dom%do_list(1:state%n_sites))
         dom%do_list = .true.
      end if
      call get_time(time1)
      !   Parallel neighbors list build
      call comm_bcast(comm, nl%rebuild_neighbors_list)

      if (nl%rebuild_neighbors_list) then
         deallocate (nl%rjs, nl%xyz, nl%thetas, nl%phis, nl%neighbor_species)
         deallocate (nl%neighbors_list, nl%n_neigh)
         deallocate (nl%n_neigh_local)
      end if
      if ((params%do_nested_sampling .and. .not. params%do_mc) .and. &
          (params%do_md .and. (loop%md_istep == params%md_nsteps .or. loop%exit_loop))) then
         deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, state%species_supercell, &
                     dom%do_list)
         if (allocated(state%velocities)) deallocate (state%velocities)
      end if
      if (params%do_mc .and. params%do_md) then
         if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop)) then
            deallocate (state%positions, state%xyz_species, state%xyz_species_supercell, state%species, &
                        state%species_supercell, dom%do_list)
            if (allocated(state%velocities)) deallocate (state%velocities)
         end if
      end if

      if ((params%do_md .and. .not. params%do_mc) .and. &
          (loop%md_istep == params%md_nsteps .or. loop%exit_loop) .and. rank == 0) then
         deallocate (state%positions_prev, state%forces_prev)
      end if
      if (params%do_mc .and. (loop%mc_istep == params%mc_nsteps .or. loop%exit_loop) .and. rank == 0) then
         if (allocated(state%forces_prev)) deallocate (state%forces_prev)
         if (allocated(state%positions_prev)) deallocate (state%positions_prev)
      end if

      if (params%exp_forces .and. (loop%md_istep == params%md_nsteps .or.&
           & loop%mc_istep == params%mc_nsteps .or. loop%exit_loop)) then
         do i = 1, params%n_exp
            if (allocated(params%exp_data(i)%x)) deallocate (params%exp_data(i)%x)
            if (allocated(params%exp_data(i)%y)) deallocate (params%exp_data(i)%y)
            if (allocated(params%exp_data(i)%y_pred)) deallocate (params%exp_data(i)%y_pred)
         end do
      end if

      call ir_step_end(ir, params, loop)

      if (.not. params%do_mc) loop%n_sites_prev = state%n_sites
      dom%n_atom_pairs_by_rank_prev = dom%n_atom_pairs_by_rank(rank + 1)

      call comm_bcast(comm, loop%exit_loop)
      if (loop%exit_loop) exit
      ! End of loop through structures in the xyz file or MD steps
   end do

!  i-PI has said EXIT, or something else ended the loop. Close the socket
!  before the reports below, so that i-PI sees the driver leave cleanly
!  rather than timing out on a half-open connection.
   if (mode == "ipi") call ipi_driver_close(rank)

   call ir_finish(ir, params, comm, time)

   if (params%do_md .or. params%do_prediction .or. params%do_mc) then
      call get_time(time2)
      if (rank == 0) then
         if (params%do_md .and. .not. params%do_nested_sampling) then
            write (*, *) '                                       |'
            write (*, '(I8,A,F13.3,A)') loop%md_istep, ' MD steps:', time2 - time3, ' seconds |'
         end if
         if (params%do_mc) then
            write (*, *)
            write (*, *) '                                       |'
            write (*, '(I8,A,F13.3,A)') loop%mc_istep, ' MC steps:', time2 - time3, ' seconds |'
         end if

         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *          Setup:', time%setup(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - input+pot.:', time%read_input(3), ' seconds |'
         if (comm_with_mpi) then
            write (*, '(A,F13.3,A)') '     -  MPI setup:', time%mpi_setup(3), ' seconds |'
         end if
         write (*, '(A,F13.3,A)') ' * Read XYZ files:', time%read_xyz(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' * Neighbor lists:', time%neigh(3), ' seconds |'
         write (*, '(A,F13.3,A)') ' *  GAP desc/pred:', time%gap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     - soap_turbo:', time%soap(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         2b:', time%gap_2b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -         3b:', time%gap_3b(3), ' seconds |'
         write (*, '(A,F13.3,A)') '     -   core_pot:', time%gap_core_pot(3), ' seconds |'
!       vdw is a parent, not one of the GAP children above: compute_vdw runs
!       outside the time%gap region and sum_times adds it in its own right.
!       Printing it indented under GAP said otherwise.
         if (params%vdw_type /= "none") then
            write (*, '(A,F13.3,A)') ' *            vdw:', time%vdw(3), ' seconds |'
         end if
         if (model%valid_xps .or. params%do_pair_distribution .or. params&
              &%do_structure_factor .or. params%do_xrd .or. params%do_nd) write (*, '(A&
              &,F13.3,A)') ' *  Exp. pred.   :', time%pdf(3) + time%sf(3) + time%xrd(3) + time%nd(3), ' seconds&
              & |'
         if (model%valid_xps) write (*, '(A,F13.3,A)') '     -        xps:',&
              & time%xps(3), ' seconds |'
         if (params%do_pair_distribution) write (*, '(A,F13.3,A)') '     -        pdf:', time%pdf(3), ' seconds |'
         if (params%do_structure_factor) write (*, '(A,F13.3,A)') '     -         sf:', time%sf(3), ' seconds |'
         if (params%do_xrd) write (*, '(A,F13.3,A)') '     -        xrd:', time%xrd(3), ' seconds |'
         if (params%do_nd) write (*, '(A,F13.3,A)') '     -         nd:', time%nd(3), ' seconds |'

         call ir_report(ir, params, time)

         if (do_electrostatics) then
            write (*, '(A,F13.3,A)') ' * Electrostatics:', time%estat(3), ' seconds |'
         end if
         if (params%do_md) then
            write (*, '(A,F13.3,A)') ' *  MD algorithms:', time%md(3), ' seconds |'
         end if
         if (params%do_mc) then
            write (*, '(A,F13.3,A)') ' *  MC algorithms:', time%mc(3), ' seconds |'
         end if

         if (comm_with_mpi) then
            write (*, '(A,F13.3,A)') ' *  MPI comms.   :', time%mpi(3) + time%mpi_positions(3) + time%mpi_ef(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -  pos & vel:', time%mpi_positions(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     - E & F brc.:', time%mpi_ef(3), ' seconds |'
            write (*, '(A,F13.3,A)') '     -  MPI misc.:', time%mpi(3), ' seconds |'
         end if
!       Miscellaneous is what the parent buckets do not account for.  It used
!       to be written out here as one long subtraction, which is how it came to
!       subtract time%gap and the mpi_ef reduce nested inside it and print a
!       negative number.  sum_times owns the list now (src/timing.f90), so the
!       set summed here and the set declared as parents there cannot disagree.
!
!       Accounted-for is printed beside it so the arithmetic is visible: the
!       three numbers below have to add up, and a reader can see at a glance how
!       much of the run the buckets actually name.  A large Miscellaneous is a
!       statement that something real is not being measured -- which is how the
!       setup bucket above came to exist.
         time%total(3) = time2 - time3
         write (*, *) '                                       |'
         write (*, '(A,F13.3,A)') ' *  Accounted for:', sum_times(time), ' seconds |'
         write (*, '(A,F13.3,A)') ' *  Miscellaneous:', time%total(3) - sum_times(time), ' seconds |'
         write (*, '(A,F13.3,A)') ' *     Total time:', time%total(3), ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end if

#ifdef _GPU
   do i = 1, model%n_soap_turbo
      if (.not. model%soap_turbo_hypers(i)%recompute_basis) then
         call gpu_free_async(model%soap_turbo_hypers(i)%W_d, gpu_stream)
         call gpu_free_async(model%soap_turbo_hypers(i)%S_d, gpu_stream)
         call gpu_free_async(model%soap_turbo_hypers(i)%multiplicity_array_d, gpu_stream)
      end if
   end do
#endif
   if (allocated(state%fix_atom)) deallocate (state%fix_atom)
   if (allocated(state%positions)) deallocate (state%positions)
   if (allocated(state%velocities)) deallocate (state%velocities)
   if (allocated(state%positions_diff)) deallocate (state%positions_diff)

   if (allocated(res%energies)) deallocate (res%energies)
   if (allocated(res%local_dipoles)) deallocate (res%local_dipoles, res%this_local_dipoles)
   if (allocated(res%energies_dipole)) deallocate (res%energies_dipole, res%this_energies_dipole)
   if (allocated(res%energies_soap)) deallocate (res%energies_soap)
   if (allocated(res%energies_2b)) deallocate (res%energies_2b)
   if (allocated(res%energies_3b)) deallocate (res%energies_3b)
   if (allocated(res%energies_core_pot)) deallocate (res%energies_core_pot)
   if (allocated(res%energies_vdw)) deallocate (res%energies_vdw)
   if (allocated(res%energies_exp)) deallocate (res%energies_exp)
   if (allocated(res%energies_lp)) deallocate (res%energies_lp)
   if (allocated(res%energies_pdf)) deallocate (res%energies_pdf)
   if (allocated(res%energies_sf)) deallocate (res%energies_sf)
   if (allocated(res%energies_xrd)) deallocate (res%energies_xrd)
   if (allocated(res%energies_nd)) deallocate (res%energies_nd)

   if (allocated(res%this_energies)) deallocate (res%this_energies)
   if (allocated(res%this_energies_vdw)) deallocate (res%this_energies_vdw)
   if (allocated(res%this_energies_lp)) deallocate (res%this_energies_lp)
   if (allocated(res%this_energies_pdf)) deallocate (res%this_energies_pdf)
   if (allocated(res%this_energies_sf)) deallocate (res%this_energies_sf)
   if (allocated(res%this_energies_xrd)) deallocate (res%this_energies_xrd)
   if (allocated(res%this_energies_nd)) deallocate (res%this_energies_nd)

   if (allocated(res%forces)) deallocate (res%forces)
   if (allocated(res%forces_soap)) deallocate (res%forces_soap)
   if (allocated(res%forces_2b)) deallocate (res%forces_2b)
   if (allocated(res%forces_3b)) deallocate (res%forces_3b)
   if (allocated(res%forces_core_pot)) deallocate (res%forces_core_pot)
   if (allocated(res%forces_vdw)) deallocate (res%forces_vdw)
   if (allocated(res%forces_lp)) deallocate (res%forces_lp)
   if (allocated(res%forces_pdf)) deallocate (res%forces_pdf)
   if (allocated(res%forces_sf)) deallocate (res%forces_sf)
   if (allocated(res%forces_xrd)) deallocate (res%forces_xrd)
   if (allocated(res%forces_nd)) deallocate (res%forces_nd)

   if (allocated(res%this_forces)) deallocate (res%this_forces)
   if (allocated(res%this_forces_vdw)) deallocate (res%this_forces_vdw)
   if (allocated(res%this_forces_lp)) deallocate (res%this_forces_lp)
   if (allocated(res%this_forces_pdf)) deallocate (res%this_forces_pdf)
   if (allocated(res%this_forces_sf)) deallocate (res%this_forces_sf)
   if (allocated(res%this_forces_xrd)) deallocate (res%this_forces_xrd)
   if (allocated(res%this_forces_nd)) deallocate (res%this_forces_nd)

   if (allocated(res%local_properties)) deallocate (res%local_properties)
   if (allocated(res%local_properties_cart_der)) deallocate (res%local_properties_cart_der)
   if (allocated(res%this_local_properties)) deallocate (res%this_local_properties)
   if (allocated(res%this_local_properties_cart_der)) deallocate (res%this_local_properties_cart_der)

   if (allocated(model%soap_turbo_hypers)) deallocate (model%soap_turbo_hypers)
   if (allocated(model%distance_2b_hypers)) deallocate (model%distance_2b_hypers)
   if (allocated(model%angle_3b_hypers)) deallocate (model%angle_3b_hypers)
   if (allocated(model%core_pot_hypers)) deallocate (model%core_pot_hypers)

   deallocate (dom%n_atom_pairs_by_rank)
   if (allocated(model%n_local_properties_mpi)) deallocate (model%n_local_properties_mpi)
   if (allocated(model%local_properties_n_sparse_mpi_soap_turbo)) deallocate (model%local_properties_n_sparse_mpi_soap_turbo)
   if (allocated(model%local_properties_dim_mpi_soap_turbo)) deallocate (model%local_properties_dim_mpi_soap_turbo)
   if (allocated(model%has_local_properties_mpi)) deallocate (model%has_local_properties_mpi)

   if (allocated(model%local_property_labels)) deallocate (model%local_property_labels)
   if (allocated(model%local_property_indexes)) deallocate (model%local_property_indexes)
   if (allocated(dom%do_list)) deallocate (dom%do_list)
   if (allocated(params%write_local_properties)) deallocate (params%write_local_properties)

   if (params%vdw_type == "ts+mbd") then
      if (rank == 0) then
         open (unit=30, file="mbd_ts_scaling.dat", status="unknown")
         do i = 1, state%n_sites
#ifdef _MPIF90
            write (30, *) res%this_mbd_ts_scaling(i)
#else
            write (30, *) res%mbd_ts_scaling(i)
#endif
         end do
         close (30)
      end if
   end if

   if (rank == 0) then
      write (*, *) '                                       |'
      write (*, *) 'End of execution                       |'
      write (*, *) '_______________________________________/'
   end if

#ifdef _GPU
!  The high-water mark, which is the number that sizes the next run.
!
!  Before gpu_context_finalize, which calls hipDeviceReset and takes the whole
!  context down -- after it there is nothing left to ask. Printed unconditionally
!  and to stderr: it costs one line, and "what did that actually use" is the
!  first question asked after any run that was close to the limit.
   if (rank == 0) call gpu_memory_report("end of run")
#endif
   call comm_finalize(comm)

   call gpu_context_finalize(params, n_omp)

end program turbogap
