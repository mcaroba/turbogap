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
   use turbogap_structure, only: state_t, structure_acquire, structure_free
   use turbogap_domain, only: domain_t, neighbors_t, domain_sync_state, domain_build, &
                              domain_complete_sites, domain_complete_e0, &
                              domain_complete_contributions, domain_sync_after_md, &
                              domain_sync_after_ipi, domain_free
   use turbogap_results, only: results_t, results_prepare, results_free
   use turbogap_soap, only: compute_soap, soap_free_device
   use turbogap_sampling, only: sampling_t, mc_prepare_step, nested_step, mc_step
   use turbogap_ir, only: ir_run_t, ir_init, ir_step_begin, ir_before_evaluate, ir_push_frame, &
                          ir_after_forces, ir_step_end, ir_finish, ir_report
   use turbogap_loop, only: loop_t, loop_init, loop_continues, loop_begin_step, creturn
   use turbogap_output, only: print_banner, print_options, print_single_point_energies, &
                              write_debug_forces, write_single_point, print_nothing_to_do, &
                              print_timing_report
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

   type(perform_t) :: perform
   integer :: n_omp = 1

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
            call print_single_point_energies(comm, res, params, model, perform)
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

            call write_debug_forces(comm, res, state, dom, params, model, perform, loop)

         end if

         if (params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc) then
            call write_single_point(comm, res, state, dyn, params, model, loop)
         end if
      else
         call print_nothing_to_do(comm)
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

      call mc_step(smp, state, res, dyn, nl, model, params, perform, loop, comm, time)

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

   call print_timing_report(comm, params, model, loop, ir, do_electrostatics, time3, time)

   call soap_free_device(model)
   call structure_free(state)
   call results_free(res)
   call model_free(model)
   call domain_free(dom)
   if (allocated(params%write_local_properties)) deallocate (params%write_local_properties)

   call vdw_write_ts_scaling(res, state, params, comm)

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
