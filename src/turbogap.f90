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
   use turbogap_structure, only: state_t, structure_acquire, structure_update_volume, structure_free
   use turbogap_domain, only: domain_t, neighbors_t, domain_sync_state, domain_build, &
                              domain_complete_sites, domain_complete_e0, &
                              domain_complete_contributions, domain_sync_after_md, &
                              domain_sync_after_ipi, domain_free, domain_end_step
   use turbogap_results, only: results_t, results_prepare, results_free
   use turbogap_soap, only: soap_free_device
   use turbogap_evaluate, only: evaluate
   use turbogap_sampling, only: sampling_t, sampling_init, mc_prepare_step, nested_step, mc_step
   use turbogap_ir, only: ir_run_t, ir_init, ir_step_begin, ir_before_evaluate, ir_push_frame, &
                          ir_after_forces, ir_step_end, ir_finish, ir_report
   use turbogap_loop, only: loop_t, loop_init, loop_continues, loop_begin_step, loop_sync, &
                            loop_end_step, creturn
   use turbogap_output, only: handle_help_request, read_run_mode, print_end, print_banner, print_options, &
                              print_single_point_energies, &
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
   type(loop_t) :: loop
   type(dynamics_t) :: dyn
   type(ir_run_t) :: ir
   type(sampling_t) :: smp
   integer :: n_pos

   type(perform_t) :: perform
   integer :: n_omp = 1

   character*1024 :: filename
   character*1024 :: string
   character*1024 :: temp_string
   character*1024 :: temp_string2

   ! This is the mode in which we run TurboGAP
   character*16 :: mode = "none"

   ! Here we store the input parameters
   type(input_parameters) :: params

   logical :: do_electrostatics = .true.
! Persistent ts+mbd correction state, owned by turbogap_vdw
   type(vdw_state) :: vdw_ws

   character*32 :: implemented_exp_observables(1:5)
#ifdef _GPU
#endif

   call handle_help_request()

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

   call read_run_mode(mode)

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
   call sampling_init(smp)

   call loop_init(loop, params, comm)

   ! This checks if we need to do the SOAP calculation more than once, if there are several concatenated
   ! structures in the xyz file provided or we're doing molecular dynamics

   call exp_decide(perform, params, model)

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
      call structure_update_volume(state)

      !   If we are doing prediction, we run this chunk of code
      if (params%do_prediction .or. params%write_soap .or. params%write_derivatives) then
         call evaluate(res, state, nl, dom, model, params, perform, loop, dyn, ir, vdw_ws, comm, &
                       do_electrostatics, n_omp, time3, time)
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
         call md_step(dyn, state, res, nl, model, params, loop, smp%i_image, smp%i_nested, comm, time)
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

      call loop_sync(loop, params, comm, time)
      call domain_sync_state(dom, comm, state, params, time)
      call domain_end_step(dom, nl, comm, state, params, loop)

      call exp_end_run(params, loop)

      call ir_step_end(ir, params, loop)

      call loop_end_step(loop, state%n_sites, params, comm)
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

   call print_end(comm)

!  The high-water mark, which is the number that sizes the next run.
!
!  Before gpu_context_finalize, which calls hipDeviceReset and takes the whole
!  context down -- after it there is nothing left to ask. Printed unconditionally
!  and to stderr: it costs one line, and "what did that actually use" is the
!  first question asked after any run that was close to the limit.
   if (rank == 0) call gpu_memory_report("end of run")
   call comm_finalize(comm)

   call gpu_context_finalize(params, n_omp)

end program turbogap
