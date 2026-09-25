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

   use kinds, only: dp
   use timing, only: times_t, get_time, time_start, time_end
   use types, only: input_parameters, perform_t
   use turbogap_comm, only: comm_t, comm_init, comm_finalize
   use threads, only: threads_init
   use gpu_context, only: gpu_context_init, gpu_context_finalize, gpu_memory_budget_init, gpu_memory_report
   use gap_backend, only: gap_backend_init
   use turbogap_setup, only: model_t, read_input_and_gap_files, model_free
   use turbogap_output, only: handle_help_request, read_run_mode, print_banner, print_options, &
                              print_nothing_to_do, print_timing_report, print_end
   use turbogap_loop, only: loop_t, loop_init, loop_continues, loop_begin_step, loop_sync, loop_end_step
   use turbogap_structure, only: state_t, structure_acquire, structure_update_volume, structure_free
   use turbogap_domain, only: domain_t, neighbors_t, domain_init, domain_sync_state, domain_build, &
                              domain_sync_after_md, domain_sync_after_ipi, domain_end_step, domain_free
   use turbogap_results, only: results_t, results_free
   use turbogap_evaluate, only: evaluate
   use turbogap_md, only: dynamics_t, md_prepare_velocities, md_step
   use ipi_driver, only: ipi_driver_open, ipi_driver_exchange, ipi_driver_close
   use turbogap_sampling, only: sampling_t, sampling_init, mc_prepare_step, nested_step, mc_step
   use turbogap_exp, only: exp_decide, exp_end_run
   use turbogap_ir, only: ir_run_t, ir_init, ir_step_begin, ir_step_end, ir_finish
   use turbogap_vdw, only: vdw_state, vdw_write_ts_scaling
   use turbogap_soap, only: soap_free_device

   implicit none

   type(comm_t) :: comm
   type(input_parameters) :: params
   type(model_t), target :: model
   type(state_t) :: state
   type(domain_t) :: dom
   type(neighbors_t) :: nl
   type(results_t), target :: res
   type(loop_t) :: loop
   type(dynamics_t) :: dyn
   type(ir_run_t) :: ir
   type(sampling_t) :: smp
   type(perform_t) :: perform
   type(vdw_state) :: vdw_ws
   type(times_t) :: time
   real(dp) :: time1
   real(dp) :: time3
   integer :: n_omp = 1
   character*16 :: mode = "none"
   logical :: do_electrostatics = .true.

   call handle_help_request()

!  Device streams and handles; empty in the host build.
   call time_start(time%create_streams)
!  Before comm_init, so every rank passes 0: slurm gives each rank one device.
   call gpu_context_init(params, comm%rank, n_omp)
   call gap_backend_init()
   call time_end(time%create_streams)

   call get_time(time1)
   time3 = time1
   call time_start(time%setup)

   call comm_init(comm)
!  Needs the communicator: how many threads a rank may take depends on how many
!  ranks share its node.
   call threads_init()
   call domain_init(dom, comm)

   call read_run_mode(mode)

   call print_banner(comm)

   call read_input_and_gap_files(mode, comm%rank, comm%size, params, &
                                 model%soap_turbo_hypers, model%distance_2b_hypers, model%angle_3b_hypers, model%core_pot_hypers, &
                                 model%n_soap_turbo, model%n_distance_2b, model%n_angle_3b, model%n_core_pot, model%n_species, &
                                 model%rcut_max, &
                                 model%valid_xps, model%xps_idx, model%vdw_lp_index, model%core_be_lp_index, &
                                 model%valid_estat_charges, model%charge_lp_index, &
                                 model%local_property_labels, model%local_property_indexes, model%n_local_properties_mpi, &
                                 model%has_local_properties_mpi, model%local_properties_n_sparse_mpi_soap_turbo, &
                                 model%local_properties_dim_mpi_soap_turbo, dyn%nrows, dyn%allelstopdata, &
                                 dyn%ephbeta, dyn%ephfdm, dyn%ephlsc, time)

!  After the input: it reads mem_fraction and writes max_Gbytes_per_process.
   call gpu_memory_budget_init(params, comm%rank, comm%size)

   call print_options(comm, params, model)

   model%xps_idx = params%xps_idx
   call sampling_init(smp)

   call loop_init(loop, params, comm)

   call exp_decide(perform, params, model)

   call ir_init(ir, params, comm)

   call time_end(time%setup)

!  Connect now, so a missing i-PI server is reported before the first force call.
   if (mode == "ipi") call ipi_driver_open(params%ipi_address, comm%rank)

   do while (loop_continues(loop, params))
      call ir_step_begin(ir, params)

      call loop_begin_step(loop, params, comm)

      call structure_acquire(state, nl%rebuild_neighbors_list, loop, params, model, comm, smp%mc_file, time)

      call md_prepare_velocities(dyn, state, params, loop, comm)
      call mc_prepare_step(smp, dyn, state, params, loop, comm)
      call domain_sync_state(dom, comm, state, params, time)
      call domain_build(dom, nl, comm, state, params, model, loop, smp%mc_file, time)
      call structure_update_volume(state)

      if (params%do_prediction .or. params%write_soap .or. params%write_derivatives) then
         call evaluate(res, state, nl, dom, model, params, perform, loop, dyn, ir, vdw_ws, comm, &
                       do_electrostatics, n_omp, time3, time)
      else
         call print_nothing_to_do(comm)
      end if

!     Under i-PI the forces go out over the socket and the next positions come back.
      if (mode == "ipi") then
         call ipi_driver_exchange(comm%rank, state%n_sites, state%positions, state%positions_prev, state%positions_diff, &
                                  state%velocities, state%a_box, state%b_box, state%c_box, state%indices, params%neighbors_buffer, &
                                  res%forces, res%energy, res%virial, loop%exit_loop, nl%rebuild_neighbors_list)
         call domain_sync_after_ipi(dom, comm, state, nl, loop)
      else
         call md_step(dyn, state, res, nl, model, params, loop, smp%i_image, smp%i_nested, comm, time)
         call domain_sync_after_md(dom, comm, state, nl, params, time)
      end if

      call nested_step(smp, state, res, dyn, nl, params, loop, comm)

      call mc_step(smp, state, res, dyn, nl, model, params, perform, loop, comm, time)

      call loop_sync(loop, params, comm, time)
      call domain_sync_state(dom, comm, state, params, time)
      call domain_end_step(dom, nl, comm, state, params, loop)

      call exp_end_run(params, loop)

      call ir_step_end(ir, params, loop)

      call loop_end_step(loop, state%n_sites, params, comm)
      if (loop%exit_loop) exit
   end do

!  Close before the reports, so i-PI sees a clean exit rather than a timeout.
   if (mode == "ipi") call ipi_driver_close(comm%rank)

   call ir_finish(ir, params, comm, time)

   call print_timing_report(comm, params, model, loop, ir, do_electrostatics, time3, time)

   call vdw_write_ts_scaling(res, state, params, comm)

   call soap_free_device(model)
   call structure_free(state)
   call results_free(res)
   call model_free(model)
   call domain_free(dom)
   if (allocated(params%write_local_properties)) deallocate (params%write_local_properties)

   call print_end(comm)

!  Before gpu_context_finalize, which resets the device.
   if (comm%rank == 0) call gpu_memory_report("end of run")
   call comm_finalize(comm)

   call gpu_context_finalize(params, n_omp)

end program turbogap
