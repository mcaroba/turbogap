! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_evaluate.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  One evaluation of the potential on the current structure: every energy,
!  force and virial term in the order they are summed, the per-site properties
!  on the way, the IR bias, and what a single-point prediction writes.
module turbogap_evaluate

   use kinds, only: dp
   use types, only: input_parameters, perform_t, any_has_vdw, any_has_core_electron_be
   use timing, only: times_t, time_start, time_end
   use turbogap_comm, only: comm_t
   use turbogap_structure, only: state_t
   use turbogap_domain, only: domain_t, neighbors_t, domain_complete_e0, domain_complete_sites, domain_total_dipole, &
                              domain_n_total, domain_decomposed, domain_site_values, domain_site_ids, domain_exp_view, &
                              domain_complete_contributions, exp_on_domains
   use turbogap_setup, only: model_t
   use turbogap_results, only: results_t, results_prepare
   use turbogap_loop, only: loop_t
   use turbogap_md, only: dynamics_t
   use turbogap_ir, only: ir_run_t, ir_before_evaluate, ir_push_frame, ir_after_forces
   use turbogap_soap, only: compute_soap
   use turbogap_estat, only: compute_estat
   use turbogap_vdw, only: vdw_state, compute_vdw, vdw_read_ts_scaling
   use turbogap_exp, only: compute_exp_xps, compute_exp_spectra
   use turbogap_output, only: print_single_point_energies, write_debug_forces, write_single_point
   use gap_backend, only: gap_backend_begin, gap_backend_end, add_2b_contribution, &
                          add_core_pot_contribution, add_3b_contribution
   use exp_utils, only: calculate_exp_interpolation, build_exp_weights
   use exp_interface, only: get_write_condition, get_overwrite_condition, preprocess_exp_data
   use read_files, only: write_exp_data
   use soap_turbo_functions, only: cross_product
   use neighbors, only: cut_neighbors_list

   implicit none

   private
   public :: evaluate

!  Round-off of adding neighbors_buffer to the cutoff and taking it off again.
   real(dp), parameter :: RCUT_ROUNDOFF = 1.d-8

contains

   subroutine evaluate(res, state, nl, dom, model, params, perform, loop, dyn, ir, vdw_ws, comm, &
                       do_electrostatics, n_omp, time3, time)
      type(results_t), target, intent(inout) :: res
      type(state_t), intent(inout) :: state
      type(neighbors_t), intent(inout), target :: nl
      type(domain_t), intent(inout) :: dom
      type(model_t), target, intent(inout) :: model
      type(input_parameters), intent(inout) :: params
      type(perform_t), intent(in) :: perform
      type(loop_t), intent(inout) :: loop
      type(dynamics_t), intent(inout) :: dyn
      type(ir_run_t), intent(inout) :: ir
      type(vdw_state), intent(inout) :: vdw_ws
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: do_electrostatics
      integer, intent(inout) :: n_omp
      real(dp), intent(in) :: time3
      type(times_t), intent(inout) :: time
!     Output scheduling that can carry from one call into the next, as it did
!     when these were variables of the main program.
      logical, save :: write_condition = .false.
      logical, save :: overwrite_condition = .false.
      character*32, save :: exp_output = "none"
      logical :: resized
      integer, allocatable :: i_beg_list(:)
      integer, allocatable :: i_end_list(:)
      integer, allocatable :: j_beg_list(:)
      integer, allocatable :: j_end_list(:)
      integer :: this_i_beg
      integer :: this_i_end
      integer :: this_j_beg
      integer :: this_j_end
      integer :: ierr
      integer :: i
      integer :: j
!     The list the GAP terms walk: nl, or the pairs of it within their
!     cutoff when a pdf, XRD, ND or vdW term made it reach further. gap_pair is
!     then each pair's index in nl, and gj_end the last pair. A pair term with
!     a shorter cutoff walks tnl, cut from it the same way.
      type(neighbors_t), target :: gap_nl
      type(neighbors_t), pointer :: gnl
      integer, allocatable :: gap_pair(:)
      integer :: gj_end
      real(dp) :: rcut_gap
      type(neighbors_t), target :: term_nl
      type(neighbors_t), pointer :: tnl
      integer, allocatable :: term_pair(:)
      integer :: tj_end
      character*1024 :: filename
      character*1024 :: temp_string
      real(dp), allocatable :: be_all(:)
      integer, allocatable :: site_ids(:)
#ifdef _GPU
      integer :: omp_task
      integer :: n_pairs_temp
      integer :: n_sites_temp
      character*8, allocatable, target :: species_types_actual(:)
#endif

      call results_prepare(res, state, params, model, perform, loop, &
                           dom%n_atom_pairs_by_rank(comm%rank + 1), dom%n_atom_pairs_by_rank_prev, resized)
      if (resized) call vdw_read_ts_scaling(res, state, params, loop, comm)
      if (params%do_forces) call ir_before_evaluate(ir, params, state, res, loop, comm, dom)

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

      rcut_gap = gap_cutoff()
      if (model%rcut_max > rcut_gap + params%neighbors_buffer + RCUT_ROUNDOFF .and. &
          .not. params%write_derivatives) then
         call cut_neighbors_list(rcut_gap, dom%i_beg, dom%i_end, nl%n_neigh, nl%rjs, nl%thetas, nl%phis, nl%xyz, &
                                 nl%neighbor_species, nl%neighbors_list, gap_nl%n_neigh, gap_nl%rjs, gap_nl%thetas, &
                                 gap_nl%phis, gap_nl%xyz, gap_nl%neighbor_species, gap_nl%neighbors_list, gap_pair)
         gnl => gap_nl
         gj_end = sum(gap_nl%n_neigh(dom%i_beg:dom%i_end))
      else
         gnl => nl
         gj_end = dom%j_end
         allocate (gap_pair(0))
      end if

      call compute_soap(res, state, gnl, gj_end, gap_pair, dom, model, params, loop, comm, time)

      call domain_complete_sites(dom, comm, res, state, params, model, time)

      if (params%do_dipole) call domain_total_dipole(dom, comm, res, state)

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
      call compute_estat(params, do_electrostatics, model%valid_estat_charges, model%charge_lp_index, &
                         state%n_sites, nl%n_neigh, nl%neighbors_list, state%species, nl%neighbor_species, nl%rjs, nl%xyz, &
                         res%local_properties, res%local_properties_cart_der, &
                         dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, comm, n_omp, &
                         res%this_energies_estat, res%this_forces_estat, res%this_virial_estat, time, &
                         domain_n_total(dom, state), domain_decomposed(dom))

      call compute_vdw(params, any_has_vdw(model%soap_turbo_hypers), state%n_sites, &
                       nl%n_neigh, nl%neighbors_list, nl%neighbor_species, nl%rjs, nl%xyz, &
                       res%local_properties, res%local_properties_cart_der, model%vdw_lp_index, &
                       dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, dom%n_atom_pairs_by_rank, dom%site_in_rank, &
                       state%indices, comm%rank, comm%size, loop%md_istep, vdw_ws, &
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
                    & domain_n_total(dom, state), dot_product(cross_product(state%a_box,&
                    & state%b_box), state%c_box)/(dfloat(state%indices(1)*state%indices(2) &
                    &*state%indices(3))), params%exp_data(i)%input, exp_output, .true.)

               if (params%write_exp .and. .not. params&
                    &%exp_data(i)%wrote_exp .and. comm%rank == 0 .and. write_condition) then

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
                 &%exp_data(i)%wrote_exp .and. comm%rank == 0 .and. write_condition) then

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
      if (any_has_core_electron_be(model%soap_turbo_hypers) .and. params%do_prediction .and. &
          (model%valid_xps .or. params%do_xps)) then
         call domain_site_values(dom, comm, state, res%local_properties(:, model%core_be_lp_index), be_all)
      else
         allocate (be_all(0))
      end if
      call domain_site_ids(dom, state, site_ids)
      call compute_exp_xps(params, state%n_sites, loop%n_xyz, nl%xyz, nl%neighbors_list, nl%n_neigh, &
                           res%local_properties, res%local_properties_cart_der, model%soap_turbo_hypers, &
                           state%a_box, state%b_box, state%c_box, state%indices, dom%i_beg, dom%i_end, dom%j_beg, &
                           dom%j_end, comm, &
                           loop%md_istep, loop%mc_istep, model%valid_xps, model%xps_idx, model%core_be_lp_index, &
                           write_condition, overwrite_condition, exp_output, &
                           res%this_energies_lp, res%this_forces_lp, res%this_virial_lp, time, be_all, site_ids)

      !###---   (Partial) Pair distribution functions and XRD   ---###!
      !
      ! Partial sums into the this_ arrays; exp_interface allocates them.
!     Under a domain decomposition these are computed on the whole structure,
!     which they are a property of; the potential stays decomposed.
      if (domain_decomposed(dom) .and. (params%do_pair_distribution .or. params%do_structure_factor .or. &
                                        params%do_xrd .or. params%do_nd) .and. .not. exp_on_domains(params)) then
         call domain_exp_view(dom, comm, state, params, model)
         call exp_spectra(dom%view%state, dom%view%nl, dom%view%i_beg, dom%view%i_end, dom%view%j_beg, &
                          dom%view%j_end, dom%view%state%n_sites)
      else
         call exp_spectra(state, nl, dom%i_beg, dom%i_end, dom%j_beg, dom%j_end, domain_n_total(dom, state))
      end if

      if (params%do_prediction) then
         !       Two-body, core-potential and three-body contributions, via the
         !       gap_backend seam. The CPU implementation is in
         !       src/gap_backend_cpu.f90; the GPU branch provides the same three
         !       names from src/gap_backend_gpu.f90 and the Makefile picks one.
         call time_start(time%gap)

         if (model%n_distance_2b > 0) then
            call term_list(twob_cutoff())
            call gap_backend_begin(params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                   tnl%neighbors_list, dom%i_beg, dom%i_end, dom%j_beg, tj_end)
            call add_2b_contribution(model%n_distance_2b, model%distance_2b_hypers, &
                                     params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, tj_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%energies_2b, res%forces_2b, res%virial_2b, time)
            call gap_backend_end()
         end if

         if (model%n_core_pot > 0) then
            call term_list(core_pot_cutoff())
            call gap_backend_begin(params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                   tnl%neighbors_list, dom%i_beg, dom%i_end, dom%j_beg, tj_end)
            call add_core_pot_contribution(model%n_core_pot, model%core_pot_hypers, &
                                           params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                           dom%i_beg, dom%i_end, dom%j_beg, tj_end, res%this_energies, res%this_forces, &
                                           res%this_virial, &
                                           res%energies_core_pot, res%forces_core_pot, res%virial_core_pot, time)
            call gap_backend_end()
         end if

         if (model%n_angle_3b > 0) then
            call term_list(angle_cutoff())
            call gap_backend_begin(params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                   tnl%neighbors_list, dom%i_beg, dom%i_end, dom%j_beg, tj_end)
            call add_3b_contribution(model%n_angle_3b, model%angle_3b_hypers, tnl%neighbors_list, &
                                     params, tnl%rjs, tnl%xyz, tnl%n_neigh, state%species, tnl%neighbor_species, &
                                     dom%i_beg, dom%i_end, dom%j_beg, tj_end, res%this_energies, res%this_forces, &
                                     res%this_virial, &
                                     res%forces, res%energies_3b, res%forces_3b, res%virial_3b, time)
            call gap_backend_end()
         end if

         call time_end(time%gap)
         !       Communicate all energies and forces here for all
         !       terms
         call domain_complete_contributions(dom, comm, res, state, params, model, perform, time)

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
!        Replicated data sums the families before the reduce, in this order.
         if (.not. res%forces_summed) then
            res%forces = res%forces_soap + res%forces_2b + res%forces_3b + res%forces_core_pot + res%forces_vdw
         end if
         res%virial = res%virial_soap + res%virial_2b + res%virial_3b + res%virial_core_pot + res%virial_vdw

!           MAD IR bias. The dipole of this configuration joins the ensemble,
!           the spectrum is compared with the experiment, and the gradient of
!           the mismatch with respect to THIS configuration is added to the
!           forces. Nothing is applied until the ensemble is full, because a
!           partly filled one has a resolution that changes step to step.
!
!           No virial: the bias is a function of the dipole, not of the cell,
!           and a stress from it would be wrong rather than merely missing.
         call ir_after_forces(ir, params, state, res, loop, comm, dyn%md_time, time3, write_condition, time, dom)

         if (.not. res%forces_summed) then
            if (model%valid_estat_charges) res%forces = res%forces + res%forces_estat
            if (perform%xps_forces) res%forces = res%forces + res%forces_lp
            if (perform%pdf_forces) res%forces = res%forces + res%forces_pdf
            if (perform%sf_forces) res%forces = res%forces + res%forces_sf
            if (perform%xrd_forces) res%forces = res%forces + res%forces_xrd
            if (perform%nd_forces) res%forces = res%forces + res%forces_nd
         end if

         if (model%valid_estat_charges) res%virial = res%virial + res%virial_estat
         if (perform%xps_forces) res%virial = res%virial + res%virial_lp
         if (perform%pdf_forces) res%virial = res%virial + res%virial_pdf
         if (perform%sf_forces) res%virial = res%virial + res%virial_sf
         if (perform%xrd_forces) res%virial = res%virial + res%virial_xrd
         if (perform%nd_forces) res%virial = res%virial + res%virial_nd

         call write_debug_forces(comm, res, state, dom, params, model, perform, loop)

      end if

      if (params%do_prediction .and. .not. params%do_md .and. .not. params%do_mc) then
         call write_single_point(comm, res, state, dyn, params, model, loop)
      end if
   contains

!     The pair distribution, structure factor, XRD and ND, on the structure and
!     neighbour list given: this rank's own, or the whole-structure view.
!     tnl: gnl, or gnl cut at rcut when that is shorter than the GAP terms'
!     reach. tj_end is its last pair.
      subroutine term_list(rcut)
         real(dp), intent(in) :: rcut
         if (rcut < rcut_gap - RCUT_ROUNDOFF) then
            call cut_neighbors_list(rcut, dom%i_beg, dom%i_end, gnl%n_neigh, gnl%rjs, gnl%thetas, gnl%phis, gnl%xyz, &
                                    gnl%neighbor_species, gnl%neighbors_list, term_nl%n_neigh, term_nl%rjs, &
                                    term_nl%thetas, term_nl%phis, term_nl%xyz, term_nl%neighbor_species, &
                                    term_nl%neighbors_list, term_pair)
            tnl => term_nl
            tj_end = sum(term_nl%n_neigh(dom%i_beg:dom%i_end))
         else
            tnl => gnl
            tj_end = gj_end
         end if
      end subroutine term_list

      real(dp) function twob_cutoff()
         integer :: k
         twob_cutoff = 0.d0
         do k = 1, model%n_distance_2b
            twob_cutoff = max(twob_cutoff, model%distance_2b_hypers(k)%rcut)
         end do
      end function twob_cutoff

      real(dp) function core_pot_cutoff()
         integer :: k
         core_pot_cutoff = 0.d0
         do k = 1, model%n_core_pot
            core_pot_cutoff = max(core_pot_cutoff, maxval(model%core_pot_hypers(k)%x))
         end do
      end function core_pot_cutoff

      real(dp) function angle_cutoff()
         integer :: k
         angle_cutoff = 0.d0
         do k = 1, model%n_angle_3b
            angle_cutoff = max(angle_cutoff, model%angle_3b_hypers(k)%rcut)
         end do
      end function angle_cutoff

!     The furthest any GAP term reaches: SOAP, 2b, 3b, core potential.
      real(dp) function gap_cutoff()
         integer :: k
         gap_cutoff = 0.d0
         do k = 1, model%n_soap_turbo
            gap_cutoff = max(gap_cutoff, model%soap_turbo_hypers(k)%rcut_max)
         end do
         do k = 1, model%n_distance_2b
            gap_cutoff = max(gap_cutoff, model%distance_2b_hypers(k)%rcut)
         end do
         do k = 1, model%n_angle_3b
            gap_cutoff = max(gap_cutoff, model%angle_3b_hypers(k)%rcut)
         end do
         do k = 1, model%n_core_pot
            gap_cutoff = max(gap_cutoff, maxval(model%core_pot_hypers(k)%x))
         end do
      end function gap_cutoff

      subroutine exp_spectra(es, enl, ib, ie, jb, je, n_tot)
         type(state_t), intent(inout) :: es
         type(neighbors_t), intent(in) :: enl
         integer, intent(in) :: ib
         integer, intent(in) :: ie
         integer, intent(in) :: jb
         integer, intent(in) :: je
!        The atoms in the system; es holds only a domain's when decomposed.
         integer, intent(in) :: n_tot

#ifdef _GPU
         call compute_exp_spectra(params, es%n_sites, es%species, es%positions, enl%rjs, enl%xyz, enl%neighbors_list, &
                                  enl%n_neigh, enl%neighbor_species, es%indices, es%a_box, es%b_box, es%c_box, &
                                  ib, ie, jb, &
                                  je, comm%rank, comm%size, ierr, loop%md_istep, loop%mc_istep, res%this_energies_pdf, &
                                  res%this_forces_pdf, res%this_virial_pdf, res%this_energies_sf, &
                                  res%this_forces_sf, res%this_virial_sf, res%this_energies_xrd, res%this_forces_xrd, &
                                  res%this_virial_xrd, res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, time, &
                                  i_beg_list, i_end_list, j_beg_list, &
                                  j_end_list, n_omp, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, &
                                  n_sites_temp, n_pairs_temp, write_condition, overwrite_condition, &
                                  temp_string, species_types_actual, es%v_uc, n_tot)
#else
         call compute_exp_spectra(params, es%n_sites, es%species, es%positions, enl%rjs, enl%xyz, enl%neighbors_list, &
                                  enl%n_neigh, enl%neighbor_species, es%indices, es%a_box, es%b_box, es%c_box, &
                                  ib, ie, jb, je, comm%rank, comm%size, ierr, loop%md_istep, loop%mc_istep, &
                                  res%this_energies_pdf, res%this_forces_pdf, res%this_virial_pdf, &
                                  res%this_energies_sf, res%this_forces_sf, res%this_virial_sf, &
                                  res%this_energies_xrd, res%this_forces_xrd, res%this_virial_xrd, &
                                  res%this_energies_nd, res%this_forces_nd, res%this_virial_nd, &
                                  time, n_tot)
#endif
      end subroutine exp_spectra

   end subroutine evaluate

end module turbogap_evaluate
