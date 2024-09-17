! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, electrostatics.f90, is copyright (c) 2023, Max Veit
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

module electrostatics

    use neighbors
    use types


    ! Both of these from the NIST website, references 2018 CODATA values
    real(dp), parameter :: HARTREE_EV = 27.2113862460_dp
    real(dp), parameter :: BOHR_ANG = 0.5291772109_dp
    ! TODO this needs to be checked with the turbogap unit system
    ! it will probably just be Hartree * Bohr, since the charges are
    ! expressed in atomic units but the lengths in angstroms
    real(dp), parameter :: COUL_CONSTANT = HARTREE_EV * BOHR_ANG ! = q_e^2 / 4πε_0
    real(dp), parameter :: FLOAT_ZERO = 10*EPSILON(1.0_dp)
    real(dp), parameter :: PI = dacos(-1.0_dp)
    real(dp), parameter :: TWO_OVER_SQRT_PI = 2.0_dp / sqrt(PI)

    contains

    ! Outer product for fixed-size vectors. Doesn't work for allocatables
    ! TODO move this to some math utils file
    function outer_prod(vec_i, vec_j)

        implicit none
        real(dp), dimension(:), intent(in) :: vec_i
        real(dp), dimension(:), intent(in) :: vec_j
        real(dp), dimension(size(vec_i), size(vec_j)) :: outer_prod
        integer :: ii, jj

        do ii = 1, size(vec_i)
            do jj = 1, size(vec_j)
                outer_prod(ii, jj) = vec_i(ii) * vec_j(jj)
            end do
        end do
    end function

    ! This is just the Coulomb 1/r form without any prefactors
    function pair_energy_direct(rij)
        implicit none
        real(dp), intent(in) :: rij
        real(dp) :: pair_energy_direct

        pair_energy_direct = 1.0_dp / rij
    end function

    ! Derivative of the above
    ! Multiply by the r_ij _unit vector_ to get the force
    ! (so divide the full vector by rij)
    function der_pair_energy_direct(rij)
        implicit none
        real(dp), intent(in) :: rij
        real(dp) :: der_pair_energy_direct

        der_pair_energy_direct = -1.0_dp / rij**2
    end function

    ! Ok, this is only the Wolf pair energy.  DSF itself does some shifting
    ! on top of this, but it's more convenient to keep just this function separate
    function pair_energy_dsf(rij, alpha)
        implicit none
        real(dp), intent(in) :: rij, alpha
        real(dp) :: pair_energy_dsf

        pair_energy_dsf = erfc(alpha*rij) / rij
    end function

    ! Again - Wolf, not DSF, to make shifting easier later on.
    function der_pair_energy_dsf(rij, alpha)
        implicit none
        real(dp), intent(in) :: rij, alpha
        real(dp) :: der_pair_energy_dsf

        der_pair_energy_dsf = -1. * erfc(alpha*rij) / rij**2 - &
            TWO_OVER_SQRT_PI * alpha * exp(-1. * alpha**2 * rij**2) / rij
    end function

    ! Compute electrostatic energies, forces, and virials via a direct
    ! summation of the Coulomb law.  This should _not_ be used for any
    ! serious applications, since it dos not converge with increasing
    ! cutoff!
    subroutine compute_coulomb_direct(charges, charge_gradients, &
                                      n_neigh, neighbors_list, &
                                      rcut, rcut_in, rcin_width, rjs, xyz, &
                                      neighbor_charges, do_gradients, &
                                      local_energies, forces, virial)
        implicit none
        real(dp), dimension(:), intent(in) :: charges, neighbor_charges, rjs
        real(dp), dimension(:,:), intent(in) :: charge_gradients, xyz
        integer, dimension(:), intent(in) :: n_neigh, neighbors_list
        real(dp), intent(in) :: rcut, rcut_in, rcin_width
        logical, intent(in) :: do_gradients

        ! inout because they are initialized outside this procedure and filled with zeros
        real(dp), intent(inout), dimension(:) :: local_energies
        real(dp), intent(inout), dimension(:,:) :: forces
        real(dp), intent(inout), dimension(3,3) :: virial

        integer :: center_i, neigh_id, soap_neigh_id, neigh_seq, soap_neigh_seq
        integer :: n_sites_global, n_sites_this, pair_counter, soap_pair_counter
        real(dp) :: rij, center_term, pair_energy, neigh_charge, inner_damp_ij
        real(dp), dimension(3) :: rij_vec, fij_vec, fki_vec
        real(dp), dimension(:), allocatable :: vc_grad_prefactor

        n_sites_this = size(n_neigh)
        n_sites_global = size(forces, 2)
        allocate(vc_grad_prefactor(1:n_sites_this))
        vc_grad_prefactor = 0.0_dp

        pair_counter = 0
        soap_pair_counter = 0
        do center_i = 1, n_sites_this
            pair_counter = pair_counter + 1
            !soap_pair_counter = soap_pair_counter + 1 ! No, because we include the center as a SOAP neighbour
            ! First we precompute q_i/4πε_0
            ! TODO this is where we add an effective dielectric constant to scale the interaction
            center_term = charges(center_i) * COUL_CONSTANT
            do neigh_seq = 2, n_neigh(center_i)
                pair_counter = pair_counter + 1 ! ???
                ! so... the neighbourlist is indexed by pair ID. This needs to be
                ! counted and updated MANUALLY (ughhhh) while iterating through pairs.
                ! (there are so, so many ways that this could go wrong....)
                ! Oh, and the neighbour ID can exceed the number of atoms...?
                ! Maybe these are periodic replicas (see below).
                neigh_id = modulo(neighbors_list(pair_counter) - 1, n_sites_global) + 1
                rij = rjs(pair_counter)
                rij_vec = xyz(:, pair_counter)
                neigh_charge = neighbor_charges(pair_counter)
                ! Sharp cutoff -- for smooth cutoff, use the DSF method instead
                if (rij > rcut) then
                    cycle
                end if
                ! Inner cutoff -- damp singularity at r->0
                if (rij < rcut_in) then
                    if (rij < (rcut_in - rcin_width)) then
                        cycle
                    else ! we are in the transition region
                        inner_damp_ij = damping_function_cosine(rij, rcut_in - rcin_width, rcut_in)
                    end if
                else
                    inner_damp_ij = 1.0_dp
                end if
                ! Note the "pair energy" does NOT include damping
                pair_energy =  center_term * neigh_charge * pair_energy_direct(rij)
                ! We use half the pair energy here, since we double-count
                local_energies(center_i) = local_energies(center_i) + &
                    0.5_dp * pair_energy * inner_damp_ij
                if (do_gradients) then
                    ! ...but we don't double-count the centers (?)
                    fij_vec = rij_vec * inner_damp_ij * der_pair_energy_direct(rij) / rij
                    if (rij < rcut_in) then
                        fij_vec = fij_vec + der_damping_function_cosine(&
                                rij, rcut_in - rcin_width, rcut_in)&
                            * pair_energy * rij_vec / rij
                    end if
                    ! TODO omit forces on periodic replicas? They should cancel in any case.
                    forces(:, center_i) = forces(:, center_i) + fij_vec
                    ! TODO check virial sign convention
                    ! This convention aligns with more positive virials indicating greater internal pressure
                    virial = virial - outer_prod(fij_vec, rij_vec)
                    ! Accumulate ij-pair prefactor for variable-charge gradient term
                    ! (the inner-damping factor is 1.0 outside the inner cutoff)
                    ! Note we need to divide by the center charge later, since it's
                    ! included in the pair_energy
                    vc_grad_prefactor(center_i) = vc_grad_prefactor(center_i) + &
                            inner_damp_ij * pair_energy
                end if
            end do
            ! Now add the variable-charge gradient contribution
            ! This time they are added to the neighbours, not the centers
            if (do_gradients) then
                ! Now we iterate over SOAP neighbours only, but _including_ the center atom
                do soap_neigh_seq = 1, n_neigh(center_i)
                    soap_pair_counter = soap_pair_counter + 1
                    ! ???
                    soap_neigh_id = modulo(neighbors_list(soap_pair_counter) - 1, n_sites_global) + 1
                    ! Avoid computing grad contribution for neighbours that are not in SOAP cutoff
                    if(all(abs(charge_gradients(:, soap_pair_counter)) < FLOAT_ZERO) ) then
                        continue
                    end if
                    ! uses gradient _of_ charge i (center) w.r.t. atom k (soap neighbour)
                    fki_vec = -1.0_dp * vc_grad_prefactor(center_i) / charges(center_i) &
                                      * charge_gradients(:, soap_pair_counter)
                    forces(:, soap_neigh_id) = forces(:, soap_neigh_id) + fki_vec
                    ! Different sign than above because the position vector is reversed
                    ! (f_ki versus r_ik)
                    virial = virial + outer_prod(fki_vec, xyz(:, soap_pair_counter))
                end do
            end if
        end do
        ! Symmetrize the viral (is this necessary?)
        virial = 0.5_dp * (virial + transpose(virial))

        if(allocated(vc_grad_prefactor)) deallocate(vc_grad_prefactor)

    end subroutine

    ! TODO this is basically copy-pasted from the above routine.  We could
    !      definitely avoid a lot of reuse with some refactoring, but I'm not
    !      yet sure how to do it without incurring a potentially large
    !      performance penalty.
    subroutine compute_coulomb_dsf(charges, charge_gradients, &
                                   n_neigh, neighbors_list, &
                                   dsf_alpha, rcut, rcut_in, rcin_width, rjs, xyz, &
                                   neighbor_charges, do_gradients, &
                                   local_energies, forces, virial)
        implicit none
        real(dp), dimension(:), intent(in) :: charges, neighbor_charges, rjs
        real(dp), dimension(:,:), intent(in) :: charge_gradients, xyz
        integer, dimension(:), intent(in) :: n_neigh, neighbors_list
        real(dp), intent(in) :: rcut, rcut_in, rcin_width, dsf_alpha
        logical, intent(in) :: do_gradients

        ! inout because they are initialized outside this procedure and filled with zeros
        real(dp), intent(inout), dimension(:) :: local_energies
        real(dp), intent(inout), dimension(:,:) :: forces
        real(dp), intent(inout), dimension(3,3) :: virial

        integer :: center_i, neigh_id, soap_neigh_id, neigh_seq, soap_neigh_seq
        integer :: n_sites_global, n_sites_this, pair_counter, soap_pair_counter
        real(dp) :: rij, center_term, pair_energy, neigh_charge, inner_damp_ij
        real(dp) :: pair_energy_rcut, der_pair_energy_rcut
        real(dp), dimension(3) :: rij_vec, fij_vec, fki_vec
        real(dp), dimension(:), allocatable :: vc_grad_prefactor

        n_sites_this = size(n_neigh)
        n_sites_global = size(forces, 2)
        allocate(vc_grad_prefactor(1:n_sites_this))
        vc_grad_prefactor = 0.0_dp
        pair_energy_rcut = pair_energy_dsf(rcut, alpha)
        der_pair_energy_rcut = der_pair_energy_dsf(rcut, alpha)

        pair_counter = 0
        soap_pair_counter = 0
        do center_i = 1, n_sites_this
            pair_counter = pair_counter + 1
            !soap_pair_counter = soap_pair_counter + 1 ! No, because we include the center as a SOAP neighbour
            ! First we precompute q_i/4πε_0
            ! TODO this is where we add an effective dielectric constant to scale the interaction
            center_term = charges(center_i) * COUL_CONSTANT
            do neigh_seq = 2, n_neigh(center_i)
                pair_counter = pair_counter + 1
                neigh_id = modulo(neighbors_list(pair_counter) - 1, n_sites_global) + 1
                rij = rjs(pair_counter)
                rij_vec = xyz(:, pair_counter)
                neigh_charge = neighbor_charges(pair_counter)
                if (rij > rcut) then
                    cycle
                end if
                ! Inner cutoff -- damp singularity at r->0
                if (rij < rcut_in) then
                    if (rij < (rcut_in - rcin_width)) then
                        cycle
                    else ! we are in the transition region
                        inner_damp_ij = damping_function_cosine(rij, rcut_in - rcin_width, rcut_in)
                    end if
                else
                    inner_damp_ij = 1.0_dp
                end if
                ! Note the "pair energy" does NOT include damping
                pair_energy =  center_term * neigh_charge * &
                    (pair_energy_dsf(rij, dsf_alpha) - pair_energy_rcut &
                     - der_pair_energy_rcut * (rij - rcut))
                ! We use half the pair energy here, since we double-count
                local_energies(center_i) = local_energies(center_i) + &
                    0.5_dp * pair_energy * inner_damp_ij
                if (do_gradients) then
                    fij_vec = rij_vec / rij * inner_damp_ij * &
                        (der_pair_energy_dsf(rij, dsf_alpha) - der_pair_energy_rcut)
                    if (rij < rcut_in) then
                        fij_vec = fij_vec + der_damping_function_cosine(&
                                rij, rcut_in - rcin_width, rcut_in)&
                            * pair_energy * rij_vec / rij
                    end if
                    forces(:, center_i) = forces(:, center_i) + fij_vec
                    ! This sign convention aligns with more positive virials indicating greater internal pressure
                    virial = virial - outer_prod(fij_vec, rij_vec)
                    ! Accumulate ij-pair prefactor for variable-charge gradient term
                    ! (the inner-damping factor is 1.0 outside the inner cutoff)
                    ! Note we need to divide by the center charge later, since it's
                    ! included in the pair_energy
                    vc_grad_prefactor(center_i) = vc_grad_prefactor(center_i) + &
                            inner_damp_ij * pair_energy
                end if
            end do
            ! Now add the variable-charge gradient contribution
            ! This time they are added to the neighbours, not the centers
            if (do_gradients) then
                ! Now we iterate over SOAP neighbours only, but _including_ the center atom
                do soap_neigh_seq = 1, n_neigh(center_i)
                    soap_pair_counter = soap_pair_counter + 1
                    ! ???
                    soap_neigh_id = modulo(neighbors_list(soap_pair_counter) - 1, n_sites_global) + 1
                    ! Avoid computing grad contribution for neighbours that are not in SOAP cutoff
                    if(all(abs(charge_gradients(:, soap_pair_counter)) < FLOAT_ZERO) ) then
                        continue
                    end if
                    ! uses gradient _of_ charge i (center) w.r.t. atom k (soap neighbour)
                    fki_vec = -1.0_dp * vc_grad_prefactor(center_i) / charges(center_i) &
                                      * charge_gradients(:, soap_pair_counter)
                    forces(:, soap_neigh_id) = forces(:, soap_neigh_id) + fki_vec
                    ! Different sign than above because the position vector is reversed
                    ! (f_ki versus r_ik)
                    virial = virial + outer_prod(fki_vec, xyz(:, soap_pair_counter))
                end do
            end if
        end do
        ! Symmetrize the viral (is this necessary?)
        virial = 0.5_dp * (virial + transpose(virial))

        if(allocated(vc_grad_prefactor)) deallocate(vc_grad_prefactor)

    end subroutine

    ! A simple half-cosine damping function - designed to keep the short-range
    ! electrostatic terms from blowing up and to be computationally convenient,
    ! not to be physically realistic.  The short-range energy should always be
    ! corrected by a GAP afterwards.
    function damping_function_cosine(distance, r_inner, r_outer)
        real(dp), intent(in) :: distance, r_inner, r_outer
        real(dp) :: damping_function_cosine

        if (distance < r_inner) then
            damping_function_cosine = 0
        else if (distance > r_outer) then
            damping_function_cosine = 1
        else
            damping_function_cosine = 0.5 - 0.5*dcos((distance - r_inner) * PI &
                                                     / (r_outer - r_inner))
        end if
    end function

    ! The r-derivative of above function
    function der_damping_function_cosine(distance, r_inner, r_outer)
        real(dp), intent(in) :: distance, r_inner, r_outer
        real(dp) :: der_damping_function_cosine

        if (distance < r_inner) then
            der_damping_function_cosine = 0
        else if (distance > r_outer) then
            der_damping_function_cosine = 0
        else
            der_damping_function_cosine = 0.5 / (r_outer - r_inner) * dsin(&
                (distance - r_inner) * PI / (r_outer - r_inner))
        end if
    end function

end module
