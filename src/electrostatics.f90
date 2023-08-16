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
    ! it will probably just be Hartree / Bohr, since the charges are
    ! expressed in atomic units but the lengths in angstroms
    real(dp), parameter :: COUL_CONSTANT = HARTREE_EV / BOHR_ANG

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

        real(dp), intent(out), dimension(:) :: local_energies
        real(dp), intent(out), dimension(:,:) :: forces
        real(dp), intent(out), dimension(3,3) :: virial

        integer :: center_i, neigh_id, neigh_offset, neigh_seq, n_sites
        real :: rij, center_term, pair_energy
        real(dp), dimension(3) :: rij_vec, fij_vec, chg_grad_force

        n_sites = size(n_neigh)
        neigh_offset = 0

        do center_i = 1, n_sites
            ! First we precompute q_i/4πε_0
            center_term = charges(center_i) * COUL_CONSTANT
            local_energies(center_i) = 0.0_dp
            do neigh_seq = 1, n_neigh(center_i)
                ! I _think_ this is how the neighbourlist indexing works...
                neigh_id = neigh_seq + neigh_offset
                rij = rjs(neigh_id)
                rij_vec = xyz(:, neigh_id)
                ! Sharp cutoff -- for smooth cutoff, use the DSF method instead
                if (rij > rcut) then
                    continue
                end if
                ! Technically half the pair energy, since we double-count
                pair_energy = 0.5 * center_term * neighbor_charges(neigh_id) / rij
                local_energies(center_i) = local_energies(center_i) + pair_energy
                if (do_gradients) then
                    fij_vec = -2.0 * pair_energy * rij_vec / rij**2
                    ! Add the extra term due to the charge gradients
                    ! chg_grad_force = neigh_charges(neigh_id) * charge_gradients(1:3,center_i) / rij
                    ! TODO we actually need to loop through a third index, representing
                    ! the SOAP neighbors of the central atom, and compute the effect of
                    ! the central atom's charge gradient _w.r.t. that neighbor_
                    ! on the pair interaction.
                    ! Each SOAP neighbor's effects should be additive via the chain rule
                    ! Review previous notes on this!
                    forces(:,center_i) = forces(:,center_i) + fij_vec
                    ! TODO symmetrize like it's done elsewhere in the code?
                    virial = virial + outer_prod(fij_vec, rij_vec)
                end if
            end do
            neigh_offset = neigh_offset + n_neigh(center_i)
        end do

    end subroutine

end module
