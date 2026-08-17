! A synthetic system and neighbour list for driving get_soap directly.
!
! Non-periodic on purpose: it keeps the caller's supercell folding out of the
! picture, and edge atoms with few neighbours are a useful case in their own
! right. The neighbour list follows the convention get_soap expects -- for each
! site, neighbour 1 is the site itself at rj = 0, then every atom inside
! rcut_max, in a deterministic order so that pair indices survive the small
! displacements the finite-difference test makes.
module gharness
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   integer, parameter :: dp = real64

! ---- system
   integer :: n_atoms, n_species
   real(dp), allocatable :: pos(:, :)          ! (3, n_atoms)
   integer, allocatable :: sp_of_atom(:)      ! 1..n_species

! ---- neighbour list, rebuilt from pos by build_neighbours()
   integer :: n_sites, n_atom_pairs
   integer, allocatable :: n_neigh(:), nl_atom(:)   ! nl_atom(k) = which atom pair k points at
   real(dp), allocatable :: rjs(:), thetas(:), phis(:)
   logical, allocatable :: mask(:, :)
   integer, allocatable :: species(:, :), species_multiplicity(:)

! ---- descriptor hyperparameters
   integer :: l_max, n_max, n_soap_uncompressed
   integer, allocatable :: alpha_max(:)
   real(dp), allocatable :: rcut_hard(:), rcut_soft(:), nf(:), global_scaling(:)
   real(dp), allocatable :: atom_sigma_r(:), atom_sigma_r_scaling(:)
   real(dp), allocatable :: atom_sigma_t(:), atom_sigma_t_scaling(:)
   real(dp), allocatable :: amplitude_scaling(:), central_weight(:)
   integer :: radial_enhancement
   character(len=64) :: basis, scaling_mode
   real(dp) :: rcut_max

! ---- trivial compression
   integer :: compress_P_nonzero, n_soap_compressed
   integer, allocatable :: compress_P_i(:), compress_P_j(:)
   real(dp), allocatable :: compress_P_el(:)

contains

   subroutine setup_system(na, seed, box)
      integer, intent(in) :: na, seed
      real(dp), intent(in) :: box
      integer :: i, sz
      integer, allocatable :: sd(:)
      real(dp) :: u(3)

      n_atoms = na
      n_species = 2
      call random_seed(size=sz); allocate (sd(sz)); sd = seed
      call random_seed(put=sd); deallocate (sd)

      if (allocated(pos)) deallocate (pos, sp_of_atom)
      allocate (pos(3, n_atoms), sp_of_atom(n_atoms))
      do i = 1, n_atoms
         call random_number(u)
         pos(:, i) = u*box
         sp_of_atom(i) = 1 + mod(i, 2)
      end do

      l_max = 8
      if (allocated(alpha_max)) deallocate (alpha_max, rcut_hard, rcut_soft, nf, global_scaling, &
                                            atom_sigma_r, atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, &
                                            amplitude_scaling, central_weight)
      allocate (alpha_max(n_species), rcut_hard(n_species), rcut_soft(n_species), nf(n_species), &
                global_scaling(n_species), atom_sigma_r(n_species), atom_sigma_r_scaling(n_species), &
                atom_sigma_t(n_species), atom_sigma_t_scaling(n_species), &
                amplitude_scaling(n_species), central_weight(n_species))
      alpha_max = 8
      rcut_hard = 4.5d0
      rcut_soft = 4.0d0
      nf = 4.d0
      global_scaling = 1.d0
      atom_sigma_r = 0.5d0
      atom_sigma_r_scaling = 0.05d0
      atom_sigma_t = 0.2d0
      atom_sigma_t_scaling = 0.1d0
      amplitude_scaling = 1.d0
      central_weight = 1.d0
      radial_enhancement = 1
      basis = "poly3gauss"
      scaling_mode = "polynomial"
      rcut_max = maxval(rcut_hard)
      n_max = sum(alpha_max)
      n_soap_uncompressed = n_max*(n_max + 1)/2*(l_max + 1)

      call build_compression()
   end subroutine setup_system

! The "trivial" recipe from soap_turbo_compress.f90: keep every (n,m,l) whose
! n or m is the first radial index of some species.
   subroutine build_compression()
      integer :: i, n, m, l, k, counter
      integer, allocatable :: pivot(:)
      allocate (pivot(n_species))
      pivot(1) = 1
      do i = 1, n_species - 1
         pivot(i + 1) = pivot(i) + alpha_max(i)
      end do
      counter = 0
      do n = 1, n_max
         do m = n, n_max
            do l = 0, l_max
               if (any(n == pivot) .or. any(m == pivot)) counter = counter + 1
            end do
         end do
      end do
      compress_P_nonzero = counter
      n_soap_compressed = counter
      if (allocated(compress_P_i)) deallocate (compress_P_i, compress_P_j, compress_P_el)
      allocate (compress_P_i(counter), compress_P_j(counter), compress_P_el(counter))
      counter = 0; k = 1
      do n = 1, n_max
         do m = n, n_max
            do l = 0, l_max
               if (any(n == pivot) .or. any(m == pivot)) then
                  counter = counter + 1
                  compress_P_i(counter) = counter
                  compress_P_j(counter) = k
                  compress_P_el(counter) = 1.d0
               end if
               k = k + 1
            end do
         end do
      end do
      deallocate (pivot)
   end subroutine build_compression

! Rebuild rjs / thetas / phis / mask from pos. Deterministic order.
   subroutine build_neighbours()
      integer :: i, j, k, cnt
      real(dp) :: d(3), r

      n_sites = n_atoms
      if (allocated(n_neigh)) deallocate (n_neigh)
      allocate (n_neigh(n_sites))

!   pass 1: count
      n_atom_pairs = 0
      do i = 1, n_sites
         cnt = 1                                   ! the site itself
         do j = 1, n_atoms
            if (j == i) cycle
            d = pos(:, j) - pos(:, i)
            if (dsqrt(sum(d*d)) < rcut_max) cnt = cnt + 1
         end do
         n_neigh(i) = cnt
         n_atom_pairs = n_atom_pairs + cnt
      end do

      if (allocated(rjs)) deallocate (rjs, thetas, phis, mask, nl_atom)
      allocate (rjs(n_atom_pairs), thetas(n_atom_pairs), phis(n_atom_pairs))
      allocate (mask(n_atom_pairs, n_species), nl_atom(n_atom_pairs))
      mask = .false.

!   pass 2: fill
      k = 0
      do i = 1, n_sites
         k = k + 1
         rjs(k) = 0.d0; thetas(k) = 0.d0; phis(k) = 0.d0
         nl_atom(k) = i
         mask(k, sp_of_atom(i)) = .true.
         do j = 1, n_atoms
            if (j == i) cycle
            d = pos(:, j) - pos(:, i)
            r = dsqrt(sum(d*d))
            if (r >= rcut_max) cycle
            k = k + 1
            rjs(k) = r
            thetas(k) = dacos(d(3)/r)
            phis(k) = datan2(d(2), d(1))
            nl_atom(k) = j
            mask(k, sp_of_atom(j)) = .true.
         end do
      end do

      if (allocated(species)) deallocate (species, species_multiplicity)
      allocate (species(1, n_sites), species_multiplicity(n_sites))
      do i = 1, n_sites
         species(1, i) = sp_of_atom(i)
         species_multiplicity(i) = 1
      end do
   end subroutine build_neighbours

end module gharness
