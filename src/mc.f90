! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mc.f90, is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   This file was authored by Tigany Zarrouk
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

module mc

  use neighbors
  use md


  contains

!   These are the routines for calculating the monte-carlo acceptance criteria for different move types
  subroutine monte_carlo_insertion(p_accept, e_new, e_prev, temp, mu, m, volume, volume_bias, N_exch)
    implicit none

    real*8, intent(in) :: e_new, e_prev, temp, mu, m, volume, volume_bias
    integer, intent(in) :: N_exch
    real*8 :: lam
    real*8 :: kB = 8.617333262e-5, hbar = 6.582119569e-1, pi=3.1415926535
    real*8, intent(out) :: p_accept
    ! mass has units eV*fs^2/A^2
    ! hbar in eV.fs
    ! lam is thermal debroglie wavelength
    lam = sqrt( ( 2.0 * pi * hbar*hbar ) / ( m * kB * temp ) )

    p_accept =  (volume*volume_bias) / (lam**3.0 * (N_exch + 1)) * &
                exp( -( e_new - e_prev - mu  ) / (kB * temp) )

  end subroutine monte_carlo_insertion


  subroutine monte_carlo_removal(p_accept, e_new, e_prev, temp, mu, m, volume, volume_bias, N_exch)
    implicit none

    real*8, intent(in) :: e_new, e_prev, temp, mu, m, volume, volume_bias
    integer, intent(in) :: N_exch
    real*8 :: lam
    real*8 :: kB = 8.617333262e-5, hbar = 6.582119569e-1, pi=3.1415926535
    real*8, intent(out) :: p_accept
    ! mass has units eV*fs^2/A^2
    ! hbar in eV.fs
    ! lam is thermal debroglie wavelength
    lam = sqrt( ( 2.0 * pi * hbar*hbar ) / ( m * kB * temp ) )

    p_accept = (lam**3.0 * N_exch ) / (volume*volume_bias)  * &
                exp( -( e_new - e_prev + mu  ) / (kB * temp) )

  end subroutine monte_carlo_removal

  subroutine monte_carlo_move(p_accept, e_new, e_prev, temp)
    implicit none

    real*8, intent(in) :: e_new, e_prev, temp
    real*8 :: kB = 8.617333262e-5
    real*8, intent(out) :: p_accept
    p_accept = exp( -( e_new - e_prev ) / (kB * temp) )

  end subroutine monte_carlo_move

  subroutine monte_carlo_volume(p_accept, e_new, e_prev, temp, V_new, V_prev, P, N_exch)
    implicit none

    real*8, intent(in) :: e_new, e_prev, temp, V_new, V_prev, P
    real*8 :: kB = 8.617333262e-5, beta
    integer, intent(in) :: N_exch
    real*8, intent(out) :: p_accept

    beta = (1./(kB * temp))
    p_accept = exp( - beta * ( (e_new - e_prev) + &
                      P * ( V_new-V_prev ) + &
                     -(N_exch+1) * log( V_new/V_prev )/ beta ) )
  end subroutine monte_carlo_volume

  ! Not implemented this yet as it is a little complicated as one has
  ! to recalculate the neighbours but with the mc_min_dist as the
  ! cutoff!

  ! subroutine get_volume_bias(volume_bias, mc_min_dist, positions, &
  !      n_neigh, neighbor_list, n_atom_pairs)
  !   implicit none

  !   real*8, intent(in) :: mc_min_dist, positions(:,:)
  !   real*8, intent(out) :: volume_bias
  !   integer :: n_neigh(:), neighbor_list(:), n_atom_pairs

  !   ! is the number of neighbors for a given atom index including itself


  ! end subroutine get_volume_bias



    

  subroutine get_mc_acceptance(mc_move, p_accept, energy, energy_prev, temp, &
       mu, n_mc_species, v_uc, v_uc_prev, mass, pressure)
    implicit none

    character*32, intent(in) :: mc_move
    real*8, intent(out) :: p_accept
    real*8, intent(in) ::  energy, energy_prev, temp, &
         mu, v_uc, v_uc_prev, mass, pressure
    integer, intent(in) :: n_mc_species



    if (mc_move == "move")then
       call monte_carlo_move(p_accept, energy, energy_prev, temp)
       !     Not implemented the volume bias yet
    else if (mc_move == "insertion")then
       call monte_carlo_insertion(p_accept, energy, energy_prev, temp, mu, &
            mass, v_uc, 1.0d0, n_mc_species)
    else if (mc_move == "removal")then
       call monte_carlo_removal(p_accept, energy, energy_prev, temp, mu, &
            mass, v_uc, 1.0d0, n_mc_species)
    else if (mc_move == "volume")then
       call monte_carlo_volume(p_accept, energy, energy_prev, temp, &
            v_uc, v_uc_prev, pressure, n_mc_species)
    end if

  end subroutine get_mc_acceptance

  subroutine mc_get_atom_disp(n_sites, move_max, idx, disp, d_disp)
    implicit none


    real*8, intent(out) :: disp(1:3), d_disp
    integer, intent(out) :: idx
    real*8, intent(in) :: move_max
    integer, intent(in) :: n_sites
    real*8 :: ranf, ranv(1:3)

!          Choose a random atom and displace it
    call random_number(ranf)
    idx = floor( ranf * n_sites ) + 1
    call random_number(ranv)
!          displacement
    disp = move_max * ( ranv - 0.5d0 )
    d_disp = norm2(disp)
!          Now pick a random index and displace

  end subroutine mc_get_atom_disp


  subroutine get_mc_move(n_mc_species, mc_types, mc_move)
    implicit none

    integer, intent(in) :: n_mc_species
    integer :: n_mc
    character*32, intent(in) ::  mc_types(:)
    character*32, intent(out) :: mc_move
    real*8 :: ranf
    logical :: invalid_move, cant_remove

    invalid_move = .true.
    do while( invalid_move )
       call random_number(ranf)
       n_mc = floor( size( mc_types,1 ) * ranf ) + 1
       mc_move = mc_types(n_mc)

       cant_remove =  (n_mc_species == 0 .and. mc_move == "removal" )

       if(cant_remove)then
          write(*,*) "Cannot remove any mc species! Choosing another move!"
          invalid_move = .true.
       else
          invalid_move = .false.
       end if
    end do

  end subroutine get_mc_move


  subroutine mc_insert_site(mc_species, positions, ref_positions, idx, n_sites, a_box, b_box, c_box, &
       indices, species, ref_species, xyz_species, ref_xyz_species, min_dist )

    implicit none

    real*8, intent(inout) :: positions(:,:)
    real*8, intent(out) :: ref_positions(:,:)
    real*8, intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3), min_dist
    real*8 :: ranv(1:3)
    integer, intent(in) :: idx, n_sites, indices(1:3), ref_species(:)
    integer, intent(inout) :: species(:)
    character*8, intent(in) :: ref_xyz_species(:)
    character*8, intent(inout) :: xyz_species(:)
    character*32, intent(in) :: mc_species
    logical :: too_close

    positions(1:3,1:n_sites-1) = ref_positions(1:3,1:n_sites-1)
    positions(1:3,1:n_sites-1) = ref_positions(1:3,1:n_sites-1)
    xyz_species(1:n_sites-1)= ref_xyz_species(1:n_sites-1)
    species(1:n_sites-1)= ref_species(1:n_sites-1)

    xyz_species(n_sites) = mc_species

    too_close = .true.
    do while ( too_close )

       call random_number(ranv)
       positions(1,n_sites) = ranv(1)*norm2(a_box)
       positions(2,n_sites) = ranv(2)*norm2(b_box)
       positions(3,n_sites) = ranv(3)*norm2(c_box)

       call wrap_pbc(positions(1:3,1:n_sites), a_box/dfloat(indices(1)), &
            b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))

       call check_if_atoms_too_close(positions(1:3,n_sites), ref_positions, n_sites-1, &
            a_box, b_box, c_box, min_dist, too_close)

    end do

    write(*,*) "Found insertion position at ", positions(1:3,n_sites)

  end subroutine mc_insert_site


  subroutine check_if_atoms_too_close(position, ref_positions, n_sites, &
       a_box, b_box, c_box, min_dist, too_close)
    implicit none
    integer, intent(in) :: n_sites
    integer :: i, i_shift(1:3)
    real*8, intent(in) :: position(:), ref_positions(:,:), min_dist
    real*8, intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8 :: d, dist(1:3)
    logical, intent(inout) :: too_close

    too_close = .false.
    do i = 1, n_sites
       call get_distance(position(1:3), ref_positions(1:3,i), &
            a_box, b_box, c_box, (/ .true., .true., .true. /), dist, d, i_shift)

       if( d < min_dist )then
          write(*,*) "Insertion site too close, finding another one"
          too_close=.true.
          exit
       end if
    end do

  end subroutine check_if_atoms_too_close


end module mc
