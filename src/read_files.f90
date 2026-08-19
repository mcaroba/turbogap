! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_files.f90, is copyright (c) 2019-2023, Miguel A. Caro and
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

module read_files

   use printing, only: print_parameter, print_parameters, print_error, &
                       print_warning, print_note, print_message
   use read_utils, only: read_parameters, check_iostatus, check_file_exists
   use error, only: turbogap_abort
   use kinds

   use neighbors
   use types
   use splines
   use vdw
   use soap_turbo_compress_module
   use xyz_module
   use md

!  Constant lists of what the code implements. Module parameters rather than
!  locals rebuilt on every call, so the keyword-family subroutines can validate
!  against them without being handed the lists.
   character*32, parameter :: implemented_thermostats(1:3) = &
                              [character*32 :: "none", "berendsen", "bussi"]
   character*32, parameter :: implemented_barostats(1:2) = &
                              [character*32 :: "none", "berendsen"]
   character*32, parameter :: implemented_mc_types(1:8) = &
                              [character*32 :: "none", "move", "insertion", "removal", "relax", "md", &
                                                "swap", "volume"]

contains

!**************************************************************************
! This subroutine reads in the XYZ file
!
! WE NEED TO WRITE A PROPER EXTXYZ READER THAT CAN IDENTIFY WHICH COLUMN CONTAINS
! EACH PROPERTY. THIS SUBROUTINE CAN ONLY READ IN FILES WITH THE FOLLOWING CONVENTION:
!
! SPECIES X Y Z (VX VY VZ (FIXX FIXY FIXZ))
!
   subroutine read_xyz(filename, ase_format, all_atoms, do_timing, n_species, species_types, &
                       repeat_xyz, rcut_max, which_atom, positions, &
                       do_md, velocities, masses_types, masses, xyz_species, xyz_species_supercell, &
                       species, species_supercell, indices, a_box, b_box, c_box, n_sites, &
                       supercell_check_only, fix_atom, t_beg, write_masses, recalculate_supercell, &
                       randomize_velocities)

      implicit none

!   Input variables
      real(dp), intent(in) :: rcut_max
      real(dp), intent(in) :: masses_types(:)
      real(dp), intent(in) :: t_beg
      integer, intent(in) :: which_atom
      integer, intent(in) :: n_species
      character*8, intent(in) :: species_types(:)
      character*1024, intent(in) :: filename
      logical, intent(in) :: ase_format
      logical, intent(in) :: all_atoms
      logical, intent(in) :: do_timing
      logical, intent(in) :: do_md
      logical, intent(in) :: supercell_check_only
      logical, intent(in) :: recalculate_supercell

!   In and out variables
      real(dp), allocatable, intent(inout) :: positions(:, :)
      real(dp), allocatable, intent(inout) :: velocities(:, :)
      real(dp), allocatable, intent(inout) :: masses(:)
      real(dp), intent(inout) :: a_box(1:3)
      real(dp), intent(inout) :: b_box(1:3)
      real(dp), intent(inout) :: c_box(1:3)
      integer, allocatable, intent(inout) :: species(:)
      integer, allocatable, intent(inout) :: species_supercell(:)
      integer, intent(inout) :: n_sites
      integer, intent(inout) :: indices(1:3)
      character*8, allocatable, intent(inout) :: xyz_species(:)
      character*8, allocatable, intent(inout) :: xyz_species_supercell(:)
      logical, intent(inout) :: repeat_xyz
      logical, intent(inout) :: write_masses
      logical, allocatable, intent(inout) :: fix_atom(:, :)
      logical, intent(in) :: randomize_velocities

!   Internal variables
      real(dp), allocatable :: positions_supercell(:, :)
      real(dp), allocatable :: velocities_supercell(:, :)
      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: dist(1:3)
      real(dp) :: read_time
      real(dp) :: E_kinetic
      real(dp) :: instant_temp
      real(dp) :: kB = 8.6173303d-5
      real(dp) :: rjunk(1:3)
      real(dp) :: rjunk1d
      integer :: i
      integer :: iostatus
      integer :: j
      integer :: n_sites_supercell
      integer :: counter
      integer :: ijunk
      integer :: k2
      integer :: i2
      integer :: j2
      integer :: indices_prev(1:3)
      character*8 :: i_char
      character*128 :: cjunk
      character*128 :: cjunk_array(1:100)
      character*1024 :: cjunk1024
      character*1024 :: properties
      character*12800 :: cjunk_array_flat
      logical :: masses_from_xyz
      logical :: has_velocities
      logical :: ljunk(1:3)

      indices_prev = indices

      if (do_timing) then
         call cpu_time(time1)
      end if

      if (.not. supercell_check_only) then
         inquire (file=filename, number=i)
         if (i /= 11) then
            open (unit=11, file=filename, status="old")
         end if
         read (11, *) n_sites
         if (ase_format) then
            read (11, fmt='(A)') cjunk_array_flat
            cjunk_array = ""
            read (cjunk_array_flat, *, iostat=iostatus) cjunk_array(:)
!     Read in lattice vectors
            i = 0
            do
               i = i + 1
               cjunk = cjunk_array(i)
               call upper_to_lower_case(cjunk)
               if (cjunk(1:7) == "lattice") then
                  read (cjunk(10:), *) a_box(1)
                  read (cjunk_array(i + 1), *) a_box(2)
                  read (cjunk_array(i + 2), *) a_box(3)
                  read (cjunk_array(i + 3), *) b_box(1)
                  read (cjunk_array(i + 4), *) b_box(2)
                  read (cjunk_array(i + 5), *) b_box(3)
                  read (cjunk_array(i + 6), *) c_box(1)
                  read (cjunk_array(i + 7), *) c_box(2)
                  cjunk = adjustr(cjunk_array(i + 8))
                  read (cjunk(1:127), *) c_box(3)
                  exit
               end if
            end do
!     Read in properties string
            i = 0
            do
               i = i + 1
               cjunk = cjunk_array_flat(i:i + 9)
               call upper_to_lower_case(cjunk)
               if (cjunk == "properties") then
                  j = i + 9
                  do
                     j = j + 1
                     if (cjunk_array_flat(j:j) == " ") then
                        properties = cjunk_array_flat(i:j - 1)
                        call upper_to_lower_case(properties)
                        exit
                     end if
                  end do
                  exit
               end if
            end do
         else
            write (*, *) "You must use ASE's XYZ format so that I can read the lattice vectors <-- ERROR"
            stop
!    else
!      a_box = 0.d0
!      b_box = 0.d0
!      c_box = 0.d0
!      read(10, *) cjunk, cjunk, a_box(1), junk, junk, junk, b_box(2), junk, junk, junk, c_box(3)
         end if
         if (.not. recalculate_supercell) then
            if (allocated(positions)) deallocate (positions)
            if (allocated(xyz_species)) deallocate (xyz_species)
            if (allocated(species)) deallocate (species)
            allocate (positions(1:3, 1:n_sites))
            allocate (xyz_species(1:n_sites))
            allocate (species(1:n_sites))
            xyz_species = ""
            species = 0
!   We need to comment this out here for nested sampling
!    if( do_md )then
            if (.true.) then
               if (allocated(velocities)) deallocate (velocities)
               if (allocated(masses)) deallocate (masses)
               if (allocated(fix_atom)) deallocate (fix_atom)
               allocate (velocities(1:3, 1:n_sites))
               velocities = 0.d0
               allocate (masses(1:n_sites))
               masses = 0.d0
               masses_from_xyz = .false.
               allocate (fix_atom(1:3, 1:n_sites))
               fix_atom = .false.
            end if
!
!   I should raise a warning if the user is specifying more species than exist in the simulation,
!   since that would lead to arrays that are bigger than necessary                                  <------- Fix this
!    allocate( species(1:max_species_multiplicity, 1:n_sites) )
!    allocate( species_multiplicity(1:n_sites) )
!    species = 0
!    species_multiplicity = 0
            do i = 1, n_sites
               read (11, '(A)') cjunk1024
               if (do_md) then
                  call read_xyz_line(properties, cjunk1024, i_char, positions(1:3, i), velocities(1:3, i), fix_atom(1:3, i), &
                                     has_velocities, masses(i), masses_from_xyz)
                  if (masses_from_xyz) then
                     masses(i) = masses(i)*103.6426965268d0
                     write_masses = .true.
                  end if
               else
                  call read_xyz_line(properties, cjunk1024, i_char, positions(1:3, i), rjunk(1:3), ljunk(1:3), has_velocities, &
                                     rjunk1d, masses_from_xyz)
               end if
               do j = 1, n_species
                  if (trim(i_char) == trim(species_types(j))) then
!          species_multiplicity(i) = species_multiplicity(i) + 1
!          species(species_multiplicity(i), i) = j
                     xyz_species(i) = species_types(j)
                     species(i) = j
!         This is commented out because we also need masses with nested sampling when used in combination with MD
!          if( do_md .and. .not. masses_from_xyz )then
                     if (.not. masses_from_xyz) then
                        masses(i) = masses_types(j)
                     end if
!          exit
                  end if
               end do
!      if( species(1, i) == 0 )then
               if (xyz_species(i) == "") then
                  write (*, *) '                                       |'
                  write (*, *) 'ERROR: atom', i, 'has no known species |  <-- ERROR'
                  write (*, *) '                                       |'
                  write (*, *) '.......................................|'
                  stop
               end if
            end do
!   Randomize velocities if velocities are not provided
            if ((do_md .and. .not. has_velocities) .or. randomize_velocities) then
               write (*, *) '                                       |'
               write (*, *) 'WARNING: you have not provided initial |  <-- WARNING'
               write (*, *) 'velocities. I am randomizing them so   |'
               write (*, *) 'that they match your initial target    |'
               write (*, *) 'temperature:                           |'
               write (*, *) '                                       |'
               write (*, '(A, F16.4, A)') ' t_beg = ', t_beg, ' K             |'
               write (*, *) '                                       |'
               write (*, *) '.......................................|'

               if (.not. allocated(velocities)) then
                  allocate (velocities(1:3, 1:n_sites), source=0.0d0)
               end if

               if (allocated(velocities)) then
                  if (size(velocities, 2) /= n_sites) then
                     deallocate (velocities)
                     allocate (velocities(1:3, 1:n_sites), source=0.0d0)
                  end if
               end if

               call random_number(velocities)
               call remove_cm_vel(velocities(1:3, 1:n_sites), masses(1:n_sites))
               E_kinetic = 0.d0
               do i = 1, n_sites
                  E_kinetic = E_kinetic + 0.5d0*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
               end do
               instant_temp = 2.d0/3.d0/dfloat(n_sites - 1)/kB*E_kinetic
               velocities = velocities*dsqrt(t_beg/instant_temp)

               do i = 1, n_sites
                  E_kinetic = E_kinetic + 0.5d0*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
               end do
            end if
!   Check if there are more structures in the xyz file
            read (11, *, iostat=iostatus) cjunk
            if (iostatus == 0) then
               backspace (11)
               repeat_xyz = .true.
            else
               close (11)
               repeat_xyz = .false.
            end if
            indices_prev = 1
         end if
      end if

!   Now we construct a supercell of the required size to accommodate the given rcut_max
!   This needs to be done when the simulation box cannot accommodate one cutoff sphere

      a_box = a_box/dfloat(indices_prev(1))
      b_box = b_box/dfloat(indices_prev(2))
      c_box = c_box/dfloat(indices_prev(3))
      call number_of_unit_cells_for_given_cutoff(a_box, b_box, c_box, rcut_max, [.true., .true., .true.], indices)

      if (.not. supercell_check_only .or. (supercell_check_only .and. any(indices /= indices_prev)) &
          .or. recalculate_supercell) then
      if (indices(1) > 1 .or. indices(2) > 1 .or. indices(3) > 1) then
         n_sites_supercell = n_sites*indices(1)*indices(2)*indices(3)
         allocate (positions_supercell(1:3, 1:n_sites_supercell))
!     We need to comment this out here for nested sampling
!      if( do_md )then
         if (allocated(velocities)) then
            allocate (velocities_supercell(1:3, 1:n_sites_supercell))
         end if
!      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
!      species_supercell = 0
         if (allocated(xyz_species_supercell)) deallocate (xyz_species_supercell)
         if (allocated(species_supercell)) deallocate (species_supercell)
         allocate (xyz_species_supercell(1:n_sites_supercell))
         allocate (species_supercell(1:n_sites_supercell))
         xyz_species_supercell = ""
         species_supercell = 0
         counter = 0
         do i2 = 1, indices(1)
            do j2 = 1, indices(2)
               do k2 = 1, indices(3)
                  do i = 1, n_sites
                     counter = counter + 1
                     positions_supercell(1:3, counter) = positions(1:3, i) + dfloat(i2 - 1)*a_box(1:3) &
                                                         + dfloat(j2 - 1)*b_box(1:3) &
                                                         + dfloat(k2 - 1)*c_box(1:3)
!             We need to comment this out here for nested sampling
!              if( do_md )then
                     if (allocated(velocities)) then
                        velocities_supercell(1:3, counter) = velocities(1:3, i)
                     end if
!              species_supercell(:, counter) = species(:, i)
                     xyz_species_supercell(counter) = xyz_species(i)
                     species_supercell(counter) = species(i)
                  end do
               end do
            end do
         end do
         deallocate (positions)
         allocate (positions(1:3, 1:n_sites_supercell))
         positions(1:3, 1:n_sites_supercell) = positions_supercell(1:3, 1:n_sites_supercell)
         deallocate (positions_supercell)
!     We need to comment this out here for nested sampling
!      if( do_md )then
         if (allocated(velocities)) then
            deallocate (velocities)
            allocate (velocities(1:3, 1:n_sites_supercell))
            velocities(1:3, 1:n_sites_supercell) = velocities_supercell(1:3, 1:n_sites_supercell)
            deallocate (velocities_supercell)
         end if
         a_box = dfloat(indices(1))*a_box
         b_box = dfloat(indices(2))*b_box
         c_box = dfloat(indices(3))*c_box
      else
         n_sites_supercell = n_sites
         if (allocated(xyz_species_supercell)) deallocate (xyz_species_supercell)
         if (allocated(species_supercell)) deallocate (species_supercell)
         allocate (xyz_species_supercell(1:n_sites_supercell))
         allocate (species_supercell(1:n_sites_supercell))
         xyz_species_supercell = xyz_species
         species_supercell = species
!      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
!      species_supercell = species
         if (supercell_check_only) then
            allocate (positions_supercell(1:3, 1:n_sites_supercell))
            positions_supercell = positions(1:3, 1:n_sites_supercell)
            deallocate (positions)
            allocate (positions(1:3, 1:n_sites_supercell))
            positions(1:3, 1:n_sites_supercell) = positions_supercell(1:3, 1:n_sites_supercell)
            deallocate (positions_supercell)
!       We need to comment this out here for nested sampling
!        if( do_md )then
            if (allocated(velocities)) then
               allocate (velocities_supercell(1:3, 1:n_sites_supercell))
               velocities_supercell = velocities(1:3, 1:n_sites_supercell)
               deallocate (velocities)
               allocate (velocities(1:3, 1:n_sites_supercell))
               velocities(1:3, 1:n_sites_supercell) = velocities_supercell(1:3, 1:n_sites_supercell)
               deallocate (velocities_supercell)
            end if
         end if
      end if

!   This is perhaps not the most efficient way to select only one atom, fix in the future <----- FIX THIS
      if (.not. all_atoms) then
         n_sites = 1
         dist(1:3) = positions(1:3, 1)
         positions(1:3, 1) = positions(1:3, which_atom)
         positions(1:3, which_atom) = dist(1:3)

!      deallocate( species )
!      allocate( species(1:max_species_multiplicity, 1:n_sites) )
!      species(1:max_species_multiplicity, 1) = species_supercell(1:max_species_multiplicity, which_atom)
!      species_supercell(1:max_species_multiplicity, which_atom) = species_supercell(1:max_species_multiplicity, 1)
!      species_supercell(1:max_species_multiplicity, 1) = species(1:max_species_multiplicity, 1)

!      ijunk = species_multiplicity(which_atom)
!      deallocate( species_multiplicity )
!      allocate( species_multiplicity(1:n_sites) )
!      species_multiplicity(1) = ijunk
      end if
      else
      a_box = a_box*dfloat(indices_prev(1))
      b_box = b_box*dfloat(indices_prev(2))
      c_box = c_box*dfloat(indices_prev(3))
      end if

      if (do_timing) then
         call cpu_time(time2)
         read_time = time2 - time1
         write (*, *) '                                       |'
         write (*, *) 'Atoms timings (read):                  |'
         write (*, *) '                                       |'
         write (*, '(A, F11.3, A)') '  *) Read xyz file: ', read_time, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine
!**************************************************************************

!**************************************************************************
   subroutine read_exp_data(file_data, n_points, data)

      implicit none

!   Input variables
      character*1024, intent(in) :: file_data
!   Output variables
      real(dp), allocatable, intent(out) :: data(:, :)
      integer, intent(out) :: n_points

!   Internal variables
      integer :: i
      integer :: j
      integer :: iostatus
      integer :: dim
      integer :: unit_number

      ! if the file_data == none then we allocate and exit
      if (trim(file_data) == "none") then
         n_points = 1
         allocate (data(1:2, 1:n_points))
      else

         !   Read data file to figure out data file size
         open (newunit=unit_number, file=file_data, status="old")
         iostatus = 0
         n_points = -1
         do while (iostatus == 0)
            read (unit_number, *, iostat=iostatus)
            n_points = n_points + 1
         end do
         close (unit_number)

         allocate (data(1:2, 1:n_points))
         !     Read local_property data
         open (newunit=unit_number, file=file_data, status="old")
         do i = 1, n_points
            read (unit_number, *) data(1, i), data(2, i)
         end do
         close (unit_number)
      end if

   end subroutine read_exp_data

!**************************************************************************
!  Read one rigid molecule from a plain xyz file, for grand-canonical moves
!  that exchange whole molecules rather than single atoms.
!
!  The file is the ordinary two-header-lines-then-symbol-and-xyz format; the
!  comment line is ignored, and so is any cell it declares -- a molecule to be
!  inserted has no cell of its own. Positions are shifted to the centre of mass
!  on the way in, so that placing a copy is a rotation about the origin
!  followed by a translation of the centre of mass.
!
!  Every symbol has to appear in `species`, because the potential is indexed by
!  species and an insertion of an atom the potential does not know about would
!  be evaluated as whichever species happened to sit at index zero.
   subroutine read_mc_molecule(mol, file, species_types, masses_types, e0, n_species, rank)

      implicit none

      type(mc_molecule), intent(inout) :: mol
      character(len=*), intent(in) :: file
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp), allocatable, intent(in) :: masses_types(:)
      real(dp), allocatable, intent(in) :: e0(:)
      integer, intent(in) :: n_species
      integer, intent(in) :: rank
      character*128 :: cjunk
      character*8 :: symbol
      real(dp) :: com(1:3)
      real(dp) :: d
      integer :: iostatus
      integer :: i
      integer :: j

      call check_file_exists(file)

      open (unit=71, file=file, status="old")
      read (71, *, iostat=iostatus) mol%n_atoms
      if (iostatus /= 0 .or. mol%n_atoms < 1) then
         if (rank == 0) write (*, *) "ERROR -> could not read an atom count from mc molecule file ", trim(file)
         call turbogap_abort()
      end if
      read (71, "(A)", iostat=iostatus) cjunk

      if (allocated(mol%positions)) deallocate (mol%positions, mol%masses, mol%species, mol%xyz_species)
      allocate (mol%positions(1:3, 1:mol%n_atoms))
      allocate (mol%masses(1:mol%n_atoms))
      allocate (mol%species(1:mol%n_atoms))
      allocate (mol%xyz_species(1:mol%n_atoms))

      do i = 1, mol%n_atoms
         read (71, *, iostat=iostatus) symbol, mol%positions(1:3, i)
         if (iostatus /= 0) then
            if (rank == 0) write (*, *) "ERROR -> mc molecule file ", trim(file), &
               " ended before its declared atom count"
            call turbogap_abort()
         end if
         mol%xyz_species(i) = symbol
         mol%species(i) = 0
         do j = 1, n_species
            if (trim(adjustl(species_types(j))) == trim(adjustl(symbol))) mol%species(i) = j
         end do
         if (mol%species(i) == 0) then
            if (rank == 0) then
               write (*, *) "ERROR -> mc molecule file ", trim(file), " contains species ", trim(symbol)
               write (*, *) "which is not in the species list:"
               write (*, *) species_types
            end if
            call turbogap_abort()
         end if
         mol%masses(i) = masses_types(mol%species(i))
      end do
      close (71)

      mol%total_mass = sum(mol%masses(1:mol%n_atoms))
      mol%e0_total = 0.d0
      do i = 1, mol%n_atoms
         mol%e0_total = mol%e0_total + e0(mol%species(i))
      end do

      com = 0.d0
      do i = 1, mol%n_atoms
         com(1:3) = com(1:3) + mol%masses(i)*mol%positions(1:3, i)
      end do
      com = com/mol%total_mass
      do i = 1, mol%n_atoms
         mol%positions(1:3, i) = mol%positions(1:3, i) - com(1:3)
      end do

      mol%radius = 0.d0
      do i = 1, mol%n_atoms
         d = dsqrt(dot_product(mol%positions(1:3, i), mol%positions(1:3, i)))
         if (d > mol%radius) mol%radius = d
      end do

      mol%is_molecule = .true.
      mol%file = file

!     masses_types is already in the internal eV fs^2 / A^2, so total_mass is
!     too; the report converts back for the reader's benefit only.
      if (rank == 0) then
         write (*, '(1X,A,1X,A,1X,A,1X,I0,1X,A,1X,F0.4,1X,A)') "! mc molecule", trim(file), "read:", &
            mol%n_atoms, "atoms, mass", mol%total_mass/103.6426965268d0, "amu"
      end if

   end subroutine read_mc_molecule
!**************************************************************************

   subroutine write_exp_data(x, y, overwrite, filename, label)

      implicit none

!   Input variables
      character(len=*), intent(in) :: filename, label
!   Output variables
      real(dp), allocatable, intent(in) :: x(:)
      real(dp), allocatable, intent(in) :: y(:)
      logical, intent(in) :: overwrite
!   Internal variables
      integer :: i

      if (overwrite) then
         open (unit=200, file=filename, status="unknown")
         write (200, '(A,1X,A)') '# ', trim(label)
      else
         open (unit=200, file=filename, status="old", position="append")
         write (200, *) ' '
      end if

      do i = 1, size(x)
         write (200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
      end do
      close (200)

   end subroutine write_exp_data

   subroutine write_exp_datan(x, y, overwrite, filename, label)

      implicit none

!   Input variables
      character(len=*), intent(in) :: filename, label
!   Output variables
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      logical, intent(in) :: overwrite
!   Internal variables
      integer :: i

      if (overwrite) then
         open (unit=200, file=filename, status="unknown")
         write (200, '(A,1X,A)') '# ', trim(label)
      else
         open (unit=200, file=filename, status="old", position="append")
         write (200, *) ' '
      end if

      do i = 1, size(x)
         write (200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
      end do
      close (200)

   end subroutine write_exp_datan

!**************************************************************************

!**************************************************************************
   subroutine read_alphas_and_descriptors(file_desc, file_alphas, n_sparse, descriptor_type, alphas, Qs, cutoff)

      implicit none

!   Input variables
      character*1024, intent(in) :: file_desc
      character*1024, intent(in) :: file_alphas
      character(len=*), intent(in) :: descriptor_type

!   Output variables
      real(dp), allocatable, intent(out) :: alphas(:)
      real(dp), allocatable, intent(out) :: Qs(:, :)
      real(dp), allocatable, intent(out) :: cutoff(:)
      integer, intent(out) :: n_sparse

!   Internal variables
      integer :: i
      integer :: j
      integer :: iostatus
      integer :: dim
      integer :: unit_number

!   Read alphas to figure out sparse set size
      open (newunit=unit_number, file=file_alphas, status="old")
      iostatus = 0
      n_sparse = -1
      do while (iostatus == 0)
         read (unit_number, *, iostat=iostatus)
         n_sparse = n_sparse + 1
      end do
      close (unit_number)

!   Read descriptor vectors in spare set
      open (newunit=unit_number, file=file_desc, status="old")
      iostatus = 0
      i = -1
      do while (iostatus == 0)
         read (unit_number, *, iostat=iostatus)
         i = i + 1
      end do
      dim = i/n_sparse
      close (unit_number)

!   We do things differently for each descriptor
      if (descriptor_type == "soap_turbo") then
!     Allocate stuff
         allocate (alphas(1:n_sparse))
         allocate (Qs(1:dim, 1:n_sparse))
!     Read alphas SOAP
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i)
         end do
         close (unit_number)
!     Read sparse set descriptors
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            do j = 1, dim
               read (unit_number, *) Qs(j, i)
            end do
         end do
         close (unit_number)
      else if (descriptor_type == "distance_2b") then
         if (dim /= 1) then
            write (*, *) "ERROR: Bad 2b descriptor/alphas file(s), dimensions/n_sparse don't match number of data entries"
            stop
         end if
!     Allocate stuff
         allocate (alphas(1:n_sparse))
         allocate (cutoff(1:n_sparse))
         allocate (Qs(1:n_sparse, 1:1))
!     Read alphas 2b and cutoff
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i), cutoff(i)
         end do
         close (unit_number)
!     Read soap vectors in spare set
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            read (unit_number, *) Qs(i, 1)
         end do
         close (unit_number)
      else if (descriptor_type == "angle_3b") then
         if (dim /= 3) then
            write (*, *) "ERROR: Bad 3b descriptor/alphas file(s), dimensions/n_sparse don't match number of data entries"
            stop
         end if
!     Allocate stuff. NOTE: the array indices are the opposite of SOAP convention, this makes execution faster
         allocate (alphas(1:n_sparse))
         allocate (cutoff(1:n_sparse))
         allocate (Qs(1:n_sparse, 1:3))
!     Read alphas 3b and cutoff
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i), cutoff(i)
         end do
         close (unit_number)
!     Read soap vectors in spare set
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            do j = 1, 3
               read (unit_number, *) Qs(i, j)
            end do
         end do
         close (unit_number)
      end if

   end subroutine
!**************************************************************************

!**************************************************************************
! This reads the input file
   subroutine read_input_file(n_species, mode, params, rank)

      implicit none

!   Input variables
      integer, intent(in) :: n_species
      integer, intent(in) :: rank
      character(len=*) :: mode

!   Output variables
      type(input_parameters), intent(out) :: params

!   Internal variables
      real(dp) :: c6_ref
      real(dp) :: r0_ref
      real(dp) :: alpha0_ref
      real(dp) :: bsf
      real(dp) :: k
      integer :: iostatus
      integer :: i
      integer :: j
      integer :: i2
      integer :: nw
      integer :: iostatus2
      character*1024 :: long_line
      character*128, allocatable :: long_line_items(:)
      character*64 :: keyword
      character*64 :: cjunk
      character*64 :: keyword_notrim
      character*32 :: implemented_exp_observables(1:5)
      character*2 :: element
      character*1 :: keyword_first
      logical :: are_vdw_refs_read(1:3)
      logical :: valid_choice
      logical :: masses_in_input_file = .false.
      logical :: keyword_found

      implemented_exp_observables(1) = "xps"
      implemented_exp_observables(2) = "xrd"
      implemented_exp_observables(3) = "saxs"
      implemented_exp_observables(4) = "pair_distribution"
      implemented_exp_observables(5) = "structure_factor"

      k = 0.d0

!   Some defaults before reading the input file (the values in the input file will override them)
      if (mode == "md") then
         params%do_md = .true.
         params%do_prediction = .true.
         params%do_forces = .true.
         params%do_derivatives = .true.
      else if (mode == "mc") then
         params%do_mc = .true.
         params%do_prediction = .true.
         params%do_forces = .true.
         params%do_derivatives = .true.
      else if (mode == "soap") then
         params%write_soap = .true.
      else if (mode == "predict") then
         params%do_prediction = .true.
         params%do_forces = .true.
      end if

!   Let's allocate some arrays:
      allocate (params%species_types(1:n_species))
      allocate (params%masses_types(1:n_species))
      allocate (params%radii(1:n_species))
      allocate (params%e0(1:n_species))
      allocate (params%vdw_c6_ref(1:n_species))
      allocate (params%vdw_r0_ref(1:n_species))
      allocate (params%vdw_alpha0_ref(1:n_species))
!   Some defaults before reading from file
      params%masses_types = 0.d0
      params%radii = 0.5d0
      params%e0 = 0.d0
      params%vdw_c6_ref = 0.d0
      params%vdw_r0_ref = 0.d0
      params%vdw_alpha0_ref = 0.d0
      are_vdw_refs_read = .false.

!   Read the input file now
      iostatus = 0
      i2 = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         call upper_to_lower_case(keyword)
         keyword = trim(keyword)
         keyword_notrim = keyword
         keyword_notrim = adjustr(keyword_notrim)
!      i2 = len(trim(keyword))
         i2 = len(keyword_notrim)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!' .or. keyword == 'pot_file' &
             .or. keyword == 'n_species') cycle

         ! Offer the line to each family in turn. The first that recognises the
         ! keyword sets keyword_found; falling off the end means nothing claimed it.
         keyword_found = .false.
         call read_options_general(10, iostatus, rank, keyword, mode, params, n_species, keyword_found, masses_in_input_file)
         if (keyword_found) cycle
         call read_options_control(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_md(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_nested(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_mc(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_vdw(10, iostatus, rank, keyword, mode, params, n_species, keyword_found, are_vdw_refs_read)
         if (keyword_found) cycle
         call read_options_estat(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_exp(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_stopping(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_local_properties(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle
         call read_options_output(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle

!      Accepted here and ignored: these tune the GPU implementation and a CPU
!      build has nothing to tune. They are read anyway so that one input deck
!      runs on both branches -- without this, any deck written for the GPU
!      build aborts here with "I do not recognize the input file keyword", and
!      the two suites cannot share a case.
         call read_options_gpu(10, iostatus, rank, keyword, mode, params, n_species, keyword_found)
         if (keyword_found) cycle

         write (*, *) "ERROR: I do not recognize the input file keyword ", trim(keyword)
         call turbogap_abort()
      end do

!   Do some checks
      if (params%write_xyz == 0) then
         params%write_xyz = params%md_nsteps
      end if

      if (params%do_md) then
         params%do_prediction = .true.
         params%do_forces = .true.
         params%which_atom = 0
      end if

      if (params%do_forces) then
         params%do_derivatives = .true.
      end if

      if (params%which_atom /= 0) then
         params%all_atoms = .false.
      else
         params%all_atoms = .true.
      end if

      do i = 1, n_species
         call get_vdw_ref_params(params%species_types(i), c6_ref, r0_ref, alpha0_ref, rank)
         if (.not. are_vdw_refs_read(1)) then
            params%vdw_c6_ref(i) = c6_ref
         end if
         if (.not. are_vdw_refs_read(2)) then
            params%vdw_r0_ref(i) = r0_ref
         end if
         if (.not. are_vdw_refs_read(3)) then
            params%vdw_alpha0_ref(i) = alpha0_ref
         end if
      end do

!   If we don't use van der Waals, then unset the default cutoff
      if (params%vdw_type == "none") then
         params%vdw_rcut = 0.d0
      else
!     If van der Waals is enabled, make sure the inner and outer cutoff regions do not overlap
!     and other sanity checks
         if (params%vdw_rcut - params%vdw_buffer < params%vdw_rcut_inner + params%vdw_buffer_inner) then
            write (*, *) "ERROR: vdW inner and outer cutoff regions can't overlap. Check your vdw_* definitions"
            stop
         end if
      end if

!   Nested sampling checks
      if (params%do_nested_sampling) then
         if (params%thermostat /= "none") then
            write (*, *) '                                       |'
            write (*, *) 'WARNING: Nested sampling only works    |  <-- WARNING'
            write (*, *) '(currently) in combination with total  |'
            write (*, *) 'energy MD. The selected thermostat has |'
            write (*, *) 'been disabled.                         |'
         end if
!     Prepare directory where we create the latest version of the walkers
         call system("rm -rf walkers/")
         call system("mkdir -p walkers/")
      end if

!   An IR spectrum is a time-series observable: it comes from the dipole
!   autocorrelation over an ensemble of configurations, not from one frame's
!   structure. It uses the exp_* keywords -- exp_labels to declare it,
!   exp_energy_scales for its weight and ramp, exp_forces and exp_energies as
!   the gates -- but it must not be routed through the per-frame prediction
!   pipeline below, which exists to turn a single structure into a pair
!   distribution, a structure factor or a diffraction pattern. If IR is the
!   only observable, there is nothing for that pipeline to do.
!   A do_ir run exists to produce a spectrum, so asking for one and getting no
!   file would be a surprise. write_ir still turns it off for anyone who wants
!   the ensemble filled and nothing written.
      if (params%do_ir .and. .not. params%write_ir) params%write_ir = .true.

      params%do_exp_structural = .false.
      do i = 1, params%n_exp
         if (trim(params%exp_data(i)%label) /= "ir" .and. &
             trim(params%exp_data(i)%label) /= "none") params%do_exp_structural = .true.
      end do
      if (params%do_exp .and. .not. params%do_exp_structural) then
         params%do_exp = .false.
         if (rank == 0 .and. params%valid_ir) then
            write (*, *) '                                       |'
            write (*, *) ' IR is the only experimental observable |'
            write (*, *) ' asked for, so the per-frame prediction |'
            write (*, *) ' pipeline stays off. exp_energy_scales, |'
            write (*, *) ' exp_forces and exp_energies still      |'
            write (*, *) ' apply to it.                           |'
         end if
      end if

!   Experimental prediction checks
      if (params%do_exp) then
         if (rank == 0) write (*, *) '                                       |'
         if (rank == 0) write (*, *) ' Experimental prediction mode          |'
         do i = 1, params%n_exp
            if (trim(params%exp_data(i)%label) == "ir") cycle
            ! check if a user range has been submitted
            write (*, *) '                                       |'

            if (params%exp_data(i)%user_range) then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') 'User exp. range specified for:', trim(params%exp_data(i)%label), '     |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' WARNING!! This feature is obselete    |'
               if (rank == 0) write (*, *) '                                       |'
            else
               if (rank == 0) write (*, '(A,1X,A,1X,A)') 'Exp data range will be used for:', trim(params%exp_data(i)%label), ' |'
               if (rank == 0) write (*, '(A,1X,A,1X,A)') ' from the file:', trim(params%exp_data(i)%file_data), ' |'

            end if

            if (params%exp_data(i)%range_min == 0.d0 .and. params%exp_data(i)%range_max == 1.d0) then
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) 'WARNING: Data range being used for exp.|'
               if (rank == 0) write (*, *) ' observable is the default (0.0, 1.0)! |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' To modify specify:                    |'
          if (rank == 0) write (*, '(A,1X,A,1X,A)') '  `range_', trim(params%exp_data(i)%label), ' = {lower_bound} {upper_bound}` |'
               if (rank == 0) write (*, *) ' in the input file.                    |'
               if (rank == 0) write (*, *) '                                       |'
            end if

            if (trim(params%exp_data(i)%label) == 'pair_distribution') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting r_range_min/max ', '     |'
               ! Note: for consistency with the implementation, we can
               ! change the value of r_min/r_max such that the x_i
               ! generated
               ! by the bin_edges of the pair_distribution function
               ! match those
               ! of the actual experimental data

               params%do_pair_distribution = .true.

               params%r_range_min = params%exp_data(i)%range_min - &
                    & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
                    & (dfloat(2*(params%exp_data(i)%n_samples - 1)))

               params%r_range_max = params%exp_data(i)%range_min + &
                    & (dfloat(2*params%exp_data(i)%n_samples - 1)* &
                    & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
                    & (dfloat(2*(params%exp_data(i)%n_samples - 1))))

               params%pair_distribution_n_samples = params%exp_data(i)%n_samples
            elseif (trim(params%exp_data(i)%label) == 'xrd') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_xrd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               ! params%q_units = 'twotheta'
               params%xrd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples

            elseif (trim(params%exp_data(i)%label) == 'nd') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_nd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               ! params%q_units = 'twotheta'
               params%nd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples

            elseif (trim(params%exp_data(i)%label) == 'saxs') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params&
                    &%exp_data(i)%label), ' found, setting q_range_min&
                    &/max with q_units = "q"', ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_xrd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               params%q_units = 'q'
               params%xrd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples
            elseif (trim(params%exp_data(i)%label) == 'structure_factor') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units =&
                    & "q"', ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               params%q_units = 'q'
               params%structure_factor_n_samples = params%exp_data(i)%n_samples
               params%xrd_n_samples = params%exp_data(i)%n_samples
            end if

            if (rank == 0) write (*, '(A,1X,F12.6,1X,A,F12.6,1X,A)') ' min =', params&
                 &%exp_data(i)%range_min, ' max =', params%exp_data(i)&
                 &%range_max, ' |'

            if (rank == 0) write (*, '(A,1X,I8,1X,A)') ' n_samples   =', params%exp_data(i)%n_samples, '                |'
            if (rank == 0) write (*, '(A,1X,L4,1X,A)') ' compute_exp =', params%exp_data(i)%compute_exp, '                    |'

            if (.not. allocated(params%exp_energy_scales) .and. (params%exp_forces .or. params%mc_optimize_exp)) then
               if (rank == 0) write (*, *) 'WARNING: No energy scales set for exp .|'
               if (rank == 0) write (*, *) ' optimisation by forces / MC!          |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' To modify specify:                    |'
               if (rank == 0) write (*, '(A)') '  `exp_energy_scales = {E1} {E2}`  |'
               if (rank == 0) write (*, *) ' In the input file.                    |'
               if (rank == 0) write (*, *) ' (example above is for n_exp = 2)      |'
               if (rank == 0) write (*, *) '                                       |'
            end if

         end do
      end if

!   The Debye route goes from the positions straight to I(q), so it needs
!   neither the pair distribution nor the structure factor. Both get switched
!   on as a side effect of asking for XRD/ND -- by the do_xrd/do_nd keywords
!   and again by the exp_labels block above -- and leaving them on would make
!   xrd_debye a pure slowdown, computing a whole pdf and sf that nothing
!   reads. A deck that asked for either in as many words keeps it: wanting the
!   pair distribution written out is a legitimate reason to compute one.
!
!   This has to run after the exp_labels loop, not inside the keyword
!   handlers, because the deck may set xrd_debye either side of do_xrd.
      if (params%xrd_debye .and. (params%do_xrd .or. params%do_nd)) then
         if (params%do_pair_distribution .and. .not. params%do_pair_distribution_explicit) then
            params%do_pair_distribution = .false.
            if (rank == 0) write (*, *) 'NOTE: xrd_debye, so the pdf is not |'
            if (rank == 0) write (*, *) 'computed. Set do_pair_distribution |'
            if (rank == 0) write (*, *) '= .true. to get it anyway.         |'
         end if
         if (params%do_structure_factor .and. .not. params%do_structure_factor_explicit) then
            params%do_structure_factor = .false.
            if (rank == 0) write (*, *) 'NOTE: xrd_debye, so the structure  |'
            if (rank == 0) write (*, *) 'factor is not computed. Set        |'
            if (rank == 0) write (*, *) 'do_structure_factor = .true. to    |'
            if (rank == 0) write (*, *) 'get it anyway.                     |'
         end if
      end if

!   The exchangeable objects, resolved once the species table, the masses, e0
!   and mc_species are all in. Doing it here rather than in the mc_species
!   handler is what lets a deck give the keywords in any order.
      if (params%do_mc .and. params%n_mc_mu > 0) then
         do i = 1, params%n_mc_mu
            if (trim(params%mc_molecule_files(i)) /= "none") then
               call read_mc_molecule(params%mc_molecules(i), trim(params%mc_molecule_files(i)), &
                                     params%species_types, params%masses_types, params%e0, &
                                     n_species, rank)
               params%mc_exchange_mass(i) = params%mc_molecules(i)%total_mass
               params%mc_exchange_e0(i) = params%mc_molecules(i)%e0_total
            else
               do j = 1, n_species
                  if (trim(adjustl(params%species_types(j))) == trim(adjustl(params%mc_species(i)))) then
                     params%mc_exchange_mass(i) = params%masses_types(j)
                     params%mc_exchange_e0(i) = params%e0(j)
                  end if
               end do
               if (params%mc_exchange_mass(i) == 0.d0) then
                  if (rank == 0) then
                     write (*, *) "ERROR -> mc_species ", trim(params%mc_species(i)), &
                        " is not in the species list:"
                     write (*, *) params%species_types
                     write (*, *) "Give a molecule file in mc_molecule_files if it is not an atom."
                  end if
                  call turbogap_abort()
               end if
            end if
         end do
      end if

!   Monte-carlo checks
      if (params%do_mc) then
         do i = 1, params%n_mc_types
            if (params%mc_types(i) == "md") then
               if (params%thermostat == "none") then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: You need to specify a         |  <-- WARNING'
                  if (rank == 0) write (*, *) 'thermostat when using md type mc steps!|'
               end if
            end if

            if (params%mc_types(i) == "relax") then
               if (params%optimize == "none") then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: You need to specify an        |  <-- WARNING'
                  if (rank == 0) write (*, *) 'optimizer when using relax type mc     |'
                  if (rank == 0) write (*, *) 'steps!!                                |'
               end if
            end if

            if (params%mc_types(i) == "volume") then
               if (params%p_beg == 1.0d0) then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: p_beg is the default          |  <-- WARNING'
                  if (rank == 0) write (*, *) 'value of 1.0 bar. For MC volume moves  |'
                  if (rank == 0) write (*, *) 'please make sure this is specified!!   |'
               end if
               if (params%mc_lnvol_max == 0.01d0) then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: mc_lnvol_max is the default   |  <-- WARNING'
                  if (rank == 0) write (*, *) 'value of 0.01. For MC volume moves     |'
                  if (rank == 0) write (*, *) 'please make sure this is specified!!   |'
               end if

            end if
         end do

         do i = 1, n_species
            if (params%accessible_volume .and. (params%radii(i) == 0.5d0)) then
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) 'WARNING: radii for accessible volume   |  <-- WARNING'
               if (rank == 0) write (*, *) 'is the default value of 0.5A.          |'
               if (rank == 0) write (*, *) 'please make sure this correct!!        |'
            end if
         end do

      end if

!   Set the writeouts
      if (.not. params%do_md) then
!     Do not write temperature
         params%write_property(3) = .false.
!     Do not write pressure
         params%write_property(4) = .false.
!     Do not write time step
         params%write_property(5) = .false.
!     Do not write time
         params%write_property(6) = .false.
!     Do not write MD step
         params%write_property(11) = .false.
!     Do not write velocities
         params%write_array_property(3) = .false.
!     Do not write masses
         params%write_array_property(6) = .false.
!     Do not write fixes
         params%write_array_property(8) = .false.
      end if
      if (.not. params%do_forces) then
!     Do not write pressure
         params%write_property(4) = .false.
!     Do not write virial
         params%write_property(8) = .false.
!     Do not write stress
         params%write_property(9) = .false.
!     Do not write forces
         params%write_array_property(4) = .false.
      end if
      if (params%vdw_type == "none") then
!     Do not write Hirshfeld volume
         params%write_array_property(7) = .false.
      end if
!   Now individual flags
      if (.not. params%write_velocities) then
         params%write_array_property(3) = .false.
      end if
      if (.not. params%write_forces) then
         params%write_array_property(4) = .false.
      end if
      if (.not. params%write_local_energies) then
         params%write_array_property(5) = .false.
      end if
      if (.not. params%write_masses) then
         params%write_array_property(6) = .false.
      end if
      if (.not. params%write_fixes) then
         params%write_array_property(8) = .false.
      end if
      if (.not. params%write_pressure) then
         params%write_property(7) = .false.
      end if
      if (.not. params%write_virial) then
         params%write_property(8) = .false.
      end if
      if (.not. params%write_stress) then
         params%write_property(9) = .false.
      end if

!   Get masses from database
      if ((params%do_md .or. params%do_mc) .and. .not. masses_in_input_file) then
         if (rank == 0) then
            write (*, *) '                                       |'
            write (*, *) 'WARNING: you have not provided masses  |  <-- WARNING'
            write (*, *) 'in your input file. I am attempting to |'
            write (*, *) 'read them from a database. If you have |'
            write (*, *) 'provided masses in your XYZ file these |'
            write (*, *) 'values will be overwritten and you can |'
            write (*, *) 'safely disregard any further warnings  |'
            write (*, *) 'printed below if a given element is not|'
            write (*, *) 'in the database (usually because you   |'
            write (*, *) 'provided a non-standard name; note that|'
            write (*, *) 'element names are case sensitive).     |'
            write (*, *) '                                       |'
            write (*, *) '               Element      Mass (amu) |'
         end if
         do i = 1, n_species
            call get_atomic_mass(params%species_types(i), params%masses_types(i), valid_choice)
            if (rank == 0) then
               write (*, *) '                                       |'
               if (valid_choice) then
            write (*, '(A, A8, A, F15.6, A)') ' ', adjustr(params%species_types(i)), ' (in database) ', params%masses_types(i), ' |'
               else
           write (*, '(A, A8, A, F11.6, A)') ' ', adjustr(params%species_types(i)), ' (not in database) ', params%masses_types(i), &
                     ' |  <-- WARNING'
               end if
            end if
         end do
!     We convert the masses in amu to eV*fs^2/A^2
         params%masses_types = params%masses_types*103.6426965268d0
         if (rank == 0) then
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if

   end subroutine

!**************************************************************************
!  Input keywords for the system: files, species, masses, and global limits.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_general(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found, masses_in_input_file)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: masses_in_input_file

      character*64 :: cjunk

      !> @kw atoms_file, input_file
      !> Path to the extended-XYZ file holding the atomic configuration, or the trajectory of
      !> configurations to loop over in predict mode. Required.
      if (keyword == 'atoms_file' .or. keyword == 'input_file') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%atoms_file
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("atoms_file", params%atoms_file)
         !> @kw e0
         !> Constant energy offset per species, added to every atom of that species. One value per
         !> entry in species, in the same order. The GAP predicts energies relative to these, so they
         !> set the zero of the energy scale.
         !> @units eV
      else if (keyword == 'e0') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%e0(1:n_species)
         if (rank == 0) call print_parameters("e0", params%e0)
         !> @kw ir_lag_factor
         !> Ratio of the stored ensemble length to the longest correlation lag. The
         !> autocorrelation at lag tau is averaged over n_window - tau pairs, so a factor of 1
         !> would leave the longest lag estimated from a single pair. 2 or more.
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_resolution
      else if (keyword == 'ir_lag_factor') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_lag_factor
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_lag_factor", params%ir_lag_factor)
         !> @kw ir_match_scale
         !> Fit an overall scale factor between the computed and experimental spectra before
         !> comparing them. The computed spectrum is in arbitrary units, so with this off the
         !> loss compares two things on different scales and its gradient is meaningless. Turn
         !> it off only if the experimental spectrum has already been put on the same scale.
         !> @modes md
         !> @see exp_labels exp_energy_scales
      else if (keyword == 'ir_match_scale') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_match_scale
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_match_scale", params%ir_match_scale)
         !> @kw ir_nu_max
         !> Highest wavenumber to fit. This is what fixes the sampling interval: nothing above
         !> the Nyquist limit of the stored series can be represented, and power above it folds
         !> back into the fitted range, so a combination of ir_stride and md_step that
         !> cannot reach this value is refused rather than aliased.
         !> @units cm^-1
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_stride
      else if (keyword == 'ir_nu_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_nu_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_nu_max", params%ir_nu_max)
         !> @kw ir_nu_min
         !> Lowest wavenumber to fit. Experimental points outside [ir_nu_min, ir_nu_max]
         !> are dropped and take no part in the loss.
         !> @units cm^-1
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_nu_max
      else if (keyword == 'ir_nu_min') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_nu_min
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_nu_min", params%ir_nu_min)
         !> @kw ir_nu_power
         !> Exponent of the wavenumber prefactor in I(nu) = nu^p * FT[C(tau)]. The classical
         !> lineshape is p = 2; p = 0 compares the bare Fourier transform of the dipole
         !> autocorrelation instead.
         !> @modes md
         !> @see exp_labels exp_energy_scales
      else if (keyword == 'ir_nu_power') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_nu_power
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_nu_power", params%ir_nu_power)
         !> @kw ir_resolution
         !> Frequency resolution wanted, which fixes the longest correlation lag as
         !> 33356.41/(resolution * dt) with dt the interval between stored configurations. The
         !> ensemble is ir_lag_factor times that. Asking for fine resolution is expensive in
         !> memory and in how long the run must go before the first force is applied.
         !> @units cm^-1
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_lag_factor ir_stride
      else if (keyword == 'ir_resolution') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_resolution
         params%ir_resolution_set = .true.
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_resolution", params%ir_resolution)
         !> @kw do_ir
         !> Predict an IR spectrum from the trajectory. The total dipole is accumulated over the
         !> whole run and transformed into ir_spectrum.dat, whose header records the sampling
         !> interval, the Nyquist limit, the resolution and hence the band over which the result
         !> can be read. Needs a dipole model in the potential file, but no experimental spectrum
         !> and no exp_* keywords: nothing is fitted and no force is added. Implies write_ir, so
         !> ir_prediction.dat carries the same spectrum over the whole trajectory. The ensemble is
         !> the entire run rather than a rolling window, so the resolution follows from md_nsteps
         !> unless ir_resolution asks for one the run is long enough to give. Naming "ir" in
         !> exp_labels instead is the other thing -- that biases the trajectory towards an
         !> experiment.
         !> @modes md
         !> @see ir_stride ir_resolution ir_nu_min ir_nu_max ir_lag_factor ir_n_samples exp_labels
      else if (keyword == 'do_ir') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_ir
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_ir", params%do_ir)
         !> @kw ir_n_samples
         !> Number of points on the predicted spectrum's grid. The default, zero, puts one point
         !> per resolution element, which is all the transform can carry; a larger number draws the
         !> same information as a smoother curve. Prediction runs only -- a fit uses the
         !> experimental grid.
         !> @modes md
         !> @see do_ir ir_resolution
      else if (keyword == 'ir_n_samples') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_n_samples", params%ir_n_samples)
         !> @kw ir_restart_file
         !> Where the dipole history is written and read back. The ensemble is the expensive
         !> part of a MAD IR run -- resolving 4 cm^-1 at 1 fs sampling is 8 ps of trajectory --
         !> so without this a restart spends that long refilling before any force is applied. A
         !> file written with different sizing is refused rather than adopted, since it
         !> describes a different spectrum.
         !> @modes md
         !> @see exp_labels exp_energy_scales
      else if (keyword == 'ir_restart_file') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_restart_file
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_restart_file", params%ir_restart_file)
         !> @kw ir_stride
         !> MD steps between stored configurations. The sampling interval is this times md_step,
         !> and that interval is what sets the Nyquist limit, so a large stride is what makes
         !> ir_nu_max unreachable.
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_nu_max md_step
      else if (keyword == 'ir_stride') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_stride
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_stride", params%ir_stride)
         !> @kw ir_window
         !> Lag window applied to the autocorrelation before the cosine transform:
         !> "hann" (default), "bartlett", "lorch", "welch" or "none". This is not optional in
         !> substance -- "none" is a boxcar whose kernel is a sinc, with side lobes at -13.3 dB
         !> and ringing to -21.7%, which drives an absorption spectrum negative. The choice is
         !> which taper, not whether. Hann has the lowest far-field leakage of the set (13 dB
         !> below the next best, because its kernel falls as 1/f^3) and is the only one whose
         !> main lobe FWHM equals the resolution ir_resolution advertises. Bartlett is the one
         !> shape whose kernel (Fejer) is non-negative everywhere, so with ir_estimator =
         !> "biased" it makes a non-negative spectrum a theorem rather than an observation, at
         !> the cost of 5 dB more leakage. Lorch is sinc(pi tau / n_lag), the same modification
         !> function exp_utils applies to g(r) under structure_factor_window, offered so the IR
         !> path can make the approximation the XRD path makes; it buys 13% narrower bands for
         !> 5 dB more leakage and 1.8x more ringing. "none" is for demonstrating the ringing,
         !> not for production. See ana/formalism_windows.py for the measured kernels.
         !> @modes md
         !> @see exp_labels exp_energy_scales ir_resolution ir_estimator structure_factor_window
      else if (keyword == 'ir_window') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_window
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_window", params%ir_window)
         !> @kw ir_subtract_mean
         !> Correlate the dipole FLUCTUATION mu - <mu> rather than mu itself (default .true.).
         !> Linear response gives the spectrum as the transform of <dM(0).dM(t)> with
         !> dM = M - <M>; a static dipole does not absorb, and the mean drops out of the
         !> derivation for that reason. The usual defence for skipping it -- "a constant only
         !> puts a delta at nu = 0, which the nu^2 prefactor kills" -- is false, because only
         !> the |<mu>|^2 piece of the un-centred correlation is constant. The cross term is a
         !> partial sum of a mean-zero series divided by its own length, i.e. a random walk in
         !> tau whose amplitude GROWS with lag. On a 100-molecule water buffer at 2 fs,
         !> |<mu>|^2 was 86% of C(0) and past ~300 fs of lag the cross term was 17x the real
         !> correlation, shifting the spectrum by 19% at the O-H stretch. Set .false. only to
         !> reproduce results from before this was fixed.
         !> @modes md
         !> @see ir_estimator ir_taper_partial exp_labels
      else if (keyword == 'ir_subtract_mean') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_subtract_mean
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_subtract_mean", params%ir_subtract_mean)
         !> @kw ir_estimator
         !> How the autocorrelation is normalised: "biased" (default) divides C(tau) by N,
         !> "unbiased" divides by the number of pairs actually averaged, N - tau. The unbiased
         !> estimator is unbiased and is also not a positive semi-definite sequence, so its
         !> cosine transform can be negative -- and an absorption spectrum cannot be. The biased
         !> one IS positive semi-definite: its transform is a periodogram (Percival &
         !> Walden ch. 6; Numerical Recipes 13.4), at the price of scaling C(tau) by
         !> (1 - tau/N). Since the lag window already tapers far harder than that, nothing is
         !> lost. Note that guarantee covers the UNWINDOWED estimate: windowing convolves the
         !> spectrum with the window kernel, so a kernel that dips negative can still drive the
         !> result negative. Measured on one trajectory with every other fix on, ir_window =
         !> "none" gave 124 negative points of 400 while hann, lorch and bartlett gave zero.
         !> Combine with ir_window = "bartlett", whose Fejer kernel is non-negative everywhere,
         !> for an unconditional guarantee.
         !> @modes md
         !> @see ir_window ir_subtract_mean ir_lag_factor
      else if (keyword == 'ir_estimator') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_estimator
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_estimator", params%ir_estimator)
         !> @kw ir_taper_partial
         !> Rebuild the lag window over the lags actually available while the ensemble is still
         !> filling (default .true.). Only affects do_ir prediction runs, where the buffer fills
         !> over the whole trajectory. The window is sized for n_lag, but early on only
         !> n_stored - 1 lags exist, so with this off the autocorrelation is truncated at a
         !> point where the window is still near 1 -- a boxcar cut, which convolves the spectrum
         !> with a sinc and makes it ring. How much it matters depends on ir_estimator: dividing
         !> C(tau) by N rather than N - tau multiplies it by (1 - tau/N), which is itself a
         !> triangular taper, so the biased estimator applies an implicit Bartlett window that
         !> has already fallen to 1/N at the truncation point. Measured with the window sized
         !> for n_lag = 417 and only 51 frames stored, ir_estimator = unbiased gave 160 negative
         !> points of 400 and ir_estimator = biased gave 0. With the default biased estimator
         !> this switch therefore does not change the ringing; what it still buys is a header
         !> that quotes the resolution the transformed lags actually support, and correct
         !> behaviour if the unbiased estimator is selected.
         !> @modes md
         !> @see do_ir ir_window ir_lag_factor write_xyz
      else if (keyword == 'ir_taper_partial') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_taper_partial
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_taper_partial", params%ir_taper_partial)
         !> @kw ir_weight_by_spacing
         !> Weight each experimental point by the local grid spacing (default .true.), so that
         !> the loss is an integral over frequency rather than a sum over however the
         !> experimental file happened to be sampled. The fit grid is the experimental grid --
         !> interpolating the experiment would invent structure between its points and then fit
         !> to it -- but experimental grids are rarely uniform. The Downing & Williams water
         !> data shipped with tests/mad_ir is sampled at 23.9 cm^-1 across the O-H stretch and
         !> ~55 cm^-1 elsewhere, so 45% of an unweighted loss falls in the top quarter of the
         !> range: a weighting chosen by whoever digitised the paper, not by the physics.
         !> Weights are normalised to mean 1, so exp_energy_scales keeps its magnitude.
         !> @modes md
         !> @see exp_energy_scales ir_nu_min ir_nu_max exp_data_files
      else if (keyword == 'ir_weight_by_spacing') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_weight_by_spacing
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_weight_by_spacing", params%ir_weight_by_spacing)
         !> @kw ir_match_offset
         !> Fit an additive baseline b alongside the scale s, comparing s*I_calc + b with the
         !> experiment (default .false.; requires ir_match_scale). Needed when the experimental
         !> file has had a baseline removed such that its minimum is exactly zero, as
         !> tests/mad_ir's process_data.py does: the model spectrum is strictly positive, so
         !> without an offset the residual in the transparency window can never vanish, and a
         !> quadratic loss responds by shrinking the whole prediction -- fighting the very band
         !> intensity the bias is trying to build. Both s and b sit at their own least-squares
         !> optimum, so the gradient picks up no extra term from either.
         !> CAUTION: the two-parameter solve is unconstrained and can return a NEGATIVE scale
         !> when the predicted shape does not resemble the experiment -- on this potential,
         !> whose dipole model produces no vibrational bands, it returned s = -9.7e-9 with
         !> b = +2.70. Since dL/dI carries a factor of s, that would reverse the bias and drive
         !> the model away from the bands it is being asked to grow. The code detects s <= 0 and
         !> falls back to the scale-only fit with b = 0; the header reports the offset actually
         !> used, so a run that took the fallback is identifiable. If your fits keep falling
         !> back, the prediction and the experiment disagree in shape, not just in baseline.
         !> @modes md
         !> @see ir_match_scale exp_energy_scales exp_data_files
      else if (keyword == 'ir_match_offset') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_match_offset
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("ir_match_offset", params%ir_match_offset)
         !> @kw ir_write_spectrum
         !> The old name for write_ir, kept because inputs use it. Sets write_ir.
         !> @modes md
         !> @see write_ir exp_labels exp_energy_scales
      else if (keyword == 'ir_write_spectrum') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%ir_write_spectrum
         call check_iostatus(iostatus, keyword)
         if (params%ir_write_spectrum) params%write_ir = .true.
         if (rank == 0) call print_parameter("ir_write_spectrum", params%ir_write_spectrum)
         !> @kw masses
         !> Atomic mass of each species, one value per entry in species and in the same order. Read in
         !> amu and converted internally to eV fs^2 / A^2. If absent they are taken from the XYZ file
         !> instead.
         !> @units amu
         !> @see species
      else if (keyword == 'masses') then
         backspace (unit)
         call read_parameters(unit, iostatus, n_species, params%masses_types)
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameters("masses", params%masses_types)
!       We convert the masses in amu to eV*fs^2/A^2
         params%masses_types = params%masses_types*103.6426965268d0
         masses_in_input_file = .true.
         !> @kw max_gbytes_per_process
         !> Memory an MPI rank may use for the SOAP descriptor and its derivatives; the rank splits
         !> its atoms into as many batches as it takes to stay under this. Naming it in the input file
         !> also stops the run sizing the budget from the node automatically, so a value chosen
         !> deliberately is never overwritten.
         !> @units GB
         !> @see mem_fraction
      else if (keyword == 'max_gbytes_per_process') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%max_Gbytes_per_process
         call check_iostatus(iostatus, keyword)
!       Remember that the input named it, so gpu_memory_budget_init leaves it
!       alone. Someone who sets this keyword has a reason the node's memory
!       cannot know about.
         params%max_Gbytes_set = .true.
         if (rank == 0) call print_parameter("max_Gbytes_per_process", params%max_Gbytes_per_process)
         !> @kw mem_fraction
         !> Fraction of the node's memory to divide between the ranks on it when
         !> max_Gbytes_per_process was not given. Only consulted for that automatic budget, and
         !> ignored once max_Gbytes_per_process appears in the deck.
         !> @see max_gbytes_per_process
      else if (keyword == 'mem_fraction') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mem_fraction
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mem_fraction", params%mem_fraction)
         !> @kw neighbors_buffer
         !> Extra distance added to every cutoff when the neighbour lists are built, so that a list
         !> stays valid for several steps as atoms move. Larger values cost memory and neighbour-loop
         !> time but rebuild less often.
         !> @units A
      else if (keyword == 'neighbors_buffer') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%neighbors_buffer
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("neighbors_buffer", params%neighbors_buffer)
         !> @kw radii
         !> Per-species radius used by the Monte-Carlo insertion and accessible-volume tests, one
         !> value per entry in species. Not a physical parameter of the potential.
         !> @units A
         !> @see accessible_volume, mc_types
      else if (keyword == 'radii') then
         backspace (unit)
         call read_parameters(unit, iostatus, n_species, params%radii)
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameters("radii", params%radii)
         !> @kw random_seed
         !> Seed for the intrinsic pseudo-random number generator. Zero, the default, leaves the
         !> compiler's own sequence alone; any other value makes a run repeatable, which is what the
         !> regression comparisons rely on when the initial velocities are randomized.
         !> @see randomize_velocities
      else if (keyword == 'random_seed') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%random_seed_value
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("random_seed", params%random_seed_value)
         !> @kw species
         !> Chemical symbols of the species in the system, in the order every other per-species list
         !> is written in. n_species entries.
      else if (keyword == 'species') then
         backspace (unit)
         call read_parameters(unit, iostatus, n_species, params%species_types)
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameters("species", params%species_types)
         if (iostatus > 0) then
            write (*, *) '                                       |'
            write (*, *) 'ERROR: your "species" keyword is wrong |  <-- ERROR'
            stop
         end if
         !> @kw timing
         !> Print a breakdown of where wall time went, by phase, at the end of the run.
      else if (keyword == 'timing') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_timing
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("timing", params%do_timing)
         !> @kw verbosity, verb
         !> How much diagnostic output to print. 0 is quiet; the MC and MD drivers gate extra
         !> reporting on thresholds up to about 50.
      else if (keyword == 'verbosity' .or. keyword == 'verb') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%verb
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("verbosity", params%verb)
         !> @kw which_atom
         !> Restrict the prediction to a single atom index, for debugging a local environment. 0, the
         !> default, means every atom.
         !> @modes predict
      else if (keyword == 'which_atom') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%which_atom
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("which_atom", params%which_atom)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_general

!**************************************************************************
!  Input keywords for what the run does: prediction, forces, optimisation.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_control(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk

      !> @kw do_derivatives
      !> Compute the derivatives of the SOAP descriptor with respect to atomic positions. Needed
      !> for forces, and switched on automatically in md and mc mode. Turning it off in predict
      !> mode is what makes an energy-only pass cheap.
      if (keyword == 'do_derivatives') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_derivatives
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_derivatives", params%do_derivatives)
         !> @kw do_derivatives_fd
         !> Compute the descriptor derivatives by finite differences instead of analytically. A
         !> correctness check on the analytic route, far slower, and not meant for production.
         !> @see do_derivatives
      else if (keyword == 'do_derivatives_fd') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_derivatives_fd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_derivatives_fd", params%do_derivatives_fd)
         !> @kw do_forces
         !> Compute forces and the virial. Set automatically in md, mc and predict mode; the only
         !> reason to name it is to switch it off.
      else if (keyword == 'do_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_forces", params%do_forces)
         !> @kw do_mc
         !> Run a Monte-Carlo walk. Set automatically by `turbogap mc`.
         !> @modes mc
      else if (keyword == 'do_mc') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_mc
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_mc", params%do_mc)
         !> @kw do_md
         !> Run molecular dynamics. Set automatically by `turbogap md`.
         !> @modes md
      else if (keyword == 'do_md') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_md
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_md", params%do_md)
         !> @kw do_prediction
         !> Evaluate the GAP at all. Set automatically in every mode except soap; switching it off
         !> leaves only the descriptor calculation.
      else if (keyword == 'do_prediction') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_prediction
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_prediction", params%do_prediction)
         !> @kw e_tol
         !> Convergence threshold on the energy change between successive relaxation steps. The
         !> relaxation stops when the energy, force and pressure criteria are all met.
         !> @units eV
         !> @needs optimize
         !> @see f_tol, p_tol
      else if (keyword == 'e_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%e_tol
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("e_tol", params%e_tol)
         !> @kw f_tol
         !> Convergence threshold on the largest force component during a relaxation.
         !> @units eV/A
         !> @needs optimize
         !> @see e_tol, p_tol
      else if (keyword == 'f_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%f_tol
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("f_tol", params%f_tol)
         !> @kw gamma0
         !> The gradient-descent prefactor this once scaled. The step is now taken from max_opt_step
         !> and the largest force, so the value is read and discarded.
         !> @see max_opt_step
         !> @type real
         !> @status ignored
      else if (keyword == "gamma0") then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw gd_box_weight
         !> How heavily the lattice degrees of freedom count against the atomic ones in a gd-box
         !> relaxation. The lattice is optimized in the scaled variable
         !> gd_box_weight*sqrt(n_sites)*A*diag(1/|a0|,1/|b0|,1/|c0|), which is what lets one step
         !> length serve both blocks. Larger values make the cell move more slowly than the atoms.
         !> @needs optimize
         !> @see optimize, max_opt_step
      else if (keyword == "gd_box_weight") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gd_box_weight
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gd_box_weight", params%gd_box_weight)
         !> @kw max_opt_step
         !> Largest distance any atom may move in one gradient-descent step. The step size is chosen
         !> so that the biggest displacement equals this.
         !> @units A
         !> @needs optimize
      else if (keyword == "max_opt_step") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%max_opt_step
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("max_opt_step", params%max_opt_step)
         !> @kw max_opt_step_eps
         !> Largest strain any cell vector could take in one step of the old alternating cell
         !> relaxation. gd-box now relaxes positions and lattice together in one preconditioned
         !> descent with a single step length, taken from max_opt_step, so this is read and
         !> discarded. Use gd_box_weight to change how far the cell moves per step.
         !> @see max_opt_step, gd_box_weight
         !> @type real
         !> @status ignored
      else if (keyword == "max_opt_step_eps") then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw optimize
         !> How the geometry is driven: "vv" for velocity-Verlet dynamics, "gd" for gradient descent
         !> on the positions, "gd-box" to relax the cell as well, and "gd-box-ortho" to relax the cell
         !> keeping it orthorhombic. Anything else aborts the run.
         !> @modes md mc
      else if (keyword == "optimize") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%optimize
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("optimize", params%optimize)
         if (params%optimize == "vv" .or. params%optimize == "gd" .or. params%optimize == "gd-box" .or. &
             params%optimize == "gd-box-ortho") then
            continue
         else
            write (*, *) "ERROR: optimize algorithm not implemented:", params%optimize
            stop
         end if
         !> @kw p_tol
         !> Convergence threshold on the pressure during a cell relaxation.
         !> @units GPa
         !> @needs optimize
         !> @see e_tol, f_tol
      else if (keyword == 'p_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%p_tol
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("p_tol", params%p_tol)
         !> @kw soap_radial_legacy_filter
         !> Keep the pre-2026 seed of the SOAP radial filter recursion. That seed drops the
         !> surface term at rcut_hard, which leaves the radial derivatives disagreeing with a
         !> finite difference of the coefficients by about exp(-nf^2/2) -- 5e-6 relative near
         !> the hard cutoff for the default nf. Setting this to .false. restores the term, which
         !> makes the derivatives finite-difference exact and the expansion slightly faster, at
         !> the cost of moving energies by ~3e-8 relative and forces by ~3e-6 eV/A, so existing
         !> baselines have to be regenerated. Second radial derivatives require .false.
         !> @see do_derivatives
      else if (keyword == 'soap_radial_legacy_filter') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%soap_radial_legacy_filter
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("soap_radial_legacy_filter", params%soap_radial_legacy_filter)
         !> @kw target_pos_step
         !> Displacement the variable time step aims for: the step is rescaled so that the fastest
         !> atom moves about this far. Naming this keyword is what enables the variable time step.
         !> @units A
         !> @modes md
         !> @see tau_dt
      else if (keyword == 'target_pos_step') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%target_pos_step
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("target_pos_step", params%target_pos_step)
         params%variable_time_step = .true.
         !> @kw tau_dt
         !> Relaxation time over which the variable time step is allowed to change, so that the step
         !> size follows target_pos_step smoothly instead of jumping.
         !> @units fs
         !> @modes md
         !> @needs target_pos_step
      else if (keyword == 'tau_dt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%tau_dt
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("tau_dt", params%tau_dt)
         !> @kw write_derivatives
         !> Write the SOAP descriptor derivatives to file alongside the descriptors. Large output; a
         !> debugging aid.
         !> @see write_soap
      else if (keyword == 'write_derivatives') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_derivatives
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_derivatives", params%write_derivatives)
         !> @kw write_soap
         !> Write the SOAP descriptor vectors to file. Set automatically by `turbogap soap`.
      else if (keyword == 'write_soap') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_soap
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_soap", params%write_soap)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_control

!**************************************************************************
!  Input keywords for molecular dynamics: thermostat, barostat, step and duration.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_md(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: i
      integer :: nw
      integer :: iostatus2
      real(dp) :: bsf
      character*1024 :: long_line
      character*128, allocatable :: long_line_items(:)
      logical :: valid_choice

      !> @kw barostat
      !> Pressure coupling: "none" or "berendsen". Anything else aborts the run.
      !> @modes md mc
      !> @see p_beg, p_end, tau_p, barostat_sym
      if (keyword == 'barostat') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%barostat
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("barostat", params%barostat)
         call upper_to_lower_case(params%barostat)
         valid_choice = .false.
         do i = 1, size(implemented_barostats)
            if (trim(params%barostat) == trim(implemented_barostats(i))) then
               valid_choice = .true.
            end if
         end do
         if (.not. valid_choice) then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid barostat keyword:", params%barostat
               write (*, *) "This is a list of valid options:"
               write (*, *) implemented_barostats
            end if
            stop
         end if
         !> @kw barostat_sym
         !> Which components of the cell the barostat is allowed to change: "isotropic" scales all
         !> three axes together, and the anisotropic settings let them move independently.
         !> @modes md mc
         !> @needs barostat
      else if (keyword == 'barostat_sym') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%barostat_sym
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("barostat_sym", params%barostat_sym)
         !> @kw box_scaling_factor
         !> A 3x3 matrix the cell is multiplied by, applied once at the start of the run. Nine
         !> numbers, read column by column. Used to strain a cell without editing the XYZ file.
         !> @needs scale_box
      else if (keyword == 'box_scaling_factor') then
         backspace (unit)
         read (10, '(A)', iostat=iostatus) long_line
         allocate (long_line_items(1:9))
         do i = 1, 9
            read (long_line, *, iostat=iostatus2) cjunk, cjunk, long_line_items(1:i)
            if (iostatus2 == -1) exit
         end do
         i = i - 1
         if (i == 1) then
            read (long_line_items(1), *) bsf
            params%box_scaling_factor(1, 1) = bsf
            params%box_scaling_factor(2, 2) = bsf
            params%box_scaling_factor(3, 3) = bsf
         else if (i == 3) then
            read (long_line_items(1), *) bsf
            params%box_scaling_factor(1, 1) = bsf
            read (long_line_items(2), *) bsf
            params%box_scaling_factor(2, 2) = bsf
            read (long_line_items(3), *) bsf
            params%box_scaling_factor(3, 3) = bsf
         else if (i == 9) then
            read (long_line_items(1), *) bsf
            params%box_scaling_factor(1, 1) = bsf
            read (long_line_items(2), *) bsf
            params%box_scaling_factor(1, 2) = bsf
            read (long_line_items(3), *) bsf
            params%box_scaling_factor(1, 3) = bsf
            read (long_line_items(4), *) bsf
            params%box_scaling_factor(2, 1) = bsf
            read (long_line_items(5), *) bsf
            params%box_scaling_factor(2, 2) = bsf
            read (long_line_items(6), *) bsf
            params%box_scaling_factor(2, 3) = bsf
            read (long_line_items(7), *) bsf
            params%box_scaling_factor(3, 1) = bsf
            read (long_line_items(8), *) bsf
            params%box_scaling_factor(3, 2) = bsf
            read (long_line_items(9), *) bsf
            params%box_scaling_factor(3, 3) = bsf
         else
            write (*, *) "ERROR: the box_scaling_factor must be given by 1, 3 or 9 numbers"
            stop
         end if
         deallocate (long_line_items)
         if (rank == 0) call print_parameters("box_scaling_factor", reshape(params%box_scaling_factor, [9]))
         !> @kw gamma_p
         !> Damping of the Berendsen barostat's cell response, on top of tau_p.
         !> @modes md mc
         !> @needs barostat
      else if (keyword == 'gamma_p') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gamma_p
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gamma_p", params%gamma_p)
         !> @kw md_nsteps
         !> Number of molecular-dynamics steps to take.
         !> @modes md
      else if (keyword == 'md_nsteps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%md_nsteps
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("md_nsteps", params%md_nsteps)
         !> @kw md_step
         !> Time step. With a variable time step this is the starting value.
         !> @units fs
         !> @modes md
         !> @see target_pos_step
      else if (keyword == 'md_step') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%md_step
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("md_step", params%md_step)
         !> @kw n_t_hold
         !> Number of entries in the t_hold list, which must be given before it.
         !> @modes md
         !> @see t_hold
      else if (keyword == 'n_t_hold') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_t_hold
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_t_hold", params%n_t_hold)
         allocate (params%t_hold(1:params%n_t_hold*3))
         !> @kw p_beg
         !> Target pressure at the start of the run. With p_end it defines a linear ramp over the run.
         !> @units GPa
         !> @modes md mc
         !> @needs barostat
      else if (keyword == 'p_beg') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%p_beg
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("p_beg", params%p_beg)
         !> @kw p_end
         !> Target pressure at the end of the run.
         !> @units GPa
         !> @modes md mc
         !> @needs barostat
         !> @see p_beg
      else if (keyword == 'p_end') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%p_end
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("p_end", params%p_end)
         !> @kw randomize_velocities
         !> Draw fresh velocities at t_beg instead of taking them from the XYZ file. Also what the
         !> "md" Monte Carlo move does before each burst. velocity_distribution chooses the draw.
         !> Use random_seed to make it repeatable.
         !> @modes md mc
         !> @see velocity_distribution, random_seed, t_beg
      else if (keyword == 'randomize_velocities') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%randomize_velocities
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("randomize_velocities", params%randomize_velocities)
         !> @kw velocity_distribution
         !> Which distribution randomize_velocities draws from: "maxwell" for Maxwell-Boltzmann at
         !> t_beg, or "uniform" for the historical draw, whose components are uniform on [0,1) and
         !> are then scaled so the kinetic energy is exactly (3/2)(N-1)kT. Only "maxwell" samples
         !> the canonical ensemble, and mc_hamiltonian needs it for detailed balance; "uniform" is
         !> the default so that existing trajectories reproduce. Anything else aborts the run.
         !> @modes md mc
         !> @see randomize_velocities, mc_hamiltonian, t_beg
      else if (keyword == 'velocity_distribution') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%velocity_distribution
         call check_iostatus(iostatus, keyword)
         call upper_to_lower_case(params%velocity_distribution)
         if (rank == 0) call print_parameter("velocity_distribution", params%velocity_distribution)
         if (params%velocity_distribution /= "maxwell" .and. params%velocity_distribution /= "uniform") then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid velocity_distribution keyword:", params%velocity_distribution
               write (*, *) "This is a list of valid options:"
               write (*, *) "maxwell uniform"
            end if
            call turbogap_abort()
         end if
         !> @kw scale_box
         !> Apply box_scaling_factor to the cell at the start of the run.
         !> @see box_scaling_factor
      else if (keyword == 'scale_box') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%scale_box
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("scale_box", params%scale_box)
         !> @kw t_beg
         !> Target temperature at the start of the run. With t_end it defines a linear ramp over the
         !> run; give only t_beg for a constant-temperature run.
         !> @units K
         !> @modes md mc
         !> @see thermostat, t_end, t_hold
      else if (keyword == 't_beg') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%t_beg
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("t_beg", params%t_beg)
         !> @kw t_end
         !> Target temperature at the end of the run.
         !> @units K
         !> @modes md mc
         !> @see t_beg
      else if (keyword == 't_end') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%t_end
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("t_end", params%t_end)
         !> @kw t_hold
         !> Segments of the temperature schedule, as a flat list read n_t_hold entries at a time: a
         !> piecewise ramp-and-hold instead of the single linear ramp that t_beg and t_end give.
         !> n_t_hold must appear first.
         !> @units K
         !> @modes md
         !> @needs n_t_hold
      else if (keyword == 't_hold') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%t_hold(nw), nw=1, params%n_t_hold*3)
         if (rank == 0) call print_parameters("t_hold", params%t_hold)
         !> @kw tau_p
         !> Relaxation time of the barostat: how quickly the cell is pushed towards the target
         !> pressure. Larger is gentler.
         !> @units fs
         !> @modes md mc
         !> @needs barostat
      else if (keyword == 'tau_p') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%tau_p
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("tau_p", params%tau_p)
         !> @kw tau_t
         !> Relaxation time of the thermostat: how quickly the kinetic energy is pushed towards the
         !> target temperature. Larger is gentler.
         !> @units fs
         !> @modes md mc
         !> @needs thermostat
      else if (keyword == 'tau_t') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%tau_t
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("tau_t", params%tau_t)
         !> @kw thermostat
         !> Temperature coupling: "none", "berendsen" or "bussi". Bussi is the stochastic velocity
         !> rescaling that samples the canonical ensemble properly; Berendsen does not. Anything else
         !> aborts the run.
         !> @modes md mc
         !> @see t_beg, t_end, tau_t
      else if (keyword == 'thermostat') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%thermostat
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("thermostat", params%thermostat)
         call upper_to_lower_case(params%thermostat)
         valid_choice = .false.
         do i = 1, size(implemented_thermostats)
            if (trim(params%thermostat) == trim(implemented_thermostats(i))) then
               valid_choice = .true.
            end if
         end do
         if (.not. valid_choice) then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid thermostat keyword:", params%thermostat
               write (*, *) "This is a list of valid options:"
               write (*, *) implemented_thermostats
            end if
            stop
         end if
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_md

!**************************************************************************
!  Input keywords for nested sampling.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_nested(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk

      !> @kw n_nested
      !> Number of nested-sampling iterations. Naming it is what enables nested sampling.
      if (keyword == 'n_nested') then
         if (mode /= "predict") then
            write (*, *) 'ERROR: the "n_nested" option for nested sampling can only be used with "turbogap predict"'
            stop
         end if
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_nested
         if (rank == 0) call print_parameter("n_nested", params%n_nested)
         if (params%n_nested > 0) then
            params%do_nested_sampling = .true.
         end if
         !> @kw nested_max_strain
         !> Largest strain a nested-sampling cell-shape move may apply.
      else if (keyword == 'nested_max_strain') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%nested_max_strain
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nested_max_strain", params%nested_max_strain)
         !> @kw nested_max_volume_change
         !> Largest fractional volume change a nested-sampling volume move may make.
         !> @needs scale_box_nested
      else if (keyword == 'nested_max_volume_change') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%nested_max_volume_change
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nested_max_volume_change", params%nested_max_volume_change)
         !> @kw p_nested
         !> External pressure entering the nested-sampling enthalpy.
         !> @units GPa
         !> @needs n_nested
      else if (keyword == 'p_nested') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%p_nested
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("p_nested", params%p_nested)
         !> @kw scale_box_nested
         !> Let nested sampling change the cell as well as the positions.
         !> @needs n_nested
         !> @see nested_max_volume_change, nested_max_strain
      else if (keyword == 'scale_box_nested') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%scale_box_nested
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("scale_box_nested", params%scale_box_nested)
         !> @kw t_extra
         !> Temperature offset from an earlier nested-sampling design; nothing consults it.
         !> @type real
         !> @status ignored
      else if (keyword == 't_extra') then
         call ignored_keyword(unit, iostatus, rank, keyword)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_nested

!**************************************************************************
!  Input keywords for Monte Carlo: move types, limits and acceptance.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_mc(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: i
      integer :: j
!     Normalisation accumulator for the acceptance weight lists. It was an
!     integer, so every partial sum of the real weights was truncated:
!     "mc_acceptance = 0.3 0.7" summed to zero and the normalisation divided
!     by it. Integer weights survived only because their partial sums happen
!     to be exact.
      real(dp) :: k
      integer :: nw
      logical :: valid_choice

      !> @kw accessible_volume
      !> Use the volume actually accessible to a new atom, with the volume of the existing atoms'
      !> spheres removed, rather than the cell volume, in the insertion and removal acceptance
      !> ratios. The spheres are sized by radii.
      !> @modes mc
      !> @see radii
      if (keyword == 'accessible_volume') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%accessible_volume
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("accessible_volume", params%accessible_volume)

         !> @kw mc_acceptance
         !> Relative probability of picking each move type, one weight per entry in mc_types and in
         !> the same order. Normalised internally, so the weights need not sum to one. Its length sets
         !> n_mc_types.
         !> @modes mc
         !> @see mc_types
      else if (keyword == 'mc_acceptance') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_acceptance(nw), nw=1, params%n_mc_types)
         ! The acceptance probability is based on this sum and normalised
         k = 0.d0
         do i = 1, params%n_mc_types
            k = k + params%mc_acceptance(i)
         end do

         do i = 1, params%n_mc_types
            params%mc_acceptance(i) = params%mc_acceptance(i)/k
         end do

         if (rank == 0) call print_parameters("mc_acceptance", params%mc_acceptance)
         !> @kw mc_hamiltonian
         !> Include the kinetic energy in the Metropolis test, so that the walk samples configurations
         !> and momenta together rather than configurations alone. Also what makes the "md" move type
         !> meaningful.
         !> @modes mc
         !> @see mc_types
      else if (keyword == 'mc_hamiltonian') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_hamiltonian
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_hamiltonian", params%mc_hamiltonian)
         !> @kw mc_hybrid_opt
         !> How the short MD burst of a hybrid "md" move is propagated: same choices as optimize.
         !> @modes mc
         !> @needs mc_types
         !> @see optimize
      else if (keyword == 'mc_hybrid_opt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_hybrid_opt
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_hybrid_opt", params%mc_hybrid_opt)
         !> @kw mc_lnvol_max
         !> Largest change in ln(V) a volume move may make, so the trial volume is drawn from
         !> V*exp(+-mc_lnvol_max). Working in the logarithm keeps the move symmetric in volume.
         !> @modes mc
         !> @needs mc_types
      else if (keyword == 'mc_lnvol_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_lnvol_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_lnvol_max", params%mc_lnvol_max)
         !> @kw mc_max_dist
         !> Largest distance from an existing atom at which a new atom may be inserted. Together with
         !> mc_min_dist it bounds where an insertion trial may land.
         !> @units A
         !> @modes mc
         !> @see mc_min_dist, mc_max_insertion_trials
      else if (keyword == 'mc_max_dist') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_max_dist
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_max_dist", params%mc_max_dist)
         !> @kw mc_max_dist_to_planes
         !> Distance from each plane within which a moved or inserted atom must stay, one value per
         !> plane. mc_n_planes must appear first.
         !> @units A
         !> @modes mc
         !> @needs mc_n_planes
         !> @see mc_planes
      else if (keyword == 'mc_max_dist_to_planes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_max_dist_to_planes(nw), nw=1, params%mc_n_planes)

         if (rank == 0) call print_parameters("mc_max_dist_to_planes", params%mc_max_dist_to_planes, "A")
         !> @kw mc_max_insertion_trials
         !> How many random positions an insertion move may try before giving up. The run reports the
         !> failure and says which of mc_min_dist, mc_max_dist or this to change.
         !> @modes mc
         !> @see mc_min_dist, mc_max_dist
      else if (keyword == 'mc_max_insertion_trials') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_max_insertion_trials
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_max_insertion_trials", params%mc_max_insertion_trials)
         !> @kw mc_min_dist
         !> Closest an inserted atom may come to an existing one. Trials nearer than this are rejected
         !> without evaluating the potential, which is what keeps an insertion from landing on top of
         !> an atom.
         !> @units A
         !> @modes mc
         !> @see mc_max_dist
      else if (keyword == 'mc_min_dist') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_min_dist
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_min_dist", params%mc_min_dist)
         !> @kw mc_move_max
         !> Largest displacement of a single-atom "move". The trial displacement is drawn uniformly up
         !> to this.
         !> @units A
         !> @modes mc
         !> @needs mc_types
      else if (keyword == 'mc_move_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_move_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_move_max", params%mc_move_max)
         !> @kw mc_mu
         !> Chemical potential of each species that may be inserted or removed, one value per entry in
         !> mc_species. Its length sets n_mc_mu.
         !> @units eV
         !> @modes mc
         !> @see mc_species, mc_types
      else if (keyword == 'mc_mu') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_mu(nw), nw=1, params%n_mc_mu)
         if (rank == 0) call print_parameters("mc_mu", params%mc_mu, "eV")
         !> @kw mc_molecule_files
         !> Insert and remove a whole rigid molecule instead of a single atom, one xyz file per
         !> entry in mc_species, or "none" for an entry that really is an atom. n_mc_mu must appear
         !> first. Every species in the file has to be in species. The molecule is placed with a
         !> uniformly random orientation and its centre of mass drawn like an atomic insertion, and
         !> a removal takes all of its atoms together, so mc_mu is the chemical potential of the
         !> molecule.
         !> @modes mc
         !> @needs n_mc_mu
         !> @see mc_species, mc_mu, mc_min_dist
      else if (keyword == 'mc_molecule_files') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_molecule_files(nw), nw=1, params%n_mc_mu)
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameters("mc_molecule_files", params%mc_molecule_files)
         !> @kw mc_mu_reference
         !> What mc_mu is measured from. "absolute" compares it against the total energy change of
         !> the exchange, which carries the e0 of whatever was inserted or removed. "e0" adds that
         !> e0 back into the acceptance ratio, so mc_mu is quoted relative to the isolated-species
         !> reference energy instead. For a molecule the e0 of all its atoms is summed. Anything
         !> else aborts the run.
         !> @modes mc
         !> @see mc_mu, e0
      else if (keyword == 'mc_mu_reference') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_mu_reference
         call check_iostatus(iostatus, keyword)
         call upper_to_lower_case(params%mc_mu_reference)
         if (rank == 0) call print_parameter("mc_mu_reference", params%mc_mu_reference)
         if (params%mc_mu_reference /= "absolute" .and. params%mc_mu_reference /= "e0") then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid mc_mu_reference keyword:", params%mc_mu_reference
               write (*, *) "This is a list of valid options:"
               write (*, *) "absolute e0"
            end if
            call turbogap_abort()
         end if
         !> @kw mc_mu_acceptance
         !> Relative probability of picking each species for a grand-canonical move, one weight per
         !> entry in mc_species.
         !> @modes mc
         !> @see mc_mu, mc_species
      else if (keyword == 'mc_mu_acceptance') then
         backspace (unit)
!        One weight per chemical potential, not per move type: mc_mu_acceptance
!        is allocated 1:n_mc_mu, so reading n_mc_types of them wrote past the
!        end whenever a deck had more move types than species.
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_mu_acceptance(nw), nw=1, params%n_mc_mu)
         ! The acceptance probability is based on this sum and normalised
         k = 0.d0
         do i = 1, params%n_mc_mu
            k = k + params%mc_mu_acceptance(i)
         end do

         do i = 1, params%n_mc_mu
            params%mc_mu_acceptance(i) = params%mc_mu_acceptance(i)/k
         end do

         if (rank == 0) call print_parameters("mc_mu_acceptance", params%mc_mu_acceptance)
         !> @kw mc_n_planes
         !> Number of planes bounding the region moves are restricted to. Must appear before mc_planes
         !> and mc_max_dist_to_planes, which it allocates.
         !> @modes mc
         !> @see mc_planes
      else if (keyword == 'mc_n_planes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_n_planes
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_n_planes", params%mc_n_planes)

         ! Allocate 4 for a b c d, in ax + by + cz = d
         if (params%mc_n_planes > 0) then
            allocate (params%mc_planes(4*params%mc_n_planes))
            allocate (params%mc_max_dist_to_planes(params%mc_n_planes))
         end if

         !> @kw mc_nrelax
         !> Number of relaxation steps to take after an accepted move when mc_relax is on.
         !> @modes mc
         !> @needs mc_relax
      else if (keyword == 'mc_nrelax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_nrelax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_nrelax", params%mc_nrelax)
         !> @kw mc_nsteps
         !> Number of Monte-Carlo steps to take.
         !> @modes mc
      else if (keyword == 'mc_nsteps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_nsteps
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_nsteps", params%mc_nsteps)
         !> @kw mc_optimize_exp
         !> Declares that this is a reverse Monte-Carlo run, so that a missing exp_energy_scales is
         !> reported. It does not itself put the experimental penalty into the Metropolis test:
         !> the acceptance ratio reads the total energy, and what folds the penalty into that is
         !> exp_energies. Set both.
         !> @modes mc
         !> @needs do_exp
         !> @see exp_energies, exp_energy_scales
      else if (keyword == 'mc_optimize_exp') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_optimize_exp
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_optimize_exp", params%mc_optimize_exp)
         !> @kw mc_planes
         !> The planes bounding the region moves are restricted to, four numbers per plane giving its
         !> normal and offset. mc_n_planes must appear first.
         !> @modes mc
         !> @needs mc_n_planes
         !> @see mc_max_dist_to_planes, mc_planes_restrict_to_polyhedron
      else if (keyword == 'mc_planes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_planes(nw), nw=1, 4*params%mc_n_planes)

         if (rank == 0) call print_parameters("mc_planes", params%mc_planes)
         !> @kw mc_planes_restrict_to_polyhedron
         !> Require an atom to be inside the polyhedron the planes enclose, rather than merely within
         !> mc_max_dist_to_planes of each plane.
         !> @modes mc
         !> @needs mc_planes
      else if (keyword == 'mc_planes_restrict_to_polyhedron') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_planes_restrict_to_polyhedron
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_planes_restrict_to_polyhedron", params%mc_planes_restrict_to_polyhedron)

         !> @kw mc_relax
         !> Relax the geometry after each accepted move, so the walk samples relaxed configurations.
         !> @modes mc
         !> @see mc_nrelax, mc_relax_opt, mc_relax_after
      else if (keyword == 'mc_relax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_relax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_relax", params%mc_relax)
         !> @kw mc_relax_after
         !> Which move types trigger a relaxation, by name. n_mc_relax_after must appear first.
         !> Without it, mc_relax applies to every accepted move.
         !> @modes mc
         !> @needs n_mc_relax_after
         !> @see mc_relax
      else if (keyword == 'mc_relax_after') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_relax_after(nw), nw=1, params%n_mc_relax_after)
         if (rank == 0) call print_parameters("mc_relax_after", params%mc_relax_after)
         !> @kw mc_relax_opt
         !> How the post-move relaxation is driven: same choices as optimize.
         !> @modes mc
         !> @needs mc_relax
         !> @see optimize
      else if (keyword == 'mc_relax_opt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_relax_opt
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_relax_opt", params%mc_relax_opt)
         !> @kw mc_reverse
         !> Reverse-Monte-Carlo switch from an earlier design; nothing consults it. Use
         !> mc_optimize_exp.
         !> @see mc_optimize_exp
         !> @type logical
         !> @modes mc
         !> @status ignored
      else if (keyword == 'mc_reverse') then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw mc_reverse_lambda
         !> Mixing weight for the same earlier design; nothing consults it. The weight that does exist
         !> is exp_energy_scales.
         !> @see exp_energy_scales
         !> @type real
         !> @modes mc
         !> @status ignored
      else if (keyword == 'mc_reverse_lambda') then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw mc_species
         !> Species that may be inserted or removed by grand-canonical moves. Its length sets n_mc_mu.
         !> @modes mc
         !> @see mc_mu, mc_types
      else if (keyword == 'mc_species') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_species(nw), nw=1, params%n_mc_mu)
         if (rank == 0) call print_parameters("mc_species", params%mc_species)
         !> @kw mc_swaps
         !> Species pairs that a "swap" move may exchange, given as a flat list of symbols read two at
         !> a time. Its length sets n_mc_swaps, and the symbols are resolved against species.
         !> @modes mc
         !> @see mc_types, species
      else if (keyword == 'mc_swaps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_swaps(nw), nw=1, 2*params%n_mc_swaps)
         !       Need the check the implemented types
         valid_choice = .false.
         do j = 1, 2*params%n_mc_swaps
            valid_choice = .false.
            do i = 1, n_species
               if (trim(params%species_types(i)) == trim(params%mc_swaps(j))) then
!                 Indexed by swap slot, not by species: count_swap_species reads
!                 mc_swaps_id(1:2*n_mc_swaps) as "the species in slot j". Writing
!                 it at the species index went out of bounds as soon as a swap
!                 named a species whose index exceeds 2*n_mc_swaps, and left the
!                 slot it should have filled undefined.
                  params%mc_swaps_id(j) = i
                  valid_choice = .true.
               end if
            end do
            if (.not. valid_choice) then
               if (rank == 0) then
                  write (*, *) "ERROR -> Invalid mc_swaps species keyword:", params%mc_swaps(j)
                  write (*, *) "This is a list of valid options:"
                  write (*, *) params%species_types
               end if
               stop
            end if
         end do

         if (rank == 0) call print_parameters("mc_swaps", params%mc_swaps)
         !> @kw mc_types
         !> Which moves the walk may make, from "move", "insertion", "removal", "relax", "md", "swap",
         !> "volume" and "none". Any other name aborts the run. Its length sets n_mc_types, and
         !> mc_acceptance weights them.
         !> @modes mc
         !> @see mc_acceptance
      else if (keyword == 'mc_types') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%mc_types(nw), nw=1, params%n_mc_types)
         !       Need the check the implemented types
         valid_choice = .false.
         do j = 1, params%n_mc_types
            call upper_to_lower_case(params%mc_types(j))
            valid_choice = .false.
            do i = 1, size(implemented_mc_types)
               if (trim(params%mc_types(j)) == trim(implemented_mc_types(i))) then
                  valid_choice = .true.
               end if
            end do
            if (.not. valid_choice) then
               if (rank == 0) then
                  write (*, *) "ERROR -> Invalid mc_type keyword:", params%mc_types(j)
                  write (*, *) "This is a list of valid options:"
                  write (*, *) implemented_mc_types
               end if
               stop
            end if
         end do
         if (rank == 0) call print_parameters("mc_types", params%mc_types)
         !> @kw mc_write_xyz
         !> Write the configuration after every Monte-Carlo step, not only the accepted ones.
         !> @modes mc
         !> @see write_xyz
      else if (keyword == 'mc_write_xyz') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mc_write_xyz
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mc_write_xyz", params%mc_write_xyz)
         !> @kw n_mc_mu
         !> Number of species in the grand-canonical lists. Give it before mc_species, mc_mu and
         !> mc_mu_acceptance, which it allocates; giving those lists instead sets it implicitly.
         !> @modes mc
         !> @see mc_species
      else if (keyword == 'n_mc_mu') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_mc_mu
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_mc_mu", params%n_mc_mu)
         allocate (params%mc_mu(1:params%n_mc_mu))
         allocate (params%mc_species(1:params%n_mc_mu))
         allocate (params%mc_mu_acceptance(1:params%n_mc_mu))
         params%mc_mu_acceptance = 1.d0/dfloat(params%n_mc_mu)
         allocate (params%mc_molecule_files(1:params%n_mc_mu))
         params%mc_molecule_files = "none"
         allocate (params%mc_molecules(1:params%n_mc_mu))
         allocate (params%mc_exchange_mass(1:params%n_mc_mu))
         allocate (params%mc_exchange_e0(1:params%n_mc_mu))
         params%mc_exchange_mass = 0.d0
         params%mc_exchange_e0 = 0.d0

         !> @kw n_mc_relax_after
         !> Number of entries in mc_relax_after, which it allocates and must precede.
         !> @modes mc
         !> @see mc_relax_after
      else if (keyword == 'n_mc_relax_after') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_mc_relax_after
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_mc_relax_after", params%n_mc_relax_after)
         allocate (params%mc_relax_after(1:params%n_mc_relax_after))
         !> @kw n_mc_swaps
         !> Number of swap pairs, which is half the length of mc_swaps. Give it before mc_swaps, which
         !> it allocates.
         !> @modes mc
         !> @see mc_swaps
      else if (keyword == 'n_mc_swaps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_mc_swaps
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_mc_swaps", params%n_mc_swaps)
         allocate (params%mc_swaps(1:2*params%n_mc_swaps))
         allocate (params%mc_swaps_id(1:2*params%n_mc_swaps))
         !> @kw n_mc_types
         !> Number of move types. Give it before mc_types and mc_acceptance, which it allocates;
         !> giving those lists instead sets it implicitly.
         !> @modes mc
         !> @see mc_types
      else if (keyword == 'n_mc_types') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_mc_types
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_mc_types", params%n_mc_types)
         allocate (params%mc_types(1:params%n_mc_types))
         allocate (params%mc_acceptance(1:params%n_mc_types))
         params%mc_acceptance = 1.d0/dfloat(params%n_mc_types)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_mc

!**************************************************************************
!  Input keywords for van der Waals: TS, MBD, and the Hirshfeld model.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_vdw(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found, are_vdw_refs_read)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: are_vdw_refs_read(1:3)

      character*64 :: cjunk

      !> @kw do_nnls
      !> Solve the Hirshfeld-volume fit with non-negative least squares, which keeps the fitted
      !> volumes positive where an unconstrained solve would not.
      !> @needs vdw_type
      if (keyword == "do_nnls") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_nnls
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_nnls", params%do_nnls)
         !> @kw mbd_correction_freq
         !> Recompute the many-body dispersion every this many MD steps, reusing the previous
         !> correction in between. MBD is the expensive part of ts+mbd, so this trades accuracy for
         !> time.
         !> @modes md mc
         !> @needs vdw_type
      else if (keyword == 'mbd_correction_freq') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%mbd_correction_freq
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("mbd_correction_freq", params%mbd_correction_freq)
! NEW VDW STUFF ENDS HERE
         !> @kw poly_cut_xmax
         !> Upper end of the range over which the polynomial cutoff turns the dispersion off.
         !> @units A
         !> @needs vdw_polynomial
         !> @see poly_cut_xmin
      else if (keyword == "poly_cut_xmax") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%poly_cut_xmax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("poly_cut_xmax", params%poly_cut_xmax)

         !> @kw poly_cut_xmin
         !> Lower end of that range: below it the interaction is at full strength.
         !> @units A
         !> @needs vdw_polynomial
         !> @see poly_cut_xmax
      else if (keyword == "poly_cut_xmin") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%poly_cut_xmin
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("poly_cut_xmin", params%poly_cut_xmin)
         !> @kw print_estat_forces
         !> Write the electrostatic contribution to the forces separately, for checking how much of
         !> the total it is.
         !> @see estat_method
      else if (keyword == 'print_estat_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%print_estat_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("print_estat_forces", params%print_estat_forces)
         !> @kw print_vdw_forces
         !> Write the dispersion contribution to the forces separately.
         !> @see vdw_type
      else if (keyword == 'print_vdw_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%print_vdw_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("print_vdw_forces", params%print_vdw_forces)
         !> @kw vdw_2b_rcut
         !> Cutoff of the pairwise part of the many-body dispersion.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_2b_rcut2
      else if (keyword == "vdw_2b_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_2b_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_2b_rcut", params%vdw_2b_rcut)
         !> @kw vdw_2b_rcut2
         !> Inner radius at which that pairwise part starts to be turned off; it is at full strength
         !> inside this and zero beyond vdw_2b_rcut.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_2b_rcut
      else if (keyword == "vdw_2b_rcut2") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_2b_rcut2
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_2b_rcut2", params%vdw_2b_rcut2)
         !> @kw vdw_alpha0_ref
         !> Free-atom static polarizability of each species, one value per entry in species. Together
         !> with vdw_c6_ref and vdw_r0_ref these are the reference data the Hirshfeld rescaling starts
         !> from.
         !> @units bohr^3
         !> @needs vdw_type
         !> @see vdw_c6_ref, vdw_r0_ref
      else if (keyword == "vdw_alpha0_ref") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_alpha0_ref(1:n_species)
         are_vdw_refs_read(3) = .true.
! NEW VDW STUFF HERE
         if (rank == 0) call print_parameters("vdw_alpha0_ref", params%vdw_alpha0_ref)
         !> @kw vdw_buffer
         !> Width of the smooth turn-off applied at vdw_rcut, so the dispersion energy and its
         !> gradient reach zero together.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_rcut
      else if (keyword == "vdw_buffer") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_buffer
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_buffer", params%vdw_buffer)
         !> @kw vdw_buffer_inner
         !> Width of the smooth turn-on applied at vdw_rcut_inner.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_rcut_inner
      else if (keyword == "vdw_buffer_inner") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_buffer_inner
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_buffer_inner", params%vdw_buffer_inner)
         !> @kw vdw_c6_ref
         !> Free-atom C6 coefficient of each species, one value per entry in species.
         !> @units Ha bohr^6
         !> @needs vdw_type
         !> @see vdw_alpha0_ref, vdw_r0_ref
      else if (keyword == "vdw_c6_ref") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_c6_ref(1:n_species)
         are_vdw_refs_read(1) = .true.
         if (rank == 0) call print_parameters("vdw_c6_ref", params%vdw_c6_ref)
         !> @kw vdw_d
         !> Steepness of the Fermi damping function in the Tkatchenko-Scheffler energy. Larger makes
         !> the switch at vdw_sr*R0 sharper.
         !> @needs vdw_type
         !> @see vdw_sr
      else if (keyword == "vdw_d") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_d
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_d", params%vdw_d)
         !> @kw vdw_d_mbd
         !> The same steepness for the many-body dispersion.
         !> @needs vdw_type
         !> @see vdw_sr_mbd
      else if (keyword == "vdw_d_mbd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_d_mbd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_d_mbd", params%vdw_d_mbd)
         !> @kw vdw_hirsh_grad
         !> Include the derivative of the Hirshfeld volumes with respect to positions in the
         !> dispersion forces. Switching it off gives a cheaper but inconsistent force.
         !> @needs vdw_type
      else if (keyword == "vdw_hirsh_grad") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_hirsh_grad
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_hirsh_grad", params%vdw_hirsh_grad)
         !> @kw vdw_loc_rcut
         !> Radius of the local environment whose atoms enter the screened polarizability of a given
         !> atom.
         !> @units A
         !> @needs vdw_type
      else if (keyword == "vdw_loc_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_loc_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_loc_rcut", params%vdw_loc_rcut)
         !> @kw vdw_mbd_cent_appr
         !> Evaluate the many-body dispersion about the central atom only, rather than over every atom
         !> of its environment. The cheaper approximation.
         !> @needs vdw_type
      else if (keyword == "vdw_mbd_cent_appr") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_cent_appr
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_cent_appr", params%vdw_mbd_cent_appr)
         !> @kw vdw_mbd_grad
         !> Compute the gradient of the many-body dispersion energy. Off gives MBD energies with
         !> TS-only forces.
         !> @needs vdw_type
      else if (keyword == "vdw_mbd_grad") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_grad
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_grad", params%vdw_mbd_grad)
         !> @kw vdw_mbd_nfreq
         !> Number of imaginary-frequency quadrature points used for the dynamic polarizability. More
         !> is more accurate and proportionally slower.
         !> @needs vdw_type
      else if (keyword == "vdw_mbd_nfreq") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_nfreq
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_nfreq", params%vdw_mbd_nfreq)
         !> @kw vdw_mbd_norder
         !> Order at which the many-body expansion is truncated.
         !> @needs vdw_type
      else if (keyword == "vdw_mbd_norder") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_norder
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_norder", params%vdw_mbd_norder)
         !> @kw vdw_mbd_rcut
         !> Cutoff of the many-body dispersion term.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_mbd_rcut2
      else if (keyword == "vdw_mbd_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_rcut", params%vdw_mbd_rcut)
         !> @kw vdw_mbd_rcut2
         !> Inner radius at which the many-body term starts to be turned off.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_mbd_rcut
      else if (keyword == "vdw_mbd_rcut2") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_rcut2
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_mbd_rcut2", params%vdw_mbd_rcut2)
         !> @kw vdw_omega_ref
         !> Reference characteristic excitation frequency of the oscillator model.
         !> @units Ha
         !> @needs vdw_type
      else if (keyword == "vdw_omega_ref") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_omega_ref
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_omega_ref", params%vdw_omega_ref)
         !> @kw vdw_polynomial
         !> Use a polynomial cutoff for the dispersion rather than the default form.
         !> @needs vdw_type
         !> @see poly_cut_xmin, poly_cut_xmax
      else if (keyword == "vdw_polynomial") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_polynomial
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_polynomial", params%vdw_polynomial)
         !> @kw vdw_r0_ref
         !> Free-atom van der Waals radius of each species, one value per entry in species.
         !> @units bohr
         !> @needs vdw_type
         !> @see vdw_c6_ref, vdw_alpha0_ref
      else if (keyword == "vdw_r0_ref") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_r0_ref(1:n_species)
         are_vdw_refs_read(2) = .true.
         if (rank == 0) call print_parameters("vdw_r0_ref", params%vdw_r0_ref)
         !> @kw vdw_rcut
         !> Cutoff of the pairwise Tkatchenko-Scheffler dispersion.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_buffer
      else if (keyword == "vdw_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_rcut", params%vdw_rcut)
         !> @kw vdw_rcut_inner
         !> Radius below which the dispersion is smoothly switched off, so that it does not
         !> double-count the short-range interaction the GAP already describes.
         !> @units A
         !> @needs vdw_type
         !> @see vdw_buffer_inner
      else if (keyword == "vdw_rcut_inner") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_rcut_inner
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_rcut_inner", params%vdw_rcut_inner)
         !> @kw vdw_scs_rcut
         !> Cutoff of the self-consistent screening that turns free-atom polarizabilities into
         !> screened ones.
         !> @units A
         !> @needs vdw_type
      else if (keyword == "vdw_scs_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_scs_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_scs_rcut", params%vdw_scs_rcut)
         !> @kw vdw_sr
         !> Position of the Fermi damping function in the Tkatchenko-Scheffler energy, as a fraction
         !> of the sum of the two van der Waals radii. Fitted per exchange-correlation functional.
         !> @needs vdw_type
         !> @see vdw_d
      else if (keyword == "vdw_sr") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_sr
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_sr", params%vdw_sr)
         !> @kw vdw_sr_mbd
         !> The same fraction for the many-body dispersion.
         !> @needs vdw_type
         !> @see vdw_d_mbd
      else if (keyword == "vdw_sr_mbd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_sr_mbd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_sr_mbd", params%vdw_sr_mbd)
         !> @kw vdw_type
         !> Which dispersion correction to add: "none", "ts" for Tkatchenko-Scheffler, "mbd" for
         !> many-body dispersion, or "ts+mbd" for both. Any other value aborts the run. The correction
         !> needs Hirshfeld volumes, so the GAP has to provide them.
      else if (keyword == "vdw_type") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%vdw_type
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("vdw_type", params%vdw_type)
         call upper_to_lower_case(params%vdw_type)
         if (params%vdw_type == "ts") then
            continue
         else if (params%vdw_type == "none") then
            continue
         else if (params%vdw_type == "mbd") then
            continue
         else if (params%vdw_type == "ts+mbd") then
            continue
         else
            write (*, *) "ERROR: I do not recognize the vdw_type keyword ", params%vdw_type
            stop
         end if
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_vdw

!**************************************************************************
!  Input keywords for electrostatics.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_estat(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk

      !> @kw estat_damped
      !> Damp the Coulomb interaction at the cutoff rather than truncating it, which removes the
      !> jump in the energy an atom crossing the cutoff would otherwise cause.
      !> @needs estat_method
      !> @type logical
      if (keyword == "estat_damped") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%damped
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_damped", params%estat_options%damped)
         !> @kw estat_damped_cosine
         !> Use a cosine taper for that damping instead of the default shifted form.
         !> @needs estat_method
         !> @type logical
      else if (keyword == "estat_damped_cosine") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%damped_cosine
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_damped_cosine", params%estat_options%damped_cosine)
         !> @kw estat_dsf_alpha
         !> Damping parameter of the damped shifted force. A negative value, the default, means the
         !> run picks one from the cutoff.
         !> @units 1/A
         !> @needs estat_method
      else if (keyword == "estat_dsf_alpha") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_dsf_alpha
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_dsf_alpha", params%estat_dsf_alpha)
         !> @kw estat_gsf
         !> Use the gradient-shifted force, which subtracts the force at the cutoff so that the force
         !> goes to zero there.
         !> @needs estat_method
         !> @type logical
      else if (keyword == "estat_gsf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%gsf
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_gsf", params%estat_options%gsf)
         !> @kw estat_inner_width
         !> Width over which the inner cutoff turns the electrostatics on.
         !> @units A
         !> @needs estat_method
         !> @see estat_rcut_inner
      else if (keyword == "estat_inner_width") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_inner_width
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_inner_width", params%estat_inner_width)
         !> @kw estat_method
         !> Which electrostatic sum to use: "none", "direct" for a bare cutoff sum, "dsf" for the
         !> damped shifted force, or "gsf" for the gradient-shifted force. An unknown name warns and
         !> computes nothing. The GAP has to provide charges for any of this to run.
      else if (keyword == "estat_method") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_method
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_method", params%estat_method)
         !> @kw estat_rcut
         !> Cutoff of the electrostatic sum.
         !> @units A
         !> @needs estat_method
      else if (keyword == "estat_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_rcut", params%estat_rcut)
         !> @kw estat_rcut_inner
         !> Radius below which the electrostatics is smoothly switched off, so it does not
         !> double-count what the GAP already describes at short range.
         !> @units A
         !> @needs estat_method
         !> @see estat_inner_width
      else if (keyword == "estat_rcut_inner") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_rcut_inner
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_rcut_inner", params%estat_rcut_inner)
         !> @kw estat_self_energy_correction
         !> Subtract each charge's interaction with its own periodic images, the self term of the
         !> shifted sum.
         !> @needs estat_method
         !> @type logical
      else if (keyword == "estat_self_energy_correction") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%self_energy_correction
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_self_energy_correction", params%estat_options%self_energy_correction)
         !> @kw estat_sp
         !> Use the shifted-potential form, which shifts the energy to zero at the cutoff but leaves a
         !> discontinuous force.
         !> @needs estat_method
         !> @type logical
      else if (keyword == "estat_sp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%sp
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_sp", params%estat_options%sp)
         !> @kw estat_tsf
         !> Use the truncated shifted force. On by default.
         !> @needs estat_method
         !> @type logical
      else if (keyword == "estat_tsf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estat_options%tsf
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estat_tsf", params%estat_options%tsf)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_estat

!**************************************************************************
!  Input keywords for experimental observables: pdf, structure factor, XRD, ND, XPS.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_exp(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: i2
      integer :: nw
      logical :: matched_label
      character*64 :: keyword_notrim

!     The suffix dispatch below ("*_range", "*_file_data", "*_n_samples")
!     matches on the tail of the keyword, so it needs the name right-justified
!     in a fixed-width buffer and the index of its last character. When this
!     family was split out of read_input_file both were left behind in the
!     parent and re-declared here without ever being assigned, so the substring
!     bounds were whatever the stack held: the branch was entered at random,
!     and its backspace(unit) with no matching read re-presented the same
!     record to the caller's loop for ever. Recomputed here from the keyword,
!     exactly as read_input_file does it.
      keyword_notrim = keyword
      keyword_notrim = adjustr(keyword_notrim)
      i2 = len(keyword_notrim)

      !> @kw do_exp
      !> Compare the prediction against experimental data and add the mismatch to the energy. The
      !> individual exp_ keywords switch this on for you.
      !> @see exp_data_files, exp_energy_scales
      if (keyword == 'do_exp') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_exp
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_exp", params%do_exp)

         !> @kw do_nd
         !> Compute a neutron diffraction pattern. Needs the partial pair distributions and the
         !> structure factors it is built from, so it switches both on.
         !> @see do_pair_distribution, do_structure_factor, nd_output
      else if (keyword == 'do_nd') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_nd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_nd", params%do_nd)

         if (params%do_nd) then
            params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
         end if

         !> @kw do_pair_distribution
         !> Compute the pair distribution function g(r). Also records that the deck asked for it in as
         !> many words, which is what lets xrd_debye tell a request it must honour from a side effect
         !> it may undo.
         !> @see pair_distribution_rcut, pair_distribution_n_samples
      else if (keyword == 'do_pair_distribution') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_pair_distribution
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_pair_distribution", params%do_pair_distribution)
         params%do_pair_distribution_explicit = params%do_pair_distribution

         !> @kw do_structure_factor
         !> Compute the structure factor S(q). Needs the pair distributions when
         !> structure_factor_from_pdf is on, so it switches those on too.
         !> @see structure_factor_from_pdf, q_range_min, q_range_max
      else if (keyword == 'do_structure_factor') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_structure_factor
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_structure_factor", params%do_structure_factor)
         params%do_structure_factor_explicit = params%do_structure_factor
         if (params%do_structure_factor) then
            params%do_pair_distribution = .true.
         end if

         !> @kw do_xps
         !> Compute an X-ray photoelectron spectrum from the predicted core-electron binding energies.
         !> The GAP has to provide those.
         !> @see xps_sigma, xps_e_min, xps_e_max
      else if (keyword == 'do_xps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_xps
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_xps", params%do_xps)

         !> @kw do_xrd
         !> Compute an X-ray diffraction pattern. Needs the partial pair distributions and the
         !> structure factors it is built from, so it switches both on, unless xrd_debye takes the
         !> direct route instead.
         !> @see xrd_wavelength, xrd_debye, xrd_output
      else if (keyword == 'do_xrd') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%do_xrd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("do_xrd", params%do_xrd)

         if (params%do_xrd) then
            params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
         end if

         !> @kw exp_data_files
         !> Files holding the measured data, one per observable, in the order of exp_labels. Two
         !> columns: abscissa and intensity.
         !> @see exp_labels, n_exp
         !> @type string list
      else if (keyword == "exp_data_files") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%file_data, nw=1, params%n_exp)
         if (rank == 0) then
            do nw = 1, params%n_exp
               call print_parameter("exp_data_files", trim(params%exp_data(nw)%file_data))
            end do
         end if

         do nw = 1, params%n_exp
            if (trim(params%exp_data(nw)%file_data) == "none") then
               ! Make sure that no type of exp data is written
               params%exp_data(nw)%compute_exp = .false.
               params%exp_data(nw)%compute_similarity = .false.
               ! If the compute exp is false, then a user range must be specified
               params%exp_data(nw)%wrote_exp = .true.
            else

               call read_exp_data( &
                  params%exp_data(nw)%file_data, &
                  params%exp_data(nw)%n_data, &
                  params%exp_data(nw)%data)

               params%exp_data(nw)%compute_exp = .true.
               params%exp_data(nw)%compute_similarity = .true.
               params%exp_data(nw)%range_min = params%exp_data(nw)%data(1, 1)
               params%exp_data(nw)%range_max = params&
                    &%exp_data(nw)%data(1, params%exp_data(nw)%n_data)
            end if
         end do

         !> @kw exp_energies
         !> Add the data mismatch to the energy. Naming it enables do_exp.
         !> @see exp_energy_scales
      else if (keyword == 'exp_energies') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_energies
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("exp_energies", params%exp_energies)
! do experimental
         params%do_exp = .true.

         !> @kw exp_energy_scales
         !> Weight of each observable's mismatch in the total energy, one value per observable. This
         !> is what sets how hard the data pulls against the potential. Also seeds the initial and
         !> final values of a ramp.
         !> @units eV
         !> @see exp_energy_scales_final, n_exp
      else if (keyword == 'exp_energy_scales' .or. keyword&
           &== 'exp_energy_scales_initial' .or. keyword&
           &== 'exp_energy_scales_beg') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params&
              &%exp_energy_scales(nw), nw=1, params&
              &%n_exp)
         if (rank == 0) call print_parameters("exp_energy_scales", params%exp_energy_scales, "eV")

         ! Set the final gamma to the initial in case
         do nw = 1, params%n_exp
            params%exp_energy_scales_initial(nw) = params%exp_energy_scales(nw)
            params%exp_energy_scales_final(nw) = params%exp_energy_scales(nw)
         end do

         !> @kw exp_energy_scales_final, exp_energy_scales_end
         !> End point of a linear ramp in the weights over the run, so the data can be brought in
         !> gradually. Without it the weights stay at exp_energy_scales.
         !> @units eV
         !> @see exp_energy_scales
         !> @type real list
         !> @default same as exp_energy_scales
      else if (keyword == 'exp_energy_scales_final' .or. keyword == 'exp_energy_scales_end') then
         backspace (unit)
         if (params%n_moments > 0) then
            read (unit, *, iostat=iostatus) cjunk, cjunk, (params&
             &%exp_energy_scales_final(nw), nw=1, params&
             &%n_moments)
         else
            read (unit, *, iostat=iostatus) cjunk, cjunk, (params&
                 &%exp_energy_scales_final(nw), nw=1, params&
                 &%n_exp)
         end if

         if (rank == 0) call print_parameters("exp_energy_scales_final", params%exp_energy_scales_final, "eV")
         !> @kw exp_forces
         !> Compute the derivative of the data mismatch with respect to positions and add it to the
         !> forces. Naming it enables do_exp. Without it the mismatch contributes to the energy alone.
         !> @see exp_energies
      else if (keyword == 'exp_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("exp_forces", params%exp_forces)
! do experimental
         params%do_exp = .true.

         !> @kw exp_input_type
         !> What the corresponding data file contains, one per observable, when it is not the
         !> observable's default representation: for instance "D(r)" for a pair distribution supplied
         !> as the reduced form.
         !> @see exp_labels, pair_distribution_output
         !> @type string list
      else if (keyword == 'exp_input_type') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%input, nw=1, params%n_exp)
         if (rank == 0) then
            do nw = 1, params%n_exp
               call print_parameter("exp_input_type", trim(params%exp_data(nw)%input))
            end do
         end if

         !> @kw exp_labels
         !> Which observable each data file is, one per file, from "xps", "xrd", "saxs",
         !> "pair_distribution" and "structure_factor". Matching a label to a file is what marks that
         !> observable valid, and an observable that was asked for but has no data stays switched off.
         !> @see exp_data_files, n_exp
         !> @type string list
      else if (keyword == "exp_labels") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%label, nw=1, params%n_exp)
         if (rank == 0) then
            do nw = 1, params%n_exp
               call print_parameter("exp_labels", trim(params%exp_data(nw)%label))
            end do
         end if
         do nw = 1, params%n_exp
            call upper_to_lower_case(params%exp_data(nw)%label)
            if (trim(params%exp_data(nw)%label) == "xps") then
               params%xps_idx = nw
               if (rank == 0) write (*, *) ' - Valid exp. XPS found                |'

            else if (trim(params%exp_data(nw)%label) == "xrd") then
               params%xrd_idx = nw
               params%valid_xrd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if (trim(params%exp_data(nw)%label) == "nd") then
               params%nd_idx = nw
               params%valid_nd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. ND found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.

            else if (trim(params%exp_data(nw)%label) == "saxs") then
               params%valid_xrd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if (trim(params%exp_data(nw)%label) == "pair_distribution") then
               params%pdf_idx = nw
               params%valid_pdf = .true.
               if (rank == 0) write (*, *) ' - Valid exp. pair distribution found  |'
            else if (trim(params%exp_data(nw)%label) == "structure_factor") then
               params%sf_idx = nw
               params%valid_sf = .true.
               if (rank == 0) write (*, *) ' - Valid exp. structure factor found   |'
            else if (trim(params%exp_data(nw)%label) == "ir") then
               params%ir_idx = nw
               params%valid_ir = .true.
               if (rank == 0) write (*, *) ' - Valid exp. IR spectrum found        |'
            end if
         end do
         !> @kw exp_n_samples
         !> Number of points to predict for each observable, one per observable, overriding that
         !> observable's own n_samples keyword.
         !> @see exp_labels
         !> @type integer list
      else if (keyword == "exp_n_samples") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%n_samples, nw=1, params%n_exp)
         if (rank == 0) then
            do nw = 1, params%n_exp
               call print_parameter("exp_n_samples", params%exp_data(nw)%n_samples)
            end do
         end if

         !> @kw exp_similarity_type
         !> How the mismatch between prediction and data is measured: "squared_diff" for a sum of
         !> squared residuals, or "similarity"/"overlap" for a normalised overlap.
         !> @see exp_energy_scales
      else if (keyword == 'exp_similarity_type') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_similarity_type
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("exp_similarity_type", params%exp_similarity_type)
         !> @kw n_exp
         !> Number of experimental observables. Give it before exp_data_files, exp_labels and
         !> exp_energy_scales, which it allocates.
         !> @see exp_labels
      else if (keyword == 'n_exp') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_exp
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_exp", params%n_exp)
         allocate (params%exp_data(1:params%n_exp))
         allocate (params%exp_energy_scales(1:params%n_exp))
         allocate (params%exp_energy_scales_initial(1:params%n_exp))
         allocate (params%exp_energy_scales_final(1:params%n_exp))

         ! Turning on exp prediction
         params%do_exp = .true.

         !> @kw nd_n_samples
         !> Number of points in the predicted neutron pattern. Writes the same sample count the XRD
         !> and structure-factor routines use, since the three share one q grid.
         !> @see xrd_n_samples
      else if (keyword == "nd_n_samples") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nd_n_samples", params%xrd_n_samples)
         params%structure_factor_n_samples = params%nd_n_samples

         !> @kw nd_output
         !> Abscissa of the neutron pattern, in the same choices as xrd_output. Writes xrd_output as
         !> well, the two patterns sharing a grid.
         !> @see xrd_output
      else if (keyword == 'nd_output') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%nd_output
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nd_output", params%nd_output)

         ! else if(keyword=='xrd_input')then
         !    backspace(10)
         !    read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_output

         !> @kw nd_rcut
         !> Real-space cutoff of the pair sum entering the neutron pattern.
         !> @units A
         !> @see pair_distribution_rcut
      else if (keyword == "nd_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%nd_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nd_rcut", params%nd_rcut)

         !> @kw nd_wavelength
         !> Neutron wavelength. Nothing consults it: the pattern is built on the grid xrd_wavelength
         !> defines.
         !> @see xrd_wavelength
         !> @type real
         !> @status ignored
      else if (keyword == 'nd_wavelength') then
         call ignored_keyword(unit, iostatus, rank, keyword)

         !> @kw pair_distribution_kde_sigma
         !> Width of the Gaussian each pair distance is smeared with when binning g(r). Zero, the
         !> default, means plain histogram binning. A non-zero width makes g(r) differentiable, which
         !> is what the forces need.
         !> @units A
         !> @see exp_forces
      else if (keyword == "pair_distribution_kde_sigma") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_kde_sigma
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("pair_distribution_kde_sigma", params%pair_distribution_kde_sigma)

         !> @kw pair_distribution_n_samples
         !> Number of r points in the predicted g(r).
         !> @see r_range_min, r_range_max
      else if (keyword == "pair_distribution_n_samples") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("pair_distribution_n_samples", params%pair_distribution_n_samples)

         !> @kw pair_distribution_output
         !> What g(r) is reported as: "pdf" for g(r) itself, or "D(r)" for the reduced pair
         !> distribution 4 pi r rho (g(r) - 1). The energy and its gradient follow the choice, so it
         !> changes the fit and not only the plot.
         !> @see exp_input_type
      else if (keyword == 'pair_distribution_output') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_output
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("pair_distribution_output", params%pair_distribution_output)

         !> @kw pair_distribution_partial
         !> Resolve g(r) into per-species-pair partials. The structure factor is built from these, so
         !> it is on by default.
         !> @see do_structure_factor
      else if (keyword == "pair_distribution_partial") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_partial
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("pair_distribution_partial", params%pair_distribution_partial)

         !> @kw pair_distribution_rcut
         !> Largest pair distance binned into g(r). A hard cut, not a smooth one, which is why a
         !> finite-difference check of the virial needs a small strain.
         !> @units A
         !> @see r_range_max
      else if (keyword == "pair_distribution_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("pair_distribution_rcut", params%pair_distribution_rcut)

         !> @kw q_range_max
         !> Upper end of the q range the structure factor and the diffraction patterns are predicted
         !> on.
         !> @units 1/A
         !> @see q_range_min, q_units
      else if (keyword == "q_range_max") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%q_range_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("q_range_max", params%q_range_max)

         !> @kw q_range_min
         !> Lower end of that range.
         !> @units 1/A
         !> @see q_range_max
      else if (keyword == "q_range_min") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%q_range_min
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("q_range_min", params%q_range_min)

         !> @kw q_units
         !> What the abscissa of the supplied and predicted patterns means: "q" or "saxs" for a
         !> scattering vector, "xrd" or "twotheta" for a diffraction angle. The conversion between
         !> them uses xrd_wavelength.
         !> @see xrd_wavelength
      else if (keyword == "q_units") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%q_units
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("q_units", params%q_units)

         !> @kw r_range_max
         !> Upper end of the r range g(r) is reported on. Points beyond pair_distribution_rcut are
         !> empty.
         !> @units A
         !> @see pair_distribution_rcut
      else if (keyword == "r_range_max") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%r_range_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("r_range_max", params%r_range_max)

         !> @kw r_range_min
         !> Lower end of that range.
         !> @units A
         !> @see r_range_max
      else if (keyword == "r_range_min") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%r_range_min
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("r_range_min", params%r_range_min)

         !> @kw sf_output
         !> Abscissa of the structure factor, in the same choices as q_units.
         !> @see q_units
      else if (keyword == 'sf_output') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%sf_output
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("sf_output", params%sf_output)

         !> @kw structure_factor_from_pdf
         !> Build S(q) by Fourier-transforming the binned partial pair distributions, rather than from
         !> the positions directly. The transform is far cheaper and is what the analytic forces are
         !> written for.
         !> @see pair_distribution_partial
      else if (keyword == "structure_factor_from_pdf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_from_pdf
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("structure_factor_from_pdf", params%structure_factor_from_pdf)
         !> @kw structure_factor_matrix
         !> Assemble the transform as a matrix product over the r grid instead of looping the pairs.
         !> Much faster, at the cost of holding the kernel matrix.
         !> @see structure_factor_matrix_forces
      else if (keyword == "structure_factor_matrix") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_matrix
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("structure_factor_matrix", params%structure_factor_matrix)
         !> @kw structure_factor_matrix_forces
         !> Use the same matrix assembly for the derivatives. Falls back to the loop when the run
         !> cannot batch.
         !> @see structure_factor_matrix
      else if (keyword == "structure_factor_matrix_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_matrix_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("structure_factor_matrix_forces", params%structure_factor_matrix_forces)

         !> @kw structure_factor_n_samples
         !> Number of q points in the predicted S(q).
         !> @see q_range_min, q_range_max
      else if (keyword == "structure_factor_n_samples") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("structure_factor_n_samples", params%structure_factor_n_samples)

         !> @kw structure_factor_window
         !> Apply a Lorch window to the truncated pair distribution before transforming, which
         !> suppresses the ringing the hard cut at pair_distribution_rcut would otherwise put into
         !> S(q).
         !> @see pair_distribution_rcut
      else if (keyword == 'structure_factor_window') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_window
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("structure_factor_window", params%structure_factor_window)

         !> @kw xps_e_max
         !> Upper end of the binding-energy range the XPS spectrum is predicted on.
         !> @units eV
         !> @needs do_xps
         !> @see xps_e_min
      else if (keyword == 'xps_e_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xps_e_max
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xps_e_max", params%xps_e_max)

         !> @kw xps_e_min
         !> Lower end of that range.
         !> @units eV
         !> @needs do_xps
         !> @see xps_e_max
      else if (keyword == 'xps_e_min') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xps_e_min
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xps_e_min", params%xps_e_min)

         !> @kw xps_force_type
         !> Which XPS force expression to use. Nothing consults it; the spectrum is differentiated one
         !> way only.
         !> @type string
         !> @status ignored
      else if (keyword == 'xps_force_type') then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw xps_n_samples
         !> Number of points in the predicted XPS spectrum.
         !> @needs do_xps
      else if (keyword == 'xps_n_samples') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xps_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xps_n_samples", params%xps_n_samples)

         !> @kw xps_sigma
         !> Width of the Gaussian each core-electron binding energy is broadened by to make a
         !> spectrum.
         !> @units eV
         !> @needs do_xps
      else if (keyword == 'xps_sigma') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xps_sigma
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xps_sigma", params%xps_sigma)
         !> @kw xrd_alpha
         !> Exponent of the sharpening applied to the diffraction pattern.
         !> @needs do_xrd
         !> @see xrd_damping
      else if (keyword == 'xrd_alpha') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_alpha
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_alpha", params%xrd_alpha)
         !> @kw xrd_damping
         !> Damping applied to the pair sum with distance, which suppresses the truncation ripple from
         !> the finite cutoff. Zero, the default, is no damping.
         !> @needs do_xrd
         !> @see xrd_alpha
      else if (keyword == 'xrd_damping') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_damping
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_damping", params%xrd_damping)
         !> @kw xrd_debye
         !> Compute the pattern from the N^2 Debye sum over positions instead of transforming the
         !> partial pair distributions. The two routes answer the same question by different
         !> approximations: the pdf route bins distances and truncates at pair_distribution_rcut, the
         !> Debye route sums every pair exactly. Turning it on makes the pdf and sf calculations
         !> unnecessary, and whichever of them was switched on only to feed the XRD is switched off
         !> again.
         !> @needs do_xrd
         !> @see pair_distribution_rcut
      else if (keyword == 'xrd_debye') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_debye
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_debye", params%xrd_debye)
         !> @kw xrd_lorentz_polarization
         !> Multiply the predicted pattern by the powder Lorentz-polarization factor. Only the Debye
         !> route applies it, where the energy and the gradient stay consistent because it is one
         !> multiplicative weight per q.
         !> @needs xrd_debye
         !> @see xrd_lp_polarization
      else if (keyword == 'xrd_lorentz_polarization') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_lorentz_polarization
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_lorentz_polarization", params%xrd_lorentz_polarization)
         !> @kw xrd_lp_polarization
         !> Degree of polarization P of the incident beam in the factor (1 + P cos^2 2theta) / (sin^2
         !> theta cos theta). P = 1 is an unpolarized source; a graphite monochromator at 2theta_M
         !> gives P = cos^2(2theta_M).
         !> @needs xrd_lorentz_polarization
      else if (keyword == 'xrd_lp_polarization') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_lp_polarization
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_lp_polarization", params%xrd_lp_polarization)
         !> @kw xrd_lp_sin_theta_min
         !> Below this |sin theta| the Lorentz factor 1/(sin^2 theta cos theta) is unusable, so the
         !> whole factor is set to zero rather than blown up.
         !> @needs xrd_lorentz_polarization
      else if (keyword == 'xrd_lp_sin_theta_min') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_lp_sin_theta_min
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_lp_sin_theta_min", params%xrd_lp_sin_theta_min)
         !> @kw xrd_iwasa
         !> Iwasa correction switch from an earlier design; the correction is applied unconditionally
         !> now.
         !> @status deprecated
      else if (keyword == 'xrd_iwasa') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_iwasa
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_iwasa", params%xrd_iwasa)
         !> @kw xrd_method
         !> Which scattering factors the pattern is built from. "xrd" uses the X-ray form factors.
         !> @needs do_xrd
      else if (keyword == 'xrd_method') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_method
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_method", params%xrd_method)

         !> @kw xrd_n_samples
         !> Number of points in the predicted X-ray pattern. Sets the structure-factor sample count
         !> too, the two sharing a q grid.
         !> @needs do_xrd
      else if (keyword == "xrd_n_samples") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_n_samples", params%xrd_n_samples)
         params%structure_factor_n_samples = params%xrd_n_samples

         !> @kw xrd_output
         !> Abscissa the pattern is reported on, in the same choices as q_units.
         !> @needs do_xrd
         !> @see q_units
      else if (keyword == 'xrd_output') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_output
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_output", params%xrd_output)
         !> @kw xrd_rcut
         !> Real-space cutoff of the pair sum entering the X-ray pattern.
         !> @units A
         !> @needs do_xrd
         !> @see pair_distribution_rcut
      else if (keyword == "xrd_rcut") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_rcut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_rcut", params%xrd_rcut)

         !> @kw xrd_wavelength
         !> Wavelength of the incident X-rays. Fixes the map between q and 2theta, so it matters
         !> whenever q_units or xrd_output name an angle.
         !> @units A
         !> @see q_units
      else if (keyword == 'xrd_wavelength') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%xrd_wavelength
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("xrd_wavelength", params%xrd_wavelength)
!     Per-observable keywords, dispatched on the tail of the name:
!     <label>_range, <label>_file_data and <label>_n_samples, where <label> is
!     an entry in exp_labels. Last in the chain on purpose -- xrd_n_samples and
!     its four siblings have branches of their own above, and those set the
!     global sample counts rather than one observable's. This catches the names
!     no dedicated branch claims.
      else if (keyword_notrim(i2 - 4:i2) == "range" .or. keyword_notrim(i2 - 8:i2) ==&
            & "file_data" .or. keyword_notrim(i2 - 8:i2) == "n_samples") then
!        Ending in one of those three suffixes is necessary but not sufficient:
!        the stem has to name an observable this deck actually declared. Give
!        the line back untouched when it does not, so the next family gets its
!        turn. Claiming it and finding nothing to do would leave the
!        backspace(unit) below unmatched by a read, and the caller's loop would
!        be handed the same record for ever.
         matched_label = .false.
         do nw = 1, params%n_exp
            if (keyword == trim(params%exp_data(nw)%label)//"_range" .or. &
                keyword == trim(params%exp_data(nw)%label)//"_file_data" .or. &
                keyword == trim(params%exp_data(nw)%label)//"_n_samples") matched_label = .true.
         end do
         if (.not. matched_label) return

         backspace (unit)
         ! Check if experimental range or data files are specified
         do nw = 1, params%n_exp
            ! See if the keyword matches any exp observables
            if (keyword == trim(params%exp_data(nw)%label)//"_range") then
               ! Expect two values which are in order of lower higher for the range to do the prediction
               params%exp_data(nw)%user_range = .true.
               read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%range_min, params%exp_data(nw)%range_max
            elseif (keyword == trim(params%exp_data(nw)%label)//"_file_data") then

               read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%file_data
               if (trim(params%exp_data(nw)%file_data) /= "none") then

                  call read_exp_data( &
                     params%exp_data(nw)%file_data, &
                     params%exp_data(nw)%n_data, &
                     params%exp_data(nw)%data)

                  params%exp_data(nw)%wrote_exp = .false.
                  params%exp_data(nw)%compute_exp = .true.
                  params%exp_data(nw)%compute_similarity = .true.
                  params%exp_data(nw)%range_min = params%exp_data(nw)%data(1, 1)
                  params%exp_data(nw)%range_max = params&
                       &%exp_data(nw)%data(1, params%exp_data(nw)%n_data)
               elseif (trim(params%exp_data(nw)%file_data) == "none") then
                  ! Make sure that no type of exp data is written
                  params%exp_data(nw)%compute_exp = .false.
                  params%exp_data(nw)%compute_similarity = .false.
                  ! If the compute exp is false, then a user range must be specified
                  params%exp_data(nw)%wrote_exp = .true.

               end if
            elseif (keyword == trim(params%exp_data(nw)%label)//"_n_samples") then
               read (unit, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%n_samples
            end if
         end do

      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_exp

!**************************************************************************
!  Input keywords for electronic stopping and the electron-phonon model.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_stopping(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: i

      !> @kw adapt_emax
      !> Largest energy change tolerated over one step; the adaptive time step shrinks until the
      !> change is under this.
      !> @units eV
      !> @modes md
      !> @needs adaptive_time
      if (keyword == 'adapt_emax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adapt_emax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adapt_emax", params%adapt_emax)

        !! --------------------------                        ******** until here for adaptive time

        !! ------- option for radiation cascade simulation with electronic stopping

         !> @kw adapt_tmax
         !> Upper bound on the adaptive time step.
         !> @units fs
         !> @modes md
         !> @needs adaptive_time
      else if (keyword == 'adapt_tmax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adapt_tmax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adapt_tmax", params%adapt_tmax)
         !> @kw adapt_tmin
         !> Lower bound on the adaptive time step. Reaching it means the dynamics cannot be resolved
         !> and the run should be reconsidered.
         !> @units fs
         !> @modes md
         !> @needs adaptive_time
      else if (keyword == 'adapt_tmin') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adapt_tmin
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adapt_tmin", params%adapt_tmin)
         !> @kw adapt_tstep_interval
         !> How often, in steps, the adaptive time step is re-examined.
         !> @modes md
         !> @needs adaptive_time
      else if (keyword == 'adapt_tstep_interval') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adapt_tstep_interval
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adapt_tstep_interval", params%adapt_tstep_interval)
         if (params%adapt_tstep_interval <= 0) then
            write (*, *) "ERROR: Interval of timesteps in adaptive time-step must be positive."
            stop
         end if
         !> @kw adapt_xmax
         !> Largest displacement of any atom tolerated over one step.
         !> @units A
         !> @modes md
         !> @needs adaptive_time
      else if (keyword == 'adapt_xmax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adapt_xmax
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adapt_xmax", params%adapt_xmax)
         !> @kw adaptive_time
         !> Choose the time step from how fast the fastest atom is moving, rather than holding md_step
         !> fixed. Written for radiation cascades, where the first femtoseconds need a far smaller
         !> step than the rest.
         !> @modes md
         !> @see adapt_tmin, adapt_tmax, adapt_xmax, adapt_emax
      else if (keyword == 'adaptive_time') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%adaptive_time
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("adaptive_time", params%adaptive_time)
         !> @kw eel_cut
         !> Kinetic energy below which an atom is too slow to lose energy to the electrons, so
         !> electronic stopping is not applied to it.
         !> @units eV
         !> @modes md
         !> @needs electronic_stopping
      else if (keyword == 'eel_cut') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eel_cut
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eel_cut", params%eel_cut)
         if (params%eel_cut <= 0) then
            write (*, *) "ERROR: Cut off energy for electronic stopping should be positive, few tens of eV!"
            stop
         end if
         !> @kw eel_freq_out
         !> How often, in steps, to report the energy lost to electronic stopping.
         !> @modes md
         !> @needs electronic_stopping
      else if (keyword == 'eel_freq_out') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eel_freq_out
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eel_freq_out", params%eel_freq_out)
         !> @kw electronic_stopping
         !> Drain energy from fast atoms into the electronic system, the dominant energy loss at the
         !> start of a radiation cascade. Reads its stopping powers from estop_filename.
         !> @modes md
         !> @see estop_filename, eel_cut
      else if (keyword == 'electronic_stopping') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%electronic_stopping
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("electronic_stopping", params%electronic_stopping)
         !> @kw eph_betafile
         !> File holding the beta(rho) coupling data the EPH model needs, one entry per species.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_betafile') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_betafile
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_betafile", params%eph_betafile)
         !> @kw eph_box_limits
         !> Extent of the electronic grid, as xmin xmax ymin ymax zmin zmax. Also written into the six
         !> individual bounds the solver reads.
         !> @units A
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_box_limits') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params%eph_box_limits(i), i=1, 6)
         params%in_x0 = params%eph_box_limits(1); params%in_x1 = params%eph_box_limits(2)
         params%in_y0 = params%eph_box_limits(3); params%in_y1 = params%eph_box_limits(4)
         params%in_z0 = params%eph_box_limits(5); params%in_z1 = params%eph_box_limits(6)
         if (rank == 0) call print_parameters("eph_box_limits", params%eph_box_limits)
         !> @kw eph_c_e
         !> Electronic heat capacity per unit volume, used when the electronic temperature is evolved
         !> rather than held fixed.
         !> @units eV/A^3/K
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_c_e') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_c_e
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_c_e", params%eph_c_e)
         !> @kw eph_e_prev_time
         !> Electronic energy carried over from a previous run, for restarting a cascade.
         !> @units eV
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_e_prev_time') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_E_prev_time
         if (rank == 0) call print_parameter("eph_e_prev_time", params%eph_e_prev_time)
         !> @kw eph_fdm_option
         !> How the electronic temperature grid is initialised: from eph_Tinfile, or uniformly at
         !> eph_Ti_e.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_fdm_option') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_fdm_option
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_fdm_option", params%eph_fdm_option)
         !> @kw eph_fdm_steps
         !> Number of electronic diffusion sub-steps taken per MD step, the electronic system being
         !> much faster than the ionic one.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_fdm_steps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_fdm_steps
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_fdm_steps", params%eph_fdm_steps)
         !> @kw eph_freq_mesh_tout
         !> How often, in steps, to write the whole electronic temperature grid.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_freq_mesh_tout') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_freq_mesh_Tout
         if (rank == 0) call print_parameter("eph_freq_mesh_tout", params%eph_freq_mesh_tout)
         !> @kw eph_freq_tout
         !> How often, in steps, to write the electronic temperature summary.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_freq_tout') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_freq_Tout
         if (rank == 0) call print_parameter("eph_freq_tout", params%eph_freq_Tout)
         !> @kw eph_friction_option
         !> Which friction term of the electron-phonon coupling to apply, or none.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_friction_option') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_friction_option
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_friction_option", params%eph_friction_option)
         !> @kw eph_gsx
         !> Number of electronic grid cells along x.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_gsx') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_gsx
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_gsx", params%eph_gsx)
         !> @kw eph_gsy
         !> Number of electronic grid cells along y.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_gsy') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_gsy
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_gsy", params%eph_gsy)
         !> @kw eph_gsz
         !> Number of electronic grid cells along z.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_gsz') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_gsz
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_gsz", params%eph_gsz)
         !> @kw eph_kappa_e
         !> Electronic thermal conductivity entering the diffusion of the electronic temperature.
         !> @units eV/A/fs/K
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_kappa_e') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_kappa_e
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_kappa_e", params%eph_kappa_e)
         !> @kw eph_md_last_step
         !> Step number the previous run finished at, for restarting a cascade.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_md_last_step') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_md_last_step
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_md_last_step", params%eph_md_last_step)
         !> @kw eph_md_prev_time
         !> Simulated time the previous run finished at, for restarting a cascade.
         !> @units fs
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_md_prev_time') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_md_prev_time
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_md_prev_time", params%eph_md_prev_time)
         !> @kw eph_random_option
         !> Which random term of the electron-phonon coupling to apply, or none. Paired with the
         !> friction term it is what makes the coupling a thermostat rather than a drain.
         !> @modes md
         !> @needs nonadiabatic_processes
         !> @see eph_friction_option
      else if (keyword == 'eph_random_option') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_random_option
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_random_option", params%eph_random_option)
         !> @kw eph_rho_e
         !> Electronic density entering the coupling.
         !> @units 1/A^3
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_rho_e') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_rho_e
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("eph_rho_e", params%eph_rho_e)
         !> @kw eph_ti_e
         !> Initial electronic temperature, when the grid is started uniform rather than read from a
         !> file.
         !> @units K
         !> @modes md
         !> @needs nonadiabatic_processes
         !> @see eph_fdm_option
      else if (keyword == 'eph_ti_e') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_Ti_e

        !! --------------------                        ******** until here for electronic stopping based on EPH model

         if (rank == 0) call print_parameter("eph_ti_e", params%eph_Ti_e)
         !> @kw eph_tinfile
         !> File holding an initial electronic temperature grid.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_tinfile') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_Tinfile
         if (rank == 0) call print_parameter("eph_tinfile", trim(params%eph_Tinfile))
         !> @kw eph_toutfile
         !> File the electronic temperature grid is written to.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'eph_toutfile') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%eph_Toutfile
         if (rank == 0) call print_parameter("eph_toutfile", trim(params%eph_Toutfile))
         !> @kw estop_filename
         !> File holding the electronic stopping power per species as a function of velocity.
         !> @modes md
         !> @needs electronic_stopping
      else if (keyword == 'estop_filename') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%estop_filename
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("estop_filename", params%estop_filename)

        !! -------------------------------                ******** until here for electronic stopping

        !! ------- option for radiation cascade simulation with EPH model

         !> @kw model_eph
         !> Which electron-phonon model to evaluate.
         !> @modes md
         !> @needs nonadiabatic_processes
      else if (keyword == 'model_eph') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%model_eph
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("model_eph", params%model_eph)
         !> @kw nonadiabatic_processes
         !> Exchange energy between the ions and an explicit electronic subsystem through the EPH
         !> model, rather than the one-way drain of electronic stopping.
         !> @modes md
         !> @see model_eph, eph_friction_option, eph_random_option
      else if (keyword == 'nonadiabatic_processes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%nonadiabatic_processes
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("nonadiabatic_processes", params%nonadiabatic_processes)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_stopping

!**************************************************************************
!  Input keywords for local properties and the core potential.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_local_properties(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: nw

      !> @kw compute_local_properties
      !> Which local properties to predict, by name, one per property. The GAP has to carry a model
      !> for each. n_local_properties must appear first.
      !> @needs n_local_properties
      !> @type string list
      if (keyword == 'compute_local_properties') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params&
              &%compute_local_properties(nw), nw=1, params&
             &%n_local_properties)
         if (rank == 0) call print_parameters("compute_local_properties", params%compute_local_properties)

         !> @kw core_pot_buffer
         !> Width over which the core potential is smoothly turned off approaching core_pot_cutoff, so
         !> the tabulated potential and its gradient reach zero together.
         !> @units A
         !> @see core_pot_cutoff
      else if (keyword == "core_pot_buffer") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%core_pot_buffer
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("core_pot_buffer", params%core_pot_buffer)
         !> @kw core_pot_cutoff
         !> Distance beyond which the tabulated core potential is dropped. The default is effectively
         !> infinite, meaning the whole table is used.
         !> @units A
         !> @see core_pot_buffer
      else if (keyword == "core_pot_cutoff") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%core_pot_cutoff
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("core_pot_cutoff", params%core_pot_cutoff)
         !> @kw n_local_properties
         !> Number of local properties to predict. Give it before compute_local_properties and
         !> write_local_properties, which it allocates.
         !> @see compute_local_properties
      else if (keyword == 'n_local_properties') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_local_properties
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_local_properties", params%n_local_properties)
         allocate (params%write_local_properties(1:params%n_local_properties))
         allocate (params%compute_local_properties(1:params%n_local_properties))
         params%write_local_properties = .true.
         !> @kw print_lp_forces
         !> Local-property force reporting switch; nothing consults it.
         !> @type logical
         !> @status ignored
      else if (keyword == 'print_lp_forces') then
         call ignored_keyword(unit, iostatus, rank, keyword)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_local_properties

!**************************************************************************
!  Input keywords for what gets written, and how often.
!
!  Sets keyword_found when it recognises the keyword and leaves it alone
!  otherwise, so read_input_file can offer the line to the next family.
   subroutine read_options_output(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk
      integer :: nw

      !> @kw print_progress
      !> Print a progress line as the run advances.
      if (keyword == 'print_progress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%print_progress
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("print_progress", params%print_progress)
         !> @kw write_exp
         !> Write the predicted experimental observables to file.
         !> @needs do_exp
      else if (keyword == "write_exp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_exp
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_exp", params%write_exp)

         !> @kw write_fixes
         !> Include the per-atom fix flags as a column in the XYZ output.
      else if (keyword == 'write_fixes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_fixes
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_fixes", params%write_fixes)
         !> @kw write_forces
         !> Include forces as columns in the XYZ output.
      else if (keyword == 'write_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_forces
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_forces", params%write_forces)
         !> @kw write_ir
         !> Write the predicted IR spectrum. ir_spectrum.dat holds the current one, with the
         !> sampling interval, Nyquist limit, resolution and trustworthy band in its header, and
         !> with the experiment beside it when there is one. On every write_xyz step a two-column
         !> block is also appended to ir_prediction.dat, giving the prediction over the whole
         !> trajectory; the blocks are separated by one blank line, the same format
         !> xrd_prediction.dat uses, so gnuplot reaches block i with "every :::i::i" rather than
         !> with "index". A fitting run additionally writes ir_exp.dat, the experiment restricted
         !> to the fitted range. Needs either do_ir, which implies this switch, or "ir" in
         !> exp_labels. ir_write_spectrum is the old name for it.
         !> @modes md
         !> @see do_ir exp_labels write_xyz ir_stride
      else if (keyword == 'write_ir') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_ir
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_ir", params%write_ir)
         !> @kw write_hirshfeld_v
         !> Hirshfeld-volume output switch; the volumes follow the local-property output instead, so
         !> nothing consults it.
         !> @see write_local_properties
         !> @type logical
         !> @status ignored
      else if (keyword == 'write_hirshfeld_v') then
         call ignored_keyword(unit, iostatus, rank, keyword)
         !> @kw write_local_energies
         !> Include the per-atom energy decomposition as a column in the XYZ output.
      else if (keyword == 'write_local_energies') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_local_energies
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_local_energies", params%write_local_energies)
         !> @kw write_local_properties
         !> Which predicted local properties to write, one flag per property in the order of
         !> compute_local_properties.
         !> @needs n_local_properties
         !> @see compute_local_properties
         !> @type logical list
      else if (keyword == 'write_local_properties') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, (params&
              &%write_local_properties(nw), nw=1, params&
              &%n_local_properties)
         if (rank == 0) call print_parameters("write_local_properties", params%write_local_properties)
         !> @kw write_lv
         !> Write the lattice vectors on every frame.
      else if (keyword == 'write_lv') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_lv
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_lv", params%write_lv)

        !! ------- option for doing simulation with adaptive time step

         !> @kw write_masses
         !> Include atomic masses as a column in the XYZ output.
      else if (keyword == 'write_masses') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_masses
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_masses", params%write_masses)
         !> @kw write_nd
         !> Write the predicted neutron diffraction pattern to file.
         !> @needs do_nd
      else if (keyword == "write_nd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_nd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_nd", params%write_nd)

         !> @kw write_pair_distribution
         !> Write the predicted pair distribution to file.
         !> @needs do_pair_distribution
      else if (keyword == "write_pair_distribution") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_pair_distribution
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_pair_distribution", params%write_pair_distribution)
         !> @kw write_pressure
         !> Include the pressure in the frame comment line.
      else if (keyword == 'write_pressure') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_pressure
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_pressure", params%write_pressure)
         !> @kw write_stress
         !> Include the stress tensor in the frame comment line.
         !> @see write_virial
      else if (keyword == 'write_stress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_stress
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_stress", params%write_stress)
         !> @kw write_structure_factor
         !> Write the predicted structure factor to file.
         !> @needs do_structure_factor
      else if (keyword == "write_structure_factor") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_structure_factor
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_structure_factor", params%write_structure_factor)

         !> @kw write_thermo
         !> Write a thermodynamic summary line every this many steps. Zero switches it off.
         !> @modes md mc
      else if (keyword == 'write_thermo') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_thermo
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_thermo", params%write_thermo)
         !> @kw write_velocities
         !> Include velocities as columns in the XYZ output.
      else if (keyword == 'write_velocities') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_velocities
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_velocities", params%write_velocities)
         !> @kw write_virial
         !> Include the virial tensor in the frame comment line.
         !> @see write_stress
      else if (keyword == 'write_virial') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_virial
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_virial", params%write_virial)
         !> @kw write_xrd
         !> Write the predicted X-ray diffraction pattern to file.
         !> @needs do_xrd
      else if (keyword == "write_xrd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_xrd
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_xrd", params%write_xrd)

         !> @kw write_xyz
         !> Write a trajectory frame every this many steps. Zero, the default, writes only the last.
         !> @modes md mc
      else if (keyword == 'write_xyz') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%write_xyz
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("write_xyz", params%write_xyz)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_output

   subroutine read_options_gpu(unit, iostatus, rank, keyword, mode, params, n_species, keyword_found)

      implicit none

      integer, intent(in) :: unit, rank, n_species
      integer, intent(inout) :: iostatus
      character*64, intent(in) :: keyword
      character(len=*), intent(in) :: mode
      type(input_parameters), intent(inout) :: params
      logical, intent(inout) :: keyword_found

      character*64 :: cjunk

      !> @kw gpu_batched
      !> Evaluate the experimental observables in batches on the device, staging one batch of pairs
      !> at a time. Switching it off takes the unbatched route, which holds everything at once and
      !> is the reference the batched route is checked against.
      if (keyword == 'gpu_batched') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gpu_batched
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gpu_batched", params%gpu_batched)
         !> @kw gpu_low_memory
         !> Free each batch's device arrays as soon as its contribution has been accumulated, rather
         !> than holding them for the whole calculation. Slower and much smaller.
         !> @see gpu_batched
      else if (keyword == 'gpu_low_memory') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gpu_low_memory
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gpu_low_memory", params%gpu_low_memory)
         !> @kw gpu_max_batch_size
         !> Largest device allocation a single batch may make. The batch count follows from it.
         !> @units GB
         !> @see gpu_n_batches
      else if (keyword == 'gpu_max_batch_size') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gpu_max_batch_size
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gpu_max_batch_size", params%gpu_max_batch_size)
         !> @kw gpu_n_batches
         !> Split the work into exactly this many batches, instead of choosing the count from the
         !> memory budget.
         !> @see gpu_max_batch_size
      else if (keyword == 'gpu_n_batches') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%gpu_n_batches
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("gpu_n_batches", params%gpu_n_batches)
         !> @kw n_batches
         !> Split the SOAP descriptor calculation into exactly this many batches, instead of choosing
         !> the count from max_Gbytes_per_process.
         !> @see max_gbytes_per_process
      else if (keyword == 'n_batches') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, params%n_batches
         call check_iostatus(iostatus, keyword)
         if (rank == 0) call print_parameter("n_batches", params%n_batches)
      else
         return
      end if

      keyword_found = .true.

   end subroutine read_options_gpu

!**************************************************************************

   ! Deprecated Variable Check
   subroutine check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
      ! Input variables
      integer, intent(in) :: n_deprecated
      character*64, intent(in) :: deprecated_keywords(:)
      character*64, intent(in) :: updated_keywords(:)
      character*64, intent(in) :: keyword
      ! Internal variables
      integer :: i

      do i = 1, n_deprecated
         if (trim(keyword) == trim(deprecated_keywords(i))) then
            call print_deprecation_message(keyword, updated_keywords(i))
         end if
      end do
   end subroutine check_deprecated

   subroutine print_deprecation_message(keyword, updated_keyword)
      character*64 :: keyword
      character*64 :: updated_keyword
      integer :: length

      length = len_trim(keyword)

      write (*, *) '.......................................|'
      write (*, *) '                                       |'
      write (*, *) 'WARNING: Found deprecated keyword      |  <-- WARNING'
      write (*, '(A41)') trim(keyword)//' |'
      write (*, *) 'Please replace this keyword with       |'
      write (*, '(A41)') trim(updated_keyword)//' |'
      write (*, *) '                                       |'

   end subroutine print_deprecation_message

!**************************************************************************
!  A keyword that is still parsed so that existing input files keep running,
!  but whose value nothing consults. The field it used to write has been
!  removed from input_parameters rather than left there looking meaningful.
!
!  Deleting the keyword outright is not an option: read_input_file treats an
!  unrecognised keyword as fatal, so every deck that sets one of these would
!  abort. Warning here is the visible half of the same decision -- the
!  generated keyword reference lists these as "accepted and ignored".
   subroutine ignored_keyword(unit, iostatus, rank, keyword)

      implicit none

      integer, intent(in) :: unit, rank
      integer, intent(inout) :: iostatus
      character(len=*), intent(in) :: keyword

      character*64 :: cjunk

      backspace (unit)
      read (unit, *, iostat=iostatus) cjunk, cjunk, cjunk
      call check_iostatus(iostatus, keyword)
      if (rank == 0) write (*, '(1X,A)') 'WARNING: "'//trim(keyword)// &
         '" is accepted for backwards compatibility and ignored'

   end subroutine ignored_keyword

!**************************************************************************
   subroutine read_gap_hypers(file_gap, &
                              n_soap_turbo, soap_turbo_hypers, &
                              n_distance_2b, distance_2b_hypers, &
                              n_angle_3b, angle_3b_hypers, &
                              n_core_pot, core_pot_hypers, &
                              rcut_max, do_prediction, params)

      implicit none

!   Input variables
      type(input_parameters), intent(in) :: params
      logical, intent(in) :: do_prediction
      character(len=*), intent(in) :: file_gap

!   Output variables
      real(dp), intent(out) :: rcut_max
      integer, intent(out) :: n_soap_turbo
      integer, intent(out) :: n_distance_2b
      integer, intent(out) :: n_angle_3b
      integer, intent(out) :: n_core_pot
      integer :: nw
      type(soap_turbo), allocatable, intent(out) :: soap_turbo_hypers(:)
      type(distance_2b), allocatable, intent(out) :: distance_2b_hypers(:)
      type(angle_3b), allocatable, intent(out) :: angle_3b_hypers(:)
      type(core_pot), allocatable, intent(out) :: core_pot_hypers(:)
!   Internal variables
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: V(:)
      real(dp) :: sig
      real(dp) :: p
      real(dp) :: qn
      real(dp) :: un
      integer :: iostatus
      integer :: i
      integer :: counter
      integer :: n_species
      integer :: n_sparse
      integer :: ijunk
      integer :: n
      integer :: n_nonzero
      integer :: j
      character*64 :: keyword
      character*64 :: cjunk
      character*64 :: compress_string
      character*1 :: keyword_first
      integer, parameter :: n_deprecated = 6
      character*64 :: deprecated_keywords(n_deprecated)
      character*64 :: updated_keywords(n_deprecated)

      deprecated_keywords(1) = "has_vdw"
      deprecated_keywords(2) = "vdw_qs"
      deprecated_keywords(3) = "vdw_alphas"
      deprecated_keywords(4) = "vdw_zeta"
      deprecated_keywords(5) = "vdw_delta"
      deprecated_keywords(6) = "vdw_v0"

      updated_keywords(1) = "has_local_properties"
      updated_keywords(2) = "local_property_qs"
      updated_keywords(3) = "local_property_alphas"
      updated_keywords(4) = "local_property_zetas"
      updated_keywords(5) = "local_property_deltas"
      updated_keywords(6) = "local_property_v0s"

      open (unit=10, file=file_gap, status="old", iostat=iostatus)
!   Look for the number of instances of each GAP
      n_soap_turbo = 0
      n_distance_2b = 0
      n_angle_3b = 0
      n_core_pot = 0
      rcut_max = 0.d0
      iostatus = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'gap_beg') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, keyword
            if (keyword == "soap_turbo") then
               n_soap_turbo = n_soap_turbo + 1
            else if (keyword == "distance_2b") then
               n_distance_2b = n_distance_2b + 1
            else if (keyword == "angle_3b") then
               n_angle_3b = n_angle_3b + 1
            else if (keyword == "core_pot") then
               n_core_pot = n_core_pot + 1
            end if
         end if
      end do
!   Allocate the variables
      if (n_soap_turbo > 0) then
         allocate (soap_turbo_hypers(1:n_soap_turbo))
      end if
      if (n_distance_2b > 0) then
         allocate (distance_2b_hypers(1:n_distance_2b))
      end if
      if (n_angle_3b > 0) then
         allocate (angle_3b_hypers(1:n_angle_3b))
      end if
      if (n_core_pot > 0) then
         allocate (core_pot_hypers(1:n_core_pot))
      end if

!   Now record the hypers
      rewind (10)
      n_soap_turbo = 0
      n_distance_2b = 0
      n_angle_3b = 0
      n_core_pot = 0
      iostatus = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'gap_beg') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, keyword
!       soap_turbo definitions here
            if (keyword == "soap_turbo") then
               n_soap_turbo = n_soap_turbo + 1
               counter = 0
               do while (iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  counter = counter + 1
                  !> @kw n_species
                  !> Number of species appearing in the neighbour environments this descriptor
                  !> sees. It must come first in the block: every per-species list below is read
                  !> n_species values wide, and the block is pre-scanned for this keyword alone
                  !> before anything else is parsed so that those arrays can be allocated.
                  !> @see species
                  if (keyword == "n_species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, n_species
                     soap_turbo_hypers(n_soap_turbo)%n_species = n_species
                     allocate (soap_turbo_hypers(n_soap_turbo)%nf(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%nf = 4.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft = 0.5d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling = 0.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling = 0.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%amplitude_scaling = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%central_weight = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%global_scaling = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species))
                     do i = 1, counter
                        backspace (10)
                     end do
                     exit
                  end if
               end do
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  !> @kw nf
                  !> Steepness of the radial filter applied to each species' density, one per
                  !> species. Larger values make the transition sharper.
                  !> @see rcut, buffer
                  if (keyword == "nf") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%nf(1:n_species)
                     !> @kw rcut
                     !> Hard cutoff, one per species: beyond it a neighbour contributes nothing. The
                     !> largest over every descriptor in the file becomes the cutoff the neighbour
                     !> lists are built to.
                     !> @units A
                     !> @see buffer
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species)
                     soap_turbo_hypers(n_soap_turbo)%rcut_max = 0.d0
                     do i = 1, n_species
                        if (soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > soap_turbo_hypers(n_soap_turbo)%rcut_max) then
                           soap_turbo_hypers(n_soap_turbo)%rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                        end if
                     end do
                     !> @kw buffer
                     !> Width of the buffer zone inside the hard cutoff over which a neighbour's
                     !> contribution is smoothly taken to zero, one per species. Read as a width and
                     !> converted to the soft cutoff rcut - buffer once the block closes; zero means
                     !> no buffer, and the density is cut abruptly.
                     !> @units A
                     !> @see rcut
                  else if (keyword == "buffer") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species)
                     !> @kw atom_sigma_r
                     !> Radial width of the Gaussian each neighbour is smeared into, one per species.
                     !> Small values resolve fine radial structure and need a larger n_max to
                     !> represent.
                     !> @units A
                     !> @see atom_sigma_r_scaling, n_max
                  else if (keyword == "atom_sigma_r") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species)
                     !> @kw atom_sigma_t
                     !> Angular width of that Gaussian, measured as an arc length rather than an
                     !> angle, one per species. Divided by the neighbour distance to give the angular
                     !> spread, so a fixed value keeps the angular resolution roughly constant with
                     !> distance.
                     !> @units A
                     !> @see atom_sigma_t_scaling, l_max
                  else if (keyword == "atom_sigma_t") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species)
                     !> @kw atom_sigma_r_scaling
                     !> How the radial width grows with neighbour distance, one per species: the width
                     !> used is atom_sigma_r * (1 + scaling * r). Letting distant neighbours blur is
                     !> what keeps the descriptor from resolving structure the fit cannot support.
                     !> @see atom_sigma_r
                  else if (keyword == "atom_sigma_r_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species)
                     !> @kw atom_sigma_t_scaling
                     !> The same linear growth for the angular width.
                     !> @see atom_sigma_t
                  else if (keyword == "atom_sigma_t_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species)
                     !> @kw amplitude_scaling
                     !> How the amplitude of a neighbour's Gaussian falls with distance, one per
                     !> species, in the form set by scaling_mode. Larger values weight the near
                     !> neighbours more heavily.
                     !> @see scaling_mode
                  else if (keyword == "amplitude_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species)
                     !> @kw central_weight
                     !> Weight given to the central atom's own contribution to the density, one per
                     !> species. Zero leaves the centre out of its own descriptor.
                  else if (keyword == "central_weight") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species)
                     !> @kw global_scaling
                     !> Overall multiplier on each species' contribution to the density, one per
                     !> species. This is where a species is made to count for more or less than its
                     !> neighbours.
                  else if (keyword == "global_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species)
                     !> @kw n_max
                     !> Number of radial basis functions per species. The descriptor's radial
                     !> resolution: it has to be large enough to represent features of width
                     !> atom_sigma_r out to the cutoff. The block's total radial basis size is the sum
                     !> over species.
                     !> @see l_max, atom_sigma_r
                  else if (keyword == "n_max") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species)
!             But we can actually use n_max to referred to the "total" radial basis (the sum of orthogonal bases
!             for different species)
                     soap_turbo_hypers(n_soap_turbo)%n_max = 0
                     do i = 1, n_species
                        soap_turbo_hypers(n_soap_turbo)%n_max = soap_turbo_hypers(n_soap_turbo)%n_max + &
                                                                soap_turbo_hypers(n_soap_turbo)%alpha_max(i)
                     end do
                     !> @kw species
                     !> Chemical symbols of the neighbour species, in the order every per-species list
                     !> in this block is written. n_species entries.
                     !> @see n_species, central_species
                  else if (keyword == "species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species)
                     !> @kw l_max
                     !> Highest spherical-harmonic degree kept in the angular expansion. The
                     !> descriptor's angular resolution, and the counterpart of n_max: the cost of the
                     !> descriptor grows with both.
                     !> @see n_max, atom_sigma_t
                  else if (keyword == "l_max") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%l_max
                     !> @kw radial_enhancement
                     !> Multiply the radial density by r raised to this power before expanding it,
                     !> which weights the outer shells more heavily. 0, 1 and 2 are the usual choices.
                  else if (keyword == "radial_enhancement") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%radial_enhancement
                     if (soap_turbo_hypers(n_soap_turbo)%radial_enhancement < 0 .or. &
                         soap_turbo_hypers(n_soap_turbo)%radial_enhancement > 2) then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: radial_enhancement must be    |  <-- WARNING'
                        write (*, *) 'and 0 <= n <= 2. I am defaulting to 0! |'
                        write (*, *) '                                       |'
                        write (*, *) '.......................................|'
                        soap_turbo_hypers(n_soap_turbo)%radial_enhancement = 0
                     end if
                     !> @kw compress_soap
                     !> Project the descriptor onto a smaller set of components before the kernel is
                     !> evaluated. The SOAP vector grows as n_max^2 l_max, and compression is what
                     !> keeps a many-species descriptor affordable.
                     !> @see compress_mode, file_compress_soap
                  else if (keyword == "compress_soap") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_soap
                     !> @kw file_compress_soap
                     !> File giving the compression explicitly, as the list of components to keep or a
                     !> transformation matrix. The alternative to naming a compress_mode.
                     !> @needs compress_soap
                     !> @see compress_mode
                  else if (keyword == "file_compress_soap") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_compress
                     !> @kw compress_mode
                     !> Named compression scheme to use, resolved into the projection internally. Use
                     !> this or file_compress_soap, not both.
                     !> @needs compress_soap
                     !> @see file_compress_soap
                  else if (keyword == "compress_mode") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_mode
                     !> @kw dipole_model
                     !> This descriptor's fitted scalar is not an energy but a potential whose
                     !> gradient with respect to the central atom's own position is the local dipole.
                     !> Its energy, forces and virial are meaningless and are not summed into the
                     !> totals; only the dipole is taken.
                     !> @see delta
                  else if (keyword == "dipole_model") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%is_dipole_model
                     !> @kw zeta
                     !> Exponent of the SOAP kernel, (q . q')^zeta. Raising it sharpens the kernel and
                     !> makes the fit more local in descriptor space.
                     !> @see delta
                  else if (keyword == "zeta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%zeta
                     !> @kw delta
                     !> Energy scale of this GAP: the standard deviation of the prior on the energy it
                     !> contributes. It sets how much of the total energy this descriptor is allowed
                     !> to account for, and setting it to zero switches the descriptor's contribution
                     !> off without removing the block.
                     !> @units eV
                     !> @see zeta
                  else if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%delta
                     !> @kw central_species
                     !> Which entry of species this descriptor is centred on, as a 1-based index. One
                     !> soap_turbo block per central species is the usual arrangement, each seeing all
                     !> of them as neighbours.
                     !> @see species
                  else if (keyword == "central_species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_species
                     !> @kw basis
                     !> Radial basis to expand the density in. "poly3" is the polynomial basis, and
                     !> "poly3gauss" adds a Gaussian for the central atom.
                     !> @see n_max
                  else if (keyword == "basis") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%basis
                     if (soap_turbo_hypers(n_soap_turbo)%basis /= "poly3" .and. &
                         soap_turbo_hypers(n_soap_turbo)%basis /= "poly3operator" .and. &
                         soap_turbo_hypers(n_soap_turbo)%basis /= "poly3gauss") then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: I didn''t understand your      |  <-- WARNING'
                        write (*, *) 'keywork for basis; defaulting to       |'
                        write (*, *) '"poly3"                                |'
                        soap_turbo_hypers(n_soap_turbo)%basis = "poly3"
                     end if
                     !> @kw scaling_mode
                     !> Functional form amplitude_scaling is applied through. "polynomial" is the
                     !> default and the only form the GPU kernels implement.
                     !> @see amplitude_scaling
                  else if (keyword == "scaling_mode") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%scaling_mode
                     if (soap_turbo_hypers(n_soap_turbo)%scaling_mode /= "polynomial") then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: I didn''t understand your     |  <-- WARNING'
                        write (*, *) 'keywork for scaling_mode; defaulting   |'
                        write (*, *) 'to "polynomial"                        |'
                        soap_turbo_hypers(n_soap_turbo)%scaling_mode = "polynomial"
                     end if
                     !> @kw alphas_sparse
                     !> File holding the fitted sparse-set coefficients, one per sparse point. Read
                     !> together with desc_sparse, which must have the same number of rows.
                     !> @see desc_sparse
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_alphas
                     !> @kw desc_sparse
                     !> File holding the sparse-set descriptors, one row per sparse point. Its row
                     !> count sets n_sparse and so the cost of every prediction.
                     !> @see alphas_sparse
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_desc
                     !> @kw has_vdw
                     !> Superseded by has_local_properties. A Hirshfeld-volume model is now one local
                     !> property among several rather than a special case, and this keyword still sets
                     !> up that one property so that older potential files load.
                     !> @see has_local_properties
                     !> @status deprecated
                  else if (keyword == "has_vdw") then
                     backspace (10)
                     !               read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_vdw
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_local_properties
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)

                     write (*, *) '---------------------------------------|'
                     write (*, *) '--------------vdW Notice---------------|'
                     write (*, *) '---------------------------------------|'
                     write (*, *) '                                       |'
                     write (*, *) 'When upgrading from deprecated vdW,    |'
                     write (*, *) 'set in *.gap file (for just one model) |'
                     write (*, *) '`n_local_properties = 1`               |'
                     write (*, *) '`local_property_labels = "hirshfeld_v"`|'
                     write (*, *) '                                       |'
                     write (*, *) '---------------------------------------|'
                     write (*, *) '                                       |'
                     write (*, *) 'WARNING: Defaulting to just 1          |  <-- WARNING'
                     write (*, *) 'local property due to deprecated       |'
                     write (*, *) 'keyword. Other loc models will not run |'
                     write (*, *) '                                       |'
                     write (*, *) '---------------------------------------|'

                     soap_turbo_hypers(n_soap_turbo)%n_local_properties = 1

                     allocate (soap_turbo_hypers(n_soap_turbo)%local_property_models( &
                               1:soap_turbo_hypers(n_soap_turbo)%n_local_properties))

                     soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%label = 'hirshfeld_v'

                     soap_turbo_hypers(n_soap_turbo)%has_vdw = .true.
                     soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%do_derivatives = .true.
                     soap_turbo_hypers(n_soap_turbo)%vdw_index = 1

                     !> @kw vdw_qs
                     !> Superseded by local_property_qs.
                     !> @see local_property_qs
                     !> @status deprecated
                  else if (keyword == "vdw_qs") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !             read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_desc
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%file_desc
                     !> @kw vdw_alphas
                     !> Superseded by local_property_alphas.
                     !> @see local_property_alphas
                     !> @status deprecated
                  else if (keyword == "vdw_alphas") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas
                    read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%file_alphas
                     !> @kw vdw_zeta
                     !> Superseded by local_property_zetas.
                     !> @see local_property_zetas
                     !> @status deprecated
                  else if (keyword == "vdw_zeta") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_zeta
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%zeta
                     !> @kw vdw_delta
                     !> Superseded by local_property_deltas.
                     !> @see local_property_deltas
                     !> @status deprecated
                  else if (keyword == "vdw_delta") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_delta
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%delta
                     !> @kw vdw_v0
                     !> Superseded by local_property_v0s.
                     !> @see local_property_v0s
                     !> @status deprecated
                  else if (keyword == "vdw_v0") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_v0
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%V0
                     !> @kw has_local_properties
                     !> This descriptor carries per-atom quantities fitted alongside the energy --
                     !> Hirshfeld volumes, charges, core-electron binding energies. Set implicitly by
                     !> n_local_properties.
                     !> @see n_local_properties, local_property_labels
                  else if (keyword == "has_local_properties") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_local_properties
                     !> @kw n_local_properties
                     !> How many local properties this descriptor carries. Give it before the
                     !> local_property_* lists, which it allocates and which are all read this many
                     !> values wide.
                     !> @see local_property_labels
                  else if (keyword == "n_local_properties") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                     ! Now allocate the local_property_soap_turbo object in the soap_turbo_hypers
                     allocate (soap_turbo_hypers(n_soap_turbo)%local_property_models( &
                               1:soap_turbo_hypers(n_soap_turbo)%n_local_properties))

                     !> @kw local_property_labels
                     !> What each local property is, one name per property. The names are meaningful,
                     !> not decorative: "hirshfeld_v" is what the van der Waals correction looks for,
                     !> and "core_electron_be" is what the XPS spectrum is built from. Naming either
                     !> one switches on the machinery that consumes it and records which slot it
                     !> occupies.
                     !> @needs n_local_properties
                     !> @see local_property_qs, local_property_alphas
                  else if (keyword == "local_property_labels") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%label, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                     do nw = 1, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                        if (trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                             &%label) == "hirshfeld_v") then
                           soap_turbo_hypers(n_soap_turbo)%has_vdw = .true.
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)%do_derivatives = .true.
                           soap_turbo_hypers(n_soap_turbo)%vdw_index = nw
                        end if

                        if (trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                             &%label) == "core_electron_be") then
                           soap_turbo_hypers(n_soap_turbo)%has_core_electron_be = .true.
                           soap_turbo_hypers(n_soap_turbo)%core_electron_be_index = nw
                        end if

                     end do

                     !> @kw local_property_qs
                     !> Sparse-set descriptor file for each local property, one per property.
                     !> @needs n_local_properties
                     !> @see local_property_alphas
                  else if (keyword == "local_property_qs") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%file_desc, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                     !> @kw local_property_alphas
                     !> Sparse-set coefficient file for each local property, one per property.
                     !> @needs n_local_properties
                     !> @see local_property_qs
                  else if (keyword == "local_property_alphas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%file_alphas, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                     !> @kw local_property_zetas
                     !> Kernel exponent for each local property, one per property. The local
                     !> properties are fitted with their own kernels and need not share the energy's
                     !> zeta.
                     !> @needs n_local_properties
                     !> @see zeta
                  else if (keyword == "local_property_zetas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk,&
                          & (soap_turbo_hypers(n_soap_turbo)&
                          &%local_property_models(nw)%zeta, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                     !> @kw local_property_deltas
                     !> Prior scale for each local property, one per property, in whatever units that
                     !> property has.
                     !> @needs n_local_properties
                     !> @see delta
                  else if (keyword == "local_property_deltas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%delta, nw=1, soap_turbo_hypers(n_soap_turbo)&
                          &%n_local_properties)
                     !> @kw local_property_v0s
                     !> Baseline added to each local property's prediction, one per property. For
                     !> Hirshfeld volumes this is the free-atom volume the fit is a correction to.
                     !> @needs n_local_properties
                  else if (keyword == "local_property_v0s") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%V0, nw=1, soap_turbo_hypers(n_soap_turbo)&
                          &%n_local_properties)
                  end if
               end do
!         We actually read in the "buffer" zone width, so transform to rcut_soft:
               do i = 1, n_species
                  if (soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) == 0.d0) then
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                  else
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) - &
                                                                    soap_turbo_hypers(n_soap_turbo)%rcut_soft(i)
                  end if
               end do
!         Read the sparse set information
               if (do_prediction) then
                  call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_desc, &
                                                   soap_turbo_hypers(n_soap_turbo)%file_alphas, &
                                                   soap_turbo_hypers(n_soap_turbo)%n_sparse, &
                                                   "soap_turbo", soap_turbo_hypers(n_soap_turbo)%alphas, &
                                                   soap_turbo_hypers(n_soap_turbo)%Qs, &
                                                   soap_turbo_hypers(n_soap_turbo)%cutoff)
                  ! Commenting this out as it will be subsumed into local property prediction
                  ! if( soap_turbo_hypers(n_soap_turbo)%has_vdw )then
                  !    call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_vdw_desc, &
                  !         soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_n_sparse, &
                  !         "soap_turbo", soap_turbo_hypers(n_soap_turbo)%vdw_alphas, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_Qs, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_cutoff)

                  ! end if

                  if (soap_turbo_hypers(n_soap_turbo)%has_local_properties) then
                     do j = 1, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                        call read_alphas_and_descriptors( &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_desc, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_alphas, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%n_sparse, &
                           "soap_turbo", &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%alphas, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%Qs, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%cutoff)

                        ! Really, this could actually just not be associated
                        ! with the soap turbo type as the same data might be
                        ! reread into separate soap turbo descriptors when
                        ! only one is needed, and further this is
                        ! broadcasted. But this way, all the files are
                        ! specified in the .gap file rather than in the
                        ! input file.

                        ! soap_turbo_hypers(n_soap_turbo)&
                        !      &%local_property_models(j)%dim =&
                        !      & size(soap_turbo_hypers(n_soap_turbo)&
                        !      &%local_property_models(j)%Qs,1)

                     end do
                  end if

               end if
               do i = 1, n_species
                  if (soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > rcut_max) then
                     rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                  end if
               end do
!         Handle SOAP compression here
!         Here we read in the compression information from a file (compress_file) or rely on a keyword provided
!         by the user (compress_mode) which leads to a predefined recipe to compress the soap_turbo descriptor
!         The file always takes precedence over the keyword.
               if (soap_turbo_hypers(n_soap_turbo)%compress_soap) then
!           A compress file takes priority over compress mode
                  if (soap_turbo_hypers(n_soap_turbo)%file_compress /= "none") then
                     open (unit=20, file=soap_turbo_hypers(n_soap_turbo)%file_compress, status="old")
                     read (20, *) (ijunk, i=1, n_species), ijunk, soap_turbo_hypers(n_soap_turbo)%dim
!             This enables definition of arbitrary compression transformations via a file
                     read (20, '(A)') compress_string
                     if (compress_string == "P_transformation") then
                        n_nonzero = -1
                        do while (compress_string /= "end_transformation")
                           read (20, '(A)') compress_string
                           n_nonzero = n_nonzero + 1
                        end do
                        soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = n_nonzero
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:n_nonzero))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:n_nonzero))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:n_nonzero))
                        do i = 1, n_nonzero + 1
                           backspace (20)
                        end do
                        do i = 1, n_nonzero
                           read (20, *) soap_turbo_hypers(n_soap_turbo)%compress_P_i(i), &
                              soap_turbo_hypers(n_soap_turbo)%compress_P_j(i), &
                              soap_turbo_hypers(n_soap_turbo)%compress_P_el(i)
                        end do
                     else
!               Old way to handle compression for backcompatibility
                        backspace (20)
                        soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = soap_turbo_hypers(n_soap_turbo)%dim
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        do i = 1, soap_turbo_hypers(n_soap_turbo)%dim
                           read (20, *) soap_turbo_hypers(n_soap_turbo)%compress_P_j(i)
                           soap_turbo_hypers(n_soap_turbo)%compress_P_i(i) = i
                           soap_turbo_hypers(n_soap_turbo)%compress_P_el(i) = 1.d0
                        end do
                     end if
                     close (20)
                  else if (soap_turbo_hypers(n_soap_turbo)%compress_mode /= "none") then
                     call get_compress_indices(soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                               soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                               soap_turbo_hypers(n_soap_turbo)%l_max, &
                                               soap_turbo_hypers(n_soap_turbo)%dim, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                               "get_dim")
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     call get_compress_indices(soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                               soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                               soap_turbo_hypers(n_soap_turbo)%l_max, &
                                               soap_turbo_hypers(n_soap_turbo)%dim, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                               "set_indices")
                  else
                     write (*, *) "ERROR: you're trying to use compression but neither a file_compress_soap nor", &
                        "compress_mode are defined!"
                     stop
                  end if
               else
                  soap_turbo_hypers(n_soap_turbo)%dim = soap_turbo_hypers(n_soap_turbo)%n_max* &
                                                        (soap_turbo_hypers(n_soap_turbo)%n_max + 1)/2* &
                                                        (soap_turbo_hypers(n_soap_turbo)%l_max + 1)
               end if
!       distance_2b definitions here
            else if (keyword == "distance_2b") then
               n_distance_2b = n_distance_2b + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  !> @kw delta
                  !> Energy scale of this two-body GAP, the standard deviation of the prior on the
                  !> energy it contributes. Zero switches the block off without removing it.
                  !> @units eV
                  !> @see sigma
                  if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%delta
                     !> @kw sigma
                     !> Width of the kernel in the descriptor, which here is the interatomic distance
                     !> itself. It sets how far apart two distances have to be before the fit treats
                     !> them as different.
                     !> @units A
                     !> @see delta
                  else if (keyword == "sigma") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%sigma
                     !> @kw rcut
                     !> Cutoff of the pair interaction: beyond it the pair contributes nothing.
                     !> @units A
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%rcut
                     !> @kw Z1, z1, species1
                     !> First species of the pair this block describes, as a chemical symbol. With Z2
                     !> it decides which pairs the descriptor is evaluated on.
                     !> @see Z2
                  else if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species1
                     !> @kw Z2, z2, species2
                     !> Second species of the pair. A block with Z1 and Z2 swapped would be the same
                     !> interaction, so only one of the two orders is given.
                     !> @see Z1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species2
                     !> @kw desc_sparse
                     !> File holding the sparse-set descriptors, one distance per row. Its row count
                     !> sets the number of sparse points.
                     !> @see alphas_sparse
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_desc
                     !> @kw alphas_sparse
                     !> File holding the fitted sparse-set coefficients, one per sparse point.
                     !> @see desc_sparse
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_alphas
                  end if
               end do
!         Read the sparse set information
               call read_alphas_and_descriptors(distance_2b_hypers(n_distance_2b)%file_desc, &
                                                distance_2b_hypers(n_distance_2b)%file_alphas, &
                                                distance_2b_hypers(n_distance_2b)%n_sparse, &
                                                "distance_2b", distance_2b_hypers(n_distance_2b)%alphas, &
                                                distance_2b_hypers(n_distance_2b)%Qs, &
                                                distance_2b_hypers(n_distance_2b)%cutoff)
               if (distance_2b_hypers(n_distance_2b)%rcut > rcut_max) then
                  rcut_max = distance_2b_hypers(n_distance_2b)%rcut
               end if
!       angle_3b definitions here
            else if (keyword == "angle_3b") then
               n_angle_3b = n_angle_3b + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  !> @kw delta
                  !> Energy scale of this three-body GAP. Zero switches the block off without
                  !> removing it.
                  !> @units eV
                  !> @see sigma
                  if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%delta
                     !> @kw sigma
                     !> Kernel widths of the three-body descriptor, three values: one for each of the
                     !> two distances and one for the coordinate built from the angle between them.
                     !> @see delta, kernel_type
                  else if (keyword == "sigma") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%sigma(1:3)
                     !> @kw rcut
                     !> Cutoff on each of the two distances from the centre. Both neighbours must be
                     !> inside it for the triplet to contribute.
                     !> @units A
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%rcut
                     !> @kw Z1, z1, species1
                     !> Species of the first neighbour in the triplet.
                     !> @see Z_center, Z2
                  else if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species1
                     !> @kw Z2, z2, species2
                     !> Species of the second neighbour. Swapping Z1 and Z2 gives the same triplet, so
                     !> only one order is given and the descriptor is symmetrised internally.
                     !> @see Z1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species2
                     !> @kw Z_center, z_center, species_center
                     !> Species at the vertex of the triplet, as a chemical symbol. The three-body
                     !> term is centred on this atom with Z1 and Z2 as its two neighbours.
                     !> @see Z1, Z2
                  else if (keyword == "Z_center" .or. keyword == "z_center" .or. keyword == "species_center") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species_center
                     !> @kw kernel_type
                     !> Form of the three-body kernel. "exp" is the squared-exponential the fits use.
                     !> @see sigma
                  else if (keyword == "kernel_type") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%kernel_type
                     !> @kw desc_sparse
                     !> File holding the sparse-set descriptors, one triplet per row, three components
                     !> each.
                     !> @see alphas_sparse
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_desc
                     !> @kw alphas_sparse
                     !> File holding the fitted sparse-set coefficients, one per sparse point.
                     !> @see desc_sparse
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_alphas
                  end if
               end do
!         Read the sparse set information
               call read_alphas_and_descriptors(angle_3b_hypers(n_angle_3b)%file_desc, &
                                                angle_3b_hypers(n_angle_3b)%file_alphas, &
                                                angle_3b_hypers(n_angle_3b)%n_sparse, &
                                                "angle_3b", angle_3b_hypers(n_angle_3b)%alphas, &
                                                angle_3b_hypers(n_angle_3b)%Qs, &
                                                angle_3b_hypers(n_angle_3b)%cutoff)
               if (angle_3b_hypers(n_angle_3b)%rcut > rcut_max) then
                  rcut_max = angle_3b_hypers(n_angle_3b)%rcut
               end if
!       core_pot definitions here
            else if (keyword == "core_pot") then
               n_core_pot = n_core_pot + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
!            if( keyword == "n" .or. keyword == "N" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%n
!            else if( keyword == "yp1" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%yp1
!            else if( keyword == "ypn" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%ypn
!            else if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
                  !> @kw Z1, z1, species1
                  !> First species of the pair this core potential acts on, as a chemical symbol.
                  !> @see Z2
                  if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species1
                     !> @kw Z2, z2, species2
                     !> Second species of the pair.
                     !> @see Z1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species2
                     !> @kw core_pot_file
                     !> File holding the tabulated potential, as distance and energy columns. It is
                     !> splined on read, and truncated at core_pot_cutoff with a taper of width
                     !> core_pot_buffer, both set in the input file rather than here. This is the
                     !> short-range repulsion the GAP is not fitted to describe, so the table usually
                     !> only covers distances the training data never visited.
                     !> @see core_pot_cutoff, core_pot_buffer
                  else if (keyword == "core_pot_file") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%core_pot_file
                  end if
               end do
!         Allocate some arrays, read in potential, etc.
               open (20, file=core_pot_hypers(n_core_pot)%core_pot_file, status="unknown")
               read (20, *) core_pot_hypers(n_core_pot)%n, core_pot_hypers(n_core_pot)%yp1, core_pot_hypers(n_core_pot)%ypn
               n = core_pot_hypers(n_core_pot)%n

               allocate (V(1:n))
               allocate (x(1:n))
               counter = 0
               do i = 1, n
                  read (20, *) x(i), V(i)
                  if (x(i) <= params%core_pot_cutoff) then
                     counter = counter + 1
                     x(counter) = x(i)
                     if (x(i) <= params%core_pot_cutoff - params%core_pot_buffer) then
                        V(counter) = V(i)
                     else
                        V(counter) = V(i)*0.5d0*(dcos(dacos(-1.d0)/params%core_pot_buffer*(x(i) - params%core_pot_cutoff &
                                                                                           + params%core_pot_buffer)) + 1.d0)
                     end if
                  end if
               end do
               close (20)
               n = counter
               core_pot_hypers(n_core_pot)%n = n
               allocate (core_pot_hypers(n_core_pot)%V(1:n))
               core_pot_hypers(n_core_pot)%V(1:n) = V(1:n)
               deallocate (V)
               allocate (core_pot_hypers(n_core_pot)%x(1:n))
               core_pot_hypers(n_core_pot)%x(1:n) = x(1:n)
               deallocate (x)
               allocate (core_pot_hypers(n_core_pot)%dVdx2(1:n))
!         This code below for spline second derivative is more or less copy-pasted from QUIP. It's the
!         easiest way to make sure both interpolations give the same numbers
               allocate (u(1:n))
               if (core_pot_hypers(n_core_pot)%yp1 > 0.99d30) then
                  core_pot_hypers(n_core_pot)%dVdx2(1) = 0.d0
                  u(1) = 0.d0
               else
                  core_pot_hypers(n_core_pot)%dVdx2(1) = -0.5d0

                  u(1) = (3.d0/(core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1)))* &
                         ((core_pot_hypers(n_core_pot)%V(2) - core_pot_hypers(n_core_pot)%V(1))/ &
                          (core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1)) - core_pot_hypers(n_core_pot)%yp1)
               end if
               do i = 2, n - 1
                  sig = (core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i - 1))/ &
                        (core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i - 1))
                  p = sig*core_pot_hypers(n_core_pot)%dVdx2(i - 1) + 2.d0
                  core_pot_hypers(n_core_pot)%dVdx2(i) = (sig - 1.d0)/p
                  u(i) = (core_pot_hypers(n_core_pot)%V(i + 1) - core_pot_hypers(n_core_pot)%V(i))/ &
                         (core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i)) - &
                         (core_pot_hypers(n_core_pot)%V(i) - core_pot_hypers(n_core_pot)%V(i - 1))/ &
                         (core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i - 1))
                  u(i) = (6.d0*u(i)/(core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i - 1)) &
                          - sig*u(i - 1))/p
               end do
               if (core_pot_hypers(n_core_pot)%ypn > 0.99d30) then
                  qn = 0.d0
                  un = 0.d0
               else
                  qn = 0.5d0
                  un = (3.d0/(core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n - 1)))* &
                      (core_pot_hypers(n_core_pot)%ypn - (core_pot_hypers(n_core_pot)%V(n) - core_pot_hypers(n_core_pot)%V(n - 1)) &
                        /(core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n - 1)))
               end if
               core_pot_hypers(n_core_pot)%dVdx2(n) = (un - qn*u(n - 1))/(qn*core_pot_hypers(n_core_pot)%dVdx2(n - 1) + 1.d0)
               do i = n - 1, 1, -1
                  core_pot_hypers(n_core_pot)%dVdx2(i) = core_pot_hypers(n_core_pot)%dVdx2(i)* &
                                                         core_pot_hypers(n_core_pot)%dVdx2(i + 1) + u(i)
               end do
               if (core_pot_hypers(n_core_pot)%yp1 > 0.99d30) then
                  u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                      core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                      core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(1:1), 1.d100)
                  core_pot_hypers(n_core_pot)%yp1 = u(1)
               end if
               if (core_pot_hypers(n_core_pot)%ypn > 0.99d30) then
                  u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                      core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                      core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(n:n), 1.d100)
                  core_pot_hypers(n_core_pot)%ypn = u(1)
               end if
               deallocate (u)
               if (maxval(core_pot_hypers(n_core_pot)%x) > rcut_max) then
                  rcut_max = maxval(core_pot_hypers(n_core_pot)%x)
               end if
            end if
         end if
      end do
      close (10)

   end subroutine
!**************************************************************************

!**************************************************************************
! ------- option for radiation cascade simulation with electronic stopping

   subroutine read_electronic_stopping_file(n_species, species_types, estopfilename, nrows, allelstopdata)
! read the given electronic stopping file
! send the data for required calculations
! also give error messages if the data in the file is not in proper format
      implicit none

      character*1024, intent(in) :: estopfilename
      integer, intent(in) :: n_species
      character*8, intent(in) :: species_types(n_species)
      integer, intent(out) :: nrows
      real(dp), allocatable :: allelstopdata(:)
      character*8, allocatable :: infoline(:)
      integer :: i
      integer :: ncols
      integer :: ndata

      open (unit=1000, file=estopfilename)
! first line gives information
! second line gives number of energy-stopping data points, i.e no. of rows of data
      read (1000, *)
      read (1000, *) nrows
      if (nrows <= 0) then
         write (*, *) "ERROR: Number of data rows in stopping file is 0 or less."
         stop
      end if
      ncols = n_species + 1
      allocate (infoline(ncols))
! third line gives energy units, names of elements in order of the atom species types in input file
      read (1000, *) (infoline(i), i=1, ncols)
      do i = 2, ncols
         if (trim(infoline(i)) /= trim(species_types(i - 1))) then
            write (*, *) "ERROR: Stopping powers for Elements are not given in order."
            stop
         end if
      end do
      ndata = nrows*ncols
      allocate (allelstopdata(ndata))
      read (1000, *) (allelstopdata(i), i=1, ndata)

      close (unit=1000)
   end subroutine read_electronic_stopping_file
!**************************************************************************

!**************************************************************************
   subroutine upper_to_lower_case(string)

      implicit none

      character(len=*), intent(inout) :: string
      character*1 :: upper_case_dict(1:26)
      character*1 :: lower_case_dict(1:26)
      integer :: i
      integer :: j
      integer :: n

      upper_case_dict = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
                         "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
                         "U", "V", "W", "X", "Y", "Z"]
      lower_case_dict = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", &
                         "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", &
                         "u", "v", "w", "x", "y", "z"]

      do i = 1, len(string)
         do j = 1, size(upper_case_dict)
            if (string(i:i) == upper_case_dict(j)) then
               string(i:i) = lower_case_dict(j)
            end if
         end do
      end do

   end subroutine
!**************************************************************************

   subroutine get_irreducible_local_properties(params, n_local_properties_tot, n_soap_turbo, soap_turbo_hypers, &
                           local_property_labels, local_property_labels_temp, local_property_labels_temp2, local_property_indexes, &
                                valid_vdw, vdw_lp_index, valid_estat_charges, charge_lp_index, core_be_lp_index, valid_xps, xps_idx)
      implicit none
      type(input_parameters), intent(inout) :: params
      integer, intent(in) :: n_soap_turbo
      integer, intent(inout) :: n_local_properties_tot
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)
      character*1024, allocatable, intent(inout) :: local_property_labels(:)
      character*1024, allocatable, intent(inout) :: local_property_labels_temp(:)
      character*1024, allocatable, intent(inout) :: local_property_labels_temp2(:)
      integer, allocatable, intent(inout) :: local_property_indexes(:)
      integer, intent(inout) :: vdw_lp_index
      integer, intent(inout) :: core_be_lp_index
      integer, intent(inout) :: charge_lp_index
      integer, intent(inout) :: xps_idx
      logical, intent(inout) :: valid_vdw
      logical, intent(inout) :: valid_xps
      logical, intent(inout) :: valid_estat_charges
      logical :: label_in_list = .false.
      integer :: i
      integer :: j
      integer :: i2
      integer :: j2
      integer :: k
      integer :: k2
      integer :: nprop
      integer :: length

      n_local_properties_tot = 0
      i2 = 1 ! using this as a counter for the labels
      do j = 1, n_soap_turbo
         if (soap_turbo_hypers(j)%has_local_properties) then
            ! This property has the labels of the quantities to
            ! compute. We must specify the number of local properties, for the sake of coding simplicity

            n_local_properties_tot = n_local_properties_tot + soap_turbo_hypers(j)%n_local_properties

            if (.not. allocated(local_property_labels)) then
               allocate (local_property_labels(1:n_local_properties_tot))
               do i = 1, n_local_properties_tot
                  local_property_labels(i) = soap_turbo_hypers(j)%local_property_models(i)%label
                  length = len_trim(soap_turbo_hypers(j)%local_property_models(i)%label)
                  write (*, *) ' Local property found                  |'
                  write (*, '(A,1X,I8,1X,A20)') ' Descriptor', j,&
                       & trim(soap_turbo_hypers(j)&
                       &%local_property_models(i)%label)//'|'
               end do
            else
               ! Allocate temporary array which is of the size before
               allocate (local_property_labels_temp(1:n_local_properties_tot - soap_turbo_hypers(j)%n_local_properties))
               local_property_labels_temp = local_property_labels
               deallocate (local_property_labels)
               allocate (local_property_labels(1:n_local_properties_tot))

               nprop = soap_turbo_hypers(j)%n_local_properties
               do i = 1, n_local_properties_tot - nprop
                  local_property_labels(i) = local_property_labels_temp(i)
               end do

               deallocate (local_property_labels_temp)

               do i = 1, nprop
                  local_property_labels(i + n_local_properties_tot -&
                       & nprop) = soap_turbo_hypers(j)&
                       &%local_property_models(i)%label
                  write (*, *) ' Local property found                  |'
                  write (*, '(A,1X,I8,1X,A20)') ' Descriptor ', j,&
                       & trim(soap_turbo_hypers(j)&
                       &%local_property_models(i)%label)//'|'

               end do
            end if
         end if
      end do

      ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties

      ! Now we create an irreducible list of the labels
      i2 = 0
      if (n_local_properties_tot > 0) then
         allocate (local_property_labels_temp(1:1))
         local_property_labels_temp(1) = local_property_labels(1)
         i2 = 1
         if (n_local_properties_tot > 1) then
            do i = 2, n_local_properties_tot
               label_in_list = .false.
               ! Iterate through irreducible list to see if there is a mismatch
               do j = 1, size(local_property_labels_temp, 1)
                  if (trim(local_property_labels_temp(j)) == trim(local_property_labels(i))) label_in_list = .true.
               end do
               if (.not. label_in_list) then
                  i2 = i2 + 1
                  allocate (local_property_labels_temp2(1:i2))
                  local_property_labels_temp2(1:i2 - 1) = local_property_labels_temp(1:i2 - 1)
                  local_property_labels_temp2(i2) = local_property_labels(i)
                  deallocate (local_property_labels_temp)
                  allocate (local_property_labels_temp(1:i2))
                  local_property_labels_temp(1:i2) = local_property_labels_temp2(1:i2)
                  deallocate (local_property_labels_temp2)
               end if
            end do
         end if

         params%n_local_properties = i2

         ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties
         !                local_property_labels_temp( 1:params%n_local_properties ) has irreducible labels of local properties

         ! Now we can have an array which has a soap turbo index as an input and it can give us the corresponding label
         allocate (local_property_indexes(1:n_local_properties_tot))
         i2 = 1
         do i = 1, params%n_local_properties
            do j = 1, n_local_properties_tot
               if (trim(local_property_labels(j)) == trim(local_property_labels_temp(i))) then

                  local_property_indexes(j) = i

                  if (trim(local_property_labels(j)) == "hirshfeld_v") then
                     vdw_lp_index = i
                     valid_vdw = .true.
                     do k2 = 1, n_soap_turbo
                        do k = 1, soap_turbo_hypers(k2)%n_local_properties
                           if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == "hirshfeld_v") then
                              if (params%do_derivatives .or. params%do_forces) then
                                 soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
                              else
                                 soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
                              end if

                           end if
                        end do
                     end do
                  end if

                  if (trim(local_property_labels(j)) == "atomic_charge") then
                     charge_lp_index = i
                     valid_estat_charges = .true.
                     do k2 = 1, n_soap_turbo
                        do k = 1, soap_turbo_hypers(k2)%n_local_properties
                           if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == "atomic_charge") then
                              soap_turbo_hypers(k2)%has_charges = .true.
                              if (params%do_derivatives .or. params%do_forces) then
                                 soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
                              else
                                 soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
                              end if
                              ! This is important -- no truncating the charges; it doesn't make sense here!
                           end if
                        end do
                     end do
                  end if

                  if (trim(local_property_labels(j)) == "core_electron_be") then
                     core_be_lp_index = i

                     ! Check if there is experimental data for one to do xps fitting
                     do i2 = 1, params%n_exp
                        if ((trim(params%exp_data(i2)%label) == "xps" .and. &
                             .not. (trim(params%exp_data(i2)%file_data) == "none"))) then
                           valid_xps = .true.
                           xps_idx = i2
                           do k2 = 1, n_soap_turbo
                              do k = 1, soap_turbo_hypers(k2)%n_local_properties
                                 if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == "core_electron_be") then
                                    soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
                                    if (params%exp_forces .and. params%do_derivatives) then
                                       soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
                                    end if
                                 end if
                              end do
                           end do
                        end if
                     end do
                  end if

               end if
            end do
         end do

         deallocate (local_property_labels)
         allocate (local_property_labels(1:size(local_property_labels_temp, 1)))
         local_property_labels = local_property_labels_temp
         deallocate (local_property_labels_temp)
      end if

   end subroutine get_irreducible_local_properties

!**************************************************************************
   subroutine init_random_seed(seed_value, rank)

      implicit none

      integer, intent(in) :: seed_value
      integer, intent(in) :: rank
      integer, allocatable :: seed(:)
      integer :: n
      integer :: i

      if (seed_value == 0) return

      call random_seed(size=n)
      allocate (seed(1:n))
      do i = 1, n
         seed(i) = seed_value + 37*(i - 1) + 7919*rank
      end do
      call random_seed(put=seed)
      deallocate (seed)

   end subroutine init_random_seed
!**************************************************************************

end module
