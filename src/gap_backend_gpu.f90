! Copyright (c) 2020-2026 by Albert Bartok and Miguel Caro
!
! The two-body, core-potential and three-body contributions: the GPU
! implementation of the gap_backend seam. The CPU implementation is
! gap_backend_cpu.f90 on the master branch. Same module name, same three
! contribution procedures, same argument lists; the Makefile compiles one.
!
! The point of this file is what is NOT in its interface. These three
! computations were internal procedures of turbogap.f90 reaching ten device
! buffers by host association -- rjs_d, xyz_d, species_d, neighbors_list_d,
! neighbor_species_d, n_neigh_d, alphas_d, qs_d, cutoff_d and the stream. The
! driver allocated them, the procedures used them, the driver freed them, and
! none of that could be expressed in an interface a CPU build could also
! satisfy. They are module state here, so the driver never names a device
! pointer and the argument lists match the CPU side exactly.
!
! gap_backend_begin uploads the neighbour data once for all three calls, which
! is what the driver used to do inline; gap_backend_end frees it. On the CPU
! both are empty.
!
! alphas_d is declared here rather than shared with the driver on purpose. The
! driver reuses one alphas_d for the SOAP path and the 2b/3b path -- separate
! allocate/upload/free cycles holding different data through one name. That is
! the reuse pattern that hid bug 5 (handoff section 4); the backend owning its
! own removes one instance of it.
module gap_backend

  use kinds

  use types
  use gap
  use timing
  use F_B_C
  use iso_c_binding

  implicit none

  private
 public :: gap_backend_gpu_init
 public :: gap_backend_begin
 public :: gap_backend_end
 public :: add_2b_contribution
 public :: add_core_pot_contribution
 public :: add_3b_contribution

!   Device state owned by this module. Named as in the driver so the lifted
!   bodies below are unchanged.
 type(c_ptr) :: n_neigh_d
 type(c_ptr) :: species_d
 type(c_ptr) :: neighbor_species_d
 type(c_ptr) :: rjs_d
 type(c_ptr) :: xyz_d
 type(c_ptr) :: neighbors_list_d
 type(c_ptr) :: alphas_d
 type(c_ptr) :: cutoff_d
 type(c_ptr) :: qs_d
 type(c_ptr) :: gpu_stream

contains

!**************************************************************************
! Hand the backend the stream it should launch on. GPU-only, and called from
! the driver's existing GPU start-up block, which is #ifdef-guarded anyway --
! the shared contribution interface stays free of device handles.
  subroutine gap_backend_gpu_init( stream )
    implicit none
 type(c_ptr), intent(in) :: stream
    gpu_stream = stream
  end subroutine gap_backend_gpu_init

!**************************************************************************
! Upload the neighbour data once for all three contribution calls.
  subroutine gap_backend_begin( params, rjs, xyz, n_neigh, species, neighbor_species, &
       neighbors_list, i_beg, i_end, j_beg, j_end )

    implicit none
 type(input_parameters), intent(inout) :: params
 real(dp), intent(in), allocatable, target :: rjs(:)
 real(dp), intent(in), allocatable, target :: xyz(:,:)
 integer, intent(in), allocatable, target :: n_neigh(:)
 integer, intent(in), allocatable, target :: species(:)
 integer, intent(in), allocatable, target :: neighbor_species(:)
 integer, intent(in), allocatable, target :: neighbors_list(:)
 integer, intent(in) :: i_beg
 integer, intent(in) :: i_end
 integer, intent(in) :: j_beg
 integer, intent(in) :: j_end

 integer(c_size_t) :: st_n_sites_int
 integer(c_size_t) :: st_n_atom_pairs_int
 integer(c_size_t) :: st_n_atom_pairs_double
 integer(c_size_t) :: size_maxnp_bytes
 integer :: n_sites

    n_sites = i_end - i_beg + 1

    st_n_sites_int = n_sites*sizeof(n_neigh(1)) ! (i_end - i_beg + 1)
    !        print *, rank, " >> Allocating 2b on gpu"        
    call gpu_malloc_all(n_neigh_d,st_n_sites_int,gpu_stream)
    call cpy_htod(c_loc(n_neigh),n_neigh_d, st_n_sites_int,gpu_stream)
    call gpu_malloc_all(species_d,st_n_sites_int,gpu_stream)
    call cpy_htod(c_loc(species),species_d, st_n_sites_int,gpu_stream)
    ! Changed n_atom_pairs to n_atom_pairs_by_rank
    st_n_atom_pairs_int = j_end * sizeof(neighbor_species(1))
    call gpu_malloc_all(neighbor_species_d,st_n_atom_pairs_int,gpu_stream)
    call cpy_htod(c_loc(neighbor_species),neighbor_species_d, st_n_atom_pairs_int,gpu_stream)
    st_n_atom_pairs_double = j_end*sizeof(rjs(1))
    call gpu_malloc_all(rjs_d,st_n_atom_pairs_double,gpu_stream)
    call cpy_htod(c_loc(rjs),rjs_d, st_n_atom_pairs_double,gpu_stream)
    call gpu_malloc_all(xyz_d,3*st_n_atom_pairs_double,gpu_stream)
    call cpy_htod(c_loc(xyz),xyz_d,3*st_n_atom_pairs_double,gpu_stream)

    size_maxnp_bytes = size(neighbors_list) * c_int 
    call gpu_malloc_all(neighbors_list_d,size_maxnp_bytes, gpu_stream)
    call cpy_htod(c_loc(neighbors_list), neighbors_list_d, size_maxnp_bytes, gpu_stream )

  end subroutine gap_backend_begin

!**************************************************************************
  subroutine gap_backend_end()
    implicit none

    call gpu_free_async(species_d,gpu_stream)
    call gpu_free_async(n_neigh_d,gpu_stream)
    call gpu_free_async(neighbor_species_d,gpu_stream)
    call gpu_free_async(rjs_d,gpu_stream)
    call gpu_free_async(xyz_d,gpu_stream)
    call gpu_free(neighbors_list_d)        

  end subroutine gap_backend_end

!**************************************************************************
! Accumulate the two-body energies, forces and virial on the GPU.
  subroutine add_2b_contribution( n_distance_2b, distance_2b_hypers, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       energies_2b, forces_2b, virial_2b, time_2b )

    implicit none
 integer, intent(in) :: n_distance_2b
 type(distance_2b), allocatable, intent(inout), target :: distance_2b_hypers(:)
 type(input_parameters), intent(inout) :: params
 real(dp), intent(in), allocatable, target :: rjs(:)
 real(dp), intent(in), allocatable, target :: xyz(:,:)
 integer, intent(in), allocatable, target :: n_neigh(:)
 integer, intent(in), allocatable, target :: species(:)
 integer, intent(in), allocatable, target :: neighbor_species(:)
 integer, intent(in) :: i_beg
 integer, intent(in) :: i_end
 integer, intent(in) :: j_beg
 integer, intent(in) :: j_end
 real(dp), intent(inout), allocatable, target :: this_energies(:)
 real(dp), intent(inout), allocatable, target :: this_forces(:,:)
 real(dp), intent(inout), target :: this_virial(1:3,1:3)
 real(dp), intent(inout), allocatable, target :: energies_2b(:)
 real(dp), intent(inout), allocatable, target :: forces_2b(:,:)
 real(dp), intent(inout), target :: virial_2b(1:3,1:3)
 real(dp), intent(inout) :: time_2b(1:3)
 type(c_ptr) :: energies_2b_d
 type(c_ptr) :: forces_2b_d
 type(c_ptr) :: virial_2b_d
 integer :: i
 integer :: j
 integer :: sp1
 integer :: sp2
 integer :: n_sparse
 integer(c_size_t) :: st_n_sites_double
 integer(c_size_t) :: st_n_sparse_double
 integer(c_size_t) :: st_virial
 logical(c_bool) :: c_do_forces
 integer :: n_sites
 real(dp) :: t1
 real(dp) :: t2

    n_sites = i_end - i_beg + 1
    c_do_forces = logical( params%do_forces, kind=c_bool )
!   These two were set in the 2b procedure and read by the other two through
!   host association -- an undocumented ordering dependency between the three,
!   which is exactly what a shared host scope hides. Each computes its own now.
    st_n_sites_double = n_sites * c_double
    st_virial = 9 * c_double

    if( n_distance_2b == 0 )return

    st_n_sites_double=n_sites*sizeof(energies_2b(1))
    call gpu_malloc_all(energies_2b_d,st_n_sites_double,gpu_stream)
    call cpy_htod(c_loc(energies_2b),energies_2b_d, st_n_sites_double,gpu_stream)
    call gpu_malloc_all(forces_2b_d,3*st_n_sites_double,gpu_stream)
    call cpy_htod(c_loc(forces_2b),forces_2b_d, 3*st_n_sites_double,gpu_stream)
    st_virial=9*sizeof(virial_2b(1,1))
    call gpu_malloc_all(virial_2b_d,st_virial,gpu_stream)
    call cpy_htod(c_loc(virial_2b),virial_2b_d, st_virial,gpu_stream)



    do i = 1, n_distance_2b
    !! time_2b(1)=MPI_Wtime()
    !call get_time( ! time_2b(1) )

    call get_time( time_2b(1) )


    !                The kernel filters neighbours by species index, so sp1/sp2
    !                must be the position of the descriptor's species within
    !                params%species_types (j), not the descriptor index (i).
    do j = 1, size(params%species_types)
    if( distance_2b_hypers(i)%species1 == params%species_types(j) )then
    sp1 = j
    exit
    end if
    end do
    do j = 1, size(params%species_types)
    if( distance_2b_hypers(i)%species2 == params%species_types(j) )then
    sp2 = j
    exit
    end if
    end do



    n_sparse = distance_2b_hypers(i)%n_sparse
    st_n_sparse_double=n_sparse*sizeof( distance_2b_hypers(i)%alphas(1))
    call gpu_malloc_all(alphas_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(distance_2b_hypers(i)%alphas),alphas_d,st_n_sparse_double,gpu_stream)
    call gpu_malloc_all(cutoff_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(distance_2b_hypers(i)%cutoff),cutoff_d,st_n_sparse_double,gpu_stream)
    call gpu_malloc_all(qs_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(distance_2b_hypers(i)%Qs(:,1)),qs_d,st_n_sparse_double,gpu_stream)

    call get_time( t1 )

    call gpu_get_2b_forces_energies(i_beg, i_end,&
    & n_sparse, energies_2b_d, 0.0d0, n_neigh_d,&
    & c_do_forces, forces_2b_d, virial_2b_d,&
    & rjs_d, distance_2b_hypers(i)%rcut,&
    & species_d, neighbor_species_d, sp1, sp2,&
    & 0.5d0, distance_2b_hypers(i)%delta,&
    & cutoff_d, qs_d, distance_2b_hypers(i)&
    &%sigma, alphas_d, xyz_d, gpu_stream)

    !time_2b(2) = MPI_wtime()
    call get_time( time_2b(2)  )

    !          print *, rank, " >>--- Finished 2b energies forces on gpu ---"
    call gpu_free_async(alphas_d,gpu_stream)
    call gpu_free_async(cutoff_d,gpu_stream)
    call gpu_free_async(qs_d,gpu_stream)

    call get_time( time_2b(2) )

    time_2b(3) = time_2b(3) + time_2b(2) - time_2b(1)
    end do

    !              The kernels accumulate into the device buffers, so after the
    !              loop these already hold the sum over every 2b descriptor.
    !              Reading them back inside the loop and adding to the host arrays
    !              counted each descriptor once more for every later iteration.
    call cpy_dtoh(energies_2b_d, c_loc(this_energies),st_n_sites_double, gpu_stream)
    call cpy_dtoh(forces_2b_d, c_loc(this_forces),3*st_n_sites_double, gpu_stream)
    call cpy_dtoh(virial_2b_d, c_loc(this_virial),st_virial, gpu_stream)

    call gpu_stream_sync( gpu_stream )
    energies_2b = energies_2b + this_energies
    if( params%do_forces )then
    forces_2b = forces_2b + this_forces
    virial_2b = virial_2b + this_virial
    end if

    call gpu_free_async(energies_2b_d,gpu_stream)
    call gpu_free_async(forces_2b_d,gpu_stream)
    call gpu_free_async(virial_2b_d,gpu_stream)
    !        print *, rank, " >>~~~ Finished freeing 2b energies forces on gpu ~~~<<"


  end subroutine add_2b_contribution

!**************************************************************************
! Accumulate the core-potential energies, forces and virial on the GPU.
  subroutine add_core_pot_contribution( n_core_pot, core_pot_hypers, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       energies_core_pot, forces_core_pot, virial_core_pot, time_core_pot )

    implicit none
 integer, intent(in) :: n_core_pot
 type(core_pot), allocatable, intent(inout), target :: core_pot_hypers(:)
 type(input_parameters), intent(inout) :: params
 real(dp), intent(in), allocatable, target :: rjs(:)
 real(dp), intent(in), allocatable, target :: xyz(:,:)
 integer, intent(in), allocatable, target :: n_neigh(:)
 integer, intent(in), allocatable, target :: species(:)
 integer, intent(in), allocatable, target :: neighbor_species(:)
 integer, intent(in) :: i_beg
 integer, intent(in) :: i_end
 integer, intent(in) :: j_beg
 integer, intent(in) :: j_end
 real(dp), intent(inout), allocatable, target :: this_energies(:)
 real(dp), intent(inout), allocatable, target :: this_forces(:,:)
 real(dp), intent(inout), target :: this_virial(1:3,1:3)
 real(dp), intent(inout), allocatable, target :: energies_core_pot(:)
 real(dp), intent(inout), allocatable, target :: forces_core_pot(:,:)
 real(dp), intent(inout), target :: virial_core_pot(1:3,1:3)
 real(dp), intent(inout) :: time_core_pot(1:3)
 type(c_ptr) :: energies_core_pot_d
 type(c_ptr) :: forces_core_pot_d
 type(c_ptr) :: virial_core_pot_d
 type(c_ptr) :: x_d
 type(c_ptr) :: V_d
 type(c_ptr) :: dVdx2_d
 integer :: i
 integer :: j
 integer :: sp1
 integer :: sp2
 integer :: n_sparse
 integer(c_size_t) :: st_n_sites_double
 integer(c_size_t) :: st_n_sparse_double
 integer(c_size_t) :: st_virial
 logical(c_bool) :: c_do_forces
 integer :: n_sites
 real(dp) :: t1
 real(dp) :: t2

    n_sites = i_end - i_beg + 1
    c_do_forces = logical( params%do_forces, kind=c_bool )
!   These two were set in the 2b procedure and read by the other two through
!   host association -- an undocumented ordering dependency between the three,
!   which is exactly what a shared host scope hides. Each computes its own now.
    st_n_sites_double = n_sites * c_double
    st_virial = 9 * c_double

    if( n_core_pot == 0 )return


    !        print *, rank, " > Allocating core_pot on gpu "
    st_n_sites_double=n_sites*sizeof(energies_core_pot(1))
    call gpu_malloc_all(energies_core_pot_d,st_n_sites_double,gpu_stream)
    call cpy_htod(c_loc(energies_core_pot),energies_core_pot_d, st_n_sites_double,gpu_stream)
    call gpu_malloc_all(forces_core_pot_d,3*st_n_sites_double,gpu_stream)
    call cpy_htod(c_loc(forces_core_pot),forces_core_pot_d, 3*st_n_sites_double,gpu_stream)
    st_virial=9*sizeof(virial_core_pot(1,1))
    call gpu_malloc_all(virial_core_pot_d,st_virial,gpu_stream)
    call cpy_htod(c_loc(virial_core_pot),virial_core_pot_d, st_virial,gpu_stream)

    !       Loop through core_pot descriptors
    do i = 1, n_core_pot

    !           print *, " > Getting core potential"
    call get_time( time_core_pot(1) )

    n_sparse = core_pot_hypers(i)%n
    st_n_sparse_double=n_sparse*sizeof( core_pot_hypers(i)%x(1))
    call gpu_malloc_all(x_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(core_pot_hypers(i)%x),x_d,st_n_sparse_double,gpu_stream)
    call gpu_malloc_all(V_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(core_pot_hypers(i)%V),V_d,st_n_sparse_double,gpu_stream)
    call gpu_malloc_all(dVdx2_d,st_n_sparse_double,gpu_stream)
    call cpy_htod(c_loc(core_pot_hypers(i)%dVdx2),dVdx2_d,st_n_sparse_double,gpu_stream)


    ! print *, rank, " >>--- Finished allocating core_pot on gpu ---"
    ! print *, rank, " > Starting core_pot energies forces on gpu "                            

!   The kernel filters neighbours by species index, so sp1/sp2 must be the
!   position of this descriptor's species within params%species_types. They
!   were never set here: the routine passed whatever add_2b_contribution had
!   left in the enclosing scope, which is the same defect as bug 3 in section 4
!   and was fixed there but not here.
    do j = 1, size(params%species_types)
      if( core_pot_hypers(i)%species1 == params%species_types(j) )then
        sp1 = j
        exit
      end if
    end do
    do j = 1, size(params%species_types)
      if( core_pot_hypers(i)%species2 == params%species_types(j) )then
        sp2 = j
        exit
      end if
    end do

    call gpu_get_core_pot_energy_and_forces(i_beg, i_end, c_do_forces, species_d, sp1, sp2, n_neigh_d, neighbor_species_d,&
    rjs_d, n_sparse, x_d, V_d, dVdx2_d, core_pot_hypers(i)%yp1, core_pot_hypers(i)%ypn,&
    xyz_d, forces_core_pot_d, virial_core_pot_d, energies_core_pot_d, gpu_stream)
    !          print *, rank, " >>--- Finished core_pot energies forces on gpu ---<<"                                      
    call gpu_free_async(x_d,gpu_stream)
    call gpu_free_async(V_d,gpu_stream)
    call gpu_free_async(dVdx2_d,gpu_stream)

    call get_time( time_core_pot(2) )

    time_core_pot(3) = time_core_pot(3) + time_core_pot(2) - time_core_pot(1)
    end do

    !              As for the 2b and 3b terms, the kernel accumulates into the
    !              device buffers, so the totals are read back only once the loop
    !              over descriptors has finished.
    call cpy_dtoh(energies_core_pot_d, c_loc(this_energies),st_n_sites_double, gpu_stream)
    call cpy_dtoh(forces_core_pot_d, c_loc(this_forces),3*st_n_sites_double, gpu_stream)
    call cpy_dtoh(virial_core_pot_d, c_loc(this_virial),st_virial, gpu_stream)

    call gpu_stream_sync( gpu_stream )

    energies_core_pot = energies_core_pot + this_energies
    if( params%do_forces )then
    forces_core_pot = forces_core_pot + this_forces
    virial_core_pot = virial_core_pot + this_virial
    end if

    call gpu_free_async(energies_core_pot_d,gpu_stream)
    call gpu_free_async(forces_core_pot_d,gpu_stream)
    call gpu_free_async(virial_core_pot_d,gpu_stream)


  end subroutine add_core_pot_contribution

!**************************************************************************
! Accumulate the three-body energies, forces and virial on the GPU.
  subroutine add_3b_contribution( n_angle_3b, angle_3b_hypers, neighbors_list, &
       params, rjs, xyz, n_neigh, species, neighbor_species, &
       i_beg, i_end, j_beg, j_end, this_energies, this_forces, this_virial, &
       forces, energies_3b, forces_3b, virial_3b, time_3b )

    implicit none
 integer, intent(in) :: n_angle_3b
 type(angle_3b), allocatable, intent(inout), target :: angle_3b_hypers(:)
 integer, intent(in), allocatable, target :: neighbors_list(:)
 type(input_parameters), intent(inout) :: params
 real(dp), intent(in), allocatable, target :: rjs(:)
 real(dp), intent(in), allocatable, target :: xyz(:,:)
 integer, intent(in), allocatable, target :: n_neigh(:)
 integer, intent(in), allocatable, target :: species(:)
 integer, intent(in), allocatable, target :: neighbor_species(:)
 integer, intent(in) :: i_beg
 integer, intent(in) :: i_end
 integer, intent(in) :: j_beg
 integer, intent(in) :: j_end
 real(dp), intent(inout), allocatable, target :: this_energies(:)
 real(dp), intent(inout), allocatable, target :: this_forces(:,:)
 real(dp), intent(inout), target :: this_virial(1:3,1:3)
!   Read only for its extent; the 3b kernel sizes a device buffer from it.
 real(dp), intent(in), allocatable :: forces(:,:)
 real(dp), intent(inout), allocatable, target :: energies_3b(:)
 real(dp), intent(inout), allocatable, target :: forces_3b(:,:)
 real(dp), intent(inout), target :: virial_3b(1:3,1:3)
 real(dp), intent(inout) :: time_3b(1:3)
 type(c_ptr) :: energies_3b_d
 type(c_ptr) :: forces_3b_d
 type(c_ptr) :: virial_3b_d
 integer :: i
 integer :: j
 integer :: sp1
 integer :: sp2
 integer :: n_sparse
 integer(c_size_t) :: st_n_sites_double
 integer(c_size_t) :: st_n_sparse_double
 integer(c_size_t) :: st_virial
 logical(c_bool) :: c_do_forces
 integer :: n_sites
 real(dp) :: t1
 real(dp) :: t2
 integer :: k
 integer :: max_np
 integer :: sp0_3b
 integer :: sp1_3b
 integer :: sp2_3b
 integer, allocatable, target :: kappas(:)
 type(c_ptr) :: kappas_array_d
    character(kind=c_char,len=4) :: c_name_3b
 integer(c_size_t) :: size_maxnp_bytes
 integer(c_size_t) :: size_maxnp_qs_bytes
 integer(c_size_t) :: size_alphas_bytes
 integer(c_size_t) :: size_energy3b
 integer(c_size_t) :: size_forces3b
 integer(c_size_t) :: size_virial3b
 type(c_ptr) :: sigma_d

    n_sites = i_end - i_beg + 1
    c_do_forces = logical( params%do_forces, kind=c_bool )
!   These two were set in the 2b procedure and read by the other two through
!   host association -- an undocumented ordering dependency between the three,
!   which is exactly what a shared host scope hides. Each computes its own now.
    st_n_sites_double = n_sites * c_double
    st_virial = 9 * c_double

    if( n_angle_3b == 0 )return

    size_energy3b = size(n_neigh)*c_double
    call gpu_malloc_all(energies_3b_d,size_energy3b, gpu_stream)
    call gpu_memset_async (energies_3b_d, 0, size_energy3b,gpu_stream)
    size_forces3b = size(forces,2) * 3 * c_double
    call gpu_malloc_all(forces_3b_d,size_forces3b, gpu_stream)
    call gpu_memset_async (forces_3b_d, 0, size_forces3b,gpu_stream)
    size_virial3b = 9 * c_double
    call gpu_malloc_all(virial_3b_d,size_virial3b, gpu_stream)
    call gpu_memset_async (virial_3b_d, 0, size_virial3b,gpu_stream)

    size_maxnp_bytes = size(n_neigh)* c_int
    call gpu_malloc_all(kappas_array_d,size_maxnp_bytes,gpu_stream)



    allocate(kappas(1:n_sites))


    k = 0
    do i = i_beg, i_end 
    kappas( i ) = k            
    do j = 1, n_neigh(i)
    k = k + 1
    end do
    end do


    ! do i = 1, size(n_neigh)
    !    if( i == 1 )then
    !       kappas( i ) = 0
    !    else
    !       kappas( i ) = n_neigh( i-1 ) + kappas( i-1 )
    !    end if
    ! end do

    call cpy_htod(c_loc(kappas), kappas_array_d, size_maxnp_bytes, gpu_stream )
    call gpu_stream_sync( gpu_stream )
    deallocate(kappas)

    !        call gpu_create_kappas(kappas_array_d, c_loc(n_neigh),gpu_stream, size(n_neigh))



    max_np = 0
    do i = 1, n_angle_3b
    if (angle_3b_hypers(i)%n_sparse > max_np)then
    max_np = angle_3b_hypers(i)%n_sparse
    end if
    end do


    !write(0,*) "max np is: ",max_np
    size_maxnp_bytes = max_np* c_double
    call gpu_malloc_all(cutoff_d,size_maxnp_bytes,gpu_stream)
    call gpu_malloc_all(alphas_d,size_maxnp_bytes,gpu_stream)
    size_maxnp_qs_bytes = 3*size_maxnp_bytes
    call gpu_malloc_all(qs_d,size_maxnp_qs_bytes,gpu_stream)
    size_alphas_bytes = 3 * c_double
    call gpu_malloc_all(sigma_d,size_alphas_bytes,gpu_stream)


    !       Loop through angle_3b descriptors
    do i = 1, n_angle_3b
    call get_time( time_3b(1) )

    call cpy_htod(c_loc(angle_3b_hypers(i)%cutoff),cutoff_d,size_maxnp_bytes,gpu_stream)
    call cpy_htod(c_loc(angle_3b_hypers(i)%sigma),sigma_d,size_alphas_bytes,gpu_stream)
    call cpy_htod(c_loc(angle_3b_hypers(i)%qs),qs_d,size_maxnp_qs_bytes,gpu_stream)
    call cpy_htod(c_loc(angle_3b_hypers(i)%alphas),alphas_d,size_maxnp_bytes,gpu_stream)

    ! print *, rank, " >> Finished allocation 3b on gpu "
    ! print *, rank, " >> Starting setup 3b on gpu "                                                
    call setup_3b_gpu(angle_3b_hypers(i)%kernel_type,angle_3b_hypers(i)%species_center,&
    angle_3b_hypers(i)%species1, angle_3b_hypers(i)%species2, params%species_types,&
    c_name_3b, sp0_3b, sp1_3b, sp2_3b)

    call gpu_3b(size(angle_3b_hypers(i)%alphas),  i_end-i_beg+1, size(rjs), size(forces,2), &
    sp0_3b, sp1_3b, sp2_3b, alphas_d, angle_3b_hypers(i)%delta, 0.d0, cutoff_d, gpu_stream,&
    rjs_d, xyz_d, n_neigh_d,species_d,neighbors_list_d,neighbor_species_d, &
    c_do_forces, angle_3b_hypers(i)%rcut, 0.5d0, sigma_d, qs_d, c_name_3b,&
    i_beg, i_end, energies_3b_d, forces_3b_d, virial_3b_d, kappas_array_d)

    call get_time( time_3b(2) )

    time_3b(3) = time_3b(3) + time_3b(2) - time_3b(1)

    ! print *, rank, " >>--- Finished 3b on gpu ---<< took ",&
    !   time_3b(2) - time_3b(1), " with total being ", time_3b(3)

    end do

    !              As for the 2b term, the kernels accumulate into the device
    !              buffers, so the totals are only read back once the loop over
    !              descriptors has finished.
    call cpy_dtoh(energies_3b_d, c_loc(this_energies),st_n_sites_double, gpu_stream)
    call cpy_dtoh(forces_3b_d, c_loc(this_forces), 3*st_n_sites_double, gpu_stream)
    call cpy_dtoh(virial_3b_d, c_loc(this_virial),st_virial, gpu_stream)

    call gpu_stream_sync( gpu_stream )

    energies_3b = energies_3b + this_energies
    if( params%do_forces )then
    forces_3b = forces_3b + this_forces
    virial_3b = virial_3b + this_virial
    end if

    !              All of the buffers below are allocated once before the descriptor
    !              loop (the parameter buffers are sized for max_np so they can be
    !              reused by every descriptor). They must therefore only be released
    !              after the loop has finished -- freeing them per iteration left
    !              dangling device pointers that the next iteration copied into.
    call gpu_free_async(energies_3b_d,gpu_stream)
    call gpu_free_async(forces_3b_d,gpu_stream)
    call gpu_free_async(kappas_array_d, gpu_stream)
    call gpu_free_async(virial_3b_d, gpu_stream)

    call gpu_free_async(cutoff_d,gpu_stream)
    call gpu_free_async(sigma_d,gpu_stream)
    call gpu_free_async(qs_d,gpu_stream)
    call gpu_free_async(alphas_d,gpu_stream)


  end subroutine add_3b_contribution

end module gap_backend
