program main
  !! Test program which splits the mpi_comm_world into two groups, and launches
  !! two instances of turbogap_predict at once, one with each mpi group.
  !! The input filenames are different, each instance is computing a different
  !! atomic structure.
  !!
  use mpi
  use turbogap_wrap, only: turbogap_predict
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  integer :: ierr
  integer :: world_me, world_np
  integer :: group_me, group_np
  integer :: mygroup, group_comm
  integer :: turbo_err
  real( c_double ) :: energy
  real( c_double ), allocatable :: force(:,:)

  call mpi_init(ierr)
  ! get world info
  call mpi_comm_size( mpi_comm_world, world_np, ierr )
  call mpi_comm_rank( mpi_comm_world, world_me, ierr )

  ! split the world into 2
  if( world_me==0) then
     mygroup=0
  else
     mygroup=1
  end if
  call mpi_comm_split( mpi_comm_world, mygroup, world_me, group_comm, ierr )
  ! get group info
  call mpi_comm_size( group_comm, group_np, ierr )
  call mpi_comm_rank( group_comm, group_me, ierr )


  if( mygroup==0 ) then

     ! pass input file input0
     call turbogap_predict( group_comm, "input0", turbo_err, output_energy=energy, output_forces=force )

  elseif( mygroup==1 )then

     ! pass input file input1
     call turbogap_predict( group_comm, "input1", turbo_err, output_energy=energy, output_forces=force )

  end if


  if( group_me == 0) then
     write(*,*) "group", mygroup," got turbo_err:", turbo_err
     write(*,*) "group",mygroup, "etot:",etot
  end if


  call mpi_finalize( ierr )
end program main

