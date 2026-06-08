cmake_minimum_required(VERSION 3.21)

project(turbogap
  DESCRIPTION ""
  VERSION 0.0.1
  LANGUAGES Fortran
  )

# source files
set( SRC
  src/splines.f90
  src/types.f90
  src/neighbors.f90
  src/gap.f90
  src/vdw.f90
  src/local_properties.f90
  src/exp_utils.f90
  src/xyz.f90
  src/md.f90
  src/mc.f90
  src/read_files.f90
  src/gap_interface.f90
  src/mpi.f90
  src/exp_interface.f90
  src/turbogap_routine.f90
  # SRC_TP_BT
  src/third_party/bussi_thermostat/resamplekin.f90
  # SRC_ST
  src/soap_turbo/src/soap_turbo_functions.f90
  src/soap_turbo/src/soap_turbo_radial.f90
  src/soap_turbo/src/soap_turbo_angular.f90
  src/soap_turbo/src/soap_turbo.f90
  src/soap_turbo/src/soap_turbo_compress.f90
  # SRC_STOP
  src/stopping/adaptive_time.f90
  src/stopping/electronic_stopping.f90
  src/stopping/eph_beta.f90
  src/stopping/eph_fdm.f90
  src/stopping/eph_electronic_stopping.f90
  )

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# define library from source files
add_library( libturbogap ${SRC} )
# library properties
set_target_properties( libturbogap
  PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
  OUTPUT_NAME "turbogap"
  )
# set public include directories for consumer
target_include_directories( libturbogap
  PUBLIC
  "${CMAKE_CURRENT_BINARY_DIR}/include"
  )
# compile flags
target_compile_options( libturbogap
  PRIVATE
  "-O3"
  )
# compile definitions
target_compile_definitions( libturbogap
  PRIVATE
  _MPIF90
  )
# link to blas,lapack, mpi
target_link_libraries( libturbogap
  PUBLIC
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  MPI::MPI_Fortran
  )
# alias name for libturbogap
add_library( turbogap ALIAS libturbogap )


# define program
add_executable( binturbogap src/turbogap.f90 )
set_target_properties( binturbogap
  PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
  EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/bin"
  OUTPUT_NAME "turbogap"
  )
# compile flags
target_compile_options( binturbogap
  PRIVATE
  "-O3"
  )
# compile definitions
target_compile_definitions( binturbogap
  PRIVATE
  _MPIF90
  )
# link to libturbogap
target_link_libraries( binturbogap PRIVATE libturbogap )


# prepare exports
install(TARGETS
  libturbogap
  EXPORT
  "${PROJECT_NAME}-targets"
  )
# export during build, for `include( <#> NAMES turbogap-Build )`
export(EXPORT
  "${PROJECT_NAME}-targets"
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}-BuildConfig.cmake"
  NAMESPACE "${PROJECT_NAME}::"
  )
