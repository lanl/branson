# Unit tests
cmake_minimum_required (VERSION 3.11)

# force mpiexec to be srun in cray environments
message( STATUS "Looking to see if we are building in a Cray Environment...")
if( DEFINED ENV{CRAYPE_VERSION} )
  set(MPIEXEC_EXECUTABLE srun)
  message("In Cray Environment, MPI executable is: ${MPIEXEC_EXECUTABLE}")
endif()

# allow run as root
if( "$ENV{GITLAB_CI}" STREQUAL "true" OR "$ENV{TRAVIS}" STREQUAL "true")
  message("In CI mode, allow running mpiexec as root")
  set(runasroot "--allow-run-as-root --oversubscribe")
endif()

# set extra flags for mpiexec in testing
set( MPIEXEC_POSTFLAGS "-bind-to none ${runasroot}" CACHE STRING
  "extra mpirun flags (list)." FORCE)

string( REPLACE " " ";" mpiexec_postflags_list "${MPIEXEC_POSTFLAGS}" )

#------------------------------------------------------------------------------#
# Helper for building and registering tests.
#
# add_branson_test(
#   SOURCE  "foo.cc"
#   PE_LIST "1;4" )

function( add_branson_test )
  cmake_parse_arguments(
    abt
    ""
    "SOURCE"
    "PE_LIST"
    ${ARGV} )
  string( REPLACE ".cc" "_exe" build_target ${abt_SOURCE} )
  string( REPLACE ".cc" ""     ctest_name   ${abt_SOURCE} )

  # Build the test
  add_executable( ${build_target} ${abt_SOURCE} )
  target_include_directories( ${build_target} PRIVATE
    $<BUILD_INTERFACE:${BRANSON_BINARY_DIR}> )
  target_link_libraries( ${build_target} PUBLIC ${branson_deps} )

  # Register the test command with ctest.
  foreach( np ${abt_PE_LIST} )
    add_test(
      NAME    ${ctest_name}_${np}pe
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${mpiexec_postflags_list}
              $<TARGET_FILE:${build_target}> )
    set_tests_properties( ${ctest_name}_${np}pe PROPERTIES PROCESSORS ${np} )
  endforeach()
endfunction()

#------------------------------------------------------------------------------#
# 1 PE tests

set(branson_test_sources_1pe
  test_cell.cc
  test_proto_cell.cc
  test_counter_rng.cc
  test_input.cc
  test_mesh.cc
  test_mpi_types.cc
  test_photon.cc
  test_work_packet.cc
 )
foreach( btest ${branson_test_sources_1pe} )
  add_branson_test( SOURCE  ${btest} PE_LIST "1" )
endforeach()

#------------------------------------------------------------------------------#
# Other tests

add_branson_test( SOURCE test_imc_state.cc     PE_LIST "2" )

#------------------------------------------------------------------------------#
# copy these input files for Input, IMC_State, Mesh and write_silo tests

set( inputfiles
  simple_input.xml
  large_particle_input.xml
  three_region_mesh_input.xml )
foreach( ifile ${inputfiles} )
  configure_file(${ifile} ${CMAKE_CURRENT_BINARY_DIR}/${ifile} COPYONLY)
endforeach()

#------------------------------------------------------------------------------#
# write_silo function test
if (VIZ_LIBRARIES_FOUND)
   add_branson_test( SOURCE test_write_silo.cc "1")
endif()

#------------------------------------------------------------------------------#
# End test/CMakeLists.txt
#------------------------------------------------------------------------------#
