# Unit tests
cmake_minimum_required (VERSION 3.11)

# if (${SITENAME} STREQUAL "Trinitite")
#   set(MPI_EXEC srun)
# endif ()

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
  target_link_libraries( ${build_target} PUBLIC
    ParMETIS::parmetis;MPI::MPI_CXX;MPI::MPI_C;METIS::metis )

  # Register the test command with ctest.
  foreach( np ${abt_PE_LIST} )
    add_test(
      NAME    ${ctest_name}_${np}pe
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              $<TARGET_FILE:${build_target}> )
    set_tests_properties( ${ctest_name}_${np}pe PROPERTIES PROCESSORS ${np} )
  endforeach()
endfunction()

#------------------------------------------------------------------------------#
# 1 PE tests

set(branson_test_sources_1pe
  test_buffer.cc
  test_cell.cc
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

add_branson_test( SOURCE test_remap_census.cc  PE_LIST "8" )
add_branson_test( SOURCE test_load_balance.cc  PE_LIST "4" )
add_branson_test( SOURCE test_tally_manager.cc PE_LIST "8" )
add_branson_test( SOURCE test_parmetis.cc      PE_LIST "1;4" )
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

# # write_silo function test
# if (VIZ_LIBRARIES_FOUND)
#   add_executable(test_write_silo_exe test_write_silo.cc)
#   target_link_libraries(test_write_silo_exe parmetis)
#   target_link_libraries(test_write_silo_exe metis)
#   target_link_libraries(test_write_silo_exe hdf5)
#   target_link_libraries(test_write_silo_exe siloh5)
#   add_test(NAME test_write_silo COMMAND ${MPI_EXEC} -n 1 ./test_write_silo_exe)
# endif ()

#------------------------------------------------------------------------------#
# End test/CMakeLists.txt
#------------------------------------------------------------------------------#
