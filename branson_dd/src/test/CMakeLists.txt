# Unit tests
cmake_minimum_required (VERSION 2.8)
# Set Intel MPI Compiler
#set(CMAKE_CXX_COMPILER mpiicpc)
# Set OpenMPI compiler
set(CMAKE_CXX_COMPILER mpic++)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# Cell class test
add_executable(test_cell_exe test_cell.cc)
target_link_libraries(test_cell_exe boost_mpi-mt)
target_link_libraries(test_cell_exe boost_serialization-mt)
target_link_libraries(test_cell_exe boost_wserialization-mt)
add_test(NAME test_cell COMMAND test_cell_exe)

# Photon class test
add_executable(test_photon_exe test_photon.cc)
target_link_libraries(test_photon_exe boost_serialization-mt)
target_link_libraries(test_photon_exe boost_wserialization-mt)
add_test(NAME test_photon COMMAND test_photon_exe)

# Buffer class test
add_executable(test_buffer_exe test_buffer.cc)
add_test(NAME test_buffer COMMAND test_buffer_exe)
