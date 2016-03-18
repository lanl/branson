# Unit tests
cmake_minimum_required (VERSION 2.8)
# Set Intel MPI Compiler
#set(CMAKE_CXX_COMPILER mpiicpc)
# Set OpenMPI compiler
set(CMAKE_CXX_COMPILER mpic++)

include_directories(${Boost_INCLUDE_DIR})

# Cell class test
add_executable(test_cell_exe test_cell.cc)
add_test(NAME test_cell COMMAND test_cell_exe)

# Photon class test
add_executable(test_photon_exe test_photon.cc)
add_test(NAME test_photon COMMAND test_photon_exe)

# Buffer class test
add_executable(test_buffer_exe test_buffer.cc)
add_test(NAME test_buffer COMMAND test_buffer_exe)

# copy these input files for Input and IMC_State class tests
configure_file(simple_input.xml ${CMAKE_CURRENT_BINARY_DIR}/simple_input.xml COPY_ONLY)
configure_file(large_particle_input.xml ${CMAKE_CURRENT_BINARY_DIR}/large_particle_input.xml COPY_ONLY)

# Input class test
add_executable(test_input_exe test_input.cc)
add_test(NAME test_input COMMAND test_input_exe)

# IMC_State class test
add_executable(test_imc_state_exe test_imc_state.cc)
add_test(NAME test_imc_state COMMAND mpirun -np 2 test_imc_state_exe)
