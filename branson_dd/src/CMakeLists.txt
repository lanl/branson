cmake_minimum_required (VERSION 2.8)
# Set Intel MPI Compiler
#set(CMAKE_CXX_COMPILER mpiicpc)
# Set OpenMPI compiler
set(CMAKE_CXX_COMPILER mpic++)
project (BRANSON)

#set compiler
message ("Compiler: " ${CMAKE_CXX_COMPILER} )

# add compiler flags
set( DEBUG_FLAGS -gdwarf-3)
set( DEBUG_FLAGS "${DEBUG_FLAGS} -Wall -std=c++0x")
set( DEBUG_FLAGS "${DEBUG_FLAGS} -O3")
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS}" )
message ("Compiler Flags: " ${CMAKE_CXX_FLAGS} )

#set version number
set(BRANSON_VERSION_MAJOR 0)
set(BRANSON_VERSION_MINOR 1)

###############################################################################
# Find Parmetis and link libraries
# (user_contrib module sets environment variables pointing to correct ParMetis)
###############################################################################
# CCS MACHINES
#set(METIS_INCLUDE_DIR  /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/include)
#set(PARMETIS_INCLUDE_DIR  /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/include)
#set(PARMETIS_LIBRARY_DIR /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/lib )

# MOONLIGHT (DRACO) 
#set(PARMETIS_INCLUDE_DIRS /usr/include/)
#set(METIS_INCLUDE_DIR /usr/projects/draco/vendors/ParMetis-4.0.3/moonlight-openmpi-1.6.5-intel-15.0.3/include )
#set(PARMETIS_LIBRARY_DIR /usr/projects/draco/vendors/ParMetis-4.0.3/moonlight-openmpi-1.6.5-intel-15.0.3/lib)

message("ParMetis from: " $ENV{ParMETIS_ROOT_DIR})
include_directories( $ENV{METIS_INC_DIR})
link_directories($ENV{PARMETIS_LIB_DIR})

###############################################################################
# Add Boost library
###############################################################################
#set(BOOST_ROOT /projects/opt/boost/1.58.0/draco/vendors/boost-1.58.0)

# MOONLIGHT 
#set(BOOST_ROOT /usr/projects/draco/vendors/boost-1.59.0/moonlight-openmpi-1.6.5-intel-15.0.3)
#set(BOOST_ROOT /usr/projects/draco/vendors/boost-1.60.0/moonlight-openmpi-1.10.1-intel-15.0.3)

# CCS
#set(BOOST_ROOT /ccs/codes/radtran/vendors/boost-1.58.0/openmpi-1.6.3-gcc-4.8.1/)

message("Boost from: " $ENV{BOOST_ROOT})
set(Boost_INCLUDE_DIR $ENV{BOOST_ROOT}/include)
set(Boost_LIBRARY_DIR $ENV{BOOST_ROOT}/lib)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

###############################################################################
# Configure executable
###############################################################################
add_executable(BRANSON main.cc)
target_link_libraries(BRANSON parmetis)
target_link_libraries(BRANSON metis)
target_link_libraries(BRANSON boost_mpi-mt)
target_link_libraries(BRANSON boost_serialization-mt)
target_link_libraries(BRANSON boost_wserialization-mt)

enable_testing()

add_subdirectory(test)
# Define variables needed to compile and link unit test executables
get_directory_property(Boost_INCLUDE_DIR DIRECTORY test VARIABLES Boost_INCLUDE_DIR)
get_directory_property(Boost_LIB_DIR DIRECTORY test VARIABLES Boost_LIB_DIR)

#install(TARGETS BRANSON DESTINATION bin)
