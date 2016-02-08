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

# set the PREFIX_PATH to look for the CMake Trillinos libraries


# For Moonlight
#set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}; /usr/projects/draco/vendors/trilinos-12.0.1/moonlight-openmpi-1.6.5-intel-15/lib/cmake/Trilinos/)

# For CCS machines
#${CMAKE_PREFIX_PATH};/var/tmp/vendors/trilinos-11.10.2/openmpi-1.6.3-gcc-4.8.1/lib/cmake/Trilinos/)
# For Darwin
#set(${CMAKE_PREFIX_PATH}  ${CMAKE_PREFIX_PATH};/projects/groups/chicoma/packages/zoltan-3.8_intel-15.0.1_impi-5.0.2)
#set(Zoltan_base /projects/groups/chicoma/packages/zoltan-3.8_intel-15.0.1_impi-5.0.2)

###############################################################################
# Find Parmetis and link libraries 
###############################################################################
# CCS MACHINES
#set(METIS_INCLUDE_DIR  /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/include)
#set(PARMETIS_INCLUDE_DIR  /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/include)
#set(PARMETIS_LIBRARY_DIR /var/tmp/vendors/ParMetis-4.0.3/openmpi-1.10.0-gcc-5.2.0/lib )

# MOONLIGHT (DRACO) 
set(PARMETIS_INCLUDE_DIRS /usr/include/)
set(METIS_INCLUDE_DIR /usr/projects/draco/vendors/ParMetis-4.0.3/moonlight-openmpi-1.6.5-intel-15.0.3/include )
set(PARMETIS_LIBRARY_DIR /usr/projects/draco/vendors/ParMetis-4.0.3/moonlight-openmpi-1.6.5-intel-15.0.3/lib)

# MOONLIGHT (HOME DIRECTORY PARMETIS BUILT WITH INTEL MPI)
#set(METIS_INCLUDE_DIR /users/along/parmetis/install/include)
#set(PARMETIS_LIBRARY_DIR /users/along/parmetis/install/lib)

include_directories(${METIS_INCLUDE_DIR})
link_directories(${PARMETIS_LIBRARY_DIR})

#include_directories(${Zoltan_base}/include)
#link_directories(${Zoltan_base}/lib)
#message ("Zoltan include path: " ${Zoltan_INCLUDE_DIRS})
#message ("Zoltan lib path: " ${Zoltan_base}/lib)


###############################################################################
# Add Boost library
###############################################################################
#set(BOOST_ROOT /projects/opt/boost/1.58.0/draco/vendors/boost-1.58.0)
# MOONLIGHT 
set(BOOST_ROOT /usr/projects/draco/vendors/boost-1.59.0/moonlight-openmpi-1.6.5-intel-15.0.3)

# MOONLIGHT (HOME DIRECTORY COMPILED WITH INTEL MPI)
#set(BOOST_ROOT /users/along/boost/boost_install)

# CCS
#set(BOOST_ROOT /ccs/codes/radtran/vendors/boost-1.58.0/openmpi-1.6.3-gcc-4.8.1/)

#set(BOOST_ROOT /usr/lib/i386-linux-gnu)
#set(BOOST_ROOT $ENV{BOOST_ROOT})
#message ("Boost Root: " ${BOOST_ROOT})
#set(Boost_NO_SYSTEM_PATH ON)
#find_package(Boost 1.58.0 COMPONENTS mpi REQUIRED)
set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
set(Boost_LIBRARY_DIR ${BOOST_ROOT}/lib)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
#message ("Boost Flags: " ${Boost_INCLUDE_DIRS} " " ${Boost_VERSION} )
message ("Boost Libraries: " ${Boost_LIBRARIES})

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


add_executable(BRANSON main.cc)
target_link_libraries(BRANSON parmetis)
target_link_libraries(BRANSON metis)
target_link_libraries(BRANSON boost_mpi-mt)
target_link_libraries(BRANSON boost_serialization-mt)
target_link_libraries(BRANSON boost_wserialization-mt)

enable_testing()

add_subdirectory(test)
get_directory_property(Boost_INCLUDE_DIR DIRECTORY test VARIABLES Boost_INCLUDE_DIR)
get_directory_property(Boost_LIB_DIR DIRECTORY test VARIABLES Boost_LIB_DIR)

#install(TARGETS BRANSON DESTINATION bin)
