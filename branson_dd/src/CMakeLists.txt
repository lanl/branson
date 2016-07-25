###############################################################################
# NOTES ON COMPILATION
###############################################################################
# Load the modules for each libary and the environment variables will be set 
# correctly

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

# Intel C++ 11 flags
set( DEBUG_FLAGS "${DEBUG_FLAGS} -Wall -std=c++0x")

# Cray flag for C++11
#set (DEBUG_FLAGS "${DEBUG_FLAGS} -h std=c++11")

# append optimization flags
set( DEBUG_FLAGS "${DEBUG_FLAGS} -O3")

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS}" )
message ("Compiler Flags: " ${CMAKE_CXX_FLAGS} )

#set version number
set(BRANSON_VERSION_MAJOR 0)
set(BRANSON_VERSION_MINOR 8)


###############################################################################
# Random123
# Download and install Random123
###############################################################################
#ExternalProject_Add(random123 URL https://www.deshawresearch.com/downloads/download_random123.cgi/Random123-1.09.tar.gz)

###############################################################################
# metis and parmetis
# Load modules for metis and parmetis to get correct environment variables
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
include_directories( $ENV{PARMETIS_INC_DIR})
link_directories($ENV{PARMETIS_LIB_DIR})

message("Metis from: " $ENV{METIS_ROOT_DIR})
include_directories( $ENV{METIS_INC_DIR})
link_directories($ENV{METIS_LIB_DIR})

###############################################################################
# boost (headers only)
# Load boost module to get correct environement variables
###############################################################################
#set(BOOST_ROOT /projects/opt/boost/1.58.0/draco/vendors/boost-1.58.0)

# MOONLIGHT 
#set(BOOST_ROOT /usr/projects/draco/vendors/boost-1.59.0/moonlight-openmpi-1.6.5-intel-15.0.3)
#set(BOOST_ROOT /usr/projects/draco/vendors/boost-1.60.0/moonlight-openmpi-1.10.1-intel-15.0.3)

# CCS
#set(BOOST_ROOT /ccs/codes/radtran/vendors/boost-1.58.0/openmpi-1.6.3-gcc-4.8.1/)

message("Boost from: " $ENV{BOOST_ROOT})
set(Boost_INCLUDE_DIR $ENV{BOOST_ROOT}/include)
#set(Boost_LIBRARY_DIR $ENV{BOOST_ROOT}/lib)

include_directories(${Boost_INCLUDE_DIR})
#link_directories(${Boost_LIBRARY_DIR})


###############################################################################
# Silo and HDF5 libraries
# Load modules for hdf5 and solo to get correct environment variables
# use find package
###############################################################################
message("Silo from: $ENV{SILO_ROOT} " )
include_directories($ENV{SILO_ROOT}/include)
link_directories($ENV{SILO_ROOT}/lib)
if ("$ENV{SILO_ROOT}x" STREQUAL "x")
  set(SILO_FOUND FALSE)
else ()
  set(SILO_FOUND TRUE)
endif ()
message("SILO_FOUND = ${SILO_FOUND}")

find_package(HDF5)
message("HDF5 from: $ENV{HDF5_ROOT} " )
message("HDF5_FOUND = ${HDF5_FOUND}")
message("HDF from: $ENV{HDF5_ROOT} ")
link_directories($ENV{HDF5_ROOT}/lib)

if (HDF5_FOUND AND SILO_FOUND)
  set(VIZ_LIBRARIES_FOUND TRUE)
endif ()

###############################################################################
# Set up libraries for Cray RMA routines
###############################################################################
SET(DMAPP_DYNAMIC -Wl,--whole-archive,-ldmapp,--no-whole-archive)

###############################################################################
# Create config file
###############################################################################
configure_file(config.h.in ${PROJECT_BINARY_DIR}/config.h)
include_directories(${PROJECT_BINARY_DIR})

###############################################################################
# Configure executable and link libraries
###############################################################################
add_executable(BRANSON main.cc)
target_link_libraries(BRANSON parmetis)
target_link_libraries(BRANSON metis)
if (VIZ_LIBRARIES_FOUND)
  target_link_libraries(BRANSON hdf5)
  target_link_libraries(BRANSON siloh5)
endif ()

# these lines link the Cray DMAPP library
#first, create the dynamic version of dmapp library
#target_link_libraries(BRANSON ${DMAPP_DYNAMIC})
#then link to it dynamically
#target_link_libraries(BRANSON dmapp)

enable_testing()

add_subdirectory(test)
# Define variables needed to compile and link unit test executables
get_directory_property(Boost_INCLUDE_DIR DIRECTORY test VARIABLES Boost_INCLUDE_DIR)
get_directory_property(VIZ_LIBRARIES_FOUND test VARIABLES VIZ_LIBRARIES_FOUND)

#install(TARGETS BRANSON DESTINATION bin)
