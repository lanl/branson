###############################################################################
# NOTES ON COMPILATION
###############################################################################
# Load the modules for each libary and the environment variables will be set 
# correctly

cmake_minimum_required (VERSION 2.8)
project (BRANSON)

###############################################################################
site_name( SITENAME )
string( REGEX REPLACE "([A-z0-9]+).*" "\\1" SITENAME ${SITENAME} )
if( ${SITENAME} MATCHES "c[it]" )
  set( SITENAME "Cielito" )
elseif( ${SITENAME} MATCHES "ml" OR ${SITENAME} MATCHES "lu" )
  set( SITENAME "Moonlight" )
elseif( ${SITENAME} MATCHES "tt") #" -login[0-9]+" OR ${SITENAME} MATCHES "tt-fey[0-9]+" )
  set( SITENAME "Trinitite" )
elseif( ${SITENAME} MATCHES "tr") #" -login[0-9]+" OR ${SITENAME} MATCHES "tr-fe[0-9]+" )
  set( SITENAME "Trinity" )
elseif( ${SITENAME} MATCHES "ccscs[0-9]+" )
  # do nothing (keep the fullname)
endif()
set( SITENAME ${SITENAME} CACHE "STRING" "Name of the current machine" FORCE)
message("Machine name: ${SITENAME}")


# set compiler
if ("${SITENAME}" STREQUAL "Moonlight")
  set(CMAKE_CXX_COMPILER mpic++)
elseif ("${SITENAME}" STREQUAL "Trinitite")
  set(CMAKE_CXX_COMPILER cc)
endif ()
message ("Compiler: " ${CMAKE_CXX_COMPILER} )

# add compiler flags
set( DEBUG_FLAGS -gdwarf-3)

# Intel C++ 11 flags
set( DEBUG_FLAGS "${DEBUG_FLAGS} -Wall -std=c++0x")

# Cray flag for C++11
#set (DEBUG_FLAGS "${DEBUG_FLAGS} -h std=c++11")

# append optimization flags
#set( DEBUG_FLAGS "${DEBUG_FLAGS} -O3")

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
message("Boost from: " $ENV{BOOST_ROOT})
set(Boost_INCLUDE_DIR $ENV{BOOST_ROOT}/include)

include_directories(${Boost_INCLUDE_DIR})

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
message("HDF5_FOUND = ${HDF5_FOUND}")
message("HDF from: $ENV{HDF5_ROOT} ")
link_directories($ENV{HDF5_ROOT}/lib)

if (HDF5_FOUND AND SILO_FOUND)
  set(VIZ_LIBRARIES_FOUND TRUE)
else ()
  message("Visualization libraries not loaded...skipping")
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
get_directory_property(SITENAME test VARIABLES SITENAME)

#install(TARGETS BRANSON DESTINATION bin)
