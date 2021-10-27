#-----------------------------*-cmake-*----------------------------------------#
# file   src/config/find_tpls.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Aug 14, 2018, 15:24 pm
# brief  Look for third party libraries like metis
# note   Copyright (C) 2018 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include( FeatureSummary )

macro(setupTPLs)

  ##############################################################################
  # MPI
  ##############################################################################
  if( NOT TARGET MPI::MPI_C )

    message(STATUS "Looking for MPI...")
    find_package(MPI QUIET REQUIRED)
    if( ${MPI_C_FOUND} )
      message(STATUS "Looking for MPI...${MPIEXEC}")
    else()
      message(STATUS "Looking for MPI...not found")
    endif()

    set_package_properties( MPI PROPERTIES
      URL "http://www.open-mpi.org/"
      DESCRIPTION "A High Performance Message Passing Library"
      TYPE REQUIRED
      PURPOSE "A parallel communication library is required in BRANSON.")

  endif()

  ##############################################################################
  # Caliper
  ##############################################################################
  #include(ExternalProject)
  #set(CALIPER_INSTALL_DIR ${CMAKE_BINARY_DIR}/caliper)
  #ExternalProject_Add(caliper_local SOURCE_DIR "${CMAKE_SOURCE_DIR}/Caliper" PREFIX ${CALIPER_INSTALL_DIR} CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CALIPER_INSTALL_DIR}")
  #include_directories(${CALIPER_INSTALL_DIR}/include/caliper)
  #link_directories(${CALIPER_INSTALL_DIR}/lib)

endmacro()

#------------------------------------------------------------------------------#
# End find_tpls.cmake
#------------------------------------------------------------------------------#
