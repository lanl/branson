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

  ##############################################################################
  # metis
  # Load modules for metis to get correct environment variables
  ##############################################################################
  if( NOT TARGET METIS::metis )

    message( STATUS "Looking for METIS..." )
    find_package( METIS QUIET REQUIRED )
    if( METIS_FOUND )
      message( STATUS "Looking for METIS.....found ${METIS_LIBRARY}" )
    else()
      message( STATUS "Looking for METIS.....not found" )
    endif()

    set_package_properties( METIS PROPERTIES
      DESCRIPTION "METIS"
      TYPE REQUIRED
      URL "http://glaros.dtc.umn.edu/gkhome/metis/metis/overview"
      PURPOSE "METIS is a set of serial programs for partitioning graphs,
   partitioning finite element meshes, and producing fill reducing orderings for
   sparse matrices.")

  endif()

  ##############################################################################
  # Silo and HDF5 libraries
  # Load modules for hdf5 and solo to get correct environment variables
  # use find package
  ##############################################################################

  if( NOT HDF5_FOUND )

    message( STATUS "Looking for HDF5..." )
    find_package( HDF5 QUIET )
    if( HDF5_FOUND )
      list(GET HDF5_LIBRARIES 0 hdf5lib)
      message( STATUS "Looking for HDF5..found ${hdf5lib}" )
      unset(hdf5lib)
    else()
      message( STATUS "Looking for HDF5..not found" )
    endif()

    set_package_properties( HDF5 PROPERTIES
      DESCRIPTION "HDF5 is a data model, library, and file format for storing
   and managing data. It supports an unlimited variety of datatypes, and is
   designed for flexible and efficient I/O and for high volume and complex
   data."
      TYPE OPTIONAL
      URL "https://support.hdfgroup.org/HDF5/"
      PURPOSE "Provides optional visualization support for Branson." )

  endif()

  if( HDF5_FOUND AND NOT TARGET Silo::silo )

    message( STATUS "Looking for Silo..." )
    find_package( Silo QUIET )
    if( Silo_FOUND )
      message( STATUS "Looking for Silo..found ${Silo_LIBRARY}" )
    else()
      message( STATUS "Looking for Silo..not found" )
    endif()

    set_package_properties( Silo PROPERTIES
      DESCRIPTION "Silo is a library for reading and writing a wide variety of
   scientific data to binary, disk files."
      TYPE OPTIONAL
      URL "http://wci.llnl.gov/simulation/computer-codes/silo"
      PURPOSE "Provides optional visualization support for Branson.")

  endif()

  if (HDF5_FOUND AND Silo_FOUND)
    set(VIZ_LIBRARIES_FOUND TRUE)
  else ()
    message(STATUS "Optional visualization libraries not loaded...skipping")
  endif ()

endmacro()

#------------------------------------------------------------------------------#
# End find_tpls.cmake
#------------------------------------------------------------------------------#
