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
  # OpenMP
  ##############################################################################
  message(STATUS "Looking for Threads...")
  find_package(Threads QUIET)
  if(Threads_FOUND)
    message(STATUS "Looking for Threads...found")
  else()
    message(STATUS "Looking for Threads...not found")
  endif()

  message(STATUS "Looking for OpenMP...")
  if(DEFINED USE_OPENMP)
    # no-op (use defined value, -DUSE_OPENMP=<OFF|ON>,  instead of attempting to guess)
  else()
    # Assume we want to use it if it is found.
    set(USE_OPENMP ON)
  endif()
  set(USE_OPENMP
      ${USE_OPENMP}
      CACHE BOOL "Enable OpenMP threading support if detected." FORCE)

  # Find package if desired:
  if(USE_OPENMP)
    find_package(OpenMP)
  else()
    set(OpenMP_FOUND FALSE)
  endif()

  if(OpenMP_CXX_FLAGS)
    # [2022-10-27 KT] cmake/3.22 doesn't report OpenMP_C_VERSION for nvc++. Fake it for now.
    if("${OpenMP_C_VERSION}x" STREQUAL "x" AND CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
      set(OpenMP_C_VERSION
          "5.0"
          CACHE BOOL "OpenMP version." FORCE)
      set(OpenMP_FOUND TRUE)
    endif()
    message(STATUS "Looking for OpenMP... ${OpenMP_C_FLAGS} (supporting the ${OpenMP_C_VERSION} "
                   "standard)")
    # if(OpenMP_C_VERSION VERSION_LESS 3.0)
    #   message(STATUS "OpenMP standard support is too old (< 3.0). Disabling OpenMP build features.")
    #   set(OpenMP_FOUND FALSE)
    #   set(OpenMP_C_FLAGS
    #       ""
    #       CACHE BOOL "OpenMP disabled (too old)." FORCE)
    # endif()
    set(OpenMP_FOUND
        ${OpenMP_FOUND}
        CACHE BOOL "Is OpenMP available?" FORCE)
  else()
    if(USE_OPENMP)
      # Not detected, though desired.
      message(STATUS "Looking for OpenMP... not found")
    else()
      # Detected, but not desired.
      message(STATUS "Looking for OpenMP... found, but disabled for this build")
    endif()
  endif()
  message("CXX OPEN MP FLAGS: ${OpenMP_CXX_FLAGS}")

  ##############################################################################
  # Caliper
  ##############################################################################
  if( NOT TARGET caliper )
    #=============================================================================
    # If the user has provided ``CALIPER_ROOT_DIR``, use it!  Choose items found
    # at this location over system locations.
    if( EXISTS "$ENV{CALIPER_ROOT_DIR}" )
      file( TO_CMAKE_PATH "$ENV{CALIPER_ROOT_DIR}" CALIPER_ROOT_DIR )
      set( CALIPER_ROOT_DIR "${CALIPER_ROOT_DIR}" CACHE PATH
        "Prefix for Caliper installation." )
    endif()

    message( STATUS "Looking for caliper..." )
    find_package( caliper QUIET)
    if( caliper_FOUND )
      message( STATUS "Looking for caliper.....found ${CALIPER_LIBRARY}" )
    else()
      message( STATUS "Looking for caliper.....not found" )
    endif()

    set_package_properties( caliper PROPERTIES
      DESCRIPTION "CALIPER"
      TYPE OPTIONAL
      URL "https://software.llnl.gov/Caliper/"
      PURPOSE "Code instrumentation for performance analysis"
   )

  endif()

  ##############################################################################
  # metis
  # Load modules for metis to get correct environment variables
  ##############################################################################
  if( NOT TARGET METIS::metis )

    message( STATUS "Looking for METIS..." )
    find_package( METIS QUIET)
    if( METIS_FOUND )
      message( STATUS "Looking for METIS.....found ${METIS_LIBRARY}" )
    else()
      message( STATUS "Looking for METIS.....not found" )
    endif()

    set_package_properties( METIS PROPERTIES
      DESCRIPTION "METIS"
      TYPE OPTIONAL
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
