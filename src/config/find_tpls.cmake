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
  # metis and parmetis
  # Load modules for metis and parmetis to get correct environment variables
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

  if( NOT TARGET ParMETIS::parmetis )

    message( STATUS "Looking for ParMETIS..." )
    find_package( ParMETIS QUIET REQUIRED )
    if( ParMETIS_FOUND )
      message( STATUS "Looking for ParMETIS..found ${ParMETIS_LIBRARY}" )
    else()
      message( STATUS "Looking for ParMETIS..not found" )
    endif()

    set_package_properties( ParMETIS PROPERTIES
      DESCRIPTION "MPI Parallel METIS"
      TYPE REQUIRED
      URL "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview"
      PURPOSE "ParMETIS is an MPI-based parallel library that implements a
   variety of algorithms for partitioning unstructured graphs, meshes, and for
   computing fill-reducing orderings of sparse matrices." )

endif()

  # message("ParMetis from: " $ENV{ParMETIS_ROOT_DIR})
  # include_directories( $ENV{PARMETIS_INC_DIR})
  # link_directories($ENV{PARMETIS_LIB_DIR})

  # message("Metis from: " $ENV{METIS_ROOT_DIR})
  # include_directories( $ENV{METIS_INC_DIR})
  # link_directories($ENV{METIS_LIB_DIR})

  ###############################################################################
  # boost (headers only)
  # Load boost module to get correct environement variables
  ###############################################################################
  if( NOT TARGET Boost::boost )

    message( STATUS "Looking for Boost..." )
    find_package( Boost QUIET REQUIRED )
    if( Boost_FOUND )
      message( STATUS "Looking for Boost..found ${Boost_INCLUDE_DIRS}" )
    else()
      message( STATUS "Looking for Boost..not found" )
    endif()

    set_package_properties( Boost PROPERTIES
      DESCRIPTION "Boost provides free peer-reviewed portable C++ source libraries."
      TYPE REQUIRED
      URL "http://www.boost.org"
      PURPOSE "Boost provides free peer-reviewed portable C++ source libraries,
   emphasizing libraries that work well with the C++ Standard Library.  Boost
   libraries are intended to be widely useful, and usable across a broad
   spectrum of applications. The Boost license encourages both commercial and
   non-commercial use.")

  endif()

  # message("Boost from: " $ENV{BOOST_ROOT})
  # set(Boost_INCLUDE_DIR $ENV{BOOST_ROOT}/include)
  # include_directories(${Boost_INCLUDE_DIR})

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

  # message("OPTIONAL: Silo from: $ENV{SILO_ROOT_DIR} " )
  # include_directories($ENV{SILO_ROOT_DIR}/include)
  # link_directories($ENV{SILO_ROOT_DIR}/lib)
  # if ("$ENV{SILO_ROOT_DIR}x" STREQUAL "x")
  #   set(SILO_FOUND FALSE)
  # else ()
  #   set(SILO_FOUND TRUE)
  # endif ()
  # message("OPTIONAL: SILO_FOUND = ${SILO_FOUND}")

  # find_package(HDF5)
  # message("OPTIONAL: HDF5_FOUND = ${HDF5_FOUND}")
  # if (${HDF5_FOUND} STREQUAL "TRUE")
  #   message("OPTIONAL: HDF5 from: $ENV{HDF5_ROOT_DIR}" )
  #   include_directories("$ENV{HDF5_ROOT_DIR}/include")
  #   link_directories("$ENV{HDF5_ROOT_DIR}/lib")
  # endif ()

  if (HDF5_FOUND AND Silo_FOUND)
    set(VIZ_LIBRARIES_FOUND TRUE)
  else ()
    message(STATUS "Optional visualization libraries not loaded...skipping")
  endif ()


#add_dependencies(BRANSON caliper_local)

# these lines link the Cray DMAPP library
#first, create the dynamic version of dmapp library
#target_link_libraries(BRANSON ${DMAPP_DYNAMIC})
#then link to it dynamically
#target_link_libraries(BRANSON dmapp)

#target_link_libraries(BRANSON caliper)
#target_link_libraries(BRANSON ${CALIPER_INSTALL_DIR}/lib/libcaliper.so)
# target_link_libraries(BRANSON parmetis)
# target_link_libraries(BRANSON metis)
# if (VIZ_LIBRARIES_FOUND)
#   target_link_libraries(BRANSON hdf5)
#   target_link_libraries(BRANSON siloh5)
# endif ()


endmacro()

#------------------------------------------------------------------------------#
# End find_tpls.cmake
#------------------------------------------------------------------------------#
