#.rst:
# FindSilo
# --------
#
# Find the Silo includes and libraries.
#
# Silo is a library for reading and writing a wide variety of scientific data to
# binary, disk files.
# http://wci.llnl.gov/simulation/computer-codes/silo
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If Silo is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  Silo::silo      - The main Silo library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  Silo_FOUND          - True if Silo found on the local system
#  Silo_INCLUDE_DIRS   - Location of Silo header files.
#  Silo_LIBRARIES      - The Silo libraries.
#  Silo_VERSION        - The version of the discovered Silo install.
#
# Hints
# ^^^^^
#
# Set ``Silo_ROOT_DIR`` to a directory that contains a Silo installation.
#
# This script expects to find libraries at ``$Silo_ROOT_DIR/lib`` and the Silo
# headers at ``$Silo_ROOT_DIR/include``.  The library directory may optionally
# provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# Silo installation discovered.  These variables may optionally be set to help
# this module find the correct files::
#
#  Silo_LIBRARY             - Location of the Silo library.
#  Silo_LIBRARY_DEBUG       - Location of the debug Silo library (if any).
#

#=============================================================================
# Copyright 2018 Kelly Thompson <kgt@lanl.gov>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# Prerequisite:
# Silo requires that HDF5 be found.
find_package( HDF5 QUIET )

#=============================================================================
# If the user has provided ``Silo_ROOT_DIR``, use it!  Choose items found at
# this location over system locations.
if( EXISTS "$ENV{Silo_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{Silo_ROOT_DIR}" Silo_ROOT_DIR )
  set( Silo_ROOT_DIR "${Silo_ROOT_DIR}" CACHE PATH
    "Prefix for Silo installation." )
endif()

#=============================================================================
# Set Silo_INCLUDE_DIRS and Silo_LIBRARIES. Try to find the libraries at
# $Silo_ROOT_DIR (if provided) or in standard system locations.  These
# find_library and find_path calls will prefer custom locations over standard
# locations (HINTS).  If the requested file is not found at the HINTS location,
# standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path( Silo_INCLUDE_DIR
  NAMES silo.h
  HINTS ${Silo_ROOT_DIR}/include ${Silo_INCLUDEDIR}
  PATH_SUFFIXES Release Debug
)
find_library( Silo_LIBRARY
  NAMES siloh5
  HINTS ${Silo_ROOT_DIR}/lib ${Silo_LIBDIR}
  PATH_SUFFIXES Release Debug
)
# Do we also have debug versions?
find_library( Silo_LIBRARY_DEBUG
  NAMES silo
  HINTS ${Silo_ROOT_DIR}/lib ${Silo_LIBDIR}
  PATH_SUFFIXES Debug
)
set( Silo_INCLUDE_DIRS ${Silo_INCLUDE_DIR} )
set( Silo_LIBRARIES ${Silo_LIBRARY} )

# Try to find the version.
if( NOT Silo_VERSION )
  if( EXISTS "${Silo_INCLUDE_DIRS}/silo.h" )
    file( STRINGS "${Silo_INCLUDE_DIRS}/silo.h" silo_h_version
        REGEX "define SILO_VERS_TAG" )
    string( REGEX REPLACE ".*([0-9]+)_([0-9]+)_([0-9]+)" "\\1" Silo_MAJOR
      ${silo_h_version} )
    string( REGEX REPLACE ".*([0-9]+)_([0-9]+)_([0-9]+)" "\\2" Silo_MINOR
      ${silo_h_version} )
    string( REGEX REPLACE ".*([0-9]+)_([0-9]+)_([0-9]+)" "\\3" Silo_SUBMINOR
      ${silo_h_version} )
    set( Silo_VERSION "${Silo_MAJOR}.${Silo_MINOR}.${Silo_SUBMINOR}")
  endif()
  # We might also try scraping the directory name for a regex match
  # "silo-X.X.X"
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set Silo_FOUND to TRUE if
# all listed variables are TRUE.
find_package_handle_standard_args( Silo
  FOUND_VAR
    Silo_FOUND
  REQUIRED_VARS
    Silo_INCLUDE_DIR
    Silo_LIBRARY
    HDF5_FOUND
  VERSION_VAR
    Silo_VERSION
    )

mark_as_advanced( Silo_ROOT_DIR Silo_VERSION Silo_LIBRARY
  Silo_INCLUDE_DIR Silo_LIBRARY_DEBUG Silo_USE_PKGCONFIG
  Silo_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" Silo_LIBRARY_DLL
    "${Silo_LIBRARY}" )
  string( REPLACE ".lib" ".dll" Silo_LIBRARY_DEBUG_DLL
    "${Silo_LIBRARY_DEBUG}" )
endif()

if( Silo_FOUND AND NOT TARGET Silo::silo )
  if( EXISTS "${Silo_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( Silo::silo SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
      set_target_properties( Silo::silo PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${Silo_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${Silo_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${Silo_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          "${HDF5_LIBRARYS}" )

    # If we have both Debug and Release libraries
    if( EXISTS "${Silo_LIBRARY_DEBUG_DLL}" )
      set_property( TARGET Silo::silo APPEND PROPERTY
        IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( Silo::silo PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${Silo_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${Silo_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( Silo::silo UNKNOWN IMPORTED )
    set_target_properties( Silo::silo PROPERTIES
      IMPORTED_LOCATION                 "${Silo_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${Silo_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      INTERFACE_LINK_LIBRARIES          "${HDF5_LIBRARYS}" )
  endif()
endif()
