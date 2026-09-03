################################################################################
# NOTES ON COMPILATION
################################################################################

# Load the modules for each library and the environment variables will be set
# correctly (load cmake, compiler, mpi, metis, and optionally
# hdf5 and silo)

# Enable debug mode by passing -DCMAKE_BUILD_TYPE=Debug to CMake, default is
# Release

cmake_minimum_required (VERSION 3.24)

option(USE_UMPIRE "Use Umpire" FALSE)
option(USE_CUDA "Use CUDA" FALSE)
option(USE_HIP "Use HIP" FALSE)
option(USE_GPU "user-set flag to compile in GPU code" FALSE)
option(USE_CALIPER "Use Caliper" FALSE)
option(USE_OPENMP "Use OpenMP" FALSE)
option(BRANSON_USE_CLANG_CUDA "Use Clang instead of NVCC for CUDA compilation" FALSE)

if(USE_CUDA)
    set(GPU_DBS_STRING "CUDA" CACHE STRING "If CUDA is available, this variable is 'CUDA'")
    # Set the CUDA compiler based on the toggle option BEFORE enabling the language
    if(BRANSON_USE_CLANG_CUDA)
      set(CMAKE_CUDA_COMPILER "clang" CACHE STRING "CUDA compiler" FORCE)
      message(STATUS "Configuring CUDA build to use Clang")
    else()
      set(CMAKE_CUDA_COMPILER "nvcc" CACHE STRING "CUDA compiler" FORCE)
      message(STATUS "Configuring CUDA build to use NVCC")
    endif()
    enable_language(CUDA)
elseif(USE_HIP)
    set(GPU_DBS_STRING "HIP" CACHE STRING "If HIP is available, this variable is 'HIP'")
    enable_language(HIP)
elseif(DEFINED ENV{CUDADIR}
   OR DEFINED ENV{CUDACXX}
   OR DEFINED ENV{CUDA_HOME})
  if(USE_GPU)
    set(GPU_DBS_STRING "CUDA" CACHE STRING "If CUDA is available, this variable is 'CUDA'")
    if(BRANSON_USE_CLANG_CUDA)
      set(CMAKE_CUDA_COMPILER "clang" CACHE STRING "CUDA compiler" FORCE)
    else()
      set(CMAKE_CUDA_COMPILER "nvcc" CACHE STRING "CUDA compiler" FORCE)
    endif()
    enable_language(CUDA)
  endif()
elseif(
  "$ENV{LOADEDMODULES}" MATCHES "rocmcc"
  OR DEFINED ENV{HIPCXX}
  OR DEFINED CMAKE_HIP_COMPILER
  OR DEFINED ENV{ROCM_PATH})
  if(USE_GPU)
    set(GPU_DBS_STRING "HIP" CACHE STRING "If HIP is available, this variable is 'HIP'")
    #if(CMAKE_HIP_COMPILER MATCHES "CC")
    #  set(CMAKE_HIP_FLAGS "-x hip")
    #endif()
    enable_language(HIP)
  endif()
endif()

project (BRANSON
  VERSION 0.8
  DESCRIPTION "Branson can be used to study different algorithms for parallel Monte Carlo transport. Currently it contains particle passing and mesh passing methods for domain decomposition."
#  HOMEPAGE URL "https://github.com/lanl/branson" # needs cmake 3.12+
  LANGUAGES CXX C ${GPU_DBS_STRING})

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)

find_package(MPI REQUIRED C CXX)

set(branson_mpi_include_dirs
  ${MPI_C_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_DIRS}
  ${MPI_C_ADDITIONAL_INCLUDE_DIRS}
  ${MPI_CXX_ADDITIONAL_INCLUDE_DIRS}
  ${MPI_C_COMPILER_INCLUDE_DIRS}
  ${MPI_CXX_COMPILER_INCLUDE_DIRS})
list(REMOVE_DUPLICATES branson_mpi_include_dirs)

list(APPEND branson_deps
  MPI::MPI_C
  MPI::MPI_CXX)

message("Languages: ${_LANGUAGES_}")
if((_LANGUAGES_ MATCHES CUDA OR _LANGUAGES_ MATCHES HIP) AND USE_GPU)
  message("CUDA/HIP module found (CUDA/HIP environment variables set) and USE_GPU is on, making GPU BRANSON")
  add_compile_definitions(HAS_GPU)
elseif(USE_GPU)
  message(FATAL_ERROR "CUDA/HIP module NOT found (CUDA/HIP environment variables set) but USE_GPU is on, reconfigure with USE_GPU off or fix modules")
elseif(NOT USE_GPU)
  message("GPU mode not requested, making CPU only BRANSON")
else()
  message("CUDA/HIP module not found and GPU mode not requested, making CPU only BRANSON")
endif()


# Build system support files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})

message( "
This is ${CMAKE_PROJECT_NAME} v. ${PROJECT_VERSION}.
")

#------------------------------------------------------------------------------#
# User configurable options

set (BRANSON_N_GROUPS ${N_GROUPS})
if (NOT BRANSON_N_GROUPS)
  message(STATUS
    "HARDCODED NUMBER OF GROUPS NOT SET, DEFAULTING TO ONE GROUP (GRAY)
     use '-DN_GROUPS=<int>' when running cmake to select multigroup.\n")
  set(BRANSON_N_GROUPS "1")
else ()
  message(STATUS
    "HARDCODED NUMBER OF GROUPS SET TO ${BRANSON_N_GROUPS} (FAUX MULTIGROUP)\n")
endif()
set(BRANSON_N_GROUPS ${BRANSON_N_GROUPS} CACHE STRING "Number of energy groups." FORCE)

set(ENABLE_VERBOSE_GPU_TRANSPORT ${ENABLE_VERBOSE_GPU_TRANSPORT} CACHE STRING "Flag for verbose GPU output." FORCE)

#------------------------------------------------------------------------------#
# Record the site name

site_name( SITENAME )
string( REGEX REPLACE "([A-z0-9]+).*" "\\1" SITENAME ${SITENAME} )
include (lanl-setup)
set( SITENAME ${SITENAME} CACHE "STRING" "Name of the current machine" FORCE)

#------------------------------------------------------------------------------#
# Setup compiler options
set( CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 17 )
set(CMAKE_CUDA_ARCHITECTURES 75 80 86 90)
# Do not enable extensions (e.g.: --std=gnu++11)
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_C_EXTENSIONS   OFF )

#------------------------------------------------------------------------------#
# Build type and custom compiler flags

if ( "${CMAKE_BUILD_TYPE}notset" STREQUAL "notset" AND
    NOT DEFINED CMAKE_CONFIGURATION_TYPES )
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release|Debug|RelWithDebInfo"
    FORCE)
endif ()

# GCC options for address or undefined sanitizing
#set(GCC_SANITIZE_COMPILE_FLAGS "-fsanitize=address")
#set(GCC_SANITIZE_LINK_FLAGS "-fsanitize=address")
#set(GCC_SANITIZE_COMPILE_FLAGS "-fsanitize=undefined")
#set(GCC_SANITIZE_LINK_FLAGS "-fsanitize=undefined")

# add compiler flags
if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-contract=off")
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0")
  #set( CMAKE_CXX_FLAGS_RELEASE "-march=native  -O3 -funroll-loops -fno-var-tracking-assignments")
  set( CMAKE_CXX_FLAGS_RELEASE "-g -O3 -funroll-loops -fno-var-tracking-assignments")
  # optionally add sanitize flags
  #set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_SANITIZE_COMPILE_FLAGS}")
  #set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_SANITIZE_LINK_FLAGS}")
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0")
  # note that the "fp-model=precise" flag significantly limits vectorization and slows down
  # branson but it makes it the hardware counter data match Jayenne more closely (Jayenne uses
  # this flag)
  #set( CMAKE_CXX_FLAGS_RELEASE "-g -O3 -fp-speculation=fast -xhost -qno-opt-dynamic-align")
       set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -xhost -qopenmp -qopt-zmm-usage=high -fp-speculation=fast -qopt-report=3 -qopt-report-file=runAVX512.optrpt -no-vec") # "-no-vec" to disable vectorization
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -O0 -xhost -qopenmp")
       #set(CMAKE_CXX_FLAGS_RELEASE "-g ") # "-no-vec" to disable vectorization
        set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -xhost -qopenmp -qopt-zmm-usage=high -fp-speculation=fast -qopt-report=3 -qopt-report-file=runAVX512.optrpt") # "-no-vec" to disable vectorization

  elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Cray" )
  set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -g -DR123_USE_GNU_UINT128=0")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DR123_USE_GNU_UINT128=0")
endif()

##############################################################################
# umpire
##############################################################################
if(USE_CUDA OR USE_HIP)
  if(USE_UMPIRE)
    message(STATUS "Looking for umpire...")
    find_package(umpire REQUIRED)
    message(STATUS "Looking for umpire......found ${umpire_DIR}")
    list(APPEND branson_deps umpire)
  endif()
endif()


#------------------------------------------------------------------------------#
# Look for Third Party Libraries (metis, etc.)
add_subdirectory(pugixml)

include(find_tpls)
setupTPLs()

#------------------------------------------------------------------------------#
# Report build configuration

# Summary of TPLs
include(FeatureSummary)
message(" ")
feature_summary( WHAT ALL
  INCLUDE_QUIET_PACKAGES
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  QUIET_ON_EMPTY )

# Summary of user-selectable build options
message( "\nBuild Summary:\n")
message( STATUS "Machine name    : ${SITENAME}")
report_lanl_hpc_features()
if( ${BRANSON_N_GROUPS} GREATER 1 )
  message( STATUS "Energy groups   : ${BRANSON_N_GROUPS}")
else()
  message( STATUS "Energy groups   : Gray (1-group)")
endif()
message( STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message( STATUS "Compiler        : ${CMAKE_CXX_COMPILER}")
if( CMAKE_CONFIGURATION_TYPES )
  message( STATUS "Compiler Flags (All)    : ${CMAKE_CXX_FLAGS}")
  message( STATUS "Compiler Flags (Debug)  : ${CMAKE_CXX_FLAGS_DEBUG}")
  message( STATUS "Compiler Flags (Release): ${CMAKE_CXX_FLAGS_RELEASE}")
else()
  string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )
  message( STATUS "Compiler Flags  : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
endif()
message("\n")

##############################################################################
# caliper/adiak
##############################################################################
if(USE_CALIPER)
  message(STATUS "Looking for caliper...")
  find_package(caliper REQUIRED)
  if(caliper_FOUND)
    message(STATUS "Looking for caliper......found ${caliper_DIR}")
    list(APPEND branson_deps
      caliper)
  endif()

  message(STATUS "Looking for adiak...")
  find_package(adiak REQUIRED)
  if(adiak_FOUND)
    message(STATUS "Looking for adiak......found ${adiak_DIR}")
    list(APPEND branson_deps
      adiak::adiak)
  endif()
else()
  message(STATUS "Disabling caliper...")
endif()

#------------------------------------------------------------------------------#
# Set up libraries for Cray RMA routines

# set(DMAPP_DYNAMIC -Wl,--whole-archive,-ldmapp,--no-whole-archive)


#------------------------------------------------------------------------------#
# Targets
file(GLOB headers *.h)
add_executable(BRANSON main.cc ${headers})
target_include_directories( BRANSON PRIVATE
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  ${PROJECT_SOURCE_DIR}/pugixml/src/
  ${PROJECT_SOURCE_DIR}/random123/
  ${branson_mpi_include_dirs} )

if("${GPU_DBS_STRING}" STREQUAL "CUDA" )
    find_package(CUDAToolkit REQUIRED)
    list(APPEND branson_deps
      CUDA::cudart)

    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  if(CUDA_ARCH)
    set(CMAKE_CUDA_ARCHITECTURES "${CUDA_ARCH}")
  else()

    # Set CUDA compute capabilities to use host architecture (read: must configure on GPU node!)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      # This value is automatically added to all CUDA targets. See cmake/component_macros.cmake
      set(CMAKE_CUDA_ARCHITECTURES
          "native"
          CACHE STRING "target architecture for gpu code.")
    endif()
    message("CMAKE_CUDA_ARCHITECTURES = ${CUDA_ARCHITECTURES}")

  endif()

  # Conditional flags depending on whether we are using Clang or NVCC for CUDA
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    string(APPEND CMAKE_CUDA_FLAGS " -g -fno-gpu-rdc")
    message("Making GPU(CUDA) BRANSON using Clang with architecture(s) ${CMAKE_CUDA_ARCHITECTURES}")
  else()
    string(APPEND CMAKE_CUDA_FLAGS " -g --expt-relaxed-constexpr --fmad=false -lineinfo --expt-extended-lambda")
    message("Making GPU(CUDA) BRANSON using NVCC with architecture(s) ${CMAKE_CUDA_ARCHITECTURES}")
  endif()

  set_source_files_properties("main.cc" PROPERTIES LANGUAGE CUDA)
  target_compile_definitions(BRANSON PRIVATE USE_GPU)

elseif("${GPU_DBS_STRING}" STREQUAL "HIP" )
  if(USE_HIP)
    if(NOT ROCM_PATH)
      if(DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
      else()
        message(FATAL_ERROR "Unable to find rocm. Please use "
                            "ROCM_PATH to provide the rocm installation.")
      endif()
    endif()

    set(CMAKE_HIP_COMPILER "${ROCM_PATH}/llvm/bin/clang++")
    list(APPEND CMAKE_PREFIX_PATH "${ROCM_PATH}/lib/cmake")

    find_package(hip REQUIRED)
    list(APPEND branson_deps
      hip::host
      hip::device)

    set(CMAKE_HIP_HOST_COMPILER ${CMAKE_CXX_COMPILER})

    set(CMAKE_HIP_STANDARD ${CMAKE_CXX_STANDARD})
    set(CMAKE_HIP_STANDARD_REQUIRED ON)
    set(CMAKE_HIP_ARCHITECTURES "${HIP_ARCH}")
    string(APPEND CMAKE_HIP_FLAGS "-g")
    set_source_files_properties("main.cc" PROPERTIES LANGUAGE HIP)
    set(CMAKE_HIP_FLAGS "-g ${CMAKE_HIP_FLAGS}")
    target_compile_definitions(BRANSON PRIVATE USE_GPU)
    message("Making GPU(HIP) BRANSON")
  else()
    message("Setting HIP compiler options")
    if(CMAKE_HIP_COMPILER MATCHES "CC")
      set(CMAKE_HIP_FLAGS "-x hip")
    endif()
    set_target_properties(BRANSON PROPERTIES HIP_ARCHITECTURES "gfx942") # MI300
    set_target_properties(BRANSON PROPERTIES HIP_STANDARD 17)
    string(APPEND CMAKE_HIP_FLAGS " -g -mllvm=--disable-peephole")
    #string(APPEND CMAKE_HIP_FLAGS " --expt-extended-lambda"  )
    set_source_files_properties("main.cc" PROPERTIES LANGUAGE HIP)
  endif()
else()
  message("GPU Options: Not a GPU build or GPU_DBS_STRING not recognized")
endif()

list(APPEND branson_deps pugixml)

if(OpenMP_FOUND)
  set(branson_deps "OpenMP::OpenMP_CXX;${branson_deps}")
elseif(USE_OPENMP AND CMAKE_CXX_COMPILER_WRAPPER STREQUAL "CrayPrgEnv")
  target_compile_options(BRANSON PRIVATE -fopenmp)
  target_link_options(BRANSON PRIVATE -fopenmp)
elseif(USE_OPENMP)
  message(FATAL_ERROR "USE_OPENMP=ON but CMake could not find OpenMP")
endif()

if(METIS_FOUND)
  set( branson_deps "METIS::metis;${branson_deps}")
endif()

if( VIZ_LIBRARIES_FOUND )
  set( branson_deps "Silo::silo;${HDF5_LIBRARIES};${branson_deps}" )
endif()

#------------------------------------------------------------------------------#
# Generate config.h

configure_file(config.h.in ${PROJECT_BINARY_DIR}/config.h)
#------------------------------------------------------------------------------#

target_link_libraries( BRANSON PRIVATE ${branson_deps} )

#------------------------------------------------------------------------------#
# Testing

option( BUILD_TESTING "Should we compile the tests?" ON )
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
  message(STATUS
    "Building tests enabled (default), disable tests with BUILD_TESTING=FALSE")
else()
  message(STATUS
    "Building tests disabled, set BUILD_TESTING=TRUE or don't set BUILD_TESTING to enable test builds")
endif()


#------------------------------------------------------------------------------#
# Targets for installation

install(TARGETS BRANSON DESTINATION ${CMAKE_INSTALL_BINDIR})

#------------------------------------------------------------------------------#
# End src/CMakeLists.txt
#------------------------------------------------------------------------------#
