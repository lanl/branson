#!/bin/bash -l
##---------------------------------------------------------------------------##
## File  : scripts/travis-run-tests.sh
## Date  : Tuesday, Jan 17, 2017, 15:55 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017-2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# .travis.yml calls this script to build draco and run the tests.

# preliminaries and environment
set -e

cd ${SOURCE_DIR:-/home/travis/branson}

#------------------------------------------------------------------------------#
# Helper functions

function die () { echo " "; echo "FATAL ERROR: $1"; exit 1;}

# Echo a command and then run it.
function run ()
{
  echo "==> $1"; if test ${dry_run:-no} = "no"; then eval $1; fi
}

# Return 0 if provided name is a bash function.
function fn_exists ()
{
  type $1 2>/dev/null | grep -c 'is a function'
}

# Return integer > 0 if 'develop' branch is found.
function find_dev_branch
{
  set -f
  git branch -a | grep -c develop
  set +f
}

#------------------------------------------------------------------------------#
# The main script
#------------------------------------------------------------------------------#

if [[ ${STYLE} ]]; then

  #----------------------------------------------------------------------------#
  # Code style conformance
  #----------------------------------------------------------------------------#

  echo "checking style conformance..."

  # Ensure the 'develop' branch is available.  In some cases (merge a branch
  # that lives at github.com/lanl), the develop branch is missing in the
  # travis checkout. Since we only test files that are modified when comapred to
  # the 'develop' branch, the develop branch must be available locally.
  num_dev_branches_found=`find_dev_branch`
  if [[ $num_dev_branches_found == 0 ]]; then
    echo "no develop branches found."
    # Register the develop branch in draco/.git/config
    run "git config --local remote.origin.fetch +refs/heads/develop:refs/remotes/origin/develop"
    # Download the meta-data for the 'develop' branch
    run "git fetch"
    # Create a local tracking branch
    run "git branch -t develop origin/develop"
  fi

  # clang-format is installed at /usr/bin.
  export PATH=$PATH:/usr/bin
  # extract the TPL list from the Dockerfile
  export CLANG_FORMAT_VER="`grep \"ENV CLANG_FORMAT_VER\" scripts/Dockerfile | sed -e 's/.*=//' | sed -e 's/\"//g'`"
  scripts/check_style.sh -t

else

  #----------------------------------------------------------------------------#
  # Try to build the code and run the tests
  #----------------------------------------------------------------------------#

  echo "checking build and test..."

  # extract the TPL list from the Dockerfile
  export DRACO_TPL="`grep \"ENV DRACO_TPL\" scripts/Dockerfile | sed -e 's/.*=//' | sed -e 's/\"//g'`"

  # Environment setup for the build...
  spack load ${DRACO_TPL} boost hdf5

  if [[ $GCCVER ]]; then
    export CXX=`which g++-${GCCVER}`
    export CC=`which gcc-${GCCVER}`
    export FC=`which gfortran-${GCCVER}`
    export GCOV=`which gcov-${GCCVER}`
  else
    export CXX=`which g++`
    export CC=`which gcc`
    export FC=`which gfortran`
    export GCOV=`which gcov`
  fi
  echo "GCCVER = ${GCCVER}"
  echo "CXX    = ${CXX}"
  echo "FC     = ${FC}"
#  ls -1 /usr/bin/g*
#  if ! [[ $FC ]]; then export FC=gfortran; fi
  # ${FC} --version

  export OMP_NUM_THREADS=2
  if [[ ${WERROR} ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" -Werror\"
    done
  fi
  if [[ ${COVERAGE} ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" --coverage\"
    done
  fi

  # echo " "
  # echo "========== printenv =========="
  # printenv
  # echo " "

  if ! [[ $BUILD_DIR ]]; then die "BUILD_DIR not set by environment."; fi
  run "mkdir -p ${BUILD_DIR}"
  run "cd ${BUILD_DIR}"

  echo " "
  if [[ -f CMakeCache.txt ]]; then
    echo "===== CMakeCache.txt ====="
    run "cat CMakeCache.txt"
  fi
  echo "========"
  if [[ ${DEBUG} == ON ]]; then
    run "cmake -DCMAKE_BUILD_TYPE=Debug ${SOURCE_DIR}/src"
  else
    run "cmake ${SOURCE_DIR}/src"
  fi
  echo "========"
  run "make -j 2"
  echo "========"
  run "ctest -j 2 --output-on-failure"
  cd -
  if [[ ${COVERAGE} ]]; then
    echo "========"
    #which codecov
    #run "codecov --gcov-exec $GCOV"
    # https://docs.codecov.io/docs/testing-with-docker
    /bin/bash <(curl -s https://codecov.io/bash)
  fi
  echo "======== end scripts/travis-run-tests.sh =========="
fi

#------------------------------------------------------------------------------#
# End scripts/travis-run-tests.sh
#------------------------------------------------------------------------------#
