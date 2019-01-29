Branson
----------------

[![Linux Build Status](https://travis-ci.org/lanl/branson.svg?branch=develop)](https://travis-ci.org/lanl/branson)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/xxxx/branch/develop?svg=true)](https://ci.appveyor.com/project/lanl/branson)
[![codecov.io](https://codecov.io/github/lanl/branson/coverage.svg?branch=develop)](https://codecov.io/github/lanl/branson/branch/develop)
[![Latest Version](https://img.shields.io/github/release/lanl/branson.svg?style=flat-square)](https://github.com/lanl/branson/releases)
[![PyPI](https://img.shields.io/pypi/l/Django.svg)](https://github.com/lanl/branson/blob/develop/LICENSE.md)

## Introduction

So you've decided to use Branson...

Here are some things to know:

- Branson is not an acronym.
- The point of Branson is to study different algorithms for parallel Monte Carlo
  transport. Currently it contains particle passing and mesh passing methods for
  domain decomposition.
- Many of the parameters that impact parallel performance can be set in the
  input file.
- The input file format is ugly right now and I don't have a good list of
  parameters and what they do. Sorry.
- Input files are in XML, which makes them easy to generate and change in
  python.
- Input files are complicated when you want to have multiple spatial regions but
  are pretty simple otherwise.

## How to install

Accessing the sources

- Download a tarball from Github, or
- Fork and clone the git repository.
```
# after forking the repository
git clone git@github.com:[github-username]/branson
cd branson
git remote add upstream git@github.com:lanl/branson
git checkout -b develop upstream/develop
```

Installing Branson:

- Build requirements:
  - C/C++ compiler(s) with support for C11 and C++14.
  - [CMake 3.9+](https://cmake.org/download/)
  - MPI 3.0+ ([OpenMPI 1.10+](https://www.open-mpi.org/software/ompi/),
    [mpich](http://www.mpich.org), etc.)
  - [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
- Optional dependencies (needed for visualization)
  - [HDF5](https://support.hdfgroup.org/HDF5/)
  - [Silo](http://wci.llnl.gov/simulation/computer-codes/silo)
  - If these tools aren't readily avilable on your target machine, consider
    using a package manager like [spack](https://github.com/spack/spack) to help
    you install these tools.
- There is only one CMake user option right now: `CMAKE_BUILD_TYPE` which can be
  set on the command line with `-DCMAKE_BUILD_TYPE=<Debug|Release>` and the
  default is Release.
- If cmake has trouble finding your installed TPLs, you can try
  - appending their locations to `CMAKE_PREFIX_PATH`,
  - Setting helper variables like `HDF5_ROOT` (refer to the
    [cmake
    documentation](https://cmake.org/cmake/help/latest/module/FindHDF5.html?highlight=findhdf5)
    or the `FindXXX.cmake` scripts in Branson's `src/config` directory for a
    list of variables), or
  - try running `ccmake .` from the build directory and changing the values of
    build system variables related to TPL locations.
```
EXPORT CXX=`which g++`
cd $build_dir
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<install-location> ${branson_source_dir}/src
make -j
```

Testing the build:

```
cd $build_dir
ctest -j 32
```

## Special builds

### Fake Multigrop Branson:

- The clean_mg branch has a new capability that mimics the data flow in real
  multigroup transport.
- This branch sets the number of groups at the configure step in CMake. I know
  this is messy from a usability standpoint but it makes the data layout very
  easy to control when the number of groups is known at compile time. Use the
  `N_GROUPS` CMake variable to set the number of groups (e.g. `cmake
  -DN_GROUPS=10 ../path/to/CMakeLists.txt`).
- The code will still produce gray results! The physical value in each group is
  the same and it's still set via the input deck.
- This branch samples a group with a uniform PDF (it does not weight the opacity
  with a Planckian spectrum).

## Authors

- Branson was written by Alex R. Long
- [Random123](http://www.deshawresearch.com/resources_random123.html)
  is by D. E. Shaw Research, Copyright 2010-2011
- [PugiXML](https://github.com/zeux/pugixml)
  is by zeux (Arseny Kapoulkine), MIT License
- RNG.h uses code from [Draco](https://github.com/lanl/Draco)

## Release

Branson is released under the BSD 3-Clause License. For more details see the
[LICENSE.md file](https://github.com/lanl/branson/blob/develop/LICENSE.md).

LANL code designation: `C17048`
