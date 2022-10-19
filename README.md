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
  transport. Currently it contains a particle passing method for domain decomposition (two mesh
  passing methods for domain decomposition are available in older versions)
- Many of the parameters that impact parallel performance can be set in the input file.
- Input files are in XML, which makes them easy to generate and change in python.
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

## Parameters ##

- In the `common` block of the XML input file you can set several parameters related to parallel
 performance:
  - `batch_size`: how many particles are processed by an MPI rank before it stops to process
    particle messages (small means more processing of MPI messages and possibly more interleaving
    of off rank and on rank work)
  - `particle_message_size`: the size of the communication buffer for particle data (small means
     more message but possibly more interleaving of on rank and off rank work)
  - `mesh_decomposition`: Can be `METIS` or `CUBE`, generally use Metis unless you're trying to run
    a very large problem (Metis is serial and ParMetis can't be used due to licensing restrictions).
    For a cube decomposition, the number of ranks must be perfect cubes (x^(1/3) is an interger)

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

## Running Branson on performance problems

- There are two performance problems of interest in the `inputs` folder, they are both simplified
 3D hohlraums and should be run with a 30 group build of Branson (see Special builds section above).
- The `3D_hohlaum_single_node.xml` problem is meant to be run on a full node. It is in replicated
 mode which means there is very little MPI communication (end of cycle reductions). It is run with:
```
mpirun -n <procs_on_node> <path/to/branson> 3D_hohlaum_single_node.xml
```
- The `3D_hohlaum_multi_node.xml` problem is meant to be run on many nodes. It is in domain
  decomposed mode which means particles must be communicated between spatial domains. It is run
  with:
```
mpirun -n <n_ranks> <path/to/branson> 3D_hohlaum_multi_node.xml
```
 Note that Branson does not currently have any threading capability so `n_ranks` should usually be
 `n_nodes*n_ranks_per_node`. This problem is meant to consume about 30\% of a 128 GB node.

## Authors

- Branson was written by Alex R. Long
- [Random123](http://www.deshawresearch.com/resources_random123.html)
  is by D. E. Shaw Research, Copyright 2010-2011
- [PugiXML](https://github.com/zeux/pugixml)
  is by zeux (Arseny Kapoulkine), [MIT License](https://github.com/zeux/pugixml/blob/master/LICENSE.md)
- RNG.h uses code from [Draco](https://github.com/lanl/Draco), [BSD-3 Clause License](https://github.com/lanl/Draco/blob/develop/LICENSE.md)

## Release

- Branson is released under the BSD 3-Clause License. For more details see the
[LICENSE.md file](https://github.com/lanl/branson/blob/develop/LICENSE.md).

- Branson development follows the development model outlined in the CS-memo-2020-11-18.pdf file,
  available in the top-level of this repo

- LANL code designation: `C17048
