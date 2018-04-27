
Branson
----------------

So you've decided to use Branson...

Here are some things to know:

- Branson is not an acronym.
- Branson requires MPI 3.0+, Boost (headers only), Metis and ParMetis.
- The point of Branson is to study different algorithms for parallel
 Monte Carlo transport. Currently it contains particle passing and
 mesh passing methods for domain decomposition.
- Many of the parameters that impact parallel performance can be set in the
 input file.
- The input file format is ugly right now and I don't have a good list of
 parameters and what they do. Sorry.
- Input files are in XML, which makes them easy to generate and change in
 python.
- Input files are complicated when you want to have multiple spatial regions
 but are pretty simple otherwise.

Installing Branson:

- Branson uses CMake. The main CMakelist file is in the src directory.
- There is only one CMake user option right now: "CMAKE_BUILD_TYPE" which
 can be set on the command line with "DCMAKE_BUILD_TYPE=<Debug|Release>" and
 the default is release.
- Visualization is optional and will be turned on if HDF5 and SILO libraries
are found by CMake

Fake Multigrop Branson:

- The clean_mg branch has a new capability that mimics the data flow in real
 multigroup transport.
- This branch sets the number of groups at the configure step in CMake. I know 
  this is  messy from a usability standpoint but it makes the data layout very 
  easy to control when the number of groups is known at compile time. Use
  the "N_GROUPS" CMake variable to set the number of groups (e.g. cmake
  -DN_GROUPS=10 ../path/to/CMakeLists.txt).
- The code will still produce gray results! The physical value in each group 
  is the same and it's still set via the input deck.
- This branch samples a group with a uniform PDF (it does not weight the
  opacity with a Planckian spectrum).

Authors
----------------
- Branson was written by Alex R. Long
- random123 is by D. E. Shaw Research, Copyright 2010-2011
- RNG.h uses code from [Draco](https://github.com/lanl/Draco)

Release
----------------
Branson is released under the BSD 3-Clause License. For more details see the
LICENSE.txt file.

LANL code designation: `C17048`
