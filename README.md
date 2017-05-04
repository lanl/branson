
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
- There is only one cmake user option right now: "CMAKE_BUILD_TYPE" which
 can be set on the command line with "DCMAKE_BUILD_TYPE=<Debug|Release>" and
 the default is release.
- Visualization is optional and will be turned on if HDF5 and SILO libraries
are found by CMake


Authors
----------------
- Branson was written by Alex R. Long
- random123 is by D. E. Shaw Research, Copyright 2010-2011
- RNG.h uses code from [Draco](https://github.com/lanl/Draco)

Release
----------------
Branson is released under the BSD 3-Clause License. For more details see the
LICENSE.txt file.

LANL code designation: `LA-CC-17-048`
