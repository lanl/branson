//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_mpi_types.cc
 * \author Alex Long
 * \date   May 12 2016
 * \brief  Test custom MPI types for consistency with their datatypes
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include "../cell.h"
#include "../mpi_types.h"
#include "../photon.h"
#include "../work_packet.h"
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  int nfail = 0;

  // test get and size functions
  {
    MPI_Types mpi_types;

    bool size_functions_pass = true;

    // test proto_cell size
    MPI_Datatype MPI_Proto_Cell = mpi_types.get_proto_cell_type();

    int proto_cell_size;
    MPI_Type_size(MPI_Proto_Cell, &proto_cell_size);

    // copy should be the same size as size recorded in class
    if (proto_cell_size != mpi_types.get_proto_cell_size())
      size_functions_pass = false;

    // copy should be the same size as actual Proto_Cell class
    if (proto_cell_size != sizeof(Proto_Cell))
      size_functions_pass = false;

    cout << "Proto_Cell object size :" << sizeof(Proto_Cell) << endl;
    cout << "MPI Proto_Cell object size :" << proto_cell_size << endl;

    if (size_functions_pass)
      cout << "TEST PASSED: MPI_Types size functions " << endl;
    else {
      cout << "TEST FAILED: MPI_Types size functions" << endl;
      nfail++;
    }
  }

  MPI_Finalize();

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_mpi_types.cc
//---------------------------------------------------------------------------//
