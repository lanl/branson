//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_write_silo.cc
 * \author Alex Long
 * \date   April 15 2016
 * \brief  Test ability to write SILO files
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>

#include "../decompose_mesh.h"
#include "../mesh.h"
#include "../mpi_types.h"
#include "../write_silo.h"
#include "testing_functions.h"

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  const Info mpi_info;

  using std::cout;
  using std::endl;
  using std::string;

  int nfail = 0;

  // Test that a simple mesh (one division in each dimension) is constructed
  // correctly from the input file (simple_input.xml) and that the write
  // silo function exits without error
  {
    // put MPI_Types in scope so it's destructor is called before MPI_Finalize
    MPI_Types mpi_types;
    string filename("simple_input.xml");
    Input input(filename, mpi_types);
    IMC_Parameters imc_p(input);

    Mesh mesh(input, mpi_types, mpi_info, imc_p);

    bool silo_write_pass = true;

    // get fake vector of mesh requests
    std::vector<uint32_t> n_requests(mesh.get_n_global_cells(), 0);

    double time = 0.0;
    int step = 0;
    double transport_runtime = 10.0;
    double mpi_time = 5.0;
    write_silo(mesh, time, step, transport_runtime, mpi_time, rank, n_rank);

    if (silo_write_pass)
      cout << "TEST PASSED: writing simple mesh silo file" << endl;
    else {
      cout << "TEST FAILED:  writing silo file" << endl;
      nfail++;
    }
  }

  // Test that a simple mesh (one division in each dimension) is constructed
  // correctly from the input file (simple_input.xml) and that each cell
  // is assigned the correct region
  {
    // put MPI_Types in scope so it's destructor is called before MPI_Finalize
    MPI_Types mpi_types;
    string filename("three_region_mesh_input.xml");
    Input input(filename, mpi_types);
    IMC_Parameters imc_p(input);

    Mesh mesh(input, mpi_types, mpi_info, imc_p);

    bool three_reg_silo_write_pass = true;

    double time = 2.0;
    int step = 1;
    double transport_runtime = 7.0;
    double mpi_time = 2.0;
    write_silo(mesh, time, step, transport_runtime, mpi_time, rank, n_rank);

    if (three_reg_silo_write_pass) {
      cout << "TEST PASSED: writing three region mesh silo file" << endl;
    } else {
      cout << "TEST FAILED:  writing silo file" << endl;
      nfail++;
    }
  }

  MPI_Finalize();

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_write_silo.cc
//---------------------------------------------------------------------------//
