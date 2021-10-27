//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   main.cc
 * \author Alex Long
 * \date   July 24 2014
 * \brief  Reads input file, sets up mesh and runs transport
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <mpi.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <vector>

#include "constants.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "input.h"
#include "mesh.h"
#include "replicated_driver.h"
#include "timer.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int n_rank = 0;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  if(n_rank > 1) {
    cout << "Cycle counting can only be run with one rank!"<<endl;
    exit(EXIT_FAILURE);
  }

  // check to see if number of arguments is correct
  if (argc != 2) {
    cout << "Usage: BRANSON <path_to_input_file>" << endl;
    exit(EXIT_FAILURE);
  }


  // wrap main loop scope so objcts are destroyed before mpi_finalize is called
  {
    // get MPI parmeters and set them in mpi_info
    cout << "----- Branson LITE (with cycle counting), a massively parallel proxy app for Implicit "
            "Monte Carlo ----"
         << endl;
    cout << "-------- Author: Alex Long (along@lanl.gov) "
            "------------------------------------"
         << endl;
    cout << "-------- Version: 0.81"
            "----------------------------------------------------------"
         << endl
         << endl;

    // get input object from filename
    std::string filename(argv[1]);
    Input input(filename);
    input.print_problem_info();

    // IMC paramters setup
    IMC_Parameters imc_p(input);

    // IMC state setup
    IMC_State imc_state(input);

    // timing
    Timer timers;

    // make mesh from input object
    timers.start_timer("Total setup");

    Mesh mesh(input, imc_p);
    mesh.initialize_physical_properties(input);

    timers.stop_timer("Total setup");

    MPI_Barrier(MPI_COMM_WORLD);

    //--------------------------------------------------------------------------//
    // TRT PHYSICS CALCULATION
    //--------------------------------------------------------------------------//

    timers.start_timer("Total transport");

    imc_replicated_driver(mesh, imc_state, imc_p);

    timers.stop_timer("Total transport");

    cout << "****************************************";
    cout << "****************************************" << endl;
    timers.print_timers();

  } // end main loop scope, objects destroyed here

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}
//---------------------------------------------------------------------------//
// end of main.cc
//---------------------------------------------------------------------------//
