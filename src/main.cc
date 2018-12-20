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
#include "info.h"
#include "input.h"
#include "mesh.h"
#include "mesh_pass_driver.h"
#include "mpi_types.h"
#include "particle_pass_driver.h"
#include "replicated_driver.h"
#include "rma_mesh_pass_driver.h"
#include "timer.h"

using Constants::CELL_PASS;
using Constants::CELL_PASS_RMA;
using Constants::PARTICLE_PASS;
using Constants::REPLICATED;
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  
  // check to see if number of arguments is correct
  if (argc != 2) {
    cout << "Usage: BRANSON <path_to_input_file>" << endl;
    exit(EXIT_FAILURE);
  }

  // wrap main loop scope so objcts are destroyed before mpi_finalize is called
  {
    // get MPI parmeters and set them in mpi_info
    const Info mpi_info;
    if (mpi_info.get_rank() == 0) {
      cout << "-------- Branson, a massively parallel proxy app for Implicit Monte Carlo ------"<<endl;
      cout << "-------- Author: Alex Long (along@lanl.gov) ------------------------------------"<<endl;
      cout << "-------- Version: 0.8 ----------------------------------------------------------"<<endl<<endl;
      cout << " Branson compiled on: " << mpi_info.get_machine_name() << endl;
    }

    // make MPI types object
    MPI_Types mpi_types;

    // get input object from filename
    std::string filename(argv[1]);
    Input input(filename, mpi_types);
    if (mpi_info.get_rank() == 0)
      input.print_problem_info();

    // IMC paramters setup
    IMC_Parameters imc_p(input);

    // IMC state setup
    IMC_State imc_state(input, mpi_info.get_rank());

    // timing
    Timer timers;

    // make mesh from input object
    timers.start_timer("Total setup");

    Mesh mesh(input, mpi_types, mpi_info, imc_p);
    mesh.initialize_physical_properties(input);

    timers.stop_timer("Total setup");

    MPI_Barrier(MPI_COMM_WORLD);
    // print_MPI_out(mesh, rank, n_rank);

    //--------------------------------------------------------------------------//
    // TRT PHYSICS CALCULATION
    //--------------------------------------------------------------------------//

    timers.start_timer("Total transport");

    if (input.get_dd_mode() == PARTICLE_PASS)
      imc_particle_pass_driver(mesh, imc_state, imc_p, mpi_types, mpi_info);
    else if (input.get_dd_mode() == CELL_PASS)
      imc_mesh_pass_driver(mesh, imc_state, imc_p, mpi_types, mpi_info);
    else if (input.get_dd_mode() == CELL_PASS_RMA)
      imc_rma_mesh_pass_driver(mesh, imc_state, imc_p, mpi_types, mpi_info);
    else if (input.get_dd_mode() == REPLICATED)
      imc_replicated_driver(mesh, imc_state, imc_p, mpi_types, mpi_info);
    else {
      cout << "Driver for DD transport method currently not supported" << endl;
      exit(EXIT_FAILURE);
    }

    timers.stop_timer("Total transport");

    if (mpi_info.get_rank() == 0) {
      cout << "****************************************";
      cout << "****************************************" << endl;
      imc_state.print_simulation_footer(input.get_dd_mode());
      timers.print_timers();
    }

  } // end main loop scope, objects destroyed here

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}
//---------------------------------------------------------------------------//
// end of main.cc
//---------------------------------------------------------------------------//
