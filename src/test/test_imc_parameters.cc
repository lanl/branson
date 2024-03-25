//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_imc_parameters.cc
 * \author Alex Long
 * \date   February 11 2016
 * \brief  Test Input class for correct reading of XML imc_parameters files
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>

#include "../constants.h"
#include "../input.h"
#include "../imc_parameters.h"
#include "../region.h"
#include "testing_functions.h"

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  using std::cout;
  using std::endl;
  using std::string;

  int nfail = 0;

  // scope for MPI_Types
  {
    MPI_Types mpi_types;

    // test the get functions to make sure correct values are set from the imc_parameters class,
    // right now the imc_parameters are just a subset of the input file se we construct an input
    // file first
    {
      using Constants::REPLICATED;
      // test simple imc_parameters file (one division in each dimension and one region)
      string filename("simple_input.xml");
      Input input(filename, mpi_types);
      IMC_Parameters imc_parameters(input);
      bool simple_imc_parameters_pass = true;

      if (imc_parameters.get_use_comb_flag() != false)
        simple_imc_parameters_pass = false;
      if (imc_parameters.get_output_frequency() != 1)
        simple_imc_parameters_pass = false;
      if (imc_parameters.get_rng_seed() != 14706)
        simple_imc_parameters_pass = false;
      if (imc_parameters.get_n_user_photons() != 10000)
        simple_imc_parameters_pass = false;
      // Even though input file says 10k, this defaults to the replicated
      // value of 100M when you're running 1 MPI rank.
      if (imc_parameters.get_batch_size() != 100000000)
        simple_imc_parameters_pass = false;
      if (imc_parameters.get_particle_message_size() != 1000)
        simple_imc_parameters_pass = false;
      // Even though imc_parameters says PARTICLE_PASS, this defaults back to replicated
      // when you're only running 1 MPI rank.
      if (imc_parameters.get_dd_mode() != REPLICATED)
        simple_imc_parameters_pass = false;
      if(imc_parameters.get_use_gpu_transporter_flag() != false)
        simple_imc_parameters_pass = false;

      if (simple_imc_parameters_pass)
        cout << "TEST PASSED: simple IMC_Parameters get functions" << endl;
      else {
        cout << "TEST FAILED: simple IMC_Parameters get functions" << endl;
        nfail++;
      }
    }

  } // end scope of MPI Types

  MPI_Finalize();

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_imc_parameters.cc
//---------------------------------------------------------------------------//

