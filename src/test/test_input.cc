//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_input.cc
 * \author Alex Long
 * \date   February 11 2016
 * \brief  Test Input class for correct reading of XML input files
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>

#include "../constants.h"
#include "../input.h"
#include "../region.h"
#include "testing_functions.h"

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  using Constants::REPLICATED;
  using Constants::REFLECT;
  using Constants::VACUUM;
  using Constants::X_NEG;
  using Constants::X_POS;
  using Constants::Y_NEG;
  using Constants::Y_POS;
  using Constants::Z_NEG;
  using Constants::Z_POS;
  using std::cout;
  using std::endl;
  using std::string;

  int nfail = 0;

  // scope for MPI_Types
  {
    MPI_Types mpi_types;

    // test the get functions to make sure correct values are set from the input
    // file (reader is working correctly) and that the get functions are working
    // these values are hardcoded in simple_input.xml
    {
      // test simple input file (one division in each dimension and one region)
      string filename("simple_input.xml");
      Input input(filename, mpi_types);

      bool simple_input_pass = true;
      if (input.get_global_n_x_cells() != 10)
        simple_input_pass = false;
      if (input.get_global_n_y_cells() != 20)
        simple_input_pass = false;
      if (input.get_global_n_z_cells() != 30)
        simple_input_pass = false;
      if (input.get_dx(0) != 1.0)
        simple_input_pass = false;
      if (input.get_dy(0) != 2.0)
        simple_input_pass = false;
      if (input.get_dz(0) != 3.0)
        simple_input_pass = false;

      if (input.get_comb_bool() != false)
        simple_input_pass = false;
      if (input.get_verbose_print_bool() != true)
        simple_input_pass = false;
      if (input.get_print_mesh_info_bool() != true)
        simple_input_pass = false;
      if (input.get_output_freq() != 1)
        simple_input_pass = false;

      if (input.get_dt() != 0.01)
        simple_input_pass = false;
      if (input.get_time_start() != 0.0)
        simple_input_pass = false;
      if (input.get_time_finish() != 0.1)
        simple_input_pass = false;
      if (input.get_time_mult() != 1.0)
        simple_input_pass = false;
      if (input.get_rng_seed() != 14706)
        simple_input_pass = false;
      if (input.get_number_photons() != 10000)
        simple_input_pass = false;
      // Even though input says 10k, this defaults to the replicated
      // value of 100M when you're running 1 MPI rank.
      if (input.get_batch_size() != 100000000)
        simple_input_pass = false;
      if (input.get_particle_message_size() != 1000)
        simple_input_pass = false;
      // Even though input says PARTICLE_PASS, this defaults back to replicated
      // when you're only running 1 MPI rank.
      if (input.get_dd_mode() != REPLICATED)
        simple_input_pass = false;
      if(input.get_use_gpu_transporter_bool() != false)
        simple_input_pass = false;

      if (input.get_bc(X_NEG) != REFLECT)
        simple_input_pass = false;
      if (input.get_bc(X_POS) != REFLECT)
        simple_input_pass = false;
      if (input.get_bc(Y_NEG) != VACUUM)
        simple_input_pass = false;
      if (input.get_bc(Y_POS) != VACUUM)
        simple_input_pass = false;
      if (input.get_bc(Z_NEG) != VACUUM)
        simple_input_pass = false;
      if (input.get_bc(Z_POS) != REFLECT)
        simple_input_pass = false;

      // test region functionality
      uint32_t region_index = input.get_region_index(0, 0, 0);
      std::vector<Region> regions = input.get_regions();
      Region region = regions[region_index];

      if (input.get_n_regions() != 1)
        simple_input_pass = false;
      if (region_index != 0)
        simple_input_pass = false;

      if (region.get_ID() != 6)
        simple_input_pass = false;
      if (region.get_cV() != 2.0)
        simple_input_pass = false;
      if (region.get_rho() != 1.0)
        simple_input_pass = false;
      if (region.get_opac_A() != 3.0)
        simple_input_pass = false;
      if (region.get_opac_B() != 1.5)
        simple_input_pass = false;
      if (region.get_opac_C() != 0.1)
        simple_input_pass = false;
      if (region.get_scattering_opacity() != 5.0)
        simple_input_pass = false;
      if (region.get_T_e() != 1.0)
        simple_input_pass = false;
      if (region.get_T_r() != 1.1)
        simple_input_pass = false;
      if (region.get_T_s() != 0.0)
        simple_input_pass = false;

      if (simple_input_pass)
        cout << "TEST PASSED: simple input get functions" << endl;
      else {
        cout << "TEST FAILED: simple input get functions" << endl;
        nfail++;
      }
    }

    // test the get functions to make sure correct values are set from the input
    // file with a more complicated mesh (reader is working correctly) these
    // values are hardcoded in three_region_mesh_input.xml
    {
      // test simple input file (one division in each dimension and one region)
      string filename("three_region_mesh_input.xml");
      Input input(filename, mpi_types);

      bool three_region_pass = true;
      if (input.get_global_n_x_cells() != 21)
        three_region_pass = false;
      if (input.get_global_n_y_cells() != 10)
        three_region_pass = false;
      if (input.get_global_n_z_cells() != 1)
        three_region_pass = false;

      if (input.get_n_x_divisions() != 3)
        three_region_pass = false;
      if (input.get_dx(0) != 1.0)
        three_region_pass = false;
      if (input.get_x_division_cells(0) != 4)
        three_region_pass = false;
      if (input.get_dx(1) != 2.0)
        three_region_pass = false;
      if (input.get_x_division_cells(1) != 2)
        three_region_pass = false;
      if (input.get_dx(2) != 2.0 / 15.0)
        three_region_pass = false;
      if (input.get_x_division_cells(2) != 15)
        three_region_pass = false;

      if (input.get_n_y_divisions() != 1)
        three_region_pass = false;
      if (input.get_dy(0) != 3.0)
        three_region_pass = false;
      if (input.get_y_division_cells(0) != 10)
        three_region_pass = false;

      if (input.get_n_z_divisions() != 1)
        three_region_pass = false;
      if (input.get_dz(0) != 1.0)
        three_region_pass = false;
      if (input.get_z_division_cells(0) != 1)
        three_region_pass = false;

      if (input.get_comb_bool() != true)
        three_region_pass = false;
      if (input.get_verbose_print_bool() != false)
        three_region_pass = false;
      if (input.get_print_mesh_info_bool() != false)
        three_region_pass = false;
      if (input.get_output_freq() != 1)
        three_region_pass = false;

      if (input.get_dt() != 0.01)
        three_region_pass = false;
      if (input.get_time_start() != 0.0)
        three_region_pass = false;
      if (input.get_time_finish() != 0.1)
        three_region_pass = false;
      if (input.get_time_mult() != 1.0)
        three_region_pass = false;
      if (input.get_time_mult() != 1.0)
        three_region_pass = false;
      if (input.get_rng_seed() != 14706)
        three_region_pass = false;
      if (input.get_number_photons() != 10000)
        three_region_pass = false;
      // Even though input says 10k, this defaults to the replicated
      // value of 100M when you're running 1 MPI rank.
      if (input.get_batch_size() != 100000000)
        three_region_pass = false;
      if (input.get_particle_message_size() != 1000)
        three_region_pass = false;
      // Even though input says PARTICLE_PASS, this defaults back to replicated
      // when you're only running 1 MPI rank.
      if (input.get_dd_mode() != REPLICATED)
        three_region_pass = false;

      if (input.get_bc(X_NEG) != REFLECT)
        three_region_pass = false;
      if (input.get_bc(X_POS) != REFLECT)
        three_region_pass = false;
      if (input.get_bc(Y_NEG) != VACUUM)
        three_region_pass = false;
      if (input.get_bc(Y_POS) != VACUUM)
        three_region_pass = false;
      if (input.get_bc(Z_NEG) != VACUUM)
        three_region_pass = false;
      if (input.get_bc(Z_POS) != REFLECT)
        three_region_pass = false;

      // test region functionality
      uint32_t region_index;
      Region region;
      std::vector<Region> regions = input.get_regions();

      if (input.get_n_regions() != 3)
        three_region_pass = false;

      region_index = input.get_region_index(0, 0, 0);
      if (region_index != 0)
        three_region_pass = false;
      region_index = input.get_region_index(1, 0, 0);
      if (region_index != 1)
        three_region_pass = false;
      region_index = input.get_region_index(2, 0, 0);
      if (region_index != 2)
        three_region_pass = false;

      // region 230
      region = regions[0];
      if (region.get_ID() != 230)
        three_region_pass = false;
      if (region.get_cV() != 2.0)
        three_region_pass = false;
      if (region.get_rho() != 1.0)
        three_region_pass = false;
      if (region.get_opac_A() != 3.0)
        three_region_pass = false;
      if (region.get_opac_B() != 1.5)
        three_region_pass = false;
      if (region.get_opac_C() != 0.1)
        three_region_pass = false;
      if (region.get_scattering_opacity() != 5.0)
        three_region_pass = false;
      if (region.get_T_e() != 1.0)
        three_region_pass = false;
      if (region.get_T_r() != 1.1)
        three_region_pass = false;
      if (region.get_T_s() != 0.0)
        three_region_pass = false;

      // region 177
      region = regions[1];
      if (region.get_ID() != 177)
        three_region_pass = false;
      if (region.get_cV() != 0.99)
        three_region_pass = false;
      if (region.get_rho() != 5.0)
        three_region_pass = false;
      if (region.get_opac_A() != 101.0)
        three_region_pass = false;
      if (region.get_opac_B() != 10.5)
        three_region_pass = false;
      if (region.get_opac_C() != 0.3)
        three_region_pass = false;
      if (region.get_scattering_opacity() != 0.01)
        three_region_pass = false;
      if (region.get_T_e() != 0.01)
        three_region_pass = false;
      if (region.get_T_r() != 0.1)
        three_region_pass = false;
      if (region.get_T_s() != 0.0)
        three_region_pass = false;

      // region 11
      region = regions[2];
      if (region.get_ID() != 11)
        three_region_pass = false;
      if (region.get_cV() != 5.0)
        three_region_pass = false;
      if (region.get_rho() != 100.0)
        three_region_pass = false;
      if (region.get_opac_A() != 0.001)
        three_region_pass = false;
      if (region.get_opac_B() != 0.01)
        three_region_pass = false;
      if (region.get_opac_C() != 4.8)
        three_region_pass = false;
      if (region.get_scattering_opacity() != 100.0)
        three_region_pass = false;
      if (region.get_T_e() != 1.2)
        three_region_pass = false;
      if (region.get_T_r() != 0.0)
        three_region_pass = false;
      if (region.get_T_s() != 0.0)
        three_region_pass = false;

      if (three_region_pass)
        cout << "TEST PASSED: three region input" << endl;
      else {
        cout << "TEST FAILED: three region input" << endl;
        nfail++;
      }
    }

    // test assigning a larger number than uint32_t to the number of photons and
    // make sure it's recognized
    {
      bool large_particle_pass = true;
      // first test large particle input file
      std::string large_filename("large_particle_input.xml");
      Input large_input(large_filename, mpi_types);

      if (large_input.get_number_photons() != 6000000000)
        large_particle_pass = false;

      cout << "particle count = " << large_input.get_number_photons() << endl;
      if (large_particle_pass)
        cout << "TEST PASSED: 64 bit particle count" << endl;
      else {
        cout << "TEST FAILED: 64 bit particle count" << endl;
        nfail++;
      }
    }
  } // end scope of MPI Types

  MPI_Finalize();

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_input.cc
//---------------------------------------------------------------------------//
