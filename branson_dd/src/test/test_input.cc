/*
  Author: Alex Long
  Date: 2/11/2016
  Name: test_input.cc
*/

#include <iostream>
#include <string>
#include "../input.h"
#include "../constants.h"
#include "testing_functions.h"

int main (void) {

  using std::cout;
  using std::endl;
  using std::string;
  using Constants::PARTICLE_PASS;
  using Constants::CELL_PASS;
  using Constants::VACUUM; 
  using Constants::REFLECT;
  using Constants::X_NEG;
  using Constants::X_POS;
  using Constants::Y_NEG;
  using Constants::Y_POS;
  using Constants::Z_NEG;
  using Constants::Z_POS; 

  int nfail = 0;

  // test the get functions to make sure correct values are set from the input
  // file (reader is working correctly) and that the get functions are working 
  // these values are hardcoded in simple_input.xml
  {
    // test simple input file
    string filename("simple_input.xml");
    Input input(filename);

    bool get_functions_pass = true;
    if (input.get_n_x_cells() != 10) get_functions_pass =false;
    if (input.get_n_y_cells() != 20) get_functions_pass =false;
    if (input.get_n_z_cells() != 30) get_functions_pass =false;
    if (input.get_dx() != 1.0) get_functions_pass =false;
    if (input.get_dy() != 2.0) get_functions_pass =false;
    if (input.get_dz() != 3.0) get_functions_pass =false;

    if (input.get_initial_Tm() != 1.0) get_functions_pass =false;
    if (input.get_initial_Tr() != 1.1) get_functions_pass =false;


    if (input.get_tilt_bool() != false) get_functions_pass =false;
    if (input.get_comb_bool() != true) get_functions_pass =false;
    if (input.get_stratified_bool() != false) get_functions_pass =false;
    if (input.get_ghost_cell_bool() != false) get_functions_pass =false;
    if (input.get_verbose_print_bool() != false) get_functions_pass =false;
    if (input.get_print_mesh_info_bool() != false) get_functions_pass =false;
    if (input.get_output_freq() != 1) get_functions_pass =false;

    if (input.get_dt() != 0.01) get_functions_pass =false;
    if (input.get_time_start() != 0.0 ) get_functions_pass =false;
    if (input.get_time_finish() != 0.1 ) get_functions_pass =false;
    if (input.get_time_mult() != 1.0 ) get_functions_pass =false;
    if (input.get_time_mult() != 1.0 ) get_functions_pass =false;
    if (input.get_rng_seed() != 14706) get_functions_pass =false;
    if (input.get_number_photons() != 10000) get_functions_pass =false;
    if (input.get_batch_size() != 10000) get_functions_pass=false; 
    if (input.get_particle_message_size() != 1000) get_functions_pass=false; 
    if (input.get_map_size() != 50000) get_functions_pass=false; 
    if (input.get_dd_mode() != PARTICLE_PASS) get_functions_pass=false; 

    if (input.get_source_T() != 0.0) get_functions_pass=false; 
    if (input.get_CV() != 2.0) get_functions_pass=false; 
    if (input.get_rho() != 1.0) get_functions_pass=false; 
    if (input.get_opacity_A() != 3.0) get_functions_pass=false; 
    if (input.get_opacity_B() != 0.0) get_functions_pass=false; 
    if (input.get_opacity_C() != 0.0) get_functions_pass=false; 
    if (input.get_opacity_S() != 5.0) get_functions_pass=false; 

    if (input.get_bc(X_NEG) != REFLECT) get_functions_pass=false;
    if (input.get_bc(X_POS) != REFLECT) get_functions_pass=false;
    if (input.get_bc(Y_NEG) != VACUUM) get_functions_pass=false;
    if (input.get_bc(Y_POS) != VACUUM) get_functions_pass=false;
    if (input.get_bc(Z_NEG) != VACUUM) get_functions_pass=false;
    if (input.get_bc(Z_POS) != REFLECT) get_functions_pass=false;

    if (get_functions_pass) cout<<"TEST PASSED: get functions"<<endl;
    else { 
      cout<<"TEST FAILED: get functions"<<endl; 
      nfail++;
    }
  }

  // test assigning a larger number than uint32_t to the number of photons and
  // make sure it's recognized
  {
    bool large_particle_pass = true;
    // first test large particle input file
    std::string large_filename("large_particle_input.xml");
    Input large_input(large_filename);
    
    if (large_input.get_number_photons() != 6000000000) large_particle_pass =false;

    cout<<"particle count = "<<large_input.get_number_photons()<<endl;
    if (large_particle_pass) cout<<"TEST PASSED: 64 bit particle count"<<endl;
    else { 
      cout<<"TEST FAILED: 64 bit particle count"<<endl; 
      nfail++;
    }
  }

  return nfail;
}
