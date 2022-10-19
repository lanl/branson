//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   test_sampling_functions.cc
 * \author Alex Long
 * \date   January 14 2016
 * \brief  Test construction, fill and get functions
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#include <iostream>
#include <mpi.h>
#include <vector>

#include "testing_functions.h"
#include "../sampling_functions.h"
#include "../RNG.h"
#include "../cell.h"

int main(void) {

  using std::cout;
  using std::endl;
  using std::vector;

  int nfail = 0;

  // test position sampling functions
  {
    bool test_sampling_functions_position = true;
    RNG rng;
    rng.set_seed(72412);

    //setup cell
    Cell cell;
    // simple cube of size 1.0
    double x_low = 0.0;
    double x_high = 1.0;
    double y_low = 0.0;
    double y_high = 1.0;
    double z_low = 0.0;
    double z_high = 1.0;
    cell.set_coor(x_low, x_high, y_low, y_high, z_low, z_high);

    // sample the position to makes it matches
    std::array<double, 3> avg_pos{0.0,0.0,0.0};
    constexpr int n_samples = 10000;

    for (int i = 0;i<n_samples;++i) {
      auto uniform_pos = get_uniform_position_in_cell(cell, &rng);
      avg_pos[0] += uniform_pos[0];
      avg_pos[1] += uniform_pos[1];
      avg_pos[2] += uniform_pos[2];
    }
    avg_pos[0] = avg_pos[0]/static_cast<double>(n_samples);
    avg_pos[1] = avg_pos[1]/static_cast<double>(n_samples);
    avg_pos[2] = avg_pos[2]/static_cast<double>(n_samples);

    std::cout<<"get_uniform_position_in_cell average position with unit cell: ";
    std::cout<<avg_pos[0]<<" "<<avg_pos[1]<<" "<<avg_pos[2]<<std::endl;

    constexpr double tolerance = 3.0e-3;
    if(!soft_equiv(avg_pos[0], 0.5, tolerance))
      test_sampling_functions_position = false;
    if(!soft_equiv(avg_pos[1], 0.5, tolerance))
      test_sampling_functions_position = false;
    if(!soft_equiv(avg_pos[2], 0.5, tolerance))
      test_sampling_functions_position = false;

    if (test_sampling_functions_position)
      cout << "TEST PASSED: Sampling Functions--position " << endl;
    else {
      cout << "TEST FAILED: Sampling Functions--position" << endl;
      nfail++;
    }
  }

  // test angle sampling functions
  {
    bool test_sampling_functions_angle = true;
    RNG rng;
    rng.set_seed(72412);

    // sample the angle to makes it matches
    std::array<double, 3> avg_angle{0.0,0.0,0.0};
    constexpr int n_samples = 10000;

    for (int i = 0;i<n_samples;++i) {
      auto uniform_angle = get_uniform_angle(&rng);
      avg_angle[0] += uniform_angle[0];
      avg_angle[1] += uniform_angle[1];
      avg_angle[2] += uniform_angle[2];
    }
    avg_angle[0] = avg_angle[0]/static_cast<double>(n_samples);
    avg_angle[1] = avg_angle[1]/static_cast<double>(n_samples);
    avg_angle[2] = avg_angle[2]/static_cast<double>(n_samples);

    std::cout<<"get_uniform_angle average angle: ";
    std::cout<<avg_angle[0]<<" "<<avg_angle[1]<<" "<<avg_angle[2]<<std::endl;

    constexpr double tolerance = 2.0e-2;
    if(!soft_equiv(avg_angle[0], 0.0, tolerance))
      test_sampling_functions_angle = false;
    if(!soft_equiv(avg_angle[1], 0.0, tolerance))
      test_sampling_functions_angle = false;
    if(!soft_equiv(avg_angle[2], 0.0, tolerance))
      test_sampling_functions_angle = false;

    if (test_sampling_functions_angle)
      cout << "TEST PASSED: Sampling Functions--angle " << endl;
    else {
      cout << "TEST FAILED: Sampling Functions--angle" << endl;
      nfail++;
    }
  }

  // test group sampling function (reminder, group is just used to give more random access patterns
  // to otherwise gray branson)
  {
    bool test_sampling_functions_group = true;
    RNG rng;
    rng.set_seed(72412);

    //setup cell
    Cell cell;
    // simple cube of size 1.0
    double x_low = 0.0;
    double x_high = 1.0;
    double y_low = 0.0;
    double y_high = 1.0;
    double z_low = 0.0;
    double z_high = 1.0;
    cell.set_coor(x_low, x_high, y_low, y_high, z_low, z_high);
    cell.set_op_a(7.0);
    cell.set_op_s(0.12);

    // sample the group to makes it matches
    double avg_group{0.0};
    constexpr int n_samples = 10000;

    for (int i = 0;i<n_samples;++i) {
      auto uniform_group = sample_emission_group(&rng, cell);
      avg_group += uniform_group;
    }
    avg_group = avg_group/static_cast<double>(n_samples);

    std::cout<<"sample_emission_group average group: "<<avg_group<<std::endl;

    constexpr double tolerance = 1.0e-1;
    constexpr double half_group = static_cast<double>(BRANSON_N_GROUPS-1)/2.0;
    std::cout<<"half of (n_groups-1):" <<half_group<<std::endl;
    if(!soft_equiv(avg_group, half_group, tolerance))
      test_sampling_functions_group = false;

    if (test_sampling_functions_group)
      cout << "TEST PASSED: Sampling Functions--group " << endl;
    else {
      cout << "TEST FAILED: Sampling Functions--group" << endl;
      nfail++;
    }
  }


  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_sampling_functions.cc
//---------------------------------------------------------------------------//

