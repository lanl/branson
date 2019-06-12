//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_proto_cell.cc
 * \author Alex Long
 * \date   October 29 2018
 * \brief  Test proto cell construction, get and set functions
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include "../constants.h"
#include "../proto_cell.h"
#include "testing_functions.h"
#include <iostream>
#include <vector>

int main(void) {

  using Constants::bc_type;
  using Constants::dir_type;
  using std::cout;
  using std::endl;
  using std::vector;
  int nfail = 0;

  // get/set functions
  {
    bool get_set_functions_pass = true;

    const vector<dir_type> dirs{Constants::X_NEG, Constants::X_POS,
                                Constants::Y_NEG, Constants::Y_POS,
                                Constants::Z_NEG, Constants::Z_POS};

    //setup cell
    Proto_Cell proto_cell;

    uint32_t cell_ID = 3271733928; // 64-bit cell ID
    uint32_t region_ID = 12;
    uint32_t silo_index = 1231;

    vector<Constants::bc_type> bcs{Constants::REFLECT, Constants::VACUUM,
                                   Constants::ELEMENT, Constants::PROCESSOR,
                                   Constants::ELEMENT, Constants::VACUUM};
    vector<uint32_t> neighbors{3500000000, 3500000001, 3500000002,
                               3500000003, 3500000004, 3500000005};
    vector<uint32_t> grip_neighbors{2500000000, 2500000001, 2500000002,
                                    2500000003, 2500000004, 2500000005};

    // simple cube of size 1.0
    double x_low = 0.0;
    double x_high = 1.0;
    double y_low = 0.0;
    double y_high = 1.0;
    double z_low = 0.0;
    double z_high = 1.0;
    vector<double> coords = {x_low, x_high, y_low, y_high, z_low, z_high};

    // set values
    proto_cell.set_coor(x_low, x_high, y_low, y_high, z_low, z_high);
    proto_cell.set_ID(cell_ID);
    proto_cell.set_grip_ID(cell_ID);
    proto_cell.set_region_ID(region_ID);
    proto_cell.set_silo_index(silo_index);
    for (auto i : dirs) {
      proto_cell.set_bc(i, bcs[i]);
      proto_cell.set_neighbor(i, neighbors[i]);
      proto_cell.set_grip_neighbor(i, grip_neighbors[i]);
    }

    // test the get methods
    const double *cell_coords = proto_cell.get_node_array();
    if (proto_cell.get_ID() != cell_ID)
      get_set_functions_pass = false;
    if (proto_cell.get_grip_ID() != cell_ID)
      get_set_functions_pass = false;
    if (proto_cell.get_region_ID() != region_ID)
      get_set_functions_pass = false;
    if (proto_cell.get_silo_index() != silo_index)
      get_set_functions_pass = false;
    for (int i = 0; i < 6; ++i) {
      if (proto_cell.get_next_cell(i) != neighbors[i])
        get_set_functions_pass = false;
      if (proto_cell.get_next_grip(i) != grip_neighbors[i])
        get_set_functions_pass = false;
      if (proto_cell.get_bc(i) != bcs[i])
        get_set_functions_pass = false;
      if (cell_coords[i] != coords[i])
        get_set_functions_pass = false;
    }

    if (get_set_functions_pass)
      cout << "TEST PASSED: Proto_Cell get/set functions" << endl;
    else {
      cout << "TEST FAILED: Proto_Cell get/set functions" << endl;
      nfail++;
    }
  }

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_proto_cell.cc
//---------------------------------------------------------------------------//
