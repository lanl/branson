//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_cell.cc
 * \author Alex Long
 * \date   January 11 2016
 * \brief  Test cell check_in_cell and distance_to_boundary functions
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include "../cell.h"
#include "../constants.h"
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

  // test the check_in_cell function
  {
    double true_pos_1[3] = {0.5, 0.5, 0.5};
    double true_pos_2[3] = {0.9, 0.9, 0.9};
    double true_pos_3[3] = {0.1, 0.1, 0.1};
    double true_pos_4[3] = {0.1, 0.9, 0.1};
    double false_pos_1[3] = {-1.0, -1.0, -1.0};
    double false_pos_2[3] = {1.0, -1.0, -1.0};
    double false_pos_3[3] = {1.0, 1.0, -1.0};
    double false_pos_4[3] = {-1.0, 1.0, -1.0};

    bool check_in_cell_pass = true;
    //positions in cell
    if (!cell.check_in_cell(true_pos_1))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_2))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_3))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_4))
      check_in_cell_pass = false;
    //positions out of cell
    if (cell.check_in_cell(false_pos_1))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_2))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_3))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_4))
      check_in_cell_pass = false;

    if (check_in_cell_pass)
      cout << "TEST PASSED: check_in_cell function" << endl;
    else {
      cout << "TEST FAILED: check_in_cell function" << endl;
      nfail++;
    }
  }

  // test distance to boundary function
  {
    double tolerance = 1.0e-8;

    double pos[3] = {0.5, 0.5, 0.5};

    double angle_1[3] = {0.999, 0.031614, 0.031614};
    double angle_2[3] = {0.031614, 0.999, 0.031614};
    double angle_3[3] = {0.031614, 0.31614, 0.999};

    double angle_4[3] = {-0.999, 0.031614, 0.031614};
    double angle_5[3] = {0.031614, -0.999, 0.031614};
    double angle_6[3] = {0.031614, 0.31614, -0.999};

    unsigned int surface_cross;
    bool distance_to_bound_pass = true;

    // tests for true
    double distance_1 =
        cell.get_distance_to_boundary(pos, angle_1, surface_cross);
    if (!soft_equiv(distance_1, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_2 =
        cell.get_distance_to_boundary(pos, angle_2, surface_cross);
    if (!soft_equiv(distance_2, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_3 =
        cell.get_distance_to_boundary(pos, angle_3, surface_cross);
    if (!soft_equiv(distance_3, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_4 =
        cell.get_distance_to_boundary(pos, angle_4, surface_cross);
    if (!soft_equiv(distance_4, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_5 =
        cell.get_distance_to_boundary(pos, angle_5, surface_cross);
    if (!soft_equiv(distance_5, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_6 =
        cell.get_distance_to_boundary(pos, angle_6, surface_cross);
    if (!soft_equiv(distance_6, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    if (distance_to_bound_pass)
      cout << "TEST PASSED: distance_to_boundary function" << endl;
    else {
      cout << "TEST FAILED: distance_to_boundary function" << endl;
      nfail++;
    }
  }

  // test get_volume function
  {
    bool get_volume_pass = true;

    double tolerance = 1.0e-8;

    if (!soft_equiv(cell.get_volume(), 1.0, tolerance))
      get_volume_pass = false;

    //make an oblong cell
    Cell oblong_cell;
    oblong_cell.set_coor(0.01, 0.02, 0.0, 10.0, -0.1, 0.1);

    if (!soft_equiv(oblong_cell.get_volume(), 0.02, tolerance))
      get_volume_pass = false;

    if (get_volume_pass)
      cout << "TEST PASSED: get_volume function" << endl;
    else {
      cout << "TEST FAILED: get_volume function" << endl;
      nfail++;
    }
  }

  // test construction from Proto_Cell
  {
    bool proto_construction_pass = true;

    const vector<dir_type> dirs{Constants::X_NEG, Constants::X_POS,
                                Constants::Y_NEG, Constants::Y_POS,
                                Constants::Z_NEG, Constants::Z_POS};

    //setup cell
    Proto_Cell proto_cell;

    uint32_t cell_ID = 3271733928; // 64-bit cell ID
    uint32_t grip_ID = 3271733920; // 64-bit cell ID
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
    proto_cell.set_grip_ID(grip_ID);
    proto_cell.set_region_ID(region_ID);
    proto_cell.set_silo_index(silo_index);
    for (auto i : dirs) {
      proto_cell.set_bc(i, bcs[i]);
      proto_cell.set_neighbor(i, neighbors[i]);
      proto_cell.set_grip_neighbor(i, grip_neighbors[i]);
    }

    Cell from_proto = Cell(proto_cell);

    // test the get methods
    const double *cell_coords = from_proto.get_node_array();
    if (from_proto.get_ID() != cell_ID)
      proto_construction_pass = false;
    if (from_proto.get_grip_ID() != grip_ID)
      proto_construction_pass = false;
    if (from_proto.get_region_ID() != region_ID)
      proto_construction_pass = false;
    if (from_proto.get_silo_index() != silo_index)
      proto_construction_pass = false;
    for (int i = 0; i < 6; ++i) {
      if (from_proto.get_next_cell(i) != neighbors[i])
        proto_construction_pass = false;
      if (from_proto.get_next_grip(i) != grip_neighbors[i])
        proto_construction_pass = false;
      if (from_proto.get_bc(i) != bcs[i])
        proto_construction_pass = false;
      if (cell_coords[i] != coords[i])
        proto_construction_pass = false;
    }

    if (proto_construction_pass)
      cout << "TEST PASSED: construction from Proto_Cell" << endl;
    else {
      cout << "TEST FAILED: construction from proto cell" << endl;
      nfail++;
    }
  }

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_cell.cc
//---------------------------------------------------------------------------//
