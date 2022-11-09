//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_photon.cc
 * \author Alex Long
 * \date   January 11 2016
 * \brief  Test photon construction and move functionality
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "../photon.h"
#include "../mpi_types.h"
#include "testing_functions.h"

using std::cout;
using std::endl;

int simple_photon_tests() {
  int n_fail = 0;

  // test construction, get and set functions
  {
    bool test_photon = true;
    Photon photon;

    // set position matches get position
    std::array<double,3> pos{0.1, 0.2, 0.3};
    photon.set_position(pos);
    const auto pos_from_get = photon.get_position();

    if (pos[0] != pos_from_get[0])
      test_photon = false;
    if (pos[1] != pos_from_get[1])
      test_photon = false;
    if (pos[2] != pos_from_get[2])
      test_photon = false;

    // set angle matches get angle
    std::array<double,3> angle{0.57735, 0.37735, 0.52427};
    photon.set_angle(angle);
    const auto angle_from_get = photon.get_angle();

    if (angle[0] != angle_from_get[0])
      test_photon = false;
    if (angle[1] != angle_from_get[1])
      test_photon = false;
    if (angle[2] != angle_from_get[2])
      test_photon = false;

    // set cell matches get cell
    uint32_t cell = 129120;
    photon.set_cell(cell);

    if (photon.get_cell() != cell)
      test_photon = false;

    if (test_photon)
      cout << "TEST PASSED: Photon construction, get and set functions" << endl;
    else {
      cout << "TEST FAILED: Photon construction, get and set functions function"
           << endl;
      n_fail++;
    }
  }

  // test move function
  {
    bool test_photon_move = true;

    double tolerance = 1.0e-8;
    Photon photon;

    // first move
    {
      std::array<double,3> pos{0.0, 0.0, 0.0};
      std::array<double,3> angle{1.0, 0.0, 0.0};
      photon.set_position(pos);
      photon.set_angle(angle);
      photon.set_distance_to_census(10.0);

      photon.move(7.5);
      const auto moved_position_1 = photon.get_position();

      if (!soft_equiv(moved_position_1[0], 7.5, tolerance))
        test_photon_move = false;
      if (!soft_equiv(moved_position_1[1], 0.0, tolerance))
        test_photon_move = false;
      if (!soft_equiv(moved_position_1[2], 0.0, tolerance))
        test_photon_move = false;
      //distance remaining is distance to census - move distance
      if (!soft_equiv(photon.get_distance_remaining(), 2.5, tolerance))
        test_photon_move = false;
    }

    // second move
    {
      std::array<double,3> pos{1.64, 6.40, -5.64};
      std::array<double,3> angle{0.57735, -0.57735, 0.57735};
      photon.set_position(pos);
      photon.set_angle(angle);
      photon.set_distance_to_census(10.0);

      photon.move(1.001648);
      const auto moved_position_2 = photon.get_position();

      if (!soft_equiv(moved_position_2[0], 2.2183014728, tolerance))
        test_photon_move = false;
      if (!soft_equiv(moved_position_2[1], 5.8216985272, tolerance))
        test_photon_move = false;
      if (!soft_equiv(moved_position_2[2], -5.0616985272, tolerance))
        test_photon_move = false;
      //distance remaining is distance to census - move distance
      if (!soft_equiv(photon.get_distance_remaining(), 8.998352, tolerance))
        test_photon_move = false;
    }

    if (test_photon_move)
      cout << "TEST PASSED: Photon move function" << endl;
    else {
      cout << "TEST FAILED: Photon move function" << endl;
      n_fail++;
    }
  }

  return n_fail;
}


int test_photon_send() {

  constexpr int photon_tag{11};
  constexpr int rank_0{0};
  constexpr int rank_1{1};
  constexpr uint32_t seed{341U};
  MPI_Types mpi_types;
  MPI_Datatype MPI_Particle = mpi_types.get_particle_type();

  bool test_photon_send = true;

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  if (rank == 0) {
    std::array<double,3> pos{7.0, 8.0, 9.0};
    std::array<double,3> angle{0.574,0.574, 0.574};
    RNG rng(seed, 1000);

    Photon photon_0;
    photon_0.set_position(pos);
    photon_0.set_angle(angle);
    photon_0.set_distance_to_census(10.0);
    photon_0.set_rng(rng);

    Photon photon_1;
    MPI_Send(&photon_0, 1, MPI_Particle, rank_1, photon_tag, MPI_COMM_WORLD);
    MPI_Recv(&photon_1, 1, MPI_Particle, rank_1, photon_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // first check all the expected parts of photon 0
    auto pos_1 = photon_1.get_position();
    if(!soft_equiv(pos_1[0], 2.0, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(pos_1[1], 3.0, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(pos_1[2], 4.0, 1.0e-16))
      test_photon_send = false;
    auto angle_1 = photon_1.get_angle();
    if(!soft_equiv(angle_1[0], -0.574, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(angle_1[1], -0.574, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(angle_1[2], -0.574, 1.0e-16))
      test_photon_send = false;

    if (!soft_equiv(photon_1.get_distance_remaining(), 777.0, 1.0e-16))
      test_photon_send = false;

    // now get next three random numbers, make sure they're the expected ones
    RNG test_rng_1(seed, 2000);
    RNG &rng_1 = photon_1.get_rng();

    if(!soft_equiv(test_rng_1.generate_random_number(), rng_1.generate_random_number()))
      test_photon_send = false;
    if(!soft_equiv(test_rng_1.generate_random_number(), rng_1.generate_random_number()))
      test_photon_send = false;
    if(!soft_equiv(test_rng_1.generate_random_number(), rng_1.generate_random_number()))
      test_photon_send = false;

    if (!test_photon_send)
      std::cout<<"Failure on rank 0"<<std::endl;
  }
  if(rank == 1) {
    std::array<double,3> pos{2.0, 3.0, 4.0};
    std::array<double,3> angle{-0.574, -0.574, -0.574};
    RNG rng(seed, 2000);

    Photon photon_1;
    photon_1.set_position(pos);
    photon_1.set_angle(angle);
    photon_1.set_distance_to_census(777.0);
    photon_1.set_rng(rng);

    Photon photon_0;
    MPI_Send(&photon_1, 1, MPI_Particle, rank_0, photon_tag, MPI_COMM_WORLD);
    MPI_Recv(&photon_0, 1, MPI_Particle, rank_0, photon_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // first check all the expected parts of photon 0
    auto pos_0 = photon_0.get_position();
    if(!soft_equiv(pos_0[0], 7.0, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(pos_0[1], 8.0, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(pos_0[2], 9.0, 1.0e-16))
      test_photon_send = false;
    auto angle_0 = photon_0.get_angle();
    if(!soft_equiv(angle_0[0], 0.574, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(angle_0[1], 0.574, 1.0e-16))
      test_photon_send = false;
    if(!soft_equiv(angle_0[2], 0.574, 1.0e-16))
      test_photon_send = false;

    if (!soft_equiv(photon_0.get_distance_remaining(), 10.0, 1.0e-16))
      test_photon_send = false;

    // now get next three random numbers, make sure they're the expected ones
    RNG test_rng_0(seed, 1000);
    RNG &rng_0 = photon_0.get_rng();

    if(!soft_equiv(test_rng_0.generate_random_number(), rng_0.generate_random_number()))
      test_photon_send = false;
    if(!soft_equiv(test_rng_0.generate_random_number(), rng_0.generate_random_number()))
      test_photon_send = false;
    if(!soft_equiv(test_rng_0.generate_random_number(), rng_0.generate_random_number()))
      test_photon_send = false;

    if (!test_photon_send)
      std::cout<<"Failure on rank 1"<<std::endl;
  }

  // reduce with global AND to make sure both ranks report test_photon_send as true
  MPI_Allreduce(MPI_IN_PLACE, &test_photon_send, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

  if (test_photon_send)
    cout << "TEST PASSED: Photon send function" << endl;
  else {
    cout << "TEST FAILED: Photon send function" << endl;
  }

  return test_photon_send ? 0 : 1;
}


int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n_fail{0};

  if (rank ==0)
    n_fail += simple_photon_tests();

  MPI_Barrier(MPI_COMM_WORLD);

  n_fail += test_photon_send();

  MPI_Finalize();

  return n_fail;
}
//---------------------------------------------------------------------------//
// end of test_photon.cc
//---------------------------------------------------------------------------//
