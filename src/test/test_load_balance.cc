//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_load_balance.cc
 * \author Alex Long
 * \date   Novemeber 9, 2016
 * \brief  Test load balance function
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "../constants.h"
#include "../info.h"
#include "../load_balance.h"
#include "testing_functions.h"

int main (int argc, char *argv[]) {

  using std::cout;
  using std::endl;
  using std::vector;

  MPI_Init(&argc, &argv);

  MPI_Types *mpi_types = new MPI_Types();
  const Info mpi_info;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  int nfail = 0;

  //test load balance function with standard off-rank photons
  {

    bool test_load_balance = true;

    uint64_t n_particles_on_rank = 0;
    // make photons and work on rank 0
    vector<Photon> census;
    uint64_t n_total_census_photons = 10000;
    if (rank ==0) {
      for (uint32_t i=0;i<n_total_census_photons;++i) {
        Photon temp_photon;
        double pos[3] =   {0.0, 0.0, 0.0};
        double angle[3] = {1.0, 0.0, 0.0};
        temp_photon.set_position(pos);
        temp_photon.set_angle(angle);
        temp_photon.set_distance_to_census(10.0);
        temp_photon.set_cell(10);
        census.push_back(temp_photon);
        n_particles_on_rank++;
      }
    }

    // make work on rank 0
    vector<Work_Packet> work;
    uint32_t n_work_packets = 100;
    uint32_t n_particles = 100;
    uint32_t g_cell_ID = 655;
    uint32_t g_grip_ID = 540;
    double packet_E = 10.0;
    if (rank ==0) {
      for (uint32_t i=0;i<n_work_packets;++i) {

        Work_Packet temp_work;
        temp_work.set_global_cell_ID(g_cell_ID);
        temp_work.set_global_grip_ID(g_grip_ID);
        temp_work.attach_creation_work(packet_E, n_particles);
        temp_work.set_source_type(Constants::EMISSION);
        work.push_back(temp_work);
        n_particles_on_rank+=n_particles;
      }
    }

    load_balance(work, census, n_particles_on_rank, mpi_types, mpi_info);

    uint64_t n_post_balanced_particles = 0;

    for (auto w_itr = work.begin(); w_itr != work.end(); ++w_itr)
      n_post_balanced_particles += w_itr->get_n_particles();
    n_post_balanced_particles += census.size();

    if (work.size() != 0) test_load_balance = false;
    if (census.size() != 0) test_load_balance = false;
    if (work.front().get_source_type() != Constants::INITIAL_CENSUS)
      test_load_balance=false;

    // balanced particles should be n_particles_on_rank for 0 divided
    // by n_rank
    if (n_post_balanced_particles != 20000) test_load_balance=false;

    if (test_load_balance) {
      cout<<"TEST PASSED: load_balance all work on rank 0"<<endl;
    }
    else {
      cout<<"TEST FAILED: load_balance all work on rank 0"<<endl;
      nfail++;
    }
  }

  delete mpi_types;

  MPI_Finalize();

  return nfail;
}
//---------------------------------------------------------------------------//
// end of test_load_balance.cc
//---------------------------------------------------------------------------//
