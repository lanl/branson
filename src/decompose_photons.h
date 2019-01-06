//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   decompose_photons.h
 * \author Alex Long
 * \date   June 17 2015
 * \brief  Load balance census photons after a step in mesh passing DD mode
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef decompose_photons_h_
#define decompose_photons_h_

#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <unordered_map>
#include <vector>

#include "mpi_types.h"
#include "photon.h"
#include "timer.h"

void print_MPI_photons(const std::vector<Photon> &phtn_vec,
                       const uint32_t &rank, const uint32_t &size) {

  using std::cout;

  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank < size; ++p_rank) {
    if (rank == p_rank) {
      for (uint32_t i = 0; i < phtn_vec.size(); ++i)
        phtn_vec[i].print_info(rank);
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
  usleep(100);
  cout.flush();
  usleep(100);
}

//! Rebalance a census that only include off-rank census photons
std::vector<Photon> rebalance_census(std::vector<Photon> &off_rank_census,
                                     const Mesh &mesh,
                                     const MPI_Types &mpi_types) {
  using std::sort;
  using std::unordered_map;
  using std::vector;

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);
  uint32_t n_off_rank = n_rank - 1;

  MPI_Datatype MPI_Particle = mpi_types.get_particle_type();

  // make off processor map
  unordered_map<int, int> proc_map;
  for (int i = 0; i < int(n_off_rank); ++i) {
    int r_index = i + int(i >= rank);
    proc_map[i] = r_index;
  }

  // sort the census vector by cell ID (global ID)
  sort(off_rank_census.begin(), off_rank_census.end());

  // size of census list
  uint32_t n_census = off_rank_census.size();

  // count the photons belonging to each rank and the start index of each
  // count the ranks that you will send to, add them to a vector
  vector<uint32_t> rank_count(n_rank, 0);
  vector<uint32_t> rank_start(n_rank + 1, 0);
  vector<bool> rank_found(n_rank, false);
  uint32_t r;
  for (uint32_t i = 0; i < n_census; ++i) {
    r = mesh.get_rank(off_rank_census[i].get_cell());
    rank_count[r]++;
    if (rank_found[r] == false) {
      rank_found[r] = true;
      rank_start[r] = i;
    }
  }

  // end of rank count is the total number of census photons
  rank_start[n_rank] = n_census;

  // make requests for non-blocking communication
  MPI_Request *reqs = new MPI_Request[n_off_rank * 2];

  // make n_off_rank receive buffers
  vector<vector<Photon>> recv_photons;
  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    vector<Photon> empty_vec;
    recv_photons.push_back(empty_vec);
  }

  // get the number of photons received from each rank
  vector<int> recv_from_rank(n_off_rank, 0);

  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    int off_rank = proc_map[ir];
    MPI_Isend(&rank_count[off_rank], 1, MPI_UNSIGNED, off_rank, 0,
              MPI_COMM_WORLD, &reqs[ir]);
    MPI_Irecv(&recv_from_rank[ir], 1, MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD,
              &reqs[ir + n_off_rank]);
  }

  MPI_Waitall(n_off_rank * 2, reqs, MPI_STATUS_IGNORE);

  // now send the buffers and post receives
  // resize receive buffers with recv_from_rank
  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    int off_rank = proc_map[ir];
    int start_copy = rank_start[off_rank];
    MPI_Isend(&off_rank_census[start_copy], rank_count[off_rank], MPI_Particle,
              off_rank, 0, MPI_COMM_WORLD, &reqs[ir]);
    recv_photons[ir].resize(recv_from_rank[ir]);
    MPI_Irecv(&recv_photons[ir][0], recv_from_rank[ir], MPI_Particle, off_rank,
              0, MPI_COMM_WORLD, &reqs[ir + n_off_rank]);
  }

  MPI_Waitall(n_off_rank * 2, reqs, MPI_STATUS_IGNORE);

  // free memory from off rank census list
  off_rank_census.clear();

  // copy received census photons to a new census list
  vector<Photon> new_on_rank_census;
  for (uint32_t ir = 0; ir < uint32_t(n_rank - 1); ++ir) {
    new_on_rank_census.insert(new_on_rank_census.end(),
                              recv_photons[ir].begin(), recv_photons[ir].end());
  }

  // explicitly delete the MPI requests
  delete[] reqs;

  return new_on_rank_census;
}

//! Rebalance a census that contains photons living on mesh owned by any rank
std::vector<Photon> rebalance_raw_census(std::vector<Photon> &census,
                                         const Mesh &mesh,
                                         const MPI_Types &mpi_types) {
  using std::sort;
  using std::unordered_map;
  using std::vector;

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  MPI_Datatype MPI_Particle = mpi_types.get_particle_type();

  // sort ceneus by cell ID
  sort(census.begin(), census.end());

  // make a list of off rank photons
  std::map<int32_t, uint64_t> rank_start;
  vector<uint64_t> n_census_on_rank(n_rank, 0);
  vector<uint64_t> n_census_on_rank_global(n_rank, 0);
  int census_p_rank;
  uint64_t ip = 0;
  for (auto &p : census) {
    census_p_rank = mesh.get_rank(p.get_cell());
    if (rank_start.find(census_p_rank) == rank_start.end()) {
      rank_start[census_p_rank] = ip;
    }
    n_census_on_rank[census_p_rank]++;
    ip++;
  }

  // if you have no census on your rank make sure it's recorded
  if (rank_start.find(rank) == rank_start.end())
    rank_start[rank] = 0;

  MPI_Allreduce(&n_census_on_rank[0], &n_census_on_rank_global[0], n_rank,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  uint64_t avg_census = std::accumulate(n_census_on_rank_global.begin(),
                                        n_census_on_rank_global.end(), 0ul) /
                        n_rank;
  uint64_t max_census = *std::max_element(n_census_on_rank_global.begin(),
                                          n_census_on_rank_global.end());

  if (rank == 0) {
    std::cout << "Max census on rank: " << max_census;
    std::cout << ", average census: " << avg_census << std::endl;
  }

  // check to see how many we should send back
  // remove from map if larger than avg_census
  std::unordered_map<int32_t, uint64_t> n_to_send;
  std::set<int32_t> acceptor_ranks;
  for (auto ir : rank_start) {
    if (n_census_on_rank_global[ir.first] < avg_census && ir.first != rank) {
      acceptor_ranks.insert(ir.first);
    }
  }

  int n_global_acceptors = acceptor_ranks.size();
  MPI_Allreduce(MPI_IN_PLACE, &n_global_acceptors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  // get the number of ranks that will send us a message
  const int one = 1;
  int *n_donors_win;
  int n_donors(0);

  // create the RMA memory windows for each value
  MPI_Win win;
  MPI_Request req;
  MPI_Win_allocate(1 * sizeof(int), 1 * sizeof(int), MPI_INFO_NULL,
                   MPI_COMM_WORLD, &n_donors_win, &win);
  n_donors_win[0] = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  int assert = MPI_MODE_NOCHECK; // no conflicting locks on this window
  MPI_Win_lock_all(assert, win);
  for (auto ir : acceptor_ranks) {
    // increment the remote number of receives
    MPI_Raccumulate(&one, 1, MPI_INT, ir, 0, 1, MPI_INT, MPI_SUM, win, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }
  MPI_Win_unlock_all(win);
  MPI_Barrier(MPI_COMM_WORLD);
  n_donors = n_donors_win[0];
  MPI_Win_free(&win);

  int n_global_donors = n_donors;
  MPI_Allreduce(MPI_IN_PLACE, &n_global_donors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "Total acceptors: " << n_global_acceptors;
    std::cout << ", total donors: " << n_global_donors << std::endl;
  }

  // make requests for non-blocking communication
  MPI_Request *s_reqs = new MPI_Request[acceptor_ranks.size()];
  MPI_Request *r_reqs = new MPI_Request[n_donors];
  MPI_Status *r_status = new MPI_Status[n_donors];

  // make n_donors receive buffers
  vector<vector<Photon>> recv_photons(n_donors);
  vector<int> recv_from_rank(n_donors, 0);

  // send the number of photons to your acceptor ranks
  int index = 0;
  for (auto ir : acceptor_ranks) {
    MPI_Isend(&n_census_on_rank[ir], 1, MPI_UNSIGNED, ir, 0, MPI_COMM_WORLD,
              &s_reqs[index]);
    index++;
  }

  // get the number of photons received from donor ranks
  for (int i = 0; i < n_donors; ++i) {
    MPI_Irecv(&recv_from_rank[i], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &r_reqs[i]);
  }

  MPI_Waitall(acceptor_ranks.size(), s_reqs, MPI_STATUS_IGNORE);
  MPI_Waitall(n_donors, r_reqs, r_status);
  MPI_Barrier(MPI_COMM_WORLD);

  std::map<int, int> donor_rank_to_n_recv;
  for (int i = 0; i < n_donors; ++i) {
    donor_rank_to_n_recv[r_status[i].MPI_SOURCE] = recv_from_rank[i];
  }

  // now send the buffers and post receives
  // resize receive buffers with recv_from_rank
  index = 0;
  for (auto &ir : donor_rank_to_n_recv) {
    recv_photons[index].resize(ir.second);
    MPI_Irecv(&recv_photons[index][0], ir.second, MPI_Particle, ir.first, 0,
              MPI_COMM_WORLD, &r_reqs[index]);
    index++;
  }

  index = 0;
  for (auto &ir : acceptor_ranks) {
    MPI_Isend(&census[rank_start[ir]], n_census_on_rank[ir], MPI_Particle, ir,
              0, MPI_COMM_WORLD, &s_reqs[index]);
    index++;
  }

  MPI_Waitall(acceptor_ranks.size(), s_reqs, MPI_STATUS_IGNORE);
  MPI_Waitall(n_donors, r_reqs, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

  // erase the parts of the census that were sent, reverse iterate over the
  // starting index so rank_start indices are valid after erasing
  for (auto ir = rank_start.crbegin(); ir != rank_start.crend(); ir++) {
    if (acceptor_ranks.count(ir->first)) {
      census.erase(census.begin() + ir->second,
                   census.begin() + ir->second + n_census_on_rank[ir->first]);
    }
  }

  // copy received census photons to the end of the census list
  vector<Photon> new_on_rank_census;
  for (int i = 0; i < n_donors; ++i) {
    census.insert(census.end(), recv_photons[i].begin(), recv_photons[i].end());
  }

  census.shrink_to_fit();

  // explicitly delete the MPI requests
  delete[] r_status;
  delete[] r_reqs;
  delete[] s_reqs;

  return census;
}

/*
  // determine how many receives to post with an RMA call
  const int one = 1;
  int num_recv(0); // return value

  // Create the RMA memory windows for each value
  MPI_Win win;
  MPI_Win_create(&num_recv, 1 * sizeof(int), sizeof(int), MPI_INFO_NULL,
                 MPI_COMM_WORLD, &win);

  // Assertion value for fences.  Currently, we effectively don't set
  // anything (zero).
  const int fence_assert = 0;

  // Accumulate the local and remote data values
  MPI_Win_fence(fence_assert, win);
  for (auto ir : receiver_ranks) {
      // ...increment the remote number of receives
      MPI_Accumulate(&one, 1, MPI_Traits<int>::element_type(), ir, 0, 1,
                     MPI_Traits<int>::element_type(), MPI_SUM, win);
  }
  MPI_Win_fence(fence_assert, win);
  MPI_Win_free(&win);

  // make requests for non-blocking communication
  MPI_Request* reqs = new MPI_Request[num_recv];
*/

#endif // decompose_photons_h
//---------------------------------------------------------------------------//
// end of decompose_photons.h
//---------------------------------------------------------------------------//
