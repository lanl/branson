//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   load_balance.h
 * \author Alex Long
 * \date   May 2 2016
 * \brief  Function for determining balance and communicating work
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef load_balance_h_
#define load_balance_h_

#include <iostream>
#include <mpi.h>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "buffer.h"
#include "constants.h"
#include "info.h"
#include "mpi_types.h"
#include "photon.h"
#include "remap_census.h"
#include "work_packet.h"

//! Load balance the work on all ranks given an array of work packets and
// census particles
void load_balance(std::vector<Work_Packet> &work,
                  std::vector<Photon> &census_list,
                  const uint64_t n_particle_on_rank, const MPI_Types &mpi_types,
                  const Info &mpi_info) {
  using std::unordered_map;
  using std::vector;
  using Constants::work_tag;
  using Constants::n_work_tag;
  using Constants::photon_tag;
  using Constants::n_photon_tag;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // get MPI datatypes
  MPI_Datatype MPI_Particle = mpi_types.get_particle_type();
  MPI_Datatype MPI_WPacket = mpi_types.get_work_packet_type();

  //--------------------------------------------------------------------------//
  // Calculate load imbalance for each rank
  //--------------------------------------------------------------------------//

  vector<int64_t> domain_particles(n_rank, 0);
  vector<int64_t> domain_delta_p(n_rank, 0);

  domain_particles[rank] = n_particle_on_rank;

  MPI_Allreduce(MPI_IN_PLACE, &domain_particles[0], n_rank, MPI_UNSIGNED_LONG,
                MPI_SUM, MPI_COMM_WORLD);

  int64_t n_global_particles =
      std::accumulate(domain_particles.begin(), domain_particles.end(), 0l);

  int64_t n_balance_particles = n_global_particles / n_rank;

  for (uint32_t ir = 0; ir < uint32_t(n_rank); ++ir)
    domain_delta_p[ir] = domain_particles[ir] - n_balance_particles;

  //--------------------------------------------------------------------------//
  // Determine which ranks to send or reiceve work from
  //--------------------------------------------------------------------------//

  // first, ignore particle imbalance on the order of 5%
  for (uint32_t ir = 0; ir < uint32_t(n_rank); ++ir) {
    // ignore particle imbalance on the order of 5%
    if (domain_delta_p[ir] < 0 &&
        domain_particles[ir] > 0.95 * n_balance_particles)
      domain_delta_p[ir] = 0;
  }

  unordered_map<int64_t, int64_t> n_send_rank;
  vector<uint64_t> work_donor_ranks;

  int64_t max_particles = 0;
  int64_t rank_delta_p, send_to_nbr;
  uint32_t max_rank = 0;

  double overload_factor = 0.25;

  int left_node, right_node;
  for (uint32_t ir = 0; ir < uint32_t(n_rank); ++ir) {
    if (domain_particles[ir] > max_particles) {
      max_particles = domain_particles[ir];
      max_rank = ir;
    }
    rank_delta_p = domain_delta_p[ir];

    // allow overload on donor ranks by effectively reducing the number of
    // particles they can donate (signed integer math)
    if (rank_delta_p > 0) {
      rank_delta_p =
          domain_particles[ir] - (1.0 + overload_factor) * n_balance_particles;
    }

    left_node = ir - 1;
    right_node = ir + 1;
    while (rank_delta_p > 0 && (left_node >= 0 || right_node < n_rank)) {
      if (left_node >= 0) {
        if (domain_delta_p[left_node] < 0 && rank_delta_p > 0) {
          send_to_nbr = abs(domain_delta_p[left_node]);
          if (send_to_nbr > rank_delta_p)
            send_to_nbr = rank_delta_p;

          // if this is your rank, remember it in your map
          if (uint32_t(rank) == ir)
            n_send_rank[left_node] = send_to_nbr;

          // if you are going to receive work, remember donor rank
          if (rank == left_node)
            work_donor_ranks.push_back(ir);

          domain_delta_p[left_node] += send_to_nbr;
          rank_delta_p -= send_to_nbr;
        }
        left_node--;
      } // if left_node >= 0
      if (right_node < n_rank) {
        if (domain_delta_p[right_node] < 0 && rank_delta_p > 0) {
          send_to_nbr = abs(domain_delta_p[right_node]);
          if (send_to_nbr > rank_delta_p)
            send_to_nbr = rank_delta_p;

          // if this is your rank, remember it in your map
          if (uint32_t(rank) == ir)
            n_send_rank[right_node] = send_to_nbr;
          // if you are going to receive work, remember donor rank

          if (rank == right_node)
            work_donor_ranks.push_back(ir);

          domain_delta_p[right_node] += send_to_nbr;
          rank_delta_p -= send_to_nbr;
        }
        right_node++;
      } // if ir+dist < n_rank
    }   // end while
  }     // end for ir in n_rank

  //--------------------------------------------------------------------------//
  // Receive work from donating ranks
  //--------------------------------------------------------------------------//

  // make MPI requests for the number of work sources
  uint32_t n_donors = work_donor_ranks.size();
  if (n_donors) {
    std::vector<MPI_Request> work_recv_request(n_donors);
    std::vector<MPI_Request> phtn_recv_request(n_donors);
    std::vector<MPI_Request> work_recv_size_request(n_donors);
    std::vector<MPI_Request> phtn_recv_size_request(n_donors);
    vector<int64_t> n_work_recv(n_donors);
    vector<int64_t> n_phtn_recv(n_donors);
    vector<Buffer<Photon>> photon_buffer(n_donors);
    vector<Buffer<Work_Packet>> work_buffer(n_donors);

    // post receives for sizes from ranks that are sending you work
    for (uint32_t i = 0; i < n_donors; ++i) {
      // post number of work packet receive
      MPI_Irecv(&n_work_recv[i], 1, MPI_INT, work_donor_ranks[i], n_work_tag,
                MPI_COMM_WORLD, &work_recv_size_request[i]);

      // post number of particles receive
      MPI_Irecv(&n_phtn_recv[i], 1, MPI_INT, work_donor_ranks[i], n_photon_tag,
                MPI_COMM_WORLD, &phtn_recv_size_request[i]);
    }

    // wait on all requests
    MPI_Waitall(n_donors, &work_recv_size_request[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(n_donors, &phtn_recv_size_request[0], MPI_STATUSES_IGNORE);

    // post receives from ranks that are sending you work
    for (uint32_t i = 0; i < n_donors; ++i) {

      // size work packet buffer
      work_buffer[i].resize(n_work_recv[i]);

      // post work packet receives
      MPI_Irecv(work_buffer[i].get_buffer(), n_work_recv[i], MPI_WPacket,
                work_donor_ranks[i], work_tag, MPI_COMM_WORLD,
                &work_recv_request[i]);

      // size photon buffer
      photon_buffer[i].resize(n_phtn_recv[i]);

      // post particle receives
      MPI_Irecv(photon_buffer[i].get_buffer(), n_phtn_recv[i], MPI_Particle,
                work_donor_ranks[i], photon_tag, MPI_COMM_WORLD,
                &phtn_recv_request[i]);
    }

    // wait on all requests
    MPI_Waitall(n_donors, &work_recv_request[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(n_donors, &phtn_recv_request[0], MPI_STATUSES_IGNORE);

    for (uint32_t i = 0; i < n_donors; ++i) {
      // add received work to your work
      const vector<Work_Packet> &temp_work = work_buffer[i].get_object();
      work.insert(work.begin(), temp_work.begin(), temp_work.end());

      // add received census photons
      const vector<Photon> &temp_photons = photon_buffer[i].get_object();
      census_list.insert(census_list.begin(), temp_photons.begin(),
                         temp_photons.end());
    }
  } // end if(n_donors)

  //--------------------------------------------------------------------------//
  // Send work to accepting ranks
  //--------------------------------------------------------------------------//

  uint32_t n_acceptors = n_send_rank.size();
  if (n_acceptors) {
    // sort work to put larger work packets near the end of the vector
    sort(work.begin(), work.end());

    MPI_Request *work_send_size_request = new MPI_Request[n_acceptors];
    MPI_Request *phtn_send_size_request = new MPI_Request[n_acceptors];
    MPI_Request *work_send_request = new MPI_Request[n_acceptors];
    MPI_Request *phtn_send_request = new MPI_Request[n_acceptors];
    uint32_t dest_rank, temp_n_send;
    uint32_t ireq = 0;            //! Request index
    uint32_t start_cut_index = 0; //! Begin slice of census list
    uint32_t n_census_remain = census_list.size(); //! Where to cut census list

    // iterate over unordered map and prepare work for other ranks
    for (auto const &imap : n_send_rank) {
      dest_rank = imap.first;
      temp_n_send = imap.second;

      vector<Work_Packet> work_to_send;
      Work_Packet temp_packet, leftover_packet;
      // first, try to send work packets to the other ranks
      while (!work.empty() && temp_n_send > 0) {
        temp_packet = work.back();
        // split work if it's larger than the number needed
        if (temp_packet.get_n_particles() > temp_n_send) {
          leftover_packet = temp_packet.split(temp_n_send);
          work[work.size() - 1] = leftover_packet;
        }
        // otherwise, pop the temp work packet off the stack
        else
          work.pop_back();

        // add packet to send list
        work_to_send.push_back(temp_packet);
        // subtract particle in packet from temp_n_send
        temp_n_send -= temp_packet.get_n_particles();
      }

      // send census particles instead (not preferred, these photons are likely
      // to travel farther and thus require more memory)
      uint32_t n_send_census = 0;
      if (temp_n_send > 0 && n_census_remain > temp_n_send) {
        n_send_census = temp_n_send;
        start_cut_index = n_census_remain - n_send_census;
        n_census_remain -= n_send_census;
      }

      int32_t n_work_send = work_to_send.size();

      // send sizes for work and photon vectors
      MPI_Isend(&n_work_send, 1, MPI_INT, dest_rank, n_work_tag, MPI_COMM_WORLD,
                &work_send_size_request[ireq]);
      MPI_Isend(&n_send_census, 1, MPI_INT, dest_rank, n_photon_tag,
                MPI_COMM_WORLD, &phtn_send_size_request[ireq]);

      // send both work and photon vectors, even if they're empty
      MPI_Isend(&work_to_send[0], n_work_send, MPI_WPacket, dest_rank, work_tag,
                MPI_COMM_WORLD, &work_send_request[ireq]);
      MPI_Isend(&census_list[start_cut_index], n_send_census, MPI_Particle,
                dest_rank, photon_tag, MPI_COMM_WORLD,
                &phtn_send_request[ireq]);

      // wait on these requests so buffers are not invalidated when going
      // out of loop scope
      MPI_Wait(&work_send_size_request[ireq], MPI_STATUS_IGNORE);
      MPI_Wait(&phtn_send_size_request[ireq], MPI_STATUS_IGNORE);
      MPI_Wait(&work_send_request[ireq], MPI_STATUS_IGNORE);
      MPI_Wait(&phtn_send_request[ireq], MPI_STATUS_IGNORE);

      // increment the request index counter
      ireq++;
    }

    // remove census photons that were sent off
    census_list.erase(census_list.begin() + n_census_remain, census_list.end());

    // clean up
    delete[] work_send_size_request;
    delete[] phtn_send_size_request;
    delete[] work_send_request;
    delete[] phtn_send_request;

  } // end if(n_acceptors)

  if (rank == 0) {
    std::cout << "Load balancing " << n_global_particles << " particles"
              << std::endl;
    std::cout << "Max particles: " << max_particles << " on rank " << max_rank
              << ", ";
    std::cout << 100.0 * double(max_particles) / double(n_global_particles);
    std::cout << " percent of total work" << std::endl;
  }
}

//! Use a binary tree approach to load balance work and census particles
void bt_load_balance(std::vector<Work_Packet> &work,
                     std::vector<Photon> &census_list,
                     const uint64_t n_particle_on_rank, MPI_Types *mpi_types,
                     const Info &mpi_info) {
  using std::vector;
  using std::unordered_map;
  using Constants::photon_tag;
  using Constants::n_photon_tag;
  using Constants::work_tag;
  using Constants::n_work_tag;

  const int n_tag(100);

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  MPI_Datatype MPI_Particle = mpi_types->get_particle_type();
  MPI_Datatype MPI_WPacket = mpi_types->get_work_packet_type();

  // split the communicator based on the node number
  MPI_Comm local_comm;
  MPI_Comm_split(MPI_COMM_WORLD, mpi_info.get_color(), rank, &local_comm);
  int n_rank_local;
  MPI_Comm_size(local_comm, &n_rank_local);

  // sort the census vector by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // sort work to put larger work packets near the end of the vector
  sort(work.begin(), work.end());

  // do the binary tree pattern to the nearest log2(rank), rounding up
  int32_t n_levels = int32_t(ceil(log2(n_rank)));

  int32_t r_partner;

  vector<MPI_Request> n_p_reqs(2);
  vector<MPI_Request> donor_send_reqs(2);
  vector<MPI_Request> acceptor_recv_reqs(2);

  // start with the current number of particles on this rank
  uint64_t balanced_rank_particles = n_particle_on_rank;
  uint64_t partner_rank_particles, avg_particles;
  int64_t n_send_census;
  int64_t temp_n_send;
  bool balanced;
  uint64_t start_cut_index = 0; //! Begin slice of census list

  vector<Work_Packet> work_to_send;
  vector<Photon> census_to_send;

  Buffer<Photon> recv_photon_buffer;
  Buffer<Work_Packet> recv_work_buffer;

  uint64_t n_census_remain = census_list.size(); //! Where to cut census list

  // run the binary tree comm pattern going from smallest communication
  // distance to largest. Don't communicate if your partner doesn't exist
  // (this will occur if log_2(n_rank) is not an integer)
  for (int32_t k = 0; k < n_levels; ++k) {

    // get current communication partner
    r_partner = get_pairing(rank, n_rank, k);

    // only participate if the partner rank exists
    if (r_partner <= n_rank) {
      // send and receive number of photons on your rank
      MPI_Isend(&balanced_rank_particles, 1, MPI_UNSIGNED_LONG, r_partner,
                n_tag, MPI_COMM_WORLD, &n_p_reqs[0]);
      MPI_Irecv(&partner_rank_particles, 1, MPI_UNSIGNED_LONG, r_partner, n_tag,
                MPI_COMM_WORLD, &n_p_reqs[1]);
      MPI_Waitall(2, &n_p_reqs[0], MPI_STATUSES_IGNORE);

      // check for balance with unsigned integer math
      avg_particles = 0.5 * (partner_rank_particles + balanced_rank_particles);

      // allow imbalance of 5%
      if (balanced_rank_particles > 1.05 * avg_particles ||
          balanced_rank_particles < 0.95 * avg_particles) {
        balanced = false;
      } else
        balanced = true;

      // if ranks are not balanced, exchange work
      if (!balanced) {

        // logic for donor ranks
        if (balanced_rank_particles > partner_rank_particles) {
          temp_n_send =
              int64_t(balanced_rank_particles) - int64_t(avg_particles);

          Work_Packet temp_packet, leftover_packet;

          // first, try to send work packets to the other ranks
          while (!work.empty() && temp_n_send > 0) {
            temp_packet = work.back();
            // split work if it's larger than the number needed
            if (temp_packet.get_n_particles() > temp_n_send) {
              leftover_packet = temp_packet.split(temp_n_send);
              work[work.size() - 1] = leftover_packet;
            }
            // otherwise, pop the temp work packet off the stack
            else
              work.pop_back();

            // add packet to send list
            work_to_send.push_back(temp_packet);

            // subtract particle in packet from temp_n_send and on rank
            // particles
            temp_n_send -= temp_packet.get_n_particles();
            balanced_rank_particles -= temp_packet.get_n_particles();
          }

          // send census particles instead (not preferred, these photons are
          // likely
          // to travel farther and thus require more memory)
          n_send_census = 0;
          if (temp_n_send > 0 && n_census_remain > temp_n_send) {
            n_send_census = temp_n_send;
            start_cut_index = n_census_remain - n_send_census;
            n_census_remain -= n_send_census;
          }

          // reduce the number of particle on rank by the size of the census
          balanced_rank_particles -= n_send_census;

          // send the size of each message
          int32_t n_send_work = work_to_send.size();
          MPI_Isend(&n_send_work, 1, MPI_INT, r_partner, n_work_tag,
                    MPI_COMM_WORLD, &donor_send_reqs[0]);
          MPI_Isend(&n_send_census, 1, MPI_INT, r_partner, n_photon_tag,
                    MPI_COMM_WORLD, &donor_send_reqs[1]);
          MPI_Waitall(2, &donor_send_reqs[0], MPI_STATUSES_IGNORE);

          // send both work and photon vectors, even if they're empty
          MPI_Isend(&work_to_send[0], work_to_send.size(), MPI_WPacket,
                    r_partner, work_tag, MPI_COMM_WORLD, &donor_send_reqs[0]);
          MPI_Isend(&census_list[start_cut_index], n_send_census, MPI_Particle,
                    r_partner, photon_tag, MPI_COMM_WORLD, &donor_send_reqs[1]);

          MPI_Waitall(2, &donor_send_reqs[0], MPI_STATUSES_IGNORE);

          // clear buffer for sending work
          work_to_send.clear();

          // remove census photons that were sent off
          census_list.erase(census_list.begin() + n_census_remain,
                            census_list.end());
          n_census_remain = census_list.size();
        }

        // logic for acceptor rank
        else {

          int32_t n_work_recv, n_census_recv;
          MPI_Irecv(&n_work_recv, 1, MPI_INT, r_partner, n_work_tag,
                    MPI_COMM_WORLD, &acceptor_recv_reqs[0]);
          MPI_Irecv(&n_census_recv, 1, MPI_INT, r_partner, n_photon_tag,
                    MPI_COMM_WORLD, &acceptor_recv_reqs[1]);
          MPI_Waitall(2, &acceptor_recv_reqs[0], MPI_STATUSES_IGNORE);

          // resize buffers
          recv_work_buffer.resize(n_work_recv);
          recv_photon_buffer.resize(n_census_recv);

          // post work packet receives
          MPI_Irecv(recv_work_buffer.get_buffer(), n_work_recv, MPI_WPacket,
                    r_partner, work_tag, MPI_COMM_WORLD,
                    &acceptor_recv_reqs[0]);

          // post particle receives
          MPI_Irecv(recv_photon_buffer.get_buffer(), n_census_recv,
                    MPI_Particle, r_partner, photon_tag, MPI_COMM_WORLD,
                    &acceptor_recv_reqs[1]);

          MPI_Waitall(2, &acceptor_recv_reqs[0], MPI_STATUSES_IGNORE);

          // add received work to your work
          const vector<Work_Packet> &temp_work = recv_work_buffer.get_object();
          work.insert(work.begin(), temp_work.begin(),
                      temp_work.begin() + n_work_recv);

          // sort, if work was received
          if (!temp_work.empty())
            sort(work.begin(), work.end());

          // add received census photons
          const vector<Photon> &temp_photons = recv_photon_buffer.get_object();
          census_list.insert(census_list.begin(), temp_photons.begin(),
                             temp_photons.begin() + n_census_recv);

          // sort, if census photons were received
          if (!temp_photons.empty())
            sort(census_list.begin(), census_list.end());

          // get current count photon count and update
          n_census_remain = census_list.size();
          balanced_rank_particles = census_list.size();
          for (auto const &i_w : work)
            balanced_rank_particles += i_w.get_n_particles();
        }
      } // if balanced
    }   // end if r_partner < n_rank
  }     // end loop over levels
}

#endif // load_balance_h_

//----------------------------------------------------------------------------//
// end of load_balance.h
//----------------------------------------------------------------------------//
