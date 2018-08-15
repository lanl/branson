//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   particle_pass_transport.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  IMC transport with particle passing method
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_particle_pass_h_
#define transport_particle_pass_h_

#include <algorithm>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <stack>
#include <unordered_map>
#include <vector>

#include "RNG.h"
#include "buffer.h"
#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "photon.h"
#include "sampling_functions.h"

Constants::event_type transport_photon_particle_pass(
    Photon &phtn, Mesh *mesh, RNG *rng, double &next_dt, double &exit_E,
    double &census_E, std::vector<double> &rank_abs_E,
    std::vector<double> &rank_track_E) {
  using Constants::ELEMENT;
  using Constants::PROCESSOR;
  using Constants::REFLECT;
  using Constants::VACUUM;
  // events
  using Constants::bc_type;
  using Constants::c;
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::PASS;
  using std::min;

  uint32_t cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;
  Cell cell;

  uint32_t surface_cross = 0;
  double cutoff_fraction = 0.01; // note: get this from IMC_state

  cell_id = phtn.get_cell();
  cell = mesh->get_on_rank_cell(cell_id);
  bool active = true;

  // transport this photon
  while (active) {
    group = phtn.get_group();
    sigma_a = cell.get_op_a(group);
    sigma_s = cell.get_op_s(group);
    f = cell.get_f();

    // get distance to event
    dist_to_scatter =
        -log(rng->generate_random_number()) / ((1.0 - f) * sigma_a + sigma_s);

    dist_to_boundary = cell.get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);

    rank_track_E[cell_id] += absorbed_E / (sigma_a * f);
    rank_abs_E[cell_id] += absorbed_E;

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      rank_abs_E[cell_id] += phtn.get_E();
      active = false;
      event = KILL;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell.get_bc(surface_cross);
        if (boundary_event == ELEMENT) {
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          cell_id = next_cell;
          cell = mesh->get_on_rank_cell(cell_id);
        } else if (boundary_event == PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          event = PASS;
        } else if (boundary_event == VACUUM) {
          exit_E += phtn.get_E();
          active = false;
          event = EXIT;
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        phtn.set_distance_to_census(c * next_dt);
        active = false;
        event = CENSUS;
        census_E += phtn.get_E();
      }
    } // end event loop
  }   // end while alive
  return event;
}

std::vector<Photon>
particle_pass_transport(Source &source, Mesh *mesh, IMC_State *imc_state,
                        IMC_Parameters *imc_parameters, MPI_Types *mpi_types,
                        Message_Counter &mctr, std::vector<double> &rank_abs_E,
                        std::vector<double> &rank_track_E,
                        const Info &mpi_info) {
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::PASS;
  using Constants::photon_tag;
  using Constants::WAIT;
  using std::cout;
  using std::endl;
  using std::stack;
  using std::unordered_map;
  using std::vector;

  double census_E = 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt(); //! Set for census photons
  double dt = imc_state->get_next_dt();      //! For making current photons

  RNG *rng = imc_state->get_rng();

  // timing
  Timer t_transport;
  Timer t_mpi;
  t_transport.start_timer("timestep transport");

  // Number of particles to run between MPI communication
  const uint32_t batch_size = imc_parameters->get_batch_size();

  // Preferred size of MPI message
  const uint32_t max_buffer_size = imc_parameters->get_particle_message_size();

  MPI_Datatype MPI_Particle = mpi_types->get_particle_type();

  // get global photon count
  uint64_t n_local = source.get_n_photon();
  uint64_t n_global;
  uint64_t last_global_complete_count;

  MPI_Allreduce(&n_local, &n_global, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);

  // This flag indicates that send processing is needed for target rank
  vector<vector<Photon>> send_list;

  // Completion count request made flag
  bool req_made = false;
  int recv_allreduce_flag;

  // get adjacent processor map (off_rank_id -> adjacent_proc_number)
  auto adjacent_procs = mesh->get_proc_adjacency_list();
  uint32_t n_adjacent = adjacent_procs.size();
  // messsage requests for photon sends and receives
  MPI_Request *phtn_recv_request = new MPI_Request[n_adjacent];
  MPI_Request *phtn_send_request = new MPI_Request[n_adjacent];
  // message request for non-blocking allreduce
  MPI_Request completion_request;
  // make a send/receive particle buffer for each adjacent processor
  vector<Buffer<Photon>> phtn_recv_buffer(n_adjacent);
  vector<Buffer<Photon>> phtn_send_buffer(n_adjacent);

  // Post receives for photons from adjacent sub-domains
  {
    uint32_t i_b; // buffer index
    int adj_rank; // adjacent rank
    for (auto const &it : adjacent_procs) {
      adj_rank = it.first;
      i_b = it.second;
      // push back send and receive lists
      vector<Photon> empty_phtn_vec;
      send_list.push_back(empty_phtn_vec);
      // make receive buffer the appropriate size
      phtn_recv_buffer[i_b].resize(max_buffer_size);
      MPI_Irecv(phtn_recv_buffer[i_b].get_buffer(), max_buffer_size,
                MPI_Particle, adj_rank, photon_tag, MPI_COMM_WORLD,
                &phtn_recv_request[i_b]);
      mctr.n_receives_posted++;
      phtn_recv_buffer[i_b].set_awaiting();
    } // end loop over adjacent processors
  }

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  vector<Photon> census_list;    //!< End of timestep census list
  stack<Photon> phtn_recv_stack; //!< Stack of received photons

  int send_rank;
  uint64_t n_complete = 0; //!< Completed histories, regardless of origin
  //! Send and receive buffers for complete count
  uint64_t s_global_complete, r_global_complete;
  uint64_t n_local_sourced = 0; //!< Photons pulled from source object
  bool from_receive_stack = false;
  int send_req_flag;
  Photon phtn;
  event_type event;

  while (last_global_complete_count != n_global) {

    uint32_t n = batch_size;

    //------------------------------------------------------------------------//
    // Transport photons from source and received list
    //------------------------------------------------------------------------//
    // first, try to transport photons from the received list
    while (n && (!phtn_recv_stack.empty() || (n_local_sourced < n_local))) {

      if (!phtn_recv_stack.empty()) {
        phtn = phtn_recv_stack.top();
        from_receive_stack = true;
      } else {
        phtn = source.get_photon(rng, dt);
        n_local_sourced++;
        from_receive_stack = false;
      }

      event = transport_photon_particle_pass(
          phtn, mesh, rng, next_dt, exit_E, census_E, rank_abs_E, rank_track_E);
      switch (event) {
      // this case should never be reached
      case WAIT:
        break;
      case KILL:
        n_complete++;
        break;
      case EXIT:
        n_complete++;
        break;
      case CENSUS:
        census_list.push_back(phtn);
        n_complete++;
        break;
      case PASS:
        t_mpi.start_timer("timestep mpi");
        send_rank = mesh->get_rank(phtn.get_cell());
        int i_b = adjacent_procs[send_rank];
        send_list[i_b].push_back(phtn);

        // test completion of send buffer
        if (phtn_send_buffer[i_b].sent()) {
          MPI_Test(&phtn_send_request[i_b], &send_req_flag, MPI_STATUS_IGNORE);
          if (send_req_flag) {
            phtn_send_buffer[i_b].reset();
            mctr.n_sends_completed++;
          }
        }

        // send full photon buffers if phtn_send_buffer is empty (send has
        // completed)
        if ((phtn_send_buffer[i_b].empty() &&
             send_list[i_b].size() >= max_buffer_size)) {
          vector<Photon>::iterator copy_start = send_list[i_b].begin();
          vector<Photon>::iterator copy_end =
              send_list[i_b].begin() + max_buffer_size;
          vector<Photon> send_now_list(copy_start, copy_end);
          send_list[i_b].erase(copy_start, copy_end);
          phtn_send_buffer[i_b].fill(send_now_list);
          MPI_Isend(phtn_send_buffer[i_b].get_buffer(), max_buffer_size,
                    MPI_Particle, send_rank, photon_tag, MPI_COMM_WORLD,
                    &phtn_send_request[i_b]);
          phtn_send_buffer[i_b].set_sent();
          // update counters
          mctr.n_particles_sent += max_buffer_size;
          mctr.n_sends_posted++;
          mctr.n_particle_messages++;
        }
        t_mpi.stop_timer("timestep mpi");
        break;
      }
      n--;
      if (from_receive_stack)
        phtn_recv_stack.pop();
    }

    //------------------------------------------------------------------------//
    // process photon send and receives
    //------------------------------------------------------------------------//
    t_mpi.start_timer("timestep mpi");
    {
      int send_req_flag;
      int recv_req_flag;
      int recv_count; // recieve count is 32 bit

      MPI_Status recv_status;
      uint32_t i_b; // buffer index
      int adj_rank; // adjacent rank
      for (auto const &it : adjacent_procs) {
        adj_rank = it.first;
        i_b = it.second;

        // sends are handled above, unless we're done with local work and banks
        // are empty
        if (n_local_sourced == n_local && phtn_recv_stack.empty()) {
          // test completion of send buffer
          if (phtn_send_buffer[i_b].sent()) {
            MPI_Test(&phtn_send_request[i_b], &send_req_flag,
                     MPI_STATUS_IGNORE);
            if (send_req_flag) {
              phtn_send_buffer[i_b].reset();
              mctr.n_sends_completed++;
            }
          }

          // send full photon buffers if send_buffer is empty and send_list has
          // some photons in it
          if (phtn_send_buffer[i_b].empty() && !send_list[i_b].empty()) {
            uint32_t n_photons_to_send = max_buffer_size;
            if (send_list[i_b].size() < max_buffer_size)
              n_photons_to_send = send_list[i_b].size();
            vector<Photon>::iterator copy_start = send_list[i_b].begin();
            vector<Photon>::iterator copy_end =
                send_list[i_b].begin() + n_photons_to_send;
            vector<Photon> send_now_list(copy_start, copy_end);
            send_list[i_b].erase(copy_start, copy_end);
            phtn_send_buffer[i_b].fill(send_now_list);
            MPI_Isend(phtn_send_buffer[i_b].get_buffer(), n_photons_to_send,
                      MPI_Particle, adj_rank, photon_tag, MPI_COMM_WORLD,
                      &phtn_send_request[i_b]);
            phtn_send_buffer[i_b].set_sent();
            // update counters
            mctr.n_particles_sent += n_photons_to_send;
            mctr.n_sends_posted++;
            mctr.n_particle_messages++;
          }
        }

        // process receive buffer
        if (phtn_recv_buffer[i_b].awaiting()) {
          MPI_Test(&phtn_recv_request[i_b], &recv_req_flag, &recv_status);
          if (recv_req_flag) {
            vector<Photon> &receive_list = phtn_recv_buffer[i_b].get_object();
            // only push the number of received photons onto the recv_stack
            MPI_Get_count(&recv_status, MPI_Particle, &recv_count);
            for (uint32_t i = 0; i < uint32_t(recv_count); ++i)
              phtn_recv_stack.push(receive_list[i]);
            phtn_recv_buffer[i_b].reset();
            // post receive again, don't resize--it's already set to maximum
            MPI_Irecv(phtn_recv_buffer[i_b].get_buffer(), max_buffer_size,
                      MPI_Particle, adj_rank, photon_tag, MPI_COMM_WORLD,
                      &phtn_recv_request[i_b]);
            phtn_recv_buffer[i_b].set_awaiting();
            mctr.n_receives_completed++;
            mctr.n_receives_posted++;
          }
        }
      } // end loop over adjacent processors
    }   // end scope of particle passing

    //------------------------------------------------------------------------//
    // binary tree completion communication
    //------------------------------------------------------------------------//

    if (!req_made) {
      s_global_complete = n_complete;
      MPI_Iallreduce(&s_global_complete, &r_global_complete, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD,
                     &completion_request);
      req_made = true;
    } else {
      MPI_Test(&completion_request, &recv_allreduce_flag, MPI_STATUS_IGNORE);
      if (recv_allreduce_flag) {
        last_global_complete_count = r_global_complete;
        s_global_complete = n_complete;
        if (last_global_complete_count != n_global) {
          MPI_Iallreduce(&s_global_complete, &r_global_complete, 1,
                         MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD,
                         &completion_request);
        }
      }
    }

    t_mpi.stop_timer("timestep mpi");

  } // end while

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish then send empty photon messages
  // do this because it's possible for a rank to receive the empty message
  // while it's still in the transport loop. In that case, it will post a
  // receive again, which will never have a matching send
  MPI_Barrier(MPI_COMM_WORLD);

  // finish off posted photon receives
  {
    vector<Photon> one_photon(1);
    uint32_t i_b; // buffer index
    int adj_rank; // adjacent rank
    for (auto const &it : adjacent_procs) {
      adj_rank = it.first;
      i_b = it.second;
      // wait for completion of previous sends
      if (phtn_send_buffer[i_b].sent())
        MPI_Wait(&phtn_send_request[i_b], MPI_STATUS_IGNORE);
      // send one photon vector to finish off receives, these photons will not
      // be processed by the receiving ranks (all ranks are out of transport)
      MPI_Isend(&one_photon[0], 1, MPI_Particle, adj_rank, photon_tag,
                MPI_COMM_WORLD, &phtn_send_request[i_b]);
      mctr.n_sends_posted++;
      MPI_Wait(&phtn_send_request[i_b], MPI_STATUS_IGNORE);
      mctr.n_sends_completed++;
    } // end loop over adjacent processors
  }

  // wait for receive requests
  for (uint32_t i_b = 0; i_b < n_adjacent; ++i_b) {
    MPI_Wait(&phtn_recv_request[i_b], MPI_STATUS_IGNORE);
    mctr.n_receives_completed++;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());

  // all ranks have now finished transport
  delete[] phtn_recv_request;
  delete[] phtn_send_request;

  // set diagnostic quantities
  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_census_size(census_list.size());
  imc_state->set_network_message_counts(mctr);
  imc_state->set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));
  imc_state->set_rank_mpi_time(t_mpi.get_time("timestep mpi"));

  return census_list;
}

#endif // def transport_particle_pass_h_
//---------------------------------------------------------------------------//
// end of transport_particle_pass.h
//---------------------------------------------------------------------------//
