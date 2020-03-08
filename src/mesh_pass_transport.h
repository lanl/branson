//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_pass_transport.h
 * \author Alex Long
 * \date   June 6 2015
 * \brief  Transport routine using two sided messaging and mesh-passing DD
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_mesh_pass_h_
#define transport_mesh_pass_h_

#include <algorithm>
#include <mpi.h>
#include <numeric>
#include <queue>
#include <vector>

#include "RNG.h"
#include "comb_photons.h"
#include "constants.h"
#include "decompose_photons.h"
#include "info.h"
#include "mesh.h"
#include "mesh_request_manager.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "sampling_functions.h"
#include "source.h"
#include "tally_manager_rma.h"
#include "timer.h"
#include "transport_photon.h"

//! Transport photons from a source object using the mesh-passing algorithm
// and two-sided messaging to fulfill requests for mesh data
std::vector<Photon> mesh_pass_transport(
    Source &source, Mesh &mesh, IMC_State &imc_state,
    const IMC_Parameters &imc_parameters, Mesh_Request_Manager &req_manager,
    Tally_Manager &tally_manager, Message_Counter &mctr,
    std::vector<double> &rank_abs_E, std::vector<double> &rank_track_E,
    const MPI_Types &mpi_types, const Info &mpi_info) {
  using Constants::event_type;
  using std::queue;
  using std::vector;
  // events
  using Constants::CENSUS;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::WAIT;

  uint32_t n_local = source.get_n_photon();
  uint32_t n_local_sourced = 0;

  uint32_t cell_id;

  double census_E = 0.0;
  double exit_E = 0.0;
  double dt = imc_state.get_next_dt();      //! For making current photons
  double next_dt = imc_state.get_next_dt(); //! For census photons

  RNG *rng = imc_state.get_rng();
  Photon phtn;

  // timing
  Timer t_transport;
  Timer t_rebalance_census;
  t_transport.start_timer("timestep transport");

  // Number of particles to run between MPI communication
  const uint32_t batch_size = imc_parameters.get_batch_size();

  event_type event;
  uint32_t wait_list_size;

  // for tallying off rank data
  std::unordered_map<uint32_t, double> off_rank_abs_E;

  //--------------------------------------------------------------------------//
  // main loop over photons
  //--------------------------------------------------------------------------//
  // size the census list to the maximum size
  vector<Photon> census_list;
  census_list.resize(source.get_n_photon());
  uint64_t census_list_size = 0;
  // vector<Photon> off_rank_census_list; //! Off rank end of timestep census
  // list
  queue<Photon> wait_list; //! Photons waiting for mesh data
  while (n_local_sourced < n_local) {

    uint32_t n = batch_size;

    while (n && n_local_sourced < n_local) {

      phtn = source.get_photon(rng, dt);
      n_local_sourced++;

      // get start cell, this only changea with cell crossing event
      cell_id = phtn.get_cell();

      // if mesh available, transport and process, otherwise put on the
      // waiting list
      if (mesh.mesh_available(cell_id)) {
        event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                           census_E, rank_abs_E, rank_track_E,
                                           off_rank_abs_E);
        cell_id = phtn.get_cell();
      } else
        event = WAIT;

      if (event == CENSUS) {
        census_list[census_list_size++] = phtn;
      } else if (event == WAIT) {
        req_manager.request_cell(phtn.get_grip(), mctr);
        wait_list.push(phtn);
      }
      n--;
    } // end batch transport

    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    req_manager.process_mesh_requests(mctr);
    if (req_manager.get_n_new_cells() > 0)
      mesh.add_non_local_mesh_cells(req_manager.get_receive_buffers(),
                                    req_manager.get_n_new_cells());
    // if data was received, try to transport photons on waiting list
    if (req_manager.get_n_new_cells()) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh.mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS) {
          census_list[census_list_size++] = phtn;
        } else if (event == WAIT) {
          req_manager.request_cell(phtn.get_grip(), mctr);
          wait_list.push(phtn);
        }
      } // end wp in wait_list
    }
  } // end while (n_local_source < n_local)

  //--------------------------------------------------------------------------//
  // Main transport loop finished, transport photons waiting for data
  //--------------------------------------------------------------------------//
  while (!wait_list.empty()) {

    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    req_manager.process_mesh_requests(mctr);
    if (req_manager.get_n_new_cells() > 0)
      mesh.add_non_local_mesh_cells(req_manager.get_receive_buffers(),
                                    req_manager.get_n_new_cells());

    // if new data received, transport waiting list
    if (req_manager.get_n_new_cells() || req_manager.no_active_requests()) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh.mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS) {
          census_list[census_list_size++] = phtn;
        } else if (event == WAIT) {
          req_manager.request_cell(phtn.get_grip(), mctr);
          wait_list.push(phtn);
        }
      }
    }
  } // end while wait_list not empty

  // process off rank tally data and force send (force send requires completion
  // of all tallies before continuing
  bool force_send = true;
  tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // start non-blocking allreduce, when it's finished all ranks are done
  MPI_Request completion_request;
  int send_int = 1;
  int recv_int;
  MPI_Iallreduce(&send_int, &recv_int, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD,
                 &completion_request);
  int finished = false;

  //--------------------------------------------------------------------------//
  // While waiting for other ranks to finish, check for other messages
  //--------------------------------------------------------------------------//
  while (!finished) {
    req_manager.process_mesh_requests(mctr);
    MPI_Test(&completion_request, &finished, MPI_STATUS_IGNORE);
  } // end while

  // wait for all ranks to finish transport to finish off cell and cell id
  // requests and sends
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_info.get_rank() == 0)
    std::cout << "Transport complete" << std::endl;

  // append remote tally to the current tally
  tally_manager.add_remote_tally(rank_abs_E);

  // trim census back to actual size
  census_list.erase(census_list.begin() + census_list_size, census_list.end());

  // set the preffered census size to 10% of the user photon number and comb
  uint64_t max_census_photons = 0.1 * imc_parameters.get_n_user_photon();
  comb_photons(census_list, max_census_photons, rng);

  // all ranks have now finished transport, set diagnostic quantities
  imc_state.set_exit_E(exit_E);
  imc_state.set_post_census_E(census_E);
  imc_state.set_network_message_counts(mctr);
  imc_state.set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  // send the off-rank census back to ranks that own the mesh its on and receive
  // census particles that are on your mesh

  t_rebalance_census.start_timer("timestep rebalance_census");
  vector<Photon> rebalanced_census =
      rebalance_raw_census(census_list, mesh, mpi_types);
  t_rebalance_census.stop_timer("timestep rebalance_census");

  imc_state.set_rank_rebalance_time(
      t_rebalance_census.get_time("timestep rebalance_census"));

  // use this only if the off rank census is separate
  // census_list.insert(census_list.end(), rebalanced_census.begin(),
  //  rebalanced_census.end());

  // sort on census vectors by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // set post census size after sorting and merging
  imc_state.set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_mesh_pass_h_
//----------------------------------------------------------------------------//
// end of transport_mesh_pass.h
//----------------------------------------------------------------------------//
