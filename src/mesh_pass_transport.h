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

//! Transport a single photon until it has a terminating event (kill, exit,
// wait for data, census)
Constants::event_type
transport_photon_mesh_pass(Photon &phtn, Mesh *mesh, RNG *rng, double &next_dt,
                           double &exit_E, double &census_E,
                           std::vector<double> &rank_abs_E,
                           std::vector<double> &rank_track_E,
                           std::unordered_map<uint32_t, double> &off_rank_abs_E)

{
  using Constants::VACUUM;
  using Constants::REFLECT;
  using Constants::ELEMENT;
  using Constants::PROCESSOR;
  // events
  using Constants::WAIT;
  using Constants::CENSUS;
  using Constants::KILL;
  using Constants::EXIT;
  using Constants::bc_type;
  using Constants::event_type;
  using Constants::c;
  using std::min;

  uint32_t cell_id, next_cell, next_grip;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;
  Cell cell;

  uint32_t surface_cross = 0;
  uint32_t on_rank_start = mesh->get_offset();
  const double cutoff_fraction = 0.01; // note: get this from IMC_state

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

    // process on rank tallies as usual
    if (mesh->on_processor(cell_id)) {
      rank_track_E[cell_id - on_rank_start] += absorbed_E / (sigma_a * f);
      rank_abs_E[cell_id - on_rank_start] += absorbed_E;
    } else
      off_rank_abs_E[cell_id] += absorbed_E;

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      if (mesh->on_processor(cell_id))
        rank_abs_E[cell_id - on_rank_start] += phtn.get_E();
      else
        off_rank_abs_E[cell_id] += phtn.get_E();
      active = false;
      event = KILL;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
        // if effective scatter change frequency
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell.get_bc(surface_cross);
        if (boundary_event == ELEMENT || boundary_event == PROCESSOR) {
          next_cell = cell.get_next_cell(surface_cross);
          next_grip = cell.get_next_grip(surface_cross);
          phtn.set_cell(next_cell);
          phtn.set_grip(next_grip);
          cell_id = next_cell;
          // look for this cell, if it's not there transport later
          if (mesh->mesh_available(cell_id))
            cell = mesh->get_on_rank_cell(cell_id);
          else {
            event = WAIT;
            active = false;
          }
        } else if (boundary_event == VACUUM) {
          active = false;
          exit_E += phtn.get_E();
          event = EXIT;
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        phtn.set_distance_to_census(c * next_dt);
        active = false;
        census_E += phtn.get_E();
        event = CENSUS;
      }
    } // end event loop
  }   // end while alive
  return event;
}

//! Transport photons from a source object using the mesh-passing algorithm
// and two-sided messaging to fulfill requests for mesh data
std::vector<Photon> mesh_pass_transport(
    Source &source, Mesh *mesh, IMC_State *imc_state,
    IMC_Parameters *imc_parameters, Mesh_Request_Manager &req_manager,
    Tally_Manager &tally_manager, Message_Counter &mctr,
    std::vector<double> &rank_abs_E, std::vector<double> &rank_track_E,
    MPI_Types *mpi_types, const Info &mpi_info) {
  using std::queue;
  using std::vector;
  using Constants::event_type;
  // events
  using Constants::WAIT;
  using Constants::CENSUS;
  using Constants::KILL;
  using Constants::EXIT;

  uint32_t n_local = source.get_n_photon();
  uint32_t n_local_sourced = 0;

  uint32_t cell_id;

  double census_E = 0.0;
  double exit_E = 0.0;
  double dt = imc_state->get_next_dt();      //! For making current photons
  double next_dt = imc_state->get_next_dt(); //! For census photons

  RNG *rng = imc_state->get_rng();
  Photon phtn;

  // timing
  Timer t_transport;
  Timer t_mpi;
  Timer t_rebalance_census;
  t_transport.start_timer("timestep transport");

  // New data flag is initially false
  bool new_data = false;
  std::vector<Cell> new_cells; // New cells from completed RMA requests

  // Number of particles to run between MPI communication
  const uint32_t batch_size = imc_parameters->get_batch_size();

  event_type event;
  uint32_t wait_list_size;

  // for tallying off rank data
  std::unordered_map<uint32_t, double> off_rank_abs_E;

  //--------------------------------------------------------------------------//
  // main loop over photons
  //--------------------------------------------------------------------------//
  vector<Photon> census_list; //! Local end of timestep census list
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
      if (mesh->mesh_available(cell_id)) {
        event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                           census_E, rank_abs_E, rank_track_E,
                                           off_rank_abs_E);
        cell_id = phtn.get_cell();
      } else
        event = WAIT;

      if (event == CENSUS) {
        census_list.push_back(phtn);
      } else if (event == WAIT) {
        t_mpi.start_timer("timestep mpi");
        req_manager.request_cell(phtn.get_grip(), mctr);
        t_mpi.stop_timer("timestep mpi");
        wait_list.push(phtn);
      }
      n--;
    } // end batch transport

    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    t_mpi.start_timer("timestep mpi");
    new_cells = req_manager.process_mesh_requests(mctr);
    new_data = !new_cells.empty();
    if (new_data)
      mesh->add_non_local_mesh_cells(new_cells);
    t_mpi.stop_timer("timestep mpi");
    // if data was received, try to transport photons on waiting list
    if (new_data) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS) {
          census_list.push_back(phtn);
        } else if (event == WAIT) {
          t_mpi.start_timer("timestep mpi");
          req_manager.request_cell(phtn.get_grip(), mctr);
          t_mpi.stop_timer("timestep mpi");
          wait_list.push(phtn);
        }
      } // end wp in wait_list
    }
  } // end while (n_local_source < n_local)

  //--------------------------------------------------------------------------//
  // Main transport loop finished, transport photons waiting for data
  //--------------------------------------------------------------------------//
  while (!wait_list.empty()) {

    t_mpi.start_timer("timestep mpi");
    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    new_cells = req_manager.process_mesh_requests(mctr);
    new_data = !new_cells.empty();
    if (new_data)
      mesh->add_non_local_mesh_cells(new_cells);

    t_mpi.stop_timer("timestep mpi");

    // if new data received, transport waiting list
    if (new_data || req_manager.no_active_requests()) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS) {
          census_list.push_back(phtn);
        } else if (event == WAIT) {
          t_mpi.start_timer("timestep mpi");
          req_manager.request_cell(phtn.get_grip(), mctr);
          t_mpi.stop_timer("timestep mpi");
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

  // set the preffered census size to 10% of the user photon number and comb
  uint64_t max_census_photons = 0.1 * imc_parameters->get_n_user_photon();
  comb_photons(census_list, max_census_photons, rng);

  // all ranks have now finished transport, set diagnostic quantities
  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_network_message_counts(mctr);
  imc_state->set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  // send the off-rank census back to ranks that own the mesh its on and receive
  // census particles that are on your mesh

  t_rebalance_census.start_timer("timestep rebalance_census");
  vector<Photon> rebalanced_census =
      rebalance_raw_census(census_list, mesh, mpi_types);
  t_rebalance_census.stop_timer("timestep rebalance_census");

  imc_state->set_rank_mpi_time(t_mpi.get_time("timestep mpi"));
  imc_state->set_rank_rebalance_time(
      t_rebalance_census.get_time("timestep rebalance_census"));

  // use this only if the off rank census is separate
  // census_list.insert(census_list.end(), rebalanced_census.begin(),
  //  rebalanced_census.end());

  // sort on census vectors by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // set post census size after sorting and merging
  imc_state->set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_mesh_pass_h_
//----------------------------------------------------------------------------//
// end of transport_mesh_pass.h
//----------------------------------------------------------------------------//
