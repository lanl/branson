/*
  Author: Alex Long
  Date: 6/29/2014
  Name: transport_mesh_pass.h
*/

#ifndef transport_rma_mesh_pass_h_
#define transport_rma_mesh_pass_h_

#include <algorithm>
#include <vector>
#include <numeric>
#include <queue>
#include <mpi.h>

#include "constants.h"
#include "decompose_photons.h"
#include "mesh.h"
#include "mesh_rma_manager.h"
#include "RNG.h"
#include "sampling_functions.h"
#include "transport_mesh_pass.h"

std::vector<Photon> transport_rma_mesh_pass(Source& source,
                                            Mesh* mesh,
                                            IMC_State* imc_state,
                                            IMC_Parameters* imc_parameters,
                                            RMA_Manager* rma_manager,
                                            std::vector<double>& rank_abs_E)
{
  using Constants::finish_tag;
  using std::queue;
  using std::vector;
  using Constants::proc_null;
  using Constants::event_type;
  //events
  using Constants::WAIT; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;

  uint32_t n_local = source.get_n_photon();
  uint32_t n_local_sourced = 0;

  uint32_t cell_id;
  double census_E = 0.0;
  double exit_E = 0.0;
  double dt = imc_state->get_next_dt(); //<! For making current photons
  double next_dt = imc_state->get_next_dt(); //<! For census photons

  RNG *rng = imc_state->get_rng();
  Photon phtn;

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // parallel event counters
  uint32_t n_cell_messages=0; //! Number of cell messages
  uint32_t n_cells_sent=0; //! Number of cells passed
  uint32_t n_sends_posted=0; //! Number of sent messages posted
  uint32_t n_sends_completed=0; //! Number of sent messages completed
  uint32_t n_receives_posted=0; //! Number of received messages completed
  uint32_t n_receives_completed=0; //! Number of received messages completed

  
  bool new_data = false; //! New data flag is initially false
  std::vector<Cell> new_cells; // New cells from completed RMA requests
  // Number of particles to run between MPI communication 
  const uint32_t batch_size = imc_parameters->get_batch_size();

  event_type event;
  uint32_t wait_list_size;

  ////////////////////////////////////////////////////////////////////////
  // main loop over photons
  ////////////////////////////////////////////////////////////////////////
  vector<Photon> census_list; //!< Local end of timestep census list
  vector<Photon> off_rank_census_list; //!< Off rank end of timestep census list
  queue<Photon> wait_list; //!< Photons waiting for mesh data 
  while ( n_local_sourced < n_local) {
    
    uint32_t n = batch_size;

    while (n && n_local_sourced < n_local) {

      phtn =source.get_photon(rng, dt); 
      n_local_sourced++;

      //get start cell, this only change with cell crossing event
      cell_id=phtn.get_cell();

      // if mesh available, transport and process, otherwise put on the
      // waiting list
      if (mesh->mesh_available(cell_id)) {
        event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                        census_E, rank_abs_E);
        cell_id = phtn.get_cell();
      }
      else event = WAIT;

      if (event==CENSUS) { 
        if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
        else off_rank_census_list.push_back(phtn);
      }
      else if (event==WAIT) {
        rma_manager->request_cell_rma(cell_id, n_receives_posted);
        wait_list.push(phtn);
      }
      n--;
    } // end batch transport

    //process mesh requests
    new_cells = rma_manager->process_rma_mesh_requests(n_receives_completed);
    new_data = !new_cells.empty();
    if (new_data) mesh->add_non_local_mesh_cells(new_cells);
    // if data was received, try to transport photons on waiting list
    if (new_data) {
      wait_list_size = wait_list.size();
      for (uint32_t wp =0; wp<wait_list_size; wp++) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id=phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                          census_E, rank_abs_E);
          cell_id = phtn.get_cell();
          if (event==CENSUS) { 
            if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
            else off_rank_census_list.push_back(phtn);
          }
          else if (event==WAIT) {
            rma_manager->request_cell_rma(cell_id, n_receives_posted);
            wait_list.push(phtn);
          }
        }
        else {
          rma_manager->request_cell_rma(cell_id, n_receives_posted);
          wait_list.push(phtn);
        }
      } // end wp in wait_list
    }

  } //end while (n_local_source < n_local)

  ////////////////////////////////////////////////////////////////////////
  // Main transport loop finished, transport photons waiting for data
  ////////////////////////////////////////////////////////////////////////
  wait_list_size = wait_list.size();
  while (!wait_list.empty()) {
    new_cells = rma_manager->process_rma_mesh_requests(n_receives_completed);
    new_data = !new_cells.empty();
    if (new_data) mesh->add_non_local_mesh_cells(new_cells);
    // if new data received or there are no active mesh requests, try to 
    // transport waiting list (it could be that there are no active memory
    // requests because the request queue was full at the time)
    if (new_data || rma_manager->no_active_requests()) {
      wait_list_size = wait_list.size();
      for (uint32_t wp =0; wp<wait_list_size; wp++) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id=phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                              census_E, rank_abs_E);
          cell_id = phtn.get_cell();
          if (event==CENSUS) { 
            if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
            else off_rank_census_list.push_back(phtn);
          }
          else if (event==WAIT) {
            rma_manager->request_cell_rma(cell_id, n_receives_posted);
            wait_list.push(phtn);
          }
        }
        else {
          wait_list.push(phtn);
          rma_manager->request_cell_rma(cell_id, n_receives_posted);
        }
      }
    }
  } //end while wait_list not empty

  MPI_Barrier(MPI_COMM_WORLD);

  //All ranks have now finished transport

  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_step_cell_messages(n_receives_posted);
  imc_state->set_step_cells_sent(n_receives_completed);
  imc_state->set_step_sends_posted(n_sends_posted);
  imc_state->set_step_sends_completed(n_sends_completed);
  imc_state->set_step_receives_posted(n_receives_posted);
  imc_state->set_step_receives_completed(n_receives_completed);

  //send the off-rank census back to ranks that own the mesh its on.
  //receive census particles that are on your mesh
  vector<Photon> rebalanced_census = rebalance_census(off_rank_census_list,
                                                      mesh);
  census_list.insert(census_list.end(), 
    rebalanced_census.begin(), 
    rebalanced_census.end());

  //sort on census vectors by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // set post census size after sorting and merging
  imc_state->set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_rma_mesh_pass_h_
