//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rma_mesh_pass_driver.h
 * \author Alex Long
 * \date   March 3 2017
 * \brief  Functions to run IMC with one-sided mesh passing
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
#ifndef rma_mesh_pass_driver_h_
#define rma_mesh_pass_driver_h_

#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "census_creation.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "info.h"
#include "load_balance.h"
#include "mesh.h"
#include "mesh_rma_manager.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "pretransport_requests.h"
#include "rma_mesh_pass_transport.h"
#include "source.h"
#include "timer.h"
#include "write_silo.h"

void imc_rma_mesh_pass_driver(Mesh *mesh, IMC_State *imc_state,
                              IMC_Parameters *imc_parameters,
                              MPI_Types *mpi_types, const Info &mpi_info) {
  using std::vector;
  vector<double> abs_E(mesh->get_n_local_cells(), 0.0);
  vector<double> track_E(mesh->get_n_local_cells(), 0.0);
  vector<Photon> census_photons;
  vector<uint32_t> needed_grip_ids; //! Grips needed after load balance
  Message_Counter mctr;
  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // make object that handles RMA mesh requests and start access
  RMA_Manager rma_manager(mesh->get_off_rank_bounds(),
                          mesh->get_max_grip_size(), mpi_types,
                          mesh->get_mesh_window_ref());
  rma_manager.start_access();

  // make object that handles tally data
  Tally_Manager tally_manager(rank, mesh->get_off_rank_bounds(),
                              mesh->get_n_local_cells());

  while (!imc_state->finished()) {
    if (rank == 0)
      imc_state->print_timestep_header();

    mctr.reset_counters();

    // set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    // all reduce to get total source energy to make correct number of
    // particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    // this will be zero on first time step, source construction
    // handles initial census
    imc_state->set_pre_census_E(get_photon_list_E(census_photons));

    // setup source and load balance, time load balance
    Source source(mesh, imc_state, imc_parameters->get_n_user_photon(),
                  global_source_energy, census_photons);
    Timer t_lb;
    t_lb.start_timer("load balance");
    load_balance(source.get_work_vector(), census_photons,
                 source.get_n_photon(), mpi_types, mpi_info);

    // get new particle count after load balance, group particle work by cell
    source.post_lb_prepare_source();

    // prerequest data for work packets and census particles not on your rank
    pretransport_requests(source.get_work_vector(), census_photons, mesh,
                          rma_manager, mctr);

    // cell properties are set in calculate_photon_energy--make sure
    // everybody gets here together so that windows are not changing
    // when transport starts
    MPI_Barrier(MPI_COMM_WORLD);

    vector<Cell> new_cells = rma_manager.process_rma_mesh_requests(mctr);
    if (!new_cells.empty())
      mesh->add_non_local_mesh_cells(new_cells);

    t_lb.stop_timer("load balance");
    imc_state->set_load_balance_time(t_lb.get_time("load balance"));

    // transport photons
    census_photons = rma_mesh_pass_transport(
        source, mesh, imc_state, imc_parameters, rma_manager, tally_manager,
        mctr, abs_E, track_E, mpi_types, mpi_info);

    imc_state->set_transported_particles(source.get_n_photon());

    // cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, track_E, imc_state);

    imc_state->print_conservation(imc_parameters->get_dd_mode());

    // purge the working mesh, it will be updated by other ranks and is now
    // invalid
    mesh->purge_working_mesh();

    if (imc_parameters->get_write_silo_flag()) {
      // write SILO file
      // don't plot the n_requests vector
      double fake_mpi_runtime = 0.0;
      vector<uint32_t> n_requests(mesh->get_n_local_cells(), 0);
      write_silo(mesh, imc_state->get_time(), imc_state->get_step(),
                 imc_state->get_rank_transport_runtime(),
                 fake_mpi_runtime, rank, n_rank, n_requests);
    }
    // reset rma_manager object for next timestep
    rma_manager.end_timestep();
    tally_manager.end_timestep();

    // update time for next step
    imc_state->next_time_step();
  }

  // close access to MPI windows in RMA_Manger object
  rma_manager.end_access();
}

#endif // rma_mesh_pass_driver_h_

//---------------------------------------------------------------------------//
// end of rma_mesh_pass_driver.h
//---------------------------------------------------------------------------//
