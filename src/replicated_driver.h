//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   replicated_driver.h
 * \author Alex Long
 * \date   March 3 2017
 * \brief  Functions to run IMC with a replicated domain
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef replicated_driver_h_
#define replicated_driver_h_

#include <iostream>
#include <mpi.h>
#include <functional>
#include <vector>

#include "census_creation.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "info.h"
#include "imc_state.h"
#include "imc_parameters.h"
#include "mesh.h"
#include "source.h"
#include "timer.h"
#include "replicated_transport.h"
#include "write_silo.h"

void imc_replicated_driver(Mesh *mesh,
                              IMC_State *imc_state,
                              IMC_Parameters *imc_parameters,
                              MPI_Types * mpi_types,
                              const Info &mpi_info)
{
  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;
  Message_Counter mctr;
  int rank = mpi_info.get_rank();

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    mctr.reset_counters();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

    imc_state->set_pre_census_E(get_photon_list_E(census_photons));

    // setup source
    Source source(mesh, imc_state, imc_parameters->get_n_user_photon(),
      global_source_energy, census_photons);
    // no load balancing in particle passing method, just call the method
    // to get accurate count and map census to work correctly
    source.post_lb_prepare_source();

    imc_state->set_transported_particles(source.get_n_photon());

    census_photons = replicated_transport(source, mesh, imc_state,
      imc_parameters, mpi_types, mctr, abs_E, mpi_info);

    // using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI_Allreduce(MPI_IN_PLACE, &abs_E[0], mesh->get_global_num_cells(),
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mesh->update_temperature(abs_E, imc_state);

    // for replicated, just let root do conservation
    if (rank) {
      imc_state->set_absorbed_E(0.0);
      imc_state->set_pre_mat_E(0.0);
      imc_state->set_post_mat_E(0.0);
    }

    imc_state->print_conservation(imc_parameters->get_dd_mode());
    // update time for next step
    imc_state->next_time_step();
  }
}

#endif // replicated_driver_h_

//---------------------------------------------------------------------------//
// end of replicated_driver.h
//---------------------------------------------------------------------------//
