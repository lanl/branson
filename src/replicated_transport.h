//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   replicated_transport.h
 * \author Alex Long
 * \date   March 3 2017
 * \brief  IMC transport on a replicated domain
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_replicated_h_
#define transport_replicated_h_

#include <algorithm>
#include <mpi.h>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "particle_pass_transport.h"
#include "photon.h"
#include "sampling_functions.h"
#include "RNG.h"

std::vector<Photon> replicated_transport( Source& source,
                                          Mesh* mesh,
                                          IMC_State* imc_state,
                                          IMC_Parameters* imc_parameters,
                                          MPI_Types* mpi_types,
                                          Message_Counter& mctr,
                                          std::vector<double>& rank_abs_E,
                                          const Info& mpi_info)
{
  using Constants::event_type;
  using Constants::PASS; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;
  using Constants::WAIT;
  using std::vector;
  using std::cout;
  using std::endl;

  double census_E= 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt(); //! Set for census photons
  double dt = imc_state->get_next_dt(); //! For making current photons

  RNG *rng = imc_state->get_rng();

  // timing
  Timer t_transport;
  Timer t_mpi;
  t_transport.start_timer("timestep transport");

  // get global photon count
  uint64_t n_local = source.get_n_photon();
  uint64_t n_global;

  MPI_Allreduce(&n_local, &n_global, 1, MPI_UNSIGNED_LONG, MPI_SUM,
    MPI_COMM_WORLD);

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  vector<Photon> census_list; //! End of timestep census list

  uint64_t n_complete = 0; //! Completed histories, regardless of origin
  uint64_t n_local_sourced = 0; //! Photons pulled from source object
  Photon phtn;
  event_type event;

  //------------------------------------------------------------------------//
  // Transport photons from source
  //------------------------------------------------------------------------//
  while (n_local_sourced < n_local) {

    phtn =source.get_photon(rng, dt);
    n_local_sourced++;

    event = transport_photon_particle_pass(phtn, mesh, rng, next_dt, exit_E,
                                          census_E, rank_abs_E);
    switch(event) {
      // this case should never be reached
      case WAIT:
        break;
      // this case should never be reached
      case PASS:
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
    }
  } // end while

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());

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

#endif // def transport_replicated_h_
//---------------------------------------------------------------------------//
// end of transport_replicated.h
//---------------------------------------------------------------------------//
