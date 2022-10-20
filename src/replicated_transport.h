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
#include <functional>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "RNG.h"
#include "constants.h"
#include "gpu_setup.h"
#include "info.h"
#include "mesh.h"
#include "message_counter.h"
#include "particle_pass_transport.h"
#include "photon.h"
#include "sampling_functions.h"

std::vector<Photon> replicated_transport(
    const Mesh &mesh, const GPU_Setup &gpu_setup, IMC_State &imc_state, std::vector<double> &rank_abs_E, std::vector<double> &rank_track_E, std::vector<Photon> all_photons) {
  using std::cout;
  using std::endl;
  using std::vector;

  double census_E = 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state.get_next_dt(); //! Set for census photons
  double dt = imc_state.get_next_dt();      //! For making current photons

  RNG *rng = imc_state.get_rng();

  // timing
  Timer t_transport;
  t_transport.start_timer("timestep transport");

  // replicated transport does not require the global photon count
  uint64_t n_local = all_photons.size();

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  vector<Photon> census_list;   //! End of timestep census list

  vector<Cell_Tally> cell_tallies(mesh.get_n_cells());;
  uint32_t rank_cell_offset{0}; // no offset in replicated mesh
  if(gpu_setup.use_gpu_transporter())
    gpu_transport(rank_cell_offset, all_photons, gpu_setup.get_device_cells_ptr(), rng, cell_tallies);
  else
    cpu_transport(rank_cell_offset, all_photons, mesh.get_cells() rng, cell_tallies);

  // post process photons, account for escaped energy and add particles to census
  post_process_photons(all_photons, census_list, census_E, exit_E);

  // copy cell tallies back out to rank_abs_E and rank_track_E
  for (size_t i = 0; i<cell_tallies.size();++i) {
    rank_abs_E[i] = cell_tallies[i].get_abs_E();
    rank_track_E[i] = cell_tallies[i].get_track_E();
  }

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());

  // set diagnostic quantities
  imc_state.set_exit_E(exit_E);
  imc_state.set_post_census_E(census_E);
  imc_state.set_census_size(census_list.size());
  imc_state.set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  return census_list;
}

#endif // def transport_replicated_h_
//---------------------------------------------------------------------------//
// end of transport_replicated.h
//---------------------------------------------------------------------------//
