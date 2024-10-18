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
#include "transport_photon.h"
#include "history_based_transport.h"
#include "event_based_transport.h"
#include "photon.h"

template <typename Census_T>
Census_T replicated_transport(
    const Mesh &mesh, const GPU_Setup &gpu_setup, IMC_State &imc_state,
    std::vector<double> &rank_abs_E, std::vector<double> &rank_track_E, 
    Census_T  &all_photons, const int n_omp_threads, const uint32_t batch_size,
    const int transport_algorithm) {
  using std::cout;
  using std::endl;
  using std::vector;

  // is the GPU even available?
  #ifdef USE_CUDA
  constexpr bool gpu_available = true;
  #else
  constexpr bool gpu_available = false;
  #endif

  double census_E = 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state.get_next_dt(); //! Set for census photons
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // print warning message if GPU transport is requested but not available
  if(rank==0 && gpu_setup.use_gpu_transporter() && !gpu_available) {
    std::cout<<"WARNING: use_gpu_transporter set to true but GPU kernel not available,";
    std::cout<<" running transport on CPU"<<std::endl;
  }

  // timing
  Timer t_transport;
  t_transport.start_timer("timestep transport");

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  Census_T census_list;   //! End of timestep census list
  vector<Cell_Tally> cell_tallies(mesh.get_n_local_cells());
  uint32_t rank_cell_offset{0}; // no offset in replicated mesh

  // GPU mode
  if(gpu_setup.use_gpu_transporter() && gpu_available ) {
    if constexpr(std::is_same_v<Census_T, std::vector<Photon>>) {
      t_transport.start_timer("gpu transport");
      gpu_transport_photons(rank_cell_offset, all_photons, gpu_setup.get_device_cells_ptr(), cell_tallies);
      t_transport.stop_timer("gpu transport");
      std::cout<<"gpu transport time: "<<t_transport.get_time("gpu transport")<<std::endl;
    }
    else {
      std::cout<<"No GPU kernel for SOA particle data structures yet"<<std::endl; 
      exit(EXIT_FAILURE);
    }
  } // CPU mode
  else {
    if (transport_algorithm == Constants::HISTORY) { 
      std::cout<<"about to start history based"<<std::endl;
      history_cpu_transport_photons(rank_cell_offset, all_photons, mesh.get_cells(), cell_tallies, n_omp_threads);
      auto batch_complete = post_process_photons(next_dt, all_photons, census_list, census_E, exit_E);
    }
    else if(transport_algorithm ==Constants::EVENT) {
      if constexpr(std::is_same_v<Census_T, std::vector<Photon>>) { 
        std::vector<EmissionGroupData> emission_groups(mesh.get_n_local_cells());
        for (size_t i = 0; i < mesh.get_n_local_cells(); ++i) {
          emission_groups[i] = precompute_emission_group_data(mesh.get_cells()[i]);
        }
        // ARL: try batching this?
        std::cout<<"about to start aos event- based"<<std::endl;
        cpu_event_transport_photons(rank_cell_offset, all_photons, mesh.get_cells(), cell_tallies, n_omp_threads, emission_groups);
      }
      else {
        std::vector<EmissionGroupData> emission_groups(mesh.get_n_local_cells());
        for (size_t i = 0; i < mesh.get_n_local_cells(); ++i) {
          emission_groups[i] = precompute_emission_group_data(mesh.get_cells()[i]);
        }

        for (size_t batch_start = 0; batch_start < all_photons.size(); batch_start += batch_size) {
          size_t batch_end = std::min(batch_start + batch_size, all_photons.size());

          std::cout<<"about to start soa event-based batch"<<batch_start<<" "<<batch_end<<std::endl;
          auto batch_photons = all_photons.get_sub_batch(batch_start, batch_end);
          cpu_event_transport_photons(rank_cell_offset, batch_photons, mesh.get_cells(), cell_tallies, n_omp_threads, emission_groups);

          // post process photons, account for escaped energy and add particles to census
          auto batch_complete = post_process_photons(next_dt, batch_photons, census_list, census_E, exit_E);
        }
      }
    }
    else {
      std::cout<<"Transport algorithm type not supported... exiting"<<std::endl;
    }
  }

  // copy cell tallies back out to rank_abs_E and rank_track_E
  double total_abs = 0;
  for (size_t i = 0; i<cell_tallies.size();++i) {
    total_abs+=cell_tallies[i].get_abs_E();
    rank_abs_E[i] = cell_tallies[i].get_abs_E();
    rank_track_E[i] = cell_tallies[i].get_track_E();
  }

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);

  //std::sort(census_list.begin(), census_list.end());

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
