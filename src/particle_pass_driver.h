//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   particle_pass_driver.h
 * \author Alex Long
 * \date   March 3 2017
 * \brief  Functions to run IMC with particle passing
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef particle_pass_driver_h_
#define particle_pass_driver_h_

#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "census_creation.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "info.h"
#include "mesh.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "particle_pass_transport.h"
#include "source.h"
#include "timer.h"
#include "write_silo.h"

void imc_particle_pass_driver(Mesh &mesh, IMC_State &imc_state,
                              const IMC_Parameters &imc_parameters,
                              const MPI_Types &mpi_types,
                              const Info &mpi_info) {
  using std::vector;
  vector<double> abs_E(mesh.get_n_local_cells(), 0.0);
  vector<double> track_E(mesh.get_n_local_cells(), 0.0);
  vector<Photon> census_photons;
  auto n_user_photons = imc_parameters.get_n_user_photons();
  Message_Counter mctr;
  const int rank = mpi_info.get_rank();
  const int n_ranks = mpi_info.get_n_rank();

  const uint32_t seed = imc_parameters.get_rng_seed();
  
#ifdef USE_MEMORY_RECORD
  MemoryRecorder mem_record = MemoryRecorder();
#endif
  std::string start_string = "StartIter_";
  std::string end_string = "EndIter_";
  std::string meminfo_string;

  while (!imc_state.finished()) {
    if (rank == 0)
      imc_state.print_timestep_header();

    meminfo_string = start_string + std::to_string(imc_state.get_step());

#ifdef USE_MEMORY_RECORD
    mem_record.read_meminfo(meminfo_string);
#endif

    mctr.reset_counters();

    //set opacity, Fleck factor, all energy to source
    wrapped_cali_mark_begin("calculate photon energy");
    mesh.calculate_photon_energy(imc_state, n_user_photons);
    wrapped_cali_mark_end("calculate photon energy");

    // all reduce to get total source energy to make correct number of
    // particles on each rank
    double global_source_energy = mesh.get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    wrapped_cali_mark_begin("calculate pre census energy");
    imc_state.set_pre_census_E(get_photon_list_E(census_photons));
    wrapped_cali_mark_end("calculate pre census energy");

    // make gpu setup object, may want to source on GPU later so make it before sourcing here
    GPU_Setup gpu_setup(rank, n_ranks, imc_parameters.get_use_gpu_transporter_flag(), mesh.get_cells());

    // setup source
    if (imc_state.get_step() == 1) {
      wrapped_cali_mark_begin("make initial census photons");
      census_photons = make_initial_census_photons(imc_state.get_dt(), mesh, rank, seed, n_user_photons, global_source_energy);
      wrapped_cali_mark_end("make initial census photons");
    }

    imc_state.set_pre_census_E(get_photon_list_E(census_photons));
    MPI_Barrier(MPI_COMM_WORLD);
    // make emission and source photons
    wrapped_cali_mark_begin("make emission and source photons");
    auto all_photons = make_photons(imc_state.get_dt(), mesh, rank, imc_state.get_step(), seed, n_user_photons, global_source_energy);
    wrapped_cali_mark_end("make emission and source photons");
    // add the census photons
    wrapped_cali_mark_begin("add the census photons");
    all_photons.insert(all_photons.end(), census_photons.begin(), census_photons.end());
    wrapped_cali_mark_end("add the census photons");

    wrapped_cali_mark_begin("set transported particles");
    imc_state.set_transported_particles(all_photons.size());
    wrapped_cali_mark_end("set transported particles");
    
    imc_state.print_memory_estimate(rank, n_ranks, mesh.get_n_local_cells(), all_photons.size());

    // add barrier here to make sure the transport timer starts at roughly the same time
    MPI_Barrier(MPI_COMM_WORLD);
    wrapped_cali_mark_begin("transport photons");
    census_photons = particle_pass_transport(mesh, gpu_setup, imc_parameters, mpi_info, mpi_types, imc_state, mctr, abs_E, track_E, all_photons, imc_parameters.get_n_omp_threads());
    wrapped_cali_mark_end("transport photons");

    wrapped_cali_mark_begin("update temperature");
    mesh.update_temperature(abs_E, track_E, imc_state);
    wrapped_cali_mark_end("update temperature");

    // update time for next step
    imc_state.print_conservation(imc_parameters.get_dd_mode());

    // write SILO file if it's enabled and it's the right cycle
    if (imc_parameters.get_write_silo_flag() &&
        !(imc_state.get_step() % imc_parameters.get_output_frequency())) {
      // write SILO file
      constexpr bool replicated_flag = false;
      double fake_mpi_runtime = 0.0;
      write_silo(mesh, imc_state.get_time(), imc_state.get_step(),
                 imc_state.get_rank_transport_runtime(), fake_mpi_runtime, rank,
                 n_ranks, replicated_flag);
    }

    meminfo_string = end_string + std::to_string(imc_state.get_step());

#ifdef USE_MEMORY_RECORD
    mem_record.read_meminfo(meminfo_string);
#endif

    wrapped_cali_mark_begin("update time for next step");
    imc_state.next_time_step();
    wrapped_cali_mark_end("update time for next step");
    
  }

#ifdef USE_MEMORY_RECORD
  mem_record.write_rss();
  mem_record.write_meminfo();
#endif
}

#endif // particle_pass_driver_h_

//---------------------------------------------------------------------------//
// end of particle_pass_driver.h
//---------------------------------------------------------------------------//
