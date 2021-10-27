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

#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "census_creation.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "mesh.h"
#include "replicated_transport.h"
#include "timer.h"

void imc_replicated_driver(Mesh &mesh, IMC_State &imc_state, const IMC_Parameters &imc_parameters) {
  using std::vector;
  vector<double> abs_E(mesh.get_global_num_cells(), 0.0);
  vector<double> track_E(mesh.get_global_num_cells(), 0.0);
  vector<Photon> census_photons;
  auto n_user_photons = imc_parameters.get_n_user_photon();

  // set the max census size to be 10% more than the user requested photons.
  // Note that this does not need to be divided by the number of ranks
  // because the comb get the global census energy with an MPI_Allreduce
  uint64_t max_census_size = static_cast<uint64_t>(1.1 * n_user_photons);

  while (!imc_state.finished()) {
    imc_state.print_timestep_header();

    //--------------------------------------------------------------------------------------------//
    // work done on the host processor, don't count cycles
    //--------------------------------------------------------------------------------------------//

    // set opacity, Fleck factor, all energy to source
    mesh.calculate_photon_energy(imc_state);

    // all reduce to get total source energy to make correct number of particles on each rank
    double global_source_energy = mesh.get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // setup source
    if (imc_state.get_step() == 1)
      census_photons = make_initial_census_photons(imc_state.get_dt(), mesh, n_user_photons, global_source_energy, imc_state.get_rng());
    imc_state.set_pre_census_E(get_photon_list_E(census_photons));
    // make emission and source photons
    auto all_photons = make_photons(imc_state.get_dt(), mesh, n_user_photons, global_source_energy, imc_state.get_rng());
    // add the census photons
    all_photons.insert(all_photons.end(), census_photons.begin(), census_photons.end());

    imc_state.set_transported_particles(all_photons.size());

    // source all the particles, move them into the "particles" vector in each cel
    for (auto &p : all_photons) {
      auto &cell = mesh.get_cell_ref(p.get_cell());
      cell.push_particle(p);
    }

    //--------------------------------------------------------------------------------------------//
    // work sent to APEs count cycles to estimate performance
    //--------------------------------------------------------------------------------------------//

    // transport photons, return the particles that reached census
    census_photons = replicated_transport(mesh, imc_state, abs_E,
                                          track_E, max_census_size);

    mesh.update_temperature(abs_E, track_E, imc_state);

    imc_state.print_conservation();

    // update time for next step
    imc_state.next_time_step();
  }
}

#endif // replicated_driver_h_

//---------------------------------------------------------------------------//
// end of replicated_driver.h
//---------------------------------------------------------------------------//
