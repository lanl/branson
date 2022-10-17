//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   source.h
 * \author Alex Long
 * \date   December 2 2015
 * \brief  Allows transport function to create particles when needed
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef source_h_
#define source_h_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "constants.h"
#include "cell.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"


//! Set input photon to the next emission photon
inline Photon get_emission_photon(const Cell &cell, const double &phtn_E, const double &dt, RNG *rng) {
  using Constants::c;
  Photon emission_photon;
  emission_photon.set_source_type(2);
  emission_photon.set_position(get_uniform_position_in_cell(cell, rng));
  emission_photon.set_angle(get_uniform_angle(rng));
  emission_photon.set_E0(phtn_E);
  emission_photon.set_distance_to_census(rng->generate_random_number() * c * dt);
  emission_photon.set_cell(cell.get_ID());
  emission_photon.set_group(std::floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
  return emission_photon;
}

//! Set input photon to the next emission photon
inline Photon get_boundary_source_photon(const Cell &cell, const double &phtn_E, const double &dt, RNG *rng,
                                  int face) {
  Photon source_photon;
  using Constants::c;
  source_photon.set_source_type(1);
  source_photon.set_position(get_uniform_position_on_face(cell, rng, face));
  source_photon.set_angle(get_source_angle_on_face(rng, face));
  source_photon.set_E0(phtn_E);
  source_photon.set_distance_to_census(rng->generate_random_number() * c * dt);
  source_photon.set_cell(cell.get_ID());
  source_photon.set_group(std::floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
  return source_photon;
}

//! Set input photon to the next intiial census photon
inline Photon get_initial_census_photon(const Cell &cell, const double &phtn_E, const double &dt, RNG *rng) {
  using Constants::c;
  Photon census_photon;
  census_photon.set_source_type(0);
  census_photon.set_position(get_uniform_position_in_cell(cell, rng));
  census_photon.set_angle(get_uniform_angle(rng));
  census_photon.set_E0(phtn_E);
  // initial census particles born at the beginning of timestep
  census_photon.set_distance_to_census(c * dt);
  census_photon.set_cell(cell.get_ID());
  census_photon.set_group(std::floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
  return census_photon;
}


//! Make the census photons on cycle 0
std::vector<Photon> make_initial_census_photons(const double dt, const Mesh &mesh, const uint64_t user_photons, const double total_E, RNG *rng) {
  std::vector<Photon> initial_census_photons;
  auto E_cell_census = mesh.get_census_E();
  for (auto const &cell : mesh) {
    int i = cell.get_ID();
    if (E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0) t_num_census = 1;
      // keep track of census energy for conservation check
      const double photon_census_E = E_cell_census[i] / t_num_census;
      for (int p=0; p<t_num_census;++p)
        initial_census_photons.push_back(get_initial_census_photon(cell, photon_census_E, dt, rng));
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return initial_census_photons;
}

std::vector<Photon> make_photons(const double dt, const Mesh &mesh, const uint64_t user_photons, const double total_E, RNG *rng) {

  auto E_cell_emission = mesh.get_emission_E();
  auto E_cell_source = mesh.get_source_E();

  double total_census_E = 0.0;

  // make work packets
  // current cell pointer
  std::vector<Photon> all_photons;
  for (auto const &cell : mesh) {
    int i = cell.get_ID();
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission =
          int(user_photons * E_cell_emission[i] / total_E);
      // make at least one photon to represent emission energy
      if (t_num_emission == 0)
        t_num_emission = 1;
      const double photon_emission_E = E_cell_emission[i] / t_num_emission;
      for (int p=0; p<t_num_emission;++p)
        all_photons.push_back(get_emission_photon(cell, photon_emission_E, dt, rng));

    }
    if (E_cell_source[i] > 0.0) {
      // boundary source
        uint32_t t_num_source =
            int(user_photons * E_cell_source[i] / total_E);
        // make at least one photon to represent emission energy
        if (t_num_source == 0)
          t_num_source = 1;
        const double photon_source_E = E_cell_source[i] / t_num_source;
        const int face = cell.get_source_face();
        for (int p=0; p<t_num_source;++p)
          all_photons.push_back(get_boundary_source_photon(cell, photon_source_E, dt, rng, face));
    }
  }
  return all_photons;
}


#endif // source_h_
//----------------------------------------------------------------------------//
// end of source.h
//----------------------------------------------------------------------------//
