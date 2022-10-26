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
inline Photon get_emission_photon(const Cell &cell, const double &phtn_E, const double &dt, const uint64_t seed) {
  using Constants::c;
  Photon emission_photon;
  RNG rng;
  rng.set_seed(seed);
  emission_photon.set_source_type(2);
  emission_photon.set_position(get_uniform_position_in_cell(cell, rng));
  emission_photon.set_angle(get_uniform_angle(rng));
  emission_photon.set_E0(phtn_E);
  emission_photon.set_distance_to_census(rng.generate_random_number() * c * dt);
  emission_photon.set_cell(cell.get_global_index());
  emission_photon.set_group(std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS)));
  emission_photon.set_rng(rng);
  return emission_photon;
}

//! Set input photon to the next emission photon
inline Photon get_boundary_source_photon(const Cell &cell, const double phtn_E, const double dt,
                                  const uint64_t seed, const int face) {
  using Constants::c;
  Photon source_photon;
  RNG rng;
  rng.set_seed(seed);
  source_photon.set_source_type(1);
  source_photon.set_position(get_uniform_position_on_face(cell, rng, face));
  source_photon.set_angle(get_source_angle_on_face(rng, face));
  source_photon.set_E0(phtn_E);
  source_photon.set_distance_to_census(rng.generate_random_number() * c * dt);
  source_photon.set_cell(cell.get_global_index());
  source_photon.set_group(std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS)));
  source_photon.set_rng(rng);
  return source_photon;
}

//! Set input photon to the next intiial census photon
inline Photon get_initial_census_photon(const Cell &cell, const double &phtn_E, const double &dt, const uint64_t seed) {
  using Constants::c;
  Photon census_photon;
  RNG rng;
  rng.set_seed(seed);
  census_photon.set_source_type(0);
  census_photon.set_position(get_uniform_position_in_cell(cell, rng));
  census_photon.set_angle(get_uniform_angle(rng));
  census_photon.set_E0(phtn_E);
  // initial census particles born at the beginning of timestep
  census_photon.set_distance_to_census(c * dt);
  census_photon.set_cell(cell.get_global_index());
  census_photon.set_group(std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS)));
  census_photon.set_rng(rng);
  return census_photon;
}


//! Make the census photons on cycle 0
std::vector<Photon> make_initial_census_photons(const double dt, const Mesh &mesh, const int rank,  const uint64_t n_user_photons, const double total_E) {
  std::vector<Photon> initial_census_photons;
  auto E_cell_census = mesh.get_census_E();
  const uint64_t rank_seed_offset{n_user_photons * rank};
  uint64_t ith_census = 0;
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    if (E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0) t_num_census = 1;
      // keep track of census energy for conservation check
      const double photon_census_E = E_cell_census[i] / t_num_census;
      for (uint32_t p=0; p<t_num_census;++p)
        initial_census_photons.push_back(get_initial_census_photon(cell, photon_census_E, dt, rank_seed_offset + ith_census));
    }
    ith_census++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return initial_census_photons;
}

std::vector<Photon> make_photons(const double dt, const Mesh &mesh, const int rank, const uint32_t cycle, const uint64_t n_user_photons, const double total_E) {

  auto E_cell_emission = mesh.get_emission_E();
  auto E_cell_source = mesh.get_source_E();
  // for RNG offsets, each cycle allows for one hundred million particles across one hundred
  // thousand ranks, increment the ten trillon place for the next cycle, using cycle plus one for
  // the cycle offset gives the initial census their own space
  const uint64_t cycle_seed_offset{10000000000000UL * static_cast<uint64_t>(cycle+1)};
  const uint64_t rank_seed_offset{n_user_photons * static_cast<uint64_t>(rank)};

  // figure out how many to make to size all_photons vector
  uint64_t photons_to_make = 0;
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission =
          int(n_user_photons * E_cell_emission[i] / total_E);
      // make at least one photon to represent emission energy
      if (t_num_emission == 0)
        t_num_emission = 1;
      photons_to_make+=t_num_emission;
    }
    if (E_cell_source[i] > 0.0) {
      // boundary source
        uint32_t t_num_source =
            int(n_user_photons * E_cell_source[i] / total_E);
        // make at least one photon to represent source energy
        if (t_num_source == 0)
          t_num_source = 1;
      photons_to_make+=t_num_source;
    }
  }

  std::vector<Photon> all_photons;
  all_photons.reserve(photons_to_make);

  // use this to increment the seed for each particle
  uint64_t ith_photon{0UL};

  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission =
          int(n_user_photons * E_cell_emission[i] / total_E);
      // make at least one photon to represent emission energy
      if (t_num_emission == 0)
        t_num_emission = 1;
      const double photon_emission_E = E_cell_emission[i] / t_num_emission;
      for (uint32_t p=0; p<t_num_emission;++p)
        all_photons.push_back(get_emission_photon(cell, photon_emission_E, dt, (cycle_seed_offset + rank_seed_offset+ith_photon)));

    }
    if (E_cell_source[i] > 0.0) {
      // boundary source
        uint32_t t_num_source =
            int(n_user_photons * E_cell_source[i] / total_E);
        // make at least one photon to represent source energy
        if (t_num_source == 0)
          t_num_source = 1;
        const double photon_source_E = E_cell_source[i] / t_num_source;
        const int face = cell.get_source_face();
        for (uint32_t p=0; p<t_num_source;++p)
          all_photons.push_back(get_boundary_source_photon(cell, photon_source_E, dt, (cycle_seed_offset + rank_seed_offset+ith_photon), face));
    }
    ith_photon++;
  }
  return all_photons;
}


#endif // source_h_
//----------------------------------------------------------------------------//
// end of source.h
//----------------------------------------------------------------------------//
