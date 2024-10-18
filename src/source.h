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
inline void set_emission_photon(PhotonArray &photon_array, size_t index, const Cell &cell, const double &phtn_E, const double &dt, const uint32_t seed, const uint64_t stream_num) {
  using Constants::c;
  // Photon emission_photon;
  // Photon emission_photon(photon_array, index);
  RNG rng(seed, stream_num);
  photon_array.source_type[index] = 2;
  photon_array.pos[index] = get_uniform_position_in_cell(cell, rng);
  photon_array.angle[index] = get_uniform_angle(rng);
  photon_array.E[index] = phtn_E;
  photon_array.E0[index] = phtn_E;
  photon_array.life_dx[index] = rng.generate_random_number() * c * dt;
  photon_array.cell_ID[index] = cell.get_global_index();
  photon_array.group[index] = std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS));
  photon_array.rng[index] = rng;
  photon_array.descriptors[index] = static_cast<unsigned char>(Constants::PASS);
  // return emission_photon;
}

//! Set input photon to the next emission photon
inline void set_boundary_source_photon(PhotonArray &photon_array, size_t index, const Cell &cell, const double phtn_E, const double dt, const uint32_t seed, const uint64_t stream_num, const int face) {
  using Constants::c;
  // Photon source_photon;
  // Photon source_photon(photon_array, index);
  RNG rng(seed, stream_num);
  photon_array.source_type[index] = 1;
  photon_array.pos[index] = get_uniform_position_on_face(cell, rng, face);
  photon_array.angle[index] = get_source_angle_on_face(rng, face);
  photon_array.E[index] = phtn_E;
  photon_array.E0[index] = phtn_E;
  photon_array.life_dx[index] = rng.generate_random_number() * c * dt;
  photon_array.cell_ID[index] = cell.get_global_index();
  photon_array.group[index] = std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS));
  photon_array.rng[index] = rng;
  photon_array.descriptors[index] = static_cast<unsigned char>(Constants::PASS);
  // return source_photon;
}

//! Set input photon to the next intiial census photon
inline void set_initial_census_photon(PhotonArray &photon_array, size_t index, const Cell &cell, const double &phtn_E, const double &dt, const uint32_t seed, const uint64_t stream_num) {
  using Constants::c;
  // Photon census_photon;
  // Photon census_photon(photon_array, index);
  RNG rng(seed, stream_num);
  photon_array.source_type[index] = 0;
  photon_array.pos[index] = get_uniform_position_in_cell(cell, rng);
  photon_array.angle[index] = get_uniform_angle(rng);
  photon_array.E[index] = phtn_E;
  photon_array.E0[index] = phtn_E;
  photon_array.life_dx[index] = c * dt;
  photon_array.cell_ID[index] = cell.get_global_index();
  photon_array.group[index] = std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS));
  photon_array.rng[index] = rng;
  photon_array.descriptors[index] = static_cast<unsigned char>(Constants::PASS);
  // return census_photon;
}

//! Set input photon to the next emission photon
inline Photon get_emission_photon(const Cell &cell, const double &phtn_E, const double &dt, const uint32_t seed, const uint64_t stream_num) {
  using Constants::c;
  Photon emission_photon;
  RNG rng(seed, stream_num);
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
                                  const uint32_t seed, const uint64_t stream_num, const int face) {
  using Constants::c;
  Photon source_photon;
  RNG rng(seed, stream_num);
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
inline Photon get_initial_census_photon(const Cell &cell, const double &phtn_E, const double &dt, const uint32_t seed, const uint64_t stream_num) {
  using Constants::c;
  Photon census_photon;
  RNG rng(seed, stream_num);
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
template <typename Census_T>
Census_T make_initial_census_photons(const double dt, const Mesh &mesh, const int rank, const uint32_t seed, const uint64_t n_user_photons, const double total_E);


template <>
PhotonArray make_initial_census_photons<PhotonArray>(const double dt, const Mesh &mesh, const int rank, const uint32_t seed, const uint64_t n_user_photons, const double total_E) {

  PhotonArray initial_census_photons;

  auto E_cell_census = mesh.get_census_E();
  const uint64_t rank_stream_num_offset{n_user_photons * rank};
  uint64_t ith_census = 0;
  size_t total_photons = 0;
  
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    if (E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      if (t_num_census == 0) 
        t_num_census = 1;
      total_photons += t_num_census;
    }
  }

  initial_census_photons.resize(total_photons);
  size_t current_index = 0;

  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    if (E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0) 
        t_num_census = 1;
      // keep track of census energy for conservation check
      const double photon_census_E = E_cell_census[i] / t_num_census;
      for (uint32_t p=0; p<t_num_census;++p) {
        set_initial_census_photon(initial_census_photons, current_index, cell, photon_census_E, dt, seed, rank_stream_num_offset + ith_census);
        // initial_census_photons.push_back(get_initial_census_photon(cell, photon_census_E, dt, seed, rank_stream_num_offset + ith_census));
        ith_census++;
        current_index++;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return initial_census_photons;
}

template <>
std::vector<Photon> make_initial_census_photons<std::vector<Photon>>(const double dt, const Mesh &mesh, const int rank, const uint32_t seed, const uint64_t n_user_photons, const double total_E) {
  std::vector<Photon> initial_census_photons;
  auto E_cell_census = mesh.get_census_E();
  const uint64_t rank_stream_num_offset{n_user_photons * rank};
  uint64_t ith_census = 0;
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    if (E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0) t_num_census = 1;
      // keep track of census energy for conservation check
      const double photon_census_E = E_cell_census[i] / t_num_census;
      for (uint32_t p=0; p<t_num_census;++p) {
        initial_census_photons.push_back(get_initial_census_photon(cell, photon_census_E, dt, seed, rank_stream_num_offset + ith_census));
        ith_census++;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return initial_census_photons;
}

// template function for making photons
template <typename Census_T>
Census_T make_photons(const double dt, const Mesh &mesh, const int rank, const uint32_t cycle,
                    const uint32_t seed, const uint64_t n_user_photons, const double total_E);

// explicit template specializations for PhotonArray and std::vector<Photon> 
template <>
PhotonArray make_photons<PhotonArray>(const double dt, const Mesh &mesh, const int rank, const uint32_t cycle,
                      const uint32_t seed, const uint64_t n_user_photons, const double total_E) {

  auto E_cell_emission = mesh.get_emission_E();
  auto E_cell_source = mesh.get_source_E();
  // for RNG offsets, each cycle allows for one hundred million particles across one hundred
  // thousand ranks, increment the ten trillon place for the next cycle, using cycle plus one for
  // the cycle offset gives the initial census their own space
  const uint64_t cycle_stream_num_offset{10000000000000UL * static_cast<uint64_t>(cycle)};
  const uint64_t rank_stream_num_offset{n_user_photons * static_cast<uint64_t>(rank)};

  // figure out how many to make to size all_photons vector
  uint64_t photons_to_make = 0;
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission = int(n_user_photons * E_cell_emission[i] / total_E);
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

  // std::vector<Photon> all_photons;
  PhotonArray all_photons;
  all_photons.resize(photons_to_make);
  size_t current_index = 0;
  // all_photons.reserve(photons_to_make);

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
      for (uint32_t p=0; p<t_num_emission;++p) {
        set_emission_photon(all_photons, current_index, cell, photon_emission_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon));
        // all_photons.push_back(get_emission_photon(cell, photon_emission_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon)));
        ith_photon++;
        current_index++;
      }
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
      for (uint32_t p=0; p<t_num_source;++p) {
        set_boundary_source_photon(all_photons, current_index, cell, photon_source_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon), face);
        // all_photons.push_back(get_boundary_source_photon(cell, photon_source_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon), face));
        current_index++;
        ith_photon++;
      }
    }
  }
  return all_photons;
}

template <>
std::vector<Photon> make_photons<std::vector<Photon>>(const double dt, const Mesh &mesh, const int rank, 
    const uint32_t cycle, const uint32_t seed, const uint64_t n_user_photons, const double total_E) {

  auto E_cell_emission = mesh.get_emission_E();
  auto E_cell_source = mesh.get_source_E();
  // for RNG offsets, each cycle allows for one hundred million particles across one hundred
  // thousand ranks, increment the ten trillon place for the next cycle, using cycle plus one for
  // the cycle offset gives the initial census their own space
  const uint64_t cycle_stream_num_offset{10000000000000UL * static_cast<uint64_t>(cycle)};
  const uint64_t rank_stream_num_offset{n_user_photons * static_cast<uint64_t>(rank)};

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
      for (uint32_t p=0; p<t_num_emission;++p) {
        all_photons.push_back(get_emission_photon(cell, photon_emission_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon)));
        ith_photon++;
      }
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
      for (uint32_t p=0; p<t_num_source;++p) {
        all_photons.push_back(get_boundary_source_photon(cell, photon_source_E, dt, seed, (cycle_stream_num_offset + rank_stream_num_offset+ith_photon), face));
        ith_photon++;
      }
    }
  }
  return all_photons;
}


#endif // source_h_
//----------------------------------------------------------------------------//
// end of source.h
//----------------------------------------------------------------------------//
