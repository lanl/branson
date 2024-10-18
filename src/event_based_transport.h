//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   event_based_transport.h
 * \author Joseph Farmer
 * \date   Auhust 1 2024
 * \brief  IMC transport with an event-based algorithm 
 * \note   Copyright (C) 2024 Triad National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
#ifndef event_based_transport_h_
#define event_based_transport_h_

#include <vector>

#include "config.h"
#include "RNG.h"
#include "cell_tally.h"
#include "constants.h"
#include "photon.h"
#include "sampling_functions.h"

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
enum EventType {SCATTER, BOUNDARY, CENSUS, KILLED};
struct Event {
  size_t photon_index;
  EventType type;
  double distance;
};

inline void precompute_data(const uint32_t rank_cell_offset, const PhotonArray &photon_array, const Cell *cells, size_t *active_photons, size_t active_count, std::vector<double> &sigma_s, std::vector<double> &sigma_a, std::vector<double> &f, std::vector<double> &total_sigma_s, std::vector<uint32_t> &local_cell_indices) {
//
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    size_t photon_index = active_photons[i];
    local_cell_indices[i] = photon_array.cell_ID[photon_index] - rank_cell_offset;
    const Cell *cell = &cells[local_cell_indices[i]];
    // cell_ptrs[i] = &cells[local_cell_indices[i]];
    sigma_s[i] = cell->get_op_s(photon_array.group[photon_index]);
    sigma_a[i] = cell->get_op_a(photon_array.group[photon_index]);
    f[i] = cell->get_f();
    total_sigma_s[i] = (1.0 - f[i]) * sigma_a[i] + sigma_s[i];
    // random_numbers[i] = photon_array.rng[photon_index].generate_random_number();
  }
}

inline void precompute_data(const uint32_t rank_cell_offset, const std::vector<Photon> &photon_array, const Cell *cells, size_t *active_photons, size_t active_count, std::vector<double> &sigma_s, std::vector<double> &sigma_a, std::vector<double> &f, std::vector<double> &total_sigma_s, std::vector<uint32_t> &local_cell_indices) {
//
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    size_t photon_index = active_photons[i];
    const auto &phtn = photon_array[photon_index];
    local_cell_indices[i] = phtn.get_cell() - rank_cell_offset;
    const Cell *cell = &cells[local_cell_indices[i]];
    // cell_ptrs[i] = &cells[local_cell_indices[i]];
    sigma_s[i] = cell->get_op_s(phtn.get_group());
    sigma_a[i] = cell->get_op_a(phtn.get_group());
    f[i] = cell->get_f();
    total_sigma_s[i] = (1.0 - f[i]) * sigma_a[i] + sigma_s[i];
    // random_numbers[i] = photon_array.rng[photon_index].generate_random_number();
  }
}

inline void calculate_distances(PhotonArray &photon_array, const Cell *cells, const size_t *active_photons, size_t active_count, const std::vector<double> &total_sigma_s, const std::vector<uint32_t> &local_cell_indices, std::vector<Event> &events) {
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    size_t photon_index = active_photons[i];
    Cell const *cell = &cells[local_cell_indices[i]];
    uint32_t surface_cross = 0;
    const double dist_to_scatter = (total_sigma_s[i] > 0.0) ? -log(photon_array.rng[photon_index].generate_random_number()) / total_sigma_s[i] : 1e100;
    const double dist_to_boundary = cell->get_distance_to_boundary(photon_array.pos[photon_index], photon_array.angle[photon_index], surface_cross);
    const double dist_to_census = photon_array.life_dx[photon_index];
    events[i].distance = std::min(dist_to_scatter, std::min(dist_to_boundary, dist_to_census));
    events[i].photon_index = photon_index;

    if (events[i].distance == dist_to_scatter) 
      events[i].type = SCATTER;
    else if (events[i].distance == dist_to_boundary)
      events[i].type = BOUNDARY;
    else 
      events[i].type = CENSUS;
  }
}

inline void calculate_distances(std::vector<Photon> &photon_array, const Cell *cells, const size_t *active_photons, size_t active_count, const std::vector<double> &total_sigma_s, const std::vector<uint32_t> &local_cell_indices, std::vector<Event> &events) {
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    size_t photon_index = active_photons[i];
    auto &phtn = photon_array[photon_index];
    Cell const *cell = &cells[local_cell_indices[i]];
    auto &rng = photon_array[i].get_rng();
    uint32_t surface_cross = 0;
    const double dist_to_scatter = (total_sigma_s[i] > 0.0) ? -log(rng.generate_random_number()) / total_sigma_s[i] : 1e100;
    const double dist_to_boundary = cell->get_distance_to_boundary(phtn.get_position(), phtn.get_angle(), surface_cross);
    const double dist_to_census = phtn.get_distance_remaining();
    events[i].distance = std::min(dist_to_scatter, std::min(dist_to_boundary, dist_to_census));
    events[i].photon_index = photon_index;

    if (events[i].distance == dist_to_scatter) 
      events[i].type = SCATTER;
    else if (events[i].distance == dist_to_boundary)
      events[i].type = BOUNDARY;
    else 
      events[i].type = CENSUS;
  }
}

inline void process_scatter_events(PhotonArray &photon_array, const std::vector<double> &sigma_s, const std::vector<double> &total_sigma_s, const std::vector<uint32_t> &local_cell_indices, const Cell *cells, const size_t *scatter_events, size_t scatter_count, size_t &counter, const std::vector<EmissionGroupData> &emission_groups) {
#pragma omp simd
  for (size_t i = 0; i < scatter_count; ++i) {
    RNG &rng = photon_array.rng[scatter_events[i]];
    photon_array.angle[scatter_events[i]] = get_uniform_angle(rng);
    photon_array.descriptors[scatter_events[i]] = static_cast<unsigned char>(Constants::SCATTER);
//   }
// #pragma omp simd
//   for (size_t i = 0; i < scatter_count; ++i) {
//     RNG &rng = photon_array.rng[scatter_events[i]];
    if (rng.generate_random_number() > (sigma_s[i] / total_sigma_s[i])) {
      // photon_array.group[scatter_events[i]] = sample_emission_group(rng, cells[local_cell_indices[scatter_events[i]]]);
      photon_array.group[scatter_events[i]] = sample_emission_group(rng, emission_groups[local_cell_indices[scatter_events[i]]]);
      // counter++;
    }
  }
}

inline void process_scatter_events(std::vector<Photon> &photon_array, const std::vector<double> &sigma_s, const std::vector<double> &total_sigma_s, const std::vector<uint32_t> &local_cell_indices, const Cell *cells, const size_t *scatter_events, size_t scatter_count, size_t &counter, const std::vector<EmissionGroupData> &emission_groups) {
#pragma omp simd
  for (size_t i = 0; i < scatter_count; ++i) {
    auto &phtn = photon_array[scatter_events[i]];
    RNG &rng = phtn.get_rng();
    phtn.set_angle(get_uniform_angle(rng));
    phtn.set_descriptor(Constants::SCATTER);
//   }
// #pragma omp simd
//   for (size_t i = 0; i < scatter_count; ++i) {
//     RNG &rng = photon_array.rng[scatter_events[i]];
    if (rng.generate_random_number() > (sigma_s[i] / total_sigma_s[i])) {
      // photon_array.group[scatter_events[i]] = sample_emission_group(rng, cells[local_cell_indices[scatter_events[i]]]);
      phtn.set_group(sample_emission_group(rng, emission_groups[local_cell_indices[scatter_events[i]]]));
      // counter++;
    }
  }
}

inline void process_boundary_events(PhotonArray &photon_array, const size_t *boundary_events, const Cell *cells, uint32_t rank_cell_offset, size_t boundary_count) {
#pragma omp simd 
  for (size_t i = 0; i < boundary_count; ++i) {
    // size_t photon_index = boundary_events[i];
    uint32_t local_cell_index = photon_array.cell_ID[boundary_events[i]] - rank_cell_offset;
    Cell const *cell = &cells[local_cell_index];
    uint32_t surface_cross = 0;
    cell->get_distance_to_boundary(photon_array.pos[boundary_events[i]], photon_array.angle[boundary_events[i]], surface_cross);
    auto boundary_event = cell->get_bc(surface_cross);
    if (boundary_event == Constants::ELEMENT) {
      photon_array.cell_ID[boundary_events[i]] = cell->get_next_cell(surface_cross);
      photon_array.descriptors[boundary_events[i]] = static_cast<unsigned char>(Constants::BOUND);
    } else if (boundary_event == Constants::PROCESSOR) {
      photon_array.cell_ID[boundary_events[i]] = cell->get_next_cell(surface_cross);
      photon_array.descriptors[boundary_events[i]] = static_cast<unsigned char>(Constants::PASS);
    } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
      photon_array.descriptors[boundary_events[i]] = static_cast<unsigned char>(Constants::EXIT);
    } else {
      int reflect_angle = surface_cross / 2;
      photon_array.angle[boundary_events[i]][reflect_angle] = -photon_array.angle[boundary_events[i]][reflect_angle];
      photon_array.descriptors[boundary_events[i]] = static_cast<unsigned char>(Constants::BOUND);
    }
  }
}

inline void process_boundary_events(std::vector<Photon> &photon_array, const size_t *boundary_events, const Cell *cells, uint32_t rank_cell_offset, size_t boundary_count) {
#pragma omp simd 
  for (size_t i = 0; i < boundary_count; ++i) {
    auto & phtn = photon_array[boundary_events[i]];
    uint32_t local_cell_index = phtn.get_cell() - rank_cell_offset;
    Cell const *cell = &cells[local_cell_index];
    uint32_t surface_cross = 0;
    cell->get_distance_to_boundary(phtn.get_position(), phtn.get_angle(), surface_cross);
    auto boundary_event = cell->get_bc(surface_cross);
    if (boundary_event == Constants::ELEMENT) {
      phtn.set_cell(cell->get_next_cell(surface_cross));
      phtn.set_descriptor(Constants::BOUND);
    } else if (boundary_event == Constants::PROCESSOR) {
      phtn.set_cell(cell->get_next_cell(surface_cross));
      phtn.set_descriptor(Constants::PASS);
    } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
      phtn.set_descriptor(Constants::EXIT);
    } else {
      phtn.reflect(surface_cross);
      phtn.set_descriptor(Constants::BOUND);
    }
  }
}

void process_census_events(PhotonArray &photon_array, const size_t *census_events, size_t census_count) {
// #pragma omp simd
  for (size_t i = 0; i < census_count; ++i) {
    photon_array.descriptors[census_events[i]] = static_cast<unsigned char>(Constants::CENSUS);
  }
}

void process_census_events(std::vector<Photon> &photon_array, const size_t *census_events, size_t census_count) {
// #pragma omp simd
  for (size_t i = 0; i < census_count; ++i) {
    photon_array[census_events[i]].set_descriptor(Constants::CENSUS);
  }
}

void process_killed_events(PhotonArray &photon_array, const size_t *killed_events, const Cell *cells, Cell_Tally *cell_tallies, const uint32_t rank_cell_offset, size_t killed_count) {
// #pragma omp simd
  for (size_t i = 0; i < killed_count; ++i) {
    uint32_t local_cell_index = photon_array.cell_ID[killed_events[i]] - rank_cell_offset;
    // photon_array.descriptors[photon_index][0] = static_cast<unsigned char>(Constants::KILLED);
    cell_tallies[local_cell_index].accumulate_absorbed_E(photon_array.E[killed_events[i]]);
  }
}

void process_killed_events(std::vector<Photon> &photon_array, const size_t *killed_events, const Cell *cells, Cell_Tally *cell_tallies, const uint32_t rank_cell_offset, size_t killed_count) {
// #pragma omp simd
  for (size_t i = 0; i < killed_count; ++i) {
    auto &phtn = photon_array[killed_events[i]];
    uint32_t local_cell_index = phtn.get_cell() - rank_cell_offset;
    cell_tallies[local_cell_index].accumulate_absorbed_E(phtn.get_E());
  }
}

inline void update_photon_state(PhotonArray &photon_array, Cell_Tally *cell_tallies, const std::vector<Event> &events, const std::vector<double> &sigma_a, const std::vector<double> &f, const std::vector<uint32_t> &local_cell_indices, size_t active_count, std::vector<size_t> &scatter_events, std::vector<size_t> &boundary_events, std::vector<size_t> &census_events, std::vector<size_t> &killed_events, size_t &scatter_count, size_t &boundary_count, size_t &census_count, size_t &killed_count, size_t cellSize) {

  // double tmpAbse[cellSize];
  // double tmpTracke[cellSize];
  // for(size_t i = 0 ; i < cellSize; ++i){
  //   tmpAbse[i] = 0.;
  //   tmpTracke[i] = 0.;
  // }
  std::vector<double> absorbed_Es(active_count);
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    absorbed_Es[i] = photon_array.E[events[i].photon_index] * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
  }
// #pragma omp simd reduction(+ : tmpAbse, tmpTracke)
// #pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    // const double absorbed_E = photon_array.E[events[i].photon_index] * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
    cell_tallies[local_cell_indices[i]].accumulate_absorbed_E(absorbed_Es[i]);
    cell_tallies[local_cell_indices[i]].accumulate_track_E(absorbed_Es[i]/ (sigma_a[i] * f[i]));
  }
    // tmpAbsE[local_cell_indices[i]] += absorbed_E;
    // tmpTrackE[local_cell_indices[i]] += absorbed_E / (sigma_a[i] * f[i]);
#pragma omp simd
for (size_t i = 0; i < active_count; ++i) {
    // const double absorbed_E = photon_array.E[events[i].photon_index] * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
    photon_array.E[events[i].photon_index] -= absorbed_Es[i];
    photon_array.pos[events[i].photon_index][0] += photon_array.angle[events[i].photon_index][0] * events[i].distance;
    photon_array.pos[events[i].photon_index][1] += photon_array.angle[events[i].photon_index][1] * events[i].distance;
    photon_array.pos[events[i].photon_index][2] += photon_array.angle[events[i].photon_index][2] * events[i].distance;
    photon_array.life_dx[events[i].photon_index] -= events[i].distance;
  }

#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    const Event &event = events[i]; 
    size_t photon_index  = event.photon_index;
    if (photon_array.E[photon_index] / photon_array.E0[photon_index] < Constants::cutoff_fraction) {
      photon_array.descriptors[photon_index] = static_cast<unsigned char>(Constants::KILLED);
      killed_events[killed_count++] = photon_index;
    } else {
      switch(event.type) {
      case SCATTER:
        scatter_events[scatter_count++]  = photon_index;
        break;
      case BOUNDARY: 
        boundary_events[boundary_count++]  = photon_index;
        break;
      case CENSUS:
        census_events[census_count++]   = photon_index;
        break;
      default:
        break;
      }
    }
  }
}

inline void update_photon_state(std::vector<Photon> &photon_array, Cell_Tally *cell_tallies, const std::vector<Event> &events, const std::vector<double> &sigma_a, const std::vector<double> &f, const std::vector<uint32_t> &local_cell_indices, size_t active_count, std::vector<size_t> &scatter_events, std::vector<size_t> &boundary_events, std::vector<size_t> &census_events, std::vector<size_t> &killed_events, size_t &scatter_count, size_t &boundary_count, size_t &census_count, size_t &killed_count, size_t cellSize) {

  // double tmpAbse[cellSize];
  // double tmpTracke[cellSize];
  // for(size_t i = 0 ; i < cellSize; ++i){
  //   tmpAbse[i] = 0.;
  //   tmpTracke[i] = 0.;
  // }
  std::vector<double> absorbed_Es(active_count);
#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    absorbed_Es[i] = photon_array[events[i].photon_index].get_E() * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
  }
// #pragma omp simd reduction(+ : tmpAbse, tmpTracke)
// #pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    // const double absorbed_E = photon_array.E[events[i].photon_index] * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
    cell_tallies[local_cell_indices[i]].accumulate_absorbed_E(absorbed_Es[i]);
    cell_tallies[local_cell_indices[i]].accumulate_track_E(absorbed_Es[i]/ (sigma_a[i] * f[i]));
  }
    // tmpAbsE[local_cell_indices[i]] += absorbed_E;
    // tmpTrackE[local_cell_indices[i]] += absorbed_E / (sigma_a[i] * f[i]);
#pragma omp simd
for (size_t i = 0; i < active_count; ++i) {
    // const double absorbed_E = photon_array.E[events[i].photon_index] * (1.0 - exp(-sigma_a[i] * f[i] * events[i].distance));
    photon_array[events[i].photon_index].set_E(photon_array[events[i].photon_index].get_E() - absorbed_Es[i]);
    photon_array[events[i].photon_index].move(events[i].distance);
  }

#pragma omp simd
  for (size_t i = 0; i < active_count; ++i) {
    const Event &event = events[i]; 
    size_t photon_index  = event.photon_index;
    if (photon_array[photon_index].below_cutoff(Constants::cutoff_fraction)) {
      photon_array[photon_index].set_descriptor(Constants::KILLED);
      killed_events[killed_count++] = photon_index;
    } else {
      switch(event.type) {
      case SCATTER:
        scatter_events[scatter_count++]  = photon_index;
        break;
      case BOUNDARY: 
        boundary_events[boundary_count++]  = photon_index;
        break;
      case CENSUS:
        census_events[census_count++]   = photon_index;
        break;
      default:
        break;
      }
    }
  }
}

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void event_transport_photon(const uint32_t rank_cell_offset, PhotonArray &photon_array, const Cell *cells, Cell_Tally *cell_tallies, size_t cellSize, const std::vector<EmissionGroupData> &emission_groups) {

  using Constants::bc_type;
  using Constants::c;
  // events
  using std::min;
  const size_t maxPhotons = photon_array.cell_ID.size();

  std::vector<size_t> scatter_events(maxPhotons), boundary_events(maxPhotons), census_events(maxPhotons), killed_events(maxPhotons), active_photons(maxPhotons);
  std::vector<Event> events(maxPhotons);
  std::vector<uint32_t> local_cell_indices(maxPhotons);
  std::vector<double> sigma_s(maxPhotons), sigma_a(maxPhotons), f(maxPhotons), total_sigma_s(maxPhotons);
  
  std::iota(active_photons.begin(), active_photons.end(), 0); 

  size_t active_count = maxPhotons;
  size_t counter = 0;

  while (active_count > 0) {
    size_t scatter_count = 0, boundary_count = 0, census_count = 0, killed_count = 0;

    precompute_data(rank_cell_offset, photon_array, cells, active_photons.data(), active_count, sigma_s, sigma_a, f, total_sigma_s, local_cell_indices);
    calculate_distances(photon_array, cells, active_photons.data(), active_count, total_sigma_s, local_cell_indices, events);
    update_photon_state(photon_array, cell_tallies, events, sigma_a, f , local_cell_indices, active_count, scatter_events, boundary_events, census_events, killed_events, scatter_count, boundary_count, census_count, killed_count, cellSize);

    process_scatter_events(photon_array, sigma_s, total_sigma_s, local_cell_indices, cells, scatter_events.data(), scatter_count, counter, emission_groups);
    process_census_events(photon_array, census_events.data(), census_count);
    process_boundary_events(photon_array, boundary_events.data(), cells, rank_cell_offset, boundary_count);
    process_killed_events(photon_array, killed_events.data(), cells,  cell_tallies, rank_cell_offset, killed_count);

    active_photons.erase(std::remove_if(active_photons.begin(), active_photons.end(), [&photon_array](size_t index){
      return photon_array.descriptors[index] == static_cast<unsigned char>(Constants::KILLED) ||
      photon_array.descriptors[index] == static_cast<unsigned char>(Constants::CENSUS) ||
      photon_array.descriptors[index] == static_cast<unsigned char>(Constants::PASS) ||
      photon_array.descriptors[index] == static_cast<unsigned char>(Constants::EXIT);
      }), active_photons.end());
    active_count  = active_photons.size();
    // std::cout << " tmp track  " << tmpAbsE[7] << std::endl;
  }
  // std::cout << " counter  " << counter << std::endl;
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void event_transport_photon(const uint32_t rank_cell_offset, std::vector<Photon> &photon_array, const Cell *cells, Cell_Tally *cell_tallies, size_t cellSize, const std::vector<EmissionGroupData> &emission_groups) {

  using Constants::bc_type;
  using Constants::c;
  // events
  using std::min;
  const size_t maxPhotons = photon_array.size();

  std::vector<size_t> scatter_events(maxPhotons), boundary_events(maxPhotons), census_events(maxPhotons), killed_events(maxPhotons), active_photons(maxPhotons);
  std::vector<Event> events(maxPhotons);
  std::vector<uint32_t> local_cell_indices(maxPhotons);
  std::vector<double> sigma_s(maxPhotons), sigma_a(maxPhotons), f(maxPhotons), total_sigma_s(maxPhotons);
  
  std::iota(active_photons.begin(), active_photons.end(), 0); 

  size_t active_count = maxPhotons;
  size_t counter = 0;

  while (active_count > 0) {
    size_t scatter_count = 0, boundary_count = 0, census_count = 0, killed_count = 0;

    precompute_data(rank_cell_offset, photon_array, cells, active_photons.data(), active_count, sigma_s, sigma_a, f, total_sigma_s, local_cell_indices);
    calculate_distances(photon_array, cells, active_photons.data(), active_count, total_sigma_s, local_cell_indices, events);
    update_photon_state(photon_array, cell_tallies, events, sigma_a, f , local_cell_indices, active_count, scatter_events, boundary_events, census_events, killed_events, scatter_count, boundary_count, census_count, killed_count, cellSize);

    process_scatter_events(photon_array, sigma_s, total_sigma_s, local_cell_indices, cells, scatter_events.data(), scatter_count, counter, emission_groups);
    process_census_events(photon_array, census_events.data(), census_count);
    process_boundary_events(photon_array, boundary_events.data(), cells, rank_cell_offset, boundary_count);
    process_killed_events(photon_array, killed_events.data(), cells,  cell_tallies, rank_cell_offset, killed_count);

    active_photons.erase(std::remove_if(active_photons.begin(), active_photons.end(), [&photon_array](size_t index){
      return photon_array[index].get_descriptor() == static_cast<unsigned char>(Constants::KILLED) ||
      photon_array[index].get_descriptor() == static_cast<unsigned char>(Constants::CENSUS) ||
      photon_array[index].get_descriptor() == static_cast<unsigned char>(Constants::PASS) ||
      photon_array[index].get_descriptor() == static_cast<unsigned char>(Constants::EXIT);
      }), active_photons.end());
    active_count  = active_photons.size();
    // std::cout << " tmp track  " << tmpAbsE[7] << std::endl;
  }
  // std::cout << " counter  " << counter << std::endl;
}
//----------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
void cpu_event_transport_photons(const uint32_t rank_cell_offset,
    PhotonArray &photon_array, const std::vector<Cell> &cells, std::vector<Cell_Tally> &cell_tallies, int n_omp_threads, const std::vector<EmissionGroupData> &emission_groups) {

  auto cpu_cells_ptr{cells.data()};
  const auto n_cells = cell_tallies.size();
  
  event_transport_photon(rank_cell_offset, photon_array, cells.data(), cell_tallies.data(), cells.size(), emission_groups);
}

//------------------------------------------------------------------------------------------------//
void cpu_event_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> &photon_array, const std::vector<Cell> &cells, std::vector<Cell_Tally> &cell_tallies, int n_omp_threads, const std::vector<EmissionGroupData> &emission_groups) {

  auto cpu_cells_ptr{cells.data()};
  const auto n_cells = cell_tallies.size();
  
  event_transport_photon(rank_cell_offset, photon_array, cells.data(), cell_tallies.data(), cells.size(), emission_groups);
}
#endif // def event_based_transport_h_
//----------------------------------------------------------------------------//
// end of event_based_transport.h
//----------------------------------------------------------------------------//
