//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   particle_pass_transport.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  IMC transport with particle passing method
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_photon_h_
#define transport_photon_h_

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "RNG.h"
#include "cell_tally.h"
#include "constants.h"
#include "photon.h"
#include "sampling_functions.h"

void post_process_photons(const double next_dt, std::vector<Photon> &all_photons, std::vector<Photon> &census_list, double &census_E, double &exit_E) {
  for ( auto & phtn : all_photons) {
    auto descriptor{phtn.get_descriptor()};
    switch (descriptor) {
    case Constants::event_type::PASS:
      // handle in other function
      break;
    case Constants::event_type::KILL:
      // note: for now killed particles go into the material so separate conservation issues here
      break;
    case Constants::event_type::EXIT:
      exit_E+=phtn.get_E();
      break;
    case Constants::event_type::CENSUS:
      phtn.set_distance_to_census(Constants::c*next_dt);
      census_list.push_back(phtn);
      census_E+=phtn.get_E();
      break;
    } //switch(descriptor)
  } // phtn : all_photons
}


GPU_KERNEL
void gpu_no_accel_transport(const uint32_t rank_cell_offset,
    Photon &phtn, const Cell *cells, RNG *rng, Cell_Tally *cell_tallies) {

#ifdef USE_CUDA
  __shared__ Particle block_particles[rtt_tracking::n_threads_per_block];
  int32_t particle_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (particle_id < n_batch_particles) {
    block_particles[threadIdx.x] = particles[index_set[particle_id]];

    transport_photon(rank_cell_offest, block_particles[threadIdx.x], cells, rng, cell_tallies);

    particles[index_set[particle_id]] = block_particles[threadIdx.x];
  }

#endif
}

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void transport_photon(const uint32_t rank_cell_offset,
    Photon &phtn, const Cell *cells, RNG *rng, Cell_Tally *cell_tallies) {

  using Constants::bc_type;
  using Constants::c;
  // events
  using std::min;

  uint32_t next_cell;
  bc_type boundary_event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  int group;

  uint32_t surface_cross = 0;
  double cutoff_fraction = 0.01; // note: get this from IMC_state

  uint32_t global_cell_index = phtn.get_cell();
  uint32_t local_cell_index =  global_cell_index - rank_cell_offset;
  Cell const * cell = &cells[local_cell_index];
  bool active = true;

  // transport this photon
  while (active) {
    group = phtn.get_group();
    sigma_a = cell->get_op_a(group);
    sigma_s = cell->get_op_s(group);
    f = cell->get_f();

    // get distance to event
    const double total_sigma_s = ((1.0 - f) * sigma_a + sigma_s);;
    dist_to_scatter = (total_sigma_s > 0.0) ?
      -log(rng->generate_random_number()) / total_sigma_s : 1.0e100;

    dist_to_boundary = cell->get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);

    cell_tallies[local_cell_index].accumulate_absorbed_E(absorbed_E);
    cell_tallies[local_cell_index].accumulate_track_E(absorbed_E / (sigma_a * f));

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      cell_tallies[local_cell_index].accumulate_absorbed_E(phtn.get_E());
      active = false;
      phtn.set_descriptor(Constants::KILL);
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {

        phtn.set_angle(get_uniform_angle(rng));
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, *cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell->get_bc(surface_cross);
        if (boundary_event == Constants::ELEMENT) {
          next_cell = cell->get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          global_cell_index = next_cell;
          local_cell_index =  global_cell_index - rank_cell_offset;
          cell = &cells[local_cell_index]; // note: only for on rank mesh data
        } else if (boundary_event == Constants::PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          next_cell = cell->get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          phtn.set_descriptor(Constants::PASS);
        } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
          active = false;
          phtn.set_descriptor(Constants::EXIT);
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        active = false;
        phtn.set_descriptor(Constants::CENSUS);
      }
    } // end event loop
  }   // end while alive
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
void cpu_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> &photons, const std::vector<Cell> &cells, RNG *rng, std::vector<Cell_Tally> &cell_tallies) {

  auto cpu_cells_ptr{cells.data()};
  auto cpu_cell_tallies_ptr{cell_tallies.data()};
  for (auto &photon : photons) {
    transport_photon(rank_cell_offset, photon, cpu_cells_ptr, rng, cpu_cell_tallies_ptr);
  }
}

void gpu_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> cpu_photons, const Cell *device_cells_ptr, RNG *rng, std::vector<Cell_Tally> cpu_cell_tallies) {

#ifdef USE_CUDA
  size_t n_active_photons = cpu_photons.size();
  // allocate and copy photons
  Photon *device_photons_ptr;
  err = cudaMalloc((void **)&device_photons_ptr, sizeof(Photon) * cpu_photons.size());
  Insist(!err, "CUDA error in allocating photons data");
  err = cudaMemcpy(device_photons_ptr, cpu_photons.data(), sizeof(Photon) * cpu_photons.size(),
                   cudaMemcpyHostToDevice);
  Insist(!err, "CUDA error in copying photons data");

  // allocate and copy cell tally object
  Cell_Tally *device_cell_tallies_ptr;
  err = cudaMalloc((void **)&device_cell_tallies_ptr, sizeof(Cell_Tally) * cpu_cell_tallies.size());
  Insist(!err, "CUDA error in allocating cell tallies data");
  err = cudaMemcpy(device_cell_tallies_ptr, cpu_cell_tallies.data(), sizeof(Cell_Tally) * cpu_cell_tallies.size(),
                   cudaMemcpyHostToDevice);
  Insist(!err, "CUDA error in copying cell tallies data");

  // kernel settings
  int n_blocks = (n_active_photons + Constants::n_threads_per_block - 1) /
                 Constants::n_threads_per_block;

  cudaDeviceSynchronize();

  std::cout << "Launching with " << n_blocks << " blocks and ";
  std::cout << n_active_photons << " photons" << std::endl;
  gpu_no_accel_transport<<<n_blocks, Constants::n_threads_per_block>>>(
      rank_cell_offset, n_active_photons, device_photons_ptr, device_cells_ptr, RNG, device_cell_tallies_ptr);

  Insist(!(cudaGetLastError()), "CUDA error in transport kernel launch");
  cudaDeviceSynchronize();

  // copy particles back to host
  err = cudaMemcpy(cpu_photons.data(), device_photons_ptr, n_active_photons * sizeof(Photon),
                   cudaMemcpyDeviceToHost);
  Insist(!err, "CUDA error in copying photons back to host");

  // copy cell tallies back to host
  err = cudaMemcpy(cpu_cell_tallies.data(), device_cell_tallies_ptr, sizeof(Cell_Tally)*cpu_cell_tallies.size(),,
                   cudaMemcpyDeviceToHost);
  Insist(!err, "CUDA error in copying cell tallies back to host");

#endif
}

#endif // def transport_photon_h_
//----------------------------------------------------------------------------//
// end of transport_photon.h
//----------------------------------------------------------------------------//
