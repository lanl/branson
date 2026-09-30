//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   particle_pass_transport.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  IMC transport with particle passing method
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 *         Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 */
//---------------------------------------------------------------------------//

#ifndef transport_photon_h_
#define transport_photon_h_

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "config.h"
#include "RNG.h"
#include "cell_tally.h"
#include "constants.h"
#include "gpu_setup.h"
#include "photon.h"
#include "sampling_functions.h"

void post_process_photons(const double next_dt, std::vector<Photon> &all_photons, std::vector<Photon> &census_list, double &census_E, double &exit_E) {
  for ( auto & phtn : all_photons) {
    auto descriptor{phtn.get_descriptor()};
    switch (descriptor) {
    case Constants::event_type::PASS:
      // handle in other function
      break;
    case Constants::event_type::KILLED:
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

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void transport_photon(const uint32_t rank_cell_offset,
    Photon &phtn, const Cell *cells, Cell_Tally *cell_tallies) {

  using std::min;

  auto &rng = phtn.get_rng();
  uint32_t surface_cross = 0;
  uint32_t local_cell_index =  phtn.get_cell() - rank_cell_offset;
  Cell const * cell = &cells[local_cell_index];

  // Keep the frequently changed photon fields local until transport finishes.
  double E = phtn.get_E();
  const double cutoff = Constants::cutoff_fraction * phtn.get_E0();
  uint32_t group = phtn.get_group();
  uint32_t photon_cell = phtn.get_cell();
  auto angle = phtn.get_angle();
  auto pos = phtn.get_position();
  double life_dx = phtn.get_distance_remaining();
  Constants::event_type descriptor = Constants::CENSUS;

  double absorbed_tally = 0.0;
  double track_tally = 0.0;
  double sigma_s = 0.0;
  double sigma_a = 0.0;
  double f = 0.0;
  bool refresh_opacity = true;
  bool active = true;

  // transport this photon
  while (active) {
    if (refresh_opacity) {
      sigma_s = cell->get_op_s(group);
      sigma_a = cell->get_op_a(group);
      f = cell->get_f();
      refresh_opacity = false;
    }
    const double total_sigma_s = (1.0 - f) * sigma_a + sigma_s;
    const double inv_total_sigma_s = (total_sigma_s > 0.0) ? 1.0 / total_sigma_s : 0.0;

    // get distance to event
    const double dist_to_scatter = (total_sigma_s > 0.0) ?
      -log(rng.generate_random_number()) * inv_total_sigma_s : 1.0e100;

    const double dist_to_boundary = cell->get_distance_to_boundary(
        pos, angle, surface_cross);
    const double dist_to_census = life_dx;

    // select minimum distance event
    const double dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    const double absorbed_E = E * (1.0 - exp(-sigma_a * f * dist_to_event));
    absorbed_tally += absorbed_E;
    track_tally += (sigma_a * f > 0.0) ? absorbed_E / (sigma_a * f)
                                        : E * dist_to_event;
    E -= absorbed_E;

    // update position
    phtn.move(dist_to_event, angle, pos, life_dx);

    // apply variance/runtime reduction
    if (E < cutoff) {
      absorbed_tally += E;
      active = false;
      descriptor = Constants::KILLED;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        angle = get_uniform_angle(rng);
        if (rng.generate_random_number() > sigma_s * inv_total_sigma_s) {
          group = sample_emission_group(rng, *cell);
          refresh_opacity = true;
        }
        descriptor = Constants::SCATTER;
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        auto boundary_event = cell->get_bc(surface_cross);
        if (boundary_event == Constants::ELEMENT) {
          // Commit this cell's tally before advancing to the next cell.
          cell_tallies[local_cell_index].accumulate_absorbed_E(absorbed_tally);
          cell_tallies[local_cell_index].accumulate_track_E(track_tally);
          absorbed_tally = 0.0;
          track_tally = 0.0;
          photon_cell = cell->get_next_cell(surface_cross);
          local_cell_index = photon_cell - rank_cell_offset;
          cell = &cells[local_cell_index]; // note: only for on rank mesh data
          refresh_opacity = true;
          descriptor = Constants::BOUND;
        } else if (boundary_event == Constants::PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          photon_cell = cell->get_next_cell(surface_cross);
          descriptor = Constants::PASS;
        } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
          active = false;
          descriptor = Constants::EXIT;
        } else {
          phtn.reflect(surface_cross, angle);
          descriptor = Constants::BOUND;
        }
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        active = false;
        descriptor = Constants::CENSUS;
      }
    } // end event loop
  } // end while alive

  cell_tallies[local_cell_index].accumulate_absorbed_E(absorbed_tally);
  cell_tallies[local_cell_index].accumulate_track_E(track_tally);

  phtn.set_E(E);
  phtn.set_descriptor(descriptor);
  phtn.set_group(group);
  phtn.set_angle(angle);
  phtn.set_position(pos);
  phtn.set_distance_to_census(life_dx);
  phtn.set_cell(photon_cell);
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
GPU_KERNEL
void gpu_no_accel_transport(const uint32_t rank_cell_offset,
    Photon *all_photons, const Cell *cells, Cell_Tally *cell_tallies, const uint32_t n_batch_particles) {

#ifdef USE_CUDA
  int32_t particle_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (particle_id < n_batch_particles) {
    transport_photon(rank_cell_offset, all_photons[particle_id], cells, cell_tallies);
  } // if particle id is valid
  __syncthreads();

#endif
}
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
void cpu_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> &photons, const std::vector<Cell> &cells, std::vector<Cell_Tally> &cell_tallies, int n_omp_threads) {

  auto cpu_cells_ptr{cells.data()};
  const auto n_cells = cell_tallies.size();
#ifdef USE_OPENMP
  // this is set earlier based on input variable
  std::vector<std::vector<Cell_Tally>> thread_tallies(n_omp_threads);
#pragma omp parallel
  {
    thread_tallies[omp_get_thread_num()].resize(n_cells);
    auto thread_tally_ptr = thread_tallies[omp_get_thread_num()].data();
#pragma omp for schedule(guided)
    for (int i=0; i<photons.size(); ++i) {
      transport_photon(rank_cell_offset, photons[i], cpu_cells_ptr, thread_tally_ptr);
    }
  } // end parallel region

  // reduce tallies if using openmp
  for(size_t cell=0; cell<n_cells; ++cell) {
    auto &cell_tally{cell_tallies[cell]};
    for(int thread =0; thread<n_omp_threads;++thread)
      cell_tally.merge_in_tally(thread_tallies[thread][cell]);
  }
#else
  // normal serial version
  for (auto &photon : photons)
    transport_photon(rank_cell_offset, photon, cpu_cells_ptr, cell_tallies.data());
#endif
}
//------------------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------------------//
void gpu_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> &cpu_photons, GPU_Setup &gpu_setup, std::vector<Cell_Tally> &cpu_cell_tallies) {

#ifdef USE_CUDA
  uint32_t n_batch_photons = static_cast<uint32_t>(cpu_photons.size());
  if (n_batch_photons == 0)
    return;
#ifdef ENABLE_VERBOSE_GPU_TRANSPORT
  int my_bus_id = 0;
  int my_device = 0;
  cudaGetDevice(&my_device);
  cudaDeviceGetAttribute(&my_bus_id, cudaDevAttrPciBusId, my_device);
  std::cout << "--GPU data-- device: " << my_device << ", busID: " << my_bus_id << ", ";
  std::cout << "particle bytes for this batch: ";
  std::cout << n_batch_photons * (sizeof(Photon) + sizeof(int)) << std::endl;
#endif

  // use this vector to track active indices in the particle vector
  std::vector<int> active_indices(n_batch_photons);
  std::iota(active_indices.begin(), active_indices.end(), 0);

  gpu_setup.reserve_transport_buffers(n_batch_photons, cpu_cell_tallies.size());
  Photon *device_photons_ptr = gpu_setup.get_device_photons_ptr();
  Cell_Tally *device_cell_tallies_ptr = gpu_setup.get_device_cell_tallies_ptr();

  // Refresh batch input in the reusable device buffers.
  cudaError_t err = cudaMemcpy(device_photons_ptr, cpu_photons.data(), sizeof(Photon) * n_batch_photons,
                   cudaMemcpyHostToDevice);
  Insist(!err, "CUDA error in copying photons data");

  err = cudaMemcpy(device_cell_tallies_ptr, cpu_cell_tallies.data(), sizeof(Cell_Tally) * cpu_cell_tallies.size(),
                   cudaMemcpyHostToDevice);
  Insist(!err, "CUDA error in copying cell tallies data");

  // kernel settings
  int n_blocks = (n_batch_photons + Constants::n_threads_per_block - 1) /
                 Constants::n_threads_per_block;

  cudaDeviceSynchronize();

  std::cout << "Launching with " << n_blocks << " blocks and ";
  std::cout << n_batch_photons << " photons" << std::endl;
  gpu_no_accel_transport<<<n_blocks, Constants::n_threads_per_block>>>(
      rank_cell_offset, device_photons_ptr, gpu_setup.get_device_cells_ptr(), device_cell_tallies_ptr, n_batch_photons);


  Insist(!(cudaGetLastError()), "CUDA error in transport kernel launch");
  cudaDeviceSynchronize();

  // copy particles back to host
  err = cudaMemcpy(cpu_photons.data(), device_photons_ptr, n_batch_photons * sizeof(Photon),
                   cudaMemcpyDeviceToHost);
  Insist(!err, "CUDA error in copying photons back to host");

  // copy cell tallies back to host
  err = cudaMemcpy(cpu_cell_tallies.data(), device_cell_tallies_ptr, sizeof(Cell_Tally)*cpu_cell_tallies.size(),
                   cudaMemcpyDeviceToHost);
  Insist(!err, "CUDA error in copying cell tallies back to host");

#endif
}

#endif // def transport_photon_h_
//----------------------------------------------------------------------------//
// end of transport_photon.h
//----------------------------------------------------------------------------//
