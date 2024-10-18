
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   history_based_transport.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  IMC transport with particle passing method
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef history_based_transport_h_
#define history_based_transport_h_

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "config.h"
#include "RNG.h"
#include "cell_tally.h"
#include "constants.h"
#include "photon.h"
#include "photon_array.h"
#include "sampling_functions.h"

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void transport_photon(const uint32_t rank_cell_offset,
    Photon &phtn, const Cell *cells, Cell_Tally *cell_tallies) {

  using Constants::bc_type;
  using Constants::c;
  // events
  using std::min;

  auto &rng = phtn.get_rng();

  uint32_t surface_cross = 0;

  uint32_t local_cell_index =  phtn.get_cell() - rank_cell_offset;
  Cell const * cell = &cells[local_cell_index];
  bool active = true;

  // keep a thread local copy of these tallies and do the atomic add when the particle leaves the
  // cell or it's otherwise terminated (try to reduce atomic contention with post-move tally)
  double thread_absorbed_E{0.0};
  double thread_track_E{0.0};

  // transport this photon
  while (active) {
    const double sigma_s = cell->get_op_s(phtn.get_group());
    const double sigma_a = cell->get_op_a(phtn.get_group());
    const double f = cell->get_f();
    const double total_sigma_s = (1.0 - f) * sigma_a + sigma_s;

    // get distance to event
    const double dist_to_scatter = (total_sigma_s > 0.0) ?
      -log(rng.generate_random_number()) / total_sigma_s : 1.0e100;

    const double dist_to_boundary = cell->get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    const double dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    const double dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    const double absorbed_E = phtn.get_E() * (1.0 - exp(-sigma_a * f * dist_to_event));

    thread_absorbed_E += absorbed_E;
    thread_track_E += absorbed_E / (sigma_a * f);

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(Constants::cutoff_fraction)) {
      thread_absorbed_E += phtn.get_E();
      cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
      cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
      active = false;
      phtn.set_descriptor(Constants::KILLED);
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        phtn.set_angle(get_uniform_angle(rng));
        if (rng.generate_random_number() > (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, *cell));
        phtn.set_descriptor(Constants::SCATTER);
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        auto boundary_event = cell->get_bc(surface_cross);
        if (boundary_event == Constants::ELEMENT) {
          // dump thread energy into this cell's indexi before updating it
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
          // update photon's cell index
          phtn.set_cell(cell->get_next_cell(surface_cross));
          local_cell_index =  phtn.get_cell() - rank_cell_offset;
          cell = &cells[local_cell_index]; // note: only for on rank mesh data
          phtn.set_descriptor(Constants::BOUND);
          thread_absorbed_E = 0.0;
          thread_track_E = 0.0;
        } else if (boundary_event == Constants::PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          phtn.set_cell(cell->get_next_cell(surface_cross));
          phtn.set_descriptor(Constants::PASS);
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
        } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
          active = false;
          phtn.set_descriptor(Constants::EXIT);
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
        } else {
          phtn.reflect(surface_cross);
          phtn.set_descriptor(Constants::BOUND);
        }
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        active = false;
        phtn.set_descriptor(Constants::CENSUS);
        cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
        cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
      }
    } // end event loop
  } // end while alive
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
GPU_HOST_DEVICE
void transport_photon(const uint32_t rank_cell_offset,
    PhotonArray &phtns, const size_t i, const Cell *cells, Cell_Tally *cell_tallies) {

  using Constants::bc_type;
  using Constants::c;
  // events
  using std::min;

  auto &rng = phtns.rng[i];

  uint32_t surface_cross = 0;

  uint32_t local_cell_index =  phtns.cell_ID[i] - rank_cell_offset;
  Cell const * cell = &cells[local_cell_index];
  bool active = true;

  // keep a thread local copy of these tallies and do the atomic add when the particle leaves the
  // cell or it's otherwise terminated (try to reduce atomic contention with post-move tally)
  double thread_absorbed_E{0.0};
  double thread_track_E{0.0};

  // transport this photon
  while (active) {
    const double sigma_s = cell->get_op_s(phtns.group[i]);
    const double sigma_a = cell->get_op_a(phtns.group[i]);
    const double f = cell->get_f();
    const double total_sigma_s = (1.0 - f) * sigma_a + sigma_s;

    // get distance to event
    const double dist_to_scatter = (total_sigma_s > 0.0) ?
      -log(rng.generate_random_number()) / total_sigma_s : 1.0e100;

    const double dist_to_boundary = cell->get_distance_to_boundary(
        phtns.pos[i], phtns.angle[i], surface_cross);
    const double dist_to_census = phtns.life_dx[i];

    // select minimum distance event
    const double dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    const double absorbed_E = phtns.E[i] * (1.0 - exp(-sigma_a * f * dist_to_event));

    thread_absorbed_E += absorbed_E;
    thread_track_E += absorbed_E / (sigma_a * f);

    phtns.E[i] = (phtns.E[i] - absorbed_E);

    // update position
    phtns.pos[i][0] += phtns.angle[i][0] * dist_to_event; 
    phtns.pos[i][1] += phtns.angle[i][1] * dist_to_event; 
    phtns.pos[i][2] += phtns.angle[i][2] * dist_to_event; 
    phtns.life_dx[i] -= dist_to_event; 

    // apply runtime reduction
    if (phtns.E[i] / phtns.E0[i] < Constants::cutoff_fraction) {
      thread_absorbed_E += phtns.E[i];
      cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
      cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
      active = false;
      phtns.descriptors[i] = Constants::KILLED;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        phtns.angle[i] = get_uniform_angle(rng);
        if (rng.generate_random_number() > (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtns.group[i] = sample_emission_group(rng, *cell);
        phtns.descriptors[i] = Constants::SCATTER;
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        auto boundary_event = cell->get_bc(surface_cross);
        if (boundary_event == Constants::ELEMENT) {
          // dump thread energy into this cell's indexi before updating it
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
          // update photon's cell index
          phtns.cell_ID[i] = cell->get_next_cell(surface_cross);
          local_cell_index =  phtns.cell_ID[i] - rank_cell_offset;
          cell = &cells[local_cell_index]; // note: only for on rank mesh data
          phtns.descriptors[i] = Constants::BOUND;
          thread_absorbed_E = 0.0;
          thread_track_E = 0.0;
        } else if (boundary_event == Constants::PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          phtns.cell_ID[i] = cell->get_next_cell(surface_cross);
          phtns.descriptors[i] = Constants::PASS;
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
        } else if (boundary_event == Constants::VACUUM || boundary_event == Constants::SOURCE) {
          active = false;
          phtns.descriptors[i] = Constants::EXIT;
          cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
          cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
        } else {
          int reflect_angle = surface_cross/2; // X -> 0, Y->1, Z->2
          phtns.angle[i][reflect_angle] = -phtns.angle[i][reflect_angle];
          phtns.descriptors[i] = Constants::BOUND;
        }
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        active = false;
        phtns.descriptors[i] = Constants::CENSUS; 
        cell_tallies[local_cell_index].accumulate_absorbed_E(thread_absorbed_E);
        cell_tallies[local_cell_index].accumulate_track_E(thread_track_E);
      }
    } // end event loop
  } // end while alive
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
void history_cpu_transport_photons(const uint32_t rank_cell_offset,
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
void history_cpu_transport_photons(const uint32_t rank_cell_offset,
    PhotonArray &photons, const std::vector<Cell> &cells, std::vector<Cell_Tally> &cell_tallies, int n_omp_threads) {

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
      transport_photon(rank_cell_offset, photons, i, cpu_cells_ptr, thread_tally_ptr);
    }
  } // end parallel region

  // reduce tallies if using openmp
  for(size_t cell=0; cell<n_cells; ++cell) {
    auto &cell_tally{cell_tallies[cell]};
    for(int thread =0; thread<n_omp_threads;++thread) {
      cell_tally.merge_in_tally(thread_tallies[thread][cell]);
    }
  }
#else
  // normal serial version
  for (size_t i=0; photons.size(); ++i) {
    transport_photon(rank_cell_offset, photons, i, cpu_cells_ptr, cell_tallies.data());
  }
#endif
}
//------------------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------------------//
void gpu_transport_photons(const uint32_t rank_cell_offset,
    std::vector<Photon> &cpu_photons, const Cell *device_cells_ptr, std::vector<Cell_Tally> &cpu_cell_tallies) {

#ifdef USE_CUDA
  uint32_t n_batch_photons = static_cast<uint32_t>(cpu_photons.size());
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

  // allocate and copy photons
  Photon *device_photons_ptr;
  cudaError_t err = cudaMalloc((void **)&device_photons_ptr, sizeof(Photon) * cpu_photons.size());
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
  int n_blocks = (n_batch_photons + Constants::n_threads_per_block - 1) /
                 Constants::n_threads_per_block;

  cudaDeviceSynchronize();

  std::cout << "Launching with " << n_blocks << " blocks and ";
  std::cout << n_batch_photons << " photons" << std::endl;
  gpu_no_accel_transport<<<n_blocks, Constants::n_threads_per_block>>>(
      rank_cell_offset, device_photons_ptr, device_cells_ptr, device_cell_tallies_ptr, n_batch_photons);


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

  // free device pointers for photons and cell tallies
  cudaFree(device_photons_ptr);
  cudaFree(device_cell_tallies_ptr);

#endif
}

#endif // def history_based_transport_h_
//----------------------------------------------------------------------------//
// end of history_based_transport.h
//----------------------------------------------------------------------------//
