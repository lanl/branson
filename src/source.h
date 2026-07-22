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

#include "timer.h"
#include "gpu_setup.h"
#include "constants.h"
#include "cell.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"

#ifdef caliper_FOUND
#include <caliper/cali.h>
#include <adiak.hpp>
#endif

GPU_KERNEL void make_source_photons( Cell  const * const cells,  const double dt, const uint32_t seed, uint64_t const * const photon_stream_numbers, double const * const photon_E, int const * const photon_type,  int const * const photon_source_face,  uint32_t const * const photon_cell_index, const uint64_t n_photons, Photon * const all_photons) {
  using Constants::c;
#ifdef USE_GPU
  int32_t i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n_photons) {
#else
  for (size_t i=0;i<n_photons;++i) {
#endif
    auto source_type = photon_type[i];
    RNG rng(seed, photon_stream_numbers[i]);
    // census and emission get uniform position and isotropic angle, source gets position on face and
    // angle surface normal distribution
    auto &cell = cells[photon_cell_index[i]];
    auto pos = (source_type != 2) ? get_uniform_position_in_cell(cell, rng) : get_uniform_position_on_face(cell, rng, photon_source_face[i]);
    auto angle = (source_type != 2) ? get_uniform_angle(rng) : get_source_angle_on_face(rng, photon_source_face[i]);
    double distance_to_census = (source_type ==0) ? c*dt : rng.generate_random_number() * c * dt;
    uint32_t group = std::floor(rng.generate_random_number() * double(BRANSON_N_GROUPS));

    all_photons[i] = Photon(cell.get_global_index(), group, source_type, Constants::event_type::BORN_SOURCE,  pos, angle, photon_E[i], distance_to_census, rng);
  }
#ifdef USE_GPU
  __syncthreads();
#endif
}

// for SoA data structures, cell ID, E, E0 and source_type are already set, use them to sample
// position and angle
GPU_KERNEL void set_source_photons( Cell  const * const cells,  const double dt, const uint32_t seed, uint64_t const * const photon_stream_numbers, int const * const photon_type, int const * const photon_source_face,  const uint64_t n_photons, uint32_t *photon_cell_index, RNG *rng, std::array< double,3 > *pos, std::array<double,3> *angle, double *life_dx, uint32_t *group) {
  using Constants::c;
#ifdef USE_GPU
  int32_t i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n_photons) {
#else
  for (size_t i=0;i<n_photons;++i) {
#endif
    auto source_type = photon_type[i];
    rng[i] = RNG(seed, photon_stream_numbers[i]);
    // census and emission get uniform position and isotropic angle, source gets position on face and
    // angle surface normal distribution
    auto &cell = cells[photon_cell_index[i]];
    pos[i] = (source_type != 2) ? get_uniform_position_in_cell(cell, rng[i]) : get_uniform_position_on_face(cell, rng[i], photon_source_face[i]);
    angle[i] = (source_type != 2) ? get_uniform_angle(rng[i]) : get_source_angle_on_face(rng[i], photon_source_face[i]);
    life_dx[i] = (source_type ==0) ? c*dt : rng[i].generate_random_number() * c * dt;
    group[i] = std::floor(rng[i].generate_random_number() * double(BRANSON_N_GROUPS));
    // cell index comes in as local, is then set to global
    photon_cell_index[i] = cell.get_global_index();
  }
#ifdef USE_GPU
  __syncthreads();
#endif
}

// template function for making photons
template <typename Census_T>
void make_photons(const double dt, const Mesh &mesh, const int rank, const uint32_t cycle,
                    const uint32_t seed, const uint64_t n_user_photons, const double total_E, GPU_Setup<Census_T> &gpu_setup) {
#ifdef caliper_FOUND
  CALI_MARK_BEGIN(make_photons);
#endif
  bool make_initial_census_flag{cycle==1};
  auto E_cell_census = mesh.get_census_E();
  auto E_cell_emission = mesh.get_emission_E();
  auto E_cell_source = mesh.get_source_E();
  // for RNG offsets, each cycle allows for one hundred million photons across one hundred
  // thousand ranks, increment the ten trillon place for the next cycle, using cycle plus one for
  // the cycle offset gives the initial census their own space
  const uint64_t census_rank_stream_num_offset{n_user_photons * rank};
  const uint64_t cycle_stream_num_offset{10000000000000UL * static_cast<uint64_t>(cycle)};
  const uint64_t rank_stream_num_offset{n_user_photons * static_cast<uint64_t>(rank)};

  // figure out how many to make to size all_photons vector
  uint64_t n_photons = 0;
  for (auto const &cell : mesh) {
    int i = mesh.get_local_index(cell.get_global_index());
    // initial census
    if (make_initial_census_flag && E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0)
        t_num_census = 1;
      n_photons+=t_num_census;
    }
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission =
          int(n_user_photons * E_cell_emission[i] / total_E);
      // make at least one photon to represent emission energy
      if (t_num_emission == 0)
        t_num_emission = 1;
      n_photons+=t_num_emission;
    }
    if (E_cell_source[i] > 0.0) {
      // boundary source
        uint32_t t_num_source =
            int(n_user_photons * E_cell_source[i] / total_E);
        // make at least one photon to represent source energy
        if (t_num_source == 0)
          t_num_source = 1;
      n_photons+=t_num_source;
    }
  }

  std::vector<uint64_t> photon_stream_nums(n_photons);
  std::vector<double> photon_E(n_photons);
  std::vector<int> photon_type(n_photons);
  std::vector<int> photon_source_face(n_photons);
  std::vector<uint32_t> photon_cell_index(n_photons);

  // in serial loop through and set the seed for each photon
  // use this to increment the seed for each photon
  uint64_t ith_photon{0UL};

  for (auto const &cell : mesh) {
    uint32_t i = mesh.get_local_index(cell.get_global_index());
    // initial census
    if (make_initial_census_flag && E_cell_census[i] > 0.0) {
      uint32_t t_num_census = int(n_user_photons * E_cell_census[i] / total_E);
      // make at least one photon to represent census energy
      if (t_num_census == 0) t_num_census = 1;
      const double photon_census_E = E_cell_census[i] / t_num_census;
      for (uint32_t p=0; p<t_num_census;++p) {
        photon_stream_nums[ith_photon] = census_rank_stream_num_offset+ith_photon;
        photon_type[ith_photon] = 0; // census type
        photon_E[ith_photon] = photon_census_E;
        photon_cell_index[ith_photon] = i;
        ith_photon++;
      }
    }
    // emission
    if (E_cell_emission[i] > 0.0) {
      uint32_t t_num_emission =
          int(n_user_photons * E_cell_emission[i] / total_E);
      // make at least one photon to represent emission energy
      if (t_num_emission == 0) t_num_emission = 1;
      const double photon_emission_E = E_cell_emission[i] / t_num_emission;
      for (uint32_t p=0; p<t_num_emission;++p) {
        photon_stream_nums[ith_photon] = cycle_stream_num_offset + rank_stream_num_offset+ith_photon;
        photon_type[ith_photon] = 2;
        photon_E[ith_photon] = photon_emission_E;
        photon_cell_index[ith_photon] = i;
        ith_photon++;
      }
    }
    if (E_cell_source[i] > 0.0) {
      // boundary source
      uint32_t t_num_source =
          int(n_user_photons * E_cell_source[i] / total_E);
      // make at least one photon to represent source energy
      if (t_num_source == 1) t_num_source = 1;
      const double photon_source_E = E_cell_source[i] / t_num_source;
      const int face = cell.get_source_face();
      for (uint32_t p=0; p<t_num_source;++p) {
        photon_stream_nums[ith_photon] = cycle_stream_num_offset + rank_stream_num_offset+ith_photon;
        photon_type[ith_photon] = 1;
        photon_source_face[ith_photon] = face;
        photon_E[ith_photon] = photon_source_E;
        photon_cell_index[ith_photon] = i;
        ith_photon++;
      }
    }
  }

  #ifdef USE_GPU
  // Copy input data to device
  uint64_t *device_photon_stream_nums_ptr;
  auto alloc_err = cudaMalloc((void **)&device_photon_stream_nums_ptr, sizeof(uint64_t) * n_photons);
  Insist(!alloc_err, "CUDA/HIP error allocating photon seeds");
  auto copy_err = cudaMemcpy(device_photon_stream_nums_ptr, photon_stream_nums.data(), sizeof(uint64_t) * n_photons, cudaMemcpyHostToDevice);
  Insist(!copy_err, "CUDA/HIP error copying photon seeds to device");

  double *device_photon_E_ptr;
  alloc_err = cudaMalloc((void **)&device_photon_E_ptr, sizeof(double) * n_photons);
  Insist(!alloc_err, "CUDA/HIP error allocating photon E");
  copy_err = cudaMemcpy(device_photon_E_ptr, photon_E.data(), sizeof(double) * n_photons, cudaMemcpyHostToDevice);
  Insist(!copy_err, "CUDA/HIP error copying photon E to device");

  int  *device_source_type_ptr;
  alloc_err = cudaMalloc((void **)&device_source_type_ptr, sizeof(int) * n_photons);
  Insist(!alloc_err, "CUDA/HIP error allocating photon types");
  copy_err = cudaMemcpy(device_source_type_ptr, photon_type.data(), sizeof(int) * n_photons, cudaMemcpyHostToDevice);
  Insist(!copy_err, "CUDA/HIP error copying photon_types to device");

  int  *device_photon_source_face_ptr;
  alloc_err = cudaMalloc((void **)&device_photon_source_face_ptr, sizeof(int) * n_photons);
  Insist(!alloc_err, "CUDA/HIP error allocating photon source face");
  copy_err = cudaMemcpy(device_photon_source_face_ptr, photon_source_face.data(), sizeof(int) * n_photons, cudaMemcpyHostToDevice);
  Insist(!copy_err, "CUDA/HIP error copying photon source face to device");

  uint32_t  *device_cell_index_ptr;
  alloc_err = cudaMalloc((void **)&device_cell_index_ptr, sizeof(uint32_t) * n_photons);
  Insist(!alloc_err, "CUDA/HIP error allocating photon cell index");
  copy_err = cudaMemcpy(device_cell_index_ptr, photon_cell_index.data(), sizeof(uint32_t) * n_photons, cudaMemcpyHostToDevice);
  Insist(!copy_err, "CUDA/HIP error copying photon cell index to device");
  #else
  // use device pointers with host side data to share code below
  uint64_t *device_photon_stream_nums_ptr = photon_stream_nums.data();
  double *device_photon_E_ptr  = photon_E.data();
  int  *device_source_type_ptr = photon_type.data();
  int  *device_photon_source_face_ptr = photon_source_face.data();
  uint32_t  *device_cell_index_ptr = photon_cell_index.data();
  #endif
  // where Aos and SoA diverge--input data above is used for both, but AoS makes its array of
  // photons to be populated by GPU kernel now, SoA can use some of these input structs directly

  if constexpr(std::is_same_v<Census_T, std::vector<Photon>>) {
    #ifdef USE_GPU
    // Allocate photon data
    Photon *device_photon_ptr;
    auto alloc_err = cudaMalloc((void **)&device_photon_ptr, sizeof(Photon) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating photons");

    #else
    std::vector<Photon> &census_photons = gpu_setup.get_census_photons();
    auto n_census_photons = census_photons.size();
    census_photons.resize(n_census_photons + n_photons);
    Photon *device_photon_ptr = census_photons.data() + n_census_photons;
    #endif

    #ifdef USE_GPU
    // Kernel settings
    int n_threads = Constants::n_threads_per_block;
    int n_blocks = (n_photons + n_threads - 1) / n_threads;
    make_source_photons<<<n_blocks, n_threads>>>(gpu_setup.get_device_cells_ptr(), dt, seed, device_photon_stream_nums_ptr, device_photon_E_ptr, device_source_type_ptr, device_photon_source_face_ptr, device_cell_index_ptr, n_photons, device_photon_ptr);

    auto kernel_err = cudaGetLastError();
    Insist(!kernel_err, "CUDA/HIP error in source kernel launch");
    auto sync_err = cudaDeviceSynchronize();
    Insist(!sync_err, "CUDA/HIP error synchronizing after source kernel");

    std::vector<Photon> &census_photons = gpu_setup.get_census_photons();
    auto n_census_photons = census_photons.size();
    census_photons.resize(n_census_photons + n_photons);

    // Copy photons back to host
    copy_err = cudaMemcpy(census_photons.data() + n_census_photons, device_photon_ptr, n_photons * sizeof(Photon), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying photons back to host");

    // Free device memory specific to AoS
    auto free_err = cudaFree(device_photon_ptr);
    Insist(!free_err, "error freeing device_photon_ptr");
    #else

    make_source_photons(mesh.get_const_cells_ptr(), dt, seed, device_photon_stream_nums_ptr, device_photon_E_ptr, device_source_type_ptr, device_photon_source_face_ptr, device_cell_index_ptr, n_photons, device_photon_ptr);

    #endif
  }
  else {
    #ifdef USE_GPU
    // Allocate photon data
    RNG *device_rng_ptr;
    auto alloc_err = cudaMalloc((void **)&device_rng_ptr, sizeof(RNG) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating rng");

    std::array<double,3> *device_pos_ptr;
    alloc_err = cudaMalloc((void **)&device_pos_ptr, sizeof(std::array<double,3>) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating position");

    std::array<double,3> *device_angle_ptr;
    alloc_err = cudaMalloc((void **)&device_angle_ptr, sizeof(std::array<double,3>) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating angle");

    double *device_life_dx_ptr;
    alloc_err = cudaMalloc((void **)&device_life_dx_ptr, sizeof(double) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating life_dx");

    uint32_t *device_group_ptr;
    alloc_err = cudaMalloc((void **)&device_group_ptr,  sizeof(uint32_t) * n_photons);
    Insist(!alloc_err, "CUDA/HIP error allocating group");
    #else
    PhotonArray &census_photons = gpu_setup.get_census_photons();
    auto n_census_photons = census_photons.size();
    census_photons.resize(n_census_photons + n_photons);
    RNG *device_rng_ptr = census_photons.rng.data() + n_census_photons;
    std::array<double,3> *device_pos_ptr = census_photons.pos.data() + n_census_photons;
    std::array<double,3> *device_angle_ptr = census_photons.angle.data() + n_census_photons;
    double *device_life_dx_ptr = census_photons.life_dx.data() + n_census_photons;
    uint32_t *device_group_ptr = census_photons.group.data() + n_census_photons;
    #endif

    #ifdef USE_GPU
    // Kernel settings
    int n_threads = Constants::n_threads_per_block;
    int n_blocks = (n_photons + n_threads - 1) / n_threads;
    set_source_photons<<<n_blocks, n_threads>>>(gpu_setup.get_device_cells_ptr(), dt, seed,
      device_photon_stream_nums_ptr, device_source_type_ptr, device_photon_source_face_ptr,
      n_photons, device_cell_index_ptr, device_rng_ptr, device_pos_ptr, device_angle_ptr,
      device_life_dx_ptr, device_group_ptr);

    auto kernel_err = cudaGetLastError();
    Insist(!kernel_err, "CUDA/HIP error in source kernel launch");
    auto sync_err = cudaDeviceSynchronize();
    Insist(!sync_err, "CUDA/HIP error synchronizing after source kernel");

    PhotonArray &census_photons = gpu_setup.get_census_photons();
    auto n_census_photons = census_photons.size();
    census_photons.resize(n_census_photons + n_photons);

    // Copy photon arrays back to host photon array object
    // CPU -> CPU copies
    std::copy( photon_E.begin(), photon_E.end(),  census_photons.E.begin() + n_census_photons);
    std::copy(photon_E.begin(), photon_E.end(), census_photons.E0.begin() + n_census_photons);
    std::copy( photon_type.begin(), photon_type.end(), census_photons.source_type.begin() + n_census_photons);
    std::fill(census_photons.descriptors.begin() + n_census_photons, census_photons.descriptors.end(), Constants::event_type::BORN_SOURCE);
    // GPU -> CPU copies
    copy_err = cudaMemcpy(census_photons.cell_ID.data() + n_census_photons, device_cell_index_ptr, n_photons * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying cell indices back to host");
    copy_err = cudaMemcpy(census_photons.group.data() + n_census_photons, device_group_ptr, n_photons * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying group back to host");
    copy_err = cudaMemcpy(census_photons.pos.data() + n_census_photons, device_pos_ptr, n_photons * sizeof(std::array<double,3>), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying pos to host");
    copy_err = cudaMemcpy(census_photons.angle.data() + n_census_photons, device_angle_ptr, n_photons * sizeof(std::array<double,3>), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying angle to host");
    copy_err = cudaMemcpy(census_photons.life_dx.data() + n_census_photons, device_life_dx_ptr, n_photons * sizeof(double), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying life_dx to host");
    copy_err = cudaMemcpy(census_photons.rng.data() + n_census_photons, device_rng_ptr, n_photons * sizeof(RNG), cudaMemcpyDeviceToHost);
    Insist(!copy_err, "CUDA/HIP error copying rng back to host");

    // Free device memory specific to SoA
    auto free_err = cudaFree(device_group_ptr );
    Insist(!free_err, "error freeing device_group_ptr");
    free_err = cudaFree(device_pos_ptr );
    Insist(!free_err, "error freeing device_pos_ptr");
    free_err = cudaFree(device_angle_ptr );
    Insist(!free_err, "error freeing device_angle_ptr");
    free_err = cudaFree(device_life_dx_ptr );
    Insist(!free_err, "error freeing device_lift_dx_ptr");
    free_err = cudaFree(device_rng_ptr );
    Insist(!free_err, "error freeing device_rng_ptr");
  #else
    set_source_photons(mesh.get_const_cells_ptr(), dt, seed,
      device_photon_stream_nums_ptr, device_source_type_ptr, device_photon_source_face_ptr,
      n_photons, device_cell_index_ptr, device_rng_ptr, device_pos_ptr, device_angle_ptr,
      device_life_dx_ptr, device_group_ptr);

    std::copy( photon_cell_index.begin() , photon_cell_index.end(), census_photons.cell_ID.begin() + n_census_photons);
    std::copy( photon_E.begin(), photon_E.end(),  census_photons.E.begin() + n_census_photons);
    std::copy(photon_E.begin(), photon_E.end(), census_photons.E0.begin() + n_census_photons);
    std::copy( photon_type.begin(), photon_type.end(), census_photons.source_type.begin() + n_census_photons);
    std::fill(census_photons.descriptors.begin() + n_census_photons, census_photons.descriptors.end(), Constants::event_type::BORN_SOURCE);
  #endif
  }

  #ifdef USE_GPU
  // Free device memory
  auto free_err = cudaFree(device_photon_stream_nums_ptr);
  Insist(!free_err, "error freeing device_photon_seed_ptr");
  free_err = cudaFree(device_photon_E_ptr);
  Insist(!free_err, "error freeing device_photon_E_ptr");
  free_err = cudaFree(device_source_type_ptr);
  Insist(!free_err, "error freeing device_source_type_ptr");
  free_err = cudaFree(device_photon_source_face_ptr);
  Insist(!free_err, "error freeing device_source_face_ptr");
  free_err = cudaFree(device_cell_index_ptr);
  Insist(!free_err, "error freeing device_cell_index_ptr");
  #endif
}


#endif // source_h_
//----------------------------------------------------------------------------//
// end of source.h
//----------------------------------------------------------------------------//
