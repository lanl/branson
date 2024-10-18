//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   photon.h
 * \author Alex Long
 * \date   July 18 2014
 * \brief  Holds values and functions needed for transporting photon
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef photon_array_h_
#define photon_array_h_

#include <cmath>
#include <iostream>
#include <vector>
#include <array>

#include "constants.h"
#include "config.h"
#include "RNG.h"

//Strucutre of arrrays to store photon attributes
class PhotonArray {
public:
  std::vector<uint32_t> cell_ID;
  std::vector<uint32_t> group;
  std::vector<uint32_t> source_type;
  std::vector<unsigned char> descriptors;
  std::vector<std::array<double, 3>> pos;
  std::vector<std::array<double, 3>> angle;
  std::vector<double> E;
  std::vector<double> E0;
  std::vector<double> life_dx;
  std::vector<RNG> rng;

  void push_back(const Photon& photon) {
    cell_ID.push_back(photon.get_cell());
    group.push_back(photon.get_group());
    source_type.push_back(photon.get_source_type());
    descriptors.push_back(photon.get_descriptor());
    pos.push_back(photon.get_position());
    angle.push_back(photon.get_angle());
    E.push_back(photon.get_E());
    E0.push_back(photon.get_E0());
    life_dx.push_back(photon.get_distance_remaining());
    rng.push_back(photon.get_rng());
  }

  void resize(size_t size)
  {
    cell_ID.resize(size);
    group.resize(size);
    source_type.resize(size);
    descriptors.resize(size);
    pos.resize(size);
    angle.resize(size);
    E.resize(size);
    E0.resize(size);
    life_dx.resize(size);
    rng.resize(size);
  }

  void add_photon(const PhotonArray &source, size_t index) {
    size_t new_index = cell_ID.size();
    cell_ID.push_back(source.cell_ID[index]);
    group.push_back(source.group[index]);
    source_type.push_back(source.source_type[index]);
    descriptors.push_back(source.descriptors[index]);
    pos.push_back(source.pos[index]);
    angle.push_back(source.angle[index]);
    E.push_back(source.E[index]);
    E0.push_back(source.E0[index]);
    life_dx.push_back(source.life_dx[index]);
    rng.push_back(source.rng[index]);
  }

  bool empty() const {return cell_ID.empty();}

  size_t size() const {return cell_ID.size();}

  Photon operator [](size_t i) {
    Photon return_photon;
    return_photon.set_source_type(source_type[i]);
    return_photon.set_position(pos[i]);
    return_photon.set_angle(angle[i]);
    return_photon.set_E0(E0[i]);
    return_photon.set_distance_to_census(life_dx[i]);
    return_photon.set_cell(cell_ID[i]);
    return_photon.set_group(group[i]);
    return_photon.set_rng(rng[i]);
    return return_photon;
  }

  void insert(const PhotonArray &photons_to_add) {
     size_t original_size = size();
     size_t new_size = original_size + photons_to_add.size();
     resize(new_size);

    for (size_t i = 0; i < photons_to_add.cell_ID.size(); ++i) {
      cell_ID[original_size+i] = photons_to_add.cell_ID[i];
      group[original_size+i] = photons_to_add.group[i];
      source_type[original_size+i] = photons_to_add.source_type[i];
      descriptors[original_size+i] = photons_to_add.descriptors[i];
      pos[original_size+i] = photons_to_add.pos[i];
      angle[original_size+i] = photons_to_add.angle[i];
      E[original_size+i] = photons_to_add.E[i];
      E0[original_size+i] = photons_to_add.E0[i];
      life_dx[original_size+i] = photons_to_add.life_dx[i];
      rng[original_size+i] = photons_to_add.rng[i];
    }
  }
 
  // this copies out photon data 
  PhotonArray get_sub_batch(const size_t batch_start, const size_t batch_end) {
    size_t batch_size = batch_end - batch_start;
    PhotonArray batch_photons;
    batch_photons.resize(batch_size);
    for (size_t i = 0; i < batch_size; ++i) {
      size_t idx = batch_start + i;
      batch_photons.cell_ID[i] = cell_ID[idx];
      batch_photons.group[i] = group[idx];
      batch_photons.source_type[i] = source_type[idx];
      batch_photons.descriptors[i] = descriptors[idx];
      batch_photons.pos[i] = pos[idx];
      batch_photons.angle[i] = angle[idx];
      batch_photons.E[i] = E[idx];
      batch_photons.E0[i] = E0[idx];
      batch_photons.life_dx[i] = life_dx[idx];
      batch_photons.rng[i] = rng[idx];
    }
    return batch_photons;
  }
};

//==============================================================================
/*!
 * \class Photon
 * \brief Contains position, direction, cell ID and energy for transport.
 *
 * Holds all of the internal state of a photon and provides functions for
 * sorting photons based on census and global cell ID.
 */
//==============================================================================
/*
class Photon {
public:
  //! Constructor
  Photon(PhotonArray &ray, size_t index) : ray(ray), index(index) {}

  //! Destructor
  ~Photon() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Check to see if photon energy weight is below cutoff fraction
  GPU_HOST_DEVICE
  bool below_cutoff(const double cutoff_fraction) const {
    return (ray.E[index] / ray.E0[index] < cutoff_fraction);
  }

  GPU_HOST_DEVICE
  inline double get_fraction() const {return ray.E[index] / ray.E0[index];}

  //! Return global cell ID
  GPU_HOST_DEVICE
  inline uint32_t get_cell(void) const { return ray.cell_ID[index]; }

  //! Return photon group
  GPU_HOST_DEVICE
  inline uint32_t get_group(void) const { return ray.group[index]; }

  //! Return a constant pointer to the start of the particle position array
  GPU_HOST_DEVICE
  inline std::array<double,3> get_position(void) const { return ray.pos[index]; }

  //! Return a constant pointer to the start of the particle direction array
  GPU_HOST_DEVICE
  inline std::array<double,3> get_angle(void) const { return ray.angle[index]; }

  //! Get the particle's energy-weight
  GPU_HOST_DEVICE
  inline double get_E(void) const { return ray.E[index]; }

  //! Get the partice's initial energy-weight
  GPU_HOST_DEVICE
  inline double get_E0(void) const { return ray.E0[index]; }

  //! Get the distance to census (cm)
  GPU_HOST_DEVICE
  inline double get_distance_remaining(void) const { return ray.life_dx[index]; }

  //! Print particle information
  void print_info(const uint32_t &rank) const {
    using std::cout;
    using std::endl;
    cout << "----Photon Info----\n";
    cout << rank << " " << ray.pos[index][0] << " " << ray.pos[index][1] << " " << ray.pos[index][2]
         << endl;
    cout << "angle: " << ray.angle[index][0] << " " << ray.angle[index][1] << " " << ray.angle[index][2]
         << endl;
    cout << "Energy: " << ray.E[index] << " , Initial energy: " << ray.E0[index] << endl;
    cout << "Cell ID: " << ray.cell_ID[index] << endl;
  }

  //! Override great than operator to sort
  bool operator<(const Photon &compare) const {
    return ray.cell_ID[index] < compare.get_cell();
  }

  GPU_HOST_DEVICE
  Constants::event_type get_descriptor() {return static_cast<Constants::event_type>(ray.descriptors[index][0]);}

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Update particle position by moving it a distance
  GPU_HOST_DEVICE
  inline void move(const double distance) {
    for (int i = 0; i < 3; ++i)
    {
      ray.pos[index][i] += ray.angle[index][i] * distance;
    }
    ray.life_dx[index] -= distance;
  }

  inline uint32_t get_source_type() const {return ray.source_type[index];}
  inline void set_source_type(uint32_t source_type_in) {ray.source_type[index] = source_type_in;}

  //! Set the global cell ID
  GPU_HOST_DEVICE
  inline void set_cell(const uint32_t new_cell) { ray.cell_ID[index] = new_cell; }

  //! Set the group of the photon
  GPU_HOST_DEVICE
  inline void set_group(const uint32_t new_group) { ray.group[index] = new_group; }

  //! Set the initial energy-weight
  inline void set_E0(const double E) {
    ray.E0[index] = E;
    ray.E[index] = E;
  }

  //! Set the current energy-weight
  GPU_HOST_DEVICE
  inline void set_E(const double E) { ray.E[index] = E; }

  //! Set the distance to census (cm)
  GPU_HOST_DEVICE
  inline void set_distance_to_census(const double dist_remain) { ray.life_dx[index] = dist_remain; }

  //! Set the angle of the photon
  GPU_HOST_DEVICE
  inline void set_angle(const std::array<double,3> &new_angle) { ray.angle[index] = new_angle;}

  //! Set the spatial position of the photon
  GPU_HOST_DEVICE
  inline void set_position(const std::array<double, 3> &new_pos) { ray.pos[index] = new_pos;}

  //! Reflect a photon about a plane aligned with the X, Y, or Z axes
  GPU_HOST_DEVICE
  inline void reflect(const uint32_t surface_cross) {
    // reflect the photon over the surface it was crossing
    int reflect_angle = surface_cross/2; // X -> 0, Y->1, Z->2
    ray.angle[index][reflect_angle] = -ray.angle[index][reflect_angle];
  }

  GPU_HOST_DEVICE
  void set_descriptor(const Constants::event_type descriptor) { ray.descriptors[index][0] = static_cast<unsigned char>(descriptor);}

  GPU_HOST_DEVICE
  RNG &get_rng() {return ray.rng[index];}

  void set_rng(const RNG &rng) { ray.rng[index] = rng;}

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  PhotonArray ray;
  size_t index;

};
*/
#endif // photon_array_h_
