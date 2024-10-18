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

#include "config.h"
#include "RNG.h"
#include "cell_tally.h"
#include "constants.h"
#include "photon.h"
#include "sampling_functions.h"


template <typename Census_T>
uint64_t post_process_photons(const double next_dt, Census_T &all_photons, Census_T &census_list, double &census_E, double &exit_E);

template <>
uint64_t post_process_photons<std::vector<Photon>>(const double next_dt, std::vector<Photon> &all_photons, std::vector<Photon> &census_list, double &census_E, double &exit_E) {
 
  uint64_t n_complete = 0;
  for (auto &phtn : all_photons) {
    auto descriptor{phtn.get_descriptor()};
    switch (descriptor) {
    case Constants::event_type::PASS:
      // handle in other function
      break;
    case Constants::event_type::KILLED:
      // note: for now killed particles go into the material so separate conservation issues here
      n_complete++;
      break;
    case Constants::event_type::EXIT:
      exit_E+=phtn.get_E();
      n_complete++;
      break;
    case Constants::event_type::CENSUS:
      phtn.set_distance_to_census(Constants::c*next_dt);
      census_list.push_back(phtn);
      census_E+=phtn.get_E();
      n_complete++;
      break;
    } //switch(descriptor)
  } // phtn : all_photons
  return n_complete;
}

template <>
uint64_t post_process_photons<PhotonArray>(const double next_dt, PhotonArray &all_photons, PhotonArray &census_list, double &census_E, double &exit_E) {
  uint64_t n_complete = 0;
  size_t census_count = 0;
  for (size_t i=0; i<all_photons.size();++i) {
    auto descriptor{all_photons.descriptors[i]};
    switch (descriptor) {
    case Constants::event_type::PASS:
      // handle in other function
      break;
    case Constants::event_type::KILLED:
      // note: for now killed particles go into the material so separate conservation issues here
      n_complete++;
      break;
    case Constants::event_type::EXIT:
      exit_E+=all_photons.E[i];
      n_complete++;
      break;
    case Constants::event_type::CENSUS:
      all_photons.life_dx[i] = Constants::c*next_dt;
      census_E+=all_photons.E[i];
      census_count++;
      n_complete++;
      break;
    } //switch(descriptor)
  }

  // append census photons on to census list
  auto old_census_size = census_list.size();
  auto new_census_size = old_census_size + census_count; 
  census_list.resize(new_census_size);
  auto census_index = old_census_size; 
  for ( size_t i = 0; i < all_photons.cell_ID.size(); ++i) {
    if (static_cast<Constants::event_type>(all_photons.descriptors[i]) == Constants::event_type::CENSUS){
      census_list.cell_ID[census_index] = all_photons.cell_ID[i];
      census_list.group[census_index] = all_photons.group[i];
      census_list.source_type[census_index] = all_photons.source_type[i];
      census_list.descriptors[census_index] = all_photons.descriptors[i];
      census_list.pos[census_index] = all_photons.pos[i];
      census_list.angle[census_index] = all_photons.angle[i];
      census_list.E[census_index] = all_photons.E[i];
      census_list.E0[census_index] = all_photons.E0[i];
      census_list.life_dx[census_index] = all_photons.life_dx[i];
      census_list.rng[census_index] = all_photons.rng[i];
      census_index++;
    }
  }
  return n_complete;
}

#endif // def transport_photon_h_
//----------------------------------------------------------------------------//
// end of transport_photon.h
//----------------------------------------------------------------------------//
