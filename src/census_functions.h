//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   census_functions.h
 * \author Alex Long
 * \date   June 21 2018
 * \brief  Comb census routine, yes, I'm finally putting it in
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef census_functions_h_
#define census_functions_h_

#include "photon_array.h"
#include "photon.h"
#include <unordered_map>


template <typename Census_T>
void join_photon_arrays(Census_T &original, Census_T &to_append) {

  if constexpr(std::is_same_v<Census_T, std::vector<Photon>>) {
    original.insert(original.end(), to_append.begin(), to_append.end());
  }
  else {
    original.insert(to_append);
  }
}

double get_photon_list_E(const PhotonArray &photon_array) {
  double total_E = 0.0;
  for (size_t i = 0; i < photon_array.E.size(); ++i)
  {
    total_E += photon_array.E[i];
  }
  return total_E;
}

double get_photon_list_E(const std::vector<Photon> &photons) {
  double total_E = 0.0;
  for (auto const &iphtn : photons) {
    total_E += iphtn.get_E();
  }
  return total_E;
}

void comb_photons(std::vector<Photon> &census_photons,
                  int64_t max_census_photons, RNG *rng) {
  std::unordered_map<uint32_t, uint32_t> cell_census_count;
  std::unordered_map<uint32_t, double> cell_census_E;
  std::unordered_map<uint32_t, double> cell_corrected_E;

  std::vector<Photon> post_comb_photons;

  double census_total_E = get_photon_list_E(census_photons);
  double global_census_E = census_total_E;
  MPI_Allreduce(MPI_IN_PLACE, &global_census_E, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  double comb_photon_E = global_census_E / double(max_census_photons);

  for (auto &p : census_photons) {
    cell_census_count[p.get_cell()]++;
    cell_census_E[p.get_cell()] += p.get_E();
  }

  for (auto ip = census_photons.begin(); ip != census_photons.end(); ++ip) {
    uint32_t cell = ip->get_cell();
    double p_kill = 1.0 - ip->get_E() / comb_photon_E;
    double rand_check = rng->generate_random_number();
    // save the photon probabilistically or if there is only one photon
    // remaining in the cell
    if (rand_check > p_kill || cell_census_count[cell] == 1) {
      ip->set_E(comb_photon_E);
      cell_corrected_E[cell] += comb_photon_E;
      post_comb_photons.push_back(*ip);
    } else {
      cell_census_count[cell]--;
    }
  }

  // correct energy for conservation
  double new_E;
  uint32_t cell;
  for (auto &p : post_comb_photons) {
    cell = p.get_cell();
    new_E = p.get_E() + (cell_census_E[cell] - cell_corrected_E[cell]) /
                            double(cell_census_count[cell]);
    p.set_E(new_E);
  }

  census_photons = post_comb_photons;
}

#endif // census_functions_h_
