//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comb_photons.h
 * \author Alex Long
 * \date   June 21 2018
 * \brief  Comb census routine, yes, I'm finally putting it in
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef comb_photons_h_
#define comb_photons_h_

#include <unordered_map>

void comb_photons(std::vector<Photon> &census_photons,
                  int64_t max_census_photons, RNG *rng) {
  std::unordered_map<uint32_t, uint32_t> cell_census_count;
  std::unordered_map<uint32_t, double> cell_census_E;
  std::unordered_map<uint32_t, double> cell_corrected_E;

  // only comb if the photon population is too large
  uint64_t n_global_census = census_photons.size();
  MPI_Allreduce(MPI_IN_PLACE, &n_global_census, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  if (n_global_census > max_census_photons) {
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
}

#endif // comb_photons_h_
