/*
  Author: Alex Long
  Date: 12/1/2015
  Name: transport.h
*/

#ifndef transport_h_
#define transport_h_

#include <vector>

#include "photon.h"

double get_photon_list_E(std::vector<Photon> photons) {
  double total_E = 0.0;
  for (std::vector<Photon>::iterator iphtn=photons.begin(); 
      iphtn<photons.end(); 
      iphtn++)
    total_E += iphtn->get_E();
  return total_E;
}


std::vector<Photon> make_initial_census_photons(Mesh* mesh, 
                                          IMC_State* imc_state, 
                                          const double& total_E, 
                                          const uint32_t& n_user_photon) {
  using std::vector;
  using Constants::c;

  vector<double> census_E = mesh->get_census_E_ref();

  uint32_t n_cell =mesh->get_number_of_objects();

  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  vector<Photon> census_photons;

  // make the photons
  Cell cell;
  double pos[3]; double angle[3];
  uint32_t g_ID;
  uint32_t p_index =0;
  for (uint32_t icell = 0; icell<n_cell; icell++) {
    cell = mesh->get_cell(icell);
    g_ID = cell.get_ID();

    //Census photons scope
    if (census_E[icell] > 0.0) { 
      uint32_t t_num_census = int(n_user_photon*census_E[icell]/total_E);
      if (t_num_census == 0) t_num_census =1;
      double census_phtn_E = census_E[icell] / t_num_census;
      for (uint32_t iphtn = 0; iphtn< t_num_census; ++iphtn) {
        cell.uniform_position_in_cell(rng, pos);
        get_uniform_angle(angle, rng);
        Photon census_photon;
        census_photon.set_position(pos);
        census_photon.set_angle(angle);
        census_photon.set_E0(census_phtn_E);
        census_photon.set_distance_to_census(c*delta_t);
        census_photon.set_cell(g_ID);
        census_photons.push_back(census_photon);
        p_index++;
      }
    }
  } //end cell loop
  return census_photons;
}

#endif // def transport_h_
