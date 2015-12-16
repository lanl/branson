/*
  Author: Alex Long
  Date: 12/1/2015
  Name: transport.h
*/

#ifndef transport_h_
#define transport_h_

#include <vector>

#include "photon.h"

double get_photon_list_energy(Photon* photon_list, const unsigned int& n_photon) {
  double list_E = 0.0;
  for (unsigned int iphtn=0; iphtn<n_photon; iphtn++)
    list_E += photon_list[iphtn].get_E();
  return list_E;
}

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
                                          const unsigned int& n_user_photon) {
  using std::vector;
  using Constants::c;

  vector<double> census_E = mesh->get_census_E_ref();

  unsigned int n_cell =mesh->get_number_of_objects();

  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  vector<Photon> census_photons;

  // make the photons
  Cell cell;
  double pos[3]; double angle[3];
  unsigned int g_ID;
  unsigned int p_index =0;
  for (unsigned int icell = 0; icell<n_cell; icell++) {
    cell = mesh->get_cell(icell);
    g_ID = cell.get_ID();

    //Census photons scope
    if (census_E[icell] > 0.0) { 
      unsigned int t_num_census = int(n_user_photon*census_E[icell]/total_E);
      if (t_num_census == 0) t_num_census =1;
      double census_phtn_E = census_E[icell] / t_num_census;
      for (unsigned int iphtn = 0; iphtn< t_num_census; ++iphtn) {
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



void make_photons(Mesh* mesh, IMC_State* imc_state, Photon*& phtn_vec, unsigned int& n_photon, const double& total_E)
{
  using std::vector;
  using Constants::c;
  unsigned int n_cell =mesh->get_number_of_objects();

  vector<double> census_E = mesh->get_census_E_ref();
  vector<double> emission_E = mesh->get_emission_E_ref();
  vector<double> source_E = mesh->get_source_E_ref();

  unsigned int total_phtns = imc_state->get_total_step_photons();
  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  //reset the number of photons before counting
  n_photon = 0;

  // count the total number of photons 
  for (unsigned int icell = 0; icell<n_cell; icell++) {
    if (census_E[icell] > 0.0) { 
      unsigned int t_num_census = int(total_phtns*census_E[icell]/total_E);
      if (t_num_census == 0) t_num_census =1;
      n_photon+= t_num_census;
    }
    if (emission_E[icell] > 0.0) {
      unsigned int t_num_emission = int(total_phtns*emission_E[icell]/total_E);
      if (t_num_emission == 0) t_num_emission =1;
      n_photon+=t_num_emission;
    }
    if (source_E[icell] > 0.0) {
      unsigned int t_num_source = int(total_phtns*source_E[icell]/total_E);
      if (t_num_source == 0) t_num_source =1;
      n_photon+=t_num_source; 
    }
  }

  // make the photon vector
  phtn_vec = new Photon[n_photon];

  Cell cell;
  double pos[3]; double angle[3];
  unsigned int g_ID;
  unsigned int p_index =0;
  for (unsigned int icell = 0; icell<n_cell; icell++) {
    cell = mesh->get_cell(icell);
    g_ID = cell.get_ID();

    //Census photons scope
    {
      if (census_E[icell] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*census_E[icell]/total_E);
        if (t_num_census == 0) t_num_census =1;
        double census_phtn_E = census_E[icell] / t_num_census;
        for (unsigned int iphtn = 0; iphtn< t_num_census; ++iphtn) {
          cell.uniform_position_in_cell(rng, pos);
          get_uniform_angle(angle, rng);
          Photon& census_photon = phtn_vec[p_index]; 
          census_photon.set_position(pos);
          census_photon.set_angle(angle);
          census_photon.set_E0(census_phtn_E);
          census_photon.set_distance_to_census(c*delta_t);
          census_photon.set_cell(g_ID);
          p_index++;
        }
      } 
    } //end census photons scope 

    //Emission photons scope
    { 
      if (emission_E[icell] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*emission_E[icell]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        double emission_phtn_E = emission_E[icell] / t_num_emission;
        for (unsigned int iphtn = 0; iphtn< t_num_emission; ++iphtn) {
          cell.uniform_position_in_cell(rng, pos);
          get_uniform_angle(angle, rng);
          Photon& emission_photon = phtn_vec[p_index];
          emission_photon.set_position(pos);
          emission_photon.set_angle(angle);
          emission_photon.set_E0(emission_phtn_E);
          emission_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          emission_photon.set_cell(g_ID);
          p_index++;
        }
      }
    } //end scope of emission photons

    //Source Photons scope
    {
      if (source_E[icell] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[icell]/total_E);
        if (t_num_source == 0) t_num_source =1;
        double source_phtn_E = source_E[icell] / t_num_source;
        for (unsigned int iphtn = 0; iphtn<t_num_source; ++iphtn) {
          pos[0] = 0.0; pos[1] = 0.0;  pos[2] = 0.0;
          get_uniform_angle(angle, rng);
          Photon& source_photon = phtn_vec[p_index];
          source_photon.set_position(pos);
          source_photon.set_angle(angle);
          source_photon.set_E0(source_phtn_E);
          source_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          source_photon.set_cell(g_ID);
          p_index++;
        }
      }
    } //end source photons scope

  } //end cell loop
}


void make_stratified_photons(Mesh* mesh, IMC_State* imc_state, Photon*& phtn_vec, unsigned int& n_photon, const double& total_E)
{
  using Constants::c;
  using std::vector;
  unsigned int n_cell =mesh->get_number_of_objects();

  vector<double> census_E = mesh->get_census_E_ref();
  vector<double> emission_E = mesh->get_emission_E_ref();
  vector<double> source_E = mesh->get_source_E_ref();

  unsigned int total_phtns = imc_state->get_total_step_photons();
  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  //reset the number of photons before counting
  n_photon = 0;

  // count the total number of photons 
  for (unsigned int icell = 0; icell<n_cell; icell++) {
    if (census_E[icell] > 0.0) { 
      unsigned int t_num_census = int(total_phtns*census_E[icell]/total_E);
      if (t_num_census == 0) t_num_census =1;
      n_photon+= t_num_census;
    }
    if (emission_E[icell] > 0.0) {
      unsigned int t_num_emission = int(total_phtns*emission_E[icell]/total_E);
      if (t_num_emission == 0) t_num_emission =1;
      n_photon+=t_num_emission;
    }
    if (source_E[icell] > 0.0) {
      unsigned int t_num_source = int(total_phtns*source_E[icell]/total_E);
      if (t_num_source == 0) t_num_source =1;
      n_photon+=t_num_source; 
    }
  }

  // make the photon vector
  phtn_vec = new Photon[n_photon];

  Cell cell;
  double pos[3]; double angle[3];
  unsigned int g_ID;
  unsigned int p_index =0;
  for (unsigned int icell = 0; icell<n_cell; icell++) {
    cell = mesh->get_cell(icell);
    g_ID = cell.get_ID();

    //Census photons scope
    {
      if (census_E[icell] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*census_E[icell]/total_E);
        if (t_num_census == 0) t_num_census =1;
        double census_phtn_E = census_E[icell] / t_num_census;
        bool stratify = false;
        if (t_num_census >=8) stratify=true;
        for (unsigned int iphtn = 0; iphtn< t_num_census; ++iphtn) {
          cell.uniform_position_in_cell(rng, pos);
          if (stratify) get_stratified_angle(angle, rng, iphtn, t_num_census);
          else get_uniform_angle(angle, rng);
          Photon& census_photon = phtn_vec[p_index]; 
          census_photon.set_position(pos);
          census_photon.set_angle(angle);
          census_photon.set_E0(census_phtn_E);
          census_photon.set_distance_to_census(c*delta_t);
          census_photon.set_cell(g_ID);
          p_index++;
        }
      } 
    } //end census photons scope 

    //Emission photons scope
    { 
      if (emission_E[icell] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*emission_E[icell]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        double emission_phtn_E = emission_E[icell] / t_num_emission;
        bool stratify = false;
        if (t_num_emission >=8) stratify=true;
        for (unsigned int iphtn = 0; iphtn< t_num_emission; ++iphtn) {
          cell.uniform_position_in_cell(rng, pos);
          if (stratify) get_stratified_angle(angle, rng, iphtn, t_num_emission);
          else get_uniform_angle(angle, rng);
          Photon& emission_photon = phtn_vec[p_index];
          emission_photon.set_position(pos);
          emission_photon.set_angle(angle);
          emission_photon.set_E0(emission_phtn_E);
          emission_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          emission_photon.set_cell(g_ID);
          p_index++;
        }
      }
    } //end scope of emission photons

    //Source Photons scope
    {
      if (source_E[icell] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[icell]/total_E);
        if (t_num_source == 0) t_num_source =1;
        double source_phtn_E = source_E[icell] / t_num_source;
        for (unsigned int iphtn = 0; iphtn<t_num_source; ++iphtn) {
          pos[0] = 0.0; pos[1] = 0.0;  pos[2] = 0.0;
          get_uniform_angle(angle, rng);
          Photon& source_photon = phtn_vec[p_index];
          source_photon.set_position(pos);
          source_photon.set_angle(angle);
          source_photon.set_E0(source_phtn_E);
          source_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          source_photon.set_cell(g_ID);
          p_index++;
        }
      }
    } //end source photons scope

  } //end cell loop
}


#endif // def transport_h_
