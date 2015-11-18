/*
  Author: Alex Long
  Date: 6/29/2014
  Name: transport.h
*/

#ifndef transport_h_
#define transport_h_

#include <vector>
#include <queue>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <boost/mpi.hpp>

#include "constants.h"
#include "mesh.h"
#include "sampling_functions.h"
#include "RNG.h"
#include "request.h"

using Constants::c;
using Constants::bc_type;
using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT;
using Constants::finish_tag;

using std::vector;
using std::queue;
using std::min;
using std::sort;
using std::accumulate;

namespace mpi = boost::mpi;

bool transport_single_photon( Photon* iphtn,
                              Mesh* mesh,
                              RNG* rng,
                              double& next_dt,
                              double& exit_E,
                              double& census_E,
                              unsigned int& census_count,
                              vector<double>& rank_abs_E)
                              
{
  unsigned int elem_id, next_element;
  bc_type boundary_event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E;
  double angle[3];
  Element elem;

  unsigned int surface_cross = 0;
  double cutoff_fraction = 0.001; //note: get this from IMC_state

  elem_id=iphtn->get_element();
  elem = mesh->get_on_rank_element(elem_id);
  bool active = true;
  bool wait_flag = false;
  //transport this photon
  while(active) {
    sigma_a = elem.get_op_a();
    sigma_s = elem.get_op_s();
    f = elem.get_f();

    //get distance to event
    dist_to_scatter = -log(rng->generate_random_number())/((1.0-f)*sigma_a + sigma_s);
    dist_to_boundary = elem.get_distance_to_boundary(iphtn->get_position(),
                                                      iphtn->get_angle(),
                                                      surface_cross);
    dist_to_census = iphtn->get_distance_remaining();

    //select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    //Calculate energy absorbed by material, update photon and material energy
    absorbed_E = iphtn->get_E()*(1.0 - exp(-sigma_a*f*dist_to_event));
    iphtn->set_E(iphtn->get_E() - absorbed_E);

    rank_abs_E[elem_id] += absorbed_E;
    
    //update position
    iphtn->move(dist_to_event);

    //Apply variance/runtime reduction
    if (iphtn->below_cutoff(cutoff_fraction)) {
      rank_abs_E[elem_id] += iphtn->get_E();
      iphtn->set_dead();
      active=false;
    }
    // or apply event
    else {
      //Apply event
      //EVENT TYPE: SCATTER
      if(dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        iphtn->set_angle(angle);
      }
      //EVENT TYPE: BOUNDARY CROSS
      else if(dist_to_event == dist_to_boundary) {
        boundary_event = elem.get_bc(surface_cross);
        if(boundary_event == ELEMENT || boundary_event == PROCESSOR) {
          next_element = elem.get_next_element(surface_cross);
          iphtn->set_element(next_element);
          elem_id=next_element;
          //look for this element, if it's not there transport later
          if (mesh->mesh_available(elem_id))
            elem = mesh->get_on_rank_element(elem_id);
          else {
            mesh->request_element(elem_id);
            wait_flag = true;
            active=false;
          }
        }
        else if(boundary_event == VACUUM) {active=false; exit_E+=iphtn->get_E();}
        else iphtn->reflect(surface_cross); 
      }
      //EVENT TYPE: REACH CENSUS
      else if(dist_to_event == dist_to_census) {
        iphtn->set_census_flag(true);
        iphtn->set_distance_to_census(c*next_dt);
        active=false;
        census_count++;
        census_E+=iphtn->get_E();
      }
    } //end event loop
  } // end while alive
  return wait_flag;
}



void transport_photons(Photon*& photon_vec, 
                        unsigned int n_photon,
                        Mesh* mesh,
                        IMC_State* imc_state,
                        vector<double>& rank_abs_E,
                        Photon*& census_list,
                        int chk_freq,  
                        mpi::communicator world)
{
  unsigned int elem_id;
  unsigned int census_count = 0;
  double census_E=0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt();

  RNG *rng = imc_state->get_rng();
  Photon* iphtn;

  int n_rank =world.size();
  int rank   =world.rank();

  bool new_data = false;

  vector<vector<bool> > r_finished;
  vector<bool> b_r_finished(n_rank-1, false);
  for (unsigned int ir=0; ir<n_rank-1; ir++) {
    vector<bool> empty_bool_elem;
    r_finished.push_back(empty_bool_elem); 
  }

  mpi::request *s_finished_reqs = new mpi::request[ (n_rank-1)];
  mpi::request *r_finished_reqs = new mpi::request[ (n_rank-1)];

  // Post the receive calls for finished message
  for (unsigned int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      unsigned int r_index = ir - (ir>rank);
      r_finished_reqs[r_index] = world.irecv(ir, finish_tag, r_finished[r_index]);
    }
  }

  queue<Photon*> wait_list; //! Photons waiting for mesh data
  bool wait_flag = false;

  ////////////////////////////////////////////////////////////////////////
  // main loop over photons
  ////////////////////////////////////////////////////////////////////////
  for ( unsigned int i=0;i<n_photon; i++) {
    iphtn = &photon_vec[i];
    //get start element, only change when element crossing event
    elem_id=iphtn->get_element();

    //should always return true for photons before transport
    // but in the future we'll want to check this
    if (mesh->mesh_available(elem_id)) {
      wait_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                          census_E, census_count, rank_abs_E);
      if (wait_flag) wait_list.push(iphtn);
    } 
    else {
      mesh->request_element(elem_id);
      wait_list.push(iphtn);
    }

    // with some frequency, check for requests and try to transport the
    // waiting list
    if (i%chk_freq == 0) {
      new_data = mesh->process_mesh_requests(world);
      // if data was received, try to transport photons on waiting list
      if (new_data) {
        for (unsigned int wp =0; wp<wait_list.size(); wp++) {
          iphtn = wait_list.front();
          wait_list.pop();
          elem_id=iphtn->get_element();
          if (mesh->mesh_available(elem_id)) {
            wait_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                                census_E, census_count, rank_abs_E);
            if (wait_flag) wait_list.push(iphtn);
          }
          else wait_list.push(iphtn);
        } // end wp in wait_list
      }
    } // end if !n_photon%check_frequency

  } // end for iphtn

  ////////////////////////////////////////////////////////////////////////
  // Main transport loop finished, transport photons waiting for data
  ////////////////////////////////////////////////////////////////////////
  while (!wait_list.empty()) {
    new_data = mesh->process_mesh_requests(world);
    // if new data received, transport waiting list 
    for (unsigned int wp =0; wp<wait_list.size(); wp++) {
      iphtn = wait_list.front();
      wait_list.pop();
      elem_id=iphtn->get_element();
      if (mesh->mesh_available(elem_id)) {
        wait_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                            census_E, census_count, rank_abs_E);
        if (wait_flag) wait_list.push(iphtn);
      }
      else {
        wait_list.push(iphtn);
        mesh->request_element(elem_id);
      }
    }
  } //end while wait_list not empty

  // This rank is finished transporting, post finished
  vector<bool> s_bool(1,true);
  for (unsigned int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      unsigned int r_index = ir - (ir>rank);
      s_finished_reqs[r_index] = world.isend(ir, finish_tag, s_bool);
    }
  }

  // While waiting for other ranks to finish, check for other messages
  bool all_finished = false;
  while (!all_finished) {
    mesh->process_mesh_requests(world);
    for (unsigned int ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        unsigned int r_index = ir - (ir>rank);
        //don't check message if it's been received
        if (!b_r_finished[r_index]) {
          if (r_finished_reqs[r_index].test()) {
            all_finished = r_finished[r_index][0];
            b_r_finished[r_index] = true;
          }
          else {
            all_finished = false;
            break;
          }
        }
        else all_finished = true;
      } //end if
      else all_finished = true;
    } //end for 
  } //end while

  sort(photon_vec, photon_vec+n_photon, Photon::census_flag_compare);
  //make the census list
  census_list = new Photon[census_count];
  unsigned int num_bytes = sizeof(Photon)*census_count;
  memcpy(census_list, photon_vec, num_bytes);

  MPI::COMM_WORLD.Barrier();

  delete[] photon_vec;
  //All ranks have now finished transport
  delete[] s_finished_reqs;
  delete[] r_finished_reqs;

  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_census_size(census_count);


}



void make_photons(Mesh* mesh, IMC_State* imc_state, Photon*& phtn_vec, unsigned int& n_photon, const double& total_E)
{
  unsigned int n_element =mesh->get_number_of_objects();

  vector<double> census_E = mesh->get_census_E_ref();
  vector<double> emission_E = mesh->get_emission_E_ref();
  vector<double> source_E = mesh->get_source_E_ref();

  unsigned int total_phtns = imc_state->get_total_step_photons();
  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  //reset the number of photons before counting
  n_photon = 0;

  // count the total number of photons 
  for (unsigned int ielem = 0; ielem<n_element; ielem++) {
    if (census_E[ielem] > 0.0) { 
      unsigned int t_num_census = int(total_phtns*census_E[ielem]/total_E);
      if (t_num_census == 0) t_num_census =1;
      n_photon+= t_num_census;
    }
    if (emission_E[ielem] > 0.0) {
      unsigned int t_num_emission = int(total_phtns*emission_E[ielem]/total_E);
      if (t_num_emission == 0) t_num_emission =1;
      n_photon+=t_num_emission;
    }
    if (source_E[ielem] > 0.0) {
      unsigned int t_num_source = int(total_phtns*source_E[ielem]/total_E);
      if (t_num_source == 0) t_num_source =1;
      n_photon+=t_num_source; 
    }
  }

  // make the photon vector
  phtn_vec = new Photon[n_photon];

  Element elem;
  double pos[3]; double angle[3];
  unsigned int g_ID;
  unsigned int p_index =0;
  for (unsigned int ielem = 0; ielem<n_element; ielem++) {
    elem = mesh->get_elem(ielem);
    g_ID = elem.get_ID();

    //Census photons scope
    {
      if (census_E[ielem] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*census_E[ielem]/total_E);
        if (t_num_census == 0) t_num_census =1;
        double census_phtn_E = census_E[ielem] / t_num_census;
        for (unsigned int iphtn = 0; iphtn< t_num_census; ++iphtn) {
          elem.uniform_position_in_elem(rng, pos);
          get_uniform_angle(angle, rng);
          Photon& census_photon = phtn_vec[p_index]; 
          census_photon.set_position(pos);
          census_photon.set_angle(angle);
          census_photon.set_E0(census_phtn_E);
          census_photon.set_distance_to_census(c*delta_t);
          census_photon.set_element(g_ID);
          p_index++;
        }
      } 
    } //end census photons scope 

    //Emission photons scope
    { 
      if (emission_E[ielem] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*emission_E[ielem]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        double emission_phtn_E = emission_E[ielem] / t_num_emission;
        for (unsigned int iphtn = 0; iphtn< t_num_emission; ++iphtn) {
          elem.uniform_position_in_elem(rng, pos);
          get_uniform_angle(angle, rng);
          Photon& emission_photon = phtn_vec[p_index];
          emission_photon.set_position(pos);
          emission_photon.set_angle(angle);
          emission_photon.set_E0(emission_phtn_E);
          emission_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          emission_photon.set_element(g_ID);
          p_index++;
        }
      }
    } //end scope of emission photons

    //Source Photons scope
    {
      if (source_E[ielem] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[ielem]/total_E);
        if (t_num_source == 0) t_num_source =1;
        double source_phtn_E = source_E[ielem] / t_num_source;
        for (unsigned int iphtn = 0; iphtn<t_num_source; ++iphtn) {
          pos[0] = 0.0; pos[1] = 0.0;  pos[2] = 0.0;
          get_uniform_angle(angle, rng);
          Photon& source_photon = phtn_vec[p_index];
          source_photon.set_position(pos);
          source_photon.set_angle(angle);
          source_photon.set_E0(source_phtn_E);
          source_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          source_photon.set_element(g_ID);
          p_index++;
        }
      }
    } //end source photons scope

  } //end element loop
}


void make_stratified_photons(Mesh* mesh, IMC_State* imc_state, Photon*& phtn_vec, unsigned int& n_photon, const double& total_E)
{
  unsigned int n_element =mesh->get_number_of_objects();

  vector<double> census_E = mesh->get_census_E_ref();
  vector<double> emission_E = mesh->get_emission_E_ref();
  vector<double> source_E = mesh->get_source_E_ref();

  unsigned int total_phtns = imc_state->get_total_step_photons();
  double delta_t = imc_state->get_dt();

  RNG *rng = imc_state->get_rng(); 

  //reset the number of photons before counting
  n_photon = 0;

  // count the total number of photons 
  for (unsigned int ielem = 0; ielem<n_element; ielem++) {
    if (census_E[ielem] > 0.0) { 
      unsigned int t_num_census = int(total_phtns*census_E[ielem]/total_E);
      if (t_num_census == 0) t_num_census =1;
      n_photon+= t_num_census;
    }
    if (emission_E[ielem] > 0.0) {
      unsigned int t_num_emission = int(total_phtns*emission_E[ielem]/total_E);
      if (t_num_emission == 0) t_num_emission =1;
      n_photon+=t_num_emission;
    }
    if (source_E[ielem] > 0.0) {
      unsigned int t_num_source = int(total_phtns*source_E[ielem]/total_E);
      if (t_num_source == 0) t_num_source =1;
      n_photon+=t_num_source; 
    }
  }

  // make the photon vector
  phtn_vec = new Photon[n_photon];

  Element elem;
  double pos[3]; double angle[3];
  unsigned int g_ID;
  unsigned int p_index =0;
  for (unsigned int ielem = 0; ielem<n_element; ielem++) {
    elem = mesh->get_elem(ielem);
    g_ID = elem.get_ID();

    //Census photons scope
    {
      if (census_E[ielem] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*census_E[ielem]/total_E);
        if (t_num_census == 0) t_num_census =1;
        double census_phtn_E = census_E[ielem] / t_num_census;
        bool stratify = false;
        if (t_num_census >=8) stratify=true;
        for (unsigned int iphtn = 0; iphtn< t_num_census; ++iphtn) {
          elem.uniform_position_in_elem(rng, pos);
          if (stratify) get_stratified_angle(angle, rng, iphtn, t_num_census);
          else get_uniform_angle(angle, rng);
          Photon& census_photon = phtn_vec[p_index]; 
          census_photon.set_position(pos);
          census_photon.set_angle(angle);
          census_photon.set_E0(census_phtn_E);
          census_photon.set_distance_to_census(c*delta_t);
          census_photon.set_element(g_ID);
          p_index++;
        }
      } 
    } //end census photons scope 

    //Emission photons scope
    { 
      if (emission_E[ielem] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*emission_E[ielem]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        double emission_phtn_E = emission_E[ielem] / t_num_emission;
        bool stratify = false;
        if (t_num_emission >=8) stratify=true;
        for (unsigned int iphtn = 0; iphtn< t_num_emission; ++iphtn) {
          elem.uniform_position_in_elem(rng, pos);
          if (stratify) get_stratified_angle(angle, rng, iphtn, t_num_emission);
          else get_uniform_angle(angle, rng);
          Photon& emission_photon = phtn_vec[p_index];
          emission_photon.set_position(pos);
          emission_photon.set_angle(angle);
          emission_photon.set_E0(emission_phtn_E);
          emission_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          emission_photon.set_element(g_ID);
          p_index++;
        }
      }
    } //end scope of emission photons

    //Source Photons scope
    {
      if (source_E[ielem] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[ielem]/total_E);
        if (t_num_source == 0) t_num_source =1;
        double source_phtn_E = source_E[ielem] / t_num_source;
        for (unsigned int iphtn = 0; iphtn<t_num_source; ++iphtn) {
          pos[0] = 0.0; pos[1] = 0.0;  pos[2] = 0.0;
          get_uniform_angle(angle, rng);
          Photon& source_photon = phtn_vec[p_index];
          source_photon.set_position(pos);
          source_photon.set_angle(angle);
          source_photon.set_E0(source_phtn_E);
          source_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          source_photon.set_element(g_ID);
          p_index++;
        }
      }
    } //end source photons scope

  } //end element loop
}


double get_photon_list_energy(Photon* photon_list, const unsigned int& n_photon) {
  double list_E = 0.0;
  for (unsigned int iphtn=0; iphtn<n_photon; iphtn++)
    list_E += photon_list[iphtn].get_E();
  return list_E;
}

#endif
