/*
  Author: Alex Long
  Date: 6/29/2014
  Name: transport.h
*/

#ifndef transport_h_
#define transport_h_

#include <vector>
#include <omp.h>
#include <numeric>
#include <algorithm>

#include "constants.h"
#include "mesh.h"
#include "sampling_functions.h"
#include "RNG.h"

//Requires C++11
#include "parallel_stable_sort.h"

using Constants::c;
using Constants::bc_type;
using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT;

using std::vector;
using std::min;
using std::sort;
using std::accumulate;


void make_census_list(Photon *photon_list, 
                      unsigned int num_photons,
                      double dt, 
                      vector<Photon>& census_list) {
  for ( unsigned int i=0; i < num_photons; i++) {
    Photon iphtn = photon_list[i];
    if (iphtn.get_census_flag()) {
      iphtn.set_census_flag(false);
      iphtn.set_distance_to_census(c*dt);
      census_list.push_back(iphtn);
    }
  }
}



void transport_photons(Photon *photon_list, 
                        unsigned int num_photons,
                        Mesh* mesh,
                        RNG *rng,
                        vector<double>& thread_abs_E,
                        double& exit_E,
                        unsigned int& thread_census_count,
                        vector<unsigned int>& element_census,
                        double& census_E )
{
  bool active;
  unsigned int elem_id, next_element;
  bc_type boundary_event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma, f, absorbed_E;
  double angle[3];
  Element* elem;

  unsigned int surface_cross = 0;

  double cutoff_fraction = 0.0001; //note: get this from IMC_state

  //transport photons
  for ( unsigned int i=0; i<num_photons; i++) {
    active=true;
    Photon *iphtn = &photon_list[i];
    while(active) {
      
      //set properties
      elem_id=iphtn->get_element();
      elem = mesh->get_element(elem_id);
      sigma = elem->get_op();   
      f = elem->get_f();

      //get distance to event
      dist_to_scatter = -log(rng->generate_random_number())/((1.0-f)*sigma);
      dist_to_boundary = elem->get_distance_to_boundary(iphtn->get_position(),
                                                        iphtn->get_angle(),
                                                        surface_cross);
      dist_to_census = iphtn->get_distance_remaining();

      //select minimum distance event
      dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

      //Calculate energy absorbed by material, update photon and material energy
      absorbed_E = iphtn->get_E()*(1.0 - exp(-sigma*f*dist_to_event));
      iphtn->set_E(iphtn->get_E() - absorbed_E);

      thread_abs_E[elem_id] += absorbed_E;

      //update position
      iphtn->move(dist_to_event);

      //Apply variance/runtime reduction
      if (iphtn->below_cutoff(cutoff_fraction)) {
        thread_abs_E[elem_id] += iphtn->get_E();
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
          boundary_event = elem->get_bc(surface_cross);
          if(boundary_event == ELEMENT) {
            next_element = elem->get_next_element(surface_cross);
            iphtn->set_element(next_element);
          }
          else if(boundary_event == VACUUM) {active=false; exit_E+=iphtn->get_E();}
          else iphtn->reflect(surface_cross); 
        }
        //EVENT TYPE: REACH CENSUS
        else if(dist_to_event == dist_to_census) {
          iphtn->set_census_flag(true);
          active=false;
          thread_census_count++;
          census_E+=iphtn->get_E();
          element_census[elem_id]++;
        }
      } //end event loop
    } // end while alive
  } // end for iphtn
}



void  make_photons(Mesh* mesh, 
                   IMC_State* imc_state,
                   vector<double>& reduced_abs_E,
                   vector<unsigned int>& element_census_count,
                   Photon*& census_list)
{
  unsigned int n_element =mesh->get_num_elems();
  unsigned int n_threads = omp_get_max_threads();

  unsigned int *thread_census_count_array = new unsigned int[n_threads];
  unsigned int *prefix_census_count_array = new unsigned int[n_threads];

  //count the number of census particles in each element
  vector<unsigned int> prefix_element_census(n_element,0);
  std::partial_sum(element_census_count.begin(), 
                   element_census_count.end(), 
                   prefix_element_census.begin());


  unsigned int total_census_count = 0;
  double tot_census_E = 0.0;
  #pragma omp parallel reduction(+:tot_census_E)
  {
    unsigned int thread_ID = omp_get_thread_num();

    unsigned int thread_total_census = 0;

    double total_E = mesh->get_total_photon_E();
    vector<double> census_E = mesh->get_census_E_vector();
    vector<double> emission_E = mesh->get_emission_E_vector();
    vector<double> source_E = mesh->get_source_E_vector();

    vector<unsigned int> thread_element_census(n_element, 0);

    unsigned int total_phtns = imc_state->get_total_step_photons();
    double delta_t = imc_state->get_dt();
    RNG *rng = new RNG();
    rng->set_seed(1241512*omp_get_thread_num()+1);
    Element* elem = 0;
    vector<double> thread_abs_E(n_element, 0.0);
    vector<Photon> thread_census_list;
    double exit_E = 0.0;
    double pos[3];
    double angle[3];
    #pragma omp for schedule(guided) 
    for (unsigned int elem_ID = 0; elem_ID<n_element; elem_ID++) {
      elem = mesh->get_element(elem_ID);

      unsigned int total_count = 0;
      unsigned int num_census=0;
      unsigned int num_emission=0;
      unsigned int num_source=0;

      if (census_E[elem_ID] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*census_E[elem_ID]/total_E);
        if (t_num_census == 0) t_num_census =1;
        total_count+= t_num_census;
        num_census = t_num_census;
      }
      // add in the census photons from last timestep
      total_count+=element_census_count[elem_ID];
      num_census+=element_census_count[elem_ID];

      if (emission_E[elem_ID] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*emission_E[elem_ID]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        total_count+=t_num_emission;
        num_emission = t_num_emission;
      }

      if (source_E[elem_ID] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[elem_ID]/total_E);
        if (t_num_source == 0) t_num_source =1;
        total_count+=t_num_source; 
        num_source = t_num_source;
      }

      Photon *phtn_vec = new Photon[total_count];
      //Census photons scope
      {
        if (imc_state->get_step()==1) {
          unsigned int elem_census_phtns = num_census;
          double census_phtn_E = census_E[elem_ID] / elem_census_phtns;
          for (unsigned int iphtn = 0; iphtn< num_census; ++iphtn) {
            elem->uniform_position_in_elem(rng, pos);
            get_uniform_angle(angle, rng);
            Photon& census_photon = phtn_vec[iphtn];
            census_photon.set_position(pos);
            census_photon.set_angle(angle);
            census_photon.set_E0(census_phtn_E);
            census_photon.set_distance_to_census(c*delta_t);
            census_photon.set_element(elem_ID);
          }
        }
        else {
          //copy census photons for this element
          unsigned int start_index = prefix_element_census[elem_ID]-num_census;
          unsigned int num_bytes = sizeof(Photon)*num_census;
          memcpy( &phtn_vec[0], &census_list[start_index], num_bytes);
          //set the timestep for census photons
          for (unsigned int iphtn = 0; iphtn< num_census; ++iphtn) {
            Photon& census_photon = phtn_vec[iphtn];
            census_photon.set_distance_to_census(c*delta_t);
          }
          
        }
      } //end census photons scope

      //Emission photons scope
      { 
        unsigned int elem_emission_phtns = num_emission;
        double emission_phtn_E = emission_E[elem_ID] / elem_emission_phtns;
        for (unsigned int iphtn = num_census; iphtn< num_census+num_emission; ++iphtn) {
          elem->uniform_position_in_elem(rng, pos);
          get_uniform_angle(angle, rng);
          Photon& emission_photon = phtn_vec[iphtn];
          emission_photon.set_position(pos);
          emission_photon.set_angle(angle);
          emission_photon.set_E0(emission_phtn_E);
          emission_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          emission_photon.set_element(elem_ID);
        }
      } //end scope of emission photons

      //Source Photons scope
      {
        unsigned int elem_source_phtns = num_source;
        double source_phtn_E = source_E[elem_ID] / elem_source_phtns;
        for (unsigned int iphtn = num_census+num_emission; iphtn<total_count; ++iphtn) {
          pos[0] = 0.0; pos[1] = 0.0;  pos[2] = 0.0;
          get_uniform_angle(angle, rng);
          Photon& source_photon = phtn_vec[iphtn];
          source_photon.set_position(pos);
          source_photon.set_angle(angle);
          source_photon.set_E0(source_phtn_E);
          source_photon.set_distance_to_census(rng->generate_random_number()*c*delta_t);
          source_photon.set_element(elem_ID);
        }
      } //end source photons scope

      //transport photons
      transport_photons(phtn_vec, total_count, mesh, rng, thread_abs_E, 
                        exit_E, thread_total_census, thread_element_census, 
                        tot_census_E);
      make_census_list(phtn_vec, total_count, delta_t, thread_census_list);
      delete[] phtn_vec;
    } //end element loop
  
    //set total census photons for each thread 
    thread_census_count_array[thread_ID] = thread_total_census;
    #pragma omp barrier

    //do prefix sums on census count array
    //reset element_census list before reduction
    #pragma omp single
    {
      std::partial_sum(thread_census_count_array, 
                        thread_census_count_array+n_threads, 
                        prefix_census_count_array);
      total_census_count = prefix_census_count_array[n_threads-1];
      census_list = new Photon[total_census_count];
      for (unsigned int i=0; i<n_element; i++) element_census_count[i]=0; 
    }
    #pragma omp barrier

    // each thread writes to the census_list
    unsigned int start_census_index=prefix_census_count_array[thread_ID]-thread_total_census;
    unsigned int num_bytes = sizeof(Photon)*thread_total_census;
    memcpy( &census_list[start_census_index], &thread_census_list[0], num_bytes);

    //do reduction on absorption vector and census particles in element
    #pragma omp critical
    {
      for (unsigned int i=0; i<n_element; i++) {
        reduced_abs_E[i]+=thread_abs_E[i];
        element_census_count[i]+=thread_element_census[i];
      }
    } //end OMP critical

  } //end OMP parallel

  //call the OpenMP sort routine
  pss::parallel_stable_sort( census_list, census_list+total_census_count, census_list[0]);

  imc_state->set_census_photons(total_census_count);
  imc_state->set_post_census_E(tot_census_E);
  //imc_state->set_absorbed_E(accumulate(reduced_abs_E.begin(), reduced_abs_E.end(), 0.0));

  delete[] thread_census_count_array;
  delete[] prefix_census_count_array;
}


double get_photon_list_energy(vector<Photon> photon_list) {
  double list_E = 0.0;
  for (vector<Photon>::iterator iphtn=photon_list.begin(); iphtn!=photon_list.end(); iphtn++)
    list_E += iphtn->get_E();
  return list_E;
}


#endif
