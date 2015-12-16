/*
  Author: Alex Long
  Date: 6/29/2014
  Name: transport_mesh_pass.h
*/

#ifndef transport_mesh_pass_h_
#define transport_mesh_pass_h_

#include <algorithm>
#include <vector>
#include <numeric>
#include <queue>
#include <boost/mpi.hpp>

#include "constants.h"
#include "mesh.h"
#include "sampling_functions.h"
#include "RNG.h"
#include "request.h"

namespace mpi = boost::mpi;

bool transport_single_photon( Photon* iphtn,
                              Mesh* mesh,
                              RNG* rng,
                              double& next_dt,
                              double& exit_E,
                              double& census_E,
                              unsigned int& census_count,
                              std::vector<double>& rank_abs_E)

{
  using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT; 
  using Constants::PROCESSOR;
  using Constants::bc_type;
  using Constants::c;
  using std::min;

  unsigned int cell_id, next_cell;
  bc_type boundary_event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E;
  double angle[3];
  Cell cell;

  unsigned int surface_cross = 0;
  double cutoff_fraction = 0.001; //note: get this from IMC_state

  cell_id=iphtn->get_cell();
  cell = mesh->get_on_rank_cell(cell_id);
  bool active = true;
  bool wait_flag = false;
  //transport this photon
  while(active) {
    sigma_a = cell.get_op_a();
    sigma_s = cell.get_op_s();
    f = cell.get_f();

    //get distance to event
    dist_to_scatter = -log(rng->generate_random_number())/((1.0-f)*sigma_a + sigma_s);
    dist_to_boundary = cell.get_distance_to_boundary(iphtn->get_position(),
                                                      iphtn->get_angle(),
                                                      surface_cross);
    dist_to_census = iphtn->get_distance_remaining();

    //select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    //Calculate energy absorbed by material, update photon and material energy
    absorbed_E = iphtn->get_E()*(1.0 - exp(-sigma_a*f*dist_to_event));
    iphtn->set_E(iphtn->get_E() - absorbed_E);

    rank_abs_E[cell_id] += absorbed_E;
    
    //update position
    iphtn->move(dist_to_event);

    //Apply variance/runtime reduction
    if (iphtn->below_cutoff(cutoff_fraction)) {
      rank_abs_E[cell_id] += iphtn->get_E();
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
        boundary_event = cell.get_bc(surface_cross);
        if(boundary_event == ELEMENT || boundary_event == PROCESSOR) {
          next_cell = cell.get_next_cell(surface_cross);
          iphtn->set_cell(next_cell);
          cell_id=next_cell;
          //look for this cell, if it's not there transport later
          if (mesh->mesh_available(cell_id))
            cell = mesh->get_on_rank_cell(cell_id);
          else {
            mesh->request_cell(cell_id);
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
                        std::vector<double>& rank_abs_E,
                        Photon*& census_list,
                        int chk_freq,
                        mpi::communicator world)
{
  using Constants::finish_tag;
  using std::queue;
  using std::vector;

  unsigned int cell_id;
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
  for (int ir=0; ir<n_rank-1; ir++) {
    vector<bool> empty_bool_cell;
    r_finished.push_back(empty_bool_cell); 
  }

  mpi::request *s_finished_reqs = new mpi::request[ (n_rank-1)];
  mpi::request *r_finished_reqs = new mpi::request[ (n_rank-1)];

  // Post the receive calls for finished message
  for (int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      int r_index = ir - (ir>rank);
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
    //get start cell, only change when cell crossing event
    cell_id=iphtn->get_cell();

    //should always return true for photons before transport
    // but in the future we'll want to check this
    if (mesh->mesh_available(cell_id)) {
      wait_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                          census_E, census_count, rank_abs_E);
      if (wait_flag) wait_list.push(iphtn);
    } 
    else {
      mesh->request_cell(cell_id);
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
          cell_id=iphtn->get_cell();
          if (mesh->mesh_available(cell_id)) {
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
      cell_id=iphtn->get_cell();
      if (mesh->mesh_available(cell_id)) {
        wait_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                            census_E, census_count, rank_abs_E);
        if (wait_flag) wait_list.push(iphtn);
      }
      else {
        wait_list.push(iphtn);
        mesh->request_cell(cell_id);
      }
    }
  } //end while wait_list not empty

  // This rank is finished transporting, post finished
  vector<bool> s_bool(1,true);
  for (int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      int r_index = ir - (ir>rank);
      s_finished_reqs[r_index] = world.isend(ir, finish_tag, s_bool);
    }
  }

  // While waiting for other ranks to finish, check for other messages
  bool all_finished = false;
  while (!all_finished) {
    mesh->process_mesh_requests(world);
    for (int ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        int r_index = ir - (ir>rank);
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

  std::sort(photon_vec, photon_vec+n_photon, Photon::census_flag_compare);
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

#endif // def transport_mesh_pass_h_
