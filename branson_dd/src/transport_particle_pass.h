/*
  Author: Alex Long
  Date: 12/1/2015
  Name: transport_mesh_pass.h
*/

#ifndef transport_particle_pass_h_
#define transport_particle_pass_h_

#include <algorithm>
#include <vector>
#include <numeric>
#include <queue>
#include <boost/mpi.hpp>

#include "constants.h"
#include "mesh.h"
#include "mesh_cell_pass.h"
#include "sampling_functions.h"
#include "RNG.h"
#include "request.h"

namespace mpi = boost::mpi;

bool transport_single_photon( Photon* iphtn,
                              Mesh_Cell_Pass* mesh,
                              RNG* rng,
                              double& next_dt,
                              double& exit_E,
                              double& census_E,
                              unsigned int& census_count,
                              unsigned int& complete_count,
                              std::vector<double>& rank_abs_E)

{
  using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT; 
  using Constants::PROCESSOR;
  using Constants::finish_tag;
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
      complete_count++;
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
        if(boundary_event == ELEMENT )
          next_cell = cell.get_next_cell(surface_cross);
          iphtn->set_cell(next_cell);
          cell_id=next_cell;
        else if(boundary_event == PROCESSOR) {
          pass_flag = true;
          active=false;
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
        complete_count++;
        census_E+=iphtn->get_E();
      }
    } //end event loop
  } // end while alive
  return pass_flag;
}



void transport_photons(Photon*& photon_vec,
                        unsigned int n_photon,
                        Mesh_Particle_Pass* mesh,
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
  unsigned int complete_count = 0;
  double census_E=0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt();

  RNG *rng = imc_state->get_rng();
  Photon* iphtn;

  int n_rank =world.size();
  int rank   =world.rank();

  int parent = (rank + 1) / 2 - 1;
  int child1 = rank * 2 + 1;
  int child2 = child1 + 1;

  // set missing nodes to Constants::proc_null
  { 
    if (!rank)
        parent = proc_null;

    // maximum valid node id
    const int last_node = n_rank - 1;

    if (child1 > last_node)
    {
        child1 = proc_null;
        child2 = proc_null;
    }
    else if (child1 == last_node)
        child2 = proc_null;
  }

  vector<vector<Photon> > pass_list;
  vector<vector<Photon> >receive_list;
  unsigned int c1_count, c2_count, p_count;

  //Message requests for finished photon counts
  Request *r_c1_count = new Request();
  Request *r_c2_count = new Request();
  Request *r_p_count  = new Request();
  Request *s_c1_count = new Request();
  Request *s_c2_count = new Request();
  Request *s_p_count  = new Request();
  Request *r_photon   = new Request[n_rank-1];
  Request *s_photon   = new Request[n_rank-1];

  //Post receives for photon counts from parent and each child
  r_c1_count->request(world.irecv(child1, count_tag, c1_count));
  r_c2_count->request(world.irecv(child2, count_tag, c2_count));
  r_p_count->request(world.irecv(parent, count_tag, p_count));

  //Post receives for photons from other ranks
  //NOTE: This should just involve adjacent sub-domains
  // I can precalculate those
  for (unsigned int ir=0; ir<n_rank; i++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      unsigned int r_index = ir - (ir>rank);
      photon_recv_req[r_index]->request(world.irecv(ir, photon_tag, receive_list[r_index]));
    }
  }
 
  bool pass_flag = false;
  unsigned int n = imc_state->get_batch_size();
  
  vector<Photon> census_list;

  ////////////////////////////////////////////////////////////////////////
  // main transport loop
  ////////////////////////////////////////////////////////////////////////
  while (!finished) {


    while (n && received_bank.size()) {
      
    }   
    for (unsigned int i=p_start;i<p_end; i++) {
      iphtn = &photon_vec[i];
      pass_flag = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                          census_E, census_count, 
                                          complete_count, rank_abs_E);
      
      if (wait_flag) pass_list.push(iphtn);
    } // end for iphtn

    // with some frequency, check for requests and try to transport the
    // waiting list
    for (unsigned int ir=0; ir<n_rank; ir++) {
      unsigned int r_index = ir - (ir>rank);
      if (request[r_index]) {
        transport_list.insert(receive_list[r_index]);
        //post receive again
        photon_recv_req[r_index] = world.irecv(ir, photon_tag, receive_list[r_index]);
      }
      // send photons to other ranks
      if (send_photon[r_index]) {
        // if last send has been received, post new send
        if (photon_send_req[r_index])
          photon_send_req[r_index] = world.isend(ir, photon_tag, send_list[r_index]);
      }
    } // end loop over ranks

    //send number of completed to parent, children
    {
      if (r_c1_count.test()) c1 = true
      if (r_c2_count.test()) c2 = true

      if (c1 && c2) {
        tree_count = complete_count + c1_count + c2_count;
        c1 = false;
        c2 = false;
        r_c1_count->request(world.irecv(child1, count_tag, c1_count));
        r_c2_count->request(world.irecv(child2, count_tag, c2_count));
      }
      //if you're finished, send count down
      if (tree_count == global_photon_count) {
        finished = true;
        s_c1_count->request(world.isend(child1, count_tag, tree_count));
        s_c2_count->request(world.isend(child2, count_tag, tree_count));
      }

      if (rank) {
        //send to parent, if unsent
        if (!s_p_count.valid()) {
          s_c1_count->request(world.isend(parent, count_tag, tree_count));
        }
        //check for parent receive completion
        if (r_p_count.test()) {
          if (p_count == global_photon_count) {
            finished= true;
            s_c1_count->request(world.isend(child1, count_tag, tree_count));
            s_c2_count->request(world.isend(child2, count_tag, tree_count));
          }
        }
      } // if  non-root ranks
    }
    
    // set next loop start and end
    p_start=p_end;
    p_end = p_end + imc_state->get_batch_size();
    if (p_end > n_photon) p_end = n_photon;

  } // end while

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

#endif // def transport_particle_pass_h_
