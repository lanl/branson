/*
  Author: Alex Long
  Date: 12/1/2015
  Name: transport_mesh_pass.h
*/

#ifndef transport_particle_pass_h_
#define transport_particle_pass_h_

#include <algorithm>
#include <vector>
#include <stack>
#include <numeric>
#include <queue>
#include <boost/mpi.hpp>

#include "constants.h"
#include "buffer.h"
#include "mesh.h"
#include "sampling_functions.h"
#include "RNG.h"
#include "request.h"
#include "photon.h"

namespace mpi = boost::mpi;

Constants::event_type transport_single_photon(Photon& phtn,
                                              Mesh* mesh,
                                              RNG* rng,
                                              double& next_dt,
                                              double& exit_E,
                                              double& census_E,
                                              std::vector<double>& rank_abs_E)
{
  using Constants::VACUUM; using Constants::REFLECT; 
  using Constants::ELEMENT; using Constants::PROCESSOR;
  using Constants::PASS; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;
  using Constants::bc_type;
  using Constants::event_type;
  using Constants::c;
  using std::min;

  unsigned int cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E;
  double angle[3];
  Cell cell;

  unsigned int surface_cross = 0;
  double cutoff_fraction = 0.001; //note: get this from IMC_state

  cell_id=phtn.get_cell();
  cell = mesh->get_on_rank_cell(cell_id);
  bool active = true;

  //transport this photon
  while(active) {
    sigma_a = cell.get_op_a();
    sigma_s = cell.get_op_s();
    f = cell.get_f();

    //get distance to event
    dist_to_scatter = -log(rng->generate_random_number())/((1.0-f)*sigma_a + sigma_s);
    dist_to_boundary = cell.get_distance_to_boundary(phtn.get_position(),
                                                      phtn.get_angle(),
                                                      surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    //select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    //Calculate energy absorbed by material, update photon and material energy
    absorbed_E = phtn.get_E()*(1.0 - exp(-sigma_a*f*dist_to_event));
    phtn.set_E(phtn.get_E() - absorbed_E);

    rank_abs_E[cell_id] += absorbed_E;
    
    //update position
    phtn.move(dist_to_event);

    //Apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      rank_abs_E[cell_id] += phtn.get_E();
      phtn.set_dead();
      active=false;
      event=KILL;
    }
    // or apply event
    else {
      //Apply event
      //EVENT TYPE: SCATTER
      if(dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
      }
      //EVENT TYPE: BOUNDARY CROSS
      else if(dist_to_event == dist_to_boundary) {
        boundary_event = cell.get_bc(surface_cross);
        if(boundary_event == ELEMENT ) {
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          cell_id=next_cell;
          cell = mesh->get_on_rank_cell(cell_id);
        }
        else if(boundary_event == PROCESSOR) {
          active=false;
          //set correct cell index with global cell ID
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          event=PASS;
        }
        else if(boundary_event == VACUUM) {
          exit_E+=phtn.get_E();
          active=false; 
          event = EXIT;
        }
        else phtn.reflect(surface_cross); 
      }
      //EVENT TYPE: REACH CENSUS
      else if(dist_to_event == dist_to_census) {
        phtn.set_census_flag(true);
        phtn.set_distance_to_census(c*next_dt);
        active=false;
        event=CENSUS;
        census_E+=phtn.get_E();
      }
    } //end event loop
  } // end while alive
  return event;
}



std::vector<Photon> transport_photons(Source& source,
                                      Mesh* mesh,
                                      IMC_State* imc_state,
                                      std::vector<double>& rank_abs_E,
                                      int chk_freq,
                                      mpi::communicator world)
{
  using Constants::event_type;
  using Constants::PASS; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;
  using Constants::photon_tag;
  using std::queue;
  using std::vector;
  using std::stack;
  using Constants::proc_null;
  using Constants::count_tag;

  double census_E=0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt(); //!< Set for census photons
  double dt = imc_state->get_next_dt(); //<! For making current photons

  RNG *rng = imc_state->get_rng();

  int n_rank = world.size();
  int rank   = world.rank();

  //get global photon count
  unsigned int n_local = source.get_n_photon();
  unsigned int n_global;
  MPI::COMM_WORLD.Allreduce(&n_local, 
                            &n_global, 
                           1, 
                           MPI_UNSIGNED, 
                           MPI_SUM);


  int parent = (rank + 1) / 2 - 1;
  int child1 = rank * 2 + 1;
  int child2 = child1 + 1;

  // set missing nodes to proc_null
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

  // This flag indicates that send processing is needed for target rank
  vector<bool> process_send_flag(n_rank-1,false);
  vector<vector<Photon> > send_list;

  //Message requests for finished photon counts
  mpi::request c1_recv_request;
  mpi::request c2_recv_request;
  mpi::request p_recv_request;
  mpi::request c1_send_request;
  mpi::request c2_send_request; 
  mpi::request p_send_request;
  mpi::request *phtn_recv_request   = new mpi::request[n_rank-1];
  mpi::request *phtn_send_request   = new mpi::request[n_rank-1];

  //Buffers
  Buffer<unsigned int> c1_recv_buffer;
  Buffer<unsigned int> c2_recv_buffer;
  Buffer<unsigned int> p_recv_buffer;
  Buffer<unsigned int> c1_send_buffer;
  Buffer<unsigned int> c2_send_buffer;
  Buffer<unsigned int> p_send_buffer;
  vector<Buffer<Photon> > phtn_recv_buffer(n_rank-1);
  vector<Buffer<Photon> > phtn_send_buffer(n_rank-1);

  //Post receives for photon counts from parent and each child
  if (child1 != proc_null) {
    c1_recv_request = world.irecv(child1, count_tag, c1_recv_buffer.get_buffer());  
    c1_recv_buffer.set_awaiting();
  }
  if (child2 != proc_null) {
    c2_recv_request = world.irecv(child2, count_tag, c2_recv_buffer.get_buffer());
    c2_recv_buffer.set_awaiting();
  }
  if (parent != proc_null) {
    p_recv_request = world.irecv(parent, count_tag, p_recv_buffer.get_buffer() );
    p_recv_buffer.set_awaiting();
  }

  //Post receives for photons from other ranks
  //NOTE: This should just involve adjacent sub-domains
  // I can precalculate those
  for (int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //push back send and receive lists
      vector<Photon> empty_phtn_vec;
      send_list.push_back(empty_phtn_vec);
      //get correct index into requests and vectors 
      int r_index = ir - (ir>rank);
      phtn_recv_request[r_index] = world.irecv(ir, photon_tag, phtn_recv_buffer[r_index].get_buffer());
      phtn_recv_buffer[r_index].set_awaiting();
    }
  }
 
  ////////////////////////////////////////////////////////////////////////
  // main transport loop
  ////////////////////////////////////////////////////////////////////////

  vector<Photon> census_list; //end of timestep census list
  stack<Photon> phtn_recv_stack; //stack of received photons

  int send_rank;
  unsigned int n_complete = 0; //!< Completed histories, regardless of origin
  unsigned int n_local_sourced = 0; //!< Photons pulled from source object
  unsigned int tree_count = 0; //!< Total for this node and all children
  bool finished = false;
  bool from_receive_stack = false;
  Photon iphtn;
  event_type event;

  while (!finished) {

    //unsigned int n = imc_state->get_batch_size();
    //hardcoded for now, fix later
    unsigned int n = 100;

    //first, try to transport photons from the received list
    while (n && (!phtn_recv_stack.empty() || (n_local_sourced < n_local))) {
      
      if (!phtn_recv_stack.empty()) {
        iphtn = phtn_recv_stack.top();
        from_receive_stack=true;
      }
      else {
        iphtn =source.get_photon(rng, dt); 
        n_local_sourced++;
        from_receive_stack=false;
      }

      event = transport_single_photon(iphtn, mesh, rng, next_dt, exit_E,
                                          census_E, rank_abs_E);
      switch(event) {
        case KILL: 
          n_complete++;
          break;
        case EXIT:
          n_complete++;
          break;
        case CENSUS:
          census_list.push_back(iphtn);
          n_complete++;
          break;
        case PASS:
          send_rank = mesh->get_rank(iphtn.get_cell());
          int r_index = send_rank - (send_rank>rank);
          send_list[r_index].push_back(iphtn);
          break;
      }
      n--;
      if (from_receive_stack) phtn_recv_stack.pop();
    }

    // process photon send and receives 
    for (int ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        int r_index = ir - (ir>rank);
       
        // process send buffer
        if (phtn_send_buffer[r_index].sent()) {
          if (phtn_send_request[r_index].test()) {
            phtn_send_buffer[r_index].reset();
          } 
        }

        if (phtn_send_buffer[r_index].ready()) {
          phtn_send_request[r_index] = world.isend(ir, photon_tag, phtn_send_buffer[r_index].get_buffer());
          phtn_send_buffer[r_index].set_sent();
        }

        if (phtn_send_buffer[r_index].empty() &&  !send_list[r_index].empty()) {
          phtn_send_buffer[r_index].fill(send_list[r_index]);
          send_list[r_index].clear();
        }

        //process receive buffer
        if (phtn_recv_buffer[r_index].awaiting()) {
          if (phtn_recv_request[r_index].test()) {
            vector<Photon> receive_list = phtn_recv_buffer[r_index].get_buffer();
            for (unsigned int i=0; i<receive_list.size(); i++) 
              phtn_recv_stack.push(receive_list[i]);
            phtn_recv_buffer[r_index].reset();
            //post receive again
            phtn_recv_request[r_index] = world.irecv(ir, photon_tag, phtn_recv_buffer[r_index].get_buffer());
            phtn_recv_buffer[r_index].set_awaiting();
          }
        }
      }
    } // end loop over ranks


    //send number of completed to parent, children
    {
      if (c1_recv_buffer.awaiting()) 
        if (c1_recv_request.test()) 
          c1_recv_buffer.set_received();
      if (c2_recv_buffer.awaiting()) 
        if (c2_recv_request.test()) 
          c2_recv_buffer.set_received();

      //if recieved from both children, compute tree count
      if ( ( c1_recv_buffer.received() || child1==proc_null) 
        && ( c2_recv_buffer.received() || child2==proc_null)) {
        unsigned int c1_count = 0; unsigned int c2_count = 0;
        if (child1 == proc_null) c1_count = 0;
        else c1_count = c1_recv_buffer.get_buffer()[0]; 
        if (child2 == proc_null) c2_count = 0;
        else c2_count = c2_recv_buffer.get_buffer()[0]; 
        tree_count = n_complete + c1_count + c2_count;
      }
      else tree_count = n_complete;

      //if received from parent, overwrite tree count
      if (p_recv_buffer.awaiting()) {
        if (p_recv_request.test()) {
          p_recv_buffer.set_received();
          tree_count = p_recv_buffer.get_buffer()[0];
        }
      }

      // If finished, set flag. Otherwise, post receives from children.
      if (tree_count == n_global) finished = true;
      // Otherwise, send messages up tree and receive from down tree
      else {
        if (child1 != proc_null &&  c1_recv_buffer.received())
          c1_recv_request = world.irecv(child1, count_tag, c1_recv_buffer.get_buffer());
        if (child2 != proc_null &&  c2_recv_buffer.received()) 
          c2_recv_request = world.irecv(child2, count_tag, c2_recv_buffer.get_buffer());
        if (rank) {
          //send to parent, if unsent
          if (p_send_buffer.sent()) {
            if (p_send_request.test()) p_send_buffer.reset();
          }
          if (p_send_buffer.empty() ) {
            p_send_buffer.fill(vector<unsigned int> (1,tree_count));
            p_send_request = world.isend(parent, count_tag, p_send_buffer.get_buffer());
            p_send_buffer.set_sent();
          }
        } // if  non-root ranks
      }

    } //end send block

  } // end while

  //once finished, send finished count down
  if (child1 != proc_null) {
    c1_send_buffer.fill(vector<unsigned int> (1,tree_count));
    c1_send_request = world.isend(child1, count_tag, c1_send_buffer.get_buffer());
  }
  if (child2 != proc_null)  {
    c2_send_buffer.fill(vector<unsigned int> (1,tree_count));
    c2_send_request = world.isend(child2, count_tag, c2_send_buffer.get_buffer());
  }

  /*
  for (int ir=0; ir<n_rank; ir++) {
    if (ir != rank) {
      //get correct index into requests and vectors 
      int r_index = ir - (ir>rank);
      send_buffer[r_index].clear();
      phtn_send_request[r_index] = world.isend(ir, photon_tag, send_buffer[r_index]);
    }
  }
  */  

  MPI::COMM_WORLD.Barrier();

  std::sort(census_list.begin(), census_list.end(), Photon::census_flag_compare);
  //All ranks have now finished transport
  delete[] phtn_recv_request;
  delete[] phtn_send_request;

  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_particle_pass_h_
