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
#include "decompose_photons.h"
#include "mesh.h"
#include "RNG.h"
#include "sampling_functions.h"

namespace mpi = boost::mpi;

Constants::event_type transport_photon_mesh_pass(Photon& phtn,
                              Mesh* mesh,
                              RNG* rng,
                              double& next_dt,
                              double& exit_E,
                              double& census_E,
                              std::vector<double>& rank_abs_E)

{
  using Constants::VACUUM; using Constants::REFLECT; 
  using Constants::ELEMENT; using Constants::PROCESSOR;
  //events
  using Constants::WAIT; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;
  using Constants::bc_type;
  using Constants::event_type;
  using Constants::c;
  using std::min;

  uint32_t cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E;
  double angle[3];
  Cell cell;

  uint32_t surface_cross = 0;
  const double cutoff_fraction = 0.01; //note: get this from IMC_state

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
        if(boundary_event == ELEMENT || boundary_event == PROCESSOR) {
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          cell_id=next_cell;
          //look for this cell, if it's not there transport later
          if (mesh->mesh_available(cell_id))
            cell = mesh->get_on_rank_cell(cell_id);
          else {
            event= WAIT;
            active=false;
          }
        }
        else if(boundary_event == VACUUM) {
          active=false; 
          exit_E+=phtn.get_E();
          event=EXIT;
        }
        else phtn.reflect(surface_cross); 
      }
      //EVENT TYPE: REACH CENSUS
      else if(dist_to_event == dist_to_census) {
        phtn.set_census_flag(true);
        phtn.set_distance_to_census(c*next_dt);
        active=false;
        census_E+=phtn.get_E();
        event=CENSUS;
      }
    } //end event loop
  } // end while alive
  return event;
}



std::vector<Photon> transport_mesh_pass(Source& source,
                                        Mesh* mesh,
                                        IMC_State* imc_state,
                                        IMC_Parameters* imc_parameters,
                                        std::vector<double>& rank_abs_E,
                                        mpi::communicator world)
{
  using Constants::finish_tag;
  using std::queue;
  using std::vector;
  using Constants::proc_null;
  using Constants::event_type;
  //events
  using Constants::WAIT; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;

  uint32_t n_local = source.get_n_photon();
  uint32_t n_local_sourced = 0;

  uint32_t cell_id;
  double census_E = 0.0;
  double exit_E = 0.0;
  double dt = imc_state->get_next_dt(); //<! For making current photons
  double next_dt = imc_state->get_next_dt(); //<! For census photons

  RNG *rng = imc_state->get_rng();
  Photon phtn;

  int n_rank =world.size();
  int rank   =world.rank();

  // parallel event counters
  uint32_t n_cell_messages=0; //! Number of cell messages
  uint32_t n_cells_sent=0; //! Number of cells passed
  uint32_t n_sends_posted=0; //! Number of sent messages posted
  uint32_t n_sends_completed=0; //! Number of sent messages completed
  uint32_t n_receives_posted=0; //! Number of received messages completed
  uint32_t n_receives_completed=0; //! Number of received messages completed

  // set parent and children for binary tree communication of finish
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


  //Message requests for finished flags
  mpi::request c1_recv_request;
  mpi::request c2_recv_request;
  mpi::request p_recv_request;
  mpi::request c1_send_request;
  mpi::request c2_send_request; 
  mpi::request p_send_request;

  //Buffers for finished flags
  Buffer<bool> c1_recv_buffer;
  Buffer<bool> c2_recv_buffer;
  Buffer<bool> p_recv_buffer;
  Buffer<bool> c1_send_buffer;
  Buffer<bool> c2_send_buffer;
  Buffer<bool> p_send_buffer;

  // finished message propagates from children to parent, then back down the
  // tree

  //post finish message receives from children and parents
  if (parent != proc_null) {
    p_recv_request = world.irecv(parent, finish_tag, p_recv_buffer.get_buffer() );
    n_receives_posted++;
    p_recv_buffer.set_awaiting();
  }
  if (child1 != proc_null) {
    c1_recv_request = world.irecv(child1, finish_tag, c1_recv_buffer.get_buffer() );
    n_receives_posted++;
    c1_recv_buffer.set_awaiting();
  }
  if (child2 != proc_null) {
    c2_recv_request = world.irecv(child2, finish_tag, c2_recv_buffer.get_buffer() );
    n_receives_posted++;
    c2_recv_buffer.set_awaiting();
  }

  // New data flag is initially false
  bool new_data = false;
  // Number of particles to run between MPI communication 
  const uint32_t batch_size = imc_parameters->get_batch_size();

  event_type event;
  uint32_t wait_list_size;

  ////////////////////////////////////////////////////////////////////////
  // main loop over photons
  ////////////////////////////////////////////////////////////////////////
  vector<Photon> census_list; //!< Local end of timestep census list
  vector<Photon> off_rank_census_list; //!< Off rank end of timestep census list
  queue<Photon> wait_list; //!< Photons waiting for mesh data 
  while ( n_local_sourced < n_local) {
    
    uint32_t n = batch_size;

    while (n && n_local_sourced < n_local) {

      phtn =source.get_photon(rng, dt); 
      n_local_sourced++;

      //get start cell, this only changea with cell crossing event
      cell_id=phtn.get_cell();

      // if mesh available, transport and process, otherwise put on the
      // waiting list
      if (mesh->mesh_available(cell_id)) {
        event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                        census_E, rank_abs_E);
        cell_id = phtn.get_cell();
      }
      else event = WAIT;

      if (event==CENSUS) { 
        if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
        else off_rank_census_list.push_back(phtn);
      }
      else if (event==WAIT) {
        mesh->request_cell(cell_id);
        wait_list.push(phtn);
      }
      n--;
    } // end batch transport

    //process mesh requests
    new_data = mesh->process_mesh_requests(world, n_cell_messages, n_cells_sent,
                                            n_sends_posted, n_sends_completed,
                                            n_receives_posted, n_receives_completed);
    // if data was received, try to transport photons on waiting list
    if (new_data) {
      wait_list_size = wait_list.size();
      for (uint32_t wp =0; wp<wait_list_size; wp++) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id=phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                          census_E, rank_abs_E);
          cell_id = phtn.get_cell();
          if (event==CENSUS) { 
            if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
            else off_rank_census_list.push_back(phtn);
          }
          else if (event==WAIT) {
            mesh->request_cell(cell_id);
            wait_list.push(phtn);
          }
        }
        else wait_list.push(phtn);
      } // end wp in wait_list
    }

  } //end while (n_local_source < n_local)

  ////////////////////////////////////////////////////////////////////////
  // Main transport loop finished, transport photons waiting for data
  ////////////////////////////////////////////////////////////////////////
  while (!wait_list.empty()) {
    new_data = mesh->process_mesh_requests(world, n_cell_messages, n_cells_sent,
                                            n_sends_posted, n_sends_completed,
                                            n_receives_posted, n_receives_completed);
    // if new data received, transport waiting list 
    if (new_data) {
      wait_list_size = wait_list.size();
      for (uint32_t wp =0; wp<wait_list_size; wp++) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id=phtn.get_cell();
        if (mesh->mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                              census_E, rank_abs_E);
          cell_id = phtn.get_cell();
          if (event==CENSUS) { 
            if (mesh->on_processor(cell_id)) census_list.push_back(phtn);
            else off_rank_census_list.push_back(phtn);
          }
          else if (event==WAIT) {
            mesh->request_cell(cell_id);
            wait_list.push(phtn);
          }
        }
        else {
          wait_list.push(phtn);
          mesh->request_cell(cell_id);
        }
      }
    }
  } //end while wait_list not empty

  vector<bool> s_bool(1,true);

  ////////////////////////////////////////////////////////////////////////
  // While waiting for other ranks to finish, check for other messages
  ////////////////////////////////////////////////////////////////////////
  bool c1_finished = false || child1==proc_null;
  bool c2_finished = false || child2==proc_null;
  bool p_finished  = false || parent==proc_null;

  while (!((c1_finished && c2_finished) && p_finished)) {
    if (c1_recv_buffer.awaiting()) {
      if (c1_recv_request.test()) {
        n_receives_completed++;
        c1_recv_buffer.set_received();
        c1_finished = c1_recv_buffer.get_buffer()[0];
      }
    }
    if (c2_recv_buffer.awaiting()) {
      if (c2_recv_request.test()) {
        n_receives_completed++;
        c2_recv_buffer.set_received();
        c2_finished = c2_recv_buffer.get_buffer()[0];
      }
    }
    if (p_recv_buffer.awaiting()) {
      if (p_recv_request.test()) {
        n_receives_completed++;
        p_recv_buffer.set_received();
        p_finished = p_recv_buffer.get_buffer()[0];
      }
    }
    
    // if children complete or terminating branch, send to parent (only once)
    if ((c1_finished && c2_finished) && 
        (parent != proc_null && !p_send_buffer.sent()) ) {
      p_send_buffer.fill(s_bool);
      p_send_request = world.isend(parent, finish_tag, p_send_buffer.get_buffer());
      n_sends_posted++;
      p_send_buffer.set_sent();
    }
    
    mesh->process_mesh_requests(world, n_cell_messages, n_cells_sent,
                                n_sends_posted, n_sends_completed,
                                n_receives_posted, n_receives_completed);
  } //end while

  // wait for completion of send
  if (p_send_buffer.sent()) {
    p_send_request.wait();
    n_sends_completed++;
  }

  //send complete message down the tree
  if (child1 != proc_null) { 
    c1_send_buffer.fill(s_bool);
    c1_send_request = world.isend(child1, finish_tag, c1_send_buffer.get_buffer());
    n_sends_posted++;
    c1_send_request.wait();
    n_sends_completed++;
  }
  if (child2 != proc_null)  {
    c2_send_buffer.fill(s_bool);
    c2_send_request = world.isend(child2, finish_tag, c2_send_buffer.get_buffer());
    n_sends_posted++;
    c2_send_request.wait();
    n_sends_completed++;
  }

  //wait for all ranks to finish transport to finish off cell and cell id
  //requests and sends
  MPI::COMM_WORLD.Barrier();
  mesh->finish_mesh_pass_messages(world, n_sends_posted, n_sends_completed,
                                  n_receives_posted, n_receives_completed);

  MPI::COMM_WORLD.Barrier();

  //All ranks have now finished transport

  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_step_cell_messages(n_cell_messages);
  imc_state->set_step_cells_sent(n_cells_sent);
  imc_state->set_step_sends_posted(n_sends_posted);
  imc_state->set_step_sends_completed(n_sends_completed);
  imc_state->set_step_receives_posted(n_receives_posted);
  imc_state->set_step_receives_completed(n_receives_completed);

  //send the off-rank census back to ranks that own the mesh its on.
  //receive census particles that are on your mesh
  vector<Photon> rebalanced_census = rebalance_census(off_rank_census_list,
                                                      mesh,
                                                      world);
  census_list.insert(census_list.end(), 
    rebalanced_census.begin(), 
    rebalanced_census.end());
  
  //sort on census vectors by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // set post census size after sorting and merging
  imc_state->set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_mesh_pass_h_
