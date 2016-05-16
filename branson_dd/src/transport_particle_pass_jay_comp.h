//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   transport_particle_pass_jay_comp.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  Holds parameters needed in IMC simulation
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef transport_particle_pass_jay_comp_h_
#define transport_particle_pass_jay_comp_h_

#include <algorithm>
#include <mpi.h>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <map>
#include <stack>
#include <vector>

#include "completion_manager.h"
#include "constants.h"
#include "buffer.h"
#include "mesh.h"
#include "sampling_functions.h"
#include "transport_particle_pass.h"
#include "RNG.h"
#include "photon.h"


std::vector<Photon> jay_comp_transport_particle_pass(Source& source,
                                                      Mesh* mesh,
                                                      IMC_State* imc_state,
                                                      IMC_Parameters* imc_parameters,
                                                      Completion_Manager* comp,
                                                      std::vector<double>& rank_abs_E)
{
  using Constants::event_type;
  using Constants::PASS; using Constants::CENSUS;
  using Constants::KILL; using Constants::EXIT;
  using Constants::WAIT;
  using Constants::photon_tag;
  using std::queue;
  using std::vector;
  using std::stack;
  using std::map;
  using Constants::proc_null;
  using Constants::count_tag;
  using std::cout;
  using std::endl;

  double census_E=0.0;
  double exit_E = 0.0;
  double next_dt = imc_state->get_next_dt(); //!< Set for census photons
  double dt = imc_state->get_next_dt(); //<! For making current photons

  RNG *rng = imc_state->get_rng();

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // parallel event counters
  uint32_t n_photon_messages=0; //! Number of photon messages
  uint32_t n_photons_sent=0; //! Number of photons passed
  uint32_t n_sends_posted=0; //! Number of sent messages posted
  uint32_t n_sends_completed=0; //! Number of sent messages completed
  uint32_t n_receives_posted=0; //! Number of received messages completed
  uint32_t n_receives_completed=0; //! Number of received messages completed
  uint64_t n_complete_messages=0; //! Number of complete messages
  uint64_t n_photon_tests=0; //! Number of receive tests for photon lists

  // Number of particles to run between MPI communication 
  const uint32_t batch_size = imc_parameters->get_batch_size();

  // Preferred size of MPI message
  const uint32_t max_buffer_size 
    = imc_parameters->get_particle_message_size();

  // MPI Particle type
  //Make the MPI_Particle type
  const int entry_count = 2 ; 
  // 7 uint32_t, 6 int, 13 double
  int array_of_block_length[3] = { 2, 9};
  // Displacements of each type in the cell
  MPI_Aint array_of_block_displace[2] = 
    {0, 2*sizeof(uint32_t)};
  //Type of each memory block
  MPI_Datatype array_of_types[2] = {MPI_UNSIGNED, MPI_DOUBLE};

  MPI_Datatype MPI_Particle;
  MPI_Type_create_struct(entry_count, array_of_block_length, 
    array_of_block_displace, array_of_types, &MPI_Particle);

  // Commit the type to MPI so it recognizes it in communication calls
  MPI_Type_commit(&MPI_Particle);

  //get global photon count
  uint64_t n_local = source.get_n_photon();
  uint64_t n_global;

  MPI_Allreduce(&n_local, &n_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
    MPI_COMM_WORLD);
 
  // post receives to children and parent for completion messages 
  comp->start_timestep(n_receives_posted);

  //set global particle count in completion object
  comp->set_timestep_global_particles(n_global);

  // This flag indicates that send processing is needed for target rank
  vector<vector<Photon> > send_list;

  //Get adjacent processor map (off_rank_id -> adjacent_proc_number)
  map<uint32_t, uint32_t> adjacent_procs = mesh->get_proc_adjacency_list();
  uint32_t n_adjacent = adjacent_procs.size();
  //Messsage requests for photon sends and receives
  MPI_Request *phtn_recv_request   = new MPI_Request[n_adjacent];
  MPI_Request *phtn_send_request   = new MPI_Request[n_adjacent];
  // make a send/receive particle buffer for each adjacent processor
  vector<Buffer<Photon> > phtn_recv_buffer(n_adjacent);
  vector<Buffer<Photon> > phtn_send_buffer(n_adjacent);

  //Post receives for photons from adjacent sub-domains
  {
    uint32_t i_b; // buffer index
    int adj_rank; // adjacent rank
    for ( std::map<uint32_t, uint32_t>::iterator it=adjacent_procs.begin(); 
      it != adjacent_procs.end(); ++it) {
      adj_rank = it->first;
      i_b = it->second;
      //push back send and receive lists
      vector<Photon> empty_phtn_vec;
      send_list.push_back(empty_phtn_vec);
      //make receive buffer the appropriate size
      phtn_recv_buffer[i_b].resize(max_buffer_size);
      MPI_Irecv(phtn_recv_buffer[i_b].get_buffer(),
        max_buffer_size,
        MPI_Particle,
        adj_rank,
        photon_tag,
        MPI_COMM_WORLD,
        &phtn_recv_request[i_b]);
      n_receives_posted++;
      phtn_recv_buffer[i_b].set_awaiting();
    } // end loop over adjacent processors
  }

  ////////////////////////////////////////////////////////////////////////
  // main transport loop
  ////////////////////////////////////////////////////////////////////////

  vector<Photon> census_list; //!< End of timestep census list
  stack<Photon> phtn_recv_stack; //!< Stack of received photons

  int send_rank;
  uint64_t n_complete = 0; //!< Completed histories, regardless of origin
  uint64_t n_local_sourced = 0; //!< Photons pulled from source object
  bool finished = false;
  bool from_receive_stack = false;
  Photon phtn;
  event_type event;

  while (!finished) {

    uint32_t n = batch_size;
    
    ////////////////////////////////////////////////////////////////////////////
    // Transport photons from source and received list
    ////////////////////////////////////////////////////////////////////////////
    //first, try to transport photons from the received list
    while (n && (!phtn_recv_stack.empty() || (n_local_sourced < n_local))) {
      
      if (!phtn_recv_stack.empty()) {
        phtn = phtn_recv_stack.top();
        from_receive_stack=true;
      }
      else {
        phtn =source.get_photon(rng, dt); 
        n_local_sourced++;
        from_receive_stack=false;
      }

      event = transport_photon_particle_pass(phtn, mesh, rng, next_dt, exit_E,
                                            census_E, rank_abs_E);
      switch(event) {
        // this case should never be reached
        case WAIT:
          break;
        case KILL: 
          n_complete++;
          break;
        case EXIT:
          n_complete++;
          break;
        case CENSUS:
          census_list.push_back(phtn);
          n_complete++;
          break;
        case PASS:
          send_rank = mesh->get_rank(phtn.get_cell());
          int i_b = adjacent_procs[send_rank];
          send_list[i_b].push_back(phtn);
          break;
      }
      n--;
      if (from_receive_stack) phtn_recv_stack.pop();
    }

    ////////////////////////////////////////////////////////////////////////////
    // process photon send and receives 
    ////////////////////////////////////////////////////////////////////////////
    {
      int send_req_flag;
      int recv_req_flag;
      int recv_count; // recieve count is 32 bit

      MPI_Status recv_status;
      uint32_t i_b; // buffer index
      int adj_rank; // adjacent rank
      for ( std::map<uint32_t, uint32_t>::iterator it=adjacent_procs.begin(); 
        it != adjacent_procs.end(); ++it) {
        adj_rank = it->first;
        i_b = it->second;
        // test completion of send buffer
        if (phtn_send_buffer[i_b].sent()) {
          MPI_Test(&phtn_send_request[i_b], &send_req_flag, MPI_STATUS_IGNORE);
          if (send_req_flag) {
            phtn_send_buffer[i_b].reset();
            n_sends_completed++;
          } 
        }

        // send full photon buffers if send_buffer is empty and send_list has
        // some photons in it
        if ( (phtn_send_buffer[i_b].empty() && !send_list[i_b].empty()) &&
          (send_list[i_b].size() >= max_buffer_size || n_local_sourced == n_local) ) {
          uint32_t n_photons_to_send = max_buffer_size;
          if ( send_list[i_b].size() < max_buffer_size) 
            n_photons_to_send = send_list[i_b].size();
          vector<Photon>::iterator copy_start = send_list[i_b].begin();
          vector<Photon>::iterator copy_end = send_list[i_b].begin()+n_photons_to_send;
          vector<Photon> send_now_list(copy_start, copy_end);
          send_list[i_b].erase(copy_start,copy_end); 
          phtn_send_buffer[i_b].fill(send_now_list);
          n_photons_sent += n_photons_to_send;
          MPI_Isend(phtn_send_buffer[i_b].get_buffer(), 
            n_photons_to_send,  
            MPI_Particle, 
            adj_rank, 
            photon_tag, 
            MPI_COMM_WORLD, 
            &phtn_send_request[i_b]);
          n_sends_posted++;
          phtn_send_buffer[i_b].set_sent();
          n_photon_messages++;
        }

        //process receive buffer
        if (phtn_recv_buffer[i_b].awaiting()) {
          n_photon_tests++;
          MPI_Test(&phtn_recv_request[i_b], &recv_req_flag, &recv_status);
          if (recv_req_flag) {
            n_receives_completed++;
            vector<Photon>& receive_list = phtn_recv_buffer[i_b].get_object();
            // only push the number of received photons onto the recv_stack
            MPI_Get_count(&recv_status, MPI_Particle, &recv_count);
            for (uint32_t i=0; i<recv_count; i++) 
              phtn_recv_stack.push(receive_list[i]);
            phtn_recv_buffer[i_b].reset();
            //post receive again, don't resize--it's already set to maximum
            MPI_Irecv(phtn_recv_buffer[i_b].get_buffer(),
              max_buffer_size,
              MPI_Particle,
              adj_rank,
              photon_tag,
              MPI_COMM_WORLD,
              &phtn_recv_request[i_b]);
            n_receives_posted++;
            phtn_recv_buffer[i_b].set_awaiting();
          }
        }
      } // end loop over adjacent processors
    } //end scope of particle passing

    ////////////////////////////////////////////////////////////////////////////
    // binary tree completion communication
    ////////////////////////////////////////////////////////////////////////////
    comp->check_messages(n_complete, n_receives_posted, n_receives_completed, 
      n_sends_posted);
    if ( n_complete  &&  (n_local==n_local_sourced && 
      phtn_recv_stack.empty()) ) {
      comp->send_parent_n_tree_complete(n_complete, n_sends_posted);
    }
    finished = comp->is_finished(n_complete);
  } // end while

  comp->end_timestep(n_sends_posted, n_sends_completed, n_receives_posted,
    n_receives_completed);

  // wait for all ranks to finish then send empty photon messages.
  // Do this because it's possible for a rank to receive the empty message
  // while it's still in the transport loop. In that case, it will post a 
  // receive again, which will never have a matching send
  MPI_Barrier(MPI_COMM_WORLD);

  //finish off posted photon receives
  {
    vector<Photon> one_photon(1);
    uint32_t i_b; // buffer index
    int adj_rank; // adjacent rank
    for ( std::map<uint32_t, uint32_t>::iterator it=adjacent_procs.begin(); 
      it != adjacent_procs.end(); ++it) {
      adj_rank = it->first;
      i_b = it->second;
      //wait for completion of previous sends
      if (phtn_send_buffer[i_b].sent()) 
        MPI_Wait(&phtn_send_request[i_b], MPI_STATUS_IGNORE);
      //send one photon vector to finish off receives, these photons will not 
      //be processed by the receiving ranks (all ranks are out of transport)
      MPI_Isend(&one_photon[0], 1, MPI_Particle, adj_rank,
        photon_tag, MPI_COMM_WORLD, &phtn_send_request[i_b]);
      n_sends_posted++;
      MPI_Wait(&phtn_send_request[i_b], MPI_STATUS_IGNORE);
      n_sends_completed++;
    } // end loop over adjacent processors
  }

  // Wait for receive requests
  for (uint32_t i_b=0; i_b<n_adjacent; i_b++) {
    MPI_Wait(&phtn_recv_request[i_b], MPI_STATUS_IGNORE);
    n_receives_completed++;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Reduce diagnostic quantities
  uint64_t g_complete_messages=0;
  uint64_t g_photon_tests=0;
  MPI_Allreduce(&n_complete_messages, &g_complete_messages, 1, MPI_UNSIGNED_LONG,
    MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&n_photon_tests, &g_photon_tests, 1, MPI_UNSIGNED_LONG,
    MPI_SUM, MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());
  //All ranks have now finished transport
  delete[] phtn_recv_request;
  delete[] phtn_send_request;

  imc_state->set_exit_E(exit_E);
  imc_state->set_post_census_E(census_E);
  imc_state->set_census_size(census_list.size());
  //set diagnostic
  imc_state->set_step_particle_messages(n_photon_messages);
  imc_state->set_step_particles_sent(n_photons_sent);
  imc_state->set_step_sends_posted(n_sends_posted);
  imc_state->set_step_sends_completed(n_sends_completed);
  imc_state->set_step_receives_posted(n_receives_posted);
  imc_state->set_step_receives_completed(n_receives_completed);

  return census_list;
}

#endif // def transport_particle_pass_jay_comp_h_
//---------------------------------------------------------------------------//
// end of transport_particle_pass_jay_comp.h
//---------------------------------------------------------------------------//
