/*
  completion_manager.h
  by Alex Long
  3/31/2016
  The number of completed particles is sent up the chain and then reset.
  This allows us to send the completed count up the tree without
  trying to synchronize completion from both children. The root
  never resets the tree count
*/

#ifndef completion_manager_h_
#define completion_manager_h_

#include <iostream>
#include <mpi.h>
#include <vector>

#include "constants.h"

class Completion_Manager
{

  public:
  Completion_Manager(const int& rank,
                     const int& n_rank)
    : n_complete_c1(0),
      n_complete_c2(0),
      n_complete_p(0),
      n_particle_global(0),
      finished(false)
  {
    using Constants::proc_null;
    //set up binary tree rank structure
    parent = (rank + 1) / 2 - 1;
    child1 = rank * 2 + 1;
    child2 = child1 + 1;
    // set missing nodes to proc_null
    if (!rank) parent = proc_null;

    // maximum valid node id
    const int last_node = n_rank - 1;

    if (child1 > last_node) {
      child1 = proc_null;
      child2 = proc_null;
    }
    else if (child1 == last_node) child2 = proc_null;

    c1_recv_buffer.resize(1);
    c2_recv_buffer.resize(1);
    p_recv_buffer.resize(1); 
  }
  ~Completion_Manager() {}

  //non-const functions
  void set_timestep_global_particles(uint64_t _n_particle_global) {
    n_particle_global = _n_particle_global;
  }

  void start_timestep(uint32_t& n_receives_posted) {
    using Constants::proc_null;
    using Constants::count_tag;
    // Messages are sent up the tree whenever a rank has completed its local work
    // or received an updated particle complete count from its child
    // Messages are sent down the tree only after completion and starting at the 
    // root node. 
    // Post receives for photon counts from children and parent now
    if (child1!=proc_null) {
      MPI_Irecv(c1_recv_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child1, 
        count_tag, MPI_COMM_WORLD, &c1_recv_req);
      n_receives_posted++;
      c1_recv_buffer.set_awaiting();
    }
    if (child2!=proc_null) {
      MPI_Irecv(c2_recv_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child2,
        count_tag, MPI_COMM_WORLD, &c2_recv_req);
      n_receives_posted++;
      c2_recv_buffer.set_awaiting();
    }
    if (parent != proc_null) {
      MPI_Irecv(p_recv_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, parent,
        count_tag, MPI_COMM_WORLD, &p_recv_req);
      n_receives_posted++;
      p_recv_buffer.set_awaiting();
    }
  }



  //! Check for completed particle counts from children and parent.
  //!  Add children to current tree count. If parent count is received, 
  //!  it will be the global problem particle count, indicating completion
  void check_messages(uint64_t& n_complete_tree,
                      uint32_t& n_receives_posted, 
                      uint32_t& n_receives_completed, 
                      uint32_t& n_sends_completed) 
  { 
    using Constants::proc_null;
    using Constants::count_tag;

    //test receives from children and add work to tree count
    if (c1_recv_buffer.awaiting()) {
      MPI_Test(&c1_recv_req, &flag_c1, MPI_STATUS_IGNORE);
      if (flag_c1) {
        n_receives_completed++;
        c1_recv_buffer.set_received();
        n_complete_c1 = c1_recv_buffer.get_object()[0];
        //update tree count 
        n_complete_tree+=n_complete_c1;
        //post receive again
        c1_recv_buffer.reset();
        MPI_Irecv(c1_recv_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child1, 
          count_tag, MPI_COMM_WORLD, &c1_recv_req);
        n_receives_posted++;
        c1_recv_buffer.set_awaiting();
      }
    }
    if (c2_recv_buffer.awaiting()) {
      MPI_Test(&c2_recv_req, &flag_c2, MPI_STATUS_IGNORE);
      if (flag_c2) {
        n_receives_completed++;
        c2_recv_buffer.set_received();
        n_complete_c2 = c2_recv_buffer.get_object()[0];
        //update tree count 
        n_complete_tree+=n_complete_c2;
        //post receive again
        c2_recv_buffer.reset();
        MPI_Irecv(c2_recv_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child2, 
          count_tag, MPI_COMM_WORLD, &c2_recv_req);
        n_receives_posted++;
        c2_recv_buffer.set_awaiting();
      }
    }

    // test receive from parent (indicates completion of time step)
    if (p_recv_buffer.awaiting()) {
      MPI_Test(&p_recv_req, &flag_p, MPI_STATUS_IGNORE);
      if (flag_p) {
        n_receives_completed++;
        p_recv_buffer.set_received();
        n_complete_p = p_recv_buffer.get_object()[0];
      }
    }

    // test sends to parent
    if (p_send_buffer.sent() ) {
      MPI_Test(&p_send_req, &flag_p, MPI_STATUS_IGNORE);
      if (flag_p) {
        n_sends_completed++;
        p_send_buffer.reset();
      }
    }
  }

  void send_parent_n_tree_complete(uint64_t& n_complete_tree, 
                                    uint32_t& n_sends_posted)
  {
    using std::vector;
    using Constants::count_tag;
    using Constants::proc_null;

    if (parent!=proc_null) {
      p_send_buffer.fill(vector<uint64_t> (1,n_complete_tree));
      MPI_Isend(p_send_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, parent,
        count_tag, MPI_COMM_WORLD, &p_send_req);
      n_sends_posted++;
      //n_complete_messages++;
      p_send_buffer.set_sent();
      //reset tree count so work is not double counted
      n_complete_tree =0;
    }
  }

  bool is_finished(const uint64_t& n_complete_tree ) {
    finished  = (n_complete_tree == n_particle_global || 
      n_complete_p == n_particle_global);
    return finished;
  }

  void end_timestep(uint32_t& n_sends_posted,
                    uint32_t& n_sends_completed, 
                    uint32_t& n_receives_posted, 
                    uint32_t& n_receives_completed)
  {
    using std::vector;
    using Constants::proc_null;
    using Constants::count_tag;

    // finish off sends and send empty messages to complete awaiting receives
    //send finished count down tree to children and wait for completion
    if (child1 != proc_null) {
      if (c1_send_buffer.sent()) {
        MPI_Wait(&c1_send_req, MPI_STATUS_IGNORE);
        n_sends_completed++;
      }
      c1_send_buffer.fill(vector<uint64_t> (1,n_particle_global));
      MPI_Isend(c1_send_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child1, 
        count_tag, MPI_COMM_WORLD, &c1_send_req);
      n_sends_posted++;
      MPI_Wait(&c1_send_req, MPI_STATUS_IGNORE);
      n_sends_completed++;
    }
    if (child2 != proc_null)  {
      if (c2_send_buffer.sent()) {
        MPI_Wait(&c2_send_req, MPI_STATUS_IGNORE);
        n_sends_completed++;
      }
      c2_send_buffer.fill(vector<uint64_t> (1,n_particle_global));
      MPI_Isend(c2_send_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, child2,
        count_tag, MPI_COMM_WORLD, &c2_send_req);
      n_sends_posted++;
      MPI_Wait(&c2_send_req, MPI_STATUS_IGNORE);
      n_sends_completed++;
    }

    // wait for parent send to complete, if sent then finish off
    // parent's receive calls 
    if (parent!=proc_null) {
      if (p_send_buffer.sent()) {
        MPI_Wait(&p_send_req, MPI_STATUS_IGNORE);
        n_sends_completed++;
      }
      p_send_buffer.fill(vector<uint64_t> (1,1));
      MPI_Isend(p_send_buffer.get_buffer(), 1, MPI_UNSIGNED_LONG, parent,
        count_tag, MPI_COMM_WORLD, &p_send_req);
      n_sends_posted++;
      MPI_Wait(&p_send_req, MPI_STATUS_IGNORE);
      n_sends_completed++;
    }

    if (child1 != proc_null) {
      MPI_Wait(&c1_recv_req, MPI_STATUS_IGNORE);
      n_receives_completed++;
    }
    if (child2 != proc_null) {
      MPI_Wait(&c2_recv_req, MPI_STATUS_IGNORE);
      n_receives_completed++;
    }

    //reset tree counts
    n_particle_global = 0;
    n_complete_c1 = 0;
    n_complete_c2 = 0;
    n_complete_p = 0;
  }

  private:
  uint64_t n_complete_c1;
  uint64_t n_complete_c2; 
  uint64_t n_complete_p;
  uint64_t n_particle_global;
  bool finished;
  Buffer<uint64_t> c1_recv_buffer;
  Buffer<uint64_t> c2_recv_buffer;
  Buffer<uint64_t> p_recv_buffer;
  Buffer<uint64_t> c1_send_buffer;
  Buffer<uint64_t> c2_send_buffer;
  Buffer<uint64_t> p_send_buffer;
  bool c1_req_flag;
  bool c2_req_flag;
  bool p_req_flag;
  uint64_t child1;
  uint64_t child2;
  uint64_t parent;
  int flag_c1;
  int flag_c2;
  int flag_p;
  MPI_Request p_recv_req;
  MPI_Request c1_recv_req;
  MPI_Request c2_recv_req;
  MPI_Request p_send_req;
  MPI_Request c1_send_req;
  MPI_Request c2_send_req;
};


#endif // def completion_manager_h_
