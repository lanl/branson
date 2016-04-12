/*
  completion_manager_rma.h
  by Alex Long
  3/10/2016
*/

#ifndef binary_tree_rma_h_
#define binary_tree_rma_h_

#include <iostream>
#include <mpi.h>

#include "constants.h"

class Completion_Manager_RMA
{

  public:
  Completion_Manager_RMA(const int& rank, const int& n_rank)
    : n_complete_tree(0),
      n_complete_c1(0),
      n_complete_c2(0),
      n_complete_p(0),
      n_particle_global(0),
      buffer_c1(0),
      buffer_c2(0),
      buffer_p(0),
      c1_req_flag(false),
      c2_req_flag(false),
      p_req_flag(false)
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

    // Get the size of MPI_UNSIGNED_LONG
    int size_mpi_uint64;
    MPI_Type_size(MPI_UNSIGNED_LONG, &size_mpi_uint64);
    
    //MPI_Win_set_attr(completion_window, MPI_WIN_MODEL, MPI_WIN_UNIFIED);

    // Make the MPI window for the number of particles completed on this
    // sub-tree, which includes this ranks completed particles
    //MPI_Win_create(&n_complete_tree, size_mpi_uint64, size_mpi_uint64, MPI_INFO_NULL,
    //  MPI_COMM_WORLD, &completion_window);

    MPI_Win_allocate(size_mpi_uint64, size_mpi_uint64, MPI_INFO_NULL,
      MPI_COMM_WORLD, &n_complete_tree, &completion_window);

    int flag;
    MPI_Win_get_attr(completion_window, MPI_WIN_MODEL, &memory_model, &flag); 
  }
  ~Completion_Manager_RMA() {MPI_Win_free(&completion_window);}


  //const functions
  uint64_t get_n_complete_tree(void) const {return *n_complete_tree;}

  int get_mpi_window_memory_type(void) const {
    return *memory_model;
  }
  //non-const functions

  void set_timestep_global_particles(uint64_t _n_particle_global) {
    n_particle_global = _n_particle_global;
  }

  void end_timestep(void) {
    //reset tree counts
    *n_complete_tree = 0;
    n_complete_c1 = 0;
    n_complete_c2 = 0;
    n_complete_p = 0;
    buffer_c1 = 0;
    buffer_c2 = 0;
    buffer_p = 0;
    //wait on outstanding requests, parent has already completed
    if (c1_req_flag) {
      MPI_Wait(&req_c1, MPI_STATUS_IGNORE);
      c1_req_flag = false;
    }
    if (c2_req_flag) {
      MPI_Wait(&req_c2, MPI_STATUS_IGNORE);
      c2_req_flag = false;
    }
  }

  void start_access(void) {
    int assert =0;
    MPI_Win_lock_all(assert,completion_window);
  }

  void end_access(void) {
    MPI_Win_unlock_all(completion_window);
  }

  bool binary_tree_rma(const uint64_t& n_complete_rank) {
    using Constants::proc_null;
    // Return value
    bool finished = false;
    // Test for completion of non-blocking RMA requests
    // If child requests were completed, add to tree complete count and
    // make a new request

    // child 1
    if (child1!= proc_null) {
      if (c1_req_flag) {
        MPI_Test(&req_c1, &flag_c1, MPI_STATUS_IGNORE);
        if (flag_c1) {
          n_complete_c1 = buffer_c1;
          c1_req_flag = false;
        }
      }
      if (!c1_req_flag) {
        MPI_Rget(&buffer_c1, 1, MPI_UNSIGNED_LONG, child1, 0,
          1, MPI_UNSIGNED_LONG, completion_window, &req_c1);
        c1_req_flag=true;
      }
    }

    // child 2
    if (child2!= proc_null) {
      if (c2_req_flag) {
        MPI_Test(&req_c2, &flag_c2, MPI_STATUS_IGNORE);
        if (flag_c2) {
          n_complete_c2 = buffer_c2;
          c2_req_flag = false;
        }
      }
      if (!c2_req_flag) {
        MPI_Rget(&buffer_c2, 1, MPI_UNSIGNED_LONG, child2, 0,
          1, MPI_UNSIGNED_LONG, completion_window, &req_c2);
        c2_req_flag=true;
      }
    }

    // If parent is complete, test for overall completion
    if (parent != proc_null) {
      if (p_req_flag) { 
        MPI_Test(&req_p, &flag_p, MPI_STATUS_IGNORE);
        if (flag_p) {
          n_complete_p = buffer_p;
          if (n_complete_p == n_particle_global) {
            finished=true;
            *n_complete_tree = n_complete_p;
          }
          p_req_flag =false;
        }
      }
      if (!p_req_flag && !finished) {
        MPI_Rget(&buffer_p, 1, MPI_UNSIGNED_LONG, parent, 0,
          1, MPI_UNSIGNED_LONG, completion_window, &req_p);
        p_req_flag=true;
      }
    }
    else {
      if (*n_complete_tree == n_particle_global) finished=true;
    }
    // Update tree count
    if (!finished) *n_complete_tree = n_complete_rank + n_complete_c1 + n_complete_c2;
    //MPI_Win_sync(completion_window);
    return finished;
  }

  private:
  uint64_t *n_complete_tree;
  uint64_t n_complete_c1;
  uint64_t n_complete_c2;
  uint64_t n_complete_p;
  uint64_t n_particle_global;
  uint64_t buffer_c1;
  uint64_t buffer_c2;
  uint64_t buffer_p;
  bool c1_req_flag;
  bool c2_req_flag;
  bool p_req_flag;
  uint64_t child1;
  uint64_t child2;
  uint64_t parent;
  MPI_Win completion_window;
  int flag_c1;
  int flag_c2;
  int flag_p;
  MPI_Request req_c1;
  MPI_Request req_c2;
  MPI_Request req_p;
  int *memory_model;
};

#endif // def transport_particle_pass_h_
