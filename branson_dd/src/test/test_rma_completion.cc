/*
  Author: Alex Long
  Date: 3/18/2016
  Name: test_rma_completion.cc
*/

#include <iostream>
#include <mpi.h>

#include "../constants.h"
#include "../binary_tree_rma.h"
#include "testing_functions.h"

using std::cout;
using std::endl;
using Constants::proc_null;

int main (int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  
  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  int nfail = 0;
  
  // set up place in binary tree
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
  
  // test completion routine
  {
    bool completion_routine_pass = true;

    uint64_t rank_particles = 10000;
    uint64_t global_count = n_rank*rank_particles;
    Completion comp(child1, child2, parent);

    //begin access epoch
    comp.start_access();

    //set global particle count
    comp.set_timestep_global_particles(global_count);

    uint64_t rank_complete = 0;
    double work;
    bool finished = false;

    while(!finished ) {
      if (rank_complete <rank_particles) rank_complete++;
      work+=exp(-rank_complete%8);
      finished = comp.binary_tree_rma(rank_complete);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //test to make sure all particles were tallied over all ranks
    if (comp.get_n_complete_tree() !=  global_count) completion_routine_pass = false;

    //reset for next timestep and test value again
    comp.end_timestep();
    if (comp.get_n_complete_tree() !=  0) completion_routine_pass = false;

    //end access epoch 
    comp.end_access();

    if (completion_routine_pass) cout<<"TEST PASSED: Completion routine with "
      <<n_rank<<" ranks"<<endl;
    else {
      cout<<"TEST FAILED: Completion routine with "<<n_rank<<" ranks"<<endl;
      nfail++;
    }
  }

  MPI_Finalize();

  return nfail;
}
