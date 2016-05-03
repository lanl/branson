/*
  rma_manager.h
  by Alex Long
  3/22/2016
*/

#ifndef rma_manager_h_
#define rma_manager_h_

#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <unordered_set>
#include <vector>

#include "constants.h"

class RMA_Manager
{

  public:
  RMA_Manager(const int& _rank, 
    const std::vector<uint32_t>& _rank_bounds,
    const uint32_t n_global_cell,
    MPI_Win& _mesh_window)
  : rank(_rank),
    rank_bounds(_rank_bounds),
    mesh_window(_mesh_window)
  {
    using std::vector;
    n_max_requests = 1000;
    max_active_index = 0;
    count=0;
    requests = std::vector<MPI_Request> (n_max_requests);
    cells_buffer = vector<Cell> (n_max_requests);
    complete_indices = vector<int> (n_max_requests);

    n_requests= std::vector<uint32_t> (n_global_cell,0);

    // Make MPI_Cell datatype
    // Three type entries in the class
    const int entry_count = 3 ; 
    // 7 uint32_t, 6 int, 13 double
    int array_of_block_length[4] = {10, 6, 14};
    // Displacements of each type in the cell
    MPI_Aint array_of_block_displace[3] = 
      {0, 10*sizeof(uint32_t),  10*sizeof(uint32_t)+6*sizeof(int)};
    //Type of each memory block
    MPI_Datatype array_of_types[3] = {MPI_UNSIGNED, MPI_INT, MPI_DOUBLE}; 

    MPI_Type_create_struct(entry_count, array_of_block_length, 
      array_of_block_displace, array_of_types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);

    MPI_Type_size(MPI_Cell, &mpi_cell_size);

    int flag;
    MPI_Win_get_attr(mesh_window, MPI_WIN_MODEL, &memory_model, &flag); 
  }
  ~RMA_Manager() {}

  int get_mpi_window_memory_type(void) const {
    return *memory_model;
  }

  //non-const functions
  void request_cell_rma(const uint32_t& g_index, uint32_t& n_requests) {
    if ( index_in_use.size() != n_max_requests &&
      mesh_requested.find(g_index) == mesh_requested.end()) 
    {
      //get local index of global index
      uint32_t off_rank_id = get_off_rank_id(g_index);
      uint32_t off_rank_local = g_index - rank_bounds[off_rank_id];
      //get correct index into received photon vector
      uint32_t r_index = off_rank_id - (off_rank_id>rank);

      //check to see if request at count is in use
      while(index_in_use.find(count) != index_in_use.end() ) {
        count++;
        if (count==n_max_requests) count=0;
      }
      max_active_index = std::max(count,max_active_index);
      //record this index as in use
      index_in_use.insert(count);
      //make request at index of count in both arrays
      MPI_Rget(&cells_buffer[count], 1, MPI_Cell, off_rank_id, off_rank_local, 1, MPI_Cell,
        mesh_window, &requests[count]);
      n_requests++;
      mesh_requested.insert(g_index);
    }
  }

  void start_access(void) {
    int assert =0;
    MPI_Win_lock_all(assert,mesh_window);
  }

  void end_access(void) {
    MPI_Win_unlock_all(mesh_window);
  }

  uint32_t get_off_rank_id(const uint32_t& g_index) const {
    //find rank of index
    bool found = false;
    uint32_t min_i = 0;
    uint32_t max_i = rank_bounds.size()-1;
    uint32_t s_i; //search index
    while(!found) {
      s_i =(max_i + min_i)/2;
      if (s_i == max_i || s_i == min_i) found = true;
      else if (g_index >= rank_bounds[s_i]) min_i = s_i;
      else max_i = s_i;
    }
    return s_i;
  }

  std::vector<Cell> process_rma_mesh_requests(uint32_t& tally_n_reqs_completed) {
    new_cells.clear();
    if (!index_in_use.empty()) {
      MPI_Testsome(max_active_index+1, &requests[0], &n_req_complete, &complete_indices[0],
        MPI_STATUSES_IGNORE);

      int comp_index;
      uint32_t g_index;
      //pull completed requests out and return them (they are then passed to 
      // the mesh)

      if (n_req_complete == MPI_UNDEFINED) std::cout<<"This is abad"<<std::endl;
      tally_n_reqs_completed +=n_req_complete;
      new_cells.resize(n_req_complete);
      for (int i = 0;i<n_req_complete;i++){
        comp_index = complete_indices[i];
        new_cells[i] = cells_buffer[comp_index];
        g_index =cells_buffer[comp_index].get_ID();
        // remove request index from index_in_use set
        index_in_use.erase(comp_index);
        // remove global cell ID from mseh_requesed set
        mesh_requested.erase(g_index);
        // increment request count (for plotting)
        n_requests[g_index]++;
      }
      if (new_cells.size() != n_req_complete) std::cout<<"This is bad"<<std::endl;
    }
    last_return_count = new_cells.size();
    return new_cells;
  }

  void end_timestep(void) {
    max_active_index=0;
    for(uint32_t i=0; i<n_requests.size();i++) n_requests[i]=0;
  }

  bool mesh_is_requested(uint32_t g_index) {
    return mesh_requested.find(g_index) != mesh_requested.end();
  }
  
  bool no_active_requests(void) {
    return index_in_use.empty();
  }

  std::vector<uint32_t> get_n_request_vec(void) {return n_requests;}

  private:
  int rank;
  std::vector<uint32_t> rank_bounds;
  std::vector<uint32_t> n_requests;
  MPI_Win mesh_window;
  uint32_t n_max_requests;
  uint32_t max_active_index;
  uint32_t count;
  uint32_t last_return_count;
  int n_req_complete;
  std::vector<Cell> new_cells;
  std::vector<MPI_Request> requests;
  std::vector<Cell> cells_buffer;
  std::vector<int> complete_indices;
  std::unordered_set<int> index_in_use;
  std::unordered_set<uint32_t> mesh_requested;
  MPI_Datatype MPI_Cell;
  int mpi_cell_size;
  int *memory_model;
};

#endif // def rma_manager_h_
