//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_rma_manager.h
 * \author Alex Long
 * \date   March 22 2016
 * \brief  Manages requests for remote mesh data
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rma_manager_h_
#define rma_manager_h_

#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <unordered_set>
#include <vector>

#include "constants.h"
#include "mpi_types.h"



//==============================================================================
/*!
 * \class RMA_Manager
 * \brief Makes and processes requests for remote mesh data
 *
 * \example no test yet
 */
//==============================================================================

class RMA_Manager
{

  public:
  //! constructor
  RMA_Manager(const int& _rank, 
    const std::vector<uint32_t>& _rank_bounds,
    const uint32_t n_global_cell,
    MPI_Types * mpi_types,
    MPI_Win& _mesh_window)
  : rank(_rank),
    rank_bounds(_rank_bounds),
    MPI_Cell(mpi_types->get_cell_type()),
    mesh_window(_mesh_window)
  {
    using std::vector;
    n_max_requests = 1000;
    max_active_index = 0;
    count=0;
    requests = std::vector<MPI_Request> (n_max_requests);
    cells_buffer = vector<Cell> (n_max_requests);
    complete_indices = vector<int> (n_max_requests);

    n_requests_vec= std::vector<uint32_t> (n_global_cell,0);

    int flag;
    MPI_Win_get_attr(mesh_window, MPI_WIN_MODEL, &memory_model, &flag); 
  }

  //! destructor
  ~RMA_Manager() {}


  //////////////////////////////////////////////////////////////////////////////
  // const functions                                                          //
  //////////////////////////////////////////////////////////////////////////////

  //! Get the type of memory model used by the MPI implementation on this system
  int get_mpi_window_memory_type(void) const {return *memory_model;}

  //////////////////////////////////////////////////////////////////////////////
  // non-const functions                                                      //
  //////////////////////////////////////////////////////////////////////////////

  //! Make one-sided request for a remote mesh cell
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
      MPI_Rget(&cells_buffer[count], 1, MPI_Cell, off_rank_id, off_rank_local, 
        1, MPI_Cell, mesh_window, &requests[count]);
      n_requests++;
      mesh_requested.insert(g_index);
    }
  }

  //! Open the MPI window associated with the cell array
  void start_access(void) {
    int assert =0;
    MPI_Win_lock_all(assert,mesh_window);
  }

  //! Close access to the MPI window associated with the cell array
  void end_access(void) { MPI_Win_unlock_all(mesh_window);}

  //! Search rank bounds to get remote mesh owner's rank
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

  //! Check for completion of all active RMA requests and return new cells
  std::vector<Cell> process_rma_mesh_requests(uint32_t& tally_n_reqs_completed)
  {
    new_cells.clear();
    if (!index_in_use.empty()) {
      MPI_Testsome(max_active_index+1, &requests[0], &n_req_complete,
        &complete_indices[0], MPI_STATUSES_IGNORE);

      int comp_index;
      uint32_t g_index;
      //pull completed requests out and return them (they are then passed to 
      // the mesh)

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
        n_requests_vec[g_index]++;
      }
    }
    return new_cells;
  }

  //! End timestep and set the number of requests to zero
  void end_timestep(void) {
    max_active_index=0;
    for(uint32_t i=0; i<n_requests_vec.size();i++) n_requests_vec[i]=0;
  }

  //! Check to see if mesh has already been requestsed
  bool mesh_is_requested(uint32_t g_index) {
    return mesh_requested.find(g_index) != mesh_requested.end();
  }
 
  //! Check to see if any MPI requests are active
  bool no_active_requests(void) {return index_in_use.empty();}

  //! Return vector of requets counts (used only for plotting SILO file)
  std::vector<uint32_t> get_n_request_vec(void) {return n_requests_vec;}

  private:
  int rank; //! MPI rank
  std::vector<uint32_t> rank_bounds; //! Global cell ID bounds on each rank
  std::vector<uint32_t> n_requests_vec; //! Number of times a cell was requested 
  MPI_Datatype MPI_Cell; //! Custom MPI datatype for cells
  MPI_Win mesh_window; //! Shared memory window of cell objects
  uint32_t n_max_requests; //! Maximum concurrent MPI requests (parameter)
  uint32_t max_active_index; //! Last index that contains an active MPI request
  uint32_t count; //! Last used index in the MPI request array
  int n_req_complete; //! Number of completed requests after MPI_Testsome
  std::vector<Cell> new_cells; //! New cells after MPI_Testsome
  std::vector<MPI_Request> requests; //! Array of MPI requests
  std::vector<Cell> cells_buffer; //! Buffer for receiving cell data
  
  //! Returned from MPI_Testsome, indicates completed requests at index
  std::vector<int> complete_indices;

  std::unordered_set<int> index_in_use; //! Stores indices of active requests

  //! Stores global IDs of requested cells
  std::unordered_set<uint32_t> mesh_requested;

  int *memory_model; //! Memory model of MPI window (implementation dependent)
};

#endif // def rma_manager_h_

//---------------------------------------------------------------------------//
// end of mesh_rma_manager.h
//---------------------------------------------------------------------------//
