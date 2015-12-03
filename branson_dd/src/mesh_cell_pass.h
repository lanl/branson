/*
  Author: Alex Long
  Date: 12/1/2015
  Name: mesh_cell_pass.h
*/

#ifndef mesh_cell_pass_h_
#define mesh_cell_pass_h_

#include <vector>
#include <string>
#include <algorithm>

#include "mesh.h"
#include "input.h"

class Mesh_Cell_Pass :: Mesh {

  public:
  Mesh_Cell_Pass(Input input, unsigned int _rank, unsigned int _nrank) 
    : Mesh(input, _rank, _nrank)
  {

    // JAYENNE code used to generate numbering in binary tree
    // compute parent and child node ids, ignoring missing nodes
    parent = (rank + 1) / 2 - 1;
    child1 = rank * 2 + 1;
    child2 = child1 + 1;

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
  }
  ~Mesh_Cell_Pass(void) {
    delete[] r_cell_reqs;
    delete[] r_cell_ids_reqs;
    delete[] s_cell_reqs;
    delete[] s_cell_ids_reqs;
  }

  bool on_processor(const unsigned int& index) const { 
    return  (index>=on_rank_start) && (index<=on_rank_end) ; 
  }

  void print_map(void) {
    for ( std::map<unsigned int,Cell>::iterator map_i =
      stored_cells.begin();
      map_i!=stored_cells.end(); map_i++)
      (map_i->second).print();
  }

  bool mesh_available(const unsigned int& index) const {
    if (on_processor(index)) return true;
    else if (stored_cells.find(index) != stored_cells.end())
      return true;
    else 
      return false; 
  } 

  void request_cell(const unsigned int& index) {
    //get local index of global index
    unsigned int off_rank_id = get_off_rank_id(index);
    unsigned int off_rank_local = index - off_rank_bounds[off_rank_id];
    //get correct index into received photon vector
    unsigned int r_index = off_rank_id - (off_rank_id>rank);
    //add to the request list if not already requested (use global index)
    if (ids_requested.find(index) == ids_requested.end()) {
      ids_requested.insert(index);
      ids_needed[r_index].push_back(off_rank_local);
      //need data from off_rank_id that has not already been requested
      need_data[r_index] = true;
      //increment number of RMA requests
      off_rank_reads++;
    }
  }



  bool process_mesh_requests(mpi::communicator world) {
    using Constants::cell_tag;
    using Constants::cell_id_tag;

    bool new_data = false;
    for (unsigned int ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        unsigned int r_index = ir - (ir>rank);
        
        ////////////////////////////////////////////////////////////////////////
        // if you need data from this rank, process request
        ////////////////////////////////////////////////////////////////////////
        if (need_data[r_index]) {
          //if you haven't requested cells, do that
          if (!b_s_cell_ids_reqs[r_index]) {
            //copy needed cells to the send buffer
            s_cell_ids[r_index].assign(
              ids_needed[r_index].begin(), ids_needed[r_index].end());
            //cout<<"Number of cells needed by "<<rank<<" from "<<ir;
            //cout<<" is " <<ids_needed[r_index].size()<<endl;
            //cout<<"Total IDS requested by rank "<<rank<<" is: ";
            //cout<<off_rank_reads<<endl;
            s_cell_ids_reqs[r_index].request( 
              world.isend(ir, cell_id_tag, s_cell_ids[r_index]));
            b_s_cell_ids_reqs[r_index] = true;
            // clear requested cell ids 
            // (requested cells will not be requested again)
            // because they are still stored in requested_ids set
            ids_needed[r_index].clear();
          }
          //otherwise, check to see if send completed, then reset buffer
          //post receive for this rank, if not done
          if (!b_r_cell_reqs[r_index]) {
            r_cell_reqs[r_index].request( 
              world.irecv(ir, cell_tag, r_cells[r_index]));
            b_r_cell_reqs[r_index] = true;
          }
          else {
            //check for completion
            if (r_cell_reqs[r_index].test()) {
              //reset send and receive flags and need_data flag 
              b_s_cell_ids_reqs[r_index] = false;
              b_r_cell_reqs[r_index] = false;
              if (ids_needed[r_index].empty()) need_data[r_index] = false;
              new_data = true;

              //add received cells to working mesh
              for (unsigned int i=0; i<r_cells[r_index].size();i++) {
                unsigned int index = r_cells[r_index][i].get_ID();
                //add this cell to the map, if possible, otherwise manage map
                if (stored_cells.size() < max_map_size) 
                  stored_cells[index] = r_cells[r_index][i];
                else {
                  //remove from map and from requests so it can be 
                  // reqeusted again if needed
                  unsigned int removed_id = (stored_cells.begin())->first ;
                  ids_requested.erase(removed_id);
                  stored_cells.erase(stored_cells.begin());
                  stored_cells[index] = r_cells[r_index][i];
                }
              } // for i in r_cells[r_index]
              r_cells[r_index].clear();
              s_cell_ids[r_index].clear();
            } // if r_cells_reqs[r_index].test()
          }
        } // if need_data[r_index]


        ////////////////////////////////////////////////////////////////////////
        // receiving cell ids needed by other ranks (post receives to all)
        ////////////////////////////////////////////////////////////////////////
        // check to see if receive call made
        if (!b_r_cell_ids_reqs[r_index]) {
          r_cell_ids_reqs[r_index].request(
            world.irecv( ir, cell_id_tag, r_cell_ids[r_index]));
          b_r_cell_ids_reqs[r_index] = true;
        }
        // add cell ids to requested for a rank
        else if (!send_data[r_index]) {
          if (r_cell_ids_reqs[r_index].test())  {
            //send data to this rank
            send_data[r_index] = true;
            //make cell send list for this rank
            for (unsigned int i=0; i<r_cell_ids[r_index].size();i++)
              s_cells[r_index].push_back(cells[r_cell_ids[r_index][i]]);
          } // if (r_cell_ids[r_index].test() )
        }
   
        //send cells needed by other ranks
        //check to see if this rank need your data
        if (send_data[r_index]) {
          if(!b_s_cell_reqs[r_index]) {
            //cout<<"Number of cells sent by "<<rank<<" to "<<ir<<" is " ;
            //cout<<s_cells[r_index].size()<<endl;
            s_cell_reqs[r_index].request(
              world.isend( ir, cell_tag, s_cells[r_index]));
            b_s_cell_reqs[r_index] = true;
          } 
          else {
            //check for completion of send message
            if (s_cell_reqs[r_index].test()) {
              //data is no longer needed by this rank, buffers can be reused
              // and receive and send messages can be posted again
              send_data[r_index] = false;
              b_r_cell_ids_reqs[r_index] = false;
              b_s_cell_reqs[r_index] = false;
              r_cell_ids[r_index].clear();
              s_cells[r_index].clear();
            }
          }
        } // if send_data[r_index]
      } //if (ir != rank)
    } //for ir in rank
    return new_data;
  }


  void purge_working_mesh(void) {
    stored_cells.clear(); ghost_cells.clear();
    ids_requested.clear();  
  }

  
  private:
  unsigned int max_map_size; //!< Maximum size of map object
  unsigned int off_rank_reads; //!< Number of off rank reads

  //send and receive buffers
  std::vector<std::vector<Cell> > r_cells; //!< Receive buffer for cells
  std::vector<std::vector<unsigned int> > r_cell_ids; //!< Receive cell ids needed by other ranks
  std::vector<std::vector<Cell> > s_cells; //!< Cells to send to each rank
  std::vector<std::vector<unsigned int> > s_cell_ids; //!< Send buffer for cell ids needed by this rank

  //Data bools needed from off rank and other ranks that need data here
  std::vector<bool> need_data; //!< Vector of size nrank-1, flag for needed data from rank
  std::vector<bool> send_data; //!< Vector of size nrank-1, flag to send data to rank

  // MPI requests for non-blocking communication and the 
  // bools for if the request has been made
  //receive requests and bool flags
  Request* r_cell_reqs; //!< Received cell requests
  std::vector<bool>  b_r_cell_reqs; //!< Bool for received call requests
  Request* r_cell_ids_reqs; //!< Received cell id requests
  std::vector<bool>  b_r_cell_ids_reqs; //!< Bool for received call requests
  //send requests and bool flags
  Request* s_cell_reqs; //!< Sent cell requests
  std::vector<bool>  b_s_cell_reqs; //!< Bool for received call requests
  Request* s_cell_ids_reqs; //!< Sent cell id requests
  std::vector<bool>  b_s_cell_ids_reqs; //!< Bool for received call requests

  std::map<unsigned int, Cell> stored_cells; //!< Cells that have been accessed off rank
  std::map<unsigned int, Cell> ghost_cells; //!< Static list of off-rank cells next to boundary

  std::vector<std::vector<unsigned int> > ids_needed ; //!< Cell needed by this rank
  std::set<unsigned int> ids_requested; //!< IDs that have been requested
};

#endif // mesh_cell_pass_h_
