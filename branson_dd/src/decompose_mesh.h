/* decompose_mesh.h
 * This file includes functions to decompose the 
 * mesh with Zoltan.
 * 6/17/2015 by Alex Long
*/

#ifndef decompose_mesh_h_
#define decompose_mesh_h_

#include <mpi.h>
#include <parmetis.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include "mesh.h"
#include "mpi_types.h"
#include "buffer.h"


void print_MPI_out(Mesh *mesh, uint32_t rank, uint32_t size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      mesh->post_decomp_print();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}


void print_MPI_maps(Mesh *mesh, uint32_t rank, uint32_t size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      mesh->print_map();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}

void decompose_mesh(Mesh* mesh, MPI_Types* mpi_types) {

  using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
  using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
  using std::vector;
  using std::map;
  using std::partial_sum;
  
  int rank, nrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  // make off processor map
  int n_off_rank = nrank -1;
  map<int,int> proc_map;
  for (uint32_t i=0; i<n_off_rank; i++) {
    int r_index = i + int(i>=rank);
    proc_map[i] = r_index;
  }

  MPI_Datatype MPI_Cell = mpi_types->get_cell_type();

  //begin PARMETIS routines

  //get the number of cells on each processor
  //vtxdist has number of vertices on each rank, same for all ranks
  vector<int> start_ncells(nrank, 0);
  vector<int> vtxdist(nrank, 0);

  int ncell_on_rank = mesh->get_n_local_cells();
  start_ncells[rank] = ncell_on_rank;
  MPI_Allreduce(MPI_IN_PLACE, 
                &start_ncells[0], 
                nrank, 
                MPI_INT, 
                MPI_SUM, 
                MPI_COMM_WORLD);
  partial_sum(start_ncells.begin(), 
              start_ncells.end(), 
              vtxdist.begin());
  vtxdist.insert(vtxdist.begin(), 0); 

  // build adjacency list needed for ParMetis call for each rank
  // also get coordinates of cell centers for geometry based partitioning
  float *xyz = new float[ncell_on_rank*3];
  vector<int> xadj;
  vector<int> adjncy;
  int adjncy_ctr = 0;
  Cell cell;
  uint32_t g_ID; //! Global ID
  for (uint32_t i=0; i<mesh->get_n_local_cells();i++) {
    cell = mesh->get_pre_renumber_cell(i);
    g_ID = cell.get_ID();
    cell.get_center(&xyz[i*3]);
    uint32_t xm_neighbor =cell.get_next_cell(X_NEG);
    uint32_t xp_neighbor =cell.get_next_cell(X_POS);
    uint32_t ym_neighbor =cell.get_next_cell(Y_NEG);
    uint32_t yp_neighbor =cell.get_next_cell(Y_POS);
    uint32_t zm_neighbor =cell.get_next_cell(Z_NEG);
    uint32_t zp_neighbor =cell.get_next_cell(Z_POS);
    
    xadj.push_back(adjncy_ctr); //starting index in xadj for this cell's nodes
    if (xm_neighbor != g_ID) {adjncy.push_back(xm_neighbor); adjncy_ctr++;}
    if (xp_neighbor != g_ID) {adjncy.push_back(xp_neighbor); adjncy_ctr++;}
    if (ym_neighbor != g_ID) {adjncy.push_back(ym_neighbor); adjncy_ctr++;}
    if (yp_neighbor != g_ID) {adjncy.push_back(yp_neighbor); adjncy_ctr++;}
    if (zm_neighbor != g_ID) {adjncy.push_back(zm_neighbor); adjncy_ctr++;}
    if (zp_neighbor != g_ID) {adjncy.push_back(zp_neighbor); adjncy_ctr++;}
  }
  xadj.push_back(adjncy_ctr);  


  int ndim = 3;
  int wgtflag = 0; //no weights for cells
  int numflag = 0; //C-style numbering
  int ncon = 1; 
  int nparts = nrank; //sub-domains = nrank

  float *tpwgts = new float[nparts];
  for (int i=0; i<nparts; i++) tpwgts[i]=1.0/float(nparts);

  float *ubvec = new float[ncon];
  for (int i=0; i<ncon; i++) ubvec[i]=1.05;

  int options[3];
  options[0] = 1; // 0--use default values, 1--use the values in 1 and 2
  options[1] = 3; //output level
  options[2] = 1242; //random number seed
  
  int edgecut =0; 
  int *part = new int[ncell_on_rank];
 
  ParMETIS_V3_PartGeomKway( &vtxdist[0],   // array describing how cells are distributed
                        &xadj[0],   // how cells are stored locally
                        &adjncy[0], // how cells are stored loccaly
                        NULL,       // weight of vertices
                        NULL,       // weight of edges
                        &wgtflag,   // 0 means no weights for node or edges
                        &numflag,   // numbering style, 0 for C-style
                        &ndim,      // n dimensions
                        xyz,        // coorindates of vertices
                        &ncon,      // weights per vertex
                        &nparts,    // number of sub-domains
                        tpwgts,     // weight per sub-domain
                        ubvec,      // unbalance in vertex weight
                        options,    // options array
                        &edgecut,   // OUTPUT: Number of edgecuts
                        part,       // OUTPUT: partition of each vertex
                        &comm); // MPI communicator

  //if edgecuts are made (edgecut > 0) send cells to other processors
  //otherwise mesh is already partitioned


  vector<int> recv_from_rank(n_off_rank,0);
  vector<int> send_to_rank(n_off_rank, 0);

  MPI_Request *reqs = new MPI_Request[n_off_rank*2];
  vector<Buffer<Cell> > send_cell(n_off_rank);
  vector<Buffer<Cell> > recv_cell(n_off_rank);

  if (edgecut) {
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      //sends
      int off_rank = proc_map[ir];
      // make list of cells to send to off_rank
      vector<Cell> send_list;
      for (uint32_t i=0; i<ncell_on_rank; i++) {
        if(part[i] == off_rank)
          send_list.push_back(mesh->get_pre_renumber_cell(i));
      }
      send_to_rank[ir] = send_list.size();
      send_cell[ir].fill(send_list);
      MPI_Isend(&send_to_rank[ir], 1,  MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD, &reqs[ir]);
      MPI_Irecv(&recv_from_rank[ir], 1, MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD, &reqs[ir+n_off_rank]);

      // erase sent cells from the mesh
      for (uint32_t i=0; i<ncell_on_rank; i++) {
        if(part[i] == off_rank) mesh->remove_cell(i);
      }
    }

    MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE); 

    //now send the buffers and post receives  
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      int off_rank = proc_map[ir];
      MPI_Isend(send_cell[ir].get_buffer(), 
                send_to_rank[ir],  
                MPI_Cell, 
                off_rank, 
                0, 
                MPI_COMM_WORLD, 
                &reqs[ir]);
      recv_cell[ir].resize(recv_from_rank[ir]);
      MPI_Irecv(recv_cell[ir].get_buffer(),
                recv_from_rank[ir], 
                MPI_Cell, 
                off_rank, 
                0, 
                MPI_COMM_WORLD, 
                &reqs[ir+n_off_rank]);
    }
 
    MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE); 

    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      vector<Cell> new_cells = recv_cell[ir].get_object();
      for (uint32_t i = 0; i< new_cells.size(); i++) {
        mesh->add_mesh_cell(new_cells[i]);
      }
    } 
  }

  //update the cell list on each processor
  mesh->update_mesh();

  // Gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh->get_n_local_cells();
  vector<uint32_t> out_cells_proc(nrank, 0);
  out_cells_proc[rank] = mesh->get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, 
                &out_cells_proc[0], 
                nrank, 
                MPI_INT, 
                MPI_SUM, 
                MPI_COMM_WORLD);

  // Prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(nrank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(), prefix_cells_proc.begin());

  // Set global numbering
  uint32_t g_start = prefix_cells_proc[rank]-n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank]-1;
  mesh->set_global_bound(g_start, g_end);

  // Append zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert( prefix_cells_proc.begin(), 0);
  mesh->set_off_rank_bounds(prefix_cells_proc);

  //make sure each index is remapped ONLY ONCE!
  vector< vector<bool> > remap_flag;
  for (uint32_t i=0; i<n_cell_post_decomp; i++) remap_flag.push_back(vector<bool> (6,false));

  //change global indices to match a simple number system for easy sorting,
  //this involves sending maps to each processor to get new indicies
  map<uint32_t, uint32_t> local_map = mesh->get_map();
  vector<uint32_t> packed_map(n_cell_post_decomp*2);
  vector<Buffer<uint32_t> > recv_packed_maps(n_off_rank);
  
  uint32_t i_packed = 0;
  for(map<uint32_t, uint32_t>::iterator map_i=local_map.begin(); 
    map_i!=local_map.end(); map_i++) {
    packed_map[i_packed] = map_i->first;
    i_packed++;
    packed_map[i_packed] = map_i->second;
    i_packed++;
  }

  // Send and receive packed maps for remapping boundaries
  for (uint32_t ir=0; ir<n_off_rank; ir++) {
    int off_rank = proc_map[ir];
    // Send your packed map
    MPI_Isend(&packed_map[0],
              n_cell_post_decomp*2,  
              MPI_UNSIGNED,
              off_rank, 
              0, 
              MPI_COMM_WORLD, 
              &reqs[ir]);
    // Receive other packed maps
    recv_packed_maps[ir].resize(out_cells_proc[off_rank]*2);
    MPI_Irecv(recv_packed_maps[ir].get_buffer(),
              out_cells_proc[off_rank]*2 , 
              MPI_UNSIGNED, 
              off_rank, 
              0, 
              MPI_COMM_WORLD, 
              &reqs[ir+n_off_rank]);
  }
  
  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE); 
  
  for (uint32_t ir=0; ir<n_off_rank; ir++) {
    vector<uint32_t> off_packed_map = recv_packed_maps[ir].get_object();
    map<uint32_t, uint32_t> off_map;
    for (uint32_t m=0; m<off_packed_map.size(); m++) {
      off_map[off_packed_map[m]] =off_packed_map[m+1];
      m++;
    }
    mesh->set_indices(off_map, remap_flag);
  }

  //now update the indices of local IDs
  mesh->set_local_indices(local_map);

  //reallocate mesh data in new MPI window
  mesh->make_MPI_window();

  // clean up dynamically allocated memory
  delete[] tpwgts;
  delete[] ubvec;
  delete[] part;
  delete[] reqs;
}

#endif // decompose_mesh_h
