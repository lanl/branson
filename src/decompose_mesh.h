//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   decompose_mesh.h
 * \author Alex Long
 * \date   June 17 2015
 * \brief  Functions to decompose mesh with ParMetis
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef decompose_mesh_h_
#define decompose_mesh_h_

#include <mpi.h>
#include <parmetis.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <vector>

#include "mesh.h"
#include "mpi_types.h"
#include "buffer.h"
#include "timer.h"

//----------------------------------------------------------------------------//
//! Print the mesh information for each rank, one at a time
void print_MPI_out(Mesh *mesh, uint32_t rank, uint32_t size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank<size; ++p_rank) {
    if (rank == p_rank) {
      mesh->post_decomp_print();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//! Print the remapping information for each rank, one at a time
void print_MPI_maps(Mesh *mesh, uint32_t rank, uint32_t size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank<size; ++p_rank) {
    if (rank == p_rank) {
      mesh->print_map();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! partition a mesh with parmetis
std::vector<int> parmetis_partition(Mesh* mesh, int &edgecut, const int rank, 
  const int n_rank) {

  using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
  using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
  using std::vector;

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  // get the number of cells on each processor
  // vtxdist has number of vertices on each rank, same for all ranks
  vector<int> start_ncells(n_rank, 0);
  vector<int> vtxdist(n_rank, 0);

  uint32_t ncell_on_rank = mesh->get_n_local_cells();
  start_ncells[rank] = ncell_on_rank;

  MPI_Allreduce(MPI_IN_PLACE, &start_ncells[0], n_rank, MPI_INT, MPI_SUM,
    MPI_COMM_WORLD);
  partial_sum(start_ncells.begin(), start_ncells.end(), vtxdist.begin());
  vtxdist.insert(vtxdist.begin(), 0);

  // build adjacency list needed for ParMetis call for each rank also get
  // coordinates of cell centers for geometry based partitioning
  float *xyz = new float[ncell_on_rank*3];
  vector<int> xadj;
  vector<int> adjncy;
  int adjncy_ctr = 0;
  Cell cell;
  uint32_t g_ID; //! Global ID
  for (uint32_t i=0; i<ncell_on_rank;++i) {
    cell = mesh->get_pre_window_allocation_cell(i);
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

  int wgtflag = 0; // no weights for cells
  int numflag = 0; // C-style numbering
  int ncon = 1;
  int nparts = n_rank; // sub-domains = n_rank

  float *tpwgts = new float[nparts];
  for (int i=0; i<nparts; ++i) tpwgts[i]=1.0/float(nparts);

  float *ubvec = new float[ncon];
  for (int i=0; i<ncon; ++i) ubvec[i]=1.05;

  int options[3];
  options[0] = 1; // 0--use default values, 1--use the values in 1 and 2
  options[1] = 3; //output level
  options[2] = 1242; //random number seed

  std::vector<int> part(ncell_on_rank);

  ParMETIS_V3_PartKway( &vtxdist[0],   // array describing how cells are distributed
                        &xadj[0],   // how cells are stored locally
                        &adjncy[0], // how cells are stored loccaly
                        NULL,       // weight of vertices
                        NULL,       // weight of edges
                        &wgtflag,   // 0 means no weights for node or edges
                        &numflag,   // numbering style, 0 for C-style
                        &ncon,      // weights per vertex
                        &nparts,    // number of sub-domains
                        tpwgts,     // weight per sub-domain
                        ubvec,      // unbalance in vertex weight
                        options,    // options array
                        &edgecut,   // OUTPUT: Number of edgecuts
                        &part[0],   // OUTPUT: partition of each vertex
                        &comm); // MPI communicator

  delete[] ubvec;
  delete[] tpwgts;
  delete[] xyz;

  return part;
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
std::vector<int> cube_partition(Mesh* mesh, const int rank, const int n_rank) {

  uint32_t nx = mesh->get_global_n_x();
  uint32_t ny = mesh->get_global_n_y();
  uint32_t nz = mesh->get_global_n_z();

  // make sure this mesh is a cube
  int rank_side = std::floor(std::pow(n_rank, 1.0/3.0)) + 1;
  if (rank_side*rank_side*rank_side == n_rank) {
    if (rank == 0)
      std::cout<<"Cube decomposition OK, ranks per edge = "<<rank_side<<std::endl;
  }
  else if ( (rank_side-1)* (rank_side-1)*(rank_side-1) == n_rank) {
    rank_side--;
    if (rank ==0)
      std::cout<<"Cube decomposition OK, ranks per edge = "<<rank_side<<std::endl;
  }
  else {
    if (rank ==0) {
      std::cout<<"Can't use simple cube decomposition for this number of ranks";
      std::cout<<" Exiting..."<<std::endl;
    }
    exit(EXIT_FAILURE);
  }

  uint32_t ncell_on_rank = mesh->get_n_local_cells();
  std::vector<int> part(ncell_on_rank, rank);
  std::unordered_set<int> acceptor_ranks;

  for (uint32_t i=0; i<ncell_on_rank;++i) {
    const Cell &cell = mesh->get_pre_window_allocation_cell(i);
    // get the correct rank given the cell ID
    uint32_t index = cell.get_ID();
    uint32_t z = index / (nx*ny);
    uint32_t y = (index - z*(nx*ny))/nx;
    uint32_t x = index - z*(nx*ny) - y*nx;
    uint32_t rx = floor(double(x)/double(nx) * rank_side);
    uint32_t ry = floor(double(y)/double(ny) * rank_side);
    uint32_t rz = floor(double(z)/double(nz) * rank_side);
    part[i] = rz*rank_side*rank_side + ry*rank_side + rx;
    if (part[i] != rank) 
      acceptor_ranks.insert(part[i]);
  }
  return part;
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//! Send given partitioning scheme
void exchange_cells_post_partitioning(const int rank, MPI_Types* mpi_types, Mesh* mesh, const std::vector<int> &part) {
  using std::vector;
  MPI_Datatype MPI_Cell = mpi_types->get_cell_type();
  uint32_t ncell_on_rank = mesh->get_n_local_cells();

  std::unordered_set<int> acceptor_ranks;
  // use the partition to build the acceptor list
  for (uint32_t i=0;i<ncell_on_rank;++i) {
    if (part[i] != rank)
      acceptor_ranks.insert(part[i]);
  }

  // get the number of ranks that will send us a message
  const int one = 1;
  int *n_donors_win;
  int n_donors(0);

  // create the RMA memory windows for each value
  MPI_Win win;
  MPI_Request req;
  MPI_Win_allocate(1 * sizeof(int), 1 * sizeof(int), MPI_INFO_NULL,
    MPI_COMM_WORLD, &n_donors_win, &win);
  n_donors_win[0] = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  int assert = MPI_MODE_NOCHECK; // no conflicting locks on this window
  MPI_Win_lock_all(assert,win);
  for (auto ir : acceptor_ranks) {
    // increment the remote number of receives
    MPI_Raccumulate(&one, 1, MPI_INT, ir, 0, 1, MPI_INT, MPI_SUM, win, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }
  MPI_Win_unlock_all(win);
  MPI_Barrier(MPI_COMM_WORLD);
  n_donors = n_donors_win[0];
  MPI_Win_free(&win);

  int n_acceptors = acceptor_ranks.size();
  vector<Buffer<Cell> > send_cell(n_acceptors);
  vector<Buffer<Cell> > recv_cell(n_donors);
  vector<int> recv_from_rank(n_donors,0);
  vector<int> send_to_rank(n_acceptors, 0);
  vector<MPI_Request> reqs(n_donors + n_acceptors);
  vector<MPI_Status> status(n_donors + n_acceptors);
  std::unordered_map<int, int> donor_rank_size;

  // send sizes
  int isend = 0;
  for (auto& ir : acceptor_ranks) {
    // make list of cells to send to off_rank
    vector<Cell> send_list;
    for (uint32_t i=0; i<ncell_on_rank; ++i) {
      if(part[i] == ir)
        send_list.push_back(mesh->get_pre_window_allocation_cell(i));
    }
    send_to_rank[isend] = send_list.size();
    send_cell[isend].fill(send_list);

    MPI_Isend(&send_to_rank[isend], 1,  MPI_UNSIGNED, ir, 0,
      MPI_COMM_WORLD, &reqs[isend]);
    isend++;
  }

  // receive sizes, ranks
  for (uint32_t i = 0; i<n_donors;++i) {
    MPI_Irecv(&recv_from_rank[i], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
      MPI_COMM_WORLD, &reqs[n_acceptors + i]);
  }

  MPI_Waitall(reqs.size(), &reqs[0], &status[0]);
  MPI_Barrier(MPI_COMM_WORLD);

  // map donor rank to message size
  for (uint32_t i =0; i<n_donors;++i)
    donor_rank_size[status[n_acceptors + i].MPI_SOURCE] = recv_from_rank[i];

  // now send the buffers and post receives
  isend = 0;
  for (auto& ir : acceptor_ranks) {
    MPI_Isend(send_cell[isend].get_buffer(), send_to_rank[isend], MPI_Cell,
      ir, 0, MPI_COMM_WORLD, &reqs[isend]);
    isend++;
  }
  int ireceive = 0;
  for (auto& ir : donor_rank_size) {
    recv_cell[ireceive].resize(ir.second);
    MPI_Irecv(recv_cell[ireceive].get_buffer(), ir.second, MPI_Cell,
      ir.first, 0, MPI_COMM_WORLD, &reqs[n_acceptors+ireceive]);
    ireceive++;
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t i=0; i<n_donors; ++i) {
    vector<Cell> new_cells = recv_cell[i].get_object();
    for (uint32_t i = 0; i< new_cells.size(); ++i) {
      mesh->add_mesh_cell(new_cells[i]);
    }
  }
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Use metis to decompose on-rank mesh into chunks for better prefetching
void overdecompose_mesh(Mesh *mesh, const uint32_t grip_size) {
  using std::unordered_map;
  using std::vector;
  using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
  using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
  // get post-decomposition number of cells on this rank
  uint32_t n_cell_on_rank = mesh->get_n_local_cells();

  // make map of global IDs to local indices
  unordered_map<uint32_t, uint32_t> mesh_cell_ids;
  for (uint32_t i=0;i<n_cell_on_rank;++i) {
    mesh_cell_ids[mesh->get_pre_window_allocation_cell(i).get_ID()] = i;
  }
  unordered_map<uint32_t, uint32_t>::const_iterator end = mesh_cell_ids.end();

  vector<int> xadj;
  vector<int> adjncy;
  int adjncy_ctr = 0;
  Cell cell;
  uint32_t g_ID; // global ID
  for (uint32_t i=0; i<mesh->get_n_local_cells();++i) {
    cell = mesh->get_pre_window_allocation_cell(i);
    g_ID = cell.get_ID();
    uint32_t xm_neighbor =cell.get_next_cell(X_NEG);
    uint32_t xp_neighbor =cell.get_next_cell(X_POS);
    uint32_t ym_neighbor =cell.get_next_cell(Y_NEG);
    uint32_t yp_neighbor =cell.get_next_cell(Y_POS);
    uint32_t zm_neighbor =cell.get_next_cell(Z_NEG);
    uint32_t zp_neighbor =cell.get_next_cell(Z_POS);

     // starting index in xadj for this cell's nodes
    xadj.push_back(adjncy_ctr);
    // use local indices for the CSV to pass to METIS
    if (xm_neighbor != g_ID && mesh_cell_ids.find(xm_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[xm_neighbor]);
      adjncy_ctr++;
    }
    if (xp_neighbor != g_ID && mesh_cell_ids.find(xp_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[xp_neighbor]);
      adjncy_ctr++;
    }
    if (ym_neighbor != g_ID && mesh_cell_ids.find(ym_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[ym_neighbor]);
      adjncy_ctr++;
    }
    if (yp_neighbor != g_ID && mesh_cell_ids.find(yp_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[yp_neighbor]);
      adjncy_ctr++;
    }
    if (zm_neighbor != g_ID && mesh_cell_ids.find(zm_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[zm_neighbor]);
      adjncy_ctr++;
    }
    if (zp_neighbor != g_ID  && mesh_cell_ids.find(zp_neighbor) != end ) {
      adjncy.push_back(mesh_cell_ids[zp_neighbor]);
      adjncy_ctr++;
    }
  } // end loop over cells

  xadj.push_back(adjncy_ctr);

  int ncon = 1;

  // get the desired number of grips on a mesh, round down
  int n_grips = n_cell_on_rank/grip_size;
  if (!n_grips) n_grips=1;

  int rank_options[METIS_NOPTIONS];

  METIS_SetDefaultOptions(rank_options);
  rank_options[METIS_OPTION_NUMBERING] = 0; // C-style numbering

  int objval;
  int *grip_index = new int[n_cell_on_rank];
  int metis_return;

  // make a signed version of this for METIS call
  int signed_n_cell_on_rank = int(n_cell_on_rank);

  // Metis does not seem to like partitioning with one part, so skip in
  // that case
  if (n_grips > 1 && grip_size > 1) {

    metis_return =
      METIS_PartGraphKway( &signed_n_cell_on_rank, // number of on-rank vertices
                            &ncon,  // weights per vertex
                            &xadj[0], // how cells are stored locally
                            &adjncy[0], // how cells are stored loccaly
                            NULL, // weight of vertices
                            NULL, // size of vertices for comm volume
                            NULL, // weight of the edges
                            &n_grips,  // number of mesh grips on this mesh
                            NULL, // tpwgts (NULL = equal weight domains)
                            NULL, // unbalance in v-weight (NULL=1.001)
                            rank_options, // options array
                            &objval,  // OUTPUT: Number of edgecuts
                            grip_index);  // OUTPUT: grip ID of each vertex

    if (metis_return==METIS_ERROR_INPUT)
      std::cout<<"METIS: Input error"<<std::endl;
    else if (metis_return==METIS_ERROR_MEMORY)
      std::cout<<"METIS: Memory error"<<std::endl;
    else if (metis_return==METIS_ERROR)
      std::cout<<"METIS: Input error"<<std::endl;

    // set the grip index for each cell
    for (uint32_t i =0; i<uint32_t(signed_n_cell_on_rank);++i) {
      Cell& cell = mesh->get_pre_window_allocation_cell_ref(i);
      cell.set_grip_ID(grip_index[i]);
    }
    // sort cells by grip ID
    mesh->sort_cells_by_grip_ID();

    // delete dynamically allocated data
    delete[] grip_index;
  } // end if n_grips > 1

  else if (n_grips == 1) {
    // one grip on this rank, all cells have the same grip index
    // set the grip index for each cell
    for (uint32_t i =0; i<uint32_t(signed_n_cell_on_rank);++i) {
      Cell& cell = mesh->get_pre_window_allocation_cell_ref(i);
      cell.set_grip_ID(0);
    }
  }

  else if (grip_size ==1) {
    for (uint32_t i =0; i<uint32_t(signed_n_cell_on_rank);++i) {
      Cell& cell = mesh->get_pre_window_allocation_cell_ref(i);
      cell.set_grip_ID(i);
    }
    mesh->sort_cells_by_grip_ID();
  }
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//! Use one-sided communication to share new cell indices in mesh connectivity
void remap_cell_and_grip_indices_rma(Mesh *mesh, const int rank, const int n_rank) {
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh->get_n_local_cells();
  vector<uint32_t> out_cells_proc(n_rank, 0);
  out_cells_proc[rank] = mesh->get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, &out_cells_proc[0], n_rank, MPI_UNSIGNED,
    MPI_SUM, MPI_COMM_WORLD);

  // prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(n_rank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(), prefix_cells_proc.begin());

  // set global numbering
  uint32_t g_start = prefix_cells_proc[rank]-n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank]-1;
  mesh->set_global_bound(g_start, g_end);

  // set the grip ID to be the cells at the center of grips using global cell
  // indices
  mesh->set_grip_ID_using_cell_index();

  // get the maximum grip size for correct parallel operations
  uint32_t max_grip_size = mesh->get_max_grip_size();
  uint32_t global_max_grip_size=0;
  uint32_t global_min_grip_size=10000000;
  MPI_Allreduce(&max_grip_size, &global_max_grip_size, 1, MPI_UNSIGNED,
    MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_grip_size, &global_min_grip_size, 1, MPI_UNSIGNED,
    MPI_MIN, MPI_COMM_WORLD);
  mesh->set_max_grip_size(global_max_grip_size);

  if (rank==0) {
    std::cout<<"Minimum/Maximum grip size: "<<global_min_grip_size;
    std::cout<<" / "<<global_max_grip_size<<std::endl;
  }

  // prepend zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert(prefix_cells_proc.begin(), 0);
  mesh->set_off_rank_bounds(prefix_cells_proc); 

  // change global indices to match a simple number system for easy sorting,
  // this involves sending maps to each processor to get new indicies
  unordered_map<uint32_t, uint32_t> local_boundary_map = mesh->get_boundary_nodes();
  unordered_map<uint32_t, uint32_t> local_boundary_grip_map = mesh->get_boundary_grips();
  unordered_map<uint32_t, uint32_t> local_map = mesh->get_new_global_index_map();
  unordered_map<uint32_t, uint32_t> local_grip_map = mesh->get_grip_map();
  unordered_set<uint32_t> boundary_indices = mesh->get_boundary_neighbors();
  unordered_map<uint32_t,uint32_t> boundary_map;
  unordered_map<uint32_t,uint32_t> boundary_grip_map;

  // create the RMA memory windows to allow all ranks to write the new index
  // of a cell at window[old_index]
  MPI_Win index_win, grip_win;
  uint32_t *new_index;
  uint32_t *new_grip;
  MPI_Win_allocate(n_cell_post_decomp*sizeof(uint32_t), sizeof(uint32_t), MPI_INFO_NULL,
    MPI_COMM_WORLD, &new_index, &index_win);
  MPI_Win_allocate(n_cell_post_decomp*sizeof(uint32_t), sizeof(uint32_t), MPI_INFO_NULL,
    MPI_COMM_WORLD, &new_grip, &grip_win);
  MPI_Barrier(MPI_COMM_WORLD);
  //assert = MPI_MODE_NOCHECK; // no conflicting locks on this window
  int assert = 0;
  MPI_Win_lock_all(assert,index_win);
  MPI_Win_lock_all(assert,grip_win);

  uint32_t remapped_index, remapped_grip, local_index;
  int target_rank;
  for (auto& imap : local_boundary_map) {
    remapped_index = imap.second;
    target_rank = mesh->get_rank(imap.first);
    local_index = imap.first - prefix_cells_proc[target_rank];
    MPI_Put(&remapped_index, 1, MPI_UNSIGNED, target_rank, local_index, 1, 
      MPI_UNSIGNED, index_win);
  }

  for (auto& imap : local_boundary_grip_map) {
    remapped_grip = imap.second;
    target_rank = mesh->get_rank(imap.first);
    local_index = imap.first - prefix_cells_proc[target_rank];
    MPI_Put(&remapped_grip, 1, MPI_UNSIGNED, target_rank, local_index, 1, 
      MPI_UNSIGNED, grip_win);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  // make the memory visible to all ranks
  MPI_Win_flush_all(index_win);
  MPI_Win_flush_all(grip_win);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_sync(index_win);
  MPI_Win_sync(grip_win);
  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<uint32_t> new_boundary_indices(boundary_indices.size());
  std::vector<uint32_t> new_boundary_grips(boundary_indices.size());
  std::vector<MPI_Request> i_reqs(boundary_indices.size());
  std::vector<MPI_Request> g_reqs(boundary_indices.size());

  int ib = 0;
  for (auto& iset : boundary_indices) {
    target_rank = mesh->get_rank(iset);
    local_index = iset - prefix_cells_proc[target_rank];
    MPI_Rget(&new_boundary_indices[ib], 1, MPI_UNSIGNED, target_rank, local_index, 1, MPI_UNSIGNED, index_win,&i_reqs[ib]);
    MPI_Rget(&new_boundary_grips[ib], 1, MPI_UNSIGNED, target_rank, local_index, 1, MPI_UNSIGNED, grip_win, &g_reqs[ib]);
    ib++;
  }

  MPI_Waitall(i_reqs.size(), &i_reqs[0], MPI_STATUS_IGNORE);
  MPI_Waitall(g_reqs.size(), &g_reqs[0], MPI_STATUS_IGNORE);
  
  ib=0;
  for (auto& iset : boundary_indices) {
    boundary_map[iset] = new_boundary_indices[ib];
    boundary_grip_map[iset] = new_boundary_grips[ib];
    ib++;
  }

  MPI_Win_unlock_all(index_win);
  MPI_Win_unlock_all(grip_win);
  MPI_Win_free(&index_win);
  MPI_Win_free(&grip_win);


  // combine local index map with boundary map from communication
  local_map.insert(boundary_map.begin(), boundary_map.end());
  local_grip_map.insert(boundary_grip_map.begin(), boundary_grip_map.end());

  // now update the indices of local IDs
  mesh->renumber_local_cell_indices(local_map, local_grip_map);
}
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//! Use two-sided communication to share new cell indices in mesh connectivity
void remap_cell_and_grip_indices(Mesh *mesh, const int rank, const int n_rank) {
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;

  // make off processor map
  uint32_t n_off_rank = n_rank -1; // implicit conversion from int to uint32_t
  unordered_map<int,int> proc_map;
  for (uint32_t i=0; i<n_off_rank; ++i) {
    int r_index = i + int(int(i)>=rank);
    proc_map[i] = r_index;
  }

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh->get_n_local_cells();
  vector<uint32_t> out_cells_proc(n_rank, 0);
  out_cells_proc[rank] = mesh->get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, &out_cells_proc[0], n_rank, MPI_UNSIGNED,
    MPI_SUM, MPI_COMM_WORLD);

  // prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(n_rank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(), prefix_cells_proc.begin());

  // set global numbering
  uint32_t g_start = prefix_cells_proc[rank]-n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank]-1;
  mesh->set_global_bound(g_start, g_end);

  // set the grip ID to be the cells at the center of grips using global cell
  // indices
  mesh->set_grip_ID_using_cell_index();

  // get the maximum grip size for correct parallel operations
  uint32_t max_grip_size = mesh->get_max_grip_size();
  uint32_t global_max_grip_size=0;
  uint32_t global_min_grip_size=10000000;
  MPI_Allreduce(&max_grip_size, &global_max_grip_size, 1, MPI_UNSIGNED,
    MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_grip_size, &global_min_grip_size, 1, MPI_UNSIGNED,
    MPI_MIN, MPI_COMM_WORLD);
  mesh->set_max_grip_size(global_max_grip_size);

  if (rank==0) {
    std::cout<<"Minimum/Maximum grip size: "<<global_min_grip_size;
    std::cout<<" / "<<global_max_grip_size<<std::endl;
  }

  // prepend zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert(prefix_cells_proc.begin(), 0);
  mesh->set_off_rank_bounds(prefix_cells_proc); 
  // change global indices to match a simple number system for easy sorting,
  // this involves sending maps to each processor to get new indicies
  unordered_map<uint32_t, uint32_t> local_boundary_map = mesh->get_boundary_nodes();
  unordered_map<uint32_t, uint32_t> local_boundary_grip_map = mesh->get_boundary_grips();
  
  int n_boundary = local_boundary_map.size();
  // get the size of boundary cells on each rank
  std::vector<int> boundary_cell_size(n_rank, 0);
  boundary_cell_size[rank] = n_boundary;
  MPI_Allreduce(MPI_IN_PLACE, &boundary_cell_size[0], n_rank, MPI_INT, MPI_SUM, 
    MPI_COMM_WORLD);
  
  
  vector<uint32_t> packed_map(n_boundary*2);
  vector<uint32_t> packed_grip_map(n_boundary*2);
  vector<Buffer<uint32_t>> recv_packed_maps(n_off_rank);
  vector<Buffer<uint32_t>> recv_packed_grip_maps(n_off_rank);
                                                                                 
  MPI_Request *reqs = new MPI_Request[n_off_rank*2];                        
  MPI_Request *grip_reqs = new MPI_Request[n_off_rank*2];                        
                                                                                 
  // pack up the index map                                                       
  uint32_t i_packed = 0;                                                         
  for(auto map_i : local_boundary_map) {                                                  
    packed_map[i_packed] = map_i.first;                                          
    i_packed++;                                                                  
    packed_map[i_packed] = map_i.second;                                         
    i_packed++;                                                                  
  }                                                                              
                                                                                 
  // pack up the grip map                                                        
  uint32_t i_grip_packed = 0;                                                    
  for(auto map_i : local_boundary_grip_map) {                                             
    packed_grip_map[i_grip_packed] = map_i.first;                                
    i_grip_packed++;                                                             
    packed_grip_map[i_grip_packed] = map_i.second;                               
    i_grip_packed++;                                                             
  }                                                                              
                                                                                 
  // Send and receive packed maps for remapping boundaries                       
  for (uint32_t ir=0; ir<n_off_rank; ++ir) {                                     
    int off_rank = proc_map[ir];                                                 
    // Send your packed index map                                                
    MPI_Isend(&packed_map[0], n_boundary*2, MPI_UNSIGNED, off_rank,         
      0, MPI_COMM_WORLD, &reqs[ir]);                                             
                                                                                 
    // Send your packed grip map                                                 
    MPI_Isend(&packed_grip_map[0], n_boundary*2, MPI_UNSIGNED, off_rank, 
      1, MPI_COMM_WORLD, &grip_reqs[ir]);                                        
                                                                                 
    // Receive other packed index maps                                           
    recv_packed_maps[ir].resize(boundary_cell_size[off_rank]*2);                     
    MPI_Irecv(recv_packed_maps[ir].get_buffer(), boundary_cell_size[off_rank]*2,        
      MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD, &reqs[ir+n_off_rank]);          
                                                                                 
    // Receive other packed grip maps                                            
    recv_packed_grip_maps[ir].resize(boundary_cell_size[off_rank]*2);
    MPI_Irecv(recv_packed_grip_maps[ir].get_buffer(),                            
      boundary_cell_size[off_rank]*2, MPI_UNSIGNED, off_rank, 1, MPI_COMM_WORLD,        
      &grip_reqs[ir+n_off_rank]);                                                
  }
                                                                                 
  // wait for completion of all sends/receives                                   
  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE);                            
  MPI_Waitall(n_off_rank*2, grip_reqs, MPI_STATUS_IGNORE);                       
  
  // we just need to build a map for the non-local indices, everything else can
  // be remapped locally
  std::unordered_set<uint32_t> boundary_indices = mesh->get_boundary_neighbors();
  std::unordered_map<uint32_t,uint32_t> boundary_map;
  std::unordered_map<uint32_t,uint32_t> boundary_grip_map;
                                                                                 
  // remake the map objects from the packed maps and take boundary data from
  // them
  for (uint32_t ir=0; ir<n_off_rank; ++ir) {                                     
    vector<uint32_t> off_packed_map = recv_packed_maps[ir].get_object();         
    vector<uint32_t> off_packed_grip_map =                                       
      recv_packed_grip_maps[ir].get_object();                                    
    unordered_map<uint32_t, uint32_t> off_map;                                   
    unordered_map<uint32_t, uint32_t> off_grip_map;                              

    for (uint32_t m=0; m<off_packed_map.size(); m+=2) {                          
      off_map[off_packed_map[m]] =off_packed_map[m+1];
      off_grip_map[off_packed_grip_map[m]] =off_packed_grip_map[m+1];
    }
    for (auto& iset : boundary_indices) {
      if (off_map.find(iset) != off_map.end())
        boundary_map[iset] = off_map[iset];
      if (off_grip_map.find(iset) != off_grip_map.end())
        boundary_grip_map[iset] = off_grip_map[iset];
    }
  }                                                                              

  // combine local index map with boundary map from communication
  unordered_map<uint32_t, uint32_t> local_map = mesh->get_new_global_index_map();
  unordered_map<uint32_t, uint32_t> local_grip_map = mesh->get_grip_map();
  local_map.insert(boundary_map.begin(), boundary_map.end());
  local_grip_map.insert(boundary_grip_map.begin(), boundary_grip_map.end());

  // now update the indices of local IDs
  mesh->renumber_local_cell_indices(local_map, local_grip_map);

  // clean up dynamically allocated memory
  delete[] grip_reqs;
  delete[] reqs;
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//! Use two-sided communication to share new cell indices in mesh connectivity
void remap_cell_and_grip_indices_staged(Mesh *mesh, const int rank, const int n_rank) {
  using std::vector;
  using std::unordered_map;

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh->get_n_local_cells();
  vector<uint32_t> out_cells_proc(n_rank, 0);
  out_cells_proc[rank] = mesh->get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, &out_cells_proc[0], n_rank, MPI_UNSIGNED,
    MPI_SUM, MPI_COMM_WORLD);

  // prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(n_rank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(), prefix_cells_proc.begin());

  // set global numbering
  uint32_t g_start = prefix_cells_proc[rank]-n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank]-1;
  mesh->set_global_bound(g_start, g_end);

  // set the grip ID to be the cells at the center of grips using global cell
  // indices
  mesh->set_grip_ID_using_cell_index();

  // get the maximum grip size for correct parallel operations
  uint32_t max_grip_size = mesh->get_max_grip_size();
  uint32_t global_max_grip_size=0;
  uint32_t global_min_grip_size=10000000;
  MPI_Allreduce(&max_grip_size, &global_max_grip_size, 1, MPI_UNSIGNED,
    MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_grip_size, &global_min_grip_size, 1, MPI_UNSIGNED,
    MPI_MIN, MPI_COMM_WORLD);
  mesh->set_max_grip_size(global_max_grip_size);

  if (rank==0) {
    std::cout<<"Minimum/Maximum grip size: "<<global_min_grip_size;
    std::cout<<" / "<<global_max_grip_size<<std::endl;
  }

  // prepend zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert(prefix_cells_proc.begin(), 0);
  mesh->set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting,
  // this involves sending maps to each processor to get new indicies
  unordered_map<uint32_t, uint32_t> local_boundary_map = mesh->get_boundary_nodes();
  unordered_map<uint32_t, uint32_t> local_boundary_grip_map = mesh->get_boundary_grips();
 
  MPI_Barrier(MPI_COMM_WORLD);
 
  int n_boundary = local_boundary_map.size();
  // get the size of boundary cells on each rank
  std::vector<int> boundary_cell_size(n_rank, 0);
  boundary_cell_size[rank] = n_boundary;
  MPI_Allreduce(MPI_IN_PLACE, &boundary_cell_size[0], n_rank, MPI_INT, MPI_SUM, 
    MPI_COMM_WORLD);
  
  vector<uint32_t> packed_map(n_boundary*2);
  vector<uint32_t> packed_grip_map(n_boundary*2);

  // pack up the index map
  uint32_t i_packed = 0;
  for(auto map_i : local_boundary_map) {
    packed_map[i_packed] = map_i.first;
    i_packed++;
    packed_map[i_packed] = map_i.second;
    i_packed++;
  }

  // pack up the grip map
  uint32_t i_grip_packed = 0;
  for(auto map_i : local_boundary_grip_map) {
    packed_grip_map[i_grip_packed] = map_i.first;
    i_grip_packed++;
    packed_grip_map[i_grip_packed] = map_i.second;
    i_grip_packed++;
  }

  // we just need to build a map for the non-local indices, everything else can
  // be remapped locally
  std::unordered_set<uint32_t> boundary_indices = mesh->get_boundary_neighbors();

  // these store the needed boundary indices new map
  std::unordered_map<uint32_t,uint32_t> boundary_map;
  std::unordered_map<uint32_t,uint32_t> boundary_grip_map;

  // ok, need to communicate with a dumb, memory preserving algorithm
  // 1/n steps at a time, receive remaps from a fraction of the total ranks 
  for (uint32_t n=0; n<8;++n) {
    int send_ranks_start = n*(n_rank/8);
    int send_ranks_end = (n+1)*(n_rank/8);
    if (n==7)
      send_ranks_end == n_rank;
    int n_senders = send_ranks_end - send_ranks_start;
    
    std::vector<MPI_Request> id_reqs(n_senders);
    std::vector<MPI_Request> grip_reqs(n_senders);
    if (rank >= send_ranks_start && rank < send_ranks_end) {
      id_reqs.resize(n_rank + n_senders );
      grip_reqs.resize(n_rank + n_senders);
    }

    vector<Buffer<uint32_t>> recv_packed_maps(n_senders);
    vector<Buffer<uint32_t>> recv_packed_grip_maps(n_senders);

    // ranks in the send block send to all
    if (rank >= send_ranks_start && rank < send_ranks_end) {
      for (uint32_t i=0; i<n_rank;++i) {
          // Send your packed index map
          MPI_Isend(&packed_map[0], n_boundary*2, MPI_UNSIGNED, i,
            0, MPI_COMM_WORLD, &id_reqs[n_senders+i]);

          // Send your packed grip map
          MPI_Isend(&packed_grip_map[0], n_boundary*2, MPI_UNSIGNED, i,
            1, MPI_COMM_WORLD, &grip_reqs[n_senders+i]);
      }
    }

    // receive block (for all ranks)
    for (uint32_t i=0; i<n_senders;++i) { 
      int message_size = boundary_cell_size[send_ranks_start+i]*2;
      int source_rank = send_ranks_start+i;
      // Receive other packed index maps
      recv_packed_maps[i].resize(message_size);
      MPI_Irecv(recv_packed_maps[i].get_buffer(), message_size,
        MPI_UNSIGNED, source_rank, 0, MPI_COMM_WORLD, &id_reqs[i]);

      // Receive other packed grip maps
      recv_packed_grip_maps[i].resize(message_size);
      MPI_Irecv(recv_packed_grip_maps[i].get_buffer(), message_size,
        MPI_UNSIGNED, source_rank, 1, MPI_COMM_WORLD, &grip_reqs[i]);
    }

    // wait for completion of all sends and receives
    MPI_Waitall(id_reqs.size(), &id_reqs[0], MPI_STATUS_IGNORE);
    MPI_Waitall(grip_reqs.size(), &grip_reqs[0], MPI_STATUS_IGNORE);

    // remake the map objects from the packed maps and take boundary data from
    // them
    for (uint32_t i=0; i<n_senders; ++i) {
      if (send_ranks_start + i == rank)
        continue;
      vector<uint32_t> off_packed_map = recv_packed_maps[i].get_object();
      vector<uint32_t> off_packed_grip_map =
        recv_packed_grip_maps[i].get_object();
      unordered_map<uint32_t, uint32_t> off_map;
      unordered_map<uint32_t, uint32_t> off_grip_map;

      for (uint32_t m=0; m<off_packed_map.size(); m+=2) {                          
        off_map[off_packed_map[m]] =off_packed_map[m+1];
        off_grip_map[off_packed_grip_map[m]] =off_packed_grip_map[m+1];
      }
      for (auto& iset : boundary_indices) {
        if (off_map.find(iset) != off_map.end())
          boundary_map[iset] = off_map[iset];
        if (off_grip_map.find(iset) != off_grip_map.end())
          boundary_grip_map[iset] = off_grip_map[iset];
      }
    }
  } // end loop over send/receive steps (8)

  // combine local index map with boundary map from communication
  unordered_map<uint32_t, uint32_t> local_map = mesh->get_new_global_index_map();
  unordered_map<uint32_t, uint32_t> local_grip_map = mesh->get_grip_map();
  local_map.insert(boundary_map.begin(), boundary_map.end());
  local_grip_map.insert(boundary_grip_map.begin(), boundary_grip_map.end());

  // now update the indices of local IDs
  mesh->renumber_local_cell_indices(local_map, local_grip_map);

}
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//! Generate new partitioning with ParMetis, send and receive cells, renumber
// mesh and communicate renumbering
void decompose_mesh(Mesh* mesh, MPI_Types* mpi_types, const Info& mpi_info,
  const uint32_t& grip_size, const int decomposition_type)
{
  using std::unordered_map;
  using std::unordered_set;
  using Constants::CUBE;
  using Constants::PARMETIS;
  Timer t_partition;
  Timer t_overdecomp;
  Timer t_remap;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // parmetis sets this, if it's zero no changes were made to the default
  // partition (cube decomposition always changes paritioning)
  int edgecut = 0;
  if (rank==0) std::cout<<"partitioning..."<<std::endl;
  t_partition.start_timer("partition");
  // decomposition methods return a partition vector which is the rank of each
  // cell
  std::vector<int> part;
  if (decomposition_type == CUBE) {
    part = cube_partition(mesh, rank, n_rank);
    edgecut = 1;
  }
  else if (decomposition_type == PARMETIS)
    part = parmetis_partition(mesh, edgecut, rank, n_rank);
  else {
    if (rank ==0) {
      std::cout<<"Decomposition type not recognized.";
      std::cout<<" Exiting..."<<std::endl;
    }
    exit(EXIT_FAILURE);
  }

  // if edgecuts are made (edgecut > 0) send cells to other processors
  // otherwise mesh is already partitioned (sets "new_cells" in mesh object)
  if (edgecut)
    exchange_cells_post_partitioning(rank, mpi_types, mesh, part);
  t_partition.stop_timer("partition");
  
  // update the cell list on each processor
  mesh->set_post_decomposition_mesh_cells(part);

  if (rank==0) std::cout<<"overdecomposing..."<<std::endl;
  t_overdecomp.start_timer("overdecomp");
  // if using grips of cell data, overdecompose mesh on rank to smaller chunks
  overdecompose_mesh(mesh, grip_size);
  t_overdecomp.stop_timer("overdecomp");


  if (rank==0) std::cout<<"remapping..."<<std::endl;
  t_remap.start_timer("remap");
  remap_cell_and_grip_indices_staged(mesh, rank, n_rank);
  t_remap.stop_timer("remap");

  // reallocate mesh data in new MPI window and delete the old vector object
  mesh->make_MPI_window();

  if (rank==0) {
    std::cout<<"Partition: "<<t_partition.get_time("partition")<<std::endl;
    std::cout<<"Overdecomp: "<<t_overdecomp.get_time("overdecomp")<<std::endl;
    std::cout<<"Remap: "<<t_remap.get_time("remap")<<std::endl;
  }
}
//----------------------------------------------------------------------------//

//! Create replicated mesh by giving all cells to all other processors, renumber
// mesh and communicate renumbering
void replicate_mesh(Mesh* mesh, MPI_Types* mpi_types, const Info& mpi_info,
  const uint32_t& grip_size)
{
  using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
  using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
  using std::vector;
  using std::unordered_map;
  using std::unordered_set;
  using std::partial_sum;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // make off processor map
  uint32_t n_off_rank = n_rank -1; // implicit conversion from int to uint32_t
  unordered_map<int,int> proc_map;
  for (uint32_t i=0; i<n_off_rank; ++i) {
    int r_index = i + int(int(i)>=rank);
    proc_map[i] = r_index;
  }

  MPI_Datatype MPI_Cell = mpi_types->get_cell_type();

  // get the number of cells on each processor
  vector<int> start_ncells(n_rank, 0);

  uint32_t ncell_on_rank = mesh->get_n_local_cells();
  start_ncells[rank] = ncell_on_rank;

  MPI_Allreduce(MPI_IN_PLACE, &start_ncells[0], n_rank, MPI_INT, MPI_SUM,
    MPI_COMM_WORLD);

  vector<int> recv_from_rank(n_off_rank,0);
  vector<int> send_to_rank(n_off_rank, 0);

  MPI_Request *reqs = new MPI_Request[n_off_rank*2];
  vector<Buffer<Cell> > send_cell(n_off_rank);
  vector<Buffer<Cell> > recv_cell(n_off_rank);

  for (uint32_t ir=0; ir<n_off_rank; ++ir) {
    //sends
    int off_rank = proc_map[ir];
    // make list of cells to send to off_rank
    vector<Cell> send_list;
    for (uint32_t i=0; i<ncell_on_rank; ++i)
      send_list.push_back(mesh->get_pre_window_allocation_cell(i));
    send_to_rank[ir] = send_list.size();
    send_cell[ir].fill(send_list);

    MPI_Isend(&send_to_rank[ir], 1,  MPI_UNSIGNED, off_rank, 0,
      MPI_COMM_WORLD, &reqs[ir]);

    MPI_Irecv(&recv_from_rank[ir], 1, MPI_UNSIGNED, off_rank, 0,
      MPI_COMM_WORLD, &reqs[ir+n_off_rank]);
  }

  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE);

  // now send the buffers and post receives
  for (uint32_t ir=0; ir<n_off_rank; ++ir) {
    int off_rank = proc_map[ir];
    MPI_Isend(send_cell[ir].get_buffer(), send_to_rank[ir], MPI_Cell,
      off_rank, 0, MPI_COMM_WORLD, &reqs[ir]);

    recv_cell[ir].resize(recv_from_rank[ir]);

    MPI_Irecv(recv_cell[ir].get_buffer(), recv_from_rank[ir], MPI_Cell,
      off_rank, 0, MPI_COMM_WORLD, &reqs[ir+n_off_rank]);
  }

  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE);

  send_cell.clear();

  for (uint32_t ir=0; ir<n_off_rank; ++ir) {
    vector<Cell> new_cells = recv_cell[ir].get_object();
    for (uint32_t i = 0; i< new_cells.size(); ++i) {
      mesh->add_mesh_cell(new_cells[i]);
    }
  }

  recv_cell.clear();

  // update the cell list on each processor (use an identity mapped partition
  // vector)
  std::vector<int> part(ncell_on_rank, rank);
  mesh->set_post_decomposition_mesh_cells(part);

  // if using grips of cell data, overdecompose mesh on rank to smaller chunks
  overdecompose_mesh(mesh, grip_size);

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh->get_n_local_cells();

  mesh->set_global_bound(0, n_cell_post_decomp);

  // set the grip ID to be the cells at the center of grips using global cell
  // indices
  mesh->set_grip_ID_using_cell_index();

  // get the maximum grip size (should not be used in replicated mode)
  uint32_t max_grip_size = mesh->get_max_grip_size();
  mesh->set_max_grip_size(max_grip_size);

  // prepend zero to the prefix array to make it a standard bounds array
  vector<uint32_t> prefix_cells_proc;
  prefix_cells_proc.push_back(0);
  prefix_cells_proc.push_back(n_cell_post_decomp);
  mesh->set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting
  unordered_map<uint32_t, uint32_t> local_grip_map = mesh->get_grip_map();
  unordered_map<uint32_t, uint32_t> local_map = mesh->get_new_global_index_map();

  // now update the indices of local IDs
  mesh->renumber_local_cell_indices(local_map, local_grip_map);

  // reallocate mesh data in new MPI window and delete the old vector object
  bool rep_flag = true;
  mesh->make_MPI_window(rep_flag);

  // clean up dynamically allocated memory
  delete[] reqs;
}

#endif // decompose_mesh_h
//---------------------------------------------------------------------------//
// end of decompose_mesh.h
//---------------------------------------------------------------------------//
