//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   decompose_mesh.h
 * \author Alex Long
 * \date   June 17 2015
 * \brief  Functions to decompose mesh with Metis or cubes
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef decompose_mesh_h_
#define decompose_mesh_h_

#include <algorithm>
#include <iostream>
#include <metis.h>
#include <mpi.h>
#include <numeric>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "buffer.h"
#include "mpi_types.h"
#include "proto_mesh.h"
#include "timer.h"

//----------------------------------------------------------------------------//
//! Print the mesh information for each rank, one at a time
void print_MPI_out(const Proto_Mesh &mesh, const uint32_t rank,
                   const uint32_t size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank < size; ++p_rank) {
    if (rank == p_rank) {
      mesh.print();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! partition a mesh with metis
std::vector<int> metis_partition(Proto_Mesh &mesh, int &edgecut, const int rank,
                                 const int n_rank, const MPI_Types &mpi_types) {

  using Constants::X_NEG;
  using Constants::X_POS;
  using Constants::Y_NEG;
  using Constants::Y_POS;
  using Constants::Z_NEG;
  using Constants::Z_POS;
  using std::vector;

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  constexpr int cell_tag = 10100;
  constexpr int part_tag = 21023;
  MPI_Datatype MPI_Proto_Cell = mpi_types.get_proto_cell_type();

  // get the number of cells on each processor
  // vtxdist has number of vertices on each rank, same for all ranks
  vector<int> start_ncells(n_rank, 0);
  vector<int> vtxdist(n_rank, 0);

  uint32_t ncell_on_rank = mesh.get_n_local_cells();
  start_ncells[rank] = ncell_on_rank;

  MPI_Allreduce(MPI_IN_PLACE, &start_ncells[0], n_rank, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  partial_sum(start_ncells.begin(), start_ncells.end(), vtxdist.begin());
  vtxdist.insert(vtxdist.begin(), 0);

  uint32_t n_global_cells = vtxdist.back();

  // return partition of each cell on rank
  std::vector<int> part(ncell_on_rank);

  // non-root ranks gather send cells to root
  if (rank != 0) {
    const std::vector<Proto_Cell> send_cells =
        mesh.get_pre_window_allocation_cells();
    MPI_Send(send_cells.data(), ncell_on_rank, MPI_Proto_Cell, 0, cell_tag,
             MPI_COMM_WORLD);
    MPI_Recv(part.data(), ncell_on_rank, MPI_INT, 0, part_tag, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    edgecut = 1;
  }
  // root rank gathers all cells and uses METIS to setup partitions
  else {
    std::vector<Proto_Cell> all_cells(n_global_cells);
    const std::vector<Proto_Cell> send_cells =
        mesh.get_pre_window_allocation_cells();
    std::copy(send_cells.begin(), send_cells.end(), all_cells.begin());
    for (int irank = 1; irank < n_rank; ++irank) {
      MPI_Recv(&all_cells[vtxdist[irank]], start_ncells[irank], MPI_Proto_Cell,
               irank, cell_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // do partitioning
    // build adjacency list needed for Metis call for each rank
    vector<int> xadj;
    vector<int> adjncy;
    int adjncy_ctr = 0;
    uint32_t global_index;
    for (auto &icell : all_cells) {
      global_index = icell.get_global_index();
      uint32_t xm_neighbor = icell.get_next_cell(X_NEG);
      uint32_t xp_neighbor = icell.get_next_cell(X_POS);
      uint32_t ym_neighbor = icell.get_next_cell(Y_NEG);
      uint32_t yp_neighbor = icell.get_next_cell(Y_POS);
      uint32_t zm_neighbor = icell.get_next_cell(Z_NEG);
      uint32_t zp_neighbor = icell.get_next_cell(Z_POS);

      xadj.push_back(
          adjncy_ctr); // starting index in xadj for this cell's nodes
      if (xm_neighbor != global_index) {
        adjncy.push_back(xm_neighbor);
        adjncy_ctr++;
      }
      if (xp_neighbor != global_index) {
        adjncy.push_back(xp_neighbor);
        adjncy_ctr++;
      }
      if (ym_neighbor != global_index) {
        adjncy.push_back(ym_neighbor);
        adjncy_ctr++;
      }
      if (yp_neighbor != global_index) {
        adjncy.push_back(yp_neighbor);
        adjncy_ctr++;
      }
      if (zm_neighbor != global_index) {
        adjncy.push_back(zm_neighbor);
        adjncy_ctr++;
      }
      if (zp_neighbor != global_index) {
        adjncy.push_back(zp_neighbor);
        adjncy_ctr++;
      }
    }
    xadj.push_back(adjncy_ctr);

    int ncon = 1;
    int n_parts = n_rank; // sub-domains = n_rank

    int rank_options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(rank_options);
    rank_options[METIS_OPTION_NUMBERING] = 0; // C-style numbering
    //rank_options[1] = 3;    // output level
    //rank_options[2] = 1242; // random number seed

    int signed_n_global_cells = n_global_cells;

    std::vector<int> global_part(n_global_cells);

    int metis_return = METIS_PartGraphKway(
        &signed_n_global_cells, // number of on-rank vertices
        &ncon,                  // weight of vertices
        &xadj[0],               // how cells are stored locally
        &adjncy[0],             // how cells are stored loccaly
        NULL,                   // weight of vertices
        NULL,                   // size of vertices for comm volume
        NULL,                   // weight of the edges
        &n_parts,               // number of ranks (partitions)
        NULL,                   // tpwgts (NULL = equal weight domains)
        NULL,                   // unbalance in v-weight (NULL=1.001)
        rank_options,           // options array
        &edgecut,               // OUTPUT: Number of edgecuts
        &global_part[0]);       // OUTPUT: rank of each cell

    if (metis_return == METIS_ERROR_INPUT)
      std::cout << "METIS: Input error" << std::endl;
    else if (metis_return == METIS_ERROR_MEMORY)
      std::cout << "METIS: Memory error" << std::endl;
    else if (metis_return == METIS_ERROR)
      std::cout << "METIS: Input error" << std::endl;

    if (!edgecut) {
      std::cout << "ERROR: No partitioning occured, run with more cells or in";
      std::cout << "replicated mode. Exiting..." << std::endl;
      exit(EXIT_FAILURE);
    }

    // send partitioning to other ranks
    for (int irank = 1; irank < n_rank; ++irank) {
      MPI_Send(&global_part[vtxdist[irank]], start_ncells[irank], MPI_INT,
               irank, part_tag, MPI_COMM_WORLD);
    }

    // copy out root ranks partitioning
    std::copy(global_part.begin(), global_part.begin() + ncell_on_rank,
              part.begin());
  }

  return part;
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
std::vector<int> cube_partition(Proto_Mesh &mesh, const int rank,
                                const int n_rank) {

  uint32_t nx = mesh.get_global_n_x();
  uint32_t ny = mesh.get_global_n_y();
  uint32_t nz = mesh.get_global_n_z();

  // make sure this mesh is a cube
  int rank_side = std::floor(std::pow(n_rank, 1.0 / 3.0)) + 1;
  if (rank_side * rank_side * rank_side == n_rank) {
    if (rank == 0)
      std::cout << "Cube decomposition OK, ranks per edge = " << rank_side
                << std::endl;
  } else if ((rank_side - 1) * (rank_side - 1) * (rank_side - 1) == n_rank) {
    rank_side--;
    if (rank == 0)
      std::cout << "Cube decomposition OK, ranks per edge = " << rank_side
                << std::endl;
  } else {
    if (rank == 0) {
      std::cout
          << "Can't use simple cube decomposition for this number of ranks";
      std::cout << " Exiting..." << std::endl;
    }
    exit(EXIT_FAILURE);
  }

  uint32_t ncell_on_rank = mesh.get_n_local_cells();
  std::vector<int> part(ncell_on_rank, rank);
  std::unordered_set<int> acceptor_ranks;

  for (uint32_t i = 0; i < ncell_on_rank; ++i) {
    const Proto_Cell &cell = mesh.get_pre_window_allocation_cell(i);
    // get the correct rank given the cell global_index
    uint32_t index = cell.get_global_index();
    uint32_t z = index / (nx * ny);
    uint32_t y = (index - z * (nx * ny)) / nx;
    uint32_t x = index - z * (nx * ny) - y * nx;
    uint32_t rx = floor(double(x) / double(nx) * rank_side);
    uint32_t ry = floor(double(y) / double(ny) * rank_side);
    uint32_t rz = floor(double(z) / double(nz) * rank_side);
    part[i] = rz * rank_side * rank_side + ry * rank_side + rx;
    if (part[i] != rank)
      acceptor_ranks.insert(part[i]);
  }
  return part;
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Send given partitioning scheme
void exchange_cells_post_partitioning(const int rank,
                                      const MPI_Types &mpi_types,
                                      Proto_Mesh &mesh,
                                      const std::vector<int> &part) {
  using std::vector;
  MPI_Datatype MPI_Proto_Cell = mpi_types.get_proto_cell_type();
  uint32_t ncell_on_rank = mesh.get_n_local_cells();

  std::unordered_set<int> acceptor_ranks;
  // use the partition to build the acceptor list
  for (uint32_t i = 0; i < ncell_on_rank; ++i) {
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
  MPI_Win_lock_all(assert, win);
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
  vector<Buffer<Proto_Cell>> send_cell(n_acceptors);
  vector<Buffer<Proto_Cell>> recv_cell(n_donors);
  vector<int> recv_from_rank(n_donors, 0);
  vector<int> send_to_rank(n_acceptors, 0);
  vector<MPI_Request> reqs(n_donors + n_acceptors);
  vector<MPI_Status> status(n_donors + n_acceptors);
  std::unordered_map<int, int> donor_rank_size;

  // send sizes
  int isend = 0;
  for (auto &ir : acceptor_ranks) {
    // make list of cells to send to off_rank
    vector<Proto_Cell> send_list;
    for (uint32_t i = 0; i < ncell_on_rank; ++i) {
      if (part[i] == ir)
        send_list.push_back(mesh.get_pre_window_allocation_cell(i));
    }
    send_to_rank[isend] = send_list.size();
    send_cell[isend].fill(send_list);

    MPI_Isend(&send_to_rank[isend], 1, MPI_UNSIGNED, ir, 0, MPI_COMM_WORLD,
              &reqs[isend]);
    isend++;
  }

  // receive sizes, ranks
  for (int i = 0; i < n_donors; ++i) {
    MPI_Irecv(&recv_from_rank[i], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &reqs[n_acceptors + i]);
  }

  MPI_Waitall(reqs.size(), &reqs[0], &status[0]);
  MPI_Barrier(MPI_COMM_WORLD);

  // map donor rank to message size
  for (int i = 0; i < n_donors; ++i)
    donor_rank_size[status[n_acceptors + i].MPI_SOURCE] = recv_from_rank[i];

  // now send the buffers and post receives
  isend = 0;
  for (auto &ir : acceptor_ranks) {
    MPI_Isend(send_cell[isend].get_buffer(), send_to_rank[isend],
              MPI_Proto_Cell, ir, 0, MPI_COMM_WORLD, &reqs[isend]);
    isend++;
  }
  int ireceive = 0;
  for (auto &ir : donor_rank_size) {
    recv_cell[ireceive].resize(ir.second);
    MPI_Irecv(recv_cell[ireceive].get_buffer(), ir.second, MPI_Proto_Cell,
              ir.first, 0, MPI_COMM_WORLD, &reqs[n_acceptors + ireceive]);
    ireceive++;
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < n_donors; ++i) {
    vector<Proto_Cell> new_cells = recv_cell[i].get_object();
    for (uint32_t i = 0; i < new_cells.size(); ++i) {
      mesh.add_mesh_cell(new_cells[i]);
    }
  }
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Use one-sided communication to share new cell indices in mesh connectivity
void remap_cell_and_grip_indices_rma(Proto_Mesh &mesh, const int rank,
                                     const int n_rank) {
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;
  if (rank == 0)
    std::cout << "remapping with rma..." << std::endl;

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh.get_n_local_cells();
  vector<uint32_t> out_cells_proc(n_rank, 0);
  out_cells_proc[rank] = mesh.get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, &out_cells_proc[0], n_rank, MPI_UNSIGNED, MPI_SUM,
                MPI_COMM_WORLD);

  // prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(n_rank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(),
              prefix_cells_proc.begin());

  // set global numbering
  uint32_t g_start = prefix_cells_proc[rank] - n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank] - 1;
  mesh.set_global_bound(g_start, g_end);

  // prepend zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert(prefix_cells_proc.begin(), 0);
  mesh.set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting,
  // this involves sending maps to each processor to get new indicies
  unordered_map<uint32_t, uint32_t> local_boundary_map =
      mesh.get_boundary_nodes();
  unordered_map<uint32_t, uint32_t> local_map = mesh.get_new_global_index_map();
  unordered_set<uint32_t> boundary_indices = mesh.get_boundary_neighbors();
  unordered_map<uint32_t, uint32_t> boundary_map;

  // create the RMA memory windows to allow all ranks to write the new index
  // of a cell at window[old_index]
  MPI_Win index_win;
  uint32_t *new_index;
  MPI_Info decomp_info;
  MPI_Info_create(&decomp_info);
  MPI_Info_set(decomp_info, "same_disp_unit", "true");
  MPI_Win_allocate(n_cell_post_decomp * sizeof(uint32_t), sizeof(uint32_t),
                   decomp_info, MPI_COMM_WORLD, &new_index, &index_win);

  // initialize to max uint32_t to check for any errors
  for (uint32_t i = 0; i < n_cell_post_decomp; ++i) {
    new_index[i] = UINT32_MAX;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  int assert = MPI_MODE_NOCHECK; // no conflicting locks on this window
  //int assert = 0;
  MPI_Win_lock_all(assert, index_win);

  uint32_t remapped_index, local_index;
  int target_rank;
  for (auto &imap : local_boundary_map) {
    remapped_index = imap.second;
    target_rank = mesh.get_rank(imap.first);
    local_index = imap.first - prefix_cells_proc[target_rank];
    MPI_Put(&remapped_index, 1, MPI_UNSIGNED, target_rank, local_index, 1,
            MPI_UNSIGNED, index_win);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  // make the memory visible to all ranks
  MPI_Win_flush_all(index_win);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_sync(index_win);
  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<uint32_t> new_boundary_indices(boundary_indices.size());
  std::vector<MPI_Request> i_reqs(boundary_indices.size());
  std::vector<MPI_Request> g_reqs(boundary_indices.size());

  int ib = 0;
  for (auto &iset : boundary_indices) {
    target_rank = mesh.get_rank(iset);
    local_index = iset - prefix_cells_proc[target_rank];
    MPI_Rget(&new_boundary_indices[ib], 1, MPI_UNSIGNED, target_rank,
             local_index, 1, MPI_UNSIGNED, index_win, &i_reqs[ib]);
    ib++;
  }

  MPI_Waitall(i_reqs.size(), &i_reqs[0], MPI_STATUS_IGNORE);
  MPI_Waitall(g_reqs.size(), &g_reqs[0], MPI_STATUS_IGNORE);

  // make sure that all new boundary indices were set
  // correctly (not equal to the dummy initial value of UINT32_MAX)
  for (auto &iset : new_boundary_indices) {
    if (iset == UINT32_MAX) {
      std::cout << "This is bad, RMA remapping failed" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  ib = 0;
  for (auto &iset : boundary_indices) {
    boundary_map[iset] = new_boundary_indices[ib];
    ib++;
  }

  MPI_Win_unlock_all(index_win);
  MPI_Win_free(&index_win);

  // combine local index map with boundary map from communication
  local_map.insert(boundary_map.begin(), boundary_map.end());

  // now update local indices
  mesh.renumber_local_cell_indices(local_map);
}
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! Use two-sided communication and allreduces to share new cell indices in mesh
// connectivity
void remap_cell_and_grip_indices_allreduce(Proto_Mesh &mesh, const int rank,
                                           const int n_rank) {
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;

  if (rank == 0)
    std::cout << "remapping with allreduce..." << std::endl;

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh.get_n_local_cells();
  vector<uint32_t> out_cells_proc(n_rank, 0);
  out_cells_proc[rank] = mesh.get_n_local_cells();
  MPI_Allreduce(MPI_IN_PLACE, &out_cells_proc[0], n_rank, MPI_UNSIGNED, MPI_SUM,
                MPI_COMM_WORLD);

  // prefix sum on out_cells to get global numbering
  vector<uint32_t> prefix_cells_proc(n_rank, 0);
  partial_sum(out_cells_proc.begin(), out_cells_proc.end(),
              prefix_cells_proc.begin());

  // set global numbering
  uint32_t g_start = prefix_cells_proc[rank] - n_cell_post_decomp;
  uint32_t g_end = prefix_cells_proc[rank] - 1;
  mesh.set_global_bound(g_start, g_end);

  // prepend zero to the prefix array to make it a standard bounds array
  prefix_cells_proc.insert(prefix_cells_proc.begin(), 0);
  mesh.set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting,
  unordered_map<uint32_t, uint32_t> local_boundary_map =
      mesh.get_boundary_nodes();

  uint32_t n_boundary = local_boundary_map.size();
  // get the size of boundary cells on each rank, use a prefix sum to get this
  // vectors offsets
  vector<uint32_t> out_bcells_proc(n_rank, 0);
  vector<uint32_t> prefix_bcells_proc(n_rank, 0);
  out_bcells_proc[rank] = n_boundary;
  MPI_Allreduce(MPI_IN_PLACE, &out_bcells_proc[0], n_rank, MPI_UNSIGNED,
                MPI_SUM, MPI_COMM_WORLD);

  partial_sum(out_bcells_proc.begin(), out_bcells_proc.end(),
              prefix_bcells_proc.begin());
  // insert zero at the start to get write location
  prefix_bcells_proc.insert(prefix_bcells_proc.begin(), 0);

  // write the original boundary indices to a global array
  uint32_t n_global_bcells = prefix_bcells_proc.back();
  vector<uint32_t> original_b_indices(n_global_bcells);
  uint32_t start = prefix_bcells_proc[rank];
  for (auto &imap : local_boundary_map)
    original_b_indices[start++] = imap.first;
  MPI_Allreduce(MPI_IN_PLACE, &original_b_indices[0], n_global_bcells,
                MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  // find the indices in this global array for your boundary cells
  auto b_nodes = mesh.get_boundary_neighbors();
  std::map<uint32_t, uint32_t> b_indices; // map indices to original indices
  for (uint32_t i = 0; i < n_global_bcells; ++i) {
    if (b_nodes.count(original_b_indices[i]))
      b_indices[i] = original_b_indices[i];
  }

  // now fill the global array with the real new indices (in the same location
  // as the old index)
  start = prefix_bcells_proc[rank];
  std::fill(original_b_indices.begin(), original_b_indices.end(), 0);
  for (auto &imap : local_boundary_map)
    original_b_indices[start++] = imap.second;
  MPI_Allreduce(MPI_IN_PLACE, &original_b_indices[0], n_global_bcells,
                MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  // iterate through the indices required by your rank
  std::unordered_map<uint32_t, uint32_t> id_old_to_new;
  for (auto &imap : b_indices)
    id_old_to_new[imap.second] = original_b_indices[imap.first];

  unordered_map<uint32_t, uint32_t> local_map = mesh.get_new_global_index_map();
  local_map.insert(id_old_to_new.begin(), id_old_to_new.end());

  // now update local indices
  mesh.renumber_local_cell_indices(local_map);
}

//----------------------------------------------------------------------------//
//! Generate new partitioning with Metis, replications or cubes, send and
// receive cells, renumber mesh and communicate renumbering
void decompose_mesh(Proto_Mesh &mesh, const MPI_Types &mpi_types,
                    const Info &mpi_info, const uint32_t grip_size,
                    const int decomposition_type) {
  using Constants::CUBE;
  using Constants::METIS;
  using std::unordered_map;
  using std::unordered_set;
  Timer t_partition;
  Timer t_remap;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // metis sets this, if it's zero no changes were made to the default
  // partition (cube decomposition always changes paritioning)
  int edgecut = 0;
  if (rank == 0)
    std::cout << "partitioning..." << std::endl;
  t_partition.start_timer("partition");
  // decomposition methods return a partition vector which is the rank of each
  // cell
  std::vector<int> part;
  if (decomposition_type == CUBE) {
    part = cube_partition(mesh, rank, n_rank);
    edgecut = 1;
  } else if (decomposition_type == METIS)
    if (n_rank > 1)
      part = metis_partition(mesh, edgecut, rank, n_rank, mpi_types);
    else {
      part = std::vector<int>(mesh.get_n_local_cells(), 0);
      edgecut = 0;
    }
  else {
    if (rank == 0) {
      std::cout << "Decomposition type not recognized.";
      std::cout << " Exiting..." << std::endl;
    }
    exit(EXIT_FAILURE);
  }

  // if edgecuts are made (edgecut > 0) send cells to other processors
  // otherwise mesh is already partitioned (sets "new_cells" in mesh object)
  if (edgecut)
    exchange_cells_post_partitioning(rank, mpi_types, mesh, part);
  t_partition.stop_timer("partition");

  // update the cell list on each processor
  mesh.set_post_decomposition_mesh_cells(part);

  t_remap.start_timer("remap");
  remap_cell_and_grip_indices_allreduce(mesh, rank, n_rank);
  t_remap.stop_timer("remap");

  if (rank == 0) {
    std::cout << "Partition: " << t_partition.get_time("partition")<<" seconds"
              << std::endl;
    std::cout << "Remap: " << t_remap.get_time("remap") <<" seconds"<<std::endl;
  }
}

//! Create replicated mesh by giving all cells to all other processors, renumber
// mesh and communicate renumbering
void replicate_mesh(Proto_Mesh &mesh, const MPI_Types &mpi_types,
                    const Info &mpi_info, const uint32_t &grip_size) {
  using Constants::X_NEG;
  using Constants::X_POS;
  using Constants::Y_NEG;
  using Constants::Y_POS;
  using Constants::Z_NEG;
  using Constants::Z_POS;
  using std::partial_sum;
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;

  int rank = mpi_info.get_rank();
  int n_rank = mpi_info.get_n_rank();

  // make off processor map
  uint32_t n_off_rank = n_rank - 1; // implicit conversion from int to uint32_t
  unordered_map<int, int> proc_map;
  for (uint32_t i = 0; i < n_off_rank; ++i) {
    int r_index = i + int(int(i) >= rank);
    proc_map[i] = r_index;
  }

  MPI_Datatype MPI_Proto_Cell = mpi_types.get_proto_cell_type();

  // get the number of cells on each processor
  vector<int> start_ncells(n_rank, 0);

  uint32_t ncell_on_rank = mesh.get_n_local_cells();
  start_ncells[rank] = ncell_on_rank;

  MPI_Allreduce(MPI_IN_PLACE, &start_ncells[0], n_rank, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  vector<int> recv_from_rank(n_off_rank, 0);
  vector<int> send_to_rank(n_off_rank, 0);

  MPI_Request *reqs = new MPI_Request[n_off_rank * 2];
  vector<Buffer<Proto_Cell>> send_cell(n_off_rank);
  vector<Buffer<Proto_Cell>> recv_cell(n_off_rank);

  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    // sends
    int off_rank = proc_map[ir];
    // make list of cells to send to off_rank
    vector<Proto_Cell> send_list;
    for (uint32_t i = 0; i < ncell_on_rank; ++i)
      send_list.push_back(mesh.get_pre_window_allocation_cell(i));
    send_to_rank[ir] = send_list.size();
    send_cell[ir].fill(send_list);

    MPI_Isend(&send_to_rank[ir], 1, MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD,
              &reqs[ir]);

    MPI_Irecv(&recv_from_rank[ir], 1, MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD,
              &reqs[ir + n_off_rank]);
  }

  MPI_Waitall(n_off_rank * 2, reqs, MPI_STATUS_IGNORE);

  // now send the buffers and post receives
  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    int off_rank = proc_map[ir];
    MPI_Isend(send_cell[ir].get_buffer(), send_to_rank[ir], MPI_Proto_Cell,
              off_rank, 0, MPI_COMM_WORLD, &reqs[ir]);

    recv_cell[ir].resize(recv_from_rank[ir]);

    MPI_Irecv(recv_cell[ir].get_buffer(), recv_from_rank[ir], MPI_Proto_Cell,
              off_rank, 0, MPI_COMM_WORLD, &reqs[ir + n_off_rank]);
  }

  MPI_Waitall(n_off_rank * 2, reqs, MPI_STATUS_IGNORE);

  send_cell.clear();

  for (uint32_t ir = 0; ir < n_off_rank; ++ir) {
    vector<Proto_Cell> new_cells = recv_cell[ir].get_object();
    for (uint32_t i = 0; i < new_cells.size(); ++i) {
      mesh.add_mesh_cell(new_cells[i]);
    }
  }

  recv_cell.clear();

  // update the cell list on each processor (use an identity mapped partition
  // vector)
  std::vector<int> part(ncell_on_rank, rank);
  mesh.set_post_decomposition_mesh_cells(part);

  // gather the number of cells on each processor
  uint32_t n_cell_post_decomp = mesh.get_n_local_cells();

  mesh.set_global_bound(0, n_cell_post_decomp);

  // prepend zero to the prefix array to make it a standard bounds array
  vector<uint32_t> prefix_cells_proc;
  prefix_cells_proc.push_back(0);
  prefix_cells_proc.push_back(n_cell_post_decomp);
  mesh.set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting
  unordered_map<uint32_t, uint32_t> local_map = mesh.get_new_global_index_map();

  // now update local indices
  mesh.renumber_local_cell_indices(local_map);

  // clean up dynamically allocated memory
  delete[] reqs;
}

#endif // decompose_mesh_h
//---------------------------------------------------------------------------//
// end of decompose_mesh.h
//---------------------------------------------------------------------------//
