//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   replicate_mesh.h
 * \author Alex Long
 * \date   June 17 2015
 * \brief  Function to replicate mesh on all ranks
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef replicate_mesh_h_
#define replicate_mesh_h_

#include <algorithm>
#include <iostream>
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

//! Create replicated mesh by giving all cells to all other processors, renumber
// mesh and communicate renumbering
void replicate_mesh(Proto_Mesh &mesh, const MPI_Types &mpi_types,
                    const Info &mpi_info) {
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

  // set the grip ID to be the cells at the center of grips using global cell
  // indices
  mesh.set_grip_ID_using_cell_index();

  // prepend zero to the prefix array to make it a standard bounds array
  vector<uint32_t> prefix_cells_proc;
  prefix_cells_proc.push_back(0);
  prefix_cells_proc.push_back(n_cell_post_decomp);
  mesh.set_off_rank_bounds(prefix_cells_proc);

  // change global indices to match a simple number system for easy sorting
  unordered_map<uint32_t, uint32_t> local_grip_map = mesh.get_grip_map();
  unordered_map<uint32_t, uint32_t> local_map = mesh.get_new_global_index_map();

  // now update the indices of local IDs
  mesh.renumber_local_cell_indices(local_map, local_grip_map);

  // clean up dynamically allocated memory
  delete[] reqs;
}

#endif
