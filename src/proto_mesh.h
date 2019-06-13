//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   proto_mesh.h
 * \author Alex Long
 * \date   November 5 2018
 * \brief  Object that holds connectivity and decomposition data for mesh
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef proto_mesh_h_
#define proto_mesh_h_

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "constants.h"
#include "imc_state.h"
#include "info.h"
#include "input.h"
#include "mpi_types.h"
#include "proto_cell.h"

//==============================================================================
/*!
 * \class Proto_Mesh
 * \brief Manages data access and decomposition for primitive mesh
 *
 * Using an Input class, make the mesh with the correct material properties
 * for each region. The mesh numbering and mapping between global IDs and local
 * indices are all determined with the aid of Metis in the decompose_mesh
 * function. The Proto_Mesh does not hold any physical data (i.e. opacity, heat
 * capacity, temperature).
 *
 */
//==============================================================================
class Proto_Mesh {

public:
  //! constructor
  Proto_Mesh(const Input &input, const MPI_Types &mpi_types,
             const Info &mpi_info)
      : ngx(input.get_global_n_x_cells()), ngy(input.get_global_n_y_cells()),
        ngz(input.get_global_n_z_cells()), rank(mpi_info.get_rank()),
        n_rank(mpi_info.get_n_rank()), n_off_rank(n_rank - 1),
        silo_x(input.get_silo_x()), silo_y(input.get_silo_y()),
        silo_z(input.get_silo_z()) {
    using Constants::bc_type;
    using Constants::ELEMENT;
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    using Constants::Z_NEG;
    using Constants::Z_POS;
    using std::vector;

    double dx, dy, dz;

    regions = input.get_regions();
    // map region IDs to index in the region
    for (uint32_t i = 0; i < regions.size(); i++) {
      region_ID_to_index[regions[i].get_ID()] = i;
    }

    vector<bc_type> bc(6);
    bc[X_POS] = input.get_bc(X_POS);
    bc[X_NEG] = input.get_bc(X_NEG);
    bc[Y_POS] = input.get_bc(Y_POS);
    bc[Y_NEG] = input.get_bc(Y_NEG);
    bc[Z_POS] = input.get_bc(Z_POS);
    bc[Z_NEG] = input.get_bc(Z_NEG);

    uint32_t global_count = 0; // global cell count

    // this rank's cells
    n_global = ngx * ngy * ngz;
    uint32_t cell_id_begin = floor(rank * double(n_global) / double(n_rank));
    uint32_t cell_id_end =
        floor((rank + 1) * double(n_global) / double(n_rank));

    uint32_t on_rank_count = 0;

    uint32_t n_x_div = input.get_n_x_divisions();
    uint32_t n_y_div = input.get_n_y_divisions();
    uint32_t n_z_div = input.get_n_z_divisions();

    Region region;
    uint32_t region_index, nx, ny, nz;
    double x_start, y_start, z_start;
    double x_cell_end;
    double y_cell_end;
    double z_cell_end;

    uint32_t g_i = 0; //! Global x index
    uint32_t g_j = 0; //! Global y index
    uint32_t g_k = 0; //! Global z index

    for (uint32_t iz_div = 0; iz_div < n_z_div; iz_div++) {
      dz = input.get_dz(iz_div);
      nz = input.get_z_division_cells(iz_div);
      z_start = input.get_z_start(iz_div);
      for (uint32_t k = 0; k < nz; k++) {
        g_j = 0;
        for (uint32_t iy_div = 0; iy_div < n_y_div; iy_div++) {
          dy = input.get_dy(iy_div);
          ny = input.get_y_division_cells(iy_div);
          y_start = input.get_y_start(iy_div);
          for (uint32_t j = 0; j < ny; j++) {
            g_i = 0;
            for (uint32_t ix_div = 0; ix_div < n_x_div; ix_div++) {
              dx = input.get_dx(ix_div);
              nx = input.get_x_division_cells(ix_div);
              x_start = input.get_x_start(ix_div);
              for (uint32_t i = 0; i < nx; i++) {
                if (global_count >= cell_id_begin &&
                    global_count < cell_id_end) {
                  Proto_Cell e;
                  // find the region for this cell
                  region_index = input.get_region_index(ix_div, iy_div, iz_div);
                  region = regions[region_index];

                  // set ending coordinates explicity to match the start of
                  // the next division to avoid weird roundoff errors
                  if (i == nx - 1 && ix_div != n_x_div - 1)
                    x_cell_end = input.get_x_start(ix_div + 1);
                  else
                    x_cell_end = x_start + (i + 1) * dx;

                  if (j == ny - 1 && iy_div != n_y_div - 1)
                    y_cell_end = input.get_y_start(iy_div + 1);
                  else
                    y_cell_end = y_start + (j + 1) * dy;

                  if (k == nz - 1 && iz_div != n_z_div - 1)
                    z_cell_end = input.get_z_start(iz_div + 1);
                  else
                    z_cell_end = z_start + (k + 1) * dz;

                  e.set_coor(x_start + i * dx, x_cell_end, y_start + j * dy,
                             y_cell_end, z_start + k * dz, z_cell_end);
                  e.set_ID(global_count);
                  e.set_region_ID(region.get_ID());

                  // set the global index for SILO plotting--this will always
                  // be the current global count (g_i +g_j*ngx + g_k*(ngy_*ngz))
                  e.set_silo_index(global_count);

                  // set neighbors in x direction
                  if (g_i < (ngx - 1)) {
                    e.set_neighbor(X_POS, global_count + 1);
                    e.set_bc(X_POS, ELEMENT);
                  } else {
                    e.set_neighbor(X_POS, global_count);
                    e.set_bc(X_POS, bc[X_POS]);
                  }
                  if (g_i > 0) {
                    e.set_neighbor(X_NEG, global_count - 1);
                    e.set_bc(X_NEG, ELEMENT);
                  } else {
                    e.set_neighbor(X_NEG, global_count);
                    e.set_bc(X_NEG, bc[X_NEG]);
                  }

                  // set neighbors in y direction
                  if (g_j < (ngy - 1)) {
                    e.set_neighbor(Y_POS, global_count + ngx);
                    e.set_bc(Y_POS, ELEMENT);
                  } else {
                    e.set_neighbor(Y_POS, global_count);
                    e.set_bc(Y_POS, bc[Y_POS]);
                  }
                  if (g_j > 0) {
                    e.set_neighbor(Y_NEG, global_count - ngx);
                    e.set_bc(Y_NEG, ELEMENT);
                  } else {
                    e.set_neighbor(Y_NEG, global_count);
                    e.set_bc(Y_NEG, bc[Y_NEG]);
                  }

                  // set neighbors in z direction
                  if (g_k < (ngz - 1)) {
                    e.set_neighbor(Z_POS, global_count + ngx * ngy);
                    e.set_bc(Z_POS, ELEMENT);
                  } else {
                    e.set_neighbor(Z_POS, global_count);
                    e.set_bc(Z_POS, bc[Z_POS]);
                  }
                  if (g_k > 0) {
                    e.set_neighbor(Z_NEG, global_count - ngx * ngy);
                    e.set_bc(Z_NEG, ELEMENT);
                  } else {
                    e.set_neighbor(Z_NEG, global_count);
                    e.set_bc(Z_NEG, bc[Z_NEG]);
                  }

                  // add cell to mesh
                  cell_list.push_back(e);
                  // increment on rank count
                  on_rank_count++;
                } // end if on processor check
                global_count++;
                g_i++;
              } // end i loop
            }   // end x division loop
            g_j++;
          } // end j loop
        }   // end y division loop
        g_k++;
      } // end k loop
    }   // end z division loop
    n_cell = on_rank_count;
  }

  // destructor
  ~Proto_Mesh() { /* empty */
  }

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//
  uint32_t get_n_local_cells(void) const { return n_cell; }
  uint32_t get_my_rank(void) const { return rank; }
  uint32_t get_offset(void) const { return on_rank_start; }
  uint32_t get_global_num_cells(void) const { return n_global; }

  std::vector<uint32_t> get_off_rank_bounds(void) { return off_rank_bounds; }

  void print(void) const {
    for (uint32_t i = 0; i < n_cell; i++)
      cell_list[i].print();
  }

  //! returns a mapping of old cell indices to new simple global indices
  std::unordered_map<uint32_t, uint32_t> get_new_global_index_map(void) const {
    std::unordered_map<uint32_t, uint32_t> local_map;
    uint32_t g_ID;
    for (uint32_t i = 0; i < n_cell; i++) {
      g_ID = cell_list[i].get_ID();
      local_map[g_ID] = i + on_rank_start;
    }
    return local_map;
  }

  //! returns a mapping of old cell indices to new simple global indices
  std::unordered_map<uint32_t, uint32_t> get_grip_map(void) const {
    std::unordered_map<uint32_t, uint32_t> local_grip_map;
    uint32_t g_ID;
    for (uint32_t i = 0; i < n_cell; i++) {
      g_ID = cell_list[i].get_ID();
      local_grip_map[g_ID] = cell_list[i].get_grip_ID();
    }
    return local_grip_map;
  }

  //! Gets cell from vector list of cells before it's deleted
  Proto_Cell get_pre_window_allocation_cell(const uint32_t &local_ID) const {
    return cell_list[local_ID];
  }

  uint32_t get_local_ID(const uint32_t &index) const {
    return index - on_rank_start;
  }

  uint32_t get_global_ID(const uint32_t &local_index) const {
    return on_rank_start + local_index;
  }

  bool on_processor(const uint32_t &index) const {
    return (index >= on_rank_start) && (index <= on_rank_end);
  }

  uint32_t get_global_n_x_faces(void) const { return ngx + 1; }
  uint32_t get_global_n_y_faces(void) const { return ngy + 1; }
  uint32_t get_global_n_z_faces(void) const { return ngz + 1; }
  uint32_t get_global_n_x(void) const { return ngx; }
  uint32_t get_global_n_y(void) const { return ngy; }
  uint32_t get_global_n_z(void) const { return ngz; }

  const std::vector<float> &get_silo_x(void) const { return silo_x; }
  const std::vector<float> &get_silo_y(void) const { return silo_y; }
  const std::vector<float> &get_silo_z(void) const { return silo_z; }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! set the grip ID to be the global index of the cell at the center of the
  // grip
  void set_grip_ID_using_cell_index(void) {
    using std::max;
    using std::unordered_map;

    uint32_t new_grip_ID, grip_end_index;
    // start by looking at the first grip
    uint32_t current_grip_ID = cell_list.front().get_grip_ID();
    uint32_t grip_start_index = 0;
    uint32_t grip_count = 0;

    unordered_map<uint32_t, uint32_t> start_index_to_count;

    // map the starting index of cells with the same grip to the number in
    // that grip
    for (uint32_t i = 0; i < n_cell; ++i) {
      Proto_Cell &cell = cell_list[i];
      if (cell.get_grip_ID() != current_grip_ID) {
        grip_start_index = i;
        current_grip_ID = cell.get_grip_ID();
      }

      // if in grip, incerement count...
      if (start_index_to_count.find(grip_start_index) !=
          start_index_to_count.end()) {
        start_index_to_count[grip_start_index]++;
      }
      // otherwise initialize count to 1
      else {
        start_index_to_count[grip_start_index] = 1;
      }
    }

    // set grip ID using the map of start indices and number of cells
    for (auto const &map_i : start_index_to_count) {
      grip_start_index = map_i.first;
      grip_count = map_i.second;
      // set the new grip ID to be the index of the cell at the center of
      // this grip for odd grip sizes and one above center for even grip
      // sizes (for convenience in parallel comm)
      new_grip_ID = on_rank_start + grip_start_index + grip_count / 2;
      // update max grip size
      // loop over cells in grip and set new ID
      grip_end_index = grip_start_index + grip_count;
      for (uint32_t j = grip_start_index; j < grip_end_index; ++j)
        cell_list[j].set_grip_ID(new_grip_ID);
    }
  }

  //! set the global ID of the start and end cell on this rank
  void set_global_bound(uint32_t _on_rank_start, uint32_t _on_rank_end) {
    on_rank_start = _on_rank_start;
    on_rank_end = _on_rank_end;
  }

  //! set the global ID starting indices for all ranks
  void set_off_rank_bounds(std::vector<uint32_t> _off_rank_bounds) {
    off_rank_bounds = _off_rank_bounds;
  }

  //! Renumber the local cell IDs and connectivity of local cells after
  // decomposition using simple global numbering
  void renumber_local_cell_indices(
      std::unordered_map<uint32_t, uint32_t> local_map,
      std::unordered_map<uint32_t, uint32_t> local_grip_map) {

    using Constants::dir_type;
    using std::unordered_map;

    uint32_t next_index;
    // grip index is already set for cells, neighbors are not set!
    // renumber global cell index
    for (uint32_t i = 0; i < n_cell; ++i) {
      Proto_Cell &cell = cell_list[i];
      cell.set_ID(i + on_rank_start);
      for (uint32_t d = 0; d < 6; ++d) {
        // get the un-remapped next index
        next_index = cell.get_next_cell(d);
        // remap it
        cell.set_neighbor(dir_type(d), local_map[next_index]);
        cell.set_grip_neighbor(dir_type(d), local_grip_map[next_index]);
        if (local_map[next_index] > off_rank_bounds.back() ||
            local_grip_map[next_index] > off_rank_bounds.back())
          std::cout << "this is bad, g > global bounds!" << std::endl;
      } // end direction
    }   // end cell
  }

  //! Remove old mesh cells after decomposition and communication of new cells
  void set_post_decomposition_mesh_cells(const std::vector<int> &partition) {
    using std::vector;
    vector<Proto_Cell> new_mesh;

    // all cells that were assigned to this rank are still part of the mesh
    for (uint32_t i = 0; i < cell_list.size(); ++i) {
      if (partition[i] == rank)
        new_mesh.push_back(cell_list[i]);
    }

    // add cells received by other processors
    new_mesh.insert(new_mesh.end(), new_cell_list.begin(), new_cell_list.end());

    // reassigned the cell_list vector and clear other vectors
    cell_list = new_mesh;
    n_cell = cell_list.size();
    new_cell_list.clear();

    // sort based on global cell ID
    sort(cell_list.begin(), cell_list.end());
  }

  //! sort pre-winodw allocation cell vector based on the grip ID of each cell
  void sort_cells_by_grip_ID(void) {
    using std::sort;
    // sort based on global cell ID
    sort(cell_list.begin(), cell_list.end(), Proto_Cell::sort_grip_ID);
  }

  //! Add mesh cell (used during decomposition, not parallel communication)
  void add_mesh_cell(const Proto_Cell new_cell) {
    new_cell_list.push_back(new_cell);
  }

  const std::vector<Proto_Cell> &get_cell_list(void) const { return cell_list; }

  //--------------------------------------------------------------------------//
  // member variables
  //--------------------------------------------------------------------------//
private:
  uint32_t ngx; //!< Number of global x sizes
  uint32_t ngy; //!< Number of global y sizes
  uint32_t ngz; //!< Number of global z sizes

  int32_t rank;       //!< MPI rank of this mesh
  int32_t n_rank;     //!< Number of global ranks
  int32_t n_off_rank; //!< Number of other ranks

  const std::vector<float> &silo_x; //!< Global array x face locations for SILO
  const std::vector<float> &silo_y; //!< Global array y face locations for SILO
  const std::vector<float> &silo_z; //!< Global array z face locations for SILO

  uint32_t n_cell;   //!< Number of local cells
  uint32_t n_global; //!< Nuber of global cells

  uint32_t on_rank_start; //!< Start of global index on rank
  uint32_t on_rank_end;   //!< End of global index on rank

  std::vector<Proto_Cell> cell_list;     //!< On processor proto-cells
  std::vector<Proto_Cell> new_cell_list; //!< New received proto-cells
  std::vector<uint32_t>
      off_rank_bounds; //!< Ending value of global ID for each rank

  std::vector<Region> regions; //!< Vector of regions in the problem
  std::unordered_map<uint32_t, uint32_t>
      region_ID_to_index; //!< Maps region ID to index
};

#endif // proto_mesh_h_
//---------------------------------------------------------------------------//
// end of proto_mesh.h
//---------------------------------------------------------------------------//
