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
#include "input.h"
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
  Proto_Mesh(const Input &input)
      : ngx(input.get_global_n_x_cells()), ngy(input.get_global_n_y_cells()),
        ngz(input.get_global_n_z_cells()) {

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
    uint32_t cell_id_begin = 0;
    uint32_t cell_id_end = n_global;


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
  uint32_t get_global_num_cells(void) const { return n_global; }

  void print(void) const {
    for (uint32_t i = 0; i < n_cell; i++)
      cell_list[i].print();
  }

  uint32_t get_local_ID(const uint32_t &index) const {
    return index ;
  }

  uint32_t get_global_ID(const uint32_t &local_index) const {
    return local_index;
  }

  Proto_Cell get_cell(const uint32_t local_index) {
    return cell_list[local_index];
  }

  uint32_t get_global_n_x_faces(void) const { return ngx + 1; }
  uint32_t get_global_n_y_faces(void) const { return ngy + 1; }
  uint32_t get_global_n_z_faces(void) const { return ngz + 1; }
  uint32_t get_global_n_x(void) const { return ngx; }
  uint32_t get_global_n_y(void) const { return ngy; }
  uint32_t get_global_n_z(void) const { return ngz; }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  const std::vector<Proto_Cell> &get_cell_list(void) const { return cell_list; }

  //--------------------------------------------------------------------------//
  // member variables
  //--------------------------------------------------------------------------//
private:
  uint32_t ngx; //!< Number of global x sizes
  uint32_t ngy; //!< Number of global y sizes
  uint32_t ngz; //!< Number of global z sizes

  uint32_t n_cell;   //!< Number of local cells
  uint32_t n_global; //!< Nuber of global cells

  std::vector<Proto_Cell> cell_list;     //!< On processor proto-cells

  std::vector<Region> regions; //!< Vector of regions in the problem
  std::unordered_map<uint32_t, uint32_t>
      region_ID_to_index; //!< Maps region ID to index
};

#endif // proto_mesh_h_
//---------------------------------------------------------------------------//
// end of proto_mesh.h
//---------------------------------------------------------------------------//
