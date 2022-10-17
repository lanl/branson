//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   proto_cell.h
 * \author Alex Long
 * \date   October 29 2018
 * \brief  Holds just connectivity and location data for cells
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef proto_cell_h_
#define proto_cell_h_

#include <iostream>
#include <array>
#include <mpi.h>

#include "RNG.h"
#include "config.h"
#include "constants.h"

//==============================================================================
/*!
 * \class Proto_Cell
 * \brief Connectivity and location data needed for partitioning
 *
 * A cartesian mesh cell without the physical data (opacity, temperature, etc).
 * This class is specifically meant to keep the memory high water mark of mesh
 * decomposition low.
 */
//==============================================================================

class Proto_Cell {

public:
  Proto_Cell(void) {}

  ~Proto_Cell(void) {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//
  //! Set input array to center of cell (for mesh decomposition only)
  inline void get_center(float xyz[3]) const {
    xyz[0] = 0.5 * (nodes[0] + nodes[1]);
    xyz[1] = 0.5 * (nodes[2] + nodes[3]);
    xyz[2] = 0.5 * (nodes[4] + nodes[5]);
  }

  //! Get boundary condition type in this direction
  inline Constants::bc_type get_bc(const uint32_t &dir) const {
    return bc[dir];
  }

  //! Get global ID of cell in next direction
  inline uint32_t get_next_cell(const uint32_t &dir) const {
    return e_next[dir];
  }

  //! Return node array (for setting up work packets)
  inline const double *get_node_array(void) const { return nodes.data(); }

  //! Return SILO index (for plotting only)
  inline uint32_t get_silo_index(void) const { return silo_index; }

  //! Return cell volume
  inline double get_volume(void) const {
    return (nodes[1] - nodes[0]) * (nodes[3] - nodes[2]) *
           (nodes[5] - nodes[4]);
  }

  // Return global ID
  inline uint32_t get_ID(void) const { return g_ID; }

  // Return region ID
  inline uint32_t get_region_ID(void) const { return region_ID; }

  //! Override great than operator to sort
  bool operator<(const Proto_Cell &compare) const {
    return g_ID < compare.get_ID();
  }

  //! Print cell data (diagnostic only)
  void print(void) const {
    using Constants::PROCESSOR;
    using std::cout;
    using std::endl;
    int32_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    bool boundary = false;
    for (uint32_t i = 0; i < 6; i++) {
      if (bc[i] == PROCESSOR)
        boundary = true;
    }

    // cout<<g_ID<<" "<<boundary;
    cout<<nodes[0]<<" "<<nodes[1]<<" "<<nodes[2]<<" "<<nodes[3]<<" "<<nodes[4]<<" "<<nodes[5]<<std::endl;
    cout << "Rank: " << my_rank << " Global ID: " << g_ID << endl;
    cout << nodes[0] << " " << nodes[2] << " " << nodes[4];
    cout << " Processor bound: " << boundary << endl;
  }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Set neighbor in a given direction by global cell ID
  void set_neighbor(Constants::dir_type neighbor_dir, uint32_t nbr_g_ID) {
    e_next[neighbor_dir] = nbr_g_ID;
  }

  //! Set boundary conditions for cell in a given direction
  void set_bc(Constants::dir_type direction, Constants::bc_type _bc) {
    bc[direction] = _bc;
  }

  //! Set global ID
  void set_ID(uint32_t _id) { g_ID = _id; }

  //! Set region ID
  void set_region_ID(uint32_t _region_ID) { region_ID = _region_ID; }

  //! Set node loactions
  void set_coor(double x_low, double x_high, double y_low, double y_high,
                double z_low, double z_high) {
    nodes[0] = x_low;
    nodes[1] = x_high;
    nodes[2] = y_low;
    nodes[3] = y_high;
    nodes[4] = z_low;
    nodes[5] = z_high;
  }

  //! Set SILO index (for plotting)
  void set_silo_index(uint32_t _silo_index) { silo_index = _silo_index; }

  std::array<Constants::bc_type, 6> get_bc() const {return bc;}
  std::array<uint32_t, 6> get_e_next() const {return e_next;}
  std::array<double, 6> get_nodes() const {return nodes;}

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:

  uint32_t g_ID; //!< Global ID, valid across all ranks
  uint32_t region_ID; //!< region cell is in (for setting physical properties)
  uint32_t silo_index;      //!< Global index not remappaed, for SILO plotting
  std::array<uint32_t, 6> e_next; //!< Bordering cell, given as global ID

  std::array<Constants::bc_type, 6> bc; //!< Boundary conditions for each face
  std::array<double, 6> nodes;          //!< x_low, x_high, y_low, y_high, z_low, z_high
};

#endif // cell_h_
//---------------------------------------------------------------------------//
// end of cell.h
//---------------------------------------------------------------------------//
