//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cell.h
 * \author Alex Long
 * \date   March 3 2015
 * \brief  Holds cells data and provides basic sampling functions
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef cell_h_
#define cell_h_

#include <iostream>
#include <mpi.h>

#include "RNG.h"
#include "config.h"
#include "constants.h"
#include "proto_cell.h"

template <typename T> int sgn(T val) { return (T(0) < val); }

//==============================================================================
/*!
 * \class Cell
 * \brief Basic geometry unit, holds physical data that is read only during
 * transport.
 *
 * A cartesian mesh cell. Holds location of each node, boundary information for
 * each face, opacity data and temperature data. The temperature data is stored
 * here but is not essential information for transport
 */
//==============================================================================

class Cell {

public:
  Cell(void) {
    op_a = 0.0;
    op_s = 0.0;
    f = 0.0;
    g_ID = 0; // won't show up in finds
  }

  explicit Cell(const Proto_Cell &proto_cell)
   : g_ID(proto_cell.get_ID()), region_ID(proto_cell.get_region_ID()),
    grip_ID(proto_cell.get_grip_ID()), silo_index(proto_cell.get_silo_index()),
    e_next(proto_cell.get_e_next()),
    grip_next(proto_cell.get_grip_next()),
    nodes(proto_cell.get_nodes()),
    bc(proto_cell.get_bc()),
    abs_groups(),
    sct_groups(),
    op_a(0.0), op_s(0.0), f(0.0), rho(0.0), T_e(0.0), T_r(0.0), T_s(0.0)
  {
    abs_groups.fill(0.0);
    sct_groups.fill(0.0);
  }

  ~Cell(void) {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//
  //! Get boundary condition type in this direction
  inline Constants::bc_type get_bc(const uint32_t &dir) const {
    return bc[dir];
  }

  //! Get global ID of cell in next direction
  inline uint32_t get_next_cell(const uint32_t &dir) const {
    return e_next[dir];
  }

  //! Get grip ID of cell in next direction
  inline uint32_t get_next_grip(const uint32_t &dir) const {
    return grip_next[dir];
  }

  //! Return a distance to boundary and set surface crossing given
  // position and angle
  inline double get_distance_to_boundary(const double *pos, const double *angle,
                                         uint32_t &surface_cross) const {
    double min_dist = 1.0e16;
    double dist = 0.0;
    uint32_t index;
    // only check the positive or negative surface
    for (uint32_t i = 0; i < 3; i++) {
      index = 2 * i + sgn(angle[i]);
      dist = (nodes[index] - pos[i]) / angle[i];
      if (dist < min_dist) {
        min_dist = dist;
        surface_cross = index;
      }
    }
    return min_dist;
  }

  //! Set position array given an RNG
  inline void uniform_position_in_cell(RNG *rng, double *pos) const {
    pos[0] = nodes[0] + rng->generate_random_number() * (nodes[1] - nodes[0]);
    pos[1] = nodes[2] + rng->generate_random_number() * (nodes[3] - nodes[2]);
    pos[2] = nodes[4] + rng->generate_random_number() * (nodes[5] - nodes[4]);
  }

  //! Determine if position is inside a cell (diagnostic only)
  bool check_in_cell(const double *pos) const {
    bool in_cell = true;
    if (pos[0] < nodes[0] || pos[0] > nodes[1])
      in_cell = false;
    if (pos[1] < nodes[2] || pos[1] > nodes[3])
      in_cell = false;
    if (pos[2] < nodes[4] || pos[2] > nodes[5])
      in_cell = false;
    return in_cell;
  }

  //! Return node array (for setting up work packets)
  inline const double *get_node_array(void) const { return nodes.data(); }

  //! Return SILO index (for plotting only)
  inline uint32_t get_silo_index(void) const { return silo_index; }

  //! Return heat capacity
  inline double get_cV(void) const { return cV; }

  //! Retrun absorption opacity
  inline double get_op_a(void) const { return op_a; }

  //! Return multigroup absorption opacity
  inline double get_op_a(uint32_t g) const { return abs_groups[g]; }

  //! Retrun scattering opacity
  inline double get_op_s(void) const { return op_s; }

  //! Return multigroup scattering opacity
  inline double get_op_s(uint32_t g) const { return sct_groups[g]; }

  //! Retrun fleck factor
  inline double get_f(void) const { return f; }

  //! Return density
  inline double get_rho(void) const { return rho; }

  //! Return cell volume
  inline double get_volume(void) const {
    return (nodes[1] - nodes[0]) * (nodes[3] - nodes[2]) *
           (nodes[5] - nodes[4]);
  }

  //! Return electron temperature
  inline double get_T_e(void) const { return T_e; }

  //! Return radiation temperature
  inline double get_T_r(void) const { return T_r; }

  //! Return source temperature
  inline double get_T_s(void) const { return T_s; }

  // Return global ID
  inline uint32_t get_ID(void) const { return g_ID; }

  // Return global grip ID
  inline uint32_t get_grip_ID(void) const { return grip_ID; }

  // Return region ID
  inline uint32_t get_region_ID(void) const { return region_ID; }

  //! Set input array to center of cell (for mesh decomposition only)
  inline void get_center(float xyz[3]) const {
    xyz[0] = 0.5 * (nodes[0] + nodes[1]);
    xyz[1] = 0.5 * (nodes[2] + nodes[3]);
    xyz[2] = 0.5 * (nodes[4] + nodes[5]);
  }

  //! Override great than operator to sort
  bool operator<(const Cell &compare) const { return g_ID < compare.get_ID(); }

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

    cout<<"rank: "<<my_rank<<" global index: "<<g_ID<<" boundary cell: "<<boundary<<std::endl;
    cout<<nodes[0]<<" "<<nodes[1]<<" "<<nodes[2]<<" "<<nodes[3]<<" "<<nodes[4]<<" "<<nodes[5]<<endl;

    cout<<"Temperatures: "<<T_e<<" "<<T_r<<" "<<T_s<<endl;
    cout<<"Density: "<<rho<<" cV: "<<cV<<" f: "<<f<<endl;
    for(int i=0;i<BRANSON_N_GROUPS;++i)
      cout<<"group: "<<i<<"/"<<BRANSON_N_GROUPS<<" abs: "<<abs_groups[i]<<" sct: "<<sct_groups[i]<<std::endl;
  }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Provide static function for sorting based on grip ID
  static bool sort_grip_ID(const Cell &compare_1, const Cell &compare_2) {
    return compare_1.get_grip_ID() < compare_2.get_grip_ID();
  }

  //! Set neighbor in a given direction by global cell ID
  void set_neighbor(Constants::dir_type neighbor_dir, uint32_t nbr_g_ID) {
    e_next[neighbor_dir] = nbr_g_ID;
  }

  //! Set grip neighbor in a given direction by global grip ID
  void set_grip_neighbor(Constants::dir_type neighbor_dir,
                         uint32_t nbr_grip_ID) {
    grip_next[neighbor_dir] = nbr_grip_ID;
  }

  //! Set boundary conditions for cell in a given direction
  void set_bc(Constants::dir_type direction, Constants::bc_type _bc) {
    bc[direction] = _bc;
  }

  //! Set absorption opacity
  void set_op_a(double _op_a) {
    op_a = _op_a;
    // set "multigroup" data
    for (uint32_t i = 0; i < BRANSON_N_GROUPS; ++i) {
      abs_groups[i] = _op_a;
    }
  }

  //! Set scattering opacity
  void set_op_s(double _op_s) {
    op_s = _op_s;
    // set "multigroup" data
    for (uint32_t i = 0; i < BRANSON_N_GROUPS; ++i) {
      sct_groups[i] = _op_s;
    }
  }

  //! Set fleck factor
  void set_f(double _f) { f = _f; }

  //! Set heat capacity
  void set_cV(double _cV) { cV = _cV; }

  //! Set electron temperature
  void set_T_e(double _T_e) { T_e = _T_e; }

  //! Set density
  void set_rho(double _rho) { rho = _rho; }

  //! Set radiation temperature
  void set_T_r(double _T_r) { T_r = _T_r; }

  //! Set source temperature
  void set_T_s(double _T_s) { T_s = _T_s; }

  //! Set global ID
  void set_ID(uint32_t _id) { g_ID = _id; }

  //! Set global grip ID
  void set_grip_ID(uint32_t _grip_id) { grip_ID = _grip_id; }

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

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  uint32_t g_ID; //!< Global ID, valid across all ranks

  //! Global ID of cell at the center of grip, valid across all ranks
  uint32_t grip_ID;

  uint32_t region_ID; //!< region cell is in (for setting physical properties)
  uint32_t silo_index;      //!< Global index not remappaed, for SILO plotting
  std::array<uint32_t, 6> e_next; //!< Bordering cell, given as global ID
  std::array<uint32_t, 6> grip_next;    //!< Bordering grip, given as global cell ID
  std::array<Constants::bc_type, 6>  bc; //!< Boundary conditions for each face
  std::array<double, 6> nodes;          //!< x_low, x_high, y_low, y_high, z_low, z_high
  std::array<double, BRANSON_N_GROUPS> abs_groups; //!< Absorption groups
  std::array<double, BRANSON_N_GROUPS> sct_groups; //!< Scattering groups

  double cV;   //!< Heat capacity  GJ/g/KeV
  double op_a; //!< Absorption opacity  (1/cm)
  double op_s; //!< Physical scattering opacity (1/cm)
  double f;    //!< Fleck factor
  double rho;  //!< Density (g/cc)
  double T_e;  //!< Material temperature
  double T_r;  //!< Radiation temperature
  double T_s;  //!< Source temperature
};

#endif // cell_h_
//---------------------------------------------------------------------------//
// end of cell.h
//---------------------------------------------------------------------------//
