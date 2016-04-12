/*
  Author: Alex Long
  Date: 3/16/2015
  Name: cell.h
*/

#ifndef cell_h_
#define cell_h_

#include <iostream>
#include <mpi.h>

#include "RNG.h"
#include "constants.h"

template <typename T> int sgn(T val) {
  return (T(0) < val);
}

//==============================================================================
/*!
 * \class Cell
 * \brief Basic geometry unit, holds physical data that is read only during 
 * transport.
 *
 * A cartesian mesh cell. Holds location of each node, boundary information for
 * each face, opacity data and temperature data. The temperature data probably
 * stored in an cell because it's not essential information for transport.
 */
//==============================================================================
class Cell
{

  public:
  
  Cell(void) {
    op_a = 0.0;
    op_s = 0.0;
    f = 0.0;
  }

  ~Cell(void) {}

/*****************************************************************************/
/* const functions                                                           */
/*****************************************************************************/
  Constants::bc_type get_bc(const uint32_t& dir) const {return bc[dir];}

  uint32_t get_next_cell(const uint32_t& dir) const 
  {
    return e_next[dir];
  } 

  double get_distance_to_boundary(const double *pos, 
                                  const double *angle, 
                                  uint32_t& surface_cross) const 
  {
    double min_dist = 1.0e16;
    double dist = 0.0;
    uint32_t index;
    //only check the positive or negative surface
    for (uint32_t i = 0; i<3; i++) {
      index = 2*i + sgn(angle[i]);
      dist = (nodes[index] - pos[i])/angle[i];
      if (dist < min_dist) {
        min_dist = dist;
        surface_cross = index;
      }
    }
    return min_dist;
  }

  void uniform_position_in_cell(RNG* rng, double* pos) const {
    pos[0]= nodes[0] + rng->generate_random_number()*(nodes[1]-nodes[0]);
    pos[1]= nodes[2] + rng->generate_random_number()*(nodes[3]-nodes[2]);
    pos[2]= nodes[4] + rng->generate_random_number()*(nodes[5]-nodes[4]);
  }

  bool check_in_cell(const double * pos) const {
    bool in_cell = true;
    if (pos[0] < nodes[0] || pos[0] > nodes[1]) in_cell = false;
    if (pos[1] < nodes[2] || pos[1] > nodes[3]) in_cell = false;
    if (pos[2] < nodes[4] || pos[2] > nodes[5]) in_cell = false;
    return in_cell;
  }

  double get_x_low(void) const {return nodes[0];}
  double get_x_high(void) const {return nodes[1];}
  double get_y_low(void) const {return nodes[2];}
  double get_y_high(void) const {return nodes[3];}
  double get_z_low(void) const {return nodes[4];}
  double get_z_high(void) const {return nodes[5];}

  double get_cV(void) const {return cV;}
  double get_op_a(void) const {return op_a;}
  double get_op_s(void) const {return op_s;}
  double get_f(void) const {return f;}
  double get_rho(void) const {return rho;}
  double get_volume(void) const 
  {
    return (nodes[1]-nodes[0])*(nodes[3]-nodes[2])*(nodes[5]-nodes[4]);
  }
  double get_T_e(void) const {return T_e;}
  double get_T_r(void) const {return T_r;}
  double get_T_s(void) const {return T_s;}
  uint32_t get_ID(void) const {return g_ID;}
  uint32_t get_region_ID(void) const {return region_ID;}

  //override great than operator to sort
  bool operator <(const Cell& compare) const {
    return g_ID < compare.get_ID();
  }
  
  void print(void) const {
    using Constants::PROCESSOR;
    using std::cout;
    using std::endl;
    int32_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    bool boundary = false;
    for (uint32_t i=0;i<6;i++) {
      if (bc[i] == PROCESSOR) boundary = true;
    }
    
    //cout<<g_ID<<" "<<boundary;
    //cout<<nodes[0]<<" "<<nodes[2]<<" "<<nodes[4]<<endl;
    
    cout<<"Rank: "<<my_rank<<" Global ID: "<<g_ID<<endl;
    cout<<nodes[0]<<" "<<nodes[2]<<" "<<nodes[4];
    cout<<" Processor bound: "<<boundary<<endl;
    //cout<<"Temperatures: "<<T_e<<" "<<T_r<<" "<<T_s<<endl;
    //cout<<"Density: "<<rho<<" cV: "<<cV<<" f: "<<f<<endl;
   
  }

/*****************************************************************************/
/* non-const functions (set)                                                 */
/*****************************************************************************/
  void set_neighbor(Constants::dir_type neighbor_dir, uint32_t index) {
    e_next[neighbor_dir] = index;
  }

  void set_bc(Constants::dir_type direction, Constants::bc_type _bc) {
    bc[direction] = _bc; 
  }

  void set_op_a(double _op_a) {op_a = _op_a;}
  void set_op_s(double _op_s) {op_s = _op_s;}
  void set_f(double _f) {f = _f;}
  void set_cV(double _cV) {cV = _cV;}
  void set_T_e(double _T_e) {T_e = _T_e;}
  void set_rho(double _rho) {rho = _rho;}
  void set_T_r(double _T_r) {T_r = _T_r;}
  void set_T_s(double _T_s) {T_s = _T_s;}
  void set_ID(double _id) {g_ID = _id;}
  void set_region_ID(uint32_t _region_ID) { region_ID=_region_ID;}
  void set_coor(double x_low, double x_high, double y_low, 
                double y_high, double z_low, double z_high) 
  {
    nodes[0] = x_low;
    nodes[1] = x_high;
    nodes[2] = y_low;
    nodes[3] = y_high;
    nodes[4] = z_low;
    nodes[5] = z_high;
  }


/*****************************************************************************/
/* member variables and private functions                                    */
/*****************************************************************************/
  private:
  uint32_t g_ID; //! Global ID, valid across all ranks
  uint32_t region_ID; //! region cell is in (for setting physical properties)
  uint32_t e_next[6]; //! Bordering cell, given as global ID
  Constants::bc_type bc[6];   //! Boundary conditions for each face 
  double nodes[6]; //! x_low, x_high, y_low, y_high, z_low, z_high
  
  double cV;    //! Heat capacity  GJ/g/KeV
  double op_a;  //! Absorption opacity  (1/cm)
  double op_s;  //! Physical scattering opacity (1/cm)
  double f;     //! Fleck factor
  double rho;   //! Density (g/cc)
  double T_e;   //! Material temperature
  double T_r;   //! Radiation temperature
  double T_s;   //! Source temperature
};

#endif // cell_h_
