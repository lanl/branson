/*
  Author: Alex Long
  Date: 3/16/2015
  Name: element.h
*/

#ifndef element_h_
#define element_h_

#include <vector>

#include "RNG.h"
#include "constants.h"

using std::vector;
using Constants::dir_type;
using Constants::bc_type;
using Constants::large_distance;
using Constants::normals;


template <typename T> int sgn(T val) {
  return (T(0) < val);
}

class Element
{

  public:
  
  Element(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high) {
    nodes[0] = x_low;
    nodes[1] = x_high;
    nodes[2] = y_low;
    nodes[3] = y_high;
    nodes[4] = z_low;
    nodes[5] = z_high;

    temp_abs_E = 0.0;
    op = 0.0;
    f = 0.0;
  }

  ~Element(void) {}

  // const functions
  bc_type get_bc(const unsigned int& dir) const {return bc[dir];}
  unsigned int get_next_element(const unsigned int& dir) const {return e_next[dir];} 

  double get_distance_to_boundary(const double *pos, const double *angle, unsigned int& surface_cross) const {

    double min_dist = 1.0e16;
    double dist = 0.0;
    unsigned int index;
    //only check the positive or negative surface
    for (unsigned int i = 0; i<3; i++) {
      index = 2*i + sgn(angle[i]);
      dist = (nodes[index] - pos[i])/angle[i];
      if (dist < min_dist) {
        min_dist = dist;
        surface_cross = index;
      }
    }
    return min_dist;
  }

  void uniform_position_in_elem(RNG* rng, double* pos) const {
    pos[0]= nodes[0] + rng->generate_random_number()*(nodes[1]-nodes[0]);
    pos[1]= nodes[2] + rng->generate_random_number()*(nodes[3]-nodes[2]);
    pos[2]= nodes[4] + rng->generate_random_number()*(nodes[5]-nodes[4]);
  }

  bool check_in_element(const double * pos) const {
    bool in_element = true;
    if (pos[0] < nodes[0] || pos[0] > nodes[1]) in_element = false;
    if (pos[1] < nodes[2] || pos[1] > nodes[3]) in_element = false;
    if (pos[2] < nodes[4] || pos[2] > nodes[5]) in_element = false;
    return in_element;
  }

  double get_op(void) const {return op;}
  double get_f(void) const {return f;}
  double get_volume(void) const {return (nodes[1]-nodes[0])*(nodes[3]-nodes[2])*(nodes[5]-nodes[4]);}

  //non const
  void set_neighbor(dir_type neighbor_dir, unsigned int index) { e_next[neighbor_dir] = index;}
  void set_bc(dir_type direction, bc_type _bc) { bc[direction] = _bc; }


  void set_op(double _op) {op = _op;}
  void set_f(double _f) {f = _f;}

  private:
  double nodes[6];  //!< x_low, x_high, y_low, y_high, z_low, z_high
  bc_type bc[6];
  unsigned int e_next[6];
  
  double temp_abs_E;
  double op;
  double f;
};

#endif // element_h_
