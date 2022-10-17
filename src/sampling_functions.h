//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sampling_functions.h
 * \author Alex Long
 * \date   September 17 2014
 * \brief  Angle sampling for isotropic and surface sources
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef sampling_functions_h_
#define sampling_functions_h_

#include <cmath>

#include "RNG.h"
#include "constants.h"
#include "cell.h"

//! Set an input array to a random position within a cell
inline std::array<double, 3> get_uniform_position_in_cell(const Cell &cell, RNG *rng)  {
  auto nodes = cell.get_node_array();
  std::array<double, 3> pos{0.0, 0.0, 0.0};
  pos[0] = nodes[0] + rng->generate_random_number() * (nodes[1] - nodes[0]);
  pos[1] = nodes[2] + rng->generate_random_number() * (nodes[3] - nodes[2]);
  pos[2] = nodes[4] + rng->generate_random_number() * (nodes[5] - nodes[4]);
  return pos;
}

//! Set an input array to a random position within a cell
inline std::array<double, 3>  get_uniform_position_on_face(const Cell &cell, RNG *rng, int face) {
  auto nodes = cell.get_node_array();
  std::array<double, 3> face_pos{0.0, 0.0, 0.0};
  if (face ==0 || face ==1) {
    face_pos[0] = (face == 0) ? nodes[0] : nodes[1];
    face_pos[1] = nodes[2] + rng->generate_random_number() * (nodes[3] - nodes[2]);
    face_pos[2] = nodes[4] + rng->generate_random_number() * (nodes[5] - nodes[4]);
  }
  else if (face ==2 || face ==3) {
    face_pos[0] = nodes[0] + rng->generate_random_number() * (nodes[1] - nodes[0]);
    face_pos[1] = (face == 2) ? nodes[2] : nodes[3];
    face_pos[2] = nodes[4] + rng->generate_random_number() * (nodes[5] - nodes[4]);
  }
  else // face == 4 || face ==5)
  {
    face_pos[0] = nodes[0] + rng->generate_random_number() * (nodes[1] - nodes[0]);
    face_pos[1] = nodes[2] + rng->generate_random_number() * (nodes[3] - nodes[2]);
    face_pos[2] = (face == 4) ? nodes[4] : nodes[5];
  }
  return face_pos;
}


//! Set angle given input array and RNG
inline std::array<double, 3> get_uniform_angle(RNG *rng) {
  using Constants::pi;
  using std::cos;
  using std::sin;
  using std::sqrt;
  std::array<double, 3> angle{0.0, 0.0, 0.0};
  double mu = rng->generate_random_number() * 2.0 - 1.0;
  double phi = rng->generate_random_number() * 2.0 * pi;
  double sin_theta = sqrt(1.0 - mu * mu);
  angle[0] = sin_theta * cos(phi);
  angle[1] = sin_theta * sin(phi);
  angle[2] = mu;
  return angle;
}

//! Set angle given input array, RNG and strata
inline std::array<double,3> get_stratified_angle( RNG *rng, uint32_t isample,
                                 uint32_t nsample) {
  using Constants::pi;
  using std::cos;
  using std::sin;
  using std::sqrt;
  std::array<double, 3> angle{0.0, 0.0, 0.0};
  //stratify by octant--two polar, four azimuthal
  double frac = double(isample) / nsample;
  int imu = int(frac > 0.5);  // 0 or 1
  int iphi = int(frac * 4.0); // 0 through 3
  double mu = 0.5 * (imu + rng->generate_random_number()) * 2.0 - 1.0;
  double phi = 0.25 * (iphi + rng->generate_random_number()) * 2.0 * pi;
  double sin_theta = sqrt(1.0 - mu * mu);
  angle[0] = sin_theta * cos(phi);
  angle[1] = sin_theta * sin(phi);
  angle[2] = mu;
  return angle;
}

//! Set angle on face given input array and RNG
inline std::array<double, 3> get_source_angle_on_face( RNG *rng, int face) {
  using Constants::pi;
  using std::cos;
  using std::sin;
  using std::sqrt;

  std::array<double, 3> angle{0.0, 0.0, 0.0};
  double theta = acos(sqrt(rng->generate_random_number()));
  double phi = rng->generate_random_number() * 2.0 * pi;
  double sign = (face % 2) ? -1.0 : 1.0;
  if( face == 0 || face ==1) {
    angle[0] = cos(theta) * sign;
    angle[1] = sin(theta) * sin(phi);
    angle[2] = sin(theta) * cos(phi);
  }
  else if( face == 2 || face ==3) {
    angle[0] = sin(theta) * sin(phi);
    angle[1] = cos(theta);
    angle[2] = sin(theta) * cos(phi);
  }
  else // face == 4 || face ==5)
  {
    angle[0] = sin(theta) * cos(phi);
    angle[1] = sin(theta) * sin(phi);
    angle[2] = cos(theta);
  }
  return angle;
}


//! Sample the group after an effective scattering event
inline int sample_emission_group(RNG *rng, const Cell &cell_data) {
  // Sample a new group from a uniform CDF (but mimc non-uniform CDF algorithm)
  double cdf_value = rng->generate_random_number();
  int new_group = -1;
  // normalizes the PDF (opacity is uniform, not weighting with spectrum so
  // this is very simple)
  double norm_factor = 1.0 / (cell_data.get_op_a(0) * BRANSON_N_GROUPS);
  while (cdf_value > 0) {
    new_group++;
    cdf_value -= cell_data.get_op_a(new_group) * norm_factor;
  }
  return new_group;
}

#endif
//---------------------------------------------------------------------------//
// end of sampling_functions.h
//---------------------------------------------------------------------------//
