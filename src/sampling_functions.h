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

#include <stdlib.h>

#include "RNG.h"
#include "constants.h"

//! Set angle given input array and RNG
inline void get_uniform_angle(double* angle, RNG* rng) {
  using std::sqrt;
  using std::sin;
  using std::cos;
  using Constants::pi;
  double mu =rng->generate_random_number()*2.0-1.0; 
  double phi = rng->generate_random_number()*2.0*pi;
  double sin_theta = sqrt(1.0 - mu*mu);
  angle[0] = sin_theta*cos(phi);
  angle[1] = sin_theta*sin(phi);
  angle[2] = mu;
}

//! Set angle given input array, RNG and strata
inline void get_stratified_angle(double* angle, RNG* rng, uint32_t isample, uint32_t nsample) {
  using std::sqrt;
  using std::sin;
  using std::cos;
  using Constants::pi;
  //stratify by octant--two polar, four azimuthal
  double frac =double(isample)/nsample;
  int imu = int(frac > 0.5) ; // 0 or 1
  int iphi = int(frac*4.0); // 0 through 3
  double mu = 0.5*(imu + rng->generate_random_number())*2.0-1.0;
  double phi = 0.25*(iphi + rng->generate_random_number())*2.0*pi;
  double sin_theta = sqrt(1.0 - mu*mu);
  angle[0] = sin_theta*cos(phi);
  angle[1] = sin_theta*sin(phi);
  angle[2] = mu;
}

//! Set angle from face source given input array, RNG and strata
inline void get_source_angle(double* angle, RNG* rng) {
  using std::sqrt;
  using std::sin;
  using std::cos;
  using Constants::pi;
  double mu =sqrt(rng->generate_random_number());
  double phi = rng->generate_random_number()*2.0*pi;
  double sin_theta = sqrt(1.0 - mu*mu);
  angle[0] = sin_theta*cos(phi);
  angle[1] = sin_theta*sin(phi);
  angle[2] = mu;
}


//! Sample the group after an effective scattering event 
inline int sample_emission_group(RNG* rng, const Cell& cell_data) {
  // Sample a new group from a uniform CDF (but mimc non-uniform CDF algorithm)
  double cdf_value = rng->generate_random_number();
  int new_group = -1;
  // normalizes the PDF (opacity is uniform, not weighting with spectrum so
  // this is very simple)
  double norm_factor = 1.0/(cell_data.get_op_a(0)*BRANSON_N_GROUPS);
  while(cdf_value > 0) {
    new_group++;
    cdf_value -= cell_data.get_op_a(new_group) * norm_factor;
  }
  return new_group;
}



#endif
//---------------------------------------------------------------------------//
// end of sampling_functions.h
//---------------------------------------------------------------------------//
