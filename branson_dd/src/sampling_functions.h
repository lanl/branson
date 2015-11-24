/*
  Author: Alex Long
  Date: 9/17/2014
  Name: sampling_functions.h
*/

#ifndef sampling_functions_h_
#define sampling_functions_h_

#include <stdlib.h>

#include "RNG.h"
#include "constants.h"


double get_uniform_angle(RNG* rng) {
  return rng->generate_random_number()*2.0-1.0;
}

double get_source_angle(RNG* rng) {
  return sqrt(rng->generate_random_number());
}

void get_uniform_angle(double* angle, RNG* rng) {
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

void get_stratified_angle(double* angle, RNG* rng, unsigned int isample, unsigned int nsample) {
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

void get_source_angle(double* angle, RNG* rng) {
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

#endif
