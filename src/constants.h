//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   constants.h
 * \author Alex Long
 * \date   July 18 2014
 * \brief  Physical constants, custom enumerations and MPI tags
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef constants_h_
#define constants_h_

namespace Constants {

constexpr int LOG_CYCLES = 20;
constexpr int EXP_CYCLES = 20;
constexpr int DIVIDE_CYCLES = 10;
constexpr int RNG_CYCLES = 5;
constexpr int PARTICLE_COMM_CYCLES = 10; // send 10 words for particle data

constexpr double pi(3.1415926535897932384626433832795); //!< Pi
constexpr double c(299.792458); //!< speed of light in cm/shake
constexpr double c_SO(1.0);     //!< speed of light for Su-Olson problem
constexpr double h(6.62606957e-34 * 1.0e-9 /
               1.0e-8);      //!< Planck's constant in GJ/sh
constexpr double k(1.60219e-31); //!< energy conversion constant GJ/keV
constexpr double a(0.01372);     //!< Boltzmann constant in GJ/cm^3/keV^4
constexpr double a_SO(1.0);      //!< Boltzmann constant for SO problems

enum bc_type { REFLECT, VACUUM, ELEMENT, SOURCE }; //!< Boundary conditions
enum dir_type { X_NEG, X_POS, Y_NEG, Y_POS, Z_NEG, Z_POS }; //!< Directions
enum event_type { KILL, EXIT, CENSUS, WAIT, BOUNDARY };               //!< Events
enum { EMISSION, INITIAL_CENSUS }; //!< Particle type for work packets

}; // namespace Constants

#endif // constants_h_
//---------------------------------------------------------------------------//
// end of constants.h
//---------------------------------------------------------------------------//
