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
const double pi(3.1415926535897932384626433832795); //!< Pi
const double c(299.792458); //!< speed of light in cm/shake
const double c_SO(1.0);     //!< speed of light for Su-Olson problem
const double h(6.62606957e-34 * 1.0e-9 /
               1.0e-8);      //!< Planck's constant in GJ/sh
const double k(1.60219e-31); //!< energy conversion constant GJ/keV
const double a(0.01372);     //!< Boltzmann constant in GJ/cm^3/keV^4
const double a_SO(1.0);      //!< Boltzmann constant for SO problems

enum bc_type { REFLECT, VACUUM, ELEMENT }; //!< Boundary conditions
enum dir_type { X_NEG, X_POS, Y_NEG, Y_POS, Z_NEG, Z_POS }; //!< Directions
enum event_type { KILL, EXIT, CENSUS, WAIT };               //!< Events
enum { EMISSION, INITIAL_CENSUS }; //!< Particle type for work packets

}; // namespace Constants

#endif // constants_h_
//---------------------------------------------------------------------------//
// end of constants.h
//---------------------------------------------------------------------------//
