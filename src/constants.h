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
constexpr double pi(3.1415926535897932384626433832795); //!< Pi
constexpr double c(299.792458); //!< speed of light in cm/shake
constexpr double c_SO(1.0);     //!< speed of light for Su-Olson problem
constexpr double h(6.62606957e-34 * 1.0e-9 /
               1.0e-8);      //!< Planck's constant in GJ/sh
constexpr double k(1.60219e-31); //!< energy conversion constant GJ/keV
constexpr double a(0.01372);     //!< Boltzmann constant in GJ/cm^3/keV^4
constexpr double a_SO(1.0);      //!< Boltzmann constant for SO problems

enum bc_type { REFLECT, VACUUM, ELEMENT, SOURCE, PROCESSOR }; //!< Boundary conditions
enum dir_type { X_NEG, X_POS, Y_NEG, Y_POS, Z_NEG, Z_POS }; //!< Directions
enum event_type { KILL, EXIT, PASS, CENSUS, WAIT };         //!< Events
enum {
  PARTICLE_PASS,
  REPLICATED
};                                 //!< Parallel types
enum { NO_DECOMP, METIS, CUBE };   //!< Mesh decomposition method
constexpr int grip_id_tag(1);          //!< MPI tag for grip ID messages
constexpr int cell_id_tag(2);          //!< MPI tag for requested cell ID messages
constexpr int count_tag(3);            //!< MPI tag for completion count messages
constexpr int photon_tag(4);           //!< MPI tag for photon messages
constexpr int n_photon_tag(5);         //!< MPI tag for work packet messages
constexpr int work_tag(6);             //!< MPI tag for work packet messages
constexpr int n_work_tag(7);           //!< MPI tag for work packet messages
constexpr int tally_tag(8);            //!< MPI tag for tally messages
constexpr int n_tally_tag(9);          //!< MPI tag for number of tally messages

//! MPI tag for cell messages NOTE: the number of grips in the message will
// added to the end of this number
constexpr int cell_tag(10);

constexpr int n_threads_per_block = 512;
}; // namespace Constants

#endif // constants_h_
//---------------------------------------------------------------------------//
// end of constants.h
//---------------------------------------------------------------------------//
