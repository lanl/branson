//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cell_tally.h
 * \author Alex Long
 * \date   March 3 2015
 * \brief  Holds cells data and provides basic sampling functions
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef cell_tally_h_
#define cell_tally_h_

#include <iostream>
#include <mpi.h>

#include "config.h"

//==============================================================================
/*!
 * \class Cell_Tally
 * \brief  Simple class that holds tally data needed in each cell on the mesh
 *
 * This cell is mostly to abstract away GPU atomics
 */
//==============================================================================

class Cell_Tally {

public:
  Cell_Tally()
    :abs_E{0.0}, track_E{0.0}
  {}

  GPU_HOST_DEVICE
  inline void accumulate_absorbed_E(const double delta_abs_E) {
    accumulate(abs_E, delta_abs_E);
  }

  GPU_HOST_DEVICE
  inline void accumulate_track_E(const double delta_track_E) {
    accumulate(track_E, delta_track_E);
  }

  double get_abs_E() const {return abs_E;}
  double get_track_E() const {return track_E;}

  double abs_E;  //!< Absorbed energy in jerks
  double track_E;  //!< Track energy used for estimate of radiation temperature
};

#endif // cell_tally_h_
//---------------------------------------------------------------------------//
// end of cell_tally.h
//---------------------------------------------------------------------------//
