//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc_parameters.h
 * \author Alex Long
 * \date   December 3 2016
 * \brief  Holds parameters needed in IMC simulation
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef imc_parameters_h_
#define imc_parameters_h_

#include "input.h"

//==============================================================================
/*!
 * \class IMC_Parameters
 * \brief Holds parameters used in IMC simulation
 * 
 * Initialized with the input class and then data members are invariant
 * \example no test yet
 */
//==============================================================================

class IMC_Parameters
{
  public:
  //! constructor
  IMC_Parameters(Input *input)
    : n_user_photon(input->get_number_photons()),
      map_size(input->get_map_size()),
      dd_mode(input->get_dd_mode()) ,
      use_ghost_cells(input->get_ghost_cell_bool()),
      batch_size(input->get_batch_size()),
      particle_message_size(input->get_particle_message_size())
    {}

  // destructor
  ~IMC_Parameters() {}

  //////////////////////////////////////////////////////////////////////////////
  // const functions                                                          //
  //////////////////////////////////////////////////////////////////////////////

  //! Return total photons specified by the user
  uint64_t get_n_user_photon(void) const {return n_user_photon;}

  //! Return maximum size of stored remote mesh
  uint32_t get_map_size(void) const {return map_size;}

  //! Return domain decomposition algorithm
  uint32_t get_dd_mode(void) const {return dd_mode;}

  //! Return the ghost cell flag used in mesh passing
  bool get_ghost_cell_bool(void) const {return use_ghost_cells;}

  //! Get the number of particles to run between MPI message processing
  uint32_t get_batch_size(void) const {return batch_size;}

  //! Get the desired number of particles in messages (particle passing only)
  uint32_t get_particle_message_size(void) const {return particle_message_size;}

  //////////////////////////////////////////////////////////////////////////////
  // member data                                                              //
  //////////////////////////////////////////////////////////////////////////////
  private:

  uint64_t n_user_photon; //!< User requested number of photons per timestep
  uint32_t map_size; //!< Size of stored off-rank mesh cells
  uint32_t dd_mode; //!< Mode of domain decomposed transport algorithm
  bool use_ghost_cells; //!< Always keep first ghost cells
  uint32_t batch_size; //!< How often to check for MPI passed data
  uint32_t particle_message_size; //!< Preferred number of particles in MPI sends
};

#endif // imc_parameters_h_

//---------------------------------------------------------------------------//
// end of imc_parameters.h
//---------------------------------------------------------------------------//
