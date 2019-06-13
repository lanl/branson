//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc_parameters.h
 * \author Alex Long
 * \date   December 3 2015
 * \brief  Holds parameters needed in IMC simulation
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
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

class IMC_Parameters {
public:
  //! constructor
  IMC_Parameters(const Input &input)
      : n_user_photon(input.get_number_photons()),
        write_silo_flag(input.get_write_silo_bool()),
        output_frequency(input.get_output_freq()) {}

  //! destructor
  ~IMC_Parameters() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Return total photons specified by the user
  uint64_t get_n_user_photon(void) const { return n_user_photon; }

  //! Get SILO write flag
  bool get_write_silo_flag(void) const { return write_silo_flag; }

  //! Get SILO output frequency
  int32_t get_output_frequency(void) const { return output_frequency; }

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  uint64_t n_user_photon;   //!< User requested number of photons per timestep
  int32_t output_frequency; //!< How often to write silo (every n cycles)
  bool write_silo_flag;     //!< Write SILO output files flag
};

#endif // imc_parameters_h_
//----------------------------------------------------------------------------//
// end of imc_parameters.h
//----------------------------------------------------------------------------//
