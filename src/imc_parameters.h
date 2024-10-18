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
      : n_user_photons(input.get_number_photons()),
        seed(input.get_rng_seed()),
        dd_mode(input.get_dd_mode()), batch_size(input.get_batch_size()),
        particle_message_size(input.get_particle_message_size()),
        output_frequency(input.get_output_freq()),
        transport_algorithm(input.get_particle_algorithm()),
        n_omp_threads(input.get_n_omp_threads()),
        write_silo_flag(input.get_write_silo_bool()),
        use_gpu_transporter_flag(input.get_use_gpu_transporter_bool()),
        use_comb_flag(input.get_comb_bool()) {}

  //! destructor
  ~IMC_Parameters() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Return total photons specified by the user
  uint64_t get_n_user_photons() const { return n_user_photons; }

  //! Return the user-set RNG seed
  uint32_t get_rng_seed() const {return seed;}

  //! Return domain decomposition algorithm
  uint32_t get_dd_mode() const { return dd_mode; }

  //! Get the number of particles to run between MPI message processing
  uint32_t get_batch_size() const { return batch_size; }

  //! Get the desired number of particles in messages (particle passing only)
  uint32_t get_particle_message_size() const {
    return particle_message_size;
  }

  //! Get SILO write flag
  bool get_write_silo_flag() const { return write_silo_flag; }

  //! Get the GPU transporter flag
  bool get_use_gpu_transporter_flag() const {return use_gpu_transporter_flag;}

  //! Get the combing flag
  bool get_use_comb_flag() const {return use_comb_flag;}

  //! Get output frequency (print when cycle % frequency == 0)
  uint32_t get_output_frequency() const { return output_frequency; }

  //! Get number of OpenMP threads to use (set by user in input)
  uint32_t get_n_omp_threads() const { return n_omp_threads; }

  //! Get user-input transport algorithm (e.g., event or history)
  uint32_t get_transport_algorithm() const { return transport_algorithm; }

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  uint64_t n_user_photons; //!< User requested number of photons per timestep
  uint32_t seed;       //!< Random number seed
  uint32_t dd_mode;    //!< Mode of domain decomposed transport algorithm
  uint32_t batch_size; //!< How often to check for MPI passed data
  uint32_t
      particle_message_size; //!< Preferred number of particles in MPI sends
  uint32_t output_frequency; //!< Frequency to dump output files
  uint32_t transport_algorithm; //!< How to process transport (e.g., event or history) 
  uint32_t n_omp_threads; //!< Number of OpenMP threads, set by user
  bool write_silo_flag;      //!< Write SILO output files flag
  bool use_gpu_transporter_flag;      //!< Write SILO output files flag
  bool use_comb_flag;                 //!< Comb the census if great than  n_user_photon after cycle 
};

#endif // imc_parameters_h_
//----------------------------------------------------------------------------//
// end of imc_parameters.h
//----------------------------------------------------------------------------//
