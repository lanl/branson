/*
  Author: Alex Long
  Date: 12/3/2015
  Name: imc_parameters.h
*/
#ifndef imc_parameters_h_
#define imc_parameters_h_

#include "input.h"

class IMC_Parameters
{
  public:
  IMC_Parameters(Input *input)
    : n_user_photon(input->get_number_photons()),
      map_size(input->get_map_size()),
      dd_mode(input->get_dd_mode()) ,
      use_ghost_cells(input->get_ghost_cell_bool()),
      batch_size(input->get_batch_size()),
      particle_message_size(input->get_particle_message_size())
    {}
  ~IMC_Parameters() {}

/*****************************************************************************/
/* const functions                                                           */
/*****************************************************************************/
  uint64_t get_n_user_photon(void) const {return n_user_photon;}
  uint32_t get_map_size(void) const {return map_size;}
  uint32_t get_dd_mode(void) const {return dd_mode;}
  bool get_ghost_cell_bool(void) const {return use_ghost_cells;}
  uint32_t  get_batch_size(void) const {return batch_size;}
  uint32_t  get_particle_message_size(void) const {return particle_message_size;}

/*****************************************************************************/
/* member variables and private functions                                    */
/*****************************************************************************/
  private:

  //parallel performance parameters
  uint64_t n_user_photon; //!< User requested number of photons per timestep
  uint32_t map_size; //!< Size of stored off-rank mesh cells
  uint32_t dd_mode; //!< Mode of domain decomposed transport algorithm
  bool use_ghost_cells; //!< Always keep first ghost cells
  uint32_t batch_size; //!< How often to check for MPI passed data
  uint32_t particle_message_size; //!< Preferred number of particles in MPI sends
};

#endif // #ifdef imc_parameters_h_
