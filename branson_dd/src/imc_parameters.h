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
    : map_size(input->get_map_size()),
      dd_mode(input->get_dd_mode()) ,
      use_ghost_cells(input->get_ghost_cell_bool()),
      check_frequency(input->get_check_frequency())
    {}
  ~IMC_Parameters() {}

/*****************************************************************************/
/* const functions                                                           */
/*****************************************************************************/
  unsigned int get_map_size(void) const {return map_size;}
  unsigned int get_dd_mode(void) const {return dd_mode;}
  bool get_ghost_cell_bool(void) const {return use_ghost_cells;}
  int  get_check_frequency(void) {return check_frequency;}

/*****************************************************************************/
/* member variables and private functions                                    */
/*****************************************************************************/
  private:

  //parallel performance parameters
  unsigned int map_size; //!< Size of stored off-rank mesh cells
  unsigned int dd_mode; //!< Mode of domain decomposed transport algorithm
  bool use_ghost_cells; //!< Always keep first ghost cells
  int check_frequency; //!< How often to check for MPI passed data
};

#endif // #ifdef imc_parameters_h_
