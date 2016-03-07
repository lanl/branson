
#include <iostream>
#include <vector>
#include <string>
#include <boost/mpi.hpp>
#include <time.h>
#include <sys/time.h>

#include "constants.h"
#include "decompose_mesh.h"
#include "imc_drivers.h"
#include "imc_state.h"
#include "imc_parameters.h"
#include "input.h"
#include "mesh.h"
#include "timing_functions.h"

using std::vector;
using std::endl;
using std::cout;
using std::string;
using Constants::PARTICLE_PASS;
using Constants::CELL_PASS;

namespace mpi = boost::mpi;

int main(int argc, char *argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;
  uint32_t rank = world.rank();
  uint32_t size = world.size();

  //check to see if number of arguments is correct
  if (argc != 2) {
    cout<<"Usage: BRANSON <path_to_input_file>"<<endl;
    exit(EXIT_FAILURE); 
  }

  //get input object from filename
  std::string filename( argv[1]);
  Input *input;
  input = new Input(filename);
  if(rank ==0) input->print_problem_info();

  //IMC paramters setup
  IMC_Parameters *imc_p;
  imc_p = new IMC_Parameters(input);

  //IMC state setup
  IMC_State *imc_state;
  imc_state = new IMC_State(input, rank);

  // make mesh from input object
  Mesh *mesh = new Mesh(input, rank, size);

  //timing 
  struct timeval start,end;
  struct timezone tzp;
  gettimeofday(&start, &tzp); 

  // decompose mesh with ParMETIS
  decompose_mesh(mesh);

  world.barrier();
  //print_MPI_out(mesh, rank, size);


  /****************************************************************************/ 
  // TRT PHYSICS CALCULATION
  /****************************************************************************/

  if (input->get_dd_mode() == PARTICLE_PASS)
    imc_particle_pass_driver(rank, mesh, imc_state, imc_p, world);
  else if (input->get_dd_mode() == CELL_PASS)
    imc_cell_pass_driver(rank, mesh, imc_state, imc_p, world);
  
  if (rank==0) {
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
    gettimeofday(&end, &tzp);
    print_elapsed_inside("runtime:",&start, &end);
    imc_state->print_simulation_footer(input->get_dd_mode());
  }
  world.barrier();

  delete mesh;
  delete imc_state;
  delete imc_p;
}
