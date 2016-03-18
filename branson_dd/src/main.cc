
#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>
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

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

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
  Mesh *mesh = new Mesh(input, rank, n_rank);

  //timing 
  struct timeval start,end;
  struct timezone tzp;
  gettimeofday(&start, &tzp); 

  // decompose mesh with ParMETIS
  decompose_mesh(mesh);

  MPI_Barrier(MPI_COMM_WORLD);
  //print_MPI_out(mesh, rank, n_rank);

  /****************************************************************************/ 
  // TRT PHYSICS CALCULATION
  /****************************************************************************/

  if (input->get_dd_mode() == PARTICLE_PASS)
    imc_particle_pass_driver(rank, mesh, imc_state, imc_p);
  else if (input->get_dd_mode() == CELL_PASS)
    imc_cell_pass_driver(rank, mesh, imc_state, imc_p);
  
  if (rank==0) {
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
    gettimeofday(&end, &tzp);
    print_elapsed_inside("runtime:",&start, &end);
    imc_state->print_simulation_footer(input->get_dd_mode());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  delete mesh;
  delete imc_state;
  delete imc_p;

  MPI_Finalize();
}
