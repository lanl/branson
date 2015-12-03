
#include <iostream>
#include <vector>
#include <string>
#include <boost/mpi.hpp>
#include <time.h>
#include <sys/time.h>

#include "constants.h"
#include "mesh.h"
#include "mesh_cell_pass.h"
#include "mesh_particle_pass.h"
#include "imc_state.h"
#include "timing_functions.h"
#include "input.h"
#include "decompose_mesh.h"

using std::vector;
using std::endl;
using std::cout;
using std::string;
using Constants::PARTICLE_PASS;
using Constants::CELL_PASS;

namespace mpi = boost::mpi;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);

  mpi::environment env(argc, argv);
  mpi::communicator world;
  unsigned int rank = world.rank();
  unsigned int size = world.size();

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

  //IMC state setup
  IMC_State *imc_state;
  imc_state = new IMC_State(input);

  // make mesh from input object
  Mesh *mesh;
  if (input->get_dd_mode() == PARTICLE_PASS)
    Mesh *mesh = new Mesh_Particle_Pass(input, rank, size);
  else if (input->get_dd_mode() == CELL_PASS)
    Mesh *mesh = new Mesh_Cell_Pass(input, rank, size);
  else 
    cout<<"Error: DD type not recognized"

  // decompose mesh with ParMETIS and Boost MPI
  decompose_mesh(mesh, world, argc, argv);

  MPI::COMM_WORLD.Barrier();
  //print_MPI_out(mesh, rank, size);  

  //timing 
  struct timeval start,end;
  struct timezone tzp;
  gettimeofday(&start, &tzp); 

/******************************************************************************/ 
// TRT PHYSICS CALCULATION
/******************************************************************************/
  
  if (input->get_dd_mode() == PARTICLE_PASS)
    imc_particle_pass(mesh);
  else if (input->get_dd_mode() == CELL_PASS)
    imc_cell_pass(mesh);


  if (rank==0) {
    gettimeofday(&end, &tzp);
    print_elapsed_inside("\nruntime:",&start, &end);
    cout<<"Total RMA Operations: "<<imc_state->get_total_RMA()<<endl;
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
  }

  MPI::COMM_WORLD.Barrier();
  delete mesh;
  MPI::Finalize();
}
