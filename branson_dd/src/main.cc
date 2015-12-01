
#include <iostream>
#include <vector>
#include <string>
#include <boost/mpi.hpp>
#include <time.h>
#include <sys/time.h>

#include "mesh.h"
#include "imc_state.h"
#include "timing_functions.h"
#include "input.h"
#include "decompose_mesh.h"
#include "decompose_photons.h"
#include "transport.h"

using std::vector;
using std::endl;
using std::cout;
using std::string;

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
  Mesh *mesh = new Mesh(input);

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
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  Photon* photon_vec;
  Photon* census_list;
  unsigned int n_photon;

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, 
                              &global_source_energy, 
                              1, 
                              MPI_DOUBLE, 
                              MPI_SUM);

    //make photons on mesh owned by rank
    if (!input->get_stratified_bool()) make_photons(mesh, 
                                                    imc_state, 
                                                    photon_vec, n_photon, 
                                                    global_source_energy);
    else make_stratified_photons( mesh, 
                                  imc_state, 
                                  photon_vec, 
                                  n_photon, 
                                  global_source_energy);

    imc_state->set_transported_photons(n_photon);
    //append census list to photon vector, rebalance
    if (imc_state->get_step() > 1) {
      unsigned int n_census = imc_state->get_census_size();
      imc_state->set_pre_census_E(get_photon_list_energy(census_list, n_census));
      imc_state->set_transported_photons(n_photon + n_census);
      //rebalances census photons, adds new census to photon vector
      on_rank_rebalance_photons(photon_vec, n_photon, census_list, n_census, mesh, world);
    }

    proto_load_balance_photons( photon_vec, 
                                n_photon, 
                                mesh, 
                                world);

    //cout<<"Rank: "<<rank<<" about to transport ";
    //cout<<n_photon<<" particles."<<endl;

    //cell properties are set in calculate_photon_energy. 
    //make sure everybody gets here together so that windows are not changing 
    //when transport starts
    MPI::COMM_WORLD.Barrier();

    //transport photons
    transport_photons(photon_vec, n_photon, mesh, imc_state, abs_E, 
                      census_list, input->get_check_frequency(), world);

    //using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &abs_E[0], mesh->get_global_num_cells(), MPI_DOUBLE, MPI_SUM);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    imc_state->print_conservation();

    //purge the working mesh, it will be updated by other ranks and is now 
    //invalid
    mesh->purge_working_mesh();

    //update time for next step
    imc_state->next_time_step();
  }

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
