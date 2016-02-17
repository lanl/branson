/*
  Author: Alex Long
  Date: 12/1/2015
  Name: imc_drivers.h
*/

#include <boost/mpi.hpp>
#include <vector>

#include "imc_state.h"
#include "imc_parameters.h"
#include "mesh.h"
#include "RNG.h"
#include "source.h"
#include "transport.h"
#include "transport_particle_pass.h"
#include "transport_mesh_pass.h"

namespace mpi = boost::mpi;

void imc_cell_pass_driver(const int& rank, 
                          Mesh *mesh, 
                          IMC_State *imc_state,
                          IMC_Parameters *imc_parameters,
                          mpi::communicator world) {
  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, 
                              &global_source_energy, 
                              1, 
                              MPI_DOUBLE, 
                              MPI_SUM);

    //make initial census photons
    //on subsequent timesteps, this comes from transport
    if (imc_state->get_step() == 1) {
      census_photons = make_initial_census_photons(mesh, 
                                                  imc_state, 
                                                  global_source_energy,
                                                  imc_parameters->get_n_user_photon());
    }

    imc_state->set_pre_census_E(get_photon_list_E(census_photons));

    Source source(mesh, imc_parameters, global_source_energy, census_photons);

    imc_state->set_transported_particles(source.get_n_photon());
    //cout<<"Rank: "<<rank<<" about to transport ";
    //cout<<n_photon<<" particles."<<endl;

    //cell properties are set in calculate_photon_energy. 
    //make sure everybody gets here together so that windows are not changing 
    //when transport starts
    MPI::COMM_WORLD.Barrier();

    //transport photons
    census_photons = transport_mesh_pass(source, 
                                          mesh, 
                                          imc_state, 
                                          imc_parameters, 
                                          abs_E, 
                                          world);

    //using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, 
                              &abs_E[0], 
                              mesh->get_global_num_cells(), 
                              MPI_DOUBLE, 
                              MPI_SUM);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    imc_state->print_conservation(imc_parameters->get_dd_mode());

    //purge the working mesh, it will be updated by other ranks and is now 
    //invalid
    mesh->purge_working_mesh();

    //update time for next step
    imc_state->next_time_step();
  }
}


void imc_particle_pass_driver(const int& rank, 
                              Mesh *mesh, 
                              IMC_State *imc_state,
                              IMC_Parameters *imc_parameters,
                              mpi::communicator world) {

  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, 
                              &global_source_energy, 
                              1, 
                              MPI_DOUBLE, 
                              MPI_SUM);

    //make initial census photons
    //on subsequent timesteps, this comes from transport
    if (imc_state->get_step() ==1)
      census_photons = make_initial_census_photons(mesh, 
                                                  imc_state, 
                                                  global_source_energy,
                                                  imc_parameters->get_n_user_photon());

    imc_state->set_pre_census_E(get_photon_list_E(census_photons)); 

    Source source(mesh, imc_parameters, global_source_energy, census_photons);

    imc_state->set_transported_particles(source.get_n_photon());

    census_photons = transport_particle_pass( source, 
                                              mesh, 
                                              imc_state, 
                                              imc_parameters,
                                              abs_E, 
                                              world);
          
    mesh->update_temperature(abs_E, imc_state);
    //update time for next step
    
    imc_state->print_conservation(imc_parameters->get_dd_mode());
    imc_state->next_time_step();
  }
}
