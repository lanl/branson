/*
  Author: Alex Long
  Date: 12/1/2015
  Name: imc_drivers.h
*/

#include <boost/mpi.hpp>
#include <vector>

#include "decompose_photons.h"
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
                          IMC_Parameters *imc_p,
                          mpi::communicator world) {
  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  Photon* photon_vec;
  Photon* census_list;
  unsigned int n_photon;

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

    make_photons(mesh, imc_state, photon_vec, n_photon, global_source_energy);

    imc_state->set_transported_photons(n_photon);
    //append census list to photon vector, rebalance
    if (imc_state->get_step() > 1) {
      unsigned int n_census = imc_state->get_census_size();
      imc_state->set_pre_census_E(get_photon_list_energy(census_list, n_census));
      imc_state->set_transported_photons(n_photon + n_census);
      //rebalances census photons, adds new census to photon vector
      on_rank_rebalance_photons(photon_vec, 
                                n_photon, 
                                census_list, 
                                n_census, 
                                mesh, 
                                world);
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
                      census_list, imc_p->get_check_frequency(), world);

    //using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, 
                              &abs_E[0], 
                              mesh->get_global_num_cells(), 
                              MPI_DOUBLE, 
                              MPI_SUM);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    imc_state->print_conservation();

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
                              IMC_Parameters *imc_p,
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
                                                  imc_p->get_n_user_photon());

    imc_state->set_pre_census_E(get_photon_list_E(census_photons)); 

    Source source(mesh, imc_state, global_source_energy, census_photons);

    imc_state->set_transported_photons(source.get_n_photon());

    census_photons = transport_photons( source, 
                                        mesh, 
                                        imc_state, 
                                        abs_E, 
                                        imc_p->get_check_frequency(),
                                        world);
    
    mesh->update_temperature(abs_E, imc_state);
    //update time for next step
    
    imc_state->print_conservation();
    imc_state->next_time_step();
  }
}
