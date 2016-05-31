//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc_drivers.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  Functions to run each domain decomposed IMC method
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef imc_drivers_h_
#define imc_drivers_h_

#include <iostream>
#include <mpi.h>
#include <functional>
#include <vector>

#include "completion_manager_milagro.h"
#include "completion_manager_rma.h"
#include "mesh_rma_manager.h"
#include "mpi_types.h"
#include "imc_state.h"
#include "imc_parameters.h"
#include "load_balance.h"
#include "mesh.h"
#include "RNG.h"
#include "source.h"
#include "census_creation.h"
#include "transport_particle_pass.h"
#include "transport_mesh_pass.h"
#include "transport_mesh_pass_rma.h"
#include "write_silo.h"

void imc_cell_pass_driver(const int& rank, 
                          const int& n_rank,
                          Mesh *mesh, 
                          IMC_State *imc_state,
                          IMC_Parameters *imc_parameters,
                          MPI_Types *mpi_types) {
  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;

  // make object that handles completion messages, completion
  // object set up the binary tree structure
  Completion_Manager *comp;
  if (imc_parameters->get_completion_method() == Constants::RMA_COMPLETION)
    comp = new Completion_Manager_RMA(rank, n_rank);
  else
    comp = new Completion_Manager_Milagro(rank, n_rank);

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

    //make initial census photons
    //on subsequent timesteps, this comes from transport
    if (imc_state->get_step() == 1) {
      census_photons = make_initial_census_photons(mesh, imc_state, 
        global_source_energy, imc_parameters->get_n_user_photon());
    }

    imc_state->set_pre_census_E(get_photon_list_E(census_photons));

    Source source(mesh, imc_parameters, global_source_energy, census_photons);
    load_balance(rank, n_rank, source.get_n_photon(), source.get_work_vector(),
      census_photons, mpi_types);
    // get new particle count after load balance. Group particle work by cell
    source.post_lb_prepare_source();

    imc_state->set_transported_particles(source.get_n_photon());
    //cout<<"Rank: "<<rank<<" about to transport ";
    //cout<<n_photon<<" particles."<<endl;

    //cell properties are set in calculate_photon_energy. 
    //make sure everybody gets here together so that windows are not changing 
    //when transport starts
    MPI_Barrier(MPI_COMM_WORLD);

    //transport photons
    census_photons = transport_mesh_pass(source, 
                                          mesh, 
                                          imc_state, 
                                          imc_parameters,
                                          comp,
                                          abs_E,
                                          mpi_types);

    //using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI_Allreduce(MPI_IN_PLACE, 
                  &abs_E[0], 
                  mesh->get_global_num_cells(), 
                  MPI_DOUBLE, 
                  MPI_SUM, 
                  MPI_COMM_WORLD);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    imc_state->print_conservation(imc_parameters->get_dd_mode());

    //purge the working mesh, it will be updated by other ranks and is now 
    //invalid
    mesh->purge_working_mesh();

    //update time for next step
    imc_state->next_time_step();
  }

  // delete completion manager (closes and free window in RMA version)
  delete comp;
}


void imc_rma_cell_pass_driver(const int& rank, 
                              const int& n_rank,
                              Mesh *mesh, 
                              IMC_State *imc_state,
                              IMC_Parameters *imc_parameters,
                              MPI_Types *mpi_types) {
  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;

  //make object that handles RMA mesh requests and start access
  RMA_Manager *rma_manager = new RMA_Manager(rank, mesh->get_off_rank_bounds(),
    mesh->get_global_num_cells(), mpi_types, mesh->get_mesh_window_ref());
  rma_manager->start_access();

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

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
    load_balance(rank, n_rank, source.get_n_photon(), source.get_work_vector(),
      census_photons, mpi_types);
    // get new particle count after load balance, group particle work by cell
    source.post_lb_prepare_source();

    imc_state->set_transported_particles(source.get_n_photon());
    //cout<<"Rank: "<<rank<<" about to transport ";
    //cout<<n_photon<<" particles."<<endl;

    //cell properties are set in calculate_photon_energy. 
    //make sure everybody gets here together so that windows are not changing 
    //when transport starts
    MPI_Barrier(MPI_COMM_WORLD);

    //transport photons
    census_photons = transport_mesh_pass_rma( source, 
                                              mesh, 
                                              imc_state, 
                                              imc_parameters,
                                              rma_manager,
                                              abs_E,
                                              mpi_types);

    //using MPI_IN_PLACE allows the same vector to send and be overwritten
    MPI_Allreduce(MPI_IN_PLACE, 
                  &abs_E[0], 
                  mesh->get_global_num_cells(), 
                  MPI_DOUBLE, 
                  MPI_SUM, 
                  MPI_COMM_WORLD);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    imc_state->print_conservation(imc_parameters->get_dd_mode());

    //purge the working mesh, it will be updated by other ranks and is now 
    //invalid
    mesh->purge_working_mesh();

    if (imc_parameters->get_write_silo_flag()) {
      // write SILO file
      vector<uint32_t> n_requests = rma_manager->get_n_request_vec();
      write_silo(mesh, imc_state->get_time(), imc_state->get_step(), 
        imc_state->get_rank_transport_runtime(), rank, n_rank, n_requests);
    }
    //reset rma_manager object for next timestep
    rma_manager->end_timestep();

    //update time for next step
    imc_state->next_time_step();    
  }
  //close access to MPI windows in RMA_Manger object and delete
  rma_manager->end_access();
  delete rma_manager;
}

void imc_particle_pass_driver(const int& rank, 
                              const int& n_rank,
                              Mesh *mesh, 
                              IMC_State *imc_state,
                              IMC_Parameters *imc_parameters,
                              MPI_Types * mpi_types) {

  using std::vector;
  vector<double> abs_E(mesh->get_global_num_cells(), 0.0);
  vector<Photon> census_photons;

  // make object that handles completion messages, completion
  // object set up the binary tree structure
  Completion_Manager *comp;
  if (imc_parameters->get_completion_method() == Constants::RMA_COMPLETION)
    comp = new Completion_Manager_RMA(rank, n_rank);
  else
    comp = new Completion_Manager_Milagro(rank, n_rank);

  while (!imc_state->finished())
  {
    if (rank==0) imc_state->print_timestep_header();

    //set opacity, Fleck factor, all energy to source
    mesh->calculate_photon_energy(imc_state);

    //all reduce to get total source energy to make correct number of
    //particles on each rank
    double global_source_energy = mesh->get_total_photon_E();
    MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

    //make initial census photons
    //on subsequent timesteps, this comes from transport
    if (imc_state->get_step() ==1)
      census_photons = make_initial_census_photons(mesh, 
                                                  imc_state, 
                                                  global_source_energy,
                                                  imc_parameters->get_n_user_photon());

    imc_state->set_pre_census_E(get_photon_list_E(census_photons)); 

    Source source(mesh, imc_parameters, global_source_energy, census_photons);
    // no load balancing in particle passing method, just call the method
    // to get accurate count
    source.post_lb_prepare_source();

    imc_state->set_transported_particles(source.get_n_photon());

    census_photons = transport_particle_pass( source, 
                                              mesh, 
                                              imc_state, 
                                              imc_parameters,
                                              mpi_types,
                                              comp,
                                              abs_E);
          
    mesh->update_temperature(abs_E, imc_state);

    //update time for next step    
    imc_state->print_conservation(imc_parameters->get_dd_mode());
    imc_state->next_time_step();
  }

  // delete completion manager (closes and free window in RMA version)
  delete comp;
}

#endif // imc_drivers_h_
//---------------------------------------------------------------------------//
// end of imc_drivers.h
//---------------------------------------------------------------------------//
