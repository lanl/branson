/*
  Author: Alex Long
  Date: 7/24/2014
  Name: imc_state.h
*/
#ifndef imc_state_h_
#define imc_state_h_


#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

#include "input.h"
#include "photon.h"
#include "constants.h"
#include "RNG.h"

class IMC_State
{
  public:
  IMC_State(Input *input, uint32_t _rank)
    : rank(_rank),
      m_dt(input->get_dt()),
      m_time(input->get_time_start()),
      m_time_stop(input->get_time_finish()),
      m_step(1),
      m_dt_mult(input->get_time_mult()),
      m_dt_max(input->get_dt_max())
    {
      pre_census_E = 0.0;
      post_census_E = 0.0;
      pre_mat_E = 0.0;
      post_mat_E = 0.0;
      emission_E = 0.0;
      exit_E = 0.0;
      absorbed_E = 0.0;
      source_E = 0.0;

      //64 bit
      trans_particles = 0;
      census_size = 0;
      step_particles_sent = 0;
      total_particles_sent=0;
      //32 bit
      total_cells_requested=0;
      total_cells_sent=0;
      total_cell_messages=0;
      total_particle_messages=0;

      step_particle_messages=0;
      step_cells_requested=0;
      step_cell_messages=0;
      step_cells_sent=0;

      step_sends_posted=0;
      step_sends_completed=0;
      step_receives_posted=0;
      step_receives_completed=0;

      m_RNG = new RNG();
      m_RNG->set_seed(input->get_rng_seed()+rank*4106);
    }

  ~IMC_State() { delete m_RNG;}

/*****************************************************************************/
/* const functions                                                           */
/*****************************************************************************/
  double get_dt(void) const {return m_dt;}
  uint32_t get_step(void) const {return m_step;}
  uint64_t get_transported_particles(void) const {return trans_particles;}
  uint64_t get_census_size(void) const {return census_size;}
  uint64_t get_step_particles_sent(void) const {return step_particles_sent;}
  uint64_t get_total_particles_sent(void) const {return total_particles_sent;}
  double get_pre_census_E(void) {return pre_census_E;}
  double get_emission_E(void) {return emission_E;}
  double get_next_dt(void) const {
    double next_dt;
    //multiply by dt_mult if it's going to be less than dt_max
    //other set to dt_max
    if (m_dt*m_dt_mult < m_dt_max) next_dt = m_dt*m_dt_mult;
    else next_dt = m_dt_max;
    //don't overrun the end time
    if (m_time + next_dt > m_time_stop) next_dt = m_time_stop - m_time;

    return next_dt;
  }  

  bool finished(void) const {
    using std::abs;
    if ( abs(m_time-m_time_stop) < 1.0e-8) return true;
    else return false;
  }

  void print_timestep_header(void) const  {
    using std::cout;
    using std::endl;
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
    cout<<"Step: "<<m_step<<"  Start Time: "<<m_time<<"  End Time: "<<m_time +m_dt<<"  dt: "<<m_dt<<endl;
  }

/*****************************************************************************/
/* non-const functions                                                       */
/*****************************************************************************/
  void print_conservation(uint32_t dd_type) {
    using std::cout;
    using std::endl;
    using std::plus;
    using Constants::PARTICLE_PASS;
    using Constants::CELL_PASS;

    //define global value
    double g_absorbed_E=0.0; 
    double g_emission_E=0.0;
    double g_pre_census_E=0.0;
    double g_pre_mat_E=0.0;
    double g_post_census_E=0.0;
    double g_post_mat_E=0.0;
    double g_exit_E = 0.0;
    // 64 bit global integers
    uint64_t g_census_size=0;
    uint64_t g_trans_particles=0;
    uint64_t g_step_particles_sent=0;
    // 32 bit global integers
    uint32_t g_step_particle_messages=0;
    uint32_t g_step_cells_requested=0;
    uint32_t g_step_cell_messages=0;
    uint32_t g_step_cells_sent=0;
    uint32_t g_step_sends_posted=0;
    uint32_t g_step_sends_completed=0;
    uint32_t g_step_receives_posted=0;
    uint32_t g_step_receives_completed=0;

    //reduce energy conservation values (double)
    MPI_Allreduce(&absorbed_E, &g_absorbed_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&emission_E, &g_emission_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(&pre_census_E, &g_pre_census_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pre_mat_E, &g_pre_mat_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&post_census_E, &g_post_census_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&post_mat_E, &g_post_mat_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&exit_E, &g_exit_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // reduce diagnostic values
    // 64 bit integer reductions
    MPI_Allreduce(&trans_particles, &g_trans_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_particles_sent, &g_step_particles_sent, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&census_size, &g_census_size, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    // 32 bit integer reductions
    MPI_Allreduce(&step_cells_requested, &g_step_cells_requested, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_particle_messages, &g_step_particle_messages, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_cell_messages, &g_step_cell_messages, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_cells_sent, &g_step_cells_sent, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_sends_posted, &g_step_sends_posted, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_sends_completed, &g_step_sends_completed, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_receives_posted, &g_step_receives_posted, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&step_receives_completed, &g_step_receives_completed, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    double rad_conservation = (g_absorbed_E + g_post_census_E + g_exit_E) - 
      (g_pre_census_E + g_emission_E + source_E);

    double mat_conservation = g_post_mat_E - (g_pre_mat_E + g_absorbed_E -
      g_emission_E);

    // update total simulation counters
    total_cells_requested+=g_step_cells_requested;
    total_cells_sent+=g_step_cells_sent;
    total_cell_messages+=g_step_cell_messages;
    total_particles_sent+=g_step_particles_sent;
    total_particle_messages+=g_step_particle_messages;
    

    if (rank == 0) {
      cout<<"Total Photons transported: "<<g_trans_particles<<endl;
      cout<<"Emission E: "<<g_emission_E<<", Absorption E: "<<g_absorbed_E;
      cout<<", Exit E: "<<g_exit_E<<endl;
      cout<<"Pre census E: "<<g_pre_census_E<<" Post census E: ";
      cout<<g_post_census_E<<" Post census Size: "<< g_census_size<<endl;
      cout<<"Pre mat E: "<<g_pre_mat_E<<" Post mat E: "<<g_post_mat_E<<endl;
      cout<<"Radiation conservation: "<<rad_conservation<<endl;
      cout<<"Material conservation: "<<mat_conservation<<endl;
      cout<<"Sends posted: "<<g_step_sends_posted;
      cout<<", sends completed: "<<g_step_sends_completed<<endl;
      cout<<"Receives posted: "<<g_step_receives_posted;
      cout<<", receives completed: "<<g_step_receives_completed<<endl;
      if (dd_type == PARTICLE_PASS) {
        cout<<"Step particles messages sent: "<<g_step_particle_messages;
        cout<<", Step particles sent: "<<g_step_particles_sent<<endl;
      }
      else {
        cout<<"Step cell messages sent: "<<g_step_cell_messages;
        cout<<", Step cells sent: "<<g_step_cells_sent<<endl;
        cout<<"Step cells requested: "<<g_step_cells_requested<<endl;
      }
    } // if rank==0
  }

  void print_simulation_footer(uint32_t dd_type ) {
    using std::cout;
    using std::endl;
    using Constants::PARTICLE_PASS;
    using Constants::CELL_PASS;
    if (dd_type == PARTICLE_PASS) {
      cout<<"Total particles sent: "<<total_particles_sent<<endl;
      cout<<"Total particle messages: "<<total_particle_messages<<endl;
    }
    else {
      cout<<"Total cells requested: "<<total_cells_requested<<endl;
      cout<<"Total cells sent: "<<total_cells_sent<<endl;
      cout<<"Total cell messages: "<<total_cell_messages<<endl;
    }
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
  }



  RNG* get_rng(void) const { return m_RNG;}

  //non const functions
  void next_time_step(void) {
    //print_timestep_header();
    m_time += m_dt;
    m_dt = get_next_dt();
    m_step++;
  }

  void set_pre_census_E(double _pre_census_E) {pre_census_E = _pre_census_E;}
  void set_post_census_E(double _post_census_E) {
    post_census_E = _post_census_E;
  }
  void set_pre_mat_E(double _pre_mat_E) {pre_mat_E = _pre_mat_E;}
  void set_post_mat_E(double _post_mat_E) {post_mat_E = _post_mat_E;}
  void set_emission_E(double _emission_E) {emission_E = _emission_E;}
  void set_source_E(double _source_E) {source_E = _source_E;}
  void set_absorbed_E(double _absorbed_E) {absorbed_E = _absorbed_E;}
  void set_exit_E(double _exit_E) {exit_E = _exit_E;}
  //set diagnostic values (64 bit)
  void set_transported_particles(uint64_t _trans_particles) {
    trans_particles = _trans_particles;
  }
  void set_census_size(uint64_t _census_size) {census_size = _census_size;}
  void set_step_particles_sent(uint64_t _step_particles_sent) {
    step_particles_sent=_step_particles_sent;
  }
  // set diagnostic values (32 bit)
  void set_step_cells_requested(uint32_t _step_cells_requested) {
    step_cells_requested = _step_cells_requested;
  }
  void set_step_particle_messages(uint32_t _step_particle_messages) {
    step_particle_messages=_step_particle_messages;
  }
  void set_step_cell_messages(uint32_t _step_cell_messages) {
    step_cell_messages=_step_cell_messages;
  }
  void set_step_cells_sent(uint32_t _step_cells_sent) {
    step_cells_sent=_step_cells_sent;
  }
  void set_step_sends_posted(uint32_t _step_sends_posted) {
    step_sends_posted=_step_sends_posted;
  }
  void set_step_sends_completed(uint32_t _step_sends_completed) {
    step_sends_completed=_step_sends_completed;
  }
  void set_step_receives_posted(uint32_t _step_receives_posted) {
    step_receives_posted=_step_receives_posted;
  }
  void set_step_receives_completed(uint32_t _step_receives_completed) { 
    step_receives_completed=_step_receives_completed;
  }

/*****************************************************************************/
/* member variables and private functions                                    */
/*****************************************************************************/
  private:
  uint32_t rank; //! Rank owning this states object
  //time
  double m_dt; //! Current time step size (sh)
  double m_time; //! Current time (sh)
  double m_time_stop; //! End time (sh)
  uint32_t m_step; //! Time step (start at 1)
  double m_dt_mult; //! Time step multiplier
  double m_dt_max;  //! Max time step

  //conservation
  double pre_census_E; //! Census energy at the beginning of the timestep
  double post_census_E; //! Census energy at the end of the timestep
  double pre_mat_E; //! Material energy at the beginning of timestep
  double post_mat_E; //! Material energy at the end of timestep
  double emission_E; //! Energy emitted this timestep
  double exit_E; //! Energy exiting problem
  double absorbed_E; //! Total absorbed energy 
  double source_E; //! Sourced energy

  //diagnostic 64 bit integers relating to particle counts
  uint64_t trans_particles; //! Particles transported
  uint64_t census_size; //! Number of particles in census

  uint64_t step_particles_sent; //! Number of particles passed

  //! Total number of particles sent for simulation
  uint64_t total_particles_sent; 

  //diagnostic 32 bit integers relating to messages and cell counts
  //! Total number of cells requested for simulation
  uint32_t total_cells_requested; 

  //! Total number of cells sent for simulation
  uint32_t total_cells_sent;

  //! Total number of cell message for simulation
  uint32_t total_cell_messages;

  //! Total number of particle messages sent for simulation
  uint32_t total_particle_messages;
 
  uint32_t step_particle_messages; //! Number of particle messages
  uint32_t step_cells_requested; //! Number of cells requested by this rank
  uint32_t step_cell_messages; //! Number of cell messages
  uint32_t step_cells_sent; //! Number of cells passed
  uint32_t step_sends_posted; //! Number of sent messages posted
  uint32_t step_sends_completed; //! Number of sent messages completed
  uint32_t step_receives_posted; //! Number of received messages completed
  uint32_t step_receives_completed;  //! Number of received messages completed

  //RNG
  RNG* m_RNG;
};

#endif
