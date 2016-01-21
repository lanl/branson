/*
  Author: Alex Long
  Date: 7/24/2014
  Name: imc_state.h
*/
#ifndef imc_state_h_
#define imc_state_h_


#include <iostream>
#include <vector>
#include <cmath>

#include "input.h"
#include "photon.h"
#include "constants.h"
#include "RNG.h"

class IMC_State
{
  public:
  IMC_State(Input *input)
    : m_dt(input->get_dt()),
      m_time(input->get_time_start()),
      m_time_stop(input->get_time_finish()),
      m_step(1),
      m_dt_mult(input->get_time_mult()),
      m_dt_max(input->get_dt_max()),
      m_total_photons(input->get_number_photons()),
      m_trans_photons(0)
    {
      m_RNG = new RNG();
      m_RNG->set_seed(input->get_rng_seed()+MPI::COMM_WORLD.Get_rank()*4106);
      pre_census_E = 0.0;
      post_census_E = 0.0;      
      pre_mat_E = 0.0;
      post_mat_E = 0.0;
      exit_E = 0.0;
      absorbed_E = 0.0;
      total_off_rank_reads=0;
    }

  ~IMC_State() { delete m_RNG;}

/*****************************************************************************/
/* const functions                                                           */
/*****************************************************************************/
  double get_dt(void) const {return m_dt;}
  unsigned int get_step(void) const {return m_step;}
  unsigned int get_total_step_photons(void) const {return m_total_photons;}
  unsigned int get_census_size(void) const {return census_size;}
  unsigned int get_total_RMA(void) const {return total_off_rank_reads;}
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
  void print_conservation(unsigned int dd_type ) {
    using std::cout;
    using std::endl;
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
    unsigned int g_census_size=0;
    unsigned int g_off_rank_reads=0;
    unsigned int g_trans_photons=0;
    unsigned int g_n_photon_messages=0;
    unsigned int g_n_photons_sent=0;
    unsigned int g_n_cell_messages=0;
    unsigned int g_n_cells_sent=0;
    unsigned int g_n_sends_posted=0;
    unsigned int g_n_sends_completed=0;
    unsigned int g_n_receives_posted=0;
    unsigned int g_n_receives_completed=0;

    //reduce energy conservation values
    MPI::COMM_WORLD.Allreduce(&absorbed_E, &g_absorbed_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&emission_E, &g_emission_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&pre_census_E, &g_pre_census_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&pre_mat_E, &g_pre_mat_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&post_census_E, &g_post_census_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&post_mat_E, &g_post_mat_E, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&exit_E, &g_exit_E, 1, MPI_DOUBLE, MPI_SUM);

    // reduce diagnostic values
    MPI::COMM_WORLD.Allreduce(&census_size, &g_census_size, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&m_trans_photons, &g_trans_photons, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&off_rank_reads, &g_off_rank_reads, 1, MPI_UNSIGNED, MPI_SUM);

    MPI::COMM_WORLD.Allreduce(&n_photon_messages, &g_n_photon_messages, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_photons_sent, &g_n_photons_sent, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_cell_messages, &g_n_cell_messages, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_cells_sent, &g_n_cells_sent, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_sends_posted, &g_n_sends_posted, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_sends_completed, &g_n_sends_completed, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_receives_posted, &g_n_receives_posted, 1, MPI_UNSIGNED, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&n_receives_completed, &g_n_receives_completed, 1, MPI_UNSIGNED, MPI_SUM);
    
    double rad_conservation = (g_absorbed_E + g_post_census_E + g_exit_E) - 
                            (g_pre_census_E + g_emission_E + source_E);

    double mat_conservation = g_post_mat_E - (g_pre_mat_E + g_absorbed_E - g_emission_E);

    if (MPI::COMM_WORLD.Get_rank() == 0) {
      cout<<"Total Photons transported: "<<g_trans_photons<<endl;
      cout<<"Emission E: "<<g_emission_E<<", Absorption E: "<<g_absorbed_E;
      cout<<", Exit E: "<<g_exit_E<<endl;
      cout<<"Pre Census E: "<<g_pre_census_E<<" Post census E: "<<g_post_census_E<<" Post Census Size: "<< g_census_size<<endl;
      cout<<"Pre mat E: "<<g_pre_mat_E<<" Post mat E: "<<g_post_mat_E<<endl;
      cout<<"Radiation Conservation: "<<rad_conservation<<endl;
      cout<<"Material Conservation: "<<mat_conservation<<endl;
      cout<<"Sends posted: "<<g_n_sends_posted;
      cout<<", sends completed: "<<g_n_sends_completed<<endl;
      cout<<"Receives posted: "<<g_n_receives_posted;
      cout<<", receives completed: "<<g_n_receives_completed<<endl;
      if (dd_type == PARTICLE_PASS) {
        cout<<"Photons messages sent: "<<g_n_photon_messages;
        cout<<", Total photons sent: "<<g_n_photons_sent<<endl;
      }
      else {
        cout<<"Cell messages sent: "<<g_n_cell_messages;
        cout<<", Total cells sent: "<<g_n_cells_sent<<endl;
        cout<<"Total RMA requests: "<<g_off_rank_reads<<endl;
      }
      total_off_rank_reads+=g_off_rank_reads;
    } // if rank==0
  }

  RNG* get_rng(void) const { return m_RNG;}

  //non const functions
  void next_time_step(void) {
    //print_timestep_header();
    m_time += m_dt;
    m_dt = get_next_dt();
    m_step++;
  }

  void set_transported_photons(unsigned int _trans_photons) {m_trans_photons = _trans_photons;}
  void set_pre_census_E(double _pre_census_E) {pre_census_E = _pre_census_E;}
  void set_post_census_E(double _post_census_E) {post_census_E = _post_census_E;}
  void set_pre_mat_E(double _pre_mat_E) {pre_mat_E = _pre_mat_E;}
  void set_post_mat_E(double _post_mat_E) {post_mat_E = _post_mat_E;}
  void set_emission_E(double _emission_E) {emission_E = _emission_E;}
  void set_source_E(double _source_E) {source_E = _source_E;}
  void set_absorbed_E(double _absorbed_E) {absorbed_E = _absorbed_E;}
  void set_exit_E(double _exit_E) {exit_E = _exit_E;}
  //set diagnostic values
  void set_census_size(unsigned int _census_size) {census_size = _census_size;}
  void set_off_rank_read(unsigned int _off_rank_reads) {off_rank_reads = _off_rank_reads;}

  void set_n_photon_messages(unsigned int _n_photon_messages) {
    n_photon_messages=_n_photon_messages;
  }
  void set_n_photons_sent(unsigned int _n_photons_sent) {
    n_photons_sent=_n_photons_sent;
  }
  void set_n_cell_messages(unsigned int _n_cell_messages) {
    n_cell_messages=_n_cell_messages;
  }
  void set_n_cells_sent(unsigned int _n_cells_sent) {
    n_cells_sent=_n_cells_sent;
  }
  void set_n_sends_posted(unsigned int _n_sends_posted) {
    n_sends_posted=_n_sends_posted;
  }
  void set_n_sends_completed(unsigned int _n_sends_completed) {
    n_sends_completed=_n_sends_completed;
  }
  void set_n_receives_posted(unsigned int _n_receives_posted) {
    n_receives_posted=_n_receives_posted;
  }
  void set_n_receives_completed(unsigned int _n_receives_completed) { 
    n_receives_completed=_n_receives_completed;
  }

/*****************************************************************************/
/* member variables and private functions                                    */
/*****************************************************************************/
  private:
  //time
  double m_dt; //! Current time step size (sh)
  double m_time; //! Current time (sh)
  double m_time_stop; //! End time (sh)
  unsigned int m_step; //! Time step (start at 1)
  double m_dt_mult; //! Time step multiplier
  double m_dt_max;  //! Max time step
  //photons
  unsigned int m_total_photons; //! Photons desired
  unsigned int m_trans_photons; //! Photons transported

  //conservation
  double pre_census_E; //! Census energy at the beginning of the timestep
  double post_census_E; //! Census energy at the end of the timestep
  double pre_mat_E; //! Material energy at the beginning of timestep
  double post_mat_E; //! Material energy at the end of timestep
  double emission_E; //! Energy emitted this timestep
  double exit_E; //! Energy exiting problem
  double absorbed_E; //! Total absorbed energy 
  double source_E; //! Sourced energy

  //diagnostic
  unsigned int census_size; //! Number of particles in census
  unsigned int off_rank_reads; //! Number of off rank reads by this rank
  unsigned int total_off_rank_reads; //! Number of RMA reads total for simulation
  unsigned int n_photon_messages; //! Number of photon messages
  unsigned int n_photons_sent; //! Number of photons passed
  unsigned int n_cell_messages; //! Number of cell messages
  unsigned int n_cells_sent; //! Number of cells passed
  unsigned int n_sends_posted; //! Number of sent messages posted
  unsigned int n_sends_completed; //! Number of sent messages completed
  unsigned int n_receives_posted; //! Number of received messages completed
  unsigned int n_receives_completed; //! Number of received messages completed

  //RNG
  RNG* m_RNG;
};

#endif
