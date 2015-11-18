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
#include <omp.h>

#include "input.h"
#include "RNG.h"

using std::pow;
using std::cout;
using std::endl;
using std::abs;
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
      m_total_photons(input->get_number_photons())
    {
      for (int i=0; i<omp_get_max_threads(); i++) {
        m_RNG.push_back(new RNG());
        m_RNG[i]->set_seed(input->get_rng_seed()+i);
      }
      pre_census_E = 0.0;
      post_census_E = 0.0;      
      pre_mat_E = 0.0;
      post_mat_E = 0.0;
      exit_E = 0.0;
      absorbed_E = 0.0;

      m_census_photons=0;
    }

  ~IMC_State() { for (int i=0; i<omp_get_max_threads(); i++) delete m_RNG[i];}

  //const functions
  double get_dt(void) const {return m_dt;}
  unsigned int get_step(void) const {return m_step;}
  unsigned int get_total_step_photons(void) const {return m_total_photons;}
  double get_pre_census_E(void) {return pre_census_E;}
  double get_emission_E(void) {return emission_E;}
  unsigned int get_census_photons(void) {return m_census_photons;}

  bool finished(void) const {
    if ( abs(m_time-m_time_stop) < 1.0e-8) return true;
    else return false;
  }

  void print_timestep_header(void) const  {
    cout<<"****************************************";
    cout<<"****************************************"<<endl;
    cout<<"Step: "<<m_step<<"  Start Time: "<<m_time<<"  End Time: "<<m_time +m_dt<<"  dt: "<<m_dt<<endl;
  }

  void print_conservation(void) {
    double rad_conservation = (absorbed_E + post_census_E + exit_E) - 
                            (pre_census_E + emission_E + source_E);

    double mat_conservation = post_mat_E - (pre_mat_E + absorbed_E - emission_E);

    cout<<"Emission E: "<<emission_E<<" Exit E: "<<exit_E<<endl;
    cout<<"Pre Census E: "<<pre_census_E<<" Post census E: "<<post_census_E<<" Post Census Size: "<<m_census_photons<<endl;
    cout<<"Pre mat E: "<<pre_mat_E<<" Post mat E: "<<post_mat_E<<endl;
    cout<<"Radiation Conservation: "<<rad_conservation<<endl;
    cout<<"Material Conservation: "<<mat_conservation<<endl;
  }

  RNG* get_rng(void) const { return m_RNG[0];}
  RNG* get_thread_rng(int thread) const { return m_RNG[thread];}

  //non const functions
  void next_time_step(void) {
    m_time += m_dt;
    if (m_dt*m_dt_mult < m_dt_max) m_dt*=m_dt_mult;
    if (m_time + m_dt > m_time_stop) m_dt = m_time_stop - m_time;
    m_step++;
    
    //make post-step quantities the pre-step quantities
    pre_census_E = post_census_E;
  }

  void set_pre_census_E(double _pre_census_E) {pre_census_E = _pre_census_E;}
  void set_post_census_E(double _post_census_E) {post_census_E = _post_census_E;}
  void set_pre_mat_E(double _pre_mat_E) {pre_mat_E = _pre_mat_E;}
  void set_post_mat_E(double _post_mat_E) {post_mat_E = _post_mat_E;}
  void set_emission_E(double _emission_E) {emission_E = _emission_E;}
  void set_source_E(double _source_E) {source_E = _source_E;}
  void set_absorbed_E(double _absorbed_E) {absorbed_E = _absorbed_E;}
  void set_exit_E(double _exit_E) {exit_E = _exit_E;}
  void set_census_photons(unsigned int _census_photons) {m_census_photons = _census_photons;}

  private:
  //time
  double m_dt; //!< Current time step size (sh)
  double m_time; //!< Current time (sh)
  double m_time_stop; //!< End time (sh)
  unsigned int m_step; //!< Time step (start at 1)
  double m_dt_mult; //!< Time step multiplier
  double m_dt_max;  //!< Max time step

  //conservation
  double pre_census_E; // !< Census energy at the beginning of the timestep
  double post_census_E; //!< Census energy at the end of the timestep
  double pre_mat_E; //!< Material energy at the beginning of timestep
  double post_mat_E; //!< Material energy at the end of timestep
  double emission_E; //!< Energy emitted this timestep
  double exit_E; //!< Energy exiting problem
  double absorbed_E; //!< Total absorbed energy 
  double source_E; //!< Sourced energy

  //photons
  unsigned int m_total_photons;
  unsigned int m_census_photons;

  //RNG
  vector<RNG*> m_RNG;
};

#endif
