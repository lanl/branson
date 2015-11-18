/*
  Author: Alex Long
  Date: 9/18/2014
  Name: mesh.h
*/

#ifndef mesh_h_
#define mesh_h_

#include <iostream>
#include <vector>

#include "constants.h"
#include "imc_state.h"
#include "input.h"
#include "RNG.h"

using std::cout;
using std::endl;
using std::pow;
using std::vector;

using Constants::c;
using Constants::a;

class Mesh
{
  public:
  Mesh(unsigned int elements, double x_start, double dx)
    : m_elements(elements),
      m_x_start(x_start),
      m_dx(dx)
    {}
  ~Mesh(void) {}

  //const functions
  virtual const vector<double>& get_source_E_vector(void) const = 0;
  virtual const vector<double>& get_census_E_vector(void) const = 0;
  virtual const vector<double>& get_emission_E_vector(void) const = 0;
  virtual unsigned int get_num_elems(void) const = 0;
  virtual double get_total_mat_E(void) const = 0;
  virtual double get_total_photon_E(void) const = 0;
  virtual void print_mesh_info(void) const = 0;
  virtual double uniform_position_in_elem(unsigned int elem_id, RNG *rng) const = 0;

  //non const functions
  virtual void calculate_photon_energy(IMC_State* imc_state) = 0;
  //virtual double get_uniform_angle(void) = 0;

  protected:
  unsigned int m_elements;
  unsigned int m_source_elem;
  double m_CV;
  double m_opA;
  double m_opB;
  double m_opC;
  double m_dx;
  double m_x_start;
  vector<double> m_T;
  vector<double> m_Tr;
  vector<double> m_Ts;
  vector<double> m_source_E;
  vector<double> m_census_E;
  vector<double> m_emission_E;
  vector<double> op;
  vector<double> f;
};


class PC_Mesh : public Mesh
{
  public:
  PC_Mesh(Input *input)  
    : Mesh(input->get_n_elements(), input->get_x_start(), input->get_dx())
  {
    m_T = vector<double>(m_elements, input->get_initial_Tm());
    m_Tr = vector<double>(m_elements, input->get_initial_Tr());
    m_Ts = vector<double>(m_elements, 0.0);
    op = vector<double>(m_elements, 0.0); 
    f = vector<double>(m_elements, 0.0); 

    //source setting
    m_source_elem = input->get_source_element();
    m_Ts[m_source_elem] = input->get_source_T();

    m_source_E = vector<double>(m_elements, 0.0);
    m_census_E = vector<double>(m_elements, 0.0);
    m_emission_E = vector<double>(m_elements, 0.0);

    m_CV = input->get_CV();
    m_opA = input->get_opacity_A();
    m_opB = input->get_opacity_B();
    m_opC = input->get_opacity_C();
  }
  ~PC_Mesh(void) {}
  
  virtual const vector<double>& get_source_E_vector(void) const { return m_source_E;}
  virtual const vector<double>& get_emission_E_vector(void) const {return m_source_E;}
  virtual const vector<double>& get_census_E_vector(void) const {return m_census_E;}
  virtual unsigned int get_num_elems(void) const {return m_elements;}

  virtual double get_total_mat_E(void) const {
    double tot_E = 0.0;
    for(unsigned int i=0; i<m_elements; ++i)
      tot_E += m_T[i]*m_CV*m_dx;
    return tot_E;
  }
  
  virtual double get_total_photon_E(void) const {
    double tot_photon_E = 0.0;
    for(unsigned int i=0; i<m_elements; ++i)
      tot_photon_E += (m_source_E[i] + m_census_E[i] + m_emission_E[i]);
    return tot_photon_E;
  }
  
  virtual void print_mesh_info(void) const {
    cout<<"---- Piecewise Constant Mesh ----"<<endl;
    cout<<"Element  x   T   op   emission   census   source"<<endl;
    for (unsigned int i=0; i<Mesh::m_elements; ++i) {
      cout<<i<<" "<<Mesh::m_x_start + m_dx*i+m_dx*0.5<<" "<<m_T[i]<<" "
          <<op[i]<<" "<<m_emission_E[i]<<" "<<m_census_E[i]<<" "
          <<m_source_E[i]<<endl;
    }
  }

  virtual double uniform_position_in_elem(unsigned int elem_id, RNG* rng) const {
    return elem_id*m_dx + m_dx*rng->generate_random_number();
  } 


  virtual void calculate_photon_energy(IMC_State* imc_s) {
    double dt = imc_s->get_dt();
    unsigned int step = imc_s->get_step();
    for (unsigned int i=0; i<m_elements;++i) {
      op[i] = m_opA + m_opB*pow(m_T[i], m_opC);
      f[i] = 1.0/(1.0 + dt*op[i]*c*(4.0*a*pow(m_T[i],3)/m_CV));
      m_emission_E[i] = dt*m_dx*f[i]*op[i]*a*c*pow(m_T[i],4);
      if (step > 1) m_census_E[i] = 0.0;  
      else m_census_E[i] =m_dx*a*pow(m_Tr[i],4); 
      m_source_E[i] = (1.0/4.0)*dt*op[i]*a*c*pow(m_Ts[i],4);
    }
  }


};

#endif

