/*
  Author: Alex Long
  Date: 7/18/2014
  Name: mesh.h
*/

#ifndef mesh_h_
#define mesh_h_

#include <iostream>
#include <vector>
#include <omp.h>
#include <numeric> 

#include "constants.h"
#include "imc_state.h"
#include "input.h"
#include "RNG.h"
#include "element.h"

using std::cout;
using std::endl;
using std::pow;
using std::vector;
using std::accumulate;

using Constants::c;
using Constants::a;
using Constants::dir_type;
using Constants::bc_type;
using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT;

class Mesh
{
  public:
  Mesh( Input *input)
    : m_x_elem(input->get_n_x_elements()),
      m_y_elem(input->get_n_y_elements()),
      m_z_elem(input->get_n_z_elements()),
      m_elem(m_x_elem*m_y_elem * m_z_elem),
      m_dx(input->get_dx()),
      m_dy(input->get_dy()),
      m_dz(input->get_dz())
  {
    //source setting
    //m_source_elem = input->get_source_element();
    //m_Ts[m_source_elem] = input->get_source_T();
    bc[X_POS] = input->get_bc(X_POS); bc[X_NEG] = input->get_bc(X_NEG);
    bc[Y_POS] = input->get_bc(Y_POS); bc[Y_NEG] = input->get_bc(Y_NEG);
    bc[Z_POS] = input->get_bc(Z_POS); bc[Z_NEG] = input->get_bc(Z_NEG);

    double x_c = 0.0; double y_c = 0.0; double z_c = 0.0;
    for (unsigned int k = 0; k<m_z_elem; k++) {
      y_c = 0.0;
      for (unsigned int j = 0; j<m_y_elem; j++) {
        x_c = 0.0;
        for (unsigned int i = 0; i<m_x_elem; i++) {
          element_list.push_back(new Element(x_c, x_c + m_dx, y_c, y_c + m_dy, z_c, z_c + m_dz));
          x_c += m_dx;
        }
        y_c += m_dy;
      }
      z_c += m_dz;
    }

    //set up connectivity of mesh
    Element *e; 
    unsigned int count = 0;
    for (unsigned int k = 0; k<m_z_elem; k++) {
      for (unsigned int j = 0; j<m_y_elem; j++) {
        for (unsigned int i = 0; i<m_x_elem; i++) { 
          e = element_list[count];
          if (i<(m_x_elem-1)) {e->set_neighbor( X_POS, count+1); e->set_bc(X_POS, ELEMENT);}
          else                {e->set_neighbor( X_POS, count); e->set_bc(X_POS, bc[X_POS]);} 

          if (i>0)            {e->set_neighbor( X_NEG, count-1); e->set_bc(X_NEG, ELEMENT);}
          else                {e->set_neighbor( X_NEG, count); e->set_bc(X_NEG, bc[X_NEG]);}

          if (j<(m_y_elem-1)) {e->set_neighbor(Y_POS, count+m_y_elem); e->set_bc(Y_POS, ELEMENT);}
          else                {e->set_neighbor(Y_POS, count); e->set_bc(Y_POS, bc[Y_POS]);}

          if (j>0)            {e->set_neighbor(Y_NEG, count-m_y_elem); e->set_bc(Y_NEG, ELEMENT);}
          else                {e->set_neighbor(Y_NEG, count); e->set_bc(Y_NEG, bc[Y_NEG]);}

          if (k<(m_z_elem-1)) {e->set_neighbor(Z_POS, count+m_x_elem*m_y_elem); e->set_bc(Z_POS, ELEMENT);}
          else                {e->set_neighbor(Z_POS, count); e->set_bc(Z_POS, bc[Z_POS]);}  

          if (k>0)            {e->set_neighbor(Z_NEG, count-m_x_elem*m_y_elem); e->set_bc(Z_NEG, ELEMENT);}
          else                {e->set_neighbor(Z_NEG, count); e->set_bc(Z_NEG, bc[Z_NEG]);}
          count++;
        }
      }
    }

    m_source_E = vector<double>(m_elem, 0.0);
    m_census_E = vector<double>(m_elem, 0.0);
    m_emission_E = vector<double>(m_elem, 0.0);
    m_T = vector<double>(m_elem, input->get_initial_Tm() );
    m_Tr = vector<double>(m_elem, input->get_initial_Tr() );
    m_Ts = vector<double>(m_elem, 0.0 );

    m_CV = input->get_CV();
    m_rho = input->get_rho();
    m_opA = input->get_opacity_A();
    m_opB = input->get_opacity_B();
    m_opC = input->get_opacity_C();
  }
  ~Mesh(void) {}

  //const functions
  vector<double> get_source_E_vector(void) const { return m_source_E;}
  vector<double> get_emission_E_vector(void) const {return m_emission_E;}
  vector<double> get_census_E_vector(void) const {return m_census_E;}

  double get_total_source_E(void) const {return accumulate(m_source_E.begin(),m_source_E.end(),0);}
  double get_total_emission_E(void) const {return accumulate(m_emission_E.begin(),m_emission_E.end(),0);}
  double get_total_census_E(void) const {return accumulate(m_census_E.begin(),m_census_E.end(),0);}
  unsigned int get_num_elems(void) const {return m_elem;}
  
  double get_total_photon_E(void) const {
    double tot_photon_E = 0.0;
    for(unsigned int i=0; i<m_elem; ++i)
      tot_photon_E += (m_source_E[i] + m_census_E[i] + m_emission_E[i]);
    return tot_photon_E;
  }

  void print_mesh_info(void) const {
    cout<<"---- Piecewise Constant Mesh ----"<<endl;
    cout<<"Element  T   op   emission   census   source"<<endl;
    for (unsigned int i=0; i<m_elem; ++i) {
      const Element *e = element_list[i];
      cout<<i<<" "<<m_T[i]<<" "
          <<e->get_op()<<" "<<m_emission_E[i]<<" "<<m_census_E[i]<<" "
          <<m_source_E[i]<<endl;
    }
  }

  //non-const functions
  void calculate_photon_energy(IMC_State* imc_s) {
    double dt = imc_s->get_dt();
    double op, f;
    unsigned int step = imc_s->get_step();
    double vol =m_dx*m_dy*m_dz;
    double tot_census_E = 0.0;
    double tot_emission_E = 0.0;
    double tot_source_E = 0.0;
    double pre_mat_E = 0.0;
    for (unsigned int i=0; i<m_elem;++i) {
      Element *e = element_list[i];
      op = m_opA + m_opB*pow(m_T[i], m_opC);
      f =1.0/(1.0 + dt*op*c*(4.0*a*pow(m_T[i],3)/m_CV));

      e->set_op(op);
      e->set_f(f);

      m_emission_E[i] = dt*vol*f*op*a*c*pow(m_T[i],4);
      if (step > 1) m_census_E[i] = 0.0;  
      else m_census_E[i] =vol*a*pow(m_Tr[i],4); 
      m_source_E[i] = (1.0/4.0)*dt*op*a*c*pow(m_Ts[i],4);

      pre_mat_E+=m_T[i]*m_CV*vol;
      tot_emission_E+=m_emission_E[i];
      tot_census_E  +=m_census_E[i];
      tot_source_E  +=m_source_E[i];
    }

    //set conservation
    imc_s->set_pre_mat_E(pre_mat_E);
    imc_s->set_emission_E(tot_emission_E);
    imc_s->set_source_E(tot_source_E);
    if(imc_s->get_step() == 1) imc_s->set_pre_census_E(tot_census_E);

  }


  void update_temperature(vector<double>& abs_E, IMC_State* imc_s) {
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol =m_dx*m_dy*m_dz;
    for (unsigned int i=0; i<m_elem;++i) {
      m_T[i] = m_T[i] + (abs_E[i] - m_emission_E[i])/(m_CV*vol*m_rho);
      total_abs_E+=abs_E[i];
      total_post_mat_E+= m_T[i]*m_CV *vol;
 
      abs_E[i] = 0.0;
    }
    imc_s->set_absorbed_E(total_abs_E);
    imc_s->set_post_mat_E(total_post_mat_E);
  }

  Element* get_element(unsigned int elem_ID) {return element_list[elem_ID];}

  /*
  double get_total_mat_E(void) const
  double get_total_photon_E(void) const
  void print_mesh_info(void) const
  double uniform_position_in_elem(unsigned int elem_id, RNG *rng) const

  //non const functions
  void calculate_photon_energy(IMC_State* imc_state)
  double tally_scatter_event(unsigned int elem_ID)
  double tally_exit_event(unsigned int elem_ID)
  */

  //void tally_abs_E(const double& E_abs, const unsigned int& elem_id) {m_abs_E[elem_id] += E_abs;}

  private:

  // Geometry variables
  unsigned int m_x_elem;
  unsigned int m_y_elem;
  unsigned int m_z_elem;
  unsigned int m_elem;
  double m_dx;
  double m_dy;
  double m_dz;
  vector<Element*> element_list;
  bc_type bc[6];

  // Material variables
  double m_rho;
  double m_CV;
  double m_opA;
  double m_opB;
  double m_opC;

  vector<double> m_T;
  vector<double> m_Tr;
  vector<double> m_Ts; //< Source temperature
  vector<double> m_source_E;
  vector<double> m_census_E;
  vector<double> m_emission_E;
};

#endif
