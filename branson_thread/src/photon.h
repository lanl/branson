/*
  Author: Alex Long
  Date: 7/18/2014
  Name: photon.h
*/

#ifndef photon_h_
#define photon_h_


#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "constants.h"

using std::pow;
using std::cout;
using std::endl;
using Constants::dir_type;
using Constants::X_POS; using Constants::X_NEG; 
using Constants::Y_POS; using Constants::Y_NEG; 
using Constants::Z_POS; using Constants::Z_NEG; 

class Photon
{
  public:
  Photon() 
    {}

  ~Photon(void) {}

  void move(const double& distance) { 
    m_pos[0] += m_angle[0]*distance;
    m_pos[1] += m_angle[1]*distance;
    m_pos[2] += m_angle[2]*distance;
    m_life_dx -=distance;
  }

  bool below_cutoff(const double& cutoff_fraction) const {
    return (m_E / m_E0 < cutoff_fraction); 
  }
  
  void set_element(const unsigned int& new_elem) { m_elem_ID = new_elem;}

  void set_E0(const double& E) { 
    m_E0 = E;
    m_E = E;
  }
  void set_E(const double& E) {m_E = E;}

  void set_census_flag(const bool& census_flag) {m_census_flag = census_flag;}
  void set_distance_to_census(const double& dist_remain) {m_life_dx = dist_remain;}
  void set_angle(double *angle) { m_angle[0] = angle[0]; m_angle[1] = angle[1]; m_angle[2] = angle[2];}
  void set_position(double *pos) { m_pos[0] = pos[0]; m_pos[1] = pos[1]; m_pos[2] = pos[2];}
  void set_dead(void) { m_census_flag = false;}

  void reflect(const unsigned int& surface_cross) {
    //reflect the photon over the surface it was crossing
    if (surface_cross == X_POS || surface_cross == X_NEG) m_angle[0] = -m_angle[0];
    else if (surface_cross == Y_POS || surface_cross == Y_NEG) m_angle[1] = -m_angle[1]; 
    else m_angle[2] = -m_angle[2]; 
  }

  //const functions
  unsigned int get_element(void) const { return m_elem_ID; }
  const double* get_position(void) const { return m_pos; }
  const double* get_angle(void) const { return m_angle; }
  double get_E(void) const { return m_E;}
  double get_E0(void) const { return m_E0;}
  bool get_census_flag(void) const {return m_census_flag;}
  double get_distance_remaining(void) const {return m_life_dx;}  

  void print_info(void ) const {
    cout<<"----Photon Info----\n";
    cout<<"position: "<<m_pos<<endl;//" "<<m_pos[1]<<" "<<m_pos[2]<<endl;
    cout<<"angle: "<<m_angle<<endl;//" "<<m_angle[1]<<" "<<m_angle[2]<<endl;
    cout<<"Energy: "<<m_E<<" , Initial energy: "<<m_E0<<endl;
    cout<<"Element ID: "<<m_elem_ID<<" , Census Flag: "<<m_census_flag<<endl;
  }

  //override great than operator to sort
  bool operator <(const Photon& compare) const {
    return m_elem_ID < compare.get_element();
  }
  
  bool operator()(const Photon& compare1, const Photon& compare2) const {
    return compare1.get_element() < compare2.get_element();
  }


  //member variables  
  private:
  double m_pos[3]; //!< photon position
  double m_angle[3]; //!< photon angle array

  double m_E; //!< current photon energy
  double m_E0; //!< photon energy at creation

  double m_elem_ID; //!< Element ID
  bool m_census_flag; //!< Flag for census, true if photon reaches census
  double m_life_dx; //!< Distance remaining this time step

  //private member functions
  private:

};

#endif
