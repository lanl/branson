/*
  Author: Alex Long
  Date: 1/15/2016
  Name: work_packet.h
*/
#ifndef work_packet_h_
#define work_packet_h_

#include "RNG.h"

//==============================================================================
/*!
 * \class Work_Packet
 * \brief Holds data to make particles 
 * 
 * A way to pass work between processors in mesh passing mode that avoids 
 * constructing the emission particles
 */
//==============================================================================
class Work_Packet {

  public:
  Work_Packet()
    : n_particles(0),
      g_cell_ID(0),
      n_census(0),
      census_index_start(0),
      emission_E(0.0)
      {}
  ~Work_Packet() {}

  // constant funtions
  uint32_t get_global_cell_ID(void) const {return g_cell_ID;}
  uint32_t get_n_particles(void) const {return n_particles;}
  uint32_t get_n_census(void) const {return n_census;}
  uint32_t get_census_index(void) const {return census_index_start;}
  double get_emission_E(void) const {return emission_E;}
  double get_photon_E(void) const {return emission_E/n_particles;}

  void uniform_position_in_cell(RNG* rng, double* pos) const {
    pos[0]= nodes[0] + rng->generate_random_number()*(nodes[1]-nodes[0]);
    pos[1]= nodes[2] + rng->generate_random_number()*(nodes[3]-nodes[2]);
    pos[2]= nodes[4] + rng->generate_random_number()*(nodes[5]-nodes[4]);
  }

  const double* get_node_array(void) const {return nodes;}

  // non-const functions
  void set_global_cell_ID(const uint32_t& _global_cell_ID) {
    g_cell_ID = _global_cell_ID;
  }
  void attach_emission_work(const double& _emission_E, 
    const uint32_t& _n_particles) 
  {
    emission_E = _emission_E;
    n_particles = _n_particles;
  }

  void attach_census_work(const uint32_t& _census_index_start, 
    const uint32_t& _n_census) 
  {
    census_index_start=_census_index_start;
    n_census = _n_census;
  }

  void set_coor(const double *cell_nodes) {
    nodes[0] = cell_nodes[0];
    nodes[1] = cell_nodes[1];
    nodes[2] = cell_nodes[2];
    nodes[3] = cell_nodes[3];
    nodes[4] = cell_nodes[4];
    nodes[5] = cell_nodes[5];
  }

  Work_Packet split(const uint32_t& n_remain) {
    Work_Packet return_work;
    
    // calculate properies of return packet
    uint32_t n_return = n_particles - n_remain;
    double return_E = emission_E*double(n_return)/double(n_particles);
    
    // set properties of this work packet
    emission_E = emission_E - return_E;
    n_particles = n_remain;

    // set properties of return work packet
    return_work.set_global_cell_ID(g_cell_ID);
    return_work.set_coor(nodes);
    return_work.attach_emission_work(return_E, n_return);
    return return_work;
  }

  private:
  uint32_t n_particles; //!< Total number of particles in work packet
  uint32_t g_cell_ID; //!< Global index of cell containing this work
  uint32_t n_census; //!< Census photons in this cell
  uint32_t census_index_start; //!< Start index in census vector
  double emission_E; //!< Emission energy in work packet
  double nodes[6]; //!< Nodes forming 3D cell
};

#endif // work_packet_h_
