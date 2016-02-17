/*
  Author: Alex Long
  Date: 1/15/2016
  Name: work_packet.h
*/
#ifndef work_packet_h_
#define work_packet_h_

#include <vector>
#include "Photon.h"

//==============================================================================
/*!
 * \class Work_Packet
 * \brief Holds data to make photons and run census photons 
 * 
 * A way to pass work between processors in mesh passing mode that avoids 
 * constructing the emission and source particles
 */
//==============================================================================
class Work_Packet {

  public:
  Work_Packet() {}
  ~Work_Packet() {}

  void set_cell(uint32_t _global_cell_ID) {
    global_cell_ID = _global_cell_ID;
  }
 
  void attach_census_photons(vector<Photon> _census_photons {
    census_photons = _census_photons;
    n_census_photons= census_photons.size();
    n_photon += n_census_photons;
  }

  void attach_emission_work(double _emission_E, uint32_t _n_photons) {
    emission_E = _emission_E;
    _n_photons = n_photons;
  }

  private:
  uint32_t n_photon; //!< Total number of particles in work packet
  uint32_t n_census_photon; //!< Number of census particles in work packet
  uint32_t n_emission_photon; //!< Number of emission particles in work packet
  uint32_t n_source_photon; //!< Number of emission particles in work packet
  uint32_t global_cell_ID; //!< Index of cell containing this work
  int rank; //!< Rank that owns the mesh data for this work packet 
  double emission_E; //!< Emission energy in work packet
  vector<Photon> census_photons; //!< Census photons in this work packet

};

#endif // work_packet_h_
