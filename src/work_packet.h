//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   work_packet.h
 * \author Alex Long
 * \date   January 15 2016
 * \brief  Class that describes work in a single cell
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//---------------------------------------------------------------------------//

#ifndef work_packet_h_
#define work_packet_h_

#include "constants.h"
#include "RNG.h"

//==============================================================================
/*!
 * \class Work_Packet
 * \brief Holds data to make particles
 *
 * A way to pass work between processors in mesh passing mode that avoids
 * constructing the emission and initial census particles. In this class I
 * use the word "creation" to describe generic particles that can be created
 * at a later time (emission, initial census or source).
 */
//==============================================================================
class Work_Packet {

  public:
  Work_Packet()
    : n_particles(0),
      g_cell_ID(0),
      g_grip_ID(0),
      n_create(0),
      n_census(0),
      source_type(0),
      census_index_start(0),
      create_E(0.0)
      {}
  ~Work_Packet() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Return global cell ID where particles will be created
  uint32_t get_global_cell_ID(void) const {return g_cell_ID;}

  //! Return global grip ID where particles will be created
  uint32_t get_global_grip_ID(void) const {return g_grip_ID;}

  //! Return the number of particles in this work packet
  uint32_t get_n_particles(void) const {return n_particles;}

  //! Return the number of to be created or attached census particles
  uint32_t get_n_census(void) const {return n_census;}

  //! Return the total number of to be created particles in this work packet
  uint32_t get_n_create(void) const {return n_create;}

  //! Return the starting index of photons in this work packet
  uint32_t get_census_index(void) const {return census_index_start;}

  //! Return the total energy to be created (emission or initial census)
  double get_create_E(void) const {return create_E;}

  //! Return the energy of an emission/initial census photon in this work packet
  double get_photon_E(void) const {return create_E/n_create;}

  //! Return the source type (either EMISSION or INITIAL_CENSUS)
  uint32_t get_source_type(void) const {return source_type;}

  //! Set an input array to a random position within a cell
  void uniform_position_in_cell(RNG* rng, double* pos) const {
    pos[0]= nodes[0] + rng->generate_random_number()*(nodes[1]-nodes[0]);
    pos[1]= nodes[2] + rng->generate_random_number()*(nodes[3]-nodes[2]);
    pos[2]= nodes[4] + rng->generate_random_number()*(nodes[5]-nodes[4]);
  }

  //! Return a pointer to the start of the node array
  const double* get_node_array(void) const {return nodes;}

  //! Override greater than operator to sort work packets by the number of
  // particles in this
  bool operator <(const Work_Packet& compare) const {
    return  n_particles <  compare.get_n_particles() ;
  }

  // non-const functions

  //! Set the global cell ID
  void set_global_cell_ID(const uint32_t& _global_cell_ID) {
    g_cell_ID = _global_cell_ID;
  }

  //! Set the global grip ID
  void set_global_grip_ID(const uint32_t& _global_grip_ID) {
    g_grip_ID = _global_grip_ID;
  }

  //! Attach emission or initial census energy and particles to this work packet
  void attach_creation_work(const double _create_E,
    const uint32_t _n_create)
  {
    create_E = _create_E;
    n_create = _n_create;
    n_particles += _n_create;
  }

  //! Attach existing work in the form of index into large census particle array
  void attach_census_work(const uint32_t _census_index_start,
    const uint32_t _n_census)
  {
    census_index_start=_census_index_start;
    n_census = _n_census;
    n_particles += _n_census;
  }

  //! Set the particle type to be created by this work packet
  void set_source_type(uint32_t _source_type) {
    source_type=_source_type;
  }

  //! Set the coordinates of the work packet
  void set_coor(const double *cell_nodes) {
    nodes[0] = cell_nodes[0];
    nodes[1] = cell_nodes[1];
    nodes[2] = cell_nodes[2];
    nodes[3] = cell_nodes[3];
    nodes[4] = cell_nodes[4];
    nodes[5] = cell_nodes[5];
  }

  //! Split the work packet and return the split off work packet
 
  //! The returned Work_Packet that was split off will have 
  //! n_particles - n_remain in it
  Work_Packet split(const uint32_t& n_remain) {
    Work_Packet return_work;

    // calculate properties of return packet
    uint32_t n_return = n_particles - n_remain;
    double return_E = create_E*double(n_return)/double(n_particles);

    // set properties of this work packet
    create_E = create_E - return_E;

    // split work packets have no census particles yet
    n_particles = n_remain;
    n_create = n_remain;

    // set properties of return work packet
    return_work.set_global_cell_ID(g_cell_ID);
    return_work.set_global_grip_ID(g_grip_ID);
    return_work.set_coor(nodes);
    return_work.attach_creation_work(return_E, n_return);
    return return_work;
  }

  private:
  uint32_t n_particles; //!< Total number of particles in work packet
  uint32_t g_cell_ID; //!< Global index of cell containing this work
  uint32_t g_grip_ID; //!< Global index of this cell's grip
  uint32_t n_create; //!< Photons to create in this cell
  uint32_t n_census; //!< Census photons in this cell
  uint32_t source_type; //!< Type of particles, census or emission
  uint32_t census_index_start; //!< Start index in census vector
  double create_E; //!< Emission or intial census energy in work packet
  double nodes[6]; //!< Nodes forming 3D cell
};

#endif // work_packet_h_
//---------------------------------------------------------------------------//
// end of work_packet.h
//---------------------------------------------------------------------------//
