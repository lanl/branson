/*
  Author: Alex Long
  Date: 12/1/2015
  Name: mesh_particle_pass.h
*/

#ifndef mesh_particle_pass_h_
#define mesh_particle_pass_h_

#include <vector>
#include <string>
#include <algorithm>

#include "mesh.h"
#include "input.h"


class Mesh_Particle_Pass :: Mesh {

  public:
  Mesh_Particle_Pass(Input input, unsigned int _rank, unsigned int _nrank) 
    : Mesh(input, _rank, _nrank)
  {

    // JAYENNE code used to generate numbering in binary tree
    // compute parent and child node ids, ignoring missing nodes
  }
  ~Mesh_Particle_Pass(void) {}
  
  private:
  int parent; //!< Rank ID of parent in binary tree
  int child1; //!< Rank ID of first child, possibly null
  int child2; //!< Rank ID of second child, possibly null
  unsigned int n_complete; //!< Completed particles, includes child nodes
  unsigned int d_n_pass_messages; //!< Number of particles messages sent (diagnostic)
  unsigned int d_n_pass; //!< Number of particles sent (diagnostic)

  std::vector<Photon> recv_particles; //!< Receive particle buffers
  std::vector<Photon> send_particle; //!< Particle send buffer
  bool all_finished; //!< Finished flag for all ranks
};

#endif // mesh_particle_pass_h_
