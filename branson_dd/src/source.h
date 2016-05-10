/*
  Author: Alex Long
  Date: 12/2/2015
  Name: source.h
*/
#ifndef source_h_
#define source_h_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "constants.h"
#include "imc_parameters.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"
#include "work_packet.h"

class Source {

  public:
  Source( Mesh *_mesh, IMC_Parameters *imc_parameters, const double& total_E,
    std::vector<Photon>& _census_photons)
    : mesh(_mesh),
      census_photons(_census_photons),
      iphoton(0)
  {
    using std::vector;

    uint32_t n_cell =mesh->get_n_local_cells();
    E_cell_emission = mesh->get_emission_E();

    //user desired photons
    uint64_t user_photons = imc_parameters->get_n_user_photon();

    //reset the total number of photons before counting
    n_photon = 0;

    uint32_t global_index;

    //make work packets
    Work_Packet temp_cell_work;
    // current cell pointer
    const Cell* cell_ptr;
    for (uint32_t i = 0; i<n_cell; i++) {
      cell_ptr = mesh->get_cell_ptr(i);
      //emission 
      if (E_cell_emission[i] > 0.0) {
        uint32_t t_num_emission = 
          int(user_photons*E_cell_emission[i]/total_E);
        // make at least one photon to represent emission energy
        if (t_num_emission == 0) t_num_emission =1; 
        n_photon+=t_num_emission;
        // make work packet and add to vector
        global_index =mesh->get_global_ID(i); 
        temp_cell_work.set_global_cell_ID(global_index);
        temp_cell_work.attach_emission_work(E_cell_emission[i], t_num_emission);
        temp_cell_work.set_coor(cell_ptr->get_node_array());
        work.push_back(temp_cell_work);
      }
    }

    // add local census size to n_photon
    n_photon+=census_photons.size();

    // local and total are the same

    // total photons and local photons are the same 
    // now, count the census photons in each cell
    // increment the n_p_in_cell count
    
    /*
    for (vector<Photon>::iterator iphtn =census_photons.begin(); 
      iphtn!=census_photons.end();iphtn++) {
      global_index = iphtn->get_cell();
      n_p_in_cell[global_index]++;
    } 
    */
  }
  ~Source() {}


  uint32_t get_n_photon(void) const {return n_photon;}


  void post_lb_prepare_source(void) {
    using std::vector;
    using std::unordered_map;
    
    //recount number of photons
    n_photon = 0;

    unordered_map<uint32_t, Work_Packet*> work_map;
    // map global cell ID to work packets to determine if census photons 
    // can be attached to a work packet
    for (vector<Work_Packet>::iterator work_itr =work.begin(); 
      work_itr!=work.end();work_itr++) 
    {
      work_map[work_itr->get_global_cell_ID()] = &(*work_itr);
      n_photon+=work_itr->get_n_particles();
    }

    Work_Packet temp_packet;

    if (!census_photons.empty()) {
      uint32_t last_g_ID = census_photons.front().get_cell();
      uint32_t current_g_ID=last_g_ID;
      uint32_t census_in_cell=0;
      uint32_t i=0; //! Used to mark the loop index
      uint32_t start_index = 0;
      for (vector<Photon>::iterator iphtn =census_photons.begin(); 
        iphtn!=census_photons.end();iphtn++) 
      {
        current_g_ID = iphtn->get_cell();
        // if new global ID, try to attach to work packet to the last global index
        if (last_g_ID != current_g_ID) {
          // if work_packet exists for the last cell, attach
          if (work_map.find(last_g_ID) != work_map.end())
            work_map[last_g_ID]->attach_census_work(start_index, census_in_cell);
          // otherwise make a work packet and append it
          else {
            temp_packet.set_global_cell_ID(last_g_ID);
            temp_packet.attach_census_work(start_index, census_in_cell);
            // add this work to the total work vector
            work.push_back(temp_packet);
          }
          // update indices, reset count
          last_g_ID = current_g_ID;
          start_index = i;
          n_photon+=census_in_cell;
          census_in_cell = 0;
        }
        //increment counters
        census_in_cell++;
        i++;
      }

      //process work from last group of census particles
      if (work_map.find(current_g_ID) != work_map.end())
        work_map[last_g_ID]->attach_census_work(start_index, census_in_cell);
      // otherwise make a work packet and append it
      else {
        temp_packet.set_global_cell_ID(last_g_ID);
        temp_packet.attach_census_work(start_index, census_in_cell);
        // add this work to the total work vector
        work.push_back(temp_packet);
      }
      // add last group of census particles to count
      n_photon+=census_in_cell;
    } // end if census not empty

    // set initial parameters and iterators
    iwork = work.begin();
    n_emission = iwork->get_n_particles();
    phtn_E = iwork->get_photon_E();
    n_in_packet = iwork->get_n_particles() + iwork->get_n_census();
    census_index = iwork->get_census_index();
    iphoton = 0;
  }


  Photon get_photon(RNG *rng, const double& dt) {
    Photon return_photon;
    
    // get emission photon
    if (iphoton < n_emission) {
      get_emission_photon(return_photon, *iwork, phtn_E, dt, rng);
      iphoton++;
    }
    // get census photon
    else {
      return_photon = census_photons[census_index];
      iphoton++;
      census_index++;
    }
 
    // if work packet is done, increment the work packet and reset quantities
    if (iphoton == n_in_packet) {
      iwork++;
      n_emission = iwork->get_n_particles();
      phtn_E = iwork->get_photon_E();
      n_in_packet = iwork->get_n_particles() + iwork->get_n_census();
      census_index = iwork->get_census_index();
      iphoton = 0;
    }

    return return_photon;
  }

  void get_emission_photon( Photon& emission_photon,
                              Work_Packet& work, 
                              const double& phtn_E, 
                              const double& dt, 
                              RNG *rng) 
  {
    using Constants::c;
    double pos[3];
    double angle[3];
    work.uniform_position_in_cell(rng, pos);
    get_uniform_angle(angle, rng);
    emission_photon.set_position(pos);
    emission_photon.set_angle(angle);
    emission_photon.set_E0(phtn_E);
    emission_photon.set_distance_to_census(rng->generate_random_number()*c*dt);
    emission_photon.set_cell(work.get_global_cell_ID());
  }

  std::vector<Work_Packet>& get_work_vector(void) {return work;}
 
  private:
  const Mesh * const mesh; //!< Pointer to mesh (source cannot change Mesh)
  std::vector<Photon>& census_photons; //!< Reference to census photons on rank
  uint32_t n_emission; //!< Number of emission particles created in this packet
  double phtn_E; //!< Photon emission energy in this packet
  uint32_t n_in_packet; //!< Number of total particles in this packet
  uint32_t census_index;
  uint32_t iphoton;  //!< Local photon counter
  std::vector<Work_Packet> work; //!< Work packets
  std::vector<Work_Packet>::iterator iwork; //!< Work iterator
  uint32_t n_photon;  //!< Total photons in this source
  std::vector<double> E_cell_emission;
};

#endif // source_h_
