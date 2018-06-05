//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   source.h
 * \author Alex Long
 * \date   December 2 2015
 * \brief  Allows transport function to create particles when needed
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef source_h_
#define source_h_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "constants.h"
#include "imc_state.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"
#include "work_packet.h"

//==============================================================================
/*!
 * \class Work_Local_Map
 * \brief A class that's used as a functor in sorting local vs. non-local work
 * packets
 *
 * \example no test yet
 */
//==============================================================================
class Work_Local_Map {

  public:

  //! Constructor
  Work_Local_Map(std::unordered_map<uint32_t, bool> & _work_map)
    : work_map(_work_map)
    {}
  //! Destructor
  ~Work_Local_Map() {}

  //! Parnthesis operator for sorting work packets based on global cell ID
  bool operator() (const Work_Packet& a, const Work_Packet& b) {
    return work_map[a.get_global_cell_ID()] < work_map[b.get_global_cell_ID()];
  }

  // member data
  //! Map global cell ID to on processor bool
  std::unordered_map<uint32_t, bool>& work_map; 
};

//==============================================================================
/*!
 * \class Source
 * \brief Describes how to create emssion particles and return census particles
 *
 * \example no test yet
 */
//==============================================================================
class Source {

  public:

  //! constructor
  Source( Mesh *_mesh, IMC_State *imc_s, const uint64_t user_photons, const double total_E,
    std::vector<Photon>& _census_photons)
    : mesh(_mesh),
      census_photons(_census_photons),
      iphoton(0)
  {
    using std::vector;

    uint32_t n_cell =mesh->get_n_local_cells();
    E_cell_emission = mesh->get_emission_E();
    E_cell_census = mesh->get_census_E();

    uint32_t step = imc_s->get_step();

    //reset the total number of photons before counting
    n_photon = 0;
    n_sourced = 0;

    double total_census_E = 0.0;

    //make work packets
    // current cell pointer
    const Cell* cell_ptr;
    for (uint32_t i = 0; i<n_cell; ++i) {
      Work_Packet temp_cell_work;
      cell_ptr = mesh->get_cell_ptr(i);
      //emission
      if (E_cell_emission[i] > 0.0) {
        uint32_t t_num_emission =
          int(user_photons*E_cell_emission[i]/total_E);
        // make at least one photon to represent emission energy
        if (t_num_emission == 0) t_num_emission =1;
        n_photon+=t_num_emission;
        // make work packet and add to vector
        temp_cell_work.set_global_cell_ID(cell_ptr->get_ID());
        temp_cell_work.set_global_grip_ID(cell_ptr->get_grip_ID());
        temp_cell_work.attach_creation_work(E_cell_emission[i], t_num_emission);
        temp_cell_work.set_coor(cell_ptr->get_node_array());
        temp_cell_work.set_source_type(Constants::EMISSION);
        work.push_back(temp_cell_work);
      }
      // initial census
      if (step ==1) {
        if (E_cell_census[i] > 0.0) {
          Work_Packet temp_cell_work;
          cell_ptr = mesh->get_cell_ptr(i);
          if (E_cell_census[i] > 0.0) {
            uint32_t t_num_census = int(user_photons*E_cell_census[i]/total_E);
            // make at least one photon to represent census energy
            if (t_num_census == 0) t_num_census =1;
            n_photon+=t_num_census;
            // make work packet and add to vector
            temp_cell_work.set_global_cell_ID(cell_ptr->get_ID());
            temp_cell_work.set_global_grip_ID(cell_ptr->get_grip_ID());
            temp_cell_work.attach_creation_work(E_cell_census[i], t_num_census);
            temp_cell_work.set_coor(cell_ptr->get_node_array());
            temp_cell_work.set_source_type(Constants::INITIAL_CENSUS);
            work.push_back(temp_cell_work);
            // keep track of census energy for conservation check
            total_census_E += E_cell_census[i];
          }
        }
        imc_s->set_pre_census_E(total_census_E);
      }
    }

    // add local census size to n_photon
    n_photon+=census_photons.size();
  }

  //! destructor
  ~Source() {}

  //! Get total photon count
  uint32_t get_n_photon(void) const {return n_photon;}

  //! Print source census balance
  void print_work_summary(const int& rank) const {
    std::cout<<"rank: "<<rank<<" emission: "<<n_photon-census_photons.size();
    std::cout<<" census: "<<census_photons.size();
    std::cout<<" fraction census: "<<census_photons.size()/double(n_photon);
    std::cout<<std::endl;
  }

  //! Link work packets to census particles and get new total
  void post_lb_prepare_source(void) {
    using std::vector;
    using std::unordered_map;

    // recount number of photons
    n_photon = 0;

    unordered_map<uint32_t, uint32_t> work_map;

    // map global cell ID to work packets to determine if census photons
    // can be attached to a work packet
    for (uint32_t i=0; i<work.size();++i)
    {
      Work_Packet& work_ref = work[i];
      work_map[work_ref.get_global_cell_ID()] = i;
      n_photon+=work_ref.get_n_particles();
    }

    uint32_t cell_ID, grip_ID;
    unordered_map<uint32_t, uint32_t> census_in_cell;
    unordered_map<uint32_t, uint32_t> census_start_index;
    unordered_map<uint32_t, uint32_t> grip_index;

    if (!census_photons.empty()) {
      uint32_t i=0;
      for (auto const &iphtn : census_photons)
      {
        cell_ID = iphtn.get_cell();
        grip_ID = iphtn.get_grip();
        // at first instance of cell ID set the starting index
        if (census_in_cell.find(cell_ID) == census_in_cell.end()) {
          census_in_cell[cell_ID] = 1;
          grip_index[cell_ID] = grip_ID;
          census_start_index[cell_ID]=i;
        }
        else {
          census_in_cell[cell_ID]++;
        }
        i++;
      }

      uint32_t work_index;
      // now attach census particles to work packets or make new ones
      for (auto const &imap : census_start_index)
      {
        cell_ID = imap.first;
        grip_ID = grip_index[cell_ID];

        // apppend count in this cell to the total
        n_photon += census_in_cell[cell_ID];

        // if a work packet for this cell exists, attach the census work
        if (work_map.find(cell_ID) != work_map.end()) {
          work_index = work_map[cell_ID];
          work[work_index].attach_census_work(imap.second,
            census_in_cell[cell_ID]);
        }
        // otherwise, make a new work packet
        else {
          Work_Packet temp_packet;
          temp_packet.set_global_cell_ID(cell_ID);
          temp_packet.set_global_grip_ID(grip_ID);
          temp_packet.attach_census_work(imap.second, census_in_cell[cell_ID]);
          // add this work to the total work vector
          work.push_back(temp_packet);
        }
      }
    } // end if !census_photons.empty()

    std::unordered_map<uint32_t, bool> on_proc_map;
    uint32_t work_cell_ID;
    for (auto const &iw : work) {
      work_cell_ID = iw.get_global_cell_ID();
      on_proc_map[work_cell_ID] = mesh->on_processor(work_cell_ID);
    }

    Work_Local_Map work_local_map(on_proc_map);
    std::sort(work.begin(), work.end(), work_local_map);

    // set initial parameters and iterators
    iwork = work.begin();
    n_create = iwork->get_n_create();
    phtn_E = iwork->get_photon_E();
    n_in_packet = iwork->get_n_particles();
    census_index = iwork->get_census_index();
    current_source = iwork->get_source_type();
    iphoton = 0;
    n_work= 1;
  }

  //! Get next particle, create it if it's an emission particle
  Photon get_photon(RNG *rng, const double& dt) {
    Photon return_photon;

    // get creation photon
    if (iphoton < n_create) {
      if (current_source == Constants::EMISSION)
        get_emission_photon(return_photon, *iwork, phtn_E, dt, rng);
      else if (current_source == Constants::INITIAL_CENSUS)
        get_initial_census_photon(return_photon, *iwork, phtn_E, dt, rng);
      iphoton++;
    }
    // get census photon
    else {
      return_photon = census_photons[census_index];
      iphoton++;
      census_index++;
    }

    // if work packet is done, increment the work packet and reset quantities
    // don't increment if this is the last photon
    if (iphoton == n_in_packet && n_work!=work.size()) {
      ++iwork;
      ++n_work;
      n_create = iwork->get_n_create();
      phtn_E = iwork->get_photon_E();
      n_in_packet = iwork->get_n_particles();
      census_index = iwork->get_census_index();
      current_source = iwork->get_source_type();
      iphoton = 0;
    }

    if (n_sourced >= n_photon) {
      std::cout<<"this is bad: can't source more than this"<<std::endl;
    }

    // increment count of sourced particles
    ++n_sourced;

    return return_photon;
  }

  //! Set input photon to the next emission photon
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
    emission_photon.set_grip(work.get_global_grip_ID());
  }

  //! Set input photon to the next intiial census photon
  void get_initial_census_photon( Photon& census_photon,
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
    census_photon.set_position(pos);
    census_photon.set_angle(angle);
    census_photon.set_E0(phtn_E);
    // initial census particles born at the beginning of timestep
    census_photon.set_distance_to_census(c*dt);
    census_photon.set_cell(work.get_global_cell_ID());
    census_photon.set_grip(work.get_global_grip_ID());
  }



  //! Return reference to work vector
  std::vector<Work_Packet>& get_work_vector(void) {return work;}

  private:
  const Mesh * const mesh; //!< Pointer to mesh (source cannot change Mesh)
  std::vector<Photon>& census_photons; //!< Reference to census photons on rank
  uint32_t n_create; //!< Number of emission particles created in this packet
  double phtn_E; //!< Photon emission energy in this packet
  uint32_t n_work; //!< Work packet counter
  uint32_t n_in_packet; //!< Number of total particles in this packet
  uint32_t census_index; //!< Index of next census particle to return
  uint32_t iphoton;  //!< Local photon counter
  uint32_t current_source; //!< Current source
  uint32_t n_sourced; //!< Number of particles returned by source
  std::vector<Work_Packet> work; //!< Vector of Work_Packet objects
  std::vector<Work_Packet>::iterator iwork; //!< Work iterator
  uint32_t n_photon;  //!< Total photons in this source
  std::vector<double> E_cell_emission; //!< Emission energy in each cell
  std::vector<double> E_cell_census; //!< Initial census energy in each cell
};

#endif // source_h_
//----------------------------------------------------------------------------//
// end of source.h
//----------------------------------------------------------------------------//
