/*
  Author: Alex Long
  Date: 12/2/2015
  Name: source.h
*/
#ifndef source_h_
#define source_h_

#include <vector>

#include "constants.h"
#include "imc_state.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"

class Source {

  public:
  Source( Mesh *_mesh, 
          IMC_State *imc_state, 
          const double& total_E, 
          std::vector<Photon>& _census_photons)
    : mesh(_mesh),
      census_photons(_census_photons),
      icell(0),
      icensus(0),
      cell(mesh->get_cell(icell)),
      n_p_in_cell(0)
  {
    using std::vector;
    n_cell =mesh->get_number_of_objects();

    E_cell_emission = mesh->get_emission_E();
    E_cell_source = mesh->get_source_E();

    n_p_cell_emission = vector<unsigned int>(n_cell,0);
    n_p_cell_census = vector<unsigned int>(n_cell,0);
    n_p_cell_source = vector<unsigned int>(n_cell,0);
    n_p_cell = vector<unsigned int>(n_cell,0);

    //user desired photons
    unsigned int user_photons = imc_state->get_total_step_photons();

    //reset the total number of photons before counting
    n_photon = 0;

    //first, count the census photons in each cell
    n_photon+=census_photons.size();
    unsigned int local_index;
    for (vector<Photon>::iterator iphtn =census_photons.begin(); 
      iphtn!=census_photons.end();iphtn++) {
      local_index = mesh->get_local_ID(iphtn->get_cell());
      n_p_cell_census[local_index]++;
    }

    // count the total number of photons 
    for (unsigned int i = 0; i<n_cell; i++) {
      //emission 
      if (E_cell_emission[i] > 0.0) {
        unsigned int t_num_emission = 
          int(user_photons*E_cell_emission[icell]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        n_photon+=t_num_emission;
        n_p_cell_emission[i] = t_num_emission;
      }
      //source
      if (E_cell_source[i] > 0.0) {
        unsigned int t_num_source = 
          int(user_photons*E_cell_source[icell]/total_E);
        if (t_num_source == 0) t_num_source =1;
        n_photon+=t_num_source; 
        n_p_cell_source[i] = t_num_source;
      }
      n_p_cell[i] = 
        n_p_cell_census[i] + n_p_cell_emission[i] + n_p_cell_source[i];
    }
    emission_phtn_E = E_cell_emission[icell]/n_p_cell_emission[icell];
  }
  ~Source() {}

  unsigned int get_n_photon(void) const {return n_photon;}

  Photon get_photon(RNG *rng, const double& dt) {
    Photon return_photon;
    //check to increment cell
    if (n_p_in_cell == n_p_cell[icell]) {
      icell++; 
      cell = mesh->get_cell(icell);
      //set new photon energy for this cell
      emission_phtn_E = E_cell_emission[icell]/n_p_cell_emission[icell];
      n_p_in_cell=0;
    }
    // get correct photon
    if (n_p_in_cell < n_p_cell_emission[icell]) {
      return_photon = get_emission_photon(cell, emission_phtn_E, dt, rng);
      n_p_in_cell++;
    }
    else {
      return_photon = census_photons[icensus];
      icensus++;
      n_p_in_cell++;
    }

    if (!mesh->on_processor(return_photon.get_cell()))
      std::cout<<"THIS IS BAD";

    return return_photon;
  }

  Photon get_emission_photon( Cell& cell, 
                              const double& phtn_E, 
                              const double& dt, 
                              RNG *rng) 
  {
    using Constants::c;
    double pos[3];
    double angle[3];
    cell.uniform_position_in_cell(rng, pos);
    get_uniform_angle(angle, rng);
    Photon emission_photon; 
    emission_photon.set_position(pos);
    emission_photon.set_angle(angle);
    emission_photon.set_E0(phtn_E);
    emission_photon.set_distance_to_census(c*dt);
    emission_photon.set_cell(cell.get_ID());
    return emission_photon;
  }
 
  private:
  const Mesh * const mesh; //!< Pointer to mesh (source cannot change Mesh)
  std::vector<Photon>&  census_photons; //!< Pointer to census photons
  unsigned int icell; //!< Index of current cell to generate photons from
  unsigned int icensus; //!< Index of current census photon to be returned
  Cell cell;  //!< Current cell
  unsigned int n_p_in_cell;
  unsigned int n_cell;
  double emission_phtn_E;  //!< Energy of emitted photon in current cell
  unsigned int n_photon;  //!< User requested number of photons to create
  std::vector<unsigned int> n_p_cell_emission;
  std::vector<unsigned int> n_p_cell_census;
  std::vector<unsigned int> n_p_cell_source;
  std::vector<unsigned int> n_p_cell;
  std::vector<double> E_cell_emission;
  std::vector<double> E_cell_source;
};

#endif // source_h_
