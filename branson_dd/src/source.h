/*
  Author: Alex Long
  Date: 12/2/2015
  Name: source.h
*/
#ifndef source_h_
#define source_h_

#include "constants.h"
#include "imc_state.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"

class Source {

  public:
  Source(Mesh *_mesh, IMC_State *imc_state, const double& total_E)
    : mesh(_mesh),
      icell(0),
      cell(mesh->get_cell(icell)),
      n_p_in_cell(0)
  {
    using std::vector;
    n_cell =mesh->get_number_of_objects();

    E_cell_census = mesh->get_census_E_ref();
    E_cell_emission = mesh->get_emission_E_ref();
    E_cell_source = mesh->get_source_E_ref();

    n_cell_emission = vector<unsigned int>(n_cell,0);
    n_cell_census = vector<unsigned int>(n_cell,0);
    n_cell_source = vector<unsigned int>(n_cell,0);
    n_cell = vector<unsigned int>(n_cell,0);

    unsigned int total_phtns = imc_state->get_total_step_photons();
    dt = imc_state->get_dt();

    //reset the number of photons before counting
    n_photon = 0;

    // count the total number of photons 
    for (unsigned int i = 0; i<n_cell; i++) {
      if (E_cell_census[i] > 0.0) { 
        unsigned int t_num_census = int(total_phtns*E_cell_census[icell]/total_E);
        if (t_num_census == 0) t_num_census =1;
        n_photon+= t_num_census;
        n_cell_census[i] = t_num_census;
      }
      if (E_cell_emission[i] > 0.0) {
        unsigned int t_num_emission = int(total_phtns*E_cell_emission[icell]/total_E);
        if (t_num_emission == 0) t_num_emission =1;
        n_photon+=t_num_emission;
        n_cell_emission[i] = t_num_emission;
      }
      if (E_cell_source[i] > 0.0) {
        unsigned int t_num_source = int(total_phtns*source_E[icell]/total_E);
        if (t_num_source == 0) t_num_source =1;
        n_photon+=t_num_source; 
        n_cell_source[i] = t_num_source;
      }
      n_cell[i] = n_cell_census[i] + n_cell_emission[i] + n_cell_source[i];
    }
    emission_phtn_E = E_cell_emission[icell]/n_cell_emission[icell];
    census_phtn_E = E_cell_census[icell]/n_cell_census[icell];
  }
  ~Source() {}

  Photon get_photon(RNG *rng, const double& dt) {
    Photon return_photon;
    if (n_p_in_cell == n_cell[icell]) {
      icell++; 
      cell = mesh->get_cell(icell);
      n_p_in_cell=0;
    }
    if (n_p_in_cell < n_cell_emission[icell]) 
      return_photon = get_emission_photon(cell, emission_phtn_E, dt, rng);
    else 
      return_photon = get_census_photon(cell, census_phtn_E, dt, rng);
    return return_photon;
  }

  Photon get_census_photon( Cell& cell, 
                            const double& phtn_E, 
                            const double* dt, 
                            RNG *rng) 
  {
    using Constants::c;
    double pos[3];
    double angle[3];
    cell.uniform_position_in_cell(rng, pos);
    get_uniform_angle(angle, rng);
    Photon census_photon; 
    census_photon.set_position(pos);
    census_photon.set_angle(angle);
    census_photon.set_E0(phtn_E);
    census_photon.set_distance_to_census(c*delta_t);
    census_photon.set_cell(cell.get_ID());
    return census_photon;
  }

  Photon get_emission_photon( Cell& cell, 
                              const double& phtn_E, 
                              const double* dt, 
                              RNG *rng) 
  {
    using Constants::c;
    double pos[3];
    double angle[3];
    cell.uniform_position_in_cell(rng, pos);
    get_uniform_angle(angle, rng);
    Photon emission_photon; 
    census_photon.set_position(pos);
    census_photon.set_angle(angle);
    census_photon.set_E0(phtn_E);
    census_photon.set_distance_to_census(c*delta_t);
    census_photon.set_cell(cell.get_ID());
    return emission_photon;
  }
 
  private:
  const Mesh const *mesh;
  unsigned int icell;
  Cell * cell; 
  unsigned int n_p_in_cell;
  unsigned int n_cell;
  double dt;
  double census_phtn_E;
  double emission_phtn_E;
  unsigned int n_photon;
  vector<unsigned int> n_cell_emission;
  vector<unsigned int> n_cell_census;
  vector<unsigned int> n_cell_source;
  vector<unsigned int> n_cell;
  vector<double> E_cell_emission;
  vector<double> E_cell_census;
  vector<double> E_cell_source;
};

#endif // source_h_

