//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh.h
 * \author Alex Long
 * \date   July 18 2014
 * \brief  Object that holds mesh and manages decomposition and communication
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef mesh_h_
#define mesh_h_

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cell.h"
#include "constants.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "info.h"
#include "input.h"
#include "mpi_types.h"
#include "proto_cell.h"
#include "proto_mesh.h"
#include "replicate_mesh.h"
#include "timer.h"

//==============================================================================
/*!
 * \class Mesh
 * \brief Manages data access, decomposition and parallel communication for mesh
 *
 * Using an Input class, make the mesh with the correct material properties
 * for each region. The mesh numbering and mapping between global IDs and local
 * indices are all determined with the aid of Metis in the decompose_mesh
 * function. The mesh class also manages two-sided messaging in the mesh-
 * passing method.
 *
 */
//==============================================================================
class Mesh {

public:
  //! constructor
  Mesh(const Input &input, const MPI_Types &mpi_types, const Info &mpi_info,
       const IMC_Parameters &imc_p)
      : ngx(input.get_global_n_x_cells()), ngy(input.get_global_n_y_cells()),
        ngz(input.get_global_n_z_cells()), n_global(ngx * ngy * ngz),
        rank(mpi_info.get_rank()), n_rank(mpi_info.get_n_rank()),
        total_photon_E(0.0), off_rank_reads(0), silo_x(input.get_silo_x()),
        silo_y(input.get_silo_y()), silo_z(input.get_silo_z()),
        regions(input.get_regions()) {
    using Constants::bc_type;
    using Constants::ELEMENT;
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    using Constants::Z_NEG;
    using Constants::Z_POS;
    using std::vector;

    Proto_Mesh proto_mesh(input, mpi_types, mpi_info);

    // mode is replicated
    replicate_mesh(proto_mesh, mpi_types, mpi_info);
    // get decomposition information from proto mesh
    off_rank_bounds = proto_mesh.get_off_rank_bounds();
    on_rank_start = off_rank_bounds.front();
    on_rank_end = off_rank_bounds.back() - 1;
    const std::vector<Proto_Cell> &proto_cell_list(proto_mesh.get_cell_list());

    // this rank's cells
    n_cell = proto_cell_list.size();
    // size physics data
    m_census_E.resize(n_cell);
    m_emission_E.resize(n_cell);
    m_source_E.resize(n_cell);
    T_r.resize(n_cell);

    // for replicated mode, set the factor that reduces the emission
    // energy and initial census energy
    replicated_factor = 1.0 / n_rank;

    // don't bother with the MPI window
    cells = new Cell[proto_cell_list.size()];
    // use the proto cells to contstruct the real cells
    int i = 0;
    for (auto icell : proto_cell_list) {
      cells[i] = Cell(icell);
      i++;
    }

    // map region IDs to index in the region
    for (uint32_t i = 0; i < regions.size(); i++)
      region_ID_to_index[regions[i].get_ID()] = i;
  }

  // destructor
  ~Mesh() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//
  uint32_t get_n_local_cells(void) const { return n_cell; }
  uint32_t get_my_rank(void) const { return rank; }
  uint32_t get_offset(void) const { return on_rank_start; }
  uint32_t get_global_num_cells(void) const { return n_global; }

  double get_total_photon_E(void) const { return total_photon_E; }

  std::vector<uint32_t> get_off_rank_bounds(void) { return off_rank_bounds; }

  uint32_t get_grip_ID_from_cell_ID(uint32_t cell_ID) const {
    uint32_t local_ID = cell_ID - on_rank_start;
    return cells[local_ID].get_grip_ID();
  }

  void print(void) const {
    for (uint32_t i = 0; i < n_cell; i++)
      cells[i].print();
  }

  Cell get_cell(const uint32_t &local_ID) const { return cells[local_ID]; }

  const Cell *get_cell_ptr(const uint32_t &local_ID) const {
    return &cells[local_ID];
  }

  const Cell *get_const_cells_ptr(void) const { return cells; }

  uint32_t get_local_ID(const uint32_t &index) const {
    return index - on_rank_start;
  }

  uint32_t get_global_ID(const uint32_t &local_index) const {
    return on_rank_start + local_index;
  }

  Cell get_on_rank_cell(const uint32_t index) const {
    // this can only be called after with valid cell index (on rank or in stored
    // cells vector
    return cells[index - on_rank_start];
  }

  bool on_processor(const uint32_t &index) const {
    return (index >= on_rank_start) && (index <= on_rank_end);
  }

  std::vector<double> get_census_E(void) const { return m_census_E; }
  std::vector<double> get_emission_E(void) const { return m_emission_E; }
  std::vector<double> get_source_E(void) const { return m_source_E; }

  //! Get the radiation temperature in a cell (for plotting/diagnostics)
  double get_T_r(const uint32_t cell_index) const { return T_r[cell_index]; }

  uint32_t get_global_n_x_faces(void) const { return ngx + 1; }
  uint32_t get_global_n_y_faces(void) const { return ngy + 1; }
  uint32_t get_global_n_z_faces(void) const { return ngz + 1; }

  const std::vector<float> &get_silo_x(void) const { return silo_x; }
  const std::vector<float> &get_silo_y(void) const { return silo_y; }
  const std::vector<float> &get_silo_z(void) const { return silo_z; }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Calculate new physical properties and emission energy for each cell on
  // the mesh
  void calculate_photon_energy(IMC_State &imc_state) {
    using Constants::a;
    using Constants::c;
    total_photon_E = 0.0;
    double dt = imc_state.get_dt();
    double op_a, op_s, f, cV, rho;
    double vol;
    double T, Tr, Ts;
    uint32_t step = imc_state.get_step();
    double tot_census_E = 0.0;
    double tot_emission_E = 0.0;
    double tot_source_E = 0.0;
    double pre_mat_E = 0.0;

    uint32_t region_ID;
    Region region;
    for (uint32_t i = 0; i < n_cell; ++i) {
      Cell &e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      Tr = e.get_T_r();
      Ts = e.get_T_s();
      rho = e.get_rho();

      region_ID = e.get_region_ID();
      region = regions[region_ID_to_index[region_ID]];

      op_a = region.get_absorption_opacity(T);
      op_s = region.get_scattering_opacity();
      f = 1.0 / (1.0 + dt * op_a * c * (4.0 * a * pow(T, 3) / (cV * rho)));
      e.set_op_a(op_a);
      e.set_op_s(op_s);
      e.set_f(f);

      m_emission_E[i] =
          replicated_factor * dt * vol * f * op_a * a * c * pow(T, 4);
      if (step > 1)
        m_census_E[i] = 0.0;
      else
        m_census_E[i] = replicated_factor * vol * a * pow(Tr, 4);
      m_source_E[i] = dt * op_a * a * c * pow(Ts, 4);

      pre_mat_E += T * cV * vol * rho;
      tot_emission_E += m_emission_E[i];
      tot_census_E += m_census_E[i];
      tot_source_E += m_source_E[i];
      total_photon_E += m_emission_E[i] + m_census_E[i] + m_source_E[i];
    }

    // set energy for conservation checks
    imc_state.set_pre_mat_E(pre_mat_E);
    imc_state.set_emission_E(tot_emission_E);
    imc_state.set_source_E(tot_source_E);
    if (imc_state.get_step() == 1)
      imc_state.set_pre_census_E(tot_census_E);
  }

  //! Use the absorbed energy and update the material temperature of each
  // cell on the mesh. Set diagnostic and conservation values.
  void update_temperature(std::vector<double> &abs_E,
                          std::vector<double> &track_E, IMC_State &imc_state) {
    using Constants::a;
    using Constants::c;
    // abs E is a global vector
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol, cV, rho, T, T_new;
    uint32_t region_ID;
    Region region;
    for (uint32_t i = 0; i < n_cell; ++i) {
      region_ID = cells[i].get_region_ID();
      region = regions[region_ID_to_index[region_ID]];
      cV = region.get_cV();
      rho = region.get_rho();
      Cell &e = cells[i];
      vol = e.get_volume();
      T = e.get_T_e();
      T_new = T + (abs_E[i] - m_emission_E[i] / replicated_factor) /
                      (cV * vol * rho);
      T_r[i] = std::pow(track_E[i] / (vol * imc_state.get_dt() * a * c), 0.25);
      e.set_T_e(T_new);
      total_abs_E += abs_E[i];
      total_post_mat_E += T_new * cV * vol * rho;
    }
    // zero out absorption tallies for all cells (global)
    abs_E.assign(abs_E.size(), 0.0);
    track_E.assign(track_E.size(), 0.0);
    imc_state.set_absorbed_E(total_abs_E);
    imc_state.set_post_mat_E(total_post_mat_E);
    off_rank_reads = 0;
  }

  //! Set the physical data for the cells on your rank
  void initialize_physical_properties(const Input &input) {
    for (uint32_t i = 0; i < n_cell; ++i) {
      int region_ID = cells[i].get_region_ID();
      // find the region for this cell
      Region region = input.get_region(region_ID);
      // set cell physical properties using region
      cells[i].set_cV(region.get_cV());
      cells[i].set_T_e(region.get_T_e());
      cells[i].set_T_r(region.get_T_r());
      cells[i].set_T_s(region.get_T_s());
      cells[i].set_rho(region.get_rho());
    }
  }

  //! Get census energy vector needed to source particles
  std::vector<double> &get_census_E_ref(void) { return m_census_E; }

  //! Get emission energy vector needed to source particles
  std::vector<double> &get_emission_E_ref(void) { return m_emission_E; }

  //! Get external source energy vector needed to source particles
  std::vector<double> &get_source_E_ref(void) { return m_source_E; }

  //--------------------------------------------------------------------------//
  // member variables
  //--------------------------------------------------------------------------//
private:
  uint32_t ngx;      //!< Number of global x sizes
  uint32_t ngy;      //!< Number of global y sizes
  uint32_t ngz;      //!< Number of global z sizes
  uint32_t n_global; //!< Nuber of global cells

  int32_t rank;   //!< MPI rank of this mesh
  int32_t n_rank; //!< Number of global ranks

  double total_photon_E;   //!< Total photon energy on the mesh
  uint32_t off_rank_reads; //!< Number of off rank reads

  const std::vector<float>
      &silo_x; //!< Global array of x face locations for SILO
  const std::vector<float>
      &silo_y; //!< Global array of y face locations for SILO
  const std::vector<float>
      &silo_z; //!< Global array of z face locations for SILO

  std::vector<Region> regions; //!< Vector of regions in the problem

  //! Factor to reduce emission and initial census in replicated mode
  double replicated_factor;

  uint32_t n_cell; //!< Number of local cells

  std::vector<double> m_census_E;   //!< Census energy vector
  std::vector<double> m_emission_E; //!< Emission energy vector
  std::vector<double> m_source_E;   //!< Source energy vector
  std::vector<double> T_r;          //!< Diagnostic quantity

  Cell *cells; //!< Cell data allocated with MPI_Alloc

  std::vector<uint32_t>
      off_rank_bounds;    //!< Ending value of global ID for each rank
  uint32_t on_rank_start; //!< Start of global index on rank
  uint32_t on_rank_end;   //!< End of global index on rank

  std::unordered_map<uint32_t, uint32_t>
      region_ID_to_index; //!< Maps region ID to index

  Cell current_cell; //!< Off rank cell found in search
};

#endif // mesh_h_
//---------------------------------------------------------------------------//
// end of mesh.h
//---------------------------------------------------------------------------//
