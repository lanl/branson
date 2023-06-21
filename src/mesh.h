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
#include <iostream>
#include <iomanip>

#include "cell.h"
#include "constants.h"
#include "decompose_mesh.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "info.h"
#include "input.h"
#include "mpi_types.h"
#include "proto_cell.h"
#include "proto_mesh.h"
#include "timer.h"

//==============================================================================
/*!
 * \class Mesh
 * \brief Manages data access, decomposition and parallel communication for mesh
 *
 * Using an Input class, make the mesh with the correct material properties
 * for each region. The mesh numbering and mapping between global indices and local
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
        ngz(input.get_global_n_z_cells()), n_global(ngz * ngy * ngx),
        rank(mpi_info.get_rank()), n_ranks(mpi_info.get_n_rank()),
        verbose_print(input.get_verbose_print_bool()), replicated(false),
        silo_x(input.get_silo_x_ptr()),
        silo_y(input.get_silo_y_ptr()), silo_z(input.get_silo_z_ptr()),
        total_photon_E(0.0), replicated_factor(1.0),
        regions(input.get_regions()) {
    using Constants::bc_type;
    using Constants::CUBE;
    using Constants::ELEMENT;
    using Constants::METIS;
    using Constants::REPLICATED;
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    using Constants::Z_NEG;
    using Constants::Z_POS;
    using std::vector;

    Proto_Mesh proto_mesh(input, mpi_types, mpi_info);

    // if mode is replicated ignore decomposition options, otherwise use
    // metis or a simple cube
    if (input.get_dd_mode() == REPLICATED) {
      replicate_mesh(proto_mesh, mpi_types, mpi_info, imc_p.get_grip_size());
      // get decomposition information from proto mesh
      off_rank_bounds = proto_mesh.get_off_rank_bounds();
      on_rank_start = off_rank_bounds.front();
      on_rank_end = off_rank_bounds.back() - 1;
      replicated_factor = 1.0 / static_cast<double>(n_ranks);
      replicated = true;
    } else if (input.get_decomposition_mode() == METIS) {
      decompose_mesh(proto_mesh, mpi_types, mpi_info, imc_p.get_grip_size(),
                     METIS);
      // get decomposition information from proto mesh
      off_rank_bounds = proto_mesh.get_off_rank_bounds();
      on_rank_start = off_rank_bounds[rank];
      on_rank_end = off_rank_bounds[rank + 1] - 1;
    } else if (input.get_decomposition_mode() == CUBE) {
      decompose_mesh(proto_mesh, mpi_types, mpi_info, imc_p.get_grip_size(),
                     CUBE);
      // get decomposition information from proto mesh
      off_rank_bounds = proto_mesh.get_off_rank_bounds();
      on_rank_start = off_rank_bounds[rank];
      on_rank_end = off_rank_bounds[rank + 1] - 1;
    } else {
      std::cout << "Method/decomposition not recognized, exiting...";
      exit(EXIT_FAILURE);
    }
    const std::vector<Proto_Cell> &proto_cell_list(proto_mesh.get_cell_list());

    // this rank's cells
    n_cell = proto_cell_list.size();
    // size physics data
    m_census_E.resize(n_cell);
    m_emission_E.resize(n_cell);
    m_source_E.resize(n_cell);
    T_r.resize(n_cell);

    // for replicated mode, set the factor that reduces the emission energy and initial census
    // energy

    // get adjacent bounds from proto mesh
    adjacent_procs = proto_mesh.get_proc_adjacency_list();
    // use the proto cells to contstruct the real cells
    for (auto icell : proto_cell_list)
      cells.push_back(Cell(icell));

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
  uint32_t get_rank(void) const { return rank; }
  uint32_t get_offset(void) const { return on_rank_start; }
  uint32_t get_n_global_cells(void) const { return n_global; }
  std::unordered_map<uint32_t, uint32_t> get_proc_adjacency_list(void) const {
    return adjacent_procs;
  }
  uint32_t get_global_num_cells(void) const { return n_global; }

  double get_total_photon_E(void) const { return total_photon_E; }

  void print(void) const {
    for (uint32_t i = 0; i < n_cell; i++)
      cells[i].print();
  }

  const Cell &get_cell_ref(const uint32_t local_index) const { return cells[local_index]; }

  uint32_t get_off_rank_id(const uint32_t &index) const {
    // find rank of index
    bool found = false;
    uint32_t min_i = 0;
    uint32_t max_i = off_rank_bounds.size() - 1;
    uint32_t s_i; // search index
    while (!found) {
      s_i = (max_i + min_i) / 2;
      if (s_i == max_i || s_i == min_i)
        found = true;
      else if (index >= off_rank_bounds[s_i])
        min_i = s_i;
      else
        max_i = s_i;
    }
    return s_i;
  }

  uint32_t get_rank_cell_offset(const int rank) const {
    return off_rank_bounds[rank];
  }

  const Cell *get_cell_ptr_global(const uint32_t global_index) const {
    return &cells[global_index - on_rank_start];
  }
  int32_t get_rank(const uint32_t &index) const {
    int32_t r_rank;
    if (on_processor(index))
      r_rank = rank;
    else
      r_rank = get_off_rank_id(index);
    return r_rank;
  }
  const Cell *get_const_cells_ptr(void) const { return &cells[0]; }

  uint32_t get_local_index(const uint32_t &global_index) const {
    if (on_processor(global_index))
      return global_index - on_rank_start;
    else {
      std::cout<<"about to seg fault probably"<<std::endl;
      return UINT32_MAX;
    }
  }

   Cell get_on_rank_cell(const uint32_t global_index) const {
     // this can only be called after with valid cell index (on rank or in stored
     // cells vector
    if (on_processor(global_index))
      return cells[global_index - on_rank_start];
    else {
      std::cout<<"about to seg fault probably"<<std::endl;
      return cells[global_index];
    }
   }

  bool on_processor(const uint32_t &index) const {
    return (index >= on_rank_start) && (index <= on_rank_end);
  }

  const std::vector<Cell> &get_cells() const {
    return cells;
  }

  std::vector<double> get_census_E(void) const { return m_census_E; }
  std::vector<double> get_emission_E(void) const { return m_emission_E; }
  std::vector<double> get_source_E(void) const { return m_source_E; }

  //! Get the radiation temperature in a cell (for plotting/diagnostics)
  double get_T_r(const uint32_t cell_index) const { return T_r[cell_index]; }

  uint32_t get_global_n_x_faces(void) const { return ngx + 1; }
  uint32_t get_global_n_y_faces(void) const { return ngy + 1; }
  uint32_t get_global_n_z_faces(void) const { return ngz + 1; }

  float *get_silo_x(void) const { return silo_x; }
  float *get_silo_y(void) const { return silo_y; }
  float *get_silo_z(void) const { return silo_z; }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Calculate new physical properties and emission energy for each cell on
  // the mesh
  void calculate_photon_energy(IMC_State &imc_state, const uint32_t n_user_photons) {
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

      // source temperature will be zero
      m_source_E[i] = replicated_factor * 0.25 * a * c *  e.get_source_area() * pow(Ts, 4) * dt;

      pre_mat_E += T * cV * vol * rho;
      tot_emission_E += m_emission_E[i];
      tot_census_E += m_census_E[i];
      tot_source_E += m_source_E[i];
      total_photon_E += m_source_E[i] + m_census_E[i] + m_emission_E[i];
    }

    // adjust the census, emission and source energies for replicated mode to avoid having multiple
    // ranks make small energy photons, recaculculate total_photon_E on this rank
    if(replicated) {
      double global_source_E{tot_emission_E + tot_census_E + tot_source_E};
      MPI_Allreduce(MPI_IN_PLACE, &global_source_E, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);

      tot_census_E = 0.0;
      tot_emission_E = 0.0;
      tot_source_E = 0.0;
      total_photon_E = 0.0;
      for (uint32_t i = 0; i < n_cell; ++i) {
        if(step ==1 && m_census_E[i] > 0.0 && int(n_user_photons*(m_census_E[i] / global_source_E)) == 0 ) {
            m_census_E[i] = (i % n_ranks == rank)  ? m_census_E[i]/replicated_factor : 0.0;
        }
        if(m_emission_E[i] > 0.0 && int(n_user_photons*(m_emission_E[i] / global_source_E)) == 0 ) {
            m_emission_E[i] = (i % n_ranks == rank) ? m_emission_E[i]/replicated_factor : 0.0;
        }
        if(m_source_E[i] > 0.0 && int(n_user_photons*(m_source_E[i] / global_source_E)) == 0 ) {
            m_source_E[i] = (i % n_ranks == rank) ? m_source_E[i]/replicated_factor : 0.0;
        }
        tot_emission_E += m_emission_E[i];
        tot_census_E += m_census_E[i];
        tot_source_E += m_source_E[i];
        total_photon_E += m_source_E[i] + m_census_E[i] + m_emission_E[i];
      } // for loop over cells
    } // if replicated

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
    using std::setiosflags;
    using std::ios;
    using std::setw;
    // abs E is a global vector
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol, cV, rho, T, T_new;
    uint32_t region_ID;
    Region region;

    // in replicated mode reduce the emission energy as some ranks may have had their emission
    // energy zeroed out for some cells to try to keep photon counts close to n_user_photons
    if(replicated)
      MPI_Allreduce(MPI_IN_PLACE, m_emission_E.data(), m_emission_E.size(),
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    if (verbose_print && rank==0) {
      std::cout.precision(8);
      std::cout<<"-------- VERBOSE PRINT BLOCK: CELL TEMPERATURE --------"<<std::endl;
      std::cout<<setiosflags(ios::right) << setw(12) << "cell"<<" ";
      std::cout<<setiosflags(ios::right) << setw(12) << "T_e"<<" ";
      std::cout<<setiosflags(ios::right) << setw(12) << "T_r"<<" ";
      std::cout<<setiosflags(ios::right) << setw(12) << "abs_E"<<" ";
      std::cout<<std::endl;
    }
    // calculate new temperatures, update global conservation quantities
    for (uint32_t i = 0; i < n_cell; ++i) {
      region_ID = cells[i].get_region_ID();
      region = regions[region_ID_to_index[region_ID]];
      cV = region.get_cV();
      rho = region.get_rho();
      Cell &e = cells[i];
      vol = e.get_volume();
      T = e.get_T_e();
      T_new = T + (abs_E[i] - m_emission_E[i]) / (cV * vol * rho);
      T_r[i] = std::pow(track_E[i] / (vol * imc_state.get_dt() * a * c), 0.25);
      e.set_T_e(T_new);
      total_abs_E += abs_E[i];
      total_post_mat_E += T_new * cV * vol * rho;
    }

    // replicated verbose print of updated fields
    if (replicated && verbose_print && rank==0) {
      for (uint32_t i = 0; i < n_cell; ++i) {
        Cell &e = cells[i];
        std::cout<<setiosflags(ios::right) << setw(12) << i <<" ";
        std::cout<<setiosflags(ios::right) << setw(12) << e.get_T_e()<<" ";
        std::cout<<setiosflags(ios::right) << setw(12) << T_r[i]<<" ";
        std::cout<<setiosflags(ios::right) << setw(12) << abs_E[i]<<" ";
        std::cout<<std::endl;
      }
    }
    // ranks take turns writing out cell data for ordered cell output in domain decomposed mode
    else if(!replicated && verbose_print) {
      for (int write_rank =0; write_rank < n_ranks; ++write_rank) {
        if(rank == write_rank) {
          for (uint32_t i = 0; i < n_cell; ++i) {
            Cell &e = cells[i];
            std::cout<<setiosflags(ios::right) << setw(12) << e.get_global_index()<<" ";
            std::cout<<setiosflags(ios::right) << setw(12) << e.get_T_e()<<" ";
            std::cout<<setiosflags(ios::right) << setw(12) << T_r[i]<<" ";
            std::cout<<setiosflags(ios::right) << setw(12) << abs_E[i]<<" ";
            std::cout<<std::endl;
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<std::flush;
      }
    }
    if (verbose_print && rank == 0)
      std::cout<<"-------------------------------------------------------"<<std::endl;
    // zero out absorption tallies for all cells (global)
    abs_E.assign(abs_E.size(), 0.0);
    track_E.assign(track_E.size(), 0.0);
    imc_state.set_absorbed_E(total_abs_E);
    imc_state.set_post_mat_E(total_post_mat_E);
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
      cells[i].set_rho(region.get_rho());
      if (cells[i].get_source_face() != -1)
        cells[i].set_T_s(input.get_source_T());
    }
  }


  std::array<int,3> get_xyz_index(int index) {
    int z = index/ngz;
    int y = (index - z*ngz)/ngy;
    int x = (index - z*ngz - y*ngy)/ngx;
    return std::array<int,3> { x,y,z};
  }

  //! Get census energy vector needed to source particles
  std::vector<double> &get_census_E_ref(void) { return m_census_E; }

  //! Get emission energy vector needed to source particles
  std::vector<double> &get_emission_E_ref(void) { return m_emission_E; }

  //! Get external source energy vector needed to source particles
  std::vector<double> &get_source_E_ref(void) { return m_source_E; }

  std::vector<Cell>::iterator begin() {return cells.begin();}
  std::vector<Cell>::iterator end() {return cells.end();}
  std::vector<Cell>::const_iterator begin() const {return cells.cbegin();}
  std::vector<Cell>::const_iterator end() const {return cells.cend();}

  //--------------------------------------------------------------------------//
  // member variables
  //--------------------------------------------------------------------------//
private:
  uint32_t ngx;      //!< Number of global x sizes
  uint32_t ngy;      //!< Number of global y sizes
  uint32_t ngz;      //!< Number of global z sizes
  uint32_t n_global; //!< Nuber of global cells

  int32_t rank;   //!< MPI rank of this mesh
  int32_t n_ranks; //!< Number of global ranks

  bool verbose_print;
  bool replicated; //!< Flag for replicated mode

  float *silo_x; //!< Global array of x face locations for SILO
  float *silo_y; //!< Global array of y face locations for SILO
  float *silo_z; //!< Global array of z face locations for SILO

  double total_photon_E;   //!< Total photon energy on the mesh
  double replicated_factor;   //!< Factor to reduce emission and initial census in replicated mode

  std::vector<Region> regions; //!< Vector of regions in the problem

  uint32_t n_cell; //!< Number of local cells

  std::vector<double> m_census_E;   //!< Census energy vector
  std::vector<double> m_emission_E; //!< Emission energy vector
  std::vector<double> m_source_E;   //!< Source energy vector
  std::vector<double> T_r;          //!< Diagnostic quantity

  std::vector<Cell> cells; //!< Cell data allocated with MPI_Alloc

  std::vector<uint32_t> off_rank_bounds;    //!< Ending value of global ID for each rank
  uint32_t on_rank_start; //!< Start of global index on rank
  uint32_t on_rank_end;   //!< End of global index on rank

  std::unordered_map<uint32_t, uint32_t> adjacent_procs; //!< List of adjacent processors
  std::unordered_map<uint32_t, uint32_t> region_ID_to_index; //!< Maps region ID to index

  Cell current_cell; //!< Off rank cell found in search
};

#endif // mesh_h_
//---------------------------------------------------------------------------//
// end of mesh.h
//---------------------------------------------------------------------------//
