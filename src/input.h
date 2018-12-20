//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   input.h
 * \author Alex Long
 * \date   July 18 2014
 * \brief  Reads data from XML input file and makes mesh skeleton
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef input_h_
#define input_h_

// Suppresses warnings found in Boost headers...
// http://wiki.services.openoffice.org/wiki/Writing_warning-free_code#When_all_else_fails
#if defined __GNUC__
#pragma GCC system_header
// Intel defines __GNUC__ by default
#ifdef __INTEL_COMPILER
#pragma warning push
#endif
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#if defined __GNUC__
#pragma GCC system_header
#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#elif defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <unordered_map>

#include <pugixml.hpp>
#include "config.h"
#include "constants.h"
#include "mpi.h"
#include "region.h"

//==============================================================================
/*!
 * \class Input
 * \brief Reads input data from an XML file and stores that data
 *
 * Boost's XML parser is used to read an XML input file. This class stores that
 * information and provides functions to access it. This class also prints the
 * problem information.
 */
//==============================================================================
class Input {
public:
  //! Constructor
  Input(std::string fileName) {
    using Constants::ELEMENT;
    using Constants::REFLECT;
    using Constants::VACUUM;
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    using Constants::Z_NEG;
    using Constants::Z_POS;
    // DD methods
    using Constants::CELL_PASS;
    using Constants::CELL_PASS_RMA;
    using Constants::CUBE;
    using Constants::PARMETIS;
    using Constants::PARTICLE_PASS;
    using Constants::REPLICATED;
    using std::cout;
    using std::endl;
    using std::vector;

    // this is just use to print wanrings for the head node
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    uint32_t x_key, y_key, z_key, key;
    // initialize nunmber of divisions in each dimension to zero
    uint32_t nx_divisions = 0;
    uint32_t ny_divisions = 0;
    uint32_t nz_divisions = 0;
    uint32_t region_ID;

    vector<float> x;
    vector<float> y;
    vector<float> z;

    pugi::xml_document doc;
    pugi::xml_parse_result load_result = doc.load_file(fileName.c_str());

    // error checking
    if (!load_result) {
      cout << load_result.description() << endl;
      cout<<"Improperly formatted xml file" << endl;
      pugi::xml_node proto_node =
        doc.child("prototype");
      cout<<"parsed regions: "<<endl;
      for (pugi::xml_node_iterator it = proto_node.begin(); it != proto_node.end(); ++it)
        cout<<it->name()<<endl;
      exit(EXIT_FAILURE);
    }

    // check each necessary branch exists and create nodes for each
    pugi::xml_node settings_node =
      doc.child("prototype").child("common");
    pugi::xml_node debug_node =
      doc.child("prototype").child("debug");
    pugi::xml_node spatial_node =
      doc.child("prototype").child("spatial");
    pugi::xml_node bc_node =
      doc.child("prototype").child("boundary");
    pugi::xml_node region_node =
      doc.child("prototype").child("regions");

    if (!settings_node) {
      std::cout << "'common' section not found!" << std::endl;           
      exit(EXIT_FAILURE);                                                          
    }
    if (!spatial_node) {
      std::cout << "'spatial' section not found!" << std::endl;           
      exit(EXIT_FAILURE);                                                          
    }
    if (!bc_node) {
      std::cout << "'boundary' section not found!" << std::endl;           
      exit(EXIT_FAILURE);                                                          
    }
    if (!region_node) {
      std::cout << "'regions' section not found!" << std::endl;           
      exit(EXIT_FAILURE);                                                          
    }

    tFinish = settings_node.child("t_stop").text().as_double();
    dt = settings_node.child("dt_start").text().as_double();
    tStart = settings_node.child("t_start").text().as_double();
    tMult =settings_node.child("t_mult").text().as_double(); 
    dtMax = settings_node.child("dt_max").text().as_double();
    unsigned long long n_photons_long = settings_node.child("photons").text().as_llong();
    if (n_photons_long > UINT64_MAX) {
      cout<<"ERROR: Can't convert "<<n_photons_long;
      cout<<" to uint64, too large"<<endl;
      exit(EXIT_FAILURE);
    }
    n_photons = static_cast<uint64_t>(n_photons_long); 
    seed = settings_node.child("seed").text().as_int();
    grip_size = settings_node.child("grip_size").text().as_int();
    map_size = settings_node.child("map_size").text().as_int();
    output_freq = settings_node.child("output_frequency").text().as_int();

    // use particle combing population control
    std::string tempString = settings_node.child_value("use_combing");
    if (tempString == "FALSE")
      use_comb = 0;
    else if (tempString == "TRUE")
      use_comb = 1;
    else {
      cout<<"\"use_combing\" not found or recognized, defaulting to TRUE"<<endl;
      use_comb =1;
    }

    // write silo flag
    tempString = settings_node.child_value("write_silo");
    if (tempString == "TRUE")
      write_silo = true;
    else
      write_silo = false;

    // number of particles to run between MPI message checks
    batch_size = settings_node.child("batch_size").text().as_int(); 

    // preferred number of particles per MPI send
    particle_message_size =settings_node.child("particle_message_size").text().as_double(); 

    // domain decomposed transport aglorithm
    tempString = settings_node.child_value("dd_transport_type"); 
    if (tempString == "CELL_PASS")
      dd_mode = CELL_PASS;
    else if (tempString == "CELL_PASS_RMA")
      dd_mode = CELL_PASS_RMA;
    else if (tempString == "PARTICLE_PASS")
      dd_mode = PARTICLE_PASS;
    else if (tempString == "REPLICATED")
      dd_mode = REPLICATED;
    else {
      if (rank == 0) {
        cout << "WARNING: Domain decomposition method not recognized... ";
        cout << "setting to PARTICLE PASSING method" << endl;
      }
      dd_mode = PARTICLE_PASS;
    }

    // domain decomposition method
    tempString =  settings_node.child("mesh_decomposition").text().as_double();
    if (tempString == "PARMETIS")
      decomp_mode = PARMETIS;
    else if (tempString == "CUBE")
      decomp_mode = CUBE;
    else {
      if (rank == 0) {
        cout << "WARNING: Mesh decomposition method not recognized... ";
        cout << "setting to PARMETIS method" << endl;
      }
      decomp_mode = PARMETIS;
    }
    if (dd_mode == REPLICATED) {
      if (rank == 0) {
        std::cout << "Replicated transport mode, mesh decomposition method";
        std::cout << " ignored" << std::endl;
      }
    }

    // debug options
    print_verbose = false;
    print_mesh_info = false;
    tempString = debug_node.child_value("print_verbose");
    if (tempString == "TRUE")
      print_verbose = true;
    tempString = debug_node.child_value("print_mesh_info");
    if (tempString == "TRUE")
      print_mesh_info = true;

    // spatial inputs
    for (pugi::xml_node_iterator it = spatial_node.begin(); it != spatial_node.end(); ++it)
    {
      std::string name_string = it->name();
      double d_x_start, d_x_end, d_y_start, d_y_end, d_z_start, d_z_end;
      uint32_t d_x_cells, d_y_cells, d_z_cells;
      if (name_string == "x_division") {
        // x information for this region
        d_x_start = it->child("x_start").text().as_double();
        d_x_end = it->child("x_end").text().as_double(); 
        d_x_cells = it->child("n_x_cells").text().as_int();
        x_start.push_back(d_x_start);
        x_end.push_back(d_x_end);
        n_x_cells.push_back(d_x_cells);
        nx_divisions++;
        // push back the master x points for silo
        for (uint32_t i = 0; i < d_x_cells; ++i)
          x.push_back(d_x_start + i * (d_x_end - d_x_start) / d_x_cells);
      }

      if (name_string == "y_division") {
        // y information for this region
        d_y_start = it->child("y_start").text().as_double();
        d_y_end = it->child("y_end").text().as_double();
        d_y_cells = it->child("n_y_cells").text().as_int();
        y_start.push_back(d_y_start);
        y_end.push_back(d_y_end);
        n_y_cells.push_back(d_y_cells);
        ny_divisions++;
        // push back the master y points for silo
        for (uint32_t i = 0; i < d_y_cells; ++i)
          y.push_back(d_y_start + i * (d_y_end - d_y_start) / d_y_cells);
      }

      if (name_string == "z_division") {
        // z information for this region
        d_z_start = it->child("z_start").text().as_double();
        d_z_end = it->child("z_end").text().as_double();
        d_z_cells = it->child("n_z_cells").text().as_int();
        z_start.push_back(d_z_start);
        z_end.push_back(d_z_end);
        n_z_cells.push_back(d_z_cells);
        nz_divisions++;
        // push back the master z points for silo
        for (uint32_t i = 0; i < d_z_cells; ++i)
          z.push_back(d_z_start + i * (d_z_end - d_z_start) / d_z_cells);
      }

      if (name_string == "region_map") {
        x_key = it->child("x_div_ID").text().as_int();
        y_key = it->child("y_div_ID").text().as_int();
        z_key = it->child("z_div_ID").text().as_int();
        region_ID = it->child("region_ID").text().as_int();
        // make a unique key using the division ID of x,y and z
        // this mapping allows for 1000 unique divisions in
        // each dimension (way too many)
        key = z_key * 1000000 + y_key * 1000 + x_key;
        region_map[key] = region_ID;
      }
    }

    // read in boundary conditions
    bool b_error = false;
    tempString = bc_node.child_value("bc_right");
    if (tempString == "REFLECT")
      bc[X_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[X_POS] = VACUUM;
    else
      b_error = true;

    tempString =bc_node.child_value("bc_left"); 
    if (tempString == "REFLECT")
      bc[X_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[X_NEG] = VACUUM;
    else
      b_error = true;

    tempString =bc_node.child_value("bc_up"); 
    if (tempString == "REFLECT")
      bc[Y_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Y_POS] = VACUUM;
    else
      b_error = true;

    tempString = bc_node.child_value("bc_down"); 
    if (tempString == "REFLECT")
      bc[Y_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Y_NEG] = VACUUM;
    else
      b_error = true;

    tempString = bc_node.child_value("bc_top"); 
    if (tempString == "REFLECT")
      bc[Z_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Z_POS] = VACUUM;
    else
      b_error = true;

    tempString =bc_node.child_value("bc_bottom");  
    if (tempString == "REFLECT")
      bc[Z_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Z_NEG] = VACUUM;
    else
      b_error = true;

    if (b_error) {
      cout << "ERROR: Boundary type not reconginzed. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    // end boundary node

    // read in region data
    for (pugi::xml_node_iterator it = region_node.begin(); it != region_node.end(); ++it)
    {
      std::string name_string = it->name();
      if (name_string == "region") {
        Region temp_region;
        temp_region.set_ID(it->child("ID").text().as_int());
        temp_region.set_cV(it->child("CV").text().as_double());
        temp_region.set_rho(it->child("density").text().as_double());
        temp_region.set_opac_A(it->child("opacA").text().as_double());
        temp_region.set_opac_B(it->child("opacB").text().as_double());
        temp_region.set_opac_C(it->child("opacC").text().as_double());
        temp_region.set_opac_S(it->child("opacS").text().as_double());
        temp_region.set_T_e(it->child("initial_T_e").text().as_double());
        // default T_r to T_e if not specified
        temp_region.set_T_r(it->child("initial_T_r").text().as_double());
        // map user defined ID to index in region vector
        region_ID_to_index[temp_region.get_ID()] = regions.size();
        // add to list of regions
        regions.push_back(temp_region);
      }
    }

    // get global cell counts
    n_global_x_cells = std::accumulate(n_x_cells.begin(), n_x_cells.end(), 0);
    n_global_y_cells = std::accumulate(n_y_cells.begin(), n_y_cells.end(), 0);
    n_global_z_cells = std::accumulate(n_z_cells.begin(), n_z_cells.end(), 0);

    // set total number of divisions
    n_divisions = nx_divisions * ny_divisions * nz_divisions;

    // make sure at least one region is specified
    if (!regions.size()) {
      cout << "ERROR: No regions were specified. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }

    // the total number of divisions  must equal the number of unique region maps
    if (n_divisions != region_map.size()) {
      cout << "ERROR: Number of total divisions must match the number of ";
      cout << "unique region maps. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }

    // append the last point values and allocate SILO array
    x.push_back(x_end.back());
    y.push_back(y_end.back());
    z.push_back(z_end.back());
    silo_x = new float[x.size()];
    silo_y = new float[y.size()];
    silo_z = new float[z.size()];
    for (uint32_t i = 0; i < x.size(); ++i)
      silo_x[i] = x[i];
    for (uint32_t j = 0; j < y.size(); ++j)
      silo_y[j] = y[j];
    for (uint32_t k = 0; k < z.size(); ++k)
      silo_z[k] = z[k];

    // batch size should be very large in replicated mode since there is no
    // need to check buffers
    if (dd_mode == REPLICATED)
      batch_size = 100000000;
  }

  //! Destructor
  ~Input() {
    delete[] silo_x;
    delete[] silo_y;
    delete[] silo_z;
  };

  //! Print the information read from the input file
  void print_problem_info(void) const {
    using Constants::a;
    using Constants::c;
    using std::cout;
    using std::endl;
    // DD methods
    using Constants::CELL_PASS;
    using Constants::CELL_PASS_RMA;
    using Constants::PARTICLE_PASS;
    using Constants::REPLICATED;
    // mesh decomposition
    using Constants::CUBE;
    using Constants::PARMETIS;

    cout << "Problem Specifications:";
    cout << "Constants -- c: " << c << " (cm/sh) , a: " << a << endl;
    cout << "Run Parameters-- Photons: " << n_photons
         << ", time finish: " << tFinish;
    cout << " (sh), time step: " << dt << " (sh) ," << endl;

    cout << " time multiplier: " << tMult << " , max dt:" << dtMax;
    cout << " (sh), Random number seed: " << seed;
    cout << " , output frequency: " << output_freq << endl;

    // cout<<"Sampling -- Emission Position: ";
    // if (use_tilt) cout<<"source tilting (x only), ";
    // else cout<<"uniform (default), ";
    // cout<<" Angle: ";
    cout << "Stratified sampling in angle: ";
    if (use_strat)
      cout << "TRUE" << endl;
    else
      cout << "FALSE" << endl;

    if (use_comb)
      cout << "Combing census enabled (default)" << endl;
    else
      cout << "No combing" << endl;

    if (print_verbose)
      cout << "Verbose printing mode enabled" << endl;
    else
      cout << "Terse printing mode (default)" << endl;

#ifdef VIZ_LIBRARIES_FOUND
    if (write_silo)
      cout << "SILO output enabled" << endl;
    else
      cout << "SILO output disabled (default)" << endl;
#else
    if (write_silo)
      cout << "NOTE: SILO libraries not linked... no visualization" << endl;
#endif
    cout << "Spatial Information -- cells x,y,z: " << n_global_x_cells << " ";
    cout << n_global_y_cells << " " << n_global_z_cells << endl;

    cout << "--Material Information--" << endl;
    for (uint32_t r = 0; r < regions.size(); ++r) {
      cout << " heat capacity: " << regions[r].get_cV();
      cout << " opacity constants: " << regions[r].get_opac_A() << " + "
           << regions[r].get_opac_B() << "^" << regions[r].get_opac_C();
      cout << ", scattering opacity: " << regions[r].get_scattering_opacity()
           << endl;
    }

    cout << "--Parallel Information--" << endl;
    cout << "DD algorithm: ";
    if (dd_mode == CELL_PASS) {
      cout << "CELL PASSING" << endl;
      cout << "grip size: " << grip_size;
      cout << ", map size: " << map_size;
      cout << ", batch size: " << batch_size;
      cout << endl;
    } else if (dd_mode == CELL_PASS_RMA) {
      cout << "CELL PASSING (with RMA on MPI windows)" << endl;
      cout << "COMPILE PARAMETER: maximum number of RMA requests" << endl;
      cout << "grip size: " << grip_size;
      cout << ", map size: " << map_size;
      cout << ", batch size: " << batch_size;
      cout << endl;
    } else if (dd_mode == PARTICLE_PASS) {
      cout << "PARTICLE PASSING" << endl;
      cout << "Batch size: " << batch_size;
      cout << ", particle message size: " << particle_message_size;
      cout << endl;
    } else if (dd_mode == REPLICATED) {
      cout << "REPLICATED" << endl;
      cout << "No parameters are needed in replicated mode";
      cout << endl;
    } else {
      cout << "ERROR: Parallel method not specific correctly";
      cout << " Exiting..." << endl;
      exit(EXIT_FAILURE);
    }

    cout << "Mesh decomposition: ";
    if (decomp_mode == PARMETIS && dd_mode != REPLICATED)
      cout << "PARMETIS" << endl;
    else if (decomp_mode == CUBE && dd_mode != REPLICATED)
      cout << "CUBE" << endl;
    else if (dd_mode == REPLICATED)
      cout << "N/A (no decomposition in replicated mode)" << std::endl;
    else {
      cout << "ERROR: Decomposition mode not specified correctly. Exiting...";
      exit(EXIT_FAILURE);
    }

    cout << endl;
  }

  //! Return the number of global cells in the x direction
  uint32_t get_global_n_x_cells(void) const { return n_global_x_cells; }
  //! Return the number of global cells in the y direction
  uint32_t get_global_n_y_cells(void) const { return n_global_y_cells; }
  //! Return the number of global cells in the z direction
  uint32_t get_global_n_z_cells(void) const { return n_global_z_cells; }

  //! Return the number of x cells in a given x division
  uint32_t get_x_division_cells(const uint32_t &div) const {
    return n_x_cells[div];
  }
  //! Return the number of y cells in a given y division
  uint32_t get_y_division_cells(const uint32_t &div) const {
    return n_y_cells[div];
  }
  //! Return the number of z cells in a given z division
  uint32_t get_z_division_cells(const uint32_t &div) const {
    return n_z_cells[div];
  }

  //! Return a pointer to the x coordinates of the mesh in SILO format
  float const *const get_silo_x_ptr(void) const { return silo_x; }
  //! Return a pointer to the y coordinates of the mesh in SILO format
  float const *const get_silo_y_ptr(void) const { return silo_y; }
  //! Return a pointer to the z coordinates of the mesh in SILO format
  float const *const get_silo_z_ptr(void) const { return silo_z; }

  //! Return the total number of x divisions in the problem
  uint32_t get_n_x_divisions(void) const { return n_x_cells.size(); }
  //! Return the total number of y divisions in the problem
  uint32_t get_n_y_divisions(void) const { return n_y_cells.size(); }
  //! Return the total number of z divisions in the problem
  uint32_t get_n_z_divisions(void) const { return n_z_cells.size(); }

  //! Return the x grid spacing in a given x division
  double get_dx(const uint32_t &div) const {
    return (x_end[div] - x_start[div]) / n_x_cells[div];
  }
  //! Return the y grid spacing in a given y division
  double get_dy(const uint32_t &div) const {
    return (y_end[div] - y_start[div]) / n_y_cells[div];
  }
  //! Return the z grid spacing in a given z division
  double get_dz(const uint32_t &div) const {
    return (z_end[div] - z_start[div]) / n_z_cells[div];
  }

  //! Return the starting x position of a given x division
  double get_x_start(const uint32_t &div) const { return x_start[div]; }
  //! Return the starting y position of a given y division
  double get_y_start(const uint32_t &div) const { return y_start[div]; }
  //! Return the starting z position of a given z division
  double get_z_start(const uint32_t &div) const { return z_start[div]; }

  //! Return the value of the use tilt option
  bool get_tilt_bool(void) const { return use_tilt; }
  //! Return the value of the use comb option
  bool get_comb_bool(void) const { return use_comb; }
  //! Return the value of the stratified option
  bool get_stratified_bool(void) const { return use_strat; }
  //! Return the value of the write SILO option
  bool get_write_silo_bool(void) const { return write_silo; }
  //! Return the value of the verbose printing option
  bool get_verbose_print_bool(void) const { return print_verbose; }
  //! Return the value of the mesh print option
  bool get_print_mesh_info_bool(void) const { return print_mesh_info; }
  //! Return the frequency of timestep summary printing
  int get_output_freq(void) const { return output_freq; }

  //! Return the timestep size (shakes)
  double get_dt(void) const { return dt; }
  //! Return the starting time (shakes)
  double get_time_start(void) const { return tStart; }
  //! Return the finish time (shakes)
  double get_time_finish(void) const { return tFinish; }
  //! Return the multiplication factor for the timestep
  double get_time_mult(void) const { return tMult; }
  //! Return the maximum timestep size (shakes)
  double get_dt_max(void) const { return dtMax; }
  //! Return the input seed for the RNG
  int get_rng_seed(void) const { return seed; }
  //! Return the number of photons set in the input file to run
  uint64_t get_number_photons(void) const { return n_photons; }
  //! Return the batch size (particles to run between parallel processing)
  uint32_t get_batch_size(void) const { return batch_size; }
  //! Return the user requested number of particles in a message
  uint32_t get_particle_message_size(void) const {
    return particle_message_size;
  }
  //! Return the user requested grip size
  uint32_t get_grip_size(void) const { return grip_size; }
  //! Return the size of the working mesh map
  uint32_t get_map_size(void) const { return map_size; }
  //! Return the domain decomposition algorithm
  uint32_t get_dd_mode(void) const { return dd_mode; }
  //! Return the domain decomposition algorithm
  uint32_t get_decomposition_mode(void) const { return decomp_mode; }

  // source functions
  //! Return the temperature of the face source
  double get_source_T(void) const { return T_source; }

  //! Return the number of material regions
  uint32_t get_n_regions(void) const { return regions.size(); }

  //! Return vector of regions
  const std::vector<Region> &get_regions(void) const { return regions; }

  //! Return region given user set ID
  const Region &get_region(const uint32_t region_ID) const {
    return regions[region_ID_to_index.at(region_ID)];
  }

  //! Return unique index given division indices

  //! Make a unique key using the division ID of x,y and z. This mapping
  //! allows for 1000 unique divisions in each dimension (way too many)
  uint32_t get_region_index(const uint32_t x_div, const uint32_t y_div,
                            const uint32_t z_div) const {
    const uint32_t key = z_div * 1000000 + y_div * 1000 + x_div;
    return region_ID_to_index.at(region_map.at(key));
  }

  //! Return boundary condition at this direction index
  Constants::bc_type get_bc(const Constants::dir_type &direction) const {
    return bc[direction];
  }

private:
  // flags
  bool write_silo;          //!< Dump SILO output files

  Constants::bc_type bc[6]; //!< Boundary condition array

  // timing
  double tStart;  //!< Starting time
  double dt;      //!< Timestep size
  double tFinish; //!< Finish time
  double tMult;   //!< Timestep multiplier
  double dtMax;   //!< Maximum timestep size

  // material
  std::vector<Region> regions; //!< Vector of regions in the problem

  //! Maps unique key to user set ID for a region
  std::unordered_map<uint32_t, uint32_t> region_map;

  //! Maps user set region ID to the index in the regions vector
  std::unordered_map<uint32_t, uint32_t> region_ID_to_index;

  // source
  double T_source; //!< Temperature of source

  // Monte Carlo parameters
  uint64_t n_photons; //!< Photons to source each timestep
  uint32_t seed;      //!< Random number seed

  // Method parameters
  bool use_tilt;        //!< Use tilting for emission sampling
  bool use_comb;        //!< Comb census photons
  bool use_strat;       //!< Use strafifed sampling
  uint32_t dd_mode;     //!< Mode of domain decomposed transport algorithm
  uint32_t decomp_mode; //!< Mode of decomposing mesh

  // Debug parameters
  int output_freq;      //!< How often to print temperature information
  bool print_verbose;   //!< Verbose printing flag
  bool print_mesh_info; //!< Mesh information printing flag

  // parallel performance parameters
  uint32_t grip_size; //!< Preferred number of cells in a parallel communication
  uint32_t map_size;  //!< Size of stored off-rank mesh cells
  uint32_t batch_size; //!< Particles to run between MPI message checks
  uint32_t
      particle_message_size; //!< Preferred number of particles in MPI sends

  // detailed mesh specifications
  uint32_t n_divisions;            //!< Number of divisions in the mesh
  std::vector<double> x_start;     //!< x starting positions for each division
  std::vector<double> x_end;       //!< x ending positions for each division
  std::vector<double> y_start;     //!< y starting positions for each division
  std::vector<double> y_end;       //!< y ending positions for each division
  std::vector<double> z_start;     //!< z starting positions for each division
  std::vector<double> z_end;       //!< z ending positions for each division
  std::vector<uint32_t> n_x_cells; //!< Number of x cells in each division
  std::vector<uint32_t> n_y_cells; //!< Number of y cells in each division
  std::vector<uint32_t> n_z_cells; //!< Number of z cells in each division
  uint32_t n_global_x_cells; //!< Total number of x cells over all divisions
  uint32_t n_global_y_cells; //!< Total number of y cells over all divisions
  uint32_t n_global_z_cells; //!< Total number of z cells over all divisions

  // arrays for silo
  float *silo_x; //!< x positions for SILO arrays (i + nx*j + nx*ny*k)
  float *silo_y; //!< y positions for SILO arrays (i + nx*j + nx*ny*k)
  float *silo_z; //!< z positions for SILO arrays (i + nx*j + nx*ny*k)
};

#endif // input_h_
//---------------------------------------------------------------------------//
// end of input.h
//---------------------------------------------------------------------------//
