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

#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <pugixml.hpp>
#include <stdlib.h>
#include <string>
#include <unordered_map>

#include "config.h"
#include "constants.h"
#include "region.h"

//==============================================================================
/*!
 * \class Input
 * \brief Reads input data from an XML file and stores that data
 *
 * Pugi's XML parser is used to read an XML input file. This class stores that
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
    using Constants::SOURCE;
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    using Constants::Z_NEG;
    using Constants::Z_POS;
    using std::cout;
    using std::endl;
    using std::vector;

    uint32_t x_key, y_key, z_key, key;
    // initialize nunmber of divisions in each dimension to zero
    uint32_t nx_divisions = 0;
    uint32_t ny_divisions = 0;
    uint32_t nz_divisions = 0;
    uint32_t region_ID;

    vector<double> x;
    vector<double> y;
    vector<double> z;

    pugi::xml_document doc;
    pugi::xml_parse_result load_result = doc.load_file(fileName.c_str());

    // error checking
    if (!load_result) {
      cout << load_result.description() << endl;
      cout << "Improperly formatted xml file" << endl;
      pugi::xml_node proto_node = doc.child("prototype");
      cout << "parsed regions: " << endl;
      for (pugi::xml_node_iterator it = proto_node.begin();
           it != proto_node.end(); ++it)
        cout << it->name() << endl;
      exit(EXIT_FAILURE);
    }

    // check each necessary branch exists and create nodes for each
    pugi::xml_node settings_node = doc.child("prototype").child("common");
    pugi::xml_node debug_node = doc.child("prototype").child("debug_options");
    // One of these spatial input types should be specified, but not both!
    pugi::xml_node spatial_node = doc.child("prototype").child("spatial");
    pugi::xml_node simple_spatial_node =
        doc.child("prototype").child("simple_spatial");
    pugi::xml_node bc_node = doc.child("prototype").child("boundary");
    pugi::xml_node region_node = doc.child("prototype").child("regions");

    if (!settings_node) {
      std::cout << "'common' section not found!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (!spatial_node && !simple_spatial_node) {
      std::cout << "'spatial' or 'simple_spatial' section not found!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    if (spatial_node && simple_spatial_node) {
      std::cout
          << "Cannot specify both 'spatial' and 'simple_spatial' sections!"
          << std::endl;
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

    t_finish = settings_node.child("t_stop").text().as_double();
    dt = settings_node.child("dt_start").text().as_double();
    t_start = settings_node.child("t_start").text().as_double();
    dt_mult = settings_node.child("t_mult").text().as_double();
    dt_max = settings_node.child("dt_max").text().as_double();
    unsigned long long n_photons_long =
        settings_node.child("photons").text().as_llong();
    if (n_photons_long > UINT64_MAX) {
      cout << "ERROR: Can't convert " << n_photons_long;
      cout << " to uint64, too large" << endl;
      exit(EXIT_FAILURE);
    }
    n_photons = static_cast<uint64_t>(n_photons_long);
    seed = settings_node.child("seed").text().as_int();
    output_freq = settings_node.child("output_frequency").text().as_int();

    // use particle combing population control
    std::string tempString = settings_node.child_value("use_combing");
    if (tempString == "FALSE")
      use_comb = 0;
    else if (tempString == "TRUE")
      use_comb = 1;
    else {
      cout << "\"use_combing\" not found or recognized, defaulting to TRUE"
           << endl;
      use_comb = 1;
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
    for (pugi::xml_node_iterator it = spatial_node.begin();
         it != spatial_node.end(); ++it) {
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

    // spatial inputs
    if (simple_spatial_node) {
      double d_x_start, d_x_end, d_y_start, d_y_end, d_z_start, d_z_end;
      uint32_t d_x_cells, d_y_cells, d_z_cells;

      nx_divisions = 1;
      ny_divisions = 1;
      nz_divisions = 1;

      d_x_start = simple_spatial_node.child("x_start").text().as_double();
      d_x_end = simple_spatial_node.child("x_end").text().as_double();
      d_x_cells = simple_spatial_node.child("n_x_cells").text().as_int();

      // push back the master x points for silo
      for (uint32_t i = 0; i < d_x_cells; ++i)
        x.push_back(d_x_start + i * (d_x_end - d_x_start) / d_x_cells);

      d_y_start = simple_spatial_node.child("y_start").text().as_double();
      d_y_end = simple_spatial_node.child("y_end").text().as_double();
      d_y_cells = simple_spatial_node.child("n_y_cells").text().as_int();

      // push back the master y points for silo
      for (uint32_t i = 0; i < d_y_cells; ++i)
        y.push_back(d_y_start + i * (d_y_end - d_y_start) / d_y_cells);

      d_z_start = simple_spatial_node.child("z_start").text().as_double();
      d_z_end = simple_spatial_node.child("z_end").text().as_double();
      d_z_cells = simple_spatial_node.child("n_z_cells").text().as_int();

      // push back the master z points for silo
      for (uint32_t i = 0; i < d_z_cells; ++i)
        z.push_back(d_z_start + i * (d_z_end - d_z_start) / d_z_cells);

      region_ID = simple_spatial_node.child("region_ID").text().as_int();

      x_start.push_back(d_x_start);
      x_end.push_back(d_x_end);
      n_x_cells.push_back(d_x_cells);
      y_start.push_back(d_y_start);
      y_end.push_back(d_y_end);
      n_y_cells.push_back(d_y_cells);
      z_start.push_back(d_z_start);
      z_end.push_back(d_z_end);
      n_z_cells.push_back(d_z_cells);

      // single region -- maps to 0:
      region_map[0] = region_ID;
    }

    // read in boundary conditions
    bool b_error = false;
    bool source_on = false;

    tempString = bc_node.child_value("bc_right");
    if (tempString == "REFLECT")
      bc[X_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[X_POS] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[X_POS] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    tempString = bc_node.child_value("bc_left");
    if (tempString == "REFLECT")
      bc[X_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[X_NEG] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[X_NEG] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    tempString = bc_node.child_value("bc_up");
    if (tempString == "REFLECT")
      bc[Y_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Y_POS] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[Y_POS] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    tempString = bc_node.child_value("bc_down");
    if (tempString == "REFLECT")
      bc[Y_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Y_NEG] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[Y_NEG] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    tempString = bc_node.child_value("bc_top");
    if (tempString == "REFLECT")
      bc[Z_POS] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Z_POS] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[Z_POS] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    tempString = bc_node.child_value("bc_bottom");
    if (tempString == "REFLECT")
      bc[Z_NEG] = REFLECT;
    else if (tempString == "VACUUM")
      bc[Z_NEG] = VACUUM;
    else if (tempString == "SOURCE") {
      bc[Z_NEG] = SOURCE;
      source_on = true;
    }
    else
      b_error = true;

    if(source_on) {
      T_source = bc_node.child("T_source").text().as_double();
    }

    if (b_error) {
      cout << "ERROR: Boundary type not reconginzed. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    // end boundary node

    // read in region data
    for (pugi::xml_node_iterator it = region_node.begin();
         it != region_node.end(); ++it) {
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

    // set total number of divisions
    n_divisions = nx_divisions * ny_divisions * nz_divisions;

    // get global cell counts
    n_global_x_cells = std::accumulate(n_x_cells.begin(), n_x_cells.end(), 0);
    n_global_y_cells = std::accumulate(n_y_cells.begin(), n_y_cells.end(), 0);
    n_global_z_cells = std::accumulate(n_z_cells.begin(), n_z_cells.end(), 0);

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
  }

  //! Destructor
  ~Input(){};

  //! Print the information read from the input file
  void print_problem_info(void) const {
    using Constants::a;
    using Constants::c;
    using Constants::SOURCE;
    using std::cout;
    using std::endl;

    cout << "Problem Specifications:";
    cout << "Constants -- c: " << c << " (cm/sh) , a: " << a << endl;
    cout << "Run Parameters-- Photons: " << n_photons
         << ", time finish: " << t_finish;
    cout << " (sh), time step: " << dt << " (sh) ," << endl;

    cout << " time multiplier: " << dt_mult << " , max dt:" << dt_max;
    cout << " (sh), Random number seed: " << seed;
    cout << " , output frequency: " << output_freq << endl;

    if (use_comb)
      cout << "Combing census enabled (default)" << endl;
    else
      cout << "No combing" << endl;

    if(bc[0] == SOURCE || bc[1] == SOURCE || bc[2] == SOURCE || bc[3] == SOURCE || bc[4] == SOURCE || bc[5] == SOURCE)
      std::cout<<"Source boundary on, T_source: "<<T_source<<std::endl;

    if (print_verbose)
      cout << "Verbose printing mode enabled" << endl;
    else
      cout << "Terse printing mode (default)" << endl;

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
    cout << "REPLICATED" << endl;
    cout << "No parameters are needed in replicated mode";
    cout << endl;

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

  //! Return the value of the use comb option
  bool get_comb_bool(void) const { return use_comb; }
  //! Return the value of the verbose printing option
  bool get_verbose_print_bool(void) const { return print_verbose; }
  //! Return the value of the mesh print option
  bool get_print_mesh_info_bool(void) const { return print_mesh_info; }
  //! Return the frequency of timestep summary printing
  uint32_t get_output_freq(void) const { return output_freq; }

  //! Return the timestep size (shakes)
  double get_dt(void) const { return dt; }
  //! Return the starting time (shakes)
  double get_time_start(void) const { return t_start; }
  //! Return the finish time (shakes)
  double get_time_finish(void) const { return t_finish; }
  //! Return the multiplication factor for the timestep
  double get_time_mult(void) const { return dt_mult; }
  //! Return the maximum timestep size (shakes)
  double get_dt_max(void) const { return dt_max; }
  //! Return the input seed for the RNG
  int get_rng_seed(void) const { return seed; }
  //! Return the number of photons set in the input file to run
  uint64_t get_number_photons(void) const { return n_photons; }

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

  Constants::bc_type bc[6]; //!< Boundary condition array

  // timing
  double t_start;  //!< Starting time
  double dt;      //!< Timestep size
  double t_finish; //!< Finish time
  double dt_mult;   //!< Timestep multiplier
  double dt_max;   //!< Maximum timestep size

  // material
  std::vector<Region> regions; //!< Vector of regions in the problem

  //! Maps unique key to user set ID for a region
  std::map<uint32_t, uint32_t> region_map;

  //! Maps user set region ID to the index in the regions vector
  std::map<uint32_t, uint32_t> region_ID_to_index;

  // source
  double T_source; //!< Temperature of source

  // Monte Carlo parameters
  uint64_t n_photons; //!< Photons to source each timestep
  uint32_t seed;      //!< Random number seed

  // Method parameters
  bool use_comb;  //!< Comb census photons

  // Debug parameters
  uint32_t output_freq; //!< How often to print temperature information
  bool print_verbose;   //!< Verbose printing flag
  bool print_mesh_info; //!< Mesh information printing flag

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
};

#endif // input_h_
//---------------------------------------------------------------------------//
// end of input.h
//---------------------------------------------------------------------------//
