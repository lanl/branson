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
#include "mpi.h"
#include "mpi_types.h"
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
  Input(std::string fileName, const MPI_Types &mpi_types) {
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
    // DD methods
    using Constants::CUBE;
    using Constants::METIS;
    using Constants::NO_DECOMP;
    using Constants::PARTICLE_PASS;
    using Constants::REPLICATED;
    using std::cout;
    using std::endl;
    using std::vector;

    // root rank reads file, prints warnings, and broadcasts to others
    int rank;
    int n_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
    if (rank == 0) {
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
      pugi::xml_node spatial_node = doc.child("prototype").child("spatial");
      pugi::xml_node bc_node = doc.child("prototype").child("boundary");
      pugi::xml_node region_node = doc.child("prototype").child("regions");

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
      tMult = settings_node.child("t_mult").text().as_double();
      dtMax = settings_node.child("dt_max").text().as_double();
      unsigned long long n_photons_long =
          settings_node.child("photons").text().as_llong();
      if (n_photons_long > UINT64_MAX) {
        cout << "ERROR: Can't convert " << n_photons_long;
        cout << " to uint64, too large" << endl;
        exit(EXIT_FAILURE);
      }
      n_photons = static_cast<uint64_t>(n_photons_long);
      seed = settings_node.child("seed").text().as_int();
      // if grip size is not found, set it to very high so no overdecomposition
      // is used
      if (!(settings_node.child("grip_size")))
        grip_size = 100000000;
      else
        grip_size = settings_node.child("grip_size").text().as_int();

      map_size = settings_node.child("map_size").text().as_int();
      output_freq = settings_node.child("output_frequency").text().as_int();

      // use gpu transporter if available
      std::string tempString = settings_node.child_value("use_gpu_transporter");
      if (tempString == "FALSE")
        use_gpu_transporter = 0;
      else if (tempString == "TRUE")
        use_gpu_transporter = 1;
      else {
        cout << "\"use_gpu_transporter\" not found or recognized, defaulting to FALSE"
             << endl;
        use_gpu_transporter = 0;
      }

      // use particle combing population control
      tempString = settings_node.child_value("use_combing");
      if (tempString == "FALSE")
        use_comb = 0;
      else if (tempString == "TRUE")
        use_comb = 1;
      else {
        cout << "\"use_combing\" not found or recognized, defaulting to TRUE"
             << endl;
        use_comb = 1;
      }

      // write silo flag
      write_silo = false;
      tempString = settings_node.child_value("write_silo");
      if (tempString == "TRUE")
        write_silo = true;

      // domain decomposed transport aglorithm
      tempString = settings_node.child_value("dd_transport_type");
      if (tempString == "PARTICLE_PASS")
        dd_mode = PARTICLE_PASS;
      else if (tempString == "REPLICATED")
        dd_mode = REPLICATED;
      else {
        cout << "WARNING: Domain decomposition method not recognized... ";
        cout << "setting to PARTICLE PASSING method" << endl;
        dd_mode = PARTICLE_PASS;
      }

      if (n_ranks == 1 && dd_mode == PARTICLE_PASS) {
        cout<<"WARNING: Domain decomposition method set to PARTICLE_PASS but there is only one";
        cout<<" rank, setting to REPLICATED"<<endl;
        dd_mode = REPLICATED;
      }

      // set number of OpenMP threads, default to 1, warn if set and no OpenMP
      if (settings_node.child("n_omp_threads")) {
        n_omp_threads = settings_node.child("n_omp_threads").text().as_int();
#ifndef USE_OPENMP
        std::cout<<"WARNING: Number of OpenMP threads set but Branson was not configured with";
        std::cout<<" OpenMP, OpenMP threads will not be used"<<std::endl;
#endif
      }
      else
        n_omp_threads = 1;

      // domain decomposition method, only do non-repliacted
      tempString = settings_node.child_value("mesh_decomposition");
      if (tempString == "METIS")
        decomp_mode = METIS;
      else if (tempString == "CUBE")
        decomp_mode = CUBE;
      else if (tempString != "" && dd_mode == REPLICATED) {
        std::cout << "Replicated transport mode, mesh decomposition method";
        std::cout << " ignored" << std::endl;
        decomp_mode = NO_DECOMP;
      }
      else if (dd_mode == REPLICATED) {
        decomp_mode = NO_DECOMP;
      }
      else {
        cout << "WARNING: Mesh decomposition method is required but not ";
        cout << "recognized... ";
        cout << "setting to METIS method" << endl;
        decomp_mode = METIS;
      }

      // preferred number of particles per MPI send
      tempString = settings_node.child_value("particle_message_size");
      if (tempString == "" && dd_mode == PARTICLE_PASS) {
        std::cout<<"particle_message_size not found in settings, defaulting to 10000"<<std::endl;
        particle_message_size = 10000;
      }
      else {
        particle_message_size =
            settings_node.child("particle_message_size").text().as_double();
      }

      // number of particles to run between MPI message checks
      tempString = settings_node.child_value("batch_size");
      if (tempString == "" && dd_mode == PARTICLE_PASS) {
        std::cout<<"batch_size not found in settings, defaulting to 10000"<<std::endl;
        batch_size = 10000;
      }
      else {
        batch_size =
            settings_node.child("batch_size").text().as_double();
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

      if (b_error) {
        cout << "ERROR: Boundary type not reconginzed. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
      // end boundary node

      // if source was turned on read in source temperature
      if(source_on) {
        T_source = bc_node.child("T_source").text().as_double();
      }

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
    } // end xml parse

    const int n_bools = 5;
    const int n_uint = 17;
    const int n_doubles = 6;
    MPI_Datatype MPI_Region = mpi_types.get_region_type();

    // root rank broadcasts read values
    if (rank == 0) {
      // some helper values
      uint32_t n_regions = regions.size();
      uint32_t n_x_div = n_x_cells.size();
      uint32_t n_y_div = n_y_cells.size();
      uint32_t n_z_div = n_z_cells.size();

      // bools
      vector<int> all_bools = {write_silo, use_comb,
                               print_verbose, print_mesh_info, use_gpu_transporter};
      MPI_Bcast(all_bools.data(), n_bools, MPI_INT, 0, MPI_COMM_WORLD);

      // bcs
      vector<int> bcast_bcs = {bc[0], bc[1], bc[2], bc[3], bc[4], bc[5]};
      MPI_Bcast(&bcast_bcs[0], 6, MPI_INT, 0, MPI_COMM_WORLD);

      // uint32
      vector<uint32_t> all_uint = {seed,
                                   dd_mode,
                                   decomp_mode,
                                   n_omp_threads,
                                   output_freq,
                                   grip_size,
                                   map_size,
                                   batch_size,
                                   particle_message_size,
                                   n_divisions,
                                   n_global_x_cells,
                                   n_global_y_cells,
                                   n_global_z_cells,
                                   n_regions,
                                   n_x_div,
                                   n_y_div,
                                   n_z_div};

      MPI_Bcast(all_uint.data(), n_uint, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      // uint64
      MPI_Bcast(&n_photons, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

      // double
      vector<double> all_doubles = {tStart, dt,    tFinish,
                                    tMult,  dtMax, T_source};
      MPI_Bcast(all_doubles.data(), n_doubles, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // region processing
      MPI_Bcast(regions.data(), n_regions, MPI_Region, 0, MPI_COMM_WORLD);
      vector<uint32_t> division_key;
      vector<uint32_t> region_at_division;
      for (auto rmap : region_map) {
        division_key.push_back(rmap.first);
        region_at_division.push_back(rmap.second);
      }
      if (division_key.size() != n_divisions ||
          region_at_division.size() != n_divisions)
        std::cout << "something went wrong in division key communication"
                  << std::endl;

      MPI_Bcast(&division_key[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&region_at_division[0], n_divisions, MPI_UNSIGNED, 0,
                MPI_COMM_WORLD);

      // mesh spacing and coordinate processing
      MPI_Bcast(&x_start[0], n_x_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&x_end[0], n_x_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_x_cells[0], n_x_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&y_start[0], n_y_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&y_end[0], n_y_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_y_cells[0], n_y_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&z_start[0], n_z_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&z_end[0], n_z_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_z_cells[0], n_z_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    } else {
      // set bools
      vector<int> all_bools(8);

      MPI_Bcast(&all_bools[0], n_bools, MPI_INT, 0, MPI_COMM_WORLD);
      write_silo = all_bools[0];
      use_comb = all_bools[1];
      print_verbose = all_bools[2];
      print_mesh_info = all_bools[3];
      use_gpu_transporter = all_bools[4];

      // set bcs
      vector<int> bcast_bcs(6);
      MPI_Bcast(&bcast_bcs[0], 6, MPI_INT, 0, MPI_COMM_WORLD);
      for (int i = 0; i < 6; ++i)
        bc[i] = Constants::bc_type(bcast_bcs[i]);

      // set uints
      vector<uint32_t> all_uint(n_uint);
      MPI_Bcast(all_uint.data(), n_uint, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      seed = all_uint[0];
      dd_mode = all_uint[1];
      decomp_mode = all_uint[2];
      n_omp_threads = all_uint[3];
      output_freq = all_uint[4];
      grip_size = all_uint[5];
      map_size = all_uint[6];
      batch_size = all_uint[7];
      particle_message_size = all_uint[8];
      n_divisions = all_uint[9];
      n_global_x_cells = all_uint[10];
      n_global_y_cells = all_uint[11];
      n_global_z_cells = all_uint[12];
      const uint32_t n_regions = all_uint[13];
      const uint32_t n_x_div = all_uint[14];
      const uint32_t n_y_div = all_uint[15];
      const uint32_t n_z_div = all_uint[16];

      // uint64
      MPI_Bcast(&n_photons, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

      vector<double> all_doubles(n_doubles);
      MPI_Bcast(&all_doubles[0], n_doubles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      tStart = all_doubles[0];
      dt = all_doubles[1];
      tFinish = all_doubles[2];
      tMult = all_doubles[3];
      dtMax = all_doubles[4];
      T_source = all_doubles[5];

      // region processing (broadcast directly into member variable)
      regions.resize(n_regions);
      MPI_Bcast(&regions[0], n_regions, MPI_Region, 0, MPI_COMM_WORLD);

      vector<uint32_t> division_key(n_divisions);
      vector<uint32_t> region_at_division(n_divisions);
      MPI_Bcast(&division_key[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&region_at_division[0], n_divisions, MPI_UNSIGNED, 0,
                MPI_COMM_WORLD);
      for (uint32_t i = 0; i < n_divisions; ++i) {
        region_map[division_key[i]] = region_at_division[i];
      }
      for (uint32_t i = 0; i < n_regions; ++i)
        region_ID_to_index[regions[i].get_ID()] = i;

      // mesh spacing and coordinate processing
      x_start.resize(n_x_div);
      x_end.resize(n_x_div);
      n_x_cells.resize(n_x_div);

      y_start.resize(n_y_div);
      y_end.resize(n_y_div);
      n_y_cells.resize(n_y_div);

      z_start.resize(n_z_div);
      z_end.resize(n_z_div);
      n_z_cells.resize(n_z_div);

      MPI_Bcast(&x_start[0], n_x_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&x_end[0], n_x_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_x_cells[0], n_x_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&y_start[0], n_y_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&y_end[0], n_y_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_y_cells[0], n_y_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&z_start[0], n_z_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&z_end[0], n_z_div, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&n_z_cells[0], n_z_div, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      silo_x = new float[n_global_x_cells];
      silo_y = new float[n_global_y_cells];
      silo_z = new float[n_global_z_cells];
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  //! Destructor
  ~Input() {
    delete[] silo_x;
    delete[] silo_y;
    delete[] silo_z;
  };

  //! Print the information read from the input file
  void print_problem_info() const {
    using Constants::a;
    using Constants::c;
    using Constants::SOURCE;
    using std::cout;
    using std::endl;
    // DD methods
    using Constants::PARTICLE_PASS;
    using Constants::REPLICATED;
    // mesh decomposition
    using Constants::CUBE;
    using Constants::METIS;

    cout << "Problem Specifications:";
    cout << "Constants -- c: " << c << " (cm/sh) , a: " << a << endl;
    cout << "Run Parameters-- Photons: " << n_photons
         << ", time finish: " << tFinish;
    cout << " (sh), time step: " << dt << " (sh) ," << endl;

    cout << " time multiplier: " << tMult << " , max dt:" << dtMax;
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
    if (dd_mode == PARTICLE_PASS) {
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
#ifdef USE_OPENMP
    std::cout<<"Branson built with OpenMP"<<std::endl;
    std::cout<<"Number of OpenMP threads per rank: "<<n_omp_threads<<std::endl;
#else
    std::cout<<"Branson built WITHOUT OpenMP, number of OpenMP threads set in the input will have";
    std::cout<<" no affect"<<std::endl;
#endif

    cout << "Mesh decomposition: ";
    if (decomp_mode == METIS && dd_mode != REPLICATED)
      cout << "METIS" << endl;
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
  uint32_t get_global_n_x_cells() const { return n_global_x_cells; }
  //! Return the number of global cells in the y direction
  uint32_t get_global_n_y_cells() const { return n_global_y_cells; }
  //! Return the number of global cells in the z direction
  uint32_t get_global_n_z_cells() const { return n_global_z_cells; }

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
  float *get_silo_x_ptr() const { return silo_x; }
  //! Return a pointer to the y coordinates of the mesh in SILO format
  float *get_silo_y_ptr() const { return silo_y; }
  //! Return a pointer to the z coordinates of the mesh in SILO format
  float *get_silo_z_ptr() const { return silo_z; }

  //! Return the total number of x divisions in the problem
  uint32_t get_n_x_divisions() const { return n_x_cells.size(); }
  //! Return the total number of y divisions in the problem
  uint32_t get_n_y_divisions() const { return n_y_cells.size(); }
  //! Return the total number of z divisions in the problem
  uint32_t get_n_z_divisions() const { return n_z_cells.size(); }

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
  bool get_comb_bool() const { return use_comb; }
  //! Return the value of the write SILO option
  bool get_write_silo_bool() const { return write_silo; }
  //! Return the value of the verbose printing option
  bool get_verbose_print_bool() const { return print_verbose; }
  //! Return the value of the mesh print option
  bool get_print_mesh_info_bool() const { return print_mesh_info; }
  //! Return the use_gpu_transporter option
  bool get_use_gpu_transporter_bool() const { return use_gpu_transporter; }
  //! Return the frequency of timestep summary printing
  uint32_t get_output_freq() const { return output_freq; }

  //! Return the timestep size (shakes)
  double get_dt() const { return dt; }
  //! Return the starting time (shakes)
  double get_time_start() const { return tStart; }
  //! Return the finish time (shakes)
  double get_time_finish() const { return tFinish; }
  //! Return the multiplication factor for the timestep
  double get_time_mult() const { return tMult; }
  //! Return the maximum timestep size (shakes)
  double get_dt_max() const { return dtMax; }
  //! Return the input seed for the RNG
  int get_rng_seed() const { return seed; }
  //! Return the number of photons set in the input file to run
  uint64_t get_number_photons() const { return n_photons; }
  //! Return the batch size (particles to run between parallel processing)
  uint32_t get_batch_size() const { return batch_size; }
  //! Return the user requested number of particles in a message
  uint32_t get_particle_message_size() const {
    return particle_message_size;
  }
  //! Return the user requested grip size
  uint32_t get_grip_size() const { return grip_size; }
  //! Return the size of the working mesh map
  uint32_t get_map_size() const { return map_size; }
  //! Return the domain decomposition algorithm
  uint32_t get_dd_mode() const { return dd_mode; }
  //! Return the domain decomposition algorithm
  uint32_t get_decomposition_mode() const { return decomp_mode; }
  //! Return the number of OpenMP threads
  int32_t get_n_omp_threads() const { return static_cast<int>(n_omp_threads);}

  // source functions
  //! Return the temperature of the face source
  double get_source_T() const { return T_source; }

  //! Return the number of material regions
  uint32_t get_n_regions() const { return regions.size(); }

  //! Return vector of regions
  const std::vector<Region> &get_regions() const { return regions; }

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
  std::map<uint32_t, uint32_t> region_map;

  //! Maps user set region ID to the index in the regions vector
  std::map<uint32_t, uint32_t> region_ID_to_index;

  // source
  double T_source; //!< Temperature of source

  // Monte Carlo parameters
  uint64_t n_photons; //!< Photons to source each timestep
  uint32_t seed;      //!< Random number seed

  // Parallel parameters
  uint32_t dd_mode;     //!< Mode of domain decomposed transport algorithm
  uint32_t decomp_mode; //!< Mode of decomposing mesh
  uint32_t n_omp_threads; //!< Number of OpenMP threads, 1 if no OpenMP

  // Debug parameters
  uint32_t output_freq; //!< How often to print temperature information

  // Bools
  bool use_comb;        //!< Comb census photons
  bool write_silo; //!< Dump SILO output files
  bool print_verbose;   //!< Verbose printing flag
  bool print_mesh_info; //!< Mesh information printing flag
  bool use_gpu_transporter; //!< Run on GPU if availabile

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
