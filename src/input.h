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
#include <map>

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
 * Boost's XML parser is used to read an XML input file. This class stores that
 * information and provides functions to access it. This class also prints the
 * problem information.
 */
//==============================================================================
class Input {
public:
  //! Constructor
  Input(std::string fileName, const MPI_Types &mpi_types)
      : using_simple_mesh(false), using_detailed_mesh(false) {
    using boost::property_tree::ptree;
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

      ptree pt;
      read_xml(fileName, pt);
      std::string tempString;
      // traverse pt
      BOOST_FOREACH (ptree::value_type const &v, pt.get_child("prototype")) {
        // read in basic problem parameters
        if (v.first == "common") {
          tFinish = v.second.get<double>("t_stop");
          dt = v.second.get<double>("dt_start");
          tStart = v.second.get<double>("t_start");
          tMult = v.second.get<double>("t_mult", 1.0);
          dtMax = v.second.get<double>("dt_max", dt);
          n_photons = v.second.get<uint64_t>("photons");
          seed = v.second.get<int>("seed");
          grip_size = v.second.get<int>("grip_size", 10);
          map_size = v.second.get<int>("map_size");
          output_freq = v.second.get<uint32_t>("output_frequency", 1);
          tempString = v.second.get<std::string>("tilt", std::string("FALSE"));
          if (tempString == "TRUE")
            use_tilt = 1;
          else
            use_tilt = 0;
          tempString =
              v.second.get<std::string>("use_combing", std::string("TRUE"));
          if (tempString == "TRUE")
            use_comb = 1;
          else
            use_comb = 0;

          // stratified sampling
          tempString = v.second.get<std::string>("stratified_sampling",
                                                 std::string("FALSE"));
          if (tempString == "TRUE")
            use_strat = true;
          else
            use_strat = false;

          // write silo flag
          tempString =
              v.second.get<std::string>("write_silo", std::string("FALSE"));
          if (tempString == "TRUE")
            write_silo = true;
          else
            write_silo = false;

          // number of particles to run between MPI message checks
          batch_size = v.second.get<uint32_t>("batch_size", 100);

          // preferred number of particles per MPI send
          particle_message_size =
              v.second.get<uint32_t>("particle_message_size", 100);

          // domain decomposed transport aglorithm
          tempString = v.second.get<std::string>("dd_transport_type",
                                                 std::string("CELL_PASS"));
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
          tempString = v.second.get<std::string>("mesh_decomposition",
                                                 std::string("PARMETIS"));
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
        } // end common

        // read in basic problem parameters
        else if (v.first == "debug_options") {
          tempString = v.second.get<std::string>("print_verbose", "FALSE");
          if (tempString == "TRUE")
            print_verbose = true;
          else
            print_verbose = false;
          tempString = v.second.get<std::string>("print_mesh_info", "FALSE");
          if (tempString == "TRUE")
            print_mesh_info = true;
          else
            print_mesh_info = false;
        } // end common

        // read in detailed spatial information
        else if (v.first == "spatial") {
          using_detailed_mesh = true;
          double d_x_start, d_x_end, d_y_start, d_y_end, d_z_start, d_z_end;
          uint32_t d_x_cells, d_y_cells, d_z_cells;
          BOOST_FOREACH (ptree::value_type const &g, v.second) {
            if (g.first == "x_division") {
              // x information for this region
              d_x_start = g.second.get<double>("x_start");
              d_x_end = g.second.get<double>("x_end");
              d_x_cells = g.second.get<uint32_t>("n_x_cells");
              x_start.push_back(d_x_start);
              x_end.push_back(d_x_end);
              n_x_cells.push_back(d_x_cells);
              nx_divisions++;
              // push back the master x points for silo
              for (uint32_t i = 0; i < d_x_cells; ++i)
                x.push_back(d_x_start + i * (d_x_end - d_x_start) / d_x_cells);
            }

            if (g.first == "y_division") {
              // y information for this region
              d_y_start = g.second.get<double>("y_start");
              d_y_end = g.second.get<double>("y_end");
              d_y_cells = g.second.get<uint32_t>("n_y_cells");
              y_start.push_back(d_y_start);
              y_end.push_back(d_y_end);
              n_y_cells.push_back(d_y_cells);
              ny_divisions++;
              // push back the master y points for silo
              for (uint32_t i = 0; i < d_y_cells; ++i)
                y.push_back(d_y_start + i * (d_y_end - d_y_start) / d_y_cells);
            }

            if (g.first == "z_division") {
              // z information for this region
              d_z_start = g.second.get<double>("z_start");
              d_z_end = g.second.get<double>("z_end");
              d_z_cells = g.second.get<uint32_t>("n_z_cells");
              z_start.push_back(d_z_start);
              z_end.push_back(d_z_end);
              n_z_cells.push_back(d_z_cells);
              nz_divisions++;
              // push back the master z points for silo
              for (uint32_t i = 0; i < d_z_cells; ++i)
                z.push_back(d_z_start + i * (d_z_end - d_z_start) / d_z_cells);
            }

            if (g.first == "region_map") {
              x_key = (g.second.get<uint32_t>("x_div_ID"));
              y_key = (g.second.get<uint32_t>("y_div_ID"));
              z_key = (g.second.get<uint32_t>("z_div_ID"));
              region_ID = (g.second.get<uint32_t>("region_ID"));
              // make a unique key using the division ID of x,y and z
              // this mapping allows for 1000 unique divisions in
              // each dimension (way too many)
              key = z_key * 1000000 + y_key * 1000 + x_key;
              region_map[key] = region_ID;
            }
          }
        }

        // read in simple spatial data
        else if (v.first == "simple_spatial") {
          using_simple_mesh = true;

          x_start.push_back(v.second.get<double>("x_start"));
          x_end.push_back(v.second.get<double>("x_end"));
          n_x_cells.push_back(v.second.get<uint32_t>("n_x_cells"));

          y_start.push_back(v.second.get<double>("y_start"));
          y_end.push_back(v.second.get<double>("y_end"));
          n_y_cells.push_back(v.second.get<uint32_t>("n_y_cells"));

          z_start.push_back(v.second.get<double>("z_start"));
          z_end.push_back(v.second.get<double>("z_end"));
          n_z_cells.push_back(v.second.get<uint32_t>("n_z_cells"));

          region_ID = (v.second.get<uint32_t>("region_ID", 0));

          // add spatial information to SILO information
          for (uint32_t i = 0; i < n_x_cells[0]; ++i)
            x.push_back(x_start[0] + i * (x_end[0] - x_start[0]) / n_x_cells[0]);
          for (uint32_t i = 0; i < n_y_cells[0]; ++i)
            y.push_back(y_start[0] + i * (y_end[0] - y_start[0]) / n_y_cells[0]);
          for (uint32_t i = 0; i < n_z_cells[0]; ++i)
            z.push_back(z_start[0] + i * (z_end[0] - z_start[0]) / n_z_cells[0]);

          // map zero key to region_ID
          region_map[0] = region_ID;
        }

        else if (v.first == "boundary") {
          // read in boundary conditions
          bool b_error = false;
          tempString = v.second.get<std::string>("bc_right");
          if (tempString == "REFLECT")
            bc[X_POS] = REFLECT;
          else if (tempString == "VACUUM")
            bc[X_POS] = VACUUM;
          else
            b_error = true;

          tempString = v.second.get<std::string>("bc_left");
          if (tempString == "REFLECT")
            bc[X_NEG] = REFLECT;
          else if (tempString == "VACUUM")
            bc[X_NEG] = VACUUM;
          else
            b_error = true;

          tempString = v.second.get<std::string>("bc_up");
          if (tempString == "REFLECT")
            bc[Y_POS] = REFLECT;
          else if (tempString == "VACUUM")
            bc[Y_POS] = VACUUM;
          else
            b_error = true;

          tempString = v.second.get<std::string>("bc_down");
          if (tempString == "REFLECT")
            bc[Y_NEG] = REFLECT;
          else if (tempString == "VACUUM")
            bc[Y_NEG] = VACUUM;
          else
            b_error = true;

          tempString = v.second.get<std::string>("bc_top");
          if (tempString == "REFLECT")
            bc[Z_POS] = REFLECT;
          else if (tempString == "VACUUM")
            bc[Z_POS] = VACUUM;
          else
            b_error = true;

          tempString = v.second.get<std::string>("bc_bottom");
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
        }

        // read in region data
        else if (v.first == "regions") {
          BOOST_FOREACH (ptree::value_type const &g, v.second) {
            if (g.first == "region") {
              Region temp_region;
              temp_region.set_ID(g.second.get<uint32_t>("ID"));
              temp_region.set_cV(g.second.get<double>("CV", 0.0));
              temp_region.set_rho(g.second.get<double>("density", 0.0));
              temp_region.set_opac_A(g.second.get<double>("opacA", 0.0));
              temp_region.set_opac_B(g.second.get<double>("opacB", 0.0));
              temp_region.set_opac_C(g.second.get<double>("opacC", 0.0));
              temp_region.set_opac_S(g.second.get<double>("opacS", 0.0));
              temp_region.set_T_e(g.second.get<double>("initial_T_e", 0.0));
              // default T_r to T_e if not specified
              temp_region.set_T_r(
                  g.second.get<double>("initial_T_r", temp_region.get_T_e()));
              // map user defined ID to index in region vector
              region_ID_to_index[temp_region.get_ID()] = regions.size();
              // add to list of regions
              regions.push_back(temp_region);
            }
          }
        }

        // if both simple_spatial and detailed_spatial are true, exit with an
        // error message
        if (using_detailed_mesh && using_simple_mesh) {
          cout << "ERROR: Spatial information cannot be specified in both";
          cout << " simple and detailed XML regions. Exiting...\n";
          exit(EXIT_FAILURE);
        }
      } // end xml parse

      // get global cell counts
      n_global_x_cells = std::accumulate(n_x_cells.begin(), n_x_cells.end(), 0);
      n_global_y_cells = std::accumulate(n_y_cells.begin(), n_y_cells.end(), 0);
      n_global_z_cells = std::accumulate(n_z_cells.begin(), n_z_cells.end(), 0);

      // set total number of divisions
      if (using_simple_mesh)
        n_divisions = 1;
      else
        n_divisions = nx_divisions * ny_divisions * nz_divisions;

      // make sure at least one region is specified
      if (!regions.size()) {
        cout << "ERROR: No regions were specified. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }

      // only one region can be specified in simple mesh mode
      if (using_simple_mesh && regions.size() != 1) {
        cout << "ERROR: Only one region may be specified in simple mesh mode. ";
        cout << " Exiting..." << endl;
        exit(EXIT_FAILURE);
      }

      // for simple spatial input, region ID must match region ID in region map
      if (using_simple_mesh) {
        if (regions[0].get_ID() != region_map[0]) {
          cout
              << "ERROR: Region ID in simple spatial blocl must match region ID ";
          cout << "in region block. Exiting..." << endl;
          exit(EXIT_FAILURE);
        }
      }

      // for detailed meshes, the total number of divisions  must equal the
      // number of unique region maps
      if (using_detailed_mesh && n_divisions != region_map.size()) {
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

    const int n_bools = 8;
    const int n_uint = 16;
    const int n_doubles = 6;
    MPI_Datatype MPI_Region = mpi_types.get_region_type();

    if (rank ==0) {
      // some helper values
      uint32_t n_regions= regions.size();
      uint32_t n_x_div = n_x_cells.size();
      uint32_t n_y_div = n_y_cells.size();
      uint32_t n_z_div = n_z_cells.size();
  
      // bools
      vector<int> all_bools = {using_simple_mesh, using_detailed_mesh, write_silo, use_tilt, use_comb, use_strat, print_verbose, print_mesh_info};
      MPI_Bcast(&all_bools[0], n_bools, MPI_INT, 0, MPI_COMM_WORLD);

      // bcs
      vector<int> bcast_bcs = {bc[0], bc[1], bc[2], bc[3], bc[4], bc[5]};
      MPI_Bcast(&bcast_bcs[0], 6, MPI_INT, 0, MPI_COMM_WORLD);

      // uint32
      vector<uint32_t> all_uint = {seed, dd_mode, decomp_mode, output_freq, grip_size, map_size, batch_size, particle_message_size, n_divisions, n_global_x_cells, n_global_y_cells, n_global_z_cells, n_regions, n_x_div, n_y_div, n_z_div};
      MPI_Bcast(&all_uint[0], n_uint, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      // uint64
      MPI_Bcast(&n_photons, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD); 

      // double
      vector<double> all_doubles = {tStart, dt, tFinish, tMult, dtMax, T_source};
      MPI_Bcast(&all_doubles[0], n_doubles, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // region processing
      MPI_Bcast(&regions[0], n_regions, MPI_Region, 0, MPI_COMM_WORLD);
      vector<uint32_t> division_key;
      vector<uint32_t> region_at_division;
      for (auto rmap : region_map) {
        division_key.push_back(rmap.first);
        region_at_division.push_back(rmap.second);
      }
      if (division_key.size() != n_divisions || region_at_division.size() != n_divisions)
        std::cout<<"something went wrong in division key communication"<<std::endl;

      MPI_Bcast(&division_key[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&region_at_division[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD); 

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
    }
    else {
      // set bools
      vector<int> all_bools(8);

      MPI_Bcast(&all_bools[0], n_bools, MPI_INT, 0, MPI_COMM_WORLD);
      using_simple_mesh = all_bools[0]; using_detailed_mesh = all_bools[1]; write_silo = all_bools[2]; use_tilt = all_bools[3]; use_comb = all_bools[4]; use_strat = all_bools[5]; print_verbose = all_bools[6]; print_mesh_info = all_bools[7];

      // set bcs
      vector<int> bcast_bcs(6);
      MPI_Bcast(&bcast_bcs[0], 6, MPI_INT, 0, MPI_COMM_WORLD);
      for (int i =0;i<6;++i) bc[i] = Constants::bc_type(bcast_bcs[i]);

      // set uints
      vector<uint32_t> all_uint(n_uint);
      MPI_Bcast(&all_uint[0], n_uint, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      seed = all_uint[0]; dd_mode = all_uint[1]; decomp_mode = all_uint[2]; output_freq = all_uint[3]; grip_size = all_uint[4]; map_size = all_uint[5]; batch_size = all_uint[6]; particle_message_size = all_uint[7]; n_divisions = all_uint[8]; n_global_x_cells = all_uint[9]; n_global_y_cells = all_uint[10]; n_global_z_cells = all_uint[11];
      const uint32_t n_regions = all_uint[12];
      const uint32_t n_x_div = all_uint[13];
      const uint32_t n_y_div = all_uint[14];
      const uint32_t n_z_div = all_uint[15];
      
      // uint64
      MPI_Bcast(&n_photons, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD); 

      vector<double> all_doubles(n_doubles);
      MPI_Bcast(&all_doubles[0], n_doubles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      tStart = all_doubles[0]; dt = all_doubles[1]; tFinish = all_doubles[2]; tMult = all_doubles[3];  dtMax =  all_doubles[4]; T_source = all_doubles[5];

      // region processing (broadcast directly into member variable)
      regions.resize(n_regions);
      MPI_Bcast(&regions[0], n_regions, MPI_Region, 0, MPI_COMM_WORLD);

      n_divisions = all_uint[8];
      vector<uint32_t> division_key(n_divisions);
      vector<uint32_t> region_at_division(n_divisions);
      MPI_Bcast(&division_key[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast(&region_at_division[0], n_divisions, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      for (uint32_t i=0;i<n_divisions;++i) {
        region_map[division_key[i]] = region_at_division[i];
      }
      for (uint32_t i = 0; i<n_regions;++i)
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
    if (using_simple_mesh) {
      cout << "dx: " << (x_end[0] - x_start[0]) / n_x_cells[0]
           << " dy: " << (y_end[0] - y_start[0]) / n_y_cells[0]
           << " dz: " << (z_end[0] - z_start[0]) / n_z_cells[0] << endl;
    }

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
  uint32_t get_output_freq(void) const { return output_freq; }

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

  void communicate_read_data(void) {

  }


private:
  // flags
  bool using_simple_mesh;   //!< Use the simple mesh specification
  bool using_detailed_mesh; //!< Use the detailed mesh specification
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
  std::map<uint32_t, uint32_t> region_map;

  //! Maps user set region ID to the index in the regions vector
  std::map<uint32_t, uint32_t> region_ID_to_index;

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
  uint32_t output_freq;      //!< How often to print temperature information
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
