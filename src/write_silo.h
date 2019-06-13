//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   write_silo.h
 * \author Alex Long
 * \date   April 11 2016
 * \brief  Writes SILO output file for data visualization
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef write_silo_h_
#define write_silo_h_

#include <array>
#include <sstream>
#include <string>
#include <vector>

#ifdef VIZ_LIBRARIES_FOUND
#include <iomanip>
#include <silo.h>
#endif

#include "config.h"
#include "constants.h"
#include "imc_state.h"

//! All ranks perform reductions to produce global arrays and rank zero
// writes the SILO file for visualization
void write_silo(const Mesh &mesh, const double &arg_time, const uint32_t &step,
                const double &r_transport_time, const double &r_mpi_time,
                const int &rank, const int &n_rank, bool rep_flag = true) {

#ifdef VIZ_LIBRARIES_FOUND
  using Constants::ELEMENT;
  using Constants::X_NEG;
  using Constants::X_POS;
  using Constants::Y_NEG;
  using Constants::Y_POS;
  using std::array;
  using std::string;
  using std::stringstream;
  using std::vector;

  // need a non-const double to pass to SILO
  double time = arg_time;

  //generate name for this silo file
  stringstream ss;
  ss.setf(std::ios::showpoint);
  ss << std::setprecision(3);
  ss << "output_" << step << ".silo";
  string file = ss.str();

  int nx = mesh.get_global_n_x_faces();
  int ny = mesh.get_global_n_y_faces();
  int nz = mesh.get_global_n_z_faces();

  // set number of dimensions
  int ndims;
  // Use a 2D mesh for one z cell (2 faces)
  if (nz == 2)
    ndims = 2;
  // Otherwise use 3D mesh for 3 or more z faces
  else
    ndims = 3;

  // generate title of plot
  stringstream tt;
  tt.setf(std::ios::showpoint);
  tt << std::setprecision(3);
  if (ndims == 2)
    tt << "2D rectangular mesh, t = " << time << " (sh)";
  else
    tt << "3D rectangular mesh, t = " << time << " (sh)";
  string title = tt.str();

  // get total cells for MPI all_reduce calls
  uint32_t n_xyz_cells;
  if (ndims == 2)
    n_xyz_cells = (nx - 1) * (ny - 1);
  else
    n_xyz_cells = (nx - 1) * (ny - 1) * (nz - 1);

  // make vectors of data for plotting
  vector<int> region_data(n_xyz_cells, -1);
  vector<double> T_e(n_xyz_cells, 0.0);
  vector<double> T_r(n_xyz_cells, 0.0);
  vector<double> transport_time(n_xyz_cells, 0.0);
  vector<double> mpi_time(n_xyz_cells, 0.0);
  vector<int> grip_ID(n_xyz_cells, 0);

  // get rank data, map values from from global ID to SILO ID
  uint32_t n_local = mesh.get_n_local_cells();
  Cell cell;
  uint32_t g_index, silo_index;
  for (uint32_t i = 0; i < n_local; i++) {
    cell = mesh.get_cell(i);
    g_index = cell.get_ID();
    silo_index = cell.get_silo_index();
    region_data[silo_index] = mesh.get_cell(i).get_region_ID();
    // set silo plot variables
    T_e[silo_index] = cell.get_T_e();
    T_r[silo_index] = mesh.get_T_r(i);
    transport_time[silo_index] = r_transport_time;
    mpi_time[silo_index] = r_mpi_time;
    grip_ID[silo_index] = cell.get_grip_ID();
  }

  // don't reduce these quantities in replicated mode
  if (!rep_flag) {
    // reduce to get rank of each cell across all ranks
    MPI_Allreduce(MPI_IN_PLACE, &region_data[0], n_xyz_cells, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    // reduce to get T_e across all ranks
    MPI_Allreduce(MPI_IN_PLACE, &T_e[0], n_xyz_cells, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    // reduce to get T_r across all ranks
    MPI_Allreduce(MPI_IN_PLACE, &T_r[0], n_xyz_cells, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    // reduce to get transport runtime from all ranks
    MPI_Allreduce(MPI_IN_PLACE, &transport_time[0], n_xyz_cells, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    // reduce to get mpi time from all ranks
    MPI_Allreduce(MPI_IN_PLACE, &mpi_time[0], n_xyz_cells, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    // reduce to get grip_ID across all ranks
    MPI_Allreduce(MPI_IN_PLACE, &grip_ID[0], n_xyz_cells, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
  }

  // First rank writes the SILO file
  if (rank == 0) {

    // write the global mesh
    auto x_vec = mesh.get_silo_x();
    auto y_vec = mesh.get_silo_y();
    auto z_vec = mesh.get_silo_z();
    float *x = &x_vec[0];
    float *y = &y_vec[0];
    float *z = &z_vec[0];

    int *dims;
    float **coords;
    int *cell_dims;
    uint32_t n_xyz_cells;

    // do 2D write
    if (ndims == 2) {
      dims = new int[2];
      dims[0] = nx;
      dims[1] = ny;
      coords = new float *[2];
      coords[0] = x;
      coords[1] = y;
      cell_dims = new int[2];
      cell_dims[0] = nx - 1;
      cell_dims[1] = ny - 1;
    }
    // do 3D write
    else {
      dims = new int[3];
      dims[0] = nx;
      dims[1] = ny;
      dims[2] = nz;
      coords = new float *[3];
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
      cell_dims = new int[3];
      cell_dims[0] = nx - 1;
      cell_dims[1] = ny - 1;
      cell_dims[2] = nz - 1;
    }

    //create SILO file for this mesh
    DBfile *dbfile = NULL;
    dbfile = DBCreate(file.c_str(), 0, DB_LOCAL, NULL, DB_PDB);

    // make the correct potion list for 2D and 3D meshes
    DBoptlist *optlist;
    if (ndims == 2) {
      optlist = DBMakeOptlist(4);
      DBAddOption(optlist, DBOPT_XLABEL, (void *)"x");
      DBAddOption(optlist, DBOPT_XUNITS, (void *)"cm");
      DBAddOption(optlist, DBOPT_YLABEL, (void *)"y");
      DBAddOption(optlist, DBOPT_YUNITS, (void *)"cm");
    }
    if (ndims == 3) {
      optlist = DBMakeOptlist(6);
      DBAddOption(optlist, DBOPT_XLABEL, (void *)"x");
      DBAddOption(optlist, DBOPT_XUNITS, (void *)"cm");
      DBAddOption(optlist, DBOPT_YLABEL, (void *)"y");
      DBAddOption(optlist, DBOPT_YUNITS, (void *)"cm");
      DBAddOption(optlist, DBOPT_ZLABEL, (void *)"z");
      DBAddOption(optlist, DBOPT_ZUNITS, (void *)"cm");
    }

    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims, DB_FLOAT,
                  DB_COLLINEAR, optlist);

    // write rank IDs xy to 1D array
    std::vector<int> region_ids = region_data;
    std::sort(region_ids.begin(), region_ids.end());
    auto last = std::unique(region_ids.begin(), region_ids.end());
    region_ids.erase(last, region_ids.end());

    DBPutMaterial(dbfile, "region_ID", "quadmesh", region_ids.size(),
                  &region_ids[0], &region_data[0], cell_dims, ndims, 0, 0, 0, 0,
                  0, DB_INT, NULL);

    // write the material temperature scalar field
    DBoptlist *Te_optlist = DBMakeOptlist(2);
    DBAddOption(Te_optlist, DBOPT_UNITS, (void *)"keV");
    DBAddOption(Te_optlist, DBOPT_DTIME, &time);
    DBPutQuadvar1(dbfile, "T_e", "quadmesh", &T_e[0], cell_dims, ndims, NULL, 0,
                  DB_DOUBLE, DB_ZONECENT, Te_optlist);

    // write the radiation temperature scalar field
    DBoptlist *Tr_optlist = DBMakeOptlist(2);
    DBAddOption(Tr_optlist, DBOPT_UNITS, (void *)"keV");
    DBAddOption(Tr_optlist, DBOPT_DTIME, &time);
    DBPutQuadvar1(dbfile, "T_r", "quadmesh", &T_r[0], cell_dims, ndims, NULL, 0,
                  DB_DOUBLE, DB_ZONECENT, Tr_optlist);

    // free option lists
    DBFreeOptlist(optlist);
    DBFreeOptlist(Te_optlist);
    DBFreeOptlist(Tr_optlist);

    // free data
    delete[] cell_dims;
    delete[] coords;
    delete[] dims;

    // close file
    DBClose(dbfile);
  } // end rank==0
#endif
}

#endif // write_silo_h_
//---------------------------------------------------------------------------//
// end of write_silo.h
//---------------------------------------------------------------------------//
