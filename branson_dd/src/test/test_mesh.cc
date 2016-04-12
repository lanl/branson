/*
  Author: Alex Long
  Date: 4/7/2016
  Name: test_mesh.cc
*/

#include <iostream>
#include <string>
#include <vector>

#include "../constants.h"
#include "../input.h"
#include "../mesh.h"
#include "testing_functions.h"

int main (int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  
  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  using std::cout;
  using std::endl;
  using std::string;

  int nfail = 0;

  // Test that a simple mesh (one division in each dimension) is constructed
  // correctly from the input file (simple_input.xml) and that each cell 
  // is assigned the correct region
  {
    string filename("simple_input.xml");
    Input *input = new Input(filename);

    Mesh mesh(input, rank, n_rank);

    bool simple_mesh_pass = true;

    uint32_t n_cell = mesh.get_n_local_cells();
    if (n_cell != 10*20*30) simple_mesh_pass =false;

    Cell cell;
    for (uint32_t i = 0; i<n_cell; i++) {
      cell = mesh.get_pre_renumber_cell(i);
      if ( cell.get_region_ID() != 6) simple_mesh_pass =false;
    }

    if (simple_mesh_pass) cout<<"TEST PASSED: simple mesh construction"<<endl;
    else { 
      cout<<"TEST FAILED: simple mesh construction"<<endl; 
      nfail++;
    }
    delete input;
  }
  
  // Test that a multi-region mesh is constructed correctly from the input file
  // (three_region_input_mesh.xml) and that each cell is assigned the correct 
  // region
  {
    bool three_region_mesh_pass = true;
    // first test large particle input file
    string three_reg_filename("three_region_mesh_input.xml");
    Input *three_reg_input = new Input(three_reg_filename);

    Mesh mesh(three_reg_input, rank, n_rank);

    uint32_t n_cell = mesh.get_n_local_cells();
    if (n_cell != 21*10) three_region_mesh_pass =false;

    Cell cell;
    double x_low;
    // check the lower x position of the cell to see if the region matches
    // the divisions set in the input file
    for (uint32_t i = 0; i<n_cell; i++) {
      cell = mesh.get_pre_renumber_cell(i);
      x_low = cell.get_x_low();
      //cells in the first region
      if (x_low < 4.0) {
        if ( cell.get_region_ID() != 230) three_region_mesh_pass =false;
      }
      else if (x_low >= 4.0 && cell.get_x_low() < 8.0) {
        if ( cell.get_region_ID() != 177) three_region_mesh_pass =false;
      }
      else if (x_low >= 8.0) {
        if ( cell.get_region_ID() != 11) three_region_mesh_pass =false;
      }
      else {
        // this should not occur, test fails
        three_region_mesh_pass = false;
      }
    }

    if (three_region_mesh_pass) cout<<"TEST PASSED: three region mesh construction"<<endl;
    else { 
      cout<<"TEST FAILED: three region mesh construction"<<endl; 
      nfail++;
    }
    delete three_reg_input;
  }

  MPI_Finalize();

  return nfail;
}
