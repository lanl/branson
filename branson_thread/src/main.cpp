/*
  Author: Alex Long
  Date: 7/18/2014
  Name: main.cpp

  This is a version of BRANSON for on-node IMC parallelism. It
  has one OpenMP loop over elements--each element's particles
  are made and transported. The census is sorted and then copied
  by each element during the next timestep.
*/

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "photon.h"
#include "mesh.h"
#include "imc_state.h"
#include "transport.h"
#include "input.h"

using std::cout;
using std::endl;


void print_elapsed_inside(const char* desc, struct timeval* start, struct timeval* end) 
{
  //prints timing statistics
  struct timeval elapsed;
  /* calculate elapsed time */
  if(start->tv_usec > end->tv_usec) 
  {
    end->tv_usec += 1000000;
    end->tv_sec--;
  }
  elapsed.tv_usec = end->tv_usec - start->tv_usec;
  elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
  printf(" %s  %f \n", desc, (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.0 );
}



int main(int argc, char *argv[] ) 
{

  if (argc != 2) {
    cout<<"Usage: BRANSON <path_to_input_file>"<<endl;
    exit(EXIT_FAILURE); 
  }  


 //timing 
 struct timeval start,end;
 struct timezone tzp;
 gettimeofday(&start, &tzp); 


  //Input
  std::string filename( argv[1]);
  Input *input;
  input = new Input(filename);
  input->print_problem_info();

  //IMC state setup
  IMC_State *imc_state;
  imc_state = new IMC_State(input);

  Mesh* mesh = new Mesh(input);
  //mesh->print_mesh_info();
  Photon *census_list;

  vector<double> abs_E(mesh->get_num_elems(), 0.0) ;
  vector<unsigned int> element_census_count(mesh->get_num_elems(), 0);

  while (!imc_state->finished())
  {
    imc_state->print_timestep_header();

    //set opacity, Fleck factor all energy to source
    mesh->calculate_photon_energy(imc_state);

    //make and transport photons in paralell 
    make_photons(mesh, imc_state, abs_E, element_census_count, census_list);

    //cout<<"updating temperature..."<<endl;
    mesh->update_temperature(abs_E, imc_state);

    //mesh->print_mesh_info();
    imc_state->print_conservation();

    //update time for next step
    imc_state->next_time_step();
  }
  cout<<"****************************************";
  cout<<"****************************************"<<endl;

  gettimeofday(&end, &tzp);
  print_elapsed_inside("\nruntime:",&start, &end);
  delete mesh;
  delete imc_state;
  delete input;

  return 0;
}
