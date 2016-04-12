/* decompose_photons.h
 * This file includes functions to load balance census photons after a step
 * in the mesh-passing DD algotithm
 * 6/17/2015 by Alex Long
*/

#ifndef decompose_photons_h_
#define decompose_photons_h_

#include <iostream>
#include <mpi.h>
#include <map>
#include <algorithm>
#include <vector>

#include "photon.h"


void print_MPI_photons( const std::vector<Photon>& phtn_vec, 
                        const uint32_t& rank, 
                        const uint32_t& size) {

  using std::cout;

  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (uint32_t p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      for(uint32_t i=0; i<phtn_vec.size();i++)
        phtn_vec[i].print_info(rank);
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
  usleep(100);
  cout.flush();
  usleep(100);
}


std::vector<Photon> rebalance_census(std::vector<Photon>& off_rank_census,
                                     Mesh* mesh)
{
  using std::map;
  using std::sort;
  using std::vector;

  int rank, n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);
  uint32_t n_off_rank = n_rank-1;

  // make the MPI Particle 
  const int entry_count = 2 ; 
  // 7 uint32_t, 6 int, 13 double
  int array_of_block_length[3] = { 2, 9};
  // Displacements of each type in the cell
  MPI_Aint array_of_block_displace[2] = 
    {0, 2*sizeof(uint32_t)};
  //Type of each memory block
  MPI_Datatype array_of_types[2] = {MPI_UNSIGNED, MPI_DOUBLE};

  MPI_Datatype MPI_Particle;
  MPI_Type_create_struct(entry_count, array_of_block_length, 
    array_of_block_displace, array_of_types, &MPI_Particle);

  // Commit the type to MPI so it recognizes it in communication calls
  MPI_Type_commit(&MPI_Particle);

  int mpi_particle_size;
  MPI_Type_size(MPI_Particle, &mpi_particle_size);
 
  int particle_size = sizeof(Photon); 

  // make off processor map
  map<int,int> proc_map;
  for (int32_t i=0; i<n_off_rank; i++) {
    int r_index = i + int(i>=rank);
    proc_map[i] = r_index;
  }

  //sort the census vector by cell ID (global ID)
  sort(off_rank_census.begin(), off_rank_census.end());
 
  //size of census list 
  uint32_t n_census = off_rank_census.size();

  //count the photons belonging to each rank and the start index of each
  //count the ranks that you will send to, add them to a vector
  vector<uint32_t> rank_count(n_rank, 0);
  vector<uint32_t> rank_start(n_rank+1, 0);
  vector<bool> rank_found(n_rank, false);
  uint32_t r;
  for (uint32_t i=0; i<n_census; i++) {
    r = mesh->get_rank(off_rank_census[i].get_cell());
    rank_count[r]++;
    if(rank_found[r]==false) {
      rank_found[r]=true;
      rank_start[r] =i;
    }
  }

  // end of rank count is the total number of census photons
  rank_start[n_rank] = n_census;

  // make requests for non-blocking communication
  MPI_Request* reqs = new MPI_Request[n_off_rank*2];

  // make n_off_rank receive buffers
  vector<vector<Photon> > recv_photons;
  for (uint32_t ir=0; ir<n_off_rank; ir++) {
    vector<Photon> empty_vec;
    recv_photons.push_back(empty_vec);
  }

  //get the number of photons received from each rank
  vector<int> recv_from_rank(n_off_rank, 0);

  for (uint32_t ir=0; ir<n_off_rank; ir++) {
    int off_rank = proc_map[ir];
    MPI_Isend(&rank_count[off_rank], 1,  MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD, &reqs[ir]);
    MPI_Irecv(&recv_from_rank[ir], 1, MPI_UNSIGNED, off_rank, 0, MPI_COMM_WORLD, &reqs[ir+n_off_rank]);
  }
  
  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE); 

  // now send the buffers and post receives  
  // resize receive buffers with recv_from_rank
  for (uint32_t ir=0; ir<n_off_rank; ir++) {
    int off_rank = proc_map[ir];
    int start_copy = rank_start[off_rank];
    MPI_Isend(&off_rank_census[start_copy], 
              rank_count[off_rank],
              MPI_Particle,
              off_rank, 
              0, 
              MPI_COMM_WORLD, 
              &reqs[ir]);
    recv_photons[ir].resize(recv_from_rank[ir]);
    MPI_Irecv(&recv_photons[ir][0],
              recv_from_rank[ir], 
              MPI_Particle, 
              off_rank, 
              0, 
              MPI_COMM_WORLD, 
              &reqs[ir+n_off_rank]);
  }

  MPI_Waitall(n_off_rank*2, reqs, MPI_STATUS_IGNORE); 

  //free memory from off rank census list
  off_rank_census.clear();

  //copy received census photons to a new census list
  vector<Photon> new_on_rank_census;
  for (uint32_t ir=0; ir<n_rank-1; ir++) {
    new_on_rank_census.insert(new_on_rank_census.end(), 
      recv_photons[ir].begin(), 
      recv_photons[ir].end());
  }

  //Explicitly delete the MPI requests
  delete[] reqs;

  return new_on_rank_census;
}


/*
void proto_load_balance_photons(Photon*& photon_vec, 
                                uint32_t& n_photon,
                                Mesh *mesh, 
                                mpi::communicator world) {

  using std::vector;

  uint32_t rank = MPI::COMM_WORLD.Get_rank();
  uint32_t size = MPI::COMM_WORLD.Get_size();


  //sort the census vector by cell ID (global ID)
  std::sort(photon_vec, photon_vec+n_photon);

  uint32_t g_n_photon = n_photon;
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &g_n_photon, 1, MPI_UNSIGNED, MPI_SUM);
  uint32_t lb_photon(double(g_n_photon)/size);
  uint32_t send_photons = 0;
  if ( n_photon > lb_photon) send_photons = n_photon - lb_photon;
  uint32_t send_per_rank(double(send_photons)/(size-1));

  //send photons to other processors
  mpi::request* reqs = new mpi::request[ (size-1)*2];
  vector<vector<Photon> > recv_photons;
  //make size-1 receive lists
  for (uint32_t ir=0; ir<size-1; ir++) {
    vector<Photon> empty_vec;
    recv_photons.push_back(empty_vec);
  }

  uint32_t icount = 0;
  for (uint32_t ir=0; ir<size; ir++) {
    if (rank != ir) {
      //build send list
      vector<Photon> send_photons(send_per_rank);
      //get correct index into received photon vector
      uint32_t r_index = ir - (ir>rank);
      uint32_t start_send = r_index*send_per_rank;
      memcpy(&send_photons[0], photon_vec+start_send, sizeof(Photon)*send_per_rank);
      reqs[icount] = world.isend(ir, 0, send_photons);
      icount++;
      reqs[icount] = world.irecv(ir,0,recv_photons[r_index]);
      icount++;
    }
  }

  mpi::wait_all(reqs, reqs+(size-1)*2);


  //make new photon list
  //first get the size of it
  uint32_t n_remain = n_photon - send_per_rank*(size-1);
  uint32_t new_photon_size = 0;
  new_photon_size +=  n_remain;
  for (uint32_t ir=0; ir<size-1; ir++) 
    new_photon_size+=recv_photons[ir].size();

  //now form the total photon list for transport
  Photon *new_photon_vec = new Photon[new_photon_size];

  //copy received photon vector on to new photon vector
  uint32_t copy_start = 0; 
  for (uint32_t ir=0; ir<size-1; ir++) {
    memcpy(new_photon_vec+copy_start, &recv_photons[ir][0], sizeof(Photon)*recv_photons[ir].size());
    copy_start+=recv_photons[ir].size();
    recv_photons[ir].clear();
  }

  //copy over old photons on to new vector
  uint32_t start_remain = send_per_rank*(size-1);
  memcpy(new_photon_vec+copy_start, photon_vec+start_remain, sizeof(Photon)*n_remain);
  delete[] photon_vec;

  //Explicitly delete the MPI requests
  delete[] reqs;

  photon_vec = new_photon_vec;
  n_photon = new_photon_size;
  //print_MPI_photons(photon_vec, rank, size);
}
*/
#endif // decompose_photons_h
