/* decompose_mesh.g
 * This file includes functions to decompose the 
 * mesh with Zoltan.
 * 6/17/2015 by Alex Long
*/

#ifndef decompose_photons_h_
#define decompose_photons_h_

#include <iostream>
#include <boost/mpi.hpp>
#include <algorithm>
#include <vector>

#include "photon.h"


namespace mpi = boost::mpi;

void print_MPI_photons(const std::vector<Photon>& phtn_vec, const unsigned int& rank, const unsigned int& size) {

  using std::cout;

  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (unsigned int p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      for(unsigned int i=0; i<phtn_vec.size();i++)
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


std::vector<Photon> rebalance_census(std::vector<Photon>& census_list,
                                     Mesh* mesh, 
                                     mpi::communicator world) 
{
  using std::vector;
  using std::sort;

  unsigned int rank = MPI::COMM_WORLD.Get_rank();
  unsigned int size = MPI::COMM_WORLD.Get_size();

  //sort the census vector by cell ID (global ID)
  sort(census_list.begin(), census_list.end());
  
  unsigned int n_census = census_list.size();

  //count the photons belonging to each rank and the start index of each
  //count the ranks that you will send to, add them to a vector
  vector<unsigned int> rank_count(size, 0);
  vector<unsigned int> rank_start(size+1, 0);
  vector<bool> rank_found(size, false);
  unsigned int r;
  for (unsigned int i=0; i<n_census; i++) {
    r = mesh->get_rank(census_list[i].get_cell());
    rank_count[r]++;
    if(rank_found[r]==false) {
      rank_found[r]=true;
      rank_start[r] =i;
    }
  }

  // end of rank count is the total number of census photons
  rank_start[size] = n_census;

  //send photons to other processors
  mpi::request* reqs = new mpi::request[ (size-1)*2];
  vector<vector<Photon> > recv_photons;
  //make size-1 receive lists
  for (unsigned int ir=0; ir<size-1; ir++) {
    vector<Photon> empty_vec;
    recv_photons.push_back(empty_vec);
  }

  unsigned int icount = 0;
  vector<Photon>::iterator start = census_list.begin();
  for (unsigned int ir=0; ir<size; ir++) {
    if (rank != ir) {
      //build send list
      vector<Photon> send_photons(rank_count[ir]);
      unsigned int start_send;
      start_send = rank_start[ir];
      memcpy(&send_photons[0], &census_list[start_send], sizeof(Photon)*rank_count[ir]);
      reqs[icount] = world.isend(ir, 0, send_photons);
      icount++;
      //get correct index into received photon vector
      unsigned int r_index = ir - (ir>rank);
      reqs[icount] = world.irecv(ir,0,recv_photons[r_index]);
      icount++;
    }
  }

  mpi::wait_all(reqs, reqs+(size-1)*2);

  //copy on rank census photons to new census list
  vector<Photon> new_census_list;
  unsigned int start_on_rank = rank_start[rank];
  new_census_list.insert(new_census_list.end(),
                         start+start_on_rank,
                         start+(start_on_rank+rank_count[rank]));

  //free memory from census list
  census_list.clear();

  for (unsigned int ir=0; ir<size-1; ir++)
    new_census_list.insert(new_census_list.end(), 
      recv_photons[ir].begin(), 
      recv_photons[ir].end());

  //Explicitly delete the MPI requests
  delete[] reqs;

  sort(new_census_list.begin(), new_census_list.end());

  return new_census_list;
}


/*
void proto_load_balance_photons(Photon*& photon_vec, 
                                unsigned int& n_photon,
                                Mesh *mesh, 
                                mpi::communicator world) {

  using std::vector;

  unsigned int rank = MPI::COMM_WORLD.Get_rank();
  unsigned int size = MPI::COMM_WORLD.Get_size();


  //sort the census vector by cell ID (global ID)
  std::sort(photon_vec, photon_vec+n_photon);

  unsigned int g_n_photon = n_photon;
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &g_n_photon, 1, MPI_UNSIGNED, MPI_SUM);
  unsigned int lb_photon(double(g_n_photon)/size);
  unsigned int send_photons = 0;
  if ( n_photon > lb_photon) send_photons = n_photon - lb_photon;
  unsigned int send_per_rank(double(send_photons)/(size-1));

  //send photons to other processors
  mpi::request* reqs = new mpi::request[ (size-1)*2];
  vector<vector<Photon> > recv_photons;
  //make size-1 receive lists
  for (unsigned int ir=0; ir<size-1; ir++) {
    vector<Photon> empty_vec;
    recv_photons.push_back(empty_vec);
  }

  unsigned int icount = 0;
  for (unsigned int ir=0; ir<size; ir++) {
    if (rank != ir) {
      //build send list
      vector<Photon> send_photons(send_per_rank);
      //get correct index into received photon vector
      unsigned int r_index = ir - (ir>rank);
      unsigned int start_send = r_index*send_per_rank;
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
  unsigned int n_remain = n_photon - send_per_rank*(size-1);
  unsigned int new_photon_size = 0;
  new_photon_size +=  n_remain;
  for (unsigned int ir=0; ir<size-1; ir++) 
    new_photon_size+=recv_photons[ir].size();

  //now form the total photon list for transport
  Photon *new_photon_vec = new Photon[new_photon_size];

  //copy received photon vector on to new photon vector
  unsigned int copy_start = 0; 
  for (unsigned int ir=0; ir<size-1; ir++) {
    memcpy(new_photon_vec+copy_start, &recv_photons[ir][0], sizeof(Photon)*recv_photons[ir].size());
    copy_start+=recv_photons[ir].size();
    recv_photons[ir].clear();
  }

  //copy over old photons on to new vector
  unsigned int start_remain = send_per_rank*(size-1);
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
