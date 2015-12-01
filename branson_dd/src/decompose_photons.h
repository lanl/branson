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

/*
static int get_number_of_photons(void *data, int *ierr) {
  vector<Photon>* photon_vec = (vector<Photon> *) data;
  *ierr = ZOLTAN_OK;
  return photon_vec->size();
}

static void get_photon_list(void *data, int sizeGID, int sizeLID, 
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
            int wgt_dim, float *obj_wgts, int *ierr) {

  unsigned int *start_id_and_total = (unsigned int*) data;
  unsigned int start_id =start_id_and_total[0];
  unsigned int total =start_id_and_total[1];
  *ierr = ZOLTAN_OK;
  for (unsigned int i=0; i<total; i++){
    globalID[i] = start_id + i; 
    localID[i] = i;
  }
}


static int get_num_geometry_photon(void *data, int *ierr)
{
  // 3 points are needed to identify points on mesh
  *ierr = ZOLTAN_OK;
  return 3;
}


static void get_geometry_list_photon(void *data, int sizeGID, int sizeLID, int num_obj,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr)
{
  vector<Photon>* photon_vec = (vector<Photon> *) data;

  if (num_dim != 3) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  *ierr = ZOLTAN_OK;

  const double* pos; 
  
  for (unsigned int i=0; i< num_obj; i++) {
    Photon& phtn = (*photon_vec)[i];
    pos = phtn.get_position();
    geom_vec[3*i] = pos[0];
    geom_vec[3*i+1] = pos[1];
    geom_vec[3*i+2] = pos[2];
  }
  return;
}
*/


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


/*
vector<Photon> decompose_photons_zoltan(vector<Photon>& photon_vec,  mpi::communicator world, int argc, char **argv) {

  //Dynamically create Zoltan object.
  Zoltan *zz = new Zoltan(MPI::COMM_WORLD);
  zz->Set_Param("DEBUG_LEVEL", "0");
  zz->Set_Param("LB_METHOD", "RCB");

  unsigned int rank = MPI::COMM_WORLD.Get_rank();
  unsigned int size = MPI::COMM_WORLD.Get_size();
  int rc;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  //Get the correct numbering scheme for the photons global IDs with an allreduce
  vector<unsigned int> out_photons_proc(world.size(), 0);
  vector<unsigned int> prefix_photons_proc(world.size(), 0);
  unsigned int n_photons = photon_vec.size();
  mpi::all_gather(world, n_photons, out_photons_proc);
  partial_sum(out_photons_proc.begin(), out_photons_proc.end(), prefix_photons_proc.begin());
  unsigned int start_id_and_total[2] = { prefix_photons_proc[rank]-n_photons, n_photons};

  //Initialize the Zoltan library with a C language call
  float version;
  rc = Zoltan_Initialize(argc, argv, &version);

  zz->Set_Num_Obj_Fn(get_number_of_photons, &photon_vec);
  zz->Set_Obj_List_Fn(get_photon_list, start_id_and_total);
  zz->Set_Num_Geom_Fn(get_num_geometry_photon, &photon_vec);
  zz->Set_Geom_Multi_Fn(get_geometry_list_photon, &photon_vec);

  rc = zz->LB_Partition(changes,        // 1 if partitioning was changed, 0 otherwise 
                        numGidEntries,  // Number of integers used for a global ID 
                        numLidEntries,  // Number of integers used for a local ID 
                        numImport,      // Number of vertices to be sent to me 
                        importGlobalGids,  // Global IDs of vertices to be sent to me 
                        importLocalGids,   // Local IDs of vertices to be sent to me 
                        importProcs,    // Process rank for source of each incoming vertex 
                        importToPart,   // New partition for each incoming vertex 
                        numExport,      // Number of vertices I must send to other processes
                        exportGlobalGids,  // Global IDs of the vertices I must send 
                        exportLocalGids,   // Local IDs of the vertices I must send 
                        exportProcs,    // Process to which I send each of the vertices 
                        exportToPart);  // Partition to which each vertex will belong 


  //send cells to other processors
  //assume all photons are staying
  vector<bool> stay_flag(n_photons, true);
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
      vector<Photon> send_photons;
      for (unsigned int i=0; i<numExport; i++) {
        if(exportProcs[i] == ir) {
          send_photons.push_back( photon_vec[exportLocalGids[i]] );
          stay_flag[exportLocalGids[i]] = false;
        }
      }
      reqs[icount] = world.isend(ir, 0, send_photons);
      icount++;
      //get correct index into received photon vector
      unsigned int r_index = ir - (ir>rank);
      reqs[icount] = world.irecv(ir,0,recv_photons[r_index]);
    }
  }
  mpi::wait_all(reqs, reqs+(size-1)*2);

  //make new photon vector
  vector<Photon> new_photon_vec;
  for (unsigned int i=0; i<n_photons; i++) {
    if (stay_flag[i]) new_photon_vec.push_back(photon_vec[i]);
  }
  //push received photon vector on to new photon vector
  for (unsigned int ir=0; ir<size-1; ir++) 
    new_photon_vec.insert(new_photon_vec.end(), recv_photons[ir].begin(), recv_photons[ir].end());
  
  //Explicitly delete the Zoltan object and MPI requests
  delete zz;
  delete[] reqs;


  //print_MPI_photons(photon_vec, rank, size);

  return new_photon_vec;
}
*/


void on_rank_rebalance_photons(Photon*& photon_vec, 
                               unsigned int& n_photon,
                               Photon*& census_list,
                               const unsigned int& n_census,
                               Mesh *mesh, 
                               mpi::communicator world) {

  using std::vector;

  unsigned int rank = MPI::COMM_WORLD.Get_rank();
  unsigned int size = MPI::COMM_WORLD.Get_size();

  //sort the census vector by cell ID (global ID)
  std::sort(census_list, census_list+ n_census);

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

  //finding exactly who you're receiving from requires an all to all 
  //communication. I don't want to do that right now. I'll allow empty
  //messages
  //
  /*
  //reduce the rank_found array to get the number of ranks to receive from
  vector<unsigned int> recv_message(world.size(),0);
  for (unsigned int ir=0; ir<size; ir++) 
    if(rank_found[ir]==true) recv_message[ir]=1;
  recv_message[rank]=0; // can't receive message from self
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &recv_message[0], size, MPI_UNSIGNED, MPI_SUM);
  */

  //end of rank count is the total number of census photons
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
  for (unsigned int ir=0; ir<size; ir++) {
    if (rank != ir) {
      //build send list
      vector<Photon> send_photons(rank_count[ir]);
      unsigned int start_send;
      start_send = rank_start[ir];
      memcpy(&send_photons[0], census_list+start_send, sizeof(Photon)*rank_count[ir]);
      reqs[icount] = world.isend(ir, 0, send_photons);
      icount++;
      //get correct index into received photon vector
      unsigned int r_index = ir - (ir>rank);
      reqs[icount] = world.irecv(ir,0,recv_photons[r_index]);
      icount++;
    }
  }

  mpi::wait_all(reqs, reqs+(size-1)*2);

  //make new census list
  //first get the size of it
  unsigned int new_census_size = 0;
  new_census_size += rank_count[rank];
  for (unsigned int ir=0; ir<size-1; ir++) 
    new_census_size+=recv_photons[ir].size();

  //keep memory use low, copy over the old, then delete the census list  
  Photon *on_rank_census = new Photon[rank_count[rank]];
  unsigned int start_on_rank = rank_start[rank];
  memcpy(on_rank_census, census_list+start_on_rank, sizeof(Photon)*rank_count[rank]);
  delete[] census_list;
  
  //now form the total photon list for transport
  //alllocate new photon vector
  Photon *new_photon_vec = new Photon[n_photon+new_census_size];
  memcpy(new_photon_vec, photon_vec, sizeof(Photon)*n_photon);
  delete[] photon_vec;

  //copy on rank census  
  memcpy(new_photon_vec+n_photon, on_rank_census, sizeof(Photon)*rank_count[rank]);
  delete[] on_rank_census;

  //push received photon vector on to new photon vector
  unsigned int copy_start = n_photon+rank_count[rank];
  for (unsigned int ir=0; ir<size-1; ir++) {
    memcpy(new_photon_vec+copy_start, &recv_photons[ir][0], sizeof(Photon)*recv_photons[ir].size());
    copy_start+=recv_photons[ir].size();
    recv_photons[ir].clear();
  }

  //Explicitly delete the MPI requests
  delete[] reqs;

  photon_vec = new_photon_vec;
  n_photon = n_photon+new_census_size;
  //print_MPI_photons(photon_vec, rank, size);
}



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

#endif // decompose_photons_h
