/* decompose_mesh.h
 * This file includes functions to decompose the 
 * mesh with Zoltan.
 * 6/17/2015 by Alex Long
*/

#ifndef decompose_mesh_h_
#define decompose_mesh_h_

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp> //!< This provides serialization for std::map
#include <parmetis.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>


namespace mpi = boost::mpi;

void print_MPI_out(Mesh *mesh, unsigned int rank, unsigned int size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (unsigned int p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      mesh->print();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}


void print_MPI_maps(Mesh *mesh, unsigned int rank, unsigned int size) {
  using std::cout;
  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  for (unsigned int p_rank = 0; p_rank<size; p_rank++) {
    if (rank == p_rank) {
      mesh->print_map();
      cout.flush();
    }
    usleep(100);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
  }
}

void decompose_mesh(Mesh* mesh, mpi::communicator world, int argc, char **argv) {

  using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
  using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
  using std::vector;
  using std::map;
  using std::partial_sum;
  
  unsigned int rank = MPI::COMM_WORLD.Get_rank();
  unsigned int nrank = MPI::COMM_WORLD.Get_size();
  MPI_Comm comm = MPI::COMM_WORLD.Dup();

  //begin PARMETIS routines

  //get the number of elements on each processor
  //vtxdist has number of vertices on each rank, same for all ranks
  vector<int> start_nelements(nrank, 0);
  vector<int> vtxdist(nrank, 0);
  int nelem_on_rank = mesh->get_number_of_objects();
  mpi::all_gather(world, nelem_on_rank, start_nelements);
  partial_sum(start_nelements.begin(), 
              start_nelements.end(), 
              vtxdist.begin());
  vtxdist.insert(vtxdist.begin(), 0); 

  //build adjacency lists for each rank
  vector<int> xadj;
  vector<int> adjncy;
  int adjncy_ctr = 0;
  Element elem;
  unsigned int g_ID; //! Global ID
  for (unsigned int i=0; i<mesh->get_number_of_objects();i++) {
    elem = mesh->get_pre_elem(i);
    g_ID = elem.get_ID();
    unsigned int xm_neighbor =elem.get_next_element(X_NEG);
    unsigned int xp_neighbor =elem.get_next_element(X_POS);
    unsigned int ym_neighbor =elem.get_next_element(Y_NEG);
    unsigned int yp_neighbor =elem.get_next_element(Y_POS);
    unsigned int zm_neighbor =elem.get_next_element(Z_NEG);
    unsigned int zp_neighbor =elem.get_next_element(Z_POS);
    
    xadj.push_back(adjncy_ctr); //starting index in xadj for this element's nodes
    if (xm_neighbor != g_ID) {adjncy.push_back(xm_neighbor); adjncy_ctr++;}
    if (xp_neighbor != g_ID) {adjncy.push_back(xp_neighbor); adjncy_ctr++;}
    if (ym_neighbor != g_ID) {adjncy.push_back(ym_neighbor); adjncy_ctr++;}
    if (yp_neighbor != g_ID) {adjncy.push_back(yp_neighbor); adjncy_ctr++;}
    if (zm_neighbor != g_ID) {adjncy.push_back(zm_neighbor); adjncy_ctr++;}
    if (zp_neighbor != g_ID) {adjncy.push_back(zp_neighbor); adjncy_ctr++;}
  }
  xadj.push_back(adjncy_ctr);  

  int wgtflag = 0; //no weights for elements
  int numflag = 0; //C-style numbering
  int ncon = 1; 
  int nparts = nrank; //sub-domains = nrank

  float *tpwgts = new float[nparts];
  for (unsigned int i=0; i<nparts; i++) tpwgts[i]=1.0/nparts;

  float *ubvec = new float[ncon];
  for (unsigned int i=0; i<ncon; i++) ubvec[i]=1.05;

  int options[3];
  options[0] = 1; // 0--use default values, 1--use the values in 1 and 2
  options[1] = 0; //output level
  options[2] = 1242; //random number seed
  
  int edgecut =0; 
  int *part = new int[nelem_on_rank];

  ParMETIS_V3_PartKway( &vtxdist[0],   // array describing how elements are distributed
                        &xadj[0],   // how elements are stored locally
                        &adjncy[0], // how elements are stored loccaly
                        NULL,       // weight of vertices
                        NULL,       // weight of edges
                        &wgtflag,   // 0 means no weights for node or edges
                        &numflag,   // numbering style, 0 for C-style
                        &ncon,      // weights per vertex
                        &nparts,    // number of sub-domains
                        tpwgts,     // weight per sub-domain
                        ubvec,      // unbalance in vertex weight
                        options,    // options array
                        &edgecut,   // OUTPUT: Number of edgecuts
                        part,       // OUTPUT: partition of each vertex
                        &comm); // MPI communicator

  //send elements to other processors
  for (unsigned int send_rank =0; send_rank<nrank; send_rank++) {
    for (unsigned int recv_rank =0; recv_rank<nrank; recv_rank++) {
      if (  (send_rank != recv_rank)  && (rank == send_rank || rank == recv_rank) ) {
        if(rank == send_rank) {
          vector<Element> send_list;
          for (unsigned int i=0; i<nelem_on_rank; i++) {
            if(part[i] == recv_rank)
              send_list.push_back(mesh->get_pre_elem(i));
          }
          world.send(recv_rank, 0, send_list);
          //Erase these elements from the mesh
          for (unsigned int i=0; i<nelem_on_rank; i++) {
            if(part[i] == recv_rank) mesh->remove_elem(i);
          }
        }
        // rank == recv_rank
        else {
          vector<Element> recv_list;
          world.recv(send_rank, 0, recv_list);
          // add these elements to the mesh
          for (unsigned int i = 0; i< recv_list.size(); i++) 
            mesh->add_mesh_elem(recv_list[i]);
        }
      } //send_rank != recv_rank
    } // loop over recv_ranks
  } // loop over send ranks

  //update the cell list on each processor
  mesh->update_mesh();

  //get the number of elements on each processor
  vector<unsigned int> out_elements_proc(nrank, 0);
  vector<unsigned int> prefix_elements_proc(nrank, 0);
  unsigned int n_elem = mesh->get_number_of_objects();
  mpi::all_gather(world, n_elem, out_elements_proc);
  partial_sum(out_elements_proc.begin(), out_elements_proc.end(), prefix_elements_proc.begin());

  unsigned int g_start = prefix_elements_proc[rank]-n_elem;
  unsigned int g_end = prefix_elements_proc[rank]-1;
  mesh->set_global_bound(g_start, g_end);
  //append zero to the prefix array to make it a standard bounds array
  prefix_elements_proc.insert( prefix_elements_proc.begin(), 0);
  mesh->set_off_rank_bounds(prefix_elements_proc);

  //make sure each index is remapped ONLY ONCE!
  vector< vector<bool> > remap_flag;
  for (unsigned int i=0; i<n_elem; i++) remap_flag.push_back(vector<bool> (6,false));

  //change global indices to match a simple number system for easy sorting,
  //this involves sending maps to each processor to get new indicies
  map<unsigned int, unsigned int> local_map = mesh->get_map();
  // Send maps
  for (unsigned int send_rank =0; send_rank<nrank; send_rank++) {
    for (unsigned int recv_rank =0; recv_rank<nrank; recv_rank++) {
      if (  (send_rank != recv_rank)  && (rank == send_rank || rank == recv_rank) ) {
        if(rank == send_rank) {
          world.send(recv_rank, 0, local_map);
        }
        // rank == recv_rank
        else {
          map<unsigned int, unsigned int> off_map;
          world.recv(send_rank, 0, off_map);
          //remap off processor indices
          mesh->set_indices(off_map, remap_flag);
        }
      } //send_rank != recv_rank
    } // loop over recv_ranks
  } // loop over send ranks

  //now update the indices of local IDs
  mesh->set_local_indices(local_map);

  //reallocate mesh data in new MPI window
  mesh->make_MPI_window();
}

#endif // decompose_mesh_h
