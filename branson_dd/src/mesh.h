
#ifndef mesh_h_
#define mesh_h_

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>

#include "mpi.h"
#include "input.h"
#include "imc_state.h"
#include "cell.h"
#include "constants.h"
#include "buffer.h"

namespace mpi = boost::mpi;

class Mesh {

  public:

  Mesh(Input* input, int _rank, int _n_rank)
  : ngx(input->get_n_x_cells()),
    ngy(input->get_n_y_cells()),
    ngz(input->get_n_z_cells()),
    rank(_rank),
    n_rank(_n_rank)
  {
    using std::vector;
    using Constants::bc_type;
    using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
    using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
    using Constants::ELEMENT;

    max_map_size = input->get_map_size();
    double dx = input->get_dx();
    double dy = input->get_dy();
    double dz = input->get_dz();

    vector<bc_type> bc(6);
    bc[X_POS] = input->get_bc(X_POS); bc[X_NEG] = input->get_bc(X_NEG);
    bc[Y_POS] = input->get_bc(Y_POS); bc[Y_NEG] = input->get_bc(Y_NEG);
    bc[Z_POS] = input->get_bc(Z_POS); bc[Z_NEG] = input->get_bc(Z_NEG);

    //initialize number of DMA requests as zero
    off_rank_reads=0;

    uint32_t g_count =0; //global count

    //this rank's cells
    n_global = ngx*ngy*ngz;
    uint32_t cell_id_begin = floor(rank*n_global/double(n_rank));
    uint32_t cell_id_end = floor((rank+1)*n_global/double(n_rank));

    uint32_t l_count =0;
    for (uint32_t k=0; k<ngz; k++) {
      for (uint32_t j=0; j<ngy; j++) {
        for (uint32_t i=0; i<ngx; i++) {
          if (g_count >= cell_id_begin && g_count < cell_id_end) {
            //global_ID.push_back(count*10 + MPI::COMM_WORLD.Get_rank()  );
            Cell e;
            e.set_coor(i*dx, (i+1)*dx, j*dy, (j+1)*dy, k*dz, (k+1)*dz);
            e.set_ID(g_count);
            e.set_cV(input->get_CV());
            e.set_T_e(input->get_initial_Tm()); 
            e.set_T_r(input->get_initial_Tr());
            e.set_T_s(0.0);
            e.set_rho(input->get_rho());

            if (i<(ngx-1)) {
              e.set_neighbor( X_POS, g_count+1); e.set_bc(X_POS, ELEMENT);
            }
            else {
              e.set_neighbor( X_POS, g_count); e.set_bc(X_POS, bc[X_POS]);
            } 

            if (i>0) {
              e.set_neighbor( X_NEG, g_count-1); e.set_bc(X_NEG, ELEMENT);
            }
            else {
              e.set_neighbor( X_NEG, g_count); e.set_bc(X_NEG, bc[X_NEG]);
            }
            if (j<(ngy-1)) {
              e.set_neighbor(Y_POS, g_count+ngx); e.set_bc(Y_POS, ELEMENT);
            }
            else {
              e.set_neighbor(Y_POS, g_count); e.set_bc(Y_POS, bc[Y_POS]);
            }
            if (j>0) {
              e.set_neighbor(Y_NEG, g_count-ngx); e.set_bc(Y_NEG, ELEMENT);
            }
            else {
              e.set_neighbor(Y_NEG, g_count); e.set_bc(Y_NEG, bc[Y_NEG]);
            }
            if (k<(ngz-1)) {
              e.set_neighbor(Z_POS, g_count+ngx*ngy); e.set_bc(Z_POS, ELEMENT);
            }
            else {
              e.set_neighbor(Z_POS, g_count); e.set_bc(Z_POS, bc[Z_POS]);
            }
            if (k>0) {
              e.set_neighbor(Z_NEG, g_count-ngx*ngy); e.set_bc(Z_NEG, ELEMENT);
            }
            else {
              e.set_neighbor(Z_NEG, g_count); e.set_bc(Z_NEG, bc[Z_NEG]);
            }
            cell_list.push_back(e);
            l_count++;
          }
          g_count++;
        } //end i loop
      } // end j loop
    } // end k loop
    n_cell = l_count;
    m_opA = input->get_opacity_A();
    m_opB = input->get_opacity_B();
    m_opC = input->get_opacity_C();
    m_opS = input->get_opacity_S();

    total_photon_E = 0.0;

    //define the MPI cell type used for MPI windows
    //make the MPI datatype for my cell class
    // Three type entries in the class
    const int entry_count = 3 ; 
    // 7 uint32_t, 6 int, 13 double
    int array_of_block_length[4] = {8, 6, 14};
    // Displacements of each type in the cell
    MPI_Aint array_of_block_displace[3] = 
      {0, 8*sizeof(uint32_t),  8*sizeof(uint32_t)+6*sizeof(int)};
    //Type of each memory block
    MPI_Datatype array_of_types[3] = {MPI_UNSIGNED, MPI_INT, MPI_DOUBLE}; 

    MPI_Datatype MPI_Cell;
    MPI_Type_create_struct(entry_count, array_of_block_length, 
      array_of_block_displace, array_of_types, &MPI_Cell);

    MPI_Type_commit(&MPI_Cell);

    MPI_Type_size(MPI_Cell, &mpi_cell_size);

    //bool flags to say if data is needed
    need_data = vector<bool>(n_rank-1, false);

    //allocate mpi requests objects
    r_cell_reqs = new mpi::request[n_rank-1];
    r_id_reqs = new mpi::request[n_rank-1];
    s_cell_reqs = new mpi::request[n_rank-1];
    s_id_reqs = new mpi::request[n_rank-1];
   
    recv_cell_buffer = vector<Buffer<Cell> >(n_rank-1);
    send_cell_buffer = vector<Buffer<Cell> >(n_rank-1);
    recv_id_buffer =  vector<Buffer<uint32_t> >(n_rank-1);
    send_id_buffer =  vector<Buffer<uint32_t> >(n_rank-1);

    // size the send and receive buffers
    // they are size n_rank -1
    for (uint32_t ir=0; ir<n_rank-1; ir++) {
      vector<uint32_t> empty_vec_ids;
      //ids needed from other ranks
      ids_needed.push_back(empty_vec_ids); 
    }
  }

  //free buffers and delete MPI allocated cell
  ~Mesh() {
    delete[] r_cell_reqs;
    delete[] r_id_reqs;
    delete[] s_cell_reqs;
    delete[] s_id_reqs;
  }

/*****************************************************************************/
  //const functions
/*****************************************************************************/
  uint32_t get_n_local_cells(void) const {return n_cell;}
  uint32_t get_my_rank(void) const {return  rank;}
  uint32_t get_offset(void) const {return on_rank_start;}
  uint32_t get_global_num_cells(void) const {return n_global;}
  std::map<uint32_t, uint32_t> get_proc_adjacency_list(void) const {
    return adjacent_procs;
  }
  double get_total_photon_E(void) const {return total_photon_E;}

  void pre_decomp_print(void) const {
    for (uint32_t i= 0; i<n_cell; i++)
      cell_list[i].print();
  }

  void post_decomp_print(void) const {
    for (uint32_t i= 0; i<n_cell; i++) 
      cells[i].print();
  }

  std::map<uint32_t, uint32_t> get_map(void) const {
    std::map<uint32_t, uint32_t> local_map;
    uint32_t g_ID;
    for (uint32_t i=0; i<n_cell; i++) {
      g_ID = cell_list[i].get_ID();
      local_map[g_ID] = i+on_rank_start;
    }
    return local_map;
  }

  Cell get_pre_renumber_cell(const uint32_t& local_ID) const 
  {
    return cell_list[local_ID];
  }

  Cell get_cell(const uint32_t& local_ID) const {
    return cells[local_ID];
  }

  uint32_t get_off_rank_id(const uint32_t& index) const {
    //find rank of index
    bool found = false;
    uint32_t min_i = 0;
    uint32_t max_i = off_rank_bounds.size()-1;
    uint32_t s_i; //search index
    while(!found) {
      s_i =(max_i + min_i)/2;
      if (s_i == max_i || s_i == min_i) found = true;
      else if (index >= off_rank_bounds[s_i]) min_i = s_i;
      else max_i = s_i;
    }
    return s_i;
  }

  uint32_t get_rank(const uint32_t& index) const {
    uint32_t r_rank;
    if (on_processor(index)) r_rank = rank;
    else  r_rank = get_off_rank_id(index);
    return r_rank;
  }

  uint32_t get_local_ID(const uint32_t& index) const {
    return index-on_rank_start;
  }

  Cell get_on_rank_cell(const uint32_t& index) {
    //this can only be called with valid on rank indexes
    if (on_processor(index)) 
      return cells[index-on_rank_start];
    else 
      return stored_cells[index];
  }

  bool on_processor(const uint32_t& index) const { 
    return  (index>=on_rank_start) && (index<=on_rank_end) ; 
  }

  void print_map(void) {
    for ( std::map<uint32_t,Cell>::iterator map_i =
      stored_cells.begin();
      map_i!=stored_cells.end(); map_i++)
      (map_i->second).print();
  }

  bool mesh_available(const uint32_t& index) const {
    if (on_processor(index)) return true;
    else if (stored_cells.find(index) != stored_cells.end())
      return true;
    else 
      return false; 
  } 



/*****************************************************************************/
  //non-const functions
/*****************************************************************************/
  void set_global_bound(uint32_t _on_rank_start, 
                        uint32_t _on_rank_end) 
  {
    on_rank_start = _on_rank_start;
    on_rank_end = _on_rank_end;
  }

  void set_off_rank_bounds(std::vector<uint32_t> _off_rank_bounds) {
    off_rank_bounds=_off_rank_bounds;
  }

  void calculate_photon_energy(IMC_State* imc_s) {
    using Constants::c;
    using Constants::a;
    total_photon_E = 0.0;
    double dt = imc_s->get_dt();
    double op_a, op_s, f, cV, rho;
    double vol;
    double T, Tr, Ts;
    uint32_t step = imc_s->get_step();
    double tot_census_E = 0.0;
    double tot_emission_E = 0.0;
    double tot_source_E = 0.0;
    double pre_mat_E = 0.0;
    for (uint32_t i=0; i<n_cell;++i) {
      Cell& e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      Tr = e.get_T_r();
      Ts = e.get_T_s();
      rho = e.get_rho();
      op_a = m_opA + m_opB*pow(T, m_opC);
      op_s = m_opS;
      f =1.0/(1.0 + dt*op_a*c*(4.0*a*pow(T,3)/(cV*rho)));
      e.set_op_a(op_a);
      e.set_op_s(op_s);
      e.set_f(f);

      m_emission_E[i] = dt*vol*f*op_a*a*c*pow(T,4);
      if (step > 1) m_census_E[i] = 0.0;  
      else m_census_E[i] =vol*a*pow(Tr,4); 
      m_source_E[i] = dt*op_a*a*c*pow(Ts,4);

      pre_mat_E+=T*cV*vol*rho;
      tot_emission_E+=m_emission_E[i];
      tot_census_E  +=m_census_E[i];
      tot_source_E  +=m_source_E[i];
      total_photon_E += m_emission_E[i] + m_census_E[i] + m_source_E[i];
    }

    //set conservation
    imc_s->set_pre_mat_E(pre_mat_E);
    imc_s->set_emission_E(tot_emission_E);
    imc_s->set_source_E(tot_source_E);
    if(imc_s->get_step() == 1) imc_s->set_pre_census_E(tot_census_E);
  }


  void set_indices( std::map<uint32_t, uint32_t> off_map, 
                    std::vector< std::vector<bool> >& remap_flag) {

    using Constants::PROCESSOR;
    using Constants::dir_type;
    using std::map;
    using std::set;

    uint32_t next_index;
    map<uint32_t, uint32_t>::iterator end = off_map.end();
    uint32_t new_index;
    //check to see if neighbors are on or off processor
    for (uint32_t i=0; i<n_cell; i++) {
      Cell& cell = cell_list[i];
      for (uint32_t d=0; d<6; d++) {
        next_index = cell.get_next_cell(d);
        map<uint32_t, uint32_t>::iterator map_i = 
          off_map.find(next_index);
        if (off_map.find(next_index) != end && remap_flag[i][d] ==false ) {
          //update index and bc type, this will always be an off processor so
          //if an index is updated it will always be at a processor bound
          remap_flag[i][d] = true;
          new_index = map_i->second;
          cell.set_neighbor( dir_type(d) , new_index );
          cell.set_bc(dir_type(d), PROCESSOR);
          boundary_cells.push_back(new_index);
          // add off-rank ID to map
          uint32_t off_rank = get_off_rank_id(new_index);
          if (adjacent_procs.find(off_rank) 
            == adjacent_procs.end()) {
            uint32_t rank_count = adjacent_procs.size();
            adjacent_procs[off_rank] = rank_count;
          } // if adjacent_proc.find(off_rank) 
        }
      }
    }
  }

  void set_local_indices(std::map<uint32_t, uint32_t> local_map) {

    using Constants::PROCESSOR;
    using Constants::bc_type;
    using Constants::dir_type;
    using std::map;

    uint32_t next_index;
    map<uint32_t, uint32_t>::iterator end = local_map.end();
    uint32_t new_index;
    bc_type current_bc;
    //check to see if neighbors are on or off processor
    for (uint32_t i=0; i<n_cell; i++) {
      Cell& cell = cell_list[i];
      cell.set_ID(i+on_rank_start);
      for (uint32_t d=0; d<6; d++) {
        current_bc = cell.get_bc(bc_type(d));
        next_index = cell.get_next_cell(d);
        map<uint32_t, uint32_t>::iterator map_i = 
          local_map.find(next_index);
        //if this index is not a processor boundary, update it
        if (local_map.find(next_index) != end && current_bc != PROCESSOR) {
          new_index = map_i->second;
          cell.set_neighbor( dir_type(d) , new_index );
        }
      } // end direction
    } // end cell
  }

  void update_mesh(void) {
    using std::vector;
    vector<Cell> new_mesh;
    for (uint32_t i =0; i< cell_list.size(); i++) {
      bool delete_flag = false;
      for (vector<uint32_t>::iterator rmv_itr= remove_cell_list.begin();
        rmv_itr != remove_cell_list.end(); 
        rmv_itr++) 
      {
        if (*rmv_itr == i)  delete_flag = true;
      }
      if (delete_flag == false) new_mesh.push_back(cell_list[i]);
    }

    for (uint32_t i =0; i< new_cell_list.size(); i++) 
      new_mesh.push_back(new_cell_list[i]);
    cell_list = new_mesh;
    n_cell = cell_list.size();
    new_cell_list.clear();
    remove_cell_list.clear();
    sort(cell_list.begin(), cell_list.end());

    //use the final number of cells to size vectors
    m_census_E = vector<double>(n_cell, 0.0);
    m_emission_E = vector<double>(n_cell, 0.0);
    m_source_E = vector<double>(n_cell, 0.0);
  }

  void make_MPI_window(void) {
    //make the MPI window with the sorted cell list
    MPI_Aint num_bytes(n_cell*mpi_cell_size);
    MPI_Alloc_mem(num_bytes, MPI_INFO_NULL, &cells);
    memcpy(cells,&cell_list[0], num_bytes);
    
    cell_list.clear();
  }

  void update_temperature(std::vector<double>& abs_E, IMC_State* imc_s) {
    //abs E is a global vector
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol,cV,rho,T, T_new;
    for (uint32_t i=0; i<n_cell;++i) {
      Cell& e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      rho = e.get_rho();
      T_new = T + (abs_E[i+on_rank_start] - m_emission_E[i])/(cV*vol*rho);
      e.set_T_e(T_new);
      total_abs_E+=abs_E[i+on_rank_start];
      total_post_mat_E+= T_new*cV*vol*rho;
    }
    //zero out absorption tallies for all cells (global) 
    for (uint32_t i=0; i<abs_E.size();++i) {
      abs_E[i] = 0.0;
    }
    imc_s->set_absorbed_E(total_abs_E);
    imc_s->set_post_mat_E(total_post_mat_E);
    imc_s->set_step_cells_requested(off_rank_reads);
    off_rank_reads = 0;
  }

  void request_cell(const uint32_t& index) {
    //get local index of global index
    uint32_t off_rank_id = get_off_rank_id(index);
    uint32_t off_rank_local = index - off_rank_bounds[off_rank_id];
    //get correct index into received photon vector
    uint32_t r_index = off_rank_id - (off_rank_id>rank);
    //add to the request list if not already requested (use global index)
    if (ids_requested.find(index) == ids_requested.end()) {
      ids_requested.insert(index);
      ids_needed[r_index].push_back(off_rank_local);
      //need data from off_rank_id that has not already been requested
      need_data[r_index] = true;
      //increment number of RMA requests
      off_rank_reads++;
    }
  }



  bool process_mesh_requests(mpi::communicator world,
                              uint32_t& n_cell_messages,
                              uint32_t& n_cells_sent,
                              uint32_t& n_sends_posted,
                              uint32_t& n_sends_completed,
                              uint32_t& n_receives_posted,
                              uint32_t& n_receives_completed) {
    using Constants::cell_tag;
    using Constants::cell_id_tag;
    using std::vector;

    bool new_data = false;
    for (uint32_t ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        uint32_t r_index = ir - (ir>rank);
        
        ////////////////////////////////////////////////////////////////////////
        // if you need data from this rank, process request
        ////////////////////////////////////////////////////////////////////////
        if (need_data[r_index]) {
  
          //test completion of send/receives
          if (recv_cell_buffer[r_index].awaiting()) {
            if (r_cell_reqs[r_index].test()) {
              n_receives_completed++;
              recv_cell_buffer[r_index].set_received();
            }
          }
          if (send_id_buffer[r_index].sent()) {
            if (s_id_reqs[r_index].test()) {
              n_sends_completed++;
              send_id_buffer[r_index].reset();
            }
          }

          // if receive buffer is empty, all receive processing has completed
          // a new send id receive cell cycle can begin
          if (recv_cell_buffer[r_index].empty() ) {
            //copy needed cells to the send buffer
            send_id_buffer[r_index].fill(ids_needed[r_index]);
            //send to ir
            s_id_reqs[r_index] = 
              world.isend(ir, cell_id_tag, send_id_buffer[r_index].get_object());
            n_sends_posted++;
            send_id_buffer[r_index].set_sent();
            // clear requested cell ids 
            // (requested cells will not be requested again)
            // because they are still stored in requested_ids set
            ids_needed[r_index].clear();
            //post receive
            r_cell_reqs[r_index] = 
              world.irecv(ir, cell_tag, recv_cell_buffer[r_index].get_object());
            n_receives_posted++;
            recv_cell_buffer[r_index].set_awaiting();
          }
          // if send id buffer has completed and cell buffer has receieved, 
          // process received data
          if (send_id_buffer[r_index].empty() && recv_cell_buffer[r_index].received()) {
            //reset need_data flag if there is no more data needed
            if (ids_needed[r_index].empty()) need_data[r_index] = false;
            new_data = true;

            //add received cells to working mesh
            vector<Cell> r_cells = recv_cell_buffer[r_index].get_object();
            for (uint32_t i=0; i<r_cells.size();i++) {
              uint32_t index = r_cells[i].get_ID();
              //add this cell to the map, if possible, otherwise manage map
              if (stored_cells.size() < max_map_size) 
                stored_cells[index] = r_cells[i];
              else {
                //remove from map and from requests so it can be 
                // reqeusted again if needed
                uint32_t removed_id = (stored_cells.begin())->first ;
                ids_requested.erase(removed_id);
                stored_cells.erase(stored_cells.begin());
                stored_cells[index] = r_cells[i];
              }
            } // for i in r_cells[r_index]
            recv_cell_buffer[r_index].reset();
          } // if recv_cell_buffer[r_index].received()
        } // if need_data[r_index]


        ////////////////////////////////////////////////////////////////////////
        // receiving cell ids needed by other ranks (post receives to all)
        ////////////////////////////////////////////////////////////////////////
        // check to see if cell ids received  and send has completed
        if ( recv_id_buffer[r_index].empty() &&
             send_cell_buffer[r_index].empty()) {
          r_id_reqs[r_index] =
            world.irecv( ir, cell_id_tag, recv_id_buffer[r_index].get_object());
          n_receives_posted++;
          recv_id_buffer[r_index].set_awaiting();
        }

        // add cell ids to requested for a rank
        if (recv_id_buffer[r_index].received()) {
          //make cell send list for this rank
          vector<uint32_t> r_ids = recv_id_buffer[r_index].get_object();
          vector<Cell> s_cell;
          for (uint32_t i=0; i<r_ids.size();i++)
            s_cell.push_back(cells[r_ids[i]]);
          send_cell_buffer[r_index].fill(s_cell);
          s_cell_reqs[r_index] = 
            world.isend( ir, cell_tag, send_cell_buffer[r_index].get_object());
          n_sends_posted++;
          n_cell_messages++;
          n_cells_sent+=s_cell.size();
          send_cell_buffer[r_index].set_sent();
          //reset receive buffer
          recv_id_buffer[r_index].reset();
        }
   
        // test receive id buffer
        if (recv_id_buffer[r_index].awaiting()) {
          if (r_id_reqs[r_index].test()) {
            n_receives_completed++;
            recv_id_buffer[r_index].set_received();
          }
        }
        //test send cell buffer
        if (send_cell_buffer[r_index].sent()) {
          if (s_cell_reqs[r_index].test()) {
            n_sends_completed++;
            send_cell_buffer[r_index].reset();
          }
        }
      } //if (ir != rank)
    } //for ir in rank
    return new_data;
  }


  void purge_working_mesh(void) {
    stored_cells.clear(); ghost_cells.clear();
    ids_requested.clear();  
  }

  void finish_mesh_pass_messages(mpi::communicator world,
                                uint32_t& n_sends_posted,
                                uint32_t& n_sends_completed,
                                uint32_t& n_receives_posted,
                                uint32_t& n_receives_completed) {
    using Constants::cell_id_tag;
    using std::vector;
    //post receives for ids to all ranks
    for (uint32_t ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        uint32_t r_index = ir - (ir>rank);
        //if receive buffer is not awaiting, post receive
        if (!recv_id_buffer[r_index].awaiting()) {
          recv_id_buffer[r_index].reset();
          r_id_reqs[r_index] = 
              world.irecv(ir, cell_id_tag, recv_id_buffer[r_index].get_object());
          recv_id_buffer[r_index].set_awaiting();
          n_receives_posted++;
        }
      }
    }
    
    //send empty cell id vector to all ranks 
    for (uint32_t ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        uint32_t r_index = ir - (ir>rank);
        //make sure sent cell id messages have completed
        if (send_cell_buffer[r_index].sent()) {
          s_cell_reqs[r_index].wait();
          n_sends_completed++;
        }
        send_cell_buffer[r_index].reset();
        //make sure sent id messages have completed
        if (send_id_buffer[r_index].sent()) {
          s_id_reqs[r_index].wait();
          n_sends_completed++;
        }
        //send empty id message to finish awaiting receives
        vector<uint32_t> empty_id_vector;
        send_id_buffer[r_index].fill(empty_id_vector);
        s_id_reqs[r_index] = 
          world.isend(ir, cell_id_tag, send_id_buffer[r_index].get_object());
        n_sends_posted++;
        //wait for message to send
        s_id_reqs[r_index].wait();
        n_sends_completed++;
        //reset send if buffer
        send_id_buffer[r_index].reset();
      }
    }

    //wait for id receives to complete 
    for (uint32_t ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        uint32_t r_index = ir - (ir>rank);
        //wait for receives to complete
        r_id_reqs[r_index].wait();
        n_receives_completed++;
        recv_id_buffer[r_index].reset();
      }
    }
  }

  void add_mesh_cell(Cell new_cell) {new_cell_list.push_back(new_cell);}
  void remove_cell(uint32_t index) {remove_cell_list.push_back(index);}

  std::vector<double>& get_census_E_ref(void) {return m_census_E;}
  std::vector<double>& get_emission_E_ref(void) {return m_emission_E;}
  std::vector<double>& get_source_E_ref(void) {return m_source_E;}
  std::vector<double> get_census_E(void) const {return m_census_E;}
  std::vector<double> get_emission_E(void) const {return m_emission_E;}
  std::vector<double> get_source_E(void) const {return m_source_E;}

/*****************************************************************************/
  //member variables
/*****************************************************************************/
  private:

  uint32_t ngx; //! Number of global x sizes
  uint32_t ngy; //! Number of global y sizes
  uint32_t ngz; //! Number of global z sizes
  uint32_t rank; //! MPI rank of this mesh
  uint32_t n_rank; //! Number of global ranks

  uint32_t n_cell; //! Number of local cells
  uint32_t n_global; //! Nuber of global cells
  
  uint32_t on_rank_start; //! Start of global index on rank
  uint32_t on_rank_end; //! End of global index on rank

  std::vector<double> m_census_E; //! Census energy vector
  std::vector<double> m_emission_E; //! Emission energy vector
  std::vector<double> m_source_E; //! Source energy vector

  Cell *cells; //! Cell data allocated with MPI_Alloc
  std::vector<Cell> cell_list; //! On processor cells
  std::vector<Cell> new_cell_list; //! New received cells
  std::vector<uint32_t> remove_cell_list; //! Cells to be removed
  std::vector<uint32_t> off_rank_bounds; //! Ending value of global ID for each rank
  std::vector<uint32_t> boundary_cells; //! Index of adjacent ghost cells
 
  std::map<uint32_t, uint32_t> adjacent_procs; //! List of adjacent processors

  double m_opA; //! Opacity coefficient A in A + B^C
  double m_opB; //! Opacity coefficient B in A + B^C
  double m_opC; //! Opacity coefficient C in A + B^C
  double m_opS; //! Scattering opacity constant

  double total_photon_E; //! Total photon energy on the mesh

  MPI_Datatype MPI_Cell; //! MPI type, allows simpler parallel communication
  int32_t mpi_cell_size; //! Size of custom MPI_Cell type

  uint32_t max_map_size; //! Maximum size of map object
  uint32_t off_rank_reads; //! Number of off rank reads

  //send and receive buffers
  std::vector<Buffer<Cell> > recv_cell_buffer; //! Receive buffer for cells
  std::vector<Buffer<Cell> > send_cell_buffer; //! Send buffer for cells
  std::vector<Buffer<uint32_t> > recv_id_buffer; //! Receive buffer for cell IDs
  std::vector<Buffer<uint32_t> > send_id_buffer; //! Receive buffer for cell IDs

  //Data bools needed from off rank 
  std::vector<bool> need_data; //! Vector of size nrank-1, flag for needed data from rank

  // MPI requests for non-blocking cell communication
  mpi::request* r_cell_reqs; //! Received cell requests
  mpi::request* s_cell_reqs; //! Send cell requests
  mpi::request* r_id_reqs; //! Received cell ID requests
  mpi::request* s_id_reqs; //! Send cell ID requests

  std::map<uint32_t, Cell> stored_cells; //! Cells that have been accessed off rank
  std::map<uint32_t, Cell> ghost_cells; //! Static list of off-rank cells next to boundary

  std::vector<std::vector<uint32_t> > ids_needed ; //! Cell needed by this rank
  std::set<uint32_t> ids_requested; //! IDs that have been requested

};

#endif //mesh_h_
