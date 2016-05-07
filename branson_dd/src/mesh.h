
#ifndef mesh_h_
#define mesh_h_

#include <algorithm>
#include <map>
#include <mpi.h>
#include <set>
#include <string>
#include <vector>

#include "buffer.h"
#include "cell.h"
#include "constants.h"
#include "input.h"
#include "imc_state.h"

class Mesh {

  public:

  Mesh(Input* input, int _rank, int _n_rank)
  : ngx(input->get_global_n_x_cells()),
    ngy(input->get_global_n_y_cells()),
    ngz(input->get_global_n_z_cells()),
    rank(_rank),
    n_rank(_n_rank),
    n_off_rank(n_rank-1),
    silo_x(input->get_silo_x_ptr()),
    silo_y(input->get_silo_y_ptr()),
    silo_z(input->get_silo_z_ptr())
  {
    using std::vector;
    using Constants::bc_type;
    using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
    using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
    using Constants::ELEMENT;

    max_map_size = input->get_map_size();
    double dx,dy,dz;

    // make off processor map
    for (uint32_t i=0; i<n_off_rank; i++) {
      int r_index = i + int(i>=rank);
      proc_map[i] = r_index;
    }

    regions = input->get_regions();
    // map region IDs to index in the region
    for (uint32_t i =0; i<regions.size(); i++) {
      region_ID_to_index[regions[i].get_ID()] = i;
    }

    vector<bc_type> bc(6);
    bc[X_POS] = input->get_bc(X_POS); bc[X_NEG] = input->get_bc(X_NEG);
    bc[Y_POS] = input->get_bc(Y_POS); bc[Y_NEG] = input->get_bc(Y_NEG);
    bc[Z_POS] = input->get_bc(Z_POS); bc[Z_NEG] = input->get_bc(Z_NEG);

    //initialize number of RMA requests as zero
    off_rank_reads=0;

    uint32_t global_count =0; //global cell count

    //this rank's cells
    n_global = ngx*ngy*ngz;
    uint32_t cell_id_begin = floor(rank*n_global/double(n_rank));
    uint32_t cell_id_end = floor((rank+1)*n_global/double(n_rank));

    uint32_t on_rank_count =0;

    uint32_t n_x_div = input->get_n_x_divisions();
    uint32_t n_y_div = input->get_n_y_divisions();
    uint32_t n_z_div = input->get_n_z_divisions();

    Region region;
    uint32_t region_index, nx, ny, nz;
    double x_start, y_start, z_start;
    double x_cell_end;
    double y_cell_end;
    double z_cell_end;

    for (uint32_t iz_div=0; iz_div<n_z_div; iz_div++) {
      dz = input->get_dz(iz_div);
      nz = input->get_z_division_cells(iz_div);
      z_start = input->get_z_start(iz_div);
      for (uint32_t k=0; k<nz; k++) {
        for (uint32_t iy_div=0; iy_div<n_y_div; iy_div++) {
          dy = input->get_dy(iy_div);
          ny = input->get_y_division_cells(iy_div);
          y_start = input->get_y_start(iy_div);
          for (uint32_t j=0; j<ny; j++) {
            for (uint32_t ix_div=0; ix_div<n_x_div; ix_div++) {
              dx = input->get_dx(ix_div);
              nx = input->get_x_division_cells(ix_div);
              x_start = input->get_x_start(ix_div);
              for (uint32_t i=0; i<nx; i++) {
                if (global_count >= cell_id_begin && global_count < cell_id_end) {
                  //find the region for this cell
                  region_index = input->get_region_index(ix_div, iy_div, iz_div);
                  region = regions[region_index];
                  Cell e;
                  // set ending coordinates explicity to match the start of
                  // the next division to avoid werid roundoff errors

                  if (i == nx-1 && ix_div != n_x_div-1) 
                    x_cell_end = input->get_x_start(ix_div+1);
                  else x_cell_end = x_start+(i+1)*dx;

                  if (i == nx-1 && ix_div != n_x_div-1) 
                    x_cell_end = input->get_x_start(ix_div+1);
                  else x_cell_end = x_start+(i+1)*dx;

                  e.set_coor(x_start+i*dx, x_start+(i+1)*dx, 
                             y_start+j*dy, y_start+(j+1)*dy, 
                             z_start+k*dz, z_start+(k+1)*dz);
                  e.set_ID(global_count);
                  e.set_region_ID(region.get_ID());

                  // set cell physical properties using region
                  e.set_cV(region.get_cV());
                  e.set_T_e(region.get_T_e()); 
                  e.set_T_r(region.get_T_r());
                  e.set_T_s(region.get_T_s());
                  e.set_rho(region.get_rho());

                  //set the global index for SILO plotting
                  e.set_silo_index(i + j*ngx + k*(ngy*ngx));

                  // set neighbors in x direction
                  if (i<(ngx-1)) {
                    e.set_neighbor( X_POS, global_count+1); 
                    e.set_bc(X_POS, ELEMENT);
                  }
                  else {
                    e.set_neighbor( X_POS, global_count); 
                    e.set_bc(X_POS, bc[X_POS]);
                  } 
                  if (i>0) {
                    e.set_neighbor( X_NEG, global_count-1); 
                    e.set_bc(X_NEG, ELEMENT);
                  }
                  else {
                    e.set_neighbor( X_NEG, global_count); 
                    e.set_bc(X_NEG, bc[X_NEG]);
                  }

                  // set neighbors in y direction
                  if (j<(ngy-1)) {
                    e.set_neighbor(Y_POS, global_count+ngx); 
                    e.set_bc(Y_POS, ELEMENT);
                  }
                  else {
                    e.set_neighbor(Y_POS, global_count); 
                    e.set_bc(Y_POS, bc[Y_POS]);
                  }
                  if (j>0) {
                    e.set_neighbor(Y_NEG, global_count-ngx);
                    e.set_bc(Y_NEG, ELEMENT);
                  }
                  else {
                    e.set_neighbor(Y_NEG, global_count);
                    e.set_bc(Y_NEG, bc[Y_NEG]);
                  }

                  // set neighbors in z direction
                  if (k<(ngz-1)) {
                    e.set_neighbor(Z_POS, global_count+ngx*ngy); 
                    e.set_bc(Z_POS, ELEMENT);
                  }
                  else {
                    e.set_neighbor(Z_POS, global_count); 
                    e.set_bc(Z_POS, bc[Z_POS]);
                  }
                  if (k>0) {
                    e.set_neighbor(Z_NEG, global_count-ngx*ngy); 
                    e.set_bc(Z_NEG, ELEMENT);
                  }
                  else {
                    e.set_neighbor(Z_NEG, global_count);
                    e.set_bc(Z_NEG, bc[Z_NEG]);
                  }

                  // add cell to mesh
                  cell_list.push_back(e);
                  //increment on rank count
                  on_rank_count++;
                } // end if on processor check
                global_count++;
              } // end i loop
            } // end x division loop
          } // end j loop
        } // end y division loop
      } // end k loop
    } // end z division loop
    n_cell = on_rank_count;

    total_photon_E = 0.0;

    //define the MPI cell type used for MPI windows
    //make the MPI datatype for my cell class
    // Three type entries in the class
    const int entry_count = 3 ; 
    // 7 uint32_t, 6 int, 13 double
    int array_of_block_length[4] = {10, 6, 14};
    // Displacements of each type in the cell
    MPI_Aint array_of_block_displace[3] = 
      {0, 10*sizeof(uint32_t),  10*sizeof(uint32_t)+6*sizeof(int)};
    //Type of each memory block
    MPI_Datatype array_of_types[3] = {MPI_UNSIGNED, MPI_INT, MPI_DOUBLE}; 

    MPI_Type_create_struct(entry_count, array_of_block_length, 
      array_of_block_displace, array_of_types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);

    MPI_Type_size(MPI_Cell, &mpi_cell_size);

    uint32_t cell_size = sizeof(Cell); 

    //bool flags to say if data is needed
    need_data = vector<bool>(n_rank-1, false);

    //allocate mpi requests objects
    r_cell_reqs = new MPI_Request[n_off_rank];
    r_id_reqs = new MPI_Request[n_off_rank];
    s_cell_reqs = new MPI_Request[n_off_rank];
    s_id_reqs = new MPI_Request[n_off_rank];
   
    recv_cell_buffer = vector<Buffer<Cell> >(n_off_rank);
    send_cell_buffer = vector<Buffer<Cell> >(n_off_rank);
    recv_id_buffer =  vector<Buffer<uint32_t> >(n_off_rank);
    send_id_buffer =  vector<Buffer<uint32_t> >(n_off_rank);

    n_send_ids = vector<int32_t>(n_off_rank);
    // size the send and receive buffers
    // they are size n_off_rank
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
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

  std::vector<uint32_t> get_off_rank_bounds(void) {return off_rank_bounds;}

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

  const Cell* get_cell_ptr(const uint32_t& local_ID) const {
    return &cells[local_ID];
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

  uint32_t get_global_ID(const uint32_t& local_index) const {
    return on_rank_start+local_index;
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

  std::vector<double> get_census_E(void) const {return m_census_E;}
  std::vector<double> get_emission_E(void) const {return m_emission_E;}
  std::vector<double> get_source_E(void) const {return m_source_E;}

  uint32_t get_global_n_x_faces(void) const {return ngx+1;}
  uint32_t get_global_n_y_faces(void) const {return ngy+1;}
  uint32_t get_global_n_z_faces(void) const {return ngz+1;}

  float * get_silo_x(void) const {return silo_x;}
  float * get_silo_y(void) const {return silo_y;}
  float * get_silo_z(void) const {return silo_z;}

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

    uint32_t region_ID;
    Region region;
    for (uint32_t i=0; i<n_cell;++i) {
      Cell& e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      Tr = e.get_T_r();
      Ts = e.get_T_s();
      rho = e.get_rho();

      region_ID = e.get_region_ID();
      region =  regions[region_ID_to_index[region_ID]];
          
      op_a = region.get_absorption_opacity(T);
      op_s = region.get_scattering_opacity();
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
    MPI_Aint n_bytes(n_cell*mpi_cell_size);
    //MPI_Alloc_mem(n_bytes, MPI_INFO_NULL, &cells);
    MPI_Win_allocate(n_bytes, mpi_cell_size, MPI_INFO_NULL,
      MPI_COMM_WORLD, &cells, &mesh_window);
    //copy the cells list data into the cells array
    memcpy(cells,&cell_list[0], n_bytes);

    cell_list.clear();
  }

  void update_temperature(std::vector<double>& abs_E, IMC_State* imc_s) {
    //abs E is a global vector
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol,cV,rho,T, T_new;
    uint32_t region_ID;
    Region region;
    for (uint32_t i=0; i<n_cell;++i) {
      region_ID = cells[i].get_region_ID();
      region = regions[region_ID_to_index[ region_ID]];
      cV = region.get_cV();
      rho = region.get_rho();
      Cell& e = cells[i];
      vol = e.get_volume();
      T = e.get_T_e();
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

  MPI_Win& get_mesh_window_ref(void) {return mesh_window;}

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

  void add_non_local_mesh_cells(std::vector<Cell> new_recv_cells) {
    for (uint32_t i=0; i<new_recv_cells.size();i++) {
      uint32_t index = new_recv_cells[i].get_ID();

      //add this cell to the map, if possible, otherwise manage map
      if (stored_cells.size() < max_map_size) 
        stored_cells[index] = new_recv_cells[i];
      else {
        //remove first cell from map, it will naturally be requested again
        stored_cells.erase(stored_cells.begin());
      }
    } // for i in new_recv_cells[ir] 
  }

  bool process_mesh_requests(uint32_t& n_cell_messages,
                              uint32_t& n_cells_sent,
                              uint32_t& n_sends_posted,
                              uint32_t& n_sends_completed,
                              uint32_t& n_receives_posted,
                              uint32_t& n_receives_completed) 
  {
    using Constants::cell_tag;
    using Constants::cell_id_tag;
    using std::vector;

    const uint32_t n_max_id_recv(5000);
    bool new_data = false;
    int off_rank;
    uint32_t n_req_cells;
    uint32_t n_send_cells;
    MPI_Status mpi_recv_ids_status;
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      off_rank = proc_map[ir]; 
      ////////////////////////////////////////////////////////////////////////
      // if you need data from this rank, process request
      ////////////////////////////////////////////////////////////////////////
      if (need_data[ir]) {

        //test completion of send/receives
        if (recv_cell_buffer[ir].awaiting()) {
          MPI_Test(&r_cell_reqs[ir], &r_cell_req_flag, MPI_STATUS_IGNORE);
          if (r_cell_req_flag) {
            n_receives_completed++;
            recv_cell_buffer[ir].set_received();
          }
        }
        if (send_id_buffer[ir].sent()) {
          MPI_Test(&s_id_reqs[ir], &s_id_req_flag, MPI_STATUS_IGNORE);
          if (s_id_req_flag) {
            n_sends_completed++;
            send_id_buffer[ir].reset();
          }
        }

        // if receive buffer is empty, all receive processing has completed
        // a new send id receive cell cycle can begin
        if (recv_cell_buffer[ir].empty() ) {
          //copy needed cells to the send buffer
          send_id_buffer[ir].fill(ids_needed[ir]);
          // get the number of cells for proper receiving
          n_req_cells = ids_needed[ir].size();
          //send to off_rank
          MPI_Isend(send_id_buffer[ir].get_buffer(), n_req_cells, MPI_UNSIGNED,
            off_rank, cell_id_tag, MPI_COMM_WORLD, &s_id_reqs[ir]);
          n_sends_posted++;
          send_id_buffer[ir].set_sent();
          // clear requested cell ids 
          // (requested cells will not be requested again)
          // because they are still stored in requested_ids set
          ids_needed[ir].clear();
          //post receive
          recv_cell_buffer[ir].resize(n_req_cells);
          MPI_Irecv(recv_cell_buffer[ir].get_buffer(), n_req_cells, MPI_Cell,
            off_rank, cell_tag, MPI_COMM_WORLD, &r_cell_reqs[ir]);
          n_receives_posted++;
          recv_cell_buffer[ir].set_awaiting();
        }

        // if send id buffer has completed and cell buffer has receieved, 
        // process received data
        if (send_id_buffer[ir].empty() && recv_cell_buffer[ir].received()) {
          //reset need_data flag if there is no more data needed
          if (ids_needed[ir].empty()) need_data[ir] = false;
          new_data = true;

          //add received cells to working mesh
          vector<Cell> r_cells = recv_cell_buffer[ir].get_object();
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
          } // for i in r_cells[ir]
          recv_cell_buffer[ir].reset();
        } // if recv_cell_buffer[ir].received()
      } // if need_data[ir]


      ////////////////////////////////////////////////////////////////////////
      // receiving cell ids needed by other ranks (post receives to all)
      ////////////////////////////////////////////////////////////////////////
      // check to see if cell ids received  and send has completed
      if ( recv_id_buffer[ir].empty() &&
           send_cell_buffer[ir].empty()) {
        recv_id_buffer[ir].resize(n_max_id_recv);
        MPI_Irecv(recv_id_buffer[ir].get_buffer(), n_max_id_recv, MPI_UNSIGNED,
          off_rank, cell_id_tag, MPI_COMM_WORLD, &r_id_reqs[ir]);
        n_receives_posted++;
        recv_id_buffer[ir].set_awaiting();
      }

      // add cell ids to requested for a rank
      if (recv_id_buffer[ir].received()) {
        n_send_cells = n_send_ids[ir];
        //make cell send list for this rank
        vector<uint32_t> r_ids = recv_id_buffer[ir].get_object();
        vector<Cell> s_cell;
        for (uint32_t i=0; i<n_send_cells;i++)
          s_cell.push_back(cells[r_ids[i]]);
        send_cell_buffer[ir].fill(s_cell);
        MPI_Isend(send_cell_buffer[ir].get_buffer(), n_send_cells, MPI_Cell,
          off_rank, cell_tag, MPI_COMM_WORLD, &s_cell_reqs[ir]);
        n_sends_posted++;
        n_cell_messages++;
        n_cells_sent+=s_cell.size();
        send_cell_buffer[ir].set_sent();
        //reset receive buffer
        recv_id_buffer[ir].reset();
      }
 
      // test receive id buffer
      if (recv_id_buffer[ir].awaiting()) {
        MPI_Test(&r_id_reqs[ir], &r_id_req_flag, &mpi_recv_ids_status);
        if (r_id_req_flag) {
          MPI_Get_count(&mpi_recv_ids_status, MPI_UNSIGNED, &n_send_ids[ir]);
          n_receives_completed++;
          recv_id_buffer[ir].set_received();
        }
      }
      //test send cell buffer
      if (send_cell_buffer[ir].sent()) {
        MPI_Test(&s_cell_reqs[ir], &s_cell_req_flag, MPI_STATUS_IGNORE);
        if (s_cell_req_flag) {
          n_sends_completed++;
          send_cell_buffer[ir].reset();
        }
      }
    } //for ir in n_off_rank
    return new_data;
  }

  void purge_working_mesh(void) {
    stored_cells.clear(); ghost_cells.clear();
    ids_requested.clear();  
  }

  void finish_mesh_pass_messages(uint32_t& n_sends_posted,
                                uint32_t& n_sends_completed,
                                uint32_t& n_receives_posted,
                                uint32_t& n_receives_completed) 
  {
    using Constants::cell_id_tag;
    using std::vector;
    //post receives for ids to all ranks
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      //get correct index into requests and vectors 
      int off_rank = proc_map[ir];
      //if receive buffer is not awaiting, post receive
      if (!recv_id_buffer[ir].awaiting()) {
        recv_id_buffer[ir].reset();
        recv_id_buffer[ir].resize(1);
        MPI_Irecv(recv_id_buffer[ir].get_buffer(), 1, MPI_UNSIGNED,
          off_rank, cell_id_tag, MPI_COMM_WORLD, &r_id_reqs[ir]);
        recv_id_buffer[ir].set_awaiting();
        n_receives_posted++;
      }
    }
    
    //send empty cell id vector to all ranks 
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      //get correct index into requests and vectors 
      int off_rank = proc_map[ir];
      //make sure sent cell id messages have completed
      if (send_cell_buffer[ir].sent()) {
        MPI_Wait(&s_cell_reqs[ir], MPI_STATUS_IGNORE);
        n_sends_completed++;
      }
      send_cell_buffer[ir].reset();
      //make sure sent id messages have completed
      if (send_id_buffer[ir].sent()) {
        MPI_Wait(&s_id_reqs[ir], MPI_STATUS_IGNORE);
        n_sends_completed++;
      }
      //send empty id message to finish awaiting receives
      vector<uint32_t> empty_id_vector;
      send_id_buffer[ir].fill(empty_id_vector);
      MPI_Isend(send_id_buffer[ir].get_buffer(), 1, MPI_UNSIGNED,
        off_rank, cell_id_tag, MPI_COMM_WORLD, &s_id_reqs[ir]);
      n_sends_posted++;
      //wait for message to send
      MPI_Wait(&s_id_reqs[ir], MPI_STATUS_IGNORE);
      n_sends_completed++;
      //reset send if buffer
      send_id_buffer[ir].reset();
    }

    //wait for id receives to complete 
    for (uint32_t ir=0; ir<n_off_rank; ir++) {
      //get correct index into requests and vectors 
      int off_rank = proc_map[ir];
      //make sure sent cell id messages have completed
      //wait for receives to complete
      MPI_Wait(&r_id_reqs[ir], MPI_STATUS_IGNORE);
      n_receives_completed++;
      recv_id_buffer[ir].reset();
    }
  }

  void add_mesh_cell(Cell new_cell) {new_cell_list.push_back(new_cell);}
  void remove_cell(uint32_t index) {remove_cell_list.push_back(index);}

  std::vector<double>& get_census_E_ref(void) {return m_census_E;}
  std::vector<double>& get_emission_E_ref(void) {return m_emission_E;}
  std::vector<double>& get_source_E_ref(void) {return m_source_E;}

  //////////////////////////////////////////////////////////////////////////////
  //member variables
  //////////////////////////////////////////////////////////////////////////////
  private:

  uint32_t ngx; //! Number of global x sizes
  uint32_t ngy; //! Number of global y sizes
  uint32_t ngz; //! Number of global z sizes
  uint32_t rank; //! MPI rank of this mesh
  uint32_t n_rank; //! Number of global ranks
  uint32_t n_off_rank; //! Number of other ranks
  float *silo_x; //! Global array of x face locations for SILO
  float *silo_y; //! Global array of y face locations for SILO
  float *silo_z; //! Global array of z face locations for SILO

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

  std::vector<Region> regions; //! Vector of regions in the problem
  std::map<uint32_t, uint32_t> region_ID_to_index; //! Maps region ID to index

  double total_photon_E; //! Total photon energy on the mesh

  MPI_Datatype MPI_Cell; //! MPI type, allows simpler parallel communication
  int32_t mpi_cell_size; //! Size of custom MPI_Cell type
  int32_t mpi_particle_size; //! Size of custom MPI_Particle type

  uint32_t max_map_size; //! Maximum size of map object
  uint32_t off_rank_reads; //! Number of off rank reads

  MPI_Win mesh_window;

  //send and receive buffers
  std::vector<Buffer<Cell> > recv_cell_buffer; //! Receive buffer for cells
  std::vector<Buffer<Cell> > send_cell_buffer; //! Send buffer for cells
  std::vector<Buffer<uint32_t> > recv_id_buffer; //! Receive buffer for cell IDs
  std::vector<Buffer<uint32_t> > send_id_buffer; //! Receive buffer for cell IDs

  //Data bools needed from off rank 
  std::vector<bool> need_data; //! Vector of size nrank-1, flag for needed data from rank

  // MPI requests for non-blocking cell communication
  MPI_Request *r_cell_reqs; //! Received cell requests
  MPI_Request *s_cell_reqs; //! Send cell requests
  MPI_Request *r_id_reqs; //! Received cell ID requests
  MPI_Request *s_id_reqs; //! Send cell ID requests

  //flags for MPI request testing
  int r_cell_req_flag; //! Flag for received cell communications
  int s_cell_req_flag; //! Flag for sent cell communications
  int r_id_req_flag;  //! Flag for received cell indices
  int s_id_req_flag; //! Flag for send cell indices

  std::map<uint32_t, Cell> stored_cells; //! Cells that have been accessed off rank
  std::map<uint32_t, Cell> ghost_cells; //! Static list of off-rank cells next to boundary

  std::map<int,int> proc_map; //! Maps number of off-rank processor to global rank

  std::vector<std::vector<uint32_t> > ids_needed; //! Cell needed by this rank
  std::vector<int32_t > n_send_ids; //! Number of IDs to send
  std::set<uint32_t> ids_requested; //! IDs that have been requested

};

#endif //mesh_h_
