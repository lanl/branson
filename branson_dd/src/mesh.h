
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
#include "request.h"

namespace mpi = boost::mpi;

class Mesh {

  public:

  Mesh(Input* input, int _rank, int _nrank)
  : ngx(input->get_n_x_cells()),
    ngy(input->get_n_y_cells()),
    ngz(input->get_n_z_cells()),
    rank(_rank),
    n_rank(_nrank)
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

    //number of ranks
    n_rank =MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();

    //initialize number of DMA requests as zero
    off_rank_reads=0;

    unsigned int g_count =0; //global count

    //this rank's cells
    n_global = ngx*ngy*ngz;
    unsigned int cell_id_begin = floor(rank*n_global/double(n_rank));
    unsigned int cell_id_end = floor((rank+1)*n_global/double(n_rank));

    unsigned int l_count =0;
    for (unsigned int k=0; k<ngz; k++) {
      for (unsigned int j=0; j<ngy; j++) {
        for (unsigned int i=0; i<ngx; i++) {
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
    // 7 unsigned int, 6 int, 13 double
    const int array_of_block_length[4] = {8, 6, 14};
    // Displacements of each type in the cell
    const MPI::Aint array_of_block_displace[3] = 
      {0, 8*sizeof(unsigned int),  8*sizeof(unsigned int)+6*sizeof(int)};
    //Type of each memory block
    MPI::Datatype array_of_types[3] = {MPI_UNSIGNED, MPI_INT, MPI_DOUBLE}; 

    MPI_Cell = MPI::Datatype::Create_struct(entry_count,
      array_of_block_length, array_of_block_displace, array_of_types);
    MPI_Cell.Commit();

    //bool flags to say if data is needed
    need_data = vector<bool>(n_rank-1, false);
    send_data = vector<bool>(n_rank-1, false);

    //allocate mpi requests objects
    r_cell_reqs = new Request[ (n_rank-1)];
    r_cell_ids_reqs = new Request[ (n_rank-1)];
    s_cell_reqs = new Request[ (n_rank-1)];
    s_cell_ids_reqs = new Request[ (n_rank-1)];
   
    b_r_cell_reqs    = vector<bool> (n_rank-1, false);
    b_r_cell_ids_reqs = vector<bool> (n_rank-1, false);
    b_s_cell_reqs    = vector<bool> (n_rank-1, false);
    b_s_cell_ids_reqs = vector<bool> (n_rank-1, false);

    // size the send and receive buffers
    // they are size n_rank -1
    for (unsigned int ir=0; ir<n_rank-1; ir++) {
      vector<Cell> empty_vec_cell;
      r_cells.push_back(empty_vec_cell);
      s_cells.push_back(empty_vec_cell);
      vector<unsigned int> empty_vec_ids;
      r_cell_ids.push_back(empty_vec_ids);
      s_cell_ids.push_back(empty_vec_ids);
      //ids needed from other ranks
      ids_needed.push_back(empty_vec_ids); 
    }
    

  }

  //free buffers and delete MPI allocated cell
  ~Mesh() { }

/*****************************************************************************/
  //const functions
/*****************************************************************************/
  unsigned int get_number_of_objects(void) const {return n_cell;}
  unsigned int get_global_ID(unsigned int index) const 
  {
    return  cell_list[index].get_ID();
  }
  unsigned int get_rank(void) const {return  rank;}
  unsigned int get_offset(void) const {return on_rank_start;}
  unsigned int get_global_num_cells(void) const {return n_global;}
  double get_total_photon_E(void) const {return total_photon_E;}

  void print(void) {
    for (unsigned int i= 0; i<n_cell; i++)
      cell_list[i].print();
  }

  std::map<unsigned int, unsigned int> get_map(void) const {
    std::map<unsigned int, unsigned int> local_map;
    unsigned int g_ID;
    for (unsigned int i=0; i<n_cell; i++) {
      g_ID = cell_list[i].get_ID();
      local_map[g_ID] = i+on_rank_start;
    }
    return local_map;
  }

  Cell get_pre_cell(const unsigned int& local_ID) const 
  {
    return cell_list[local_ID];
  } 
  Cell get_cell(const unsigned int& local_ID) const {
    return cells[local_ID];
  } 

  unsigned int get_off_rank_id(const unsigned int& index) const {
    //find rank of index
    bool found = false;
    unsigned int min_i = 0;
    unsigned int max_i = off_rank_bounds.size()-1;
    unsigned int s_i; //search index
    while(!found) {
      s_i =(max_i + min_i)/2;
      if (s_i == max_i || s_i == min_i) found = true;
      else if (index >= off_rank_bounds[s_i]) min_i = s_i;
      else max_i = s_i;
    }
    return s_i;
  }

  unsigned int get_rank(const unsigned int& index) const {
    unsigned int r_rank;
    if (on_processor(index)) r_rank = rank;
    else  r_rank = get_off_rank_id(index);
    return r_rank;
  }

  Cell get_on_rank_cell(const unsigned int& index) {
    //this can only be called with valid on rank indexes
    if (on_processor(index)) 
      return cells[index-on_rank_start];
    else 
      return stored_cells[index];
  }

/*****************************************************************************/
  //non-const functions
/*****************************************************************************/
  void set_global_bound(unsigned int _on_rank_start, 
                        unsigned int _on_rank_end) 
  {
    on_rank_start = _on_rank_start;
    on_rank_end = _on_rank_end;
  }

  void set_off_rank_bounds(std::vector<unsigned int> _off_rank_bounds) {
    off_rank_bounds=_off_rank_bounds;
  }

  void calculate_photon_energy(IMC_State* imc_s) {
    using Constants::c;
    using Constants::a;
    total_photon_E = 0.0;
    double dt = imc_s->get_dt();
    double op_a, op_s, f, cV;
    double vol;
    double T, Tr, Ts;
    unsigned int step = imc_s->get_step();
    double tot_census_E = 0.0;
    double tot_emission_E = 0.0;
    double tot_source_E = 0.0;
    double pre_mat_E = 0.0;
    for (unsigned int i=0; i<n_cell;++i) {
      Cell& e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      Tr = e.get_T_r();
      Ts = e.get_T_s();
      op_a = m_opA + m_opB*pow(T, m_opC);
      op_s = m_opS;
      f =1.0/(1.0 + dt*op_a*c*(4.0*a*pow(T,3)/cV));

      e.set_op_a(op_a);
      e.set_op_s(op_s);
      e.set_f(f);

      m_emission_E[i] = dt*vol*f*op_a*a*c*pow(T,4);
      if (step > 1) m_census_E[i] = 0.0;  
      else m_census_E[i] =vol*a*pow(Tr,4); 
      m_source_E[i] = dt*op_a*a*c*pow(Ts,4);

      pre_mat_E+=T*cV*vol;
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


  void set_indices( std::map<unsigned int, unsigned int> off_map, 
                    std::vector< std::vector<bool> >& remap_flag) {

    using Constants::PROCESSOR;
    using Constants::dir_type;
    using std::map;

    unsigned int next_index;
    map<unsigned int, unsigned int>::iterator end = off_map.end();
    unsigned int new_index;
    //check to see if neighbors are on or off processor
    for (unsigned int i=0; i<n_cell; i++) {
      Cell& cell = cell_list[i];
      for (unsigned int d=0; d<6; d++) {
        next_index = cell.get_next_cell(d);
        map<unsigned int, unsigned int>::iterator map_i = 
          off_map.find(next_index);
        if (off_map.find(next_index) != end && remap_flag[i][d] ==false ) {
          //update index and bc type, this will always be an off processor so
          //if an index is updated it will always be at a processor bound
          remap_flag[i][d] = true;
          new_index = map_i->second;
          cell.set_neighbor( dir_type(d) , new_index );
          cell.set_bc(dir_type(d), PROCESSOR);
          boundary_cells.push_back(new_index);
        }
      }
    }
  }

  void set_local_indices(std::map<unsigned int, unsigned int> local_map) {

    using Constants::PROCESSOR;
    using Constants::bc_type;
    using Constants::dir_type;
    using std::map;

    unsigned int next_index;
    map<unsigned int, unsigned int>::iterator end = local_map.end();
    unsigned int new_index;
    bc_type current_bc;
    //check to see if neighbors are on or off processor
    for (unsigned int i=0; i<n_cell; i++) {
      Cell& cell = cell_list[i];
      cell.set_ID(i+on_rank_start);
      for (unsigned int d=0; d<6; d++) {
        current_bc = cell.get_bc(bc_type(d));
        next_index = cell.get_next_cell(d);
        map<unsigned int, unsigned int>::iterator map_i = 
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
    for (unsigned int i =0; i< cell_list.size(); i++) {
      bool delete_flag = false;
      for (vector<unsigned int>::iterator rmv_itr= remove_cell_list.begin();
        rmv_itr != remove_cell_list.end(); 
        rmv_itr++) 
      {
        if (*rmv_itr == i)  delete_flag = true;
      }
      if (delete_flag == false) new_mesh.push_back(cell_list[i]);
    }

    for (unsigned int i =0; i< new_cell_list.size(); i++) 
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
    unsigned int num_bytes =n_cell*MPI_Cell.Get_size();
    cells = (Cell*) MPI::Alloc_mem(num_bytes, MPI_INFO_NULL);
    memcpy(cells,&cell_list[0], num_bytes);
    
    cell_list.clear();
  }

  void update_temperature(std::vector<double>& abs_E, IMC_State* imc_s) {
    //abs E is a global vector
    double total_abs_E = 0.0;
    double total_post_mat_E = 0.0;
    double vol,cV,rho,T, T_new;
    for (unsigned int i=0; i<n_cell;++i) {
      Cell& e = cells[i];
      vol = e.get_volume();
      cV = e.get_cV();
      T = e.get_T_e();
      rho = e.get_rho();
      T_new = T + (abs_E[i+on_rank_start] - m_emission_E[i])/(cV*vol*rho);
      e.set_T_e(T_new);
      total_abs_E+=abs_E[i+on_rank_start];
      total_post_mat_E+= T_new*cV*vol;
    }
    //zero out absorption tallies for all cells (global) 
    for (unsigned int i=0; i<abs_E.size();++i) {
      abs_E[i] = 0.0;
    }
    imc_s->set_absorbed_E(total_abs_E);
    imc_s->set_post_mat_E(total_post_mat_E);
    imc_s->set_off_rank_read(off_rank_reads);
    off_rank_reads = 0;
  }


  /*
  void clear_messages(void) {
    for (unsigned int ir=0; ir<n_rank; ir++) {
      if (ir != rank) {
        //get correct index into requests and vectors 
        unsigned int r_index = ir - (ir>rank);
        b_r_cell_ids_reqs[r_index] = false;
        b_s_cell_ids_reqs[r_index] = false;
        b_r_cell_reqs[r_index] = false;
        b_s_cell_reqs[r_index] = false;
        r_cell_ids_reqs[r_index].cancel();
        s_cell_ids_reqs[r_index].cancel();
        r_cell_reqs[r_index].cancel();
        s_cell_reqs[r_index].cancel();
      }
    }
  }
  */

  void add_mesh_cell(Cell new_cell) {new_cell_list.push_back(new_cell);}
  void remove_cell(unsigned int index) {remove_cell_list.push_back(index);}

  std::vector<double>& get_census_E_ref(void) {return m_census_E;}
  std::vector<double>& get_emission_E_ref(void) {return m_emission_E;}
  std::vector<double>& get_source_E_ref(void) {return m_source_E;}

/*****************************************************************************/
  //member variables
/*****************************************************************************/
  private:

  unsigned int ngx; //!< Number of global x sizes
  unsigned int ngy; //!< Number of global y sizes
  unsigned int ngz; //!< Number of global z sizes
  unsigned int rank; //!< MPI rank of this mesh
  unsigned int n_rank; //!< Number of global ranks

  unsigned int n_cell; //!< Number of local cells
  unsigned int n_global; //!< Nuber of global cells
  
  unsigned int on_rank_start; //!< Start of global index on rank
  unsigned int on_rank_end; //!< End of global index on rank

  std::vector<double> m_census_E; //!< Census energy vector
  std::vector<double> m_emission_E; //!< Emission energy vector
  std::vector<double> m_source_E; //!< Source energy vector

  Cell *cells; //!< Cell data allocated with MPI_Alloc
  std::vector<Cell> cell_list; //!< On processor cells
  std::vector<Cell> new_cell_list; //!< New received cells
  std::vector<unsigned int> remove_cell_list; //!< Cells to be removed
  std::vector<unsigned int> off_rank_bounds; //!< Ending value of global ID for each rank
  std::vector<unsigned int> boundary_cells; //!< Index of adjacent ghost cells

  double m_opA; //!< Opacity coefficient A in A + B^C
  double m_opB; //!< Opacity coefficient B in A + B^C
  double m_opC; //!< Opacity coefficient C in A + B^C
  double m_opS; //!< Scattering opacity constant

  double total_photon_E; //!< Total photon energy on the mesh

  MPI::Datatype MPI_Cell; //!< MPI type, allows simpler function calls
};

#endif //mesh_h_
