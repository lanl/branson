/*
  Author: Alex Long
  Date: 7/18/2014
  Name: input.h
*/

#ifndef input_h_
#define input_h_

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "constants.h"

//==============================================================================
/*!
 * \class Input
 * \brief Reads input data from an XML file and stores that data
 *
 * Boost's XML parser is used to read an XML input file. This class stores that
 * information and provides functions to access it. This class also prints the 
 * problem information.
 */
//==============================================================================
class Input
{
  public:
  Input( std::string fileName)
  {
    using boost::property_tree::ptree;
    using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT;
    using Constants::X_POS;  using Constants::Y_POS; using Constants::Z_POS;
    using Constants::X_NEG;  using Constants::Y_NEG; using Constants::Z_NEG;
    using Constants::PARTICLE_PASS;
    using Constants::CELL_PASS;
    using std::cout;  
    using std::endl;

    ptree pt;
    read_xml(fileName, pt);
    std::string tempString;
    // traverse pt
    BOOST_FOREACH( ptree::value_type const& v, pt.get_child("prototype") ) 
    {
      //read in basic problem parameters      
      if(v.first =="common")
      {
        tFinish = v.second.get<double>("t_stop");
        dt = v.second.get<double>("dt_start");
        tStart= v.second.get<double>("t_start");
        tMult = v.second.get<double>("t_mult" , 1.0);
        dtMax = v.second.get<double>("dt_max" , dt);
        Tm_initial = v.second.get<double>("Tm_initial");
        Tr_initial = v.second.get<double>("Tr_initial");
        n_photons =v.second.get<int>("photons"); 
        seed = v.second.get<int>("seed");
        map_size = v.second.get<int>("map_size");
        output_freq = v.second.get<int>("output_frequency",1);
        tempString = v.second.get<std::string>("tilt", std::string("FALSE"));
        if (tempString == "TRUE") use_tilt = 1;
        else use_tilt = 0;
        tempString = v.second.get<std::string>("use_combing", 
                                                std::string("TRUE"));
        if (tempString == "TRUE") use_comb= 1;
        else use_comb = 0;

        //stratified sampling
        tempString = v.second.get<std::string>("stratified_sampling", 
                                                std::string("FALSE"));
        if (tempString == "TRUE") use_strat= true;
        else use_strat = false;

        //ghost map flag
        tempString = v.second.get<std::string>("use_ghost_map", 
                                               std::string("FALSE"));
        if (tempString == "TRUE") use_ghost_cells= true;
        else use_ghost_cells = false;

        //number of particles to run between MPI message checks
        batch_size = v.second.get<unsigned int>("batch_size", 100);

        //preferred number of particles per MPI send
        particle_message_size = v.second.get<unsigned int>("particle_message_size", 100);

        // domain decomposed transport aglorithm
        tempString = v.second.get<std::string>("dd_transport_type", 
                                               std::string("CELL_PASS"));
        if (tempString == "CELL_PASS") dd_mode = CELL_PASS;
        else if (tempString == "PARTICLE_PASS") dd_mode = PARTICLE_PASS;
        else {
          cout<<"WARNING: Domain decomposition method not recognized ";
          cout<<"setting to PARTICLE PASSING method";
          dd_mode = PARTICLE_PASS;
        }
      } //end common

      //read in basic problem parameters      
      if(v.first =="debug_options")
      {
        tempString = v.second.get<std::string>("print_verbose", "FALSE");
        if (tempString == "TRUE") print_verbose = true;
        else print_verbose = false;
        tempString = v.second.get<std::string>("print_mesh_info", "FALSE");
        if (tempString == "TRUE") print_mesh_info = true;
        else print_mesh_info = false;
      } //end common

      //read in spatial data
      if(v.first =="spatial")
      {
        n_x_cells = v.second.get<unsigned int>("n_x_cells");
        dx = v.second.get<double>("dx");

        n_y_cells = v.second.get<unsigned int>("n_y_cells");
        dy = v.second.get<double>("dy");

        n_z_cells = v.second.get<unsigned int>("n_z_cells");
        dz = v.second.get<double>("dz");

        //read in boundary conditions
        bool b_error = false;
        tempString = v.second.get<std::string>("bc_right");
        if      (tempString == "REFLECT") bc[X_POS] = REFLECT;
        else if (tempString == "VACUUM")  bc[X_POS] = VACUUM; 
        else    b_error = true;
        
        tempString = v.second.get<std::string>("bc_left");
        if      (tempString == "REFLECT") bc[X_NEG] = REFLECT;
        else if (tempString == "VACUUM")  bc[X_NEG] = VACUUM; 
        else    b_error = true;

        tempString = v.second.get<std::string>("bc_up");
        if      (tempString == "REFLECT") bc[Y_POS] = REFLECT;
        else if (tempString == "VACUUM")  bc[Y_POS] = VACUUM; 
        else    b_error = true;

        tempString = v.second.get<std::string>("bc_down");
        if      (tempString == "REFLECT") bc[Y_NEG] = REFLECT;
        else if (tempString == "VACUUM")  bc[Y_NEG] = VACUUM; 
        else    b_error = true;

        tempString = v.second.get<std::string>("bc_top");
        if      (tempString == "REFLECT") bc[Z_POS] = REFLECT;
        else if (tempString == "VACUUM")  bc[Z_POS] = VACUUM; 
        else    b_error = true;

        tempString = v.second.get<std::string>("bc_bottom");
        if      (tempString == "REFLECT") bc[Z_NEG] = REFLECT;
        else if (tempString == "VACUUM")  bc[Z_NEG] = VACUUM; 
        else    b_error = true;

        if (b_error) {
          cout<<"Boundary type not reconginzed"<<endl;
          exit(EXIT_FAILURE);
        }
      }
      
      //read in material data
      else if(v.first == "material")
      {
        CV = v.second.get<double>("CV", 0.0);
        rho = v.second.get<double>("density", 0.0);
        opacA = v.second.get<double>("opacA", 0.0);
        opacB = v.second.get<double>("opacB", 0.0);
        opacC = v.second.get<double>("opacC", 0.0);
        opacS = v.second.get<double>("opacS", 0.0);
      } 

      //read in source data
      else if(v.first == "source")
      {
        T_source = v.second.get<double>("T_source", 0.0);
      }
    } //end xml parse

  }
  
  ~Input() {};
  
  void print_problem_info(void)
  {
    using Constants::a;
    using Constants::c;
    using std::cout;  
    using std::endl;
    using Constants::PARTICLE_PASS;
    using Constants::CELL_PASS;

    cout<<"Problem Specifications:";
    cout<<"Constants -- c: "<<c<<" (cm/sh) , a: "<<a <<endl;
    cout<<"Run Parameters-- Photons: "<<n_photons<<", time finish: "<<tFinish;
    cout<<" (sh), time step: "<<dt<<" (sh) ,"<<endl;

    cout<<" time multiplier: "<<tMult<<" , max dt:"<<dtMax;
    cout<<" (sh), Random number seed: "<<seed;
    cout<<" , output frequency: "<<output_freq<<endl;
    
    cout<<"material temperature: "<<Tm_initial;
    cout<<" (keV), radiation temperature: "<<Tr_initial<<" (keV)"<<endl;

    //cout<<"Sampling -- Emission Position: ";
    //if (use_tilt) cout<<"source tilting (x only), ";
    //else cout<<"uniform (default), ";
    //cout<<" Angle: ";
    cout<<"Stratified sampling in angle: ";
    if(use_strat) cout<<"TRUE"<<endl;
    else cout<<"FALSE"<<endl;

    if (use_comb) cout<<"Combing census enabled (default)"<<endl;
    else cout<<"No combing"<<endl;

    if (print_verbose) cout<<"Verbose printing mode enabled"<<endl;
    else cout<<"Terse printing mode (default)"<<endl;

    cout<<"Spatial Information -- cells x,y,z: "<<n_x_cells<<" ";
    cout<<n_y_cells<<" "<<n_z_cells<<endl;
    cout<<"dx: "<< dx<<" dy: "<<dy<<" dz: "<<dz<<endl;

    cout<<"Material Information -- heat capacity: "<<CV;
    cout<<" opacity constants: "<<opacA<<" + "<<opacB<<"^"<<opacC;
    cout<<", scattering opacity: "<<opacS<<endl;
    
    cout<<"Parallel Information -- DD algorithm: ";
    if (dd_mode == CELL_PASS) {
      cout<<"CELL PASSING"<<endl;
      cout<<"map size: "<<map_size;
      cout<<",  make ghost cell map: ";
      if (use_ghost_cells) cout<<"TRUE";
      else cout<<"FALSE";
      cout<<", Batch size: "<<batch_size;
      cout<<endl;
    }
    else {
      cout<<"PARTICLE PASSING"<<endl;
      cout<<"Batch size: "<<batch_size;
      cout<<", particle message size: "<<particle_message_size;
    }

    cout<<endl;
  }

  int get_n_x_cells(void) const {return n_x_cells;}
  int get_n_y_cells(void) const {return n_y_cells;}
  int get_n_z_cells(void) const {return n_z_cells;}

  double get_dx(void) const {return dx;}
  double get_dy(void) const {return dy;}
  double get_dz(void) const {return dz;}

  double get_initial_Tm(void) const {return Tm_initial;}
  double get_initial_Tr(void) const {return Tr_initial;}
  int get_output_freq(void) const {return output_freq;}

  bool get_tilt_bool(void) const {return use_tilt;}
  bool get_comb_bool(void) const {return use_comb;}
  bool get_stratified_bool(void) const {return use_strat;}
  bool get_ghost_cell_bool(void) const {return use_ghost_cells;}
  bool get_verbose_print_bool(void) const {return print_verbose;}
  bool get_print_mesh_info_bool(void) const {return print_mesh_info;}

  double get_dt(void) const {return dt;}
  double get_time_start(void) const {return tStart;}
  double get_time_finish(void) const {return tFinish;}
  double get_time_mult(void) const {return tMult;}
  double get_dt_max(void) const {return dtMax;}
  int get_number_photons(void) const {return n_photons;}
  int get_rng_seed(void) const {return seed;}
  unsigned int get_batch_size(void) const {return batch_size;}
  unsigned int get_particle_message_size(void) const {return particle_message_size;}

  //source functions
  double get_source_T(void) const {return T_source;}

  //material functions
  double get_CV(void) {return CV;}
  double get_rho(void) {return rho;}
  double get_opacity_A(void) {return opacA;}
  double get_opacity_B(void) {return opacB;}
  double get_opacity_C(void) {return opacC;}
  double get_opacity_S(void) {return opacS;}

  Constants::bc_type get_bc(const Constants::dir_type& direction) const 
  { 
    return bc[direction];
  }
  
  unsigned int get_map_size(void) const {return map_size;}
  unsigned int get_dd_mode(void) const {return dd_mode;}

  private:

  // geometry
  unsigned int n_x_cells; //!< Total x cells
  unsigned int n_y_cells; //!< Total y cells
  unsigned int n_z_cells; //!< Total z cells

  double dx; //!< x grid spacing (cm)
  double dy; //!< y grid spacing (cm)
  double dz; //!< z grid spacing (cm)

  Constants::bc_type bc[6]; //!< Boundary condition array 

  // timing
  double tStart; //!< Starting time
  double dt; //!< Timestep size
  double tFinish; //!< Finish time
  double tMult; //!< Timestep multiplier
  double dtMax; //!< Maximum timestep size

  // initial conditions
  double Tm_initial; //!< Initial material temperature 
  double Tr_initial; //!< Initial radiation temperature

  //material
  double rho; //!< Density (g/cc)
  double CV; //!< Heat capacity (jk/keV/g)
  double opacA; //!< Constant opacity
  double opacB; //!< Opacity temperature multiplier 
  double opacC; //!< Opacity temperature power
  double opacS; //!< Scattering opacity

  //source
  double T_source; //!< Temperature of source

  // Monte Carlo parameters
  unsigned int n_photons; //!< Photons to source each timestep
  unsigned int seed; //!< Random number seed

  // Method parameters
  bool use_tilt; //!< Use tilting for emission sampling
  bool use_comb; //!< Comb census photons
  bool use_strat; //!< Use strafifed sampling

  // Debug parameters
  int output_freq; //!< How often to print temperature information
  bool print_verbose; //!< Verbose printing flag
  bool print_mesh_info; //!< Mesh information printing flag

  //parallel performance parameters
  unsigned int map_size; //!< Size of stored off-rank mesh cells
  unsigned int dd_mode; //!< Mode of domain decomposed transport algorithm
  bool use_ghost_cells; //!< Always keep first ghost cells
  unsigned int batch_size; //!< Particles to run between MPI message checks
  unsigned int particle_message_size; //!< Preferred number of particles in MPI sends

};

#endif // input_h_
