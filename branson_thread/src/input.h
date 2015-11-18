/*
  Author: Alex Long
  Date: 7/18/2014
  Name: input.h
*/

#ifndef Input_h_
#define Input_h_

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <stdlib.h>
#include <string>

#include "constants.h"

using std::cout;
using std::endl;
using std::vector;
using boost::property_tree::ptree;
using Constants::a;
using Constants::c;
using Constants::bc_type;
using Constants::dir_type;
using Constants::VACUUM; using Constants::REFLECT; using Constants::ELEMENT;

class Input
{
  public:
  Input( std::string fileName)
  {
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
        n_photons = v.second.get<int>("photons");
        seed = v.second.get<int>("seed");
        output_freq = v.second.get<int>("output_frequency",1);
        tempString = v.second.get<std::string>("tilt", std::string("FALSE"));       
        if (tempString == "TRUE") use_tilt = 1;
        else use_tilt = 0;
        tempString = v.second.get<std::string>("use_combing", std::string("TRUE"));
        if (tempString == "TRUE") use_comb= 1;
        else use_comb = 0;
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
        n_x_elements = v.second.get<unsigned int>("n_x_elements");
        dx = v.second.get<double>("dx");

        n_y_elements = v.second.get<unsigned int>("n_y_elements");
        dy = v.second.get<double>("dy");

        n_z_elements = v.second.get<unsigned int>("n_z_elements");
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
      } 

      //read in source data
      else if(v.first == "source")
      {
        T_source = v.second.get<double>("T_source", 0.0);
        source_element = v.second.get<unsigned int>("source_element", 0);
      }
    } //end xml parse
  }
  
  ~Input() {};
  
  void print_problem_info(void)
  {
    cout<<"OpenMP Threads: "<<omp_get_max_threads()<<endl;
    cout<<"Problem Specifications:";
    cout<<"Constants -- c: "<<c<<" (cm/sh) , a: "<<a <<endl;
    cout<<"Run Parameters-- Photons: "<<n_photons<<", time finish: "<<tFinish<<" (sh), time step: "<<dt<<" (sh) ,"<<endl;
    cout<<" time multiplier: "<<tMult<<" , max dt:"<<dtMax<<" (sh), Random number seed: "<<seed<<" , output frequency: "<<output_freq<<endl;
    
    cout<<"material temperature: "<<Tm_initial<<" (keV), radiation temperature: "<<Tr_initial<<" (keV)"<<endl;

    cout<<"Sampling -- Emission Position: "; 
    if (use_tilt) cout<<"source tilting (x only), ";
    else cout<<"uniform (default), ";
    cout<<" Angle: ";

    if (use_comb) cout<<"Combing census enabled (default)"<<endl;
    else cout<<"No combing"<<endl;

    if (print_verbose) cout<<"Verbose printing mode enabled"<<endl;
    else cout<<"Terse printing mode (default)"<<endl;

    cout<<"Spatial Information -- elements x,y,z: "<<n_x_elements<<" "<<n_y_elements<<" "<<n_z_elements<<endl;
    cout<<"dx: "<< dx<<" dy: "<<dy<<" dz: "<<dz<<endl;

    cout<<"Material Information -- heat capacity: "<<CV<<" opacity constants: "<<opacA<<" + "<<opacB
        <<"^"<<opacC<<endl;
    
    cout<<endl;
  }

  int get_n_x_elements(void) const {return n_x_elements;}
  int get_n_y_elements(void) const {return n_y_elements;}
  int get_n_z_elements(void) const {return n_z_elements;}

  double get_dx(void) const {return dx;}
  double get_dy(void) const {return dy;}
  double get_dz(void) const {return dz;}

  double get_initial_Tm(void) const {return Tm_initial;}
  double get_initial_Tr(void) const {return Tr_initial;}
  int get_output_freq(void) const {return output_freq;}

  bool get_tilt_bool(void) const {return use_tilt;}
  bool get_comb_bool(void) const {return use_comb;}
  bool get_verbose_print_bool(void) const {return print_verbose;}
  bool get_print_mesh_info_bool(void) const {return print_mesh_info;}

  double get_dt(void) const {return dt;}
  double get_time_start(void) const {return tStart;}
  double get_time_finish(void) const {return tFinish;}
  double get_time_mult(void) const {return tMult;}
  double get_dt_max(void) const {return dtMax;}
  int get_number_photons(void) const {return n_photons;}
  int get_rng_seed(void) const {return seed;}

  //source functions
  unsigned int get_source_element(void) const {return source_element;}
  double get_source_T(void) const {return T_source;}

  //material functions
  double get_CV(void) {return CV;}
  double get_rho(void) {return rho;}
  double get_opacity_A(void) {return opacA;}
  double get_opacity_B(void) {return opacB;}
  double get_opacity_C(void) {return opacC;}

  bc_type get_bc(const dir_type& direction) const { return bc[direction];}

  private:

  // geometry
  unsigned int n_x_elements;
  unsigned int n_y_elements;
  unsigned int n_z_elements;

  double dx;
  double dy;
  double dz;

  bc_type bc[6];  

  // timing
  double tStart;
  double dt;
  double tFinish;
  double tMult;
  double dtMax;

  // initial conditions
  double Tm_initial;
  double Tr_initial;

  //material
  double rho;
  double CV;
  double opacA;
  double opacB;
  double opacC;

  //source
  double T_source;
  double source_element;

  // Monte Carlo parameters
  unsigned int n_photons;
  unsigned int seed;

  // Method parameters
  bool use_tilt;
  bool use_comb;

  // Debug paramters
  int output_freq;
  bool print_verbose;
  bool print_mesh_info;

};

#endif // Input_h_
