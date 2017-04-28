//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   info.h
 * \author Alex Long
 * \date   September 2 2016
 * \brief  Stores MPI information to compress function signatures
 * \note   ***COPYRIGHT_GOES_HERE****
 */
//----------------------------------------------------------------------------//

#ifndef info_h_
#define info_h_

#include <mpi.h>
#include <regex>
#include <string>

#include "config.h"
//==============================================================================
/*!
 * \class Info
 * \brief Stores MPI information
 *
 * This class is used to store global and node based MPI information. This
 * keeps machine specific code in a high-level place. The name MPI_Info is
 * already used by MPI so this is just called Info.
 *
 * \example no test yet
 */
//==============================================================================

class Info
{
  public:
  //! Constructor
  Info(void) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

    // default values (CCS node)
    color = 1;
    node_mem = 16000000000;
    machine_name = "CCS";

    // these ifdef macros are messy but for some reason some compilers don't
    // like the #elif directive
#ifdef SNOW_NODE
    machine_name = "snow";
    char * proc_name = new char[MPI_MAX_PROCESSOR_NAME];
    int result_length;
    MPI_Get_processor_name(proc_name, &result_length);
    std::string p_name(proc_name, result_length);
    std::regex r_ml_node("sn([0-9]*)[.]localdomain");
    std::smatch result;
    std::regex_match(p_name, result, r_ml_node);
    color = std::stoi(result.str(1));
    node_mem = 32000000000;
    delete[] proc_name;
#endif

#ifdef MOONLIGHT_NODE
    machine_name = "moonlight";
    char * proc_name = new char[MPI_MAX_PROCESSOR_NAME];
    int result_length;
    MPI_Get_processor_name(proc_name, &result_length);
    std::string p_name(proc_name, result_length);
    std::regex r_ml_node("ml([0-9]*)[.]localdomain");
    std::smatch result;
    std::regex_match(p_name, result, r_ml_node);
    color = std::stoi(result.str(1));
    node_mem = 32000000000;
    delete[] proc_name;
#endif

#ifdef TRINITITE_NODE
    machine_name = "trinitite";
    char * proc_name = new char[MPI_MAX_PROCESSOR_NAME];
    int result_length;
    MPI_Get_processor_name(proc_name, &result_length);
    std::string p_name(proc_name, result_length);
    std::regex r_tt_node("nid([0-9]*)");
    std::smatch result;
    std::regex_match(p_name, result, r_tt_node);
    color = std::stoi(result.str(1));
    node_mem = 128000000000;
    delete[] proc_name;
#endif

  }

  // Destructor
  ~Info() {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Return global MPI rank of process
  int get_rank(void) const {return rank;}

  //! Return number of total MPI ranks
  int get_n_rank(void) const {return n_rank;}

  //! Return unique color identifier for the current node
  int get_color(void) const {return color;}

  //! Return total memory on this node
  int64_t get_node_mem(void) const {return node_mem;}

  //! Return the name of the machine the branson makefile was configured on
  std::string get_machine_name(void) const {return machine_name;}

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
  private:
  int rank; //!< Global rank ID
  int n_rank; //!< Global number of ranks
  int color; //!< Unique identifier for this node
  int64_t node_mem; //!< Total memory available for this node
  std::string machine_name; //!< Name of the machine compiled on
};

#endif // info_h_
//----------------------------------------------------------------------------//
// end of info.h
//----------------------------------------------------------------------------//
