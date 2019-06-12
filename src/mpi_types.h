//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   mpi_types.h
 * \author Alex Long
 * \date   May 12 2016
 * \brief  Creates and provides access to custom MPI dataypes
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef mpi_types_h_
#define mpi_types_h_

//==============================================================================
/*!
 * \class MPI_Types
 * \brief Creates, commits and provides access to custom MPI types
 *
 * MPI allows custom types to be specified and then used directly in MPI send
 * receive and get calls. This class creates MPI types needed for particles,
 * cells and work packets
 * \example no test yet
 */
//==============================================================================

class MPI_Types {
public:
  //! constructor
  MPI_Types(void) {

    // make and commit the MPI proto cell type
    {
      MPI_Datatype og_MPI_Proto_Cell;

      // remake the MPI cell datatype from mesh
      const int cell_entry_count = 3;
      // 16 uint32_t, 6 int, 6 doubles
      int cell_array_of_block_length[3] = {16, 6, 6};
      // Displacements of each type in the cell
      MPI_Aint cell_array_of_block_displace[3] = {
          0, 16 * sizeof(uint32_t), 16 * sizeof(uint32_t) + 6 * sizeof(int)};
      //Type of each memory block
      MPI_Datatype cell_array_of_types[3] = {MPI_UNSIGNED, MPI_INT, MPI_DOUBLE};

      MPI_Type_create_struct(cell_entry_count, cell_array_of_block_length,
                             cell_array_of_block_displace, cell_array_of_types,
                             &og_MPI_Proto_Cell);

      // Commit the type to MPI so it recognizes it in communication calls
      MPI_Type_commit(&og_MPI_Proto_Cell);

      MPI_Type_size(og_MPI_Proto_Cell, &mpi_proto_cell_size);
      // Duplicate the type so it's recognized when returned out of this
      // context (I don't know why this is necessary)
      MPI_Type_dup(og_MPI_Proto_Cell, &MPI_Proto_Cell);
    }

    // make and commit the MPI region type
    {
      // make the Region
      const int region_entry_count = 2;
      // 2 uint32_t, 9 doubles
      int region_array_of_block_length[2] = {2, 9};
      // Displacements of each type in the cell
      MPI_Aint region_array_of_block_displace[2] = {0, 2 * sizeof(uint32_t)};
      //Type of each memory block
      MPI_Datatype region_array_of_types[2] = {MPI_UNSIGNED, MPI_DOUBLE};

      MPI_Datatype og_MPI_Region;
      MPI_Type_create_struct(region_entry_count, region_array_of_block_length,
                             region_array_of_block_displace,
                             region_array_of_types, &og_MPI_Region);

      // Commit the type to MPI so it recognizes it in communication calls
      MPI_Type_commit(&og_MPI_Region);

      MPI_Type_size(og_MPI_Region, &mpi_region_size);
      // Duplicate the type so it's recognized when returned out of this
      // context (I don't know why this is necessary)
      MPI_Type_dup(og_MPI_Region, &MPI_Region);
    }
  }

  // destructor
  ~MPI_Types() {
    MPI_Type_free(&MPI_Proto_Cell);
    MPI_Type_free(&MPI_Region);
  }

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Return reference to MPI_Proto_Cell for use in communication calls
  MPI_Datatype get_proto_cell_type(void) const { return MPI_Proto_Cell; }

  //! Return reference to Region for use in communication calls
  MPI_Datatype get_region_type(void) const { return MPI_Region; }

  //! Return size of MPI cell in bytes
  int get_proto_cell_size(void) const { return mpi_proto_cell_size; }

  //! Return size of MPI Region in bytes
  int get_region_size(void) const { return mpi_region_size; }

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  MPI_Datatype MPI_Proto_Cell; //!< Custom MPI datatype for proto mesh cell
  MPI_Datatype MPI_Region;     //!< Custom MPI datatype for region
  int mpi_proto_cell_size;     //!< Size of MPI_Cell datatype
  int mpi_region_size;         //!< Size of MPI_Region datatype
};

#endif // mpi_types_h_
//----------------------------------------------------------------------------//
// end of mpi_types.h
//----------------------------------------------------------------------------//
