//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   gpu_setup.h
 * \author Alex Long
 * \date   Oct 20 2022
 * \brief  Struct that holds device pointers to particles, mesh and tallies
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef gpu_setup_h_
#define gpu_setup_h_

#include <vector>

#include "cell.h"

class GPU_Setup {

public:
  //! Constructor
  Message_Counter(const bool use_gpu_transporter, std::vector<Cell> &cpu_cells)
    : m_use_gpu_transporter(use_gpu_transporter), device_cells_ptr(nullptr),
  {
#ifdef USE_CUDA
    if(use_gpu_transporter) {

      // allocate and copy cells
      err = cudaMalloc((void **)&device_cells_ptr, sizeof(Cell) * cpu_cells.size());
      Insist(!err, "CUDA error in allocating cells data");
      err = cudaMemcpy(device_cells_ptr, cpu_cells.data(), sizeof(Cell) * cpu_cells.size(),
                       cudaMemcpyHostToDevice);
      Insist(!err, "CUDA error in copying cells data");
    }
  }

  //! Destructor
  ~Message_Counter() {
#ifdef USE_CUDA
    if(use_gpu_transporter) {
      cudaFree(device_cells_ptr);
      cudaFree(device_abs_E_ptr);
      cudaFree(device_track_E_ptr);
    }
  }

  Cell *get_device_cells_ptr() const {return device_cells_ptr;}

  private:
  bool m_use_gpu_transporter;
  Cell *device_cells_ptr;
};

#endif // gpu_setup_h_
//----------------------------------------------------------------------------//
// end of gpu_setup.h
//----------------------------------------------------------------------------//

