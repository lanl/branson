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
#include "config.h"

class GPU_Setup {

public:
  //! Constructor
  GPU_Setup(const int rank, const int n_ranks, const bool use_gpu_transporter, const std::vector<Cell> &cpu_cells)
    : m_use_gpu_transporter(use_gpu_transporter), device_cells_ptr(nullptr)
  {
#ifdef USE_CUDA
    if(m_use_gpu_transporter) {
      // MPI rank to GPU mapping
      set_device_ID(rank, n_ranks);

      std::cout<<"Allocating and transferring "<<cpu_cells.size()<<" cell(s) to the GPU"<<std::endl;
      // allocate and copy cells
      cudaError_t err = cudaMalloc((void **)&device_cells_ptr, sizeof(Cell) * cpu_cells.size());
      Insist(!err, "CUDA error in allocating cells data");
      err = cudaMemcpy(device_cells_ptr, cpu_cells.data(), sizeof(Cell) * cpu_cells.size(),
                       cudaMemcpyHostToDevice);
      Insist(!err, "CUDA error in copying cells data");
    }
#endif
  }



  //! Destructor
  ~GPU_Setup() {
#ifdef USE_CUDA
    if(m_use_gpu_transporter) {
      cudaFree(device_cells_ptr);
    }
#endif
  }

  Cell *get_device_cells_ptr() const {return device_cells_ptr;}
  bool use_gpu_transporter() const {return m_use_gpu_transporter;}

private:

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Assign GPUs to MPI ranks
 *
 * If there are more ranks on a node than GPUs available on a node, multiple ranks will be allowed
 * to use the same GPU.
 *
 * \param[in] rank MPI rank of calling processor
 * \param[in] n_ranks total number of MPI ranks
 */
void set_device_ID(const int rank, const int n_ranks) {
#ifdef USE_CUDA
  // device set
  int n_devices;
  cudaGetDeviceCount(&n_devices);
  int my_device = 0;
  if (n_ranks <= n_devices)
    my_device = rank;
  else
    my_device = rank % n_devices;
  cudaSetDevice(my_device);
  int my_bus_id = 0;
  cudaDeviceGetAttribute(&my_bus_id, cudaDevAttrPciBusId, my_device);
  std::cout << "N devices: " << n_devices << std::endl;
  std::cout << "rank: " << rank << ", device: " << my_device << " bus id: ";
  std::cout << my_bus_id << std::endl;
#endif
}


  bool m_use_gpu_transporter;
  Cell *device_cells_ptr;
};

#endif // gpu_setup_h_
//----------------------------------------------------------------------------//
// end of gpu_setup.h
//----------------------------------------------------------------------------//

