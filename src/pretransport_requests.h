//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   pretransport_requests.h
 * \author Alex Long
 * \date   February 14 2017
 * \brief  Requests mesh for needed paricles before transport
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef pretransport_requests_h_
#define pretransport_requests_h_

#include "mesh.h"
#include "mesh_request_manager.h"
#include "mesh_rma_manager.h"
#include "message_counter.h"
#include "photon.h"
#include "work_packet.h"

//! One-sided request for ajdacent rank data and cells where your work is
void pretransport_requests(std::vector<Work_Packet>& work,
  std::vector<Photon>& census_list, const Mesh &mesh, RMA_Manager &rma_manager,
  Message_Counter& mctr)
{

  for (auto const &i_w : work) {
    uint32_t grip_ID = i_w.get_global_grip_ID();
    if (!mesh.on_processor(grip_ID))
      rma_manager.request_cell_rma(grip_ID, mctr);
  }

  for (auto const &i_p : census_list) {
    uint32_t grip_ID = i_p.get_grip();
    if (!mesh.on_processor(grip_ID))
      rma_manager.request_cell_rma(grip_ID, mctr);
  }

  const Cell * const cell_ptr = mesh.get_const_cells_ptr();
  for (uint32_t i = 0; i<mesh.get_n_local_cells(); ++i) {
    const Cell& cell = cell_ptr[i];
    for (uint32_t dir = 0; dir<6;++dir) {
      if (!mesh.on_processor(cell.get_next_grip(dir)))
        rma_manager.request_cell_rma(cell.get_next_grip(dir), mctr);
    }
  }

}

//! Two-sided request for ajdacent rank data and cells where your work is
void pretransport_requests(std::vector<Work_Packet>& work,
  std::vector<Photon>& census_list, const Mesh &mesh,
  Mesh_Request_Manager &req_manager, Message_Counter& mctr)
{

  for (auto const &i_w : work) {
    uint32_t grip_ID = i_w.get_global_grip_ID();
    if (!mesh.on_processor(grip_ID))
      req_manager.request_cell(grip_ID, mctr);
  }

  for (auto const &i_p : census_list) {
    uint32_t grip_ID = i_p.get_grip();
    if (!mesh.on_processor(grip_ID))
      req_manager.request_cell(grip_ID, mctr);
  }

  const Cell * const cell_ptr = mesh.get_const_cells_ptr();
  for (uint32_t i = 0; i<mesh.get_n_local_cells(); ++i) {
    const Cell& cell = cell_ptr[i];
    for (uint32_t dir = 0; dir<6;++dir) {
      if (!mesh.on_processor(cell.get_next_grip(dir)))
        req_manager.request_cell(cell.get_next_grip(dir), mctr);
    }
  }
}

#endif // pretransport_requests_h_

//----------------------------------------------------------------------------//
// end of pretransport_requests.h
//----------------------------------------------------------------------------//
