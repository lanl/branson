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

//! Request data for mesh that you particles live on right now
void pretransport_requests(std::vector<Work_Packet>& work,
  std::vector<Photon>& census_list, Mesh *mesh, RMA_Manager *rma_manager,
  Message_Counter& mctr)
{

  for (auto const &i_w : work) {
    uint32_t grip_ID = i_w.get_global_grip_ID();
    if (!mesh->on_processor(grip_ID))
      rma_manager->request_cell_rma(grip_ID, mctr);
  }

  for (auto const &i_p : census_list) {
    uint32_t grip_ID = i_p.get_grip();
    if (!mesh->on_processor(grip_ID))
      rma_manager->request_cell_rma(grip_ID, mctr);
  }
}

//! Request data for mesh that you particles live on right now
void pretransport_requests(std::vector<Work_Packet>& work,
  std::vector<Photon>& census_list, Mesh *mesh,
  Mesh_Request_Manager *req_manager, Message_Counter& mctr)
{

  for (auto const &i_w : work) {
    uint32_t grip_ID = i_w.get_global_grip_ID();
    if (!mesh->on_processor(grip_ID))
      req_manager->request_cell(grip_ID, mctr);
  }

  for (auto const &i_p : census_list) {
    uint32_t grip_ID = i_p.get_grip();
    if (!mesh->on_processor(grip_ID))
      req_manager->request_cell(grip_ID, mctr);
  }
}

#endif // pretransport_requests_h_

//----------------------------------------------------------------------------//
// end of pretransport_requests.h
//----------------------------------------------------------------------------//
