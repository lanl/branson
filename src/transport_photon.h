//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   particle_pass_transport.h
 * \author Alex Long
 * \date   December 1 2015
 * \brief  IMC transport with particle passing method
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_photon_h_
#define transport_photon_h_

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "RNG.h"
#include "constants.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"

//----------------------------------------------------------------------------//
//! Transport a photon when the mesh is always available
Constants::event_type transport_photon_particle_pass(
    Photon &phtn, const Mesh &mesh, RNG *rng, double &next_dt, double &exit_E,
    double &census_E, std::vector<double> &rank_abs_E,
    std::vector<double> &rank_track_E) {
  using Constants::ELEMENT;
  using Constants::PROCESSOR;
  using Constants::REFLECT;
  using Constants::VACUUM;
  // events
  using Constants::bc_type;
  using Constants::c;
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::PASS;
  using std::min;

  uint32_t cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;
  const Cell *cell;

  uint32_t surface_cross = 0;
  double cutoff_fraction = 0.01; // note: get this from IMC_state

  cell_id = phtn.get_cell();
  cell = mesh.get_cell_ptr(cell_id); // note: only for on rank mesh data
  bool active = true;

  // transport this photon
  while (active) {
    group = phtn.get_group();
    sigma_a = cell->get_op_a(group);
    sigma_s = cell->get_op_s(group);
    f = cell->get_f();

    // get distance to event
    dist_to_scatter =
        -log(rng->generate_random_number()) / ((1.0 - f) * sigma_a + sigma_s);

    dist_to_boundary = cell->get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);

    rank_track_E[cell_id] += absorbed_E / (sigma_a * f);
    rank_abs_E[cell_id] += absorbed_E;

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      rank_abs_E[cell_id] += phtn.get_E();
      active = false;
      event = KILL;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, *cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell->get_bc(surface_cross);
        if (boundary_event == ELEMENT) {
          next_cell = cell->get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          cell_id = next_cell;
          cell = mesh.get_cell_ptr(cell_id); // note: only for on rank mesh data
        } else if (boundary_event == PROCESSOR) {
          active = false;
          // set correct cell index with global cell ID
          next_cell = cell->get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          event = PASS;
        } else if (boundary_event == VACUUM) {
          exit_E += phtn.get_E();
          active = false;
          event = EXIT;
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        phtn.set_distance_to_census(c * next_dt);
        active = false;
        event = CENSUS;
        census_E += phtn.get_E();
      }
    } // end event loop
  }   // end while alive
  return event;
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//! Transport a single photon until it has a terminating event (kill, exit,
// wait for data, census)
Constants::event_type
transport_photon_mesh_pass(Photon &phtn, const Mesh &mesh, RNG *rng,
                           double &next_dt, double &exit_E, double &census_E,
                           std::vector<double> &rank_abs_E,
                           std::vector<double> &rank_track_E,
                           std::unordered_map<uint32_t, double> &off_rank_abs_E)

{
  using Constants::ELEMENT;
  using Constants::PROCESSOR;
  using Constants::REFLECT;
  using Constants::VACUUM;
  // events
  using Constants::bc_type;
  using Constants::c;
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::WAIT;
  using std::min;

  uint32_t cell_id, next_cell, next_grip;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;
  Cell cell;

  uint32_t surface_cross = 0;
  uint32_t on_rank_start = mesh.get_offset();
  const double cutoff_fraction = 0.01; // note: get this from IMC_state

  cell_id = phtn.get_cell();
  cell = mesh.get_on_rank_cell(cell_id);
  bool active = true;

  // transport this photon
  while (active) {
    group = phtn.get_group();
    sigma_a = cell.get_op_a(group);
    sigma_s = cell.get_op_s(group);
    f = cell.get_f();

    // get distance to event
    dist_to_scatter =
        -log(rng->generate_random_number()) / ((1.0 - f) * sigma_a + sigma_s);

    dist_to_boundary = cell.get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);

    // process on rank tallies as usual
    if (mesh.on_processor(cell_id)) {
      rank_track_E[cell_id - on_rank_start] += absorbed_E / (sigma_a * f);
      rank_abs_E[cell_id - on_rank_start] += absorbed_E;
    } else
      off_rank_abs_E[cell_id] += absorbed_E;

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      if (mesh.on_processor(cell_id))
        rank_abs_E[cell_id - on_rank_start] += phtn.get_E();
      else
        off_rank_abs_E[cell_id] += phtn.get_E();
      active = false;
      event = KILL;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
        // if effective scatter change frequency
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell.get_bc(surface_cross);
        if (boundary_event == ELEMENT || boundary_event == PROCESSOR) {
          next_cell = cell.get_next_cell(surface_cross);
          next_grip = cell.get_next_grip(surface_cross);
          phtn.set_cell(next_cell);
          phtn.set_grip(next_grip);
          cell_id = next_cell;
          // look for this cell, if it's not there transport later
          if (mesh.mesh_available(cell_id))
            cell = mesh.get_on_rank_cell(cell_id);
          else {
            event = WAIT;
            active = false;
          }
        } else if (boundary_event == VACUUM) {
          active = false;
          exit_E += phtn.get_E();
          event = EXIT;
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        phtn.set_distance_to_census(c * next_dt);
        active = false;
        census_E += phtn.get_E();
        event = CENSUS;
      }
    } // end event loop
  }   // end while alive
  return event;
}
//----------------------------------------------------------------------------//

#endif // def transport_photon_h_
//----------------------------------------------------------------------------//
// end of transport_photon.h
//----------------------------------------------------------------------------//
