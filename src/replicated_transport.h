//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   replicated_transport.h
 * \author Alex Long
 * \date   March 3 2017
 * \brief  IMC transport on a replicated domain
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#ifndef transport_replicated_h_
#define transport_replicated_h_

#include <algorithm>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "RNG.h"
#include "comb_photons.h"
#include "constants.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"
#include "source.h"
#include "low_precision_functions.h"


int comm_particles(std::vector<Cell> &cells, std::vector<Photon> &comm_particles, int max_steps) {
  // 10 CYCLES PER PARTICLE COMM
  std::vector<int> comms_for_cells(cells.size(), 0);
  for (auto &particle : comm_particles) {
    auto cell_index = particle.get_cell();
    cells[cell_index].push_particle(particle);
    comms_for_cells[cell_index]++;
  }

  // clear wait list now that everything has been given to the correct cell
  comm_particles.clear();

  // get the maximum number of comms, that's the maximum of "comms_for_cell" times the max steps
  auto max_element = std::max_element(comms_for_cells.begin(), comms_for_cells.end());
  auto max_cell = std::distance(comms_for_cells.begin(), std::max_element(comms_for_cells.begin(), comms_for_cells.end()));
  auto max_particles = *max_element;
  int max_comms = max_particles*max_steps;

  // determines number of cycles for this communication phase
  int comm_cycles = Constants::PARTICLE_COMM_CYCLES*max_comms;
  return comm_cycles;
}

Constants::event_type transport_photon(Photon &phtn, Mesh &mesh, RNG *rng,
                                       double &next_dt, double &exit_E,
                                       double &census_E,
                                       std::vector<double> &rank_abs_E,
                                       std::vector<double> &rank_track_E,
                                        uint64_t &cycle_count) {
  // cycle counting
  using Constants::LOG_CYCLES;
  using Constants::EXP_CYCLES;
  using Constants::RNG_CYCLES;
  using Constants::DIVIDE_CYCLES;

  using Constants::ELEMENT;
  using Constants::REFLECT;
  using Constants::VACUUM;
  using Constants::BOUNDARY;
  using Constants::SOURCE;
  // events
  using Constants::bc_type;
  using Constants::c;
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using std::min;

  uint32_t cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;

  uint32_t surface_cross = 0;
  constexpr double cutoff_fraction = 0.01; // treat this as constant for now (standard value)

  // initialize cell data---this should be 0 CYCLES if each core represents one cell (it should be in
  // memory and does not need to be loaded from the particle state)
  cell_id = phtn.get_cell();
  const Cell *cell = mesh.get_cell_ptr(cell_id);
  bool active = true;

  // transport this photon
  while (active) {
    // 4 LOADS, 4 CYCLES
    group = phtn.get_group();
    sigma_a = cell->get_op_a(group);
    sigma_s = cell->get_op_s(group);
    f = cell->get_f();
    cycle_count+= 4;

    // get distance to each event type

    // distance to scattering event
    dist_to_scatter =
        -log(rng->generate_random_number()) / ((1.0 - f) * sigma_a + sigma_s);
    // log operation, random number generation, division
    cycle_count += LOG_CYCLES + RNG_CYCLES + DIVIDE_CYCLES;

    // distance to cell boundary event
    dist_to_boundary = cell->get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    // loop over three dimensions, each iteration has divide, 3 LOADS and a min operation (estimating the min as 3 cycles)
    cycle_count += 3*(3 + DIVIDE_CYCLES + 3);

    // distance until end of timestep (time remaining * c)
    dist_to_census = phtn.get_distance_remaining();
    // just one LOAD of particle's m_life_dx field
    cycle_count += 1;

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));
    // assume this takes 4 cycles, two compares, two possible swaps?
    cycle_count += 4;

    // calculate energy absorbed by material, update photon and material energy
    // and update the path-length weighted tally for T_r
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);
    // load photon energy, 3 multiples, one exponential, one subtract
    cycle_count += 1 + 3 + EXP_CYCLES;

    rank_track_E[cell_id] += absorbed_E / (sigma_a * f);
    rank_abs_E[cell_id] += absorbed_E;
    // one multiply one divide, all data probably in registers
    cycle_count += 1 + DIVIDE_CYCLES;

    phtn.set_E(phtn.get_E() - absorbed_E);
    // one subtract
    cycle_count += 1;

    // update position
    phtn.move(dist_to_event);
    // update each dimension (x,y,z) independently, 6 reads, 3 multiply-adds, one subtraction for time remaining
    cycle_count += 6 + 3 + 3 + 1;

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
        get_uniform_angle(rng, angle);
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
          // in mesh-computer IMC this is now a terminating event (don't need to update cell here)
          event = BOUNDARY;
          active = false;
        } else if (boundary_event == VACUUM || boundary_event == SOURCE) {
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

    // UPDATE ME 3 cycles for below cutoff
    // UPDATE ME 20 cycles for scatter event
    // UPDATE ME 5 cycles for boundary crossing event
    // UPDATE ME 4 cycles for exiting problem
    // UPDADE ME 6 cycles for census
    cycle_count += 3 + 20 + 5 + 4 + 6;
  }   // end while alive
  return event;
}

std::vector<Photon> replicated_transport(Mesh &mesh,
                                         IMC_State &imc_state,
                                         std::vector<double> &rank_abs_E,
                                         std::vector<double> &rank_track_E,
                                         const uint64_t max_census_photons) {
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::WAIT;
  using Constants::BOUNDARY;

  // cycle counting
  using Constants::PARTICLE_COMM_CYCLES;

  using std::cout;
  using std::endl;
  using std::vector;

  double census_E = 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state.get_next_dt(); //! Set for census photons

  RNG *rng = imc_state.get_rng();

  // timing
  Timer t_transport;
  t_transport.start_timer("timestep transport");

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  vector<Photon> census_list;   //! End of timestep census list
  vector<Photon> boundary_wait_list;   //! Particles moving between cores
  event_type event;
  std::vector<Cell> &cells = mesh.get_cells();
  uint64_t global_cycle_count = 0; //! Cycle count across all cores
  uint64_t total_comm_cycle_count = 0; //! Total cycle counts from comms this step

  int sub_step = 0;
  while(std::any_of(cells.begin(), cells.end(), [] (const Cell &cell) {return cell.still_working();})) {

    uint64_t step_max_cycle_count = 0; //! Max cycles by a core this transport phase
    uint64_t step_min_cycle_count = UINT64_MAX; //! Min cycles by a core this transport phase
    std::vector<uint64_t> cell_sends(cells.size(), 0); //! Number of particles this cell has to comm
    // for each cell, run the particles until they die or reach a cell boundary
    for (auto &cell : cells) {
      std::vector<Photon> &particle_list = cell.get_particles();
      //std::cout<<"particle list size: "<<particle_list.size()<<std::endl;
      uint64_t cycle_count =0;
      for (auto &phtn : particle_list) {
        event = transport_photon(phtn, mesh, rng, next_dt, exit_E, census_E,
                               rank_abs_E, rank_track_E, cycle_count);

        switch (event) {
        // this case should never be reached
        case WAIT:
          break;
        case KILL:
          break;
        case EXIT:
          break;
        case CENSUS:
          census_list.push_back(phtn);
          break;
        case BOUNDARY:
          boundary_wait_list.push_back(phtn);
          cell_sends[cell.get_ID()]++;
        }
      }

      //std::cout<<"cell done, census size: "<<census_list.size()<<" waiting list size: "<<boundary_wait_list.size();
      //std::cout<<" cycles on this cell: "<<cycle_count<<std::endl;

      step_max_cycle_count = std::max(cycle_count, step_max_cycle_count);
      step_min_cycle_count = std::min(cycle_count, step_min_cycle_count);
      // particles should either be on comm list, census list or dead so clear this vector
      particle_list.clear();
    }

    global_cycle_count += step_max_cycle_count;

    // communicate boundary cells
    // the max distance between physically adjacent cells, the array of APEs and cells are not 1-1
    // (as would be the case for 3D mesh on 2D apes)
    size_t comm_size = boundary_wait_list.size();
    int max_steps = mesh.get_max_distance_between_cells();
    int particle_comm_cycles = comm_particles(cells, boundary_wait_list, max_steps);
    auto max_itr = std::max_element(cell_sends.begin(), cell_sends.end());
    uint64_t n_max_sends = *max_itr;
    size_t max_cell = std::distance(cell_sends.begin(), max_itr);
    global_cycle_count += particle_comm_cycles;
    total_comm_cycle_count += particle_comm_cycles;
    std::cout<<"max particle sends: "<<n_max_sends<<" , max cell: "<<max_cell<<" comm cycles: "<<particle_comm_cycles<<" comm size: "<<comm_size<<std::endl;
    std::cout<<"sub step: "<<sub_step<<" , min cycle count: "<<step_min_cycle_count<<" , max cycle count: ";
    std::cout<<step_max_cycle_count<<" , global cycle count: "<<global_cycle_count<<std::endl;
    sub_step++;
  }

  std::cout<<"Total sub-steps: "<<sub_step<<" total cycle count: "<<global_cycle_count<<" total comm cycles: "<<total_comm_cycle_count<<std::endl;

  // comb the photon population to keep it from growing unbounded. I hardcode
  // the taget comb value to be 10% of the user requested photons divided
  // by the number of replicated ranks
  comb_photons(census_list, max_census_photons, rng);

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());

  // set diagnostic quantities
  imc_state.set_exit_E(exit_E);
  imc_state.set_post_census_E(census_E);
  imc_state.set_census_size(census_list.size());
  imc_state.set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  return census_list;
}

#endif // def transport_replicated_h_
//---------------------------------------------------------------------------//
// end of transport_replicated.h
//---------------------------------------------------------------------------//
