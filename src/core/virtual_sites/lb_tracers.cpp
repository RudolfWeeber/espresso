/*
 * Copyright (C) 2010-2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_VIRTUAL_SITES_INERTIALESS_TRACERS

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/CellStructure.hpp"
#include "errorhandling.hpp"
#include "lb/Solver.hpp"
#include "lb/particle_coupling.hpp"

#include <utils/math/sqr.hpp>

static bool lb_sanity_checks(LB::Solver const &lb) {
  if (not lb.is_solver_set()) {
    runtimeErrorMsg() << "LB needs to be active for inertialess tracers.";
    return true;
  }
  return false;
}

void lb_tracers_add_particle_force_to_fluid(CellStructure &cell_structure,
                                            BoxGeometry const &box_geo,
                                            LocalBox const &local_box,
                                            LB::Solver &lb) {
  if (lb_sanity_checks(lb)) {
    return;
  }
  auto const agrid = lb.get_agrid();

  // Distribute summed-up forces from physical particles to ghosts
  cell_structure.ghosts_reset_forces();
  cell_structure.update_ghosts_and_resort_particle(Cells::DATA_PART_FORCE);

  // Keep track of ghost particles (ids) that have already been coupled
  LB::CouplingBookkeeping bookkeeping{cell_structure};
  // Apply particle forces to the LB fluid at particle positions.
  // For physical particles, also set particle velocity = fluid velocity.
  for (auto const &particle_range :
       {cell_structure.local_particles(), cell_structure.ghost_particles()}) {
    for (auto const &p : particle_range) {
      if (!LB::is_tracer(p))
        continue;
      if (bookkeeping.should_be_coupled(p)) {
        for (auto const &pos :
             positions_in_halo(p.pos(), box_geo, local_box, agrid)) {
          lb.add_force_density(pos, p.force());
        }
      }
    }
  }

  // Clear ghost forces to avoid double counting later
  cell_structure.ghosts_reset_forces();
}

void lb_tracers_propagate(CellStructure &cell_structure,
                          BoxGeometry const &box_geo, LB::Solver const &lb,
                          double time_step) {
  if (lb_sanity_checks(lb)) {
    return;
  }
  /* Same budget as @ref CellStructure::check_resort_required: both partners of
   * a pair move, so each may travel at most half the skin before the Verlet
   * list can go stale.  The displacement is measured with the (Lees-Edwards
   * aware) minimum image so that a tracer repositioned across a periodic or
   * shear boundary is not mistaken for a box-length jump. */
  auto const lim = Utils::sqr(cell_structure.get_verlet_skin() / 2.);

  // Advect particles
  for (auto &p : cell_structure.local_particles()) {
    if (!LB::is_tracer(p))
      continue;
    p.v() = lb.get_coupling_interpolated_velocity(p.pos());
    for (auto i = 0u; i < 3u; i++) {
      if (!p.is_fixed_along(i)) {
        p.pos()[i] += p.v()[i] * time_step;
      }
    }
    // Verlet list update check
    if (box_geo.get_mi_dist2(p.pos(), p.pos_at_last_verlet_update()) > lim) {
      cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
    }
  }
}
#endif // ESPRESSO_VIRTUAL_SITES_INERTIALESS_TRACERS
