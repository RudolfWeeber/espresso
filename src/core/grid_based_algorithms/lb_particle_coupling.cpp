/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "lb_particle_coupling.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"
#include "lb_interpolation.hpp"
#include "random.hpp"

#include <profiler/profiler.hpp>
#include <utils/Counter.hpp>
#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

LB_Particle_Coupling lb_particle_coupling;

void mpi_bcast_lb_particle_coupling_local() {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

REGISTER_CALLBACK(mpi_bcast_lb_particle_coupling_local)

void mpi_bcast_lb_particle_coupling() {
  mpi_call_all(mpi_bcast_lb_particle_coupling_local);
}

void lb_lbcoupling_activate() { lb_particle_coupling.couple_to_md = true; }

void lb_lbcoupling_deactivate() {
  if (lattice_switch != ActiveLB::NONE && this_node == 0 &&
      lb_particle_coupling.gamma > 0.) {
    runtimeWarningMsg()
        << "Recalculating forces, so the LB coupling forces are not "
           "included in the particle force the first time step. This "
           "only matters if it happens frequently during sampling.";
  }

  lb_particle_coupling.couple_to_md = false;
}

void lb_lbcoupling_set_gamma(double gamma) {
  lb_particle_coupling.gamma = gamma;
}

double lb_lbcoupling_get_gamma() { return lb_particle_coupling.gamma; }

bool lb_lbcoupling_is_seed_required() {
  if (lattice_switch == ActiveLB::WALBERLA) {
    return not lb_particle_coupling.rng_counter_coupling.is_initialized();
  }
  return false;
}

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling->value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_coupling_get_rng_state_cpu();
  }
  throw std::runtime_error("No LB active");
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::WALBERLA) {
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
  } else
    throw std::runtime_error("No LB active");
}

/**
 * @brief Add a force to the lattice force density.
 * @param pos Position of the force
 * @param force Force in MD units.
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  /* transform momentum transfer to lattice units
     (eq. (12) @cite ahlrichs99a) */
  auto const delta_j = -(time_step / lb_lbfluid_get_lattice_speed()) * force;
  lb_lbinterpolation_add_force_density(pos, delta_j);
}
/** @brief Calculate particle drift velocity offset due to ENGINE and
 *  ELECTROHYDRODYNAMICS */
Utils::Vector3d lb_particle_coupling_drift_vel_offset(const Particle &p) {
  Utils::Vector3d vel_offset{};
#ifdef ENGINE
  if (p.p.swim.swimming) {
    vel_offset += p.p.swim.v_swim * p.r.calc_director();
  }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  vel_offset += p.p.mu_E;
#endif
  return vel_offset;
}

/** calculate drag force on a single particle
 *
 *  Section II.C. @cite ahlrichs99a
 *
 *  @param[in] p             The coupled particle.
 *  @param vel_offset        Velocity offset to be added to interpolated LB
 * velocity before calculating the force
 *
 *  @return The viscous coupling force
 */
Utils::Vector3d lb_drag_force(Particle const &p,
                              const Utils::Vector3d &vel_offset) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const interpolated_u =
      lb_lbinterpolation_get_interpolated_velocity(p.r.p) *
      lb_lbfluid_get_lattice_speed();

  Utils::Vector3d v_drift = interpolated_u + vel_offset;
  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return -lb_lbcoupling_get_gamma() * (p.m.v - v_drift);
}

using Utils::Vector;
using Utils::Vector3d;
using Utils::Vector3i;

template <class T, size_t N> using Box = std::pair<Vector<T, N>, Vector<T, N>>;

/**
 * @brief Check if a position is in a box.
 *
 * The left boundary belong to the box, the
 * right one does not. Periodic boundaries are
 * not considered.
 *
 * @param pos Position to check
 * @param box Box to check
 *
 * @return True iff the point is inside of the box.
 */
template <class T, size_t N>
bool in_box(Vector<T, N> const &pos, Box<T, N> const &box) {
  return (pos >= box.first) and (pos < box.second);
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param pos Position to check
 * @param local_box Geometry to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
template <class T>
bool in_local_domain(Vector<T, 3> const &pos, LocalBox<T> const &local_box,
                     T const &halo = {}) {
  auto const halo_vec = Vector<T, 3>::broadcast(halo);

  return in_box(
      pos, {local_geo.my_left() - halo_vec, local_geo.my_right() + halo_vec});
}

/**
 * @brief Check if a position is within the local LB domain
 *       plus halo.
 *
 * @param pos Position to check
 *
 * @return True iff the point is inside of the domain.
 */
bool in_local_halo(Vector3d const &pos) {
  auto const halo = 0.5 * lb_lbfluid_get_agrid();

  return in_local_domain(pos, local_geo, halo);
}

#ifdef ENGINE
void add_swimmer_force(Particle const &p) {
  if (p.p.swim.swimming) {
    // calculate source position
    const double direction =
        double(p.p.swim.push_pull) * p.p.swim.dipole_length;
    auto const director = p.r.calc_director();
    auto const source_position = p.r.p + direction * director;

    if (in_local_halo(source_position)) {
      add_md_force(source_position, p.p.swim.f_swim * director);
    }
  }
}
#endif

Utils::Vector3d lb_particle_coupling_noise(bool enabled, int part_id,
                                           const OptionalCounter &rng_counter) {
  if (enabled) {
    return Random::noise_uniform<RNGSalt::PARTICLES>(rng_counter->value(), 0,
                                                     part_id);
  }
  return {};
}

void couple_particle(Particle &p, bool couple_virtual, double noise_amplitude,
                     const OptionalCounter &rng_counter) {

  if (p.p.is_virtual and not couple_virtual)
    return;

  /* Particles within one agrid of the outermost lattice point
   * of the lb domain can contribute forces to the local lb due to
   * interpolation on neighboring LB nodes. If the particle
   * IS in the local domain, we also add the opposing
   * force to the particle. */
  if (in_local_halo(p.r.p)) {
    auto const drag_force =
        lb_drag_force(p, lb_particle_coupling_drift_vel_offset(p));
    auto const random_force =
        noise_amplitude * lb_particle_coupling_noise(noise_amplitude > 0.0,
                                                     p.identity(), rng_counter);
    auto const coupling_force = drag_force + random_force;
    add_md_force(p.r.p, coupling_force);
    if (in_local_domain(p.r.p, local_geo)) {
      /* Particle is in our LB volume, so this node
       * is responsible to adding its force */
      p.f.f += coupling_force;
    }
  }
}

void lb_lbcoupling_calc_particle_lattice_ia(
    bool couple_virtual, const ParticleRange &particles,
    const ParticleRange &more_particles) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch == ActiveLB::WALBERLA) {
    if (lb_particle_coupling.couple_to_md) {
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::quadratic):
        throw std::runtime_error("The non-linear interpolation scheme is not "
                                 "implemented for the CPU LB.");
      case (InterpolationOrder::linear): {
        using Utils::sqr;
        auto const kT = lb_lbfluid_get_kT() * sqr(lb_lbfluid_get_agrid()) /
                        sqr(lb_lbfluid_get_tau());
        /* Eq. (16) @cite ahlrichs99a.
         * The factor 12 comes from the fact that we use random numbers
         * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
         * time_step comes from the discretization.
         */
        auto const noise_amplitude =
            (kT > 0.) ? std::sqrt(12. * 2. * lb_lbcoupling_get_gamma() * kT /
                                  time_step)
                      : 0.0;

        /* Couple particles ranges */
        for (auto &p : particles) {
          couple_particle(p, couple_virtual, noise_amplitude,
                          lb_particle_coupling.rng_counter_coupling);
#ifdef ENGINE
          add_swimmer_force(p);
#endif
        }

        for (auto &p : more_particles) {
          couple_particle(p, couple_virtual, noise_amplitude,
                          lb_particle_coupling.rng_counter_coupling);
#ifdef ENGINE
          add_swimmer_force(p);
#endif
        }

        break;
      }
      }
    }
  }
}

void lb_lbcoupling_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    if (lb_lbfluid_get_kT() > 0.0) {
      if (lattice_switch == ActiveLB::WALBERLA) {
        lb_particle_coupling.rng_counter_coupling->increment();
      }
    }
  }
}
