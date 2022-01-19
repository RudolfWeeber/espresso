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
#ifndef CORE_LB_INTERFACE
#define CORE_LB_INTERFACE

#include "config.hpp"

#include <utils/Vector.hpp>

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace LB {
int get_steps_per_md_step(double md_timestep);
}

/** @brief LB implementation currently active. */
enum class ActiveLB : int { NONE, WALBERLA };

/** @brief Switch determining the type of lattice dynamics. */
extern ActiveLB lattice_switch;

struct NoLBActive : public std::exception {
  const char *what() const noexcept override { return "LB not activated"; }
};

/**
 * @brief Propagate the LB fluid.
 */
void lb_lbfluid_propagate();

/**
 * @brief Perform a full initialization of the lattice-Boltzmann system.
 * All derived parameters and the fluid are reset to their default values.
 */
void lb_lbfluid_init();

/**
 * @brief Get the global variable @ref lattice_switch.
 */
ActiveLB lb_lbfluid_get_lattice_switch();

/**
 * @brief Check if tau is an integer multiple of time_step, throws if not
 */
void check_tau_time_step_consistency(double tau, double time_step);

/**
 * @brief Perform LB parameter and boundary velocity checks.
 */
void lb_lbfluid_sanity_checks(double time_step);

/**
 * @brief Set the LB fluid velocity for a single node.
 */
void lb_lbnode_set_velocity(const Utils::Vector3i &ind,
                            const Utils::Vector3d &u);

/**
 * @brief Get the LB time step.
 */
double lb_lbfluid_get_tau();

/**
 * @brief Get the LB grid spacing.
 */
double lb_lbfluid_get_agrid();

/**
 * @brief Get the thermal energy parameter of the LB fluid.
 */
double lb_lbfluid_get_kT();

/**
 * @brief Get the lattice speed (agrid/tau).
 */
double lb_lbfluid_get_lattice_speed();

/** @brief Calculate the average pressure tensor of all nodes by accumulating
 *  over all nodes and dividing by the number of nodes.
 *  Returns the lower triangle of the LB pressure tensor.
 */
const Utils::VectorXd<9> lb_lbfluid_get_pressure_tensor();

Utils::Vector3d lb_lbfluid_calc_fluid_momentum();

/**
 * @brief Calculates the interpolated fluid velocity on the head node process.
 * @param pos Position at which the velocity is to be calculated.
 * @retval interpolated fluid velocity.
 */
const Utils::Vector3d
lb_lbfluid_get_interpolated_velocity(const Utils::Vector3d &pos);

/**
 * @brief Calculates the interpolated fluid density on the head node process.
 * @param pos Position at which the density is to be calculated.
 * @retval interpolated fluid density.
 */
double lb_lbfluid_get_interpolated_density(const Utils::Vector3d &pos);

namespace Walberla {
/**
 * @brief Access the ID of the velocity field
 * @return velocity field id
 */
std::size_t get_velocity_field_id();

/**
 * @brief Access the ID of the force field
 * @return force field id
 */
std::size_t get_force_field_id();
} // namespace Walberla
#endif
