/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#ifndef GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#define GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA
#include <LBWalberlaBase.hpp>
#include <lb_walberla_init.hpp>

struct LBWalberlaParams {
  LBWalberlaParams(double agrid, double tau) : m_agrid(agrid), m_tau(tau) {}
  double get_agrid() const { return m_agrid; };
  double get_tau() const { return m_tau; };

private:
  double m_agrid;
  double m_tau;
};

/** @brief Access the per-MPI-node LBWalberla instance */
LBWalberlaBase *lb_walberla();

/** @brief Access the Walberla parameters */
LBWalberlaParams *lb_walberla_params();

void init_lb_walberla(double viscosity, double density, double agrid,
                      double tau, const Utils::Vector3i &grid_dimensions,
                      const Utils::Vector3i &node_grid, double kT,
                      unsigned int seed);

/** @brief Create the LBWalberla instance and sets the lattice switch to
 *  WALBERLA
 *
 *  @param viscosity Fluid viscosity
 *  @param density   Fluid density
 *  @param agrid     Size of one LB cell
 *  @param tau       LB time step
 *  @param box_size  Box dimensions
 *  @param kT        Temperature
 *  @param seed      LB random seed
 */
void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau, Utils::Vector3d box_size, double kT,
                          unsigned int seed);

/** @brief Destruct the LBWalberla instance and set lattice switch to NONE */
void mpi_destruct_lb_walberla();

#endif // LB_WALBERLA

#endif
