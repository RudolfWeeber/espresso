/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#pragma once

/**
 * @file
 * Color-gradient (two-component) LB extension interface.
 * Inherits the common fluid interface and adds model-specific
 * configuration, initialization, and per-component accessors.
 */

#include "LBWalberlaBase.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <optional>
#include <vector>

class LBWalberlaColorGradientBase : public virtual LBWalberlaBase {
public:
  ~LBWalberlaColorGradientBase() override = default;

  // CG-specific configuration
  virtual void set_collision_model_color_gradient(double sigma,
                                                  double beta) = 0;
  virtual void init_pdfs_from_components() = 0;

  // CG-specific observables / forces
  virtual std::vector<Utils::Vector3d>
  get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &positions) = 0;
  virtual void
  add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                              std::vector<double> const &delta_mus) = 0;

  // per-component density accessors
  virtual std::optional<std::array<double, 2>>
  get_node_component_densities(Utils::Vector3i const &node,
                               bool consider_ghosts = false) const = 0;
  virtual bool
  set_node_component_densities(Utils::Vector3i const &node,
                               std::array<double, 2> const &rho) = 0;
  virtual std::vector<double>
  get_slice_component_densities(Utils::Vector3i const &lower,
                                Utils::Vector3i const &upper) const = 0;

  // per-component viscosity
  virtual void set_component_viscosities(std::array<double, 2> const &nu) = 0;
  [[nodiscard]] virtual std::array<double, 2>
  get_component_viscosities() const = 0;
};
