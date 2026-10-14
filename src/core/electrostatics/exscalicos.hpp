/*
 * Copyright (C) 2026 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ESPRESSO_EXSCALICOS

/**
 * @file
 * @ref CoulombExscalicos provides the public interface of the
 * ExScaLiCoS bridge. It relies on type erasure to hide the
 * ExScaLiCoS implementation details from the ESPResSo core.
 * It is implemented by @ref CoulombExscalicosImpl.
 */

#include "electrostatics/actor.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Parameters of an ExScaLiCoS method.
 * Non-positive values are tuned by the library.
 */
struct ExscalicosParameters {
  /** Target RMS error of the forces. */
  double accuracy = -1.;
  /** Real-space cutoff. */
  double r_cut = -1.;
  /** Ewald splitting parameter. */
  double alpha = -1.;
  /** Charge assignment order. */
  int cao = -1;
  /** Number of mesh points. */
  Utils::Vector3i mesh = {-1, -1, -1};
};

struct CoulombExscalicos : public Coulomb::Actor<CoulombExscalicos> {
  virtual ~CoulombExscalicos() = default;

  void on_activation() {
    sanity_checks();
    update_system_params();
    tune();
  }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { update_system_params(); }
  void on_node_grid_change() const {}
  void on_periodicity_change() {
    sanity_checks_periodicity();
    update_system_params();
  }
  void on_cell_structure_change() const {}
  void init() const {}

  void sanity_checks() const {
    sanity_checks_periodicity();
    sanity_checks_charge_neutrality();
  }
  void sanity_checks_periodicity() const;

  bool is_tuned() const { return m_is_tuned; }
  void tune() {
    if (not is_tuned()) {
      tune_impl();
      m_is_tuned = true;
    }
  }

  /** @brief Collect particle data in continuous arrays. */
  virtual void update_particle_data() = 0;
  /** @brief Write forces back to particles. */
  virtual void update_particle_forces() const = 0;
  /** @brief Calculate long-range part of the energy. */
  virtual double long_range_energy() = 0;
  /** @brief Add long-range part of the forces to particles. */
  virtual void add_long_range_forces() = 0;
  /** @brief Reinitialize box shape and periodicity. */
  virtual void update_system_params() = 0;
  virtual std::string get_method() const = 0;
  /** @brief Get the current, possibly tuned, parameters. */
  virtual ExscalicosParameters get_parameters() const = 0;
  /** @brief Short-range cutoff of the ESPResSo kernels. */
  virtual double get_r_cut() const = 0;
  virtual double get_pair_force(double dist) const = 0;
  virtual double get_pair_energy(double dist) const = 0;
  virtual void set_near_field_delegation(bool delegate) = 0;
  virtual bool get_near_field_delegation() const = 0;

  static std::vector<std::string> available_methods();

  /** @brief Calculate near-field pair force. */
  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const {
    if (dist > get_r_cut())
      return {};

    return d * (q1q2 * prefactor * get_pair_force(dist));
  }

  /** @brief Calculate near-field pair energy. */
  double pair_energy(double q1q2, double dist) const {
    if (dist > get_r_cut())
      return 0.;

    return q1q2 * prefactor * get_pair_energy(dist);
  }

protected:
  virtual void tune_impl() = 0;
  bool m_is_tuned = false;
};

std::shared_ptr<CoulombExscalicos>
make_coulomb_exscalicos(std::string const &method,
                        ExscalicosParameters const &parameters);

#endif // ESPRESSO_EXSCALICOS
