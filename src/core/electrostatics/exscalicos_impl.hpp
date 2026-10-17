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

#include "electrostatics/exscalicos.hpp"

#include "communication.hpp"

#include <exscalicos/Coulomb.hpp>

#include <mpi.h>

#include <boost/range/numeric.hpp>

#include <memory>
#include <string>
#include <vector>

struct CoulombExscalicosImpl : public CoulombExscalicos {
  CoulombExscalicosImpl(MPI_Comm comm, std::string const &method,
                        ExscalicosParameters const &parameters);

  void update_particle_data() override;
  void update_particle_forces() const override;
  void update_system_params() override;

  double long_range_energy() override {
    if (not run_library(true)) {
      return 0.;
    }
    return 0.5 * prefactor * boost::inner_product(charges, potentials, 0.);
  }

  void add_long_range_forces() override {
    if (run_library(false)) {
      update_particle_forces();
    }
  }

  std::string get_method() const override { return m_bridge.get_method(); }

  ExscalicosParameters get_parameters() const override;

  double get_r_cut() const override {
    return m_bridge.get_near_field_delegation() ? m_bridge.r_cut() : 0.;
  }

  double get_pair_force(double dist) const override {
    return m_bridge.pair_force(dist);
  }

  double get_pair_energy(double dist) const override {
    return m_bridge.pair_energy(dist);
  }

  void set_near_field_delegation(bool delegate) override;

  bool get_near_field_delegation() const override {
    return m_bridge.get_near_field_delegation();
  }

private:
  /** Keep Kokkos alive until the library data structures are destroyed. */
  std::shared_ptr<KokkosHandle> m_kokkos_handle;
  ExscalicosBridge::Coulomb m_bridge;
  /** Target RMS error of the forces. */
  double m_accuracy;
  /** Whether the real-space cutoff was provided by the user. */
  bool m_fixed_r_cut;
  /** Whether the library parameters were tuned. */
  bool m_library_is_tuned = false;
  /** Accumulated wall time of the library runs, in seconds. */
  double m_library_time = 0.;
  /** Number of library runs in @ref m_library_time. */
  int m_library_runs = 0;
  /** Inputs of the charged particles */
  std::vector<double> positions, charges;
  /** Outputs of the charged particles */
  std::vector<double> fields, potentials;

  /** Largest real-space cutoff supported by the cell system. */
  double max_r_cut() const;
  /**
   * @brief Tune the library parameters.
   * @return false if there are no charged particles, true otherwise
   */
  bool tune_library();
  /**
   * @brief Balance the real-space cutoff between the ESPResSo near field
   * and the library.
   * The near-field time is the time of a force calculation minus the time
   * of the library run. It thus contains all short-range interactions,
   * as in a production run.
   */
  void balance_r_cut();
  /**
   * @brief Calculate the fields and potentials of the charged particles.
   * @param with_potentials  Calculate the potentials (only needed for energies)
   * @return false if there are no charged particles, true otherwise
   */
  bool run_library(bool with_potentials);
  void tune_impl() override;
};

#endif // ESPRESSO_EXSCALICOS
