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

/**
 * @file
 * Bridge to the ExScaLiCoS library. The library is header-only and
 * relies on Kokkos, Cabana and heFFTe; this bridge hides these
 * dependencies from the ESPResSo core with the pointer-to-implementation
 * idiom. Only the short-range kernels are exposed, since they are
 * evaluated in the ESPResSo near-field loop.
 */

#include <common/ExscalicosShortRange.hpp>

#include <mpi.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace ExscalicosBridge {

/**
 * @brief Parameters of the ExScaLiCoS P3M method.
 * Non-positive values are tuned by the library.
 */
struct Parameters {
  /** Target RMS error of the fields of a unit charge. */
  double accuracy = -1.;
  /** Real-space cutoff. */
  double r_cut = -1.;
  /** Ewald splitting parameter. */
  double alpha = -1.;
  /** Charge assignment order. */
  int cao = -1;
  /** Number of mesh points. */
  std::array<int, 3> mesh = {-1, -1, -1};
};

/** @brief Abstraction of a Coulomb method from the ExScaLiCoS library. */
class Coulomb {
public:
  Coulomb(MPI_Comm comm, std::string method, Parameters const &parameters);
  ~Coulomb();

  /** @brief Methods available in the library. */
  static std::vector<std::string> available_methods();

  std::string get_method() const { return m_method; }

  /** @brief Get the current, possibly tuned, parameters. */
  Parameters get_parameters() const;

  /** @brief Set box geometry and number of particles. */
  void set_runtime_parameters(double const *box_l, int const *periodicity,
                              int total_particles);

  /** @brief Set the target RMS error of the fields of a unit charge. */
  void set_accuracy(double accuracy);

  /** @brief Tune the parameters which were not set by the user. */
  void tune(std::vector<double> &charges, std::vector<double> &positions);

  /**
   * @brief Calculate the fields and potentials for charges.
   * @param with_potentials  Calculate the potentials, otherwise only the
   *                         fields are valid.
   */
  void run(std::vector<double> &charges, std::vector<double> &positions,
           std::vector<double> &fields, std::vector<double> &potentials,
           bool with_potentials);

  /**
   * @brief Delegate the short-range calculation.
   * By default, ESPResSo calculates the short-range forces and energies.
   * This decision can be overridden to obtain the result from a full
   * ExScaLiCoS calculation.
   * @param delegate  Delegate short-range calculation to ESPResSo if true,
   *                  or to ExScaLiCoS if false.
   */
  void set_near_field_delegation(bool delegate) {
    m_delegate_near_field = delegate;
  }

  bool get_near_field_delegation() const { return m_delegate_near_field; }

  /** @brief Get the real-space cutoff. */
  double r_cut() const { return m_r_cut; }

  /** @brief Set the real-space cutoff. Takes effect on the next tuning. */
  void set_r_cut(double r_cut);

  /**
   * @brief Provide the times per step of the near-field calculation
   * done by ESPResSo and of the library run, for balancing the
   * real-space cutoff. The library time includes all communication.
   * @param near_field_seconds  Near-field time in seconds
   * @param library_seconds     Library time in seconds
   */
  void set_step_times(double near_field_seconds, double library_seconds);

  /**
   * @brief Balance the near-field time and the library time of the
   * last run by adjusting the real-space cutoff, and re-tune the
   * remaining parameters.
   * @return true if the parameters changed, false if already balanced
   */
  bool balance_r_cut(std::vector<double> &charges,
                     std::vector<double> &positions);

  /**
   * @brief Calculate the short-range pair field magnitude
   * of a unit charge, divided by the distance.
   */
  double pair_force(double dist) const { return m_kernel.field(dist); }

  /** @brief Calculate the short-range pair potential of a unit charge. */
  double pair_energy(double dist) const { return m_kernel.potential(dist); }

private:
  struct Implementation;
  std::unique_ptr<Implementation> m_impl;
  std::string m_method;
  bool m_delegate_near_field;
  /** Cached real-space cutoff, used in the short-range kernels. */
  double m_r_cut;
  /** Short-range kernels with the cached splitting parameter. */
  ExscalicosEwaldShortRange<double> m_kernel;

  void update_cached_parameters();
};

} // namespace ExscalicosBridge
