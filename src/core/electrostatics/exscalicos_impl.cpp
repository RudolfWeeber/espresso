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

#include "config/config.hpp"

#ifdef ESPRESSO_EXSCALICOS

#include "electrostatics/exscalicos.hpp"
#include "electrostatics/exscalicos_impl.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <exscalicos/Coulomb.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

std::shared_ptr<CoulombExscalicos>
make_coulomb_exscalicos(std::string const &method,
                        ExscalicosParameters const &parameters) {
  return std::make_shared<CoulombExscalicosImpl>(comm_cart, method, parameters);
}

std::vector<std::string> CoulombExscalicos::available_methods() {
  return ExscalicosBridge::Coulomb::available_methods();
}

void CoulombExscalicos::sanity_checks_periodicity() const {
  auto const &box_geo = *get_system().box_geo;
  if (not box_geo.periodic(0) or not box_geo.periodic(1) or
      not box_geo.periodic(2)) {
    throw std::runtime_error(
        "CoulombExscalicos: requires periodicity (True, True, True)");
  }
}

static auto to_bridge_parameters(ExscalicosParameters const &parameters) {
  ExscalicosBridge::Parameters result;
  result.r_cut = parameters.r_cut;
  result.alpha = parameters.alpha;
  result.cao = parameters.cao;
  std::ranges::copy(parameters.mesh, result.mesh.begin());
  return result;
}

CoulombExscalicosImpl::CoulombExscalicosImpl(
    MPI_Comm comm, std::string const &method,
    ExscalicosParameters const &parameters)
    : m_kokkos_handle{::kokkos_handle},
      m_bridge{comm, method, to_bridge_parameters(parameters)},
      m_accuracy{parameters.accuracy}, m_fixed_r_cut{parameters.r_cut > 0.} {
  if (m_accuracy <= 0.) {
    throw std::domain_error("Parameter 'accuracy' must be > 0");
  }
}

ExscalicosParameters CoulombExscalicosImpl::get_parameters() const {
  auto const parameters = m_bridge.get_parameters();
  ExscalicosParameters result;
  result.accuracy = m_accuracy;
  result.r_cut = parameters.r_cut;
  result.alpha = parameters.alpha;
  result.cao = parameters.cao;
  std::ranges::copy(parameters.mesh, result.mesh.begin());
  return result;
}

void CoulombExscalicosImpl::update_particle_data() {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;

  positions.clear();
  charges.clear();

  for (auto const &p : cell_structure.local_particles()) {
    if (p.q() != 0.) {
      auto const pos = box_geo.folded_position(p.pos());
      positions.push_back(pos[0]);
      positions.push_back(pos[1]);
      positions.push_back(pos[2]);
      charges.push_back(p.q());
    }
  }
}

void CoulombExscalicosImpl::update_particle_forces() const {
  if (positions.empty())
    return;

  auto const &cell_structure = *get_system().cell_structure;

  auto it_fields = fields.begin();
  for (auto &p : cell_structure.local_particles()) {
    if (p.q() != 0.) {
      p.force() += prefactor * p.q() *
                   Utils::Vector3d(std::span<const double>(&*it_fields, 3ul));
      std::advance(it_fields, 3);
    }
  }

  /* Check that the particle number did not change */
  assert(it_fields == fields.end());
}

void CoulombExscalicosImpl::update_system_params() {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;
  auto const periodicity =
      Utils::Vector3i{static_cast<int>(box_geo.periodic(0u)),
                      static_cast<int>(box_geo.periodic(1u)),
                      static_cast<int>(box_geo.periodic(2u))};
  auto n_local = std::size_t{0ul};
  for (auto const &p : cell_structure.local_particles()) {
    if (p.q() != 0.) {
      ++n_local;
    }
  }
  auto const n_total =
      boost::mpi::all_reduce(comm_cart, n_local, std::plus<std::size_t>());
  m_bridge.set_runtime_parameters(box_geo.length().data(), periodicity.data(),
                                  static_cast<int>(n_total));
}

double CoulombExscalicosImpl::max_r_cut() const {
  auto const &system = get_system();
  auto const min_box_l = std::ranges::min(system.box_geo->length());
  auto const min_local_box_l = std::ranges::min(system.local_geo->length());
  auto const verlet_skin =
      std::max(0., system.cell_structure->get_verlet_skin());
  /* ExScaLiCoS limits the cutoff to 0.49 times the smallest box length */
  return std::min(min_local_box_l, 0.49 * min_box_l) - verlet_skin;
}

bool CoulombExscalicosImpl::tune_library() {
  update_system_params();
  auto const n_charges = boost::mpi::all_reduce(comm_cart, charges.size(),
                                                std::plus<std::size_t>());
  if (n_charges == 0ul) {
    return false;
  }
  m_bridge.set_accuracy(m_accuracy / prefactor);
  m_bridge.tune(charges, positions);
  if (m_bridge.get_near_field_delegation() and not m_fixed_r_cut) {
    // the cutoff tuned by the library might not fit in the cell system
    auto const r_cut_max = max_r_cut();
    if (m_bridge.r_cut() > r_cut_max) {
      m_bridge.set_r_cut(r_cut_max);
      m_bridge.tune(charges, positions);
    }
  }
  m_library_is_tuned = true;
  return true;
}

void CoulombExscalicosImpl::balance_r_cut() {
  auto &system = get_system();
  auto constexpr max_iterations = 10;
  auto constexpr timings = 10;
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    system.on_coulomb_change();
    m_library_time = 0.;
    m_library_runs = 0;
    auto const force_time = benchmark_integration_step(system, timings) / 1e3;
    auto const library_time = boost::mpi::all_reduce(
        comm_cart, m_library_time / std::max(m_library_runs, 1),
        boost::mpi::maximum<double>());
    auto const near_field_time = force_time - library_time;
    if (near_field_time <= 0. or library_time <= 0.) {
      break;
    }
    m_bridge.set_step_times(near_field_time, library_time);
    if (not m_bridge.balance_r_cut(charges, positions)) {
      break;
    }
    auto const r_cut_max = max_r_cut();
    if (m_bridge.r_cut() > r_cut_max) {
      m_bridge.set_r_cut(r_cut_max);
      m_bridge.tune(charges, positions);
      break;
    }
  }
}

void CoulombExscalicosImpl::tune_impl() {
  update_particle_data();
  if (not tune_library() and not m_fixed_r_cut) {
    throw std::runtime_error(
        "CoulombExscalicos: no charged particles in the system");
  }
  if (m_library_is_tuned and m_bridge.get_near_field_delegation() and
      not m_fixed_r_cut) {
    balance_r_cut();
  }
  // with a fixed cutoff, tuning can be deferred until charges are present,
  // e.g. when the solver is restored from a checkpoint before the particles
  get_system().on_coulomb_change();
}

bool CoulombExscalicosImpl::run_library(bool with_potentials) {
  update_particle_data();
  if (not m_library_is_tuned and not tune_library()) {
    return false;
  }
  auto const tick = MPI_Wtime();
  m_bridge.run(charges, positions, fields, potentials, with_potentials);
  m_library_time += MPI_Wtime() - tick;
  ++m_library_runs;
  return true;
}

void CoulombExscalicosImpl::set_near_field_delegation(bool delegate) {
  m_bridge.set_near_field_delegation(delegate);
  if (is_tuned()) {
    m_is_tuned = false;
    if (not m_system.expired()) {
      tune();
    }
  }
}

#endif // ESPRESSO_EXSCALICOS
