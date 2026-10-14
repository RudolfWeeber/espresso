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

#include "exscalicos/Coulomb.hpp"

#include <Exscalicos.hpp>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ExscalicosBridge {

namespace {
using ExecutionSpace = Kokkos::DefaultHostExecutionSpace;
using MemorySpace = Kokkos::HostSpace;
using Data = ExscalicosData<double, ExecutionSpace, MemorySpace>;
using P3M =
    Exscalicos<ExscalicosSolver::P3M, double, ExecutionSpace, MemorySpace>;

void handle_error(ExscalicosReturn status) {
  if (not status.check(false)) {
    throw std::runtime_error("ExScaLiCoS: " + status.errorMsg());
  }
}
} // namespace

struct Coulomb::Implementation {
  Data data;
  P3M solver;
  /** Whether the library was tuned at least once. */
  bool is_tuned = false;
  /** Parameters provided by the user. */
  Parameters parameters;
};

std::vector<std::string> Coulomb::available_methods() { return {"p3m"}; }

Coulomb::Coulomb(MPI_Comm comm, std::string method,
                 Parameters const &parameters)
    : m_impl{std::make_unique<Implementation>()}, m_method{std::move(method)},
      m_delegate_near_field{true}, m_r_cut{parameters.r_cut},
      m_alpha{parameters.alpha} {
  auto const methods = available_methods();
  if (std::ranges::find(methods, m_method) == methods.end()) {
    throw std::invalid_argument("Method '" + m_method +
                                "' is unknown or not supported in ExScaLiCoS;"
                                " currently available methods are 'p3m'");
  }
  auto &data = m_impl->data;
  auto &solver = m_impl->solver;
  m_impl->parameters = parameters;
  data.globalComm(comm);
  handle_error(solver.create(data));
  // all ranks take part in both the short-range and long-range calculations
  auto const last_rank = data.globalSize() - 1;
  handle_error(solver.setParameter(ExscalicosP3MShortRangeProcess{}, data, 0,
                                   last_rank));
  handle_error(
      solver.setParameter(ExscalicosP3MLongRangeProcess{}, data, 0, last_rank));
  if (parameters.r_cut > 0.) {
    auto r_cut = parameters.r_cut;
    handle_error(solver.setRMax(r_cut));
  }
  if (parameters.alpha > 0.) {
    auto alpha = parameters.alpha;
    handle_error(solver.setAlpha(alpha));
  }
  if (parameters.cao > 0) {
    auto cao = parameters.cao;
    handle_error(solver.setCAO(cao));
  }
  if (std::ranges::all_of(parameters.mesh, [](int n) { return n > 0; })) {
    auto mesh = parameters.mesh;
    handle_error(solver.setNGrid(std::span<int, 3>(mesh)));
  }
}

Coulomb::~Coulomb() { m_impl->solver.destroy(m_impl->data); }

Parameters Coulomb::get_parameters() const {
  auto result = m_impl->parameters;
  if (m_impl->is_tuned) {
    auto const &solver = m_impl->solver;
    result.r_cut = solver.getRMax();
    result.alpha = solver.getAlpha();
    result.cao = solver.getCAO();
    result.mesh = solver.getNGrid();
  }
  return result;
}

void Coulomb::set_runtime_parameters(double const *box_l,
                                     int const *periodicity,
                                     int total_particles) {
  auto &data = m_impl->data;
  for (int i = 0; i < 3; ++i) {
    data.sysSize(i, box_l[i]);
    data.periodicity(i, periodicity[i]);
  }
  data.nSysParticles(static_cast<long>(total_particles));
}

void Coulomb::set_accuracy(double accuracy) {
  m_impl->parameters.accuracy = accuracy;
  m_impl->data.accuracy(accuracy);
}

void Coulomb::set_r_cut(double r_cut) {
  handle_error(m_impl->solver.setRMax(r_cut));
  m_r_cut = r_cut;
}

void Coulomb::update_cached_parameters() {
  m_r_cut = m_impl->solver.getRMax();
  m_alpha = m_impl->solver.getAlpha();
}

void Coulomb::tune(std::vector<double> &charges,
                   std::vector<double> &positions) {
  auto &solver = m_impl->solver;
  handle_error(solver.tune(m_impl->data, positions, charges));
  m_impl->is_tuned = true;
  update_cached_parameters();
  // The library re-tunes when the box volume changes; keep the
  // cutoff fixed, since it determines the ESPResSo cell system.
  handle_error(solver.setRMax(m_r_cut));
}

void Coulomb::run(std::vector<double> &charges, std::vector<double> &positions,
                  std::vector<double> &fields,
                  std::vector<double> &potentials) {
  auto &solver = m_impl->solver;
  auto const n_part = charges.size();
  fields.resize(3ul * n_part);
  potentials.resize(n_part);
  handle_error(solver.setNearFieldFlag(not m_delegate_near_field));
  handle_error(
      solver.run(m_impl->data, positions, charges, potentials, fields));
  update_cached_parameters();
}

} // namespace ExscalicosBridge
