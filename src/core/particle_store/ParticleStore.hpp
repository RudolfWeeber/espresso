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

#pragma once

#include <config/config.hpp>

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>

class Particle; // attach_to_store is defined in Particle.hpp

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns per-particle quantities in a single index space: local particles
 * first (cell-traversal order), ghosts appended. Vector/quaternion fields are
 * particle-major (@ref StateVectorLayout, LayoutRight) Kokkos columns; see
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 *
 * Migration phase 2: force and torque columns (observables). Phase 3 adds the
 * STATE columns (position, image box, quaternion, position-at-last-verlet-
 * update, position-at-last-time-step, Lees-Edwards offset and flag). Rebuild
 * protocol: mark_dirty() on any topology change; the owner (CellStructure)
 * later runs begin_rebuild / assign_row-per-particle / finish_rebuild.
 * Rebuild preserves values by old row and seeds new rows from the particle's
 * migration carrier (defaults for genuinely new particles).
 * Rebuilds are purely rank-local (no MPI).
 *
 * Since the phase-3 flip, state and observable columns are the source of
 * truth for attached particles; detached particles live in their carriers.
 */
class ParticleStore {
public:
  /**
   * @brief Memory layout of the vector/quaternion state columns.
   *
   * Particle-major (LayoutRight): one particle's components are contiguous.
   *
   * DECISION (phase 3.5): interleaved same-state A/B measurements (with and
   * without AVX2) showed component-major (LayoutLeft) columns cost +10-12% on
   * gather-dominated CPU paths -- pair-kernel neighbor gathers and ghost
   * packing touch 3 cache lines per particle-vector under LayoutLeft instead
   * of 1, with zero SIMD offset (vector_length == 1 on the host). The host
   * columns are therefore particle-major. Component-major is re-evaluated at
   * phase 8 (GPU), where coalesced per-component access may flip the balance.
   */
  using StateVectorLayout = Kokkos::LayoutRight;
  using Column = Kokkos::DualView<double *[3], StateVectorLayout>;
  using IntegerColumn = Kokkos::DualView<int *[3], StateVectorLayout>;
  using QuaternionColumn = Kokkos::DualView<double *[4], StateVectorLayout>;
  using ScalarColumn = Kokkos::DualView<double *>;
  using ShortColumn = Kokkos::DualView<short *>;

  std::size_t number_of_local_particles() const {
    return m_number_of_local_particles;
  }
  std::size_t number_of_ghost_particles() const {
    return m_number_of_ghost_particles;
  }
  std::size_t number_of_particles() const {
    return m_number_of_local_particles + m_number_of_ghost_particles;
  }

  void mark_dirty() { m_dirty = true; }
  /** @brief Monotonic row-assignment generation; bumped by finish_rebuild.
   *  Consumers may cache row-derived data (e.g. ghost-communication row
   *  ranges) tagged with this value and reuse it while it is unchanged. */
  std::uint64_t generation() const { return m_generation; }
  bool is_dirty() const { return m_dirty; }

  void begin_rebuild(std::size_t number_of_local_particles,
                     std::size_t number_of_ghost_particles);
  void assign_row(Particle &particle, int row);
  void finish_rebuild();

  /** Release all Kokkos-backed columns. Must be called while the Kokkos
   *  runtime is still alive (e.g. before Kokkos::finalize); afterwards the
   *  store is empty and dirty. */
  void release_columns() {
    m_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_torque = Column{};
#endif
    m_position = Column{};
    m_image_box = IntegerColumn{};
#ifdef ESPRESSO_ROTATION
    m_quaternion = QuaternionColumn{};
#endif
    m_position_at_last_verlet_update = Column{};
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_position_last_time_step = Column{};
#endif
    m_lees_edwards_offset = ScalarColumn{};
    m_lees_edwards_flag = ShortColumn{};
    m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_old_torque = Column{};
#endif
    m_old_position = Column{};
    m_old_image_box = IntegerColumn{};
#ifdef ESPRESSO_ROTATION
    m_old_quaternion = QuaternionColumn{};
#endif
    m_old_position_at_last_verlet_update = Column{};
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_old_position_last_time_step = Column{};
#endif
    m_old_lees_edwards_offset = ScalarColumn{};
    m_old_lees_edwards_flag = ShortColumn{};
    m_number_of_local_particles = 0u;
    m_number_of_ghost_particles = 0u;
    m_old_number_of_particles = 0u;
    m_dirty = true;
  }

  // -- observable columns (phase 2) -----------------------------------------
  VectorReference force_reference(int const row) {
    return column_reference(m_force, row);
  }
  Utils::Vector3d force_value(int const row) const {
    return column_value(m_force, row);
  }
#ifdef ESPRESSO_ROTATION
  VectorReference torque_reference(int const row) {
    return column_reference(m_torque, row);
  }
  Utils::Vector3d torque_value(int const row) const {
    return column_value(m_torque, row);
  }
#endif

  // -- state columns (phase 3) ----------------------------------------------
  VectorReference position_reference(int const row) {
    return column_reference(m_position, row);
  }
  Utils::Vector3d position_value(int const row) const {
    return column_value(m_position, row);
  }
  IntegerVectorReference image_box_reference(int const row) {
    return integer_column_reference(m_image_box, row);
  }
  Utils::Vector3i image_box_value(int const row) const {
    return integer_column_value(m_image_box, row);
  }
#ifdef ESPRESSO_ROTATION
  QuaternionReference quaternion_reference(int const row) {
    return quaternion_column_reference(m_quaternion, row);
  }
  Utils::Quaternion<double> quaternion_value(int const row) const {
    return quaternion_column_value(m_quaternion, row);
  }
#endif
  VectorReference position_at_last_verlet_update_reference(int const row) {
    return column_reference(m_position_at_last_verlet_update, row);
  }
  Utils::Vector3d position_at_last_verlet_update_value(int const row) const {
    return column_value(m_position_at_last_verlet_update, row);
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  VectorReference position_last_time_step_reference(int const row) {
    return column_reference(m_position_last_time_step, row);
  }
  Utils::Vector3d position_last_time_step_value(int const row) const {
    return column_value(m_position_last_time_step, row);
  }
#endif
  double &lees_edwards_offset(int const row) {
    return scalar_reference(m_lees_edwards_offset, row);
  }
  double lees_edwards_offset(int const row) const {
    return scalar_value(m_lees_edwards_offset, row);
  }
  short &lees_edwards_flag(int const row) {
    return scalar_reference(m_lees_edwards_flag, row);
  }
  short lees_edwards_flag(int const row) const {
    return scalar_value(m_lees_edwards_flag, row);
  }

  // -- host-view getters for kernel wiring (phase 7) ------------------------
  auto force_view() { return m_force.view_host(); }
#ifdef ESPRESSO_ROTATION
  auto torque_view() { return m_torque.view_host(); }
#endif
  auto position_view() { return m_position.view_host(); }
  auto image_box_view() { return m_image_box.view_host(); }
#ifdef ESPRESSO_ROTATION
  auto quaternion_view() { return m_quaternion.view_host(); }
#endif
  auto position_at_last_verlet_update_view() {
    return m_position_at_last_verlet_update.view_host();
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto position_last_time_step_view() {
    return m_position_last_time_step.view_host();
  }
#endif
  auto lees_edwards_offset_view() { return m_lees_edwards_offset.view_host(); }
  auto lees_edwards_flag_view() { return m_lees_edwards_flag.view_host(); }

private:
  // -- proxy factories for the various column kinds -------------------------
  VectorReference column_reference(Column &column, int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto &view = column.view_host();
    // Layout-agnostic proxy over one row: address of component 0 plus the
    // component-to-component stride. stride(1) is 1 for LayoutRight (row
    // contiguous) and number-of-rows for LayoutLeft; &view(row, 0) is the
    // correct base under either layout, unlike view.data() + row (which
    // assumes LayoutLeft). stride(1) is the non-deprecated spelling of
    // stride_1() in the installed Kokkos version.
    return VectorReference(&view(row, 0), view.stride(1));
  }
  Utils::Vector3d column_value(Column const &column, int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }
  IntegerVectorReference integer_column_reference(IntegerColumn &column,
                                                  int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto &view = column.view_host();
    return IntegerVectorReference(&view(row, 0), view.stride(1));
  }
  Utils::Vector3i integer_column_value(IntegerColumn const &column,
                                       int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }
  QuaternionReference quaternion_column_reference(QuaternionColumn &column,
                                                  int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto &view = column.view_host();
    return QuaternionReference(&view(row, 0), view.stride(1));
  }
  Utils::Quaternion<double>
  quaternion_column_value(QuaternionColumn const &column, int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    Utils::Quaternion<double> q;
    q[0] = view(row, 0);
    q[1] = view(row, 1);
    q[2] = view(row, 2);
    q[3] = view(row, 3);
    return q;
  }
  template <class ScalarColumnType>
  auto scalar_reference(ScalarColumnType &column, int const row)
      -> decltype(column.view_host()(row)) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return column.view_host()(row);
  }
  template <class ScalarColumnType>
  auto scalar_value(ScalarColumnType const &column, int const row) const
      -> std::remove_reference_t<decltype(column.view_host()(row))> {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return column.view_host()(row);
  }

  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
  bool m_dirty = false;
  std::uint64_t m_generation = 0u;

  // -- current-generation columns -------------------------------------------
  Column m_force;
#ifdef ESPRESSO_ROTATION
  Column m_torque;
#endif
  Column m_position;
  IntegerColumn m_image_box;
#ifdef ESPRESSO_ROTATION
  QuaternionColumn m_quaternion;
#endif
  Column m_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_position_last_time_step;
#endif
  ScalarColumn m_lees_edwards_offset;
  ShortColumn m_lees_edwards_flag;

  // -- spare (previous-generation) columns ----------------------------------
  // Capacity-cached double buffering (phase 3.5): these are kept alive across
  // rebuilds as the swap-in write target. During a rebuild they hold the
  // just-current data (the preserve_or_seed READ source); after finish_rebuild
  // they hold the retired generation, reused as the spare on the next swap.
  Column m_old_force;
#ifdef ESPRESSO_ROTATION
  Column m_old_torque;
#endif
  Column m_old_position;
  IntegerColumn m_old_image_box;
#ifdef ESPRESSO_ROTATION
  QuaternionColumn m_old_quaternion;
#endif
  Column m_old_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_old_position_last_time_step;
#endif
  ScalarColumn m_old_lees_edwards_offset;
  ShortColumn m_old_lees_edwards_flag;

  std::size_t m_old_number_of_particles = 0u;
};
