/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "attributes.hpp"
#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <cstdint>

struct CellStructure::AoSoA_pack {
  // Particle-major (@ref ParticleStore::StateVectorLayout, LayoutRight) to
  // match the ParticleStore columns these views alias / are derived from.
  //
  // phase 3.5: position/image/director are no longer owned here. They alias
  // the authoritative ParticleStore host columns (position/image) and the
  // store-side derived director view. Kernels index them by *store row*, which
  // is obtained by translating a pack index i through row(i) (see
  // pack_index_to_store_row). For i < n_local the translation is the identity
  // (both are built in cell-traversal order); only the deduped ghost tail is
  // remapped. commit_particle no longer copies position/image/director.
  //
  // phase 4: velocity joins the store-aliased views — it aliases the
  // authoritative ParticleStore velocity host column and is likewise indexed by
  // *store row* via row(i). commit_particle no longer copies velocity, and
  // resize() no longer allocates it.
  //
  // phase 5: charge/dipm and id/type/mass ALSO join the store-aliased views —
  // they alias the authoritative ParticleStore scalar columns (q/dipm/id/type/
  // mass) and are likewise indexed by *store row* via row(i). The per-step
  // charge/dipm copies and the rebuild-time id/type/mass writes in
  // commit_particle DIE; commit_particle only sets the pack-owned exclusion
  // flags. `flags` becomes the allocation sentinel (the last pack-owned view).
  //
  // The layout of the store-aliased vector views (position/image/director/
  // velocity) MUST match the store columns' layout.
  using PositionViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using VelocityViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using DirectorViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using ImageViewType = Kokkos::View<int *[3], ParticleStore::StateVectorLayout,
                                     Kokkos::HostSpace>;
  using RowMapViewType = Kokkos::View<int const *, Kokkos::HostSpace>;
  using ChargeViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using DipmViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using IdViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using TypeViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using MassViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using FlagsViewType = Kokkos::View<uint8_t *, Kokkos::HostSpace>;

  // Store-aliased / store-derived columns (indexed by store row).
  PositionViewType position;
  ImageViewType image;
  DirectorViewType director;
  VelocityViewType velocity;
  ChargeViewType charge;
  DipmViewType dipm;
  IdViewType id;
  TypeViewType type;
  MassViewType mass;
  // Pack-index -> store-row translation. Identity on the local prefix.
  RowMapViewType row_map;

  // Pack-owned per-step column (indexed by pack index).
  FlagsViewType flags;

  AoSoA_pack() = default;

  AoSoA_pack(std::size_t num_particles) { resize(num_particles); }

  /** @brief Translate a pack index to its ParticleStore row. */
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION int
  row(std::size_t i) const {
    return row_map(i);
  }

  void resize(std::size_t num_particles) {
    // phase 5: charge/dipm/id/type/mass are store-aliased (bound in
    // bind_pack_store_views), not allocated here; `flags` is the last
    // pack-owned column and serves as the allocation sentinel.
    if (flags.extent(0) == 0) {
      // First allocation
      flags = FlagsViewType("flags", num_particles);
    } else {
      // Reallocation
      Kokkos::realloc(flags, num_particles);
    }
  }

  template <typename array_layout, typename T, std::size_t N>
  Utils::Vector<T, N> get_vector_at(
      Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> const &view,
      std::size_t i) const {
    Utils::Vector<T, N> result;
    auto const data = result.data();
#if !defined(__NVCOMPILER)
#if (defined(__GNUC__) or defined(__GNUG__)) && !defined(__clang__)
#pragma GCC unroll 8
#else
#pragma omp unroll
#endif
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      data[j] = view(i, j);
    }
    return result;
  }

  template <typename array_layout, typename T, std::size_t N>
  void
  set_vector_at(Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> &view,
                std::size_t i, Utils::Vector<T, N> const &value) {
#if !defined(__NVCOMPILER)
#if (defined(__GNUC__) or defined(__GNUG__)) && !defined(__clang__)
#pragma GCC unroll 8
#else
#pragma omp unroll
#endif
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      view(i, j) = value[j];
    }
  }

  void set_has_exclusion(std::size_t i, bool value) {
    flags(i) = value ? uint8_t{1} : uint8_t{0};
  }

  bool has_exclusion(std::size_t i) const { return flags(i) == uint8_t{1}; }
};
