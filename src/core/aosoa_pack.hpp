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
  // phase 5: id/mass join the store-aliased views — they alias the
  // authoritative ParticleStore id/mass columns and are indexed by *store row*
  // via row(i). They are read only on cold paths (bond kernels), so the strided
  // store gather is acceptable. The rebuild-time id/mass writes in
  // commit_particle DIE.
  //
  // phase-5 perf recovery: `type` is once again a PACK-OWNED contiguous array,
  // written at pack-rebuild time in commit_particle (`type(index)=p.type()`)
  // and read pack-indexed by the hot pair kernels
  // (forces/energy/pressure_cabana). The ParticleStore type column REMAINS
  // authoritative; this pack copy is a derived cache refreshed on rebuild (a
  // mid-run type change forces a rebuild via on_particle_type_change ->
  // set_resort_particles).
  //
  // Likewise `charge`/`dipm` are pack-owned contiguous arrays, refreshed per
  // step ONLY when a coulomb (resp. dipolar) actor is active (see
  // refresh_pack_charges / refresh_pack_dipm). The hot pair kernels and the P3M
  // gather/spread loops read them PACK-INDEXED (contiguous). Pure-LJ runs never
  // touch them, paying zero cost. The store q/dipm columns stay authoritative.
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
  IdViewType id;
  MassViewType mass;
  // `charge` aliases the authoritative ParticleStore q column and is read by
  // *store row* on the cold BOND path (BondedCoulomb needs charges even with NO
  // coulomb solver, so it must always be valid — hence it aliases the store
  // rather than the guarded pack-owned pair_charge below).
  ChargeViewType charge;
  // Pack-index -> store-row translation. Identity on the local prefix.
  RowMapViewType row_map;

  // Pack-owned columns (indexed by pack index).
  // `type` is written on rebuild (read by the hot pair kernels). `pair_charge`/
  // `pair_dipm` are the contiguous hot-path charge/dipm columns, refreshed per
  // step ONLY when the respective solver is active (see refresh_pack_charges /
  // refresh_pack_dipm); the hot pair kernels and P3M gather/spread read them
  // pack-indexed. `flags` is written every commit and is the allocation
  // sentinel (type/pair_charge/pair_dipm/flags are resized together).
  TypeViewType type;
  ChargeViewType pair_charge;
  DipmViewType pair_dipm;
  FlagsViewType flags;

  AoSoA_pack() = default;

  AoSoA_pack(std::size_t num_particles) { resize(num_particles); }

  /** @brief Translate a pack index to its ParticleStore row. */
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION int
  row(std::size_t i) const {
    return row_map(i);
  }

  void resize(std::size_t num_particles) {
    // phase 5: id/mass/charge are store-aliased (bound in
    // bind_pack_store_views), not allocated here.
    // `type`/`pair_charge`/`pair_dipm`/`flags` are pack-owned contiguous
    // columns allocated here; `flags` serves as the allocation sentinel (all
    // four are resized together).
    if (flags.extent(0) == 0) {
      // First allocation
      type =
          TypeViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "type"),
                       num_particles);
      pair_charge = ChargeViewType(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair_charge"),
          num_particles);
      pair_dipm = DipmViewType(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair_dipm"),
          num_particles);
      flags = FlagsViewType("flags", num_particles);
    } else {
      // Reallocation
      Kokkos::realloc(Kokkos::WithoutInitializing, type, num_particles);
      Kokkos::realloc(Kokkos::WithoutInitializing, pair_charge, num_particles);
      Kokkos::realloc(Kokkos::WithoutInitializing, pair_dipm, num_particles);
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
