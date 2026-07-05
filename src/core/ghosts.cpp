/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.hpp "ghosts.hpp"
 *
 * Note on variable naming:
 * - a "GhostCommunicator" is always named "gcr",
 * - a "GhostCommunication" is always named "ghost_comm".
 */
#include "ghosts.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/ParticleStore.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/range/numeric.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <vector>

/** Tag for ghosts communications. */
#define REQ_GHOST_SEND 100

/**
 * Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 */
class CommBuf {
public:
  /** Returns a pointer to the non-bond storage.
   */
  char *data() { return buf.data(); }
  const char *data() const { return buf.data(); }

  /** Returns the number of elements in the non-bond storage.
   */
  std::size_t size() const { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(std::size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage.
   */
  auto &bonds() { return bondbuf; }
  const auto &bonds() const { return bondbuf; }

  auto make_span() { return std::span(buf.data(), buf.size()); }

private:
  std::vector<char> buf;     ///< Buffer for everything but bonds
  std::vector<char> bondbuf; ///< Buffer for bond lists
};

/** @brief Pseudo-archive to calculate the size of the serialization buffer. */
class SerializationSizeCalculator {
  std::size_t m_size = 0;

public:
  auto size() const { return m_size; }

  template <class T> auto &operator<<(T &) {
    m_size += sizeof(T);
    return *this;
  }

  template <class T> auto &operator&(T &t) { return *this << t; }
};

/**
 * @brief Raw column pointers + strides for the mixed-parts POSITION fast path.
 *
 * Captured once per communication from the (clean) ParticleStore's host views
 * so the generic per-particle POSITION serialization can read/write the four
 * position columns directly by row instead of through per-field accessor
 * proxies. See make_position_row_context() below for how it is built and for
 * the applicability rationale. A row-r component-c element of a column lives at
 * base + r*row_stride + c*comp_stride (comp_stride==1 and one row contiguous
 * under the store's LayoutRight).
 */
struct PositionRowContext {
  double *pos; ///< base of the position column (row 0, component 0)
  int *img;    ///< base of the image-box column
  std::size_t pos_row_stride;
  std::size_t pos_comp_stride;
  std::size_t img_row_stride;
  std::size_t img_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *quat;
  std::size_t quat_row_stride;
  std::size_t quat_comp_stride;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  double *plast;
  std::size_t plast_row_stride;
  std::size_t plast_comp_stride;
#endif
};

/** @brief Raw-pointer view of the force/torque columns, hoisted once per
 *  communication so the per-part_list bulk loops pay no view-handle
 *  (refcount) or DualView bookkeeping cost per list. */
struct ForceRowContext {
  double *force;
  std::size_t force_row_stride;
  std::size_t force_comp_stride;
#ifdef ESPRESSO_ROTATION
  double *torque;
  std::size_t torque_row_stride;
  std::size_t torque_comp_stride;
#endif
};

/** @brief Type of reduction to carry out during serialization. */
enum class ReductionPolicy {
  /** @brief Reduction for domain-to-domain particle communication. */
  MOVE,
  /** @brief Reduction for cell-to-cell particle update. */
  UPDATE,
};

/** @brief Whether to save the state to or load the state from the archive. */
enum class SerializationDirection { SAVE, LOAD };

/**
 * @brief Serialize particle data, possibly with reduction.
 * The reduction can take place during the save stage, e.g. to apply
 * a ghost shift to the particle position, or during the load stage,
 * e.g. to transfer momentum between particles in two local cells.
 */
template <class Archive>
static void
serialize_and_reduce(Archive &ar, Particle &p, unsigned int data_parts,
                     ReductionPolicy policy, SerializationDirection direction,
                     BoxGeometry const &box_geo,
                     Utils::Vector3d const *ghost_shift,
                     PositionRowContext const *pos_ctx = nullptr) {
  if (data_parts & GHOSTTRANS_PROPRTS) {
    ar & p.id() & p.mol_id() & p.type() & p.propagation();
#ifdef ESPRESSO_ROTATION
    ar & p.rotation();
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    ar & p.rinertia();
#endif
#endif
#ifdef ESPRESSO_MASS
    ar & p.mass();
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    ar & p.q();
#endif
#ifdef ESPRESSO_DIPOLES
    ar & p.dipm();
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    ar & p.mu_E();
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    ar & p.vs_relative();
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    ar & p.gamma();
#ifdef ESPRESSO_ROTATION
    ar & p.gamma_rot();
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & p.fixed();
    ar & p.ext_force();
#ifdef ESPRESSO_ROTATION
    ar & p.ext_torque();
#endif
#endif
#ifdef ESPRESSO_ENGINE
    ar & p.swimming();
#endif
  }
  if (data_parts & GHOSTTRANS_POSITION) {
    /* Position has no reduction policy: it is MOVE semantics for both
     * ReductionPolicy values. Branch on direction first (the lesson from the
     * FORCE path): on LOAD we always write the received value INTO the
     * particle; on SAVE we read the value FROM the particle, optionally with
     * the ghost shift + fold applied. Never bind particle-struct memory to the
     * archive directly, so a future field-storage flip cannot serialize a
     * proxy. Wire layout is byte-identical to the previous implementation.
     *
     * When a PositionRowContext is present AND the particle is attached
     * (store_row() >= 0), read/write pos/image (and quat, pos_last_time_step)
     * directly through the raw column pointers rather than the per-field
     * accessor proxies -- this is the mixed-parts (POSITION|PROPERTIES) hot
     * path. A detached particle (resort-step fresh ghost) or an absent context
     * falls back to the accessor, which uses the migration carriers. Local
     * copies are still made for the archive (SAVE) / from the archive (LOAD),
     * so the shift+fold and the assignment semantics are exactly as before. */
    auto const row = p.store_row();
    bool const use_ctx = pos_ctx != nullptr and row >= 0;
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d position;
      Utils::Vector3i image_box;
      ar & position;
      ar & image_box;
      if (use_ctx) {
        auto *pos = pos_ctx->pos +
                    static_cast<std::size_t>(row) * pos_ctx->pos_row_stride;
        auto *img = pos_ctx->img +
                    static_cast<std::size_t>(row) * pos_ctx->img_row_stride;
        pos[0 * pos_ctx->pos_comp_stride] = position[0];
        pos[1 * pos_ctx->pos_comp_stride] = position[1];
        pos[2 * pos_ctx->pos_comp_stride] = position[2];
        img[0 * pos_ctx->img_comp_stride] = image_box[0];
        img[1 * pos_ctx->img_comp_stride] = image_box[1];
        img[2 * pos_ctx->img_comp_stride] = image_box[2];
      } else {
        p.pos() = position;
        p.image_box() = image_box;
      }
    } else {
      Utils::Vector3d position;
      Utils::Vector3i image_box;
      if (use_ctx) {
        auto const *pos = pos_ctx->pos + static_cast<std::size_t>(row) *
                                             pos_ctx->pos_row_stride;
        auto const *img = pos_ctx->img + static_cast<std::size_t>(row) *
                                             pos_ctx->img_row_stride;
        position = {pos[0 * pos_ctx->pos_comp_stride],
                    pos[1 * pos_ctx->pos_comp_stride],
                    pos[2 * pos_ctx->pos_comp_stride]};
        image_box = {img[0 * pos_ctx->img_comp_stride],
                     img[1 * pos_ctx->img_comp_stride],
                     img[2 * pos_ctx->img_comp_stride]};
      } else {
        position = p.pos();
        image_box = p.image_box();
      }
      if (ghost_shift != nullptr) {
        /* ok, this is not nice, but perhaps fast */
        position += *ghost_shift;
        box_geo.fold_position(position, image_box);
      }
      ar & position;
      ar & image_box;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Quaternion<double> quat;
      ar & quat;
      if (use_ctx) {
        auto *q = pos_ctx->quat +
                  static_cast<std::size_t>(row) * pos_ctx->quat_row_stride;
        q[0 * pos_ctx->quat_comp_stride] = quat[0];
        q[1 * pos_ctx->quat_comp_stride] = quat[1];
        q[2 * pos_ctx->quat_comp_stride] = quat[2];
        q[3 * pos_ctx->quat_comp_stride] = quat[3];
      } else {
        p.quat() = quat;
      }
    } else {
      Utils::Quaternion<double> quat;
      if (use_ctx) {
        auto const *q = pos_ctx->quat + static_cast<std::size_t>(row) *
                                            pos_ctx->quat_row_stride;
        quat[0] = q[0 * pos_ctx->quat_comp_stride];
        quat[1] = q[1 * pos_ctx->quat_comp_stride];
        quat[2] = q[2 * pos_ctx->quat_comp_stride];
        quat[3] = q[3 * pos_ctx->quat_comp_stride];
      } else {
        quat = p.quat();
      }
      ar & quat;
    }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d pos_last_time_step;
      ar & pos_last_time_step;
      if (use_ctx) {
        auto *pl = pos_ctx->plast +
                   static_cast<std::size_t>(row) * pos_ctx->plast_row_stride;
        pl[0 * pos_ctx->plast_comp_stride] = pos_last_time_step[0];
        pl[1 * pos_ctx->plast_comp_stride] = pos_last_time_step[1];
        pl[2 * pos_ctx->plast_comp_stride] = pos_last_time_step[2];
      } else {
        p.pos_last_time_step() = pos_last_time_step;
      }
    } else {
      Utils::Vector3d pos_last_time_step;
      if (use_ctx) {
        auto const *pl = pos_ctx->plast + static_cast<std::size_t>(row) *
                                              pos_ctx->plast_row_stride;
        pos_last_time_step = {pl[0 * pos_ctx->plast_comp_stride],
                              pl[1 * pos_ctx->plast_comp_stride],
                              pl[2 * pos_ctx->plast_comp_stride]};
      } else {
        pos_last_time_step = p.pos_last_time_step();
      }
      ar & pos_last_time_step;
    }
#endif
  }
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    ar & p.v();
#ifdef ESPRESSO_ROTATION
    ar & p.omega();
#endif
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d force;
      ar & force;
      if (policy == ReductionPolicy::UPDATE) {
        p.force() += force;
      } else {
        p.force() = force;
      }
    } else {
      Utils::Vector3d force = p.force();
      ar & force;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d torque;
      ar & torque;
      if (policy == ReductionPolicy::UPDATE) {
        p.torque() += torque;
      } else {
        p.torque() = torque;
      }
    } else {
      Utils::Vector3d torque = p.torque();
      ar & torque;
    }
#endif
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (data_parts & GHOSTTRANS_RATTLE) {
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d correction;
      ar & correction;
      p.rattle_correction() += correction;
    } else {
      ar & p.rattle_correction();
    }
  }
#endif
}

/* ------------------------------------------------------------------------- *
 *  Columnar fast paths for the two per-step ghost communications.
 *
 *  The two hot, per-integration-step ghost transfers carry exactly one field
 *  group each: POSITION (ghost position distribution) and FORCE (ghost force
 *  reduction). For these the ParticleStore state/observable columns can be
 *  streamed to/from the CommBuf in bulk over contiguous store-row ranges,
 *  bypassing the per-Particle proxy construction and per-field archive
 *  dispatch that the generic serialize_and_reduce path incurs.
 *
 *  Wire layout is BYTE-IDENTICAL to serialize_and_reduce (the MemcpyArchive
 *  packs each field as its raw component bytes, tightly, no padding):
 *    POSITION: pos (3 double) | image (3 int)
 *              [| quat (4 double)] [| pos_last_time_step (3 double)]
 *    FORCE:    force (3 double) [| torque (3 double)]
 *  so buffer sizes (calc_transmit_size, untouched) and the GHOST_RDCE bitwise
 *  reduction over raw doubles stay identical.
 *
 *  Applicability (see ghost_communicator): the fast path is taken only when
 *  data_parts is EXACTLY GHOSTTRANS_POSITION or GHOSTTRANS_FORCE, AND the
 *  ParticleStore is not dirty. FORCE reduction always runs on a clean store
 *  (calculate_forces syncs the store, then reduces with no intervening
 *  topology change). POSITION update runs clean on non-resort steps, but on a
 *  resort step it runs while the store is dirty (ghosts_count marks it dirty;
 *  the sync happens later, at the top of calculate_forces): freshly created
 *  ghost particles are then DETACHED (store()==nullptr) and the generic path
 *  correctly reads/writes their migration carriers rather than columns. We
 *  therefore never sync here (that would be wrong: rows must match the current
 *  columns, not a fresh generation) -- we simply fall back to the per-particle
 *  path when the store is dirty. This is the safe branch required by the task,
 *  and it is exercised every resort step.
 * ------------------------------------------------------------------------- */

/** @brief The ParticleStore backing the current cell structure, or nullptr if
 *  no system/cell structure is set up (e.g. unit-test harnesses). */
static ParticleStore *active_particle_store() {
  auto const &system = System::get_system();
  if (system.cell_structure == nullptr) {
    return nullptr;
  }
  return &system.cell_structure->particle_store();
}

/** @brief Contiguous store-row range [first_row, first_row + size) for a
 *  part_list, or std::nullopt if the list is empty or its first particle is
 *  detached (store()==nullptr / store_row() < 0).
 *
 *  Contiguity is a STRUCTURAL invariant, not something we scan for in release:
 *  ensure_particle_store_synchronized assigns rows in cell-traversal order
 *  (every local cell's particles, then every ghost cell's particles), and every
 *  ghost-comm part_list is exactly one cell's ParticleList, so each cell's
 *  particles occupy one consecutive row block. The fast path is only reachable
 *  on a clean store (columnar_eligible), where this holds by construction. We
 *  therefore derive the range from the FIRST particle alone
 *  (first_row = first->store_row(); size = list size) with no per-particle
 * scan. A detached first particle (empty list or dirty/fresh-ghost store)
 * yields nullopt, detectable BEFORE any caller mutates the store, and callers
 * then fall back to the generic per-particle communication for the whole comm.
 *
 *  In debug builds ONLY, the whole list is validated up front (before any
 * caller reads/writes a column) and a contiguity/attachment violation asserts
 * with a clear message. */
static std::optional<std::pair<int, std::size_t>>
contiguous_store_rows(ParticleList const &part_list) {
  auto const size = part_list.size();
  if (size == 0u) {
    return std::nullopt;
  }
  auto const first = part_list.begin();
  if (first->store() == nullptr) {
    return std::nullopt;
  }
  auto const first_row = first->store_row();
  if (first_row < 0) {
    return std::nullopt;
  }
#ifndef NDEBUG
  for (std::size_t i = 0u; i < size; ++i) {
    auto const &p = first[i];
    if (p.store() == nullptr or
        p.store_row() != first_row + static_cast<int>(i)) {
      assert(false and "ghost part_list store rows are not contiguous");
    }
  }
#endif
  return std::make_pair(first_row, size);
}

/** @brief Resolve the contiguous store-row range for every part_list in a comm,
 *  deciding fast-path applicability BEFORE any caller mutates the store.
 *
 *  Returns std::nullopt (fall back to the generic path for the WHOLE comm) if
 *  any non-empty part_list has a detached first particle; empty lists yield an
 *  entry with size 0. Reads only each list's first particle (no per-particle
 *  scan): contiguity within a list is a structural invariant on a clean store
 *  (see contiguous_store_rows). Because the fall-back decision is taken here,
 *  before any column is read or written, the mutating callers (FORCE `+=`,
 *  POSITION assign) can never partially apply and then fall back -- which for
 *  FORCE would double-count. */
static std::optional<std::vector<std::pair<int, std::size_t>>>
columnar_resolve_ranges(GhostCommunication const &ghost_comm) {
  std::vector<std::pair<int, std::size_t>> ranges;
  ranges.reserve(ghost_comm.part_lists.size());
  for (auto part_list : ghost_comm.part_lists) {
    auto const range = contiguous_store_rows(*part_list);
    if (not range.has_value()) {
      if (part_list->empty()) {
        ranges.emplace_back(0, 0u);
        continue;
      }
      return std::nullopt; // detached first particle: fall back before mutating
    }
    ranges.push_back(*range);
  }
  return ranges;
}

/** @brief Whether data_parts selects a columnar-eligible per-step comm on a
 *  clean store. When false, callers use the generic per-particle path. */
static bool columnar_eligible(unsigned int data_parts) {
  if (data_parts != GHOSTTRANS_POSITION and data_parts != GHOSTTRANS_FORCE) {
    return false;
  }
  auto const *store = active_particle_store();
  return store != nullptr and not store->is_dirty();
}

/* ------------------------------------------------------------------------- *
 *  Per-communication row-access context for the MIXED-parts POSITION branch.
 *
 *  The per-step ghost update slot NEVER carries a pure GHOSTTRANS_POSITION
 *  transfer: get_global_ghost_flags() always yields POSITION|PROPERTIES, so
 *  the whole-communication POSITION-only fast path above cannot engage for it.
 *  It therefore goes through the generic per-particle serialize_and_reduce,
 *  whose POSITION branch would otherwise route every pos/image/quat/pos_last
 *  field through a per-field VectorReference/QuaternionReference proxy (SAVE:
 *  proxy construct + materialize per field; LOAD: proxy construct + write per
 *  field), ~2k ghosts x both sides x per step -- the measured +49% regression
 *  in update_ghosts_and_resort_particle.
 *
 *  The context captures, ONCE per ghost_communicator/cell_cell_transfer
 *  invocation (only when the store is present and clean, WITHOUT the
 *  data_parts==POSITION restriction of columnar_eligible), the raw host base
 *  pointers + row/component strides for the four POSITION columns. Any attached
 *  particle's row-r vector then lives at base + r*row_stride, its components at
 *  component_stride apart (component_stride==1 and the row contiguous under the
 *  store's LayoutRight, matching the accessor proxies' invariants). The generic
 *  POSITION branch reads/writes these directly by row instead of through the
 *  accessor, producing a byte-identical wire layout.
 *
 *  Passed to serialize_and_reduce as an optional pointer (nullptr = accessor
 *  fallback), and consulted only when p.store_row() >= 0 (an attached
 *  particle). Detached ghosts during a resort-step update (store()==nullptr /
 *  store_row() < 0) and any dirty-store situation fall back to the accessor,
 *  which correctly reads/writes their migration carriers.
 * ------------------------------------------------------------------------- */

/** @brief Build a PositionRowContext for the current clean store, or
 *  std::nullopt if there is no store or it is dirty (accessor fallback). This
 *  intentionally omits the data_parts==POSITION gate of columnar_eligible: it
 *  is meant for the mixed-parts (POSITION|PROPERTIES) per-step update. */
static std::optional<PositionRowContext> make_position_row_context() {
  auto *store = active_particle_store();
  if (store == nullptr or store->is_dirty() or
      store->number_of_particles() == 0u) {
    return std::nullopt;
  }
  auto pos = store->position_view();
  auto img = store->image_box_view();
  PositionRowContext ctx{};
  ctx.pos = &pos(0, 0);
  ctx.img = &img(0, 0);
  ctx.pos_row_stride = pos.stride(0);
  ctx.pos_comp_stride = pos.stride(1);
  ctx.img_row_stride = img.stride(0);
  ctx.img_comp_stride = img.stride(1);
#ifdef ESPRESSO_ROTATION
  auto quat = store->quaternion_view();
  ctx.quat = &quat(0, 0);
  ctx.quat_row_stride = quat.stride(0);
  ctx.quat_comp_stride = quat.stride(1);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store->position_last_time_step_view();
  ctx.plast = &plast(0, 0);
  ctx.plast_row_stride = plast.stride(0);
  ctx.plast_comp_stride = plast.stride(1);
#endif
  return ctx;
}

/** @brief Build a ForceRowContext for the current store (caller has already
 *  established store validity/cleanliness). */
static ForceRowContext make_force_row_context(ParticleStore &store) {
  ForceRowContext ctx{};
  auto force = store.force_view();
  ctx.force = &force(0, 0);
  ctx.force_row_stride = force.stride(0);
  ctx.force_comp_stride = force.stride(1);
#ifdef ESPRESSO_ROTATION
  auto torque = store.torque_view();
  ctx.torque = &torque(0, 0);
  ctx.torque_row_stride = torque.stride(0);
  ctx.torque_comp_stride = torque.stride(1);
#endif
  return ctx;
}

/** @brief Build a PositionRowContext unconditionally from a valid store
 *  (helper for the columnar bulk paths, which have already checked
 *  cleanliness via columnar_eligible/columnar_resolve_ranges). */
static PositionRowContext
make_position_row_context_unchecked(ParticleStore &store) {
  PositionRowContext ctx{};
  auto pos = store.position_view();
  auto img = store.image_box_view();
  ctx.pos = &pos(0, 0);
  ctx.img = &img(0, 0);
  ctx.pos_row_stride = pos.stride(0);
  ctx.pos_comp_stride = pos.stride(1);
  ctx.img_row_stride = img.stride(0);
  ctx.img_comp_stride = img.stride(1);
#ifdef ESPRESSO_ROTATION
  auto quat = store.quaternion_view();
  ctx.quat = &quat(0, 0);
  ctx.quat_row_stride = quat.stride(0);
  ctx.quat_comp_stride = quat.stride(1);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store.position_last_time_step_view();
  ctx.plast = &plast(0, 0);
  ctx.plast_row_stride = plast.stride(0);
  ctx.plast_comp_stride = plast.stride(1);
#endif
  return ctx;
}

namespace {
/** @brief Byte cursor over the CommBuf that packs/unpacks scalars via memcpy,
 *  exactly like the MemcpyArchive. Using memcpy (not a reinterpret_cast +
 *  dereference) is REQUIRED: the tightly packed per-particle wire layout is not
 *  8-byte aligned (e.g. 36 B/particle for POSITION without ROTATION), so a
 *  direct double* dereference would be misaligned UB (and trip UBSan). */
template <class CharPtr> struct ByteCursor {
  CharPtr ptr;
  template <class T> void put(T const &value) {
    std::memcpy(ptr, &value, sizeof(T));
    ptr += sizeof(T);
  }
  template <class T> T get() {
    T value;
    std::memcpy(&value, ptr, sizeof(T));
    ptr += sizeof(T);
    return value;
  }
};
using WriteCursor = ByteCursor<char *>;
using ReadCursor = ByteCursor<char const *>;
} // namespace

/** @brief Bulk-pack POSITION for a contiguous row range into @p out (SAVE).
 *  Reproduces serialize_and_reduce's SAVE branch exactly: with a ghost shift,
 *  apply +shift then fold_position (writing the folded image); otherwise copy
 *  the stored pos/image. Then quat (ROTATION) and pos_last_time_step
 *  (BOND_CONSTRAINT) verbatim. */
static char *pack_position_range(char *out, PositionRowContext const &ctx,
                                 int first_row, std::size_t n,
                                 BoxGeometry const &box_geo,
                                 Utils::Vector3d const *ghost_shift) {
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto const *p = ctx.pos + row * ctx.pos_row_stride;
    auto const *im = ctx.img + row * ctx.img_row_stride;
    Utils::Vector3d position{p[0], p[ctx.pos_comp_stride],
                             p[2u * ctx.pos_comp_stride]};
    Utils::Vector3i image_box{im[0], im[ctx.img_comp_stride],
                              im[2u * ctx.img_comp_stride]};
    if (ghost_shift != nullptr) {
      position += *ghost_shift;
      box_geo.fold_position(position, image_box);
    }
    cursor.put(position[0]);
    cursor.put(position[1]);
    cursor.put(position[2]);
    cursor.put(image_box[0]);
    cursor.put(image_box[1]);
    cursor.put(image_box[2]);
#ifdef ESPRESSO_ROTATION
    auto const *q = ctx.quat + row * ctx.quat_row_stride;
    cursor.put(q[0]);
    cursor.put(q[ctx.quat_comp_stride]);
    cursor.put(q[2u * ctx.quat_comp_stride]);
    cursor.put(q[3u * ctx.quat_comp_stride]);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto const *pl = ctx.plast + row * ctx.plast_row_stride;
    cursor.put(pl[0]);
    cursor.put(pl[ctx.plast_comp_stride]);
    cursor.put(pl[2u * ctx.plast_comp_stride]);
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack POSITION for a contiguous row range from @p in (LOAD).
 *  Assigns pos/image (and quat, pos_last_time_step) directly into columns. */
static char const *unpack_position_range(char const *in,
                                         PositionRowContext const &ctx,
                                         int first_row, std::size_t n) {
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto *p = ctx.pos + row * ctx.pos_row_stride;
    auto *im = ctx.img + row * ctx.img_row_stride;
    p[0] = cursor.get<double>();
    p[ctx.pos_comp_stride] = cursor.get<double>();
    p[2u * ctx.pos_comp_stride] = cursor.get<double>();
    im[0] = cursor.get<int>();
    im[ctx.img_comp_stride] = cursor.get<int>();
    im[2u * ctx.img_comp_stride] = cursor.get<int>();
#ifdef ESPRESSO_ROTATION
    auto *q = ctx.quat + row * ctx.quat_row_stride;
    q[0] = cursor.get<double>();
    q[ctx.quat_comp_stride] = cursor.get<double>();
    q[2u * ctx.quat_comp_stride] = cursor.get<double>();
    q[3u * ctx.quat_comp_stride] = cursor.get<double>();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto *pl = ctx.plast + row * ctx.plast_row_stride;
    pl[0] = cursor.get<double>();
    pl[ctx.plast_comp_stride] = cursor.get<double>();
    pl[2u * ctx.plast_comp_stride] = cursor.get<double>();
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-pack FORCE for a contiguous row range into @p out (SAVE). */
static char *pack_force_range(char *out, ForceRowContext const &ctx,
                              int first_row, std::size_t n) {
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto const *f = ctx.force + row * ctx.force_row_stride;
    cursor.put(f[0]);
    cursor.put(f[ctx.force_comp_stride]);
    cursor.put(f[2u * ctx.force_comp_stride]);
#ifdef ESPRESSO_ROTATION
    auto const *t = ctx.torque + row * ctx.torque_row_stride;
    cursor.put(t[0]);
    cursor.put(t[ctx.torque_comp_stride]);
    cursor.put(t[2u * ctx.torque_comp_stride]);
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack FORCE for a contiguous row range from @p in.
 *  With @p accumulate (UPDATE policy) the received value is ADDED; otherwise
 *  (MOVE policy) it is assigned. Matches serialize_and_reduce's FORCE LOAD. */
static char const *unpack_force_range(char const *in,
                                      ForceRowContext const &ctx, int first_row,
                                      std::size_t n, bool accumulate) {
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = static_cast<std::size_t>(first_row) + i;
    auto *f = ctx.force + row * ctx.force_row_stride;
    if (accumulate) {
      f[0] += cursor.get<double>();
      f[ctx.force_comp_stride] += cursor.get<double>();
      f[2u * ctx.force_comp_stride] += cursor.get<double>();
    } else {
      f[0] = cursor.get<double>();
      f[ctx.force_comp_stride] = cursor.get<double>();
      f[2u * ctx.force_comp_stride] = cursor.get<double>();
    }
#ifdef ESPRESSO_ROTATION
    auto *t = ctx.torque + row * ctx.torque_row_stride;
    if (accumulate) {
      t[0] += cursor.get<double>();
      t[ctx.torque_comp_stride] += cursor.get<double>();
      t[2u * ctx.torque_comp_stride] += cursor.get<double>();
    } else {
      t[0] = cursor.get<double>();
      t[ctx.torque_comp_stride] = cursor.get<double>();
      t[2u * ctx.torque_comp_stride] = cursor.get<double>();
    }
#endif
  }
  return cursor.ptr;
}

static auto calc_transmit_size(BoxGeometry const &box_geo,
                               unsigned data_parts) {
  std::size_t force_size = 0ul;
  if (data_parts & GHOSTTRANS_FORCE) {
#ifdef ESPRESSO_ROTATION
    force_size = 6ul * sizeof(double);
#else
    force_size = 3ul * sizeof(double);
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_FORCE);
    if (data_parts == 0u) {
      return force_size;
    }
  }
  std::size_t position_size = 0ul;
  if (data_parts & GHOSTTRANS_POSITION) {
    /* Compositional, compile-time-constant size of the POSITION wire layout:
     * pos (Vector3d, 24 B) + image_box (Vector3i, 12 B)
     * [+ quat (Quaternion<double>, 32 B) under ROTATION]
     * [+ pos_last_time_step (Vector3d, 24 B) under BOND_CONSTRAINT].
     * The MemcpyArchive packs every field tightly (exactly sizeof(T), no
     * alignment padding), and SerializationSizeCalculator likewise accumulates
     * sizeof(T), so this constant matches the sizer output exactly. */
    position_size = 3ul * sizeof(double) + 3ul * sizeof(int);
#ifdef ESPRESSO_ROTATION
    position_size += sizeof(Utils::Quaternion<double>);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    position_size += 3ul * sizeof(double);
#endif
    data_parts &= ~static_cast<unsigned>(GHOSTTRANS_POSITION);
    if (data_parts == 0u) {
      return force_size + position_size;
    }
  }
  SerializationSizeCalculator sizeof_archive;
  Particle p{};
  serialize_and_reduce(sizeof_archive, p, data_parts, ReductionPolicy::MOVE,
                       SerializationDirection::SAVE, box_geo, nullptr);
  return sizeof_archive.size() + force_size + position_size;
}

static auto calc_transmit_size(GhostCommunication const &ghost_comm,
                               BoxGeometry const &box_geo,
                               unsigned int data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(unsigned int) * ghost_comm.part_lists.size();

  auto const n_part = boost::accumulate(
      ghost_comm.part_lists, std::size_t{0},
      [](std::size_t sum, auto part_list) { return sum + part_list->size(); });

  return n_part * calc_transmit_size(box_geo, data_parts);
}

/** @brief Columnar bulk fill of @p send_buffer for POSITION/FORCE on a clean
 *  store. Returns false (leaving the buffer untouched) if any part_list breaks
 *  the contiguity/attachment precondition, so the caller can fall back. */
static bool columnar_prepare_send_buffer(CommBuf &send_buffer,
                                         GhostCommunication const &ghost_comm,
                                         BoxGeometry const &box_geo,
                                         unsigned int data_parts) {
  auto const ranges = columnar_resolve_ranges(ghost_comm);
  if (not ranges.has_value()) {
    return false; // decided before touching the buffer/store
  }
  auto &store = *active_particle_store();
  auto const position_ctx =
      (data_parts == GHOSTTRANS_POSITION)
          ? std::optional<PositionRowContext>(
                make_position_row_context_unchecked(store))
          : std::nullopt;
  auto const force_ctx =
      (data_parts == GHOSTTRANS_FORCE)
          ? std::optional<ForceRowContext>(make_force_row_context(store))
          : std::nullopt;
  auto *cursor = send_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    if (data_parts == GHOSTTRANS_POSITION) {
      cursor = pack_position_range(cursor, *position_ctx, first_row, n, box_geo,
                                   &ghost_comm.shift);
    } else {
      cursor = pack_force_range(cursor, *force_ctx, first_row, n);
    }
  }
  assert(static_cast<std::size_t>(cursor - send_buffer.data()) ==
         send_buffer.size());
  return true;
}

static void prepare_send_buffer(CommBuf &send_buffer,
                                GhostCommunication const &ghost_comm,
                                BoxGeometry const &box_geo,
                                unsigned int data_parts) {

  /* reallocate send buffer */
  send_buffer.resize(calc_transmit_size(ghost_comm, box_geo, data_parts));
  send_buffer.bonds().clear();

  /* Columnar fast path for the two per-step comms on a clean store. */
  if (columnar_eligible(data_parts) and
      columnar_prepare_send_buffer(send_buffer, ghost_comm, box_geo,
                                   data_parts)) {
    return;
  }

  auto archiver = Utils::MemcpyOArchive{send_buffer.make_span()};

  /* Mixed-parts POSITION fast path: read pos/image/quat/pos_last straight from
   * the (clean) store columns by row instead of via per-field accessor proxies.
   * nullptr when POSITION is absent or the store is dirty/absent (accessor
   * fallback, per-particle). */
  auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                           ? make_position_row_context()
                           : std::nullopt;
  auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;

  /* Construct archive that pushes back to the bond buffer */
  namespace io = boost::iostreams;
  io::stream<io::back_insert_device<std::vector<char>>> os{
      io::back_inserter(send_buffer.bonds())};
  boost::archive::binary_oarchive bond_archiver{os};

  /* put in data */
  for (auto part_list : ghost_comm.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      assert(part_list->size() <= std::numeric_limits<unsigned int>::max());
      auto np = static_cast<unsigned int>(part_list->size());
      archiver << np;
    } else {
      for (auto &p : *part_list) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::SAVE, box_geo,
                             &ghost_comm.shift, pos_ctx_ptr);
        if (data_parts & GHOSTTRANS_BONDS) {
          bond_archiver << p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_written() == send_buffer.size());
}

static void prepare_recv_buffer(CommBuf &recv_buffer,
                                GhostCommunication const &ghost_comm,
                                BoxGeometry const &box_geo,
                                unsigned int data_parts) {
  /* reallocate recv buffer */
  recv_buffer.resize(calc_transmit_size(ghost_comm, box_geo, data_parts));
  /* clear bond buffer */
  recv_buffer.bonds().clear();
}

/** @brief Columnar bulk unpack of @p recv_buffer into columns for POSITION on
 *  a clean store. Returns false (untouched) on a broken precondition. */
static bool columnar_put_recv_buffer(CommBuf &recv_buffer,
                                     GhostCommunication const &ghost_comm,
                                     unsigned int data_parts) {
  assert(data_parts == GHOSTTRANS_POSITION);
  static_cast<void>(data_parts);
  auto const ranges = columnar_resolve_ranges(ghost_comm);
  if (not ranges.has_value()) {
    return false; // decided before writing any column
  }
  auto &store = *active_particle_store();
  auto const position_ctx = make_position_row_context_unchecked(store);
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_position_range(cursor, position_ctx, first_row, n);
  }
  assert(static_cast<std::size_t>(cursor - recv_buffer.data()) ==
         recv_buffer.size());
  return true;
}

static void put_recv_buffer(CommBuf &recv_buffer,
                            GhostCommunication const &ghost_comm,
                            BoxGeometry const &box_geo,
                            unsigned int data_parts) {
  /* Columnar fast path: pure POSITION on a clean store. */
  if (data_parts == GHOSTTRANS_POSITION and columnar_eligible(data_parts) and
      columnar_put_recv_buffer(recv_buffer, ghost_comm, data_parts)) {
    recv_buffer.bonds().clear();
    return;
  }

  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};

  if (data_parts & GHOSTTRANS_PARTNUM) {
    for (auto part_list : ghost_comm.part_lists) {
      unsigned int np;
      archiver >> np;
      CellParticleStorage::resize_ghost_storage(*part_list, np);
    }
  } else {
    /* Mixed-parts POSITION fast path (see prepare_send_buffer): write
     * pos/image/quat/pos_last straight into the (clean) store columns by row.
     * nullptr => accessor fallback (POSITION absent, or dirty/absent store, or
     * detached fresh ghosts on a resort step). */
    auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                             ? make_position_row_context()
                             : std::nullopt;
    auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;
    for (auto part_list : ghost_comm.part_lists) {
      for (auto &p : *part_list) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::LOAD, box_geo, nullptr,
                             pos_ctx_ptr);
      }
    }
    if (data_parts & GHOSTTRANS_BONDS) {
      namespace io = boost::iostreams;
      io::stream<io::array_source> bond_stream(io::array_source{
          recv_buffer.bonds().data(), recv_buffer.bonds().size()});
      boost::archive::binary_iarchive bond_archiver(bond_stream);

      for (auto part_list : ghost_comm.part_lists) {
        for (auto &p : *part_list) {
          bond_archiver >> p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_read() == recv_buffer.size());

  recv_buffer.bonds().clear();
}

#ifdef ESPRESSO_BOND_CONSTRAINT
static void
add_rattle_correction_from_recv_buffer(CommBuf &recv_buffer,
                                       const GhostCommunication &ghost_comm) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : *part_list) {
      ParticleRattle pr;
      archiver >> pr;
      part.rattle_params() += pr;
    }
  }
}
#endif

/** @brief Columnar bulk force accumulation from @p recv_buffer on a clean
 *  store. Returns false (untouched) on a broken precondition. */
static bool
columnar_add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                     GhostCommunication const &ghost_comm) {
  /* The fall-back decision is made BEFORE any `+=`: a mid-list fall-back after
   * a partial accumulate would double-count on the generic retry. */
  auto const ranges = columnar_resolve_ranges(ghost_comm);
  if (not ranges.has_value()) {
    return false;
  }
  auto &store = *active_particle_store();
  auto const force_ctx = make_force_row_context(store);
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_force_range(cursor, force_ctx, first_row, n,
                                /*accumulate=*/true);
  }
  assert(static_cast<std::size_t>(cursor - recv_buffer.data()) ==
         recv_buffer.size());
  return true;
}

static void add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                        const GhostCommunication &ghost_comm) {
  /* Columnar fast path: FORCE reduction always runs on a clean store. */
  if (columnar_eligible(GHOSTTRANS_FORCE) and
      columnar_add_forces_from_recv_buffer(recv_buffer, ghost_comm)) {
    return;
  }

  /* put back data */
  auto archiver = Utils::MemcpyIArchive{recv_buffer.make_span()};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : *part_list) {
      Utils::Vector3d force;
      archiver >> force;
      part.force() += force;
#ifdef ESPRESSO_ROTATION
      Utils::Vector3d torque;
      archiver >> torque;
      part.torque() += torque;
#endif
    }
  }
}

/** @brief Columnar POSITION cell-to-cell transfer: shift+fold the source rows
 *  and assign into the destination rows (quat, pos_last_time_step verbatim).
 *  Matches the generic LOCL SAVE(shift)+LOAD(assign) semantics: the shift is
 *  ALWAYS applied (the generic LOCL SAVE always passes a non-null shift), and
 * it is the zero vector for a local cell-to-cell transfer (no periodic wrap),
 * so the `+= shift` is a no-op before the fold in that common case rather than
 * a branch. */
static void locl_transfer_position(PositionRowContext const &ctx, int src_first,
                                   int dst_first, std::size_t n,
                                   BoxGeometry const &box_geo,
                                   Utils::Vector3d const &shift) {
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = static_cast<std::size_t>(src_first) + i;
    auto const d = static_cast<std::size_t>(dst_first) + i;
    auto const *ps = ctx.pos + s * ctx.pos_row_stride;
    auto const *is = ctx.img + s * ctx.img_row_stride;
    auto *pd = ctx.pos + d * ctx.pos_row_stride;
    auto *id = ctx.img + d * ctx.img_row_stride;
    Utils::Vector3d position{ps[0], ps[ctx.pos_comp_stride],
                             ps[2u * ctx.pos_comp_stride]};
    Utils::Vector3i image_box{is[0], is[ctx.img_comp_stride],
                              is[2u * ctx.img_comp_stride]};
    position += shift;
    box_geo.fold_position(position, image_box);
    pd[0] = position[0];
    pd[ctx.pos_comp_stride] = position[1];
    pd[2u * ctx.pos_comp_stride] = position[2];
    id[0] = image_box[0];
    id[ctx.img_comp_stride] = image_box[1];
    id[2u * ctx.img_comp_stride] = image_box[2];
#ifdef ESPRESSO_ROTATION
    auto const *qs = ctx.quat + s * ctx.quat_row_stride;
    auto *qd = ctx.quat + d * ctx.quat_row_stride;
    qd[0] = qs[0];
    qd[ctx.quat_comp_stride] = qs[ctx.quat_comp_stride];
    qd[2u * ctx.quat_comp_stride] = qs[2u * ctx.quat_comp_stride];
    qd[3u * ctx.quat_comp_stride] = qs[3u * ctx.quat_comp_stride];
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto const *ls = ctx.plast + s * ctx.plast_row_stride;
    auto *ld = ctx.plast + d * ctx.plast_row_stride;
    ld[0] = ls[0];
    ld[ctx.plast_comp_stride] = ls[ctx.plast_comp_stride];
    ld[2u * ctx.plast_comp_stride] = ls[2u * ctx.plast_comp_stride];
#endif
  }
}

/** @brief Columnar FORCE cell-to-cell transfer: accumulate source rows into
 *  destination rows (UPDATE policy: += force and torque). */
static void locl_transfer_force(ForceRowContext const &ctx, int src_first,
                                int dst_first, std::size_t n) {
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = static_cast<std::size_t>(src_first) + i;
    auto const d = static_cast<std::size_t>(dst_first) + i;
    auto const *fs = ctx.force + s * ctx.force_row_stride;
    auto *fd = ctx.force + d * ctx.force_row_stride;
    fd[0] += fs[0];
    fd[ctx.force_comp_stride] += fs[ctx.force_comp_stride];
    fd[2u * ctx.force_comp_stride] += fs[2u * ctx.force_comp_stride];
#ifdef ESPRESSO_ROTATION
    auto const *ts = ctx.torque + s * ctx.torque_row_stride;
    auto *td = ctx.torque + d * ctx.torque_row_stride;
    td[0] += ts[0];
    td[ctx.torque_comp_stride] += ts[ctx.torque_comp_stride];
    td[2u * ctx.torque_comp_stride] += ts[2u * ctx.torque_comp_stride];
#endif
  }
}

/** @brief Columnar LOCL fast path. Returns false on a broken precondition
 *  (so the caller falls back to the per-particle path).
 *
 *  All src/dst row ranges are resolved (reading only each list's first
 *  particle, no per-particle scan) BEFORE any transfer, so the FORCE `+=`
 *  transfer can never partially apply and then fall back (which would
 *  double-count on the generic retry). */
static bool columnar_cell_cell_transfer(GhostCommunication const &ghost_comm,
                                        BoxGeometry const &box_geo,
                                        unsigned int data_parts) {
  auto const offset = ghost_comm.part_lists.size() / 2;
  // Phase 1: resolve every src/dst range pair; bail before mutating on break.
  std::vector<std::pair<int, int>> pairs; // (src_first, dst_first); n implicit
  std::vector<std::size_t> sizes;
  pairs.reserve(offset);
  sizes.reserve(offset);
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto *src_list = ghost_comm.part_lists[pl];
    auto *dst_list = ghost_comm.part_lists[pl + offset];
    assert(src_list->size() == dst_list->size());
    if (src_list->empty()) {
      pairs.emplace_back(0, 0);
      sizes.push_back(0u);
      continue;
    }
    auto const src_range = contiguous_store_rows(*src_list);
    auto const dst_range = contiguous_store_rows(*dst_list);
    if (not src_range.has_value() or not dst_range.has_value()) {
      return false;
    }
    pairs.emplace_back(src_range->first, dst_range->first);
    sizes.push_back(src_range->second);
  }
  // Phase 2: apply the transfers.
  auto &store = *active_particle_store();
  auto const position_ctx =
      (data_parts == GHOSTTRANS_POSITION)
          ? std::optional<PositionRowContext>(
                make_position_row_context_unchecked(store))
          : std::nullopt;
  auto const force_ctx =
      (data_parts == GHOSTTRANS_FORCE)
          ? std::optional<ForceRowContext>(make_force_row_context(store))
          : std::nullopt;
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto const n = sizes[pl];
    if (n == 0u) {
      continue;
    }
    auto const [src_first, dst_first] = pairs[pl];
    if (data_parts == GHOSTTRANS_POSITION) {
      locl_transfer_position(*position_ctx, src_first, dst_first, n, box_geo,
                             ghost_comm.shift);
    } else {
      locl_transfer_force(*force_ctx, src_first, dst_first, n);
    }
  }
  return true;
}

static void cell_cell_transfer(GhostCommunication const &ghost_comm,
                               BoxGeometry const &box_geo,
                               unsigned int data_parts) {
  /* Columnar fast path for the two per-step comms on a clean store. */
  if (columnar_eligible(data_parts) and
      columnar_cell_cell_transfer(ghost_comm, box_geo, data_parts)) {
    return;
  }

  CommBuf buffer;
  if (!(data_parts & GHOSTTRANS_PARTNUM)) {
    buffer.resize(calc_transmit_size(box_geo, data_parts));
  }
  /* Mixed-parts POSITION fast path (see prepare_send_buffer): both source and
   * destination particles read/write pos/image/quat/pos_last through the
   * (clean) store columns by their own row. nullptr => accessor fallback. */
  auto const pos_ctx = (data_parts & GHOSTTRANS_POSITION)
                           ? make_position_row_context()
                           : std::nullopt;
  auto const *pos_ctx_ptr = pos_ctx.has_value() ? &(*pos_ctx) : nullptr;
  /* transfer data */
  auto const offset = ghost_comm.part_lists.size() / 2;
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto *src_list = ghost_comm.part_lists[pl];
    auto *dst_list = ghost_comm.part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      CellParticleStorage::resize_ghost_storage(*dst_list, src_list->size());
    } else {
      auto &src_part = *src_list;
      auto &dst_part = *dst_list;
      assert(src_part.size() == dst_part.size());

      for (std::size_t i = 0; i < src_part.size(); i++) {
        auto ar_out = Utils::MemcpyOArchive{buffer.make_span()};
        auto ar_in = Utils::MemcpyIArchive{buffer.make_span()};
        auto &p1 = src_part.begin()[i];
        auto &p2 = dst_part.begin()[i];
        serialize_and_reduce(ar_out, p1, data_parts, ReductionPolicy::UPDATE,
                             SerializationDirection::SAVE, box_geo,
                             &ghost_comm.shift, pos_ctx_ptr);
        serialize_and_reduce(ar_in, p2, data_parts, ReductionPolicy::UPDATE,
                             SerializationDirection::LOAD, box_geo, nullptr,
                             pos_ctx_ptr);
        if (data_parts & GHOSTTRANS_BONDS) {
          p2.bonds() = p1.bonds();
        }
      }
    }
  }
}

static bool is_send_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
          (comm_type == GHOST_BCST && node == this_node));
}

static bool is_recv_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_RECV) ||
          (comm_type == GHOST_BCST && node != this_node) ||
          (comm_type == GHOST_RDCE && node == this_node));
}

static bool is_prefetchable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const prefetch = ghost_comm.type & GHOST_PREFETCH;
  int const node = ghost_comm.node;
  return is_send_op(comm_type, node, this_node) && prefetch;
}

static bool is_poststorable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const poststore = ghost_comm.type & GHOST_PSTSTORE;
  int const node = ghost_comm.node;
  return is_recv_op(comm_type, node, this_node) && poststore;
}

void ghost_communicator(GhostCommunicator const &gcr,
                        BoxGeometry const &box_geo, unsigned int data_parts) {
  if (GHOSTTRANS_NONE == data_parts)
    return;

  static CommBuf send_buffer, recv_buffer;

  auto const &comm = gcr.mpi_comm;

  for (auto cit = gcr.communications.cbegin(); cit != gcr.communications.cend();
       ++cit) {
    auto const &ghost_comm = *cit;
    int const comm_type = ghost_comm.type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(ghost_comm, box_geo, data_parts);
      continue;
    }

    int const prefetch = ghost_comm.type & GHOST_PREFETCH;
    int const poststore = ghost_comm.type & GHOST_PSTSTORE;
    int const node = ghost_comm.node;

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node, comm.rank())) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch) {
        prepare_send_buffer(send_buffer, ghost_comm, box_geo, data_parts);
      }
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      assert(send_buffer.size() ==
             calc_transmit_size(ghost_comm, box_geo, data_parts));
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto prefetch_ghost_comm =
          std::find_if(std::next(cit), gcr.communications.cend(),
                       [this_node = comm.rank()](auto const &other_ghost_comm) {
                         return is_prefetchable(other_ghost_comm, this_node);
                       });

      if (prefetch_ghost_comm != gcr.communications.end())
        prepare_send_buffer(send_buffer, *prefetch_ghost_comm, box_geo,
                            data_parts);
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node, comm.rank()))
      prepare_recv_buffer(recv_buffer, ghost_comm, box_geo, data_parts);

    /* transfer data */
    // Use two send/recvs in order to avoid having to serialize CommBuf
    // (which consists of already serialized data).
    switch (comm_type) {
    case GHOST_RECV:
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.data(),
                static_cast<int>(recv_buffer.size()));
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.bonds());
      break;
    case GHOST_SEND:
      comm.send(node, REQ_GHOST_SEND, send_buffer.data(),
                static_cast<int>(send_buffer.size()));
      comm.send(node, REQ_GHOST_SEND, send_buffer.bonds());
      break;
    case GHOST_BCST:
      if (node == comm.rank()) {
        boost::mpi::broadcast(comm, send_buffer.data(),
                              static_cast<int>(send_buffer.size()), node);
        boost::mpi::broadcast(comm, send_buffer.bonds(), node);
      } else {
        boost::mpi::broadcast(comm, recv_buffer.data(),
                              static_cast<int>(recv_buffer.size()), node);
        boost::mpi::broadcast(comm, recv_buffer.bonds(), node);
      }
      break;
    case GHOST_RDCE:
      if (node == comm.rank())
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            reinterpret_cast<double *>(recv_buffer.data()), std::plus<double>{},
            node);
      else
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            std::plus<double>{}, node);
      break;
    }

    // recv op; write back data directly, if no PSTSTORE delay is requested.
    if (is_recv_op(comm_type, node, comm.rank())) {
      if (!poststore) {
        /* forces have to be added, the rest overwritten. Exception is RDCE,
         * where the addition is integrated into the communication. */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, ghost_comm);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          add_rattle_correction_from_recv_buffer(recv_buffer, ghost_comm);
#endif
        else
          put_recv_buffer(recv_buffer, ghost_comm, box_geo, data_parts);
      }
    } else if (poststore) {
      /* send op; write back delayed data from last recv, when this was a
       * prefetch send. */
      /* find previous action where we recv and which has PSTSTORE set */
      auto poststore_ghost_comm = std::find_if(
          std::make_reverse_iterator(cit), gcr.communications.crend(),
          [this_node = comm.rank()](auto const &other_ghost_comm) {
            return is_poststorable(other_ghost_comm, this_node);
          });

      if (poststore_ghost_comm != gcr.communications.rend()) {
        assert(recv_buffer.size() ==
               calc_transmit_size(*poststore_ghost_comm, box_geo, data_parts));
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, *poststore_ghost_comm);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          add_rattle_correction_from_recv_buffer(recv_buffer,
                                                 *poststore_ghost_comm);
#endif
        else
          put_recv_buffer(recv_buffer, *poststore_ghost_comm, box_geo,
                          data_parts);
      }
    }
  }
}
