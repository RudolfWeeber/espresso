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
                     Utils::Vector3d const *ghost_shift) {
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
     * proxy. Wire layout is byte-identical to the previous implementation. */
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d position;
      Utils::Vector3i image_box;
      ar & position;
      ar & image_box;
      p.pos() = position;
      p.image_box() = image_box;
    } else if (ghost_shift != nullptr) {
      /* ok, this is not nice, but perhaps fast */
      Utils::Vector3d position = Utils::Vector3d(p.pos()) + *ghost_shift;
      Utils::Vector3i image_box = p.image_box();
      box_geo.fold_position(position, image_box);
      ar & position;
      ar & image_box;
    } else {
      Utils::Vector3d position = p.pos();
      Utils::Vector3i image_box = p.image_box();
      ar & position;
      ar & image_box;
    }
#ifdef ESPRESSO_ROTATION
    if (direction == SerializationDirection::LOAD) {
      Utils::Quaternion<double> quat;
      ar & quat;
      p.quat() = quat;
    } else {
      Utils::Quaternion<double> quat = p.quat();
      ar & quat;
    }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d pos_last_time_step;
      ar & pos_last_time_step;
      p.pos_last_time_step() = pos_last_time_step;
    } else {
      Utils::Vector3d pos_last_time_step = p.pos_last_time_step();
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
static char *pack_position_range(char *out, ParticleStore &store, int first_row,
                                 std::size_t n, BoxGeometry const &box_geo,
                                 Utils::Vector3d const *ghost_shift) {
  auto const pos = store.position_view();
  auto const img = store.image_box_view();
#ifdef ESPRESSO_ROTATION
  auto const quat = store.quaternion_view();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto const plast = store.position_last_time_step_view();
#endif
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = first_row + static_cast<int>(i);
    Utils::Vector3d position{pos(row, 0), pos(row, 1), pos(row, 2)};
    Utils::Vector3i image_box{img(row, 0), img(row, 1), img(row, 2)};
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
    cursor.put(quat(row, 0));
    cursor.put(quat(row, 1));
    cursor.put(quat(row, 2));
    cursor.put(quat(row, 3));
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    cursor.put(plast(row, 0));
    cursor.put(plast(row, 1));
    cursor.put(plast(row, 2));
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack POSITION for a contiguous row range from @p in (LOAD).
 *  Assigns pos/image (and quat, pos_last_time_step) directly into columns. */
static char const *unpack_position_range(char const *in, ParticleStore &store,
                                         int first_row, std::size_t n) {
  auto pos = store.position_view();
  auto img = store.image_box_view();
#ifdef ESPRESSO_ROTATION
  auto quat = store.quaternion_view();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store.position_last_time_step_view();
#endif
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = first_row + static_cast<int>(i);
    pos(row, 0) = cursor.get<double>();
    pos(row, 1) = cursor.get<double>();
    pos(row, 2) = cursor.get<double>();
    img(row, 0) = cursor.get<int>();
    img(row, 1) = cursor.get<int>();
    img(row, 2) = cursor.get<int>();
#ifdef ESPRESSO_ROTATION
    quat(row, 0) = cursor.get<double>();
    quat(row, 1) = cursor.get<double>();
    quat(row, 2) = cursor.get<double>();
    quat(row, 3) = cursor.get<double>();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    plast(row, 0) = cursor.get<double>();
    plast(row, 1) = cursor.get<double>();
    plast(row, 2) = cursor.get<double>();
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-pack FORCE for a contiguous row range into @p out (SAVE). */
static char *pack_force_range(char *out, ParticleStore &store, int first_row,
                              std::size_t n) {
  auto const force = store.force_view();
#ifdef ESPRESSO_ROTATION
  auto const torque = store.torque_view();
#endif
  WriteCursor cursor{out};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = first_row + static_cast<int>(i);
    cursor.put(force(row, 0));
    cursor.put(force(row, 1));
    cursor.put(force(row, 2));
#ifdef ESPRESSO_ROTATION
    cursor.put(torque(row, 0));
    cursor.put(torque(row, 1));
    cursor.put(torque(row, 2));
#endif
  }
  return cursor.ptr;
}

/** @brief Bulk-unpack FORCE for a contiguous row range from @p in.
 *  With @p accumulate (UPDATE policy) the received value is ADDED; otherwise
 *  (MOVE policy) it is assigned. Matches serialize_and_reduce's FORCE LOAD. */
static char const *unpack_force_range(char const *in, ParticleStore &store,
                                      int first_row, std::size_t n,
                                      bool accumulate) {
  auto force = store.force_view();
#ifdef ESPRESSO_ROTATION
  auto torque = store.torque_view();
#endif
  ReadCursor cursor{in};
  for (std::size_t i = 0u; i < n; ++i) {
    auto const row = first_row + static_cast<int>(i);
    auto const fx = cursor.get<double>();
    auto const fy = cursor.get<double>();
    auto const fz = cursor.get<double>();
    if (accumulate) {
      force(row, 0) += fx;
      force(row, 1) += fy;
      force(row, 2) += fz;
    } else {
      force(row, 0) = fx;
      force(row, 1) = fy;
      force(row, 2) = fz;
    }
#ifdef ESPRESSO_ROTATION
    auto const tx = cursor.get<double>();
    auto const ty = cursor.get<double>();
    auto const tz = cursor.get<double>();
    if (accumulate) {
      torque(row, 0) += tx;
      torque(row, 1) += ty;
      torque(row, 2) += tz;
    } else {
      torque(row, 0) = tx;
      torque(row, 1) = ty;
      torque(row, 2) = tz;
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
  auto *cursor = send_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    if (data_parts == GHOSTTRANS_POSITION) {
      cursor = pack_position_range(cursor, store, first_row, n, box_geo,
                                   &ghost_comm.shift);
    } else {
      cursor = pack_force_range(cursor, store, first_row, n);
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
                             &ghost_comm.shift);
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
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_position_range(cursor, store, first_row, n);
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
    for (auto part_list : ghost_comm.part_lists) {
      for (auto &p : *part_list) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::LOAD, box_geo, nullptr);
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
  auto const *cursor = recv_buffer.data();
  for (auto const &[first_row, n] : *ranges) {
    if (n == 0u) {
      continue;
    }
    cursor = unpack_force_range(cursor, store, first_row, n,
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
static void locl_transfer_position(ParticleStore &store, int src_first,
                                   int dst_first, std::size_t n,
                                   BoxGeometry const &box_geo,
                                   Utils::Vector3d const &shift) {
  auto pos = store.position_view();
  auto img = store.image_box_view();
#ifdef ESPRESSO_ROTATION
  auto quat = store.quaternion_view();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto plast = store.position_last_time_step_view();
#endif
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = src_first + static_cast<int>(i);
    auto const d = dst_first + static_cast<int>(i);
    Utils::Vector3d position{pos(s, 0), pos(s, 1), pos(s, 2)};
    Utils::Vector3i image_box{img(s, 0), img(s, 1), img(s, 2)};
    position += shift;
    box_geo.fold_position(position, image_box);
    pos(d, 0) = position[0];
    pos(d, 1) = position[1];
    pos(d, 2) = position[2];
    img(d, 0) = image_box[0];
    img(d, 1) = image_box[1];
    img(d, 2) = image_box[2];
#ifdef ESPRESSO_ROTATION
    quat(d, 0) = quat(s, 0);
    quat(d, 1) = quat(s, 1);
    quat(d, 2) = quat(s, 2);
    quat(d, 3) = quat(s, 3);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    plast(d, 0) = plast(s, 0);
    plast(d, 1) = plast(s, 1);
    plast(d, 2) = plast(s, 2);
#endif
  }
}

/** @brief Columnar FORCE cell-to-cell transfer: accumulate source rows into
 *  destination rows (UPDATE policy: += force and torque). */
static void locl_transfer_force(ParticleStore &store, int src_first,
                                int dst_first, std::size_t n) {
  auto force = store.force_view();
#ifdef ESPRESSO_ROTATION
  auto torque = store.torque_view();
#endif
  for (std::size_t i = 0u; i < n; ++i) {
    auto const s = src_first + static_cast<int>(i);
    auto const d = dst_first + static_cast<int>(i);
    force(d, 0) += force(s, 0);
    force(d, 1) += force(s, 1);
    force(d, 2) += force(s, 2);
#ifdef ESPRESSO_ROTATION
    torque(d, 0) += torque(s, 0);
    torque(d, 1) += torque(s, 1);
    torque(d, 2) += torque(s, 2);
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
  for (std::size_t pl = 0; pl < offset; pl++) {
    auto const n = sizes[pl];
    if (n == 0u) {
      continue;
    }
    auto const [src_first, dst_first] = pairs[pl];
    if (data_parts == GHOSTTRANS_POSITION) {
      locl_transfer_position(store, src_first, dst_first, n, box_geo,
                             ghost_comm.shift);
    } else {
      locl_transfer_force(store, src_first, dst_first, n);
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
                             &ghost_comm.shift);
        serialize_and_reduce(ar_in, p2, data_parts, ReductionPolicy::UPDATE,
                             SerializationDirection::LOAD, box_geo, nullptr);
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
