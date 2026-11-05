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
 *  Reusable particle packing/unpacking for ghost communications.
 */

#include "particle_packing.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ghosts.hpp"
#include "ghosts/threaded_cell_loop.hpp"

#include <utils/attributes.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

namespace GhostComm {

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
 * Forced inline: it runs once per particle in every ghost update, and GCC
 * otherwise keeps it as an out-of-line call once the threaded loops are
 * added to the packing functions.
 */
template <class Archive>
ESPRESSO_ATTR_ALWAYS_INLINE static inline void
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
    if (direction == SerializationDirection::SAVE and ghost_shift != nullptr) {
      /* ok, this is not nice, but perhaps fast */
      auto pos = p.pos() + *ghost_shift;
      auto img = p.image_box();
      box_geo.fold_position(pos, img);
      ar & pos;
      ar & img;
    } else {
      ar & p.pos();
      ar & p.image_box();
    }
#ifdef ESPRESSO_BOND_CONSTRAINT
    ar & p.pos_last_time_step();
#endif
  }
  // Wire-symmetry: pack and unpack use the same data_parts per exchange,
  // so the layout is always symmetric; QUAT follows POSITION in the stream
  // when both are set (position push), and TORQUE follows FORCE when both
  // are set (force reduce).
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_QUAT) {
    ar & p.quat();
  }
#endif
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    ar & p.v();
#ifdef ESPRESSO_ROTATION
    ar & p.omega();
#endif
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d force;
      ar & force;
      p.force() += force;
    } else {
      ar & p.force();
    }
  }
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_TORQUE) {
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d torque;
      ar & torque;
      p.torque() += torque;
    } else {
      ar & p.torque();
    }
  }
#endif
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
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (data_parts & GHOSTTRANS_DIPFLD) {
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d dip_fld;
      ar & dip_fld;
      p.dip_fld() += dip_fld;
    } else {
      ar & p.dip_fld();
    }
  }
#endif
}

static void prepare_ghost_cell(ParticleList *cell, std::size_t size) {
  /* Adapt size */
  cell->resize(size);

  /* Mark particles as ghosts */
  for (auto &p : *cell) {
    p.set_ghost(true);
  }
}

/**
 * @brief Offsets of the cells in a buffer, in particles.
 * Entry i is the number of particles in the cells before cell i, the last
 * entry is the total number of particles.
 */
static auto cell_offsets(std::span<ParticleList *const> cells) {
  std::vector<std::size_t> offsets(cells.size() + 1u, 0u);
  for (std::size_t i = 0u; i < cells.size(); ++i) {
    offsets[i + 1u] = offsets[i] + cells[i]->size();
  }
  return offsets;
}

std::size_t calc_transmit_size(BoxGeometry const &box_geo,
                               unsigned data_parts) {
  SerializationSizeCalculator sizeof_archive;
  Particle p{};
  serialize_and_reduce(sizeof_archive, p, data_parts, ReductionPolicy::MOVE,
                       SerializationDirection::SAVE, box_geo, nullptr);
  return sizeof_archive.size();
}

std::size_t calc_transmit_size(std::span<ParticleList *const> cells,
                               BoxGeometry const &box_geo,
                               unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(unsigned int) * cells.size();

  auto const n_part = std::accumulate(
      cells.begin(), cells.end(), std::size_t{0},
      [](std::size_t sum, auto part_list) { return sum + part_list->size(); });

  return n_part * calc_transmit_size(box_geo, data_parts);
}

void pack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                unsigned data_parts) {
  /* reallocate send buffer */
  buf.resize(calc_transmit_size(cells, box_geo, data_parts));
  buf.bonds().clear();

  if (not(data_parts & (GHOSTTRANS_PARTNUM | GHOSTTRANS_BONDS)) and
      cells.size() > 1u and threads_available()) {
    auto const offsets = cell_offsets(cells);
    if (use_threaded_loop(offsets.back())) {
      // every particle takes the same number of bytes, so each cell has a
      // fixed slice of the buffer and the cells can be packed in parallel
      auto const part_size = calc_transmit_size(box_geo, data_parts);
      auto const buffer = buf.make_span();
      threaded_loop("ghost_pack_cells", cells.size(), [&](std::size_t i) {
        auto ar = Utils::MemcpyOArchive{buffer.subspan(
            offsets[i] * part_size, cells[i]->size() * part_size)};
        for (auto &p : *cells[i]) {
          serialize_and_reduce(ar, p, data_parts, ReductionPolicy::MOVE,
                               SerializationDirection::SAVE, box_geo, &shift);
        }
        assert(ar.bytes_written() == cells[i]->size() * part_size);
      });
      return;
    }
  }

  auto archiver = Utils::MemcpyOArchive{buf.make_span()};

  /* Construct archive that pushes back to the bond buffer */
  namespace io = boost::iostreams;
  io::stream<io::back_insert_device<std::vector<char>>> os{
      io::back_inserter(buf.bonds())};
  boost::archive::binary_oarchive bond_archiver{os};

  /* put in data */
  for (auto part_list : cells) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      assert(part_list->size() <= std::numeric_limits<unsigned int>::max());
      auto np = static_cast<unsigned int>(part_list->size());
      archiver << np;
    } else {
      for (auto &p : *part_list) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::SAVE, box_geo, &shift);
        if (data_parts & GHOSTTRANS_BONDS) {
          bond_archiver << p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_written() == buf.size());
}

void unpack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                  BoxGeometry const &box_geo, unsigned data_parts) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};

  if (data_parts & GHOSTTRANS_PARTNUM) {
    for (auto part_list : cells) {
      unsigned int np;
      archiver >> np;
      prepare_ghost_cell(part_list, np);
    }
  } else {
    auto unpacked = false;
    if (not(data_parts & GHOSTTRANS_BONDS) and cells.size() > 1u and
        threads_available()) {
      auto const offsets = cell_offsets(cells);
      if (use_threaded_loop(offsets.back())) {
        // each ghost cell is filled by exactly one message, so the cells can
        // be unpacked in parallel from their fixed slices of the buffer
        auto const part_size = calc_transmit_size(box_geo, data_parts);
        auto const buffer = buf.make_span();
        threaded_loop("ghost_unpack_cells", cells.size(), [&](std::size_t i) {
          auto ar = Utils::MemcpyIArchive{buffer.subspan(
              offsets[i] * part_size, cells[i]->size() * part_size)};
          for (auto &p : *cells[i]) {
            serialize_and_reduce(ar, p, data_parts, ReductionPolicy::MOVE,
                                 SerializationDirection::LOAD, box_geo,
                                 nullptr);
          }
          assert(ar.bytes_read() == cells[i]->size() * part_size);
        });
        archiver.skip(offsets.back() * part_size);
        unpacked = true;
      }
    }
    if (not unpacked) {
      for (auto part_list : cells) {
        for (auto &p : *part_list) {
          serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                               SerializationDirection::LOAD, box_geo, nullptr);
        }
      }
    }
    if (data_parts & GHOSTTRANS_BONDS) {
      namespace io = boost::iostreams;
      io::stream<io::array_source> bond_stream(
          io::array_source{buf.bonds().data(), buf.bonds().size()});
      boost::archive::binary_iarchive bond_archiver(bond_stream);

      for (auto part_list : cells) {
        for (auto &p : *part_list) {
          bond_archiver >> p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_read() == buf.size());

  buf.bonds().clear();
}

template <class Archive>
static void add_force_and_torque(Archive &ar, Particle &part,
                                 unsigned data_parts) {
  if (data_parts & GHOSTTRANS_FORCE) {
    Utils::Vector3d force;
    ar >> force;
    part.force() += force;
  }
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_TORQUE) {
    Utils::Vector3d torque;
    ar >> torque;
    part.torque() += torque;
  }
#endif
}

void add_forces(CommBuf &buf, std::span<ParticleList *const> cells,
                unsigned data_parts) {
  if (cells.size() > 1u and threads_available()) {
    auto const offsets = cell_offsets(cells);
    if (use_threaded_loop(offsets.back())) {
      // a local cell can occur several times in one buffer (images in several
      // directions): one thread adds all occurrences of a cell, in buffer
      // order, which gives the same summation order as the serial loop
      std::size_t part_size = 0u;
      if (data_parts & GHOSTTRANS_FORCE) {
        part_size += sizeof(Utils::Vector3d);
      }
#ifdef ESPRESSO_ROTATION
      if (data_parts & GHOSTTRANS_TORQUE) {
        part_size += sizeof(Utils::Vector3d);
      }
#endif
      assert(part_size * offsets.back() == buf.size());
      auto const buffer = buf.make_span();
      auto const groups = group_by_target_cell(
          cells.size(), [&cells](std::size_t i) { return cells[i]; });
      threaded_loop("ghost_add_forces", groups.size(), [&](std::size_t g) {
        for (auto k = groups.begin[g]; k < groups.begin[g + 1u]; ++k) {
          auto const i = groups.order[k];
          auto ar = Utils::MemcpyIArchive{buffer.subspan(
              offsets[i] * part_size, cells[i]->size() * part_size)};
          for (Particle &part : *cells[i]) {
            add_force_and_torque(ar, part, data_parts);
          }
        }
      });
      return;
    }
  }
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto &part_list : cells) {
    for (Particle &part : *part_list) {
      add_force_and_torque(archiver, part, data_parts);
    }
  }
}

#ifdef ESPRESSO_BOND_CONSTRAINT
void add_rattle(CommBuf &buf, std::span<ParticleList *const> cells) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto &part_list : cells) {
    for (Particle &part : *part_list) {
      ParticleRattle pr;
      archiver >> pr;
      part.rattle_params() += pr;
    }
  }
}
#endif

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void add_dip_fld(CommBuf &buf, std::span<ParticleList *const> cells) {
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto &part_list : cells) {
    for (Particle &part : *part_list) {
      Utils::Vector3d dip_fld;
      archiver >> dip_fld;
      part.dip_fld() += dip_fld;
    }
  }
}
#endif

void local_cell_copy(ParticleList &src, ParticleList &dst,
                     Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                     unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM) {
    prepare_ghost_cell(&dst, src.size());
  } else {
    assert(src.size() == dst.size());
    CommBuf buffer;
    buffer.resize(calc_transmit_size(box_geo, data_parts));

    for (std::size_t i = 0; i < src.size(); i++) {
      auto ar_out = Utils::MemcpyOArchive{buffer.make_span()};
      auto ar_in = Utils::MemcpyIArchive{buffer.make_span()};
      auto &p1 = src.begin()[i];
      auto &p2 = dst.begin()[i];
      serialize_and_reduce(ar_out, p1, data_parts, ReductionPolicy::UPDATE,
                           SerializationDirection::SAVE, box_geo, &shift);
      serialize_and_reduce(ar_in, p2, data_parts, ReductionPolicy::UPDATE,
                           SerializationDirection::LOAD, box_geo, nullptr);
      if (data_parts & GHOSTTRANS_BONDS) {
        p2.bonds() = p1.bonds();
      }
    }
  }
}

} // namespace GhostComm
