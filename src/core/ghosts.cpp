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
#include "cell_system/ParticleListOperations.hpp"
#include "system/System.hpp"

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
#include <functional>
#include <iterator>
#include <limits>
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
      Utils::Vector3d position = p.pos() + *ghost_shift;
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

static void prepare_send_buffer(CommBuf &send_buffer,
                                GhostCommunication const &ghost_comm,
                                BoxGeometry const &box_geo,
                                unsigned int data_parts) {

  /* reallocate send buffer */
  send_buffer.resize(calc_transmit_size(ghost_comm, box_geo, data_parts));
  send_buffer.bonds().clear();

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

static void put_recv_buffer(CommBuf &recv_buffer,
                            GhostCommunication const &ghost_comm,
                            BoxGeometry const &box_geo,
                            unsigned int data_parts) {
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

static void add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                        const GhostCommunication &ghost_comm) {
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

static void cell_cell_transfer(GhostCommunication const &ghost_comm,
                               BoxGeometry const &box_geo,
                               unsigned int data_parts) {
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
