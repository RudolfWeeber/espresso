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
#include "ghosts/particle_packing.hpp"

#include "BoxGeometry.hpp"
#include "system/System.hpp"

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <span>

/** Tag for ghosts communications. */
#define REQ_GHOST_SEND 100

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

  static GhostComm::CommBuf send_buffer, recv_buffer;

  auto const &comm = gcr.mpi_comm;

  for (auto cit = gcr.communications.cbegin(); cit != gcr.communications.cend();
       ++cit) {
    auto const &ghost_comm = *cit;
    int const comm_type = ghost_comm.type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      /* local cell-to-cell transfer: split part_lists into src/dst halves */
      auto const offset = ghost_comm.part_lists.size() / 2;
      for (std::size_t pl = 0; pl < offset; pl++) {
        GhostComm::local_cell_copy(*ghost_comm.part_lists[pl],
                                   *ghost_comm.part_lists[pl + offset],
                                   ghost_comm.shift, box_geo, data_parts);
      }
      continue;
    }

    int const prefetch = ghost_comm.type & GHOST_PREFETCH;
    int const poststore = ghost_comm.type & GHOST_PSTSTORE;
    int const node = ghost_comm.node;

    auto const part_lists_span =
        std::span<ParticleList *const>(ghost_comm.part_lists);

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node, comm.rank())) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch) {
        GhostComm::pack_cells(send_buffer, part_lists_span, ghost_comm.shift,
                              box_geo, data_parts);
      }
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      assert(send_buffer.size() == GhostComm::calc_transmit_size(
                                       part_lists_span, box_geo, data_parts));
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto prefetch_ghost_comm =
          std::find_if(std::next(cit), gcr.communications.cend(),
                       [this_node = comm.rank()](auto const &other_ghost_comm) {
                         return is_prefetchable(other_ghost_comm, this_node);
                       });

      if (prefetch_ghost_comm != gcr.communications.end()) {
        auto const prefetch_lists_span =
            std::span<ParticleList *const>(prefetch_ghost_comm->part_lists);
        GhostComm::pack_cells(send_buffer, prefetch_lists_span,
                              prefetch_ghost_comm->shift, box_geo, data_parts);
      }
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node, comm.rank())) {
      recv_buffer.resize(
          GhostComm::calc_transmit_size(part_lists_span, box_geo, data_parts));
      recv_buffer.bonds().clear();
    }

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
          GhostComm::add_forces(recv_buffer, part_lists_span);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          GhostComm::add_rattle(recv_buffer, part_lists_span);
#endif
        else
          GhostComm::unpack_cells(recv_buffer, part_lists_span, box_geo,
                                  data_parts);
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
        auto const poststore_lists_span =
            std::span<ParticleList *const>(poststore_ghost_comm->part_lists);
        assert(recv_buffer.size() ==
               GhostComm::calc_transmit_size(poststore_lists_span, box_geo,
                                             data_parts));
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          GhostComm::add_forces(recv_buffer, poststore_lists_span);
#ifdef ESPRESSO_BOND_CONSTRAINT
        else if (data_parts == GHOSTTRANS_RATTLE && comm_type != GHOST_RDCE)
          GhostComm::add_rattle(recv_buffer, poststore_lists_span);
#endif
        else
          GhostComm::unpack_cells(recv_buffer, poststore_lists_span, box_geo,
                                  data_parts);
      }
    }
  }
}
