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

#include "cell_system/AtomDecomposition.hpp"

#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/MigrationPack.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_to_all.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

void AtomDecomposition::configure_neighbors() {
  std::vector<Cell *> red_neighbors;
  std::vector<Cell *> black_neighbors;

  /* distribute force calculation work  */
  for (int n = 0; n < m_comm.size(); n++) {
    if (m_comm.rank() == n) {
      continue;
    }

    if (n < m_comm.rank()) {
      red_neighbors.push_back(&cells.at(n));
    } else {
      black_neighbors.push_back(&cells.at(n));
    }
  }

  local().m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
}

GhostCommunicator AtomDecomposition::prepare_comm() {
  /* no need for comm for only 1 node */
  if (m_comm.size() == 1) {
    return GhostCommunicator{m_comm, 0};
  }

  auto ghost_comm =
      GhostCommunicator{m_comm, static_cast<std::size_t>(m_comm.size())};
  /* every node has its dedicated comm step */
  for (int n = 0; n < m_comm.size(); n++) {
    ghost_comm.communications[n].part_lists.resize(1);
    ghost_comm.communications[n].part_lists[0] = &(cells.at(n));
    ghost_comm.communications[n].node = n;
  }

  return ghost_comm;
}

void AtomDecomposition::configure_comms() {
  m_exchange_ghosts_comm = prepare_comm();
  m_collect_ghost_force_comm = prepare_comm();

  if (m_comm.size() > 1) {
    for (int n = 0; n < m_comm.size(); n++) {
      /* use the prefetched send buffers. Node 0 transmits first and never
       * prefetches. */
      if (m_comm.rank() == 0 || m_comm.rank() != n) {
        m_exchange_ghosts_comm.communications[n].type = GHOST_BCST;
      } else {
        m_exchange_ghosts_comm.communications[n].type =
            GHOST_BCST | GHOST_PREFETCH;
      }
      m_collect_ghost_force_comm.communications[n].type = GHOST_RDCE;
    }
    /* first round: all nodes except the first one prefetch their send data */
    if (m_comm.rank() != 0) {
      m_exchange_ghosts_comm.communications[0].type |= GHOST_PREFETCH;
    }
  }
}

void AtomDecomposition::mark_cells() {
  m_local_cells.resize(1, std::addressof(local()));
  m_ghost_cells.clear();
  for (int n = 0; n < m_comm.size(); n++) {
    if (n != m_comm.rank()) {
      m_ghost_cells.push_back(std::addressof(cells.at(n)));
    }
  }
}

void AtomDecomposition::resort(bool global_flag,
                               std::vector<ParticleChange> &diff) {
  auto &store = local().store();
  // Fold positions of all committed rows (write-through views).
  for (auto &p : local().particles()) {
    Utils::Vector3d position = p.pos();
    Utils::Vector3i image_box = p.image_box();
    m_box.fold_position(position, image_box);
    p.pos() = position;
    p.image_box() = image_box;

    p.pos_at_last_verlet_update() = p.pos();
  }

  /* Local updates are a NoOp for this decomposition. */
  if (not global_flag) {
    return;
  }

  // Phase 7b flip: a mis-owned particle is copied out of the live store into a
  // staging-store row (stage_row); its staging-row index goes into the target
  // rank's per-rank bucket, which is then packed per-field
  // (MigrationPack::pack_rows) and exchanged as a byte buffer via all_to_all.
  assert(m_migration_staging && "migration staging store not installed");
  auto &staging = *m_migration_staging.store;
  m_migration_staging.clear();

  // Sort displaced particles into per-rank staging-row buckets (iterate the
  // committed rows by position; drop_row removes the row via swap-with-back, so
  // re-examine the swapped-in position).
  std::vector<std::vector<int>> send_rows(m_comm.size());
  for (std::size_t index = 0u; index < local().rows().size();) {
    auto const live_row = local().rows().begin()[index];
    auto const id = store.id(live_row);
    auto const target_node = id_to_rank(id);
    if (target_node != m_comm.rank()) {
      diff.emplace_back(RemovedParticle{id});
      send_rows.at(target_node)
          .push_back(m_migration_staging.stage_row(live_row));
      CellParticleStorage::drop_row(local(), index);
    } else {
      ++index;
    }
  }

  // Pack each rank's bucket into a byte buffer and exchange.
  std::vector<std::vector<char>> send_buf(m_comm.size());
  for (int n = 0; n < m_comm.size(); ++n) {
    MigrationPack::pack_rows(staging, send_rows[static_cast<std::size_t>(n)],
                             send_buf[static_cast<std::size_t>(n)]);
  }
  std::vector<std::vector<char>> recv_buf(m_comm.size());
  boost::mpi::all_to_all(m_comm, send_buf, recv_buf);

  diff.emplace_back(ModifiedList{local().rows()});

  // Unpack the received buffers (in rank order, matching the pre-flip recv_buf
  // iteration) into fresh staging rows, then stage each into the local cell --
  // the pre-flip cell staging path, so the rebuild commits them in the same
  // order.
  for (auto const &buffer : recv_buf) {
    if (buffer.empty()) {
      continue;
    }
    std::uint64_t count = 0u;
    std::memcpy(&count, buffer.data(), sizeof(count));
    if (count == 0u) {
      continue;
    }
    auto const first_row =
        m_migration_staging.reserve_rows(static_cast<int>(count));
    MigrationPack::unpack_rows(staging, first_row, buffer);
    for (int k = 0; k < static_cast<int>(count); ++k) {
      CellParticleStorage::insert_particle(
          local(), m_migration_staging.snapshot_row(first_row + k));
    }
  }

  m_migration_staging.clear();
}

AtomDecomposition::AtomDecomposition(BoxGeometry const &box_geo)
    : m_box(box_geo) {}

AtomDecomposition::AtomDecomposition(boost::mpi::communicator comm,
                                     BoxGeometry const &box_geo)
    : m_comm(std::move(comm)), cells(m_comm.size()), m_box(box_geo) {
  /* create communicators */
  configure_comms();
  /* configure neighbor relations */
  configure_neighbors();
  /* fill local and ghost cell lists */
  mark_cells();
}

Utils::Vector3d AtomDecomposition::max_cutoff() const {
  return Utils::Vector3d::broadcast(std::numeric_limits<double>::infinity());
}

Utils::Vector3d AtomDecomposition::max_range() const { return max_cutoff(); }
