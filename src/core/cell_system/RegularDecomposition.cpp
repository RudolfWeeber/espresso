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

#include "cell_system/RegularDecomposition.hpp"

#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/MigrationPack.hpp"

#include "communication.hpp"
#include "error_handling/RuntimeErrorStream.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/mpi/cart_comm.hpp>
#include <utils/mpi/sendrecv.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>
#include <boost/range/algorithm/reverse.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <utility>
#include <vector>

int RegularDecomposition::position_to_cell_index(
    Utils::Vector3d const &pos) const {
  Utils::Vector3i cpos;

  for (auto i = 0u; i < 3u; i++) {
    cpos[i] = static_cast<int>(std::floor(pos[i] * inv_cell_size[i])) + 1 -
              cell_offset[i];

    /* particles outside our box. Still take them if
       nonperiodic boundary. We also accept the particle if we are at
       the box boundary, and the particle is within the box. In this case
       the particle belongs here and could otherwise potentially be dismissed
       due to rounding errors. */
    if (cpos[i] < 1) {
      if ((!m_box.periodic(i) or (pos[i] >= m_box.length()[i])) and
          m_local_box.boundary()[2u * i])
        cpos[i] = 1;
      else
        return -1;
    } else if (cpos[i] > cell_grid[i]) {
      if ((!m_box.periodic(i) or (pos[i] < m_box.length()[i])) and
          m_local_box.boundary()[2u * i + 1u])
        cpos[i] = cell_grid[i];
      else
        return -1;
    }
  }

  return get_linear_index(cpos, ghost_cell_grid);
}

void RegularDecomposition::move_if_local(
    std::vector<int> &src, std::vector<int> &rest,
    std::vector<ParticleChange> &modified_cells) {
  auto &staging = *m_migration_staging.store;
  for (auto const staging_row : src) {
    // The particle's home cell is decided by its position, read from the
    // staging store row (no live view needed; the row is not in the main store
    // yet).
    auto target_cell = position_to_cell(staging.position_value(staging_row));

    if (target_cell) {
      // Stage the staging row into the target cell as a row reference: the next
      // store rebuild copies it into a committed row, in
      // [surviving..., staged...] order per cell. The staging row stays valid
      // until CellStructure commits it.
      CellParticleStorage::insert_staged_row(*target_cell, staging,
                                             staging_row);
      modified_cells.emplace_back(ModifiedList{*target_cell});
    } else {
      // Undeliverable this round: keep the staging row for the next round's
      // re-pack (the multi-round intermediate-rank case).
      rest.push_back(staging_row);
    }
  }

  src.clear();
}

void RegularDecomposition::move_left_or_right(std::vector<int> &src,
                                              std::vector<int> &left,
                                              std::vector<int> &right,
                                              int dir) const {
  auto const &staging = *m_migration_staging.store;
  auto const is_open_boundary_left = m_local_box.boundary()[2 * dir] != 0;
  auto const is_open_boundary_right = m_local_box.boundary()[2 * dir + 1] != 0;
  auto const can_move_left = m_box.periodic(dir) or not is_open_boundary_left;
  auto const can_move_right = m_box.periodic(dir) or not is_open_boundary_right;
  auto const my_left = m_local_box.my_left()[dir];
  auto const my_right = m_local_box.my_right()[dir];
  auto keep = src.begin();
  for (auto const staging_row : src) {
    auto const pos = staging.position_value(staging_row)[dir];
    if (m_box.get_mi_coord(pos, my_right, dir) >= 0. and can_move_right) {
      right.push_back(staging_row);
    } else if (m_box.get_mi_coord(pos, my_left, dir) < 0. and can_move_left) {
      left.push_back(staging_row);
    } else {
      // Row stays for a later direction: compact it into the surviving prefix
      // (preserving relative order).
      *keep++ = staging_row;
    }
  }
  src.erase(keep, src.end());
}

void RegularDecomposition::exchange_neighbors(
    std::vector<int> &displaced_rows,
    std::vector<ParticleChange> &modified_cells) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(m_comm);
  auto &staging = *m_migration_staging.store;
  // Per-direction staging-row buckets, packed byte buffers, and the reusable
  // scratch for the received staging rows. Static so the buffers' capacity is
  // reused across resorts.
  static std::vector<int> send_rows_l, send_rows_r, recv_rows;
  static std::vector<char> send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;

  // Unpack a received byte buffer into fresh staging rows and record their
  // indices in `recv_rows` (the received-this-direction set, delivered by
  // move_if_local below).
  auto unpack_received = [&](std::vector<char> const &buffer) {
    if (buffer.empty()) {
      return;
    }
    // unpack_rows reads the row-count header (u64) to learn the count, but the
    // caller must reserve the staging rows first. Peek the count from the
    // row-count header (u64 prefix), reserve that many rows, then unpack.
    std::uint64_t count = 0u;
    std::memcpy(&count, buffer.data(), sizeof(count));
    if (count == 0u) {
      return;
    }
    auto const first_row =
        m_migration_staging.reserve_rows(static_cast<int>(count));
    MigrationPack::unpack_rows(staging, first_row, buffer);
    for (int k = 0; k < static_cast<int>(count); ++k) {
      recv_rows.push_back(first_row + k);
    }
  };

  for (int dir = 0; dir < 3; dir++) {
    /* Single node direction, no action needed. */
    if (Utils::Mpi::cart_get<3>(m_comm).dims[dir] == 1) {
      continue;
      /* In this (common) case left and right neighbors are
         the same, and we need only one communication */
    }
    send_rows_l.clear();
    send_rows_r.clear();
    recv_rows.clear();
    if (Utils::Mpi::cart_get<3>(m_comm).dims[dir] == 2) {
      move_left_or_right(displaced_rows, send_rows_l, send_rows_l, dir);
      MigrationPack::pack_rows(staging, send_rows_l, send_buf_l);

      Utils::Mpi::sendrecv(m_comm, node_neighbors[2 * dir], 0, send_buf_l,
                           node_neighbors[2 * dir], 0, recv_buf_l);

      unpack_received(recv_buf_l);
      recv_buf_l.clear();
    } else {
      using boost::mpi::request;
      using Utils::Mpi::isendrecv;

      move_left_or_right(displaced_rows, send_rows_l, send_rows_r, dir);
      MigrationPack::pack_rows(staging, send_rows_l, send_buf_l);
      MigrationPack::pack_rows(staging, send_rows_r, send_buf_r);

      auto req_l = isendrecv(m_comm, node_neighbors[2 * dir], 0, send_buf_l,
                             node_neighbors[2 * dir], 0, recv_buf_l);
      auto req_r = isendrecv(m_comm, node_neighbors[2 * dir + 1], 0, send_buf_r,
                             node_neighbors[2 * dir + 1], 0, recv_buf_r);

      std::array<request, 4> reqs{{req_l[0], req_l[1], req_r[0], req_r[1]}};
      boost::mpi::wait_all(reqs.begin(), reqs.end());

      unpack_received(recv_buf_l);
      unpack_received(recv_buf_r);
      recv_buf_l.clear();
      recv_buf_r.clear();
    }

    // Deliver the freshly received rows: local ones are staged into their home
    // cell, undeliverable ones are appended back onto `displaced_rows` (kept
    // for the next direction / round). The non-routed rows already sitting in
    // `displaced_rows` are NOT re-examined here.
    move_if_local(recv_rows, displaced_rows, modified_cells);
  }
}

/**
 * @brief Fold coordinates to box and reset the old position.
 */
static void fold_and_reset(Particle &p, BoxGeometry const &box_geo) {
  Utils::Vector3d position = p.pos();
  Utils::Vector3i image_box = p.image_box();
  box_geo.fold_position(position, image_box);
  p.pos() = position;
  p.image_box() = image_box;

  p.pos_at_last_verlet_update() = p.pos();
}

void RegularDecomposition::resort(bool global,
                                  std::vector<ParticleChange> &diff) {
  // Both a mis-celled NON-local particle (routed through the exchange rounds)
  // and a wrong-cell-but-LOCAL particle are copied out of the live store into a
  // staging-store row (stage_row) and then dropped from their cell;
  // `displaced_rows` holds the non-local staging-row indices. The staging store
  // is NOT cleared here: it is reset by CellStructure once it has committed the
  // staged rows (ensure_particle_store_synchronized). Clearing at resort
  // start/end would recycle rows that the deferred commit=false hot path (and
  // the hybrid children sharing this staging store) still reference.
  assert(m_migration_staging && "migration staging store not installed");
  auto &staging = *m_migration_staging.store;
  std::vector<int> displaced_rows;

  for (auto &c : local_cells()) {
    // Iterate the cell's committed rows by RAW position. drop_row marks the row
    // pending-removed (the range keeps its store-row order; no swap-with-back),
    // so we advance unconditionally -- a dropped row stays in the raw range but
    // is skipped by every live iteration until the next rebuild resolves it.
    // fold_and_reset writes through the view into the store column BEFORE the
    // row is copied to staging, so the staged copy carries the folded value.
    // One cached view reused across this cell's rows, REBOUND per position via
    // attach_to_store.
    Particle view;
    auto const offset = c->offset();
    auto const count = c->count();
    for (std::size_t index = 0u; index < count; ++index) {
      auto const live_row = static_cast<int>(offset + index);
      view.attach_to_store(c->store(), live_row);
      fold_and_reset(view, m_box);

      auto target_cell = particle_to_cell(view);

      /* Particle is in place */
      if (target_cell == c) {
        // A pending-removed row reaching here is safe ONLY because removal
        // never moves a particle: fold_and_reset leaves it in the same cell,
        // so particle_to_cell returns c, and we continue without re-staging.
        // The rebuild then drops the pending row. If removal could ever move
        // a particle (changing particle_to_cell), a pending row could enter
        // the migrate/stage branch below and corrupt survivors.
        continue;
      }

      /* Particle is not local */
      if (target_cell == nullptr) {
        diff.emplace_back(RemovedParticle{view.id()});
        // Copy the (folded) live row into the staging store, then drop the row
        // from this cell. The staging row index rides `displaced_rows` through
        // the exchange rounds.
        displaced_rows.push_back(m_migration_staging.stage_row(live_row));
        CellParticleStorage::drop_row(*c, index);
        diff.emplace_back(ModifiedList{*c});
      }
      /* Particle belongs on this node but is in the wrong cell. Copy the
       * (folded) live row into the staging store, drop it from this cell, and
       * stage a reference to that staging row in the target cell -- the next
       * store rebuild copies it into a committed row (surviving-then-staged per
       * cell order contract).
       */
      else {
        auto const staging_row = m_migration_staging.stage_row(live_row);
        CellParticleStorage::drop_row(*c, index);
        diff.emplace_back(ModifiedList{*c});
        CellParticleStorage::insert_staged_row(*target_cell, staging,
                                               staging_row);
        diff.emplace_back(ModifiedList{*target_cell});
      }
    }
  }

  if (global) {
    auto const grid = Utils::Mpi::cart_get<3>(m_comm).dims;
    /* Worst case we need grid - 1 rounds per direction.
     * This correctly implies that if there is only one node,
     * no action should be taken. */
    int rounds_left = grid[0] + grid[1] + grid[2] - 3;
    for (; rounds_left > 0; rounds_left--) {
      exchange_neighbors(displaced_rows, diff);

      auto left_over = boost::mpi::all_reduce(m_comm, displaced_rows.size(),
                                              std::plus<std::size_t>());

      if (left_over == 0) {
        break;
      }
    }
  } else {
    exchange_neighbors(displaced_rows, diff);
  }

  if (not displaced_rows.empty()) {
    auto sort_cell = local_cells()[0];

    for (auto const staging_row : displaced_rows) {
      runtimeErrorMsg() << "Particle " << staging.id(staging_row) << " moved "
                        << "more than one local box length in one timestep";
      CellParticleStorage::insert_staged_row(*sort_cell, staging, staging_row);

      diff.emplace_back(ModifiedList{*sort_cell});
    }
  }

  // The staging store is NOT cleared here: any staged row references (delivered
  // above, or staged by wrong-cell-local moves) must survive until
  // CellStructure commits them (ensure_particle_store_synchronized resets the
  // staging store).
}

void RegularDecomposition::mark_cells() {
  m_local_cells.clear();
  m_ghost_cells.clear();

  int cnt_c = 0;
  for (int o = 0; o < ghost_cell_grid[2]; o++)
    for (int n = 0; n < ghost_cell_grid[1]; n++)
      for (int m = 0; m < ghost_cell_grid[0]; m++) {
        if ((m > 0 && m < ghost_cell_grid[0] - 1 && n > 0 &&
             n < ghost_cell_grid[1] - 1 && o > 0 && o < ghost_cell_grid[2] - 1))
          m_local_cells.push_back(&cells.at(cnt_c++));
        else
          m_ghost_cells.push_back(&cells.at(cnt_c++));
      }
}

void RegularDecomposition::fill_comm_cell_lists(Cell **part_lists,
                                                Utils::Vector3i const &lc,
                                                Utils::Vector3i const &hc) {
  for (int o = lc[0]; o <= hc[0]; o++)
    for (int n = lc[1]; n <= hc[1]; n++)
      for (int m = lc[2]; m <= hc[2]; m++) {
        auto const i = Utils::get_linear_index(o, n, m, ghost_cell_grid);

        *part_lists++ = &(cells.at(i));
      }
}

Utils::Vector3d RegularDecomposition::max_cutoff() const {
  auto dir_max_range = [this](unsigned int i) {
    return std::min(0.5 * m_box.length()[i], m_local_box.length()[i]);
  };

  return {dir_max_range(0u), dir_max_range(1u), dir_max_range(2u)};
}

Utils::Vector3d RegularDecomposition::max_range() const { return cell_size; }
int RegularDecomposition::calc_processor_min_num_cells() const {
  /* the minimal number of cells can be lower if there are at least two nodes
     serving a direction,
     since this also ensures that the cell size is at most half the box
     length. However, if there is only one processor for a direction, there
     have to be at least two cells for this direction. */
  return boost::accumulate(Utils::Mpi::cart_get<3>(m_comm).dims, 1,
                           [](int n_cells, int grid) {
                             return (grid == 1) ? 2 * n_cells : n_cells;
                           });
}

void RegularDecomposition::create_cell_grid(double range) {
  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);

  int n_local_cells;
  auto cell_range = Utils::Vector3d::broadcast(range);
  auto const min_num_cells = calc_processor_min_num_cells();

  if (range <= 0.) {
    /* this is the non-interacting case */
    auto const cells_per_dir =
        static_cast<int>(std::ceil(std::cbrt(min_num_cells)));

    cell_grid = Utils::Vector3i::broadcast(cells_per_dir);
    n_local_cells = Utils::product(cell_grid);
  } else {
    /* Calculate initial cell grid */
    auto const &local_box_l = m_local_box.length();
    auto const volume = Utils::product(local_box_l);
    auto const scale = std::cbrt(RegularDecomposition::max_num_cells / volume);

    for (auto i = 0u; i < 3u; i++) {
      /* this is at least 1 */
      cell_grid[i] = static_cast<int>(std::ceil(local_box_l[i] * scale));
      cell_range[i] = local_box_l[i] / static_cast<double>(cell_grid[i]);

      if (cell_range[i] < range) {
        /* ok, too many cells for this direction, set to minimum */
        cell_grid[i] = static_cast<int>(std::floor(local_box_l[i] / range));
        if (cell_grid[i] < 1) {
          runtimeErrorMsg()
              << "interaction range " << range << " in direction " << i
              << " is larger than the local box size " << local_box_l[i];
          cell_grid[i] = 1;
        }
        cell_range[i] = local_box_l[i] / static_cast<double>(cell_grid[i]);
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the
       coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that.
       */
    for (;;) {
      n_local_cells = Utils::product(cell_grid);

      /* done */
      if (n_local_cells <= RegularDecomposition::max_num_cells)
        break;

      /* find coordinate with the smallest cell range */
      auto min_ind = 0u;
      auto min_size = cell_range[0];

      for (auto i = 1u; i < 3u; ++i) {
        if (cell_grid[i] > 1 and cell_range[i] < min_size) {
          min_ind = i;
          min_size = cell_range[i];
        }
      }

      cell_grid[min_ind]--;
      cell_range[min_ind] = m_local_box.length()[min_ind] / cell_grid[min_ind];
    }

    /* sanity check */
    if (n_local_cells < min_num_cells) {
      runtimeErrorMsg() << "number of cells " << n_local_cells
                        << " is smaller than minimum " << min_num_cells
                        << ": either interaction range is too large for "
                        << "the current skin (range=" << range << ", "
                        << "half_local_box_l=[" << local_box_l / 2. << "]) "
                        << "or min_num_cells too large";
    }
  }

  if (n_local_cells > RegularDecomposition::max_num_cells) {
    runtimeErrorMsg() << "no suitable cell grid found";
  }

  auto const node_pos = cart_info.coords;

  /* now set all dependent variables */
  int new_cells = 1;
  for (auto i = 0u; i < 3u; i++) {
    ghost_cell_grid[i] = cell_grid[i] + 2;
    new_cells *= ghost_cell_grid[i];
    cell_size[i] = m_local_box.length()[i] / static_cast<double>(cell_grid[i]);
    inv_cell_size[i] = 1.0 / cell_size[i];
    cell_offset[i] = node_pos[i] * cell_grid[i];
  }

  /* allocate cell array and cell pointer arrays */
  cells.clear();
  cells.resize(static_cast<unsigned int>(new_cells));
  m_local_cells.resize(n_local_cells);
  m_ghost_cells.resize(new_cells - n_local_cells);
}

template <class K, class Comparator> auto make_flat_set(Comparator &&comp) {
  return boost::container::flat_set<K, std::remove_reference_t<Comparator>>(
      std::forward<Comparator>(comp));
}

void RegularDecomposition::init_cell_interactions() {

  // Note: the global index for physical cells is 0-based.
  // I.e., a global index of -1 refers to a ghost cell.
  auto const halo = Utils::Vector3i{1, 1, 1}; // number of ghost layers
  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);
  // 3D index of the MPI rank in the Cartesian grid of MPI ranks
  auto const &node_pos = cart_info.coords;
  // size of the Cartesian grid of MPI ranks
  auto const &node_grid = ::communicator.node_grid;
  auto const global_halo_offset = hadamard_product(node_pos, cell_grid) - halo;
  // MD cell index of lower halo layer on this MPI rank
  auto const global_size = hadamard_product(node_grid, cell_grid);

  // is a cell at the system boundary in the given coord
  auto const at_boundary = [&global_size](int coord, Utils::Vector3i cell_idx) {
    return (cell_idx[coord] == 0 or cell_idx[coord] == global_size[coord] - 1);
  };

  // For the fully connected feature (cells that don't share at least a corner)
  // only apply if one cell is a ghost cell (i.e. connections across the
  // periodic boundary.
  auto const fcb_is_inner_connection = [&global_size, this](Utils::Vector3i a,
                                                            Utils::Vector3i b) {
    if (fully_connected_boundary()) {
      auto const [fc_normal, fc_dir] = *fully_connected_boundary();
      auto const involves_ghost_cell =
          (a[fc_normal] == -1 or a[fc_normal] == global_size[fc_normal] or
           b[fc_normal] == -1 or b[fc_normal] == global_size[fc_normal]);
      if (not involves_ghost_cell) {
        // check if cells do not share at least a corner
        return std::abs((a - b)[fc_dir]) > 1;
      }
    }
    return false;
  };

  /* Translate a node local index (relative to the origin of the local grid)
   * to a global index. */
  auto global_index =
      [&](Utils::Vector3i const &local_index) -> Utils::Vector3i {
    return (global_halo_offset + local_index);
  };

  /* Linear index in the global cell grid. */
  auto folded_linear_index = [&](Utils::Vector3i const &global_index) {
    auto const folded_index = (global_index + global_size) % global_size;

    return get_linear_index(folded_index, global_size);
  };

  /* Translate a global index into a local one */
  auto local_index =
      [&](Utils::Vector3i const &global_index) -> Utils::Vector3i {
    return (global_index - global_halo_offset);
  };

  // sanity checks
  if (fully_connected_boundary()) {
    auto const [fc_normal, fc_dir] = *fully_connected_boundary();
    if (fc_normal == fc_dir) {
      throw std::domain_error("fully_connected_boundary normal and connection "
                              "coordinates need to differ.");
    }
    if (node_grid[fc_dir] != 1) {
      throw std::runtime_error(
          "The MPI nodegrid must be 1 in the fully connected direction.");
    }
    if (not m_box.periodic(fc_normal)) {
      throw std::runtime_error(
          "The fully connected boundary requires periodicity in the "
          "boundary normal direction.");
    }
  }

  /* We only consider local cells (e.g. not halo cells), which
   * span the range [(1,1,1), cell_grid) in local coordinates. */
  auto const start = global_index(Utils::Vector3i{1, 1, 1});
  auto const end = start + cell_grid;

  bool one_mpi_rank = m_comm.size() == 1;

  /* loop all local cells */
  for (int o = start[2]; o < end[2]; o++)
    for (int n = start[1]; n < end[1]; n++)
      for (int m = start[0]; m < end[0]; m++) {
        /* next-nearest neighbors in every direction */
        Utils::Vector3i lower_index = {m - 1, n - 1, o - 1};
        Utils::Vector3i upper_index = {m + 1, n + 1, o + 1};

        /* In the fully connected case, we consider all cells
        * in the direction as neighbors, not only the nearest ones.
        //         */
        if (fully_connected_boundary()) {
          auto const [fc_boundary, fc_direction] = *fully_connected_boundary();

          // Fully connected is only needed at the box surface
          if (at_boundary(fc_boundary, {m, n, o})) {
            lower_index[fc_direction] = -1;
            upper_index[fc_direction] = global_size[fc_direction];
          }
        }

        /* In non-periodic directions, the halo needs not
         * be considered. */
        for (auto i = 0u; i < 3u; i++) {
          if (not m_box.periodic(i)) {
            lower_index[i] = std::max(0, lower_index[i]);
            upper_index[i] = std::min(global_size[i] - 1, upper_index[i]);
          }
        }

        /* Unique set of neighbors, cells are compared by their linear
         * index in the global cell grid. */
        auto neighbors = make_flat_set<Utils::Vector3i>(
            [&](Utils::Vector3i const &a, Utils::Vector3i const &b) {
              return folded_linear_index(a) < folded_linear_index(b);
            });

        /* Collect neighbors */
        for (int p = lower_index[2]; p <= upper_index[2]; p++)
          for (int q = lower_index[1]; q <= upper_index[1]; q++)
            for (int r = lower_index[0]; r <= upper_index[0]; r++) {
              if (fully_connected_boundary()) {
                // Avoid fully connecting the boundary layer and the
                // next INNER layer
                if (fcb_is_inner_connection({m, n, o}, {r, q, p}))
                  continue;
              }
              neighbors.insert(Utils::Vector3i{r, q, p});
            }

        /* Red-black partition by global index. */
        auto const ind1 = folded_linear_index({m, n, o});

        std::vector<Cell *> red_neighbors;
        std::vector<Cell *> black_neighbors;

        /* If we are running on a single MPI rank, it is not necessary to use
         * ghost cells. Instead of adding a ghost cell as neighbor,
         * we directly connect to the corresponding
         * physical cell across the periodic boundary */
        for (auto &neighbor : neighbors) {
          if (one_mpi_rank) {
            for (auto coord : {0u, 1u, 2u}) {
              if (neighbor[coord] == -1) {
                neighbor[coord] += cell_grid[coord];
              } else if (neighbor[coord] == cell_grid[coord]) {
                neighbor[coord] -= cell_grid[coord];
              }
            }
          }
          auto const ind2 = folded_linear_index(neighbor);
          /* Exclude cell itself */
          if (ind1 == ind2)
            continue;

          auto cell = &cells.at(
              get_linear_index(local_index(neighbor), ghost_cell_grid));

          // Divide red and black neighbors
          if (ind2 > ind1) {
            red_neighbors.push_back(cell);
          } else {
            black_neighbors.push_back(cell);
          }
        }

        // Assign neighbors to the cell
        cells[get_linear_index(local_index({m, n, o}), ghost_cell_grid)]
            .m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

namespace {
/** Revert the order of a communicator: After calling this the
 *  communicator is working in reverted order with exchanged
 *  communication types GHOST_SEND <-> GHOST_RECV.
 */
void revert_comm_order(GhostCommunicator &comm) {
  /* revert order */
  boost::reverse(comm.communications);

  /* exchange SEND/RECV */
  for (auto &c : comm.communications) {
    if (c.type == GHOST_SEND)
      c.type = GHOST_RECV;
    else if (c.type == GHOST_RECV)
      c.type = GHOST_SEND;
    else if (c.type == GHOST_LOCL) {
      boost::reverse(c.part_lists);
    }
  }
}

/** Of every two communication rounds, set the first receivers to prefetch and
 *  poststore
 */
void assign_prefetches(GhostCommunicator &comm) {
  for (auto it = comm.communications.begin(); it != comm.communications.end();
       it += 2) {
    auto next = std::next(it);
    if (it->type == GHOST_RECV && next->type == GHOST_SEND) {
      it->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
      next->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
    }
  }
}

} // namespace

GhostCommunicator RegularDecomposition::prepare_comm() {
  // When running on a single MPI rank, the ghost cells are not needed,
  // as the corresponding physical cell is available on the same MPI rank.
  // The cell neighborships are set up such that cells across the periodic
  // boundaries are connected as neighbors directly.
  if (m_comm.size() == 1)
    return {};

  int dir, lr, i, cnt, n_comm_cells[3];
  Utils::Vector3i lc{}, hc{}, done{};

  auto const comm_info = Utils::Mpi::cart_get<3>(m_comm);
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(m_comm);

  /* calculate number of communications */
  std::size_t num = 0;
  for (dir = 0; dir < 3; dir++) {
    for (lr = 0; lr < 2; lr++) {
      /* No communication for border of non periodic direction */
      if (comm_info.dims[dir] == 1)
        num++;
      else
        num += 2;
    }
  }

  /* prepare communicator */
  auto ghost_comm = GhostCommunicator{m_comm, num};

  /* number of cells to communicate in a direction */
  n_comm_cells[0] = cell_grid[1] * cell_grid[2];
  n_comm_cells[1] = cell_grid[2] * ghost_cell_grid[0];
  n_comm_cells[2] = ghost_cell_grid[0] * ghost_cell_grid[1];

  cnt = 0;
  /* direction loop: x, y, z */
  for (dir = 0; dir < 3; dir++) {
    lc[(dir + 1) % 3] = 1 - done[(dir + 1) % 3];
    lc[(dir + 2) % 3] = 1 - done[(dir + 2) % 3];
    hc[(dir + 1) % 3] = cell_grid[(dir + 1) % 3] + done[(dir + 1) % 3];
    hc[(dir + 2) % 3] = cell_grid[(dir + 2) % 3] + done[(dir + 2) % 3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for (lr = 0; lr < 2; lr++) {
      if (comm_info.dims[dir] == 1) {
        /* just copy cells on a single node */
        ghost_comm.communications[cnt].type = GHOST_LOCL;
        ghost_comm.communications[cnt].node = m_comm.rank();

        /* Buffer has to contain Send and Recv cells -> factor 2 */
        ghost_comm.communications[cnt].part_lists.resize(2 * n_comm_cells[dir]);

        /* fill send ghost_comm cells */
        lc[dir] = hc[dir] = 1 + lr * (cell_grid[dir] - 1);

        fill_comm_cell_lists(ghost_comm.communications[cnt].part_lists.data(),
                             lc, hc);

        /* fill recv ghost_comm cells */
        lc[dir] = hc[dir] = 0 + (1 - lr) * (cell_grid[dir] + 1);

        /* place receive cells after send cells */
        fill_comm_cell_lists(
            &ghost_comm.communications[cnt].part_lists[n_comm_cells[dir]], lc,
            hc);

        cnt++;
      } else {
        /* i: send/recv loop */
        for (i = 0; i < 2; i++) {
          if ((comm_info.coords[dir] + i) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_SEND;
            ghost_comm.communications[cnt].node = node_neighbors[2 * dir + lr];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);
            /* prepare folding of ghost positions */

            lc[dir] = hc[dir] = 1 + lr * (cell_grid[dir] - 1);

            fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
          if ((comm_info.coords[dir] + (1 - i)) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_RECV;
            ghost_comm.communications[cnt].node =
                node_neighbors[2 * dir + (1 - lr)];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);

            lc[dir] = hc[dir] = (1 - lr) * (cell_grid[dir] + 1);

            fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
        }
      }
      done[dir] = 1;
    }
  }

  return ghost_comm;
}

RegularDecomposition::RegularDecomposition(
    boost::mpi::communicator comm, double range, BoxGeometry const &box_geo,
    LocalBox const &local_geo,
    std::optional<std::pair<int, int>> fully_connected)
    : m_comm(std::move(comm)), m_box(box_geo), m_local_box(local_geo),
      m_fully_connected_boundary(std::move(fully_connected)) {

  /* set up new regular decomposition cell structure */
  create_cell_grid(range);

  /* setup cell neighbors */
  init_cell_interactions();

  /* mark local and ghost cells */
  mark_cells();

  /* create communicators */
  m_exchange_ghosts_comm = prepare_comm();
  m_collect_ghost_force_comm = prepare_comm();

  /* collect forces has to be done in reverted order! */
  revert_comm_order(m_collect_ghost_force_comm);

  assign_prefetches(m_exchange_ghosts_comm);
  assign_prefetches(m_collect_ghost_force_comm);
}
