/*
 * Copyright (C) 2019-2024 The ESPResSo project
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

/**
 * @file
 * Boundary access methods for @ref walberla::LBWalberlaImpl.
 * This file is included inside the class body of LBWalberlaImpl.
 */

// ---- Boundary access methods (included inside LBWalberlaImpl) ----

public:
std::optional<Utils::Vector3d>
get_node_velocity_at_boundary(Utils::Vector3i const &node,
                              bool consider_ghosts = false) const override {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc or !m_boundary->node_is_boundary(node))
    return std::nullopt;

  return {to_vector3d(m_boundary->get_node_value_at_boundary(node))};
}

bool set_node_velocity_at_boundary(Utils::Vector3i const &node,
                                   Utils::Vector3d const &velocity) override {
  on_boundary_add();
  m_pending_ghost_comm.set(GhostComm::UBB);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc) {
    bc = get_block_and_cell(get_lattice(), node, true);
  }
  if (bc) {
    m_boundary->set_node_value_at_boundary(
        node, to_vector3<FloatType>(velocity), *bc);
  }
  return bc.has_value();
}

std::vector<std::optional<Utils::Vector3d>> get_slice_velocity_at_boundary(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const override {
  std::vector<std::optional<Utils::Vector3d>> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(ci.numCells());

        auto kernel = [&out, this](unsigned const, unsigned const local_index,
                                   Utils::Vector3i const &node) {
          if (m_boundary->node_is_boundary(node)) {
            out[local_index] =
                to_vector3d(m_boundary->get_node_value_at_boundary(node));
          } else {
            out[local_index] = std::nullopt;
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

void set_slice_velocity_at_boundary(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<std::optional<Utils::Vector3d>> const &velocity) override {
  on_boundary_add();
  m_pending_ghost_comm.set(GhostComm::UBB);
  auto const &lattice = get_lattice();
  for_each_block_in_slice(
      lattice, lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(velocity.size() == ci.numCells());

        auto kernel = [&, this](unsigned const, unsigned const local_index,
                                Utils::Vector3i const &node) {
          auto const bc = get_block_and_cell(lattice, node, false);
          assert(bc->block->getAABB() == block.getAABB());
          auto const &opt = velocity[local_index];
          if (opt) {
            m_boundary->set_node_value_at_boundary(
                node, to_vector3<FloatType>(*opt), *bc);
          } else {
            m_boundary->remove_node_from_boundary(node, *bc);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
}

std::optional<Utils::Vector3d>
get_node_boundary_force(Utils::Vector3i const &node) const override {
  auto const bc = get_block_and_cell(get_lattice(), node, true);
  if (!bc or !m_boundary->node_is_boundary(node))
    return std::nullopt;

  return get_node_last_applied_force(node, true);
}

bool remove_node_from_boundary(Utils::Vector3i const &node) override {
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc) {
    bc = get_block_and_cell(get_lattice(), node, true);
  }
  if (bc) {
    m_boundary->remove_node_from_boundary(node, *bc);
  }
  return bc.has_value();
}

std::optional<bool>
get_node_is_boundary(Utils::Vector3i const &node,
                     bool consider_ghosts = false) const override {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  return {m_boundary->node_is_boundary(node)};
}

std::vector<bool>
get_slice_is_boundary(Utils::Vector3i const &lower_corner,
                      Utils::Vector3i const &upper_corner) const override {
  std::vector<bool> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(ci.numCells());

        auto kernel = [&out, this](unsigned const, unsigned const local_index,
                                   Utils::Vector3i const &node) {
          out[local_index] = m_boundary->node_is_boundary(node);
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

void reallocate_ubb_field() override { m_boundary->boundary_update(); }

void on_boundary_add() {
  if (not m_has_boundaries) {
    m_has_boundaries = true;
    setup_streaming_communicator();
  }
  m_has_boundaries = true;
}

void clear_boundaries() override {
  reset_boundary_handling(get_lattice().get_blocks());
  m_pending_ghost_comm.set(GhostComm::UBB);
  ghost_communication();
  m_has_boundaries = false;
  setup_streaming_communicator();
}

void update_boundary_from_shape(std::vector<int> const &raster_flat,
                                std::vector<double> const &data_flat) override {
  on_boundary_add();
  m_pending_ghost_comm.set(GhostComm::UBB);
  auto const &grid_size = get_lattice().get_grid_dimensions();
  auto data = fill_3D_vector_array(data_flat, grid_size);
  set_boundary_from_grid(*m_boundary, get_lattice(), raster_flat, data);
  ghost_communication();
  reallocate_ubb_field();
}

private:
[[nodiscard]] Utils::Vector3i flat_index_to_node(int index) const {
  Utils::Vector3i node({0, 0, 0});
  auto const &grid_size = get_lattice().get_grid_dimensions();
  node[2] = index % grid_size[2];
  int tmp = index / grid_size[2];
  node[1] = tmp % grid_size[1];
  node[0] = tmp / grid_size[1];
  return node;
}

[[nodiscard]] Utils::Vector3i get_neighbor_node(Utils::Vector3i const &node,
                                                int dir) const {
  Utils::Vector3i neighbor({0, 0, 0});
  auto const &grid_size = get_lattice().get_grid_dimensions();
  auto constexpr neighbor_offset = DynamicUBB::neighborOffset;
  for (int i = 0; i < neighbor.size(); i++) {
    neighbor[i] =
        (node[i] - neighbor_offset[i][dir] + grid_size[i]) % grid_size[i];
  }
  return neighbor;
}

public:
[[nodiscard]] Utils::Vector3d get_boundary_force_from_shape(
    std::vector<int> const &raster_flat) const override {
  Utils::Vector3d force({0, 0, 0});
  auto const &grid_size = get_lattice().get_grid_dimensions();
  for (auto &block : *get_lattice().get_blocks()) {
    auto const offset = get_lattice().get_block_corner(block, true);
    auto const &force_field = m_boundary->get_force_vector(&block);
    auto const &index_field = m_boundary->get_index_vector(&block);
    for (int i = 0; i < raster_flat.size(); i++) {
      if (raster_flat[i] != 0) {
        auto node = flat_index_to_node(i);
        if (get_lattice().node_in_local_halo(node)) {
          // shift node to local frame
          node = (node - offset + grid_size) % grid_size;
          for (int j = 0; j < index_field.size(); j++) {
            auto neighbor_node = get_neighbor_node(node, index_field[j].dir);
            if (index_field[j].x == neighbor_node[0] &&
                index_field[j].y == neighbor_node[1] &&
                index_field[j].z == neighbor_node[2]) {
              force[0] += force_field[j].F_0;
              force[1] += force_field[j].F_1;
              force[2] += force_field[j].F_2;
            }
          }
        }
      }
    }
  }
  return zero_centered_to_md(force);
}

// Global boundary force
[[nodiscard]] Utils::Vector3d get_boundary_force() const override {
  Vector3<double> force(0.);
  for (auto &block : *get_lattice().get_blocks()) {
    force += m_boundary->get_total_force(&block);
  }
  return zero_centered_to_md(to_vector3d(force));
}
