/*
 * Copyright (C) 2026 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ESPRESSO_H5MD

#include "Shape.hpp"

#include "core/BoxGeometry.hpp"
#include "core/LocalBox.hpp"
#include "core/io/voxel_mask/VoxelMaskFile.hpp"
#include "core/system/System.hpp"

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Shapes {

namespace detail {

/** Contiguous run of voxels along one axis, mapped to the mask. */
struct VoxelMaskSegment {
  int region_start;
  int mask_start;
  int length;
};

/** Split the stored region along @p axis into runs of mask indices. */
inline auto voxel_mask_segments(::Shapes::VoxelMask::Geometry const &g,
                                ::Shapes::VoxelMask::Region const &region,
                                int axis) {
  std::vector<VoxelMaskSegment> segments;
  for (int i = 0; i < region.extent[axis]; ++i) {
    auto const m =
        ::Shapes::VoxelMask::mask_index(g, axis, region.lower[axis] + i);
    if (m < 0) {
      continue;
    }
    if (not segments.empty()) {
      auto &last = segments.back();
      if (last.region_start + last.length == i and
          last.mask_start + last.length == m) {
        ++last.length;
        continue;
      }
    }
    segments.push_back({i, m, 1});
  }
  return segments;
}

/** Local geometry of the last load. */
struct VoxelMaskLoadState {
  Utils::Vector3d my_left;
  Utils::Vector3d my_right;
  Utils::Vector3d box_l;
  std::array<bool, 3> periodic;
  bool operator==(VoxelMaskLoadState const &) const = default;
};

/**
 * @brief Make a loader that reads the stored region of this rank.
 * The loader is collective. It re-reads the file only if the local box
 * or the box geometry changed on any rank since the last load.
 */
inline auto make_voxel_mask_loader(std::string path, std::string dataset,
                                   boost::mpi::communicator comm) {
  auto last = std::make_shared<std::optional<VoxelMaskLoadState>>();
  return [path = std::move(path), dataset = std::move(dataset),
          comm = std::move(comm), last](::Shapes::VoxelMask &mask) {
    auto const &system = ::System::get_system();
    auto const &box_geo = *system.box_geo;
    auto const &local_geo = *system.local_geo;
    VoxelMaskLoadState const state{
        local_geo.my_left(),
        local_geo.my_right(),
        box_geo.length(),
        {{box_geo.periodic(0u), box_geo.periodic(1u), box_geo.periodic(2u)}}};
    auto const changed = not last->has_value() or * *last != state;
    if (not boost::mpi::all_reduce(comm, changed, std::logical_or<>())) {
      return;
    }
    auto g = mask.geometry();
    g.box_l = state.box_l;
    g.periodic = state.periodic;
    mask.set_geometry(g);
    last->reset();

    IO::VoxelMaskFile const file(path, dataset, comm);
    if (file.shape() != g.mask_shape) {
      throw std::runtime_error("The shape of dataset '" + dataset +
                               "' in file '" + path + "' has changed");
    }
    auto const region =
        ::Shapes::VoxelMask::stored_region(g, state.my_left, state.my_right);
    std::array<std::vector<VoxelMaskSegment>, 3> const segments{
        {voxel_mask_segments(g, region, 0), voxel_mask_segments(g, region, 1),
         voxel_mask_segments(g, region, 2)}};
    auto const &extent = region.extent;
    std::vector<std::uint8_t> voxels(static_cast<std::size_t>(extent[0]) *
                                         static_cast<std::size_t>(extent[1]) *
                                         static_cast<std::size_t>(extent[2]),
                                     std::uint8_t{0});
    for (auto const &sx : segments[0]) {
      for (auto const &sy : segments[1]) {
        for (auto const &sz : segments[2]) {
          auto const block =
              file.read({sx.mask_start, sy.mask_start, sz.mask_start},
                        {sx.length, sy.length, sz.length});
          std::size_t n = 0;
          for (int i = 0; i < sx.length; ++i) {
            for (int j = 0; j < sy.length; ++j) {
              auto const offset =
                  (static_cast<std::size_t>(sx.region_start + i) *
                       static_cast<std::size_t>(extent[1]) +
                   static_cast<std::size_t>(sy.region_start + j)) *
                      static_cast<std::size_t>(extent[2]) +
                  static_cast<std::size_t>(sz.region_start);
              std::copy_n(block.begin() + static_cast<std::ptrdiff_t>(n),
                          sz.length,
                          voxels.begin() + static_cast<std::ptrdiff_t>(offset));
              n += static_cast<std::size_t>(sz.length);
            }
          }
        }
      }
    }
    mask.set_local_data(region, voxels);
    *last = state;
  };
}

} // namespace detail

class VoxelMask : public Shape {
public:
  VoxelMask() : m_shape(std::make_shared<::Shapes::VoxelMask>()) {
    add_parameters(
        {{"file", AutoParameter::read_only, [this]() { return m_file; }},
         {"dataset", AutoParameter::read_only, [this]() { return m_dataset; }},
         {"voxel_size", AutoParameter::read_only,
          [this]() { return m_shape->geometry().voxel_size; }},
         {"origin", AutoParameter::read_only,
          [this]() { return m_shape->geometry().origin; }},
         {"halo", AutoParameter::read_only,
          [this]() { return m_shape->geometry().halo; }},
         {"mask_shape", AutoParameter::read_only,
          [this]() { return m_shape->geometry().mask_shape; }}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_shape; }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([&]() {
      m_file = get_value<std::string>(params, "file");
      m_dataset = get_value<std::string>(params, "dataset");
      ::Shapes::VoxelMask::Geometry g;
      g.voxel_size = get_value<double>(params, "voxel_size");
      g.origin = get_value_or<Utils::Vector3d>(params, "origin", {});
      g.halo = get_value_or<int>(params, "halo", 2);
      if (not(g.voxel_size > 0.)) {
        throw std::domain_error("Parameter 'voxel_size' must be > 0");
      }
      if (g.halo < 1) {
        throw std::domain_error("Parameter 'halo' must be >= 1");
      }
      auto const &comm = context()->get_comm();
      g.mask_shape = IO::VoxelMaskFile(m_file, m_dataset, comm).shape();
      auto const &box_geo = *::System::get_system().box_geo;
      g.box_l = box_geo.length();
      for (unsigned int a = 0u; a < 3u; ++a) {
        g.periodic[a] = box_geo.periodic(a);
      }
      m_shape->set_geometry(g);
      m_shape->set_loader(
          detail::make_voxel_mask_loader(m_file, m_dataset, comm));
      m_shape->on_local_geometry_change();
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    context()->parallel_try_catch(
        [this]() { m_shape->on_local_geometry_change(); });
    auto const &comm = context()->get_comm();

    if (name == "calc_distance" or name == "is_inside") {
      context()->parallel_try_catch([&]() {
        if (not params.contains("position")) {
          throw std::runtime_error("Parameter 'position' is missing");
        }
      });
      auto const &box_geo = *::System::get_system().box_geo;
      auto const pos = box_geo.folded_position(
          get_value<Utils::Vector3d>(params.at("position")));
      auto const owner = owner_rank(pos);
      // result: distance, vector, inside
      std::array<double, 5> result{};
      context()->parallel_try_catch([&]() {
        if (comm.rank() != owner) {
          return;
        }
        if (name == "calc_distance") {
          double dist;
          Utils::Vector3d vec;
          m_shape->calculate_dist(pos, dist, vec);
          result = {dist, vec[0], vec[1], vec[2], 0.};
        } else {
          result[4] = m_shape->is_inside(pos) ? 1. : 0.;
        }
      });
      if (owner != 0) {
        if (comm.rank() == owner) {
          comm.send(0, 42, result.data(), static_cast<int>(result.size()));
        } else if (comm.rank() == 0) {
          comm.recv(owner, 42, result.data(), static_cast<int>(result.size()));
        }
      }
      if (comm.rank() != 0) {
        return {};
      }
      if (name == "calc_distance") {
        return std::vector<Variant>{
            result[0], Utils::Vector3d{result[1], result[2], result[3]}};
      }
      return {result[4] != 0.};
    }

    if (name == "rasterize") {
      context()->parallel_try_catch([]() {
        throw std::runtime_error("VoxelMask cannot be rasterized globally");
      });
      return {};
    }

    if (name == "get_local_voxel_counts") {
      return gather_on_head(m_shape->local_voxel_count());
    }

    if (name == "get_local_memory_bytes") {
      return gather_on_head(m_shape->local_memory_bytes());
    }

    return Shape::do_call_method(name, params);
  }

private:
  /** @brief Rank whose local box contains a folded position. */
  int owner_rank(Utils::Vector3d const &folded_pos) const {
    auto const &comm = context()->get_comm();
    auto const &system = ::System::get_system();
    auto const &box_geo = *system.box_geo;
    auto const &local_geo = *system.local_geo;
    auto folded = folded_pos;
    auto is_local = true;
    for (unsigned int a = 0u; a < 3u; ++a) {
      if (not box_geo.periodic(a)) {
        // positions beyond a non-periodic boundary go to the boundary rank
        folded[a] = std::clamp(folded[a], 0., box_geo.length()[a]);
      }
      // boundaries are inclusive to close floating-point gaps between ranks
      is_local &= local_geo.my_left()[a] <= folded[a] and
                  folded[a] <= local_geo.my_right()[a];
    }
    auto const owner = boost::mpi::all_reduce(
        comm, is_local ? comm.rank() : comm.size(), boost::mpi::minimum<int>());
    context()->parallel_try_catch([&]() {
      if (owner == comm.size()) {
        throw std::runtime_error("No MPI rank owns the position");
      }
    });
    return owner;
  }

  Variant gather_on_head(std::size_t value) const {
    auto const &comm = context()->get_comm();
    std::vector<std::size_t> values;
    boost::mpi::gather(comm, value, values, 0);
    if (comm.rank() != 0) {
      return {};
    }
    return make_vector_of_variants(values);
  }

  std::shared_ptr<::Shapes::VoxelMask> m_shape;
  std::string m_file;
  std::string m_dataset;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif // ESPRESSO_H5MD
