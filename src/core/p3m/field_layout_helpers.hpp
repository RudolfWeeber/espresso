/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include "for_each_3d.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <iterator>
#include <memory>
#include <new>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace detail {
/** @brief Allocator that leaves trivially-constructible elements uninitialized
 *  on default construction, so a buffer that is about to be fully overwritten
 *  is not first zero-filled. Behaves like std::allocator otherwise.
 */
template <typename T, typename Base = std::allocator<T>>
struct default_init_allocator : Base {
  using traits = std::allocator_traits<Base>;
  template <typename U> struct rebind {
    using other =
        default_init_allocator<U, typename traits::template rebind_alloc<U>>;
  };
  using Base::Base;
  default_init_allocator() = default;
  template <typename U>
  default_init_allocator(default_init_allocator<U> const &) noexcept {}
  template <typename U>
  void construct(U *ptr) noexcept(std::is_nothrow_default_constructible_v<U>) {
    ::new (static_cast<void *>(ptr)) U;
  }
  template <typename U, typename... Args> void construct(U *ptr, Args &&...a) {
    traits::construct(static_cast<Base &>(*this), ptr,
                      std::forward<Args>(a)...);
  }
};
} // namespace detail

// Function to extract a 3D block from the halo field
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename Container>
auto extract_block(Container const &in_array, Utils::Vector3i const &dimensions,
                   Utils::Vector3i const &start, Utils::Vector3i const &stop) {
  // Calculate the size of the block excluding halo regions
  auto const block_dim = stop - start;
  auto const size = static_cast<std::size_t>(Utils::product(block_dim));

  // Output vector to hold the block. The block is written in full below, so
  // skip the zero-initialization (no-init allocator).
  using value_t = typename Container::value_type;
  std::vector<value_t, detail::default_init_allocator<value_t>> out_array(size);

  // Extract the block
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_src = dimensions[2] * dimensions[1];
    auto const lane_src = dimensions[2];
    auto const lane_dst = block_dim[2];
    auto const *const src = in_array.data();
    auto *const dst = out_array.data();
    auto const n_y = stop[1] - start[1];
    // Explicit destination offset per (x, y) lane lets the outer loop run in
    // parallel: each x writes a disjoint set of contiguous lanes, so the copy
    // is thread-safe and bitwise-identical to the serial version.
    _Pragma("omp parallel for") for (int x = start[0]; x < stop[0]; ++x) {
      for (int y = start[1]; y < stop[1]; ++y) {
        auto const offset_src = x * plane_src + y * lane_src + start[2];
        auto const offset_dst =
            (static_cast<std::size_t>(x - start[0]) * n_y + (y - start[1])) *
            lane_dst;
        std::copy_n(src + offset_src, lane_dst, dst + offset_dst);
      }
    }
  } else {
    for_each_3d_lin<output_memory_order>(
        start, stop, [&](Utils::Vector3i const &indices, int out_index) {
          // Compute indices for input and output arrays
          auto const in_index =
              Utils::get_linear_index<memory_order>(indices, dimensions);
          assert(out_index == Utils::get_linear_index<output_memory_order>(
                                  indices - start, block_dim));
          // Copy the value
          out_array[out_index] = in_array[in_index];
        });
  }

  return out_array;
}

/** @brief Pad a 3D matrix with zeros to restore halo regions. */
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename T>
auto pad_with_zeros_discard_imag(std::span<T> cropped_array,
                                 Utils::Vector3i const &cropped_dim,
                                 Utils::Vector3i const &pad_left,
                                 Utils::Vector3i const &pad_right) {

  using value_type = decltype(std::real(std::declval<T>()));
  auto constexpr get_real = [](auto const &v) { return std::real(v); };
  // Calculate dimensions and strides
  auto const padded_dim = cropped_dim + pad_left + pad_right;
  // Output vector to hold the padded field (initialized with zeros)
  std::vector<value_type> padded_array(Utils::product(padded_dim));

  // Fill in the original cropped field into the padded array by chunks
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_dst = padded_dim[2] * padded_dim[1];
    auto const lane_dst = padded_dim[2];
    auto const lane_src = cropped_dim[2];
    auto const base_offset_dst =
        pad_left[0] * plane_dst + pad_left[1] * lane_dst + pad_left[2];
    auto *const dst = padded_array.data();
    auto const *const src = cropped_array.data();
    // Explicit source offset per (x, y) lane lets the outer loop run in
    // parallel: each x reads a disjoint set of lanes and writes a disjoint
    // region of the padded array, so it is thread-safe and bitwise-identical.
    _Pragma("omp parallel for") for (int x = 0; x < cropped_dim[0]; ++x) {
      for (int y = 0; y < cropped_dim[1]; ++y) {
        auto const offset_dst = base_offset_dst + x * plane_dst + y * lane_dst;
        auto const offset_src =
            (static_cast<std::size_t>(x) * cropped_dim[1] + y) * lane_src;
        if constexpr (std::is_floating_point_v<T>) {
          std::copy_n(src + offset_src, lane_src, dst + offset_dst);
        } else {
          std::transform(src + offset_src, src + offset_src + lane_src,
                         dst + offset_dst, get_real);
        }
      }
    }
  } else {
    for_each_3d_lin<memory_order>(
        Utils::Vector3i{0, 0, 0}, cropped_dim,
        [&](Utils::Vector3i const &indices, int in_index) {
          // Compute indices for input and output arrays
          auto const out_index = Utils::get_linear_index<output_memory_order>(
              indices + pad_left, padded_dim);
          // Copy the value
          padded_array[out_index] = std::real(cropped_array[in_index]);
        });
  }

  return padded_array;
}
