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

/**
 * @file
 * Duplicated Kokkos ScatterView without a private replica when the
 * execution space runs on a single thread.
 *
 * A duplicated ScatterView keeps one replica of its destination view per
 * thread; kernels accumulate into the replica of their thread, and
 * @c contribute() adds all replicas to the destination. With one thread
 * the single replica only adds memory traffic: it is zeroed, accumulated
 * into, and then added to the (zeroed) destination, i.e. one more array of
 * the destination's size per step that has to stay in cache.
 *
 * @ref make_scatter_view instead lets the single replica alias the
 * destination, so kernels write into the destination directly. Kokkos
 * supports such an aliased replica: @c contribute() skips a replica that is
 * the destination, and @c reset_except(destination) does not zero it. The
 * ScatterView type is unchanged, so the choice is made at run time and
 * the kernels are not affected. The summation order per value is the same
 * as with a separate replica (0 + a_1 + a_2 + ...), so the results are
 * bitwise identical.
 *
 * Setting the replica relies on the protected replica member of the
 * duplicated ScatterView, for which Kokkos (5.0) has no public interface;
 * @c scatter_view_single_thread_test checks that the alias is set up.
 */

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <optional>
#include <type_traits>

namespace ScatterViewSingleThread {
/** @brief Duplicated ScatterView with right layout (one replica per thread
 *  in the leading dimension).
 */
template <class> struct is_duplicated_layout_right : std::false_type {};
template <class DataType, class Device, class Op, class Contribution>
struct is_duplicated_layout_right<Kokkos::Experimental::ScatterView<
    DataType, Kokkos::LayoutRight, Device, Op,
    Kokkos::Experimental::ScatterDuplicated, Contribution>> : std::true_type {
};

/** @brief Access to the protected members of a duplicated ScatterView. */
template <class ScatterView> struct Aliased : ScatterView {
  /**
   * @brief ScatterView whose only replica is @p dest, or nothing if the
   * execution space has more than one thread.
   */
  template <class View>
  static std::optional<ScatterView> make(View const &dest) {
    using internal_view_type = typename ScatterView::internal_view_type;
    static_assert(std::is_same_v<typename View::array_layout,
                                 Kokkos::LayoutRight>);
    static_assert(internal_view_type::rank() == View::rank() + 1u);
    static_assert(std::is_same_v<typename View::non_const_value_type,
                                 typename internal_view_type::value_type>);
    // the default constructor allocates no replicas
    Aliased scatter_view{};
    if (scatter_view.unique_token.size() != 1) {
      return std::nullopt;
    }
    // unmanaged view of dest with a leading replica extent of one
    if constexpr (View::rank_dynamic == 0u) {
      scatter_view.internal_view = internal_view_type(dest.data(), 1);
    } else {
      static_assert(View::rank_dynamic == 1u);
      scatter_view.internal_view =
          internal_view_type(dest.data(), 1, dest.extent(0));
    }
    return static_cast<ScatterView const &>(scatter_view);
  }

  /** @brief Whether the replicas of @p scatter_view are @p dest itself. */
  template <class View>
  static bool is_aliased(ScatterView const &scatter_view, View const &dest) {
    return (scatter_view.*(&Aliased::internal_view)).data() == dest.data();
  }
};

/**
 * @brief Create a ScatterView of type @p ScatterView for @p dest.
 *
 * If @p ScatterView is duplicated and its execution space has one thread,
 * the single replica is @p dest itself, otherwise this is the usual
 * ScatterView with its own replicas. Requirements:
 * - @p dest must be zero when the ScatterView is created, since the
 *   aliased replica is not zeroed by the constructor as a separate replica
 *   would be (allocate @p dest with initialization).
 * - In the aliased case the ScatterView, and every copy of it, holds an
 *   unmanaged view of @p dest: it must not outlive or be used after a
 *   reallocation of @p dest; create a new one instead.
 * - Reset the replicas with @c reset_except(dest) after zeroing @p dest.
 */
template <class ScatterView, class View>
ScatterView make_scatter_view(View const &dest) {
  if constexpr (is_duplicated_layout_right<ScatterView>::value) {
    if (auto scatter_view = Aliased<ScatterView>::make(dest)) {
      return *scatter_view;
    }
  }
  return ScatterView(dest);
}

/** @brief Whether @p scatter_view has no replica other than @p dest. */
template <class ScatterView, class View>
bool writes_directly(ScatterView const &scatter_view, View const &dest) {
  if constexpr (is_duplicated_layout_right<ScatterView>::value) {
    return Aliased<ScatterView>::is_aliased(scatter_view, dest);
  }
  return true;
}
} // namespace ScatterViewSingleThread
