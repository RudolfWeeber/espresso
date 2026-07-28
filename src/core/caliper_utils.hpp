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
 * @file caliper_utils.hpp
 * @brief Zero-overhead Caliper guards for the inactive (no @c CALI_CONFIG)
 * case.
 *
 * When Caliper is compiled in but @c CALI_CONFIG is not set (the common
 * production case), the standard @c CALI_MARK_BEGIN/END and
 * @c CALI_CXX_MARK_FUNCTION macros still pay the full Caliper entry cost
 * (thread-blackboard update, siglock acquisition).  On a hot path with
 * many markers per step this amounts to measurable wall-clock overhead.
 *
 * This header provides guarded replacements:
 * - @c espresso_cali_active() — checked once per process; in steady state this
 *   is a single byte load + branch.
 * - @c ESPRESSO_CALI_MARK_FUNCTION — RAII guard equivalent to
 *   @c CALI_CXX_MARK_FUNCTION but zero-cost when inactive.
 * - @c ESPRESSO_CALI_MARK_BEGIN(name) / @c ESPRESSO_CALI_MARK_END(name)
 * - @c ESPRESSO_CALI_MARK_LOOP_BEGIN / @c ESPRESSO_CALI_MARK_LOOP_ITERATION /
 *   @c ESPRESSO_CALI_MARK_LOOP_END
 *
 * All macros are no-ops when @c ESPRESSO_CALIPER is not defined.
 */

#ifdef ESPRESSO_CALIPER

#include <caliper/cali.h>
#include <cstdlib>

/**
 * @brief Return true if Caliper is configured for this process.
 *
 * Reads @c CALI_CONFIG from the environment exactly once (on the first call)
 * and caches the result.  Subsequent calls pay only one byte load + branch.
 *
 * @c inline ensures a single shared static across all translation units
 * (C++ ODR for inline functions with static locals).
 */
inline bool espresso_cali_active() noexcept {
  static const bool active = (std::getenv("CALI_CONFIG") != nullptr);
  return active;
}

/**
 * @brief RAII region guard: begin on construction, end on destruction.
 *
 * Conditionally calls @c cali_begin_region / @c cali_end_region only when
 * @c espresso_cali_active() is true.  Correctly handles all exit paths.
 */
struct EspressoCaliRegion {
  const char *name_;
  bool active_;

  EspressoCaliRegion(const char *name, bool active)
      : name_(name), active_(active) {
    if (active_)
      cali_begin_region(name_);
  }
  ~EspressoCaliRegion() {
    if (active_)
      cali_end_region(name_);
  }
  EspressoCaliRegion(const EspressoCaliRegion &) = delete;
  EspressoCaliRegion &operator=(const EspressoCaliRegion &) = delete;
};

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

/**
 * @brief Guarded drop-in replacement for @c CALI_CXX_MARK_FUNCTION.
 *
 * Captures @c __func__ at the call site so that the region name is the
 * enclosing function name, not the destructor name.
 */
#define ESPRESSO_CALI_MARK_FUNCTION                                            \
  EspressoCaliRegion CALI_CREATE_VAR_NAME(__espresso_cali_fn, __LINE__)(       \
      __func__, ::espresso_cali_active())

/**
 * @brief Guarded @c CALI_MARK_BEGIN — no-op when inactive.
 */
#define ESPRESSO_CALI_MARK_BEGIN(name)                                         \
  do {                                                                         \
    if (::espresso_cali_active())                                              \
      CALI_MARK_BEGIN(name);                                                   \
  } while (false)

/**
 * @brief Guarded @c CALI_MARK_END — no-op when inactive.
 */
#define ESPRESSO_CALI_MARK_END(name)                                           \
  do {                                                                         \
    if (::espresso_cali_active())                                              \
      CALI_MARK_END(name);                                                     \
  } while (false)

/**
 * @brief Guarded loop annotation: begin the loop region.
 *
 * Creates a per-loop RAII guard that emits the "loop" annotation only when
 * active.  The loop_id identifier must be unique within the enclosing scope.
 * End the loop region by calling @c ESPRESSO_CALI_MARK_LOOP_END(loop_id)
 * before the loop ends.
 *
 * Unlike @c CALI_CXX_MARK_LOOP_BEGIN, which creates a @c cali::Loop with
 * iteration-level tracking, this guard uses a plain region for the loop and
 * omits per-iteration attributes (iteration counts) when inactive.  When
 * active, both the region and the per-loop object are created so that
 * @c ESPRESSO_CALI_MARK_LOOP_ITERATION can use the loop object.
 */
#define ESPRESSO_CALI_MARK_LOOP_BEGIN(loop_id, name)                           \
  bool const __espresso_cali_loop_active_##loop_id = ::espresso_cali_active(); \
  CALI_CXX_MARK_LOOP_BEGIN(loop_id, name)

/**
 * @brief Guarded per-step loop-iteration annotation.
 *
 * No-op when loop was marked inactive by @c ESPRESSO_CALI_MARK_LOOP_BEGIN.
 */
#define ESPRESSO_CALI_MARK_LOOP_ITERATION(loop_id, iter)                       \
  do {                                                                         \
    if (__espresso_cali_loop_active_##loop_id)                                 \
      CALI_CXX_MARK_LOOP_ITERATION(loop_id, iter);                             \
  } while (false)

/**
 * @brief End the loop region opened by @c ESPRESSO_CALI_MARK_LOOP_BEGIN.
 */
#define ESPRESSO_CALI_MARK_LOOP_END(loop_id)                                   \
  do {                                                                         \
    if (__espresso_cali_loop_active_##loop_id)                                 \
      CALI_CXX_MARK_LOOP_END(loop_id);                                         \
  } while (false)

// NOLINTEND(cppcoreguidelines-macro-usage)

#endif // ESPRESSO_CALIPER
