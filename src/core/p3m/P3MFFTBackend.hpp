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

#include <utils/Vector.hpp>

#include <complex>
#include <type_traits>

/**
 * @brief Abstract interface for the P3M reciprocal-space FFT.
 *
 * P3M only ever touches the FFT through this handful of methods: the k-space
 * box geometry accessors and the forward/backward transforms. Concrete
 * backends (heFFTe for the general/distributed case, kokkos-fft for the
 * single-rank case) implement it, so the solver is agnostic to which one is
 * active. The transform buffers are flat, contiguous, no-halo arrays; halo
 * exchange is handled outside the FFT.
 */
template <typename FloatType, class FFTConfig> struct P3MFFTBackend {
  using ComplexType = std::complex<FloatType>;
  /** Real-space scalar type: real for an r2c transform, complex for c2c. */
  using RSpaceScalar =
      std::conditional_t<FFTConfig::use_r2c, FloatType, ComplexType>;

  virtual ~P3MFFTBackend() = default;

  /** @brief Lower-left corner of this rank's k-space box (global coords). */
  virtual Utils::Vector3i ks_local_ld_index() const = 0;
  /** @brief Upper-right corner (exclusive) of this rank's k-space box. */
  virtual Utils::Vector3i ks_local_ur_index() const = 0;
  /** @brief Extent of this rank's k-space box. */
  virtual Utils::Vector3i ks_local_size() const = 0;
  /** @brief Extent of this rank's real-space (no-halo) box. */
  virtual Utils::Vector3i rs_local_size() const = 0;
  /** @brief Forward transform (real/complex charge density to k-space). */
  virtual void forward(RSpaceScalar const *in, ComplexType *out) = 0;
  /** @brief Backward transform (k-space field to real space). */
  virtual void backward(ComplexType const *in, RSpaceScalar *out) = 0;
};
