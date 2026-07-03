/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

#include <cassert>
#include <cmath>
#include <cstddef>

/**
 * @brief Write-through reference to one particle's 3-vector stored in a
 * component-major (LayoutLeft) column.
 *
 * Component @c j of the referenced vector lives at <tt>base[j * stride]</tt>.
 * The proxy is a cheap value type (pointer + stride); copying it copies the
 * reference, not the data. Reads convert to @ref Utils::Vector3d by value.
 *
 * Note: because @ref Utils::Vector arithmetic operators are function
 * templates, a VectorReference does not participate in template argument
 * deduction — expressions like <tt>proxy - vector</tt> require an explicit
 * conversion: <tt>Utils::Vector3d(proxy) - vector</tt>.
 */
class VectorReference {
  double *m_base;
  std::size_t m_stride;

public:
  VectorReference(double *base, std::size_t stride)
      : m_base(base), m_stride(stride) {
    assert(base != nullptr);
  }

  operator Utils::Vector3d() const {
    return {m_base[0u], m_base[m_stride], m_base[2u * m_stride]};
  }

  VectorReference &operator=(Utils::Vector3d const &value) {
    m_base[0u] = value[0u];
    m_base[m_stride] = value[1u];
    m_base[2u * m_stride] = value[2u];
    return *this;
  }

  VectorReference &operator+=(Utils::Vector3d const &value) {
    m_base[0u] += value[0u];
    m_base[m_stride] += value[1u];
    m_base[2u * m_stride] += value[2u];
    return *this;
  }

  VectorReference &operator-=(Utils::Vector3d const &value) {
    m_base[0u] -= value[0u];
    m_base[m_stride] -= value[1u];
    m_base[2u * m_stride] -= value[2u];
    return *this;
  }

  VectorReference &operator*=(double const factor) {
    m_base[0u] *= factor;
    m_base[m_stride] *= factor;
    m_base[2u * m_stride] *= factor;
    return *this;
  }

  double &operator[](std::size_t const j) {
    assert(j < 3u);
    return m_base[j * m_stride];
  }
  double const &operator[](std::size_t const j) const {
    assert(j < 3u);
    return m_base[j * m_stride];
  }

  double norm2() const {
    auto const x = m_base[0u];
    auto const y = m_base[m_stride];
    auto const z = m_base[2u * m_stride];
    return x * x + y * y + z * z;
  }
  double norm() const { return std::sqrt(norm2()); }
};
