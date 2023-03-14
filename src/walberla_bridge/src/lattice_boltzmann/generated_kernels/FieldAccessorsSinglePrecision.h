// kernel generated with pystencils v1.1.1, lbmpy v1.1.1,
// lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

/*
 * Copyright (C) 2021-2023 The ESPResSo project
 * Copyright (C) 2020 The waLBerla project
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
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>

#include <field/GhostLayerField.h>
#include <stencil/D3Q19.h>

#include <array>
#include <tuple>

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla {
namespace lbm {
namespace accessor {

namespace Population {
inline std::array<float, 19u>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field, Cell const &cell) {
  float const &xyz0 = pdf_field->get(cell, uint_t{0u});
  std::array<float, 19u> pop;
  pop[0] = pdf_field->getF(&xyz0, 0);
  pop[1] = pdf_field->getF(&xyz0, 1);
  pop[2] = pdf_field->getF(&xyz0, 2);
  pop[3] = pdf_field->getF(&xyz0, 3);
  pop[4] = pdf_field->getF(&xyz0, 4);
  pop[5] = pdf_field->getF(&xyz0, 5);
  pop[6] = pdf_field->getF(&xyz0, 6);
  pop[7] = pdf_field->getF(&xyz0, 7);
  pop[8] = pdf_field->getF(&xyz0, 8);
  pop[9] = pdf_field->getF(&xyz0, 9);
  pop[10] = pdf_field->getF(&xyz0, 10);
  pop[11] = pdf_field->getF(&xyz0, 11);
  pop[12] = pdf_field->getF(&xyz0, 12);
  pop[13] = pdf_field->getF(&xyz0, 13);
  pop[14] = pdf_field->getF(&xyz0, 14);
  pop[15] = pdf_field->getF(&xyz0, 15);
  pop[16] = pdf_field->getF(&xyz0, 16);
  pop[17] = pdf_field->getF(&xyz0, 17);
  pop[18] = pdf_field->getF(&xyz0, 18);
  return pop;
}

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                std::array<float, 19u> const &pop, Cell const &cell) {
  float &xyz0 = pdf_field->get(cell, uint_t{0u});
  pdf_field->getF(&xyz0, 0) = pop[0];
  pdf_field->getF(&xyz0, 1) = pop[1];
  pdf_field->getF(&xyz0, 2) = pop[2];
  pdf_field->getF(&xyz0, 3) = pop[3];
  pdf_field->getF(&xyz0, 4) = pop[4];
  pdf_field->getF(&xyz0, 5) = pop[5];
  pdf_field->getF(&xyz0, 6) = pop[6];
  pdf_field->getF(&xyz0, 7) = pop[7];
  pdf_field->getF(&xyz0, 8) = pop[8];
  pdf_field->getF(&xyz0, 9) = pop[9];
  pdf_field->getF(&xyz0, 10) = pop[10];
  pdf_field->getF(&xyz0, 11) = pop[11];
  pdf_field->getF(&xyz0, 12) = pop[12];
  pdf_field->getF(&xyz0, 13) = pop[13];
  pdf_field->getF(&xyz0, 14) = pop[14];
  pdf_field->getF(&xyz0, 15) = pop[15];
  pdf_field->getF(&xyz0, 16) = pop[16];
  pdf_field->getF(&xyz0, 17) = pop[17];
  pdf_field->getF(&xyz0, 18) = pop[18];
}
} // namespace Population

namespace Vector {
inline Vector3<float> get(GhostLayerField<float, uint_t{3u}> const *vec_field,
                          Cell const &cell) {
  const float &xyz0 = vec_field->get(cell, 0);
  Vector3<float> vec;
  vec[0] = vec_field->getF(&xyz0, 0);
  vec[1] = vec_field->getF(&xyz0, 1);
  vec[2] = vec_field->getF(&xyz0, 2);
  return vec;
}

inline void set(GhostLayerField<float, uint_t{3u}> *vec_field,
                Vector3<float> const &vec, Cell const &cell) {
  float &xyz0 = vec_field->get(cell, 0);
  vec_field->getF(&xyz0, 0) = vec[0];
  vec_field->getF(&xyz0, 1) = vec[1];
  vec_field->getF(&xyz0, 2) = vec[2];
}

inline void add(GhostLayerField<float, uint_t{3u}> *vec_field,
                Vector3<float> const &vec, Cell const &cell) {
  float &xyz0 = vec_field->get(cell, 0);
  vec_field->getF(&xyz0, 0) += vec[0];
  vec_field->getF(&xyz0, 1) += vec[1];
  vec_field->getF(&xyz0, 2) += vec[2];
}

inline void broadcast(GhostLayerField<float, uint_t{3u}> *vec_field,
                      Vector3<float> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    float &xyz0 = vec_field->get(x, y, z, 0);
    vec_field->getF(&xyz0, 0) = vec[0];
    vec_field->getF(&xyz0, 1) = vec[1];
    vec_field->getF(&xyz0, 2) = vec[2];
  });
}

inline void add_to_all(GhostLayerField<float, uint_t{3u}> *vec_field,
                       Vector3<float> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    float &xyz0 = vec_field->get(x, y, z, 0);
    vec_field->getF(&xyz0, 0) += vec[0];
    vec_field->getF(&xyz0, 1) += vec[1];
    vec_field->getF(&xyz0, 2) += vec[2];
  });
}
} // namespace Vector

namespace EquilibriumDistribution {
inline float get(stencil::Direction const direction,
                 Vector3<float> const &u = Vector3<float>(float(0.0)),
                 float rho = float(1.0)) {

  using namespace stencil;
  switch (direction) {
  case C:
    return rho * -0.33333333333333331f * (u[0] * u[0]) +
           rho * -0.33333333333333331f * (u[1] * u[1]) +
           rho * -0.33333333333333331f * (u[2] * u[2]) +
           rho * 0.33333333333333331f;
  case N:
    return rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[1] +
           rho * 0.16666666666666666f * (u[1] * u[1]);
  case S:
    return rho * -0.16666666666666666f * u[1] +
           rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[1] * u[1]);
  case W:
    return rho * -0.16666666666666666f * u[0] +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[0] * u[0]);
  case E:
    return rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[0] +
           rho * 0.16666666666666666f * (u[0] * u[0]);
  case T:
    return rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[2] +
           rho * 0.16666666666666666f * (u[2] * u[2]);
  case B:
    return rho * -0.16666666666666666f * u[2] +
           rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[2] * u[2]);
  case NW:
    return rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[1] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]);
  case NE:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.25f * u[0] * u[1];
  case SW:
    return rho * -0.083333333333333329f * u[0] +
           rho * -0.083333333333333329f * u[1] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.25f * u[0] * u[1];
  case SE:
    return rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[0] * u[1] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]);
  case TN:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[1] * u[2];
  case TS:
    return rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[1] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case TW:
    return rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case TE:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[0] * u[2];
  case BN:
    return rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[1] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case BS:
    return rho * -0.083333333333333329f * u[1] +
           rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[1] * u[2];
  case BW:
    return rho * -0.083333333333333329f * u[0] +
           rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[0] * u[2];
  case BE:
    return rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[0] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  default:
    WALBERLA_ABORT("Invalid Direction")
  }
}
} // namespace EquilibriumDistribution

namespace Equilibrium {
inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                Vector3<float> const &u, float const rho, Cell const &cell) {

  float &xyz0 = pdf_field->get(cell, 0);
  pdf_field->getF(&xyz0, 0) = rho * -0.33333333333333331f * (u[0] * u[0]) +
                              rho * -0.33333333333333331f * (u[1] * u[1]) +
                              rho * -0.33333333333333331f * (u[2] * u[2]) +
                              rho * 0.33333333333333331f;
  pdf_field->getF(&xyz0, 1) = rho * -0.16666666666666666f * (u[0] * u[0]) +
                              rho * -0.16666666666666666f * (u[2] * u[2]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * u[1] +
                              rho * 0.16666666666666666f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 2) = rho * -0.16666666666666666f * u[1] +
                              rho * -0.16666666666666666f * (u[0] * u[0]) +
                              rho * -0.16666666666666666f * (u[2] * u[2]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 3) = rho * -0.16666666666666666f * u[0] +
                              rho * -0.16666666666666666f * (u[1] * u[1]) +
                              rho * -0.16666666666666666f * (u[2] * u[2]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * (u[0] * u[0]);
  pdf_field->getF(&xyz0, 4) = rho * -0.16666666666666666f * (u[1] * u[1]) +
                              rho * -0.16666666666666666f * (u[2] * u[2]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * u[0] +
                              rho * 0.16666666666666666f * (u[0] * u[0]);
  pdf_field->getF(&xyz0, 5) = rho * -0.16666666666666666f * (u[0] * u[0]) +
                              rho * -0.16666666666666666f * (u[1] * u[1]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * u[2] +
                              rho * 0.16666666666666666f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 6) = rho * -0.16666666666666666f * u[2] +
                              rho * -0.16666666666666666f * (u[0] * u[0]) +
                              rho * -0.16666666666666666f * (u[1] * u[1]) +
                              rho * 0.055555555555555552f +
                              rho * 0.16666666666666666f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 7) =
      rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[1] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 8) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]) + rho * 0.25f * u[0] * u[1];
  pdf_field->getF(&xyz0, 9) =
      rho * -0.083333333333333329f * u[0] +
      rho * -0.083333333333333329f * u[1] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]) + rho * 0.25f * u[0] * u[1];
  pdf_field->getF(&xyz0, 10) =
      rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[0] * u[1] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 11) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[1] * u[2];
  pdf_field->getF(&xyz0, 12) =
      rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[1] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 13) =
      rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 14) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[0] * u[2];
  pdf_field->getF(&xyz0, 15) =
      rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[1] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 16) =
      rho * -0.083333333333333329f * u[1] +
      rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[1] * u[2];
  pdf_field->getF(&xyz0, 17) =
      rho * -0.083333333333333329f * u[0] +
      rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[0] * u[2];
  pdf_field->getF(&xyz0, 18) =
      rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[0] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
inline float get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
                 Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, 0);
  const float f_0 = pdf_field->getF(&xyz0, 0);
  const float f_1 = pdf_field->getF(&xyz0, 1);
  const float f_2 = pdf_field->getF(&xyz0, 2);
  const float f_3 = pdf_field->getF(&xyz0, 3);
  const float f_4 = pdf_field->getF(&xyz0, 4);
  const float f_5 = pdf_field->getF(&xyz0, 5);
  const float f_6 = pdf_field->getF(&xyz0, 6);
  const float f_7 = pdf_field->getF(&xyz0, 7);
  const float f_8 = pdf_field->getF(&xyz0, 8);
  const float f_9 = pdf_field->getF(&xyz0, 9);
  const float f_10 = pdf_field->getF(&xyz0, 10);
  const float f_11 = pdf_field->getF(&xyz0, 11);
  const float f_12 = pdf_field->getF(&xyz0, 12);
  const float f_13 = pdf_field->getF(&xyz0, 13);
  const float f_14 = pdf_field->getF(&xyz0, 14);
  const float f_15 = pdf_field->getF(&xyz0, 15);
  const float f_16 = pdf_field->getF(&xyz0, 16);
  const float f_17 = pdf_field->getF(&xyz0, 17);
  const float f_18 = pdf_field->getF(&xyz0, 18);
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float vel2Term = f_12 + f_13 + f_5;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;
  return rho;
}
} // namespace Density

namespace DensityAndVelocity {
inline std::tuple<float, Vector3<float>>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
    GhostLayerField<float, uint_t{3u}> const *force_field, Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const float &xyz0 = pdf_field->get(cell, 0);
  const float f_0 = pdf_field->getF(&xyz0, 0);
  const float f_1 = pdf_field->getF(&xyz0, 1);
  const float f_2 = pdf_field->getF(&xyz0, 2);
  const float f_3 = pdf_field->getF(&xyz0, 3);
  const float f_4 = pdf_field->getF(&xyz0, 4);
  const float f_5 = pdf_field->getF(&xyz0, 5);
  const float f_6 = pdf_field->getF(&xyz0, 6);
  const float f_7 = pdf_field->getF(&xyz0, 7);
  const float f_8 = pdf_field->getF(&xyz0, 8);
  const float f_9 = pdf_field->getF(&xyz0, 9);
  const float f_10 = pdf_field->getF(&xyz0, 10);
  const float f_11 = pdf_field->getF(&xyz0, 11);
  const float f_12 = pdf_field->getF(&xyz0, 12);
  const float f_13 = pdf_field->getF(&xyz0, 13);
  const float f_14 = pdf_field->getF(&xyz0, 14);
  const float f_15 = pdf_field->getF(&xyz0, 15);
  const float f_16 = pdf_field->getF(&xyz0, 16);
  const float f_17 = pdf_field->getF(&xyz0, 17);
  const float f_18 = pdf_field->getF(&xyz0, 18);
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const float vel2Term = f_12 + f_13 + f_5;
  const float momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;
  const float md_0 =
      force_field->get(x, y, z, 0) * 0.50000000000000000f + momdensity_0;
  const float md_1 =
      force_field->get(x, y, z, 1) * 0.50000000000000000f + momdensity_1;
  const float md_2 =
      force_field->get(x, y, z, 2) * 0.50000000000000000f + momdensity_2;

  const float conversion = float(1) / rho;
  Vector3<float> velocity;
  velocity[0] = md_0 * conversion;
  velocity[1] = md_1 * conversion;
  velocity[2] = md_2 * conversion;

  return {rho, velocity};
}

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                GhostLayerField<float, uint_t{3u}> const *force_field,
                Vector3<float> const &u, float const rho_in, Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const float rho = rho_in;
  const float delta_rho = rho - 1;
  const float u_0 =
      -force_field->get(x, y, z, 0) * 0.50000000000000000f / rho + u[0];
  const float u_1 =
      -force_field->get(x, y, z, 1) * 0.50000000000000000f / rho + u[1];
  const float u_2 =
      -force_field->get(x, y, z, 2) * 0.50000000000000000f / rho + u[2];

  Equilibrium::set(pdf_field, Vector3<float>(u_0, u_1, u_2), rho, cell);
}
} // namespace DensityAndVelocity

namespace DensityAndMomentumDensity {
inline std::tuple<float, Vector3<float>>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
    GhostLayerField<float, uint_t{3u}> const *force_field, Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const float &xyz0 = pdf_field->get(cell, 0);
  const float f_0 = pdf_field->getF(&xyz0, 0);
  const float f_1 = pdf_field->getF(&xyz0, 1);
  const float f_2 = pdf_field->getF(&xyz0, 2);
  const float f_3 = pdf_field->getF(&xyz0, 3);
  const float f_4 = pdf_field->getF(&xyz0, 4);
  const float f_5 = pdf_field->getF(&xyz0, 5);
  const float f_6 = pdf_field->getF(&xyz0, 6);
  const float f_7 = pdf_field->getF(&xyz0, 7);
  const float f_8 = pdf_field->getF(&xyz0, 8);
  const float f_9 = pdf_field->getF(&xyz0, 9);
  const float f_10 = pdf_field->getF(&xyz0, 10);
  const float f_11 = pdf_field->getF(&xyz0, 11);
  const float f_12 = pdf_field->getF(&xyz0, 12);
  const float f_13 = pdf_field->getF(&xyz0, 13);
  const float f_14 = pdf_field->getF(&xyz0, 14);
  const float f_15 = pdf_field->getF(&xyz0, 15);
  const float f_16 = pdf_field->getF(&xyz0, 16);
  const float f_17 = pdf_field->getF(&xyz0, 17);
  const float f_18 = pdf_field->getF(&xyz0, 18);
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const float vel2Term = f_12 + f_13 + f_5;
  const float momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;
  const float md_0 =
      force_field->get(x, y, z, 0) * 0.50000000000000000f + momdensity_0;
  const float md_1 =
      force_field->get(x, y, z, 1) * 0.50000000000000000f + momdensity_1;
  const float md_2 =
      force_field->get(x, y, z, 2) * 0.50000000000000000f + momdensity_2;

  Vector3<float> momentumDensity;
  momentumDensity[0] = md_0;
  momentumDensity[1] = md_1;
  momentumDensity[2] = md_2;

  return {rho, momentumDensity};
}
} // namespace DensityAndMomentumDensity

namespace MomentumDensity {
inline Vector3<float>
reduce(GhostLayerField<float, uint_t{19u}> const *pdf_field,
       GhostLayerField<float, uint_t{3u}> const *force_field) {
  Vector3<float> momentumDensity(float{0});
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    const float &xyz0 = pdf_field->get(x, y, z, 0);
    const float f_0 = pdf_field->getF(&xyz0, 0);
    const float f_1 = pdf_field->getF(&xyz0, 1);
    const float f_2 = pdf_field->getF(&xyz0, 2);
    const float f_3 = pdf_field->getF(&xyz0, 3);
    const float f_4 = pdf_field->getF(&xyz0, 4);
    const float f_5 = pdf_field->getF(&xyz0, 5);
    const float f_6 = pdf_field->getF(&xyz0, 6);
    const float f_7 = pdf_field->getF(&xyz0, 7);
    const float f_8 = pdf_field->getF(&xyz0, 8);
    const float f_9 = pdf_field->getF(&xyz0, 9);
    const float f_10 = pdf_field->getF(&xyz0, 10);
    const float f_11 = pdf_field->getF(&xyz0, 11);
    const float f_12 = pdf_field->getF(&xyz0, 12);
    const float f_13 = pdf_field->getF(&xyz0, 13);
    const float f_14 = pdf_field->getF(&xyz0, 14);
    const float f_15 = pdf_field->getF(&xyz0, 15);
    const float f_16 = pdf_field->getF(&xyz0, 16);
    const float f_17 = pdf_field->getF(&xyz0, 17);
    const float f_18 = pdf_field->getF(&xyz0, 18);
    const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const float vel1Term = f_1 + f_11 + f_15 + f_7;
    const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const float vel2Term = f_12 + f_13 + f_5;
    const float momdensity_2 =
        f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                      vel1Term + vel2Term;
    const float md_0 =
        force_field->get(x, y, z, 0) * 0.50000000000000000f + momdensity_0;
    const float md_1 =
        force_field->get(x, y, z, 1) * 0.50000000000000000f + momdensity_1;
    const float md_2 =
        force_field->get(x, y, z, 2) * 0.50000000000000000f + momdensity_2;

    momentumDensity[0] += md_0;
    momentumDensity[1] += md_1;
    momentumDensity[2] += md_2;
  });
  return momentumDensity;
}
} // namespace MomentumDensity

namespace PressureTensor {
inline Matrix3<float> get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
                          Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, 0);
  const float f_0 = pdf_field->getF(&xyz0, 0);
  const float f_1 = pdf_field->getF(&xyz0, 1);
  const float f_2 = pdf_field->getF(&xyz0, 2);
  const float f_3 = pdf_field->getF(&xyz0, 3);
  const float f_4 = pdf_field->getF(&xyz0, 4);
  const float f_5 = pdf_field->getF(&xyz0, 5);
  const float f_6 = pdf_field->getF(&xyz0, 6);
  const float f_7 = pdf_field->getF(&xyz0, 7);
  const float f_8 = pdf_field->getF(&xyz0, 8);
  const float f_9 = pdf_field->getF(&xyz0, 9);
  const float f_10 = pdf_field->getF(&xyz0, 10);
  const float f_11 = pdf_field->getF(&xyz0, 11);
  const float f_12 = pdf_field->getF(&xyz0, 12);
  const float f_13 = pdf_field->getF(&xyz0, 13);
  const float f_14 = pdf_field->getF(&xyz0, 14);
  const float f_15 = pdf_field->getF(&xyz0, 15);
  const float f_16 = pdf_field->getF(&xyz0, 16);
  const float f_17 = pdf_field->getF(&xyz0, 17);
  const float f_18 = pdf_field->getF(&xyz0, 18);
  const float p_0 =
      f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
  const float p_1 = -f_10 - f_7 + f_8 + f_9;
  const float p_2 = -f_13 + f_14 + f_17 - f_18;
  const float p_3 = -f_10 - f_7 + f_8 + f_9;
  const float p_4 =
      f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
  const float p_5 = f_11 - f_12 - f_15 + f_16;
  const float p_6 = -f_13 + f_14 + f_17 - f_18;
  const float p_7 = f_11 - f_12 - f_15 + f_16;
  const float p_8 =
      f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;

  Matrix3<float> pressureTensor;
  pressureTensor[0] = p_0;
  pressureTensor[1] = p_1;
  pressureTensor[2] = p_2;

  pressureTensor[3] = p_3;
  pressureTensor[4] = p_4;
  pressureTensor[5] = p_5;

  pressureTensor[6] = p_6;
  pressureTensor[7] = p_7;
  pressureTensor[8] = p_8;

  return pressureTensor;
}
} // namespace PressureTensor

} // namespace accessor
} // namespace lbm
} // namespace walberla

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif