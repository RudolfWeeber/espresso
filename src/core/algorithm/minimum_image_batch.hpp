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

/** @file
 *  Batched (vectorizable) minimum-image primitives for the ORTHORHOMBIC,
 *  periodic, non-Lees-Edwards fast path.
 *
 *  These reproduce the EXACT per-coordinate arithmetic of
 *  @ref BoxGeometry::get_mi_vector / @ref BoxGeometry::get_mi_dist2 on the
 *  cuboid branch (see @c detail::get_mi_coord in BoxGeometry.hpp), so that
 *  for each reference/neighbor pair the produced minimum-image vector (resp.
 *  squared distance) is BITWISE-identical to the scalar call. Callers must fall
 *  back to the scalar path when Lees-Edwards is active
 *  (@c box_geo.type() == BoxType::LEES_EDWARDS).
 *
 *  Autovectorization (gcc 13, `-march=native`, AVX2+FMA): the per-coordinate
 *  fold is expressed as pure per-lane arithmetic with a single comparison-
 *  derived mask and no libc call, so `#pragma omp simd` (the core is built with
 *  `-fopenmp`) turns each hot loop into `vmulpd`/`vcmpltpd`/`vroundpd`/
 *  `vfnmadd*pd` over 256-bit ymm lanes.
 *
 *  Bitwise-identity subtleties:
 *   1. The scalar fold is *conditional*:
 *        `(periodic && abs(dx) > box_half) ? dx - round(dx*inv)*box : dx`.
 *      We keep the identical `>` select via a 1.0/0.0 mask; at the exact
 *      boundary `abs(dx) == box_half` the mask is 0 (matching the scalar
 *      `return dx`). Per-axis periodicity is folded into an effective half
 *      length: a non-periodic axis uses `+inf`, so `abs(dx) > +inf` is always
 *      false and the axis is never folded (matching `periodic && ...`).
 *   2. The scalar uses `std::round` (round-half-away-from-zero), which does not
 *      vectorize. `std::rint` maps to `vroundpd` but is round-to-nearest-EVEN;
 *      @c round_half_away corrects the tie cases branchlessly so the result
 *      is bitwise-identical to `std::round` for every argument the fold rounds.
 *      See @c minimum_image_batch_test for the exhaustive tie / near-tie sweep.
 *
 *  Layout: the reference position is a scalar triple; neighbor positions are
 *  passed as three restrict-qualified contiguous component arrays (SoA), which
 *  is what the pair kernel / Verlet build gather into their scratch. Plain
 *  doubles throughout; no lane intrinsics.
 */

#include "BoxGeometry.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace Algorithm {

/**
 * @brief Precomputed orthorhombic box parameters for the batch fold.
 *
 * Mirrors the members @c detail::get_mi_coord reads for each coordinate:
 * length and inverse length. Periodicity is baked into @c half_or_inf: the
 * per-axis half length for a periodic axis, or `+inf` for a non-periodic axis
 * (so the per-lane comparison `abs(dx) > half_or_inf` is always false there and
 * the axis is never folded, exactly matching the scalar `periodic && ...`).
 * Keeping periodicity out of the inner comparison is what lets gcc vectorize
 * the fold: a separate loop-invariant `periodic_mask` multiply on the masked
 * term defeats gcc 13's if-conversion.
 */
struct OrthoBoxParams {
  double length[3];
  double length_inv[3];
  double half_or_inf[3]; //< box_half if periodic, else +inf

  /** Build from a cuboid @ref BoxGeometry (LE-inactive fast path only). */
  static OrthoBoxParams from(BoxGeometry const &box_geo) {
    OrthoBoxParams p;
    auto const &l = box_geo.length();
    auto const &li = box_geo.length_inv();
    auto const &lh = box_geo.length_half();
    for (unsigned c = 0u; c < 3u; ++c) {
      p.length[c] = l[c];
      p.length_inv[c] = li[c];
      p.half_or_inf[c] =
          box_geo.periodic(c) ? lh[c] : std::numeric_limits<double>::infinity();
    }
    return p;
  }
};

/**
 * @brief Vectorizable round-half-away-from-zero, bitwise-identical to
 *        @c std::round for every argument the minimum-image fold rounds.
 *
 * `std::round` (round-half-away) is a libc call with control flow that blocks
 * autovectorization. `std::rint` lowers to `vroundpd` but is round-to-nearest-
 * EVEN, which differs from `std::round` only at exact half-integer ties
 * `y == n + 0.5`. We correct the ties that `rint` rounded TOWARD zero: at a tie
 * `diff = y - rint(y)` is exactly `+/-0.5`, and `diff` has the same sign as `y`
 * iff `rint` rounded toward zero (round-half-away wants the away neighbor
 * there). Both the tie test and the direction test are equalities/inequalities
 * that fold into vector masks, so the whole helper vectorizes.
 */
[[nodiscard]] inline double round_half_away(double y) {
  double const r = std::rint(y); // vroundpd, round to nearest even
  double const diff = y - r;     // exact, in [-0.5, 0.5]
  auto const is_tie = static_cast<double>(std::fabs(diff) == 0.5);
  auto const toward_zero = static_cast<double>(diff * y > 0.0);
  return r + is_tie * toward_zero * std::copysign(1.0, y);
}

/**
 * @brief Fold one coordinate difference, bitwise-identical to the scalar path.
 *
 * Reproduces @c detail::get_mi_coord exactly on the cuboid branch:
 *   `(periodic && abs(dx) > box_half) ? dx - round(dx*inv)*box : dx`.
 * @p half_or_inf is `box_half` (periodic axis) or `+inf` (non-periodic axis),
 * so the mask below reproduces the full `periodic && abs(dx) > box_half`
 * condition. The result is `dx - do_fold * round(dx*inv)*box`, equal to `dx`
 * when `do_fold == 0` (scalar `return dx`, including the exact boundary) and to
 * `dx - round(dx*inv)*box` when `do_fold == 1` (scalar folded arm). On folded
 * lanes `|dx*box_inv| > 0.5`, where @ref round_half_away matches `std::round`.
 */
[[nodiscard]] inline double fold_coord_branchless(double dx, double box,
                                                  double box_inv,
                                                  double half_or_inf) {
  auto const do_fold = static_cast<double>(std::abs(dx) > half_or_inf);
  return dx - do_fold * (round_half_away(dx * box_inv) * box);
}

/**
 * @brief One reference position vs N neighbor positions -> N MI vectors.
 *
 * @param[in]  box   Precomputed orthorhombic box parameters (cuboid branch).
 * @param[in]  ref   Reference position (terminal point @c a); MI vector is
 *                   `ref - neighbor` per axis, matching @c get_mi_vector(a, b).
 * @param[in]  nx    Neighbor x components (contiguous, @p n entries).
 * @param[in]  ny    Neighbor y components.
 * @param[in]  nz    Neighbor z components.
 * @param[out] out_x MI-vector x components (contiguous, @p n entries).
 * @param[out] out_y MI-vector y components.
 * @param[out] out_z MI-vector z components.
 * @param[in]  n     Number of neighbors.
 *
 * Each output element equals `box_geo.get_mi_vector(ref, {nx,ny,nz})` for the
 * cuboid branch, bitwise. Hot loop autovectorizes to ymm under `-march=native`.
 */
inline void get_mi_vector_batch(OrthoBoxParams const &box, double const ref[3],
                                double const *__restrict nx,
                                double const *__restrict ny,
                                double const *__restrict nz,
                                double *__restrict out_x,
                                double *__restrict out_y,
                                double *__restrict out_z, std::size_t n) {
  double const rx = ref[0], ry = ref[1], rz = ref[2];
  double const lx = box.length[0], ly = box.length[1], lz = box.length[2];
  double const lix = box.length_inv[0], liy = box.length_inv[1],
               liz = box.length_inv[2];
  double const hx = box.half_or_inf[0], hy = box.half_or_inf[1],
               hz = box.half_or_inf[2];
#pragma omp simd
  for (std::size_t k = 0u; k < n; ++k) {
    out_x[k] = fold_coord_branchless(rx - nx[k], lx, lix, hx);
    out_y[k] = fold_coord_branchless(ry - ny[k], ly, liy, hy);
    out_z[k] = fold_coord_branchless(rz - nz[k], lz, liz, hz);
  }
}

/**
 * @brief One reference position vs N neighbor positions -> N squared MI dists.
 *
 * @param[in]  box   Precomputed orthorhombic box parameters (cuboid branch).
 * @param[in]  ref   Reference position (terminal point @c a).
 * @param[in]  nx    Neighbor x components (contiguous, @p n entries).
 * @param[in]  ny    Neighbor y components.
 * @param[in]  nz    Neighbor z components.
 * @param[out] out   Squared MI distances (contiguous, @p n entries).
 * @param[in]  n     Number of neighbors.
 *
 * Each output element equals `box_geo.get_mi_dist2(ref, {nx,ny,nz})` for the
 * cuboid branch, bitwise (same `d0*d0 + d1*d1 + d2*d2` term order). The
 * verlet-build gate consumes @p out as the precomputed
 * `box_geo.get_mi_dist2(...)` value.
 */
inline void get_mi_dist2_batch(OrthoBoxParams const &box, double const ref[3],
                               double const *__restrict nx,
                               double const *__restrict ny,
                               double const *__restrict nz,
                               double *__restrict out, std::size_t n) {
  double const rx = ref[0], ry = ref[1], rz = ref[2];
  double const lx = box.length[0], ly = box.length[1], lz = box.length[2];
  double const lix = box.length_inv[0], liy = box.length_inv[1],
               liz = box.length_inv[2];
  double const hx = box.half_or_inf[0], hy = box.half_or_inf[1],
               hz = box.half_or_inf[2];
#pragma omp simd
  for (std::size_t k = 0u; k < n; ++k) {
    double const d0 = fold_coord_branchless(rx - nx[k], lx, lix, hx);
    double const d1 = fold_coord_branchless(ry - ny[k], ly, liy, hy);
    double const d2 = fold_coord_branchless(rz - nz[k], lz, liz, hz);
    out[k] = d0 * d0 + d1 * d1 + d2 * d2;
  }
}

} // namespace Algorithm
