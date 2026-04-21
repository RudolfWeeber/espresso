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

#include <config/config.hpp>

#include "aosoa_pack.hpp"
#include "bond_error.hpp"
#include "cell_system/LocalBondState.hpp"
#include "pressure_cabana.hpp"
#include "pressure_inline.hpp"

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <optional>
#include <variant>

struct BondsPressureKernelData {
  BondedInteractionsMap const &bonded_ias;
  BoxGeometry const &box_geo;
  Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
};

struct PairBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::PairBondlistType bond_list;
  LocalBondState::PairBondIDType bond_ids;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_f_kernel;

  PairBondsPressureKernel(
      BondsPressureKernelData data_,
      LocalBondState::PairBondlistType bond_list_,
      LocalBondState::PairBondIDType bond_ids_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)), coulomb_f_kernel(coulomb_f_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_pressure = data.local_pressure;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    auto const tid = omp_get_thread_num();

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const dx = box_geo.get_mi_vector(pos1, pos2);

    auto const pair_force =
        calc_bond_pair_force(iaparams, dx,
#ifdef ESPRESSO_ELECTROSTATICS
                             aosoa.charge(i) * aosoa.charge(j), coulomb_f_kernel
#else
                             0.0, nullptr
#endif
        );

    if (pair_force) {
      auto const pressure = Utils::tensor_product(*pair_force, dx);
      auto const bin = layout.bonded_idx(bond_id);
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          local_pressure(tid, bin + k * 3 + l) += pressure(k, l);
    } else {
      auto partner_id = aosoa.id(j);
      bond_broken_error(aosoa.id(i), {&partner_id, 1});
    }
  }
};

struct AngleBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::AngleBondlistType bond_list;
  LocalBondState::AngleBondIDType bond_ids;

  AngleBondsPressureKernel(BondsPressureKernelData data_,
                           LocalBondState::AngleBondlistType bond_list_,
                           LocalBondState::AngleBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_pressure = data.local_pressure;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    auto const tid = omp_get_thread_num();

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const k = bond_list(idx, 2);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
    auto const dx21 = -box_geo.get_mi_vector(pos1, pos2);
    auto const dx31 = box_geo.get_mi_vector(pos3, pos1);

    auto const result = calc_bonded_three_body_force(iaparams, dx21, dx31);

    if (result) {
      Utils::Vector3d force2, force3;
      std::tie(std::ignore, force2, force3) = result.value();
      auto const pressure = Utils::tensor_product(force2, dx21) +
                            Utils::tensor_product(force3, dx31);
      auto const bin = layout.bonded_idx(bond_id);
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          local_pressure(tid, bin + m * 3 + n) += pressure(m, n);
    } else {
      std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
      bond_broken_error(aosoa.id(i), {pids.data(), 2});
    }
  }
};

struct DihedralBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::DihedralBondlistType bond_list;
  LocalBondState::DihedralBondIDType bond_ids;

  DihedralBondsPressureKernel(BondsPressureKernelData data_,
                              LocalBondState::DihedralBondlistType bond_list_,
                              LocalBondState::DihedralBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto &local_pressure = data.local_pressure;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    auto const tid = omp_get_thread_num();

    auto const bin = layout.bonded_idx(bond_id);
    for (int m = 0; m < 3; ++m)
      for (int n = 0; n < 3; ++n)
        local_pressure(tid, bin + m * 3 + n) += 0.0;
  }
};
