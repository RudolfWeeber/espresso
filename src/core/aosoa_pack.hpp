/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#ifdef SHARED_MEMORY_PARALLELISM

#include <Cabana_Core.hpp>

const int vector_length = 1;
using data_types = Cabana::MemberTypes<double[3], double, int, int>; //, bool>;
using memory_space = Kokkos::HostSpace; // Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;
using AoSoA_type = Cabana::AoSoA<data_types, memory_space, vector_length>;

struct AoSoA_pack {
  AoSoA_type::member_slice_type<0> position;
  // AoSoA_type::member_slice_type<1> force;
  // AoSoA_type::member_slice_type<2> torque;
  AoSoA_type::member_slice_type<1> charge;
  AoSoA_type::member_slice_type<2> id;
  AoSoA_type::member_slice_type<3> type;
  // AoSoA_type::member_slice_type<4> ghost;

  AoSoA_pack() = default;

  AoSoA_pack(AoSoA_type &aosoa)
      : // position(Cabana::slice<0>(aosoa)), force(Cabana::slice<1>(aosoa)),
        // torque(Cabana::slice<2>(aosoa)), charge(Cabana::slice<3>(aosoa)),
        position(Cabana::slice<0>(aosoa)), charge(Cabana::slice<1>(aosoa)),
        id(Cabana::slice<2>(aosoa)), type(Cabana::slice<3>(aosoa)) {}
  // ghost(Cabana::slice<4>(aosoa)) {}
};
#endif
