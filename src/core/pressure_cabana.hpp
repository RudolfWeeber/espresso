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
#include "pressure_inline.hpp"
#include "short_range_cabana_helpers.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

struct PressureBinLayout {
  std::size_t n_bonded;
  std::size_t n_types;
  std::size_t off_bonded = 0;
  std::size_t off_nb_inter;
  std::size_t off_nb_intra;
  std::size_t off_coulomb;
  std::size_t off_dipolar;
  std::size_t off_dpd;
  std::size_t off_kinetic;
  std::size_t total;

  PressureBinLayout(std::size_t n_bonded_, std::size_t n_types_)
      : n_bonded(n_bonded_), n_types(n_types_) {
    auto const n_nb = n_types * (n_types + 1) / 2;
    off_nb_inter = off_bonded + 9 * n_bonded;
    off_nb_intra = off_nb_inter + 9 * n_nb;
    off_coulomb = off_nb_intra + 9 * n_nb;
    off_dipolar = off_coulomb + 9;
    off_dpd = off_dipolar + 9;
    off_kinetic = off_dpd + 9;
    total = off_kinetic + 9;
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t nb_inter_idx(int t1, int t2) const {
    return off_nb_inter +
           9 * Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t nb_intra_idx(int t1, int t2) const {
    return off_nb_intra +
           9 * Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  KOKKOS_INLINE_FUNCTION std::size_t coulomb_idx() const { return off_coulomb; }
  KOKKOS_INLINE_FUNCTION std::size_t dipolar_idx() const { return off_dipolar; }
  KOKKOS_INLINE_FUNCTION std::size_t dpd_idx() const { return off_dpd; }
  KOKKOS_INLINE_FUNCTION std::size_t kinetic_idx() const { return off_kinetic; }
  KOKKOS_INLINE_FUNCTION std::size_t bonded_idx(int b) const {
    return off_bonded + 9 * b;
  }
};

struct KineticPressureKernel {
  CellStructure::AoSoA_pack const &aosoa;
  Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  std::size_t n_local_particles;

  KineticPressureKernel(
      CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
      PressureBinLayout layout_, std::size_t n_local_particles_)
      : aosoa(aosoa_), local_pressure(local_pressure_),
        layout(std::move(layout_)), n_local_particles(n_local_particles_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(std::size_t idx) const {
    if (idx >= n_local_particles)
      return;
    if (aosoa.is_virtual(idx))
      return;

    auto const tid = omp_get_thread_num();
    auto const bin = layout.kinetic_idx();

    auto const mass = aosoa.mass(idx);
    auto const vx = aosoa.velocity(idx, 0);
    auto const vy = aosoa.velocity(idx, 1);
    auto const vz = aosoa.velocity(idx, 2);

    local_pressure(tid, bin + 0) += mass * vx * vx;
    local_pressure(tid, bin + 1) += mass * vx * vy;
    local_pressure(tid, bin + 2) += mass * vx * vz;
    local_pressure(tid, bin + 3) += mass * vy * vx;
    local_pressure(tid, bin + 4) += mass * vy * vy;
    local_pressure(tid, bin + 5) += mass * vy * vz;
    local_pressure(tid, bin + 6) += mass * vz * vx;
    local_pressure(tid, bin + 7) += mass * vz * vy;
    local_pressure(tid, bin + 8) += mass * vz * vz;
  }
};

struct PressureKernel {
  BondedInteractionsMap const &bonded_ias;
  InteractionsNonBonded const &nonbonded_ias;
  Coulomb::Solver const &coulomb;
  Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel;
  Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel;
  Dipoles::Solver const &dipoles;
  Thermostat::Thermostat const &thermostat;
  BoxGeometry const &box_geo;
  std::vector<Particle *> const &unique_particles;
  Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
  Kokkos::View<int *> mol_id_view;
  double system_max_cutoff;

  PressureKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::Solver const &coulomb_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel_,
      Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel_,
      Dipoles::Solver const &dipoles_,
      Thermostat::Thermostat const &thermostat_, BoxGeometry const &box_geo_,
      std::vector<Particle *> const &unique_particles_,
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
      PressureBinLayout layout_, CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<int *> mol_id_view_, double system_max_cutoff_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb(coulomb_), coulomb_f_kernel(coulomb_f_kernel_),
        coulomb_p_kernel(coulomb_p_kernel_), dipoles(dipoles_),
        thermostat(thermostat_), box_geo(box_geo_),
        unique_particles(unique_particles_), local_pressure(local_pressure_),
        layout(std::move(layout_)), aosoa(aosoa_),
        mol_id_view(std::move(mol_id_view_)),
        system_max_cutoff(system_max_cutoff_) {}

  KOKKOS_INLINE_FUNCTION
  void write_tensor(std::size_t tid, std::size_t bin,
                    Utils::Matrix<double, 3, 3> const &tensor) const {
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        local_pressure(tid, bin + k * 3 + l) += tensor(k, l);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(std::size_t i, std::size_t j) const {
    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const d = box_geo.get_mi_vector(pos1, pos2);
    auto const dist = d.norm();
    if (dist > system_max_cutoff)
      return;

    auto const t1 = aosoa.type(i);
    auto const t2 = aosoa.type(j);
    auto const &ia_params = nonbonded_ias.get_ia_param(t1, t2);

    auto const tid = omp_get_thread_num();

#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_THOLE)
    auto const flag = compute_pair_data_flags(
        dist, ia_params, coulomb_f_kernel != nullptr, aosoa, i, j);

    Particle const *p1_ptr = nullptr;
    Particle const *p2_ptr = nullptr;
    if (flag.need_particle_pointers) {
      p1_ptr = unique_particles.at(i);
      p2_ptr = unique_particles.at(j);
    }
#endif

    Utils::Vector3d f_nb{};

    if (dist <= ia_params.max_cut) {
#ifdef ESPRESSO_EXCLUSIONS
      bool skip = false;
      if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j))
        skip = not do_nonbonded(*p1_ptr, *p2_ptr);
      if (not skip)
#endif
      {
        f_nb += calc_central_radial_force(ia_params, d, dist);

#ifdef ESPRESSO_THOLE
        if (thole_active(ia_params, coulomb_f_kernel != nullptr)) {
          f_nb += thole_pair_force(*p1_ptr, *p2_ptr, ia_params, d, dist,
                                   bonded_ias, coulomb_f_kernel);
        }
#endif
#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
          auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
          auto const pf = gb_pair_force(dir1, dir2, ia_params, d, dist);
          f_nb += pf.f;
        }
#endif
      }

      auto const stress_nb = Utils::tensor_product(d, f_nb);
      auto const bin = (mol_id_view(i) == mol_id_view(j))
                           ? layout.nb_intra_idx(t1, t2)
                           : layout.nb_inter_idx(t1, t2);
      write_tensor(tid, bin, stress_nb);
    }

#ifdef ESPRESSO_ELECTROSTATICS
    if (coulomb_p_kernel != nullptr) {
      auto const q1 = aosoa.charge(i), q2 = aosoa.charge(j);
      if (q1 != 0. and q2 != 0.) {
        auto const p_coulomb = (*coulomb_p_kernel)(q1 * q2, d, dist);
        write_tensor(tid, layout.coulomb_idx(), p_coulomb);
      }
    }
#endif

#ifdef ESPRESSO_DIPOLES
    if (dipoles.impl->solver) {
      fprintf(stderr, "calculating pressure for magnetostatics which doesn't "
                      "have it implemented\n");
    }
#endif

#ifdef ESPRESSO_DPD
    if (dpd_active(ia_params, thermostat.thermo_switch)) {
      auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
      auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
      auto const v21 = box_geo.velocity_difference(pos1, pos2, vel1, vel2);

      Utils::Vector3d f_r =
          dpd_pair_force(ia_params.dpd.radial, v21, dist, Utils::Vector3d{});
      Utils::Vector3d f_t =
          dpd_pair_force(ia_params.dpd.trans, v21, dist, Utils::Vector3d{});

      auto const dist2 = dist * dist;
      auto const P = Utils::tensor_product(d / dist2, d);
      auto const f_dpd = P * (f_r - f_t) + f_t;

      auto const stress_dpd = -Utils::tensor_product(d, f_dpd);
      write_tensor(tid, layout.dpd_idx(), stress_dpd);
    }
#endif
  }
};

static void reduce_cabana_pressure(
    Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure,
    PressureBinLayout const &layout, Observable_stat &obs,
    BondedInteractionsMap const &bonded_ias, int n_types) {
  auto const nthreads = local_pressure.extent(0);
  auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, local_pressure);

  auto sum_bin = [&](std::size_t bin) {
    double s = 0.;
    for (std::size_t t = 0; t < nthreads; ++t)
      s += host(t, bin);
    return s;
  };

  for (int b = 0; b < int(layout.n_bonded); ++b) {
    auto const base = layout.bonded_idx(b);
    auto output = obs.bonded_contribution(b);
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        output[k * 3 + l] += sum_bin(base + k * 3 + l);
  }

  for (int t1 = 0; t1 < n_types; ++t1)
    for (int t2 = 0; t2 <= t1; ++t2) {
      auto const base_inter = layout.nb_inter_idx(t1, t2);
      auto output_inter = obs.non_bonded_inter_contribution(t1, t2);
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          output_inter[k * 3 + l] += sum_bin(base_inter + k * 3 + l);

      auto const base_intra = layout.nb_intra_idx(t1, t2);
      auto output_intra = obs.non_bonded_intra_contribution(t1, t2);
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          output_intra[k * 3 + l] += sum_bin(base_intra + k * 3 + l);
    }

  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l) {
      obs.coulomb[k * 3 + l] += sum_bin(layout.coulomb_idx() + k * 3 + l);
      obs.dipolar[k * 3 + l] += sum_bin(layout.dipolar_idx() + k * 3 + l);
      obs.dpd[k * 3 + l] += sum_bin(layout.dpd_idx() + k * 3 + l);
      obs.kinetic_lin[k * 3 + l] += sum_bin(layout.kinetic_idx() + k * 3 + l);
    }
}
