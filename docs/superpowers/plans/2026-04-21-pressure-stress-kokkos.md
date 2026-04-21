# Pressure/Stress Tensor Kokkos/Cabana Conversion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `short_range_loop` in `System::calculate_pressure()` with the same Kokkos/Cabana pair-loop infrastructure used by `calculate_energy()`, integrating DPD pressure into the main pair loop and switching the kinetic virial to `reduce_over_local_particles`.

**Architecture:** Two new header files (`pressure_cabana.hpp`, `bond_pressure_kokkos.hpp`) mirror `energy_cabana.hpp` / `bond_energy_kokkos.hpp` but store 9-component flat tensors per bin instead of scalars. `pressure.cpp` is restructured to match `energy.cpp`. Three existing helper functions get overloads taking AoSoA-available data (positions, charges, directors) instead of `Particle const &`.

**Tech Stack:** C++23, Kokkos (OpenMP backend), Cabana, ESPResSo CellStructure/AoSoA, `Observable_stat` with `chunk_size=9`.

**Spec:** `docs/superpowers/specs/2026-04-21-pressure-stress-kokkos-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/core/pressure_cabana.hpp` | `PressureBinLayout`, `PressureKernel`, `reduce_cabana_pressure` |
| Create | `src/core/bond_pressure_kokkos.hpp` | `BondsPressureKernelData`, `PairBondsPressureKernel`, `AngleBondsPressureKernel`, `NullDihedralPressureKernel` |
| Modify | `src/core/pressure_inline.hpp` | Add position-based overloads; remove dead `add_kinetic_virials` and `add_non_bonded_pair_virials` |
| Modify | `src/core/forces_inline.hpp` | Add director-based `calc_non_central_force` overload for Cabana Gay-Berne |
| Modify | `src/core/pressure.cpp` | Replace `short_range_loop` + `dpd_pressure_local` + kinetic loop with Cabana path |

---

## Task 1: Add position-based overloads to `pressure_inline.hpp`

The existing `calc_bonded_virial_pressure_tensor` and `calc_bonded_three_body_pressure_tensor` take `Particle const &` to get positions and charges. The bond Kokkos kernels need versions taking raw positions (and a charge product for the pair version), since the AoSoA has positions and charges but not full particle structs in device code.

**Files:**
- Modify: `src/core/pressure_inline.hpp`

- [ ] **Step 1: Open `src/core/pressure_inline.hpp` and locate the two existing virial functions**

  The functions start at approximately lines 104 and 124. Read the file to confirm exact locations.
  ```bash
  grep -n "calc_bonded_virial_pressure_tensor\|calc_bonded_three_body_pressure_tensor" src/core/pressure_inline.hpp
  ```

- [ ] **Step 2: Add position-based overload of `calc_bonded_virial_pressure_tensor` after the existing one**

  Insert after the closing `}` of `calc_bonded_virial_pressure_tensor` (the Particle-based version):
  ```cpp
  // Overload for Cabana kernels: takes positions and charge product directly
  inline std::optional<Utils::Matrix<double, 3, 3>>
  calc_bonded_virial_pressure_tensor(
      Bonded_IA_Parameters const &iaparams,
      Utils::Vector3d const &pos1, Utils::Vector3d const &pos2,
      BoxGeometry const &box_geo,
      Coulomb::ShortRangeForceKernel::kernel_type const *kernel,
      double q1q2) {
    auto const dx = box_geo.get_mi_vector(pos1, pos2);
    auto const pair_force = calc_bond_pair_force(iaparams, dx,
  #ifdef ESPRESSO_ELECTROSTATICS
                                                 q1q2, kernel
  #else
                                                 0.0, nullptr
  #endif
    );
    std::optional<Utils::Matrix<double, 3, 3>> pressure{std::nullopt};
    if (pair_force) {
      pressure = Utils::tensor_product(*pair_force, dx);
    }
    return pressure;
  }
  ```

- [ ] **Step 3: Add position-based overload of `calc_bonded_three_body_pressure_tensor` after the existing one**

  Insert after the closing `}` of `calc_bonded_three_body_pressure_tensor` (the Particle-based version):
  ```cpp
  // Overload for Cabana kernels: takes positions directly
  inline std::optional<Utils::Matrix<double, 3, 3>>
  calc_bonded_three_body_pressure_tensor(
      Bonded_IA_Parameters const &iaparams,
      Utils::Vector3d const &pos1, Utils::Vector3d const &pos2,
      Utils::Vector3d const &pos3,
      BoxGeometry const &box_geo) {
    if (std::holds_alternative<AngleHarmonicBond>(iaparams) or
        std::holds_alternative<AngleCosineBond>(iaparams) or
  #ifdef ESPRESSO_TABULATED
        std::holds_alternative<TabulatedAngleBond>(iaparams) or
  #endif
        std::holds_alternative<AngleCossquareBond>(iaparams)) {
      auto const dx21 = -box_geo.get_mi_vector(pos1, pos2);
      auto const dx31 = box_geo.get_mi_vector(pos3, pos1);

      auto const result = calc_bonded_three_body_force(iaparams, dx21, dx31);
      if (result) {
        Utils::Vector3d force2, force3;
        std::tie(std::ignore, force2, force3) = result.value();
        return Utils::tensor_product(force2, dx21) +
               Utils::tensor_product(force3, dx31);
      }
    } else {
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.index()) +
                                 " in pressure calculation.";
      return Utils::Matrix<double, 3, 3>{};
    }
    return {};
  }
  ```

- [ ] **Step 4: Build to verify compilation**

  ```bash
  cd build && make -j$(nproc) 2>&1 | grep -E "error:|warning:" | head -30
  ```
  Expected: no errors in `pressure_inline.hpp`.

- [ ] **Step 5: Commit**

  ```bash
  git add src/core/pressure_inline.hpp
  git commit -m "pressure: add position-based overloads of bonded virial helpers"
  ```

---

## Task 2: Add director-based `calc_non_central_force` overload to `forces_inline.hpp`

The energy Cabana kernel calls `gb_pair_energy(dir1, dir2, ...)` using directors from the AoSoA. For pressure we need forces, not energies. `gb_pair_force(Vector3d, Vector3d, ...)` already exists in `gay_berne.hpp`. We add a `calc_non_central_force` overload taking directors so the pressure kernel can compute Gay-Berne virials without falling back to `unique_particles`.

**Files:**
- Modify: `src/core/forces_inline.hpp`

- [ ] **Step 1: Locate the existing `calc_non_central_force(Particle, Particle, ...)` in `forces_inline.hpp`**

  ```bash
  grep -n "calc_non_central_force" src/core/forces_inline.hpp
  ```

- [ ] **Step 2: Add director-based overload immediately after the existing `calc_non_central_force`**

  ```cpp
  // Director-based overload for use in Cabana pressure kernel.
  // Returns only the force vector (not torque) — sufficient for the virial.
  inline Utils::Vector3d
  calc_non_central_force(Utils::Vector3d const &dir1,
                         Utils::Vector3d const &dir2,
                         IA_parameters const &ia_params,
                         Utils::Vector3d const &d,
                         double const dist) {
    Utils::Vector3d f{};
  #ifdef ESPRESSO_GAY_BERNE
    f += gb_pair_force(dir1, dir2, ia_params, d, dist).f;
  #endif
    return f;
  }
  ```
  Note: `gb_pair_force(Vector3d, Vector3d, ...)` already exists in
  `src/core/nonbonded_interactions/gay_berne.hpp` (line ~49). The `.f` member
  extracts the force component from `ParticleForce`.

- [ ] **Step 3: Build to verify compilation**

  ```bash
  cd build && make -j$(nproc) 2>&1 | grep -E "error:|warning:" | head -30
  ```
  Expected: no errors.

- [ ] **Step 4: Commit**

  ```bash
  git add src/core/forces_inline.hpp
  git commit -m "forces: add director-based calc_non_central_force overload for pressure Cabana kernel"
  ```

---

## Task 3: Create `src/core/pressure_cabana.hpp`

This file contains `PressureBinLayout` (bin-index arithmetic, identical to `EnergyBinLayout` but with a DPD bin and a `tensor_offset(bin, k)` helper), `PressureKernel` (pair non-bonded + Coulomb + DPD virials), and `reduce_cabana_pressure` (sums thread-local bins into `Observable_stat`).

**Files:**
- Create: `src/core/pressure_cabana.hpp`

- [ ] **Step 1: Create the file with header guard and includes**

  ```cpp
  /*
   * Copyright (C) 2026 The ESPResSo project
   * ...license boilerplate matching energy_cabana.hpp...
   */

  #pragma once

  #include <config/config.hpp>

  #include "BoxGeometry.hpp"
  #include "Observable_stat.hpp"
  #include "Particle.hpp"
  #include "aosoa_pack.hpp"
  #include "bonded_interactions/bonded_interaction_data.hpp"
  #include "electrostatics/coulomb.hpp"
  #include "exclusions.hpp"
  #include "forces_inline.hpp"
  #include "nonbonded_interactions/nonbonded_interaction_data.hpp"
  #include "pressure_inline.hpp"
  #include "short_range_cabana_helpers.hpp"

  #ifdef ESPRESSO_DPD
  #include "dpd.hpp"
  #endif

  #include <utils/Vector.hpp>
  #include <utils/math/tensor_product.hpp>

  #include <Cabana_Core.hpp>

  #include <omp.h>

  #include <cstddef>
  #include <span>
  #include <vector>
  ```

- [ ] **Step 2: Add `PressureBinLayout`**

  ```cpp
  struct PressureBinLayout {
    std::size_t n_bonded;
    std::size_t n_types;
    std::size_t off_bonded = 0;
    std::size_t off_nb_inter;
    std::size_t off_nb_intra;
    std::size_t off_coulomb;
    std::size_t off_dipolar;
    std::size_t off_dpd;
    std::size_t total;

    PressureBinLayout(std::size_t n_bonded_, std::size_t n_types_)
        : n_bonded(n_bonded_), n_types(n_types_) {
      auto const n_nb = n_types * (n_types + 1) / 2;
      off_nb_inter = off_bonded + n_bonded;
      off_nb_intra = off_nb_inter + n_nb;
      off_coulomb  = off_nb_intra + n_nb;
      off_dipolar  = off_coulomb + 1;
      off_dpd      = off_dipolar + 1;
      total        = off_dpd + 1;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t tensor_offset(std::size_t bin, std::size_t k) const {
      return bin * 9 + k;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t nb_inter_idx(int t1, int t2) const {
      return off_nb_inter +
             Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t nb_intra_idx(int t1, int t2) const {
      return off_nb_intra +
             Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
    }

    KOKKOS_INLINE_FUNCTION std::size_t coulomb_idx() const { return off_coulomb; }
    KOKKOS_INLINE_FUNCTION std::size_t dipolar_idx() const { return off_dipolar; }
    KOKKOS_INLINE_FUNCTION std::size_t dpd_idx()     const { return off_dpd; }
    KOKKOS_INLINE_FUNCTION std::size_t bonded_idx(int b) const {
      return off_bonded + static_cast<std::size_t>(b);
    }
  };
  ```

- [ ] **Step 3: Add `PressureKernel`**

  ```cpp
  struct PressureKernel {
    BondedInteractionsMap const &bonded_ias;
    InteractionsNonBonded const &nonbonded_ias;
    Coulomb::Solver const &coulomb;
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel;
    Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel;
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
        BoxGeometry const &box_geo_,
        std::vector<Particle *> const &unique_particles_,
        Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
        PressureBinLayout layout_,
        CellStructure::AoSoA_pack const &aosoa_,
        Kokkos::View<int *> mol_id_view_,
        double system_max_cutoff_)
        : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
          coulomb(coulomb_), coulomb_f_kernel(coulomb_f_kernel_),
          coulomb_p_kernel(coulomb_p_kernel_), box_geo(box_geo_),
          unique_particles(unique_particles_), local_pressure(local_pressure_),
          layout(layout_), aosoa(aosoa_),
          mol_id_view(std::move(mol_id_view_)),
          system_max_cutoff(system_max_cutoff_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(std::size_t i, std::size_t j) const {
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const d    = box_geo.get_mi_vector(pos1, pos2);
      auto const dist = d.norm();
      if (dist > system_max_cutoff)
        return;

      auto const t1 = aosoa.type(i);
      auto const t2 = aosoa.type(j);
      auto const &ia_params = nonbonded_ias.get_ia_param(t1, t2);
      auto const tid = omp_get_thread_num();

      // --- Non-bonded virial ---
      if (dist <= ia_params.max_cut) {
        bool skip = false;
  #ifdef ESPRESSO_EXCLUSIONS
        if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j)) {
          skip = not do_nonbonded(*unique_particles.at(i),
                                  *unique_particles.at(j));
        }
  #endif
        if (not skip) {
          Utils::Vector3d f = calc_central_radial_force(ia_params, d, dist);

  #ifdef ESPRESSO_GAY_BERNE
          if (gay_berne_active(dist, ia_params)) {
            auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
            auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
            f += calc_non_central_force(dir1, dir2, ia_params, d, dist);
          }
  #endif
  #ifdef ESPRESSO_THOLE
          if (thole_active(ia_params, coulomb_f_kernel != nullptr)) {
            f += thole_pair_force(*unique_particles.at(i),
                                  *unique_particles.at(j),
                                  ia_params, d, dist,
                                  bonded_ias, coulomb, coulomb_f_kernel);
          }
  #endif

          auto const stress = Utils::flatten(Utils::tensor_product(d, f));
          auto const bin = (mol_id_view(i) == mol_id_view(j))
                               ? layout.nb_intra_idx(t1, t2)
                               : layout.nb_inter_idx(t1, t2);
          for (std::size_t k = 0; k < 9; ++k)
            local_pressure(tid, layout.tensor_offset(bin, k)) += stress[k];
        }
      }

  #ifdef ESPRESSO_ELECTROSTATICS
      if (coulomb_p_kernel != nullptr) {
        auto const q1 = aosoa.charge(i), q2 = aosoa.charge(j);
        if (q1 != 0. and q2 != 0.) {
          auto const p_c = Utils::flatten((*coulomb_p_kernel)(q1 * q2, d, dist));
          for (std::size_t k = 0; k < 9; ++k)
            local_pressure(tid,
                layout.tensor_offset(layout.coulomb_idx(), k)) += p_c[k];
        }
      }
  #endif

  #ifdef ESPRESSO_DPD
      {
        auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
        auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
        auto const v21  = box_geo.velocity_difference(pos1, pos2, vel1, vel2);
        auto const dist2 = d.norm2();
        auto const f_r = dpd_pair_force(ia_params.dpd.radial, v21, dist, {});
        auto const f_t = dpd_pair_force(ia_params.dpd.trans,  v21, dist, {});
        auto const P   = Utils::tensor_product(d / dist2, d);
        auto const f_d = P * (f_r - f_t) + f_t;
        auto const s   = Utils::flatten(Utils::tensor_product(d, f_d));
        for (std::size_t k = 0; k < 9; ++k)
          local_pressure(tid,
              layout.tensor_offset(layout.dpd_idx(), k)) += s[k];
      }
  #endif
    }
  };
  ```

- [ ] **Step 4: Add `reduce_cabana_pressure`**

  ```cpp
  static void reduce_cabana_pressure(
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure,
      PressureBinLayout const &layout, Observable_stat &obs,
      BondedInteractionsMap const &bonded_ias, int n_types) {
    auto const nthreads = local_pressure.extent(0);
    auto host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, local_pressure);

    auto sum_tensor = [&](std::size_t bin, std::span<double> dest) {
      for (std::size_t k = 0; k < 9; ++k) {
        double s = 0.;
        for (std::size_t t = 0; t < nthreads; ++t)
          s += host(t, layout.tensor_offset(bin, k));
        dest[k] += s;
      }
    };

    for (int b = 0; b < int(layout.n_bonded); ++b)
      sum_tensor(layout.bonded_idx(b), obs.bonded_contribution(b));

    for (int t1 = 0; t1 < n_types; ++t1)
      for (int t2 = 0; t2 <= t1; ++t2)
        sum_tensor(layout.nb_inter_idx(t1, t2),
                   obs.non_bonded_inter_contribution(t1, t2));

    for (int t1 = 0; t1 < n_types; ++t1)
      for (int t2 = 0; t2 <= t1; ++t2)
        sum_tensor(layout.nb_intra_idx(t1, t2),
                   obs.non_bonded_intra_contribution(t1, t2));

    if (!obs.coulomb.empty())
      sum_tensor(layout.coulomb_idx(), obs.coulomb.subspan(0, 9));
    if (!obs.dipolar.empty())
      sum_tensor(layout.dipolar_idx(), obs.dipolar.subspan(0, 9));
    if (!obs.dpd.empty())
      sum_tensor(layout.dpd_idx(), obs.dpd);
  }
  ```

- [ ] **Step 5: Build to verify compilation**

  ```bash
  cd build && make -j$(nproc) 2>&1 | grep -E "error:" | head -30
  ```
  The file is not yet `#include`d anywhere so this only catches syntax errors
  if another TU includes it. That's fine — full build validation happens in Task 6.

- [ ] **Step 6: Commit**

  ```bash
  git add src/core/pressure_cabana.hpp
  git commit -m "pressure: add PressureBinLayout, PressureKernel, reduce_cabana_pressure"
  ```

---

## Task 4: Create `src/core/bond_pressure_kokkos.hpp`

Bond kernels for pair and angle bonds. Mirrors `bond_energy_kokkos.hpp` but calls the position-based virial overloads from Task 1 and stores 9-component tensors. Includes a no-op dihedral kernel struct that satisfies the `cabana_short_range` interface.

**Files:**
- Create: `src/core/bond_pressure_kokkos.hpp`

- [ ] **Step 1: Create file with header and `BondsPressureKernelData`**

  ```cpp
  /*
   * Copyright (C) 2026 The ESPResSo project
   * ...license boilerplate...
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

  #include <array>
  #include <cstddef>

  struct BondsPressureKernelData {
    BondedInteractionsMap const &bonded_ias;
    BoxGeometry const &box_geo;
    Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
    PressureBinLayout layout;
    CellStructure::AoSoA_pack const &aosoa;
  };
  ```

- [ ] **Step 2: Add `PairBondsPressureKernel`**

  ```cpp
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
      auto const &bonded_ias   = data.bonded_ias;
      auto const &box_geo      = data.box_geo;
      auto &local_pressure     = data.local_pressure;
      auto const &layout       = data.layout;
      auto const &aosoa        = data.aosoa;
      auto const bond_id       = bond_ids(idx);
      auto const thread_id     = omp_get_thread_num();

      auto const i = bond_list(idx, 0);
      auto const j = bond_list(idx, 1);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);

      auto const pressure = calc_bonded_virial_pressure_tensor(
          iaparams, pos1, pos2, box_geo, coulomb_f_kernel,
  #ifdef ESPRESSO_ELECTROSTATICS
          aosoa.charge(i) * aosoa.charge(j)
  #else
          0.0
  #endif
      );

      if (pressure) {
        auto const flat = Utils::flatten(*pressure);
        for (std::size_t k = 0; k < 9; ++k)
          local_pressure(thread_id,
              layout.tensor_offset(layout.bonded_idx(bond_id), k)) += flat[k];
      } else {
        auto partner_id = aosoa.id(j);
        bond_broken_error(aosoa.id(i), {&partner_id, 1});
      }
    }
  };
  ```

- [ ] **Step 3: Add `AngleBondsPressureKernel`**

  ```cpp
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
      auto const &bonded_ias   = data.bonded_ias;
      auto const &box_geo      = data.box_geo;
      auto &local_pressure     = data.local_pressure;
      auto const &layout       = data.layout;
      auto const &aosoa        = data.aosoa;
      auto const bond_id       = bond_ids(idx);
      auto const thread_id     = omp_get_thread_num();

      auto const i   = bond_list(idx, 0);
      auto const j   = bond_list(idx, 1);
      auto const k   = bond_list(idx, 2);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);

      auto const pressure =
          calc_bonded_three_body_pressure_tensor(iaparams, pos1, pos2, pos3, box_geo);

      if (pressure) {
        auto const flat = Utils::flatten(*pressure);
        for (std::size_t k2 = 0; k2 < 9; ++k2)
          local_pressure(thread_id,
              layout.tensor_offset(layout.bonded_idx(bond_id), k2)) += flat[k2];
      } else {
        std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
        bond_broken_error(aosoa.id(i), {pids.data(), 2});
      }
    }
  };
  ```

- [ ] **Step 4: Add `NullDihedralPressureKernel`**

  ```cpp
  struct NullDihedralPressureKernel {
    KOKKOS_INLINE_FUNCTION void operator()(std::size_t) const {}
  };
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/core/bond_pressure_kokkos.hpp
  git commit -m "pressure: add Kokkos bond pressure kernels (pair, angle, null dihedral)"
  ```

---

## Task 5: Update `pressure.cpp` to use Cabana infrastructure

Replace `short_range_loop` + `dpd_pressure_local` + plain kinetic loop with the new Cabana path. The file structure mirrors `energy.cpp`.

**Files:**
- Modify: `src/core/pressure.cpp`

- [ ] **Step 1: Update includes**

  Replace:
  ```cpp
  #include "dpd.hpp"
  #include "pressure_inline.hpp"
  #include "short_range_loop.hpp"
  ```
  With:
  ```cpp
  #include "bond_pressure_kokkos.hpp"
  #include "integrators/Propagation.hpp"
  #include "nonbonded_interactions/VerletCriterion.hpp"
  #include "particle_reduction.hpp"
  #include "pressure_cabana.hpp"
  #include "pressure_inline.hpp"
  #include "short_range_cabana.hpp"
  ```
  Keep all other existing includes (`BoxGeometry.hpp`, `Observable_stat.hpp`, etc.).

- [ ] **Step 2: Replace the entire body of `System::calculate_pressure()`**

  The new body (replace everything after the opening `{` through the closing `}`):
  ```cpp
  std::shared_ptr<Observable_stat> System::calculate_pressure() {

    auto obs_pressure_ptr = std::make_shared<Observable_stat>(
        9ul, static_cast<std::size_t>(bonded_ias->get_next_key()),
        nonbonded_ias->get_max_seen_particle_type());

    if (long_range_interactions_sanity_checks()) {
      return obs_pressure_ptr;
    }

    auto &obs_pressure = *obs_pressure_ptr;

    on_observable_calc();

    auto const volume = box_geo->volume();

    // Kinetic virial — OpenMP-parallel reduction over local particles
    auto const kinetic = reduce_over_local_particles<Utils::Matrix<double, 3, 3>>(
        *cell_structure,
        [](Utils::Matrix<double, 3, 3> &acc, Particle const &p) {
          if (!p.is_virtual())
            acc += Utils::tensor_product(p.v(), p.mass() * p.v());
        },
        [](auto &a, auto const &b) { a = a + b; });
    std::ranges::copy(Utils::flatten(kinetic), obs_pressure.kinetic_lin.begin());

    auto const coulomb_force_kernel    = coulomb.pair_force_kernel();
    auto const coulomb_pressure_kernel = coulomb.pair_pressure_kernel();

    VerletCriterion<> const verlet_criterion{*this,
                                             cell_structure->get_verlet_skin(),
                                             get_interaction_range(),
                                             coulomb.cutoff(),
                                             dipoles.cutoff(),
                                             inactive_cutoff};
    update_cabana_state(*cell_structure, verlet_criterion,
                        get_interaction_range(), propagation->integ_switch);

    PressureBinLayout layout{
        static_cast<std::size_t>(bonded_ias->get_next_key()),
        std::size_t(nonbonded_ias->get_max_seen_particle_type() + 1)};

    using exec = Kokkos::DefaultExecutionSpace;
    Kokkos::View<double **, Kokkos::LayoutRight> local_pressure(
        "local_pressure", exec().concurrency(), layout.total * 9);

    auto const &unique_particles = cell_structure->get_unique_particles();
    auto const n_particles = static_cast<int>(unique_particles.size());
    Kokkos::View<int *> mol_id_view("mol_id", n_particles);
    auto mol_id_host = Kokkos::create_mirror_view(mol_id_view);
    for (int i = 0; i < n_particles; ++i)
      mol_id_host(i) = unique_particles[i]->mol_id();
    Kokkos::deep_copy(mol_id_view, mol_id_host);

    PressureKernel pair_p_kernel{*bonded_ias,
                                 *nonbonded_ias,
                                 coulomb,
                                 get_ptr(coulomb_force_kernel),
                                 get_ptr(coulomb_pressure_kernel),
                                 *box_geo,
                                 cell_structure->get_unique_particles(),
                                 local_pressure,
                                 layout,
                                 cell_structure->get_aosoa(),
                                 mol_id_view,
                                 maximal_cutoff()};

    auto &bs = cell_structure->bond_state();
    BondsPressureKernelData bonds_p_data{*bonded_ias, *box_geo, local_pressure,
                                         layout, cell_structure->get_aosoa()};
    PairBondsPressureKernel pair_bp_kernel{bonds_p_data, bs.pair_list,
                                           bs.pair_ids,
                                           get_ptr(coulomb_force_kernel)};
    AngleBondsPressureKernel angle_bp_kernel{bonds_p_data, bs.angle_list,
                                             bs.angle_ids};
    NullDihedralPressureKernel null_dih_kernel{};

    cabana_short_range(pair_bp_kernel, angle_bp_kernel, null_dih_kernel,
                       pair_p_kernel, *cell_structure, get_interaction_range(),
                       bonded_ias->maximal_cutoff(), verlet_criterion,
                       propagation->integ_switch);

    reduce_cabana_pressure(local_pressure, layout, obs_pressure, *bonded_ias,
                           nonbonded_ias->get_max_seen_particle_type() + 1);

  #ifdef ESPRESSO_ELECTROSTATICS
    auto const coulomb_pressure = coulomb.calc_pressure_long_range();
    std::ranges::copy(coulomb_pressure, obs_pressure.coulomb.begin() + 9u);
  #endif
  #ifdef ESPRESSO_DIPOLES
    dipoles.calc_pressure_long_range();
  #endif

  #ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    if (!obs_pressure.virtual_sites.empty()) {
      auto const vs_pressure = vs_relative_pressure_tensor(*cell_structure);
      std::ranges::copy(Utils::flatten(vs_pressure),
                        obs_pressure.virtual_sites.begin());
    }
  #endif

    obs_pressure.rescale(volume);

    obs_pressure.mpi_reduce();
    return obs_pressure_ptr;
  }
  ```

- [ ] **Step 3: Build**

  ```bash
  cd build && make -j$(nproc) 2>&1 | grep -E "error:" | head -40
  ```
  Fix any compile errors before continuing. Common issues:
  - Missing include for `update_cabana_state` → add `#include "short_range_cabana.hpp"` if not present
  - Missing `inactive_cutoff` → check how `energy.cpp` accesses it (it's a member of `System`)
  - `Utils::flatten` on `Utils::Matrix<double,3,3>` → confirm it's in `<utils/math/tensor_product.hpp>`

- [ ] **Step 4: Commit**

  ```bash
  git add src/core/pressure.cpp
  git commit -m "pressure: replace short_range_loop with Cabana pair loop; integrate DPD virial"
  ```

---

## Task 6: Remove dead functions from `pressure_inline.hpp`

`add_kinetic_virials` and `add_non_bonded_pair_virials` are no longer called from anywhere in the core (`dpd_stress()` in `dpd.cpp` uses `cell_structure.non_bonded_loop()` directly, not `add_non_bonded_pair_virials`).

**Files:**
- Modify: `src/core/pressure_inline.hpp`

- [ ] **Step 1: Verify `add_kinetic_virials` is not called anywhere**

  ```bash
  grep -r "add_kinetic_virials" src/ --include="*.cpp" --include="*.hpp"
  ```
  Expected: no matches.

- [ ] **Step 2: Verify `add_non_bonded_pair_virials` is not called anywhere**

  ```bash
  grep -r "add_non_bonded_pair_virials" src/ --include="*.cpp" --include="*.hpp"
  ```
  Expected: no matches.

- [ ] **Step 3: Remove both functions from `pressure_inline.hpp`**

  Delete the `add_kinetic_virials` function (the one that takes `Particle const &p1, Observable_stat &obs_pressure` and writes to `obs_pressure.kinetic_lin`).

  Delete the `add_non_bonded_pair_virials` function (the one that takes two `Particle const &`, a distance vector, `IA_parameters`, etc.).

  Leave in place: `calc_bonded_virial_pressure_tensor` (both overloads), `calc_bonded_three_body_pressure_tensor` (both overloads), `calc_bonded_pressure_tensor`.

- [ ] **Step 4: Build to verify no remaining uses**

  ```bash
  cd build && make -j$(nproc) 2>&1 | grep -E "error:" | head -20
  ```
  Expected: clean build.

- [ ] **Step 5: Commit**

  ```bash
  git add src/core/pressure_inline.hpp
  git commit -m "pressure: remove dead add_kinetic_virials and add_non_bonded_pair_virials"
  ```

---

## Task 7: Build full binary and run tests

Verify the Cabana pressure path produces numerically identical results to the previous implementation.

**Files:**
- Test: `testsuite/python/pressure.py`
- Test: `testsuite/python/dpd.py`

- [ ] **Step 1: Full build**

  ```bash
  cd build && make -j$(nproc)
  ```
  Expected: zero errors, zero new warnings.

- [ ] **Step 2: Run the pressure test suite**

  ```bash
  cd build && ./pypresso ../testsuite/python/pressure.py 2>&1
  ```
  Expected: all tests pass. If tests print per-component pressure tensor values,
  verify they are within floating-point tolerance of pre-change values.
  If the test runner reports FAIL, read the assertion message — most pressure
  tests compare `system.analysis.pressure()` output against analytic or
  reference values.

- [ ] **Step 3: Run the DPD test**

  ```bash
  cd build && ./pypresso ../testsuite/python/dpd.py 2>&1
  ```
  Expected: all tests pass. This test exercises `system.analysis.pressure()`
  with DPD active, confirming the integrated DPD virial is correct.

- [ ] **Step 4: Run the DPD stats test**

  ```bash
  cd build && ./pypresso ../testsuite/python/dpd_stats.py 2>&1
  ```
  Expected: all tests pass.

- [ ] **Step 5: Commit (if any test fixes were needed)**

  If any numerical tolerances or test infrastructure changes were needed,
  commit them now:
  ```bash
  git add -p
  git commit -m "pressure: fix test tolerances after Cabana conversion"
  ```
  If all tests passed without changes, no additional commit needed.
