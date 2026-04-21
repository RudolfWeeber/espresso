# Pressure Tensor Kokkos/Cabana Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the pressure/stress tensor calculation from CPU `short_range_loop()` to Kokkos/Cabana, mirroring the energy calculation architecture.

**Architecture:** Create `PressureBinLayout` + `PressureKernel` + `KineticPressureKernel` in `pressure_cabana.hpp` and bonded pressure kernels in `bond_pressure_kokkos.hpp`. Replace the `short_range_loop()` call in `pressure.cpp` with `cabana_short_range()` + Kokkos kinetic loop. Each bin holds 9 consecutive elements (flattened 3x3 tensor). DPD stress is computed inside `PressureKernel` with zero noise. Virtual sites and Coulomb k-space remain CPU calls.

**Tech Stack:** C++17, Kokkos, Cabana, OpenMP (thread-local accumulation via `omp_get_thread_num()`)

**Spec:** `docs/superpowers/specs/2026-04-21-pressure-tensor-kokkos-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `src/core/aosoa_pack.hpp` | Add `is_virtual` flag to AoSoA |
| Modify | `src/core/short_range_cabana.hpp` | Commit `is_virtual` flag in `commit_particle()` |
| Create | `src/core/pressure_cabana.hpp` | `PressureBinLayout`, `PressureKernel`, `KineticPressureKernel`, `reduce_cabana_pressure()` |
| Create | `src/core/bond_pressure_kokkos.hpp` | `BondsPressureKernelData`, `PairBondsPressureKernel`, `AngleBondsPressureKernel`, `DihedralBondsPressureKernel` |
| Modify | `src/core/pressure.cpp` | Replace `short_range_loop()` with Cabana path |

### Key Reference Files (read before implementing)

| File | Why |
|------|-----|
| `src/core/energy_cabana.hpp` | Template for `PressureBinLayout` and `PressureKernel` |
| `src/core/bond_energy_kokkos.hpp` | Template for bonded pressure kernels |
| `src/core/forces_cabana.hpp` | Force computation pattern, DPD integration, Thole/GB |
| `src/core/short_range_cabana.hpp` | `cabana_short_range()` signature, `update_cabana_state()` |
| `src/core/short_range_cabana_helpers.hpp` | `dpd_active()`, `thole_active()`, `gay_berne_active()`, `compute_pair_data_flags()` |
| `src/core/pressure_inline.hpp` | Current pressure inline helpers (force-based pressure calculations) |
| `src/core/pressure.cpp` | Current `calculate_pressure()` to be modified |
| `src/core/energy.cpp` | Reference for how energy integrates Cabana path |
| `src/core/dpd.hpp` / `dpd.cpp` | DPD pair force signatures and `dpd_pair_force()` overloads |
| `src/core/electrostatics/coulomb.hpp` | `ShortRangePressureKernel` type alias |
| `src/core/electrostatics/coulomb_inline.hpp` | `pair_pressure_kernel()` implementation |

---

### Task 0: Add is_virtual flag and first_ghost_idx to AoSoA/CellStructure

The `KineticPressureKernel` needs to skip virtual particles (the CPU path uses `p1.is_virtual()`) and ghost particles (the CPU path only iterates `local_parts`, not ghosts). The AoSoA currently has no `is_virtual` flag, and there's no way to know where ghost particles start in the AoSoA.

**Files:**
- Modify: `src/core/aosoa_pack.hpp`
- Modify: `src/core/short_range_cabana.hpp`
- Modify: `src/core/cell_system/CellStructure.hpp`
- Modify: `src/core/cell_system/CellStructure.cpp`

- [ ] **Step 1: Extend flags encoding in aosoa_pack.hpp**

The current `flags` field uses `uint8_t` with values 0 or 1 (exclusion only). Extend to a bitfield:

- Bit 0: `has_exclusion` (existing)
- Bit 1: `is_virtual` (new)

Replace the flag access methods in `aosoa_pack.hpp`:

```cpp
  void set_has_exclusion(std::size_t i, bool value) {
    if (value)
      flags(i) |= uint8_t{1};
    else
      flags(i) &= ~uint8_t{1};
  }

  bool has_exclusion(std::size_t i) const { return (flags(i) & uint8_t{1}) != 0; }

  void set_is_virtual(std::size_t i, bool value) {
    if (value)
      flags(i) |= uint8_t{2};
    else
      flags(i) &= ~uint8_t{2};
  }

  bool is_virtual(std::size_t i) const { return (flags(i) & uint8_t{2}) != 0; }
```

- [ ] **Step 2: Add first_ghost_idx to CellStructure**

In `src/core/cell_system/CellStructure.hpp`, add a new member variable alongside the existing cached values (around line 197):

```cpp
  std::size_t m_first_ghost_idx = 0;
```

And add a getter (around line 448, next to `get_num_local_particles_cached()`):

```cpp
  std::size_t get_first_ghost_idx() const { return m_first_ghost_idx; }
```

- [ ] **Step 3: Set first_ghost_idx in CellStructure::set_index_map()**

In `src/core/cell_system/CellStructure.cpp`, in `set_index_map()`, add after line 280 (after the `Kokkos::fence()` following `enumerate_local_particles`) and before the ghost particle loop (line 290):

```cpp
  m_first_ghost_idx = count_local_particles();
```

This captures the index where ghost particles begin in `unique_particles`. Local particles occupy indices `[0, m_first_ghost_idx)` and ghost particles occupy `[m_first_ghost_idx, unique_particles.size())`.

- [ ] **Step 4: Commit is_virtual flag in commit_particle() in short_range_cabana.hpp**

In `commit_particle()`, after the exclusion flag setting, add the virtual flag:

```cpp
  // Always update exclusion flags (they can change during simulation)
#ifdef ESPRESSO_EXCLUSIONS
  aosoa.set_has_exclusion(index, !p.exclusions().empty());
#else
  aosoa.flags(index) = 0;
#endif

  // Always update virtual flag
  aosoa.set_is_virtual(index, p.is_virtual());
```

Note: `p.is_virtual()` is available when `ESPRESSO_VIRTUAL_SITES` is defined, and returns `false` otherwise (constexpr). So this line is always safe to call.

- [ ] **Step 5: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

- [ ] **Step 6: Run existing pressure test to verify no regression**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/pressure.py`

Expected: All tests pass. The flags change is backward-compatible because bit 0 still correctly represents exclusion. `first_ghost_idx` is a new field not yet used by existing code.

- [ ] **Step 7: Commit**

```bash
git add src/core/aosoa_pack.hpp src/core/short_range_cabana.hpp src/core/cell_system/CellStructure.hpp src/core/cell_system/CellStructure.cpp
git commit -m "feat: add is_virtual flag and first_ghost_idx for Kokkos pressure calculation"
```

---

### Task 1: Create PressureBinLayout in pressure_cabana.hpp

**Files:**
- Create: `src/core/pressure_cabana.hpp`

- [ ] **Step 1: Write PressureBinLayout struct**

Create `src/core/pressure_cabana.hpp` with the license header (copy from `energy_cabana.hpp`), includes, and `PressureBinLayout` struct.

```cpp
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
```

Key differences from `EnergyBinLayout`:
- Each category offset is multiplied by 9 (9 elements per tensor)
- DPD and kinetic bins added
- No external_fields bin
- Index functions return the *first* of 9 flat indices

- [ ] **Step 2: Verify the file compiles**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

This won't link yet (nothing includes it), but syntax errors will show.

- [ ] **Step 3: Commit**

```bash
git add src/core/pressure_cabana.hpp
git commit -m "feat: add PressureBinLayout for Kokkos pressure calculation"
```

---

### Task 2: Create KineticPressureKernel in pressure_cabana.hpp

**Files:**
- Modify: `src/core/pressure_cabana.hpp`

- [ ] **Step 1: Add KineticPressureKernel struct to pressure_cabana.hpp**

Add after `PressureBinLayout`:

```cpp
struct KineticPressureKernel {
  CellStructure::AoSoA_pack const &aosoa;
  Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  std::size_t n_local_particles;

  KineticPressureKernel(
      CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
      PressureBinLayout layout_,
      std::size_t n_local_particles_)
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
```

Note: `n_local_particles` is `cell_structure->get_first_ghost_idx()`, which marks where ghost particles start in the AoSoA. Only local particles (indices 0 to `n_local_particles-1`) contribute to kinetic pressure. Ghost particles are skipped. `aosoa.is_virtual(idx)` uses the flag added in Task 0 to skip virtual particles, matching the CPU path's `p1.is_virtual()` check.

- [ ] **Step 2: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

- [ ] **Step 3: Commit**

```bash
git add src/core/pressure_cabana.hpp
git commit -m "feat: add KineticPressureKernel for Kokkos pressure calculation"
```

---

### Task 3: Create PressureKernel (non-bonded pair pressure) in pressure_cabana.hpp

**Files:**
- Modify: `src/core/pressure_cabana.hpp`

This is the largest kernel. It mirrors `EnergyKernel` but computes `tensor_product(d, f)` instead of scalar energy, and adds DPD stress.

- [ ] **Step 1: Add PressureKernel struct to pressure_cabana.hpp**

Add after `KineticPressureKernel`. Read `forces_cabana.hpp` (the `ForcesKernel` struct) and `energy_cabana.hpp` (the `EnergyKernel` struct) for the exact pattern.

```cpp
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
      Thermostat::Thermostat const &thermostat_,
      BoxGeometry const &box_geo_,
      std::vector<Particle *> const &unique_particles_,
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
      PressureBinLayout layout_,
      CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<int *> mol_id_view_,
      double system_max_cutoff_)
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
      auto const dist2 = dist * dist;

      Utils::Vector3d f_r = dpd_pair_force(ia_params.dpd.radial, v21, dist,
                                            Utils::Vector3d{});
      Utils::Vector3d f_t = dpd_pair_force(ia_params.dpd.trans, v21, dist,
                                            Utils::Vector3d{});

      auto const P = Utils::tensor_product(d / dist2, d);
      auto const f_dpd = P * (f_r - f_t) + f_t;

      auto const stress_dpd = -Utils::tensor_product(d, f_dpd);
      write_tensor(tid, layout.dpd_idx(), stress_dpd);
    }
#endif
  }
};
```

**Implementation notes:**
- The non-bonded force computation uses `calc_central_radial_force` (same as `ForcesKernel`), not energy
- The non-bonded stress is `tensor_product(d, f_nb)` — this matches `add_non_bonded_pair_virials()` in `pressure_inline.hpp` which computes `tensor_product(d, f)` where f includes central + non-central + Thole + GB forces
- DPD uses the same projection operator pattern as `dpd.cpp:dpd_viscous_stress_local()` but with zero noise (empty `Utils::Vector3d{}`)
- The `write_tensor` helper avoids repeating the `k*3+l` indexing pattern
- Coulomb uses `coulomb_p_kernel` (returns 3x3 matrix) not the energy or force kernel
- The Thole contribution needs the force kernel (`coulomb_f_kernel`), same as `ForcesKernel`
- `dipolar` just prints a warning (matching current behavior)

- [ ] **Step 2: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

- [ ] **Step 3: Commit**

```bash
git add src/core/pressure_cabana.hpp
git commit -m "feat: add PressureKernel for non-bonded pair pressure tensor"
```

---

### Task 4: Add reduce_cabana_pressure() to pressure_cabana.hpp

**Files:**
- Modify: `src/core/pressure_cabana.hpp`

- [ ] **Step 1: Add reduce_cabana_pressure() function**

Add after `PressureKernel` at the end of `pressure_cabana.hpp`:

```cpp
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
```

This mirrors `reduce_cabana_energy()` but writes 9 elements per category into the `Observable_stat` spans. The `Observable_stat` was constructed with `chunk_size=9`, so each `bonded_contribution(b)` returns a span of 9 elements, and `non_bonded_inter_contribution(t1,t2)` returns a span of 9 elements.

- [ ] **Step 2: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

- [ ] **Step 3: Commit**

```bash
git add src/core/pressure_cabana.hpp
git commit -m "feat: add reduce_cabana_pressure for Kokkos pressure reduction"
```

---

### Task 5: Create bond_pressure_kokkos.hpp

**Files:**
- Create: `src/core/bond_pressure_kokkos.hpp`

- [ ] **Step 1: Write the full bonded pressure kernel file**

Create `src/core/bond_pressure_kokkos.hpp`:

```cpp
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
      BondsPressureKernelData data_, LocalBondState::PairBondlistType bond_list_,
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

    auto const pair_force = calc_bond_pair_force(iaparams, dx,
#ifdef ESPRESSO_ELECTROSTATICS
                                                  aosoa.charge(i) * aosoa.charge(j),
                                                  coulomb_f_kernel
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

  DihedralBondsPressureKernel(
      BondsPressureKernelData data_,
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
```

**Implementation notes:**
- `PairBondsPressureKernel` uses `calc_bond_pair_force()` (same as `calc_bonded_virial_pressure_tensor()` in `pressure_inline.hpp`) and then `tensor_product(force, dx)` — this matches the CPU path exactly
- `AngleBondsPressureKernel` uses `calc_bonded_three_body_force()` and the same tensor product decomposition as `calc_bonded_three_body_pressure_tensor()`
- `DihedralBondsPressureKernel` writes zero (matching current CPU behavior where dihedral bonds return zero matrix + warning)
- The angle kernel uses `dx21 = -box_geo.get_mi_vector(pos1, pos2)` and `dx31 = box_geo.get_mi_vector(pos3, pos1)` — matching `calc_bonded_three_body_pressure_tensor()` in `pressure_inline.hpp:135-136`

- [ ] **Step 2: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -20`

- [ ] **Step 3: Commit**

```bash
git add src/core/bond_pressure_kokkos.hpp
git commit -m "feat: add bonded pressure kernels for Kokkos pressure calculation"
```

---

### Task 6: Integrate Cabana pressure path into pressure.cpp

**Files:**
- Modify: `src/core/pressure.cpp`

This is the integration task. Replace `short_range_loop()` with the Cabana path, add kinetic Kokkos loop, remove DPD CPU call.

- [ ] **Step 1: Rewrite pressure.cpp**

Replace the entire `calculate_pressure()` function body. The new includes needed are:

```cpp
#include "bond_pressure_kokkos.hpp"
#include "energy_cabana.hpp"  // for EnergyBinLayout (not used, but for consistency)
#include "pressure_cabana.hpp"
#include "short_range_cabana.hpp"
```

Remove:
```cpp
#include "short_range_loop.hpp"
```

The new `calculate_pressure()`:

```cpp
namespace System {
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

  auto const coulomb_force_kernel = coulomb.pair_force_kernel();
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
      "local_pressure", exec().concurrency(), layout.total);

  auto const &unique_particles = cell_structure->get_unique_particles();
  auto const n_particles = static_cast<int>(unique_particles.size());

  // Kinetic pressure (Kokkos parallel loop over local particles only)
  auto const n_local_particles = cell_structure->get_first_ghost_idx();
  KineticPressureKernel kinetic_kernel{cell_structure->get_aosoa(),
                                       local_pressure, layout,
                                       n_local_particles};
  Kokkos::parallel_for("kinetic_pressure", n_local_particles, kinetic_kernel);
  Kokkos::fence();

  // Build mol_id view
  Kokkos::View<int *> mol_id_view("mol_id", n_particles);
  auto mol_id_host = Kokkos::create_mirror_view(mol_id_view);
  for (int i = 0; i < n_particles; ++i) {
    mol_id_host(i) = unique_particles[i]->mol_id();
  }
  Kokkos::deep_copy(mol_id_view, mol_id_host);

  // Non-bonded pair pressure
  PressureKernel pair_p_kernel{*bonded_ias,
                               *nonbonded_ias,
                               coulomb,
                               get_ptr(coulomb_force_kernel),
                               get_ptr(coulomb_pressure_kernel),
                               dipoles,
                               *propagation->thermostat,
                               *box_geo,
                               cell_structure->get_unique_particles(),
                               local_pressure,
                               layout,
                               cell_structure->get_aosoa(),
                               mol_id_view,
                               maximal_cutoff()};

  // Bonded pressure
  auto &bs = cell_structure->bond_state();
  BondsPressureKernelData bonds_p_data{*bonded_ias, *box_geo, local_pressure,
                                       layout, cell_structure->get_aosoa()};
  PairBondsPressureKernel pair_bp_kernel{bonds_p_data, bs.pair_list,
                                         bs.pair_ids,
                                         get_ptr(coulomb_force_kernel)};
  AngleBondsPressureKernel angle_bp_kernel{bonds_p_data, bs.angle_list,
                                           bs.angle_ids};
  DihedralBondsPressureKernel dih_bp_kernel{bonds_p_data, bs.dihedral_list,
                                            bs.dihedral_ids};

  cabana_short_range(pair_bp_kernel, angle_bp_kernel, dih_bp_kernel,
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

  // DPD is now computed inside PressureKernel, no separate CPU call needed

  obs_pressure.rescale(volume);

  obs_pressure.mpi_reduce();
  return obs_pressure_ptr;
}
} // namespace System
```

Key changes from the original `pressure.cpp`:
1. Replaced `short_range_loop()` with `cabana_short_range()` + `update_cabana_state()`
2. Added `PressureBinLayout` + `local_pressure` view allocation (mirroring `energy.cpp`)
3. Added `KineticPressureKernel` with `Kokkos::parallel_for` instead of serial kinetic loop
4. Added `PressureKernel` + bond pressure kernels + `cabana_short_range()` call
5. Added `reduce_cabana_pressure()` to write bins into `Observable_stat`
6. Removed `#include "short_range_loop.hpp"` and `#include "dpd.hpp"` (DPD is now in-kernel)
7. Added `#include "pressure_cabana.hpp"`, `#include "bond_pressure_kokkos.hpp"`, `#include "short_range_cabana.hpp"`
8. Removed `dpd_pressure_local()` CPU call
9. Kept virtual sites, Coulomb k-space, and dipolar warning as CPU calls

- [ ] **Step 2: Verify compilation**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc) 2>&1 | tail -30`

Fix any compilation errors. Common issues:
- Missing includes for `VerletCriterion`, `propagation->thermostat`
- `Utils::Vector3d{}` may need to be `Utils::Vector3d{0., 0., 0.}` for some compilers
- The `propagation->thermostat` access path may differ — check `system/System.hpp` for the exact member name

- [ ] **Step 3: Commit**

```bash
git add src/core/pressure.cpp
git commit -m "feat: integrate Cabana pressure calculation into calculate_pressure()"
```

---

### Task 7: Run existing pressure tests

**Files:**
- Test: `testsuite/python/pressure.py`
- Test: `testsuite/python/interactions_bonded.py` (bonded pressure)
- Test: `testsuite/python/interactions_non-bonded.py` (non-bonded pressure)
- Test: `testsuite/python/stress_tensor_covariance.py` (rotation covariance)

- [ ] **Step 1: Run the core pressure test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/pressure.py`

Expected: All tests pass. The kinetic, bonded, non-bonded inter/intra, and total pressure tensors must match analytical results within tolerance.

- [ ] **Step 2: Run the bonded interaction pressure test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/interactions_bonded.py -v 2>&1 | grep -i pressure`

Expected: Bonded pressure assertions pass.

- [ ] **Step 3: Run the non-bonded interaction pressure test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/interactions_non-bonded.py -v 2>&1 | grep -i pressure`

Expected: Non-bonded pressure tensor assertions pass.

- [ ] **Step 4: Run the stress tensor covariance test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/stress_tensor_covariance.py`

Expected: All tests pass.

- [ ] **Step 5: Run the virtual sites pressure test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/virtual_sites_relative.py VirtualSitesRelTest.test_zz_pressure_tensor`

Expected: Pass.

- [ ] **Step 6: If any tests fail, debug and fix**

Common failure modes:
- Kinetic pressure: virtual particle filtering might not work via `has_exclusion` — check if `aosoa.flags` correctly marks virtuals
- Bonded pressure: `calc_bond_pair_force` might differ from `calc_bonded_virial_pressure_tensor` in the Coulomb-bond handling
- Non-bonded pressure: Thole or Gay-Berne stress contribution might be missing or double-counted
- DPD stress: sign convention (the CPU path uses `-tensor_product(d, f_dpd)`)

- [ ] **Step 7: Commit any fixes**

```bash
git add -u
git commit -m "fix: resolve pressure tensor Kokkos test failures"
```

---

### Task 8: Run MPI parallel tests

**Files:**
- Test: `testsuite/python/pressure.py` (with MPI)

- [ ] **Step 1: Run pressure test with MPI + OpenMP**

Run: `cd /tikhome/weeber/es/build && OMP_NUM_THREADS=4 mpirun -n 4 ./pypresso ../testsuite/python/pressure.py`

Expected: All tests pass. This verifies that thread-local accumulation and MPI reduction work correctly.

- [ ] **Step 2: Run bonded interaction test with MPI**

Run: `cd /tikhome/weeber/es/build && OMP_NUM_THREADS=2 mpirun -n 2 ./pypresso ../testsuite/python/interactions_bonded.py -v 2>&1 | grep -i pressure`

Expected: Pass.

- [ ] **Step 3: If any tests fail, debug and fix**

Common MPI issues:
- Thread-local accumulation with `omp_get_thread_num()` is OpenMP-specific — check that Kokkos is using the OpenMP backend
- Ghost particle handling: the AoSoA includes ghost particles, but kinetic pressure should only sum over local particles. The `KineticPressureKernel` currently iterates over all particles in the AoSoA including ghosts. Check if this matches the original `local_parts` loop which only iterates local particles.

- [ ] **Step 4: Commit any fixes**

```bash
git add -u
git commit -m "fix: resolve pressure tensor Kokkos MPI test failures"
```

---

### Task 9: Run DPD-specific pressure test

**Files:**
- Test: `testsuite/python/lees_edwards.py` (contains DPD pressure tensor test)

- [ ] **Step 1: Run the DPD pressure test**

Run: `cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/lees_edwards.py LeesEdwardsDPDTest -v 2>&1 | tail -30`

Expected: DPD pressure tensor test passes. The DPD dissipative stress computed inside `PressureKernel` must match the previous CPU-only result.

- [ ] **Step 2: If test fails, debug DPD stress**

The DPD stress calculation in the kernel uses zero-noise `dpd_pair_force(params, v21, dist, Utils::Vector3d{})` which calls the inner overload. Verify:
- `ia_params.dpd.radial` and `ia_params.dpd.trans` are correctly accessed from AoSoA type info
- The projection operator `P = tensor_product(d/dist2, d)` is numerically stable for small distances
- The sign convention matches: CPU uses `-Utils::flatten(stress)` in `dpd_pressure_local()`

- [ ] **Step 3: Commit any fixes**

```bash
git add -u
git commit -m "fix: resolve DPD pressure tensor Kokkos test failures"
```

---

### Task 10: Format and final verification

**Files:**
- All modified/created files

- [ ] **Step 1: Format C++ files**

Run:
```bash
cd /tikhome/weeber/es && maintainer/format/clang-format.sh -i src/core/pressure_cabana.hpp src/core/bond_pressure_kokkos.hpp src/core/pressure.cpp
```

If the formatter complains about missing packages, stop and ask the user to load the environment.

- [ ] **Step 2: Run full build**

Run: `cd /tikhome/weeber/es/build && make -j$(nproc)`

Expected: Clean build with no warnings.

- [ ] **Step 3: Run all pressure-related tests one more time**

Run:
```bash
cd /tikhome/weeber/es/build && ./pypresso ../testsuite/python/pressure.py && ./pypresso ../testsuite/python/stress_tensor_covariance.py && ./pypresso ../testsuite/python/virtual_sites_relative.py VirtualSitesRelTest.test_zz_pressure_tensor
```

Expected: All pass.

- [ ] **Step 4: Final commit**

```bash
git add -u
git commit -m "style: format pressure tensor Kokkos files"
```
