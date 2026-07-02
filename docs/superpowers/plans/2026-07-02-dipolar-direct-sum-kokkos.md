# Kokkos + kokkos-simd Dipolar Direct Sum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the three CPU dipolar direct-sum kernels (forces, energy, dipole-field) in `src/core/magnetostatics/dipolar_direct_sum.cpp` to Kokkos `parallel_for`/`parallel_reduce` with a `Kokkos::Experimental::native_simd<double>` inner particle loop, preserving results within test tolerances.

**Architecture:** Factor the pair math into `template <class T>` kernels shared by a scalar tail path (`T = double`) and a SIMD body (`T = native_simd<double>`), living in a new header. Rewrite the three methods to gather particle data (unchanged), repack it into SoA Kokkos views, and launch host (OpenMP) Kokkos kernels. The forces kernel keeps its two-phase MPI-overlap + Newton's-third-law `ScatterView` accumulation; the energy kernel gains the same two-phase MPI overlap; the field kernel is a plain `parallel_for`. The GPU (`ESPRESSO_CUDA`) path is untouched.

**Tech Stack:** C++20, Kokkos (≥4.6, OpenMP backend), kokkos-simd (`<Kokkos_SIMD.hpp>`), `Kokkos::Experimental::ScatterView`, Boost.MPI, `Utils::Vector<T,N>`, Boost.Test (unit tests), ESPResSo Python testsuite.

**Authoritative design reference:** `docs/superpowers/specs/2026-07-02-dipolar-direct-sum-kokkos-design.md`. Section numbers (§3.1, §3.6.1, …) below refer to it. When a task says "implement per §X", the spec section is the detailed contract; this plan shows the skeleton, novel code, interfaces, and exact verification.

## Global Constraints

- Kokkos execution space is the host: `Kokkos::DefaultExecutionSpace` (OpenMP in this build). CPU-only; do **not** touch the `ESPRESSO_CUDA` GPU path (`dipolar_direct_sum.hpp:43-86`, `dipolar_direct_sum_gpu*`).
- All new code guarded by `#ifdef ESPRESSO_DIPOLES`; `dipole_field` / field kernel additionally by `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING`. Torque code needs **no** `ESPRESSO_ROTATION` guard here (`DIPOLES implies ROTATION`, `src/config/features.def:46`).
- No CMake changes to `src/core/CMakeLists.txt` for the library (Kokkos already linked at line 79). The **unit test** does need a new `espresso_unit_test(...)` registration line.
- Preserve the exact reaction-torque relation: `fji.torque += vector_product(pf.f, rn) - pf.torque` (do **not** simplify to a recomputed pair). Preserve self-interaction exclusion (primary image of the self pair only).
- Preserve method signatures in `dipolar_direct_sum.hpp:81-88` verbatim; only bodies change.
- Build with `make -j8` (never `-j$(nproc)`). Use `git -C <path>` rather than `cd <path> && git`.
- SIMD load flag: `Kokkos::Experimental::simd_flag_default` (element-aligned). SIMD alias: `Kokkos::Experimental::simd<double>` (native ABI). NOTE: `native_simd` is deprecated and gated behind `KOKKOS_ENABLE_DEPRECATED_CODE_4`, which this build turns OFF — it will not compile; use `simd<double>`. Horizontal reduce: `Kokkos::Experimental::reduce(v, std::plus<>{})` returns the scalar `T`.

## File Structure

| File | Responsibility |
| --- | --- |
| `src/core/magnetostatics/dipolar_direct_sum_kernels.hpp` | **New.** Templated `pair_force`/`pair_potential`/`dipole_field`, `PairForce<T>`, `simd_double` alias, `Utils::operator*(simd, Vector<simd,N>)` overloads. Pure math, no MPI/Kokkos-launch code. Unit-testable in isolation. |
| `src/core/unit_tests/dipolar_direct_sum_kernels_test.cpp` | **New.** Boost.Test: analytic anchor + scalar/SIMD lane-equivalence for the three kernels. |
| `src/core/unit_tests/CMakeLists.txt` | **Modify.** Register the new unit test. |
| `src/core/magnetostatics/dipolar_direct_sum.cpp` | **Modify.** SoA repack + image-shift + distance helpers; rewrite the three methods. Remove old static `double` kernels, `for_each_image`, `image_sum`. |
| `src/core/magnetostatics/dipolar_direct_sum.hpp` | Unchanged. |

Work happens on branch `dds-kokkos` (already created). Commit after each task.

---

### Task 1: Templated pair kernels header + unit test

**Files:**
- Create: `src/core/magnetostatics/dipolar_direct_sum_kernels.hpp`
- Test: `src/core/unit_tests/dipolar_direct_sum_kernels_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Utils::Vector<T,3>`, `Utils::vector_product` (`src/utils/include/utils/Vector.hpp`); kokkos-simd.
- Produces (used by Tasks 3-5):
  - `using simd_double = Kokkos::Experimental::simd<double>;`
  - `template <class T> struct PairForce { Utils::Vector<T,3> f; Utils::Vector<T,3> torque; PairForce &operator+=(PairForce const &o); };`
  - `template <class T> PairForce<T> pair_force(Utils::Vector<T,3> const &d, Utils::Vector<T,3> const &m1, Utils::Vector<T,3> const &m2);`
  - `template <class T> T pair_potential(Utils::Vector<T,3> const &d, Utils::Vector<T,3> const &m1, Utils::Vector<T,3> const &m2);`
  - `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING template <class T> Utils::Vector<T,3> dipole_field(Utils::Vector<T,3> const &d, Utils::Vector<T,3> const &m1);`
  - `namespace Utils { template <std::size_t N> Vector<simd_double,N> operator*(simd_double const &, Vector<simd_double,N> const &); /* + symmetric */ }`

- [ ] **Step 1: Write the failing unit test**

Create `src/core/unit_tests/dipolar_direct_sum_kernels_test.cpp`:

```cpp
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

#define BOOST_TEST_MODULE "dipolar direct sum kernels"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#include <cstddef>

// Kokkos must be initialized before any simd/view use.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};
BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {
// Build a native_simd<double> whose lane l holds vals[l].
simd_double make_simd(double const *vals) {
  simd_double s;
  s.copy_from(vals, Kokkos::Experimental::simd_flag_default);
  return s;
}
Utils::Vector<simd_double, 3> broadcast(Utils::Vector3d const &v) {
  return {simd_double(v[0]), simd_double(v[1]), simd_double(v[2])};
}
} // namespace

BOOST_AUTO_TEST_SUITE(suite)

// Anchor: two identical dipoles m=(1,0,0) separated by d=(1,0,0).
// Analytically f=(-6,0,0), torque=0 (parallel moments => no torque).
BOOST_AUTO_TEST_CASE(pair_force_analytic) {
  Utils::Vector3d const d{1., 0., 0.};
  Utils::Vector3d const m{1., 0., 0.};
  auto const pf = pair_force<double>(d, m, m);
  BOOST_CHECK_CLOSE(pf.f[0], -6., 1e-10);
  BOOST_CHECK_SMALL(pf.f[1], 1e-12);
  BOOST_CHECK_SMALL(pf.f[2], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[0], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[1], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[2], 1e-12);
}

// Core SIMD-correctness property: each lane of the simd kernel equals the
// scalar kernel evaluated on that lane's inputs, for distinct per-lane data.
BOOST_AUTO_TEST_CASE(pair_force_simd_matches_scalar) {
  constexpr std::size_t w = simd_double::size();
  // Distinct, non-degenerate inputs per lane.
  std::vector<Utils::Vector3d> ds(w), m1s(w), m2s(w);
  for (std::size_t l = 0; l < w; ++l) {
    auto const x = 1. + 0.37 * double(l);
    ds[l]  = {x, 0.5 - 0.1 * double(l), 0.2 + 0.05 * double(l)};
    m1s[l] = {0.3 + 0.2 * double(l), 1.1, -0.4};
    m2s[l] = {-0.7, 0.9 - 0.05 * double(l), 0.6};
  }
  // Pack per-component lane buffers into simd inputs.
  auto pack = [&](std::vector<Utils::Vector3d> const &v, int c) {
    std::vector<double> buf(w);
    for (std::size_t l = 0; l < w; ++l) buf[l] = v[l][c];
    return make_simd(buf.data());
  };
  Utils::Vector<simd_double, 3> const d_s{pack(ds,0), pack(ds,1), pack(ds,2)};
  Utils::Vector<simd_double, 3> const m1_s{pack(m1s,0), pack(m1s,1), pack(m1s,2)};
  Utils::Vector<simd_double, 3> const m2_s{pack(m2s,0), pack(m2s,1), pack(m2s,2)};

  auto const pf_s = pair_force<simd_double>(d_s, m1_s, m2_s);

  for (std::size_t l = 0; l < w; ++l) {
    auto const pf = pair_force<double>(ds[l], m1s[l], m2s[l]);
    for (int c = 0; c < 3; ++c) {
      BOOST_CHECK_CLOSE(pf_s.f[c][l], pf.f[c], 1e-10);
      BOOST_CHECK_CLOSE(pf_s.torque[c][l], pf.torque[c], 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(pair_potential_simd_matches_scalar) {
  constexpr std::size_t w = simd_double::size();
  std::vector<Utils::Vector3d> ds(w), m1s(w), m2s(w);
  for (std::size_t l = 0; l < w; ++l) {
    ds[l]  = {1. + 0.3 * double(l), 0.4, -0.2};
    m1s[l] = {0.5, -0.3 + 0.1 * double(l), 0.8};
    m2s[l] = {0.2, 0.7, -0.6 + 0.05 * double(l)};
  }
  auto pack = [&](std::vector<Utils::Vector3d> const &v, int c) {
    std::vector<double> buf(w);
    for (std::size_t l = 0; l < w; ++l) buf[l] = v[l][c];
    return make_simd(buf.data());
  };
  Utils::Vector<simd_double, 3> const d_s{pack(ds,0), pack(ds,1), pack(ds,2)};
  Utils::Vector<simd_double, 3> const m1_s{pack(m1s,0), pack(m1s,1), pack(m1s,2)};
  Utils::Vector<simd_double, 3> const m2_s{pack(m2s,0), pack(m2s,1), pack(m2s,2)};
  auto const u_s = pair_potential<simd_double>(d_s, m1_s, m2_s);
  for (std::size_t l = 0; l < w; ++l) {
    auto const u = pair_potential<double>(ds[l], m1s[l], m2s[l]);
    BOOST_CHECK_CLOSE(u_s[l], u, 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()
```

- [ ] **Step 2: Register the test in CMake**

In `src/core/unit_tests/CMakeLists.txt`, add near the other `Kokkos::kokkos` tests (e.g. after the `field_layout_helpers_test` line ~67):

```cmake
espresso_unit_test(SRC dipolar_direct_sum_kernels_test.cpp DEPENDS espresso::core
                   espresso::utils Kokkos::kokkos)
```

- [ ] **Step 3: Run the test to verify it fails (header missing)**

Configure a build dir if none exists, then build the test:

```bash
cmake -S . -B build -D ESPRESSO_BUILD_WITH_CUDA=OFF \
      -D ESPRESSO_TEST_NP=2 -D CMAKE_BUILD_TYPE=Release \
      -D ESPRESSO_BUILD_WITH_UNIT_TESTS=ON
cmake --build build --target dipolar_direct_sum_kernels_test -j8
```

Expected: **compile failure** — `magnetostatics/dipolar_direct_sum_kernels.hpp: No such file or directory`.

- [ ] **Step 4: Create the kernels header**

Create `src/core/magnetostatics/dipolar_direct_sum_kernels.hpp` (math copied verbatim from the current `dipolar_direct_sum.cpp:62-124`, templated, with ADL `sqrt`):

```cpp
/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifdef ESPRESSO_DIPOLES

#include <utils/Vector.hpp>

#include <Kokkos_SIMD.hpp>

#include <cmath>
#include <cstddef>

using simd_double = Kokkos::Experimental::simd<double>;

namespace Utils {
/** simd-scalar * Vector<simd> — the arithmetic-constrained scalar overload in
 *  Vector.hpp excludes class scalar types, so provide it here for ADL. */
template <std::size_t N>
Vector<simd_double, N> operator*(simd_double const &a,
                                 Vector<simd_double, N> const &b) {
  Vector<simd_double, N> out;
  for (std::size_t i = 0; i < N; ++i)
    out[i] = a * b[i];
  return out;
}
template <std::size_t N>
Vector<simd_double, N> operator*(Vector<simd_double, N> const &a,
                                 simd_double const &b) {
  return b * a;
}
} // namespace Utils

/** @brief Force and torque of one pair interaction, generic scalar type. */
template <class T> struct PairForce {
  Utils::Vector<T, 3> f{};
  Utils::Vector<T, 3> torque{};
  PairForce &operator+=(PairForce const &o) {
    f += o.f;
    torque += o.torque;
    return *this;
  }
};

/** @brief Pair force of two interacting dipoles (see dipolar_direct_sum.cpp). */
template <class T>
PairForce<T> pair_force(Utils::Vector<T, 3> const &d,
                        Utils::Vector<T, 3> const &m1,
                        Utils::Vector<T, 3> const &m2) {
  using std::sqrt; // ADL: std::sqrt for double, Kokkos simd sqrt for simd
  auto const pe2 = m1 * d;
  auto const pe3 = m2 * d;

  auto const r2 = d.norm2();
  auto const r = sqrt(r2);
  auto const r5 = r2 * r2 * r;
  auto const r7 = r5 * r2;

  auto const a = 3.0 * (m1 * m2) / r5;
  auto const b = -15.0 * pe2 * pe3 / r7;

  auto const f = (a + b) * d + 3.0 * (pe3 * m1 + pe2 * m2) / r5;
  auto const r3 = r2 * r;
  auto const t =
      -vector_product(m1, m2) / r3 + 3.0 * pe3 * vector_product(m1, d) / r5;

  return PairForce<T>{f, t};
}

/** @brief Pair potential for two interacting dipoles. */
template <class T>
T pair_potential(Utils::Vector<T, 3> const &d, Utils::Vector<T, 3> const &m1,
                 Utils::Vector<T, 3> const &m2) {
  using std::sqrt;
  auto const r2 = d * d;
  auto const r = sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;

  auto const pe1 = m1 * m2;
  auto const pe2 = m1 * d;
  auto const pe3 = m2 * d;

  return pe1 / r3 - 3.0 * pe2 * pe3 / r5;
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
/** @brief Dipole field contribution from a dipole @c m1 at distance @c d. */
template <class T>
Utils::Vector<T, 3> dipole_field(Utils::Vector<T, 3> const &d,
                                 Utils::Vector<T, 3> const &m1) {
  using std::sqrt;
  auto const r2 = d * d;
  auto const r = sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;
  auto const pe2 = m1 * d;

  return 3.0 * pe2 * d / r5 - m1 / r3;
}
#endif // ESPRESSO_DIPOLE_FIELD_TRACKING

#endif // ESPRESSO_DIPOLES
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build --target dipolar_direct_sum_kernels_test -j8
ctest --test-dir build -R dipolar_direct_sum_kernels --output-on-failure
```

Expected: build succeeds and **all cases PASS**. If `pair_force_simd_matches_scalar` fails only at machine epsilon, that is a real problem (the two paths must be bit-close); if it fails to compile on `simd_flag_default` or `native_simd`, fix the alias/flag spelling (spec §7.2) and note the correct spelling for later tasks.

- [ ] **Step 6: Commit**

```bash
git add src/core/magnetostatics/dipolar_direct_sum_kernels.hpp \
        src/core/unit_tests/dipolar_direct_sum_kernels_test.cpp \
        src/core/unit_tests/CMakeLists.txt
git commit -m "magnetostatics: add templated dipolar pair kernels (scalar + simd)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: SoA repack + image-shift + distance helpers

Add the file-static support machinery to `dipolar_direct_sum.cpp` needed by all three kernel rewrites. No behavior change yet (helpers unused until Tasks 3-5), so verification is: the translation unit still compiles and existing tests still pass.

**Files:**
- Modify: `src/core/magnetostatics/dipolar_direct_sum.cpp`

**Interfaces:**
- Consumes: Task 1 header (`simd_double`), `gather_particle_data` (existing, unchanged), `BoxGeometry`.
- Produces (used by Tasks 3-5):
  - `struct PosMomViews { Kokkos::View<double*[3], Kokkos::LayoutLeft, Kokkos::HostSpace> pos, m; };`
  - `static PosMomViews make_posmom_views(std::size_t n_total);`
  - `static void fill_posmom_views(PosMomViews &, std::vector<PosMom> const &, std::size_t begin, std::size_t end);`
  - `static std::vector<Utils::Vector3d> make_image_shifts(Utils::Vector3i const &ncut, Utils::Vector3d const &box_l);` — index 0 is the zero (primary) shift.

- [ ] **Step 1: Add includes and helpers**

At the top of `dipolar_direct_sum.cpp`, after the existing includes add:

```cpp
#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_SIMD.hpp>
```

Then add the helpers (near the top, after `PosMom` / `get_n_cut`). Implement per spec §3.3 and §3.4:

```cpp
struct PosMomViews {
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> pos;
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> m;
};

static PosMomViews make_posmom_views(std::size_t n_total) {
  return {decltype(PosMomViews::pos)("dds_pos", n_total),
          decltype(PosMomViews::m)("dds_m", n_total)};
}

static void fill_posmom_views(PosMomViews &views,
                              std::vector<PosMom> const &all_posmom,
                              std::size_t begin, std::size_t end) {
  for (auto i = begin; i < end; ++i) {
    for (int c = 0; c < 3; ++c) {
      views.pos(i, c) = all_posmom[i].pos[c];
      views.m(i, c) = all_posmom[i].m[c];
    }
  }
}

/** Real-space image shifts n .* box_l inside the |ncut| sphere; index 0 is the
 *  primary (zero) shift so self-interaction loops start at index 1. */
static std::vector<Utils::Vector3d>
make_image_shifts(Utils::Vector3i const &ncut, Utils::Vector3d const &box_l) {
  auto const ncut2 = ncut.norm2();
  std::vector<Utils::Vector3d> shifts;
  shifts.push_back({0., 0., 0.});
  for (int nx = -ncut[0]; nx <= ncut[0]; ++nx)
    for (int ny = -ncut[1]; ny <= ncut[1]; ++ny)
      for (int nz = -ncut[2]; nz <= ncut[2]; ++nz) {
        if (nx == 0 && ny == 0 && nz == 0)
          continue;
        if (nx * nx + ny * ny + nz * nz <= ncut2)
          shifts.push_back(
              {nx * box_l[0], ny * box_l[1], nz * box_l[2]});
      }
  return shifts;
}
```

> Note: `make_image_shifts` reproduces `for_each_image`'s sphere predicate (`dipolar_direct_sum.cpp:141-159`) but emits the zero shift first. The old `for_each_image` and `image_sum` are removed only in Tasks 3-5 once their callers are gone; leaving them temporarily is fine (they may warn as unused — acceptable mid-task, resolved by Task 5).

- [ ] **Step 2: Build the core library to verify it compiles**

```bash
cmake --build build --target espresso_core -j8
```

Expected: **compiles cleanly** (helpers may be reported unused — OK for now).

- [ ] **Step 3: Commit**

```bash
git add src/core/magnetostatics/dipolar_direct_sum.cpp
git commit -m "magnetostatics: add SoA repack and image-shift helpers for DDS

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Rewrite the forces kernel (Kokkos + SIMD, two-phase, ScatterView)

Replace `add_long_range_forces_cpu` (`dipolar_direct_sum.cpp:291-380`) with the two-phase Kokkos implementation. This is the largest task; implement **exactly** per spec §3.6.1.

**Files:**
- Modify: `src/core/magnetostatics/dipolar_direct_sum.cpp:291-380`

**Interfaces:**
- Consumes: Task 1 (`pair_force<T>`, `PairForce<T>`, `simd_double`), Task 2 (`PosMomViews`, `make_posmom_views`, `fill_posmom_views`, `make_image_shifts`), `gather_particle_data`, `get_n_cut`.
- Produces: no new external interface (method body only).

- [ ] **Step 1: Establish the regression baseline (test must pass on current code)**

```bash
cmake --build build --target pypresso espresso_pyscript -j8
ctest --test-dir build -R dipolar_direct_summation --output-on-failure
```

Expected: **PASS** on the current (pre-rewrite) code — this is the reference the rewrite must preserve. Record which cases run at this NP.

- [ ] **Step 2: Rewrite `add_long_range_forces_cpu` per spec §3.6.1**

Structure (see spec §3.6.1 steps 1-9 for the full contract). Skeleton with the required accumulator types, the SIMD body, the **per-lane scatter**, and the exact torque relation:

```cpp
void DipolarDirectSum::add_long_range_forces_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset] =
      gather_particle_data(box_geo, particles);

  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();
  constexpr std::size_t w = simd_double::size();

  // SoA views; fill local slice only (MPI overlap).
  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, offset, offset + n_local);

  using execution_space = Kokkos::DefaultExecutionSpace;
  using ForceView =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ScatterForce =
      Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight>;
  ForceView local_force("dds_force", n_local);
  ForceView local_torque("dds_torque", n_local);
  ScatterForce scatter_force(local_force);
  ScatterForce scatter_torque(local_torque);

  // Phase A: local pairs. Each i owns its own force accumulation; partner-j
  // contributions go through the ScatterView (per-lane scatter).
  Kokkos::parallel_for(
      "dds_local_pairs",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        auto const gi = offset + i;
        Utils::Vector3d const pos_i{views.pos(gi, 0), views.pos(gi, 1),
                                    views.pos(gi, 2)};
        Utils::Vector3d const m_i{views.m(gi, 0), views.m(gi, 1),
                                  views.m(gi, 2)};
        PairForce<double> fi{};

        // (a) self-images (shifts[1..], primary excluded)
        for (std::size_t s = 1; s < shifts.size(); ++s)
          fi += pair_force<double>(shifts[s], m_i, m_i);

        auto force_access = scatter_force.access();
        auto torque_access = scatter_torque.access();
        auto const m_i_s = broadcast_simd(m_i); // Utils::Vector<simd_double,3>

        // (b) SIMD body over j in (gi, offset + n_local)
        auto j = gi + 1;
        for (; j + w <= offset + n_local; j += w) {
          auto const m_j = load_simd_moment(views, j); // §3.5 helper
          auto const d0 =
              primary_distance_simd(pos_i, views, j, with_replicas, box_geo);
          PairForce<simd_double> fij{}, fji{};
          for (auto const &shift : shifts) {
            Utils::Vector<simd_double, 3> const rn{
                d0[0] + shift[0], d0[1] + shift[1], d0[2] + shift[2]};
            auto const pf = pair_force<simd_double>(rn, m_i_s, m_j);
            fij += pf;
            fji.f -= pf.f;
            /* Conservation of angular momentum mandates that
             * 0 = t_i + r_ij x F_ij + t_j */
            fji.torque += vector_product(pf.f, rn) - pf.torque;
          }
          // i-side: horizontal reduce into fi
          for (int c = 0; c < 3; ++c) {
            fi.f[c] += Kokkos::Experimental::reduce(fij.f[c], std::plus<>{});
            fi.torque[c] +=
                Kokkos::Experimental::reduce(fij.torque[c], std::plus<>{});
          }
          // j-side: per-lane scatter
          for (std::size_t l = 0; l < w; ++l) {
            auto const jl = (j + l) - offset;
            for (int c = 0; c < 3; ++c) {
              force_access(jl, c) += fji.f[c][l];
              torque_access(jl, c) += fji.torque[c][l];
            }
          }
        }
        // (c) scalar tail over remaining j in (.., offset + n_local)
        for (; j < offset + n_local; ++j) {
          Utils::Vector3d const pos_j{views.pos(j, 0), views.pos(j, 1),
                                      views.pos(j, 2)};
          Utils::Vector3d const m_j{views.m(j, 0), views.m(j, 1),
                                    views.m(j, 2)};
          auto const d0 = with_replicas
                              ? (pos_i - pos_j)
                              : box_geo.get_mi_vector(pos_i, pos_j);
          for (auto const &shift : shifts) {
            auto const rn = d0 + shift;
            auto const pf = pair_force<double>(rn, m_i, m_j);
            fi.f += pf.f;
            fi.torque += pf.torque;
            auto const jl = j - offset;
            for (int c = 0; c < 3; ++c) {
              force_access(jl, c) -= pf.f[c];
              torque_access(jl, c) +=
                  (vector_product(pf.f, rn) - pf.torque)[c];
            }
          }
        }
        // (d) write i's own total directly (unique owner, no race)
        local_particles[i]->force() += prefactor * fi.f;
        local_particles[i]->torque() += prefactor * fi.torque;
      });
  Kokkos::fence();

  // Wait for remote data, fill remote slices.
  boost::mpi::wait_all(reqs.begin(), reqs.end());
  fill_posmom_views(views, all_posmom, 0, offset);
  fill_posmom_views(views, all_posmom, offset + n_local, n_total);

  // Phase B: remote pairs (red [0,offset) + black [offset+n_local, n_total)),
  // visit-twice, no scatter — accumulate only i.
  Kokkos::parallel_for(
      "dds_remote_pairs",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        // for each remote range: SIMD body + scalar tail summing into a
        // PairForce<double> fi (only i is written), see spec §3.6.1 step 7.
        // ... (mirror of the (b)/(c) blocks above, but no scatter, and looping
        //     the two remote ranges [0, offset) and [offset+n_local, n_total))
        // local_particles[i]->force()  += prefactor * fi.f;
        // local_particles[i]->torque() += prefactor * fi.torque;
      });
  Kokkos::fence();

  // Reduce the Newton's-3rd-law contributions and add to particles.
  Kokkos::Experimental::contribute(local_force, scatter_force);
  Kokkos::Experimental::contribute(local_torque, scatter_torque);
  Kokkos::parallel_for(
      "dds_reduction",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        local_particles[i]->force() +=
            prefactor * Utils::Vector3d{local_force(i, 0), local_force(i, 1),
                                        local_force(i, 2)};
        local_particles[i]->torque() +=
            prefactor * Utils::Vector3d{local_torque(i, 0), local_torque(i, 1),
                                        local_torque(i, 2)};
      });
  Kokkos::fence();

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (not m_is_gpu) {
    dipole_field_at_part_cpu();
  }
#endif
}
```

Also add the two SIMD load/distance helpers (spec §3.5) as file-static functions before the method, plus a `broadcast_simd(Utils::Vector3d)` helper:

```cpp
static Utils::Vector<simd_double, 3> broadcast_simd(Utils::Vector3d const &v) {
  return {simd_double(v[0]), simd_double(v[1]), simd_double(v[2])};
}

static Utils::Vector<simd_double, 3>
load_simd_moment(PosMomViews const &views, std::size_t j) {
  Utils::Vector<simd_double, 3> m;
  for (int c = 0; c < 3; ++c)
    m[c].copy_from(&views.m(j, c), Kokkos::Experimental::simd_flag_default);
  return m;
}

static Utils::Vector<simd_double, 3>
primary_distance_simd(Utils::Vector3d const &pos_i, PosMomViews const &views,
                      std::size_t j, bool with_replicas,
                      BoxGeometry const &box_geo) {
  constexpr std::size_t w = simd_double::size();
  Utils::Vector<simd_double, 3> d;
  if (with_replicas) {
    for (int c = 0; c < 3; ++c) {
      simd_double pj;
      pj.copy_from(&views.pos(j, c), Kokkos::Experimental::simd_flag_default);
      d[c] = simd_double(pos_i[c]) - pj;
    }
  } else {
    // per-lane minimum image (preserves Lees-Edwards semantics)
    double buf[3][w];
    for (std::size_t l = 0; l < w; ++l) {
      Utils::Vector3d const pos_j{views.pos(j + l, 0), views.pos(j + l, 1),
                                  views.pos(j + l, 2)};
      auto const mi = box_geo.get_mi_vector(pos_i, pos_j);
      for (int c = 0; c < 3; ++c)
        buf[c][l] = mi[c];
    }
    for (int c = 0; c < 3; ++c)
      d[c].copy_from(&buf[c][0], Kokkos::Experimental::simd_flag_default);
  }
  return d;
}
```

Fully implement the Phase B body (the `...` comment) mirroring the (b) SIMD-body and (c) scalar-tail blocks for each of the two remote ranges, **without** the scatter (only `fi` accumulates), per spec §3.6.1 step 7.

- [ ] **Step 3: Build**

```bash
cmake --build build --target espresso_core pypresso -j8
```

Expected: compiles cleanly. Fix any `native_simd`/`simd_flag_default`/`reduce` spelling issues using the resolution found in Task 1.

- [ ] **Step 4: Run the forces regression tests**

```bash
ctest --test-dir build -R dipolar_direct_summation --output-on-failure
```

Expected: **PASS**, including `test_dds_cpu` (1E-12 vs stored `.npy`) and, at NP≥2, `test_inner_loop_consistency_cpu` (which directly exercises the Phase A/B split). If `test_dds_cpu` fails at ~1e-13, that is the reordering effect flagged in spec §7.1 — **stop and escalate** (do not regenerate references).

- [ ] **Step 5: Run at 1 and 2 threads to isolate threading**

```bash
OMP_NUM_THREADS=1 ctest --test-dir build -R dipolar_direct_summation --output-on-failure
OMP_NUM_THREADS=4 ctest --test-dir build -R dipolar_direct_summation --output-on-failure
```

Expected: PASS in both. A pass at 1 thread but fail at 4 indicates a ScatterView/aliasing bug in Phase A — debug there.

- [ ] **Step 6: Commit**

```bash
git add src/core/magnetostatics/dipolar_direct_sum.cpp
git commit -m "magnetostatics: port DDS forces to Kokkos + kokkos-simd

Two-phase MPI-overlap preserved; Newton's-third-law partner writes via
ScatterView with per-lane scatter; inner particle loop vectorized with
native_simd<double>.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Rewrite the energy kernel (Kokkos parallel_reduce, two-phase)

Replace `long_range_energy_cpu` (`dipolar_direct_sum.cpp:387-416`) per spec §3.6.2, adding the two-phase MPI overlap (user decision).

**Files:**
- Modify: `src/core/magnetostatics/dipolar_direct_sum.cpp:387-416`

**Interfaces:**
- Consumes: Task 1 (`pair_potential<T>`, `simd_double`), Task 2/3 helpers (`PosMomViews`, `make_posmom_views`, `fill_posmom_views`, `make_image_shifts`, `broadcast_simd`, `load_simd_moment`, `primary_distance_simd`), `gather_particle_data`, `get_n_cut`.
- Produces: none (method body only).

- [ ] **Step 1: Rewrite `long_range_energy_cpu` per spec §3.6.2**

Skeleton (two `parallel_reduce` phases; local-upper `[gi, offset+n_local)` before wait, remote-black `[offset+n_local, n_total)` after; **no red range**):

```cpp
double DipolarDirectSum::long_range_energy_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset] =
      gather_particle_data(box_geo, particles);

  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_geo.length());
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();
  constexpr std::size_t w = simd_double::size();

  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, offset, offset + n_local);

  using execution_space = Kokkos::DefaultExecutionSpace;

  // Phase A: local-upper triangular sum [gi, offset + n_local), self incl.
  double uA = 0.;
  Kokkos::parallel_reduce(
      "dds_energy_local",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i, double &u_local) {
        auto const gi = offset + i;
        // pos_i, m_i as in forces; self-images over shifts[1..] with
        // pair_potential<double>; SIMD body over j in (gi, offset+n_local)
        // with a simd_double acc reduced via Kokkos::Experimental::reduce;
        // scalar tail. Accumulate into u_local. (spec §3.6.2 step 3)
      },
      uA);

  boost::mpi::wait_all(reqs.begin(), reqs.end());
  fill_posmom_views(views, all_posmom, offset + n_local, n_total);

  // Phase B: remote-black [offset + n_local, n_total); no self term.
  double uB = 0.;
  Kokkos::parallel_reduce(
      "dds_energy_remote",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i, double &u_local) {
        // SIMD body + scalar tail over j in [offset+n_local, n_total)
        // with pair_potential; accumulate into u_local. (spec §3.6.2 step 6)
      },
      uB);
  Kokkos::fence();

  return prefactor * (uA + uB);
}
```

Fully implement the two lambda bodies. The self-image term (Phase A only): `for (s=1; s<shifts.size(); ++s) u_local += pair_potential<double>(shifts[s], m_i, m_i);`. SIMD body accumulator: `simd_double acc{0.}; ... acc += pair_potential<simd_double>(rn, m_i_s, m_j);` then `u_local += Kokkos::Experimental::reduce(acc, std::plus<>{});`.

- [ ] **Step 2: Build**

```bash
cmake --build build --target espresso_core pypresso -j8
```

Expected: compiles cleanly.

- [ ] **Step 3: Run energy + full DDS tests (1 and 2 ranks)**

```bash
ctest --test-dir build -R dipolar_direct_summation --output-on-failure
```

Expected: **PASS**. `test_dds_cpu` checks energy at 1E-12 against the stored reference; the two-phase reassociation may perturb it at ~1e-14 (spec §7.1/§7.8) — escalate if it fails marginally rather than regenerating references.

- [ ] **Step 4: Commit**

```bash
git add src/core/magnetostatics/dipolar_direct_sum.cpp
git commit -m "magnetostatics: port DDS energy to Kokkos + kokkos-simd

parallel_reduce over the primed triangular sum, vectorized inner loop,
with two-phase MPI latency-hiding (local block computed before wait).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Rewrite the dipole-field kernel + remove dead code

Replace `dipole_field_at_part_cpu` (`dipolar_direct_sum.cpp:428-457`) per spec §3.6.3, then remove the now-unused `image_sum` and `for_each_image` and the old `double`-only pair functions (superseded by the header).

**Files:**
- Modify: `src/core/magnetostatics/dipolar_direct_sum.cpp` (field kernel + dead-code removal)

**Interfaces:**
- Consumes: Task 1 (`dipole_field<T>`, `simd_double`), Task 2/3 helpers.
- Produces: none.

- [ ] **Step 1: Rewrite `dipole_field_at_part_cpu` per spec §3.6.3**

`parallel_for` over local `i` (waits first, full fill), SIMD sweep over **all** `j` with the self-index `gi` split out (ranges `[0, gi)` and `[gi+1, n_total)`) plus a scalar self-image term over `shifts[1..]`, then **assign** `dip_fld() = prefactor * u`:

```cpp
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void DipolarDirectSum::dipole_field_at_part_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset] =
      gather_particle_data(box_geo, particles);

  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_geo.length());
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();

  boost::mpi::wait_all(reqs.begin(), reqs.end());
  auto views = make_posmom_views(n_total);
  fill_posmom_views(views, all_posmom, 0, n_total);

  using execution_space = Kokkos::DefaultExecutionSpace;
  Kokkos::parallel_for(
      "dds_dipole_field",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_local),
      [=](std::size_t const i) {
        auto const gi = offset + i;
        Utils::Vector3d u{};
        // SIMD sweep [0, gi) and [gi+1, n_total) with dipole_field<simd_double>,
        // horizontal-reduce per component + scalar tails; plus scalar
        // self-image term over shifts[1..] with dipole_field<double>.
        // (spec §3.6.3)
        local_particles[i]->dip_fld() = prefactor * u;
      });
  Kokkos::fence();
}
#endif
```

Fully implement the sweep bodies.

- [ ] **Step 2: Remove dead code**

Delete from `dipolar_direct_sum.cpp`: the old static `pair_force`, `pair_potential`, `dipole_field` (lines ~62-124), `for_each_image` (~141-159), and `image_sum` (~198-221). Remove now-unused includes `<utils/cartesian_product.hpp>` and `<ranges>` if nothing else uses them (grep the file first). Keep `PosMom`, `gather_particle_data`, `get_n_cut`.

- [ ] **Step 3: Build with field tracking enabled (maxset)**

The default config lacks `DIPOLE_FIELD_TRACKING`; use maxset to compile and exercise the field kernel:

```bash
cmake -S . -B build-maxset -D ESPRESSO_BUILD_WITH_CUDA=OFF \
      -D CMAKE_BUILD_TYPE=Release -D ESPRESSO_BUILD_WITH_UNIT_TESTS=ON \
      -D ESPRESSO_TEST_NP=2 \
      -D ESPRESSO_CONFIG_FILE=$PWD/maintainer/configs/maxset.hpp
cmake --build build-maxset --target espresso_core pypresso -j8
```

Expected: compiles cleanly with all three kernels active.

- [ ] **Step 4: Run DDS tests on maxset build (field tracking active)**

```bash
ctest --test-dir build-maxset -R dipolar_direct_summation --output-on-failure
```

Expected: **PASS** (the field path is exercised via `dipolar_direct_summation.py` when `DIPOLE_FIELD_TRACKING` is compiled; check the test file for the field-specific asserts).

- [ ] **Step 5: Also rebuild + retest the default (non-maxset) build**

```bash
cmake --build build --target espresso_core pypresso -j8
ctest --test-dir build -R dipolar_direct_summation --output-on-failure
```

Expected: PASS (confirms the `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING` guards keep the field code out cleanly when the feature is off).

- [ ] **Step 6: Commit**

```bash
git add src/core/magnetostatics/dipolar_direct_sum.cpp
git commit -m "magnetostatics: port DDS dipole-field to Kokkos, drop scalar kernels

Vectorized dipole-field sweep; remove image_sum/for_each_image and the
double-only pair functions now superseded by dipolar_direct_sum_kernels.hpp.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: Full verification sweep

Confirm the whole feature across rank/thread/config combinations and run the broader dipolar test group to catch regressions in dependent code.

**Files:** none (verification only).

- [ ] **Step 1: Run the full dipolar test group (default build)**

```bash
ctest --test-dir build -R "dipol" --output-on-failure
```

Expected: all selected tests PASS. `dawaanr-and-dds-gpu.py` is skipped (GPU-gated) in this CPU build — confirm it reports as skipped, not failed.

- [ ] **Step 2: Multi-rank + multi-thread matrix**

```bash
for np in 1 2; do for th in 1 4; do
  echo "=== NP=$np OMP=$th ==="; \
  OMP_NUM_THREADS=$th ctest --test-dir build -R dipolar_direct_summation \
    --output-on-failure; \
done; done
```

Expected: PASS in all four combinations. (NP is fixed by the test's `MAX_NUM_PROC 2`; ctest runs the registered rank counts.)

- [ ] **Step 3: maxset build full dipolar group**

```bash
ctest --test-dir build-maxset -R "dipol" --output-on-failure
```

Expected: PASS (with field tracking active).

- [ ] **Step 4: Warning check on touched files**

Rebuild the two touched translation units verbose and confirm no new warnings:

```bash
cmake --build build --target espresso_core -j8 2>&1 | \
  grep -Ei "dipolar_direct_sum|warning" || echo "no warnings in DDS files"
```

Expected: no new warnings referencing `dipolar_direct_sum*`.

- [ ] **Step 5: Final review + finishing skill**

Use superpowers:requesting-code-review on the branch diff, then superpowers:finishing-a-development-branch to decide integration (PR vs merge). Report the full test matrix results.

---

## Self-Review

**Spec coverage:**
- §3.1 templated pair kernels → Task 1 (header + tests). ✓
- §3.2 `Utils::operator*` simd overloads → Task 1 (in header). ✓
- §3.3 SoA repack (`PosMomViews`, `make/fill_posmom_views`) → Task 2. ✓
- §3.4 image-shift precompute (`make_image_shifts`) → Task 2. ✓
- §3.5 distance helpers (`primary_distance_simd`, scalar counterpart) → Task 3 (defined there, reused by 4/5). ✓
- §3.6.1 forces (two-phase, ScatterView, per-lane scatter, torque relation) → Task 3. ✓
- §3.6.2 energy (two-phase reduce, MPI overlap) → Task 4. ✓
- §3.6.3 field (assign, self-split) + dead-code removal → Task 5. ✓
- §5 correctness (tail, ADL sqrt, ScatterView, MPI order, self-exclusion, guards) → threaded through Tasks 1-5 steps + verification. ✓
- §6 build/verification (make -j8, python tests, maxset for field, 1/2 rank, threads) → Tasks 3-6. ✓
- §7 risks (tolerance escalation, simd spelling, layout fallback) → called out in Task 1 Step 5, Task 3 Step 4, and thread test Task 3 Step 5. ✓

**Placeholder scan:** The `...` blocks in Tasks 3-5 skeletons are explicitly delegated to named spec sections with the novel/hard code (per-lane scatter, SIMD distance, torque relation, horizontal reduce) shown inline; each such step says "fully implement the bodies" and cites the spec contract. No vague "add error handling"-type placeholders.

**Type consistency:** `PairForce<T>` / `pair_force<T>` / `pair_potential<T>` / `dipole_field<T>` / `simd_double` / `PosMomViews` / `make_posmom_views` / `fill_posmom_views` / `make_image_shifts` / `broadcast_simd` / `load_simd_moment` / `primary_distance_simd` — names and signatures are consistent between the Interfaces blocks and the code across Tasks 1-5. `broadcast_simd` (used in Task 3) is defined in Task 3; the unit test in Task 1 has its own local `broadcast`/`make_simd` helpers (test-only, no collision).
