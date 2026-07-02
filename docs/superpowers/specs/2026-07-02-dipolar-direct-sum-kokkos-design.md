# Design: Kokkos + kokkos-simd port of the CPU dipolar direct sum

- **Date:** 2026-07-02
- **Branch / worktree:** `python` @ `/tikhome/weeber/es/.claude/worktrees/es-dds`
- **Status:** approved design, ready for implementation
- **Scope:** `src/core/magnetostatics/dipolar_direct_sum.cpp` — all three CPU kernels

## 1. Overview / goal

Convert the CPU dipolar direct sum (`DipolarDirectSum::add_long_range_forces_cpu`,
`DipolarDirectSum::long_range_energy_cpu`, and — behind
`ESPRESSO_DIPOLE_FIELD_TRACKING` — `DipolarDirectSum::dipole_field_at_part_cpu`) from
serial N-square loops to Kokkos `parallel_for`/`parallel_reduce` over the local
particle index `i`, with the inner `j`-particle loop vectorized using
`Kokkos::Experimental::native_simd<double>` (kokkos-simd). The pair math is factored
into `template <class T>` kernels shared between the scalar tail path
(`T = double`) and the SIMD body (`T = native_simd<double>`), operating on
`Utils::Vector<T, 3>`. Both existing optimizations of the forces kernel are
preserved: the MPI latency-hiding two-phase structure, and the Newton's-third-law
"visit local pairs once" trick (which under a parallel `i`-loop requires a
`Kokkos::Experimental::ScatterView` for the partner-`j` writes). The energy kernel,
which is serial-in-MPI today, additionally **gains** the same two-phase MPI
latency-hiding (compute the local triangular sum while the remote data is in
flight) as part of this change. The GPU (`ESPRESSO_CUDA`) path is untouched.
Verification is via the existing test suite.

## 2. Background: current implementation

All code below is in `src/core/magnetostatics/dipolar_direct_sum.cpp` unless noted.
The entire translation unit is guarded by `#ifdef ESPRESSO_DIPOLES` (line 24).
`DIPOLES implies ROTATION` (`src/config/features.def:46`), so `p.torque()` and
`ParticleForce::torque` (guarded by `ESPRESSO_ROTATION` in
`src/core/Particle.hpp:330-369`) are always available inside this file; the current
code already relies on this.

### 2.1 Pair kernels (currently `double`-only, file-static)

- `pair_force(d, m1, m2)` (lines 62-81): dipole-dipole force and torque on particle
  1; returns `ParticleForce{f, t}`. Uses `d.norm2()`, `std::sqrt(r2)` (qualified —
  relevant later), dot products via `Utils::Vector::operator*`, and
  `Utils::vector_product`.
- `pair_potential(d, m1, m2)` (lines 92-104): interaction energy; uses unqualified
  `sqrt(r2)`.
- `dipole_field(d, m1)` (lines 116-124, `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING`):
  dipole field contribution; unqualified `sqrt(r2)`.

All three are branch-free straight-line math.

### 2.2 Support machinery

- `for_each_image(ncut, f)` (lines 141-159): calls `f(nx, ny, nz)` for every integer
  triple in `[-ncut, ncut]^3` with `nx² + ny² + nz² <= ncut²` (sphere test), via
  `Utils::cartesian_product` over `std::views::iota` ranges.
- `PosMom` (lines 164-171): AoS record `{ Utils::Vector3d pos; Utils::Vector3d m; }`
  with Boost.MPI serialization.
- `image_sum(begin, end, it, with_replicas, ncut, box_geo, init, f)` (lines
  198-221): the "primed" pair sum for one particle `it` over a `j`-range. For each
  `jt`, the primary distance is `it->pos - jt->pos` if `with_replicas`, else
  `box_geo.get_mi_vector(...)` (minimum image). Then for each image shift the pair
  kernel is accumulated, skipping only the primary image when `it == jt`
  (self-interaction exclusion; the particle *does* interact with its own periodic
  images).
- `gather_particle_data(box_geo, particles)` (lines 223-261): collects local
  particles with `dipm() != 0` into `local_particles` (pointers, in the same order)
  and `local_posmom` (folded positions + `calc_dip()`), `all_gather`s the per-rank
  counts, then starts a **non-blocking** `Utils::Mpi::iall_gatherv` into
  `all_posmom` and returns `(local_particles, all_posmom, reqs, offset)`. With one
  rank, `all_posmom` is just swapped from `local_posmom` and `reqs` is empty. The
  local slice `[offset, offset + local_particles.size())` of `all_posmom` is valid
  immediately; the rest only after `wait_all`.
- `get_n_cut(box_geo, n_replicas)` (lines 263-267): per-axis image cutoff,
  `n_replicas` masked by periodicity flags; `with_replicas = (ncut.norm2() > 0)`.

### 2.3 The three kernels

- **Forces** (`add_long_range_forces_cpu`, lines 291-380). Two-phase structure
  documented at lines 269-290:
  1. start async gather;
  2. **Phase A** (lines 313-347): loop over local `i`; (a) interaction with own
     periodic images (`image_sum(it, next(it), it, ...)`, lines 315-319); (b) local
     pairs `j > i` (lines 322-343), each pair visited **once**, with the force and
     torque applied to both partners. The reaction torque uses the exact relation
     (lines 335-337):
     ```c++
     /* Conservation of angular momentum mandates that
      * 0 = t_i + r_ij x F_ij + t_j */
     fji.torque += vector_product(pf.f, rn) - pf.torque;
     ```
     and `fji.f -= pf.f` (Newton's third law). Both partners' accumulators are
     multiplied by `prefactor` and added to `force()`/`torque()`.
  3. `boost::mpi::wait_all` (line 350);
  4. **Phase B** (lines 356-374): loop over local `i` again, summing over the remote
     "red" range `[all_posmom.begin(), local_posmom_begin)` and "black" range
     `[local_posmom_end, all_posmom.end())`; every remote pair is visited **twice**
     (once on each owning rank), so only `i` is written — no reduction needed.
  5. If `ESPRESSO_DIPOLE_FIELD_TRACKING` and not GPU, calls
     `dipole_field_at_part_cpu()` (lines 375-379).
- **Energy** (`long_range_energy_cpu`, lines 387-416): waits for MPI **before** any
  computation (line 400, no overlap), then for each local `i` accumulates
  `image_sum(it, all_posmom.end(), it, ...)` with `pair_potential` — a triangular
  sum starting at `i` itself (self-image energy included). Returns `prefactor * u`
  (the per-rank partial; the caller reduces over MPI).
- **Field** (`dipole_field_at_part_cpu`, lines 428-457,
  `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING`): waits for MPI, then for each local `i`
  sums `dipole_field` over **all** particles (full range, self-primary excluded by
  `image_sum`) and **assigns** `(*p)->dip_fld() = prefactor * u`.

### 2.4 Existing in-repo Kokkos patterns (to follow)

- `src/core/CMakeLists.txt:79`: `espresso_core` already links
  `Cabana::Core Kokkos::kokkos OpenMP::OpenMP_CXX` — **no build-system changes**.
- Top-level `CMakeLists.txt:693-714`: `find_package(Kokkos 4.6 QUIET)` with a
  FetchContent fallback to Kokkos 5.0.2 (`Kokkos_ENABLE_OPENMP ON`,
  `Kokkos_ENABLE_DEPRECATED_CODE_4/5 OFF`). kokkos-simd (`<Kokkos_SIMD.hpp>`,
  `Kokkos::Experimental::native_simd`) is available in both. This would be the
  first use of kokkos-simd in the repo.
- `src/core/forces.cpp:209-262` (`reduce_cabana_forces_and_torques`): the canonical
  `Kokkos::Experimental::contribute(local_view, scatter_view)` +
  `Kokkos::parallel_for` reduction-into-particles + `Kokkos::fence()` pattern.
- `src/core/electrostatics/icc.cpp:134-143`: same `contribute` + reduction loop, via
  the `kokkos_parallel_range_for` helper
  (`src/core/short_range_cabana.hpp:41-53`, which falls back to a serial loop when
  `Kokkos::num_threads() == 1`).
- `src/core/cell_system/CellStructure.hpp:159-166`: view typedefs —
  `Kokkos::View<double *[3], Kokkos::LayoutRight>` for forces,
  `Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight>`,
  `memory_space = Kokkos::HostSpace`.

### 2.5 `Utils::Vector<T, N>` genericity audit (`src/utils/include/utils/Vector.hpp`)

`Utils::Vector<T, N>` (line 50) is templated on the scalar type, and most needed
operations are generic:

- dot product `operator*(Vector<T,N>, Vector<U,N>)` (lines 342-353) — generic, but
  initializes the accumulator with `R acc{};` (see §5.3);
- `norm2()` (line 158) — delegates to the dot product;
- `vector_product` (lines 373-377) — generic;
- `operator+`, `operator-`, unary `-`, `+=` (lines 260-287) — generic;
- vector/scalar `operator/` (lines 315-321) — **unconstrained** in `U`, works for a
  simd divisor;
- **scalar multiplication `operator*(U, Vector<T,N>)` / `(Vector<T,N>, U)` (lines
  289-306) is constrained with `requires(std::is_arithmetic_v<U>)`** —
  `native_simd<double>` is a class type, so expressions like `(a + b) * d` and
  `pe3 * m1` in `pair_force` do **not** compile for `T = native_simd<double>`
  without additional overloads (§3.2);
- componentwise `sqrt(Vector<T,N>)` (lines 365-371) — not used by these kernels;
- `norm()` (line 159) uses `std::sqrt` — not used by these kernels (they use
  `norm2()` + explicit `sqrt`).

## 3. Design

The port is organized into six units. Units 1-2 live in a new header; units 3-6 are
rewrites inside `dipolar_direct_sum.cpp`. Throughout, the execution space is the
host (`Kokkos::DefaultExecutionSpace`, which is OpenMP in this build — same as
`src/core/forces.cpp:227`); the target is CPU-only.

Common aliases (in the new header):

```c++
using simd_double = Kokkos::Experimental::native_simd<double>;
// width: simd_double::size() (constexpr)
```

### 3.1 Unit 1 — Templated pair kernels (new header `dipolar_direct_sum_kernels.hpp`)

**What:** the three pair kernels rewritten as `template <class T>` functions on
`Utils::Vector<T, 3>`, byte-for-byte the same math as lines 62-124 of the current
`.cpp`. Instantiated with `T = double` (scalar tail) and `T = simd_double` (vector
body). The math is branch-free, so no SIMD masking is needed anywhere.

**Interface:**

```c++
/** Force and torque of one pair interaction, generic scalar type. */
template <class T> struct PairForce {
  Utils::Vector<T, 3> f;
  Utils::Vector<T, 3> torque;
};

template <class T>
PairForce<T> pair_force(Utils::Vector<T, 3> const &d,
                        Utils::Vector<T, 3> const &m1,
                        Utils::Vector<T, 3> const &m2);

template <class T>
T pair_potential(Utils::Vector<T, 3> const &d, Utils::Vector<T, 3> const &m1,
                 Utils::Vector<T, 3> const &m2);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
template <class T>
Utils::Vector<T, 3> dipole_field(Utils::Vector<T, 3> const &d,
                                 Utils::Vector<T, 3> const &m1);
#endif
```

Notes:

- `PairForce<T>` replaces `ParticleForce` as the kernel return type because
  `ParticleForce` is `double`-only and its `torque` member is `#ifdef
  ESPRESSO_ROTATION` (`src/core/Particle.hpp:355-361`). Since `DIPOLES implies
  ROTATION`, `PairForce<T>::torque` is unconditional. Members are
  zero-initialized (`= Utils::Vector<T,3>::broadcast(T{0.})` or equivalent
  default-member-initializers) so accumulators start at zero for both `double` and
  simd instantiations (see §5.3).
- **ADL sqrt** — the one non-transparent spot: `pair_force` currently calls
  `std::sqrt(r2)` (line 68), which would never find the kokkos-simd overload. All
  three kernels must use an unqualified call after a `using`-declaration:
  ```c++
  using std::sqrt;
  auto const r = sqrt(r2);   // resolves to std::sqrt for double,
                             // to Kokkos::Experimental::sqrt via ADL for simd
  ```
- Header guards: file content behind `#ifdef ESPRESSO_DIPOLES` (after
  `#include <config/config.hpp>`), `dipole_field` additionally behind
  `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING`, mirroring the current `.cpp`.
- Required includes: `<Kokkos_SIMD.hpp>`, `<utils/Vector.hpp>`,
  `<config/config.hpp>`, `<cmath>`.

**Dependencies:** `Utils::Vector`, kokkos-simd, Unit 2 (the simd-scalar × vector
operators, which must be declared before the kernel templates are instantiated).

### 3.2 Unit 2 — simd-scalar × `Utils::Vector` operator overloads (same header)

**What:** `Utils::Vector`'s scalar multiplication requires
`std::is_arithmetic_v<U>` (`Vector.hpp:290-306`), which excludes
`native_simd<double>`. The kernel expressions `(a + b) * d`,
`3.0 * (pe3 * m1 + pe2 * m2) / r5`, `3.0 * pe3 * vector_product(m1, d) / r5` and
`3.0 * pe2 * d / r5` all form `simd * Vector<simd, 3>`. Rather than relaxing the
constraint in the core utils header (repo-wide blast radius), add narrowly-typed
overloads in the new kernels header, **in namespace `Utils`** so they are found by
ADL through the `Utils::Vector` operand:

**Interface:**

```c++
namespace Utils {

template <std::size_t N>
auto operator*(simd_double const &a, Vector<simd_double, N> const &b);
// -> Vector<simd_double, N>, componentwise a * b[i]

template <std::size_t N>
auto operator*(Vector<simd_double, N> const &a, simd_double const &b);
// symmetric form (not strictly needed by the current expressions; added for
// robustness of the shared templates)

} // namespace Utils
```

Notes:

- No overload is needed for `Vector<simd,N> / simd` — the existing unconstrained
  `operator/` (`Vector.hpp:315-321`) already handles it.
- No ambiguity with the dot product (`Vector.hpp:342-353`): that overload requires
  both arguments to be `Vector`s, and the arithmetic-constrained scalar overloads
  reject simd. `double * Vector<simd,N>` (e.g. `3.0 * (...)`) continues to use the
  existing arithmetic overload; its element type `R = decltype(double * simd)`
  resolves to `simd_double` via kokkos-simd's broadcasting hidden-friend
  operators.
- `Kokkos::Experimental::native_simd<double>` is an alias; the overloads are
  written against the alias so they track whatever concrete type
  (`basic_simd<double, Abi>`) the installed Kokkos provides (see §7.2).

**Dependencies:** `Utils::Vector`, kokkos-simd.

### 3.3 Unit 3 — SoA repack of the gathered particle data (`.cpp`)

**What:** `gather_particle_data` (lines 223-261) is kept **unchanged** (its output
order defines the `offset`/range bookkeeping and the MPI wire format). After it
returns, positions and moments are repacked from the AoS `std::vector<PosMom>` into
SoA Kokkos views so the SIMD `j`-loop can do contiguous `copy_from` loads:

**Interface:**

```c++
/** SoA copy of a slice of all_posmom, one view per field.
 *  LayoutLeft => component c is contiguous across particles, so
 *  simd.copy_from(&pos(j, c), simd_flag_default) loads pos_c for j..j+width-1. */
struct PosMomViews {
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> pos;
  Kokkos::View<double *[3], Kokkos::LayoutLeft, Kokkos::HostSpace> m;
};

/** Allocate views of size n_total (unmanaged by particles; plain allocation). */
static PosMomViews make_posmom_views(std::size_t n_total);

/** Copy all_posmom[begin, end) into the corresponding rows of the views. */
static void fill_posmom_views(PosMomViews &views,
                              std::vector<PosMom> const &all_posmom,
                              std::size_t begin, std::size_t end);
```

Notes:

- `fill_posmom_views` is called **twice** in the forces path to preserve the MPI
  overlap: once for the local slice `[offset, offset + n_local)` before Phase A,
  and once for the remainder (`[0, offset)` and `[offset + n_local, n_total)`)
  after `wait_all`. Energy and field kernels fill the whole range once, after
  `wait_all` (matching their current no-overlap structure).
- SIMD loads use
  `simd.copy_from(&views.pos(j, c), Kokkos::Experimental::simd_flag_default)`
  (element-aligned; no alignment guarantee is assumed for arbitrary `j`).
- The views are function-local (allocated per call, like `all_posmom` today). If
  profiling later shows allocation overhead, they can be cached on the actor; out
  of scope here.

**Dependencies:** `gather_particle_data` (unchanged), Kokkos views.

### 3.4 Unit 4 — Image-shift precomputation (`.cpp`)

**What:** the innermost loop is over periodic image shifts; `ncut` is typically 0
or small, so this loop **stays scalar** — vectorization is over the `j`-particle
dimension only. Instead of re-running `for_each_image` (lines 141-159, built on
`Utils::cartesian_product` + `std::views::iota`) inside every pair iteration,
precompute the shift vectors once per kernel invocation on the host:

**Interface:**

```c++
/** Real-space image shift vectors n .* box_l for all integer triples inside
 *  the |ncut| sphere. shifts[0] is the primary image (0, 0, 0), so the
 *  self-interaction loops can simply start at index 1. */
static std::vector<Utils::Vector3d>
make_image_shifts(Utils::Vector3i const &ncut, Utils::Vector3d const &box_l);
```

Notes:

- Implemented as a simple explicit triple loop over `[-ncut, ncut]^3` with the
  `nx² + ny² + nz² <= ncut²` sphere test (same predicate as `for_each_image`),
  pushing the zero shift first and the rest in the deterministic loop order.
  `for_each_image` and `Utils::cartesian_product` usage can then be removed from
  this file. (Rationale: capturing C++20 ranges machinery inside Kokkos lambdas is
  legal on host but noisy; a flat `std::vector` captured by pointer/size — or
  wrapped in an unmanaged view — is simpler and cheaper.)
- When `with_replicas == false` (`ncut == {0,0,0}`), the vector contains exactly
  the zero shift, and the kernels degenerate to the plain minimum-image sum, as
  today.

**Dependencies:** none beyond `Utils::Vector`.

### 3.5 Unit 5 — Distance computation helper: minimum image vs. raw (`.cpp`)

**What:** `image_sum` computes the primary distance as `it->pos - jt->pos` when
`with_replicas`, else `box_geo.get_mi_vector(...)` (lines 206-208).
`BoxGeometry::get_mi_vector` is a scalar host function (and also encapsulates
Lees-Edwards handling), so it cannot be templated on simd. The SIMD `j`-loop
therefore has two primary-distance paths:

- `with_replicas == true`: full SIMD — broadcast `pos_i`, `copy_from` the `j`-chunk,
  subtract: `d[c] = simd(pos_i[c]) - pos_j[c]`.
- `with_replicas == false`: per-lane scalar gather — compute
  `box_geo.get_mi_vector(pos_i, pos_j_lane)` for each of the `width` lanes into a
  small stack array `double buf[3][width]`, then `d[c].copy_from(&buf[c][0], ...)`.
  This preserves the exact minimum-image (and Lees-Edwards) semantics at the cost
  of scalarizing only the distance computation; the pair math after it is still
  SIMD.

The branch on `with_replicas` is per pair-chunk (uniform for the whole kernel), so
it is hoisted trivially by the compiler; no masking is involved.

**Interface (illustrative):**

```c++
/** Load the primary distance vectors of lanes j..j+width-1 relative to pos_i. */
static Utils::Vector<simd_double, 3>
primary_distance_simd(Utils::Vector3d const &pos_i, PosMomViews const &views,
                      std::size_t j, bool with_replicas,
                      BoxGeometry const &box_geo);
// plus the scalar counterpart for the tail loop (same logic, T = double),
// which is just the existing lines 206-208 / 324-325 expression.
```

**Dependencies:** Unit 3, `BoxGeometry`.

### 3.6 Unit 6 — The three Kokkos launch wrappers (`.cpp`)

All three methods keep their current signatures
(`src/core/magnetostatics/dipolar_direct_sum.hpp:81-88`); only their bodies change.
Common structure per method: `gather_particle_data` → `ncut`/`with_replicas` /
`make_image_shifts` → SoA fill → Kokkos launches → apply to particles. Local index
convention below: `n_local = local_particles.size()`, global index of local
particle `i` is `offset + i`, `n_total = all_posmom.size()`.

#### 3.6.1 Forces (`add_long_range_forces_cpu`, replaces lines 291-380)

Preserves the two-phase MPI overlap **ordering exactly**:

1. `gather_particle_data` (async `iall_gatherv` already in flight).
2. Allocate `PosMomViews` of size `n_total`; fill local slice
   `[offset, offset + n_local)` only.
3. Allocate accumulator views sized `n_local` (matching the CellStructure
   convention, `CellStructure.hpp:159-163`):
   ```c++
   using ForceView   = Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
   using ScatterForce = Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight>;
   ForceView local_force("dds_force", n_local), local_torque("dds_torque", n_local);
   ScatterForce scatter_force(local_force), scatter_torque(local_torque);
   ```
4. **Phase A** — `Kokkos::parallel_for("dds_local_pairs",
   RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_local), ...)`; per iteration `i`
   (global `gi = offset + i`):
   - registers `PairForce<double> fi{}` for particle `i`'s own accumulation;
   - **(a) self-images:** scalar loop over `shifts[1..]` (zero shift excluded, per
     `image_sum`'s `exclude_primary`), `fi += pair_force<double>(shift, m_i, m_i)`
     — equivalent to lines 315-319. Skipped entirely when
     `with_replicas == false` (empty shift tail);
   - **(b) local pairs, SIMD body:** for `j` from `gi + 1` in steps of
     `simd_double::size()` while a full chunk fits before `offset + n_local`:
     load `m_j` (3 × `copy_from`), compute the primary distance per §3.5,
     `PairForce<simd_double> fij{}, fji{};` then the scalar shift loop:
     ```c++
     for (auto const &shift : shifts) {          // rn = d + shift (simd)
       auto const pf = pair_force(rn, m_i_bcast, m_j);
       fij += pf;                                 // (or accumulate f/torque directly)
       fji.f -= pf.f;
       /* Conservation of angular momentum mandates that
        * 0 = t_i + r_ij x F_ij + t_j */
       fji.torque += vector_product(pf.f, rn) - pf.torque;
     }
     ```
     — the torque relation is preserved verbatim from lines 335-337, now evaluated
     per lane. After the shift loop:
     - `i`-side: horizontal-reduce each component of `fij` into `fi`
       (`Kokkos::Experimental::reduce(fij.f[c], std::plus<>{})`, likewise torque);
     - `j`-side: **per-lane scatter** — extract each lane `l` (partner
       `jl = j + l`, local index `jl - offset`) and write through the ScatterView
       access individually:
       ```c++
       auto force_access  = scatter_force.access();
       auto torque_access = scatter_torque.access();
       for (std::size_t l = 0; l < simd_double::size(); ++l)
         for (int c = 0; c < 3; ++c) {
           force_access(jl - offset, c)  += fji.f[c][l];
           torque_access(jl - offset, c) += fji.torque[c][l];
         }
       ```
       This is the standard SIMD+scatter pattern: SIMD math, scalar scatter.
   - **(c) scalar tail:** the remaining `(local_end - j) % width` partners via the
     `T = double` instantiation, identical structure, `j`-writes also through the
     ScatterView (single-lane);
   - **(d)** write `i`'s own total directly (each `i` is owned by exactly one
     iteration; scatter writes from other iterations go to duplicated buffers, so
     there is no race):
     `local_particles[i]->force() += prefactor * fi.f;`
     `local_particles[i]->torque() += prefactor * fi.torque;`
5. `boost::mpi::wait_all(reqs.begin(), reqs.end());` — unchanged position in the
   sequence (after Phase A submission **and completion**; host `parallel_for` is
   synchronous, matching the current semantics where local work fills the
   communication latency).
6. Fill the remote slices `[0, offset)` and `[offset + n_local, n_total)` of the
   SoA views.
7. **Phase B** — `Kokkos::parallel_for("dds_remote_pairs", ..., n_local)`; per `i`:
   two SIMD `j`-sweeps (red `[0, offset)`, black `[offset + n_local, n_total)`,
   mirroring lines 356-370), each = SIMD body + scalar tail, accumulating only
   `PairForce` for `i` (visit-twice — **no scatter needed**), then
   `local_particles[i]->force() += prefactor * fi.f;` etc. (direct write, one
   writer per `i`).
8. Reduce the Newton's-third-law contributions (pattern of
   `src/core/forces.cpp:215-253` and `icc.cpp:134-143`):
   ```c++
   Kokkos::Experimental::contribute(local_force, scatter_force);
   Kokkos::Experimental::contribute(local_torque, scatter_torque);
   Kokkos::parallel_for("dds_reduction", ..., n_local, /* per i: */
       // force()  += prefactor * {local_force(i,0..2)}
       // torque() += prefactor * {local_torque(i,0..2)}
   );
   Kokkos::fence();
   ```
   The `prefactor` multiplication happens here for the scattered `j`-side, and at
   steps 4(d)/7 for the `i`-side — in both cases applied exactly once per
   accumulated `PairForce`, as in the current code (lines 341-342, 345-346,
   372-373).
9. `ESPRESSO_DIPOLE_FIELD_TRACKING` trailer call unchanged (lines 375-379).

#### 3.6.2 Energy (`long_range_energy_cpu`, replaces lines 387-416)

The current code waits for MPI **before** any computation (line 400, no overlap).
This design **adds** the same two-phase MPI latency-hiding as the forces kernel
(user decision): the energy triangular sum for local particle `gi` runs over
`j ∈ [gi, n_total)`, which splits cleanly into a **local-upper** part
`[gi, offset + n_local)` — readable immediately from the local SoA slice — and a
**remote-black** part `[offset + n_local, n_total)` that needs the gathered data.
(The "red" range `[0, offset)` is *never* summed by the energy kernel — the primed
triangular sum only counts each global pair on the rank owning its lower index —
so the remote fill for energy is the black slice only.)

Sequence:

1. `gather_particle_data` (async `iall_gatherv` in flight).
2. Fill the **local** SoA slice `[offset, offset + n_local)` only.
3. **Phase A** — `Kokkos::parallel_reduce("dds_energy_local", 0, n_local, ..., uA)`
   over `j ∈ [gi, offset + n_local)`:
   ```c++
   // (a) self-images: scalar loop over shifts[1..] with pair_potential<double>
   // (b) local-upper j in (gi, offset + n_local): simd_double acc{0.} over the
   //     SIMD body (chunk loop: primary distance per §3.5, scalar shift loop,
   //     acc += pair_potential<simd_double>(rn, m_i, m_j));
   //     then uA_local += Kokkos::Experimental::reduce(acc, std::plus<>{});
   // (c) scalar tail with pair_potential<double>
   ```
4. `boost::mpi::wait_all(reqs.begin(), reqs.end());`
5. Fill the **black** SoA slice `[offset + n_local, n_total)`.
6. **Phase B** — `Kokkos::parallel_reduce("dds_energy_remote", 0, n_local, ..., uB)`
   over `j ∈ [offset + n_local, n_total)` (SIMD body + scalar tail; no self-image
   term, no primary exclusion — the range is entirely remote).
7. `return prefactor * (uA + uB);`

Together Phase A + Phase B reproduce the triangular primed sum
`image_sum(it, all_posmom.end(), it, ...)` of lines 408-413 exactly (same `j`-range
`[gi, n_total)`, same self-image handling), now with the local block computed while
the remote data transfers. The per-rank partial times `prefactor` is returned; MPI
reduction over ranks stays with the caller, as now.

#### 3.6.3 Field (`dipole_field_at_part_cpu`, replaces lines 428-457, `#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING`)

Same data staging as energy (wait first, full fill), then
`Kokkos::parallel_for("dds_dipole_field", ..., n_local)`; per `i`:

- SIMD sweep over **all** `j` in `[0, n_total)` with `dipole_field<simd_double>`,
  accumulated in a `Utils::Vector<simd_double, 3>`; the chunk containing
  `gi = offset + i` needs the self-primary exclusion — simplest faithful handling:
  process that single index scalar (skip lane `gi` in the chunk via splitting the
  range into `[0, gi)` and `[gi + 1, n_total)` plus a scalar self-image term over
  `shifts[1..]`, exactly the `exclude_primary` semantics of `image_sum`);
- horizontal-reduce per component, add the scalar tails;
- **assign** (not add): `local_particles[i]->dip_fld() = prefactor * u;` (line 454
  semantics — one writer per `i`, no scatter needed).

Note the current implementation resets `dip_fld` in `forces.cpp:146-151` before
force calculation; assignment semantics must stay.

## 4. Files created / modified

| File | Action | Content |
| --- | --- | --- |
| `src/core/magnetostatics/dipolar_direct_sum_kernels.hpp` | **new** | Units 1-2: `PairForce<T>`, templated `pair_force` / `pair_potential` / `dipole_field` (ADL `sqrt`), `simd_double` alias, `Utils::operator*` simd-scalar overloads. Guarded by `ESPRESSO_DIPOLES` / `ESPRESSO_DIPOLE_FIELD_TRACKING`. Follows the existing `*_kokkos.hpp` header convention (`bond_forces_kokkos.hpp` etc.) in spirit; named `*_kernels.hpp` since it also carries the scalar instantiation. |
| `src/core/magnetostatics/dipolar_direct_sum.cpp` | **modified** | Remove the three static `double` kernels (lines 62-124) and `for_each_image` (141-159); keep `PosMom`, `image_sum` is removed (superseded by the explicit loops), keep `gather_particle_data` and `get_n_cut` unchanged; add Units 3-5 (SoA views, `make_image_shifts`, distance helper); rewrite the three methods per §3.6. New includes: `dipolar_direct_sum_kernels.hpp`, `<Kokkos_Core.hpp>`, `<Kokkos_ScatterView.hpp>`, `<Kokkos_SIMD.hpp>`. Drop `<utils/cartesian_product.hpp>` and `<ranges>` if no longer used. |
| `src/core/magnetostatics/dipolar_direct_sum.hpp` | unchanged | Public interface identical. |
| CMake files | unchanged | Kokkos already linked (`src/core/CMakeLists.txt:79`). |

No other file changes. GPU sources (`dipolar_direct_sum_gpu*`) untouched.

## 5. Correctness considerations

1. **SIMD tail / remainder.** Every SIMD `j`-sweep has a scalar remainder of
   `count % simd_double::size()` iterations using the `T = double` instantiation of
   the *same* templated kernel — identical math, no drift between paths. Ranges
   that are shorter than one width (small N, few local particles, `n_local == 0`)
   degrade gracefully to pure scalar. Empty ranges (`n_local == 0`, e.g. a rank
   with no dipoles) must produce zero-length views and no-op launches — Kokkos
   handles 0-extent views/policies; `reqs` may still be non-empty and must still be
   waited on (as today).
2. **ADL `sqrt`.** `std::sqrt(r2)` (current line 68) is replaced by
   `using std::sqrt; ... sqrt(r2)` in all three templated kernels so
   `Kokkos::Experimental::sqrt(native_simd<double>)` is found by ADL for the simd
   instantiation while `double` keeps `std::sqrt`. This is the only spot where the
   existing kernel text is not type-transparent.
3. **Zero-initialization of simd accumulators.** `Utils::Vector`'s dot product
   initializes with `R acc{};` (`Vector.hpp:348`). For `R = native_simd<double>`
   this is value-initialization; Kokkos's simd types declare `simd() = default`
   over a trivially-copyable register member, so value-init zero-initializes —
   but this is an implementation detail of Kokkos, not a documented guarantee.
   Mitigation: all *new* accumulators in our code are explicitly zero-initialized
   (broadcast `T{0.}` / `PairForce<T>` default-member-initializers); for the dot
   product, add a compile-time/unit sanity check during implementation (e.g.
   assert `(simd_double{} == simd_double(0.)).all_of()`-style check in a debug
   path, or provide a local 3-component `dot` helper in the kernels header if the
   check fails on any supported Kokkos version).
4. **ScatterView correctness for the Newton's-third-law trick.** Phase A parallel
   iterations over `i` write to partner `j` indices owned by other iterations; all
   such writes go through `ScatterView::access()` (per-lane, scalar). On the
   OpenMP host backend the default contribution strategy is data duplication, so
   concurrent `+=` on the same `(j, c)` is safe and merged deterministically by
   `contribute` for a fixed thread count. Particle-`i`'s own accumulation is
   written directly by its unique owning iteration (no aliasing with the
   duplicated scatter buffers). Pattern matches `forces.cpp:209-253`.
5. **MPI overlap ordering.** The forces sequence must remain: async gather → local
   SoA fill → Phase A → `wait_all` → remote SoA fill → Phase B → contribute/apply.
   Reading `all_posmom` outside `[offset, offset + n_local)` before `wait_all` is
   undefined (the `iall_gatherv` buffer is being filled); the two-step
   `fill_posmom_views` (§3.3) is what enforces this. The energy kernel now uses the
   same discipline (§3.6.2): local slice → Phase A reduce → `wait_all` → black
   slice → Phase B reduce. The field kernel waits first, as today (line 440).
6. **`with_replicas` semantics.** `with_replicas == true` ⇒ raw distance
   `pos_i - pos_j` plus image shifts; `false` ⇒ `box_geo.get_mi_vector` minimum
   image (including its Lees-Edwards handling) with only the zero shift. The
   per-lane scalar gather (§3.5) keeps the minimum-image path bit-identical to the
   current per-pair computation.
7. **Self-interaction exclusion.** The primed sum excludes only the primary image
   of the self pair. Encoded as: shifts list with the zero shift at index 0;
   self-image loops iterate `shifts[1..]`; pair loops (`j > i`, red/black, field
   with `gi` split out) iterate all shifts.
8. **Torque conservation.** The reaction torque uses exactly
   `fji.torque += vector_product(pf.f, rn) - pf.torque` (current lines 335-337)
   with per-lane `rn`; do not "simplify" to a recomputed `pair_force(-rn, mj, mi)`
   — the current formulation guarantees `0 = t_i + r_ij × F_ij + t_j` exactly in
   floating point.
9. **`#ifdef` guards.** Whole TU and new header: `ESPRESSO_DIPOLES`. Field kernel
   and `dipole_field` template: `ESPRESSO_DIPOLE_FIELD_TRACKING`. Torque code needs
   no `ESPRESSO_ROTATION` guard inside this TU (`DIPOLES implies ROTATION`,
   `features.def:46`; `p.torque()` is already used unguarded at lines 342/346/373).
   `ESPRESSO_CUDA` GPU dispatch in `dipolar_direct_sum.hpp:64-86` and the GPU
   implementation files are untouched.
10. **Determinism / summation order.** Results change relative to the serial code
    (SIMD lane order, per-thread scatter duplicates), but remain deterministic for
    a fixed thread count and rank count. See §7.1 for test-tolerance implications.

## 6. Build & verification plan

- **Build:** CPU-only configuration (CUDA off), then `make -j8` in the build
  directory. No CMake changes; Kokkos ≥ 4.6 with OpenMP is already required and
  linked.
- **Python tests** (primary verification):
  - `testsuite/python/dipolar_direct_summation.py`
    - `test_dds_cpu` (line 173): CPU DDS with open boundaries vs. **stored
      reference data** (`dipolar_open_boundaries_{energy,arrays}.npy`), tolerances
      `1E-12` (energy/forces/torques) — the strictest check;
    - `test_min_image_convention_cpu` (line 252): `n_replicas = 0` minimum-image
      path vs. a python reference, `atol=1e-10` + rtol;
    - `test_inner_loop_consistency_cpu` (line 265): requires `n_nodes > 1`
      (`mpiexec -n 2` / `-n 4`), directly exercises the Phase A vs. Phase B
      (local-pair / remote-pair) split and hence the ScatterView + overlap logic,
      tolerance `1e-10`;
    - `test_gen_reference_data` / `test_dds_scafacos` as available in the config.
  - `testsuite/python/dawaanr-and-dds-gpu.py` is entirely gated on
    `@utx.skipIfMissingGPU()` + `CUDA` (lines 29-30): in the CPU-only build it is
    skipped. It remains relevant CI coverage on GPU-enabled configurations, where
    its "dawaanr" side runs the CPU code changed here against the untouched GPU
    implementation (`rtol=1e-2` forces/torques, `1e-5` energy).
  - Smoke coverage: `testsuite/python/dipolar_interface.py` (actor setup/getters).
  - Run at 1, 2, and 4 MPI ranks (the testsuite's standard parallel modes), and —
    since this introduces threading into a previously serial path — with
    `OMP_NUM_THREADS`/Kokkos threads > 1 explicitly, plus one run with 1 thread to
    isolate threading effects.
- **C++ unit tests:** there is **no** dedicated C++ unit test for
  `DipolarDirectSum` in `src/core/unit_tests` (verified by search); coverage is
  python-level only. Adding one is out of scope for this refactor (behavior is
  pinned by the stored-reference python test), but see §7.5.
- **Pass criteria:** all currently-passing tests in the above set still pass, at
  all rank/thread combinations tested; no new compiler warnings in the touched
  files; CPU-only build clean with `make -j8`.

## 7. Risks & open questions

1. **Tight `1E-12` absolute tolerances vs. reordered sums.** `test_dds_cpu`
   compares against `.npy` references generated by the current serial summation
   order; SIMD lane reduction and per-thread scatter merging change the order.
   For the N=20, O(1)-magnitude reference system the expected drift is
   ~1e-14…1e-13, comfortably inside `1E-12` — but this is the most likely test to
   flag trouble. If it fails marginally with correct physics, the options are (a)
   loosen the tolerance in the test (needs maintainer sign-off) or (b) regenerate
   the reference data via `gen_reference_data`; **do not** silently regenerate —
   escalate.
2. **`native_simd` alias availability.** The build allows either a system Kokkos
   ≥ 4.6 or fetched Kokkos 5.0.2 with `Kokkos_ENABLE_DEPRECATED_CODE_4/5 OFF`
   (`CMakeLists.txt:693-714`). Confirm during implementation that
   `Kokkos::Experimental::native_simd<double>` (rather than only
   `basic_simd<double, Abi>`) exists and is non-deprecated across both; the
   `simd_double` alias in the kernels header is the single point of change if the
   spelling differs. Same for the load-flag spelling
   (`simd_flag_default` vs. legacy `element_aligned_tag`) and
   `Kokkos::Experimental::reduce(simd, op)`.
3. **`Utils::Vector<simd, 3>` template hygiene.** `Vector` inherits from
   `Utils::Array<T, N>` and uses `std::ranges::transform`, `std::plus<T>`, etc. —
   all fine for a copyable class scalar. Watch items: the `DEVICE_QUALIFIER
   constexpr` on `Array`/`broadcast` (host-only build, `constexpr`-ness of simd
   default ctor is satisfied since it's defaulted); `boost::qvm` trait
   specializations at the bottom of `Vector.hpp` are instantiated only on demand
   and are never used for the simd instantiation. If an unexpected instantiation
   error surfaces, the fallback is a minimal local `struct Vec3<T>` in the kernels
   header — a design regression to avoid unless forced.
4. **SoA layout / `copy_from` fit.** `View<double*[3], LayoutLeft>` gives stride-1
   access in `j` for fixed component — verified against Kokkos layout rules, but if
   implementation testing shows otherwise (e.g. legacy-view padding surprises with
   `Kokkos_ENABLE_IMPL_VIEW_LEGACY ON`, `CMakeLists.txt:708`), the fallback is
   three separate `View<double*>`s per field (`pos_x/pos_y/pos_z`), which
   guarantees contiguity at the cost of clunkier repack code.
5. **ScatterView duplication cost.** Host duplication allocates
   `n_threads × n_local × 3` doubles each for force and torque per forces call.
   For typical DDS particle counts (direct sum is O(N²), so N is small-to-moderate
   by construction) this is negligible; for very high thread counts on large-N
   systems it is measurable. Acceptable for this refactor; a
   `ScatterNonDuplicated`(atomics) switch is a one-line change if ever needed.
6. **Per-lane `get_mi_vector` gather in the `n_replicas = 0` path** scalarizes the
   distance computation, capping the SIMD gain for the (common) minimum-image
   case. The pair math (sqrt, divisions — the expensive part) remains SIMD, so a
   net win is still expected. If profiling justifies it later, an orthorhombic
   fast path (`d -= box_l * round(d / box_l)` with kokkos-simd `round`, valid only
   without Lees-Edwards) can be added behind a runtime check — explicitly out of
   scope now.
7. **Threading a previously serial region.** Phase A/B `parallel_for` bodies write
   to `Particle` objects through `local_particles[i]` pointers. Each `i` is
   written by exactly one iteration, and `Particle` accessors are plain data — no
   hidden shared state. `Kokkos::fence()` after each launch (repo convention,
   `forces.cpp:253`) before touching results serially.
8. **Energy MPI latency-hiding (now in scope).** Per the user's decision, the
   energy kernel gains the same two-phase MPI overlap as forces (§3.6.2): the local
   triangular sum `[gi, offset + n_local)` is computed before `wait_all`, the
   remote-black sum `[offset + n_local, n_total)` after. This is a behavior change
   relative to the current serial-in-MPI energy path (line 400), but the summed
   result is identical up to floating-point reassociation (same caveat as §7.1 for
   the energy reference in `test_dds_cpu`). No red-range fill is needed for energy.
