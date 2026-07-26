# Async Ghost Communication — Ph0+Ph1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace ESPResSo's staggered, blocking, single-buffer ghost transport with an asynchronous, topology-agnostic engine driven by a declarative `HaloPlan`, reaching exact physics parity on all cell systems while adding Caliper instrumentation and a locked A/B baseline.

**Architecture:** A new `GhostComm::HaloPlan` (per-neighbor `{send cells+shift, recv ghost cells}` + same-rank local copies + optional n-square collective section) is emitted by each decomposition and executed by a new engine `GhostComm::halo_exchange[_start/_finish]` using non-blocking `Irecv`/`Isend`/`wait_all` with per-neighbor buffers. Migration is incremental behind the existing `CellStructure` ghost façade: a decomposition that returns a non-null `halo_plan()` is routed to the new engine; the legacy `ghost_communicator()` path stays until parity is proven, then is deleted.

**Tech Stack:** C++20, Boost.MPI (`boost::mpi::communicator::isend/irecv`, `boost::mpi::wait_all`), Boost.Test (`espresso_unit_test` CMake macro, `NUM_PROC` for multi-rank), Kokkos/Cabana (branch context), Caliper (`ESPRESSO_CALIPER`), Python testsuite via `pypresso`, `maintainer/benchmarks/{lj,p3m}.py`.

## Global Constraints

Every task's requirements implicitly include these (verbatim from the spec):

- **R5 — Preserve physics exactly.** Ghost *contents* and position convention are unchanged; **ghost positions stay unshifted** (LB particle coupling reads them assuming unshifted; shifts also entangle with Lees-Edwards). Do **not** activate meaningful `shift` values in Ph0/Ph1; plumb `shift` at its current value only.
- **R6 — Compatibility bar.** All green: regular / n-square / hybrid cell systems, Lees-Edwards, virtual sites, RATTLE, collision detection, LB coupling/tracers, bonds, `maxset` features, full testsuite, CI on the `me` remote.
- **R7 — Performance bar.** At least neutral (no regression) on LJ & P3M at 1/2/4/8 ranks; latency hiding measured via Caliper; a scaling win is desired, not a merge gate.
- **Build/workflow:** CPU-only build for this work (`-D ESPRESSO_BUILD_WITH_CUDA=OFF`). Use the `~/es-env` virtualenv for formatters/Python. Build with `make -j8` (not `-j$(nproc)`). Run `maintainer/format` on changed files before every push (CI won't start otherwise). Push finished work to the `me` remote only — **never** write to any other remote, and do not open PRs. Use `git -C <path>` rather than `cd`. Before any benchmark/profile run, confirm the machine is idle (no other users' heavy jobs).
- **No Cartesian assumptions** in `HaloPlan` or the engine — a future tree/SFC load balancer must plug into the same interface.

## Reference: current entry points (do not guess — these are verified)

- Engine to replace: `ghost_communicator()` — `src/core/ghosts.cpp:443`; data structures `src/core/ghosts.hpp:144-176`.
- Façade bodies — `src/core/cell_system/CellStructure.cpp:479-496`:
  ```cpp
  void CellStructure::ghosts_count() {
    ghost_communicator(decomposition().exchange_ghosts_comm(),
                       *get_system().box_geo, GHOSTTRANS_PARTNUM);
  }
  void CellStructure::ghosts_update(unsigned data_parts) {
    ghost_communicator(decomposition().exchange_ghosts_comm(),
                       *get_system().box_geo, map_data_parts(data_parts));
  }
  void CellStructure::ghosts_reduce_forces() {
    ghost_communicator(decomposition().collect_ghost_force_comm(),
                       *get_system().box_geo, GHOSTTRANS_FORCE);
  }
  ```
- Decomposition interface — `src/core/cell_system/ParticleDecomposition.hpp:79-83` (`exchange_ghosts_comm()`, `collect_ghost_force_comm()`).
- Regular construction to replace — `RegularDecomposition::prepare_comm()` `src/core/cell_system/RegularDecomposition.cpp:605`; neighbor stencil `init_cell_interactions` `RegularDecomposition.cpp:477-524`; `fully_connected_boundary()` sites `418,452,487,517`.
- Caliper idiom — `src/core/forces.cpp:57-59` (`#ifdef ESPRESSO_CALIPER / #include <caliper/cali.h> / #endif`), markers `CALI_CXX_MARK_FUNCTION;` and `CALI_MARK_BEGIN("x")/CALI_MARK_END("x")`.
- Ghost source registered at `src/core/CMakeLists.txt:32` (`ghosts.cpp`). Unit tests registered in `src/core/unit_tests/CMakeLists.txt` via `espresso_unit_test(SRC … DEPENDS espresso::core Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC N)`.

## File Structure (Ph0+Ph1)

- **Create** `src/core/ghosts/HaloPlan.hpp` — data model (`SendRegion`, `NeighborComm`, `LocalComm`, `CollectiveSection`, `HaloPlan`, `Direction`, `Combine`, `ExchangeOp`) in namespace `GhostComm`. Header-only.
- **Create** `src/core/ghosts/particle_packing.hpp` / `.cpp` — reusable pack/unpack extracted from `ghosts.cpp` (`CommBuf`, `serialize_and_reduce`, `calc_transmit_size`, `prepare_send_buffer`, `put_recv_buffer`, `add_forces_from_recv_buffer`, `add_rattle_correction_from_recv_buffer`, `cell_cell_transfer`). Used by both the legacy and new engine during parity.
- **Create** `src/core/ghosts/HaloExchange.hpp` / `.cpp` — the async engine (`GhostExchange` handle, `halo_exchange_start/finish`, blocking `halo_exchange`, collective handling).
- **Create** `src/core/unit_tests/HaloPlan_test.cpp`, `src/core/unit_tests/HaloExchange_test.cpp` — Boost.Test.
- **Modify** `src/core/CMakeLists.txt` — add the new `.cpp` files to the core target.
- **Modify** `src/core/unit_tests/CMakeLists.txt` — register the two new tests.
- **Modify** `src/core/cell_system/ParticleDecomposition.hpp` — add `virtual GhostComm::HaloPlan const *halo_plan() const { return nullptr; }`.
- **Modify** `src/core/cell_system/RegularDecomposition.{hpp,cpp}` — build+cache a `HaloPlan`, override `halo_plan()`.
- **Modify** `src/core/cell_system/AtomDecomposition.{hpp,cpp}`, `HybridDecomposition.{hpp,cpp}` — same, with collective section.
- **Modify** `src/core/cell_system/CellStructure.cpp` — route façade to the engine when `halo_plan()` is non-null; add Caliper markers.
- **Modify** `src/core/ghosts.cpp` / `ghosts.hpp` — Ph1 final task removes the legacy `ghost_communicator` and its helpers now living in `particle_packing`.
- **Modify** `testsuite/python/caliper.py` — extend `EXPECTED_LABELS` with ghost regions.
- **Create** `maintainer/benchmarks/ghost_ab.sh` — A/B runner + machine-idle gate.

---

## Phase 0 — Instrument the existing path and lock a baseline

No behavior change. Goal: ghost cost becomes visible and separable, and we have numbers to A/B against.

### Task 0.1: Caliper markers on the `CellStructure` ghost entry points

**Files:**
- Modify: `src/core/cell_system/CellStructure.cpp` (top-of-file includes; methods at `479-496` and `resort_particles` `512`, `update_ghosts_and_resort_particle` `593`)

**Interfaces:**
- Consumes: existing `ghost_communicator()`; nothing new.
- Produces: Caliper regions `ghosts_count`, `ghosts_update`, `ghosts_reduce_forces`, `resort_particles`, `update_ghosts_and_resort_particle` (used by the caliper test and A/B).

- [ ] **Step 1:** Add the Caliper include block after the existing includes in `src/core/cell_system/CellStructure.cpp`:

```cpp
#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif
```

- [ ] **Step 2:** Add `CALI_CXX_MARK_FUNCTION;` as the first statement of `ghosts_count()`, `ghosts_update(unsigned)`, `ghosts_reduce_forces()`, `resort_particles(bool)`, and `update_ghosts_and_resort_particle(unsigned)`. Guard each with `#ifdef ESPRESSO_CALIPER` only if the surrounding file does not already unconditionally define the macro to nothing — check `caliper/cali.h`: the macros are no-ops when Caliper is absent **only** if the header is included, so keep the include guarded and the marker calls unguarded (matches `forces.cpp`, which calls `CALI_CXX_MARK_FUNCTION;` unguarded). Example:

```cpp
void CellStructure::ghosts_update(unsigned data_parts) {
  CALI_CXX_MARK_FUNCTION;
  ghost_communicator(decomposition().exchange_ghosts_comm(),
                     *get_system().box_geo, map_data_parts(data_parts));
}
```

- [ ] **Step 3:** Configure a Caliper-enabled CPU build (out-of-source, reuse for the whole plan):

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm rev-parse --abbrev-ref HEAD   # expect worktree-comm
cmake -S /tikhome/weeber/es/.claude/worktrees/comm -B /tikhome/weeber/es/.claude/worktrees/comm/build \
  -D ESPRESSO_BUILD_WITH_CUDA=OFF -D ESPRESSO_BUILD_WITH_WALBERLA=ON \
  -D ESPRESSO_BUILD_WITH_CALIPER=ON -D ESPRESSO_BUILD_BENCHMARKS=ON \
  -D CMAKE_BUILD_TYPE=Release
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 espresso
```
Expected: builds clean.

- [ ] **Step 4:** Sanity-run a Caliper report and confirm the new labels appear:

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
CALI_CONFIG=runtime-report ./pypresso ../samples/lj_liquid.py 2>&1 | grep -iE 'ghosts_update|ghosts_reduce_forces|resort_particles'
```
Expected: the ghost region names show up in the report.

- [ ] **Step 5:** Format and commit:

```bash
~/es-env/bin/python -m pre_commit run --files src/core/cell_system/CellStructure.cpp || maintainer/format/format.sh src/core/cell_system/CellStructure.cpp
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/cell_system/CellStructure.cpp
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "core/ghosts: add Caliper markers to CellStructure ghost entry points"
```

### Task 0.2: Update the caliper test's expected label tree

**Files:**
- Modify: `testsuite/python/caliper.py` (`EXPECTED_LABELS`, lines ~29-54)

**Interfaces:**
- Consumes: the regions added in Task 0.1.
- Produces: a green `caliper.py` test that pins the new label hierarchy.

- [ ] **Step 1: Observe the actual tree.** Run the caliper child under a report and capture the label nesting the test parses:

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
CALI_CONFIG=runtime-report ./pypresso ../testsuite/python/caliper_child.py 2>&1 | sed -n '1,60p'
```
Note where `resort_particles_if_needed`/`calculate_forces` now show `ghosts_*` children.

- [ ] **Step 2: Run the test to see it fail** (mismatch between old `EXPECTED_LABELS` and the new tree):

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
ctest -R caliper --output-on-failure
```
Expected: FAIL with a label-set mismatch listing the new `ghosts_*` labels.

- [ ] **Step 3:** Edit `EXPECTED_LABELS` in `testsuite/python/caliper.py` to insert the observed ghost labels at their observed nesting (e.g. `resort_particles` under `resort_particles_if_needed`, `ghosts_update`/`ghosts_reduce_forces` under `calculate_forces` if that is where the report nests them — use the exact tree from Step 1).

- [ ] **Step 4: Run the test to verify it passes:**

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
ctest -R caliper --output-on-failure
```
Expected: PASS.

- [ ] **Step 5:** Commit:

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm add testsuite/python/caliper.py
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "testsuite: pin ghost-communication Caliper labels"
```

### Task 0.3: A/B runner with machine-idle gate + baseline capture

**Files:**
- Create: `maintainer/benchmarks/ghost_ab.sh`

**Interfaces:**
- Consumes: `build/pypresso`, `maintainer/benchmarks/{lj,p3m}.py`.
- Produces: `build/ghost_ab_<label>.csv` rows and a printed refusal if the machine is busy; a committed baseline CSV under `docs/superpowers/`.

- [ ] **Step 1:** Create `maintainer/benchmarks/ghost_ab.sh`:

```bash
#!/usr/bin/env bash
# Usage: ghost_ab.sh <label>   (run from the build dir)
# Runs lj.py and p3m.py on 1/2/4/8 ranks after checking the machine is idle.
set -euo pipefail
label="${1:?need a label, e.g. baseline or async}"
build="$(pwd)"

# --- machine-idle gate (per project constraint) ---
read -r load1 _ < <(awk '{print $1}' /proc/loadavg | tr '\n' ' ')
ncpu="$(nproc)"
others="$(who | awk '{print $1}' | sort -u | grep -v "^${USER}$" | wc -l)"
if (( $(echo "${load1} > 2.0" | bc -l) )) || (( others > 0 )); then
  echo "REFUSING to benchmark: load1=${load1}, other users=${others}, ncpu=${ncpu}." >&2
  echo "Re-run when the machine is idle." >&2
  exit 3
fi

out="${build}/ghost_ab_${label}.csv"
: > "${out}"
for ranks in 1 2 4 8; do
  mpiexec --bind-to core -n "${ranks}" ./pypresso ../maintainer/benchmarks/lj.py \
    --particles_per_core=1000 --volume_fraction=0.50 --output="${out}"
  mpiexec --bind-to core -n "${ranks}" ./pypresso ../maintainer/benchmarks/p3m.py \
    --particles_per_core=1000 --volume_fraction=0.25 --output="${out}"
done
echo "wrote ${out}"
```

- [ ] **Step 2:** Make it executable and confirm the idle gate refuses when busy (dry logic check):

```bash
chmod +x /tikhome/weeber/es/.claude/worktrees/comm/maintainer/benchmarks/ghost_ab.sh
bash -n /tikhome/weeber/es/.claude/worktrees/comm/maintainer/benchmarks/ghost_ab.sh   # syntax check
```
Expected: no output (syntax OK).

- [ ] **Step 3: Capture the baseline** (only if the idle gate passes):

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
../maintainer/benchmarks/ghost_ab.sh baseline
# per-region ghost breakdown at 4 ranks:
CALI_CONFIG=runtime-report mpiexec --bind-to core -n 4 ./pypresso ../maintainer/benchmarks/lj.py \
  --particles_per_core=1000 --volume_fraction=0.5 --output=/dev/null 2>&1 | tee ghost_regions_baseline_lj_4.txt
```
Expected: `ghost_ab_baseline.csv` populated; `ghost/…` and `ghosts_*` region times recorded.

- [ ] **Step 4:** Archive the baseline numbers in the repo for later comparison:

```bash
mkdir -p /tikhome/weeber/es/.claude/worktrees/comm/docs/superpowers/baselines
cp /tikhome/weeber/es/.claude/worktrees/comm/build/ghost_ab_baseline.csv \
   /tikhome/weeber/es/.claude/worktrees/comm/build/ghost_regions_baseline_lj_4.txt \
   /tikhome/weeber/es/.claude/worktrees/comm/docs/superpowers/baselines/
git -C /tikhome/weeber/es/.claude/worktrees/comm add maintainer/benchmarks/ghost_ab.sh docs/superpowers/baselines/
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "benchmarks: ghost A/B runner + machine-idle gate + Ph0 baseline"
```

> If the machine is busy, skip Step 3–4's measurement and record only the script; capture the baseline when idle before starting Ph1 comparisons. This is a hard gate — do not benchmark on a loaded machine.

---

## Phase 1 — Async engine + declarative `HaloPlan` at parity

### Task 1.1: `HaloPlan` data model

**Files:**
- Create: `src/core/ghosts/HaloPlan.hpp`
- Create: `src/core/unit_tests/HaloPlan_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Produces (relied on by all later tasks):
  - `namespace GhostComm { enum class Direction { Push, Reduce }; enum class Combine { Overwrite, Add }; struct ExchangeOp { Direction direction; Combine combine; }; }`
  - `struct SendRegion { ParticleList *cell; Utils::Vector3d shift; };`
  - `struct NeighborComm { int peer; std::vector<SendRegion> send; std::vector<ParticleList *> recv; };`
  - `struct LocalComm { ParticleList *src; ParticleList *dst; Utils::Vector3d shift; };`
  - `enum class CollectivePattern { None, Broadcast, ReduceSum };`
  - `struct CollectiveSection { CollectivePattern pattern; std::vector<ParticleList *> cells; };`
  - `struct HaloPlan { boost::mpi::communicator comm; std::vector<NeighborComm> neighbors; std::vector<LocalComm> local; std::optional<CollectiveSection> collective; };`

- [ ] **Step 1: Write the failing test** `src/core/unit_tests/HaloPlan_test.cpp`:

```cpp
#define BOOST_TEST_MODULE HaloPlan test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ghosts/HaloPlan.hpp"
#include "ParticleList.hpp"

BOOST_AUTO_TEST_CASE(neighbor_comm_shapes) {
  using namespace GhostComm;
  ParticleList a, b;
  NeighborComm nc{/*peer*/ 1, {SendRegion{&a, {}}}, {&b}};
  BOOST_CHECK_EQUAL(nc.peer, 1);
  BOOST_CHECK_EQUAL(nc.send.size(), nc.recv.size()); // recv[k] <-> peer.send[k]
  ExchangeOp op{Direction::Reduce, Combine::Add};
  BOOST_CHECK(op.direction == Direction::Reduce);
}
```

- [ ] **Step 2: Register the test** — add to `src/core/unit_tests/CMakeLists.txt` (near the other `espresso::core` tests):

```cmake
espresso_unit_test(SRC HaloPlan_test.cpp DEPENDS espresso::core Kokkos::kokkos
                   Boost::mpi MPI::MPI_CXX)
```

- [ ] **Step 3: Run it to verify it fails** (header missing):

```bash
cmake --build /tikhome/weeber/es/.claude/worktrees/comm/build --target HaloPlan_test 2>&1 | tail -20
```
Expected: FAIL — `ghosts/HaloPlan.hpp` not found.

- [ ] **Step 4: Create `src/core/ghosts/HaloPlan.hpp`** with the license header (copy from any core header) and:

```cpp
#pragma once
#include "ParticleList.hpp"
#include <utils/Vector.hpp>
#include <boost/mpi/communicator.hpp>
#include <optional>
#include <vector>

namespace GhostComm {
enum class Direction { Push, Reduce };   // real->ghost ; ghost->real
enum class Combine { Overwrite, Add };
struct ExchangeOp { Direction direction; Combine combine; };

struct SendRegion { ParticleList *cell; Utils::Vector3d shift; };
struct NeighborComm {
  int peer;
  std::vector<SendRegion> send;
  std::vector<ParticleList *> recv;      // recv[k] <-> peer.send[k]
};
struct LocalComm { ParticleList *src; ParticleList *dst; Utils::Vector3d shift; };

enum class CollectivePattern { None, Broadcast, ReduceSum };
struct CollectiveSection { CollectivePattern pattern; std::vector<ParticleList *> cells; };

struct HaloPlan {
  boost::mpi::communicator comm;
  std::vector<NeighborComm> neighbors;
  std::vector<LocalComm> local;
  std::optional<CollectiveSection> collective;
};
} // namespace GhostComm
```

- [ ] **Step 5: Run to verify it passes:**

```bash
cmake --build /tikhome/weeber/es/.claude/worktrees/comm/build --target HaloPlan_test -j8 && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R HaloPlan_test --output-on-failure
```
Expected: PASS.

- [ ] **Step 6: Commit:**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/ghosts/HaloPlan.hpp src/core/unit_tests/HaloPlan_test.cpp src/core/unit_tests/CMakeLists.txt
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "core/ghosts: add topology-agnostic HaloPlan data model"
```

### Task 1.2: Extract reusable particle packing (refactor, no behavior change)

**Files:**
- Create: `src/core/ghosts/particle_packing.hpp` / `.cpp`
- Modify: `src/core/ghosts.cpp` (include the new header; remove the moved definitions), `src/core/CMakeLists.txt`

**Interfaces:**
- Produces (relied on by the engine, Task 1.3):
  - `namespace GhostComm { class CommBuf { … }; }` (moved verbatim from `ghosts.cpp:64-91`)
  - `std::size_t GhostComm::calc_transmit_size(BoxGeometry const &, unsigned data_parts);`
  - `std::size_t GhostComm::calc_transmit_size(std::span<ParticleList *const>, BoxGeometry const &, unsigned data_parts);`
  - `void GhostComm::pack_cells(CommBuf &, std::span<ParticleList *const> cells, Utils::Vector3d const &shift, BoxGeometry const &, unsigned data_parts);` (from `prepare_send_buffer`)
  - `void GhostComm::unpack_cells(CommBuf &, std::span<ParticleList *const>, BoxGeometry const &, unsigned data_parts);` (from `put_recv_buffer`)
  - `void GhostComm::add_forces(CommBuf &, std::span<ParticleList *const>);` (from `add_forces_from_recv_buffer`)
  - `void GhostComm::add_rattle(CommBuf &, std::span<ParticleList *const>);` (BOND_CONSTRAINT only)
  - `void GhostComm::local_cell_copy(ParticleList &src, ParticleList &dst, Utils::Vector3d const &shift, BoxGeometry const &, unsigned data_parts);` (from `cell_cell_transfer`, one src/dst pair)
- Consumes: `serialize_and_reduce` and friends (move with them).

- [ ] **Step 1:** Move `CommBuf`, `SerializationSizeCalculator`, `ReductionPolicy`, `SerializationDirection`, `serialize_and_reduce`, both `calc_transmit_size` overloads, `prepare_send_buffer`, `prepare_ghost_cell`, `prepare_recv_buffer`, `put_recv_buffer`, `add_forces_from_recv_buffer`, `add_rattle_correction_from_recv_buffer`, and the body of `cell_cell_transfer` from `ghosts.cpp` into `particle_packing.{hpp,cpp}`, wrapped in `namespace GhostComm`. Give the public functions the signatures listed in **Interfaces** (thin wrappers over the existing statics; keep the free `ParticleList *`-vector variants as `std::span` overloads). Keep the exact serialization order — parity depends on byte-identical layout.

- [ ] **Step 2:** In `ghosts.cpp`, `#include "ghosts/particle_packing.hpp"` and rewrite `ghost_communicator()` to call the moved functions via `GhostComm::` (the legacy engine keeps working unchanged, just delegating packing). Add both new sources to `src/core/CMakeLists.txt` next to `ghosts.cpp`.

- [ ] **Step 3: Build the core and run the packing-adjacent unit tests** (parity guard — layout must be unchanged):

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 espresso Particle_serialization_test link_cell_test && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R 'Particle_serialization|link_cell' --output-on-failure
```
Expected: PASS (pure move; no semantic change).

- [ ] **Step 4: Run a quick Python smoke test** to confirm ghost physics is unchanged:

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
mpiexec --bind-to core -n 2 ./pypresso ../testsuite/python/lj.py
```
Expected: PASS.

- [ ] **Step 5:** Format changed files and commit:

```bash
maintainer/format/format.sh src/core/ghosts/particle_packing.hpp src/core/ghosts/particle_packing.cpp src/core/ghosts.cpp || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/ghosts/particle_packing.* src/core/ghosts.cpp src/core/CMakeLists.txt
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "core/ghosts: extract reusable particle packing (no behavior change)"
```

### Task 1.3: Async engine `halo_exchange` (point-to-point + local; collective stubbed)

**Files:**
- Create: `src/core/ghosts/HaloExchange.hpp` / `.cpp`
- Create: `src/core/unit_tests/HaloExchange_test.cpp`
- Modify: `src/core/CMakeLists.txt`, `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `HaloPlan` (1.1), `GhostComm::CommBuf`, `pack_cells`, `unpack_cells`, `add_forces`, `add_rattle`, `local_cell_copy` (1.2).
- Produces (relied on by the façade, Task 1.5):
  - `struct GhostComm::GhostExchange { … };` (opaque: per-neighbor `CommBuf send/recv`, `std::vector<boost::mpi::request>`, the resolved op + data_parts + a pointer to the plan)
  - `GhostExchange GhostComm::halo_exchange_start(HaloPlan const &, BoxGeometry const &, unsigned data_parts, ExchangeOp);`
  - `void GhostComm::halo_exchange_finish(GhostExchange &);`
  - `void GhostComm::halo_exchange(HaloPlan const &, BoxGeometry const &, unsigned data_parts, ExchangeOp);`

- [ ] **Step 1: Write the failing multi-rank test** `src/core/unit_tests/HaloExchange_test.cpp`. It builds a 2-rank hand-made plan where each rank sends one cell of particles to the other and checks the ghost cell is filled:

```cpp
#define BOOST_TEST_MODULE HaloExchange test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpi.hpp>

#include "ghosts/HaloExchange.hpp"
#include "ghosts/HaloPlan.hpp"
#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(push_positions_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box; box.set_length({10., 10., 10.});

  ParticleList local, ghost;
  local.resize(1);
  local.front().id() = me;                  // tag by owner rank
  local.front().pos() = {double(me), 0., 0.};
  ghost.resize(1);                          // sized as if ghosts_count ran

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_POSITION, {Direction::Push, Combine::Overwrite});

  BOOST_CHECK_CLOSE(ghost.front().pos()[0], double(other), 1e-12);
}
```

- [ ] **Step 2: Register** in `src/core/unit_tests/CMakeLists.txt`:

```cmake
espresso_unit_test(SRC HaloExchange_test.cpp DEPENDS espresso::core Kokkos::kokkos
                   Boost::mpi MPI::MPI_CXX NUM_PROC 2)
```

- [ ] **Step 3: Run to verify it fails** (engine missing):

```bash
cmake --build /tikhome/weeber/es/.claude/worktrees/comm/build --target HaloExchange_test 2>&1 | tail -20
```
Expected: FAIL — `ghosts/HaloExchange.hpp` not found.

- [ ] **Step 4: Implement `HaloExchange.hpp`** (declarations from **Interfaces**) and `HaloExchange.cpp`. Core algorithm — non-blocking, deadlock-free, sizes known a-priori from `recv` cell sizes:

```cpp
// halo_exchange_start
GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op) {
  GhostExchange st;
  st.op = op; st.data_parts = data_parts; st.box = &box; st.plan = &plan;
  if (data_parts == GHOSTTRANS_NONE) return st;

  auto const &comm = plan.comm;
  auto const n = plan.neighbors.size();
  st.send.resize(n); st.recv.resize(n);

  // For Reduce, send/recv roles swap: we send from *ghost* cells and receive
  // into *local* cells, adding on arrival.
  auto send_cells = [&](NeighborComm const &nc) {
    return op.direction == Direction::Push ? /*SendRegion cells*/ /*…*/ ;
  };
  // 1) post all receives (size = packed size of the cells we will receive into)
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    auto recv_cells = (op.direction == Direction::Push) ? span_of(nc.recv)
                                                        : span_of_send(nc.send);
    st.recv[i].resize(calc_transmit_size(recv_cells, box, data_parts));
    st.requests.push_back(
        comm.irecv(nc.peer, TAG, st.recv[i].data(), int(st.recv[i].size())));
  }
  // 2) pack and post all sends
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    if (op.direction == Direction::Push)
      pack_regions(st.send[i], nc.send, box, data_parts);        // applies per-region shift
    else
      pack_cells(st.send[i], span_of(nc.recv), {}, box, data_parts);
    st.requests.push_back(
        comm.isend(nc.peer, TAG, st.send[i].data(), int(st.send[i].size())));
  }
  return st;
}

void halo_exchange_finish(GhostExchange &st) {
  if (st.data_parts == GHOSTTRANS_NONE) return;
  // same-rank copies can run while messages are in flight
  for (auto const &lc : st.plan->local)
    local_cell_copy(*lc.src, *lc.dst, lc.shift, *st.box, st.data_parts);
  if (st.plan->collective) run_collective(*st.plan, *st.box, st.data_parts, st.op); // Task 1.6
  boost::mpi::wait_all(st.requests.begin(), st.requests.end());
  // unpack / reduce
  for (std::size_t i = 0; i < st.plan->neighbors.size(); ++i) {
    auto const &nc = st.plan->neighbors[i];
    auto dst = (st.op.direction == Direction::Push) ? span_of(nc.recv)
                                                    : span_of_send(nc.send);
    if (st.op.combine == Combine::Add && st.data_parts == GHOSTTRANS_FORCE)
      add_forces(st.recv[i], dst);
    else if (st.op.combine == Combine::Add && st.data_parts == GHOSTTRANS_RATTLE)
      add_rattle(st.recv[i], dst);
    else
      unpack_cells(st.recv[i], dst, *st.box, st.data_parts);
  }
}

void halo_exchange(HaloPlan const &p, BoxGeometry const &b, unsigned dp, ExchangeOp op) {
  auto st = halo_exchange_start(p, b, dp, op);
  halo_exchange_finish(st);
}
```
Notes for the implementer:
- `TAG` is a small `enum` per data-part kind (position/force/partnum/bonds) so future concurrent exchanges don't alias; one constant is sufficient within a single exchange since each message is matched by `(peer, tag)`.
- `span_of(std::vector<ParticleList*> const&)` → `std::span<ParticleList *const>`; `span_of_send` extracts the `.cell` pointers from a `std::vector<SendRegion>` into a scratch `std::vector<ParticleList*>` on the handle.
- `pack_regions` packs each `SendRegion` with its own `shift` (currently all zero, R5); implement as a loop calling the single-shift `pack_cells` per region into the same growing `CommBuf`, or extend `pack_cells` to take regions. Keep byte layout identical to the legacy packer.
- **Bonds:** if `data_parts & GHOSTTRANS_BONDS`, send/recv the `CommBuf::bonds()` vector as a *second* per-neighbor message (mirrors the legacy two-transfer split, `ghosts.cpp:497-506`). This is the cold resort path; folding into one message is a later optimization, out of scope for parity.
- **PARTNUM:** when `data_parts & GHOSTTRANS_PARTNUM`, recv sizes are not yet known — post a fixed-size recv of `sizeof(unsigned) * recv_cells.size()` and let `unpack_cells` resize the ghost cells (legacy `put_recv_buffer` PARTNUM path). This is why `ghosts_count()` must run before data exchanges (unchanged ordering).

- [ ] **Step 5: Run to verify it passes** (2 ranks):

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 HaloExchange_test && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R HaloExchange_test --output-on-failure
```
Expected: PASS on 2 ranks.

- [ ] **Step 6: Add a reduce (force) test case** to the same file — build a plan, put a force on the ghost cell on each rank, run `{Direction::Reduce, Combine::Add}`, and assert the owner's local force is the *sum*. Run and confirm PASS.

- [ ] **Step 7:** Format and commit:

```bash
maintainer/format/format.sh src/core/ghosts/HaloExchange.hpp src/core/ghosts/HaloExchange.cpp src/core/unit_tests/HaloExchange_test.cpp || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/ghosts/HaloExchange.* src/core/unit_tests/HaloExchange_test.cpp src/core/unit_tests/CMakeLists.txt src/core/CMakeLists.txt
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "core/ghosts: async halo_exchange engine (p2p + local, split-phase)"
```

### Task 1.4: `RegularDecomposition::make_halo_plan()` (direct-neighbor) + coverage check

**Files:**
- Modify: `src/core/cell_system/ParticleDecomposition.hpp` (add `virtual GhostComm::HaloPlan const *halo_plan() const { return nullptr; }`; forward-declare/`#include "ghosts/HaloPlan.hpp"`)
- Modify: `src/core/cell_system/RegularDecomposition.{hpp,cpp}`
- Create: `src/core/unit_tests/RegularHaloPlan_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `HaloPlan` (1.1); the existing grid/neighbor machinery (`init_cell_interactions`, `cell_grid`, `ghost_cell_grid`, `cart_neighbors`, `m_box`, `comm_cart`).
- Produces: `GhostComm::HaloPlan const *RegularDecomposition::halo_plan() const override;` returning a cached member `m_halo_plan`, rebuilt in `resort()` where `prepare_comm()` is called today (`RegularDecomposition.cpp:727`).

- [ ] **Step 1: Write the failing coverage test** `src/core/unit_tests/RegularHaloPlan_test.cpp` (4 ranks). It creates a `RegularDecomposition` on a periodic box, then asserts: (a) every ghost cell is the `recv` of **exactly one** `NeighborComm`; (b) the union of `recv` cells equals the set of `ghost_cells()`. Use the same construction pattern as `Verlet_list_test.cpp` (which already builds a `RegularDecomposition` at `NUM_PROC 4`):

```cpp
#define BOOST_TEST_MODULE RegularHaloPlan test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpi.hpp>
#include "cell_system/RegularDecomposition.hpp"
#include "ghosts/HaloPlan.hpp"
#include <set>
// … build a RegularDecomposition exactly as Verlet_list_test.cpp does …

BOOST_AUTO_TEST_CASE(every_ghost_cell_covered_once) {
  // auto dd = make_regular_decomposition(comm_cart, box, cutoff);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  std::multiset<ParticleList *> recvd;
  for (auto const &nc : plan->neighbors)
    for (auto *c : nc.recv) recvd.insert(c);
  for (auto const &lc : plan->local) recvd.insert(lc.dst);
  // (a) no ghost cell filled twice
  for (auto *c : recvd) BOOST_CHECK_EQUAL(recvd.count(c), 1u);
  // (b) coverage: every ghost cell is filled
  for (auto *g : dd.ghost_cells()) BOOST_CHECK_EQUAL(recvd.count(&g->particles()), 1u);
}
```
(Adjust `&g->particles()` to whatever `ParticleList *` the plan stores — match the type used in `make_halo_plan`.)

- [ ] **Step 2: Register** with `NUM_PROC 4` in `src/core/unit_tests/CMakeLists.txt`:

```cmake
espresso_unit_test(SRC RegularHaloPlan_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 4)
```

- [ ] **Step 3: Run to verify it fails** (`halo_plan()` returns nullptr → `BOOST_REQUIRE` trips):

```bash
cmake --build /tikhome/weeber/es/.claude/worktrees/comm/build --target RegularHaloPlan_test 2>&1 | tail -20
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R RegularHaloPlan_test --output-on-failure || true
```
Expected: FAIL (plan is null).

- [ ] **Step 4: Implement `make_halo_plan()`** in `RegularDecomposition.cpp`. Algorithm (direct-neighbor; reuses the grid the current code already computes):

  1. For each local cell `c` (grid index `m`), iterate its 3×3×3 stencil offsets `o ∈ {-1,0,1}³ \ {0}`; the neighbor cell index is `m + o`.
  2. If the neighbor index lies in the local region, skip (it's a real-real neighbor).
  3. Otherwise it maps to a ghost cell `g` on the ghost grid, owned by the peer obtained from the Cartesian shift of `o` (existing `comm_cart`/`node_grid` logic — the same `cart_neighbors` used in `prepare_comm`, but composed across all 3 axes to get the diagonal peer). Record the pair `(send: c with periodic shift for `o`, recv: g)` bucketed by `peer`.
  4. If `peer == this_node` (single rank in one or more dims, periodic self-ghost), record a `LocalComm{src: c, dst: g, shift}` instead (replaces `GHOST_LOCL`).
  5. Group buckets → one `NeighborComm{peer, send, recv}` per distinct peer; **sort** each `send`/`recv` pair list by the ghost cell's global index so both ranks agree on `recv[k] ⟷ peer.send[k]` (Section 5.3 ordering contract).
  6. Set `shift` from the existing per-offset periodic image shift the decomposition already knows (keep it at the value the legacy code used — R5: effectively the current behavior; do **not** introduce new shifts).
  7. Cache in `m_halo_plan`; call this from `resort()` right where `m_exchange_ghosts_comm = prepare_comm();` is today (`RegularDecomposition.cpp:727`). Keep `prepare_comm()` too for now (legacy path still available).
  8. Add a debug-only coverage assertion (the multiset check from the test) under `#ifdef ESPRESSO_ADDITIONAL_CHECKS`.

- [ ] **Step 5: Run to verify it passes** (4 ranks):

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 RegularHaloPlan_test && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R RegularHaloPlan_test --output-on-failure
```
Expected: PASS.

- [ ] **Step 6:** Format and commit:

```bash
maintainer/format/format.sh src/core/cell_system/ParticleDecomposition.hpp src/core/cell_system/RegularDecomposition.cpp src/core/cell_system/RegularDecomposition.hpp src/core/unit_tests/RegularHaloPlan_test.cpp || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/cell_system/ParticleDecomposition.hpp src/core/cell_system/RegularDecomposition.* src/core/unit_tests/RegularHaloPlan_test.cpp src/core/unit_tests/CMakeLists.txt
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "cell_system: RegularDecomposition emits a direct-neighbor HaloPlan"
```

### Task 1.5: Route the `CellStructure` façade to the engine (Regular only) — parity

**Files:**
- Modify: `src/core/cell_system/CellStructure.cpp` (`ghosts_count`, `ghosts_update`, `ghosts_reduce_forces`, `ghosts_reduce_rattle_correction`)

**Interfaces:**
- Consumes: `decomposition().halo_plan()` (1.4), `halo_exchange()` (1.3), `map_data_parts` (existing).
- Produces: unchanged façade signatures; behavior identical to legacy for Regular, still legacy for Atom/Hybrid (their `halo_plan()` is still null).

- [ ] **Step 1:** In each façade method, branch on the plan:

```cpp
void CellStructure::ghosts_update(unsigned data_parts) {
  CALI_CXX_MARK_FUNCTION;
  auto const parts = map_data_parts(data_parts);
  if (auto const *plan = decomposition().halo_plan()) {
    GhostComm::halo_exchange(*plan, *get_system().box_geo, parts,
                             {GhostComm::Direction::Push, GhostComm::Combine::Overwrite});
  } else {
    ghost_communicator(decomposition().exchange_ghosts_comm(),
                       *get_system().box_geo, parts);
  }
}
```
Do the same for `ghosts_count` (`GHOSTTRANS_PARTNUM`, Push/Overwrite), `ghosts_reduce_forces` (`GHOSTTRANS_FORCE`, `{Reduce, Add}`), and `ghosts_reduce_rattle_correction` (`GHOSTTRANS_RATTLE`, `{Reduce, Add}`).

- [ ] **Step 2: Build and run the full C++ unit tests:**

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -L unit_tests --output-on-failure
```
Expected: all PASS.

- [ ] **Step 3: Run the parity-critical Python tests** (regular decomposition path) on 1/2/4 ranks:

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
for t in lj p3m_fft lees_edwards virtual_sites_relative rattle collision_detection lb_tracers observable_cylindricalLB; do
  for n in 1 2 4; do
    echo "== $t on $n ranks =="; mpiexec --bind-to core -n $n ./pypresso ../testsuite/python/$t.py || echo "FAILED $t $n"
  done
done
```
Expected: PASS (or "test not present" for any name mismatch — verify actual filenames in `testsuite/python/`). LB coupling/tracers and Lees-Edwards passing here is the R5 guard.

- [ ] **Step 4: A/B measurement** (idle machine): build is `async` variant now.

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
../maintainer/benchmarks/ghost_ab.sh async
python3 - <<'PY'
import pandas as pd, glob
b=pd.read_csv('ghost_ab_baseline.csv'); a=pd.read_csv('ghost_ab_async.csv')
print(b.groupby(['script','ranks'])['mean'].mean()); print('---'); print(a.groupby(['script','ranks'])['mean'].mean())
PY
```
Expected: async mean per step ≤ baseline + CI at every rank count (R7 neutral). Record the numbers in the commit message.

- [ ] **Step 5:** Commit:

```bash
maintainer/format/format.sh src/core/cell_system/CellStructure.cpp || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/cell_system/CellStructure.cpp
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "cell_system: route regular ghost exchange through async engine (parity)"
```

- [ ] **Step 6: Push to `me` and watch CI:**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm push me HEAD:async-ghost-comm
```
Then monitor the pipeline on the `me` remote; triage any regression inline before continuing.

### Task 1.6: Migrate `AtomDecomposition` + `HybridDecomposition` (collective section)

**Files:**
- Modify: `src/core/cell_system/AtomDecomposition.{hpp,cpp}`, `HybridDecomposition.{hpp,cpp}`, `src/core/ghosts/HaloExchange.cpp` (`run_collective`)

**Interfaces:**
- Consumes: `CollectiveSection` (1.1), the current `AtomDecomposition::prepare_comm` semantics (one `GHOST_BCST` per rank for exchange, one `GHOST_RDCE` per rank for collect).
- Produces: `AtomDecomposition::halo_plan()` and `HybridDecomposition::halo_plan()` overrides; `run_collective()` in the engine.

- [ ] **Step 1:** Implement `run_collective(plan, box, data_parts, op)` in `HaloExchange.cpp`, mirroring the legacy `GHOST_BCST`/`GHOST_RDCE` loop (`ghosts.cpp:507-531`): for `Broadcast`, loop `root = 0..size-1` doing `boost::mpi::broadcast` of the packed owned cell(s); for `ReduceSum`, `boost::mpi::reduce` with `std::plus`. Reuse the packing helpers.

- [ ] **Step 2:** Implement `AtomDecomposition::make_halo_plan()`: `HaloPlan` with empty `neighbors`, a `CollectiveSection{Broadcast, {owned cell}}` for push, executed as reduce for the force direction (engine picks pattern from `ExchangeOp.direction`). Cache + override `halo_plan()`.

- [ ] **Step 3:** Implement `HybridDecomposition::make_halo_plan()`: concatenate the regular child's `neighbors`/`local` with the n-square child's `collective` into one `HaloPlan`. Cache + override.

- [ ] **Step 4: Write/extend tests** — add an n-square case to `HaloExchange_test.cpp` (2 ranks, collective broadcast+reduce round-trip) and run:

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 HaloExchange_test && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -R HaloExchange_test --output-on-failure
```
Expected: PASS.

- [ ] **Step 5: Run the n-square / hybrid Python tests:**

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
for n in 1 2 4; do mpiexec --bind-to core -n $n ./pypresso ../testsuite/python/cell_system.py; done
grep -rl "n_square\|hybrid" ../testsuite/python | while read f; do mpiexec --bind-to core -n 2 ./pypresso "$f" || echo "FAILED $f"; done
```
Expected: PASS.

- [ ] **Step 6:** Format, commit, push to `me`, watch CI:

```bash
maintainer/format/format.sh src/core/cell_system/AtomDecomposition.* src/core/cell_system/HybridDecomposition.* src/core/ghosts/HaloExchange.cpp || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add src/core/cell_system/AtomDecomposition.* src/core/cell_system/HybridDecomposition.* src/core/ghosts/HaloExchange.cpp src/core/unit_tests/HaloExchange_test.cpp
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "cell_system: n-square and hybrid emit HaloPlan collective section"
git -C /tikhome/weeber/es/.claude/worktrees/comm push me HEAD:async-ghost-comm
```

### Task 1.7: Delete the legacy engine and staggered construction

**Files:**
- Modify: `src/core/ghosts.cpp`, `src/core/ghosts.hpp`, `src/core/cell_system/CellStructure.cpp`, `ParticleDecomposition.hpp`, `RegularDecomposition.{cpp,hpp}`, `AtomDecomposition.{cpp,hpp}`, `HybridDecomposition.{cpp,hpp}`

**Interfaces:**
- Removes: `ghost_communicator()`, `GhostCommunicator`/`GhostCommunication`, `GHOST_PREFETCH`/`GHOST_PSTSTORE`, `prepare_comm()`, `revert_comm_order()`, `assign_prefetches()`, `exchange_ghosts_comm()`/`collect_ghost_force_comm()`, and the `halo_plan()` default/`nullptr` branch in the façade.
- Preserves: `map_data_parts`, `GHOSTTRANS_*` (still the payload selector consumed by the packing layer).

- [ ] **Step 1:** Make `ParticleDecomposition::halo_plan()` **pure virtual** (all three decompositions now implement it). Remove `exchange_ghosts_comm()`/`collect_ghost_force_comm()` from the interface and all overrides. Remove the `else` legacy branch from every façade method (always use the engine).

- [ ] **Step 2:** Delete `ghost_communicator()` and the `GhostCommunicator`/`GhostCommunication` types + `GHOST_*` job/flag macros from `ghosts.{cpp,hpp}` (packing already lives in `particle_packing`). Delete `prepare_comm`, `revert_comm_order`, `assign_prefetches`, and the `m_exchange_ghosts_comm`/`m_collect_ghost_force_comm` members from `RegularDecomposition`; likewise the legacy comm build in Atom/Hybrid.

- [ ] **Step 3: Build everything and run the full unit-test label:**

```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm/build -j8 && \
ctest --test-dir /tikhome/weeber/es/.claude/worktrees/comm/build -L unit_tests --output-on-failure
```
Expected: all PASS; no references to the deleted symbols remain (`grep -rn "ghost_communicator\|GHOST_PREFETCH\|prepare_comm" src/core` returns nothing but comments/history).

- [ ] **Step 4: Run the broad Python suite on 2 and 4 ranks** (regression sweep):

```bash
cd /tikhome/weeber/es/.claude/worktrees/comm/build
ctest -L python --output-on-failure -j4 2>&1 | tail -40
```
Expected: PASS (investigate any failure inline; do not proceed with failures).

- [ ] **Step 5: Final A/B** on an idle machine; confirm R7 (no regression) and record the ghost-region deltas (`ghost/wait` becomes the headline metric for later latency-hiding phases).

- [ ] **Step 6:** Format all changed files, commit, push to `me`, watch CI:

```bash
maintainer/format/format.sh $(git -C /tikhome/weeber/es/.claude/worktrees/comm diff --name-only HEAD~1) || true
git -C /tikhome/weeber/es/.claude/worktrees/comm add -A
git -C /tikhome/weeber/es/.claude/worktrees/comm commit -m "core/ghosts: remove staggered engine; async HaloPlan is the only path"
git -C /tikhome/weeber/es/.claude/worktrees/comm push me HEAD:async-ghost-comm
```

---

## Subsequent plans (authored after Ph1 lands, when interfaces are concrete)

These are separate plans (each producing working, testable software) per the writing-plans scope guidance; they depend on Ph1's finalized `HaloPlan`/engine types:

- **Ph2 — Unified validator.** Build-time + `ADDITIONAL_CHECKS` invariants across all decompositions: coverage (done as a lite check in 1.4), **neighborship match** (filled ghosts ⊇ ghosts referenced by `init_cell_interactions`), **cross-rank symmetric handshake** (`(peer, counts, cell-id hash)` exchange), reverse=mirror, interior/boundary consistency, op sanity. TDD with deliberately broken plans.
- **Ph3 — Interior/boundary tags + LE cleanup.** First-class interior/boundary cell classification; opt-in tag-based Euclidean-vs-MI pair-loop branch (measured); precise LE-sheared neighbor construction and **removal of `fully_connected_boundary`** (+ offset-drift rebuild/widen), gated by LE tests + validator + shear-boundary ghost-count assertion.
- **Ph4 — Latency hiding.** Quantify message-level overlap via `ghost/wait`; prototype one bounded compute-overlap case behind `start/finish`; keep only if measured worthwhile.
- **Ph5 — P3M CA-halo.** Port `p3m_send_mesh` onto the engine via a mesh-block pack/unpack policy (payload-generic engine), A/B on `p3m.py`.

---

## Self-Review (author checklist — completed)

**1. Spec coverage (Ph0+Ph1 slice):** Instrumentation (spec §6) → Tasks 0.1–0.2. Baseline/A-B methodology (§6) → Task 0.3. Data model (§5.1) → Task 1.1. Async engine + split-phase (§5.2) → Task 1.3. Construction/direct-neighbor (§5.3) → Task 1.4. Façade routing + parity + R5/R6 guards (§4) → Task 1.5. Collectives (§5.5) → Task 1.6. Removal of relay/PREFETCH/static buffers (§2, §9) → Task 1.7. Validator (§5.4), interior/boundary + Euclidean (§5.6), LE `fully_connected` removal (§5.3), latency-hiding (§5.2/Ph4), P3M (§3/Ph5) → explicitly deferred to subsequent plans above. No in-slice spec requirement is unassigned.

**2. Placeholder scan:** No "TBD/handle edge cases/similar to Task N". The engine `.cpp` in 1.3 gives real algorithm + exact boost::mpi calls with named helper contracts; helper signatures are defined in 1.2's Interfaces block. The one deliberate deferral (bond-message folding, PARTNUM sizing) is stated with the legacy line refs to copy from.

**3. Type consistency:** `HaloPlan`, `NeighborComm`, `SendRegion`, `LocalComm`, `CollectiveSection`, `ExchangeOp{Direction,Combine}` are defined once in 1.1 and used verbatim in 1.3–1.7. `halo_plan()` returns `HaloPlan const *` in both the interface (1.4) and the façade branch (1.5). Engine entry points `halo_exchange[_start/_finish]` are consistent across 1.3, 1.5, 1.6. `map_data_parts`/`GHOSTTRANS_*` retained through 1.7.
