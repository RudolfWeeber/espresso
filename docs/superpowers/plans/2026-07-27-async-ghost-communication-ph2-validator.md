# Async Ghost Communication — Ph2 (Unified Validator) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add a single, decomposition-agnostic validator for `GhostComm::HaloPlan` that proves every halo is communicated and that cell neighborships match the halo — satisfying spec goal #3.

**Architecture:** A pure, always-compiled function `GhostComm::validate_halo_plan(...)` returns a list of violations (empty = valid); a companion `validate_halo_plan_symmetry(...)` does a cross-rank MPI handshake. Production call sites invoke them under `#ifdef ESPRESSO_ADDITIONAL_CHECKS`; unit tests call them directly (testable in any build). This consolidates the ad-hoc asserts added in Ph1 (`RegularDecomposition.cpp:725` coverage, `HaloExchange.cpp:237` peer-uniqueness).

**Tech Stack:** C++20, Boost.Test (`espresso_unit_test`, `NUM_PROC`), Boost.MPI, the Ph1 `HaloPlan`/`HaloExchange` engine.

## Global Constraints

- Validator is a **pure function returning violations** (`std::vector<std::string>`), always compiled; only its production call sites are `#ifdef ESPRESSO_ADDITIONAL_CHECKS`. This keeps it unit-testable in the Release/maxset build (which has checks OFF).
- No behavior change to the exchange hot path. No new per-step cost in non-checked builds.
- Commit on branch `worktree-comm`; **do NOT create a new branch**. No push (controller pushes + watches CI). Format with `sh maintainer/format/clang-format.sh -i <files>`.
- Build is pre-built at `build/`; incremental `make -C build -j8 <target>` (auto-reconfigures on CMakeLists change); never re-run cmake manually.
- Verify unit tests with `make -C build -j8 check_unit_tests` (NOT `ctest -L unit_tests` — that label doesn't exist and returns a false green). Single test: `ctest --test-dir build -R <name> --output-on-failure`.
- `git -C /tikhome/weeber/es/.claude/worktrees/comm ...`; never `cd`-chain git.

## Reference (verified current signatures)

- `GhostComm::HaloPlan { boost::mpi::communicator comm; std::vector<NeighborComm> neighbors; std::vector<LocalComm> local; std::optional<CollectiveSection> collective; }`; `NeighborComm{int peer; vector<SendRegion> send; vector<ParticleList*> recv;}`; `LocalComm{ParticleList* src,*dst; Utils::Vector3d shift;}` — `src/core/ghosts/HaloPlan.hpp`.
- `ParticleDecomposition::halo_plan() const` pure virtual (`ParticleDecomposition.hpp:82`); `local_cells()`/`ghost_cells()` return `std::span<Cell *const>` (`ParticleDecomposition.hpp:93,104`).
- `Cell::neighbors()` → `Neighbors<Cell*>` with `.all()` (`cell_system/Cell.hpp:114,69`); a `Cell*`'s particle list is `Cell::particles()` (`ParticleList&`). The plan stores `ParticleList*`; map a `Cell*` to its `ParticleList*` via `&cell->particles()`.
- `RegularDecomposition::make_halo_plan()` (`RegularDecomposition.cpp:570`) has an `#ifdef ESPRESSO_ADDITIONAL_CHECKS` coverage block at `:725`; `m_halo_plan = make_halo_plan();` at `:768`. Atom (`AtomDecomposition.cpp:57`,`:86`), Hybrid (`HybridDecomposition.cpp:82,85`).
- Engine peer-uniqueness assert: `HaloExchange.cpp:237-245`.

## File Structure

- **Create** `src/core/ghosts/HaloPlanValidator.hpp` / `.cpp` — the two validator functions.
- **Create** `src/core/unit_tests/HaloPlanValidator_test.cpp` — Boost.Test (local checks, single-rank) + a `NUM_PROC 4` symmetry case.
- **Modify** `src/core/CMakeLists.txt` (add the new .cpp), `src/core/unit_tests/CMakeLists.txt` (register the test).
- **Modify** `src/core/cell_system/RegularDecomposition.cpp`, `AtomDecomposition.cpp`, `HybridDecomposition.cpp` — replace/route the ad-hoc coverage checks to call the validator under `ADDITIONAL_CHECKS`.
- **Modify** `src/core/ghosts/HaloExchange.cpp` — replace the inline peer-uniqueness assert with a call to the validator's peer-uniqueness (or leave as belt-and-suspenders; see Task 2.3), and add the op-sanity assert.

---

### Task 2.1: Local validator (coverage, neighborship-match, peer-uniqueness)

**Files:** Create `src/core/ghosts/HaloPlanValidator.{hpp,cpp}`; Create `src/core/unit_tests/HaloPlanValidator_test.cpp`; Modify `src/core/CMakeLists.txt`, `src/core/unit_tests/CMakeLists.txt`.

**Interfaces produced (Tasks 2.2/2.3 consume):**
- `std::vector<std::string> GhostComm::validate_halo_plan(HaloPlan const &plan, std::span<Cell *const> local_cells, std::span<Cell *const> ghost_cells);`
  - **Coverage:** every `ParticleList*` in `ghost_cells` (as `&c->particles()`) appears as a `recv` entry (across all `neighbors[i].recv`) or a `local[j].dst` **exactly once**; no `recv`/`dst` target is outside `ghost_cells`; no target appears twice.
  - **Neighborship-match:** for every local cell, every neighbor (`cell->neighbors().all()`) whose `ParticleList*` is a ghost (in `ghost_cells`) must be a covered `recv`/`dst` target. (Report any referenced-but-uncommunicated ghost.)
  - **Peer-uniqueness:** each `peer` appears in at most one `NeighborComm`.
  - **Send/recv shape:** each `NeighborComm` has `send.size() == recv.size()`.
  - Returns a human-readable violation string per problem; empty vector = valid.

- [ ] **Step 1: failing test** `HaloPlanValidator_test.cpp` — build small hand-made cells + a valid plan, assert `validate_halo_plan(...)` returns empty; then mutate copies to inject each defect and assert the specific violation is reported:

```cpp
#define BOOST_TEST_MODULE HaloPlanValidator test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/HaloPlan.hpp"
#include "cell_system/Cell.hpp"
#include <span>
#include <vector>

BOOST_AUTO_TEST_CASE(detects_defects) {
  using namespace GhostComm;
  // 1 local cell whose only ghost neighbor is `ghost`.
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});
  std::vector<Cell *> locals{&local}, ghosts{&ghost};

  HaloPlan good;
  good.neighbors.push_back(NeighborComm{1, {SendRegion{&local.particles(), {}}},
                                        {&ghost.particles()}});
  BOOST_CHECK(validate_halo_plan(good, locals, ghosts).empty());

  // uncovered ghost (neighborship-match + coverage failure)
  HaloPlan uncovered; // empty plan, ghost never filled
  BOOST_CHECK(!validate_halo_plan(uncovered, locals, ghosts).empty());

  // duplicate peer
  HaloPlan duppeer = good;
  duppeer.neighbors.push_back(NeighborComm{1, {SendRegion{&local.particles(), {}}},
                                           {&ghost.particles()}});
  BOOST_CHECK(!validate_halo_plan(duppeer, locals, ghosts).empty());

  // double-filled ghost
  HaloPlan dbl = good;
  dbl.neighbors.push_back(NeighborComm{2, {SendRegion{&local.particles(), {}}},
                                       {&ghost.particles()}});
  BOOST_CHECK(!validate_halo_plan(dbl, locals, ghosts).empty());

  // send/recv size mismatch
  HaloPlan badshape = good;
  badshape.neighbors[0].recv.clear();
  BOOST_CHECK(!validate_halo_plan(badshape, locals, ghosts).empty());
}
```

- [ ] **Step 2: register + confirm RED**

```cmake
espresso_unit_test(SRC HaloPlanValidator_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi MPI::MPI_CXX)
```
Run `cmake --build build --target HaloPlanValidator_test 2>&1 | tail -20` → FAIL (no header).

- [ ] **Step 3: implement `HaloPlanValidator.{hpp,cpp}`** — GPL header; the function per the Interfaces contract. Build a `std::unordered_set<ParticleList const*>` of ghost lists (`&c->particles()`); walk `plan.neighbors[*].recv` + `plan.local[*].dst` accumulating a `std::unordered_map<ParticleList const*,int>` fill-count; emit violations for count!=1, targets not in the ghost set, duplicate peers, and per-neighbor `send.size()!=recv.size()`. For neighborship-match, iterate `local_cells`, and for each `cell->neighbors().all()` entry whose `&n->particles()` is in the ghost set, verify it's a covered target. Add both new sources to `src/core/CMakeLists.txt` next to `ghosts/HaloExchange.cpp`.

- [ ] **Step 4: confirm GREEN** — `make -C build -j8 HaloPlanValidator_test && ctest --test-dir build -R HaloPlanValidator_test --output-on-failure` → PASS.

- [ ] **Step 5: format + commit** — `sh maintainer/format/clang-format.sh -i src/core/ghosts/HaloPlanValidator.{hpp,cpp} src/core/unit_tests/HaloPlanValidator_test.cpp`; commit `core/ghosts: add local HaloPlan validator (coverage, neighborship-match, peer-uniqueness)`.

---

### Task 2.2: Cross-rank symmetric handshake

**Files:** Modify `HaloPlanValidator.{hpp,cpp}`; Modify `HaloPlanValidator_test.cpp`.

**Interfaces produced:** `std::vector<std::string> GhostComm::validate_halo_plan_symmetry(HaloPlan const &plan);` — a lightweight MPI exchange: for each `NeighborComm nc`, send `nc.send.size()` to `nc.peer` and receive the count `peer` intends to `recv` from me; assert my `send`-count to P equals P's `recv`-count from me (and vice-versa). Uses non-blocking irecv/isend + wait_all over `plan.comm` (mirror the engine's deadlock-free pattern). Returns violations (empty = valid).

- [ ] **Step 1: failing multi-rank test** — extend `HaloPlanValidator_test.cpp` with a `NUM_PROC 4` case that builds a real `RegularDecomposition` (pattern: `Verlet_list_test.cpp`) and asserts `validate_halo_plan_symmetry(*dd.halo_plan()).empty()`; then a synthetic asymmetric plan (a `NeighborComm` to a peer with a `send` count the peer won't mirror) asserts a non-empty result. Register the test target with `NUM_PROC 4` (update the CMake line: `... NUM_PROC 4`).

- [ ] **Step 2: confirm RED** (function missing) — `cmake --build build --target HaloPlanValidator_test 2>&1 | tail`.

- [ ] **Step 3: implement `validate_halo_plan_symmetry`** — post an `irecv` of an `int` (the peer's send-count toward me) from each `nc.peer`, `isend` my `nc.send.size()` to each `nc.peer`, `boost::mpi::wait_all`; compare each received count to my corresponding `nc.recv.size()` (my recv from P must equal what P sends me). Report mismatches. (Collective section: skip — n-square symmetry is inherent to bcast/reduce.)

- [ ] **Step 4: confirm GREEN** (4 ranks) — `make -C build -j8 HaloPlanValidator_test && ctest --test-dir build -R HaloPlanValidator_test --output-on-failure`.

- [ ] **Step 5: format + commit** — `core/ghosts: add cross-rank HaloPlan symmetry check`.

---

### Task 2.3: Wire the validator in; consolidate Ph1 ad-hoc asserts; op-sanity

**Files:** Modify `RegularDecomposition.cpp`, `AtomDecomposition.cpp`, `HybridDecomposition.cpp`, `HaloExchange.cpp`.

**Interfaces consumed:** `validate_halo_plan`, `validate_halo_plan_symmetry` (2.1/2.2).

- [ ] **Step 1:** In each decomposition, right after `m_halo_plan = make_halo_plan();`, add (guarded):
```cpp
#ifdef ESPRESSO_ADDITIONAL_CHECKS
  {
    auto const v = GhostComm::validate_halo_plan(m_halo_plan, local_cells(), ghost_cells());
    assert(v.empty());
    auto const s = GhostComm::validate_halo_plan_symmetry(m_halo_plan);
    assert(s.empty());
  }
#endif
```
Replace the existing bespoke coverage block at `RegularDecomposition.cpp:725` with this call (delete the now-redundant hand-rolled coverage loop). For Atom (collective-only plan: `neighbors`/`local` empty) `validate_halo_plan` trivially passes; keep the call for uniformity.

- [ ] **Step 2:** In `HaloExchange.cpp`: keep the existing peer-uniqueness assert OR replace it with a comment pointing at the validator (the validator now covers peer-uniqueness at plan-build; the engine assert is redundant — remove it to avoid duplication, or keep as a cheap guard. Choose: **remove** it and rely on the plan-build validator, since every plan is validated at build under ADDITIONAL_CHECKS). Add an **op-sanity** assert at the top of `halo_exchange_start`: `assert(op.combine != Combine::Add || data_parts == GHOSTTRANS_FORCE || data_parts == GHOSTTRANS_RATTLE)` under `#ifdef ESPRESSO_ADDITIONAL_CHECKS`.

- [ ] **Step 2b:** `git grep` to confirm no dangling reference to the deleted coverage/uniqueness code.

- [ ] **Step 3: verify** — `make -C build -j8 check_unit_tests` (all pass). Then a checks-ON spot check: configure a **separate** build dir with `-D CMAKE_BUILD_TYPE=RelWithAssert` is heavyweight; instead confirm the validator is exercised by running `HaloPlanValidator_test` + `RegularHaloPlan_test` + `HaloExchange_test` (they call the validator directly / build plans). Run a 2-rank + 4-rank `mpiexec ./build/pypresso testsuite/python/{lj,lees_edwards,nsquare,hybrid_decomposition}.py` to confirm no behavior change and (in a checks-on CI build) the asserts don't trip on valid plans.

- [ ] **Step 4: format + commit** — `cell_system: validate HaloPlan at build under ADDITIONAL_CHECKS; drop ad-hoc Ph1 asserts`.

---

## Self-Review (author checklist)
- Spec §5.4 invariants: coverage (2.1) ✓, neighborship-match (2.1) ✓, cross-rank symmetry (2.2) ✓, peer-uniqueness (2.1) ✓, op-sanity (2.3) ✓. Reverse=mirror is structural (reduce is a runtime view of the same neighbors — no separate object to drift) — documented, not a runtime check. Interior/boundary consistency is deferred to Ph3 (tags don't exist yet).
- No placeholders; every step has commands/code. Types match Ph1 (`HaloPlan`, `NeighborComm`, `validate_halo_plan`) across tasks.
- Testable in Release (validator is a plain function; only call sites are `#ifdef`).
