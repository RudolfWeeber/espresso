# Async Ghost Communication — Ph3 (Interior/Boundary Tags + Euclidean-vs-MI) Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Give each local cell a first-class interior/boundary classification (foundation for Ph4 compute/comm overlap), fold it into the validator, and — as a *measured experiment* — use it to compute plain Euclidean distance for interior pairs instead of minimum-image (which the regular decomposition currently pays on every pair).

**Architecture:** A `Cell::m_is_boundary` flag, set by each decomposition after neighbor setup (a local cell is *boundary* iff any of its `neighbors().all()` is a ghost cell; else *interior*). Interior cells can only pair with local cells, so no periodic-image ghost is ever involved and plain Euclidean distance is exact. The validator gains an interior/boundary consistency check. The Cabana pair kernel optionally uses Euclidean for interior↔interior pairs; kept only if the LJ A/B shows a real gain.

## Global Constraints
- Verification build: **`/ssd/weeber/comm-build`** (`RelWithAssert` + `-Werror` + maxset, checks-on). Build/test there: `make -C /ssd/weeber/comm-build -j8 <target>`; `make -C /ssd/weeber/comm-build -j8 check_unit_tests` (NOT `ctest -L unit_tests`). Home quota is full — do NOT build on home.
- For the Ph3 T3.2 **A/B timing** a separate **Release** build on /ssd is needed (RelWithAssert distorts timing) — the controller sets this up; implementers of T3.2 measure with it.
- Commit on `worktree-comm`; NO new branch; NO push (controller pushes + watches CI). Format `sh maintainer/format/clang-format.sh -i`. `-Werror` is on in the /ssd build → fix warnings locally. `git -C /tikhome/weeber/es/.claude/worktrees/comm ...`.
- **LE `fully_connected_boundary` removal is OUT OF SCOPE** (Ph1 review established it is a Verlet-neighbor-list construct, not a ghost-comm workaround; tangential to this phase).
- R5 preserved: ghost positions stay unshifted; the Euclidean choice is a per-pair distance-function selection, not a change to stored positions.

## Reference (verified)
- `Cell` (`src/core/cell_system/Cell.hpp:96-115`): members `m_particles`, `m_neighbors`, `m_verlet_list`; `neighbors()` → `Neighbors<Cell*>` with `.all()`. Add the flag here.
- Decompositions set up neighbors in `init_cell_interactions` (RegularDecomposition.cpp:477+) / `configure_neighbors` (Atom). `local_cells()`/`ghost_cells()` return `std::span<Cell *const>`.
- Cabana pair kernel: `link_cell_kokkos` (`src/core/short_range_cabana.hpp:92-164`) computes `box_geo.get_mi_dist2(p1.pos(), p2.pos())` for the Verlet-list criterion (intra: line 114; inter: line 139). The force distance is recomputed in `nonbonded_kernel` (from `forces_cabana.hpp`).
- Validator: `src/core/ghosts/HaloPlanValidator.cpp` `validate_halo_plan(...)`.

---

### Task 3.1: Interior/boundary cell classification + validator consistency

**Files:** `src/core/cell_system/Cell.hpp` (add flag); `RegularDecomposition.cpp`, `AtomDecomposition.cpp`, `HybridDecomposition.cpp` (set flag); `src/core/ghosts/HaloPlanValidator.{hpp,cpp}` (+ consistency check); tests `src/core/unit_tests/Cell_boundary_test.cpp` (or extend an existing cell test) + a case in `HaloPlanValidator_test.cpp`.

**Interfaces produced:** `bool Cell::is_boundary() const` (and the stored `m_is_boundary`, default `false`); each decomposition sets it for local cells; the validator gains an interior/boundary consistency check (an interior local cell must have NO ghost neighbor).

- [ ] **Step 1 (failing test):** In a Boost.Test, build hand-made cells: cell A with only local neighbors → expect interior; cell B with a ghost neighbor → expect boundary; assert a small helper `mark_boundary_cells(local_cells, ghost_cells)` sets the flags correctly. (If a free helper is cleanest, declare `void GhostComm::mark_boundary_cells(std::span<Cell *const> local, std::span<Cell *const> ghost)` in a small header, or make it a decomposition detail — pick per code fit and note it.)
- [ ] **Step 2:** RED (`make -C /ssd/weeber/comm-build -j8 <test>` → fails).
- [ ] **Step 3 (implement):** Add `bool m_is_boundary = false;` + `bool is_boundary() const { return m_is_boundary; }` to `Cell`. Implement the marking: for each local cell, `m_is_boundary = any(neighbor in ghost_set for neighbor in neighbors().all())`. Call it from each decomposition right after neighbor setup + plan build (where `make_halo_plan` is called). Reset local cells to interior first (idempotent on rebuild).
- [ ] **Step 4 (validator):** Add to `validate_halo_plan` an interior/boundary consistency check: for each local cell with `is_boundary()==false`, assert none of its `neighbors().all()` is in `ghost_set` (else violation "interior cell has a ghost neighbor"). This composes with the existing checks.
- [ ] **Step 5:** GREEN — the unit test passes; and `make -C /ssd/weeber/comm-build -j8 check_unit_tests` = all pass (the wired validator now also checks boundary consistency on real decompositions under RelWithAssert — this is the guard that marking is correct on real plans, incl. 1-rank).
- [ ] **Step 6:** Format + commit `cell_system: classify interior/boundary cells; validator checks consistency`.

### Task 3.2: Euclidean-vs-MI for interior pairs (measured experiment)

**Files:** `src/core/short_range_cabana.hpp` (`link_cell_kokkos`); possibly `forces_cabana.hpp` (nonbonded distance) — only if the force kernel also needs the per-pair choice; tests + A/B.

**Interfaces consumed:** `Cell::is_boundary()` (3.1).

- [ ] **Step 1 (design + correctness):** In `link_cell_kokkos`, when BOTH cells of a pair are interior (`!cells[i]->is_boundary() && !neighbor->is_boundary()`), the pair cannot involve a periodic-image ghost, so replace `box_geo.get_mi_dist2(p1.pos(), p2.pos())` with the plain squared Euclidean distance `(p1.pos() - p2.pos()).norm2()`. Otherwise keep `get_mi_dist2`. Do the SAME in the force `nonbonded_kernel` if it recomputes the distance (verify; if it takes the Distance from the caller, thread the choice through). Add a `#ifdef ESPRESSO_ADDITIONAL_CHECKS` assert that for an interior-interior pair the Euclidean and MI distances agree (guards correctness).
- [ ] **Step 2 (correctness test):** A unit/Python parity check that forces/energies are unchanged vs the MI-only path (e.g. `lj.py`, `lees_edwards.py` still pass identically on 1/2/4 ranks on the /ssd build). This is the gate: physics MUST be identical.
- [ ] **Step 3 (A/B):** On a **Release** /ssd build (controller provides), A/B `lj.py` (and `p3m.py`) at 1/2/4/8 ranks: interior-Euclidean vs the MI-only baseline (same build, toggle via a temporary flag or two commits). Record per-step deltas + Caliper `cabana_pair_loop` region.
- [ ] **Step 4 (decision):** If a real gain (beyond CI): keep it, commit `core/short_range: Euclidean distance for interior cell pairs (A/B: <numbers>)`. If negligible/regressive or too invasive on the force kernel: **revert the kernel change, keep the tags** (they're still needed for Ph4), and document the negative result in the report + a project memory. Either way the interior/boundary tags from 3.1 remain.

---

## Self-Review (author)
- Tags (3.1) are the Ph4 foundation and are low-risk/clean; validator consistency guards correctness on real plans (incl. 1-rank) via the checks-on /ssd build.
- Euclidean-MI (3.2) is correctness-gated (physics identical) then value-gated (A/B); interior↔interior is provably wrap-free, so Euclidean is exact. Kept only if measured worthwhile.
- LE fully_connected explicitly out of scope (documented).
