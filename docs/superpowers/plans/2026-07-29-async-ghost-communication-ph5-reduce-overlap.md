# Async Ghost Communication — Ph5 (Force-Reduce / Integrator-Step-2 Overlap) Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Hide the force-reduce wait (the largest ghost cost at 10k/8: ~2.9% avg, imbalance-dominated) behind the second velocity-Verlet half-kick: `ghosts_reduce_forces_start` → `integrator_step_2` over **interior** cells → `finish` → step 2 over **boundary** cells. Plus two enablers: local transfers moved into `halo_exchange_start` (cheap latency hiding for all split-phase use), and a **robust interior definition** that also excludes periodic-wrap neighborships (user requirement: on 1 MPI rank, *interior* means the cell's neighborship does not wrap around the periodic boundary, even when the wrap partner is not a ghost cell).

**Why step 2 and not the pair loop:** step 2 is a per-particle sweep over the cell lists (`for_each_local_particle`) — no pair structure, no Verlet list, no AoSoA. Splitting by `Cell::is_boundary()` reorders only independent per-particle updates → bitwise-identical trajectories. Interior cells receive no reduce contributions by construction (reduce destinations = exported cells = boundary).

## Global Constraints
- Verification build `/ssd/weeber/comm-build` (RelWithAssert + -Werror + checks): `make -C /ssd/weeber/comm-build -j8 check_unit_tests` = **149/149**; `ctest --test-dir /ssd/weeber/comm-build -R caliper` green. Timing build `/ssd/weeber/comm-build-nocali` (Release, caliper OFF) — controller A/Bs; **machine-idle gate** before timing. Home quota FULL — build only on /ssd.
- Commit on `worktree-comm`; NO new branch; NO push (controller pushes + CI). Format `sh maintainer/format/clang-format.sh -i`. `git -C ...`.
- Physics bitwise identical on the default path; Caliper label SET unchanged (`testsuite/python/caliper.py` untouched).
- Parity gate per task: lj@2+4, lees_edwards@4, collision_detection@4, nsquare@4, hybrid_decomposition@4 on the /ssd RelWithAssert build.

## Reference (verified)
- Engine: `src/core/ghosts/HaloExchange.cpp` — `halo_exchange_start` (resolve cells → pack → irecv/isend post), `halo_exchange_finish` (GHOSTTRANS_NONE early-return → **local copies** → collective → progressive wait/unpack → send wait). `ExchangeBuffers` pool on `CellStructure::m_ghost_buffers`; `GhostExchange` is move-only (owned unique_ptr).
- Post-reduce consumers (`src/core/forces.cpp:487ff`): only `comfixed->apply`, `force_capping`, `recalc_forces=false`. Loop tail (`src/core/integrate.cpp:716-741`): LB-tracer arm (LB only) → `integrator_step_2` (integrate.cpp:476, per-particle via `for_each_local_particle`, NPT arm after the loop) → BD resort / LE offset / SHAKE / LB-EK propagation.
- Tags: `Cell::m_is_boundary` set by `GhostComm::mark_boundary_cells(local, ghost)` (`src/core/ghosts/mark_boundary_cells.hpp`), called from all 3 decompositions; validator interior/boundary consistency in `src/core/ghosts/HaloPlanValidator.cpp`.

---

### Task 5.1: Local transfers in `halo_exchange_start`

**Files:** `src/core/ghosts/HaloExchange.cpp` (+ its unit test if it asserts phase placement).

- [ ] Move the `local_cell_copy` loop from `halo_exchange_finish` to the END of `halo_exchange_start` (pool overload; after the isend posts, outside the Caliper `ghost/post` region — no new labels). Rationale: local copies then overlap in-flight messages for split-phase callers and shrink `finish`. The **collective section stays in `finish`** (blocking bcast/reduce would stall `start`).
- [ ] Correctness argument to encode as a comment: local copies never touch send buffers (Push: they write ghost cells, sends carry real-cell data; Reduce: they write real cells, sends carry ghost-cell data), and they still execute before any unpack → identical FP accumulation order → bitwise-identical results.
- [ ] GHOSTTRANS_NONE early-return in `start` must still skip them; `finish` keeps its NONE early-return.
- [ ] Verify: check_unit_tests 149/149; parity gate; caliper ctest. Commit `core/ghosts: run local transfers in halo_exchange_start (overlap with in-flight messages)`.

### Task 5.2: Robust interior classification (periodic-wrap rule) + validator

**Files:** `src/core/ghosts/mark_boundary_cells.hpp`, `src/core/cell_system/{Regular,Atom,Hybrid}Decomposition.cpp`, `src/core/ghosts/HaloPlanValidator.{hpp,cpp}`, `src/core/unit_tests/Cell_boundary_test.cpp`.

**Rule:** a local cell is **boundary** iff (a) any of `neighbors().all()` is a ghost cell **[existing]**, OR (b) any of its neighbor relations crosses a periodic box boundary — even when both cells are local (wrap without ghosts). Interior = neither.

- [ ] Extend `GhostComm::mark_boundary_cells` with an optional per-cell-pair predicate `wraps(Cell const *a, Cell const *b) -> bool` (default: never) OR an extra per-cell predicate — pick the cleanest fit; the decomposition supplies it.
- [ ] `RegularDecomposition`: supply the wrap predicate from grid knowledge — a neighbor relation wraps along axis `i` iff the axis is periodic, the local domain spans the whole box along `i` (node_grid[i]==1), and the two cells are in the first/last cell layer of that axis such that their adjacency is the folded one. (Equivalent, acceptable simplification: with node_grid[i]==1 and periodic `i`, mark ALL cells in the first and last local layer along `i` as boundary — conservative and simple; note which variant was implemented.)
- [ ] `AtomDecomposition`: all local cells boundary (no spatial locality — interior is empty; overlap degenerates gracefully to no-op interior pass). `HybridDecomposition`: conservative — all boundary (document).
- [ ] Validator: add the **overlap-safety invariant** — an interior cell must appear in NO `NeighborComm` send/recv region and NO `LocalComm` src/dst of the plan. This is the exact precondition Task 5.3 relies on; wire it into the existing `validate_halo_plan` consistency section.
- [ ] Unit test: hand-made cells — wrap-neighbor (local↔local across the boundary) ⇒ boundary even without ghost neighbors; non-wrap fully-local ⇒ interior; 1-rank RegularDecomposition (via existing real-decomp test pattern): every cell whose neighborship wraps is boundary.
- [ ] Verify: 149/149 + parity gate (validator runs under RelWithAssert on real decomps incl. 1 rank). Commit `cell_system: boundary tag accounts for periodic-wrap neighborships; validator overlap invariant`.

### Task 5.3: Split-phase force reduce + interior/boundary step-2 overlap + A/B

**Files:** `src/core/cell_system/CellStructure.{hpp,cpp}` (split-phase façade + pending state), `src/core/forces.cpp` (start instead of blocking reduce when eligible), `src/core/integrate.cpp` (interior pass → finish → boundary pass), tests.

**Interfaces:** `CellStructure::ghosts_reduce_forces_start()` (stores pending `GhostExchange` in the CellStructure; asserts none pending), `CellStructure::ghosts_reduce_forces_finish()` (finishes + clears; safe to call only when pending), `CellStructure::has_pending_ghost_reduce()`. Pending exchange MUST be resolved before any resort/decomposition change (assert in resort path + destructor).

- [ ] **Eligibility (conservative whitelist, computed in `calculate_forces`):** overlap only when ALL hold: **`comm.size() > 1`** (at 1 MPI rank the reduce is local-only — nothing to hide, skip the split entirely; user decision), `comfixed` inactive, `force_cap == 0`, integ method is VV or symplectic Euler (NOT steepest descent), NPT propagation not in use, LB-tracer arm not active (no `TRANS_LB_TRACER` with LB thermostat). Everything else takes today's blocking path unchanged. Encode as a small function with a comment per condition.
- [ ] `calculate_forces`: when eligible, call `ghosts_reduce_forces_start()` instead of the blocking reduce and skip comfixed/capping (they are inactive by eligibility — assert); when not eligible, exactly today's code.
- [ ] `integrate.cpp` after `calculate_forces`: if `has_pending_ghost_reduce()`: `integrator_step_2` restricted to **interior** cells → `ghosts_reduce_forces_finish()` → step 2 restricted to **boundary** cells. Else: unchanged full step 2. Implement the restriction with a cell-filtered helper (iterate `local_cells()`, skip by `is_boundary()` predicate — factor the existing per-particle lambda so both passes share it verbatim; the NPT arm inside `integrator_step_2` stays on the ineligible path only).
- [ ] Caliper: reuse the existing `ghosts_reduce_forces` label around BOTH start and finish episodes (same label twice per step aggregates; label set unchanged). `integrator_step_2` keeps its function marker (fires per pass — same label).
- [ ] **Safety nets:** `ADDITIONAL_CHECKS` assert after finish that no interior cell was a reduce destination (or rely on 5.2's validator invariant — state which); assert in `resort_particles`/decomposition-set path that nothing is pending.
- [ ] Verify: 149/149; parity gate (these run the eligible path for lj/nsquare, ineligible variants covered by collision/lees_edwards configs — confirm both paths execute by temporarily logging in a scratch run, not in the commit); bitwise check: `lj.py` forces/positions identical with overlap on vs off on 4 ranks (flip eligibility artificially in a scratch build for the check — document).
- [ ] Commit `core: overlap ghost force reduction with interior-cell velocity update`.
- [ ] **A/B (controller-led, idle-gated, caliper-off builds):** lj @ 1k(1/8) + 10k(4/8) vs pre-overlap HEAD binary and vs python-staggered reference. Record deltas; expectation: 10k/8 reduce-wait (~2.9%) partially hidden; NO regression at 1 rank / 1k. If no gain or regression → document negative result, keep 5.1+5.2, revert 5.3's forces/integrate changes (engine façade may stay).

---

## Self-Review (author)
- 5.1 is behavior-preserving for blocking callers (same op order) and pure win for split-phase; collective stays blocking in finish — documented.
- 5.2 implements the user-mandated wrap rule as first-class infrastructure and adds the exact validator invariant 5.3 depends on; conservative choices (Atom/Hybrid all-boundary) degrade gracefully.
- 5.3 is whitelist-gated with today's path as fallback, per-particle reorder only (bitwise identity), pending-exchange lifecycle asserted, and value-gated by A/B with an explicit revert clause.

### Task 5.4: Non-zeroing CommBuf resize (perf-mandated)

**Files:** `src/core/ghosts/particle_packing.hpp` (CommBuf), engine call sites if signatures shift.

perf @10k/8 (cycles:u, async vs staggered): `__memset_avx2` +~1.0% of total cycles — `CommBuf`'s `std::vector<char>::resize` value-initializes (memset) every grown buffer before `pack_cells` immediately overwrites it; recv buffers are likewise overwritten by MPI. Largest single fixable delta.

- [ ] Replace the zero-filling grow with an uninitialized-capacity scheme: keep `std::vector<char>` but grow via `reserve` + track a separate logical size, or switch the storage to a small buffer type with `malloc`-style uninitialized grow (pick the least invasive that keeps `data()/size()/bonds()` semantics; document the choice). Bond buffer (`std::vector<char>` serialized via boost.mpi) stays as-is (cold path).
- [ ] No behavior change: buffers are fully written by pack (send) or by MPI receive before any read — state this invariant in a comment; `ADDITIONAL_CHECKS`-guarded canary optional.
- [ ] Verify: 149/149; parity gate; A/B is covered by the Ph5 combined A/B.
- [ ] Commit `core/ghosts: grow exchange buffers without zero-fill (perf: -1% memset at 10k/8)`.

### Task 5.5: Runtime-conditional orientation payload (GHOSTTRANS_QUAT / GHOSTTRANS_TORQUE)

**Files:** `src/core/ghosts.hpp` (bits), `src/core/ghosts/particle_packing.cpp` (serialize_and_reduce keys quat/torque on the new bits), `src/core/system/System.cpp` (`get_global_ghost_flags`), `src/core/cell_system/CellStructure.cpp` (`map_data_parts`, ghosts_update/reduce paths), tests.

**Finding:** with ESPRESSO_ROTATION compiled in, every per-step position push carries `quat` (32 of 68 B/particle) and every force reduce carries `torque` (24 of 48 B/particle) keyed only on the compile-time feature — pure waste for non-rotating systems (LJ). Runtime-conditional payload halves serialization bytes in the common case; the architectural pattern already exists (MOMENTUM/BONDS bits in `get_global_ghost_flags`).

- [ ] New bits `GHOSTTRANS_QUAT`, `GHOSTTRANS_TORQUE` (under `#ifdef ESPRESSO_ROTATION`). `serialize_and_reduce`: quat moves from the POSITION arm to a QUAT arm; torque from the FORCE arm to a TORQUE arm. `calc_transmit_size` follows automatically (same function).
- [ ] **Reader audit (correctness gate):** grep every reader of `p.quat()`, `calc_director()`, `p.torque()`, `p.calc_dip()` reachable for GHOST particles (nonbonded anisotropic kernels e.g. Gay-Berne, dipolar solvers, vs_relative update/back-transfer, engine/swimming, collision detection variants, LB coupling). Derive the whitelist in `get_global_ghost_flags`: set QUAT (and TORQUE for the reduce) when ANY of: rotational propagation in use (`used_propagations & ROT_*`), dipolar solver set, vs_relative in use, anisotropic nonbonded interaction configured, engine/swimming with director coupling — document each condition next to the code like the MOMENTUM ones. When in doubt for a reader, include it (conservative).
- [ ] `ghosts_reduce_forces` (blocking + split-phase): TORQUE bit computed from the same whitelist (thread through or recompute — keep both paths consistent).
- [ ] Wire symmetry: pack and unpack use the same `data_parts` per exchange (existing invariant) → layout stays symmetric automatically; state this in a comment.
- [ ] Tests: rotation parity `testsuite/python/rotation_per_particle.py` or an existing rotating-system test @2+4 ranks (quat must still propagate when rotation used); dipolar test if cheap (`dawaanr-and-dds-gpu`? pick an existing CPU dipolar test); plus the standard parity gate (lj/lees_edwards/collision/nsquare/hybrid). `check_unit_tests` all green.
- [ ] A/B (controller-led): lj 1k/8 + 10k/4+8 — expect a visible cut of the 10k/8 residual.
- [ ] Commit `core/ghosts: send quat/torque only when orientation physics is active`.
