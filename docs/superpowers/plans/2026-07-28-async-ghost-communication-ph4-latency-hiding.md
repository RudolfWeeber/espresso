# Async Ghost Communication — Ph4 (Latency Hiding / 8-rank mitigation) Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Mitigate the 8-rank async regression (LJ +9.3%, P3M +13.2% vs the staggered baseline; neutral/better at 1/2/4) by profiling the async exchange and applying the cause-indicated latency-hiding, using the interior/boundary tags + the split-phase `start/finish` API.

**Architecture:** Profile-then-fix. First instrument the engine's phases (`ghost/pack·post·wait·unpack`) with Caliper. Profile the 8-rank LJ/P3M case to find the dominant cost. Then apply the indicated mitigation — if `ghost/wait`-dominated, overlap interior force compute (interior cells, no ghost dependency) with the in-flight position exchange via `halo_exchange_start` → interior forces → `halo_exchange_finish` → boundary forces; if pack/post/buffer-dominated, reduce per-message overhead. A/B to confirm the 8-rank regression shrinks without regressing 1/2/4.

## Global Constraints
- Correctness build: **`/ssd/weeber/comm-build`** (`RelWithAssert` + `-Werror` + maxset, caliper on). Timing/profiling build: **`/ssd/weeber/comm-build-release`** (`Release`, no asserts, caliper on) — controller sets up; RelWithAssert distorts timing.
- **Machine-idle gate before every benchmark/profile** (`maintainer/benchmarks/ghost_ab.sh`'s gate; check load / `who`).
- Commit on `worktree-comm`; NO new branch; NO push (controller pushes + CI). Format `sh maintainer/format/clang-format.sh -i`. `-Werror` on locally. `git -C ...`.
- Physics MUST stay identical (parity: lj/lees_edwards/p3m/collision/nsquare/hybrid/lb-coupling on the /ssd RelWithAssert build) — overlap reorders *when* compute happens, never the result.

## Reference (verified)
- Engine: `src/core/ghosts/HaloExchange.cpp` `halo_exchange_start` (posts irecv, packs, posts isend) / `halo_exchange_finish` (local copies, collective, `wait_all`, unpack). Blocking `halo_exchange` = start+finish.
- Façade: `CellStructure::ghosts_update` etc. call blocking `halo_exchange`.
- Force pipeline: `integrate.cpp` → `calculate_forces` (`forces.cpp`) → `cabana_short_range` (`short_range_cabana.hpp`); ghost push before, `ghosts_reduce_forces` after (`forces.cpp:424`).
- Tags: `Cell::is_boundary()` (Ph3). Caliper idiom: `#ifdef ESPRESSO_CALIPER` + `CALI_MARK_BEGIN/END`. Caliper labels test: `testsuite/python/caliper.py` `EXPECTED_LABELS`.
- Baseline numbers: `docs/superpowers/baselines/`; async-@1.5 numbers in project memory (LJ 1/2/4/8 vs staggered).

---

### Task 4.1: Caliper phase markers in the engine

**Files:** `src/core/ghosts/HaloExchange.cpp`; `testsuite/python/caliper.py` (EXPECTED_LABELS).

- [ ] **Step 1:** Wrap the phases of `halo_exchange_start`/`finish` in `#ifdef ESPRESSO_CALIPER` `CALI_MARK_BEGIN/END`: `ghost/pack` (buffer packing), `ghost/post` (irecv+isend posting), `ghost/wait` (`wait_all`), `ghost/unpack` (unpack/reduce). Keep them cheap + guarded (zero cost when caliper off).
- [ ] **Step 2:** Build on `/ssd/weeber/comm-build` (caliper on); run a short caliper report to confirm the new regions appear; update `testsuite/python/caliper.py` `EXPECTED_LABELS` to include them at observed nesting (run the caliper child, match the tree — the same observe-then-match method as Ph0). `ctest --test-dir /ssd/weeber/comm-build -R caliper` green.
- [ ] **Step 3:** `check_unit_tests` green; format; commit `core/ghosts: Caliper phase markers (pack/post/wait/unpack) in halo_exchange`.

### Task 4.2: Profile the 8-rank case (analysis)

**Controller-led measurement (like the baseline).** On `/ssd/weeber/comm-build-release` (Release, caliper on), machine idle:
- [ ] Run `CALI_CONFIG=runtime-report mpiexec --bind-to core -n 8 ./pypresso ../maintainer/benchmarks/lj.py ...` (and n=4 for comparison; and p3m). Capture per-region times: `ghost/pack·post·wait·unpack`, `ghosts_update`, `ghosts_reduce_forces`, `cabana_pair_loop`.
- [ ] **Determine the dominant cost at 8 ranks** and write it to the report: is the regression `ghost/wait` (latency — overlap will help), `ghost/post`+`pack` (per-message overhead / count — coalescing/buffer-reuse), or `unpack`? This decides Task 4.3.

### Task 4.3: Apply the cause-indicated mitigation + A/B

**Files:** depend on 4.2's finding. Most likely `forces.cpp`/`integrate.cpp` (+ `CellStructure` split-phase wrappers) for overlap; or `HaloExchange.cpp`/`make_halo_plan` for message coalescing/buffer-reuse.

- [ ] **If `ghost/wait`-dominated (latency) → compute/comm overlap:** add split-phase façade wrappers (`ghosts_update_start`/`_finish` returning/taking a handle) over `halo_exchange_start/finish`; in `calculate_forces`, start the position exchange, compute **interior** cell forces (no ghost dependency, via `Cell::is_boundary()`), then finish the exchange, then compute **boundary** cell forces. Requires partitioning the Cabana pair loop into interior/boundary passes — verify physics identical (forces summed identically regardless of order). This is the invasive path; keep it behind clear structure.
- [ ] **If `post`/`pack`-dominated (message count/overhead) → reduce overhead:** e.g. coalesce co-directional peers into fewer larger messages, or ensure per-neighbor buffers are reused across steps (no per-step realloc), or fold the bonds message. Less invasive.
- [ ] **A/B** on the Release build (machine idle), 1/2/4/8 ranks, vs the pre-Ph4 async numbers: confirm the 8-rank regression shrinks and 1/2/4 don't regress. Record deltas + the `ghost/wait` change.
- [ ] **Physics parity** on /ssd RelWithAssert (lj/lees_edwards/p3m/collision/nsquare/hybrid/lb-coupling) — identical results. Commit. If the mitigation doesn't help (or regresses), document the negative result + revert the invasive change (keep 4.1's markers).

### Task 4.4: Progressive receive-unpack in `halo_exchange_finish` (wait_any pipelining)

**Files:** `src/core/ghosts/HaloExchange.cpp` (`halo_exchange_finish` only; `halo_exchange_start` and the request layout stay as-is). No header/API change.

**Motivation (from 4.2/4.3 profiles):** at 10k/8 the update phase pays `ghost/wait` 1.9% + `ghost/unpack` 3.1% *sequentially*, the reduce phase wait 3.1% + add 1.4%. Unpacking each neighbor's buffer as soon as its receive completes hides unpack under wait (and vice versa) — ceiling ~3% of step time at 8 ranks, and it composes with any later compute/comm overlap.

**Request-vector layout (invariant, set in `halo_exchange_start`):** `bufs.requests[0..n)` = per-neighbor irecvs (neighbor order), `[n..2n)` = isends; when `GHOSTTRANS_BONDS` is set, `[2n..3n)` = bond irecvs and `[3n..4n)` = bond isends.

- [ ] **Step 1 (hot path, no bonds):** replace the single `wait_all` + unpack-all with:
  - **Overwrite combine (position/property push):** loop `boost::mpi::wait_any` over the *active* recv sub-range; on completion, immediately unpack that neighbor's buffer (arrival order — safe: the validator guarantees each ghost cell is filled by exactly one message). Swap the completed request out of the active range (`std::iter_swap` with the last active slot) and maintain a slot→neighbor index map so re-waiting a completed request is impossible.
  - **Add combine (FORCE / RATTLE reduce):** wait the receives in **fixed neighbor order** (`wait(requests[i])` → `add_forces`/`add_rattle` for i = 0..n-1). Float addition is not associative and several neighbors add into the same local cells — fixed order keeps runs bitwise reproducible while still pipelining (adding neighbor i overlaps delivery of i+1..n-1).
  - After the receive loop, `wait_all` the send sub-range `[n..2n)` (buffers must be complete before next-step reuse).
- [ ] **Step 2 (cold bonds path):** when `GHOSTTRANS_BONDS` is set, keep the existing behavior (single `wait_all` over everything, then unpack in neighbor order) — `unpack_cells` needs a neighbor's data *and* bond message together, and the path is resort-only.
- [ ] **Step 3 (Caliper):** keep the label set unchanged (`ghost/wait`, `ghost/unpack` — `testsuite/python/caliper.py` EXPECTED_LABELS must NOT change): mark each wait with `ghost/wait` BEGIN/END around the `wait_any`/`wait(i)` call and each unpack with `ghost/unpack` BEGIN/END around it (repeated sequential begin/end of the same region accumulates — valid Caliper usage; regions stay properly nested).
- [ ] **Step 4 (edge cases):** n==0 (1 rank) must remain a no-op; `wait_any` requires a non-empty range — guard the loop. Local copies + collective stay before the wait loop (unchanged).
- [ ] **Step 5 (verify):** `make -C /ssd/weeber/comm-build -j8 check_unit_tests` = 149/149; parity on the /ssd RelWithAssert build: lj@2+4, lees_edwards@4, collision_detection@4, nsquare@4, hybrid_decomposition@4; `ctest --test-dir /ssd/weeber/comm-build -R caliper` green.
- [ ] **Step 6:** format + commit `core/ghosts: progressive wait_any receive-unpack in halo_exchange_finish`.
- [ ] **A/B (controller-led):** Release build, machine idle, lj @ 1k (1/2/4/8) + 10k (4/8) vs the staggered same-build numbers; record deltas.

---

## Self-Review (author)
- 4.1 (markers) is low-risk and enables the profiling that should drive the fix — avoids blindly doing the invasive overlap.
- 4.3 is profile-driven + correctness-gated (physics identical) + value-gated (A/B). The split-phase seam + interior/boundary tags (from Ph1/Ph3) are the enablers already in place.
- The 8-rank regression's mechanism is currently a hypothesis (many-small-messages); 4.2 confirms it before 4.3 commits to a fix.
