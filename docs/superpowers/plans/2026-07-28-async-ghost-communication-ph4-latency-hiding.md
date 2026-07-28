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

---

## Self-Review (author)
- 4.1 (markers) is low-risk and enables the profiling that should drive the fix — avoids blindly doing the invasive overlap.
- 4.3 is profile-driven + correctness-gated (physics identical) + value-gated (A/B). The split-phase seam + interior/boundary tags (from Ph1/Ph3) are the enablers already in place.
- The 8-rank regression's mechanism is currently a hypothesis (many-small-messages); 4.2 confirms it before 4.3 commits to a fix.
