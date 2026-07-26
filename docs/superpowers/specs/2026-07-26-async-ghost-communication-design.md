# Async, topology-agnostic ghost communication — design

- **Date:** 2026-07-26
- **Worktree / branch:** `/tikhome/weeber/es/.claude/worktrees/comm` on `worktree-comm`
- **Status:** design approved; pending spec review → implementation plan
- **Benchmarks used for A/B:** `maintainer/benchmarks/lj.py`, `maintainer/benchmarks/p3m.py` at 1/2/4/8 MPI ranks

## 1. Goal

Extend the flexibility and performance/scalability of ESPResSo's MPI ghost communication:

1. Replace the **staggered** communication (one spatial direction after another) with a set of
   **asynchronous point-to-point** messages.
2. Enable **latency hiding** — overlap communication with communication now, and open a clean seam for
   overlapping communication with computation.
3. Provide a **clearer interface** for declaring *what* must be communicated *where*, plus **automated
   validation** that (a) all halos get communicated and (b) cell neighborships match the halo requirements.

**Forward-looking constraint:** the design must **not assume a Cartesian node grid**. Tree / space-filling-curve
(SFC) load balancing will be added later, producing irregular, per-rank, time-varying neighbor sets.

## 2. Current architecture (what we are replacing)

- **Engine** — `ghost_communicator()` (`src/core/ghosts.cpp:443`) is a dumb serial executor: it walks an ordered
  `std::vector<GhostCommunication>` front-to-back and **blocks** on each op (`comm.send`/`comm.recv`/`broadcast`/
  `reduce`). A single pair of **function-static** `send_buffer`/`recv_buffer` (`ghosts.cpp:448`) means only one
  message is ever in flight. Each op does **two** blocking transfers (payload, then a separate bonds transfer).
- **Staggering is a genuine data dependency.** `RegularDecomposition::prepare_comm` (`RegularDecomposition.cpp:605`)
  sweeps ±x, ±y, ±z; a `done[]` accumulator widens each later direction's region so ghosts received along x are
  re-shipped along y (carrying edges/corners). x→y→z is therefore a hard dependency chain. Deadlock-freedom relies
  on an even/odd parity trick (`RegularDecomposition.cpp:677,689`). The force path is the same schedule reversed
  (`revert_comm_order`).
- **"Latency hiding" is fake.** `GHOST_PREFETCH`/`GHOST_PSTSTORE` only reshuffle *serialization* around the same
  blocking calls; nothing overlaps communication with force computation.
- **No declarative halo, no validation.** What-goes-where is implicit in the ordered vector. The halo requirement is
  one hard-coded if-chain (`System::get_global_ghost_flags`, `src/core/system/System.cpp:550`). Cell neighborships
  (`init_cell_interactions`) and the comm schedule (`prepare_comm`) are derived **independently**, with nothing
  cross-checking them — so they can silently drift.
- **Fragile reduction gate.** Additive write-back is gated on `data_parts == GHOSTTRANS_FORCE` **exactly**
  (`ghosts.cpp:538`), so forces cannot be reduced bundled with any other part.
- **P3M charge-assignment halo** (`p3m_send_mesh`, `src/core/p3m/send_mesh.cpp`) is a *structurally identical*
  staggered, blocking, direction-by-direction exchange, used by **all** P3M variants. HeFFTe modernized only the FFT
  (already async point-to-point internally); the CA halo remains the un-modernized staggered part on P3M's hot path.
- **Profiling gap.** Caliper is integrated (`-D ESPRESSO_BUILD_WITH_CALIPER=ON`, source guard `ESPRESSO_CALIPER`,
  runtime `CALI_CONFIG=runtime-report`); the integrator/force path is instrumented, but **ghost communication has
  zero markers**, so its cost is only visible lumped inside `calculate_forces`/`resort`.
- **Branch context:** this is the Cabana/Kokkos branch — the force path is `cabana_short_range`; the ghost machinery
  still runs on host/CPU.

## 3. Scope

- **In scope (this effort):** rewrite the **particle ghost layer** (`ghosts.{cpp,hpp}` + the regular decomposition's
  comm construction) as an async, topology-agnostic, declarative, validated subsystem. Build the async transport as a
  **reusable, payload-generic** component.
- **Follow-on (Ph5):** port the **P3M charge-assignment halo** (`p3m_send_mesh`) onto the same engine via a
  mesh-block pack/unpack policy.
- **Out of scope:** LB/EK halo exchange (owned by waLBerla / HeFFTe); the legacy `fft.cpp` domain exchange
  (DP3M-only); replacing HeFFTe's internal FFT reshapes (already async p2p).

## 4. Requirements & constraints

- **R1 — Async point-to-point, topology-agnostic.** No Cartesian assumptions; arbitrary, per-rank neighbor sets.
- **R2 — Split-phase API now.** `start` → independent work → `finish`. Delivers message-level overlap immediately and
  a seam for compute overlap; the force-loop interior/boundary restructuring is **deferred** and evidence-gated.
- **R3 — General declarative framework (all decompositions).** Regular, N-square (Atom), and Hybrid all emit the same
  plan object; one unified validator checks any of them. Collectives are one pattern in the taxonomy.
- **R4 — Interior/boundary cell classification is first-class now.** It serves both future compute-overlap and the
  Euclidean-vs-minimum-image (MI) distance optimization.
- **R5 — Preserve physics exactly in the transport swap.** Ghost *contents* and the position convention are unchanged
  (direct routing changes only the *path*). Specifically: **ghost positions stay unshifted** — the LB particle
  coupling reads ghost positions assuming they are unshifted, and shifts would also entangle with Lees-Edwards. We do
  **not** reintroduce meaningful `shift` values in this effort.
- **R6 — Compatibility bar.** Everything stays green: all three cell systems (regular, n-square, hybrid),
  Lees-Edwards, virtual sites, RATTLE, collision detection, LB coupling/tracers, bonds, `maxset` features, the full
  testsuite, and CI on the `me` remote.
- **R7 — Performance bar.** At least **neutral** (no regression) on LJ and P3M at 1/2/4/8 ranks; latency hiding
  explored and **measured** via Caliper. A scaling improvement is the goal but **not a merge gate**.

## 5. Design

### 5.1 Data model — topology-agnostic halo declaration

One plan object per decomposition, rebuilt on resort, describing *what goes where* with no grid assumption:

```cpp
struct SendRegion {
  ParticleList *cell;
  Utils::Vector3d shift;                 // applied at pack time; preserved at current (unshifted) semantics
};
struct NeighborComm {
  int peer;                              // arbitrary rank, not a grid neighbor
  std::vector<SendRegion> send;          // local cells to pack & send
  std::vector<ParticleList *> recv;      // ghost cells to fill; recv[k] ⟷ peer.send[k]
};
struct CollectiveSection { /* n-square all-to-all: pattern + owned cell(s) */ };
struct HaloPlan {
  boost::mpi::communicator comm;
  std::vector<NeighborComm> neighbors;           // point-to-point section
  std::optional<CollectiveSection> collective;   // n-square (Atom / hybrid-nsquare)
};
enum class Direction { Push, Reduce };           // real→ghost, ghost→real
enum class Combine { Overwrite, Add };
struct ExchangeOp { Direction dir; Combine combine; };
```

- **Ordering contract** `recv[k] ⟷ peer.send[k]` is symmetric-by-construction; the validator enforces it, so no
  parity/ordering trick is needed.
- **Op is explicit, not inferred.** An exchange carries an `ExchangeOp`. This removes the
  `data_parts == GHOSTTRANS_FORCE` exact-equality gate — forces reduce because the op says `Reduce/Add`, so
  reductions compose with other parts. The **reverse (force-collect) plan is a view** of the same neighbors with
  `send`/`recv` swapped — replacing `revert_comm_order`; forward and reverse cannot drift.
- **Payload-generic.** The engine is parameterized by a pack/unpack policy (particle serialization now; contiguous
  mesh-block copy in Ph5), so the same async orchestration serves the P3M CA halo.
- **`shift` preserved, not activated.** Kept as a field and plumbed exactly as today (current value) so Ph1 is a pure
  transport swap. We do **not** make it meaningful (see R5).

### 5.2 Async engine + split-phase API

```cpp
GhostExchange ghost_exchange_start(HaloPlan const &, unsigned data_parts, ExchangeOp);
void          ghost_exchange_finish(GhostExchange &);          // Waitall + unpack/reduce
inline void   ghost_exchange(/* … */) { auto h = start(…); finish(h); }  // blocking convenience
```

- **Protocol (deadlock-free on any neighbor graph):** post all `Irecv` → pack + post all `Isend` → `Waitall` →
  unpack (or add, for `Reduce`). No parity, no relay; identical for regular and future tree/SFC.
- **Sizes known a priori.** `ghosts_count()` (PARTNUM) already runs after each resort and sizes ghost cells, so recv
  sizes are computable from the plan — **no `MPI_Probe` on the hot path**. Bonds (cold, resort-only) get a small size
  header and never burden per-step POSITION/FORCE exchanges.
- **One message per neighbor**, payload+bonds folded into one contiguous per-neighbor buffer — replaces today's two
  blocking transfers per op.
- **Buffers owned by the exchange/plan and reused across steps** (topology stable between resorts) — eliminates the
  `static` singletons that force serialization; N neighbor messages are genuinely in flight at once. One tag per
  exchange-kind (position / force / partnum / bonds) so future concurrent exchanges don't collide.
- **`GHOST_PREFETCH`/`GHOST_PSTSTORE` deleted** — real non-blocking MPI makes the hand-rolled buffer reshuffling
  obsolete.

### 5.3 Construction & neighbor discovery

- New decomposition method `HaloPlan make_halo_plan() const`, replacing `exchange_ghosts_comm()` /
  `collect_ghost_force_comm()` (reverse is a view, not a second build).
- **Regular decomposition** builds it directly: for each local cell, walk its interaction stencil; for every stencil
  cell that is a ghost, resolve the owning **peer** and the `(local send cell ⟷ ghost recv cell)` correspondence;
  group by peer → `NeighborComm`s. Corners/edges are **explicit peers** — no relay, no `done[]`.
- **Neighbor discovery is the pluggable seam:**
  - *Regular Cartesian:* peer = cart-shift of my coords by the stencil offset (periodic mod grid) — analytic, cheap.
  - *Lees-Edwards:* the shear-boundary neighbor set is offset by the LE displacement, so a boundary cell maps to
    ghosts owned by potentially several peers. The topology-agnostic model expresses this natively.
  - *Future tree/SFC:* a discovery step ("who owns the region within cutoff of this boundary cell") emitting the same
    `HaloPlan`; engine + validator unchanged.
- **Ordering contract** honored by enumerating each shared face/edge/corner region in a **canonical deterministic
  order** (e.g. by global cell index) on both ranks, so lists line up without shipping index arrays — verified by the
  build-time handshake (5.4).
- **Interior/boundary tag** computed here: a local cell is *boundary* iff any stencil neighbor is a ghost, else
  *interior*.

**Lees-Edwards `fully_connected_boundary` removal.** `fully_connected_boundary()` (threaded through
`RegularDecomposition.cpp:418, 452, 487, 517`; constructed at `:713`) exists only because the axis-relay cannot
express sheared peer connections; it over-connects the entire shear-boundary row as a correctness sledgehammer. Precise LE-aware enumeration in
`make_halo_plan()` makes it **removable** — a clarity and performance win (far fewer ghost interactions at the shear
boundary). Sequenced safely:

1. **Ph1** preserves `fully_connected` (the stencil walk subsumes it) → pure transport swap, behavior identical.
2. **Ph3** builds the exact offset-dependent sheared connection set and **deletes `fully_connected_boundary`**, gated
   by the LE tests + the validator + an assertion that shear-boundary ghost count drops.

**Open sub-task (LE offset drift):** the LE offset changes every step while the plan rebuilds only on resort. Removing
`fully_connected` surfaces this: the shear-boundary band must either rebuild when the offset crosses a cell, or be
widened to cover the maximum offset drift per resort interval. (`fully_connected` hid this.) Nailed down in Ph3.

### 5.4 Unified validator (goal #3)

Runs once per plan build (cheap) and additionally under `ESPRESSO_ADDITIONAL_CHECKS`. Decomposition-agnostic
invariants:

1. **Coverage** — every ghost cell is the `recv` of *exactly one* `NeighborComm` (no unfilled, no double-write).
2. **Neighborship match** — ghost cells filled ⊇ ghost cells referenced by `init_cell_interactions`; i.e. no
   interaction reaches into an uncommunicated ghost. (This is "cell neighborships match halo requirements" verbatim.)
3. **Cross-rank symmetric handshake** — a lightweight neighbor exchange of `(peer, counts, cell-id hash)` at build
   time; assert my `send` to P mirrors P's `recv` from me, in order.
4. **Reverse = mirror** — the reduce plan is structurally the push plan with `send`/`recv` swapped (free; a view).
5. **Interior/boundary consistency** — interior cells have zero ghosts in stencil; boundary ghosts all covered by (1).
6. **Op sanity** — `Combine::Add` only with reducible parts (FORCE/RATTLE/torque); `Push` never adds.

### 5.5 Collectives (Atom / Hybrid)

- N-square all-to-all stays `bcast`/`reduce` (no p2p benefit) but is expressed as `HaloPlan.collective`. The engine
  dispatches: point-to-point section via `Irecv/Isend/Waitall`, collective section via boost::mpi. Same façade; the
  branch lives in the engine, not in callers.
- **Hybrid** becomes a `HaloPlan` carrying *both* a `neighbors` section (regular child) and a `collective` section
  (n-square child), executed in one exchange — replacing the current "concatenate two `communications` vectors" hack
  with a structured, validated plan.

### 5.6 Interior/boundary tags & Euclidean-vs-MI

- The interior tag (5.3) lets the pair loop **skip the MI fold for interior cells**, where the naive vector
  `r_i − r_j` already equals the MI vector by construction. This is **tag-based, not shift-based**: no stored ghost
  position changes, so LB particle coupling and Lees-Edwards are untouched (R5).
- The kernel-level Euclidean branch ships as an **opt-in, measured** change in Ph3; the tag itself is delivered
  earlier so the classification exists for both this and future compute-overlap.

## 6. Instrumentation & A/B methodology

- **Caliper markers (Ph0, before any behavior change):** mark the `CellStructure` ghost entry points
  (`ghosts_update`, `ghosts_reduce_forces`, `ghosts_count`, `update_ghosts_and_resort_particle`, `resort_particles`);
  inside the engine, split-phase regions **`ghost/pack` · `ghost/post` · `ghost/wait` · `ghost/unpack`**. `ghost/wait`
  is the key latency-hiding signal (shrinks as overlap improves). Update `testsuite/python/caliper.py`
  `EXPECTED_LABELS`.
- **A/B protocol:** CPU-only Release; `ADDITIONAL_CHECKS` off for timing runs; `performance.hpp` (LJ) and a
  P3M-capable myconfig (P3M). **Pin skin and cell grid** across variants (disable `tune_skin` and lj.py's in-loop
  `retune_skin_after_steps`) so both share identical halo geometry — otherwise the A/B is invalid. **Machine-idle
  gate** before every run (check load / `who`; do not measure while other users run heavy jobs); rank pinning
  `--bind-to core`, no oversubscription, at 1/2/4/8. Metrics: `benchmarks.py` mean ± 95% CI per step **and** Caliper
  `runtime-report` per region (total ghost time and `ghost/wait`). Baseline locked on the staggered code first; each
  phase compared against it.

## 7. Phases, migration & testing (strategy: incremental behind the `CellStructure` façade)

- **Ph0 — instrument + baseline.** Caliper markers on the existing path; update the caliper test; capture the A/B
  baseline. No behavior change.
- **Ph1 — engine + regular plan.** `HaloPlan` + async split-phase engine + regular `make_halo_plan()`
  (direct-neighbor; `fully_connected` preserved). Route the façade to the new engine; keep old `ghost_communicator`
  until parity is proven, then delete it. Gate: full testsuite (all 3 cell systems, LE, virtual sites, RATTLE,
  collision, **LB coupling/tracers**, bonds) + CI on `me`. TDD: engine unit test + a Python ghost microbenchmark (via
  `espressomd.profiler.Caliper`) exercising push/reduce on 2/4/8 ranks.
- **Ph2 — unified validator.** Validator + build-time cross-rank handshake, wired for regular/atom/hybrid. TDD: feed
  deliberately broken plans (missing sender, asymmetric counts, uncovered interaction neighbor) and assert each is
  caught.
- **Ph3 — interior/boundary + LE cleanup.** Interior/boundary tags; opt-in tag-based Euclidean pair-loop branch
  (measured); precise LE-sheared construction + **remove `fully_connected_boundary`** (+ offset-drift handling), gated
  by LE tests + validator + shear-boundary ghost-count assertion.
- **Ph4 — latency hiding.** Measure message-level overlap; prototype *one* bounded compute-overlap case behind
  `start/finish`; keep only if it measures well.
- **Ph5 — follow-on.** P3M CA halo (`p3m_send_mesh`) adopts the engine via a mesh-block pack/unpack policy.

Every phase: run `maintainer/format` on changed files (CI won't start otherwise), push to `me`, watch CI, triage
regressions inline.

## 8. Risks & open sub-tasks

- **LE offset drift vs resort cadence** (5.3) — the concrete algorithm for rebuild-on-cell-cross vs band-widening is
  decided in Ph3.
- **Direct-neighbor message count** — up to 26 small messages per exchange vs 6 larger ones. Expected win from
  concurrency + no serialization dependency, but must be confirmed by the Ph1 A/B (esp. small-subdomain / high-rank
  cases). If message count dominates at 8 ranks, consider coalescing co-directional peers.
- **Parity with corner/edge handling** — the direct scheme must reproduce exactly the ghost set the relay produced;
  the validator's coverage + neighborship-match checks are the guard, backed by full-testsuite diffing.
- **Cabana/Kokkos interaction** — ghost (un)packing stays host-side; ensure ordering vs AoSoA commit is preserved.
- **boost::mpi non-blocking ergonomics** — confirm `irecv`/`isend` + `wait_all` handle the raw-buffer + bonds split
  cleanly; fall back to native `MPI_Isend/Irecv` if needed.

## 9. Acceptance criteria

- All of R6 green (all cell systems, LE, vs, RATTLE, collision, LB coupling/tracers, bonds, `maxset`, testsuite, CI
  on `me`).
- Ghost communication is async point-to-point and topology-agnostic; the axis-relay, parity trick,
  `PREFETCH`/`PSTSTORE`, and static buffers are gone.
- A single declarative `HaloPlan` per decomposition drives both neighbor stencil consistency and the comm schedule;
  the unified validator enforces coverage + neighborship-match + cross-rank symmetry.
- Interior/boundary classification exists and is validated; the Euclidean-vs-MI branch is available (opt-in) and
  measured.
- `fully_connected_boundary` removed with LE tests green.
- LJ and P3M at 1/2/4/8 ranks: no regression; ghost time and `ghost/wait` reported via Caliper; latency-hiding
  potential quantified.
