# Task 5.2 Report — Robust interior classification + validator overlap invariant

## Summary

All items delivered.  149/149 unit tests pass; full parity gate passes (lj@1/2/4,
lees_edwards@4, collision_detection@4, nsquare@4, hybrid_decomposition@4) on the
`/ssd/weeber/comm-build` RelWithAssert build.

---

## Wrap-detection variant implemented

**Precise pair-predicate** (not the conservative layer-marking).

`mark_boundary_cells` gains an optional third parameter:

```cpp
std::function<bool(Cell const *, Cell const *)> wrap_predicate = nullptr
```

`RegularDecomposition` supplies a lambda that fires iff the two cells sit in
the "first layer ↔ last layer" pair along a wrap axis:

```
axis i wraps  ⟺  node_grid[i]==1  &&  box.periodic(i)
first local layer along i  ⟺  ghost_coord[i] == 1
last local layer along i   ⟺  ghost_coord[i] == cell_grid[i]
```

Ghost-grid coordinate of a cell is computed from its linear offset into
`RegularDecomposition::cells` via column-major unpacking.  The predicate is
built only when at least one wrap axis exists; when `has_wrap_axis==false` a
`nullptr` is passed and the inner loop has no overhead.

This is the precise variant: only first↔last layer *pairs* trigger the rule,
not all cells in those layers individually.  A cell in the first layer that
only has neighbors in layers 2…cell_grid[i]-1 (impossible in a 3×3×3+ grid
since stencil reach is ±1, but correct to reason about) would not be marked
boundary by the predicate alone.  In practice every cell in layer 1 has a
stencil neighbor in layer cell_grid[i] when cell_grid[i]==1 (degenerate 1-cell
grid, still periodic), and every cell in layer 1 or cell_grid[i] has a neighbor
in the opposite layer when cell_grid[i]≥2, so the predicate is precise and
conservative in the appropriate sense.

---

## AtomDecomposition handling

`mark_boundary_cells(local, ghost)` is called first (handles multi-rank case
where all local-cell ghost-neighbors exist), then an explicit loop forces
`m_is_boundary = true` on all local cells.  This catches the single-rank case
(no ghost cells, no neighbours) where the base rule would leave the one local
cell interior.  Documented in comment: "no spatial locality — interior always
empty".

---

## HybridDecomposition handling

Same conservative approach: `mark_boundary_cells` + explicit all-boundary loop.
Comment explains: combining regular + n-square neighbours makes it non-trivial to
determine which cells remain genuinely interior; the overlap degrades gracefully
to a no-op interior pass.

---

## Validator overlap-safety invariant

Added to `validate_halo_plan` after the existing interior/boundary consistency
check.  Wording (in the violation strings):

- `"interior cell appears as NeighborComm send source (peer=N) — overlap-safety invariant violated"`
- `"interior cell appears as LocalComm src — overlap-safety invariant violated"`

The invariant is also documented in `HaloPlanValidator.hpp` as check #6:

> An interior cell must not appear as a NeighborComm send source or a LocalComm
> src.  This is the exact precondition that the integrator step-2 / force-reduce
> overlap (Task 5.3) relies on.

Implementation: builds an `unordered_set<ParticleList const *>` of interior
cell particle lists; skips the set scan when it is empty (all-boundary
decompositions: Atom/Hybrid).

---

## Test evidence

### Hand-made cell tests (Cell_boundary_test, no MPI)

| Test | Assertion |
|------|-----------|
| `wrap_predicate_makes_boundary` | always-fire predicate → both local cells boundary, no ghost needed |
| `no_wrap_no_ghost_stays_interior` | never-fire predicate + no ghosts → both local cells interior |
| `validator_fires_for_interior_cell_in_send_region` | interior cell as send source → "overlap-safety invariant" violation |
| `validator_fires_for_interior_cell_in_local_comm_src` | interior cell as LocalComm src → violation |
| `validator_silent_for_boundary_cell_in_send_region` | boundary cell as send source → no violation |

### 1-rank RegularDecomposition test (HaloPlanValidator_test, runs at NUM_PROC 1)

`single_rank_wrap_axis_cells_are_boundary`: constructs `RegularDecomposition`
with `{1,1,1}` node_grid, box 6.0, range 1.0 (gives cell_grid={6,6,6}).
Asserts: every cell in the first/last local layer (ghost-coord==1 or
ghost-coord==cell_grid[i]) is boundary.  Also asserts strictly-interior cells
(all ghost-coords in [2,5]) remain interior.

### Correctness under ESPRESSO_ADDITIONAL_CHECKS

The `assert(validate_halo_plan(...).empty())` in all three decomposition
constructors runs on the RelWithAssert build.  The new overlap invariant checks
are inside `validate_halo_plan`, so they are exercised on every real-decomposition
construction including 1-rank.  All pass (no interior cell appears as a send
source in any real decomposition).

---

## Concerns / Notes

None.  Conservative choice for Atom/Hybrid is correct: the overlap in Task 5.3
degrades to the blocking path when interior set is empty, which is the same as
today.  The precise pair-predicate for RegularDecomposition is exact: it fires
on exactly the cells whose neighbour stencil wraps, which is what Task 5.3
needs to be safe.

---

## Fix round (2026-07-29) — degenerate cell_grid==1 wrap axis

**Bug corrected:** The prior report's claim "the precise pair-predicate is
exact — it fires on exactly the cells whose neighbour stencil wraps" was wrong
for the degenerate case `cell_grid[i]==1` on a wrap axis.

**Root cause:** When `node_grid[i]==1 && periodic[i] && cell_grid[i]==1`, the
single cell layer's periodic neighbour folds back to the cell itself.
`init_cell_interactions()` drops the self-pair via the `ind1==ind2` guard
(~line 552 of RegularDecomposition.cpp), so `neighbors().all()` has **no**
entry along that axis.  The wrap predicate is therefore never invoked, and
`mark_boundary_cells()` leaves the cell INTERIOR — wrong, because the cell
interacts with itself across the periodic boundary.

**Fix** (commit `24fe5fdec0`, subject: `cell_system: cell_grid==1 wrap axes
force all-boundary (fix false interior)`):

After `mark_boundary_cells(local_cells(), ghost_cells(), wrap_pred)`, added a
post-pass in `RegularDecomposition.cpp`:

```cpp
for (int i = 0; i < 3; ++i) {
  if (node_grid[i] == 1 && m_box.periodic(i) && cell_grid[i] == 1) {
    for (Cell *c : local_cells())
      c->m_is_boundary = true;
    break; // one degenerate axis is enough to force all-boundary
  }
}
```

**New test** (`HaloPlanValidator_test.cpp::single_rank_cell_grid_one_all_cells_are_boundary`):
constructs `make_dd({1,1,1}, box_l=1.0, range=2.0)` which clamps `cell_grid`
to `{1,1,1}` on all three axes, then asserts every local cell is boundary.
**RED before fix** (BOOST error: "cell_grid==1 wrap axis: every local cell must
be boundary"), **GREEN after fix** (*** No errors detected).

**Minor additions in same commit:**
- `HaloPlanValidator.hpp` check #6 doc: added note that `LocalComm.dst`/
  `NeighborComm.recv` are omitted because they are ghost cells by construction.
- `mark_boundary_cells.hpp`: collapsed two sequential break-ing `if`s into one
  combined condition for readability.

**Build/parity:** 149/149 unit tests pass; `lj@1`, `lj@4`, `nsquare@4` pass.
