# Dynamic ghost-layer neighborships for Lees-Edwards

Date: 2026-08-19
Branch: `comm_le` (worktree `comm_le`, based on `comm`)
Status: design, awaiting review

## Problem

Lees-Edwards (LE) shear needs the short-range loop to see every pair
across the shear boundary for any position offset. Today this is bought
with the `fully_connected_boundary` attribute of the regular
decomposition: every cell on the shear-plane-normal boundary is wired to
*all* cells along the shear direction, and `LeesEdwardsBC::distance()`
applies the `±pos_offset` correction inside the force loop. The
fully-connected stencil guarantees the correct partner is always among
the candidates, at the cost of O(cells-along-shear) neighbors per
boundary cell.

Two problems follow. First, the neighbor count on the boundary is
large and offset-independent, so it never shrinks to what the geometry
actually needs. Second, `fully_connected` requires all shear-direction
cells to be local: `RegularDecomposition.cpp:550` throws
`"The MPI nodegrid must be 1 in the fully connected direction."`. You
cannot split the MPI domain along the shear direction.

## Goal

Make the physical-cell -> ghost-cell wiring (and, when the grid is split
along the shear direction, the MPI peer for each region) a function of
the current `pos_offset`, rebuilt as the offset evolves. A narrow
stencil then sees exactly the partners it needs, for any offset, without
`fully_connected_boundary`, including when the MPI grid splits along the
shear direction.

### Non-goals

- No change to the LE physics, protocols, or the `distance()` math.
- No change to the thermostat, integrator, or force kernels.
- No new user-facing API. LE parameters drive the behavior.

## Acceptance criteria

The DPD benchmark (new, see Deliverables) is the fixture. Correctness is
proven three ways.

1. Bitwise-identical trajectories, per fixed decomposition: the dynamic
   code reproduces the current `fully_connected` output bit for bit for
   the *same* MPI/OpenMP configuration. This is only defined on grids
   that `fully_connected` supports, i.e. `node_grid[shear_dir] == 1`.
   Bitwise identity across *different* decompositions is not expected
   (force-summation order differs).
2. `non_bonded_loop_trace` pair sets match across all decompositions,
   including grids that split along the shear direction. The set of
   pairs the short-range loop visits (id1, id2, vec21) is
   decomposition-invariant and equal to the fully-connected reference.
3. Time-averaged shear stress agrees across decompositions and with the
   fully-connected reference (within statistics).

## Current architecture (as-is)

Key locations on the `comm` branch.

- `LeesEdwardsBC` and `distance()`:
  `src/core/lees_edwards/LeesEdwardsBC.hpp`. `distance()` subtracts
  `pos_offset` from the shear-direction component when the pair crosses
  the shear-plane-normal boundary.
- `fully_connected_boundary`: stored as
  `std::optional<std::pair<int,int>> m_fully_connected_boundary`
  (`RegularDecomposition.hpp:84`), `{shear_normal, shear_dir}`.
- Neighbor stencil: `RegularDecomposition::init_cell_interactions()`
  (`RegularDecomposition.cpp:486`). The fully-connected expansion sets
  `lower_index[fc_dir] = -1; upper_index[fc_dir] = global_size[fc_dir]`
  at the boundary (lines 579-587). Neighbors are collected into a
  flat_set ordered by `folded_linear_index` and split red/black by that
  index (lines 600-656).
- Halo plan: `RegularDecomposition::make_halo_plan()`
  (`RegularDecomposition.cpp:660`). For each ghost cell it computes the
  mirrored global cell `ghost_global`, its owning `peer` via
  `owner_of()`, and a `recv_key`/`send` pairing keyed by `global_key()`.
  All `SendRegion.shift` are currently `{}` (line 803).
- Ghost exchange: `src/core/ghosts/HaloExchange.*` and
  `particle_packing.cpp`. `pack_regions()` applies a common per-neighbor
  shift; `HaloExchange.cpp` already documents that a Lees-Edwards
  refactor may store distinct per-region shifts.
- Force loop order: `src/core/algorithm/link_cell.hpp`. For each `p1` it
  visits same-cell later particles, then red-neighbor cells in stored
  order; pairs are counted once via the red/black split.
- Plan and interactions are built once, in the constructor
  (`RegularDecomposition.cpp:829-835`), and cached in `m_halo_plan`.
  `resort()` (`RegularDecomposition.cpp:187`) does not rebuild them.
- Resort trigger: `resort_particles_if_needed()`
  (`integrate.cpp:373`) calls `check_resort_required(offset)` with
  `offset = LeesEdwards::verlet_list_offset(...)`.
  `check_resort_required` (`CellStructure.cpp:774`) uses
  `lim = (skin/2)^2 - offset.norm2()`, so the LE offset drift eats into
  the skin budget. `m_le_pos_offset_at_last_resort` is refreshed on
  resort (`CellStructure.cpp:672`).
- `m_box` is a live `BoxGeometry const&` (`RegularDecomposition.hpp:82`),
  so the decomposition reads the current `pos_offset`,
  `shear_direction`, and `shear_plane_normal` from
  `m_box.lees_edwards_bc()`.

## Design

Keep the offset where it already lives, in `distance()`, and change only
which cell feeds each shear-crossing ghost and which cells are wired as
neighbors. Ghost positions keep carrying only periodic-image shifts,
never the LE offset, so `distance()`'s subtraction stays byte-for-byte
what it computes today. Three coordinated changes plus a rebuild hook.

Notation: `sn = shear_plane_normal`, `sd = shear_direction`,
`cs = cell_size[sd]`, `O = pos_offset`, `s = round(O / cs)` the integer
column shift, `frac = O - s*cs` the residual in `[-cs/2, cs/2]`.

### 1. Halo sourcing (`make_halo_plan`)

For a ghost cell that crosses the shear-plane-normal boundary
(`side[sn] != 0`), add an integer column shift along `sd` to the
mirrored global index before resolving owner and key:

```
ghost_global[sd] += -side[sn] * s;   // s = round(pos_offset / cs)
```

Then `owner_of(ghost_global)` and `global_key(ghost_global)` naturally
select the sheared source column and its owning rank. This is exactly
why the `node_grid[sd] == 1` restriction disappears: the sheared source
may live on another rank, and `owner_of` finds it. Any wrap of the
shifted column past the shear-direction periodic boundary becomes an
ordinary image shift stored in `SendRegion.shift` (a periodic shift,
a multiple of the box length; never the LE offset).

The dual send-cell pairing (lines 773-784) is mirrored the same way, so
`recv[k]` still lines up with the peer's `send[k]` by shared key.

When `node_grid[sd] == 1`, the sheared source is local and this
degenerates to `LocalComm` self-copies of shifted columns; one code path
covers both the local and cross-rank cases.

### 2. Neighbor stencil (`init_cell_interactions`)

Replace the fully-connected all-columns expansion at the boundary with a
narrow window along `sd`, centered on the integer-shifted position:

```
if at_boundary(sn, {m,n,o}):
    center = m - s                       // integer-shifted center
    reach  = ceil((cutoff + skin) / cs) + 1   // +1 covers |frac| <= cs/2
    lower_index[sd] = center - reach
    upper_index[sd] = center + reach
```

The window is sized from `cutoff + skin` (matching the non-LE stencil)
plus one cell for the fractional offset. `reach` is small (a few cells),
not the full row. The neighbors still go through the same
`folded_linear_index` flat_set and red/black split, so ordering is
preserved (see Correctness).

The exact minimal `reach` is confirmed by `non_bonded_loop_trace`: the
test proves no pair within cutoff is missed, and lets us tighten the
window.

### 3. Rebuild on the resort path

Rebuild the shear-crossing halo plan and the shear-boundary cell
interactions when a resort fires, using the current offset. Do not add a
separate integer-crossing trigger: the resort cadence already bounds the
offset drift (see Skin below), and since `cs >= cutoff + skin > skin/2`,
`s` changes by at most about one per resort interval, so it is subsumed.

Between resorts nothing needs rebuilding: the narrow stencil plus the
continuous offset in `distance()` stay correct, and ghost positions are
already re-communicated every step by the existing position update.

Rebuild scope: the first implementation may rebuild the whole plan and
all cell interactions on resort (simple, correct). If profiling shows
this matters, narrow it to the shear-crossing `NeighborComm`s and the
shear-boundary cells only.

### The role of skin

The skin bounds how far the offset can drift before a forced resort, and
therefore how stale the connections can get. `check_resort_required`
uses `lim = (skin/2)^2 - |offset drift|^2`; once the offset alone drifts
by `skin/2`, `lim < 0` and a resort is forced. So the offset never
drifts more than `skin/2` between resorts. Two design points follow:
the rebuild cadence is the resort cadence (already skin- and LE-aware),
and the connection window must be sized from `cutoff + skin` so it stays
valid across the permitted `skin/2` drift. This is why `reach` uses
`cutoff + skin`, not `cutoff`.

## Correctness argument

### Bitwise identity to `fully_connected` (where comparable)

For `node_grid[sd] == 1`, take any cross-boundary pair `(p1, p2)`.

- Same particle data. The sheared source column now lands in ghost slot
  `c + s` instead of slot `c`, but it holds the same particles, with the
  same real positions (only a periodic image applied, as before, never
  the LE offset). `distance()` uses the same `pos_offset`, so the
  computed distance is byte-identical.
- Same summation order. Neighbors are stored ordered by
  `folded_linear_index`; the force loop walks red-neighbors in that
  order. Because `s` is a constant integer shift, slot index `c + s` is
  monotonic in the source column `c`, so in-range partners appear in the
  same relative order as in the fully-connected list. Columns the narrow
  window omits held no in-range partner (guaranteed by the `cutoff+skin`
  window, verified by the trace), so omitting them does not change the
  order of the ones that remain.

Same pairs, same order, same values => bitwise-identical forces =>
bitwise-identical trajectory for that decomposition.

### Split along the shear direction

`fully_connected` cannot run these grids, so there is no bitwise
baseline. Correctness rests on criteria 2 and 3: the
`non_bonded_loop_trace` pair set equals the reference set (proving no
pair is missed or duplicated), and the average shear stress agrees. The
halo sourcing change routes the sheared columns from the owning ranks;
`distance()` is unchanged.

## Edge cases and risks

- Fractional offset coverage: `frac in [-cs/2, cs/2]` plus `skin/2`
  drift. The `+1` cell in `reach` covers it. Risk: an off-by-one at the
  window edge. Mitigation: the trace test is the gate; size
  conservatively first, then tighten.
- Oscillatory shear: offset drift changes sign. The window is symmetric
  (`center ± reach`), so both drift directions are covered.
- Corner ghosts (crossing `sn` and another axis at once): the `sd` shift
  applies only to the `sn`-crossing component; other crossings keep their
  normal image shift. Must be handled in `make_halo_plan` per axis.
- `sd` periodic wrap: the shifted column can wrap the box; the resulting
  image shift goes into `SendRegion.shift`. Must not be confused with the
  LE offset.
- Distinct per-region shifts: `HaloExchange.cpp` warns that per-region
  shifts must still write bonds into one shared archive. If the wrap
  shift makes shifts differ within a `NeighborComm`, honor that note.
- Boundary-cell marking: `mark_plan_cells_boundary`
  (`RegularDecomposition.cpp:925`) already anticipates LE plan-driven
  boundary cells; verify it covers the sheared sources.
- `ESPRESSO_ADDITIONAL_CHECKS` plan validation
  (`validate_halo_plan`, `RegularDecomposition.cpp:928`) must still pass
  for the sheared plan.

## Deliverables

1. `maintainer/benchmarks/dpd.py`. Modeled on `lj.py` structure
   (argparse via `benchmarks.add_common_args`, build/tune/state/run
   modes, CSV report), with the physics from `~/dpd/lin_shear.py`: DPD
   non-bonded interaction, `thermostat.set_dpd`, a `LinearShear` LE
   protocol, and `set_lees_edwards_bc`. Options for shear velocity,
   density, and node grid. It doubles as the correctness fixture.
2. The halo refactor: changes 1-3 above in
   `RegularDecomposition.cpp/.hpp`, reading LE parameters from
   `m_box.lees_edwards_bc()`, plus the resort-path rebuild hook.
3. Verification harness (script and/or testsuite addition mirroring
   `lees_edwards.py::run_lj_pair_visibility`) implementing the three
   proofs.

## Verification plan

Decomposition matrix for the DPD fixture, `shear_dir = x`,
`shear_plane_normal = y`:

- 1 rank (reference for bitwise per-config comparisons).
- Split along the neutral axis `z` (e.g. `[1,1,N]`): `fully_connected`
  valid -> bitwise vs reference.
- Split along the shear-normal `y` (e.g. `[1,N,1]`): `fully_connected`
  valid -> bitwise vs reference.
- Split along the shear direction `x` (e.g. `[N,1,1]`): new capability,
  `fully_connected` invalid -> trace + shear stress only.
- Mixed grids (e.g. `[2,2,1]`, `[2,1,2]`).
- OpenMP thread counts 1, 2, 4 for each of the above.

For each: (a) dump full trajectory, compare bitwise to the matching
`fully_connected` run where valid; (b) collect `non_bonded_loop_trace`,
compare the pair set to the reference; (c) collect the pressure tensor
shear component over time, compare the average.

## Rollout

Keep the `fully_connected` code during development as the bitwise
reference. Make the dynamic sheared remap the automatic behavior when LE
is active (no `fully_connected_boundary` needed for correctness). Once
the proofs pass, remove or deprecate `fully_connected_boundary`.

## Open questions

- Rebuild scope on first cut: whole plan vs shear-crossing only. Start
  whole, optimize if needed.
- Whether to keep `fully_connected_boundary` as a deprecated no-op alias
  or remove the API outright at the end.
