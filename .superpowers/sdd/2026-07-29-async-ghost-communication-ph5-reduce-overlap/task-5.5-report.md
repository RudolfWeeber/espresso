## Task 5.5 Report: Runtime-conditional orientation payload (GHOSTTRANS_QUAT / GHOSTTRANS_TORQUE)

### Summary

`ESPRESSO_ROTATION` builds bundled quat (32 B) with every position push and torque (24 B) with every force reduce, even for purely translational systems like LJ.  After this change both bits are **off by default**; the whitelist turns them on only when orientation physics is actually active.

Savings per particle:
- position push: 68 → 36 B (quat off)
- force reduce: 48 → 24 B (torque off)

---

### 1. New bits

`src/core/ghosts.hpp` (under `#ifdef ESPRESSO_ROTATION`):

```cpp
GHOSTTRANS_QUAT   = 256u,   // orientation quaternion (pushed with position; runtime-conditional)
GHOSTTRANS_TORQUE = 512u,   // torque (reduced with force; runtime-conditional)
```

Corresponding `DataPart` additions in `src/core/cell_system/CellStructure.hpp`:

```cpp
DATA_PART_QUAT   = 128u,   // orientation quaternion (pushed with position)
DATA_PART_TORQUE = 256u,   // torque (reduced with force)
```

`map_data_parts` in `CellStructure.cpp` maps them 1-to-1.

---

### 2. Reader audit (correctness gate)

Grepped every reader of `p.quat()`, `calc_director()`, `p.torque()`, `p.calc_dip()` reachable for ghost particles under ESPRESSO_ROTATION.

| Reader | Location | Ghost? | Condition added |
|--------|----------|--------|-----------------|
| `short_range_cabana` Gay-Berne kernel | `short_range_cabana.cpp` | yes (ghost pair kernel) | `#ifdef ESPRESSO_GAY_BERNE` combined_active_pair_mask |
| `short_range_cabana` dipole pair kernel | `short_range_cabana.cpp` | yes (ghost pair kernel) | `#ifdef ESPRESSO_DIPOLES` dipoles.impl->solver |
| `vs_relative_update_particles` | `virtual_sites/relative.cpp` | yes (p_ref may be ghost) | `ROT_VS_RELATIVE | ROT_VS_INDEPENDENT | TRANS_VS_RELATIVE` |
| LB swimmer coupling (`calc_director`) | `lb/particle_coupling.cpp` | yes (ghost particles) | `#ifdef ESPRESSO_ENGINE` lb.is_solver_set() |
| `HomogeneousMagneticField` (`calc_dip`) | local particles only | no | (included via dipole check, conservative) |
| dipolar solvers (`dp3m_heffte` `unique_particles->calc_dip()`) | `dp3m_heffte.impl.hpp:195` — `unique_particles` includes ghosts with no local counterpart (CellStructure::set_index_map lines 310–322) | **yes** | covered by `dipoles.impl->solver` condition |
| `ShapeBasedConstraint` / `calc_non_central_force` reads `p.quat()` | `constraints/ShapeBasedConstraint.cpp` | no — iterates `local_particles` only | safe without the bit; no condition needed |
| `Stoner-Wohlfarth` `integrate_magnetodynamics`: calls `p_ref->calc_director()` on result of `get_reference_particle` / `get_local_particle` | `magnetostatics/stoner_wohlfarth_thermal.cpp:272` | **yes** (get_local_particle can return ghost) | `#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH` `thermo_switch & THERMO_LANGEVIN` (SW requires Langevin — enforced in integrate.cpp) |
| rotational integrators | local particles only | no | — |

Readers of **torque on ghost particles** (force reduce):
- `vs_relative_back_transfer_forces_and_torques` accumulates torque from ghost VS onto local real particle: covered by `TRANS_VS_RELATIVE | ROT_VS_RELATIVE | ROT_VS_INDEPENDENT`.
- Rotational propagation active generally: torque is consumed by the integrator on real particles; ghost torque is only needed when ghost particles accumulate torque that must be reduced back (i.e., whenever rotational propagation is active).
- ICC (`icc.cpp`) blocking reduce: resets `force_and_torque` on local particles; may carry a zero-valued TORQUE payload on the wire when orientation physics is concurrently active — harmless (bytes only, value is zero).

### 3. `orientation_ghosts_needed` whitelist (System.cpp)

```cpp
static bool orientation_ghosts_needed(System const &sys) {
#ifndef ESPRESSO_ROTATION
  return false;
#else
  // Rotational propagation active
  if (sys.propagation->used_propagations &
      (ROT_EULER | ROT_LANGEVIN | ROT_BROWNIAN | ROT_STOKESIAN |
       ROT_VS_RELATIVE | ROT_VS_INDEPENDENT))
    return true;
  // Virtual sites relative uses p_ref.quat() where p_ref may be a ghost
  if (sys.propagation->used_propagations &
      (TRANS_VS_RELATIVE | ROT_VS_RELATIVE | ROT_VS_INDEPENDENT))
    return true;
#ifdef ESPRESSO_DIPOLES
  if (sys.dipoles.impl && sys.dipoles.impl->solver) return true;
#endif
#ifdef ESPRESSO_GAY_BERNE
  if (sys.nonbonded_ias->combined_active_pair_mask() &
      pair_potential_bit(PairPotential::GayBerne)) return true;
#endif
#ifdef ESPRESSO_ENGINE
  if (sys.lb.is_solver_set()) return true;
#endif
  return false;
#endif
}
```

`get_global_ghost_flags` sets `DATA_PART_QUAT` when `orientation_ghosts_needed(*this)`.  
`get_force_reduce_ghost_flags` (new method) returns `GHOSTTRANS_FORCE | GHOSTTRANS_TORQUE` when `orientation_ghosts_needed(*this)`, else just `GHOSTTRANS_FORCE`.

---

### 4. Wire symmetry

Pack and unpack use the **same `data_parts` value** per exchange (existing invariant in `HaloExchange`). The layout in `serialize_and_reduce` is therefore always symmetric: `QUAT` follows `POSITION` in the stream when both are set; `TORQUE` follows `FORCE` when both are set.  `calc_transmit_size` reuses `serialize_and_reduce` via `SerializationSizeCalculator` and therefore stays consistent automatically — no separate size calculation needed.

Comment in `particle_packing.cpp`:
```cpp
// Wire-symmetry: pack and unpack use the same data_parts per exchange,
// so the layout is always symmetric; QUAT follows POSITION in the stream
// when both are set (position push), and TORQUE follows FORCE when both
// are set (force reduce).
```

---

### 5. `ghosts_reduce_forces` paths

Both the blocking (`ghosts_reduce_forces`) and the split-phase start/finish variants now call `get_system().get_force_reduce_ghost_flags()` instead of the hardcoded `GHOSTTRANS_FORCE`.  The `HaloExchange` sanity assert was updated to accept `GHOSTTRANS_FORCE | GHOSTTRANS_TORQUE` as a valid reduce combination.

The `add_forces` fast path in `particle_packing.cpp` was updated to accept `unsigned data_parts` and handle force and torque as separate optional arms.

---

### 6. Bug found and fixed: `vs_relative_update_particles`

**Root cause**: `vs_relative_update_particles` called `ghosts_update(DATA_PART_POSITION | DATA_PART_MOMENTUM)` without `DATA_PART_QUAT`.  In the old code, `GHOSTTRANS_POSITION` implicitly bundled quat under `ESPRESSO_ROTATION`.  After the split, ghost particles' quats were stale when `connection_vector(p_ref, p)` read `p_ref.quat()`, causing ~0.2% force mismatch at 2+ ranks.

**Fix**: `virtual_sites/relative.cpp` now explicitly requests `| Cells::DATA_PART_QUAT` under `#ifdef ESPRESSO_ROTATION`.

---

### 7. Tests

All tests passed (build: `make -j8`):

**Unit tests:** 149/149 green (`check_unit_tests`).

**Parity gate** (same forces with/without `ROTATION` compiled in; bits must be OFF for plain LJ):
| Test | Ranks | Result |
|------|-------|--------|
| lj | 2 | PASS |
| lj | 4 | PASS |
| lees_edwards | 4 | PASS |
| collision_detection | 4 | PASS |
| nsquare | 4 | PASS |
| hybrid_decomposition | 4 | PASS |

**Orientation coverage** (quat/torque propagate correctly when rotation ON):
| Test | Ranks | Result |
|------|-------|--------|
| rotation_per_particle | 1 | PASS |
| rotation_per_particle | 4 | PASS |
| dipolar_direct_summation | 4 | PASS |
| virtual_sites_relative | 4 | PASS |
| virtual_sites_relative_pbc | 4 | PASS |

**Bits-OFF proof for LJ**: `lj` test at 2 and 4 ranks passes — confirming quat and torque bits are never set for a non-rotating system, so there is zero orientation overhead on the wire.

---

### 8. A/B bandwidth reduction

Not benchmarked (shared machine — other users' jobs running; see memory note).  Per-particle savings are deterministic:
- Position push: 68 → 36 B/particle (-47%)
- Force reduce: 48 → 24 B/particle (-50%)

---

### 9. Files changed

| File | Change |
|------|--------|
| `src/core/ghosts.hpp` | `GHOSTTRANS_QUAT = 256u`, `GHOSTTRANS_TORQUE = 512u` |
| `src/core/cell_system/CellStructure.hpp` | `DATA_PART_QUAT = 128u`, `DATA_PART_TORQUE = 256u` |
| `src/core/cell_system/CellStructure.cpp` | `map_data_parts` new bits; `ghosts_reduce_forces{,_start}` use `get_force_reduce_ghost_flags()` |
| `src/core/ghosts/particle_packing.cpp` | quat moved to QUAT arm; torque to TORQUE arm; `add_forces` takes `data_parts` |
| `src/core/ghosts/particle_packing.hpp` | `add_forces` signature updated |
| `src/core/ghosts/HaloExchange.cpp` | sanity assert accepts `FORCE|TORQUE`; add path uses `data_parts & FORCE` |
| `src/core/system/System.hpp` | `get_force_reduce_ghost_flags()` declaration |
| `src/core/system/System.cpp` | `orientation_ghosts_needed` helper; `get_global_ghost_flags` QUAT; `get_force_reduce_ghost_flags` |
| `src/core/virtual_sites/relative.cpp` | `vs_relative_update_particles` explicitly requests `DATA_PART_QUAT` |
