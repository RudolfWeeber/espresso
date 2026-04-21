# Design: Pressure/Stress Tensor Calculation via Kokkos/Cabana

**Date:** 2026-04-21
**Branch:** perf_2026_04

## Goal

Convert `System::calculate_pressure()` to use the same Kokkos/Cabana short-range loop
infrastructure as `System::calculate_energy()`, replacing the legacy `short_range_loop`
call and the separate `dpd_pressure_local()` call.

## Key Differences from the Energy Conversion

| Aspect | Energy | Pressure |
|---|---|---|
| Values per bin | 1 scalar | 9 doubles (flat 3×3 tensor) |
| Kinetic contribution | `kinetic_lin` + `kinetic_rot` | `kinetic_lin` only (no rotation analogue) |
| DPD contribution | none | yes — dissipative virial, own bin |
| Coulomb kernel | `pair_energy_kernel` | `pair_pressure_kernel` (returns 3×3) |
| Dipoles | implemented | not implemented (warning only, unchanged) |
| Dihedrals bonded | implemented | not implemented (warning, unchanged) |

## New Files

### `src/core/pressure_cabana.hpp`

Mirrors `energy_cabana.hpp`.

**`PressureBinLayout`**

Identical offset arithmetic to `EnergyBinLayout` (`off_bonded`, `off_nb_inter`,
`off_nb_intra`, `off_coulomb`, `off_dipolar`, `off_dpd`, `total`). Adds DPD bin:

```
off_dpd  = off_dipolar + 1
total    = off_dpd + 1     (bin units; multiply by 9 for flat view index)
```

Key accessor:
```cpp
KOKKOS_INLINE_FUNCTION
std::size_t tensor_offset(std::size_t bin, std::size_t k) const {
  return bin * 9 + k;
}
```

**`PressureKernel`**

Same constructor fields as `EnergyKernel`, replacing the Coulomb energy kernel with:
- `Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel`

No DPD thermostat pointer needed — the DPD dissipative force uses only `ia_params.dpd.radial`
and `ia_params.dpd.trans`, which are already accessible via `nonbonded_ias` inside the kernel.

`operator()(i, j)` computes three contributions:

1. **Non-bonded virial** (inside `ia_params.max_cut`):
   ```
   f  = calc_central_radial_force(ia_params, d, dist)
   f += calc_non_central_force(dir1, dir2, ia_params, d, dist)  // director overload
   #ifdef ESPRESSO_THOLE
     f += thole_pair_force(...)  // via unique_particles fallback
   #endif
   stress = tensor_product(d, f)  → 9 values → nb_inter or nb_intra bin
   ```

2. **Coulomb virial** (`#ifdef ESPRESSO_ELECTROSTATICS`):
   ```
   p_coulomb = (*coulomb_p_kernel)(q1*q2, d, dist)  // 3×3 matrix
   → 9 values → coulomb bin
   ```

3. **DPD virial** (`#ifdef ESPRESSO_DPD`, when `dpd != nullptr`):
   ```
   vel1 = aosoa.get_vector_at(aosoa.velocity, i)
   vel2 = aosoa.get_vector_at(aosoa.velocity, j)
   v21  = box_geo.velocity_difference(pos1, pos2, vel1, vel2)
   f_r  = dpd_pair_force(ia_params.dpd.radial, v21, dist, {})
   f_t  = dpd_pair_force(ia_params.dpd.trans,  v21, dist, {})
   P    = tensor_product(d / dist2, d)
   f    = P*(f_r - f_t) + f_t
   stress_dpd = tensor_product(d, f)  → 9 values → dpd bin
   ```

**`reduce_cabana_pressure`**

Mirrors `reduce_cabana_energy`. Copies host mirror, iterates bins × 9 components,
writes into `Observable_stat` via existing setters:
- `obs.non_bonded_inter_contribution(t1,t2)[k]`
- `obs.non_bonded_intra_contribution(t1,t2)[k]`
- `obs.bonded_contribution(b)[k]`
- `obs.coulomb[k]`
- `obs.dipolar[k]`
- `obs.dpd[k]`

### `src/core/bond_pressure_kokkos.hpp`

Mirrors `bond_energy_kokkos.hpp`.

**`BondsPressureKernelData`** — same fields as `BondsEnergyKernelData` with `local_pressure`.

**`PairBondsPressureKernel`**

Calls `calc_bonded_virial_pressure_tensor(iaparams, pos1, pos2, box_geo, coulomb_kernel)`,
flattens the resulting `Utils::Matrix<double,3,3>` into 9 values at `layout.tensor_offset(bonded_idx(bond_id), k)`.

**`AngleBondsPressureKernel`**

Calls `calc_bonded_three_body_pressure_tensor(iaparams, pos1, pos2, pos3, box_geo)`,
stores 9 values into bonded bin.

**No dihedral kernel** — matches current behavior (dihedrals silently contribute zero).

## Changes to Existing Files

### `src/core/pressure.cpp`

`calculate_pressure()` restructured to match `calculate_energy()`:

```
1. Kinetic virial
   reduce_over_local_particles<Utils::Matrix<double,3,3>>(
       *cell_structure,
       [](auto &acc, Particle const &p) {
           if (!p.is_virtual())
               acc += Utils::tensor_product(p.v(), p.mass() * p.v());
       },
       [](auto &a, auto const &b) { a = a + b; });
   → flatten into obs_pressure.kinetic_lin[0..8]

2. update_cabana_state(...)  // same verlet_criterion as energy

3. PressureBinLayout layout{...}
   Kokkos::View<double**, LayoutRight> local_pressure(
       "local_pressure", exec().concurrency(), layout.total * 9);

4. mol_id_view (same pattern as energy)

5. Construct PressureKernel, PairBondsPressureKernel, AngleBondsPressureKernel

6. cabana_short_range(pair_bp_kernel, angle_bp_kernel,
                      NullDihedralKernel{}, pair_p_kernel, ...)

7. reduce_cabana_pressure(local_pressure, layout, obs_pressure, ...)

8. Long-range k-space (coulomb, dipoles) — unchanged
9. Virtual sites — unchanged
10. rescale + mpi_reduce
```

Removed:
- `short_range_loop(...)` call
- `dpd_pressure_local(...)` call (DPD now inside `PressureKernel`)
- Plain `for (auto const &p : local_parts) { add_kinetic_virials(...) }` loop

### `src/core/pressure_inline.hpp`

- `add_kinetic_virials()` removed (dead code after above change)
- `add_non_bonded_pair_virials()` removed (dead code — `dpd_stress()` uses
  `cell_structure.non_bonded_loop()` directly, not this function)
- Add `calc_non_central_force(dir1, dir2, ia_params, d, dist)` overload taking
  directors instead of `Particle const &` to support Gay-Berne in Cabana kernel
  without `unique_particles` fallback

### `src/core/forces_inline.hpp`

Add director-based overload of `calc_non_central_force` (or move to `pressure_inline.hpp`
if Gay-Berne is the only non-central force relevant to pressure).

## Gay-Berne Handling

The energy kernel already accesses Gay-Berne via `aosoa.director` (director vector).
For pressure, `calc_non_central_force` currently takes `Particle const &` to get
`p.quat()`. A new overload accepting `(Utils::Vector3d dir1, Utils::Vector3d dir2, ...)`
avoids the `unique_particles` fallback for Gay-Berne, keeping it on the hot path.
Thole still uses the `unique_particles` fallback (same `flag.need_particle_pointers`
guard as the energy kernel).

## What Is Not Changed

- Long-range (k-space) Coulomb and dipole pressure — unchanged
- Virtual sites pressure — unchanged
- `dpd_stress()` (the full stress tensor observable, separate from `calculate_pressure`) — unchanged
- `particle_short_range_energy_contribution()` — energy only, unaffected
