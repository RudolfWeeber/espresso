# Pressure Tensor Kokkos/Cabana Migration Design

## Goal

Port the pressure/stress tensor calculation from the CPU `short_range_loop()` to Kokkos/Cabana, analogous to the existing energy calculation port. The design mirrors the energy path's architecture (bin layout + per-thread 2D view + reduce) while accounting for key differences: tensor (9-element) vs scalar output, DPD stress contribution, no rotational kinetic contribution, and no external fields contribution.

## Approach

**Flat scalar bins with 9x offset per category** (Approach A from brainstorming). The `local_pressure` view uses the same `Kokkos::View<double**, LayoutRight>` pattern as energy, but each category occupies 9 consecutive elements (one flattened 3x3 tensor). Index functions return the first of 9 elements; kernels write `bin + k*3 + l` for tensor element (k,l).

## File Structure

### New Files

| File | Analog | Purpose |
|------|--------|---------|
| `src/core/pressure_cabana.hpp` | `energy_cabana.hpp` | `PressureBinLayout`, `PressureKernel`, `KineticPressureKernel`, `reduce_cabana_pressure()` |
| `src/core/bond_pressure_kokkos.hpp` | `bond_energy_kokkos.hpp` | `BondsPressureKernelData`, `PairBondsPressureKernel`, `AngleBondsPressureKernel`, `DihedralBondsPressureKernel` |

### Modified Files

| File | Change |
|------|--------|
| `src/core/pressure.cpp` | Replace `short_range_loop()` with `cabana_short_range()` + Kokkos kinetic loop; remove DPD CPU call (now in-kernel); keep CPU calls for virtual sites, Coulomb k-space, dipolar warning |

## PressureBinLayout

Mirrors `EnergyBinLayout` with these differences:

- Each category occupies 9 consecutive flat indices (one flattened 3x3 tensor)
- **No rotational kinetic bin** (pressure has no rotational contribution)
- **DPD bin added** (energy has no DPD contribution; pressure has dissipative stress)
- **No external_fields bin** (energy has constraint energy; pressure has no constraint pressure)
- **No dipolar energy analog** — dipolar bin exists but stays zero (with stderr warning, matching current behavior)

### Flat Index Layout

```
[0 .. 9*n_bonded)              bonded bins (9 per bond_id)
[9*n_bonded .. +9*n_nb)        nb_inter bins (9 per type pair)
[.. +9*n_nb)                   nb_intra bins (9 per type pair)
[.. +9)                        coulomb real-space (1 tensor)
[.. +9)                        dipolar real-space (1 tensor, stays zero)
[.. +9)                        DPD stress (1 tensor)
[.. +9)                        kinetic linear (1 tensor)
```

Where `n_nb = n_types * (n_types + 1) / 2` (lower triangular type-pair count), same as energy.

### Index Functions

Each returns the *first* of 9 flat indices for that category:

- `bonded_idx(b)` → `off_bonded + 9 * b`
- `nb_inter_idx(t1, t2)` → `off_nb_inter + 9 * lower_triangular(max(t1,t2), min(t1,t2))`
- `nb_intra_idx(t1, t2)` → `off_nb_intra + 9 * lower_triangular(max(t1,t2), min(t1,t2))`
- `coulomb_idx()` → `off_coulomb`
- `dipolar_idx()` → `off_dipolar`
- `dpd_idx()` → `off_dpd`
- `kinetic_idx()` → `off_kinetic`

### View

`Kokkos::View<double**, Kokkos::LayoutRight> local_pressure("local_pressure", n_threads, layout.total)`

Same 2D per-thread layout as `local_energy`, but wider due to 9× bins.

## PressureKernel (Non-Bonded Pair Pressure)

`PressureKernel::operator()(i, j)` mirrors `EnergyKernel::operator()(i, j)` structure but computes the full virial tensor.

### Compute Flow

1. Get positions from AoSoA, compute `d = box_geo.get_mi_vector(pos1, pos2)`, `dist = d.norm()`, early exit if `dist > system_max_cutoff`
2. Exclusion/Thole flag check (same pattern as energy/forces)
3. If `dist <= ia_params.max_cut` (after exclusion check):
   - `f = calc_central_radial_force(ia_params, d, dist)`
   - Thole force: `f += thole_pair_force(...)` if active
   - Gay-Berne force: `f += gb_pair_force(...)` if active (includes non-central contributions)
   - `stress = Utils::tensor_product(d, f)` → 9 elements written to nb_inter/nb_intra bin based on `mol_id_view(i)` vs `mol_id_view(j)`
4. **DPD** (dissipative only, no noise): call the inner `dpd_pair_force(params, v21, dist, zero_noise)` with `zero_noise = Utils::Vector3d{}` to get deterministic dissipative force only (no stochastic term). Compute `stress_dpd = -Utils::tensor_product(d, f_dpd)` → DPD bin. Requires velocity data from AoSoA and DPD thermostat parameters.
5. **Coulomb real-space**: `kernel_pressure(q1q2, d, dist)` returns a `Utils::Matrix<double,3,3>` → write 9 elements to coulomb bin
6. **Dipolar**: stderr warning only, no tensor written

### Key Differences from EnergyKernel

| Aspect | EnergyKernel | PressureKernel |
|--------|-------------|----------------|
| Core computation | `calc_central_radial_energy()` (scalar) | `calc_central_radial_force()` → `tensor_product(d, f)` (3x3) |
| DPD | Not present | Dissipative stress tensor included |
| Coulomb kernel | `pair_energy_kernel()` (scalar) | `pair_pressure_kernel()` (3x3 matrix) |
| Bin write | 1 element per bin | 9 elements per bin |
| Rotational | Not in non-bonded | Not applicable |
| Dipolar | Energy value | Warning only |

### Constructor Parameters

Same captures as `EnergyKernel`, plus:
- `Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel` (pressure kernel, not energy)
- `Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel` (needed for Thole)
- `Dipoles::Solver const &dipoles` (for warning only)
- `Thermostat::Thermostat const &thermostat` (needed for DPD)
- `DPDThermostat const *dpd` (needed for DPD noise-free computation)
- `int thermo_switch` (needed for `dpd_active()` check)

## Bonded Pressure Kernels

### BondsPressureKernelData

```cpp
struct BondsPressureKernelData {
  BondedInteractionsMap const &bonded_ias;
  BoxGeometry const &box_geo;
  Kokkos::View<double**, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
};
```

### PairBondsPressureKernel

1. Resolve AoSoA indices from bond list, get positions, compute `dx = box_geo.get_mi_vector(pos1, pos2)`
2. Call `calc_bond_pair_force(iaparams, dx, q1q2, coulomb_force_kernel)` → `optional<Utils::Vector3d>`
3. If force exists: `pressure = Utils::tensor_product(*pair_force, dx)` → write 9 elements to `bonded_idx(bond_id)` bin
4. If bond broken: `bond_broken_error()`

### AngleBondsPressureKernel

1. Resolve 3 particle indices, compute `dx21 = -box_geo.get_mi_vector(pos1, pos2)`, `dx31 = box_geo.get_mi_vector(pos3, pos1)`
2. Call `calc_bonded_three_body_force(iaparams, dx21, dx31)` → `optional<tuple<force1, force2, force3>>`
3. If force exists: `pressure = Utils::tensor_product(force2, dx21) + Utils::tensor_product(force3, dx31)` → write 9 elements to bonded bin
4. If bond broken: `bond_broken_error()`

### DihedralBondsPressureKernel

1. Resolve 4 particle indices, get positions
2. Write zero 3x3 tensor (9 zero elements) to `bonded_idx(bond_id)` bin
3. Log runtime warning for unsupported bond type (matching current CPU behavior)

## KineticPressureKernel

New kernel that computes kinetic pressure via Kokkos parallel loop (unlike energy, which still uses a serial CPU loop for kinetic).

### KineticPressureKernel::operator()(idx)

1. Get particle index from AoSoA
2. Skip if virtual (check `aosoa.flags`)
3. Read `v[k]` for k=0,1,2 and `mass` from AoSoA
4. Compute `m * v[k] * v[l]` for all k,l → write 9 elements to `layout.kinetic_idx()` bin

Uses `Kokkos::parallel_for` over `n_particles` range. No rotational contribution (matching current pressure behavior).

## reduce_cabana_pressure()

Mirrors `reduce_cabana_energy()`:

1. `Kokkos::create_mirror_view_and_copy` to host
2. Sum across threads for each bin via `sum_bin(bin_index)` lambda
3. Write into `Observable_stat` spans (chunk_size=9):

**Bonded**: For each bond_id `b`:
```
obs_pressure.bonded_contribution(b)[k*3+l] += sum_bin(layout.bonded_idx(b) + k*3+l)
```
for k,l in 0..2.

**nb_inter/nb_intra**: For each type pair (t1,t2):
```
obs_pressure.non_bonded_inter_contribution(t1,t2)[k*3+l] += sum_bin(layout.nb_inter_idx(t1,t2) + k*3+l)
obs_pressure.non_bonded_intra_contribution(t1,t2)[k*3+l] += sum_bin(layout.nb_intra_idx(t1,t2) + k*3+l)
```

**Coulomb real-space**: `obs_pressure.coulomb[k*3+l] += sum_bin(layout.coulomb_idx() + k*3+l)`

**Dipolar**: `obs_pressure.dipolar[k*3+l] += sum_bin(layout.dipolar_idx() + k*3+l)` (will be zero)

**DPD**: `obs_pressure.dpd[k*3+l] += sum_bin(layout.dpd_idx() + k*3+l)`

**Kinetic**: `obs_pressure.kinetic_lin[k*3+l] += sum_bin(layout.kinetic_idx() + k*3+l)`

## pressure.cpp New Flow

1. Create `Observable_stat(9ul, n_bonded, max_type)` — unchanged
2. Sanity checks, `on_observable_calc()` — unchanged
3. Get Coulomb kernels: `pair_pressure_kernel()`, `pair_force_kernel()` — both needed (pressure kernel for Coulomb pressure, force kernel for Thole/bonded)
4. `update_cabana_state()` — same call as energy
5. Allocate `PressureBinLayout` + `local_pressure` view
6. Build `mol_id_view` — same pattern as energy
7. `Kokkos::parallel_for` with `KineticPressureKernel` over n_particles
8. Construct `PressureKernel`, bond pressure kernels
9. `cabana_short_range(pair_bond_p_kernel, angle_bond_p_kernel, dih_bond_p_kernel, pair_p_kernel, ...)`
10. `reduce_cabana_pressure()` into `Observable_stat`
11. Coulomb k-space: `coulomb.calc_pressure_long_range()` → `obs_pressure.coulomb[9..17]` (CPU, unchanged)
12. Dipolar warning (CPU, unchanged)
13. Virtual sites: `vs_relative_pressure_tensor()` (CPU, unchanged)
14. ~~DPD CPU call removed~~ — now computed inside `PressureKernel`
15. `obs_pressure.rescale(volume)`, `obs_pressure.mpi_reduce()` — unchanged

## Summary of Energy vs Pressure Differences

| Aspect | Energy (chunk=1) | Pressure (chunk=9) |
|--------|-----------------|-------------------|
| Bin size | 1 scalar | 9 elements (3x3 tensor) |
| Kinetic | `0.5*m*v^2` + rotational | `m*v[k]*v[l]` only, no rotational |
| Non-bonded core | `calc_central_radial_energy()` | `calc_central_radial_force()` → `tensor_product(d,f)` |
| DPD | Not present | Dissipative stress in-kernel |
| Coulomb real-space | `pair_energy_kernel()` (scalar) | `pair_pressure_kernel()` (3x3) |
| Coulomb k-space | Scalar energy | 9-element tensor (CPU, unchanged) |
| External fields | Constraint energy | Not computed (stays zero) |
| Virtual sites | Not in energy | CPU call (unchanged) |
| Dipolar | Energy value | Warning only, zero tensor |
| Dihedral bonds | Actual energy | Zero tensor + warning |
