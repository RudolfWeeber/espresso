/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#include "cells.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "constraints/Constraints.hpp"
#include "electrostatics/icc.hpp"
#include "forces_inline.hpp"
#include "galilei/ComFixed.hpp"
#include "immersed_boundary/ImmersedBoundaries.hpp"
#include "integrators/Propagation.hpp"
#include "lb/particle_coupling.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rotation.hpp"
#include "short_range_cabana.hpp"
#include "short_range_loop.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"
#include "virtual_sites/com.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

#include <Cabana_Core.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <span>
#include <variant>

#ifdef ESPRESSO_ENGINE
/** @brief Add the swimming (engine) force to a particle's force.
 *
 *  Split out of the (now vectorized) external-force column sweep: the swim
 *  term needs the quaternion-derived director, so it stays a per-row scalar
 *  update applied ONLY to active swimmers. Preserves the original combination
 *  order: force = ext_force (column copy) then force += f_swim * director,
 *  which equals the previous @c ext_force + swim for every swimmer.
 */
static void add_swim_force(Particle &p) {
  if (p.swimming().swimming and !p.swimming().is_engine_force_on_fluid) {
    auto force = p.force();
    force += p.swimming().f_swim * p.calc_director();
  }
}
#endif

/** @brief Copy an external-force component column into a force column.
 *
 *  Component-major (@ref ParticleStore::StateVectorLayout LayoutLeft) columns
 *  make each component a contiguous stride-1 run across particles, so
 *  initializing @c force[c] = @c ext_force[c] over @c [0, n_local) is a plain
 *  contiguous copy that vectorizes (packed 4-wide double moves under AVX2)
 *  instead of the per-row @ref Particle proxy write it replaces. The
 *  @c #pragma omp simd asserts independence for the auto-vectorizer.
 *
 *  @tparam DstView  destination component-major view (@c double*[3]).
 *  @tparam SrcView  source component-major view (@c double*[3]).
 */
template <class DstView, class SrcView>
static void init_component_copy(DstView const &dst, SrcView const &src,
                                std::size_t const n_local) {
  for (unsigned int c = 0u; c < 3u; ++c) {
#pragma omp simd
    for (std::size_t row = 0u; row < n_local; ++row) {
      dst(row, c) = src(row, c);
    }
  }
}

/** @brief Zero a force component column over @c [0, n_local).
 *
 *  Used when the external-force column is absent (external forces disabled at
 *  compile time): initialize @c force[c] = 0 as a contiguous memset-style
 *  vectorizable loop, matching the previous @c force = @c {} per-row write.
 */
template <class DstView>
static void init_component_zero(DstView const &dst, std::size_t const n_local) {
  for (unsigned int c = 0u; c < 3u; ++c) {
#pragma omp simd
    for (std::size_t row = 0u; row < n_local; ++row) {
      dst(row, c) = 0.;
    }
  }
}

/** @brief Apply the Langevin friction+noise contribution to one particle.
 *
 *  Per-row body of the Langevin sweep, factored into its own force-inlined
 *  helper. ESPRESSO_ATTR_ALWAYS_INLINE (NOT [[gnu::flatten]]) recovers the
 *  original single-pass lambda's inlining of the trivial helpers (Utils::Vector
 *  ctors, Utils::hadamard_product / operator+=, the body->space rotation)
 *  WITHOUT the whole-tree flatten: flattening this path under -march=native
 *  fuses pref*v into the force accumulate (an FMA contraction the original
 *  flattened lambda did not hit), which alters trajectories at the rounding
 *  level and breaks the WS3 bitwise-identity requirement. always_inline keeps
 *  the friction arithmetic in the same shape as before, so it stays
 *  bitwise-identical (verified by the same-flags Langevin identity gate).
 *  Per-particle order is unchanged: force += (friction + noise),
 *  torque += rot(friction + noise).
 */
template <class Store, class Propagation, class VelView, class IdView
#ifdef ESPRESSO_ROTATION
          ,
          class OmegaView
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
          ,
          class GammaView
#ifdef ESPRESSO_ROTATION
          ,
          class GammaRotView
#endif
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
          ,
          class QuatView
#endif
          >
ESPRESSO_ATTR_ALWAYS_INLINE inline void
apply_langevin_row(Store &store, Propagation const &propagation,
                   LangevinThermostat const &langevin, VelView const &vel_view,
                   IdView const &id_view,
#ifdef ESPRESSO_ROTATION
                   OmegaView const &omega_view,
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                   GammaView const &gamma_view,
#ifdef ESPRESSO_ROTATION
                   GammaRotView const &gamma_rot_view,
#endif
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
                   QuatView const &quat_view,
#endif
                   int const row, double const time_step, double const kT) {
  // Read the propagation bitfield once per particle (each
  // should_propagate_with call would otherwise re-read the store column).
  int const prop = store.propagation_view()(row);
  if (propagation.should_propagate_with(prop,
                                        PropagationMode::TRANS_LANGEVIN)) {
    auto const friction = friction_thermo_langevin(langevin, vel_view, id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                                                   gamma_view,
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
                                                   quat_view,
#endif
                                                   row, time_step, kT);
    auto force = store.force_reference(row);
    force += friction;
  }
#ifdef ESPRESSO_ROTATION
  if (propagation.should_propagate_with(prop, PropagationMode::ROT_LANGEVIN)) {
    Particle p;
    p.attach_to_store(store, row);
    auto torque = p.torque();
    torque += convert_vector_body_to_space(
        p, friction_thermo_langevin_rotation(langevin, omega_view, id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                                             gamma_rot_view,
#endif
                                             row, time_step, kT));
  }
#endif
}

/** Combined force initialization and Langevin noise application.
 *
 *  Restructured into contiguous component sweeps under the component-major
 *  columns (SIMD WS3):
 *   1. force/torque INIT sweep -- three (x,y,z) contiguous column copies
 *      (force = ext_force, torque = ext_torque), each trivially vectorizable,
 *      replacing the per-row proxy writes. The engine swim term needs the
 *      quaternion-derived director, so it stays a per-row scalar update applied
 *      only to active swimmers (detected by a cheap pre-scan of the swimming
 *      sidecar); non-engine builds skip it entirely.
 *   2. Langevin sweep -- the per-particle id-keyed Philox noise draw is scalar
 *      by design; the friction+noise contribution is added per row EXACTLY as
 *      before. Per-particle arithmetic order is preserved: force = ext_force
 *      (sweep 1) then force += (friction + noise) (sweep 2) == the previous
 *      force = ext_force; force += (friction + noise). No reordering within a
 *      particle's force; only the two passes over the (same) row set are split.
 *
 *  NOTE on [[gnu::flatten]]: the original single-pass version flattened the
 *  whole lambda to keep the Langevin helper tree inline. It is NOT applied to
 *  THIS outer function -- flattening it whole would (a) inhibit nothing useful
 *  and (b) fuse pref*v into the force accumulate (an FMA contraction the
 *  original flattened lambda did not hit) and break the WS3 bitwise-identity
 *  requirement. Flatten is instead scoped to @ref apply_langevin_row, which
 *  reproduces the original inlining on the friction path (bitwise-identical,
 *  verified by the same-flags Langevin identity gate) while leaving the init
 *  sweeps here un-flattened so they auto-vectorize into packed column copies.
 */
static void init_forces_and_thermostat(System::System const &system) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  auto &cell_structure = *system.cell_structure;
  auto const &propagation = *system.propagation;
  auto const &thermostat = *system.thermostat;
  auto const kT = thermostat.kT;
  auto const time_step = system.get_time_step();

  // Check if Langevin thermostat is active
  bool const langevin_active =
      thermostat.langevin &&
      (propagation.used_propagations &
       (PropagationMode::TRANS_LANGEVIN | PropagationMode::ROT_LANGEVIN));

  auto &store = cell_structure.particle_store();
  // Local particles tile store rows [0, n_local) contiguously in
  // cell-traversal order (see ensure_particle_store_synchronized), so the
  // component sweeps below run straight over [0, n_local) -- exactly the row
  // set for_each_local_particle_row visits, in the same order.
  auto const n_local = store.number_of_local_particles();

  // -- Sweep 1: contiguous force/torque initialization ----------------------
  // force[c] = ext_force[c] (or 0 when external forces are disabled).
  {
    auto force_view = store.force_view();
#ifdef ESPRESSO_EXTERNAL_FORCES
    init_component_copy(force_view, store.ext_force_view(), n_local);
#else
    init_component_zero(force_view, n_local);
#endif
#ifdef ESPRESSO_ROTATION
    auto torque_view = store.torque_view();
#if defined(ESPRESSO_EXTERNAL_FORCES)
    init_component_copy(torque_view, store.ext_torque_view(), n_local);
#else
    init_component_zero(torque_view, n_local);
#endif
#endif // ESPRESSO_ROTATION
  }

#ifdef ESPRESSO_ENGINE
  // Engine swim force: needs the per-particle director (non-vectorizable), so
  // apply it per row ONLY when at least one active swimmer exists. The pre-scan
  // is a contiguous bool read over the swimming sidecar; the benchmarks carry
  // no swimmers, so this short-circuits without touching the director path.
  {
    bool any_swimmer = false;
    for (std::size_t row = 0u; row < n_local; ++row) {
      auto const &sw = store.swimming(static_cast<int>(row));
      if (sw.swimming and !sw.is_engine_force_on_fluid) {
        any_swimmer = true;
        break;
      }
    }
    if (any_swimmer) {
      cell_structure.for_each_local_particle_row([&](int const row) {
        Particle p;
        p.attach_to_store(store, row);
        add_swim_force(p);
      });
    }
  }
#endif // ESPRESSO_ENGINE

  // -- Sweep 2: Langevin friction + noise ------------------------------------
  if (langevin_active) {
    // Hoist the Langevin column-view handles ONCE outside the parallel_for.
    // The friction is a column kernel (velocity / omega / id (+ gamma,
    // quaternion) read by row); the body->space torque rotation needs the
    // quaternion, so it stays on the per-row view path.
    auto vel_view = store.velocity_view();
    auto id_view = store.id_view();
#ifdef ESPRESSO_ROTATION
    auto omega_view = store.angular_velocity_view();
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    auto gamma_view = store.gamma_view();
#ifdef ESPRESSO_ROTATION
    auto gamma_rot_view = store.gamma_rot_view();
#endif
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
    auto quat_view = store.quaternion_view();
#endif
    auto const &langevin = *thermostat.langevin;

    cell_structure.for_each_local_particle_row([&](int const row) {
      apply_langevin_row(store, propagation, langevin, vel_view, id_view,
#ifdef ESPRESSO_ROTATION
                         omega_view,
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                         gamma_view,
#ifdef ESPRESSO_ROTATION
                         gamma_rot_view,
#endif
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
                         quat_view,
#endif
                         row, time_step, kT);
    });
  }

  cell_structure.reset_local_force_and_torque();

  // Initialize ghost forces (unchanged)
  cell_structure.ghosts_reset_forces();
}

static void force_capping(CellStructure &cell_structure, double force_cap) {
  if (force_cap > 0.) {
    auto const force_cap_sq = Utils::sqr(force_cap);
    cell_structure.for_each_local_particle(
        [&force_cap, &force_cap_sq](Particle &p) {
          auto force = p.force();
          auto const force_sq = force.norm2();
          if (force_sq > force_cap_sq) {
            force *= force_cap / std::sqrt(force_sq);
          }
        });
  }
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
static void reinit_dip_fld(CellStructure const &cell_structure) {
  cell_structure.for_each_local_particle(
      [](Particle &p) { p.dip_fld() = {0., 0., 0.}; });
}
#endif

static BondsKernelData
create_kokkos_bonds_kernel_data(System::System const &system) {
  auto scatter_force = system.cell_structure->get_scatter_force();
#ifdef ESPRESSO_NPT
  auto scatter_virial = system.cell_structure->get_scatter_virial();
#endif
  auto const &aosoa = system.cell_structure->get_aosoa();
  return /* BondsKernelData */ {*system.bonded_ias,
                                *system.bond_breakage,
                                *system.box_geo,
                                scatter_force,
#ifdef ESPRESSO_NPT
                                scatter_virial,
#endif
                                aosoa,
                                !system.bond_breakage->breakage_specs.empty()};
}

static ForcesKernel create_cabana_neighbor_kernel(
    System::System const &system, Utils::Vector3d *virial,
    auto const &elc_kernel, auto const &coulomb_kernel,
    auto const &dipoles_kernel, auto const &coulomb_u_kernel) {

  auto const &unique_particles = system.cell_structure->get_unique_particles();

  auto scatter_force = system.cell_structure->get_scatter_force();
#ifdef ESPRESSO_ROTATION
  auto scatter_torque = system.cell_structure->get_scatter_torque();
#endif
#ifdef ESPRESSO_NPT
  auto scatter_virial = system.cell_structure->get_scatter_virial();
#endif
  auto const &aosoa = system.cell_structure->get_aosoa();

  return /* ForcesKernel */ {*system.bonded_ias,
                             *system.nonbonded_ias,
                             get_ptr(coulomb_kernel),
                             get_ptr(dipoles_kernel),
                             get_ptr(elc_kernel),
                             get_ptr(coulomb_u_kernel),
                             system.coulomb,
                             *system.thermostat,
                             *system.box_geo,
                             unique_particles,
                             scatter_force,
#ifdef ESPRESSO_ROTATION
                             scatter_torque,
#endif
#ifdef ESPRESSO_NPT
                             virial,
                             scatter_virial,
#endif
                             aosoa,
                             system.maximal_cutoff()};
}

static void reduce_cabana_forces_and_torques(System::System const &system,
                                             Utils::Vector3d *virial) {

  auto const &unique_particles = system.cell_structure->get_unique_particles();
  auto &local_force = system.cell_structure->get_local_force();
  auto scatter_force = system.cell_structure->get_scatter_force();
  Kokkos::Experimental::contribute(local_force, scatter_force);
#ifdef ESPRESSO_ROTATION
  auto &local_torque = system.cell_structure->get_local_torque();
  auto scatter_torque = system.cell_structure->get_scatter_torque();
  Kokkos::Experimental::contribute(local_torque, scatter_torque);
#endif
#ifdef ESPRESSO_NPT
  auto &local_virial = system.cell_structure->get_local_virial();
  auto scatter_virial = system.cell_structure->get_scatter_virial();
  Kokkos::Experimental::contribute(local_virial, scatter_virial);
#endif

  using execution_space = Kokkos::DefaultExecutionSpace;
  Kokkos::RangePolicy<execution_space> policy(std::size_t{0},
                                              unique_particles.size());
  Kokkos::parallel_for("reduction", policy,
                       [&local_force,
#ifdef ESPRESSO_ROTATION
                        &local_torque,
#endif
                        &unique_particles](std::size_t const i) {
                         Utils::Vector3d force{};
#ifdef ESPRESSO_ROTATION
                         Utils::Vector3d torque{};
#endif
                         force[0] += local_force(i, 0);
                         force[1] += local_force(i, 1);
                         force[2] += local_force(i, 2);
#ifdef ESPRESSO_ROTATION
                         torque[0] += local_torque(i, 0);
                         torque[1] += local_torque(i, 1);
                         torque[2] += local_torque(i, 2);
#endif
                         auto *particle = unique_particles[i];
                         particle->force() += force;
#ifdef ESPRESSO_ROTATION
                         particle->torque() += torque;
#endif
                       });
  Kokkos::fence();

#ifdef ESPRESSO_NPT
  if (virial) {
    (*virial)[0] += local_virial(0);
    (*virial)[1] += local_virial(1);
    (*virial)[2] += local_virial(2);
  }
#endif
}

void System::System::calculate_forces() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  // Ensure every local/ghost particle has a valid ParticleStore row before any
  // force/torque access below. O(1) when the store is clean; rank-local.
  cell_structure->ensure_particle_store_synchronized();
#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("copy_particles_to_GPU");
#endif
    gpu->update();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("copy_particles_to_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_COLLISION_DETECTION
  collision_detection->clear_queue();
  auto const collision_detection_cutoff = collision_detection->cutoff();
#else
  auto const collision_detection_cutoff = inactive_cutoff;
#endif
  bond_breakage->clear_queue();
  auto particles = cell_structure->local_particles();
#ifdef ESPRESSO_NPT
  if (propagation->used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) {
    // reset virial part of instantaneous pressure
    npt_inst_pressure->p_vir = Utils::Vector3d{};
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  // reset dipole field
  reinit_dip_fld(*cell_structure);
#endif

  // Use combined function instead of two separate calls

  auto const elc_kernel = coulomb.pair_force_elc_kernel();
  auto const coulomb_kernel = coulomb.pair_force_kernel();
  auto const dipoles_kernel = dipoles.pair_force_kernel();
  auto const coulomb_u_kernel = coulomb.pair_energy_kernel();
  auto *const virial = get_npt_virial();

  VerletCriterion<> const verlet_criterion{*this,
                                           cell_structure->get_verlet_skin(),
                                           get_interaction_range(),
                                           coulomb.cutoff(),
                                           dipoles.cutoff(),
                                           collision_detection_cutoff};

  update_cabana_state(*cell_structure, verlet_criterion,
                      get_interaction_range(), propagation->integ_switch);
#ifdef ESPRESSO_ELECTROSTATICS
  // Refresh the pack-owned charge column once per step, ONLY when a coulomb
  // actor is active (both the P3M long-range gather below and the real-space
  // pair kernel read it contiguously). Pure-LJ runs skip it.
  if (coulomb.impl->solver) {
    refresh_pack_charges(*cell_structure);
  }
  if (coulomb.impl->extension) {
    update_icc_particles();
  }
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_DIPOLES
  // Refresh the pack-owned dipm column, guarded by an active dipolar actor.
  if (dipoles.impl->solver) {
    refresh_pack_dipm(*cell_structure);
  }
#endif // ESPRESSO_DIPOLES
  init_forces_and_thermostat(*this);
#ifdef ESPRESSO_CALIPER
  CALI_MARK_BEGIN("calc_long_range_forces");
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  coulomb.calc_long_range_force();
#endif
#ifdef ESPRESSO_DIPOLES
  dipoles.calc_long_range_force();
#endif
#ifdef ESPRESSO_CALIPER
  CALI_MARK_END("calc_long_range_forces");
#endif

#ifdef ESPRESSO_CALIPER
  CALI_MARK_BEGIN("cabana_short_range");
#endif
  auto &bs = cell_structure->bond_state();
  auto bonds_kernel_data = create_kokkos_bonds_kernel_data(*this);
  auto pair_bonds_kernel = PairBondsKernel{
      bonds_kernel_data, bs.pair_list, bs.pair_ids, get_ptr(coulomb_kernel)};
  auto angle_bonds_kernel =
      AngleBondsKernel{bonds_kernel_data, bs.angle_list, bs.angle_ids};
  auto dihedral_bonds_kernel =
      DihedralBondsKernel{bonds_kernel_data, bs.dihedral_list, bs.dihedral_ids};

  auto first_neighbor_kernel =
      create_cabana_neighbor_kernel(*this, virial, elc_kernel, coulomb_kernel,
                                    dipoles_kernel, coulomb_u_kernel);

  cabana_short_range(pair_bonds_kernel, angle_bonds_kernel,
                     dihedral_bonds_kernel, first_neighbor_kernel,
                     *cell_structure, get_interaction_range(),
                     bonded_ias->maximal_cutoff(), verlet_criterion,
                     propagation->integ_switch);

  // Force and Torque reduction
  reduce_cabana_forces_and_torques(*this, virial);

#ifdef ESPRESSO_COLLISION_DETECTION
  auto collision_kernel = [&collision_detection = *collision_detection](
                              Particle const &p1, Particle const &p2,
                              Distance const &d) {
    collision_detection.detect_collision(p1, p2, d.dist2);
  };
  if (not collision_detection->is_off()) {
    cell_structure->non_bonded_loop(collision_kernel, verlet_criterion);
  }
#endif // ESPRESSO_COLLISION_DETECTION

#ifdef ESPRESSO_CALIPER
  CALI_MARK_END("cabana_short_range");
#endif

  constraints->add_forces(particles, get_sim_time());
  oif_global->calculate_forces();

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries->volume_conservation(*cell_structure);

  if (thermostat->lb and (propagation->used_propagations &
                          PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("lb_particle_coupling");
#endif
    lb_couple_particles();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("lb_particle_coupling");
#endif
  }

#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("copy_forces_from_GPU");
#endif
    gpu->copy_forces_to_host(particles, this_node);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    gpu->copy_dip_fld_to_host(particles, this_node);
#endif

#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("copy_forces_from_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE |
       PropagationMode::ROT_VS_INDEPENDENT)) {
    vs_relative_back_transfer_forces_and_torques(*cell_structure);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_CENTER_OF_MASS)) {
    vs_com_back_transfer_forces_and_torques(*cell_structure);
  }
#endif

  // Communication step: ghost forces
  cell_structure->ghosts_reduce_forces();

  // should be pretty late, since it needs to zero out the total force
  comfixed->apply(particles);

  // Needs to be the last one to be effective
  force_capping(*cell_structure, force_cap);

  // mark that forces are now up-to-date
  propagation->recalc_forces = false;
}
