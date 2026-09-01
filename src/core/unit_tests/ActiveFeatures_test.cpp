/*
 * Copyright (C) 2026 The ESPResSo project
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

#define BOOST_TEST_MODULE ActiveFeatures test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "cell_system/CellStructure.hpp"
#include "collision_detection/ActiveProtocol.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "exclusions.hpp"
#include "ghosts.hpp"
#include "integrators/Propagation.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_node.hpp"
#include "system/ActiveFeatures.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <memory>

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    auto system = System::System::create();
    system->set_box_l(Utils::Vector3d{10., 10., 10.});
    system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(system);
  }
  ~GlobalConfig() { ::System::reset_system(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

/** Apply @p mutation on the rank that owns the particle, then mark the
 *  cached state stale on all ranks (collective). */
static void mutate_particle(int p_id, auto &&mutation) {
  auto &system = System::get_system();
  auto *p = system.cell_structure->get_local_particle(p_id);
  if (p != nullptr and not p->is_ghost()) {
    mutation(*p);
  }
  system.propagation->recalc_active_features = true;
}

/** Remove all particles between test cases. */
struct ParticleCleanup {
  ~ParticleCleanup() { ::remove_all_particles(); }
};

BOOST_AUTO_TEST_SUITE(suite)

BOOST_FIXTURE_TEST_CASE(default_particle_sets_no_bits, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not system.propagation->recalc_active_features);
  BOOST_CHECK(not active_features.particles_are_virtual());
  BOOST_CHECK(not active_features.particles_can_rotate());
  BOOST_CHECK(not active_features.particles_are_fixed());
  BOOST_CHECK(not active_features.particles_have_charge());
#ifdef ESPRESSO_EXCLUSIONS
  BOOST_CHECK(not active_features.particles_have_exclusions());
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  BOOST_CHECK(not active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK(not active_features.particles_have_dipole_moment());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK(not active_features.particles_have_ext_force());
#endif
#ifdef ESPRESSO_ENGINE
  BOOST_CHECK(not active_features.particles_are_swimmers());
#endif
}

#ifdef ESPRESSO_EXCLUSIONS
BOOST_FIXTURE_TEST_CASE(exclusion_bit_is_reduced_across_ranks,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The particle lives on exactly one rank; the bit must be visible on all.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { add_exclusion(p, 42); });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_exclusions());
  // Removing the last exclusion clears the bit again.
  mutate_particle(0, [](Particle &p) { delete_exclusion(p, 42); });
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.particles_have_exclusions());
}
#endif

#ifdef ESPRESSO_ROTATION
BOOST_FIXTURE_TEST_CASE(update_if_needed_honors_dirty_flag, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  active_features.update();
  BOOST_CHECK(not active_features.particles_can_rotate());
  // Mutate without invalidating: the stale value must survive.
  auto *p = system.cell_structure->get_local_particle(0);
  if (p != nullptr and not p->is_ghost()) {
    p->set_can_rotate_all_axes();
  }
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.particles_can_rotate());
  // After invalidation the update must pick the change up.
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_can_rotate());
}
#endif

BOOST_FIXTURE_TEST_CASE(property_bits_follow_particle_state, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  mutate_particle(0,
                  [](Particle &p) { p.stoner_wohlfarth_is_enabled() = true; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  mutate_particle(0, [](Particle &p) { p.dipm() = 1.5; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_dipole_moment());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  mutate_particle(0, [](Particle &p) {
    p.ext_force() = Utils::Vector3d{1., 0., 0.};
    p.set_fixed_along(0, true);
  });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_ext_force());
  BOOST_CHECK(active_features.particles_are_fixed());
#endif
#ifdef ESPRESSO_ENGINE
  mutate_particle(0, [](Particle &p) { p.swimming().swimming = true; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_are_swimmers());
#endif
  // Removing the carrier particle clears every bit.
  ::remove_particle(0);
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  BOOST_CHECK(not active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK(not active_features.particles_have_dipole_moment());
#endif
  BOOST_CHECK(not active_features.particles_are_fixed());
}

BOOST_FIXTURE_TEST_CASE(used_propagations_from_same_sweep, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  auto &propagation = *system.propagation;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) {
    p.propagation() = PropagationMode::TRANS_LANGEVIN;
  });
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations & PropagationMode::TRANS_LANGEVIN);
  BOOST_CHECK(
      (propagation.used_propagations & PropagationMode::SYSTEM_DEFAULT) == 0);
}

BOOST_FIXTURE_TEST_CASE(system_default_expansion_and_thermostat_invalidation,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  auto &propagation = *system.propagation;
  // A default particle carries SYSTEM_DEFAULT, which expands to
  // default_propagation, which depends on the thermostat switch.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations & PropagationMode::TRANS_NEWTON);
  BOOST_CHECK(not propagation.recalc_active_features);
  // Changing the thermostat must invalidate and change the expansion.
  system.thermostat->thermo_switch = THERMO_LANGEVIN;
  system.on_thermostat_param_change();
  BOOST_CHECK(propagation.recalc_active_features);
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations & PropagationMode::TRANS_LANGEVIN);
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(propagation.used_propagations & PropagationMode::ROT_LANGEVIN);
#endif
  // Restore for other test cases.
  system.thermostat->thermo_switch = THERMO_OFF;
  system.on_thermostat_param_change();
}

BOOST_FIXTURE_TEST_CASE(set_integ_switch_invalidates, ParticleCleanup) {
  auto &system = System::get_system();
  system.active_features->update();
  BOOST_CHECK(not system.propagation->recalc_active_features);
  system.propagation->set_integ_switch(INTEG_METHOD_NVT);
  BOOST_CHECK(system.propagation->recalc_active_features);
}

BOOST_FIXTURE_TEST_CASE(particle_creation_invalidates, ParticleCleanup) {
  auto &system = System::get_system();
  system.active_features->update();
  BOOST_CHECK(not system.propagation->recalc_active_features);
  // make_new_particle fires on_particle_change, which must set the flag.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  BOOST_CHECK(system.propagation->recalc_active_features);
}

#ifdef ESPRESSO_ROTATION
BOOST_FIXTURE_TEST_CASE(orientation_ghosts_off_for_default_particles,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // A default particle cannot rotate, so the ROT_* bits in the SYSTEM_DEFAULT
  // expansion must not trigger orientation ghost exchange.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
  BOOST_CHECK((system.get_global_ghost_flags() & Cells::DATA_PART_QUAT) == 0u);
  BOOST_CHECK_EQUAL(system.get_force_reduce_ghost_flags(), GHOSTTRANS_FORCE);
}

BOOST_FIXTURE_TEST_CASE(orientation_ghosts_on_for_rotating_particles,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { p.set_can_rotate_all_axes(); });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.orientation_ghosts_needed());
  BOOST_CHECK((system.get_global_ghost_flags() & Cells::DATA_PART_QUAT) != 0u);
  BOOST_CHECK((system.get_force_reduce_ghost_flags() & GHOSTTRANS_TORQUE) !=
              0u);
}

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_FIXTURE_TEST_CASE(orientation_ghosts_on_for_stoner_wohlfarth,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The Stoner-Wohlfarth arm must fire from the particle bit alone,
  // without any rotating particle and without the Langevin thermostat.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0,
                  [](Particle &p) { p.stoner_wohlfarth_is_enabled() = true; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.orientation_ghosts_needed());
}
#endif

#ifdef ESPRESSO_DIPOLES
BOOST_FIXTURE_TEST_CASE(dipole_moment_alone_needs_no_orientation_ghosts,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The dipole arm requires a dipolar solver AND a nonzero moment;
  // a moment without a solver must not fire it.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { p.dipm() = 1.; });
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
}
#endif

#ifdef ESPRESSO_GAY_BERNE
BOOST_FIXTURE_TEST_CASE(gay_berne_aggregate, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.has_gay_berne());
  system.nonbonded_ias->make_particle_type_exist(0);
  auto &ia_params = system.nonbonded_ias->get_ia_param(0, 0);
  ia_params.gay_berne = GayBerne_Parameters(1., 1., 2., 1., 1., 1., 1.);
  system.on_non_bonded_ia_change();
  BOOST_CHECK(active_features.has_gay_berne());
  BOOST_CHECK(active_features.orientation_ghosts_needed());
  // Reset so later test cases see a clean interaction table.
  ia_params.gay_berne = GayBerne_Parameters();
  system.on_non_bonded_ia_change();
}
#endif

#if defined(ESPRESSO_COLLISION_DETECTION) &&                                   \
    defined(ESPRESSO_VIRTUAL_SITES_RELATIVE)
BOOST_FIXTURE_TEST_CASE(collision_detection_active_needs_orientation_ghosts,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
  // BindCenters::initialize requires the bond type to exist and be a pair
  // bond.
  auto const bond_centers_id = 0;
  auto const bond = HarmonicBond(200.0, 0.3, 1.0);
  system.bonded_ias->insert(bond_centers_id,
                            std::make_shared<Bonded_IA_Parameters>(bond));
  system.collision_detection->set_protocol(
      std::make_shared<CollisionDetection::ActiveProtocol>(
          CollisionDetection::BindCenters(1., bond_centers_id)));
  BOOST_CHECK(active_features.orientation_ghosts_needed());
  // Reset so later test cases see collision detection off.
  system.collision_detection->unset_protocol();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
}
#endif
#endif // ESPRESSO_ROTATION

BOOST_AUTO_TEST_SUITE_END()
