/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "Bond breakage"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"
#include "ParticleFactory.hpp"

#include "bond_breakage/actions.hpp"
#include "bond_breakage/bond_breakage.hpp"

#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <optional>

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<System::System> system;
} // namespace espresso

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    ErrorHandling::init_error_handling(comm_cart);
    espresso::system = System::System::create();
    espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(espresso::system);
  }
  ~GlobalConfig() {
    espresso::system.reset();
    ::System::reset_system();
    ErrorHandling::deinit_error_handling();
  }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_AUTO_TEST_CASE(test_actions_equality) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3} == Action{1, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{0, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 0, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 2, 0}));
  }

  {
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3} == Action{1, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{0, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {0, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 0}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 4}, 0}));
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2} == Action{1, 2}));
    BOOST_CHECK((Action{1, 2} != Action{0, 2}));
    BOOST_CHECK((Action{1, 2} != Action{1, 0}));
  }
}

BOOST_AUTO_TEST_CASE(test_actions_hash_value) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3}.hash_value() == Action{1, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{0, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 0, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 2, 0}.hash_value()));
  }

  {
    // clang-format off
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() == Action{1, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{0, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {0, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 0}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 4}, 0}.hash_value()));
    // clang-format on
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2}.hash_value() == Action{1, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{0, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{1, 0}.hash_value()));
  }
}

/**
 * Regression test for bug-sweep finding #33: an @ref
 * BondBreakage::ActionType::NONE breakage spec applied to a pair bond must be a
 * clean no-op. On the unfixed code the NONE pair-bond queue entry falls into
 * the REVERT/angle-bond @c else branch of @c actions_for_breakage, which
 * dereferences @c bond_partners[1] -- an empty @c std::optional<int> for a pair
 * bond (the second partner is @c std::nullopt). This is undefined behaviour;
 * with libstdc++ assertions (or a sanitizer) it aborts. After the fix the
 * NONE entry returns no actions and the particle's bonds are left untouched.
 */
BOOST_FIXTURE_TEST_CASE(none_action_pair_bond_is_noop, ParticleFactory) {
  auto &system = *espresso::system;
  system.set_box_l(Utils::Vector3d::broadcast(10.));
  system.set_time_step(0.01);
  system.cell_structure->set_verlet_skin(0.4);

  auto const comm = boost::mpi::communicator();
  // single-rank scenario: place both partners on the only rank
  if (comm.size() != 1) {
    return;
  }

  auto const pid0 = 0;
  auto const pid1 = 1;
  auto const bond_type = 0;
  create_particle({1.0, 1.0, 1.0}, pid0, 0);
  create_particle({1.5, 1.0, 1.0}, pid1, 0);

  // Configure a NONE breakage spec with the default breakage_length of 0.0,
  // exactly as a default-constructed core BreakageSpec / Python
  // BreakageSpec(action_type="none") would produce.
  BondBreakage::BondBreakage bond_breakage{};
  bond_breakage.breakage_specs[bond_type] =
      std::make_shared<BondBreakage::BreakageSpec>(
          BondBreakage::BreakageSpec{0.0, BondBreakage::ActionType::NONE});

  // Queue a *pair*-bond breakage event: the second partner is std::nullopt.
  // This mirrors bond_handler() for a pair bond (partners.size() == 1).
  auto const queued = bond_breakage.check_and_handle_breakage(
      pid0, {{pid1, std::nullopt}}, bond_type, /* distance */ 0.0);
  BOOST_REQUIRE(queued);

  // On the UNFIXED code this dereferences the empty optional bond_partners[1]
  // (undefined behaviour -> abort under _GLIBCXX_ASSERTIONS or a sanitizer).
  // On the FIXED code it is a clean no-op.
  BOOST_CHECK_NO_THROW(bond_breakage.process_queue(system));

  // Even in a non-hardened build the bug is observable: dereferencing the
  // empty optional yields the value-initialized payload 0, which resolves to
  // the (non-virtual) particle pid0 and makes the REVERT/angle branch emit a
  // spurious "has to be configured for the bond on the virtual site" runtime
  // error. A correct NONE no-op must not emit any runtime error.
  BOOST_CHECK_EQUAL(check_runtime_errors_local(), 0);

  // The bonds of the involved particles must be left untouched by a NONE spec.
  BOOST_CHECK(system.cell_structure->get_local_particle(pid0)->bonds().empty());
  BOOST_CHECK(system.cell_structure->get_local_particle(pid1)->bonds().empty());
}
