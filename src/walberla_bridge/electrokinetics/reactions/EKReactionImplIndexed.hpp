/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_WALBERLA_BRIDGE_ELECTROKINETICS_REACTIONS_EKREACTIONIMPLINDEXED_HPP
#define ESPRESSO_SRC_WALBERLA_BRIDGE_ELECTROKINETICS_REACTIONS_EKREACTIONIMPLINDEXED_HPP

#include "EKReactant.hpp"
#include "EKReactionBase.hpp"
#include "LatticeWalberla.hpp"

#include <domain_decomposition/BlockDataID.h>

#include <boost/optional.hpp>

#include <memory>
#include <vector>

namespace walberla {

class EKReactionImplIndexed : public EKReactionBase {
private:
  BlockDataID m_flagfield_id;
  BlockDataID m_indexvector_id;

  bool m_pending_changes;

public:
  EKReactionImplIndexed(std::shared_ptr<LatticeWalberla> lattice,
                        std::vector<std::shared_ptr<EKReactant>> reactants,
                        double coefficient);
  ~EKReactionImplIndexed() override = default;

  using EKReactionBase::get_coefficient;
  using EKReactionBase::get_lattice;
  using EKReactionBase::get_reactants;

  void perform_reaction() override;

  void set_node_is_boundary(const Utils::Vector3i &node, bool is_boundary);
  [[nodiscard]] boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node);

  [[nodiscard]] auto get_indexvector_id() const noexcept {
    return m_indexvector_id;
  }
  [[nodiscard]] auto get_flagfield_id() const noexcept {
    return m_flagfield_id;
  }

  void boundary_update();
};

} // namespace walberla

#endif
