/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#include "LBWalberlaImpl.hpp"
#include "LBWalberlaImplColorGradient.hpp"
#include "LBWalberlaImplSingleComponent.hpp"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <memory>

std::shared_ptr<LBWalberlaBase>
new_lb_walberla_cpu(std::shared_ptr<LatticeWalberla> const &lattice,
                    std::vector<double> viscosity, double density,
                    bool single_precision) {
  if (viscosity.size() == 1u) {
    if (single_precision) {
      return std::make_shared<
          walberla::LBWalberlaImplSingleComponent<float, lbmpy::Arch::CPU>>(
          lattice, viscosity, density, false);
    }
    return std::make_shared<
        walberla::LBWalberlaImplSingleComponent<double, lbmpy::Arch::CPU>>(
        lattice, viscosity, density, false);
  }
  if (single_precision) {
    return std::make_shared<
        walberla::LBWalberlaImplColorGradient<float, lbmpy::Arch::CPU>>(
        lattice, viscosity, density, true);
  }
  return std::make_shared<
      walberla::LBWalberlaImplColorGradient<double, lbmpy::Arch::CPU>>(
      lattice, viscosity, density, true);
}

std::shared_ptr<LBWalberlaBase>
new_lb_walberla_cpu(std::shared_ptr<LatticeWalberla> const &lattice,
                    double viscosity, double density, bool single_precision) {
  return new_lb_walberla_cpu(lattice, std::vector<double>{viscosity}, density,
                             single_precision);
}
