/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLEBODYANGULARVELOCITIES_HPP
#define OBSERVABLES_PARTICLEBODYANGULARVELOCITIES_HPP

#include "PidObservable.hpp"

#include <utils/Span.hpp>

#include <cstddef>
#include <vector>

namespace Observables {

class ParticleBodyAngularVelocities : public PidObservable {
public:
  using PidObservable::PidObservable;

  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    std::vector<double> res(n_values());
#ifdef ROTATION
    for (std::size_t i = 0; i < particles.size(); i++) {
      res[3 * i + 0] = particles[i].get().m.omega[0];
      res[3 * i + 1] = particles[i].get().m.omega[1];
      res[3 * i + 2] = particles[i].get().m.omega[2];
    }
#endif
    return res;
  }

  std::vector<std::size_t> shape() const override { return {ids().size(), 3}; }
};

} // Namespace Observables
#endif
