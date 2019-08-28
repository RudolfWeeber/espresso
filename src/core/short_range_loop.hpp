/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_SHORT_RANGE_HPP
#define CORE_SHORT_RANGE_HPP

#include "algorithm/for_each_pair.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "lees_edwards.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/iterator/indirect_iterator.hpp>
#include <profiler/profiler.hpp>

#include <utility>

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  const double dist2;
};

namespace detail {
class MinimalImageDistance {
  const BoxGeometry& m_box;
  public:
    MinimalImageDistance(const BoxGeometry& box) : m_box{box} {};

  inline
  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(get_mi_vector(p1.r.p, p2.r.p, m_box));
  }
};

/**
 * @brief Functor that returns true for
 *        any arguments.
 */
struct True {
  template <class... T> bool operator()(T...) const { return true; }
};
} // namespace detail

template <class ParticleKernel, class PairKernel,
          class VerletCriterion = detail::True>
void short_range_loop(ParticleKernel &&particle_kernel,
                      PairKernel &&pair_kernel,
                      const VerletCriterion &verlet_criterion = {}) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(get_resort_particles() == Cells::RESORT_NONE);

  if (cell_structure.min_range != INACTIVE_CUTOFF) {
    auto first = boost::make_indirect_iterator(local_cells.begin());
    auto last = boost::make_indirect_iterator(local_cells.end());
    Algorithm::for_each_pair(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel),
        detail::MinimalImageDistance{box_geo}, verlet_criterion,
        cell_structure.use_verlet_list, rebuild_verletlist);

    rebuild_verletlist = 0;
#ifdef LEES_EDWARDS
    LeesEdwards::on_verlet_rebuild();
#endif    
  } else {
    for (auto &p : cell_structure.local_cells().particles()) {
      particle_kernel(p);
    }
  }
}

#endif
