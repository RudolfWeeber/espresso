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

#pragma once

#include <config/config.hpp>

#include "Particle.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"

#include <utils/index.hpp>
#include <utils/math/sqr.hpp>

struct GetNonbondedCutoff {
  GetNonbondedCutoff(System::System const &system) : m_system{system} {}
  auto operator()(int type_i, int type_j) const {
    return m_system.nonbonded_ias->get_ia_param(type_i, type_j).max_cut;
  }

private:
  System::System const &m_system;
};

/** Returns true if the particles are to be considered for short range
 *  interactions.
 */
template <typename CutoffGetter = GetNonbondedCutoff> class VerletCriterion {
  const double m_skin;
  const double m_eff_max_cut2;
  const double m_eff_coulomb_cut2 = 0.;
  const double m_eff_dipolar_cut2 = 0.;
  const double m_collision_cut2 = 0.;
  double eff_cutoff_sqr(double x) const {
    if (x == inactive_cutoff)
      return inactive_cutoff;
    return Utils::sqr(x + m_skin);
  }
  CutoffGetter get_nonbonded_cutoff;

private:
  bool check_pair(double dist2, double q1, double q2, double dipm1,
                  double dipm2, int type1, int type2) const {
    if (dist2 > m_eff_max_cut2)
      return false;
#ifdef ESPRESSO_ELECTROSTATICS
    if (dist2 <= m_eff_coulomb_cut2 and q1 != 0. and q2 != 0.)
      return true;
#endif
#ifdef ESPRESSO_DIPOLES
    if (dist2 <= m_eff_dipolar_cut2 and dipm1 != 0. and dipm2 != 0.)
      return true;
#endif
#ifdef ESPRESSO_COLLISION_DETECTION
    if (dist2 <= m_collision_cut2)
      return true;
#endif
    auto const ia_cut = get_nonbonded_cutoff(type1, type2);
    return (ia_cut != inactive_cutoff) &&
           (dist2 <= Utils::sqr(ia_cut + m_skin));
  }

public:
  VerletCriterion(System::System const &system, double skin, double max_cut,
                  double coulomb_cut = 0., double dipolar_cut = 0.,
                  double collision_detection_cutoff = 0.)
      : m_skin(skin), m_eff_max_cut2(eff_cutoff_sqr(max_cut)),
        m_eff_coulomb_cut2(eff_cutoff_sqr(coulomb_cut)),
        m_eff_dipolar_cut2(eff_cutoff_sqr(dipolar_cut)),
        m_collision_cut2(eff_cutoff_sqr(collision_detection_cutoff)),
        get_nonbonded_cutoff(system) {}

  template <typename Distance>
  bool operator()(const Particle &p1, const Particle &p2,
                  Distance const &dist) const {
    return check_pair(dist.dist2,
#ifdef ESPRESSO_ELECTROSTATICS
                      p1.q(), p2.q(),
#else
                      0., 0.,
#endif
#ifdef ESPRESSO_DIPOLES
                      p1.dipm(), p2.dipm(),
#else
                      0., 0.,
#endif
                      p1.type(), p2.type());
  }
};
