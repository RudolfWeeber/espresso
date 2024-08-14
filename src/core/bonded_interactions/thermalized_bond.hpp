/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/** \file
 *  Routines to thermalize the center of mass and distance of a particle pair.
 */

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <memory>
#include <optional>
#include <tuple>

namespace Thermostat {
class Thermostat;
}

/** Parameters for Thermalized bond */
struct ThermalizedBond {
  double temp_com;
  double gamma_com;
  double temp_distance;
  double gamma_distance;
  double r_cut;
  double pref1_com;
  double pref2_com;
  double pref1_dist;
  double pref2_dist;

  double cutoff() const { return r_cut; }

  static constexpr int num = 1;

  ThermalizedBond(double temp_com, double gamma_com, double temp_distance,
                  double gamma_distance, double r_cut) {
    this->temp_com = temp_com;
    this->gamma_com = gamma_com;
    this->temp_distance = temp_distance;
    this->gamma_distance = gamma_distance;
    this->r_cut = r_cut;

    pref1_com = -1.;
    pref2_com = -1.;
    pref1_dist = -1.;
    pref2_dist = -1.;
  }

  void recalc_prefactors(double time_step) {
    pref1_com = gamma_com;
    pref2_com = std::sqrt(24. * gamma_com / time_step * temp_com);
    pref1_dist = gamma_distance;
    pref2_dist = std::sqrt(24. * gamma_distance / time_step * temp_distance);
  }

  void set_thermostat_view(
      std::weak_ptr<Thermostat::Thermostat const> const &thermostat) {
    m_thermostat = thermostat;
  }

  void unset_thermostat_view() { m_thermostat.reset(); }

  std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d>>
  forces(Particle const &p1, Particle const &p2,
         Utils::Vector3d const &dx) const;

private:
  std::weak_ptr<Thermostat::Thermostat const> m_thermostat;
};
