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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_EXSCALICOS

#include "Actor.hpp"

#include "core/communication.hpp"
#include "core/electrostatics/exscalicos.hpp"

#include "script_interface/get_value.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Coulomb {

class CoulombExscalicos : public Actor<CoulombExscalicos, ::CoulombExscalicos> {
  std::shared_ptr<boost::mpi::environment> m_mpi_env_lock;

  static auto get_method_parameters(VariantMap const &params) {
    auto const valid_keys =
        std::set<std::string>{"accuracy", "r_cut", "alpha", "cao", "mesh"};
    for (auto const &kv : params) {
      if (not valid_keys.contains(kv.first)) {
        throw std::invalid_argument("Parameter '" + kv.first +
                                    "' is not a valid ExScaLiCoS parameter");
      }
    }
    if (not params.contains("accuracy")) {
      throw std::invalid_argument("Parameter 'accuracy' is required");
    }
    ExscalicosParameters parameters;
    parameters.accuracy = get_value<double>(params, "accuracy");
    parameters.r_cut = get_value_or<double>(params, "r_cut", -1.);
    parameters.alpha = get_value_or<double>(params, "alpha", -1.);
    parameters.cao = get_value_or<int>(params, "cao", -1);
    if (params.contains("mesh")) {
      parameters.mesh = get_value<Utils::Vector3i>(params, "mesh");
    }
    return parameters;
  }

public:
  CoulombExscalicos() {
    add_parameters({
        {"method_name", AutoParameter::read_only,
         [this]() { return actor()->get_method(); }},
        {"method_params", AutoParameter::read_only,
         [this]() {
           auto const parameters = actor()->get_parameters();
           // only report parameters that were set by the user or tuned
           auto result = VariantMap{{"accuracy", parameters.accuracy}};
           if (parameters.r_cut > 0.) {
             result["r_cut"] = parameters.r_cut;
           }
           if (parameters.alpha > 0.) {
             result["alpha"] = parameters.alpha;
           }
           if (parameters.cao > 0) {
             result["cao"] = parameters.cao;
           }
           if (parameters.mesh[0] > 0) {
             result["mesh"] = parameters.mesh;
           }
           return result;
         }},
    });
  }

  ~CoulombExscalicos() override {
    m_actor.reset();
    m_mpi_env_lock.reset();
  }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([&]() {
      auto const method_name = get_value<std::string>(params, "method_name");
      auto const method_params =
          get_method_parameters(get_value<VariantMap>(params, "method_params"));
      m_actor = make_coulomb_exscalicos(method_name, method_params);
      actor()->set_prefactor(get_value<double>(params, "prefactor"));
    });
    set_charge_neutrality_tolerance(params);
    // MPI communicator is needed to destroy the FFT plans
    m_mpi_env_lock = ::communication_environment->get_mpi_env();
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_near_field_delegation") {
      auto const delegate = get_value<bool>(params, "delegate");
      context()->parallel_try_catch(
          [&]() { actor()->set_near_field_delegation(delegate); });
      return {};
    }
    if (name == "get_near_field_delegation") {
      return actor()->get_near_field_delegation();
    }
    if (name == "get_available_methods") {
      return make_vector_of_variants(::CoulombExscalicos::available_methods());
    }
    return Actor<SIActorClass, CoreActorClass>::do_call_method(name, params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ESPRESSO_EXSCALICOS
