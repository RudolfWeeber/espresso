/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef SCRIPT_INTERFACE_REGISTRY_HPP
#define SCRIPT_INTERFACE_REGISTRY_HPP

#include "ScriptInterface.hpp"
#include "Serializer.hpp" 

namespace ScriptInterface {

template <typename ManagedType>
class ScriptObjectRegistry : public ScriptInterfaceBase {
public:
  virtual void add_in_core(std::shared_ptr<ManagedType> obj_ptr) = 0;
  virtual void remove_in_core(std::shared_ptr<ManagedType> obj_ptr) = 0;
  Variant call_method(std::string const &method,
                      VariantMap const &parameters) override {

    if (method == "add") {
      auto obj_ptr =
          get_value<std::shared_ptr<ManagedType>>(parameters.at("object"));

      add_in_core(obj_ptr);
      m_elements.push_back(obj_ptr);
    }

    if (method == "remove") {
      auto obj_ptr =
          get_value<std::shared_ptr<ManagedType>>(parameters.at("object"));

      remove_in_core(obj_ptr);
      m_elements.erase(
          std::remove(m_elements.begin(), m_elements.end(), obj_ptr),
          m_elements.end());
    }

    if (method == "get_elements") {
      std::vector<Variant> ret;
      ret.reserve(m_elements.size());

      for (auto const &e : m_elements)
        ret.emplace_back(e->id());

      return ret;
    }

    if (method == "clear") {
      for (auto const &e : m_elements) {
        remove_in_core(e);
      }

      m_elements.clear();
    }

    if (method == "size") {
      return static_cast<int>(m_elements.size());
    }

    if (method == "empty") {
      return m_elements.empty();
    }

    return none;
  }
  /**
 * @brief  Return a Variant representation of the state of the object.
 *
 * This should return the internal state of the instance, so that
 * the instance can be restored from this information.  The default
 * implementation stores all the public parameters, including object
 * parameters that are captured by calling get_state on them.
 */
  Variant get_state() const override {
  std::vector<Variant> state;

  state.reserve(m_elements.size());

  for (auto const &e : m_elements) {
    auto v = Variant{e->id()};
    state.emplace_back(boost::apply_visitor(Serializer{}, v));
  }
  return state;
}

void set_state(Variant const &state) override {
  using boost::get;
  using std::vector;

  UnSerializer u;
  
  m_elements.clear();

  for (auto &e : get<vector<Variant>>(state)) {
    auto oid = get_value<std::shared_ptr<ManagedType>>(
       boost::apply_visitor(u, e));
    m_elements.emplace_back(oid);
  }
}

private:
  std::vector<std::shared_ptr<ManagedType>> m_elements;
};
} // Namespace ScriptInterface
#endif
