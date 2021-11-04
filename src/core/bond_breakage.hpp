#pragma once

#include <unordered_map>

namespace BondBreakage {

enum class ActionType {
  NONE = 0,
  DELETE_BOND = 1,
  REVERT_BIND_AT_POINT_OF_COLLISION = 2
};

struct BreakageSpec {
  int bond_type;
  double breakage_length;
  ActionType action_type;
};
extern std::unordered_map<int, BreakageSpec> breakage_specs;
} // namespace BondBreakage
