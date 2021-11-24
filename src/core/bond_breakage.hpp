#pragma once

#include <memory>
#include <unordered_map>

namespace BondBreakage {

enum class ActionType {
  NONE = 0,
  DELETE_BOND = 1,
  REVERT_BIND_AT_POINT_OF_COLLISION = 2
};

struct BreakageSpec {
  double breakage_length;
  ActionType action_type;
};
extern std::unordered_map<int, std::shared_ptr<BreakageSpec>> breakage_specs;

/** @brief Checks if the bond between the particles shoudl break, if yes, queue
 * it */
bool check_and_handle_breakage(int particle_id, int bond_partner_id,
                               int bond_type, double distance);

void clear_queue();

void process_queue();

} // namespace BondBreakage
