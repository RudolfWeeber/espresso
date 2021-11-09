#include "bond_breakage.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "particle_data.hpp"
#include "utils/mpi/gather_buffer.hpp"

#include <boost/variant.hpp>
#include <unordered_set>
#include <vector>

namespace BondBreakage {

// Bond breakage specifications
std::unordered_map<int, BreakageSpec> breakage_specs;

// Delete Actions
struct DeleteBond {
  int particle_id;
  int bond_partner_id;
  int bond_type;
  std::size_t hash_value() const {
    std::size_t seed = 3875;
    boost::hash_combine(seed, particle_id);
    boost::hash_combine(seed, bond_partner_id);
    boost::hash_combine(seed, bond_type);
    return seed;
  };
  bool operator==(const DeleteBond &rhs) const {
    return rhs.particle_id == particle_id and
           rhs.bond_partner_id == bond_partner_id and
           rhs.bond_type == bond_type;
  };
};

struct DeleteAllBonds {
  int particle_id_1;
  int particle_id_2;
  std::size_t hash_value() const {
    std::size_t seed = 75;
    boost::hash_combine(seed, particle_id_1);
    boost::hash_combine(seed, particle_id_2);
    return seed;
  };
  bool operator==(const DeleteAllBonds &rhs) const {
    return rhs.particle_id_1 == particle_id_1 and
           rhs.particle_id_2 == particle_id_2;
  };
};
} // namespace BondBreakage

// Hash support for std::unordered_set
namespace boost {
template <> struct hash<BondBreakage::DeleteBond> {
  std::size_t operator()(const BondBreakage::DeleteBond &t) const noexcept {
    return t.hash_value();
  };
};
template <> struct hash<BondBreakage::DeleteAllBonds> {
  std::size_t operator()(const BondBreakage::DeleteAllBonds &t) const noexcept {
    return t.hash_value();
  };
};
} // namespace boost

namespace BondBreakage {
// Variant holding any of the actions
using Action = boost::variant<DeleteBond, DeleteAllBonds>;

// Set of actions
using ActionSet = std::unordered_set<Action>;

// Queue to record bond broken during a time step
struct QueueEntry {
  int particle_id;
  int bond_partner_id;
  int bond_type;

  /// Serialization for synchronization across mpi ranks
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &particle_id;
    ar &bond_partner_id;
    ar &bond_type;
  }
};
using Queue = std::vector<QueueEntry>;
Queue queue;

/** @brief Retrieve breakage specification for the bond type */
boost::optional<BreakageSpec> get_breakage_spec(int bond_type) {
  if (breakage_specs.find(bond_type) != breakage_specs.end()) {
    return {breakage_specs.at(bond_type)};
  }
  return {};
}

/** Add a particle+bond combination to the breakage queue */
void queue_breakage(int particle_id, int bond_partner_id, int bond_type) {
  queue.emplace_back(QueueEntry{particle_id, bond_partner_id, bond_type});
}

/** @brief Checks if the bond between the particles shoudl break, if yes, queue
 * it */
bool check_and_handle_breakage(int particle_id, int bond_partner_id,
                               int bond_type, double distance) {
  // Retrieve specification for this bond type
  auto spec = get_breakage_spec(bond_type);
  if (!spec)
    return false; // No breakage rule for this bond type

  // Is the bond length longer than the breakage length?
  if (distance >= (*spec).breakage_length) {
    queue_breakage(particle_id, bond_partner_id, bond_type);
    return true;
  }
  return false;
}

void clear_queue() { queue.clear(); }

/** @brief Gathers combined queue from all mpi ranks */
Queue gather_global_queue(const Queue &local_queue) {
  Queue res = local_queue;
  Utils::Mpi::gather_buffer(res, comm_cart);
  boost::mpi::broadcast(comm_cart, res, 0);
  return res;
}

/** @brief Constructs the actions to take for a breakage queue entry */
ActionSet actions_for_breakage(const QueueEntry &e) {
  // Retrieve relevant breakage spec
  auto spec = get_breakage_spec(e.bond_type);
  if (!spec)
    throw std::runtime_error("Bond breakage spec not available.");

  // Handle different action types
  if ((*spec).action_type == ActionType::DELETE_BOND)
    return {DeleteBond{e.particle_id, e.bond_partner_id, e.bond_type}};
  if ((*spec).action_type == ActionType::REVERT_BIND_AT_POINT_OF_COLLISION) {
    // We need to find the base particles for the two virtual sites
    // between which the bond broke.
    auto p1 = cell_structure.get_local_particle(e.particle_id);
    auto p2 = cell_structure.get_local_particle(e.bond_partner_id);
    if (!p1 || !p2)
      return {}; // particles not on this mpi rank

    if (!p1->p.is_virtual || !p2->p.is_virtual) {
      throw std::runtime_error(
          "The REVERT_BIND_AT_POINT_OF_COLLISION bond breakage action has to "
          "be configured for the bond on the virtual site. Encountered a "
          "particle hat is not virtual.");
    }

    return {
        // Bond between virtual sites
        DeleteBond{e.particle_id, e.bond_partner_id, e.bond_type},
        // Bond between base particles. We do not know, on which of the two
        // the bond is defined, since bonds are stored only on one partner
        DeleteAllBonds{p1->p.vs_relative.to_particle_id,
                       p2->p.vs_relative.to_particle_id},
        DeleteAllBonds{p2->p.vs_relative.to_particle_id,
                       p1->p.vs_relative.to_particle_id},
    };
  }
}

// andler for the different delete events
class execute : public boost::static_visitor<> {
public:
  void operator()(const DeleteBond &d) const {
    auto p = cell_structure.get_local_particle(d.particle_id);
    if (!p)
      return;
    local_remove_bond(*p, {d.bond_type, d.bond_partner_id});
  }
  void operator()(const DeleteAllBonds d) const {
    auto p = cell_structure.get_local_particle(d.particle_id_1);
    if (!p)
      return;
    local_remove_pair_bonds_to(*p, d.particle_id_2);
  };
};

void process_queue() {
  if (breakage_specs.empty())
    return;

  auto global_queue = gather_global_queue(queue);

  // Construc tdelete actions from breakage queue
  ActionSet actions = {};
  for (auto const &e : global_queue) {
    // Convert to merge() once we are on c++ 17
    auto to_add = actions_for_breakage(e);
    actions.insert(to_add.begin(), to_add.end());
  }

  // Execute actions
  for (auto const &a : actions) {
    boost::apply_visitor(execute{}, a);
  };
}
} // namespace BondBreakage
