
#include "correlators.hpp"
#include "integrate.hpp" 

namespace Correlators {
std::vector<std::shared_ptr<Correlators::Correlator>> auto_update_correlators;

void auto_update() {
for (auto& c : auto_update_correlators) {
  if (sim_time - c->last_update >c->dt *0.9999) {
    if (c->get_data()) {
      throw std::runtime_error("Updating of a correlator failed.");
    }
  }
}

}

}

