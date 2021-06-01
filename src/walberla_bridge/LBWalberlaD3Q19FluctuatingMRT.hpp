#include "LBWalberlaImpl.hpp"
#include "relaxation_rates.hpp"

#ifdef __AVX__
#include "generated_kernels/FluctuatingMRTLatticeModelAvx.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModelAvx
#else
#include "generated_kernels/FluctuatingMRTLatticeModel.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModel
#endif

namespace walberla {
class LBWalberlaD3Q19FluctuatingMRT : public LBWalberlaImpl<LatticeModelName> {

  using LatticeModel = LatticeModelName;

public:
  void construct_lattice_model(double viscosity, double kT, unsigned int seed) {
    const real_t omega = shear_mode_relaxation_rate(viscosity);
    const real_t omega_odd = odd_mode_relaxation_rate(omega);
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id, real_c(kT),
                     omega,     // bulk
                     omega,     // even
                     omega_odd, // odd
                     omega,     // shear
                     1, seed));
  };
  void set_viscosity(double viscosity) override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    const real_t omega = shear_mode_relaxation_rate(viscosity);
    const real_t omega_odd = odd_mode_relaxation_rate(omega);
    lm->omega_shear_ = omega;
    lm->omega_odd_ = omega_odd;
    lm->omega_even_ = omega;
    lm->omega_bulk_ = omega;
    on_lattice_model_change();
  };
  double get_viscosity() const override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return viscosity_from_shear_relaxation_rate(lm->omega_shear_);
  };
  LBWalberlaD3Q19FluctuatingMRT(double viscosity, double density,
                                const Utils::Vector3i &grid_dimensions,
                                const Utils::Vector3i &node_grid,
                                int n_ghost_layers, double kT,
                                unsigned int seed)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers) {
    m_kT = kT;
    construct_lattice_model(viscosity, kT, seed);
    setup_with_valid_lattice_model(density);
  };
  void integrate() override {
    m_time_loop->singleStep();
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    lm->time_step_ += 1;
    on_lattice_model_change();
  };
  double get_kT() const override { return m_kT; };

  uint64_t get_rng_state() const override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return lm->time_step_;
  }

  void set_rng_state(uint64_t counter) override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    lm->time_step_ = counter;
    on_lattice_model_change();
  }

private:
  double m_kT;
};

} // namespace walberla
#undef LatticeModelName

