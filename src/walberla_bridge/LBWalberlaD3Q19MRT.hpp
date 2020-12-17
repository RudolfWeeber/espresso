#include "LBWalberlaImpl.hpp"
#include "generated_kernels/MRTLatticeModel.h"

namespace walberla {
class LBWalberlaD3Q19MRT : public LBWalberlaImpl<lbm::MRTLatticeModel> {
public:
  using LatticeModel = lbm::MRTLatticeModel;

  void construct_lattice_model(double viscosity, double magic_number) {
    const real_t omega = 2 / (6 * real_c(viscosity) + 1);
    const real_t omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id,
                     omega,   // bulk
                     omega,   // even
                     omega_2, // odd
                     omega)); // shear
  };
  void construct_lattice_model(LBRelaxationRates relaxation_rates) {
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id, relaxation_rates.omega_bulk,
                     relaxation_rates.omega_even, relaxation_rates.omega_odd,
                     relaxation_rates.omega_shear));
  };
  void set_viscosity(double viscosity, double magic_number) override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    const real_t omega = 2 / (6 * real_c(viscosity) + 1);
    const real_t omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    lm->omega_shear_ = omega;
    lm->omega_odd_ = omega_2;
    lm->omega_even_ = omega;
    lm->omega_bulk_ = omega;
    on_lattice_model_change();
  };
  double get_viscosity() const override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return (2 - lm->omega_shear_) / (6 * lm->omega_shear_);
  };
  double get_bulk_viscosity() const override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return lm->omega_bulk_;
  };
  double get_shear_viscosity() const override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return lm->omega_shear_;
  };
  double get_gamma_odd() const override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return lm->omega_odd_;
  };
  double get_gamma_even() const override {
    LatticeModel *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return lm->omega_even_;
  };
  LBWalberlaD3Q19MRT(double viscosity, double magic_number, double density,
                     double agrid, double tau,
                     const Utils::Vector3d &box_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers)
      : LBWalberlaImpl(agrid, tau, box_dimensions, node_grid, n_ghost_layers) {
    construct_lattice_model(viscosity, magic_number);
    setup_with_valid_lattice_model(density);
  };
  LBWalberlaD3Q19MRT(LBRelaxationRates relaxation_rates, double density,
                     double agrid, double tau,
                     const Utils::Vector3d &box_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers)
      : LBWalberlaImpl(agrid, tau, box_dimensions, node_grid, n_ghost_layers) {
    construct_lattice_model(relaxation_rates);
    setup_with_valid_lattice_model(density);
  };
};

} // namespace walberla
