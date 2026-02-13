/*
 * Copyright (C) 2019-2024 The ESPResSo project
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

/**
 * @file
 * Out-of-class VTK writer registration definition for
 * @ref walberla::LBWalberlaImpl.
 */

// ---- VTK registration (out-of-class definition) ----

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::register_vtk_field_writers(
    walberla::vtk::VTKOutput &vtk_obj, LatticeModel::units_map const &units,
    int flag_observables) {
#if defined(__CUDACC__)
  auto const allocate_cpu_field_if_empty =
      [&]<typename Field>(auto const &blocks, std::string name,
                          std::optional<BlockDataID> &cpu_field) {
        if (not cpu_field) {
          cpu_field = field::addToStorage<Field>(
              blocks, name, FloatType{0}, field::fzyx,
              m_lattice->get_ghost_layers(), m_host_field_allocator);
        }
      };
#endif
  if (flag_observables & static_cast<int>(OutputVTK::density)) {
    auto const unit_conversion =
        FloatType_c(zero_centered_to_md(units.at("density")));
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      auto const &blocks = m_lattice->get_blocks();
      allocate_cpu_field_if_empty.template operator()<PdfFieldCpu>(
          blocks, "pdfs_cpu", m_pdf_cpu_field_id);
      vtk_obj.addBeforeFunction(gpu::fieldCpyFunctor<PdfFieldCpu, PdfField>(
          blocks, *m_pdf_cpu_field_id, m_pdf_field_id));
    }
#endif
    vtk_obj.addCellDataWriter(std::make_shared<DensityVTKWriter<float>>(
        m_pdf_field_id, "density", unit_conversion));
  }
  if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
    auto const unit_conversion = FloatType_c(units.at("velocity"));
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      auto const &blocks = m_lattice->get_blocks();
      allocate_cpu_field_if_empty.template operator()<VectorFieldCpu>(
          blocks, "vel_cpu", m_vel_cpu_field_id);
      vtk_obj.addBeforeFunction(
          gpu::fieldCpyFunctor<VectorFieldCpu, VectorField>(
              blocks, *m_vel_cpu_field_id, m_velocity_field_id));
    }
#endif
    vtk_obj.addCellDataWriter(std::make_shared<VelocityVTKWriter<float>>(
        m_velocity_field_id, "velocity_vector", unit_conversion));
  }
  if (flag_observables & static_cast<int>(OutputVTK::pressure_tensor)) {
    auto const unit_conversion =
        FloatType_c(zero_centered_to_md(units.at("pressure")));
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      auto const &blocks = m_lattice->get_blocks();
      allocate_cpu_field_if_empty.template operator()<PdfFieldCpu>(
          blocks, "pdfs_cpu", m_pdf_cpu_field_id);
      vtk_obj.addBeforeFunction(gpu::fieldCpyFunctor<PdfFieldCpu, PdfField>(
          blocks, *m_pdf_cpu_field_id, m_pdf_field_id));
    }
#endif
    vtk_obj.addCellDataWriter(std::make_shared<PressureTensorVTKWriter<float>>(
        m_pdf_field_id, "pressure_tensor", unit_conversion,
        pressure_tensor_correction_factor()));
  }
}
