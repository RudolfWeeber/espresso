/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

/**
 * @file
 * @ref walberla::LBWalberlaImpl implements the interface of the LB
 * waLBerla bridge using sweeps generated by lbmpy
 * (see <tt>maintainer/walberla_kernels</tt>).
 */

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/IBlock.h>
#include <field/GhostLayerField.h>
#include <field/vtk/FlagFieldCellFilter.h>
#include <field/vtk/VTKWriter.h>

#include <domain_decomposition/SharedSweep.h>

#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/communication/PackInfo.h>
#include <lbm/communication/PdfFieldPackInfo.h>
#include <lbm/field/AddToStorage.h>
#include <lbm/field/PdfField.h>
#include <lbm/sweeps/CellwiseSweep.h>

#include <stencil/D3Q19.h>
#include <stencil/D3Q27.h>

#include "../BoundaryHandling.hpp"
#include "ResetForce.hpp"
#include "lb_kernels.hpp"
#include "vtk_writers.hpp"

#include "walberla_bridge/BlockAndCell.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"
#include "walberla_bridge/lattice_boltzmann/InterpolateAndShiftAtBoundary.hpp"
#include "walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp"
#include "walberla_bridge/lattice_boltzmann/LeesEdwardsPack.hpp"
#include "walberla_bridge/utils/boundary_utils.hpp"
#include "walberla_bridge/utils/walberla_utils.hpp"

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>
#include <utils/math/make_lin_space.hpp>

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {

/** @brief Class that runs and controls the LB on waLBerla. */
template <typename FloatType = double>
class LBWalberlaImpl : public LBWalberlaBase {
  template <typename T> inline FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }
  using CollisionModelLeesEdwards =
      typename detail::KernelTrait<FloatType>::CollisionModelLeesEdwards;
  using CollisionModelThermalized =
      typename detail::KernelTrait<FloatType>::CollisionModelThermalized;
  using StreamSweep = typename detail::KernelTrait<FloatType>::StreamSweep;
  using InitialPDFsSetter =
      typename detail::KernelTrait<FloatType>::InitialPDFsSetter;
  using BoundaryModel = BoundaryHandling<
      Vector3<FloatType>,
      typename detail::BoundaryHandlingTrait<FloatType>::Dynamic_UBB>;

protected:
  using CollisionModel =
      boost::variant<CollisionModelThermalized, CollisionModelLeesEdwards>;

private:
  class : public boost::static_visitor<> {
  public:
    void operator()(CollisionModelThermalized &cm, IBlock *b) { cm(b); }

    void operator()(CollisionModelLeesEdwards &cm, IBlock *b) {
      cm.v_s_ = static_cast<decltype(cm.v_s_)>(
          m_lees_edwards_callbacks->get_shear_velocity());
      cm(b);
    }
    void register_lees_edwards_callbacks(
        std::shared_ptr<LeesEdwardsPack> lees_edwards_callbacks) {
      m_lees_edwards_callbacks = std::move(lees_edwards_callbacks);
    }

  private:
    std::shared_ptr<LeesEdwardsPack> m_lees_edwards_callbacks;

  } run_collide_sweep;

  FloatType shear_mode_relaxation_rate() const {
    return FloatType{2} / (FloatType{6} * m_viscosity + FloatType{1});
  }

  FloatType odd_mode_relaxation_rate(
      FloatType shear_relaxation,
      FloatType magic_number = FloatType{3} / FloatType{16}) const {
    return (FloatType{4} - FloatType{2} * shear_relaxation) /
           (FloatType{4} * magic_number * shear_relaxation + FloatType{2} -
            shear_relaxation);
  }

  void reset_boundary_handling() {
    auto const &blocks = get_lattice().get_blocks();
    m_boundary = std::make_shared<BoundaryModel>(blocks, m_pdf_field_id,
                                                 m_flag_field_id);
  }

public:
  // Type definitions
  typedef stencil::D3Q19 Stencil;
  using PdfField = GhostLayerField<FloatType, Stencil::Size>;
  using VectorField = GhostLayerField<FloatType, 3u>;
  using FlagField = typename BoundaryModel::FlagField;
  using Lattice_T = LatticeWalberla::Lattice_T;

private:
  FloatType getDensity(const BlockAndCell &bc) const {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    return lbm::accessor::Density::get(*pdf_field, bc.cell.x(), bc.cell.y(),
                                       bc.cell.z());
  }

  FloatType getDensityAndVelocity(const BlockAndCell &bc,
                                  Vector3<FloatType> &velocity) const {
    auto const pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto const force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    return getDensityAndVelocity(pdf_field, force_field, bc.cell.x(),
                                 bc.cell.y(), bc.cell.z(), velocity);
  }

  FloatType getDensityAndVelocity(const PdfField *pdf_field,
                                  const VectorField *force_field,
                                  const cell_idx_t x, const cell_idx_t y,
                                  const cell_idx_t z,
                                  Vector3<FloatType> &velocity) const {
    auto const rho = lbm::accessor::DensityAndMomentumDensity::get(
        velocity, *force_field, *pdf_field, x, y, z);
    auto const invRho = FloatType{1} / rho;
    velocity *= invRho;
    return rho;
  }

  auto get_velocity_field_ptr(const IBlock *block) const {
    return block->template uncheckedFastGetData<VectorField>(
        m_velocity_field_id);
  }

  void setDensityAndVelocity(const BlockAndCell &bc,
                             Vector3<FloatType> const &velocity,
                             FloatType rho) {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    lbm::accessor::DensityAndVelocity::set(*pdf_field, bc.cell.x(), bc.cell.y(),
                                           bc.cell.z(), *force_field, velocity,
                                           rho);
  }

  Matrix3<FloatType> getPressureTensor(const BlockAndCell &bc) const {
    Matrix3<FloatType> pressureTensor;
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    lbm::accessor::PressureTensor::get(pressureTensor, *pdf_field, bc.cell.x(),
                                       bc.cell.y(), bc.cell.z());
    return pressureTensor;
  }

  FloatType pressure_tensor_correction_factor() const {
    return m_viscosity / (m_viscosity + FloatType{1} / FloatType{6});
  }

  void pressure_tensor_correction(Matrix3<FloatType> &tensor) const {
    auto const revert_factor = pressure_tensor_correction_factor();
    for (auto const i : {1, 2, 3, 5, 6, 7}) {
      tensor[i] *= revert_factor;
    }
  }

  class interpolation_illegal_access : public std::runtime_error {
  public:
    explicit interpolation_illegal_access(std::string const &field,
                                          Utils::Vector3d const &pos,
                                          std::array<int, 3> const &node,
                                          double weight)
        : std::runtime_error("Access to LB " + field + " field failed") {
      std::cerr << "pos [" << pos << "], "
                << "node [" << Utils::Vector3i(node) << "], "
                << "weight " << weight << "\n";
    }
  };

  class vtk_runtime_error : public std::runtime_error {
  public:
    explicit vtk_runtime_error(std::string const &vtk_uid,
                               std::string const &reason)
        : std::runtime_error("VTKOutput object '" + vtk_uid + "' " + reason) {}
  };

protected:
  // Member variables
  FloatType m_viscosity;
  FloatType m_density;
  FloatType m_kT;

  // Block data access handles
  BlockDataID m_pdf_field_id;
  BlockDataID m_pdf_tmp_field_id;
  BlockDataID m_flag_field_id;

  BlockDataID m_last_applied_force_field_id;
  BlockDataID m_force_to_be_applied_id;

  BlockDataID m_velocity_field_id;
  BlockDataID m_vec_tmp_field_id;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;
  using PDFStreamingCommunicator =
      blockforest::communication::UniformBufferedScheme<
          typename stencil::D3Q19>;
  std::shared_ptr<PDFStreamingCommunicator> m_pdf_streaming_communication;

  // ResetForce sweep + external force handling
  std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  // Stream sweep
  std::shared_ptr<StreamSweep> m_stream;

  // Lees Edwards boundary interpolation
  std::shared_ptr<LeesEdwardsPack> m_lees_edwards_callbacks;
  std::shared_ptr<InterpolateAndShiftAtBoundary<PdfField, FloatType>>
      m_lees_edwards_pdf_interpol_sweep;
  std::shared_ptr<InterpolateAndShiftAtBoundary<VectorField, FloatType>>
      m_lees_edwards_vel_interpol_sweep;
  std::shared_ptr<InterpolateAndShiftAtBoundary<VectorField, FloatType>>
      m_lees_edwards_last_applied_force_interpol_sweep;
  std::shared_ptr<InterpolateAndShiftAtBoundary<VectorField, FloatType>>
      m_lees_edwards_force_to_be_applied_backwards_interpol_sweep;

  // Collision sweep
  std::shared_ptr<CollisionModel> m_collision_model;

  // boundaries
  std::shared_ptr<BoundaryModel> m_boundary;

  // lattice
  std::shared_ptr<LatticeWalberla> m_lattice;

  std::size_t stencil_size() const override {
    return static_cast<std::size_t>(Stencil::Size);
  }

  /**
   * @brief Convenience function to add a field with a custom allocator.
   *
   * When vectorization is off, let waLBerla decide which memory allocator
   * to use. When vectorization is on, the aligned memory allocator is
   * required, otherwise <tt>cpu_vectorize_info["assume_aligned"]</tt> will
   * trigger assertions. That is because for single-precision kernels the
   * waLBerla heuristic in <tt>src/field/allocation/FieldAllocator.h</tt>
   * will fall back to @c StdFieldAlloc, yet @c AllocateAligned is needed
   * for intrinsics to work.
   */
  template <typename Field> auto add_to_storage(char const *tag) {
    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();
#ifdef WITH_AVX_KERNELS
#if defined(__AVX512F__)
    constexpr uint_t alignment = 64;
#elif defined(__AVX__)
    constexpr uint_t alignment = 32;
#elif defined(__SSE__)
    constexpr uint_t alignment = 16;
#else
#error "Unsupported arch, check walberla src/field/allocation/FieldAllocator.h"
#endif
    using value_type = typename Field::value_type;
    using Allocator = field::AllocateAligned<value_type, alignment>;
    auto const allocator = std::make_shared<Allocator>();
    auto const empty_set = Set<SUID>::emptySet();
    return field::addToStorage<Field>(
        blocks, tag, field::internal::defaultSize, FloatType{0}, field::fzyx,
        n_ghost_layers, false, {}, empty_set, empty_set, allocator);
#else  // WITH_AVX_KERNELS
    return field::addToStorage<Field>(blocks, tag, FloatType{0}, field::fzyx,
                                      n_ghost_layers);
#endif // WITH_AVX_KERNELS
  }

public:
  LBWalberlaImpl(std::shared_ptr<LatticeWalberla> lattice,
                 double viscosity, double density)
      : m_viscosity(FloatType_c(viscosity)), m_density(FloatType_c(density)),
        m_kT(FloatType{0}), m_lattice(std::move(lattice)) {

    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();
    if (n_ghost_layers == 0)
      throw std::runtime_error("At least one ghost layer must be used");

    // Init and register fields
    m_pdf_field_id = add_to_storage<PdfField>("pdfs");
    m_pdf_tmp_field_id = add_to_storage<PdfField>("pdfs_tmp");
    m_last_applied_force_field_id = add_to_storage<VectorField>("force field");
    m_force_to_be_applied_id = add_to_storage<VectorField>("force field");
    m_velocity_field_id = add_to_storage<VectorField>("velocity field");
    m_vec_tmp_field_id = add_to_storage<VectorField>("velocity field");

    // Init and register pdf field
    auto pdf_setter =
        InitialPDFsSetter(m_force_to_be_applied_id, m_pdf_field_id,
                          m_velocity_field_id, m_density);
    for (auto b = blocks->begin(); b != blocks->end(); ++b) {
      pdf_setter(&*b);
    }

    // Init and register flag field (fluid/boundary)
    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        blocks, "flag field", n_ghost_layers);
    // Init boundary sweep
    reset_boundary_handling();

    // Set up the communication and register fields
    m_pdf_streaming_communication =
        std::make_shared<PDFStreamingCommunicator>(blocks);
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, n_ghost_layers));
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, n_ghost_layers));
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<FlagField>>(
            m_flag_field_id, n_ghost_layers));

    m_full_communication = std::make_shared<FullCommunicator>(blocks);
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_velocity_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<FlagField>>(
            m_flag_field_id, n_ghost_layers));

    // Instance the sweep responsible for force double buffering and
    // external forces
    m_reset_force = std::make_shared<ResetForce<PdfField, VectorField>>(
        m_last_applied_force_field_id, m_force_to_be_applied_id);

    // Prepare LB sweeps
    // Note: For now, combined collide-stream sweeps cannot be used,
    // because the collide-push variant is not supported by lbmpy.
    // The following functors are individual in-place collide and stream steps
    m_stream = std::make_shared<StreamSweep>(
        m_last_applied_force_field_id, m_pdf_field_id, m_velocity_field_id);
  }

private:
  void integrate_stream(std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_stream)(&*b);
  }

  void integrate_collide(std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      boost::apply_visitor(run_collide_sweep, *m_collision_model,
                           boost::variant<IBlock *>(&*b));
    if (auto *cm = boost::get<CollisionModelThermalized>(&*m_collision_model)) {
      cm->time_step_++;
    }
  }

  bool lees_edwards_bc() {
    return boost::get<CollisionModelLeesEdwards>(&*m_collision_model);
  }

  void apply_lees_edwards_pdf_interpolation(
      std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_lees_edwards_pdf_interpol_sweep)(&*b);
  }

  void apply_lees_edwards_vel_interpolation_and_shift(
      std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_lees_edwards_vel_interpol_sweep)(&*b);
  }

  void apply_lees_edwards_last_applied_force_interpolation(
      std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_lees_edwards_last_applied_force_interpol_sweep)(&*b);
  }

  void apply_lees_edwards_force_to_be_applied_backwards_interpolation(
      std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_lees_edwards_force_to_be_applied_backwards_interpol_sweep)(&*b);
  }

  void integrate_reset_force(std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_reset_force)(&*b);
  }

  void integrate_boundaries(std::shared_ptr<Lattice_T> const &blocks) {
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_boundary)(&*b);
  }

  void integrate_push_scheme() {
    auto const &blocks = get_lattice().get_blocks();
    // Reset force fields
    integrate_reset_force(blocks);
    // LB collide
    integrate_collide(blocks);
    (*m_pdf_streaming_communication).communicate();
    // Handle boundaries
    integrate_boundaries(blocks);
    // LB stream
    integrate_stream(blocks);
    // Refresh ghost layers
    (*m_full_communication).communicate();
  }

  void integrate_pull_scheme() {
    auto const &blocks = get_lattice().get_blocks();
    integrate_reset_force(blocks);
    // Handle boundaries
    integrate_boundaries(blocks);
    // LB stream
    integrate_stream(blocks);
    // LB collide
    integrate_collide(blocks);

    // Refresh ghost layers
    ghost_communication();
  }

protected:
  void integrate_vtk_writers() override {
    for (auto const &it : m_vtk_auto) {
      auto &vtk_handle = it.second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
  }

public:
  void integrate() override {
    reallocate_ubb_field();
    if (lees_edwards_bc()) {
      integrate_pull_scheme();
    } else {
      integrate_push_scheme();
    }
    // Handle VTK writers
    integrate_vtk_writers();
  }

  void ghost_communication() override {
    (*m_full_communication).communicate();
    if (lees_edwards_bc()) {
      auto const &blocks = get_lattice().get_blocks();
      apply_lees_edwards_pdf_interpolation(blocks);
      apply_lees_edwards_vel_interpolation_and_shift(blocks);
      apply_lees_edwards_last_applied_force_interpolation(blocks);
    }
  }

  void set_collision_model(double kT, unsigned int seed) override {
    auto const omega = shear_mode_relaxation_rate();
    auto const omega_odd = odd_mode_relaxation_rate(omega);
    m_kT = FloatType_c(kT);
    auto obj = CollisionModelThermalized(
        m_last_applied_force_field_id, m_pdf_field_id, uint32_t{0u},
        uint32_t{0u}, uint32_t{0u}, m_kT, omega, omega, omega_odd, omega, seed,
        uint32_t{0u});
    obj.block_offset_generator =
        [this](IBlock *const block, uint32_t &block_offset_0,
               uint32_t &block_offset_1, uint32_t &block_offset_2) {
          auto const &blocks = get_lattice().get_blocks();
          block_offset_0 = blocks->getBlockCellBB(*block).xMin();
          block_offset_1 = blocks->getBlockCellBB(*block).yMin();
          block_offset_2 = blocks->getBlockCellBB(*block).zMin();
        };
    m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
  }

  void set_collision_model(
      std::unique_ptr<LeesEdwardsPack> &&lees_edwards_pack) override {
    assert(m_kT == 0.);
    auto const shear_direction = lees_edwards_pack->shear_direction;
    auto const shear_plane_normal = lees_edwards_pack->shear_plane_normal;
    auto const shear_vel = FloatType_c(lees_edwards_pack->get_shear_velocity());
    auto const omega = shear_mode_relaxation_rate();
    if (shear_plane_normal != 1) {
      throw std::runtime_error(
          "Lees-Edwards LB only supports shear_plane_normal=\"y\"");
    }
    auto obj = CollisionModelLeesEdwards(
        m_last_applied_force_field_id, m_pdf_field_id,
        FloatType_c(get_lattice().get_grid_dimensions()[shear_plane_normal]),
        omega, shear_vel);
    m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
    m_lees_edwards_callbacks = std::move(lees_edwards_pack);
    run_collide_sweep.register_lees_edwards_callbacks(m_lees_edwards_callbacks);
    m_lees_edwards_pdf_interpol_sweep =
        std::make_shared<InterpolateAndShiftAtBoundary<PdfField, FloatType>>(
            get_lattice().get_blocks(), m_pdf_field_id, m_pdf_tmp_field_id,
            get_lattice().get_ghost_layers(), shear_direction,
            shear_plane_normal, m_lees_edwards_callbacks->get_pos_offset);
    m_lees_edwards_vel_interpol_sweep =
        std::make_shared<InterpolateAndShiftAtBoundary<VectorField, FloatType>>(
            get_lattice().get_blocks(), m_velocity_field_id, m_vec_tmp_field_id,
            get_lattice().get_ghost_layers(), shear_direction,
            shear_plane_normal, m_lees_edwards_callbacks->get_pos_offset,
            m_lees_edwards_callbacks->get_shear_velocity);
    m_lees_edwards_last_applied_force_interpol_sweep =
        std::make_shared<InterpolateAndShiftAtBoundary<VectorField, FloatType>>(
            get_lattice().get_blocks(), m_last_applied_force_field_id,
            m_vec_tmp_field_id, get_lattice().get_ghost_layers(),
            shear_direction, shear_plane_normal,
            m_lees_edwards_callbacks->get_pos_offset);
    m_lees_edwards_force_to_be_applied_backwards_interpol_sweep =
        std::make_shared<InterpolateAndShiftAtBoundary<VectorField, FloatType>>(
            get_lattice().get_blocks(), m_force_to_be_applied_id,
            m_vec_tmp_field_id, get_lattice().get_ghost_layers(),
            shear_direction, shear_plane_normal, [this]() {
              return -1.0 * m_lees_edwards_callbacks->get_pos_offset();
            });
  }

  void check_lebc(int shear_direction, int shear_plane_normal) const override {
    if (m_lees_edwards_callbacks) {
      if (m_lees_edwards_callbacks->shear_direction != shear_direction or
          m_lees_edwards_callbacks->shear_plane_normal != shear_plane_normal) {
        throw std::runtime_error(
            "MD and LB Lees-Edwards boundary conditions disagree");
      }
    }
  }

  void set_viscosity(double viscosity) override {
    m_viscosity = FloatType_c(viscosity);
  }

  double get_viscosity() const override {
    return numeric_cast<double>(m_viscosity);
  }

  double get_density() const override {
    return numeric_cast<double>(m_density);
  }

  // Velocity
  boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const override {
    auto const is_boundary = get_node_is_boundary(node, consider_ghosts);
    if (is_boundary)    // is info available locally
      if (*is_boundary) // is the node a boundary
        return get_node_velocity_at_boundary(node);
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {};
    auto const vel_field = get_velocity_field_ptr(bc->block);
    return Utils::Vector3d{double_c(vel_field->get((*bc).cell, uint_t(0u))),
                           double_c(vel_field->get((*bc).cell, uint_t(1u))),
                           double_c(vel_field->get((*bc).cell, uint_t(2u)))};
  }
  bool set_node_velocity(const Utils::Vector3i &node,
                         const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;
    // We have to set both, the pdf and the stored velocity field
    auto const density = getDensity(*bc);
    auto const vel = to_vector3<FloatType>(v);
    setDensityAndVelocity(*bc, vel, density);
    auto vel_field =
        (*bc).block->template getData<VectorField>(m_velocity_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      vel_field->get((*bc).cell, f) = FloatType_c(v[f]);
    }

    return true;
  }
  boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &pos,
                      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return {};
    Utils::Vector3d v{0.0, 0.0, 0.0};
    interpolate_bspline_at_pos(
        pos, [this, &v, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto const res = get_node_velocity(Utils::Vector3i(node), true);
            if (!res) {
              throw interpolation_illegal_access("velocity", pos, node, weight);
            }
            v += *res * weight;
          }
        });
    return {v};
  }

  boost::optional<double> get_interpolated_density_at_pos(
      const Utils::Vector3d &pos,
      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return {};
    double dens = 0.0;
    interpolate_bspline_at_pos(
        pos, [this, &dens, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto const res = get_node_density(Utils::Vector3i(node), true);
            if (!res) {
              throw interpolation_illegal_access("density", pos, node, weight);
            }
            dens += *res * weight;
          }
        });
    return {dens};
  }

  // Local force
  bool add_force_at_pos(const Utils::Vector3d &pos,
                        const Utils::Vector3d &force) override {
    if (!m_lattice->pos_in_local_halo(pos))
      return false;
    auto force_at_node = [this, force](const std::array<int, 3> node,
                                       double weight) {
      auto const bc =
          get_block_and_cell(get_lattice(), Utils::Vector3i(node), true);
      if (bc) {
        auto force_field =
            (*bc).block->template uncheckedFastGetData<VectorField>(
                m_force_to_be_applied_id);
        for (int i : {0, 1, 2})
          force_field->get((*bc).cell, i) += FloatType_c(force[i] * weight);
      }
    };
    interpolate_bspline_at_pos(pos, force_at_node);
    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return {};

    auto const &force_field =
        (*bc).block->template getData<VectorField>(m_force_to_be_applied_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  bool set_node_last_applied_force(Utils::Vector3i const &node,
                                   Utils::Vector3d const &force) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    auto force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      force_field->get((*bc).cell, f) = FloatType_c(force[f]);
    }

    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {};

    auto const force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  // Population
  bool set_node_pop(Utils::Vector3i const &node,
                    std::vector<double> const &population) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    assert(population.size() == Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      pdf_field->get((*bc).cell, f) = FloatType_c(population[f]);
    }

    return true;
  }

  boost::optional<std::vector<double>>
  get_node_pop(const Utils::Vector3i &node,
               bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
    std::vector<double> population(Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      population[f] = double_c(pdf_field->get((*bc).cell, f));
    }

    return {population};
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node, double density) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    Vector3<FloatType> vel;
    getDensityAndVelocity(*bc, vel);
    setDensityAndVelocity(*bc, vel, FloatType_c(density));

    return true;
  }

  boost::optional<double>
  get_node_density(const Utils::Vector3i &node,
                   bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    auto const density = getDensity(*bc);
    return {double_c(density)};
  }

  boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc or !m_boundary->node_is_boundary(node))
      return {boost::none};

    return {m_boundary->get_node_value_at_boundary(node)};
  }

  bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                     const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary->set_node_value_at_boundary(node, v, *bc);

    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc or !m_boundary->node_is_boundary(node))
      return {boost::none};

    return get_node_last_applied_force(node, true);
  }

  bool remove_node_from_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary->remove_node_from_boundary(node, *bc);

    return true;
  }

  boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary->node_is_boundary(node)};
  }

  void reallocate_ubb_field() override { m_boundary->boundary_update(); }

  void clear_boundaries() override { reset_boundary_handling(); }

  void
  update_boundary_from_shape(std::vector<int> const &raster_flat,
                             std::vector<double> const &data_flat) override {
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const data = fill_3D_vector_array(data_flat, grid_size);
    auto const field_getter = [field_id = m_pdf_field_id](auto &block) {
      return block->template getData<PdfField>(field_id);
    };
    set_boundary_from_grid(*m_boundary, field_getter, get_lattice(),
                           raster_flat, data);
  }

  // Pressure tensor
  boost::optional<Utils::VectorXd<9>>
  get_node_pressure_tensor(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return {boost::none};
    auto tensor = getPressureTensor(*bc);
    pressure_tensor_correction(tensor);
    return to_vector9d(tensor);
  }

  // Global pressure tensor
  Utils::VectorXd<9> get_pressure_tensor() const override {
    auto const &blocks = get_lattice().get_blocks();
    Matrix3<FloatType> tensor(FloatType{0});
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
        Matrix3<FloatType> local_p;
        lbm::accessor::PressureTensor::get(local_p, *pdf_field, x, y, z);
        tensor += local_p;
      });
    }
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const number_of_nodes = Utils::product(grid_size);
    pressure_tensor_correction(tensor);
    return to_vector9d(tensor) * (1. / static_cast<double>(number_of_nodes));
  }

  // Global momentum
  Utils::Vector3d get_momentum() const override {
    auto const &blocks = get_lattice().get_blocks();
    Vector3<FloatType> mom;
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      auto force_field =
          block->template getData<VectorField>(m_last_applied_force_field_id);
      Vector3<FloatType> local_v;
      WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
        FloatType local_dens =
            getDensityAndVelocity(pdf_field, force_field, x, y, z, local_v);
        mom += local_dens * local_v;
      });
    }
    return to_vector3d(mom);
  }

  // Global external force
  void set_external_force(const Utils::Vector3d &ext_force) override {
    m_reset_force->set_ext_force(ext_force);
  }
  Utils::Vector3d get_external_force() const override {
    return m_reset_force->get_ext_force();
  }

  double get_kT() const override { return numeric_cast<double>(m_kT); }

  uint64_t get_rng_state() const override {
    auto const *cm = boost::get<CollisionModelThermalized>(&*m_collision_model);
    if (!cm or m_kT == 0.)
      throw std::runtime_error("The LB does not use a random number generator");
    return cm->time_step_;
  }
  void set_rng_state(uint64_t counter) override {
    auto *cm = boost::get<CollisionModelThermalized>(&*m_collision_model);
    if (!cm or m_kT == 0.)
      throw std::runtime_error("The LB does not use a random number generator");
    cm->time_step_ = counter;
  }

  LatticeWalberla const &get_lattice() const noexcept override {
    return *m_lattice;
  }

  [[nodiscard]] std::size_t get_velocity_field_id() const override {
    return m_velocity_field_id;
  }

  [[nodiscard]] std::size_t get_force_field_id() const override {
    return m_force_to_be_applied_id;
  }

  void register_vtk_field_filters(walberla::vtk::VTKOutput &vtk_obj) override {
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_id);
    fluid_filter.addFlag(Boundary_flag);
    vtk_obj.addCellExclusionFilter(fluid_filter);
  }

  void register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                                  int flag_observables) override {
    if (flag_observables & static_cast<int>(OutputVTK::density)) {
      vtk_obj.addCellDataWriter(
          make_shared<lbm::DensityVTKWriter<LBWalberlaImpl, float>>(
              m_pdf_field_id, "DensityFromPDF"));
    }
    if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
      vtk_obj.addCellDataWriter(
          make_shared<field::VTKWriter<VectorField, float>>(
              m_velocity_field_id, "VelocityFromVelocityField"));
    }
    if (flag_observables & static_cast<int>(OutputVTK::pressure_tensor)) {
      vtk_obj.addCellDataWriter(
          make_shared<lbm::PressureTensorVTKWriter<LBWalberlaImpl, float>>(
              m_pdf_field_id, "PressureTensorFromPDF",
              pressure_tensor_correction_factor()));
    }
  }

  ~LBWalberlaImpl() override = default;
};

} // namespace walberla
