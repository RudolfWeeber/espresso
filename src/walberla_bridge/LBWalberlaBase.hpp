/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#ifndef LB_WALBERLA_BASE_HPP
#define LB_WALBERLA_BASE_HPP

/**
 * @file
 * @ref LBWalberlaBase provides the public interface of the LB
 * waLBerla bridge. It relies on type erasure to hide the waLBerla
 * implementation details from the ESPResSo core. It is implemented
 * by @ref walberla::LBWalberlaImpl.
 */

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

/** Class that runs and controls the LB on WaLBerla
 */
class LBWalberlaBase {
public:
  /** @brief Integrate LB for one time step */
  virtual void integrate() = 0;

  /** @brief perform ghost communication of PDF and applied forces */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized velocities in the PDF */
  virtual size_t stencil_size() const = 0;

  virtual boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const = 0;
  virtual bool set_node_velocity(const Utils::Vector3i &node,
                                 const Utils::Vector3d &v) = 0;
  virtual boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &position,
                      bool consider_points_in_halo = false) const = 0;
  virtual boost::optional<double> get_interpolated_density_at_pos(
      const Utils::Vector3d &position,
      bool consider_points_in_halo = false) const = 0;

  virtual bool add_force_at_pos(const Utils::Vector3d &position,
                                const Utils::Vector3d &force) = 0;
  /** @brief Get (stored) force to be applied on node in next time step */
  virtual boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const = 0;

  /** @brief Get (stored) force that was applied on node in last time step */
  virtual boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const = 0;

  virtual bool set_node_last_applied_force(const Utils::Vector3i &node,
                                           const Utils::Vector3d &force) = 0;

  // Population
  virtual boost::optional<std::vector<double>>
  get_node_pop(const Utils::Vector3i &node) const = 0;
  virtual bool set_node_pop(const Utils::Vector3i &node,
                            std::vector<double> const &population) = 0;

  // Density
  virtual bool set_node_density(const Utils::Vector3i &node,
                                double density) = 0;
  virtual boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const = 0;

  /** @brief Get velocity for node with velocity boundary condition */
  virtual boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const = 0;

  /** @brief Set a node as velocity boundary condition and assign boundary
   *  velocity
   */
  virtual bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                             const Utils::Vector3d &v) = 0;
  /** @brief Get (stored) force applied on node due to boundary condition */
  virtual boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const = 0;
  virtual bool remove_node_from_boundary(const Utils::Vector3i &node) = 0;
  virtual boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const = 0;
  virtual void clear_boundaries() = 0;

  // Pressure tensor
  virtual boost::optional<Utils::Vector6d>
  get_node_pressure_tensor(const Utils::Vector3i &node) const = 0;

  /** @brief Calculate momentum summed over all nodes on the MPI rank */
  virtual Utils::Vector3d get_momentum() const = 0;

  virtual void set_external_force(const Utils::Vector3d &ext_force) = 0;
  virtual Utils::Vector3d get_external_force() const = 0;

  // Global parameters
  virtual void set_viscosity(double viscosity) = 0;
  virtual double get_viscosity() const = 0;
  virtual double get_kT() const = 0;

  //* @brief Fet the rng counter for thermalized LBs */
  virtual uint64_t get_rng_state() const = 0;

  /** @brief set the rng state of thermalized LBs */
  virtual void set_rng_state(uint64_t counter) = 0;

  // Grid, domain, halo
  virtual int n_ghost_layers() const = 0;
  virtual Utils::Vector3i get_grid_dimensions() const = 0;
  virtual std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const = 0;
  virtual bool node_in_local_domain(const Utils::Vector3i &node) const = 0;
  virtual bool node_in_local_halo(const Utils::Vector3i &node) const = 0;
  virtual bool pos_in_local_domain(const Utils::Vector3d &pos) const = 0;
  virtual bool pos_in_local_halo(const Utils::Vector3d &pos) const = 0;

  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *         otherwise add a callback to write every @p delta_N LB steps
   *         to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *         @ref OutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  virtual void create_vtk(unsigned delta_N, unsigned initial_count,
                          unsigned flag_observables,
                          std::string const &identifier,
                          std::string const &base_folder,
                          std::string const &prefix) = 0;
  /** @brief Write a VTK observable to disk.
   *
   *  @param vtk_uid          Name of the VTK object
   */
  virtual void write_vtk(std::string const &vtk_uid) = 0;
  /** @brief Toggle a VTK observable on/off.
   *
   *  @param vtk_uid          Name of the VTK object
   *  @param status           1 to switch on, 0 to switch off
   */
  virtual void switch_vtk(std::string const &vtk_uid, int status) = 0;

  /** @brief return a pairs of global node index and node center position */
  virtual std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts) const = 0;
  virtual ~LBWalberlaBase() = default;
};

/** @brief LB statistics to write to VTK files */
enum class OutputVTK : unsigned {
  density = 1u << 0u,
  velocity_vector = 1u << 1u,
  pressure_tensor = 1u << 2u,
};

#endif // LB_WALBERLA_H

