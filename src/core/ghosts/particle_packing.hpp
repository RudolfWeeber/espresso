/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

/** \file
 *  Reusable particle packing/unpacking for ghost communications.
 *
 *  Extracted from ghosts.cpp so that both the legacy ghost_communicator()
 *  and future async engines can share the same byte-identical serialization.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "ParticleList.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <span>
#include <vector>

namespace GhostComm {

/**
 * Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 */
class CommBuf {
public:
  /** Returns a pointer to the non-bond storage. */
  char *data() { return buf.data(); }
  const char *data() const { return buf.data(); }

  /** Returns the number of elements in the non-bond storage. */
  std::size_t size() const { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(std::size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage. */
  auto &bonds() { return bondbuf; }
  const auto &bonds() const { return bondbuf; }

  auto make_span() { return std::span(buf.data(), buf.size()); }

private:
  std::vector<char> buf;     ///< Buffer for everything but bonds
  std::vector<char> bondbuf; ///< Buffer for bond lists
};

/**
 * @brief Calculate the per-particle transmit size for the given data parts.
 */
std::size_t calc_transmit_size(BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Calculate the total transmit size for a list of cells and data parts.
 *
 * When GHOSTTRANS_PARTNUM is set, returns sizeof(unsigned int) per cell.
 * Otherwise returns the total number of particles times the per-particle size.
 */
std::size_t calc_transmit_size(std::span<ParticleList *const> cells,
                               BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Pack particle data from cells into a communication buffer.
 *
 * Equivalent to prepare_send_buffer() in the original ghosts.cpp.
 * Handles GHOSTTRANS_PARTNUM (writes cell sizes) and GHOSTTRANS_BONDS
 * (separate bond buffer) special cases.
 *
 * @param buf        Buffer to pack into (resized as needed).
 * @param cells      Source particle lists.
 * @param shift      Ghost shift applied to particle positions on save.
 * @param box_geo    Box geometry for fold_position.
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void pack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                unsigned data_parts);

/**
 * @brief Unpack particle data from a communication buffer into cells.
 *
 * Equivalent to put_recv_buffer() in the original ghosts.cpp.
 * Handles GHOSTTRANS_PARTNUM (resizes ghost cells) and GHOSTTRANS_BONDS
 * (separate bond buffer) special cases.
 *
 * @param buf        Buffer to unpack from.
 * @param cells      Destination particle lists.
 * @param box_geo    Box geometry (unused here, kept for API symmetry).
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void unpack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                  BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Add forces from a communication buffer to particles in cells.
 *
 * Equivalent to add_forces_from_recv_buffer() in the original ghosts.cpp.
 */
void add_forces(CommBuf &buf, std::span<ParticleList *const> cells);

#ifdef ESPRESSO_BOND_CONSTRAINT
/**
 * @brief Add rattle corrections from a communication buffer to particles.
 *
 * Equivalent to add_rattle_correction_from_recv_buffer() in ghosts.cpp.
 */
void add_rattle(CommBuf &buf, std::span<ParticleList *const> cells);
#endif

/**
 * @brief Copy particle data from src to dst applying a ghost shift.
 *
 * Corresponds to one src/dst pair in cell_cell_transfer() from ghosts.cpp.
 * For GHOSTTRANS_PARTNUM, resizes dst to match src.
 * Otherwise serializes each particle from src and deserializes into dst,
 * applying the shift and folding the position. Bond data is copied directly.
 *
 * @param src        Source particle list.
 * @param dst        Destination particle list.
 * @param shift      Ghost shift to apply.
 * @param box_geo    Box geometry for fold_position.
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void local_cell_copy(ParticleList &src, ParticleList &dst,
                     Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                     unsigned data_parts);

} // namespace GhostComm
