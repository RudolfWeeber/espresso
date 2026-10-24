/*
 * Copyright (C) 2026 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace IO {

/**
 * @brief Parallel read access to a 3D voxel mask stored in an HDF5 dataset.
 *
 * The dataset layout is @c [x][y][z] in C order. Nonzero values denote
 * solid voxels. The file is opened read-only with the MPI-IO driver.
 * The constructor is collective; every rank throws on failure.
 * Reads are independent, so ranks may issue different numbers of reads.
 */
class VoxelMaskFile {
public:
  /**
   * @throw std::runtime_error naming the path or dataset when either is
   * missing or the dataset is not 3D.
   */
  VoxelMaskFile(std::string const &path, std::string const &dataset,
                boost::mpi::communicator const &comm);
  ~VoxelMaskFile();
  VoxelMaskFile(VoxelMaskFile const &) = delete;
  VoxelMaskFile &operator=(VoxelMaskFile const &) = delete;

  /** @brief Extent of the dataset. */
  Utils::Vector3i shape() const;

  /**
   * @brief Independent read of the box [offset, offset + count).
   * Result is row-major (x slowest), values are converted to 0/1.
   */
  std::vector<std::uint8_t> read(Utils::Vector3i const &offset,
                                 Utils::Vector3i const &count) const;

private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

} // namespace IO
