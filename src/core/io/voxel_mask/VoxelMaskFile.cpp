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

#include "io/writer/hdf5_patches.hpp" // must appear first

#include "VoxelMaskFile.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <highfive/highfive.hpp>

#include <H5Tpublic.h>
#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace IO {

struct VoxelMaskFile::Impl {
  HighFive::File file;
  HighFive::DataSet dataset;
  Utils::Vector3i shape;
  H5T_class_t type_class;
  std::size_t type_size;
  std::string description;

  Impl(std::string const &path, std::string const &name,
       boost::mpi::communicator const &comm)
      : file(open_file(path, comm)), dataset(open_dataset(path, name)),
        description("dataset '" + name + "' in file '" + path + "'") {
    auto const dims = dataset.getSpace().getDimensions();
    if (dims.size() != 3u) {
      throw std::runtime_error("The " + description + " must be 3D, got " +
                               std::to_string(dims.size()) + " dimensions");
    }
    for (int i = 0; i < 3; ++i) {
      shape[i] = static_cast<int>(dims[i]);
    }
    auto const type = dataset.getDataType();
    type_class = H5Tget_class(type.getId());
    type_size = type.getSize();
    if (type_class != H5T_INTEGER and type_class != H5T_ENUM and
        type_class != H5T_FLOAT) {
      throw std::runtime_error("The " + description +
                               " has an unsupported data type");
    }
    if (type_class == H5T_ENUM and type_size != 1u) {
      throw std::runtime_error("The " + description +
                               " is an enum with more than one byte");
    }
  }

private:
  HighFive::File open_file(std::string const &path,
                           boost::mpi::communicator const &comm) {
    HighFive::FileAccessProps fapl;
    fapl.add(HighFive::MPIOFileAccess(MPI_Comm(comm), MPI_INFO_NULL));
    try {
      return HighFive::File(path, HighFive::File::ReadOnly, fapl);
    } catch (HighFive::Exception const &) {
      throw std::runtime_error("Cannot open HDF5 file '" + path + "'");
    }
  }

  HighFive::DataSet open_dataset(std::string const &path,
                                 std::string const &name) {
    try {
      return file.getDataSet(name);
    } catch (HighFive::Exception const &) {
      throw std::runtime_error("Cannot find dataset '" + name + "' in file '" +
                               path + "'");
    }
  }
};

VoxelMaskFile::VoxelMaskFile(std::string const &path,
                             std::string const &dataset,
                             boost::mpi::communicator const &comm)
    : m_impl(std::make_unique<Impl>(path, dataset, comm)) {}

VoxelMaskFile::~VoxelMaskFile() = default;

Utils::Vector3i VoxelMaskFile::shape() const { return m_impl->shape; }

template <typename T>
static void to_mask(std::vector<T> const &raw,
                    std::vector<std::uint8_t> &mask) {
  mask.resize(raw.size());
  std::ranges::transform(raw, mask.begin(), [](T value) {
    return static_cast<std::uint8_t>(value != T{0});
  });
}

std::vector<std::uint8_t>
VoxelMaskFile::read(Utils::Vector3i const &offset,
                    Utils::Vector3i const &count) const {
  auto const &impl = *m_impl;
  std::size_t total = 1;
  for (int i = 0; i < 3; ++i) {
    if (offset[i] < 0 or count[i] < 0 or offset[i] + count[i] > impl.shape[i]) {
      throw std::runtime_error("Requested region lies outside of the " +
                               impl.description);
    }
    total *= static_cast<std::size_t>(count[i]);
  }
  std::vector<std::uint8_t> mask;
  if (total == 0u) {
    return mask;
  }
  std::vector<std::size_t> const start{static_cast<std::size_t>(offset[0]),
                                       static_cast<std::size_t>(offset[1]),
                                       static_cast<std::size_t>(offset[2])};
  std::vector<std::size_t> const extent{static_cast<std::size_t>(count[0]),
                                        static_cast<std::size_t>(count[1]),
                                        static_cast<std::size_t>(count[2])};
  auto const selection = impl.dataset.select(start, extent);
  if (impl.type_class == H5T_FLOAT) {
    std::vector<double> raw(total);
    selection.read_raw(raw.data(), HighFive::create_datatype<double>());
    to_mask(raw, mask);
  } else if (impl.type_size == 1u) {
    // one-byte integers and enums (e.g. h5py bool): the file type is used as
    // memory type, since HDF5 cannot convert enums to integers
    std::vector<std::int8_t> raw(total);
    selection.read_raw(raw.data(), impl.dataset.getDataType());
    to_mask(raw, mask);
  } else {
    std::vector<std::int64_t> raw(total);
    selection.read_raw(raw.data(), HighFive::create_datatype<std::int64_t>());
    to_mask(raw, mask);
  }
  return mask;
}

} // namespace IO
