#include "utils/Vector.hpp"
#include "utils/index.hpp"
#include <stdexcept>
#include <vector>

// Function to extract a 3D block from the halo field
template <typename T>
std::vector<T> extract_block(const std::vector<T> &in_array,
                             Utils::Vector3i dimensions, int n_halo) {
  // Extract the dimensions
  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  // Validate input
  if (nx <= 2 * n_halo || ny <= 2 * n_halo || nz <= 2 * n_halo) {
    throw std::invalid_argument(
        "n_halo is too large for the given dimensions.");
  }

  // Calculate the size of the block excluding halo regions
  int nx_block = nx - 2 * n_halo;
  int ny_block = ny - 2 * n_halo;
  int nz_block = nz - 2 * n_halo;
  Utils::Vector3i block_dimensions = {nx_block, ny_block, nz_block};

  // Output vector to hold the block
  std::vector<T> out_array(nx_block * ny_block * nz_block);

  // Extract the block
  for (int z = 0; z < nz_block; ++z) {
    for (int y = 0; y < ny_block; ++y) {
      for (int x = 0; x < nx_block; ++x) {
        // Compute indices for input and output arrays
        int in_index =
            Utils::get_linear_index(x + n_halo, y + n_halo, z + n_halo,
                                    dimensions, Utils::MemoryOrder::ROW_MAJOR);

        int out_index = Utils::get_linear_index(x, y, z, block_dimensions,
                                                Utils::MemoryOrder::ROW_MAJOR);

        // Copy the value
        out_array[out_index] = in_array[in_index];
      }
    }
  }

  return out_array;
}

// Function to pad the 3D cropped field with zeros to restore the halo regions
template <typename T>
std::vector<T> pad_with_zeros(const std::vector<T> &cropped_array,
                              Utils::Vector3i cropped_dimensions, int n_halo) {
  // Extract the dimensions of the cropped field
  int nx_block = cropped_dimensions[0];
  int ny_block = cropped_dimensions[1];
  int nz_block = cropped_dimensions[2];

  // Calculate the dimensions of the padded field (original size with halo
  // regions)
  int nx = nx_block + 2 * n_halo;
  int ny = ny_block + 2 * n_halo;
  int nz = nz_block + 2 * n_halo;

  // Output vector to hold the padded field (initialized with zeros)
  std::vector<T> padded_array(nx * ny * nz, 0.0);
  Utils::Vector3i padded_dimensions = {nx, ny, nz};

  // Fill in the original cropped field into the padded array
  for (int z = 0; z < nz_block; ++z) {
    for (int y = 0; y < ny_block; ++y) {
      for (int x = 0; x < nx_block; ++x) {
        // Compute indices for cropped and padded arrays
        int cropped_index = Utils::get_linear_index(
            x, y, z, cropped_dimensions, Utils::MemoryOrder::ROW_MAJOR);
        int padded_index = Utils::get_linear_index(
            x + n_halo, y + n_halo, z + n_halo, padded_dimensions,
            Utils::MemoryOrder::ROW_MAJOR);

        // Copy the value
        padded_array[padded_index] = cropped_array[cropped_index];
      }
    }
  }

  return padded_array;
}
