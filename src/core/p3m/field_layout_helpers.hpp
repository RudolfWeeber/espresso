#include "utils/Vector.hpp"
#include "utils/index.hpp"
#include <stdexcept>
#include <vector>

// Function to extract a 3D block from the halo field
template <typename T>
std::vector<T> extract_block(const std::vector<T> &in_array,
                             Utils::Vector3i dimensions, Utils::Vector3i start, Utils::Vector3i stop) {
  // Extract the dimensions
  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  // Validate input

  // Calculate the size of the block excluding halo regions
  Utils::Vector3i block_dim = stop-start;

  // Output vector to hold the block
  std::vector<T> out_array(block_dim[0] * block_dim[1] * block_dim[2]);

  // Extract the block
  for (int z = 0; z < block_dim[2]; ++z) {
    for (int y = 0; y < block_dim[1]; ++y) {
      for (int x = 0; x < block_dim[0]; ++x) {
        // Compute indices for input and output arrays
        int in_index =
            Utils::get_linear_index(x + start[0], y + start[1], z + start[2],
                                    dimensions, Utils::MemoryOrder::ROW_MAJOR);

        int out_index = Utils::get_linear_index(x, y, z, block_dim,
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
                              Utils::Vector3i cropped_dim, Utils::Vector3i pad_left, Utils::Vector3i pad_right) {

  // Output vector to hold the padded field (initialized with zeros)
  Utils::Vector3i padded_dim = cropped_dim+pad_left+pad_right;
  std::vector<T> padded_array(padded_dim[0]*padded_dim[1]*padded_dim[2]);

  // Fill in the original cropped field into the padded array
  for (int z = 0; z < cropped_dim[2]; ++z) {
    for (int y = 0; y < cropped_dim[1]; ++y) {
      for (int x = 0; x < cropped_dim[0]; ++x) {
        // Compute indices for cropped and padded arrays
        int cropped_index = Utils::get_linear_index(
            x, y, z, cropped_dim, Utils::MemoryOrder::ROW_MAJOR);
        int padded_index = Utils::get_linear_index(
            x + pad_left[0], y + pad_left[1], z + pad_left[2], padded_dim,
            Utils::MemoryOrder::ROW_MAJOR);

        // Copy the value
        padded_array[padded_index] = cropped_array[cropped_index];
      }
    }
  }

  return padded_array;
}
