#include "utils/Vector.hpp"
#include "utils/index.hpp"
#include <stdexcept>
#include <vector>
#include <complex>

// Function to extract a 3D block from the halo field
template <typename Container>
auto extract_block(const Container &in_array,
                             Utils::Vector3i dimensions, Utils::Vector3i start, Utils::Vector3i stop) {
  // Extract the dimensions

  // Validate input

  // Calculate the size of the block excluding halo regions
  Utils::Vector3i block_dim = stop-start;

  // Output vector to hold the block
  std::vector<typename Container::value_type> out_array(block_dim[0] * block_dim[1] * block_dim[2]);

  // Extract the block
      for (int x = 0; x < block_dim[0]; ++x) {
    for (int y = 0; y < block_dim[1]; ++y) {
  for (int z = 0; z < block_dim[2]; ++z) {
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

  return std::move(out_array);
}

// Function to pad the 3D cropped field with zeros to restore the halo regions
#include <vector>
#include <algorithm>

template <typename T>
std::vector<T> pad_with_zeros(const std::vector<T> &cropped_array,
                              Utils::Vector3i cropped_dim, Utils::Vector3i pad_left, Utils::Vector3i pad_right) {

  // Calculate dimensions and strides
  Utils::Vector3i padded_dim = cropped_dim + pad_left + pad_right;
  int cropped_xy_stride = cropped_dim[1] * cropped_dim[2];
  int padded_xy_stride = padded_dim[1] * padded_dim[2];

  // Output vector to hold the padded field (initialized with zeros)
  std::vector<T> padded_array(padded_dim[0] * padded_dim[1] * padded_dim[2]);

  // Calculate the starting position in the padded array for the inner field
  int padded_start_x = pad_left[0] * padded_xy_stride;
  int padded_start_y = pad_left[1] * padded_dim[2] + pad_left[2];

  // Fill in the original cropped field into the padded array by chunks
  for (int x = 0; x < cropped_dim[0]; ++x) {
    int cropped_x_offset = x * cropped_xy_stride;
    int padded_x_offset = padded_start_x + x * padded_xy_stride;

    for (int y = 0; y < cropped_dim[1]; ++y) {
      int cropped_y_offset = cropped_x_offset + y * cropped_dim[2];
      int padded_y_offset = padded_x_offset + y * padded_dim[2] + padded_start_y;

      // Copy a contiguous slice of the z-dimension at once
      std::copy(cropped_array.begin() + cropped_y_offset,
                cropped_array.begin() + cropped_y_offset + cropped_dim[2],
                padded_array.begin() + padded_y_offset);
    }
  }

  return std::move(padded_array);
}




template <typename T>
std::vector<T> discard_imaginary_part(std::vector<std::complex<T>>& v) {
   std::vector<T> res;
   res.reserve(v.size());
   for (int i=0;i<v.size();i++) {
    res.push_back(v[i].real());
   }
   return std::move(res);
}
