#include "solution.hpp"
#include <numeric>

uint16_t checksum(const Blob &blob) {
  // VECTORIZE_WIDTH times the element size should probably be equal to the
  // SIMD register width for the current achitecture.
  constexpr uint16_t VECTORIZE_WIDTH = 256;
  uint16_t acc[VECTORIZE_WIDTH] = {};
  static_assert((N % VECTORIZE_WIDTH) == 0, "VECTORIZE_WIDTH must divide N.");
  constexpr std::size_t NUM_WINDOWS = N / VECTORIZE_WIDTH;

  for (std::size_t window_idx = 0; window_idx < NUM_WINDOWS; ++window_idx) {
    for (uint16_t vec_idx = 0; vec_idx < VECTORIZE_WIDTH; ++vec_idx) {
      const auto value = blob[window_idx * VECTORIZE_WIDTH + vec_idx];
      acc[vec_idx] += value;
      // add carry
      acc[vec_idx] += acc[vec_idx] < value;
    }
  }
  // Accumulate the accumulator vector itself!
  uint16_t result = 0;
  for (uint16_t vec_idx = 0; vec_idx < VECTORIZE_WIDTH; ++vec_idx) {
    result += acc[vec_idx];
    // add carry
    result += result < acc[vec_idx];
  }

  return result;
}
