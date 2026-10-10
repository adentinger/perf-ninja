#include "solution.h"
#include <inttypes.h>
#include <random>

S create_entry(uint8_t first_value, uint8_t second_value) {
  S entry;

  entry.i = first_value;
  entry.s = static_cast<uint8_t>(second_value);
  entry.l = static_cast<unsigned short>(first_value * second_value);
  entry.d = static_cast<double>(first_value) / maxRandom;
  entry.b = first_value < second_value;

  return entry;
}

void init(std::vector<S> &arr) {
  std::default_random_engine generator;
  std::uniform_int_distribution<uint8_t> distribution(minRandom, maxRandom - 1);

  for (int i = 0; i < N; i++) {
    uint8_t random_int1 = distribution(generator);
    uint8_t random_int2 = distribution(generator);

    arr[i] = create_entry(random_int1, random_int2);
  }
}
