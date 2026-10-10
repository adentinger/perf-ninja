#include <vector>
#include <inttypes.h>

// Assume those constants never change
constexpr int N = 1000000;
constexpr uint8_t minRandom = 0;
constexpr uint8_t maxRandom = 100;
static_assert(
  maxRandom <= 127,
  "struct S bitfields require maxRandom <= 127!");

// FIXME: this data structure can be reduced in size
struct S {
  uint8_t i;
  uint8_t s : 7;
  bool b : 1;
  unsigned short l;
  float d;

  bool operator<(const S &s) const { return this->i < s.i; }
};

void init(std::vector<S> &arr);
S create_entry(uint8_t first_value, uint8_t second_value);
void solution(std::vector<S> &arr);
