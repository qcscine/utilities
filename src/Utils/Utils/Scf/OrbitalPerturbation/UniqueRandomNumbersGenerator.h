/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SEMO_UNIQUERANDOMNUMBERGENERATOR_H
#define SEMO_UNIQUERANDOMNUMBERGENERATOR_H

#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>
#include <vector>

namespace Scine {
namespace Utils {

/**
 * @brief Class to sample without replacement N numbers from a range bounded by min_ and max_
 *
 * The user gives a minimum number and a maximum number as the boundaries of a range of possible numbers.
 * Upon calling the UniqueRandomNumbersGenerator::generate(unsigned N) function, N values are
 * sampled without replacement from the possible values.
 */

template<class IntegerType>
class UniqueRandomNumbersGenerator {
 public:
  void setMin(IntegerType min) noexcept;
  void setMax(IntegerType max) noexcept;
  void setRange(IntegerType min, IntegerType max) noexcept;

  std::vector<IntegerType> generate(unsigned N);

  template<class Generator>
  std::vector<IntegerType> generate(Generator& generator, unsigned N);

  IntegerType min_ = 0;
  IntegerType max_ = 0;
};

template<class IntegerType>
void UniqueRandomNumbersGenerator<IntegerType>::setMin(IntegerType min) noexcept {
  min_ = min;
}

template<class IntegerType>
void UniqueRandomNumbersGenerator<IntegerType>::setRange(IntegerType min, IntegerType max) noexcept {
  min_ = min;
  max_ = max;
}

template<class IntegerType>
void UniqueRandomNumbersGenerator<IntegerType>::setMax(IntegerType max) noexcept {
  max_ = max;
}

template<class IntegerType>
std::vector<IntegerType> UniqueRandomNumbersGenerator<IntegerType>::generate(unsigned N) {
  assert(min_ + static_cast<int>(N) < max_ + 2 && "Not enough numbers to select from.");

  std::vector<IntegerType> shuffler(max_ - min_ + 1);
  std::iota(shuffler.begin(), shuffler.end(), min_);
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(shuffler.begin(), shuffler.end(), g);
  return {shuffler.begin(), shuffler.begin() + N};
}

template<class IntegerType>
template<class Generator>
std::vector<IntegerType> UniqueRandomNumbersGenerator<IntegerType>::generate(Generator& generator, unsigned N) {
  assert(min_ + static_cast<int>(N) < max_ + 2 && "Not enough numbers to select from.");

  std::vector<IntegerType> shuffler(max_ - min_ + 1);
  std::iota(shuffler.begin(), shuffler.end(), min_);
  std::shuffle(shuffler.begin(), shuffler.end(), generator);
  return {shuffler.begin(), shuffler.begin() + N};
}

} // namespace Utils
} // namespace Scine

#endif // SEMO_UNIQUERANDOMNUMBERGENERATOR_H
