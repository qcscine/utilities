/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Solvation/RandomIndexGenerator.h"

namespace Scine {
namespace Utils {

SoluteSolventComplex::RandomIndexGenerator::RandomIndexGenerator(int size, int seed) {
  generator_ = std::mt19937(seed);
  distribution_ = std::uniform_int_distribution<>(0, size - 1);
}

void SoluteSolventComplex::RandomIndexGenerator::setSeed(int seed) {
  generator_ = std::mt19937(seed);
}

void SoluteSolventComplex::RandomIndexGenerator::setSize(int size) {
  distribution_ = std::uniform_int_distribution<>(0, size - 1);
}

int SoluteSolventComplex::RandomIndexGenerator::next() {
  return distribution_(generator_);
}

} /* namespace Utils */
} /* namespace Scine */
