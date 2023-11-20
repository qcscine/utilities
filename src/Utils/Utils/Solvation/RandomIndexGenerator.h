/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_RANDOMINDEXGENERATOR_H
#define UTILSOS_RANDOMINDEXGENERATOR_H

#include <random>

namespace Scine {
namespace Utils {
namespace SoluteSolventComplex {

/**
 * @class RandomIndexGenerator RandomIndexGenerator.h
 * @brief A random number generator for indices between 0 and size - 1
 *
 * Upon calling next it generators number with a uniform distribution between 0 and size - 1
 */
class RandomIndexGenerator {
 public:
  /**
   * @brief Default empty constructor
   */
  RandomIndexGenerator() = default;
  /**
   * @brief Initialising the RandomIndexGenerator, assigning size and seed to its generator and distribution
   * @param size Size of the wanted distribution, meaning the size of the related vector.
   * @param seed Seed for the generator.
   */
  RandomIndexGenerator(int size, int seed);

  /**
   * @brief Set seed for generator
   * @param seed Seed for generator.
   */
  void setSeed(int seed);
  /**
   * @brief Set size of wanted distribution.
   * @param size Size of the wanted distribution.
   */
  void setSize(int size);
  /**
   * @brief Generate next random index from initialised RandomIndexGenerator.
   * @return Integer between 0 and size - 1 with uniform distribution.
   */
  int next();

 private:
  std::mt19937 generator_;
  std::uniform_int_distribution<> distribution_;
};

} // namespace SoluteSolventComplex
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_RANDOMINDEXGENERATOR_H
