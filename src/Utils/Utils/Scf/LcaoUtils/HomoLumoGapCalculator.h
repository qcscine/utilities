/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LcaoUtils_HOMOLUMOGAPCALCULATOR_H
#define UTILS_LcaoUtils_HOMOLUMOGAPCALCULATOR_H

#include <exception>
#include <string>
#include <utility>

namespace Scine {
namespace Utils {

class SingleParticleEnergies;
namespace LcaoUtils {
class ElectronicOccupation;

class HomoLumoGapException : public std::exception {
 public:
  explicit HomoLumoGapException(std::string errorMessage) : message_(std::move(errorMessage)) {
  }
  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

/*!
 * This class calculates a Homo-Lumo gap, given the single-particle energies and the electronic occupation.
 * In the unrestricted case, puts returns energy difference between the orbitals irrespectively of whether they
 * are alpha- or beta-polarized.
 */
class HomoLumoGapCalculator {
 public:
  static double calculate(const SingleParticleEnergies& energies, const ElectronicOccupation& occupation);

 private:
  static double calculateRestricted(const SingleParticleEnergies& energies, const ElectronicOccupation& occupation);
  static double calculateUnrestricted(const SingleParticleEnergies& energies, const ElectronicOccupation& occupation);
};

} // namespace LcaoUtils

} // namespace Utils
} // namespace Scine
#endif // UTILS_LcaoUtils_HOMOLUMOGAPCALCULATOR_H