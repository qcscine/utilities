/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONICENERGYCALCULATOR_H
#define UTILS_ELECTRONICENERGYCALCULATOR_H

namespace Scine {
namespace Utils {

/*!
 * Interface for the calculation of the electronic energy.
 */
class ElectronicEnergyCalculator {
 public:
  virtual ~ElectronicEnergyCalculator() = default;

  virtual double calculateElectronicEnergy() = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_ELECTRONICENERGYCALCULATOR_H
