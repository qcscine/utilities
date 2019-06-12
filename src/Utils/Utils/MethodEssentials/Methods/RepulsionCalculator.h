/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_REPULSIONCALCULATOR_H
#define UTILS_REPULSIONCALCULATOR_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/DerivOrderEnum.h>

namespace Scine {
namespace Utils {

/*!
 * Interface for the calculation of the repulsion energy of LCAO methods.
 */
class RepulsionCalculator {
 public:
  virtual ~RepulsionCalculator() = default;

  /*! Reinitialize after a change in the elements. */
  virtual void initialize() = 0;

  /*! Calculate the repulsion. */
  virtual void calculateRepulsion(Utils::derivOrder order) = 0;

  /*! Get the energy and derivatives. */
  virtual double getRepulsionEnergy() const = 0;
  virtual void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const = 0;
  virtual void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const = 0;
  virtual void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const = 0;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_REPULSIONCALCULATOR_H
