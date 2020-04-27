/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONICCONTRIBUTIONCALCULATOR_H
#define UTILS_ELECTRONICCONTRIBUTIONCALCULATOR_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <memory>

namespace Scine {

namespace Utils {

enum class derivOrder;
class SpinAdaptedMatrix;
class AdditiveElectronicContribution;

/**
 * Interface for the computation of the Fock matrix.
 * TODO: Separate for non-SCF and SCF specializations?
 */
class ElectronicContributionCalculator {
 public:
  virtual ~ElectronicContributionCalculator() = default;

  /*! Reinitialize after a change in the elements. */
  virtual void initialize() = 0;

  /** This function will be called only once per single-point calculation. */
  virtual void calculateDensityIndependentPart(Utils::derivOrder order) = 0;
  /** This function will be called once per SCF iteration. */
  virtual void calculateDensityDependentPart(Utils::derivOrder order) = 0;
  /** This function will be called after the last iteration has run. */
  virtual void finalize(Utils::derivOrder order) = 0;
  /**
   * This function adds an additive electronic contribution to the
   * Hamiltonian that will be evaluated each SCF iteration.
   */
  virtual void addDensityDependentElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> contribution) = 0;
  /**
   * This function adds an additive electronic contribution to the Hamiltonian
   * that will be evaluated once per single-point calculation.
   */
  virtual void addDensityIndependentElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> contribution) = 0;

  virtual SpinAdaptedMatrix getMatrix() const = 0;

  virtual double calculateElectronicEnergy() const = 0;
  virtual void
  addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const = 0;
  virtual void addDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const = 0;
  virtual void addDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_ELECTRONICCONTRIBUTIONCALCULATOR_H
