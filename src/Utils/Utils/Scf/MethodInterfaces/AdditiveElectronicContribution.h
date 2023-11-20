/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_ADDITIVEELECTRONICCONTRIBUTION_H
#define UTILSOS_ADDITIVEELECTRONICCONTRIBUTION_H

#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>

namespace Scine {
namespace Utils {

class SpinAdaptedMatrix;
class DensityMatrix;

/**
 * @class AdditiveElectronicContribution @file AdditiveElectronicContribution.h
 * This class represents an additive electronic contribution.
 * On one sides it stores the matrix with the matrix elements, on the other side
 * it is an interface for the addDerivatives<Utils::Derivative O> methods in the derived classes.
 */
class AdditiveElectronicContribution {
 public:
  /** Virtual destructor. */
  virtual ~AdditiveElectronicContribution() = default;

  /**
   * @brief Function to identify the density dependence of the electronic contribution.
   */
  virtual bool isDensityDependent() const = 0;

  /**
   * @brief Function to check if a contribution is already present. Does not check for correctness.
   */
  virtual bool isValid() const {
    return isValid_;
  };

  /**
   * @brief Function to check if a contribution in form of matrix is present.
   */
  virtual bool hasMatrixContribution() const {
    return hasMatrixContribution_;
  }
  /**
   * @brief Invalidates the contribution, which has to be set again.
   */
  void invalidate() {
    isValid_ = false;
  }
  /**
   * @brief Setter function.
   * @param contribution The matrix containing the electronic contributions.
   */
  void setElectronicContribution(const SpinAdaptedMatrix& contribution) {
    electronicContribution_ = contribution;
    isValid_ = true;
  }

  /**
   * @brief Getter for the electronic contributions.
   */
  const SpinAdaptedMatrix& getElectronicContribution() const {
    return electronicContribution_;
  }

  /**
   * @brief Returns the energy contribution. It needs the density matrix for the contraction of the matrix elements.
   */
  virtual double getElectronicEnergyContribution() const = 0;
  /**
   * @brief Calculates the energy contribution.
   * For example, calculates the dipole-field interaction energy.
   */
  virtual void calculate(const DensityMatrix& densityMatrix, DerivativeOrder order) = 0;
  /**
   * @brief Functions to calculate the first derivatives of the electronic contribution.
   * @param derivativeContainer A container for the derivatives up to the first order.
   */
  virtual void addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::First>& derivativeContainer) const = 0;
  /**
   * @brief Functions to calculate the atomic second derivatives of the electronic contribution.
   * @param derivativeContainer A container for the atomic derivatives up to the second order.
   */
  virtual void
  addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondAtomic>& derivativeContainer) const = 0;
  /**
   * @brief Functions to calculate the second derivatives of the electronic contribution.
   * @param derivativeContainer A container for the derivatives up to the second order.
   */
  virtual void
  addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondFull>& derivativeContainer) const = 0;

 private:
  SpinAdaptedMatrix electronicContribution_;
  bool isValid_{false};
  bool hasMatrixContribution_{false};
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_ADDITIVEELECTRONICCONTRIBUTION_H
