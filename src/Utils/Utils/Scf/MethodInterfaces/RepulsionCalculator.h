/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_REPULSIONCALCULATOR_H
#define UTILS_REPULSIONCALCULATOR_H

#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/First1D.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <map>

namespace Scine {
namespace Utils {

/**
 * @class RepulsionCalculator @file RepulsionCalculator.h
 * @brief Class representing a repulsion calculator between (classical) nuclei.
 */
class RepulsionCalculator {
 public:
  RepulsionCalculator(const ElementTypeCollection& elements, const PositionCollection& positions);
  virtual ~RepulsionCalculator() = default;

  /**
   * @brief Initializes a vector containing the classes computing the pairwise interactions between cores.
   */
  virtual void initialize();

  /**
   * @brief Method to calculate the pairwise core-core repulsion.
   * This functions dispatches the call to the right PairwiseCoreRepulsion class
   * computing the interaction between two cores up to the derivative order specified.
   */
  virtual void calculateRepulsion(DerivativeOrder order);

  /**
   * @brief Method to collect the energy contributions stored in the PairwiseCoreRepulsion classes.
   */
  virtual double getRepulsionEnergy() const;
  /**
   * @brief Method to collect the energy derivatives contributions stored in the PairwiseCoreRepulsion classes.
   * @{
   */
  virtual void addRepulsionDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::First>& derivatives) const;
  virtual void
  addRepulsionDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondAtomic>& derivatives) const;
  virtual void addRepulsionDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondFull>& derivatives) const;
  //! @}

 protected:
  /**
   * @brief Use AutomaticDifferentiation to get the right values for the energy derivatives.
   */
  template<DerivativeOrder order>
  AutomaticDifferentiation::Value1DType<order> calculatePairwiseCoreRepulsion(double distance, double repulsionConstant) {
    auto R = Utils::AutomaticDifferentiation::variableWithUnitDerivative<order>(distance);
    // The result of this depends on the derivative order!
    auto inverseR = 1.0 / R;
    return inverseR * repulsionConstant;
  }

  const ElementTypeCollection& elements_;
  const PositionCollection& positions_;
  double repulsionEnergy_;
  GradientCollection repulsionGradients_;
  std::map<std::pair<int, int>, AutomaticDifferentiation::Second3D> repulsionHessian_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_REPULSIONCALCULATOR_H
