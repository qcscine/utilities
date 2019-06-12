/**
 * @file NumericalHessianCalculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NUMERICALHESSIANCALCULATOR_H
#define UTILS_NUMERICALHESSIANCALCULATOR_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

namespace Scine {
namespace Core {
class Calculator;
}
namespace Utils {

class State;
class Results;

/**
 * This class calculates the Hessian and, optionally, the dipole gradient (semi-)numerically.
 * Useful if analytic second derivatives are not available.
 *
 * calculateFromEnergyDifferences uses the following, where D = delta / 2:
 * d^2E/dx^2 = 1 / delta^2 * (E(x-2D) - 2 * E(x) + E(x+2D))
 * d^2E/dxdy = 1 / delta^2 * (E(x+D,y+D) - E(x+D,y-D) - E(x-D,y+D) + E(x-D,y-D))
 *
 * calculateFromGradientDifferences uses the following, where D = delta / 2:
 * d^2E/dxdy = 1 / (2*delta) * (g_j(i+D,j) - g_j(i-D,j) + g_i(i,j+D) - g_i(i,j-D))
 *
 * The second formulation is more stable numerically and is used as default.
 *
 * In order to calculate the dipole gradient, from each displacement the dipole is calculated.
 * dmu/dx_i = 1 / (2*delta) * (mu_x(x_i + delta) - mu_x(x_i - delta))
 * Only 3N calculations are needed to fill the 3Nx3 dipole gradient matrix,
 * as the dipole comes as a x,y,z vector at each single point calculation.
 */
class NumericalHessianCalculator {
 public:
  explicit NumericalHessianCalculator(Core::Calculator& calculator);

  /**
   * @brief Calcualtes the hessian matrix and, if needed, the dipole gradient.
   * @param delta The step size for the numerical derivative.
   * @return A Results class with an hessian matrix and, if needed, a dipole gradient.
   */
  Results calculate(double delta = defaultDelta);

  /**
   * This method can ONLY calculate the hessian matrix, not the dipole gradient.
   */
  HessianMatrix calculateFromEnergyDifferences(double delta);
  /**
   * This method calculated the hessian matrix from the gradient differences and the dipole
   * gradient as dipole difference for each displacement.
   */
  Results calculateFromGradientDifferences(double delta);

  /**
   * @brief Sets whether the dipole gradient is also calculated.
   */
  void requiredDipoleGradient(bool dipoleGradient);

 private:
  double hessianElementSameFromEnergy(int i, const PositionCollection& referencePositions, double delta);
  double hessianElementDifferentFromEnergy(int i, int j, const PositionCollection& referencePositions, double delta);

  Eigen::VectorXd addGradientContribution(DipoleGradient& dipoleDiff, int i,
                                          const Utils::PositionCollection& referencePositions, double delta,
                                          std::shared_ptr<Core::Calculator> calculator, std::shared_ptr<State> state);
  Core::Calculator& calculator_;

  // A step width of 0.01 bohr is also the default setting in MoViPac (see its manual for a study of this parameter)
  static constexpr double defaultDelta = 1e-2;
  bool dipoleGradient_{false};
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_NUMERICALHESSIANCALCULATOR_H
