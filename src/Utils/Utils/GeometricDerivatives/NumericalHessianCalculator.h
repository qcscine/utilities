/**
 * @file NumericalHessianCalculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
class State;
} // namespace Core
namespace Utils {

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
   * @brief Calculates the hessian matrix and, if needed, the dipole gradient.
   * @param delta The step size for the numerical derivative.
   * @return A Results class with an hessian matrix and, if needed, a dipole gradient.
   */
  Results calculate(double delta = defaultDelta);
  /**
   * @brief Calculates the hessian matrix and, if needed, the dipole gradient from a subset of atoms.
   * @param indices A vector containing the indices of the atom from which to calculate the Hessian and dipole Gradient.
   * @param delta The step size for the numerical derivative.
   * @return A Results class with an hessian matrix and, if needed, a dipole gradient of a subset of atoms.
   */
  Results calculate(const std::vector<int>& indices, double delta = defaultDelta);

  /**
   * This method can ONLY calculate the hessian matrix, not the dipole gradient.
   */
  HessianMatrix calculateFromEnergyDifferences(double delta);
  /**
   * This method calculates the hessian matrix from the gradient differences and the dipole
   * gradient as dipole difference for each displacement.
   */
  Results calculateFromGradientDifferences(double delta);

  /**
   * This method modifies ONLY the columns/rows of the hessian matrix corresponding to the indices
   * given as argument from the gradient differences and the dipole
   * gradient as dipole difference for each displacement.
   */
  Results calculateFromGradientDifferences(double delta, const std::vector<int>& atomsToUpdate);
  /**
   * @brief Sets whether the dipole gradient is also calculated.
   */
  void requiredDipoleGradient(bool dipoleGradient);

 private:
  double hessianElementSameFromEnergy(int i, const PositionCollection& referencePositions, double delta);
  double hessianElementDifferentFromEnergy(int i, int j, const PositionCollection& referencePositions, double delta);

  Eigen::VectorXd addGradientContribution(DipoleGradient& dipoleDiff, int i, const Utils::PositionCollection& referencePositions,
                                          double delta, Core::Calculator& calculator, std::shared_ptr<Core::State> state);
  Core::Calculator& calculator_;

  // A step width of 0.01 bohr is also the default setting in MoViPac (see its manual for a study of this parameter)
  static constexpr double defaultDelta = 1e-2;
  bool dipoleGradient_{false};
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_NUMERICALHESSIANCALCULATOR_H
