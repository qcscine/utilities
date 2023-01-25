/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GRADIENTBASEDCHECK_H_
#define UTILS_GRADIENTBASEDCHECK_H_

#include "Utils/Optimizer/ConvergenceCheck.h"
#include "Utils/Settings.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {
/**
 * @brief A convergence check based on parameter, value and gradient information.
 *
 * Checks the RMS of the gradient and step, the maximum coefficient of the gradient
 * and step and also the change in value.
 *
 * For convergence the value and a set number of the other four criteria need to
 * converge.
 */
class GradientBasedCheck : public ConvergenceCheck {
 public:
  /// @brief Default constructor.
  GradientBasedCheck() = default;
  /// @brief Default destructor.
  virtual ~GradientBasedCheck() = default;

  /**
   * @brief Checks for convergence.
   *
   * @param parameter The current parameters.
   * @param value     The current value
   * @param gradient  The current gradient.
   * @return true     If converged.
   * @return false    If not converged.
   */
  virtual bool checkConvergence(const Eigen::VectorXd& parameter, double value, const Eigen::VectorXd& gradient);
  /**
   * @brief Checks if the current iteration is within the maximum of allowed iterations.
   *
   * @param currentIteration The number of the current iteration.
   * @return true            If the current iteration is smaller than the maximum allowed number of iterations.
   * @return false           If the current iteration is equal to or greater than the maximum allowed number of
   *                         iterations.
   */
  bool checkMaxIterations(unsigned int currentIteration) const;

  /**
   * @brief Sets the old parameters and value stored in the convergence checker.
   *
   * To be used to modify/init the parameters and value without performing a convergence check.
   *
   * @param parameter The current parameters.
   * @param value     The current value
   */
  void setParametersAndValue(const Eigen::VectorXd& parameter, double value);

  /// @brief See Scine::Utils::ConvergenceCheck::addSettingsDescriptors()
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final;
  /// @brief See Scine::Utils::ConvergenceCheck::applySettings()
  void applySettings(const Settings& s) final;

  /// @brief The threshold for the maximum absolute element of the last step taken.
  double stepMaxCoeff = 1.0e-4;
  /// @brief The threshold for the root mean square of the last step taken.
  double stepRMS = 5.0e-4;
  /// @brief The threshold for the maximum absolute element of the gradient.
  double gradMaxCoeff = 5.0e-5;
  /// @brief The threshold for the root mean square of the gradient.
  double gradRMS = 1.0e-5;
  /// @brief The threshold for the change in the functional value.
  double deltaValue = 1.0e-7;
  /// @brief The maximum number of iterations.
  unsigned int maxIter = 150;
  /// @brief The number of criteria that have to converge besides the value criterion.
  unsigned int requirement = 3;

 private:
  Eigen::VectorXd _oldParams;
  double _oldValue = 0.0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_GRADIENTBASEDCHECK_H_
