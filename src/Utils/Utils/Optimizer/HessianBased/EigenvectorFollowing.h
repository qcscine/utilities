/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EIGENVECTORFOLLOWING_H_
#define UTILS_EIGENVECTORFOLLOWING_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Core/Log.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of an Eigenvector following optimization algorithm.
 *
 * This algorithm is intended to find the maximum along one single eigenvector
 * and the minimum along all other eigenvectors of a given system/Hessian.
 */
class EigenvectorFollowing : public Optimizer {
 public:
  static constexpr const char* evTrustRadius = "ev_trust_radius";
  static constexpr const char* evMode = "ev_follow_mode";
  /// @brief Default constructor.
  EigenvectorFollowing() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
   *                        4. Eigen::Matrixd& the Hessian
   *                        5. bool a flag if the Hessian is to be calculated
   *
   * @param parameters The parameters to be optimized.
   * @param function   The function to be evaluated in order to get values and gradients
   *                   for a given set of parameters.
   * @param check      The ConvergenceCheck to be used in order to determine when the optimization
   *                   is finished or should stop for other reasons.
   * @return int       Returns the number of optimization cycles carried out until the conclusion
   *                   of the optimization function.
   */
  template<class UpdateFunction>
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, GradientBasedCheck& check, Core::Log& /* log */) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    if (nParams == 0) {
      throw EmptyOptimizerParametersException();
    }
    double value;
    /* Initialize gradients and Hessian. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    auto hessian = std::make_unique<Eigen::MatrixXd>(nParams, nParams);
    hessian->setZero();
    _cycle = _startCycle;
    /* Initial update */
    function(parameters, value, gradients, *hessian, true);
    this->triggerObservers(_cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    bool stop = false;
    Eigen::VectorXd prevFollowedEigenvector;
    prevFollowedEigenvector.resize(nParams);
    while (!stop) {
      _cycle++;
      Eigen::VectorXd steps = modeMaximizationWithHessian(gradients, *hessian, modeToFollow);
      // Check trustradius and take a step
      const double maxVal = steps.array().abs().maxCoeff();
      if (maxVal > trustRadius) {
        steps *= trustRadius / maxVal;
      }
      // Update and check convergence
      constrainedAdd(parameters, steps);
      function(parameters, value, gradients, *hessian, true);
      this->triggerObservers(_cycle, value, parameters);
      stop = check.checkMaxIterations(_cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients));
      // Check oscillation, if oscillating modify step and calculate new Hessian
      if (!stop && this->isOscillating(value)) {
        oscillationCorrection(steps, parameters);
        function(parameters, value, gradients, *hessian, true);
        _cycle++;
        this->triggerObservers(_cycle, value, parameters);
      }
    }
    return _cycle;
  }
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the optimizers's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor ev_trust_radius("The maximum RMS of a taken step.");
    ev_trust_radius.setMinimum(0.0);
    ev_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(EigenvectorFollowing::evTrustRadius, ev_trust_radius);
    UniversalSettings::IntDescriptor ev_follow_mode("The number of the mode that should be followed starting with 0.");
    ev_follow_mode.setMinimum(0);
    ev_follow_mode.setDefaultValue(modeToFollow);
    collection.push_back(EigenvectorFollowing::evMode, ev_follow_mode);
  };
  /**
   * @brief Updates the optimizers's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  void applySettings(const Settings& settings) final {
    trustRadius = settings.getDouble(EigenvectorFollowing::evTrustRadius);
    modeToFollow = settings.getInt(EigenvectorFollowing::evMode);
  };
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 1.0e-1;
  /// @brief The number of the mode to follow
  int modeToFollow = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_EIGENVECTORFOLLOWING_H_
