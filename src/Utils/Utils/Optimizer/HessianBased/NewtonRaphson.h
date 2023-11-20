/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NEWTONRAPHSON_H_
#define UTILS_NEWTONRAPHSON_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Log.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of a Newton-Raphson optimization algorithm.
 */
class NewtonRaphson : public Optimizer {
 public:
  /// @brief Default constructor.
  NewtonRaphson() = default;
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
    while (!stop) {
      _cycle++;
      // Decompose Hessian, use SVD in order to be less strict on the Hessian's requirements
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(*hessian, Eigen::ComputeFullV | Eigen::ComputeFullU);
      svd.setThreshold(svdThreshold);
      Eigen::VectorXd steps = svd.solve(gradients);
      // Take a step
      const double maxVal = steps.array().abs().maxCoeff();
      if (maxVal > trustRadius) {
        steps *= trustRadius / maxVal;
      }
      constrainedAdd(parameters, -steps);
      // Update and check convergence
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
   *        thus expanding it to include the Newton-Raphson's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor nr_svd_threshold("The SVD threshold for the decomopsition of the Hessian.");
    nr_svd_threshold.setMinimum(0.0);
    nr_svd_threshold.setDefaultValue(svdThreshold);
    collection.push_back(SettingsNames::Optimizations::NewtonRaphson::svdThreshold, nr_svd_threshold);
    UniversalSettings::DoubleDescriptor nr_trust_radius("The maximum RMS of a taken step.");
    nr_trust_radius.setMinimum(0.0);
    nr_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(SettingsNames::Optimizations::NewtonRaphson::trustRadius, nr_trust_radius);
  };
  /**
   * @brief Updates the Newton-Raphson's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  void applySettings(const Settings& settings) final {
    svdThreshold = settings.getDouble(SettingsNames::Optimizations::NewtonRaphson::svdThreshold);
    trustRadius = settings.getDouble(SettingsNames::Optimizations::NewtonRaphson::trustRadius);
  };
  /// @brief The SVD threshold for the decomopsition of the Hessian.
  double svdThreshold = 1.0e-12;
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.5;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NEWTONRAPHSON_H_
