/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_STEEPESTDESCENT_H_
#define UTILS_STEEPESTDESCENT_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Core/Log.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of a steepest descent optimization algorithm.
 *
 * The steepest descent algorithm is on of the most simple gradient based optimization algorithms.
 * In each step the parameters are adjusted by subtracting the scaled negative gradient.
 * The gradient is usually scaled by a factor that is smaller than 1.0 .
 */
class SteepestDescent : public Optimizer {
 public:
  /// @brief Default constructor.
  SteepestDescent() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
   * @tparam ConvergenceCheckClass The convergence check class, needs to have a check function signature
   *                               that is the same as the GradientBasedCheck.
   *
   * @param parameters The parameters to be optimized.
   * @param function   The function to be evaluated in order to get values and gradients
   *                   for a given set of parameters.
   * @param check      The ConvergenceCheck to be used in order to determine when the optimization
   *                   is finished or should stop for other reasons.
   * @return int       Returns the number of optimization cycles carried out until the conclusion
   *                   of the optimization function.
   */
  template<class UpdateFunction, class ConvergenceCheckClass>
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check,
               Core::Log& /* log */) {
    if (parameters.size() == 0) {
      throw EmptyOptimizerParametersException();
    }
    double value = 0.0;
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(parameters.size());
    /* Set "assumed previous parameters" as parametersOld. */
    Eigen::VectorXd parametersOld(parameters);
    _cycle = _startCycle;

    function(parameters, value, gradients);
    this->triggerObservers(_cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    bool stop = false;
    bool previousStepScaling = true;
    while (!stop) {
      _cycle++;
      /* Multiply factor with dynamic multiplier (1.0 per default) only if no scaling has been performed previously */
      if (!previousStepScaling) {
        factor *= dynamicMultiplier;
      }
      /* Calculate step vector with given gradients and factor */
      Eigen::VectorXd stepVector = -factor * gradients;
      /* Apply trust radius check and scaling of step if necessary */
      if (useTrustRadius) {
        double maxVal = stepVector.array().abs().maxCoeff();
        /* Scale step vector if max value exceeds the given trust radius */
        if (maxVal > trustRadius) {
          /* Use the minimum of the initial factor and the current factor scaled such as
           * the maximum step corresponds to the given trust radius. */
          stepVector.noalias() = std::min(_initFactor, factor * trustRadius / maxVal) * -gradients;
          previousStepScaling = true;
          /* Reset factor to initial factor*/
          factor = _initFactor;
        }
        else {
          previousStepScaling = false;
        }
      }
      /* Apply step on parameters */
      constrainedAdd(parameters, stepVector);
      // Calculate things of next step
      function(parameters, value, gradients);
      this->triggerObservers(_cycle, value, parameters);
      stop = check.checkMaxIterations(_cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients));

      // Check oscillation, perform new calculation if oscillating and lower factor
      if (!stop && this->isOscillating(value)) {
        oscillationCorrection(-factor * gradients, parameters);
        function(parameters, value, gradients);
        _cycle++;
        this->triggerObservers(_cycle, value, parameters);
        factor *= 0.95;
        _oscillationCounter++;
      }
      else {
        // reset factor
        factor /= std::pow(0.95, _oscillationCounter);
        _oscillationCounter = 0;
      }
    }
    return _cycle;
  }
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the steepest descent options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor sd_factor("The steepest descent scaling factor.");
    sd_factor.setDefaultValue(factor);
    collection.push_back(SettingsNames::Optimizations::SteepestDescent::factor, sd_factor);
    UniversalSettings::BoolDescriptor sd_useTrustRadius("Enable the use of a trust radius for all steps.");
    sd_useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(SettingsNames::Optimizations::SteepestDescent::useTrustRadius, sd_useTrustRadius);
    UniversalSettings::DoubleDescriptor sd_trustRadius("The maximum size (RMS) of a taken step.");
    sd_trustRadius.setMinimum(0.0);
    sd_trustRadius.setDefaultValue(trustRadius);
    collection.push_back(SettingsNames::Optimizations::SteepestDescent::trustRadius, sd_trustRadius);
    UniversalSettings::DoubleDescriptor sd_dynamicMultiplier(
        "The multiplier to increase the scaling factor after each iteration.");
    sd_dynamicMultiplier.setMinimum(1.0);
    sd_dynamicMultiplier.setDefaultValue(dynamicMultiplier);
    collection.push_back(SettingsNames::Optimizations::SteepestDescent::dynamicMultiplier, sd_dynamicMultiplier);
  };
  /**
   * @brief Updates the steepest descent's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  void applySettings(const Settings& settings) final {
    factor = settings.getDouble(SettingsNames::Optimizations::SteepestDescent::factor);
    _initFactor = factor;

    useTrustRadius = settings.getBool(SettingsNames::Optimizations::SteepestDescent::useTrustRadius);
    trustRadius = settings.getDouble(SettingsNames::Optimizations::SteepestDescent::trustRadius);
    if (!useTrustRadius && std::fabs(trustRadius - 0.1) > 1e-6) {
      throw std::logic_error("A trust radius was specified, but the trust radius was not activated. "
                             "Please also set the setting 'sd_use_trust_radius': true, if you specify a radius.");
    }

    dynamicMultiplier = settings.getDouble(SettingsNames::Optimizations::SteepestDescent::dynamicMultiplier);
    if (!useTrustRadius && std::fabs(dynamicMultiplier - 1.0) > 1e-6) {
      throw std::logic_error("A dynamic multiplier was specified, but the trust radius was not activated. "
                             "Please also set the setting 'sd_use_trust_radius': true, if you specify a multiplier.");
    }
  };
  /**
   * @brief The scaling factor alpha to be used in the steepest descent algorithm.
   *
   * The parameters \f$\{x_i\}\f$ are generated as:
   * \f[ x_{i,n+1} = x_i - \alpha * g_i \f]
   */
  double factor = 0.1;
  /// @brief Enable the use of a trust radius for all steps.
  bool useTrustRadius = false;
  /// @brief The maximum size (RMS) of a taken step.
  double trustRadius = 0.1;

  double dynamicMultiplier = 1.0;

 private:
  // The number of back-to-back oscillation detections and factor decreases
  unsigned int _oscillationCounter = 0;
  // @brief The initial scaling factor alpha stored for resetting the factor when necessary.
  double _initFactor = factor;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STEEPESTDESCENT_H_
