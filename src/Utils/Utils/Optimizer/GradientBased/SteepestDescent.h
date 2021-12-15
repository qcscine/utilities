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
  static constexpr const char* sdFactorKey = "sd_factor";

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
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check) {
    if (parameters.size() == 0) {
      throw EmptyOptimizerParametersException();
    }
    double value = 0.0;
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(parameters.size());
    int cycle = _startCycle;
    function(parameters, value, gradients);
    this->triggerObservers(cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    bool stop = false;
    while (!stop) {
      cycle++;
      constrainedAdd(parameters, -factor * gradients);
      function(parameters, value, gradients);
      this->triggerObservers(cycle, value, parameters);
      stop = check.checkMaxIterations(cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients));

      // Check oscillation, perform new calculation if oscillating and lower factor
      if (!stop && this->isOscillating(value)) {
        oscillationCorrection(-factor * gradients, parameters);
        function(parameters, value, gradients);
        cycle++;
        this->triggerObservers(cycle, value, parameters);
        factor *= 0.95;
        _oscillationCounter++;
      }
      else {
        // reset factor
        factor /= std::pow(0.95, _oscillationCounter);
        _oscillationCounter = 0;
      }
    }
    return cycle;
  }
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the steepest descent options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor sd_factor("The steepest descent scaling factor.");
    sd_factor.setDefaultValue(factor);
    collection.push_back(SteepestDescent::sdFactorKey, sd_factor);
  };
  /**
   * @brief Updates the steepest descent's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  void applySettings(const Settings& settings) final {
    factor = settings.getDouble(SteepestDescent::sdFactorKey);
  };
  /**
   * @brief The scaling factor alpha to be used in the steepest descent algorithm.
   *
   * The parameters \f$\{x_i\}\f$ are generated as:
   * \f[ x_{i,n+1} = x_i - \alpha * g_i \f]
   */
  double factor = 0.1;

 private:
  // The number of back-to-back oscillation detections and factor decreases
  unsigned int _oscillationCounter = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STEEPESTDESCENT_H_
