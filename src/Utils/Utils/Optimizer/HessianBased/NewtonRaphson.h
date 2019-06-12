/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NEWTONRAPHSON_H_
#define UTILS_NEWTONRAPHSON_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of a Newton-Raphson optimization algorithm.
 */
class NewtonRaphson : public Optimizer {
 public:
  static constexpr const char* nrTrustRadius = "nr_trust_radius";
  static constexpr const char* nrSvdThreshold = "nr_svd_threshold";
  /// @brief Default constructor.
  NewtonRaphson() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
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
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, GradientBasedCheck& check) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    double value;

    /* Initialize gradients and Hessian. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    auto hessian = std::make_unique<Eigen::MatrixXd>(nParams, nParams);
    hessian->setZero();

    /* Initial update */
    function(parameters, value, gradients, *hessian);
    bool stop = false;
    int cycle = 0;
    while (!stop) {
      cycle++;
      // Decompose Hessian, use SVD in order to be less strict on the Hessian's requirements
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(*hessian, Eigen::ComputeFullV | Eigen::ComputeFullU);
      svd.setThreshold(svdThreshold);
      Eigen::VectorXd steps = svd.solve(gradients);
      // Take a step
      double rms = sqrt(steps.squaredNorm() / steps.size());
      if (rms > trustRadius) {
        steps *= trustRadius / rms;
      }
      parameters -= steps;
      // Update and check convergence
      function(parameters, value, gradients, *hessian);
      this->triggerObservers(cycle, value, parameters);
      stop = check.checkMaxIterations(cycle);
      if (!stop) {
        stop = check.checkConvergence(parameters, value, gradients);
      }
    }
    return cycle;
  };
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the Newton-Raphson's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor nr_svd_threshold("The SVD threshold for the decomopsition of the Hessian.");
    nr_svd_threshold.setDefaultValue(svdThreshold);
    collection.push_back(NewtonRaphson::nrSvdThreshold, nr_svd_threshold);
    UniversalSettings::DoubleDescriptor nr_trust_radius("The maximum RMS of a taken step.");
    nr_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(NewtonRaphson::nrTrustRadius, nr_trust_radius);
  };
  /**
   * @brief Updates the Newton-Raphson's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  virtual void applySettings(const Settings& settings) final {
    svdThreshold = settings.getDouble(NewtonRaphson::nrSvdThreshold);
    trustRadius = settings.getDouble(NewtonRaphson::nrTrustRadius);
  };
  /// @brief The SVD threshold for the decomopsition of the Hessian.
  double svdThreshold = 1.0e-12;
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.5;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NEWTONRAPHSON_H_
