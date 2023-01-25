/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Settings.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

bool GradientBasedCheck::checkConvergence(const Eigen::VectorXd& parameter, double value, const Eigen::VectorXd& gradient) {
  if (_oldParams.size() != parameter.size()) {
    _oldParams = Eigen::VectorXd::Zero(parameter.size());
  }
  // Generate temporaries
  Eigen::VectorXd deltaParam = (parameter - _oldParams).eval();
  double deltaV = value - _oldValue;
  // Rotate stored old data
  _oldParams = parameter;
  _oldValue = value;
  // Check
  unsigned int converged = 0;
  if (gradient.cwiseAbs().maxCoeff() < gradMaxCoeff) {
    converged++;
  }
  if (deltaParam.cwiseAbs().maxCoeff() < stepMaxCoeff) {
    converged++;
  }
  if (sqrt(gradient.squaredNorm() / gradient.size()) < gradRMS) {
    converged++;
  }
  if (sqrt(deltaParam.squaredNorm() / deltaParam.size()) < stepRMS) {
    converged++;
  }
  return ((fabs(deltaV) < deltaValue) && (converged >= requirement));
}

bool GradientBasedCheck::checkMaxIterations(unsigned int currentIteration) const {
  return currentIteration >= maxIter;
}

void GradientBasedCheck::setParametersAndValue(const Eigen::VectorXd& parameter, double value) {
  _oldParams = parameter;
  _oldValue = value;
}

void GradientBasedCheck::addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const {
  UniversalSettings::DoubleDescriptor step_max_coeff(
      "Convergence threshold for step vector's maximum absolute element.");
  step_max_coeff.setMinimum(0.0);
  step_max_coeff.setDefaultValue(stepMaxCoeff);
  collection.push_back(SettingsNames::Optimizations::Convergence::stepMaxCoeff, step_max_coeff);

  UniversalSettings::DoubleDescriptor step_RMS("Convergence threshold for step vector's RMS.");
  step_RMS.setMinimum(0.0);
  step_RMS.setDefaultValue(stepRMS);
  collection.push_back(SettingsNames::Optimizations::Convergence::stepRMS, step_RMS);

  UniversalSettings::DoubleDescriptor grad_max_coeff(
      "Convergence threshold for gradient vector's maximum absolute element.");
  grad_max_coeff.setMinimum(0.0);
  grad_max_coeff.setDefaultValue(gradMaxCoeff);
  collection.push_back(SettingsNames::Optimizations::Convergence::gradMaxCoeff, grad_max_coeff);

  UniversalSettings::DoubleDescriptor grad_RMS("Convergence threshold for gradient vector's RMS.");
  grad_RMS.setMinimum(0.0);
  grad_RMS.setDefaultValue(gradRMS);
  collection.push_back(SettingsNames::Optimizations::Convergence::gradRMS, grad_RMS);

  UniversalSettings::DoubleDescriptor delta_value(
      "Convergence threshold for the absolute difference in the value between the current and the last step.");
  delta_value.setMinimum(0.0);
  delta_value.setDefaultValue(deltaValue);
  collection.push_back(SettingsNames::Optimizations::Convergence::deltaValue, delta_value);

  UniversalSettings::IntDescriptor max_iter("The maximum number of iterations.");
  max_iter.setMinimum(0.0);
  max_iter.setDefaultValue(maxIter);
  collection.push_back(SettingsNames::Optimizations::Convergence::maxIter, max_iter);

  UniversalSettings::IntDescriptor requirements("The number of thresholds besides the value one that need to "
                                                "Convergence::erge for overall Convergence::ergence.");
  requirements.setDefaultValue(requirement);
  requirements.setMaximum(4);
  requirements.setMinimum(0);
  collection.push_back(SettingsNames::Optimizations::Convergence::requirement, requirements);
}

void GradientBasedCheck::applySettings(const Settings& settings) {
  stepMaxCoeff = settings.getDouble(SettingsNames::Optimizations::Convergence::stepMaxCoeff);
  stepRMS = settings.getDouble(SettingsNames::Optimizations::Convergence::stepRMS);
  gradMaxCoeff = settings.getDouble(SettingsNames::Optimizations::Convergence::gradMaxCoeff);
  gradRMS = settings.getDouble(SettingsNames::Optimizations::Convergence::gradRMS);
  deltaValue = settings.getDouble(SettingsNames::Optimizations::Convergence::deltaValue);
  maxIter = settings.getInt(SettingsNames::Optimizations::Convergence::maxIter);
  requirement = settings.getInt(SettingsNames::Optimizations::Convergence::requirement);
}

} // namespace Utils
} // namespace Scine
