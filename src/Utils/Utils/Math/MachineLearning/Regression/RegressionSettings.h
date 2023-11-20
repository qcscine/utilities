/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_MLREGRESSIONSETTINGS_H
#define UTILSOS_MLREGRESSIONSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/OptimizationSettingsNames.h>

namespace Scine {
namespace Utils {
namespace MachineLearning {

/**
 * @brief Settings class for the iterative diagonalizers.
 */
class RegressionSettings : public Scine::Utils::Settings {
 public:
  void setRestartOptimization(UniversalSettings::DescriptorCollection& settings);
  void setNumRestarts(UniversalSettings::DescriptorCollection& settings);
  void setMaxIterations(UniversalSettings::DescriptorCollection& settings);
  void setMaxLinesearch(UniversalSettings::DescriptorCollection& settings);
  void setConvergenceTolerance(UniversalSettings::DescriptorCollection& settings);
  void setLinesearchTolerance(UniversalSettings::DescriptorCollection& settings);

  RegressionSettings() : Settings("RegressionSettings") {
    setRestartOptimization(_fields);
    setNumRestarts(_fields);
    setMaxIterations(_fields);
    setMaxLinesearch(_fields);
    setConvergenceTolerance(_fields);
    setLinesearchTolerance(_fields);
    resetToDefaults();
  }
};

inline void RegressionSettings::setRestartOptimization(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor restartOptimization(
      "Whether to restart the hyperparameter optimization from different starting points.");
  restartOptimization.setDefaultValue(true);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::restartOptimization, std::move(restartOptimization));
}

inline void RegressionSettings::setNumRestarts(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numRestarts("The number of restarts for hyperparameter optimization.");
  numRestarts.setMinimum(1);
  numRestarts.setDefaultValue(10);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::numRestarts, std::move(numRestarts));
}

inline void RegressionSettings::setMaxIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxIterations(
      "The maximum number of iterations. Note that setting this parameter to zero continues an optimization process "
      "until convergence or error.");
  maxIterations.setMinimum(1);
  maxIterations.setDefaultValue(1000);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::maxIterations, std::move(maxIterations));
}

inline void RegressionSettings::setMaxLinesearch(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxLinesearch("The maximum number of trials for the line search.");
  maxLinesearch.setMinimum(1);
  maxLinesearch.setDefaultValue(20000);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::maxLinesearch, std::move(maxLinesearch));
}

inline void RegressionSettings::setConvergenceTolerance(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor convergenceTolerance("The absolute tolerance for convergence test.");
  convergenceTolerance.setDefaultValue(1e-6);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::convergenceTolerance, std::move(convergenceTolerance));
}

inline void RegressionSettings::setLinesearchTolerance(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor ftol("A parameter to control the accuracy of the line search routine.");
  ftol.setDefaultValue(0.001);
  settings.push_back(SettingsNames::Optimizations::MachineLearning::ftol, std::move(ftol));
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_MLREGRESSIONSETTINGS_H
