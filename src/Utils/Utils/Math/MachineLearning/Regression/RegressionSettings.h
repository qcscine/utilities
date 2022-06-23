/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_MLREGRESSIONSETTINGS_H
#define UTILSOS_MLREGRESSIONSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace MachineLearning {

namespace SettingsNames {

static constexpr const char* restartOptimization = "restart_optimization";
static constexpr const char* numRestarts = "num_restarts";
static constexpr const char* maxIterations = "max_iterations";
static constexpr const char* maxLinesearch = "max_linesearch";
static constexpr const char* convergenceTolerance = "convergence_tolerance";
static constexpr const char* ftol = "linesearch_tolerance";

} // namespace SettingsNames

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
  settings.push_back(SettingsNames::restartOptimization, std::move(restartOptimization));
}

inline void RegressionSettings::setNumRestarts(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numRestarts("The number of restarts for hyperparameter optimization.");
  numRestarts.setMinimum(1);
  numRestarts.setDefaultValue(10);
  settings.push_back(SettingsNames::numRestarts, std::move(numRestarts));
}

inline void RegressionSettings::setMaxIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxIterations(
      "The maximum number of iterations. Note that setting this parameter to zero continues an optimization process "
      "until convergence or error.");
  maxIterations.setMinimum(1);
  maxIterations.setDefaultValue(1000);
  settings.push_back(SettingsNames::maxIterations, std::move(maxIterations));
}

inline void RegressionSettings::setMaxLinesearch(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxLinesearch("The maximum number of trials for the line search.");
  maxLinesearch.setMinimum(1);
  maxLinesearch.setDefaultValue(20000);
  settings.push_back(SettingsNames::maxLinesearch, std::move(maxLinesearch));
}

inline void RegressionSettings::setConvergenceTolerance(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor convergenceTolerance("The absolute tolerance for convergence test.");
  convergenceTolerance.setDefaultValue(1e-6);
  settings.push_back(SettingsNames::convergenceTolerance, std::move(convergenceTolerance));
}

inline void RegressionSettings::setLinesearchTolerance(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor ftol("A parameter to control the accuracy of the line search routine.");
  ftol.setDefaultValue(0.001);
  settings.push_back(SettingsNames::ftol, std::move(ftol));
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_MLREGRESSIONSETTINGS_H
