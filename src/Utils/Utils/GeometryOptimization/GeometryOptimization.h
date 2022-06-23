/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GEOMETRY_OPTIMIZATION_H_
#define UTILS_GEOMETRY_OPTIMIZATION_H_

#include <Core/Interfaces/Calculator.h>
#include <Utils/GeometryOptimization/AfirOptimizer.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/GeometryOptimization/IrcOptimizer.h>
#include <Utils/GeometryOptimization/NtOptimizer.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/Optimizer/GradientBased/GradientBasedCheck.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace GeometryOptimization {

/**
 * @brief Comparison of the energy convergence threshold of the energy and all optimization convergence criteria.
 *
 * @return bool whether the energy convergence is less or equal to the optimization convergence.
 */
inline bool energyThresholdIsLowEnough(const double energyThreshold, const GradientBasedCheck& check) {
  return !(energyThreshold > check.deltaValue || energyThreshold > check.gradMaxCoeff ||
           energyThreshold > check.gradRMS || energyThreshold > check.stepMaxCoeff || energyThreshold > check.stepRMS);
}

/**
 * @brief A templated check whether the convergence settings make sense with the given calculator accuracy and
 * convergence settings where the template defines a geometry optimizer
 *
 * @tparam OptimizerType Expects an Optimizer class, which can access the settings of it's calculator and
 * a GradientBasedCheck.
 *
 * @return bool whether the calculator is more or equally accurate as the convergence criteria
 */
template<class OptimizerType>
bool settingsMakeSense(const OptimizerType& optimizer) {
  const auto settings = optimizer.getCalculatorSettings();
  if (!settings || !settings->valueExists(SettingsNames::selfConsistenceCriterion)) {
    throw std::runtime_error("Necessary '" + std::string(SettingsNames::selfConsistenceCriterion) +
                             "' not implemented in the applied calculator.");
  }
  const double energyThreshold = settings->getDouble(SettingsNames::selfConsistenceCriterion);
  const GradientBasedCheck check = optimizer.getConvergenceCheck();
  return energyThresholdIsLowEnough(energyThreshold, check);
}

// @brief The Convergence of the NtOptimizer is independent of the energy convergence threshold
template<>
inline bool settingsMakeSense(const NtOptimizer& /* optimizer */) {
  return true; // convergence is independent of calculator
}

/**
 * @brief A specific sanity check for the QmmmOptimizer, because it has two different convergence checks and is a
 * templated class without a baseclass.
 *
 * @tparam OptimizerType Expects any of the Optimizer classes.
 * @return bool whether the calculator is more or equally accurate as the convergence criteria of both convergence checks
 */
template<class OptimizerType>
bool settingsMakeSense(const QmmmGeometryOptimizer<OptimizerType>& optimizer) {
  const auto settings = optimizer->getCalculatorSettings();
  if (!settings || !settings->valueExists(SettingsNames::selfConsistenceCriterion)) {
    throw std::runtime_error("Necessary '" + std::string(SettingsNames::selfConsistenceCriterion) +
                             "' not implemented in the applied calculator.");
  }
  const double energyThreshold = settings->getDouble(SettingsNames::selfConsistenceCriterion);
  const GradientBasedCheck check1 = optimizer->fullOptimizer->getConvergenceCheck();
  const GradientBasedCheck check2 = optimizer->mmOptimizer->getConvergenceCheck();
  return (energyThresholdIsLowEnough(energyThreshold, check1) && energyThresholdIsLowEnough(energyThreshold, check2));
}

} // namespace GeometryOptimization
} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRY_OPTIMIZATION_H_
