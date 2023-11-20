/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometryOptimization/AfirConvergenceCheck.h"
#include "Utils/GeometryOptimization/AfirOptimizerBase.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Settings.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

bool AfirConvergenceCheck::checkConvergence(const Eigen::VectorXd& parameter, double value, const Eigen::VectorXd& gradient) {
  bool gradientConvergence = GradientBasedCheck::checkConvergence(parameter, value, gradient);
  // If required evaluate fragment distance
  bool exceedsDistanceThreshold = false;
  if (afirUseMaxFragmentDistance) {
    Utils::PositionCollection positions;
    if (this->transformation) {
      positions = this->transformation->coordinatesToCartesian(parameter);
    }
    else {
      positions = Eigen::Map<const Utils::PositionCollection>(parameter.data(), parameter.size() / 3, 3);
    }
    exceedsDistanceThreshold = this->checkExceedsDistanceThreshold(positions);
  }
  // If either converged or distance threshold exceeded return true
  return gradientConvergence || exceedsDistanceThreshold;
}

bool AfirConvergenceCheck::checkExceedsDistanceThreshold(Utils::PositionCollection& positions) {
  return this->calculateMinFragmentDistance(positions) >= afirMaxFragmentDistance;
}

double AfirConvergenceCheck::calculateMinFragmentDistance(Utils::PositionCollection& positions) {
  double minDistance = std::numeric_limits<double>::max();
  double distance;
  if (!lhsList.empty() && !rhsList.empty()) {
    for (int i : lhsList) {
      for (int j : rhsList) {
        distance = (positions.row(i) - positions.row(j)).norm();
        // Set minimum distance to minimum of previous and current value
        minDistance = std::min(minDistance, distance);
      }
    }
  }
  // If at least one of the lists is empty this check does not make sense
  else {
    throw std::runtime_error(
        "AFIR LHS and/or RHS list empty. Cannot use a minimum fragment distance without defined fragments.");
  }
  return minDistance;
}

void AfirConvergenceCheck::addAfirSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const {
  // Add afir specific setting descriptors
  UniversalSettings::BoolDescriptor afir_use_max_fragment_distance(
      "Whether to stop the AFIR optimization when exceeding a maximum interfragment distance.");
  afir_use_max_fragment_distance.setDefaultValue(afirUseMaxFragmentDistance);
  collection.push_back(AfirConvergenceCheck::afirUseMaxFragmentDistanceKey, afir_use_max_fragment_distance);
  UniversalSettings::DoubleDescriptor afir_max_fragment_distance(
      "Interfragment distance upon exceeding which the AFIR optimization is stopped.");
  afir_max_fragment_distance.setDefaultValue(afirMaxFragmentDistance);
  collection.push_back(AfirConvergenceCheck::afirMaxFragmentDistanceKey, afir_max_fragment_distance);
}

void AfirConvergenceCheck::applyAfirSettings(const Settings& settings) {
  // Apply AFIR specific convergence settings
  afirUseMaxFragmentDistance = settings.getBool(AfirConvergenceCheck::afirUseMaxFragmentDistanceKey);
  afirMaxFragmentDistance = settings.getDouble(AfirConvergenceCheck::afirMaxFragmentDistanceKey);
}

} // namespace Utils
} // namespace Scine
