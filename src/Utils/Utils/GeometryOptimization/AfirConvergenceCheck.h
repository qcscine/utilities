/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIRCONVERGENCECHECK_H_
#define UTILS_AFIRCONVERGENCECHECK_H_

#include "Utils//Typenames.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"

namespace Scine {
namespace Utils {

/**
 * @brief A convergence check based on inter-fragment distance and a GradientBasedCheck
 *
 * Checks for convergence with a GradientBasedCheck with the additional option to
 * calculate the minimum distance between two AFIR fragments and compare it to a
 * threshold value. Convergence is signalled, if the distance exceeds the
 * threshold value and/or the gradient based convergence criterion is fulfilled.
 */

class AfirConvergenceCheck : public GradientBasedCheck {
 public:
  static constexpr const char* afirUseMaxFragmentDistanceKey = "afir_use_max_fragment_distance";
  static constexpr const char* afirMaxFragmentDistanceKey = "afir_max_fragment_distance";
  /// @brief Default constructor.
  AfirConvergenceCheck() = default;
  /// @brief Default destructor.
  virtual ~AfirConvergenceCheck() = default;

  /// @brief Whether to eventually stop the optimization based on the inter-fragment distance
  bool afirUseMaxFragmentDistance = false;
  /// @brief The threshold for the maximum inter-fragment distance not to be exceeded
  double afirMaxFragmentDistance = 8.0;

  /**
   * @brief Checks for gradient convergence and inter-fragment distance
   *
   * @param parameter The current parameters.
   * @param value     The current value.
   * @param gradient  The current gradient.
   * @return true     If converged or fragment distance triggering calculation stop.
   * @return false    If neither converged nor fragment distance triggering calculation stop.
   */
  bool checkConvergence(const Eigen::VectorXd& parameter, double value, const Eigen::VectorXd& gradient) override;

  /// @brief Adds AFIR specific settings; See Scine::Utils::ConvergenceCheck::addSettingsDescriptors()
  virtual void addAfirSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final;
  /// @brief Applies AFIR specific settings; See Scine::Utils::ConvergenceCheck::applySettings()
  virtual void applyAfirSettings(const Settings& s) final;
  /// @brief List of first fragment atom indices;
  std::vector<int> lhsList = {};
  /// @brief List of second fragment atom indices;
  std::vector<int> rhsList = {};
  /// @brief Whether the optimization uses internal coordinates;
  bool transformCoordinates = true;
  /// @brief Transformation matrix;
  std::shared_ptr<InternalCoordinates> transformation = nullptr;

 private:
  /**
   * @brief Calculates the shortest distance between the positions of two position collections
   *
   * @param positions The current positions.
   * @return Minimum distance between the fragments.
   */
  double calculateMinFragmentDistance(Utils::PositionCollection& positions);

  /**
   * @brief Checks whether the inter-fragment distance exceeds the maximum value
   *
   * @param positions The current positions.
   * @return true If he minimum fragment distance exceeds the threshold
   * @return false If the minimum fragment distance does not exceed the threshold
   */
  bool checkExceedsDistanceThreshold(Utils::PositionCollection& positions);
};
} // namespace Utils
} // namespace Scine

#endif // UTILS_AFIRCONVERGENCECHECK_H_
