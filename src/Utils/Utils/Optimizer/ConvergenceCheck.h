/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CONVERGENCECHECK_H_
#define UTILS_CONVERGENCECHECK_H_

#include "Utils/Settings.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {
/**
 * @brief The base class for all convergence checks.
 *
 * Similar to the actual Optimizer the convergence check is abstracted into
 * a class that shall allow both usage via Settings object and as a 'raw' class
 * using the public members to adjust the convergence criteria.\n
 * To this end the base class defines and enforces the presence of the functions
 * that parse the Settings into these members and vice versa.\n
 * In general, it is adviseable to derive classes that also only require the use of
 * a default constructor, after all the convergence checks are mostly defined by
 * a set of thresholds and any data should be required in the actual optimization cycles
 * and not upon creation of the derived object.
 */
class ConvergenceCheck {
 public:
  /// @brief Default Constructor.
  ConvergenceCheck() = default;
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include these ConvergenceCheck's options.
   *
   * This function has to be implemented in derived classes in order to allow for
   * all ConvergenceChecks to be used and configured at runtime by end users.
   *
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const = 0;
  /**
   * @brief Updates the ConvergenceCheck's options with those values given in the Settings.
   *
   * This function has to be implemented in derived classes in order to allow for
   * all ConvergenceChecks to be used and configured at runtime by end users.
   *
   * Note: Additional settings that are not used by the optimizer shall not generate errors,
   *       or warnings, they shall simply be ignored.
   *
   * @param settings The settings to update the option of this optimizer with.
   */
  virtual void applySettings(const Settings& settings) = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_CONVERGENCECHECK_H_
