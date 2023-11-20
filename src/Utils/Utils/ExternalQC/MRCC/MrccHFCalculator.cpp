/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "MrccHFCalculator.h"

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string MrccHFCalculator::getMethodFamily() const {
  return "HF";
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
