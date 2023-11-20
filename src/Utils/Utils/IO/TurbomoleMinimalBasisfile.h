/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_IO_TURBOMOLE_MINIMAL_BASIS_FILE_H
#define INCLUDE_UTILS_IO_TURBOMOLE_MINIMAL_BASIS_FILE_H

#include <Utils/DataStructures/AtomicGtos.h>
#include <unordered_map>

namespace Scine {
namespace Utils {

std::unordered_map<int, AtomicGtos> readTurbomoleBasisfile(const std::string& parameterFile);

} // namespace Utils
} // namespace Scine

#endif
