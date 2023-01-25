/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_STRINGS_H
#define INCLUDE_UTILS_STRINGS_H

#include <algorithm>
#include <string>

namespace Scine {
namespace Utils {

static bool caseInsensitiveEqual(const std::string& a, const std::string& b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](char a, char b) { return tolower(a) == tolower(b); });
}

} // namespace Utils
} // namespace Scine

#endif
