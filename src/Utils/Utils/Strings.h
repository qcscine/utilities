/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_STRINGS_H
#define INCLUDE_UTILS_STRINGS_H

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief Compares two strings case insensitive
 * @param a The first string
 * @param b The second string
 * @return Whether they are equal
 */
static bool caseInsensitiveEqual(const std::string& a, const std::string& b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](char a, char b) { return tolower(a) == tolower(b); });
}

/**
 * @brief Splits a given string based on spaces into a vector and filters out any blank entries in the resulting vector
 * basically a small wrapper around boost::split
 * @param line The string to be split
 * @return The vector of strings in the line, no entry can be an blank space
 */
static std::vector<std::string> splitOnSpaceWithoutResultingSpace(const std::string& line) {
  std::vector<std::string> lineSplitted;
  boost::split(lineSplitted, line, boost::is_any_of(" "), boost::token_compress_on);
  std::vector<std::string> strippedResult;
  for (const auto& value : lineSplitted) {
    if (value.find_first_of(' ') == std::string::npos && !value.empty()) {
      strippedResult.push_back(value);
    }
    else if (std::any_of(value.begin(), value.end(), [](const char v) { return v != ' '; })) {
      // entry in boost split result includes a space, but also a char that is not a space, something went wrong
      throw std::runtime_error("Boost split gave the split entry '" + value + "', which is not a valid result");
    }
  }
  return strippedResult;
}

} // namespace Utils
} // namespace Scine

#endif
