/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_FORMATTEDSTRING_H_
#define UTILS_FORMATTEDSTRING_H_

#include <stdio.h>
#include <memory>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief Generates a std::string from a format expression.
 *
 * Note that std::string arguments for '%s' expressions need to be given as
 * 'char*' e.g. 'foo' as 'foo.c_str()'.
 *
 * @tparam Args The variable types to be formatted.
 * @param format The string containing format identifiers.
 * @param args The variables to be formatted.
 * @return std::string The resulting string.
 */
template<typename... Args>
std::string format(const std::string& format, Args... args) {
  size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1;
  std::unique_ptr<char[]> buf(new char[size]);
  snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(), buf.get() + size - 1);
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_FORMATTEDSTRING_H_
