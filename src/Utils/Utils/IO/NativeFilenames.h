/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NATIVEFILENAMES_H
#define UTILS_NATIVEFILENAMES_H

#include <string>

namespace Scine {
namespace Utils {

/*!
 * @brief This class contains utility functions for storing file paths in strings, and cares for cross-platform issues.
 * Hides the use of Boost::filesystem.
 */
class NativeFilenames {
 public:
  static char getDirectorySeparatorChar();
  static std::string getDirectorySeparatorString();

  /*! concatenates path segments by adding a separator if necessary. */
  template<typename T, typename U>
  static std::string combinePathSegments(const T& l, const U& r);
  /*! Concatenates multiple path segments contained in the argument pack. */
  template<typename T, typename... Ts>
  static std::string combinePathSegments(const T& l, const Ts&... r);

  /*! remove the directory separator at the end of the string if present. */
  static std::string removeTrailingSeparator(const std::string& path);
  /*! add a directory separator at the end of the string if none is present. */
  static std::string addTrailingSeparator(const std::string& path);
  /*! @brief removes the extension from the file called filename. */
  static std::string removeExtension(const std::string& filename);

 private:
  /*! concatenates path segments by adding a separator if necessary. */
  static std::string combinePathSegmentsImpl(const std::string& left, const std::string& right);
};

template<typename T, typename U>
inline std::string NativeFilenames::combinePathSegments(const T& l, const U& r) {
  return combinePathSegmentsImpl(l, r);
}

template<typename T, typename... Ts>
std::string NativeFilenames::combinePathSegments(const T& l, const Ts&... r) {
  return combinePathSegments(l, combinePathSegments(r...));
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_NATIVEFILENAMES_H
