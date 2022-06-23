/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NativeFilenames.h"
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace Scine {
namespace Utils {

namespace boostfs = boost::filesystem;

char NativeFilenames::getDirectorySeparatorChar() {
  return boostfs::path::preferred_separator;
}

std::string NativeFilenames::getDirectorySeparatorString() {
  return std::string(1, getDirectorySeparatorChar());
}

std::string NativeFilenames::getParentDirectory(const std::string& path) {
  boostfs::path p = path;
  return p.parent_path().string();
}

std::string NativeFilenames::combinePathSegmentsImpl(const std::string& left, const std::string& right) {
  boostfs::path p = left;
  p /= right;
  return p.string();
}

std::string NativeFilenames::removeTrailingSeparator(const std::string& path) {
  boostfs::path p = path;
  p.remove_trailing_separator();
  return p.string();
}

std::string NativeFilenames::addTrailingSeparator(const std::string& path) {
  return removeTrailingSeparator(path) + getDirectorySeparatorString();
}

std::string NativeFilenames::removeExtension(const std::string& filename) {
  size_t lastdot = filename.find_last_of('.');
  if (lastdot == std::string::npos) {
    return filename;
  }
  return filename.substr(0, lastdot);
}

std::string NativeFilenames::createRandomDirectoryName(const std::string& baseDirectory) {
  boost::uuids::uuid uuid = boost::uuids::random_generator()();
  std::string uuidString = boost::uuids::to_string(uuid);

  auto directoryName = NativeFilenames::combinePathSegments(baseDirectory, uuidString);
  return NativeFilenames::addTrailingSeparator(directoryName);
}

} // namespace Utils
} // namespace Scine
