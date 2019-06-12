/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_FILESYSTEMHELPERS_H
#define UTILS_FILESYSTEMHELPERS_H

#include <string>

namespace Scine {
namespace Utils {

/**
 * @class FilesystemHelpers FilesystemHelpers.h
 * @brief This class contains utility functions for dealing with files and directories.
 *        Hides the use of Boost::filesystem.
 */
class FilesystemHelpers {
 public:
  /**
   * @brief Creates directories for fullDirectoryPath if they do not exist yet.
   */
  static void createDirectories(const std::string& fullDirectoryPath);
  /**
   * @brief Empty some given directory.
   */
  static void emptyDirectory(const std::string& directory);
  /**
   * @brief Check whether a given path is a directory.
   */
  static bool isDirectory(const std::string& path);
  /**
   * @brief Generate a temporary directory
   */
  static std::string systemTmpDirectory();
  /**
   * @brief Get the name of the current directory.
   */
  static std::string currentDirectory();
  /**
   * @brief Copy file (overwrites it if necessary)
   */
  static void copyFile(const std::string& from, const std::string& to);
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_FILESYSTEMHELPERS_H