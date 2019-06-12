/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "FilesystemHelpers.h"
#include <boost/filesystem.hpp>

namespace boostfs = boost::filesystem;

namespace Scine {
namespace Utils {

void FilesystemHelpers::createDirectories(const std::string& fullDirectoryPath) {
  try {
    boostfs::create_directories(fullDirectoryPath);
  }
  catch (std::exception&) {
    throw std::runtime_error("Could not create the directory " + fullDirectoryPath);
  }
}

void FilesystemHelpers::emptyDirectory(const std::string& directory) {
  boostfs::path directoryPath(directory);

  if (!boostfs::exists(directoryPath)) {
    throw std::runtime_error("Cannot empty directory \"" + directory + "\": it does not exist.");
  }

  for (boostfs::directory_iterator end_dir_it, it(directoryPath); it != end_dir_it; ++it) {
    boostfs::remove_all(it->path());
  }
}

bool FilesystemHelpers::isDirectory(const std::string& path) {
  return boostfs::is_directory(path);
}

std::string FilesystemHelpers::systemTmpDirectory() {
  auto tmpDirectory = boostfs::temp_directory_path();
  return tmpDirectory.string();
}

std::string FilesystemHelpers::currentDirectory() {
  boostfs::path fullPath(boostfs::current_path());
  return fullPath.string();
}

void FilesystemHelpers::copyFile(const std::string& from, const std::string& to) {
  boostfs::path fromFile{from};
  boostfs::path toFile{to};

  try {
    boostfs::copy_file(fromFile, toFile, boostfs::copy_option::overwrite_if_exists);
  }
  catch (std::exception&) {
    throw std::runtime_error("Could not copy the file \"" + fromFile.string() + "\" to \"" + toFile.string() + "\".");
  }
}

} // namespace Utils
} // namespace Scine