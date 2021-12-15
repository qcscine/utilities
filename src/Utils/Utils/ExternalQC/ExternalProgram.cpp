/**
 * @filesetWorkingDirectory
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "ExternalProgram.h"
#include "Exceptions.h"
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>
#include <boost/process.hpp>
#include <iostream>

namespace bp = boost::process;

namespace Scine {
namespace Utils {
namespace ExternalQC {

void ExternalProgram::setWorkingDirectory(const std::string& directory) {
  workingDirectory_ = NativeFilenames::addTrailingSeparator(directory);
}

const std::string& ExternalProgram::getWorkingDirectory() const {
  return workingDirectory_;
}

std::string ExternalProgram::generateFullFilename(const std::string& filename) const {
  return NativeFilenames::combinePathSegments(workingDirectory_, filename);
}

void ExternalProgram::createWorkingDirectory() {
  if (!workingDirectory_.empty()) {
    FilesystemHelpers::createDirectories(workingDirectory_);
  }
}

void ExternalProgram::executeCommand(const std::string& command) const {
  executeCommand(command, "", "");
}

void ExternalProgram::executeCommand(const std::string& command, const std::string& outputFile) const {
  executeCommand(command, "", outputFile);
}

void ExternalProgram::executeCommand(const std::string& command, const std::string& inputFile,
                                     const std::string& outputFile) const {
  int ret = executeCommandImpl(command, inputFile, outputFile);

  if (ret != 0) {
    throw UnsuccessfulSystemCommand(command, inputFile, outputFile);
  }
}

int ExternalProgram::executeCommandImpl(const std::string& command, const std::string& inputFile,
                                        const std::string& outputFile) const {
  bool hasInput = !inputFile.empty();
  bool hasOutput = !outputFile.empty();

  auto workingDirectory = bp::start_dir(workingDirectory_);
  if (workingDirectory_.empty()) {
    workingDirectory = bp::start_dir(FilesystemHelpers::currentDirectory());
  }

  if (hasInput && hasOutput) {
    return bp::system(command, bp::std_out > outputFile, bp::std_in < inputFile, workingDirectory);
  }
  if (hasInput) {
    return bp::system(command, bp::std_in < inputFile, workingDirectory);
  }
  if (hasOutput) {
    return bp::system(command, bp::std_out > outputFile);
  }
  if (hasOutput) {
    return bp::system(command, bp::std_out > outputFile, workingDirectory);
  }
  return bp::system(command, workingDirectory);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
