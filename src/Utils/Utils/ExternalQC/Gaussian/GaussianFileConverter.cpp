/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "boost/filesystem.hpp"
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/ExternalQC/Gaussian/GaussianFileConverter.h>
#include <Utils/IO/NativeFilenames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace GaussianFileConverter {

std::string generateFormattedCheckpointFile(const std::string& fileBase, const std::string& workingDirectory,
                                            const std::string& gaussianDirectory) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDirectory);
  externalProgram.createWorkingDirectory();
  std::string fullChkFile = externalProgram.generateFullFilename(fileBase + ".chk");
  std::string fullFchkFile = externalProgram.generateFullFilename(fileBase + ".fchk");
  std::string formchkExecutable = NativeFilenames::combinePathSegments(gaussianDirectory, "formchk");
  if (!boost::filesystem::exists(fullChkFile)) {
    throw GaussianCheckpointFileNotFoundException("Checkpoint file " + fullChkFile + " does not exist.");
  }
  try {
    // Chk file is not accepted as std::in in formchk
    externalProgram.executeCommand(formchkExecutable + " " + fullChkFile, fullFchkFile);
  }
  catch (const std::exception& ex) {
    throw GaussianCheckpointConversionException("Conversion of " + fullChkFile +
                                                " to a formatted checkpoint file failed: " + ex.what());
  }
  catch (...) {
    throw GaussianCheckpointConversionException("Conversion of " + fullChkFile + " to a formatted checkpoint file failed.");
  }
  return fullFchkFile;
}

std::string generateCheckpointFile(const std::string& fileBase, const std::string& workingDirectory,
                                   const std::string& gaussianDirectory) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDirectory);
  externalProgram.createWorkingDirectory();
  std::string fullFchkFile = externalProgram.generateFullFilename(fileBase + ".fchk");
  std::string fullChkFile = externalProgram.generateFullFilename(fileBase + ".chk");
  std::string unfchkExecutable = NativeFilenames::combinePathSegments(gaussianDirectory, "unfchk");
  if (!boost::filesystem::exists(fullFchkFile)) {
    throw GaussianCheckpointFileNotFoundException("Formatted checkpoint file " + fullFchkFile + " does not exist.");
  }
  try {
    externalProgram.executeCommand(unfchkExecutable + " " + fullFchkFile, fullChkFile);
  }
  catch (const std::exception& ex) {
    throw GaussianCheckpointConversionException("Conversion of " + fullFchkFile +
                                                " to a binary checkpoint file failed: " + ex.what());
  }
  catch (...) {
    throw GaussianCheckpointConversionException("Conversion of " + fullFchkFile + " to a binary checkpoint file failed.");
  }

  return fullChkFile;
}

} // namespace GaussianFileConverter
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
