/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TURBOMOLEFILES_H
#define TURBOMOLEFILES_H

#include "Utils/IO/NativeFilenames.h"
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

struct TurbomoleFiles {
  std::string filenameBase;
  std::string alphaFile;
  std::string betaFile;
  std::string controlFile;
  std::string energyFile;
  std::string hessianFile;
  std::string gradientFile;
  std::string alphaBakFile;
  std::string betaBakFile;
  std::string ridftFile;
  std::string ridftBakFile;
  std::string defineInputFile;
  std::string coordFile;
  std::string solvationInputFile;
};

inline void setCorrectTurbomoleFileNames(TurbomoleFiles& files, const std::string workingDirectory) {
  files.filenameBase = workingDirectory;
  files.coordFile = NativeFilenames::combinePathSegments(files.filenameBase, "coord");
  files.defineInputFile = NativeFilenames::combinePathSegments(files.filenameBase, "tm.input");
  files.alphaFile = NativeFilenames::combinePathSegments(files.filenameBase, "alpha");
  files.betaFile = NativeFilenames::combinePathSegments(files.filenameBase, "beta");
  files.controlFile = NativeFilenames::combinePathSegments(files.filenameBase, "control");
  files.energyFile = NativeFilenames::combinePathSegments(files.filenameBase, "energy");
  files.hessianFile = NativeFilenames::combinePathSegments(files.filenameBase, "hessian");
  files.gradientFile = NativeFilenames::combinePathSegments(files.filenameBase, "gradient");
  files.alphaBakFile = NativeFilenames::combinePathSegments(files.filenameBase, "alpha.bak");
  files.betaBakFile = NativeFilenames::combinePathSegments(files.filenameBase, "beta.bak");
  files.ridftFile = NativeFilenames::combinePathSegments(files.filenameBase, "ridft.out");
  files.ridftBakFile = NativeFilenames::combinePathSegments(files.filenameBase, "ridft_unperturbed.out");
  files.solvationInputFile = NativeFilenames::combinePathSegments(files.filenameBase, "cosmoprep.inp");
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // TURBOMOLEFILES_H