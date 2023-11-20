/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GAUSSIANORBITALWRITER_H
#define UTILS_GAUSSIANORBITALWRITER_H

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <fstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class GaussianOrbitalWriter GaussianOrbitalWriter.h
 * @brief Class to write an updated formatted Gaussian checkpoint file from orbitals and an input checkpoint file.
 *
 * @param orbitals The orbitals to be written.
 *
 */
class GaussianOrbitalWriter {
 public:
  GaussianOrbitalWriter(const MolecularOrbitals& orbitals);
  /**
   * @brief Updates the orbitals in the checkpoint file.
   *
   * @param chkFileBase The base of the file name of the chk file
   * @param workingDirectory The directory the chk file is located at
   * @param gaussianDirectory The directory the "formchk" and "unchk" executables are located in.
   *
   * @note Only the orbital coefficient block is updated.
   * All other parts of the checkpoint file are left unchanged.
   */
  void updateCheckpointFile(const std::string& chkFileBase, const std::string& workingDirectory,
                            const std::string& gaussianExecutable);

 private:
  void write();
  void openInFile(const std::string& file);
  void openOutFile(const std::string& file);
  void closeFchkFiles();
  void writeRestrictedOrbitals(const std::string& line);
  void writeAlphaOrbitals(const std::string& line);
  void writeBetaOrbitals(const std::string& line);
  void ignoreInputLines();
  void writeCoefficients(const Eigen::MatrixXd& c);
  // Converts double to string in scientific notation as used in Gaussian fchk file
  std::string convertToScientificNotation(const double number) const;

  const MolecularOrbitals& orbitals_;
  std::ifstream fchkin_;
  std::ofstream fchkout_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANORBITALWRITER
