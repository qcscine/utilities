/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GAUSSIANORBITALPARSER_H
#define UTILS_GAUSSIANORBITALPARSER_H

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <fstream>
#include <vector>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class GaussianOrbitalWriter GaussianOrbitalWriter.h
 * @brief Class to read a formatted Gaussian checkpoint file and extract the orbitals
 *
 * @param chkFileBase The base of the name of the checkpoint file.
 * @param workingDirectory The working directory with the relevant checkpoint file.
 * @param gaussianDirectory The directory the "formchk" executable is located in.
 */
class GaussianOrbitalParser {
 public:
  GaussianOrbitalParser(const std::string& chkFileBase, const std::string& workingDirectory,
                        const std::string& gaussianDirectory);

  unsigned getNumberOrbitals() const;
  unsigned getNumberAlphaElectrons() const;
  unsigned getNumberBetaElectrons() const;
  const MolecularOrbitals& getOrbitals() const;
  /**
   * @brief Gets the electronic occupation
   *
   * @note The orbitals are assumed to be filled from the bottom up
   */
  const LcaoUtils::ElectronicOccupation& getElectronicOccupation() const;

 private:
  void openFile(const std::string& file);
  void closeFile();
  void readOrbitals();
  void checkNumberOrbitalsLine(const std::string& line);
  void checkNumberAlphaElectronsLine(const std::string& line);
  void checkNumberBetaElectronsLine(const std::string& line);
  void checkAlphaOrbitals(const std::string& line);
  void checkBetaOrbitals(const std::string& line);
  void fillCoefficients(std::vector<double>& coefficients);
  void createMolecularOrbitals();
  void createOccupation();

  std::ifstream fchkFile_;
  MolecularOrbitals orbitals_;
  unsigned nOrbitals_ = 0;
  unsigned nAlpha_ = 0;
  unsigned nBeta_ = 0;
  std::vector<double> alphaCoefficients_;
  std::vector<double> betaCoefficients_;
  LcaoUtils::ElectronicOccupation occupation_;
  bool isUnrestricted_ = false;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANORBITALPARSER_H
