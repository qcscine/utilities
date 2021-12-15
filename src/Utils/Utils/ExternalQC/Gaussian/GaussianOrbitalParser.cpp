/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaussianOrbitalParser.h"
#include "GaussianFileConverter.h"
#include <cstdio>

namespace Scine {
namespace Utils {
namespace ExternalQC {

GaussianOrbitalParser::GaussianOrbitalParser(const std::string& chkFileBase, const std::string& workingDirectory,
                                             const std::string& gaussianDirectory) {
  std::string fchkFile =
      GaussianFileConverter::generateFormattedCheckpointFile(chkFileBase, workingDirectory, gaussianDirectory);
  openFile(fchkFile);
  readOrbitals();
  closeFile();
  // The fchk file could be reused to update the orbitals later, but it is not to avoid chk/fchk discrepancies
  std::remove(fchkFile.c_str());
}

void GaussianOrbitalParser::openFile(const std::string& file) {
  fchkFile_.open(file);
  if (fchkFile_.fail() | fchkFile_.bad())
    throw std::runtime_error("Could not open" + file);
}

void GaussianOrbitalParser::closeFile() {
  fchkFile_.close();
}

void GaussianOrbitalParser::readOrbitals() {
  std::string line;
  while (std::getline(fchkFile_, line)) {
    checkNumberOrbitalsLine(line);
    checkNumberAlphaElectronsLine(line);
    checkNumberBetaElectronsLine(line);
    checkAlphaOrbitals(line);
    checkBetaOrbitals(line);
  }
  isUnrestricted_ = (betaCoefficients_.size() > 0);
  createMolecularOrbitals();
  createOccupation();
}

void GaussianOrbitalParser::checkNumberOrbitalsLine(const std::string& line) {
  std::string pattern = "Number of basis functions";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    std::stringstream ss(line.substr(pattern.length()));
    char c; // 'I'
    ss >> c >> nOrbitals_;
  }
}

void GaussianOrbitalParser::checkNumberAlphaElectronsLine(const std::string& line) {
  std::string pattern = "Number of alpha electrons";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    std::stringstream ss(line.substr(pattern.length()));
    char c; // 'I'
    ss >> c >> nAlpha_;
  }
}

void GaussianOrbitalParser::checkNumberBetaElectronsLine(const std::string& line) {
  std::string pattern = "Number of beta electrons";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    std::stringstream ss(line.substr(pattern.length()));
    char c; // 'I'
    ss >> c >> nBeta_;
  }
}

void GaussianOrbitalParser::checkAlphaOrbitals(const std::string& line) {
  std::string pattern = "Alpha MO coefficients";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    fillCoefficients(alphaCoefficients_);
  }
}

void GaussianOrbitalParser::checkBetaOrbitals(const std::string& line) {
  std::string pattern = "Beta MO coefficients";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    fillCoefficients(betaCoefficients_);
  }
}

void GaussianOrbitalParser::fillCoefficients(std::vector<double>& coefficients) {
  auto numberCoefficients = nOrbitals_ * nOrbitals_;
  unsigned numberValuesPerLine = 5;
  coefficients.resize(numberCoefficients);
  for (unsigned i = 0; i < numberCoefficients; i += numberValuesPerLine) {
    std::string line;
    std::getline(fchkFile_, line);
    std::stringstream lineStream(line);
    for (unsigned j = i; j < i + numberValuesPerLine && j < numberCoefficients; ++j) {
      lineStream >> coefficients[j];
    }
  }
}

void GaussianOrbitalParser::createMolecularOrbitals() {
  auto numberCoefficients = nOrbitals_ * nOrbitals_;
  bool foundNumberOrbitals = nOrbitals_ != 0;
  bool foundOrbitalCoefficients = alphaCoefficients_.size() == numberCoefficients &&
                                  (betaCoefficients_.size() == numberCoefficients || !isUnrestricted_);
  if (!(foundNumberOrbitals && foundOrbitalCoefficients)) {
    throw std::runtime_error("Parsing orbitals from Gaussian calculation unsuccessful.");
  }
  Eigen::Map<Eigen::MatrixXd> Ca(alphaCoefficients_.data(), nOrbitals_, nOrbitals_);
  if (isUnrestricted_) {
    Eigen::Map<Eigen::MatrixXd> Cb(betaCoefficients_.data(), nOrbitals_, nOrbitals_);
    orbitals_ = MolecularOrbitals::createFromUnrestrictedCoefficients(Ca, Cb);
  }
  else {
    orbitals_ = MolecularOrbitals::createFromRestrictedCoefficients(Ca);
  }
}

void GaussianOrbitalParser::createOccupation() {
  if (isUnrestricted_) {
    occupation_.fillLowestUnrestrictedOrbitals(nAlpha_, nBeta_);
  }
  else {
    occupation_.fillLowestRestrictedOrbitalsWithElectrons(nAlpha_ + nBeta_);
  }
}

unsigned GaussianOrbitalParser::getNumberOrbitals() const {
  return nOrbitals_;
}

unsigned GaussianOrbitalParser::getNumberAlphaElectrons() const {
  return nAlpha_;
}

unsigned GaussianOrbitalParser::getNumberBetaElectrons() const {
  return nBeta_;
}

const MolecularOrbitals& GaussianOrbitalParser::getOrbitals() const {
  return orbitals_;
}

const LcaoUtils::ElectronicOccupation& GaussianOrbitalParser::getElectronicOccupation() const {
  return occupation_;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine