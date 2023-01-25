/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaussianOrbitalWriter.h"
#include "GaussianFileConverter.h"
#include <cstdio>
#include <iomanip>
#include <sstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

GaussianOrbitalWriter::GaussianOrbitalWriter(const MolecularOrbitals& orbitals) : orbitals_(orbitals) {
}

void GaussianOrbitalWriter::updateCheckpointFile(const std::string& chkFileBase, const std::string& workingDirectory,
                                                 const std::string& gaussianDirectory) {
  std::string inFchkFile =
      GaussianFileConverter::generateFormattedCheckpointFile(chkFileBase, workingDirectory, gaussianDirectory);
  openInFile(inFchkFile);
  std::string newFchkFile = inFchkFile + "_new";
  openOutFile(newFchkFile);
  write();
  closeFchkFiles();
  std::rename(newFchkFile.c_str(), inFchkFile.c_str());
  std::string outChkFile = GaussianFileConverter::generateCheckpointFile(chkFileBase, workingDirectory, gaussianDirectory);
  std::remove(inFchkFile.c_str());
}

void GaussianOrbitalWriter::openInFile(const std::string& file) {
  fchkin_.open(file);
  if (fchkin_.fail() || fchkin_.bad())
    throw std::runtime_error("Could not open " + file);
}

void GaussianOrbitalWriter::openOutFile(const std::string& file) {
  fchkout_.open(file);
  if (fchkout_.fail() || fchkout_.bad())
    throw std::runtime_error("Could not open " + file);
}

void GaussianOrbitalWriter::write() {
  std::string line;
  while (std::getline(fchkin_, line)) {
    fchkout_ << line << std::endl;
    if (orbitals_.isRestricted()) {
      writeRestrictedOrbitals(line);
    }
    else {
      writeAlphaOrbitals(line);
      writeBetaOrbitals(line);
    }
  }
}

void GaussianOrbitalWriter::writeRestrictedOrbitals(const std::string& line) {
  std::string pattern = "Alpha MO coefficients";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    ignoreInputLines();
    writeCoefficients(orbitals_.restrictedMatrix());
  }
}

void GaussianOrbitalWriter::writeAlphaOrbitals(const std::string& line) {
  std::string pattern = "Alpha MO coefficients";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    ignoreInputLines();
    writeCoefficients(orbitals_.alphaMatrix());
  }
}

void GaussianOrbitalWriter::writeBetaOrbitals(const std::string& line) {
  std::string pattern = "Beta MO coefficients";
  int cmp = line.compare(0, pattern.length(), pattern);
  if (cmp == 0) {
    ignoreInputLines();
    writeCoefficients(orbitals_.betaMatrix());
  }
}

void GaussianOrbitalWriter::ignoreInputLines() {
  unsigned numberCoefficients = orbitals_.numberOrbitals() * orbitals_.numberOrbitals();
  unsigned numberValuesPerLine = 5;
  for (unsigned i = 0; i < numberCoefficients; i += numberValuesPerLine) {
    std::string lineToIgnore;
    std::getline(fchkin_, lineToIgnore);
  }
}

void GaussianOrbitalWriter::writeCoefficients(const Eigen::MatrixXd& c) {
  unsigned numberCoefficients = orbitals_.numberOrbitals() * orbitals_.numberOrbitals();
  Eigen::Map<const Eigen::VectorXd> data(c.data(), numberCoefficients);
  unsigned numberValuesPerLine = 5;
  for (unsigned i = 0; i < numberCoefficients; i += numberValuesPerLine) {
    for (unsigned j = i; j < i + numberValuesPerLine && j < numberCoefficients; ++j) {
      fchkout_ << this->convertToScientificNotation(data[j]);
    }
    fchkout_ << std::endl;
  }
}

std::string GaussianOrbitalWriter::convertToScientificNotation(const double number) const {
  using namespace std;
  std::stringstream ss;
  ss.imbue(std::locale("C"));
  ss << scientific << setprecision(8) << setw(16) << number;
  std::string numberString = ss.str();
  numberString[numberString.size() - 4] = 'E';
  return numberString;
}

void GaussianOrbitalWriter::closeFchkFiles() {
  fchkin_.close();
  fchkout_.close();
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
