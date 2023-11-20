/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TurbomoleOrbitalFileReader.h"
#include "TurbomoleNumberStringConverter.h"

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

TurbomoleOrbitalFileReader::TurbomoleOrbitalFileReader(const std::string& file, unsigned nOrbitals)
  : nOrbitals_(nOrbitals) {
  coefficientMatrix_.resize(nOrbitals, nOrbitals);
  extractContent(file);
  readOrbitals();
}

void TurbomoleOrbitalFileReader::extractContent(const std::string& file) {
  moFile_.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  moFile_.open(file);
  // TODO: This file is never properly closed. Moreover, this function is named misleadingly;
  // it doesn't extract any content.
}

void TurbomoleOrbitalFileReader::readOrbitals() {
  readHeader();
  for (currentOrbital_ = 0; currentOrbital_ < nOrbitals_; ++currentOrbital_)
    readOneMolecularOrbital();
  readFooter();
}

void TurbomoleOrbitalFileReader::readOneMolecularOrbital() {
  std::string metaInformationString;
  getline(moFile_, metaInformationString);
  metaInformation_.addOrbitalInformation(std::move(metaInformationString));

  unsigned coefficientIndex = 0;
  while (coefficientIndex < nOrbitals_) {
    std::string coefficientsString;
    getline(moFile_, coefficientsString);
    unsigned linePosition = 0;
    while (linePosition < coefficientsString.length()) {
      std::string numberString{coefficientsString, linePosition, 20};
      double c = TurbomoleNumberStringConverter::toDouble(numberString);
      coefficientMatrix_(coefficientIndex, currentOrbital_) = c;
      linePosition += 20;
      ++coefficientIndex;
    }
  }
}

const Eigen::MatrixXd& TurbomoleOrbitalFileReader::getCoefficientMatrix() const {
  return coefficientMatrix_;
}

const TurbomoleOrbitalsMetaInformation& TurbomoleOrbitalFileReader::getMetaInformation() const {
  return metaInformation_;
}

void TurbomoleOrbitalFileReader::readHeader() {
  std::string header1, header2, header3;
  getline(moFile_, header1);
  getline(moFile_, header2);
  getline(moFile_, header3);
  metaInformation_.setHeader(std::move(header1), std::move(header2), std::move(header3));
}

void TurbomoleOrbitalFileReader::readFooter() {
  std::string footer;
  getline(moFile_, footer);
  metaInformation_.setFooter(std::move(footer));
}

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
