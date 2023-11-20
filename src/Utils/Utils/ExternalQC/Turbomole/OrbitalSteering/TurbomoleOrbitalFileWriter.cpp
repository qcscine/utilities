/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TurbomoleOrbitalFileWriter.h"
#include "TurbomoleNumberStringConverter.h"
#include "TurbomoleOrbitalsMetaInformation.h"

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

TurbomoleOrbitalFileWriter::TurbomoleOrbitalFileWriter(const Eigen::MatrixXd& coefficientMatrix,
                                                       const TurbomoleOrbitalsMetaInformation& metaInformation)
  : coefficientMatrix_(coefficientMatrix),
    metaInformation_(metaInformation),
    nOrbitals_(static_cast<unsigned>(coefficientMatrix.rows())) {
}

void TurbomoleOrbitalFileWriter::writeToFile(const std::string& file) {
  moFile_.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  moFile_.open(file);
  writeHeader();
  writeOrbitals();
  writeFooter();
  moFile_.close();
}

void TurbomoleOrbitalFileWriter::writeHeader() {
  moFile_ << metaInformation_.getFirstHeaderLine() << std::endl;
  moFile_ << metaInformation_.getSecondHeaderLine() << std::endl;
  moFile_ << metaInformation_.getThirdHeaderLine() << std::endl;
}

void TurbomoleOrbitalFileWriter::writeOrbitals() {
  for (currentOrbital_ = 0; currentOrbital_ < nOrbitals_; ++currentOrbital_)
    writeOneOrbital();
}

void TurbomoleOrbitalFileWriter::writeOneOrbital() {
  moFile_ << metaInformation_.getOrbitalInformation(currentOrbital_) << std::endl;
  unsigned coefficientIndex_ = 0;
  while (coefficientIndex_ < nOrbitals_) {
    for (unsigned i = 0; i < 4 && coefficientIndex_ < nOrbitals_; ++i) {
      double c = coefficientMatrix_(coefficientIndex_, currentOrbital_);
      moFile_ << TurbomoleNumberStringConverter::toString(c);
      ++coefficientIndex_;
    }
    moFile_ << std::endl;
  }
}

void TurbomoleOrbitalFileWriter::writeFooter() {
  moFile_ << metaInformation_.getFooter() << std::endl;
}

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
