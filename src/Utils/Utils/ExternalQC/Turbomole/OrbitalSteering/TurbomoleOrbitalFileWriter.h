/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef ORBITALPERTURBATION_TURBOMOLEORBITALFILEWRITER_H
#define ORBITALPERTURBATION_TURBOMOLEORBITALFILEWRITER_H

#include <Eigen/Core>
#include <fstream>
namespace Scine {
namespace Utils {
namespace ExternalQC {

namespace TurbomoleOrbitalPerturbation {
class TurbomoleOrbitalsMetaInformation;

/*!
 * Class to write into a single Turbomole molecular orbital file (for instance alpha, beta, mo, ...).
 * For this, it needs a matrix and an instance of TurbomoleOrbitalsMetaInformation.
 */
class TurbomoleOrbitalFileWriter {
 public:
  TurbomoleOrbitalFileWriter(const Eigen::MatrixXd& coefficientMatrix, const TurbomoleOrbitalsMetaInformation& metaInformation);
  void writeToFile(const std::string& file);

 private:
  void writeHeader();
  void writeOrbitals();
  void writeOneOrbital();
  void writeFooter();

  const Eigen::MatrixXd& coefficientMatrix_;
  const TurbomoleOrbitalsMetaInformation& metaInformation_;
  const unsigned nOrbitals_;

  unsigned currentOrbital_ = 0;
  std::ofstream moFile_;
};

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // ORBITALPERTURBATION_TURBOMOLEORBITALFILEWRITER_H
