
/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef ORBITALPERTURBATION_TURBOMOLEORBITALFILEREADER_H
#define ORBITALPERTURBATION_TURBOMOLEORBITALFILEREADER_H

#include "TurbomoleOrbitalsMetaInformation.h"
#include <Eigen/Core>
#include <fstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

/*!
 * Class to read a single Turbomole molecular orbital file (for instance alpha, beta, mo, ...).
 * The contents of the file are then available as a matrix and an instance of TurbomoleOrbitalsMetaInformation.
 */
class TurbomoleOrbitalFileReader {
 public:
  TurbomoleOrbitalFileReader(const std::string& file, unsigned nOrbitals);

  const Eigen::MatrixXd& getCoefficientMatrix() const;
  const TurbomoleOrbitalsMetaInformation& getMetaInformation() const;

 private:
  void extractContent(const std::string& file);
  void readOrbitals();
  void readOneMolecularOrbital();
  void readHeader();
  void readFooter();

  const unsigned nOrbitals_;
  unsigned currentOrbital_ = 0;
  std::ifstream moFile_;
  Eigen::MatrixXd coefficientMatrix_;
  TurbomoleOrbitalsMetaInformation metaInformation_;
};

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // ORBITALPERTURBATION_TURBOMOLEORBITALFILEREADER_H