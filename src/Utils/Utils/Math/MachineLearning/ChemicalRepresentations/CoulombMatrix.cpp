/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CoulombMatrix.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Utils {
namespace MachineLearning {

CoulombMatrix::CoulombMatrix(const AtomCollection& structure) : nAtoms_(structure.size()) {
  if (nAtoms_ == 0) {
    throw std::runtime_error("Cannot create a Coulomb matrix for an empty molecular structure!");
  }
  featureVector_ = generateCoulombMatrix(structure);
}

Eigen::VectorXd CoulombMatrix::generateCoulombMatrix(const AtomCollection& structure) {
  std::vector<double> features;

  const int N = structure.size();

  for (int i = 0; i < N; ++i) {
    const unsigned z1 = ElementInfo::Z(structure.getElement(i));
    for (int j = i; j < N; ++j) {
      const unsigned z2 = ElementInfo::Z(structure.getElement(j));
      const double value =
          (i == j) ? (0.5 * pow(z1, diagonalElementsExponent_)) : (z1 * z2 / interatomicDistance(i, j, structure));
      features.push_back(value);
    }
  }

  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(features.data(), features.size());
}

double CoulombMatrix::interatomicDistance(int i, int j, const AtomCollection& structure) {
  return (structure.getPosition(i) - structure.getPosition(j)).norm();
}

// Since the full matrix is rarely needed, it is constructed from the feature vector on the fly
Eigen::MatrixXd CoulombMatrix::getMatrix() const {
  Eigen::MatrixXd coulombMatrix(nAtoms_, nAtoms_);
  int counter = 0;
  for (int i = 0; i < nAtoms_; ++i) {
    for (int j = i; j < nAtoms_; ++j) {
      coulombMatrix(i, j) = featureVector_(counter);
      coulombMatrix(j, i) = featureVector_(counter);
      counter++;
    }
  }
  return coulombMatrix;
}

Eigen::VectorXd CoulombMatrix::getFeatures() const {
  return featureVector_;
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
