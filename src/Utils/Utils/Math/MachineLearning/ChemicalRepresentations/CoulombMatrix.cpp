/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CoulombMatrix.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Utils {
namespace MachineLearning {

CoulombMatrix::CoulombMatrix(const AtomCollection& structure) : structure_(structure), nAtoms_(structure.size()) {
  if (nAtoms_ == 0)
    throw std::runtime_error("Cannot create a Coulomb matrix for an empty molecular structure!");
  generateCoulombMatrix();
}

void CoulombMatrix::generateCoulombMatrix() {
  std::vector<double> features;

  for (int i = 0; i < nAtoms_; ++i) {
    int z1 = ElementInfo::Z(structure_.getElement(i));
    for (int j = i; j < nAtoms_; ++j) {
      int z2 = ElementInfo::Z(structure_.getElement(j));
      double value = (i == j) ? (0.5 * pow(z1, diagonalElementsExponent_)) : (z1 * z2 / interatomicDistance(i, j));
      features.push_back(value);
    }
  }

  featureVector_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(features.data(), features.size());
}

double CoulombMatrix::interatomicDistance(int i, int j) {
  auto distanceVector = structure_.getPosition(i) - structure_.getPosition(j);
  return distanceVector.norm();
}

// Since the full matrix is rarely needed, it is constructed from the feature vector on the fly
Eigen::MatrixXd CoulombMatrix::getMatrix() const {
  Eigen::MatrixXd coulombMatrix(nAtoms_, nAtoms_);
  int counter = 0;
  for (int i = 0; i < nAtoms_; ++i) {
    for (int j = i; j < nAtoms_; ++j) {
      coulombMatrix(i, j) = featureVector_(counter);
      counter++;
    }
  }
  return coulombMatrix.selfadjointView<Eigen::Upper>();
}

Eigen::VectorXd CoulombMatrix::getFeatures() const {
  return featureVector_;
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
