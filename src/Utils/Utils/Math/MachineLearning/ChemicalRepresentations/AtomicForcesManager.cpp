/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomicForcesManager.h"
#include "../PrincipalComponentAnalysis.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Utils {
namespace MachineLearning {

AtomicForcesManager::AtomicForcesManager(const AtomCollection& structure)
  : structure_(structure), currentPositions_(structure.getPositions()), nAtoms_(structure.size()) {
  if (nAtoms_ < 4) {
    throw std::runtime_error("The number of atoms in the molecular system must be larger than 3 in order to generate "
                             "internal force representations.");
  }
  transformationMatrices_.resize(nAtoms_);
  // Which atoms are close enough to be considered important, is determined only once for the original structure
  determineImportantAtoms();
  // Construct first set of transformation matrices for the original structure
  constructTransformationMatrices();
}

void AtomicForcesManager::modifyPositions(const PositionCollection& newPositions) {
  currentPositions_ = newPositions;
  constructTransformationMatrices(); // update the transformation matrices
}

void AtomicForcesManager::constructTransformationMatrices() {
  // Calculate center of charge in order to later decide on the sign of the PCA axes
  auto centerOfCharge = calculateCenterOfCharge();
  // Loop over all atoms and calculate their transformation matrices
#pragma omp parallel for
  for (int i = 0; i < nAtoms_; ++i) {
    Eigen::MatrixXd chemicalEnvironment = calculateChemicalEnvironment(i);
    MachineLearning::PrincipalComponentAnalysis pca(chemicalEnvironment);
    Eigen::MatrixXd matrix = pca.calculate(3).first;

    // Decide on the sign of the PCA axes
    // Convention is: Center of charge only has positive coordinates
    Eigen::MatrixXd inverseMatrix = matrix.inverse();
    Eigen::Vector3d internalCenterOfCharge = inverseMatrix * (centerOfCharge - currentPositions_.row(i)).transpose();
    for (int d = 0; d < 3; ++d) {
      if (internalCenterOfCharge(d) < 0.0) {
        matrix.col(d) *= -1;
      }
    }

    transformationMatrices_[i] = matrix;
  }
}

Eigen::VectorXd AtomicForcesManager::calculateFeatures(int atomIndex) {
  Eigen::MatrixXd fMat = calculateFeatureMatrix(atomIndex);
  fMat.transposeInPlace();
  Eigen::VectorXd fVec(Eigen::Map<Eigen::VectorXd>(fMat.data(), fMat.size()));
  return fVec;
}

Eigen::MatrixXd AtomicForcesManager::calculateFeatureMatrix(int atomIndex) {
  // Custom container that holds the unsorted rows of the feature matrix along with the interatomic distances
  std::vector<std::pair<Eigen::RowVector3d, double>> unsortedFeatures;
  // Loop over all atoms and fill up this matrix
  Eigen::RowVector3d pos = currentPositions_.row(atomIndex);
  for (int i : importantAtoms_.at(atomIndex)) {
    Eigen::RowVector3d row = currentPositions_.row(i) - pos;
    auto rowTransformed = toInternalRepresentation(row, atomIndex);
    auto distance = rowTransformed.norm();
    assert(std::abs(rowTransformed.norm() - row.norm()) < 1e-6);

    rowTransformed *= ElementInfo::Z(structure_.getElement(i)) / std::pow(distance, 3);
    // Add this row to the unsorted features container together with the calculated distance
    unsortedFeatures.emplace_back(std::make_pair(rowTransformed, distance));
  }
  return sortFeatureMatrix(unsortedFeatures);
}

Eigen::MatrixXd AtomicForcesManager::calculateChemicalEnvironment(int atomIndex) {
  // Initialize an empty environment matrix
  Eigen::MatrixXd chemEnv(0, 3);
  // Loop over all atoms and fill up this matrix
  Eigen::RowVector3d pos = currentPositions_.row(atomIndex);
  for (int i : importantAtoms_.at(atomIndex)) {
    Eigen::RowVector3d row = currentPositions_.row(i) - pos;
    auto distance = row.norm();
    row *= ElementInfo::Z(structure_.getElement(i)) / std::pow(distance, 3);
    // Add this row to the end of the chemical environment matrix for the atom with index 'atomIndex'
    chemEnv.conservativeResize(chemEnv.rows() + 1, chemEnv.cols());
    chemEnv.row(chemEnv.rows() - 1) = row;
  }

  if (chemEnv.size() < 3) {
    throw std::runtime_error("The atom with the index " + std::to_string(atomIndex) +
                             " is too far away from the other atoms. The feature vector could not be generated for it, "
                             "because less than 3 atoms are close enough.");
  }

  return chemEnv;
}

Eigen::RowVector3d AtomicForcesManager::toGlobalRepresentation(const Eigen::RowVector3d& force, int atomIndex) {
  return transformationMatrices_.at(atomIndex) * force.transpose();
}

Eigen::RowVector3d AtomicForcesManager::toInternalRepresentation(const Eigen::RowVector3d& force, int atomIndex) {
  return transformationMatrices_.at(atomIndex).inverse() * force.transpose();
}

Eigen::RowVector3d AtomicForcesManager::calculateCenterOfCharge() {
  Eigen::RowVector3d centerOfCharge(0.0, 0.0, 0.0);
  double totalCharge = 0.0;
  for (int i = 0; i < nAtoms_; ++i) {
    auto z = ElementInfo::Z(structure_.getElement(i));
    centerOfCharge += z * currentPositions_.row(i);
    totalCharge += z;
  }
  return centerOfCharge / totalCharge;
}

Eigen::MatrixXd AtomicForcesManager::sortFeatureMatrix(std::vector<std::pair<Eigen::RowVector3d, double>>& unsortedFeatures) {
  Eigen::MatrixXd fMat(unsortedFeatures.size(), 3);
  std::sort(unsortedFeatures.begin(), unsortedFeatures.end(),
            [](const std::pair<Eigen::RowVector3d, double>& lhs, const std::pair<Eigen::RowVector3d, double>& rhs) {
              return lhs.second < rhs.second;
            });

  const int featuresSize = unsortedFeatures.size();
  for (int i = 0; i < featuresSize; ++i) {
    fMat.row(i) = unsortedFeatures[i].first;
  }

  return fMat;
}

void AtomicForcesManager::determineImportantAtoms() {
  importantAtoms_.resize(nAtoms_);
  for (int i = 0; i < nAtoms_; ++i) {
    for (int j = i + 1; j < nAtoms_; ++j) {
      auto distance = (currentPositions_.row(i) - currentPositions_.row(j)).norm();
      if (distance <= distanceThreshold_) {
        importantAtoms_[i].push_back(j);
        importantAtoms_[j].push_back(i);
      }
    }
  }
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
