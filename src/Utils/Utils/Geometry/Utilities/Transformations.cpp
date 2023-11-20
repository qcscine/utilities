/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Transformations.h"
#include "Properties.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Math/IterativeDiagonalizer/SubspaceOrthogonalizer.h"

namespace Scine {
namespace Utils {
namespace Geometry {
namespace Transformations {

Eigen::MatrixXd positionVectorToMatrix(const Eigen::VectorXd& v) {
  assert(v.size() % 3 == 0);

  using namespace Eigen;
  using RowMajorMatrix = Matrix<double, Dynamic, Dynamic, RowMajor>;

  auto numberParticles = v.size() / 3;
  Map<const RowMajorMatrix> mapMatrix(v.data(), numberParticles, 3);
  return mapMatrix;
}

Eigen::VectorXd positionMatrixToVector(const Eigen::MatrixXd& m) {
  assert(m.cols() == 3);

  using namespace Eigen;
  using RowMajorMatrix = Matrix<double, Dynamic, Dynamic, RowMajor>;

  RowMajorMatrix m2(m);
  Map<const VectorXd> mapVector(m2.data(), m2.size());

  return mapVector;
}

Eigen::MatrixXd calculateTranslationAndRotationModes(const PositionCollection& positions, const ElementTypeCollection& elements) {
  auto nAtoms = positions.rows();
  auto masses = Properties::getMasses(elements);
  auto centerOfMass = Properties::getCenterOfMass(positions, masses);
  auto principalMoments = Properties::calculatePrincipalMoments(positions, masses, centerOfMass);

  const auto& X = principalMoments.eigenvectors;

  Eigen::MatrixXd allRotoTranslationVectors(3 * nAtoms, 6);

  for (int i = 0; i < nAtoms; ++i) {
    // Translation
    allRotoTranslationVectors.block(3 * i, 0, 3, 3) = Eigen::MatrixXd::Identity(3, 3);

    // Rotation
    Eigen::Vector3d P = X.transpose() * (positions.row(i) - centerOfMass).transpose();
    Eigen::Matrix3d d;
    d.col(0) = P.y() * X.col(2) - P.z() * X.col(1);
    d.col(1) = P.z() * X.col(0) - P.x() * X.col(2);
    d.col(2) = P.x() * X.col(1) - P.y() * X.col(0);

    allRotoTranslationVectors.block(3 * i, 3, 3, 3) = d;
  }

  // look which roto-translational modes are valid
  std::vector<int> validRotoTranslationModes;
  for (int i = 0; i < 6; ++i) {
    // The norm will be zero for "invalid" roto-translation modes (single atoms, linear molecules).
    auto squaredNorm = allRotoTranslationVectors.col(i).squaredNorm();
    if (squaredNorm > 0.1) {
      validRotoTranslationModes.push_back(i);
    }
  }

  // return only the valid modes
  auto numberRotoTranslationModes = static_cast<int>(validRotoTranslationModes.size());
  Eigen::MatrixXd rotoTranslationVectors(3 * nAtoms, numberRotoTranslationModes);
  for (int i = 0; i < numberRotoTranslationModes; ++i) {
    rotoTranslationVectors.col(i) = allRotoTranslationVectors.col(validRotoTranslationModes[i]);
  }
  rotoTranslationVectors.colwise().normalize();

  return rotoTranslationVectors;
}

Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions,
                                                     const ElementTypeCollection& elements, bool massWeighted) {
  auto rotoTranslation = calculateTranslationAndRotationModes(positions, elements);
  // If mass-weighted Hessian shall be transformed the rotation and translation modes have to be adapted accordingly
  if (massWeighted) {
    int nAtoms = elements.size();
    auto masses = Properties::getMasses(elements);
    for (int i = 0; i < nAtoms; ++i) {
      rotoTranslation.middleRows(3 * i, 3) *= std::sqrt(masses[i]);
    }
    rotoTranslation.colwise().normalize();
  }
  auto nDims = rotoTranslation.rows();
  auto numberRotoTranslationModes = rotoTranslation.cols();

  srand(42);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(nDims, nDims);

  A.leftCols(numberRotoTranslationModes) = rotoTranslation;

  // Orthogonalization
  auto n = nDims;
  auto m = nDims;
  Eigen::MatrixXd R_1(n, n);
  Eigen::MatrixXd Q_1(m, n);
  for (int i = 0; i < n; ++i) {
    R_1(i, i) = A.col(i).norm();
    Q_1.col(i) = A.col(i) / R_1(i, i);
    for (int j = i + 1; j < n; ++j) {
      R_1(i, j) = Q_1.col(i).transpose() * A.col(j);
      A.col(j) -= (Q_1.col(i) * R_1(i, j));
    }
  }

  A.colwise().normalize();

  return A.rightCols(nDims - numberRotoTranslationModes);
}

Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions, const ElementTypeCollection& elements,
                                                     const GradientCollection& gradients, bool massWeighted) {
  auto rotoTranslation = calculateTranslationAndRotationModes(positions, elements);
  rotoTranslation.conservativeResize(Eigen::NoChange, rotoTranslation.cols() + 1);
  rotoTranslation.col(rotoTranslation.cols() - 1) = Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.size());
  // If mass-weighted Hessian shall be transformed the rotation and translation modes have to be adapted accordingly
  if (massWeighted) {
    int nAtoms = elements.size();
    auto masses = Properties::getMasses(elements);
    for (int i = 0; i < nAtoms; ++i) {
      rotoTranslation.middleRows(3 * i, 3) *= std::sqrt(masses[i]);
    }
    rotoTranslation.colwise().normalize();
  }
  // Orthogonalize gradient to other components
  SubspaceOrthogonalizer::orthogonalizeToSubspace(rotoTranslation.col(rotoTranslation.cols() - 1),
                                                  rotoTranslation.leftCols(rotoTranslation.cols() - 1),
                                                  rotoTranslation.col(rotoTranslation.cols() - 1));
  auto nDims = rotoTranslation.rows();
  auto numberRotoTranslationModes = rotoTranslation.cols();

  srand(42);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(nDims, nDims);

  A.leftCols(numberRotoTranslationModes) = rotoTranslation;

  // Orthogonalization
  auto n = nDims;
  auto m = nDims;
  Eigen::MatrixXd R_1(n, n);
  Eigen::MatrixXd Q_1(m, n);
  for (int i = 0; i < n; ++i) {
    R_1(i, i) = A.col(i).norm();
    Q_1.col(i) = A.col(i) / R_1(i, i);
    for (int j = i + 1; j < n; ++j) {
      R_1(i, j) = Q_1.col(i).transpose() * A.col(j);
      A.col(j) -= (Q_1.col(i) * R_1(i, j));
    }
  }

  A.colwise().normalize();

  return A.rightCols(nDims - numberRotoTranslationModes);
}

} /* namespace Transformations */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */
