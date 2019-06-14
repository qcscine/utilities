/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Math/QuaternionFit.h"
#include <Eigen/Eigenvalues>
#include <limits>

namespace Scine {
namespace Utils {
namespace Geometry {

PositionCollection translatePositions(const PositionCollection& positions, const Displacement& translation) {
  PositionCollection pc = positions;
  translatePositions(pc, translation);
  return pc;
}

void translatePositions(PositionCollection& positions, const Displacement& translation) {
  positions.rowwise() += translation;
}

unsigned getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition) {
  assert(positions.rows() == 0 && "Cannot determine closest atom if there are no atoms!");
  double minimalDistanceSquared = std::numeric_limits<double>::max();
  unsigned closestAtom = 0;

  auto nAtoms = static_cast<unsigned int>(positions.rows());

  for (unsigned int i = 0; i < nAtoms; i++) {
    double currentSquaredDistance = (positions.row(i) - targetPosition).squaredNorm();
    if (currentSquaredDistance < minimalDistanceSquared) {
      minimalDistanceSquared = currentSquaredDistance;
      closestAtom = i;
    }
  }

  return closestAtom;
}

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

void alignPositions(const PositionCollection& reference, PositionCollection& positions) {
  QuaternionFit fit(reference, positions);
  positions = fit.getFittedData();
}

std::vector<double> getMasses(const ElementTypeCollection& elements) {
  std::vector<double> masses;
  masses.reserve(elements.size());
  std::transform(elements.begin(), elements.end(), std::back_inserter(masses), ElementInfo::mass);
  return masses;
}

Position getCenterOfMass(const PositionCollection& positions, const std::vector<double>& masses) {
  Position P;
  P.setZero();
  double totalMass = 0;
  for (int i = 0; i < positions.rows(); ++i) {
    P += masses[i] * positions.row(i);
    totalMass += masses[i];
  }
  P /= totalMass;
  return P;
}

Position getCenterOfMass(const AtomCollection& structure) {
  auto masses = getMasses(structure.getElements());
  return getCenterOfMass(structure.getPositions(), masses);
}

Position getAveragePosition(const PositionCollection& positions) {
  return positions.rowwise().sum() / positions.rows();
}

Eigen::Matrix3d calculateInertiaTensor(const PositionCollection& positions, const std::vector<double>& masses,
                                       const Position& centerOfMass) {
  double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
  for (int i = 0; i < positions.rows(); ++i) {
    double m = masses[i];
    auto dp = positions.row(i) - centerOfMass;
    double x = dp.x();
    double y = dp.y();
    double z = dp.z();
    Ixx += m * (y * y + z * z);
    Iyy += m * (x * x + z * z);
    Izz += m * (x * x + y * y);
    Ixy -= m * x * y;
    Ixz -= m * x * z;
    Iyz -= m * y * z;
  }
  Eigen::Matrix3d In;
  In << Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz;
  return In;
}

PrincipalMomentsOfInertia calculatePrincipalMoments(const PositionCollection& positions,
                                                    const std::vector<double>& masses, const Position& centerOfMass) {
  auto In = calculateInertiaTensor(positions, masses, centerOfMass);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.compute(In);

  PrincipalMomentsOfInertia pmi;
  pmi.eigenvalues = es.eigenvalues();
  pmi.eigenvectors = es.eigenvectors();
  return pmi;
}

Eigen::MatrixXd calculateRotationAndTranslationModes(const PositionCollection& positions, const ElementTypeCollection& elements) {
  auto nAtoms = positions.rows();
  auto masses = Utils::Geometry::getMasses(elements);
  auto centerOfMass = Utils::Geometry::getCenterOfMass(positions, masses);
  auto principalMoments = Utils::Geometry::calculatePrincipalMoments(positions, masses, centerOfMass);

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

Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions, const ElementTypeCollection& elements) {
  auto rotoTranslation = calculateRotationAndTranslationModes(positions, elements);
  auto nDims = rotoTranslation.rows();
  auto numberRotoTranslationModes = rotoTranslation.cols();

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

} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */
