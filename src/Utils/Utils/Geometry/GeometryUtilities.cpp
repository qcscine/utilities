/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/MolecularTrajectory.h"
#include <Eigen/Eigenvalues>
#include <chrono>
#include <limits>

namespace Scine {
namespace Utils {
namespace Geometry {

PositionCollection translatePositions(const PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& translation) {
  PositionCollection pc = positions;
  translatePositionsInPlace(pc, translation);
  return pc;
}

void translatePositionsInPlace(PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& translation) {
  positions.rowwise() += translation;
}

PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& startOrientation,
                                   const Eigen::Ref<Eigen::RowVector3d>& endOrientation,
                                   const Eigen::Ref<Eigen::RowVector3d>& pointOfRotation) {
  // set up Quaternion
  Eigen::Quaterniond rotationQuat = Eigen::Quaterniond().setFromTwoVectors(startOrientation, endOrientation);

  PositionCollection rotPos = positions;
  rotatePositionsInPlace(rotPos, rotationQuat, pointOfRotation);

  return rotPos;
}

PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& rotAxisOrientation,
                                   double angle, const Eigen::Ref<Eigen::RowVector3d>& pointOfRotation) {
  // set up AngleAxis
  Eigen::AngleAxisd rotAxis(angle, rotAxisOrientation.normalized());
  // set up Quaternion
  Eigen::Quaterniond rotationQuat(rotAxis);

  PositionCollection rotPos = positions;
  rotatePositionsInPlace(rotPos, rotationQuat, pointOfRotation);
  return rotPos;
}

void rotatePositionsInPlace(PositionCollection& pc, const Eigen::Quaterniond& rotation,
                            const Eigen::Ref<Eigen::RowVector3d>& pointOfRotation) {
  Eigen::RowVector3d pointOfRotationToOrigin = pointOfRotation * -1;
  // translate position to origin
  translatePositionsInPlace(pc, pointOfRotationToOrigin);
  for (int i = 0; i < pc.rows(); i++) {
    pc.row(i) = rotation * pc.row(i);
  }
  // translate position back to point of rotation
  translatePositionsInPlace(pc, pointOfRotation);
}

PositionCollection randomDisplacement(const PositionCollection& positions, const double maxDisplacement, const double seed) {
  srand(seed);
  return positions + (maxDisplacement * DisplacementCollection::Random(positions.rows(), positions.cols()));
}

void randomDisplacementInPlace(PositionCollection& positions, const double maxDisplacement, const double seed) {
  srand(seed);
  positions += (maxDisplacement * DisplacementCollection::Random(positions.rows(), positions.cols()));
}

MolecularTrajectory randomDisplacementTrajectory(const AtomCollection& atoms, const int numFrames, const double maxDisplacement) {
  auto result = MolecularTrajectory();
  result.setElementTypes(atoms.getElements());
  for (int i = 0; i < numFrames; ++i) {
    double seed = static_cast<double>(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
    seed *= i + 1;
    result.push_back(randomDisplacement(atoms.getPositions(), maxDisplacement, seed));
  }
  return result;
}

int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          double squaredDistanceConsideredZero) {
  assert(positions.rows() != 0 && "Cannot determine closest atom if there are no atoms!");
  double minimalDistanceSquared = std::numeric_limits<double>::max();
  int closestAtom = 0;

  auto nAtoms = static_cast<int>(positions.rows());

  for (int i = 0; i < nAtoms; i++) {
    double currentSquaredDistance = (positions.row(i) - targetPosition).squaredNorm();
    if (currentSquaredDistance <= squaredDistanceConsideredZero)
      continue;
    if (currentSquaredDistance < minimalDistanceSquared) {
      minimalDistanceSquared = currentSquaredDistance;
      closestAtom = i;
    }
  }

  return closestAtom;
}

int getIndexOfAtomInStructure(const AtomCollection& structure, const Atom& atom, double squaredDistanceConsideredZero) {
  int index = 0;
  ElementType element = atom.getElementType();
  for (const auto& a : structure) {
    if (a.getElementType() == element) { // First check whether the element types match
      double squaredDistance = (a.getPosition() - atom.getPosition()).squaredNorm();
      if (squaredDistance <= squaredDistanceConsideredZero) { // Check whether the atoms have basically the same position
        return index;
      }
    }
    index++;
  }
  throw std::runtime_error("The given atom was not found in the given structure.");
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

void alignPositions(const PositionCollection& reference, PositionCollection& positions, const ElementTypeCollection& elements) {
  QuaternionFit fit(reference, positions, elements);
  positions = fit.getFittedData();
}

std::vector<int> getListOfDivergingAtoms(const PositionCollection& reference, PositionCollection positions,
                                         double threshold, const ElementTypeCollection& elements) {
  std::vector<int> indices;
  indices.reserve(positions.rows());
  elements.empty() ? alignPositions(reference, positions) : alignPositions(reference, positions, elements);
  Eigen::VectorXd difference = (reference - positions).rowwise().norm();
  for (int i = 0; i < difference.size(); ++i) {
    if (difference[i] > threshold)
      indices.push_back(i);
  }
  return indices;
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
  return positions.colwise().sum() / positions.rows();
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

Eigen::MatrixXd calculateTranslationAndRotationModes(const PositionCollection& positions, const ElementTypeCollection& elements) {
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

Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions,
                                                     const ElementTypeCollection& elements, bool massWeighted) {
  auto rotoTranslation = calculateTranslationAndRotationModes(positions, elements);
  // If mass-weighted Hessian shall be transformed the rotation and translation modes have to be adapted accordingly
  if (massWeighted) {
    int nAtoms = elements.size();
    auto masses = Geometry::getMasses(elements);
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

} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */
