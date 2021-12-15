/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Properties.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Utils {
namespace Geometry {
namespace Properties {

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

bool operator==(const PositionCollection& lhs, const PositionCollection& rhs) {
  return lhs.isApprox(rhs);
}

bool operator!=(const PositionCollection& lhs, const PositionCollection& rhs) {
  return !(lhs == rhs);
}

} /* namespace Properties */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */
