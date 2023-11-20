/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Manipulations.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/MolecularTrajectory.h"
#include <chrono>

namespace Scine {
namespace Utils {
namespace Geometry {
namespace Manipulations {

PositionCollection translatePositions(const PositionCollection& positions, const Eigen::RowVector3d& translation) {
  PositionCollection pc = positions;
  translatePositionsInPlace(pc, translation);
  return pc;
}

void translatePositionsInPlace(PositionCollection& positions, const Eigen::RowVector3d& translation) {
  positions.rowwise() += translation;
}

PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::RowVector3d& startOrientation,
                                   const Eigen::RowVector3d& endOrientation, const Eigen::RowVector3d& pointOfRotation) {
  // set up Quaternion
  Eigen::Quaterniond rotationQuat = Eigen::Quaterniond().setFromTwoVectors(startOrientation, endOrientation);

  PositionCollection rotPos = positions;
  rotatePositionsInPlace(rotPos, rotationQuat, pointOfRotation);

  return rotPos;
}

PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::RowVector3d& rotAxisOrientation,
                                   double angle, const Eigen::RowVector3d& pointOfRotation) {
  // set up AngleAxis
  Eigen::AngleAxisd rotAxis(angle, rotAxisOrientation.normalized());
  // set up Quaternion
  Eigen::Quaterniond rotationQuat(rotAxis);

  PositionCollection rotPos = positions;
  rotatePositionsInPlace(rotPos, rotationQuat, pointOfRotation);
  return rotPos;
}

void rotatePositionsInPlace(PositionCollection& pc, const Eigen::Quaterniond& rotation,
                            const Eigen::RowVector3d& pointOfRotation) {
  Eigen::RowVector3d pointOfRotationToOrigin = pointOfRotation * -1;
  // translate position to origin
  translatePositionsInPlace(pc, pointOfRotationToOrigin);
  for (int i = 0; i < pc.rows(); i++) {
    pc.row(i) = rotation * pc.row(i);
  }
  // translate position back to point of rotation
  translatePositionsInPlace(pc, pointOfRotation);
}

PositionCollection displaceAlongModes(const PositionCollection& positions, const std::vector<NormalMode>& modes,
                                      const std::vector<double>& displacementSteps) {
  PositionCollection modifiedPositions = positions;
  displaceAlongModesInPlace(modifiedPositions, modes, displacementSteps);
  return modifiedPositions;
}

void displaceAlongModesInPlace(PositionCollection& positions, const std::vector<NormalMode>& modes,
                               const std::vector<double>& displacementSteps) {
  assert(modes.size() == displacementSteps.size());
  for (std::size_t iMode = 0; iMode < modes.size(); ++iMode) {
    DisplacementCollection modeToDisplace = modes[iMode].getMode();
    modeToDisplace *= displacementSteps[iMode];
    positions += modeToDisplace;
  }
}

PositionCollection randomDisplacement(const PositionCollection& positions, double maxDisplacement) {
  return positions + (maxDisplacement * DisplacementCollection::Random(positions.rows(), positions.cols()));
}

void randomDisplacementInPlace(PositionCollection& positions, double maxDisplacement) {
  positions += (maxDisplacement * DisplacementCollection::Random(positions.rows(), positions.cols()));
}

MolecularTrajectory randomDisplacementTrajectory(const AtomCollection& atoms, unsigned int numFrames, double maxDisplacement) {
  auto result = MolecularTrajectory();
  result.setElementTypes(atoms.getElements());
  const PositionCollection& positions = atoms.getPositions();
  for (unsigned int i = 0; i < numFrames; ++i) {
    result.push_back(randomDisplacement(positions, maxDisplacement));
  }
  return result;
}

MolecularTrajectory randomDisplacementTrajectory(const AtomCollection& atoms, unsigned int numFrames,
                                                 double maxDisplacement, unsigned int seed) {
  srand(seed);
  return randomDisplacementTrajectory(atoms, numFrames, maxDisplacement);
}

void alignPositions(const PositionCollection& reference, PositionCollection& positions) {
  QuaternionFit fit(reference, positions);
  positions = fit.getFittedData();
}

void alignPositions(const PositionCollection& reference, PositionCollection& positions, const ElementTypeCollection& elements) {
  QuaternionFit fit(reference, positions, elements);
  positions = fit.getFittedData();
}

void alignPositions(const PositionCollection& reference, PositionCollection& positions, const Eigen::VectorXd& weights,
                    Eigen::VectorXd& rmsdVector) {
  QuaternionFit fit(reference, positions, weights);
  positions = fit.getFittedData();
  rmsdVector = (positions - reference).rowwise().squaredNorm();
  rmsdVector.array() *= weights.normalized().array();
  rmsdVector = rmsdVector.cwiseSqrt();
}

} /* namespace Manipulations */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */
