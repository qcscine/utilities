/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PeriodicBoundaries.h"
#include "Utils/Technical/SpgInterface.h"
#include <utility>

namespace Scine {
namespace Utils {

PeriodicBoundaries::PeriodicBoundaries(double cubeLength, const std::string& periodicity)
  : PeriodicBoundaries(Eigen::Matrix3d::Identity() * cubeLength, periodicity) {
}

PeriodicBoundaries::PeriodicBoundaries(const PeriodicBoundaries& rhs)
  : PeriodicBoundaries(rhs.getCellMatrix(), rhs.getPeriodicityString()) {
}

PeriodicBoundaries::PeriodicBoundaries(Eigen::Matrix3d matrix, const std::string& periodicity)
  : _matrix(std::move(matrix)) {
  setMembers();
  setPeriodicity(periodicity);
}

PeriodicBoundaries::PeriodicBoundaries(const Eigen::Vector3d& lengths, const Eigen::Vector3d& angles, bool isBohr,
                                       bool isDegrees, const std::string& periodicity) {
  constructMembersFromLengthsAndAngles(lengths, angles, isBohr, isDegrees, periodicity);
}

PeriodicBoundaries::PeriodicBoundaries(std::string periodicBoundariesString, const std::string& delimiter, bool isBohr,
                                       bool isDegrees) {
  std::string input = periodicBoundariesString;
  // split string into vector
  size_t pos = 0;
  std::vector<std::string> entries;
  while ((pos = periodicBoundariesString.find(delimiter)) != std::string::npos) {
    entries.push_back(periodicBoundariesString.substr(0, pos));
    periodicBoundariesString.erase(0, pos + delimiter.length());
  }
  entries.push_back(periodicBoundariesString); // fill in remaining string
  std::string periodicity = "xyz";
  // sanity check
  if (entries.size() != 6) {
    if (entries.size() == 7) {
      periodicity = entries.back();
    }
    else {
      throw std::logic_error("The given string '" + input + "' does not contain 3 cell lengths and 3 cell angles.");
    }
  }
  // construct Eigen vectors
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  for (int i = 0; i < 3; ++i) {
    lengths[i] = std::stod(entries[i]);
    angles[i] = std::stod(entries[i + 3]);
  }
  // construct all class variables
  constructMembersFromLengthsAndAngles(lengths, angles, isBohr, isDegrees, periodicity);
}

void PeriodicBoundaries::constructMembersFromLengthsAndAngles(const Eigen::Vector3d& lengths, const Eigen::Vector3d& angles,
                                                              bool isBohr, bool isDegrees, const std::string& periodicity) {
  const double an = (isBohr) ? lengths[0] : lengths[0] * Constants::bohr_per_angstrom;
  const double bn = (isBohr) ? lengths[1] : lengths[1] * Constants::bohr_per_angstrom;
  const double cn = (isBohr) ? lengths[2] : lengths[2] * Constants::bohr_per_angstrom;
  const double alpha = (isDegrees) ? getRadians(angles[0]) : angles[0];
  const double beta = (isDegrees) ? getRadians(angles[1]) : angles[1];
  const double gamma = (isDegrees) ? getRadians(angles[2]) : angles[2];
  // basis vectors of cell in standard basis, b is init later directly via rotation
  Eigen::Vector3d a = Eigen::Vector3d::Zero();
  Eigen::Vector3d c = Eigen::Vector3d::Zero();
  // convention that a is along x-axis
  a[0] = an;
  // convention that b is in x-y axis --> rotate by gamma along z-axis
  Eigen::Matrix3d rotM;
  // clang-format off
  rotM << std::cos(gamma), -std::sin(gamma), 0.0,
          std::sin(gamma),  std::cos(gamma), 0.0,
          0.0, 0.0, 1.0;
  // clang-format on
  Eigen::Vector3d b = rotM * a;
  // scale according to norms
  b *= (bn / an);
  reduceNoise(b);
  // c0 is known because a only has entry along x axis
  c[0] = std::cos(beta) * cn;
  // for c1 we use dot product with b, because b is defined in x-y plane -> b[2] = 0
  c[1] = bn * cn * std::cos(alpha) - b[0] * c[0];
  c[1] /= b[1];
  // for c2 we can use norm
  c[2] = std::sqrt(std::pow(cn, 2) - std::pow(c[0], 2) - std::pow(c[1], 2));
  // assert that angles of created vectors are identical to input angles
  assert(std::acos(a.dot(b) / (a.norm() * b.norm())) - gamma < _eps);
  assert(std::acos(c.dot(b) / (c.norm() * b.norm())) - alpha < _eps);
  assert(std::acos(a.dot(c) / (a.norm() * c.norm())) - beta < _eps);
  _matrix = Eigen::Matrix3d::Zero();
  _matrix.row(0) = a;
  _matrix.row(1) = b;
  _matrix.row(2) = c;
  setMembers();
  setPeriodicity(periodicity);
}

void PeriodicBoundaries::reduceNoise(Eigen::Vector3d& vec) const {
  for (int i = 0; i < 3; ++i) {
    if (std::fabs(vec[i] - 0.0) < _eps) {
      vec[i] = 0.0;
    }
  }
}

void PeriodicBoundaries::reduceNoise(Eigen::Matrix3d& mat) const {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (std::fabs(mat(i, j) - 0.0) < _eps) {
        mat(i, j) = 0.0;
      }
    }
  }
}

void PeriodicBoundaries::setMembers() {
  reduceNoise(_matrix);
  auto invalid = [&](const Eigen::Matrix3d& matrix) {
    return (matrix(0, 0) < 0 || matrix(1, 1) < 0 || matrix(2, 2) < 0);
  };
  if (invalid(_matrix)) {
    // given matrix is invalid, see if we find an symmetry identical rotation that gives us a valid matrix
    bool found = false;
    auto alternatives = SpgInterface::findAlternativePbcs(_matrix, _eps);
    for (auto& alternative : alternatives) {
      reduceNoise(alternative);
      if (!invalid(alternative)) {
        found = true;
        _matrix = alternative;
        break;
      }
    }
    if (!found) {
      std::stringstream error;
      error << "Periodic boundaries received unphysical information where a unit vector extends "
               "into the opposite direction to where it should extend\n"
            << _matrix;
      throw std::runtime_error(error.str());
    }
  }
  /* vectors */
  const Eigen::Vector3d a = _matrix.row(0);
  const Eigen::Vector3d b = _matrix.row(1);
  const Eigen::Vector3d c = _matrix.row(2);
  /* lengths */
  _aNorm = a.norm();
  _bNorm = b.norm();
  _cNorm = c.norm();
  /* angles */
  _alpha = getDegrees(std::acos(c.dot(b) / (c.norm() * b.norm())));
  _beta = getDegrees(std::acos(a.dot(c) / (a.norm() * c.norm())));
  _gamma = getDegrees(std::acos(a.dot(b) / (a.norm() * b.norm())));
  /* biggest distance for sanity checks */
  Position middle;
  middle << 0.5, 0.5, 0.5;
  _biggestDistanceSquared = transform(middle).squaredNorm();
  /* used as cutoff to know if faster algorithm works */
  Eigen::Vector3d bXc = b.cross(c);
  Eigen::Vector3d cXa = c.cross(a);
  Eigen::Vector3d aXb = a.cross(b);
  std::vector<double> w;
  w.push_back(a.dot(bXc) / bXc.norm());
  w.push_back(b.dot(cXa) / cXa.norm());
  w.push_back(c.dot(aXb) / aXb.norm());
  _minPerpendicularSq = std::pow(*std::min_element(w.begin(), w.end()), 2);
  _invMatrix = _matrix.inverse();
}

PositionCollection PeriodicBoundaries::transform(const PositionCollection& positions, bool relativeToCartesian) const {
  if (relativeToCartesian) {
    return positions * _matrix;
  }
  return positions * _invMatrix;
}

Position PeriodicBoundaries::transform(const Position& position, bool relativeToCartesian) const {
  if (relativeToCartesian) {
    return position * _matrix;
  }
  return position * _invMatrix;
}

Position PeriodicBoundaries::translatePositionsIntoCell(const Position& position, const Eigen::RowVector3d& relShift) const {
  Position shiftedPosition = position;
  translatePositionsIntoCellInPlace(shiftedPosition, relShift);
  return shiftedPosition;
}

PositionCollection PeriodicBoundaries::translatePositionsIntoCell(const PositionCollection& positions,
                                                                  const Eigen::RowVector3d& relShift) const {
  PositionCollection shiftedPositions = positions;
  translatePositionsIntoCellInPlace(shiftedPositions, relShift);
  return shiftedPositions;
}

void PeriodicBoundaries::translatePositionsIntoCellInPlace(Eigen::Ref<Position> position,
                                                           const Eigen::RowVector3d& relShift) const {
  transformInPlace(position, false);
  for (int i = 0; i < 3; ++i) {
    if (_periodicity[i]) {
      // floor is 0 for [0,1), negative no problem, because - and -  leads to right addition
      position[i] -= floor(position[i]);
    }
  }
  position += relShift;
  transformInPlace(position);
}

void PeriodicBoundaries::translatePositionsIntoCellInPlace(PositionCollection& positions,
                                                           const Eigen::RowVector3d& relShift) const {
  for (int i = 0; i < positions.rows(); ++i) {
    Position pos = positions.row(i);
    translatePositionsIntoCellInPlace(pos, relShift);
    positions.row(i) = pos;
  }
}

double PeriodicBoundaries::fastMinimumImageDistanceSquared(const Position& p1, const Position& p2) const {
  const double& x = _aNorm;
  const double& y = _bNorm;
  const double& z = _cNorm;
  double distX = std::fabs(p1[0] - p2[0]);
  double distY = std::fabs(p1[1] - p2[1]);
  double distZ = std::fabs(p1[2] - p2[2]);
  // minimum image convention
  if (_periodicity[0] && distX > x * 0.5) {
    distX -= x;
  }
  if (_periodicity[1] && distY > y * 0.5) {
    distY -= y;
    distX -= _matrix(1, 0); // subtract b0 from distX for non-orthorhombic cells
  }
  if (_periodicity[2] && distZ > z * 0.5) {
    distZ -= z;
    distX -= _matrix(2, 0); // subtract c0 from distX for non-orthorhombic cells
    distY -= _matrix(2, 1); // subtract c1 from distY for non-orthorhombic cells
  }
  return std::pow(distX, 2) + std::pow(distY, 2) + std::pow(distZ, 2);
}

double PeriodicBoundaries::bruteForceMinimumImageDistanceSquared(const Position& p1, const Position& p2) const {
  std::vector<double> distances = getAllImageDistancesSquared(p1, p2);
  return *std::min_element(distances.begin(), distances.end());
}

Displacement PeriodicBoundaries::bruteForceMinimumImageDisplacementVector(const Position& p1, const Position& p2) const {
  std::vector<Displacement> displacements = getAllImageDisplacementVectors(p1, p2);
  std::vector<double> distances;
  for (auto const& disp : displacements) {
    distances.push_back(disp.squaredNorm());
  }
  int minIndex = std::min_element(distances.begin(), distances.end()) - distances.begin();
  return displacements[minIndex];
}

std::vector<double> PeriodicBoundaries::getAllImageDistancesSquared(const Position& p1, const Position& p2) const {
  std::vector<Displacement> displacements = getAllImageDisplacementVectors(p1, p2);
  std::vector<double> result;
  for (auto const& disp : displacements) {
    result.push_back(disp.squaredNorm());
  }
  return result;
}

std::vector<Displacement> PeriodicBoundaries::getAllImageDisplacementVectors(const Position& p1, Position p2) const {
  std::vector<Displacement> displacements;
  // cycle through all neighbor cell images based on allowed periodicity
  // p1 and p2 have to already be within cell
  // not the case if we are missing one or more periodicities
  assert(!_periodicity[0] || !_periodicity[1] || !_periodicity[2] || isWithinCell(p1));
  assert(!_periodicity[0] || !_periodicity[1] || !_periodicity[2] || isWithinCell(p2));
  int iRange = (_periodicity[0]) ? 1 : 0;
  int jRange = (_periodicity[1]) ? 1 : 0;
  int kRange = (_periodicity[2]) ? 1 : 0;
  for (int i = -1 * iRange; i <= iRange; ++i) {
    for (int j = -1 * jRange; j <= jRange; ++j) {
      for (int k = -1 * kRange; k <= kRange; ++k) {
        Eigen::RowVector3d shift;
        shift << i, j, k;
        transformInPlace(shift);
        p2 += shift;
        // Vector pointing from p1 to p2
        displacements.push_back((p2 - p1));
        p2 -= shift;
      }
    }
  }
  return displacements;
}

bool PeriodicBoundaries::minimumDistanceViaImage(Position p1, Position p2) const {
  translatePositionsIntoCellInPlace(p1);
  translatePositionsIntoCellInPlace(p2);
  auto minDisplacement = bruteForceMinimumImageDisplacementVector(p1, p2);
  Position minDistanceP2 = p1 + minDisplacement;
  return !isWithinCell(minDistanceP2);
}

bool PeriodicBoundaries::isWithinCell(const Position& position) const {
  Position pos = transform(position, false);
  return 1.0 > pos(0) && pos(0) >= 0.0 && 1.0 > pos(1) && pos(1) >= 0.0 && 1.0 > pos(2) && pos(2) >= 0.0;
}

bool PeriodicBoundaries::isWithinCell(const PositionCollection& positions) const {
  for (int i = 0; i < positions.rows(); ++i) {
    Position pos = positions.row(i);
    if (!isWithinCell(pos)) {
      return false;
    }
  }
  return true;
}

void PeriodicBoundaries::canonicalize() {
  auto canonicPbc = PeriodicBoundaries(getLengths(), getAngles(), true, true, getPeriodicityString());
  setCellMatrix(canonicPbc.getCellMatrix());
}

PeriodicBoundaries PeriodicBoundaries::canonicalized() const {
  auto canonicPbc = *this;
  canonicPbc.canonicalize();
  return canonicPbc;
}

Eigen::Matrix3d PeriodicBoundaries::getCanonicalizationRotationMatrix() const {
  auto canonicPbc = PeriodicBoundaries(getLengths(), getAngles(), true, true, getPeriodicityString());
  if (canonicPbc.getCellMatrix().isApprox(getCellMatrix())) {
    // we are already canonic, avoid numeric instability, just return I
    return Eigen::Matrix3d::Identity();
  }
  Eigen::Matrix3d rotationMatrix = this->_invMatrix * canonicPbc.getCellMatrix();
  return rotationMatrix;
}

PeriodicBoundaries& PeriodicBoundaries::operator=(const PeriodicBoundaries& rhs) {
  this->setCellMatrix(rhs._matrix);
  this->setPeriodicity(rhs._periodicity);
  return *this;
}

PeriodicBoundaries& PeriodicBoundaries::operator=(const Eigen::Matrix3d& rhs) {
  this->setCellMatrix(rhs);
  return *this;
}

bool PeriodicBoundaries::isApprox(const PeriodicBoundaries& rhs, double eps) const {
  if (getPeriodicity() != rhs.getPeriodicity()) {
    return false;
  }
  if (getCellMatrix().isApprox(rhs.getCellMatrix(), eps)) {
    return true;
  }
  return this->canonicalized().getCellMatrix().isApprox(rhs.canonicalized().getCellMatrix(), eps);
}

bool operator==(const PeriodicBoundaries& lhs, const PeriodicBoundaries& rhs) {
  return lhs.isApprox(rhs, 1e-12);
}

bool operator!=(const PeriodicBoundaries& lhs, const PeriodicBoundaries& rhs) {
  return !(lhs == rhs);
}

PeriodicBoundaries PeriodicBoundaries::operator*(const Eigen::Vector3d& scalingFactors) const {
  PeriodicBoundaries super = *this;
  return super *= scalingFactors;
}

PeriodicBoundaries PeriodicBoundaries::operator*(const std::vector<double>& scalingFactors) const {
  PeriodicBoundaries super = *this;
  return super *= scalingFactors;
}

PeriodicBoundaries& PeriodicBoundaries::operator*=(const std::vector<double>& scalingFactors) {
  if (scalingFactors.size() != 3) {
    throw std::runtime_error("Scaling factor of PeriodicBoundaries must be either a number or a vector of length 3");
  }
  const Eigen::Vector3d vec = Eigen::Map<const Eigen::Vector3d>(scalingFactors.data(), 3);
  this->operator*=(vec);
  return *this;
}

PeriodicBoundaries& PeriodicBoundaries::operator*=(const Eigen::Vector3d& scalingFactors) {
  Eigen::Matrix3d newMatrix = this->_matrix.array().colwise() * scalingFactors.array();
  this->setCellMatrix(newMatrix);
  return *this;
}

PeriodicBoundaries PeriodicBoundaries::operator*(double scalingFactor) const {
  PeriodicBoundaries super = *this;
  return super *= scalingFactor;
}

PeriodicBoundaries& PeriodicBoundaries::operator*=(double scalingFactor) {
  this->setCellMatrix(_matrix * scalingFactor);
  return *this;
}

PeriodicBoundaries PeriodicBoundaries::operator+(const Eigen::Matrix3d& matrix) const {
  auto super = PeriodicBoundaries(*this);
  super += matrix;
  return super;
}

PeriodicBoundaries& PeriodicBoundaries::operator+=(const Eigen::Matrix3d& matrix) {
  Eigen::Matrix3d newMatrix = this->_matrix + matrix;
  this->setCellMatrix(newMatrix);
  return *this;
}

PeriodicBoundaries PeriodicBoundaries::operator+(const PeriodicBoundaries& other) const {
  return *this + other.getCellMatrix();
}

PeriodicBoundaries& PeriodicBoundaries::operator+=(const PeriodicBoundaries& other) {
  return *this += other.getCellMatrix();
}

} // namespace Utils
} // namespace Scine
