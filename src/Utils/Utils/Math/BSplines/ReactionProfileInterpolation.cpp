/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Math/BSplines/ReactionProfileInterpolation.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Math/QuaternionFit.h"

namespace Scine {
namespace Utils {
namespace BSplines {

TrajectorySpline::TrajectorySpline(const ElementTypeCollection& e, const Eigen::VectorXd& k, const Eigen::MatrixXd& d, double ts)
  : elements(e), knots(k), data(d), tsPosition(ts) {
}

std::tuple<double, AtomCollection> TrajectorySpline::evaluate(const double& pos, const unsigned int& degree) const {
  if (pos < 0.0 || pos > 1.0) {
    throw std::runtime_error("Molecular B-Splines are only valid in the interval [0.0, 1.0].");
  }
  unsigned int posIdx = 0;
  while (pos >= knots[posIdx + 1] && posIdx < (knots.size() - 2)) {
    posIdx += 1;
  }

  // Pad and extract data
  Eigen::MatrixXd d(degree + 1, data.cols());
  for (int i = 0; i < int(degree) + 1; i++) {
    // This is the original formula given for the de Boor's algorithm
    //
    //   int idx = std::max(0, int(i + posIdx - degree));
    //   idx = std::max(0, idx);
    //   idx = std::min(int(data.rows() - 1), idx);
    //
    //   for some reason, using it results in a shift of the fitted
    //   curve, the code below takes care of removing that shift.
    //   It may very well be that there is another bug in the code
    //   below that generates this shift in the first place.
    //   Please alter this with caution.
    int idx = degree;
    if (degree != 0) {
      idx = idx - 3;
      idx = (idx < 0) ? int(double(idx) / 2.0 - 0.5) : int(double(idx) / 2.0 + 0.5);
    }
    else {
      idx = -1;
    }
    idx = i + int(posIdx) - idx - 1;
    idx = std::max(0, idx);
    idx = std::min(int(data.rows() - 1), idx);
    d.row(i) = data.row(idx);
  }
  posIdx += degree;
  Eigen::VectorXd padded(knots.size() + 2 * degree);
  padded.segment(degree, knots.size()) = knots;
  for (unsigned int i = 0; i < degree + 1; i++) {
    padded[i] = 0.0;
    padded[padded.size() - 1 - i] = 1.0;
  }

  for (unsigned int i = 1; i < degree + 1; i++) {
    for (unsigned int j = degree; j > i - 1; j--) {
      const double alpha = (pos - padded[j + posIdx - degree]) / (padded[j + 1 + posIdx - i] - padded[j + posIdx - degree]);
      d.row(j) = (1.0 - alpha) * d.row(j - 1) + alpha * d.row(j);
    }
  }
  double energy = d.row(degree)[0];
  Eigen::VectorXd posVector = d.row(degree).segment(1, this->elements.size() * 3);
  PositionCollection positions = Eigen::Map<PositionCollection>(posVector.data(), this->elements.size(), 3);
  return {energy, AtomCollection(this->elements, positions)};
}

void ReactionProfileInterpolation::clear() {
  _start = nullptr;
  _data.resize(0, 0);
  _nPoints = 0;
  _tsIdx = nullptr;
}

void ReactionProfileInterpolation::appendStructure(const Scine::Utils::AtomCollection& atoms, const double& energy, bool isTS) {
  if (this->_start != nullptr) {
    if (atoms.size() != this->_start->size()) {
      throw std::runtime_error(
          "Tried to add structure to ReactionProfileInterpolation with different number of atoms than stored ones.");
    }
    const auto& e1 = atoms.getElements();
    const auto& e2 = this->_start->getElements();
    for (int i = 0; i < atoms.size(); i++) {
      if (e1[i] != e2[i]) {
        throw std::runtime_error(
            "Tried to add structure to ReactionProfileInterpolation with different elements than stored ones.");
      }
    }
  }
  if (this->_start == nullptr) {
    _start = std::make_unique<Scine::Utils::AtomCollection>(atoms);
    _nPoints = 1;
    this->_data.resize(1, atoms.size() * 3 + 1);
  }
  else {
    _nPoints += 1;
    _data.conservativeResize(_nPoints, this->_data.cols());
  }
  if (isTS) {
    _tsIdx = std::make_unique<unsigned int>(_nPoints - 1);
  }
  _data.row(_nPoints - 1)[0] = energy;
  auto positions = atoms.getPositions();

  if (_nPoints > 1 && this->useQuaternionFit == true) {
    // Quaternion fit to last structure
    Eigen::VectorXd lasPosVector = this->_data.row(_nPoints - 2).segment(1, this->_data.cols() - 1);
    PositionCollection lastPos = Eigen::Map<PositionCollection>(lasPosVector.data(), atoms.size(), 3);
    Utils::QuaternionFit fit(lastPos, positions, this->_start->getElements());
    PositionCollection newPos = fit.getFittedData();
    Eigen::VectorXd newPosVector = Eigen::Map<Eigen::VectorXd>(newPos.data(), atoms.size() * 3);
    this->_data.row(_nPoints - 1).segment(1, this->_data.cols() - 1) = newPosVector.eval();
  }
  else {
    // Add as is
    Eigen::VectorXd posVector = Eigen::Map<Eigen::VectorXd>(positions.data(), atoms.size() * 3);
    this->_data.row(_nPoints - 1).segment(1, this->_data.cols() - 1) = posVector.eval();
  }
}

double ReactionProfileInterpolation::getCurrentTSPosition() const {
  if (_tsIdx == nullptr) {
    throw std::runtime_error("No transition state has been given to interpolation, cannot report its position.");
  }
  // Massweight coordinate part
  Eigen::MatrixXd massWeighted = this->_data.block(0, 1, _nPoints, this->_data.cols() - 1);
  const unsigned int nAtoms = this->_start->size();
  for (unsigned int i = 0; i < nAtoms; i++) {
    const auto& element = this->_start->getElements()[i];
    massWeighted.col(3 * i + 0).array() *= Utils::ElementInfo::mass(element);
    massWeighted.col(3 * i + 1).array() *= Utils::ElementInfo::mass(element);
    massWeighted.col(3 * i + 2).array() *= Utils::ElementInfo::mass(element);
  }

  // Generate distances
  Eigen::VectorXd distances = Eigen::VectorXd::Zero(_nPoints);
  for (unsigned int i = 1; i < _nPoints; i++) {
    const double norm = (massWeighted.row(i - 1) - massWeighted.row(i)).norm();
    distances[i] = distances[i - 1] + norm;
  }
  distances.array() /= distances[_nPoints - 1];
  return distances[*_tsIdx];
}

TrajectorySpline ReactionProfileInterpolation::spline(unsigned int nInterpolationPoints, unsigned int degree) const {
  // Check inpputs
  if (nInterpolationPoints < 3 && _tsIdx != nullptr) {
    throw std::runtime_error("Number of interpolation points needs to be at least 3.");
  }
  if (nInterpolationPoints < 2 && _tsIdx == nullptr) {
    throw std::runtime_error("Number of interpolation points needs to be at least 2.");
  }
  // Check TS
  if (_tsIdx != nullptr) {
    if ((*_tsIdx) == 0 || (*_tsIdx) == (_nPoints - 1)) {
      throw std::runtime_error("TS cannot be the same as start or end of reaction profile.");
    }
  }

  // Massweight coordinate part
  Eigen::MatrixXd massWeighted = this->_data.block(0, 1, _nPoints, this->_data.cols() - 1);
  const unsigned int nAtoms = this->_start->size();
  for (unsigned int i = 0; i < nAtoms; i++) {
    const auto& element = this->_start->getElements()[i];
    massWeighted.col(3 * i + 0).array() *= Utils::ElementInfo::mass(element);
    massWeighted.col(3 * i + 1).array() *= Utils::ElementInfo::mass(element);
    massWeighted.col(3 * i + 2).array() *= Utils::ElementInfo::mass(element);
  }

  // Generate distances
  Eigen::VectorXd distances = Eigen::VectorXd::Zero(_nPoints);
  for (unsigned int i = 1; i < _nPoints; i++) {
    const double norm = (massWeighted.row(i - 1) - massWeighted.row(i)).norm();
    distances[i] = distances[i - 1] + norm;
  }
  distances.array() /= distances[_nPoints - 1]; // normalize to [0,1] -> knot vector

  TrajectorySpline fullSpline(this->_start->getElements(), distances, this->_data);

  // Evenly interpolate start -> end
  Eigen::VectorXd compressedKnots(nInterpolationPoints);
  compressedKnots[0] = 0.0;
  compressedKnots[nInterpolationPoints - 1] = 1.0;
  if (_tsIdx != nullptr) {
    // in case of a TS in the trajectory make sure it is explicitly added
    const double tsPosition = distances[*_tsIdx];
    unsigned int off = 0;
    const double interval = 1.0 / (nInterpolationPoints - 2);
    for (unsigned int i = 0; i < nInterpolationPoints - 2; i++) {
      if (compressedKnots[i + off] + interval >= tsPosition && tsPosition >= compressedKnots[i + off]) {
        compressedKnots[i + 1] = tsPosition;
        off = 1;
        compressedKnots[i + 1 + off] = compressedKnots[i] + interval;
      }
      else {
        compressedKnots[i + 1 + off] = compressedKnots[i + off] + interval;
      }
    }
    compressedKnots[nInterpolationPoints - 1] = 1.0;
  }
  else {
    const double interval = 1.0 / (nInterpolationPoints - 1);
    for (unsigned int i = 0; i < nInterpolationPoints - 2; i++) {
      compressedKnots[i + 1] = compressedKnots[i] + interval;
    }
  }

  // Interpolate compressed data from full spline
  Eigen::MatrixXd compressedData(nInterpolationPoints, this->_data.cols());
  for (unsigned int i = 0; i < nInterpolationPoints; i++) {
    auto interpolated = fullSpline.evaluate(compressedKnots[i], degree);
    auto positions = std::get<1>(interpolated).getPositions();
    Eigen::VectorXd dataVector(this->_data.cols());
    dataVector[0] = std::get<0>(interpolated);
    Eigen::VectorXd posVector = Eigen::Map<Eigen::VectorXd>(positions.data(), this->_data.cols() - 1);
    dataVector.segment(1, this->_data.cols() - 1) = posVector;
    compressedData.row(i) = dataVector;
  }

  return TrajectorySpline(this->_start->getElements(), compressedKnots, compressedData,
                          (_tsIdx != nullptr) ? distances[*_tsIdx] : -1.0);
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
