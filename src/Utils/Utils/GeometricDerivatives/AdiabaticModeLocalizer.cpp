/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometricDerivatives/AdiabaticModeLocalizer.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/GeometricDerivatives/HessianUtilities.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/MolecularTrajectory.h"
#include <Eigen/Eigen>
#include <algorithm>

namespace Scine {
namespace Utils {

// Localized mode storage
void AdiabaticModesContainer::addMode(std::pair<int, int> internalCoordinate, NormalMode mode, double forceConstant) {
  // Ensure that this mode is not already there
  int storageIndex = this->getStorageIndex(internalCoordinate);
  if (storageIndex != storageOrder_.size()) {
    throw AdiabaticLocalizationExceptions::modeReadded();
  };
  // Add mode and force constant
  modesContainer_.add(mode);
  forceConstants_.insert(std::pair<std::pair<int, int>, double>(internalCoordinate, forceConstant));
  forceConstants_.insert(std::pair<std::pair<int, int>, double>(
      std::make_pair(internalCoordinate.second, internalCoordinate.first), forceConstant));
  // Store storage location
  storageOrder_.insert(
      std::pair<std::pair<int, int>, int>(std::minmax(internalCoordinate.first, internalCoordinate.second), storageIndex));
  return;
}

int AdiabaticModesContainer::getStorageIndex(std::pair<int, int> internalCoordinate) const {
  // Ensure proper sorting
  std::pair<int, int> sortedCoord = std::minmax(internalCoordinate.first, internalCoordinate.second);
  // Find value
  auto mapIterator = storageOrder_.find(sortedCoord);
  if (mapIterator != storageOrder_.end()) {
    return storageOrder_.find(sortedCoord)->second;
  }
  else
    return storageOrder_.size();
}

const DisplacementCollection& AdiabaticModesContainer::getMode(std::pair<int, int> internalCoordinate) const {
  int storageIndex = this->getStorageIndex(internalCoordinate);
  if (storageIndex == storageOrder_.size()) {
    throw AdiabaticLocalizationExceptions::InvalidModeRequested();
  };
  return modesContainer_.getMode(storageIndex);
}

MolecularTrajectory AdiabaticModesContainer::getModeAsMolecularTrajectory(std::pair<int, int> internalCoordinate,
                                                                          const Utils::AtomCollection& atoms,
                                                                          double scalingFactor) const {
  int storageIndex = this->getStorageIndex(internalCoordinate);
  if (storageIndex == storageOrder_.size()) {
    throw AdiabaticLocalizationExceptions::InvalidModeRequested();
  };
  return modesContainer_.getModeAsMolecularTrajectory(storageIndex, atoms, scalingFactor);
}

std::map<std::pair<int, int>, double> AdiabaticModesContainer::getWaveNumbers() {
  if (wavenumbers_.size() != 2 * modesContainer_.size()) {
    std::vector<double> waveVector = modesContainer_.getWaveNumbers();
    for (auto const& it : storageOrder_) {
      wavenumbers_.insert(std::pair<std::pair<int, int>, double>(std::make_pair(it.first.first, it.first.second),
                                                                 waveVector.at(it.second)));
      wavenumbers_.insert(std::pair<std::pair<int, int>, double>(std::make_pair(it.first.second, it.first.first),
                                                                 waveVector.at(it.second)));
    }
  }
  return wavenumbers_;
}

std::map<std::pair<int, int>, double> AdiabaticModesContainer::getForceConstants() const {
  return forceConstants_;
}

// Mode Localization
AdiabaticModeLocalizer::AdiabaticModeLocalizer(const Eigen::MatrixXd& hessian, const AtomCollection& atoms,
                                               const BondOrderCollection& bondOrders, const double bondingThreshold)
  : hessian_(std::move(hessian)), atoms_(std::move(atoms)) {
  if (hessian_.rows() != 3 * atoms_.size() || hessian_.cols() != 3 * atoms_.size())
    throw AdiabaticLocalizationExceptions::InvalidHessianDimensions();
  // Get pruned bond list
  Eigen::SparseMatrix<double> bondMatrix = bondOrders.getMatrix().pruned(bondingThreshold).triangularView<Eigen::Upper>();
  for (int j = 0; j < bondMatrix.outerSize(); ++j) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(bondMatrix, j); it; ++it) {
      internalCoords_.push_back(std::make_pair(it.row(), it.col()));
    }
  }
}

AdiabaticModeLocalizer::AdiabaticModeLocalizer(const Eigen::MatrixXd& hessian, const AtomCollection& atoms,
                                               const std::vector<std::pair<int, int>>& internalCoordinates)
  : hessian_(std::move(hessian)), atoms_(std::move(atoms)), internalCoords_(std::move(internalCoordinates)) {
  if (hessian_.rows() != 3 * atoms_.size() || hessian_.cols() != 3 * atoms_.size())
    throw AdiabaticLocalizationExceptions::InvalidHessianDimensions();
}

void AdiabaticModeLocalizer::calculateStretchingB() {
  B_.resize(internalCoords_.size(), 3 * atoms_.size());
  B_.setZero();
  Eigen::MatrixXd positions = atoms_.getPositions();
  for (unsigned int i = 0; i != internalCoords_.size(); ++i) {
    // Bij = del q_j / del x_i
    std::pair<int, int> const coord = internalCoords_[i];
    Eigen::RowVector3d p1 = positions.row(coord.first);
    Eigen::RowVector3d p2 = positions.row(coord.second);
    double distance = (p1 - p2).norm();
    B_.row(i).segment(3 * coord.first, 3) << p1 - p2;
    B_.row(i).segment(3 * coord.second, 3) << p2 - p1;
    B_.row(i) /= distance;
  }
}

AdiabaticModesContainer AdiabaticModeLocalizer::localizeModes() {
  AdiabaticModesContainer modesContainer;
  // Following the notation of https://doi.org/10.1002/wcms.1480
  // Get Wilson B matrix
  this->calculateStretchingB();
  HessianUtilities hessianDiagonalizer(hessian_, atoms_.getElements(), atoms_.getPositions(), false);
  Eigen::MatrixXd L = hessianDiagonalizer.getBackTransformedInternalEigenvectors(); // rototranslation free eigenvectors
  Eigen::MatrixXd K = L.transpose() * hessian_ * L;                                 // force constants/eigenvalues
  Eigen::MatrixXd KInverse = K.inverse();
  Eigen::MatrixXd D = B_ * L; // internal modes
  // Get mass matrix and adiabatic masses
  ElementTypeCollection elements = atoms_.getElements();
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> MInverse(3 * elements.size());
  MInverse.setIdentity();
  for (unsigned int i = 0; i != elements.size(); ++i) {
    MInverse.diagonal().segment(3 * i, 3) /= ElementInfo::mass(elements[i]);
  }
  Eigen::MatrixXd MAdiabaticInverse = B_ * MInverse * B_.transpose(); // Inverse adiabatic masses m_a, G matrix diagonal

  for (unsigned int i = 0; i != internalCoords_.size(); ++i) {
    Eigen::RowVectorXd dn = D.row(i);
    // Adiabatic force constant in a.u.
    double kAdiabatic = 1 / (dn * KInverse * dn.transpose());
    // Adiabatic wave number in a.u., stored as negative if imaginary
    double signConverter = kAdiabatic < 0 ? -1 : 1;
    double wavenumber = signConverter * sqrt(signConverter * kAdiabatic * MAdiabaticInverse.diagonal()[i]);
    // Conversion to cm^-1
    wavenumber *= Constants::invCentimeter_per_hartree * sqrt(Constants::u_per_electronRestMass);
    // Localized mode in cartesians
    Eigen::VectorXd modeAdiabatic = L * KInverse * dn.transpose() / (dn * KInverse * dn.transpose());

    // Transform to displacement collection
    DisplacementCollection dcAdiabatic(atoms_.size(), 3);
    for (int j = 0; j < atoms_.size(); ++j)
      dcAdiabatic.row(j) = modeAdiabatic.segment(3 * j, 3);
    NormalMode formattedMode(wavenumber, dcAdiabatic);
    modesContainer.addMode(internalCoords_.at(i), formattedMode, kAdiabatic);
  }
  return modesContainer;
}
} // namespace Utils
} // namespace Scine
