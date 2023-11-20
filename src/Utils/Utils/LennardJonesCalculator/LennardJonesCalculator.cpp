/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/LennardJonesCalculator/LennardJonesCalculator.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/Utilities/Distances.h"
#include "Utils/LennardJonesCalculator/LennardJonesCalculatorSettings.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

LennardJonesCalculator::LennardJonesCalculator() {
  settings_ = std::make_shared<LennardJonesCalculatorSettings>();
  applySettings();
}

void LennardJonesCalculator::applySettings() {
  if (!settings_->valid()) {
    settings_->throwIncorrectSettings();
  }
  // Set-up periodic boundaries
  auto pbcString = settings_->getString(Utils::SettingsNames::periodicBoundaries);
  pbc_ = (pbcString.empty()) ? nullptr : std::make_shared<PeriodicBoundaries>(pbcString);
  if (pbc_ && std::sqrt(pbc_->getSmallestPerpendicularSquared()) <=
                  2 * settings_->getDouble(Utils::SettingsNames::lennardJonesCutoff)) {
    throw Core::InitializationException("Lennard-Jones box size has to be larger than twice the cut-off distance.");
  }
  // Set-up potential
  cutoff_ = settings_->getDouble(Utils::SettingsNames::lennardJonesCutoff);
  sigma_ = settings_->getDouble(Utils::SettingsNames::lennardJonesSigma);
  // Convert epsilon/k_B in K to epsilon in E_h
  epsilon_ = settings_->getDouble(Utils::SettingsNames::lennardJonesEpsilon) * Constants::boltzmannConstant *
             Constants::hartree_per_joule;
}

bool LennardJonesCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "LENNARDJONES";
}
void LennardJonesCalculator::setStructure(const AtomCollection& structure) {
  applySettings();
  structure_ = structure;
  results_ = Results{};
}
void LennardJonesCalculator::modifyPositions(PositionCollection newPositions) {
  structure_.setPositions(newPositions);
  results_ = Results{};
}
const PositionCollection& LennardJonesCalculator::getPositions() const {
  return structure_.getPositions();
}
void LennardJonesCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}
PropertyList LennardJonesCalculator::getRequiredProperties() const {
  return requiredProperties_;
}
PropertyList LennardJonesCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::StressTensor;
}
const Results& LennardJonesCalculator::calculate(std::string description) {
  applySettings();
  auto positions = structure_.getPositions();
  // Translate positions into cell
  if (pbc_) {
    pbc_->translatePositionsIntoCellInPlace(positions);
  }

  const unsigned int nAtoms = structure_.size();
  double energy = 0.0;
  GradientCollection g(structure_.size(), 3);
  g.setZero();
  const bool calcStress = requiredProperties_.containsSubSet(Utils::Property::StressTensor);
  Eigen::Matrix3d s = Eigen::Matrix3d::Zero();

  for (unsigned int i = 0; i < nAtoms; i++) {
    Position posI = positions.row(i);
    for (unsigned int j = 0; j < i; j++) {
      /*
       * Energy
       */
      Position posJ = positions.row(j);
      const Eigen::Vector3d r = calculateDistanceVector(posJ, posI);
      const double dist = r.norm();
      if (dist < cutoff_) {
        const double lj = sigma_ / dist;
        const double lj6 = lj * lj * lj * lj * lj * lj;
        const double lj12 = lj6 * lj6;
        // compose energy contribution
        energy += 4 * epsilon_ * (lj12 - lj6);
        /*
         * Gradient
         */
        // Partial derivative of E w.r.t. r
        const double rDeriv = 24.0 * epsilon_ * (lj6 / dist - 2 * lj12 / dist);
        // Derivative of E w.r.t. x,y,z
        const double gx = (rDeriv / dist) * r[0];
        const double gy = (rDeriv / dist) * r[1];
        const double gz = (rDeriv / dist) * r[2];
        g(i, 0) += gx;
        g(i, 1) += gy;
        g(i, 2) += gz;
        g(j, 0) -= gx;
        g(j, 1) -= gy;
        g(j, 2) -= gz;

        /*
         * Stress tensor according to Virial
         */
        if (calcStress) {
          // negative because of force
          s(0, 0) -= gx * r[0]; // Fx rx
          s(1, 1) -= gy * r[1]; // Fy ry
          s(2, 2) -= gz * r[2]; // Fz rz

          s(1, 0) -= gx * r[1]; // Fx ry
          s(2, 0) -= gx * r[2]; // Fx rz

          s(0, 1) -= gx * r[1]; // Fy rx
          s(0, 2) -= gx * r[2]; // Fz rx

          s(2, 1) -= gy * r[2]; // Fy rx
          s(1, 2) -= gz * r[1]; // Fz ry
          s *= 0.5;
        }
      }
    }
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::Description>(std::move(description));
  results_.set<Property::Energy>(energy);
  results_.set<Property::Gradients>(g);
  if (calcStress) {
    results_.set<Property::StressTensor>(s);
  }

  return results_;
}
std::string LennardJonesCalculator::name() const {
  return std::string("LennardJonesCalculator");
}
const Settings& LennardJonesCalculator::settings() const {
  return *settings_;
}
Settings& LennardJonesCalculator::settings() {
  return *settings_;
}
std::shared_ptr<Core::State> LennardJonesCalculator::getState() const {
  return nullptr;
}
void LennardJonesCalculator::loadState(std::shared_ptr<Core::State> /* state */) {
}
Utils::Results& LennardJonesCalculator::results() {
  return results_;
}
const Utils::Results& LennardJonesCalculator::results() const {
  return results_;
}
std::unique_ptr<Utils::AtomCollection> LennardJonesCalculator::getStructure() const {
  return std::make_unique<AtomCollection>(structure_);
}

Displacement LennardJonesCalculator::calculateDistanceVector(const Position& p1, Position& p2) const {
  if (pbc_) {
    return pbc_->bruteForceMinimumImageDisplacementVector(p1, p2);
  }

  return p2 - p1;
}

} // namespace Utils
} // namespace Scine
