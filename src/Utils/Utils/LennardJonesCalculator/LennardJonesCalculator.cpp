/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
  const double boxlength = settings_->getDouble(Utils::SettingsNames::lennardJonesBoxsize);
  if (settings_->getBool(Utils::SettingsNames::lennardJonesUsePBCs) &&
      boxlength <= 2 * settings_->getDouble(Utils::SettingsNames::lennardJonesCutoff)) {
    throw Core::InitializationException("Lennard-Jones box size has to be larger than twice the cut-off distance.");
  }
  if (!settings_->valid()) {
    settings_->throwIncorrectSettings();
  }
  cutoff_ = settings_->getDouble(Utils::SettingsNames::lennardJonesCutoff);
  sigma_ = settings_->getDouble(Utils::SettingsNames::lennardJonesSigma);
  // Convert epsilon/k_B in K to epsilon in E_h
  epsilon_ = settings_->getDouble(Utils::SettingsNames::lennardJonesEpsilon) * Constants::boltzmannConstant *
             Constants::hartree_per_joule;
  // Set-up periodic boundaries
  if (settings_->getBool(Utils::SettingsNames::lennardJonesUsePBCs)) {
    Eigen::Vector3d lengths, angles;
    lengths << boxlength, boxlength, boxlength;
    angles << 90.0, 90.0, 90.0;
    pbc_ = std::make_shared<PeriodicBoundaries>(lengths, angles);
  }
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
  return Utils::Property::Energy | Utils::Property::Gradients;
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
        g(i, 0) += (rDeriv / dist) * r[0];
        g(i, 1) += (rDeriv / dist) * r[1];
        g(i, 2) += (rDeriv / dist) * r[2];
        g(j, 0) -= (rDeriv / dist) * r[0];
        g(j, 1) -= (rDeriv / dist) * r[1];
        g(j, 2) -= (rDeriv / dist) * r[2];
      }
    }
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::Description>(std::move(description));
  results_.set<Property::Energy>(energy);
  results_.set<Property::Gradients>(g);

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
