/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SinglePointMethod.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/UniversalSettings/Exceptions.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {

using namespace Utils::AutomaticDifferentiation;

SinglePointMethod::SinglePointMethod(Utils::derivOrder maximalOrder)
  : maximalCalculableDerivativeOrder_(maximalOrder), energy_(0) {
}

void SinglePointMethod::setAtomCollection(const Utils::AtomCollection& structure) {
  initializeStructure(structure.getElements(), structure.getPositions());
}

void SinglePointMethod::initializeStructure(const Utils::ElementTypeCollection& elements) {
  elementTypes_ = elements;
}

void SinglePointMethod::initializeStructure(const Utils::ElementTypeCollection& elements,
                                            const Utils::PositionCollection& positions) {
  assert(positions.rows() == elements.size() &&
         "Error: Position vector does not have the same size as the element vector!");
  initializeStructure(elements);
  setPositions(positions);
}

void SinglePointMethod::setPositions(Utils::PositionCollection positions) {
  assert(positions.rows() == elementTypes_.size() && "Error: Position vector does not have the same size as the "
                                                     "element vector! Call initializeStructure(...) first.");
  positions_ = std::move(positions);
}

const Utils::PositionCollection& SinglePointMethod::getPositions() const {
  return positions_;
}

const Utils::GradientCollection& SinglePointMethod::getGradients() const {
  return gradients_;
}

const Utils::AtomicSecondDerivativeCollection& SinglePointMethod::getAtomicSecondDerivatives() const {
  return secondDerivatives_;
}

const Utils::FullSecondDerivativeCollection& SinglePointMethod::getFullSecondDerivatives() const {
  return fullSecondDerivatives_;
}

const Utils::ElementTypeCollection& SinglePointMethod::getElementTypes() const {
  return elementTypes_;
}

bool SinglePointMethod::canCalculateSecondDerivatives() const {
  return maximalCalculableDerivativeOrder_ == Utils::derivOrder::two;
}

const Utils::BondOrderCollection& SinglePointMethod::getBondOrderCollection() const {
  return bondOrders_;
}

void SinglePointMethod::setBondOrderCollection(const Utils::BondOrderCollection& B) {
  bondOrders_ = B;
}

double SinglePointMethod::getEnergy() const {
  return energy_;
}

void SinglePointMethod::setEnergy(double energy) {
  energy_ = energy;
}

const std::vector<double>& SinglePointMethod::getAtomicCharges() const {
  return atomicCharges_;
};

double SinglePointMethod::getAtomicCharge(int index) const {
  return atomicCharges_[index];
};

void SinglePointMethod::setAtomicCharges(const std::vector<double>& charges) {
  atomicCharges_ = charges;
};

void SinglePointMethod::startLogger(const std::string& loggerSeverity) const {
  if (loggerSeverity == SettingsNames::LogLevels::none)
    Log::stopConsoleLogging();
  else if (loggerSeverity == SettingsNames::LogLevels::trace)
    Log::startConsoleLogging(Log::severity_level::trace);
  else if (loggerSeverity == SettingsNames::LogLevels::debug)
    Log::startConsoleLogging(Log::severity_level::debug);
  else if (loggerSeverity == SettingsNames::LogLevels::info)
    Log::startConsoleLogging(Log::severity_level::info);
  else if (loggerSeverity == SettingsNames::LogLevels::warning)
    Log::startConsoleLogging(Log::severity_level::warning);
  else if (loggerSeverity == SettingsNames::LogLevels::error)
    Log::startConsoleLogging(Log::severity_level::error);
  else if (loggerSeverity == SettingsNames::LogLevels::fatal)
    Log::startConsoleLogging(Log::severity_level::fatal);
  else
    throw UniversalSettings::OptionDoesNotExistException(loggerSeverity, SettingsNames::loggerVerbosity);
}

void SinglePointMethod::resizeRealTimeMethodMembers() {
  gradients_ = Utils::GradientCollection{getNumberAtoms(), 3};
  secondDerivatives_ = Utils::AtomicSecondDerivativeCollection(getNumberAtoms());
  fullSecondDerivatives_ = Utils::FullSecondDerivativeCollection(getNumberAtoms());
  positions_.resize(getNumberAtoms(), 3);
  bondOrders_.resize(getNumberAtoms());
  atomicCharges_.resize(getNumberAtoms());
}
} // namespace Utils
} // namespace Scine
