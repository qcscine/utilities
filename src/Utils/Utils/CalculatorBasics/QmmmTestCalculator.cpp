/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Settings.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/QmmmTestCalculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/TestCalculator.h>

namespace Scine {
namespace Utils {

// Define mock settings
class QmmmTestCalculatorSettings : public Scine::Utils::Settings {
 public:
  QmmmTestCalculatorSettings() : Settings("QmmmTestCalculatorSettings") {
    Utils::UniversalSettings::BoolDescriptor ignoreQm(
        "Whether to ignore all contributions from the QM calculation, and therefore, not performing it.");
    ignoreQm.setDefaultValue(false);
    this->_fields.push_back("ignore_qm", std::move(ignoreQm));
    Utils::UniversalSettings::IntListDescriptor qmAtoms("A list of the indices of the atoms in the QM region.");
    this->_fields.push_back("qm_atoms", std::move(qmAtoms));

    this->resetToDefaults();
  };
  ~QmmmTestCalculatorSettings() override = default;
};

QmmmTestCalculator::QmmmTestCalculator() {
  this->_settings = std::make_unique<QmmmTestCalculatorSettings>();
  this->_underlyingCalculator = std::make_shared<TestCalculator>();
  this->_underlyingCalculator->setPrecision(7);
}
void QmmmTestCalculator::setStructure(const AtomCollection& structure) {
  _underlyingCalculator->setStructure(structure);
  _structure = structure;
}
std::unique_ptr<Utils::AtomCollection> QmmmTestCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(_structure);
}
void QmmmTestCalculator::modifyPositions(PositionCollection newPositions) {
  _underlyingCalculator->modifyPositions(newPositions);
  _structure.setPositions(newPositions);
}
const PositionCollection& QmmmTestCalculator::getPositions() const {
  return _underlyingCalculator->getPositions();
}
void QmmmTestCalculator::setRequiredProperties(const PropertyList& /* requiredProperties */) {
}
PropertyList QmmmTestCalculator::getRequiredProperties() const {
  return PropertyList{};
}
PropertyList QmmmTestCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::PartialHessian;
}
const Results& QmmmTestCalculator::calculate(std::string dummy) {
  _results = _underlyingCalculator->calculate(dummy);
  return _results;
}
std::string QmmmTestCalculator::name() const {
  return "QmmmTestCalculator";
}
std::shared_ptr<Core::State> QmmmTestCalculator::getState() const {
  return nullptr;
}
void QmmmTestCalculator::loadState(std::shared_ptr<Core::State> /* state */) {
}
const Settings& QmmmTestCalculator::settings() const {
  return *_settings;
}
Settings& QmmmTestCalculator::settings() {
  return *_settings;
}
Utils::Results& QmmmTestCalculator::results() {
  return _results;
}
const Utils::Results& QmmmTestCalculator::results() const {
  return _results;
}
bool QmmmTestCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "QMMM";
}
bool QmmmTestCalculator::allowsPythonGILRelease() const {
  return true;
}
std::shared_ptr<Core::Calculator> QmmmTestCalculator::cloneImpl() const {
  std::shared_ptr<Core::Calculator> calculator;
  calculator->setStructure(_structure);

  return calculator;
}

} // namespace Utils
} // namespace Scine
