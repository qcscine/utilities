/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "QmmmEmbeddingTestCalculator.h"
#include "Core/Interfaces/Calculator.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/CalculatorBasics/TestCalculator.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/Settings.h"

namespace Scine {
namespace Utils {

// Define mock settings
class QmmmEmbeddingTestCalculatorSettings : public Scine::Utils::Settings {
 public:
  QmmmEmbeddingTestCalculatorSettings() : Settings("QmmmEmbeddingTestCalculatorSettings") {
    Utils::UniversalSettings::BoolDescriptor ignoreQm(
        "Whether to ignore all contributions from the QM calculation, and therefore, not performing it.");
    ignoreQm.setDefaultValue(false);
    this->_fields.push_back("ignore_qm", std::move(ignoreQm));
    Utils::UniversalSettings::IntListDescriptor qmAtoms("A list of the indices of the atoms in the QM region.");
    this->_fields.push_back("qm_atoms", std::move(qmAtoms));
    Utils::UniversalSettings::BoolDescriptor electrostaticEmbedding(
        "Sets whether electrostatic embedding is used in QM/MM. The alternative is applying mechanical embedding "
        "only.");
    electrostaticEmbedding.setDefaultValue(false);
    this->_fields.push_back("electrostatic_embedding", std::move(electrostaticEmbedding));
    Utils::UniversalSettings::BoolDescriptor optimizeLinks(
        "Whether to optimize the position of the link nuclei before reporting an energy.");
    optimizeLinks.setDefaultValue(false);
    this->_fields.push_back("optimize_links", std::move(optimizeLinks));

    this->resetToDefaults();
  };
  ~QmmmEmbeddingTestCalculatorSettings() override = default;
};

QmmmEmbeddingTestCalculator::QmmmEmbeddingTestCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  settings_ = std::make_unique<QmmmEmbeddingTestCalculatorSettings>();
}
QmmmEmbeddingTestCalculator::QmmmEmbeddingTestCalculator(const QmmmEmbeddingTestCalculator& rhs) : CloneInterface(rhs) {
  structure_ = rhs.structure_;
}
void QmmmEmbeddingTestCalculator::setUnderlyingCalculators(std::vector<std::shared_ptr<Calculator>> underlyingCalculators) {
  setUnderlyingCalculatorsImpl(underlyingCalculators);
}

void QmmmEmbeddingTestCalculator::setUnderlyingCalculatorsImpl(std::vector<std::shared_ptr<Calculator>> underlyingCalculators) {
  this->qmCalculator_ = underlyingCalculators.at(0)->clone();
  this->mmCalculator_ = underlyingCalculators.at(1)->clone();
}

std::vector<std::shared_ptr<Core::Calculator>> QmmmEmbeddingTestCalculator::getUnderlyingCalculators() const {
  return std::vector<std::shared_ptr<Core::Calculator>>{qmCalculator_, mmCalculator_};
}

void QmmmEmbeddingTestCalculator::addUnderlyingSettings() {
}

void QmmmEmbeddingTestCalculator::setStructure(const Utils::AtomCollection& structure) {
  setStructureImpl(structure);
}

void QmmmEmbeddingTestCalculator::setStructureImpl(const Utils::AtomCollection& structure) {
  structure_ = structure;
  mmCalculator_->setStructure(structure_);
  Utils::AtomCollection qmRegion;
  for (const auto& qmAtomIndex : settings_->getIntList("qm_atoms")) {
    qmRegion.push_back(structure.at(qmAtomIndex));
  }
  // For test purposes, the whole system is QM
  qmRegion_ = structure;

  if (results_.has<Utils::Property::SuccessfulCalculation>() &&
      qmCalculator_->getStructure()->getElements() == qmRegion_.getElements()) {
    qmCalculator_->modifyPositions(qmRegion_.getPositions());
  }
  else {
    qmCalculator_->setStructure(qmRegion_);
  }
}
std::unique_ptr<Utils::AtomCollection> QmmmEmbeddingTestCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(structure_);
}
void QmmmEmbeddingTestCalculator::modifyPositions(PositionCollection newPositions) {
  structure_.setPositions(newPositions);
  mmCalculator_->modifyPositions(std::move(newPositions));
  qmCalculator_->modifyPositions(qmRegion_.getPositions());
}
const PositionCollection& QmmmEmbeddingTestCalculator::getPositions() const {
  return structure_.getPositions();
}
void QmmmEmbeddingTestCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  // qmmm properties
  requiredProperties_ = requiredProperties;
  // qm properties
  auto qmProperties = requiredProperties_;
  // mm properties
  auto mmProperties = requiredProperties_;
  mmCalculator_->setRequiredProperties(mmProperties);
  qmCalculator_->setRequiredProperties(qmProperties);
}
PropertyList QmmmEmbeddingTestCalculator::getRequiredProperties() const {
  return requiredProperties_;
}
PropertyList QmmmEmbeddingTestCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian;
}
const Results& QmmmEmbeddingTestCalculator::calculate(std::string dummy) {
  try {
    return calculateImpl(dummy);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

const Results& QmmmEmbeddingTestCalculator::calculateImpl(std::string dummy) {
  auto qmResults = qmCalculator_->calculate();
  double qmEnergy = qmResults.get<Utils::Property::Energy>();
  auto mmResults = mmCalculator_->calculate();
  double mmEnergy = mmResults.get<Utils::Property::Energy>();

  // Total energy:
  double totalEnergy = qmEnergy + mmEnergy;

  // Assemble results
  results_.set<Utils::Property::Description>(std::move(dummy));
  results_.set<Utils::Property::Energy>(totalEnergy);
  results_.set<Utils::Property::Gradients>(mmResults.get<Utils::Property::Gradients>());
  if (requiredProperties_.containsSubSet(Utils::Property::Hessian)) {
    results_.set<Utils::Property::Hessian>(mmResults.get<Utils::Property::Hessian>());
  }
  results_.set<Property::SuccessfulCalculation>(true);
  return results_;
}
std::string QmmmEmbeddingTestCalculator::name() const {
  return "QmmmEmbeddingTestCalculator";
}
std::shared_ptr<Core::State> QmmmEmbeddingTestCalculator::getState() const {
  return nullptr;
}
void QmmmEmbeddingTestCalculator::loadState(std::shared_ptr<Core::State> /* state */) {
}
const Settings& QmmmEmbeddingTestCalculator::settings() const {
  return *settings_;
}
Settings& QmmmEmbeddingTestCalculator::settings() {
  return *settings_;
}
Utils::Results& QmmmEmbeddingTestCalculator::results() {
  return results_;
}
const Utils::Results& QmmmEmbeddingTestCalculator::results() const {
  return results_;
}
bool QmmmEmbeddingTestCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "QMMM";
}
bool QmmmEmbeddingTestCalculator::allowsPythonGILRelease() const {
  return true;
}

} // namespace Utils
} // namespace Scine
