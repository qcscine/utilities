/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/ExternalQC/MRCC/MrccCalculator.h"
#include "Utils/ExternalQC/Exceptions.h"
#include "Utils/ExternalQC/ExternalProgram.h"
#include "Utils/ExternalQC/MRCC/MrccHelper.h"
#include "Utils/ExternalQC/MRCC/MrccIO.h"
#include "Utils/ExternalQC/MRCC/MrccSettings.h"
#include "Utils/ExternalQC/MRCC/MrccState.h"
#include "Utils/IO/NativeFilenames.h"
#include "Utils/Solvation/ImplicitSolvation.h"
#include <Utils/Strings.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string MrccCalculator::name() const {
  return this->name_;
}

MrccCalculator::MrccCalculator()
  : binaryDirectory_(std::getenv(binaryEnvVariable)), settings_(std::make_unique<MrccSettings>()) {
}

MrccCalculator::MrccCalculator(const MrccCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  this->setLog(rhs.getLog());
  applySettings();
  this->setStructure(*rhs.getStructure());
  this->results() = rhs.results();
  this->binaryDirectory_ = rhs.getBinaryDirectory();
}

bool MrccCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  if (std::getenv(binaryEnvVariable)) {
    return caseInsensitiveEqual(methodFamily, this->getMethodFamily());
  }
  return false;
}

void MrccCalculator::setStructure(const AtomCollection& structure) {
  applySettings();
  atoms_ = structure;
  calculationDirectory_ = NativeFilenames::createRandomDirectoryName(baseWorkingDirectory_);
  results_ = Results{};
}

std::unique_ptr<Utils::AtomCollection> MrccCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(atoms_);
}

void MrccCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

PropertyList MrccCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

PropertyList MrccCalculator::possibleProperties() const {
  return Property::Energy | Property::SuccessfulCalculation;
}

const Results& MrccCalculator::calculate(std::string description) {
  applySettings();
  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

const Results& MrccCalculator::calculateImpl(std::string description) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(calculationDirectory_);
  externalProgram.createWorkingDirectory();

  MrccHelper helper(this->binaryDirectory_, this->calculationDirectory_);
  MrccIO io(helper.getFiles(), *this->settings_, this->getMethodFamily());
  io.writeInput(this->atoms_);
  helper.run();
  const std::string output = io.readOutput();

  results_.set<Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Property::Energy)) {
    results_.set<Property::Energy>(io.getEnergy(output));
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::ProgramName>(this->program);
  return results_;
}

const Utils::Settings& MrccCalculator::settings() const {
  return *settings_;
}

Utils::Settings& MrccCalculator::settings() {
  return *settings_;
}

const Utils::Results& MrccCalculator::results() const {
  return results_;
}

Utils::Results& MrccCalculator::results() {
  return results_;
}

void MrccCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  atoms_.setPositions(std::move(newPositions));
  results_ = Results{};
}

const Utils::PositionCollection& MrccCalculator::getPositions() const {
  return atoms_.getPositions();
}

std::shared_ptr<Core::State> MrccCalculator::getState() const {
  auto ret = std::make_shared<MrccState>(this->getCalculationDirectory());
  return ret;
}

std::string MrccCalculator::getCalculationDirectory() const {
  return calculationDirectory_;
}

std::string MrccCalculator::getBinaryDirectory() const {
  return this->binaryDirectory_;
}

void MrccCalculator::applySettings() {
  if (settings_->valid()) {
    Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels_, *settings_);
  }
  else {
    settings_->throwIncorrectSettings();
  }
  this->baseWorkingDirectory_ = settings_->getString(SettingsNames::baseWorkingDirectory);
}

void MrccCalculator::loadState(std::shared_ptr<Core::State> state) {
  auto mrccState = std::dynamic_pointer_cast<MrccState>(state);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine