/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaCalculator.h"
#include "OrcaCalculatorSettings.h"
#include "OrcaHessianOutputParser.h"
#include "OrcaInputFileCreator.h"
#include "OrcaMainOutputParser.h"
#include "OrcaStatesHandler.h"
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <iostream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string OrcaCalculator::name() const {
  return "ORCA";
}

OrcaCalculator::OrcaCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<OrcaCalculatorSettings>();
  this->statesHandler_ = std::make_unique<OrcaStatesHandler>(*this);
  applySettings();
}

OrcaCalculator::OrcaCalculator(const OrcaCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  ;
  applySettings();
  this->results() = rhs.results();
  this->setStructure(rhs.atoms_);
}

void OrcaCalculator::applySettings() {
  if (settings_->check()) {
    orcaExecutable_ = settings_->getString(SettingsNames::orcaBinaryPath);
    fileNameBase_ = settings_->getString(SettingsNames::orcaFilenameBase);
    baseWorkingDirectory_ = settings_->getString(SettingsNames::baseWorkingDirectory);
  }
  else {
    throw Core::InitializationException("Settings are invalid!");
  }
}

void OrcaCalculator::setStructure(const Utils::AtomCollection& structure) {
  atoms_ = structure;
  calculationDirectory_ = createNameForCalculationDirectory();
}

std::unique_ptr<Utils::AtomCollection> OrcaCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(atoms_);
}

void OrcaCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  atoms_.setPositions(std::move(newPositions));
}

const Utils::PositionCollection& OrcaCalculator::getPositions() const {
  return atoms_.getPositions();
}

void OrcaCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList OrcaCalculator::possibleProperties() const {
  return Property::Energy | Property::Gradients | Property::Hessian;
}

const Results& OrcaCalculator::calculate(std::string description) {
  applySettings();
  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

const Results& OrcaCalculator::calculateImpl(std::string description) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(calculationDirectory_);
  externalProgram.createWorkingDirectory();
  std::string inputFile = externalProgram.generateFullFilename(fileNameBase_ + ".inp");
  std::string outputFile = externalProgram.generateFullFilename(fileNameBase_ + ".out");

  InputFileCreator inputFileCreator;
  inputFileCreator.createInputFile(inputFile, atoms_, *settings_, requiredProperties_);

  externalProgram.executeCommand(orcaExecutable_ + " " + inputFile, outputFile);

  OrcaMainOutputParser parser(outputFile);

  results_.setDescription(description);
  if (requiredProperties_.containsSubSet(Property::Energy))
    results_.setEnergy(parser.getEnergy());
  if (requiredProperties_.containsSubSet(Property::Gradients))
    results_.setGradients(parser.getGradients());
  if (requiredProperties_.containsSubSet(Property::Hessian)) {
    std::string hessianFile = externalProgram.generateFullFilename(fileNameBase_ + ".hess");
    OrcaHessianOutputParser hessianParser(hessianFile);
    results_.setHessian(hessianParser.getHessian());
  }

  return results_;
}

std::string OrcaCalculator::createNameForCalculationDirectory() {
  boost::uuids::uuid uuid = boost::uuids::random_generator()();
  std::string uuidString = boost::uuids::to_string(uuid);

  auto directoryName = NativeFilenames::combinePathSegments(baseWorkingDirectory_, uuidString);
  return NativeFilenames::addTrailingSeparator(directoryName);
}

const Utils::Settings& OrcaCalculator::settings() const {
  return *settings_;
}

Utils::Settings& OrcaCalculator::settings() {
  return *settings_;
}

const Utils::StatesHandler& OrcaCalculator::statesHandler() const {
  return *statesHandler_;
}

Utils::StatesHandler& OrcaCalculator::statesHandler() {
  return *statesHandler_;
}

const Utils::Results& OrcaCalculator::results() const {
  return results_;
}

Utils::Results& OrcaCalculator::results() {
  return results_;
}

std::string OrcaCalculator::getFileNameBase() {
  return fileNameBase_;
}

std::string OrcaCalculator::getCalculationDirectory() {
  return calculationDirectory_;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
