/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GaussianCalculator.h"
#include "GaussianCalculatorSettings.h"
#include "GaussianInputFileCreator.h"
#include "GaussianOrbitalParser.h"
#include "GaussianOrbitalWriter.h"
#include "GaussianOutputParser.h"
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Solvation/ImplicitSolvation.h>
#include <boost/process.hpp>
#include <regex>

namespace bp = boost::process;

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string GaussianCalculator::name() const {
  return "GAUSSIAN";
}

bool GaussianCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  if (std::getenv("GAUSSIAN_BINARY_PATH")) {
    return methodFamily == "DFT" || methodFamily == "HF";
  }

  return false;
}

GaussianCalculator::GaussianCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<GaussianCalculatorSettings>();

  // Get Gaussian's binary path from an environment variable
  if (const char* pathPtr = std::getenv("GAUSSIAN_BINARY_PATH")) {
    gaussianExecutable_ = std::string(pathPtr);
    gaussianDirectory_ = NativeFilenames::getParentDirectory(gaussianExecutable_);
  }

  applySettings();
}

GaussianCalculator::GaussianCalculator(const GaussianCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  this->setLog(rhs.getLog());
  applySettings();
  atoms_ = rhs.atoms_;
  calculationDirectory_ = NativeFilenames::createRandomDirectoryName(baseWorkingDirectory_);
  results_ = rhs.results();
  this->gaussianExecutable_ = rhs.gaussianExecutable_;
  this->gaussianDirectory_ = rhs.gaussianDirectory_;
  this->binaryHasBeenChecked_ = rhs.binaryHasBeenChecked_;
}

void GaussianCalculator::applySettings() {
  if (settings_->valid()) {
    if (settings_->getDouble(Utils::SettingsNames::electronicTemperature) > 0.0) {
      throw Core::InitializationException(
          "Gaussian calculations with an electronic temperature above 0.0 K are not supported.");
    }
    if (!(settings_->getBool(SettingsNames::enforceScfCriterion)) &&
        (requiredProperties_.containsSubSet(Property::Gradients) || requiredProperties_.containsSubSet(Property::Hessian)) &&
        settings_->getDouble(Utils::SettingsNames::selfConsistenceCriterion) > 1e-8) {
      settings_->modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
      this->getLog().warning << "Warning: Energy accuracy was increased to 1e-8 to ensure valid gradients/hessian."
                             << Core::Log::nl;
    }
    fileNameBase_ = settings_->getString(SettingsNames::gaussianFilenameBase);
    baseWorkingDirectory_ = settings_->getString(SettingsNames::baseWorkingDirectory);
    /* throws error for wrong input and updates 'any' entries */
    /* information if solvation is performed is deduced from settings from InputFileCreator and therefore discarded here
     */
    Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels_, *settings_);
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

void GaussianCalculator::setStructure(const Utils::AtomCollection& structure) {
  applySettings();
  atoms_ = structure;
  calculationDirectory_ = NativeFilenames::createRandomDirectoryName(baseWorkingDirectory_);
  results_ = Results{};
}

std::unique_ptr<Utils::AtomCollection> GaussianCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(atoms_);
}

void GaussianCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  atoms_.setPositions(std::move(newPositions));
  results_ = Results{};
}

const Utils::PositionCollection& GaussianCalculator::getPositions() const {
  return atoms_.getPositions();
}

void GaussianCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

PropertyList GaussianCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList GaussianCalculator::possibleProperties() const {
  return Property::Energy | Property::Gradients | Property::AtomicCharges | Property::CoefficientMatrix |
         Property::ElectronicOccupation | Property::SuccessfulCalculation;
}

const Results& GaussianCalculator::calculate(std::string description) {
  applySettings();

  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

const Results& GaussianCalculator::calculateImpl(std::string description) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(calculationDirectory_);
  externalProgram.createWorkingDirectory();
  std::string inputFile = externalProgram.generateFullFilename(fileNameBase_ + ".inp");
  std::string checkpointFile = externalProgram.generateFullFilename(fileNameBase_ + ".chk");
  std::string outputFile = externalProgram.generateFullFilename(fileNameBase_ + ".out");

  GaussianInputFileCreator::createInputFile(inputFile, checkpointFile, atoms_, *settings_, requiredProperties_);

  // Check whether the given binary is valid
  if (!binaryIsValid()) {
    throw std::runtime_error("No or incorrect information about the binary path for Gaussian was given.");
  }
  externalProgram.executeCommand(gaussianExecutable_, inputFile, outputFile);

  GaussianOutputParser parser(outputFile);

  results_.set<Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Property::Energy)) {
    results_.set<Property::Energy>(parser.getEnergy());
  }
  if (requiredProperties_.containsSubSet(Property::Gradients)) {
    results_.set<Property::Gradients>(parser.getGradients());
  }
  if (requiredProperties_.containsSubSet(Property::AtomicCharges)) {
    results_.set<Property::AtomicCharges>(parser.getCM5Charges());
  }
  if (requiredProperties_.containsSubSet(Property::ElectronicOccupation) ||
      requiredProperties_.containsSubSet(Property::CoefficientMatrix)) {
    GaussianOrbitalParser orbitalParser(fileNameBase_, calculationDirectory_, gaussianDirectory_);
    if (requiredProperties_.containsSubSet(Property::CoefficientMatrix)) {
      results_.set<Property::CoefficientMatrix>(orbitalParser.getOrbitals());
    }
    if (requiredProperties_.containsSubSet(Property::ElectronicOccupation)) {
      results_.set<Property::ElectronicOccupation>(orbitalParser.getElectronicOccupation());
    }
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::ProgramName>("gaussian");
  if (SpinMode::Any == SpinModeInterpreter::getSpinModeFromString(settings_->getString(Utils::SettingsNames::spinMode))) {
    int multiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    auto spinMode = multiplicity == 1 ? SpinMode::Restricted : SpinMode::Unrestricted;
    settings_->modifyString(Utils::SettingsNames::spinMode, SpinModeInterpreter::getStringFromSpinMode(spinMode));
  }

  return results_;
}

void GaussianCalculator::setOrbitals(const MolecularOrbitals& mos) {
  GaussianOrbitalWriter orbitalWriter(mos);
  orbitalWriter.updateCheckpointFile(fileNameBase_, calculationDirectory_, gaussianDirectory_);
}

const Utils::Settings& GaussianCalculator::settings() const {
  return *settings_;
}

Utils::Settings& GaussianCalculator::settings() {
  return *settings_;
}

class GaussianStatesHandlingIsNotImplementedException : public std::exception {
  const char* what() const noexcept final {
    return "States handling is not available for Gaussian calculations.";
  }
};

std::shared_ptr<Core::State> GaussianCalculator::getState() const {
  throw GaussianStatesHandlingIsNotImplementedException();
}

void GaussianCalculator::loadState(std::shared_ptr<Core::State> /*state*/) {
  throw GaussianStatesHandlingIsNotImplementedException();
}

const Utils::Results& GaussianCalculator::results() const {
  return results_;
}

Utils::Results& GaussianCalculator::results() {
  return results_;
}

std::string GaussianCalculator::getFileNameBase() const {
  return fileNameBase_;
}

std::string GaussianCalculator::getCalculationDirectory() const {
  return calculationDirectory_;
}

bool GaussianCalculator::binaryIsValid() {
  if (binaryHasBeenChecked_) {
    return true;
  }
  if (gaussianExecutable_.empty()) {
    return false;
  }

  bp::ipstream out;
  std::error_code ec;
  bp::child c(gaussianExecutable_ + " non_existing_test_input_file", bp::std_out > out, bp::std_err > bp::null, ec);
  try {
    c.wait(); // Wait for the process to exit
  }
  catch (const std::exception& e) {
    // Since the g09 binary ends by Segmentation fault, old boost versions (e.g., 1.64.0) might throw an exception,
    // which is caught here.
  }

  std::string regexString = "non_existing_test_input_file\\.com";
  std::regex regex(regexString);
  std::smatch matches;

  // Get output
  std::string outputString;
  std::string line;
  while (std::getline(out, line)) {
    outputString += line;
  }

  bool isValid = std::regex_search(outputString, matches, regex);
  if (isValid) {
    binaryHasBeenChecked_ = true;
  }
  return isValid;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
