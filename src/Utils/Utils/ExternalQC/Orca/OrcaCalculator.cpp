/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/ExternalQC/Orca/OrcaCalculator.h"
#include "Utils/ExternalQC/ExternalProgram.h"
#include "Utils/ExternalQC/Orca/MoessbauerParameters.h"
#include "Utils/ExternalQC/Orca/OrcaCalculatorSettings.h"
#include "Utils/ExternalQC/Orca/OrcaHessianOutputParser.h"
#include "Utils/ExternalQC/Orca/OrcaInputFileCreator.h"
#include "Utils/ExternalQC/Orca/OrcaMainOutputParser.h"
#include "Utils/ExternalQC/Orca/OrcaPointChargesGradientsFileParser.h"
#include "Utils/ExternalQC/Orca/OrcaState.h"
#include "Utils/IO/FilesystemHelpers.h"
#include "Utils/IO/NativeFilenames.h"
#include "Utils/Properties/Thermochemistry/ThermochemistryCalculator.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include "Utils/Solvation/ImplicitSolvation.h"
#include <boost/process.hpp>
#include <regex>

namespace bp = boost::process;
namespace bfs = boost::filesystem;

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string OrcaCalculator::name() const {
  return "ORCA";
}

bool OrcaCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  if (std::getenv("ORCA_BINARY_PATH")) {
    return std::find(availableMethodFamilies_.begin(), availableMethodFamilies_.end(), methodFamily) !=
           availableMethodFamilies_.end();
  }

  return false;
}

OrcaCalculator::OrcaCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<OrcaCalculatorSettings>();

  // Get Orca's binary path from an environment variable
  if (const char* pathPtr = std::getenv("ORCA_BINARY_PATH")) {
    orcaExecutable_ = std::string(pathPtr);
  }

  applySettings();
}

OrcaCalculator::OrcaCalculator(const OrcaCalculator& rhs) : CloneInterface(rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  this->setLog(rhs.getLog());
  applySettings();
  atoms_ = rhs.atoms_;
  calculationDirectory_ = NativeFilenames::createRandomDirectoryName(baseWorkingDirectory_);
  results_ = rhs.results();
  orcaExecutable_ = rhs.orcaExecutable_;
  binaryHasBeenChecked_ = rhs.binaryHasBeenChecked_;
}

void OrcaCalculator::applySettings() {
  if (settings_->valid()) {
    if (settings_->getDouble(Utils::SettingsNames::electronicTemperature) > 0.0) {
      throw Core::InitializationException(
          "ORCA calculations with an electronic temperature above 0.0 K are not supported.");
    }
    fileNameBase_ = settings_->getString(SettingsNames::orcaFilenameBase);
    baseWorkingDirectory_ = settings_->getString(SettingsNames::baseWorkingDirectory);
    // throws error for wrong input and updates 'any' entries */
    // information if solvation is performed is deduced from settings from InputFileCreator and therefore discarded here
    Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels_, *settings_);
    if (!(settings_->getBool(SettingsNames::enforceScfCriterion)) &&
        (requiredProperties_.containsSubSet(Property::Gradients) || requiredProperties_.containsSubSet(Property::Hessian)) &&
        settings_->getDouble(Utils::SettingsNames::selfConsistenceCriterion) > 1e-8) {
      settings_->modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
      this->getLog().warning << "Warning: Energy accuracy was increased to 1e-8 to ensure valid gradients/hessian as "
                                "recommended by ORCA developers."
                             << Core::Log::nl;
    }
    if (requiredProperties_.containsSubSet(Property::Gradients) &&
        std::find(numGradMethods_.begin(), numGradMethods_.end(), settings_->getString(Utils::SettingsNames::method)) !=
            numGradMethods_.end()) {
      settings_->modifyString(SettingsNames::gradientCalculationType, "numerical");
      this->getLog().output << "Calculating gradients numerically." << Core::Log::nl;
    }
    if (requiredProperties_.containsSubSet(Property::Hessian) &&
        std::find(numFreqMethods_.begin(), numFreqMethods_.end(), settings_->getString(Utils::SettingsNames::method)) !=
            numFreqMethods_.end()) {
      settings_->modifyString(SettingsNames::hessianCalculationType, "numerical");
      this->getLog().output << "Calculating Hessian numerically." << Core::Log::nl;
    }
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

void OrcaCalculator::setStructure(const Utils::AtomCollection& structure) {
  applySettings();
  atoms_ = structure;
  calculationDirectory_ = NativeFilenames::createRandomDirectoryName(baseWorkingDirectory_);
  results_ = Results{};
}

std::unique_ptr<Utils::AtomCollection> OrcaCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(atoms_);
}

void OrcaCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  atoms_.setPositions(std::move(newPositions));
  results_ = Results{};
}

const Utils::PositionCollection& OrcaCalculator::getPositions() const {
  return atoms_.getPositions();
}

void OrcaCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
  if (requiredProperties_.containsSubSet(Property::Thermochemistry)) {
    requiredProperties_.addProperty(Property::Hessian);
  }
}

PropertyList OrcaCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList OrcaCalculator::possibleProperties() const {
  return Property::Energy | Property::Gradients | Property::Hessian | Property::BondOrderMatrix |
         Property::Thermochemistry | Property::AtomicCharges | Property::PointChargesGradients |
         Property::MoessbauerParameter | Property::OrbitalEnergies | Property::SuccessfulCalculation;
}

const Results& OrcaCalculator::calculate(std::string description) {
  applySettings();
  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    // Delete all of the .tmp files in the calculation directory if requested
    if (settings_->getBool(SettingsNames::deleteTemporaryFiles)) {
      deleteTemporaryFiles();
    }
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

const Results& OrcaCalculator::calculateImpl(std::string description) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(calculationDirectory_);
  externalProgram.createWorkingDirectory();
  std::string inputFile = externalProgram.generateFullFilename(fileNameBase_ + ".inp");
  std::string outputFile = externalProgram.generateFullFilename(fileNameBase_ + ".out");

  OrcaInputFileCreator::createInputFile(inputFile, atoms_, *settings_, requiredProperties_);

  // Check whether the given binary is valid
  if (!binaryIsValid()) {
    throw std::runtime_error("No or incorrect information about the binary path for Orca was given.");
  }

  // Delete output file before the calculation if it already exists.
  // Necessary since otherwise boost just pipes into the old output file which may lead to problems.
  bfs::remove(outputFile);

  // Execute Orca command
  externalProgram.executeCommand(orcaExecutable_ + " " + inputFile, outputFile);

  OrcaMainOutputParser parser(outputFile);
  parser.checkForErrors();

  results_.set<Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Property::Energy)) {
    results_.set<Property::Energy>(parser.getEnergy());
  }
  if (requiredProperties_.containsSubSet(Property::Gradients)) {
    results_.set<Property::Gradients>(parser.getGradients());
  }
  if (requiredProperties_.containsSubSet(Property::Hessian)) {
    std::string hessianFile = externalProgram.generateFullFilename(fileNameBase_ + ".hess");
    results_.set<Property::Hessian>(OrcaHessianOutputParser::getHessian(hessianFile));
  }
  if (requiredProperties_.containsSubSet(Property::BondOrderMatrix)) {
    results_.set<Property::BondOrderMatrix>(parser.getBondOrders());
  }
  if (requiredProperties_.containsSubSet(Property::AtomicCharges)) {
    results_.set<Property::AtomicCharges>(parser.getHirshfeldCharges());
  }
  if (requiredProperties_.containsSubSet(Property::Thermochemistry)) {
    ThermochemicalContainer thermochemicalContainer;
    thermochemicalContainer.symmetryNumber = parser.getSymmetryNumber();
    thermochemicalContainer.enthalpy = parser.getEnthalpy();
    thermochemicalContainer.entropy = parser.getEntropy();
    thermochemicalContainer.zeroPointVibrationalEnergy = parser.getZeroPointVibrationalEnergy();
    thermochemicalContainer.gibbsFreeEnergy = parser.getGibbsFreeEnergy();
    thermochemicalContainer.heatCapacityP = std::numeric_limits<double>::quiet_NaN();
    thermochemicalContainer.heatCapacityV = std::numeric_limits<double>::quiet_NaN();

    ThermochemicalComponentsContainer thermochemistry;
    thermochemistry.overall = thermochemicalContainer;
    results_.set<Property::Thermochemistry>(thermochemistry);
  }
  if (requiredProperties_.containsSubSet(Property::PointChargesGradients)) {
    std::string pcGradientsFile = externalProgram.generateFullFilename(fileNameBase_ + ".pcgrad");
    OrcaPointChargesGradientsFileParser pcParser(pcGradientsFile);
    results_.set<Property::PointChargesGradients>(pcParser.getPointChargesGradients());
  }
  if (requiredProperties_.containsSubSet(Property::MoessbauerParameter)) {
    Moessbauer::MoessbauerParameterContainer moessbauerContainer;
    int numIrons = Moessbauer::determineNumIrons(atoms_);
    moessbauerContainer.etas = parser.getMoessbauerAsymmetryParameter(numIrons);
    moessbauerContainer.quadrupoleSplittings = parser.getMoessbauerQuadrupoleSplittings(numIrons);
    moessbauerContainer.densities = parser.getMoessbauerIronElectronDensities(numIrons);
    results_.set<Property::MoessbauerParameter>(moessbauerContainer);
  }
  if (requiredProperties_.containsSubSet(Property::OrbitalEnergies)) {
    results_.set<Property::OrbitalEnergies>(parser.getOrbitalEnergies());
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::ProgramName>("orca");
  if (SpinMode::Any == SpinModeInterpreter::getSpinModeFromString(settings_->getString(Utils::SettingsNames::spinMode))) {
    int multiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    auto spinMode = multiplicity == 1 ? SpinMode::Restricted : SpinMode::Unrestricted;
    settings_->modifyString(Utils::SettingsNames::spinMode, SpinModeInterpreter::getStringFromSpinMode(spinMode));
  }

  return results_;
}

const Utils::Settings& OrcaCalculator::settings() const {
  return *settings_;
}

Utils::Settings& OrcaCalculator::settings() {
  return *settings_;
}

void OrcaCalculator::copyBackupFile(const std::string& from, const std::string& to) const {
  const auto& workingDirectory = this->getCalculationDirectory();
  auto fromFile = NativeFilenames::combinePathSegments(workingDirectory, from + ".gbw");
  auto toFile = NativeFilenames::combinePathSegments(workingDirectory, to + ".gbw");

  try {
    FilesystemHelpers::copyFile(fromFile, toFile);
  }
  catch (std::runtime_error& e) {
    throw OrcaStateSavingException();
  }
}

std::shared_ptr<Core::State> OrcaCalculator::getState() const {
  auto ret = std::make_shared<OrcaState>(this->getCalculationDirectory());
  this->copyBackupFile(this->getFileNameBase(), ret->stateIdentifier);
  return ret;
}

void OrcaCalculator::loadState(std::shared_ptr<Core::State> state) {
  auto orcaState = std::dynamic_pointer_cast<OrcaState>(state);
  this->copyBackupFile(orcaState->stateIdentifier, this->getFileNameBase());
}

const Utils::Results& OrcaCalculator::results() const {
  return results_;
}

Utils::Results& OrcaCalculator::results() {
  return results_;
}

std::string OrcaCalculator::getFileNameBase() const {
  return fileNameBase_;
}

std::string OrcaCalculator::getCalculationDirectory() const {
  return calculationDirectory_;
}

bool OrcaCalculator::binaryIsValid() {
  if (binaryHasBeenChecked_) {
    return true;
  }

  if (orcaExecutable_.empty()) {
    return false;
  }

  bp::ipstream out;
  std::error_code ec;
  bp::child c(orcaExecutable_, bp::std_out > out, bp::std_err > bp::null, ec);
  c.wait(); // Wait for the process to exit

  std::string regexString = "ORCA TEST.INP";
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

void OrcaCalculator::deleteTemporaryFiles() {
  bfs::path p(calculationDirectory_);
  if (bfs::exists(p) && bfs::is_directory(p)) {
    bfs::directory_iterator end;
    for (bfs::directory_iterator it(p); it != end; ++it) {
      try {
        if (bfs::is_regular_file(it->status()) && (it->path().extension() == ".tmp")) {
          bfs::remove(it->path());
        }
      }
      catch (const std::exception& e) {
      }
    }
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
