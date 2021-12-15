/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/ExternalQC/Cp2k/Cp2kCalculator.h"
#include "Utils/ExternalQC/Cp2k/Cp2kCalculatorSettings.h"
#include "Utils/ExternalQC/Cp2k/Cp2kInputFileCreator.h"
#include "Utils/ExternalQC/Cp2k/Cp2kMainOutputParser.h"
#include "Utils/ExternalQC/Cp2k/Cp2kState.h"
#include "Utils/ExternalQC/Exceptions.h"
#include "Utils/ExternalQC/ExternalProgram.h"
#include "Utils/IO/FilesystemHelpers.h"
#include "Utils/IO/NativeFilenames.h"
#include "Utils/Properties/Thermochemistry/ThermochemistryCalculator.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include "Utils/Solvation/ImplicitSolvation.h"
#include <boost/exception/diagnostic_information.hpp>
#include <boost/process.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <regex>

namespace bp = boost::process;
namespace bfs = boost::filesystem;

namespace Scine {
namespace Utils {
namespace ExternalQC {

std::string Cp2kCalculator::name() const {
  return "CP2K";
}

bool Cp2kCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  if (std::getenv("CP2K_BINARY_PATH")) {
    return std::find(availableMethodFamilies_.begin(), availableMethodFamilies_.end(), methodFamily) !=
           availableMethodFamilies_.end();
  }

  return false;
}

Cp2kCalculator::Cp2kCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<Cp2kCalculatorSettings>();

  // Get Cp2k's binary path from an environment variable
  if (const char* pathPtr = std::getenv("CP2K_BINARY_PATH")) {
    cp2kExecutable_ = std::string(pathPtr);
  }

  applySettings();
}

Cp2kCalculator::Cp2kCalculator(const Cp2kCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  applySettings();
  this->setStructure(rhs.atoms_);
  this->results() = rhs.results();
  this->cp2kExecutable_ = rhs.cp2kExecutable_;
  this->binaryHasBeenChecked_ = rhs.binaryHasBeenChecked_;
}

void Cp2kCalculator::applySettings() {
  if (settings_->valid()) {
    fileNameBase_ = settings_->getString(SettingsNames::cp2kFilenameBase);
    baseWorkingDirectory_ = settings_->getString(SettingsNames::baseWorkingDirectory);
    // throws error for wrong input and updates 'any' entries */
    // information if solvation is performed is deduced from settings from InputFileCreator and therefore discarded here
    Solvation::ImplicitSolvation::solvationNeededAndPossible(availableSolvationModels_, *settings_);
    if (settings_->getDouble(Utils::SettingsNames::electronicTemperature) > 0.0 &&
        settings_->getInt(SettingsNames::additionalMos) == 0) {
      throw Core::InitializationException("Specified non-zero electronic temperature, "
                                          "but no additional molecular orbitals!");
    }
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

void Cp2kCalculator::setStructure(const Utils::AtomCollection& structure) {
  applySettings();
  atoms_ = structure;
  calculationDirectory_ = createNameForCalculationDirectory();
  results_ = Results{};
}

std::unique_ptr<Utils::AtomCollection> Cp2kCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(atoms_);
}

void Cp2kCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  atoms_.setPositions(std::move(newPositions));
  results_ = Results{};
}

const Utils::PositionCollection& Cp2kCalculator::getPositions() const {
  return atoms_.getPositions();
}

void Cp2kCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
  // complete properties for thermochemistry
  if (requiredProperties.containsSubSet(Property::Thermochemistry)) {
    if (!requiredProperties.containsSubSet(Property::Hessian)) {
      requiredProperties_.addProperty(Property::Hessian);
    }
    if (!requiredProperties.containsSubSet(Property::Energy)) {
      requiredProperties_.addProperty(Property::Energy);
    }
  }
}

PropertyList Cp2kCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

PropertyList Cp2kCalculator::possibleProperties() const {
  return Property::Energy | Property::Gradients | Property::Hessian | Property::AtomicCharges |
         Property::DensityMatrix | Property::GridOccupation | Property::Thermochemistry | Property::OverlapMatrix |
         Property::AOtoAtomMapping | Property::BondOrderMatrix;
}

const Results& Cp2kCalculator::calculate(std::string description) {
  applySettings();
  // all properties that cannot be parsed from a CP2K Hessian calculation
  std::vector<Property> nonHessProperties = {Property::BondOrderMatrix, Property::DensityMatrix,
                                             Property::OverlapMatrix, Property::GridOccupation, Property::AtomicCharges};
  bool splitCalcNecessary = requiredProperties_.containsSubSet(Property::Hessian) &&
                            std::any_of(nonHessProperties.begin(), nonHessProperties.end(),
                                        [&](const auto& prop) { return requiredProperties_.containsSubSet(prop); });
  if (!splitCalcNecessary) {
    return calculateImplDeleteGuard(description);
  }
  // safe actually wanted properties
  PropertyList actualRequiredProperties = requiredProperties_;
  bool thermoToo = actualRequiredProperties.containsSubSet(Property::Thermochemistry);
  // fill in wanted properties possible without Hessian calculation
  requiredProperties_ = Property::Energy | Property::Gradients;
  for (const auto& prop : nonHessProperties) {
    if (actualRequiredProperties.containsSubSet(prop)) {
      requiredProperties_.addProperty(prop);
    }
  }
  results_ = calculateImplDeleteGuard(description);
  // fill properties with Hessian necessary
  requiredProperties_ = Property::Hessian;
  if (thermoToo) {
    requiredProperties_.addProperty(Property::Thermochemistry);
  }
  auto hessResults = calculateImplDeleteGuard(description);
  // transfer results and properties
  results_.set<Property::Hessian>(hessResults.get<Property::Hessian>());
  if (thermoToo) {
    results_.set<Property::Thermochemistry>(hessResults.get<Property::Thermochemistry>());
  }
  requiredProperties_ = actualRequiredProperties;
  return results_;
}

const Results& Cp2kCalculator::calculateImplDeleteGuard(const std::string& description) {
  try {
    return calculateImpl(description);
  }
  catch (...) {
    // Delete all of the .bak files in the calculation directory if requested
    if (settings_->getBool(SettingsNames::deleteTemporaryFiles)) {
      deleteTemporaryFiles();
    }
    throw Core::UnsuccessfulCalculationException(boost::current_exception_diagnostic_information());
  }
}

const Results& Cp2kCalculator::calculateImpl(const std::string& description) {
  ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(calculationDirectory_);
  externalProgram.createWorkingDirectory();
  std::string inputFile = externalProgram.generateFullFilename(fileNameBase_ + ".inp");
  std::string outputFile = externalProgram.generateFullFilename(fileNameBase_ + ".out");
  std::string additionalOutputFile =
      externalProgram.generateFullFilename(settings_->getString(SettingsNames::additionalOutputFile) + "-1_0.Log");

  auto inputCreator = Cp2kInputFileCreator();
  inputCreator.createInputFile(inputFile, atoms_, *settings_, requiredProperties_, fileNameBase_);

  // Check whether the given binary is valid
  if (!binaryIsValid()) {
    throw std::runtime_error("No or incorrect information about the binary path for Cp2k was given.");
  }

  // Delete output file before the calculation if it already exists.
  // Necessary since otherwise boost just pipes into the old output file which may lead to problems.
  bfs::remove(outputFile);
  bfs::remove(additionalOutputFile);

  // Execute Cp2k command
  int nCores = settings_->getInt(Utils::SettingsNames::externalProgramNProcs);
  std::string command = (nCores == 1) ? cp2kExecutable_ : "mpirun -np " + std::to_string(nCores) + " " + cp2kExecutable_;
  if (nCores > 1 && !multiProcIsPossible_) {
    this->getLog().warning << "Warning: Requested multiple cores, but 'mpirun' is not available on your system."
                           << Core::Log::nl << "Executing CP2K with a single core." << Core::Log::nl;
    settings_->modifyInt(Utils::SettingsNames::externalProgramNProcs, 1);
    command = cp2kExecutable_;
  }
  command += " -o " + outputFile + " " + inputFile;
  externalProgram.executeCommand(command);

  auto parser = (bfs::exists(additionalOutputFile)) ? Cp2kMainOutputParser(outputFile, additionalOutputFile)
                                                    : Cp2kMainOutputParser(outputFile);
  try {
    parser.checkForErrors();
  }
  catch (const ScfNotConvergedError& e) {
    if (!settings_->getBool(SettingsNames::allowUnconvergedScf)) {
      throw std::runtime_error(e.what());
    }
  }
  auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings_->getString(Utils::SettingsNames::spinMode));
  if (spinMode == SpinMode::Any) {
    int multiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    spinMode = (multiplicity == 1) ? SpinMode::Restricted : SpinMode::Unrestricted;
    settings_->modifyString(Utils::SettingsNames::spinMode, SpinModeInterpreter::getStringFromSpinMode(spinMode));
  }

  results_.set<Property::Description>(description);
  if (requiredProperties_.containsSubSet(Property::Energy)) {
    results_.set<Property::Energy>(parser.getEnergy());
  }
  if (requiredProperties_.containsSubSet(Property::Gradients)) {
    results_.set<Property::Gradients>(parser.getGradients());
  }
  if (requiredProperties_.containsSubSet(Property::AtomicCharges)) {
    results_.set<Property::AtomicCharges>(parser.getHirshfeldCharges());
  }
  if (requiredProperties_.containsSubSet(Property::BondOrderMatrix)) {
    results_.set<Property::BondOrderMatrix>(parser.getBondOrders(atoms_.getElements(), spinMode));
  }
  if (requiredProperties_.containsSubSet(Property::GridOccupation)) {
    results_.set<Property::GridOccupation>(parser.getGridCounts());
  }
  if (requiredProperties_.containsSubSet(Property::DensityMatrix)) {
    results_.set<Property::DensityMatrix>(parser.getDensityMatrix(spinMode));
  }
  if (requiredProperties_.containsSubSet(Property::OverlapMatrix)) {
    results_.set<Property::OverlapMatrix>(parser.getOverlapMatrix());
  }
  if (requiredProperties_.containsSubSet(Property::AOtoAtomMapping)) {
    results_.set<Property::AOtoAtomMapping>(parser.getAtomAoIndex(atoms_.getElements()));
  }
  if (requiredProperties_.containsSubSet(Property::Hessian)) {
    results_.set<Property::Hessian>(parser.getHessian());
  }
  if (requiredProperties_.containsSubSet(Property::Thermochemistry)) {
    auto thermoCalc = ThermochemistryCalculator(results_.get<Property::Hessian>(), atoms_,
                                                settings_->getInt(Utils::SettingsNames::spinMultiplicity),
                                                results_.get<Property::Energy>());
    thermoCalc.setMolecularSymmetryNumber(parser.getSymmetryNumber());
    thermoCalc.setTemperature(settings_->getDouble(Utils::SettingsNames::temperature));
    auto thermochemistry = thermoCalc.calculate();
    results_.set<Property::Thermochemistry>(thermochemistry);
  }
  results_.set<Property::SuccessfulCalculation>(true);
  results_.set<Property::ProgramName>("cp2k");

  return results_;
}

std::string Cp2kCalculator::createNameForCalculationDirectory() {
  boost::uuids::uuid uuid = boost::uuids::random_generator()();
  std::string uuidString = boost::uuids::to_string(uuid);

  auto directoryName = NativeFilenames::combinePathSegments(baseWorkingDirectory_, uuidString);
  return NativeFilenames::addTrailingSeparator(directoryName);
}

const Utils::Settings& Cp2kCalculator::settings() const {
  return *settings_;
}

Utils::Settings& Cp2kCalculator::settings() {
  return *settings_;
}

void Cp2kCalculator::copyBackupFile(const std::string& from, const std::string& to) const {
  const auto& workingDirectory = this->getCalculationDirectory();
  auto fromFile = NativeFilenames::combinePathSegments(workingDirectory, from + "-RESTART.wfn");
  auto toFile = NativeFilenames::combinePathSegments(workingDirectory, to + "-RESTART.wfn");

  try {
    FilesystemHelpers::copyFile(fromFile, toFile);
  }
  catch (std::runtime_error& e) {
    throw Cp2kStateSavingException();
  }
}

std::shared_ptr<Core::State> Cp2kCalculator::getState() const {
  auto ret = std::make_shared<Cp2kState>(this->getCalculationDirectory());
  this->copyBackupFile(this->getFileNameBase(), ret->stateIdentifier);
  return ret;
}

void Cp2kCalculator::loadState(std::shared_ptr<Core::State> state) {
  auto cp2kState = std::dynamic_pointer_cast<Cp2kState>(state);
  this->copyBackupFile(cp2kState->stateIdentifier, this->getFileNameBase());
}

const Utils::Results& Cp2kCalculator::results() const {
  return results_;
}

Utils::Results& Cp2kCalculator::results() {
  return results_;
}

std::string Cp2kCalculator::getFileNameBase() const {
  return fileNameBase_;
}

std::string Cp2kCalculator::getCalculationDirectory() const {
  return calculationDirectory_;
}

bool Cp2kCalculator::binaryIsValid() {
  if (binaryHasBeenChecked_) {
    return true;
  }

  if (cp2kExecutable_.empty()) {
    return false;
  }

  bp::ipstream out;
  std::error_code ec;
  bp::child c(cp2kExecutable_, bp::std_out > out, bp::std_err > bp::null, ec);
  c.wait(); // Wait for the process to exit

  std::string regexString = "The following options can be used"; // string in help message for no input file given
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
    checkMpirun();
  }
  return isValid;
}

void Cp2kCalculator::checkMpirun() {
  bp::ipstream err;
  std::error_code ec;
  bp::child c("mpirun", bp::std_out > bp::null, bp::std_err > err, ec);
  c.wait(); // Wait for the process to exit

  std::regex regex("mpirun could not find anything to do."); // string in help message if nothing after mpirun
  std::smatch matches;

  // Get output
  std::string outputString;
  std::string line;
  while (std::getline(err, line)) {
    outputString += line;
  }

  bool isValid = std::regex_search(outputString, matches, regex);
  multiProcIsPossible_ = isValid;
}

void Cp2kCalculator::deleteTemporaryFiles() {
  bfs::path p(calculationDirectory_);
  const std::regex filter(".+\\.bak.{0,}");
  if (bfs::exists(p) && bfs::is_directory(p)) {
    bfs::directory_iterator end;
    for (bfs::directory_iterator it(p); it != end; ++it) {
      try {
        if (bfs::is_regular_file(it->status())) {
          std::smatch match;
          if (std::regex_search(it->path().filename().string(), match, filter)) {
            bfs::remove(it->path());
          }
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
