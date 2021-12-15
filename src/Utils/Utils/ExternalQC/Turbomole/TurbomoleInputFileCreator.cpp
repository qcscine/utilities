/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TurbomoleInputFileCreator.h"
#include "TurbomoleCalculatorSettings.h"
#include "TurbomoleHelper.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <math.h>
#include <boost/process.hpp>
#include <fstream>

namespace bp = boost::process;

namespace Scine {
namespace Utils {
namespace ExternalQC {

TurbomoleInputFileCreator::TurbomoleInputFileCreator(std::string& calculationDirectory,
                                                     std::string& turbomoleExecutableBase, TurbomoleFiles& files)
  : calculationDirectory_(calculationDirectory), turbomoleExecutableBase_(turbomoleExecutableBase), files_(files) {
}

void TurbomoleInputFileCreator::createInputFiles(const AtomCollection& atoms, const Settings& settings) {
  writeCoordFile(atoms);
  prepareDefineSession(settings, atoms);
  runDefine();
  checkAndUpdateControlFile(settings);
}

void TurbomoleInputFileCreator::writeCoordFile(const AtomCollection& atoms) {
  std::ofstream coordStream;
  coordStream.open(files_.coordFile);
  coordStream << "$coord\n";
  for (int i = 0; i < atoms.size(); i++) {
    std::string elementName = ElementInfo::symbol(atoms.at(i).getElementType());
    std::for_each(elementName.begin(), elementName.end(), [](char& c) { c = ::tolower(c); });
    coordStream << atoms.at(i).getPosition() << " " << elementName << std::endl;
  }
  coordStream << "$end";
  coordStream.close();
}

void TurbomoleInputFileCreator::prepareDefineSession(const Settings& settings, const AtomCollection& atoms) {
  // Check this before because define yields very strange output in case that charge and multiplicity don't match.
  checkValidityOfChargeAndMultiplicity(settings, atoms);

  std::ofstream out;
  // out.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  out.open(files_.defineInputFile);

  out << "\n"
      << "\n"
      << "a coord"
      << "\n"
      << "*\nno\n";

  auto basisSet = settings.getString(Utils::SettingsNames::basisSet);
  TurbomoleHelper helper(calculationDirectory_, turbomoleExecutableBase_);
  helper.mapBasisSetToTurbomoleStringRepresentation(basisSet);

  if (basisSet.empty())
    out << "\n*\neht\n\n";
  else
    out << "\nb all " << basisSet << "\n\n\n*\neht\n\n";

  out << settings.getInt(Scine::Utils::SettingsNames::molecularCharge) << "\n";

  auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings.getString(Utils::SettingsNames::spinMode));
  auto multiplicity = settings.getInt(Scine::Utils::SettingsNames::spinMultiplicity);
  // If spin mode is not set or set to "restricted", let turbomole handle spin mode automatically (default is RHF)
  if (spinMode == SpinMode::Any) {
    out << "\n\n\n";
  }
  else if (spinMode == SpinMode::Restricted) {
    if (multiplicity != 1) {
      throw std::logic_error("Specified restricted spin for multiplicity larger than 1.");
    }
    out << "\n\n\n";
  }
  else if (spinMode == SpinMode::Unrestricted) {
    int numElectrons = multiplicity - 1;
    if (numElectrons == 0)
      out << "no\ns\n*\n\n";
    else
      out << "no\nu " << numElectrons << "\n*\n\n";
  }
  else if (spinMode == SpinMode::RestrictedOpenShell) {
    throw std::logic_error("Spin mode not implemented in Turbomole!");
  }
  else
    throw std::logic_error("Specified unknown spin mode " + SpinModeInterpreter::getStringFromSpinMode(spinMode) +
                           " in settings."); // this should have been handled by settings

  out << "ri\non\n\n";
  // Method
  auto method = settings.getString(Utils::SettingsNames::method);
  std::istringstream iss(method);
  std::vector<std::string> methods((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
  // make sure that functional is lowercase only
  std::transform(std::begin(methods[0]), std::end(methods[0]), std::begin(methods[0]),
                 [](const auto c) { return std::tolower(c); });
  out << "dft\non\nfunc " << methods[0] << "\n\n";
  // Dispersion Correction
  if (methods.size() > 1) {
    auto vdWType = methods[1];
    std::for_each(vdWType.begin(), vdWType.end(), [](char& c) { c = ::toupper(c); });
    auto it = std::find(availableD3Params_.begin(), availableD3Params_.end(), vdWType);
    if (it - availableD3Params_.begin() == 0)
      out << "dsp\non\n\n";
    else if (it - availableD3Params_.begin() == 1)
      out << "dsp\nbj\n\n";
    else
      throw std::runtime_error("Invalid dispersion correction!");
  }
  auto maxScfIterations = settings.getInt(Utils::SettingsNames::maxScfIterations);

  out << "scf\niter\n" << std::to_string(maxScfIterations) << "\n";
  out << "\n*";
  out.close();
}

void TurbomoleInputFileCreator::runDefine() {
  TurbomoleHelper helper(calculationDirectory_, turbomoleExecutableBase_);
  // empty the control file in case it is present already
  helper.emptyFile(files_.controlFile);
  helper.execute("define", files_.defineInputFile);
}

void TurbomoleInputFileCreator::checkAndUpdateControlFile(const Settings& settings) {
  bool basisSetIsCorrect = false;
  bool functionalIsCorrect = false;

  double scfConvCriterion = settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion);
  auto scfDamping = settings.getBool(Utils::SettingsNames::scfDamping);
  auto scfOrbitalShift = settings.getDouble(SettingsNames::scfOrbitalShift);
  auto solvent = settings.getString(Utils::SettingsNames::solvent);

  if (!solvent.empty() && solvent != "none")
    addSolvation(settings);

  auto method = settings.getString(Utils::SettingsNames::method);
  std::istringstream iss(method);
  std::vector<std::string> methods((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
  std::for_each(methods[0].begin(), methods[0].end(), [](char& c) { c = ::tolower(c); });

  std::ifstream in;
  std::ofstream out;
  std::string controlBackupFile = NativeFilenames::combinePathSegments(calculationDirectory_, "control.bak");
  in.open(files_.controlFile);
  out.open(controlBackupFile);
  std::string line;

  auto basisSet = settings.getString(Utils::SettingsNames::basisSet);
  TurbomoleHelper helper(calculationDirectory_, turbomoleExecutableBase_);
  helper.mapBasisSetToTurbomoleStringRepresentation(basisSet);

  while (std::getline(in, line)) {
    if ((line.find("$scfdamp") != std::string::npos) && scfDamping)
      out << "$scfdamp   start=8.500  step=0.10  min=0.10"
          << "\n";
    else if ((line.find("scforbitalshift") != std::string::npos))
      out << "$scforbitalshift closedshell=" << scfOrbitalShift << "\n";
    else if (line.find("$end") != std::string::npos)
      out << "$pop loewdin wiberg\n$end";
    else if (line.find("$scfconv") != std::string::npos)
      out << "$scfconv " << int(round(-log10(scfConvCriterion))) << "\n";
    else
      out << line << "\n";

    if (line.find(basisSet) != std::string::npos)
      basisSetIsCorrect = true;
    else if (line.find(methods[0]) != std::string::npos)
      functionalIsCorrect = true;
  }

  in.close();
  out.close();
  std::rename(controlBackupFile.c_str(), files_.controlFile.c_str());

  if (!basisSetIsCorrect)
    throw std::runtime_error("Your specified basis set " + basisSet + " is invalid. Check the spelling");
  if (!functionalIsCorrect)
    throw std::runtime_error("Your specified DFT functional " + methods[0] + " is invalid. Check the spelling");
}

void TurbomoleInputFileCreator::addSolvation(const Settings& settings) {
  auto solvent = settings.getString(Utils::SettingsNames::solvent);

  std::for_each(solvent.begin(), solvent.end(), [](char& c) { c = ::tolower(c); });

  std::ofstream out;
  out.open(files_.solvationInputFile);

  std::unordered_map<std::string, double>::iterator it = availableSolventModels_.find(solvent);
  if (it != availableSolventModels_.end()) {
    out << it->second << "\n\n\n\n\n\n\n\n\n\n\n\n"
        << "r all b"
        << "\n"
        << "*"
        << "\n\n\n";
  }
  else
    throw std::runtime_error("The solvent '" + solvent + "' is currently not supported.");
  out.close();

  // run cosmoprep
  auto directory = bp::start_dir(calculationDirectory_);
  bp::ipstream stdout, stderr;
  std::string outputFile = NativeFilenames::combinePathSegments(calculationDirectory_, "COSMO.out");
  std::string executable = NativeFilenames::combinePathSegments(turbomoleExecutableBase_, "cosmoprep");
  TurbomoleHelper helper(calculationDirectory_, turbomoleExecutableBase_);
  helper.execute("cosmoprep", files_.solvationInputFile);
}

void TurbomoleInputFileCreator::checkValidityOfChargeAndMultiplicity(const Settings& settings, const AtomCollection& atoms) {
  int numberUnpairedElectrons = settings.getInt(Utils::SettingsNames::spinMultiplicity) - 1;
  int charge = settings.getInt(Utils::SettingsNames::molecularCharge);

  int numElectrons = 0;
  for (const auto& atom : atoms) {
    numElectrons += Utils::ElementInfo::Z(atom.getElementType());
  }
  // Subtract the charge from the total number of electrons
  numElectrons -= charge;

  bool numUnpairedElecIsEven = numberUnpairedElectrons % 2 == 0;
  bool numElectronsIsEven = numElectrons % 2 == 0;
  if ((numUnpairedElecIsEven && !numElectronsIsEven) || (!numUnpairedElecIsEven && numElectronsIsEven)) {
    throw std::logic_error("Invalid charge/multiplicity pair for the given system!");
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
