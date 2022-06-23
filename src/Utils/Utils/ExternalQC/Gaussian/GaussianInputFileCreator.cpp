/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GaussianInputFileCreator.h"
#include "GaussianCalculatorSettings.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <boost/filesystem.hpp>
#include <fstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

void GaussianInputFileCreator::createInputFile(const std::string& filename, const std::string& checkpointFilename,
                                               const AtomCollection& atoms, const Settings& settings,
                                               const PropertyList& requiredProperties) {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, checkpointFilename, atoms, settings, requiredProperties);
  fout.close();
  CalculationRoutines::checkValidityOfChargeAndMultiplicity(settings.getInt(Utils::SettingsNames::molecularCharge),
                                                            settings.getInt(Utils::SettingsNames::spinMultiplicity), atoms);
}

void GaussianInputFileCreator::createInputFile(std::ostream& out, const std::string& checkpointFilename,
                                               const AtomCollection& atoms, const Settings& settings,
                                               const PropertyList& requiredProperties) {
  printCalculationType(out, checkpointFilename, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void GaussianInputFileCreator::printCalculationType(std::ostream& out, const std::string& checkpointFilename,
                                                    const Settings& settings, const PropertyList& requiredProperties) {
  out << "%NProcShared=" << settings.getInt(Utils::SettingsNames::externalProgramNProcs) << std::endl;
  out << "%Mem=" << settings.getInt(Utils::SettingsNames::externalProgramMemory) << "MB" << std::endl;
  // Store checkpoint file if molecular orbitals are requested and/or the chk file is to be reused
  std::string scfGuess = settings.getString(SettingsNames::scfGuess);
  if (requiredProperties.containsSubSet(Property::CoefficientMatrix) ||
      requiredProperties.containsSubSet(Property::ElectronicOccupation) || scfGuess == "read" || scfGuess == "(only, read)") {
    out << "%chk=" + checkpointFilename << std::endl;
  }

  auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings.getString(Utils::SettingsNames::spinMode));
  std::string spinModeMethodAddon = ""; // translated version of our spinMode to the Gaussian abbreviation
  if (spinMode == SpinMode::Restricted) {
    spinModeMethodAddon = "R";
  }
  else if (spinMode == SpinMode::Unrestricted) {
    spinModeMethodAddon = "U";
  }
  else if (spinMode == SpinMode::RestrictedOpenShell) {
    spinModeMethodAddon = "RO";
  }
  // do nothing for any
  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(
      settings.getString(Scine::Utils::SettingsNames::method));
  out << "# " << spinModeMethodAddon << methodInput.first << "/" << settings.getString(Utils::SettingsNames::basisSet)
      << " " << translateDispersion(methodInput.second);
  // Scf convergence
  double tolE = settings.getDouble(Scine::Utils::SettingsNames::selfConsistenceCriterion);
  // Gaussian only takes powers of ten -> check for valid input -> check if log is integer
  double logTolE = std::log10(tolE);
  if ((trunc(logTolE) != logTolE)) {
    throw std::logic_error("The self consistence criterion must be a power of ten for Gaussian calculations.");
  }
  out << " SCF=(Conver=" + std::to_string(static_cast<int>(-1 * logTolE)) + ")";
  // scf guess
  if (scfGuess == "read" && !boost::filesystem::exists(checkpointFilename)) {
    scfGuess = "harris";
  }
  out << " guess=" + scfGuess;
  // solvation
  auto solvent = settings.getString(Utils::SettingsNames::solvent);
  auto solvation = settings.getString(Utils::SettingsNames::solvation);
  if (!solvent.empty()) {
    out << " SCRF=(" << solvation << ",Solvent=" << solvent << ")";
  }
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    out << " Force";
  }
  if (requiredProperties.containsSubSet(Property::AtomicCharges)) {
    out << " Pop=Hirshfeld";
  }
  out << std::endl << std::endl; // Empty line needed
}

void GaussianInputFileCreator::printTitle(std::ostream& out) {
  out << "# Gaussian calculation created by SCINE" << std::endl << std::endl;
}

void GaussianInputFileCreator::printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings) {
  out << settings.getInt(Scine::Utils::SettingsNames::molecularCharge) << " "
      << settings.getInt(Scine::Utils::SettingsNames::spinMultiplicity) << std::endl;
  for (const auto& a : atoms) {
    MolecularTrajectoryIO::writeXYZLine(out, a.getElementType(), a.getPosition());
  }
  out << std::endl;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
