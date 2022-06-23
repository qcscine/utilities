/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaInputFileCreator.h"
#include "OrcaCalculatorSettings.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <fstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

void OrcaInputFileCreator::createInputFile(const std::string& filename, const AtomCollection& atoms,
                                           const Settings& settings, const PropertyList& requiredProperties) {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, atoms, settings, requiredProperties);
  fout.close();
  CalculationRoutines::checkValidityOfChargeAndMultiplicity(settings.getInt(Utils::SettingsNames::molecularCharge),
                                                            settings.getInt(Utils::SettingsNames::spinMultiplicity), atoms);
}

void OrcaInputFileCreator::createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                           const PropertyList& requiredProperties) {
  printCalculationType(out, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void OrcaInputFileCreator::printCalculationType(std::ostream& out, const Settings& settings,
                                                const PropertyList& requiredProperties) {
  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(
      settings.getString(Scine::Utils::SettingsNames::method));
  out << "! " << methodInput.first << " " << methodInput.second << " "
      << settings.getString(Utils::SettingsNames::basisSet) << std::endl;

  auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings.getString(Utils::SettingsNames::spinMode));
  if (spinMode == SpinMode::Unrestricted) {
    out << "! UHF" << std::endl; // keyword UKS also exists, but they are interchangeable, same for restricted
  }
  else if (spinMode == SpinMode::Restricted) {
    out << "! RHF AllowRHF" << std::endl; // AllowRHF enforces restricted
  }
  else if (spinMode == SpinMode::RestrictedOpenShell) {
    out << "! ROHF" << std::endl;
  }
  // do nothing if any, which is default in ORCA

  auto scfDamping = settings.getBool(Utils::SettingsNames::scfDamping);
  if (scfDamping) {
    out << "! SlowConv" << std::endl;
  }

  auto solvent = settings.getString(Utils::SettingsNames::solvent);
  if (!solvent.empty() && solvent != "none") {
    out << "! CPCM(" << solvent << ")" << std::endl;
  }
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    out << "! EnGrad TightSCF" << std::endl;
  }
  if (requiredProperties.containsSubSet(Property::Hessian)) {
    std::string freqType = (settings.getString(SettingsNames::hessianCalculationType) == "analytical") ? "AnFreq" : "NumFreq";
    out << "! " << freqType << std::endl;
  }
  auto specialOption = settings.getString(Scine::Utils::ExternalQC::SettingsNames::specialOption);
  if (!specialOption.empty()) {
    out << "! " << specialOption << std::endl;
  }

  // Number of processes and memory per core
  int numProcs = settings.getInt(Utils::SettingsNames::externalProgramNProcs);
  out << "%maxcore " << settings.getInt(Utils::SettingsNames::externalProgramMemory) / numProcs << std::endl;
  if (numProcs != 1) {
    out << "%pal\nnprocs " << numProcs << "\nend" << std::endl;
  }
  if (!solvent.empty() && solvent != "none") {
    out << "%cpcm ndiv 6" << std::endl;
    if (settings.getString(Utils::SettingsNames::solvation) == "smd") {
      out << "smd true\nSMDsolvent \"" << solvent << "\"" << std::endl;
    }
    out << "end" << std::endl;
  }

  // Print Mayer bond orders and/or Hirshfeld charges if required
  bool bondOrdersRequired = requiredProperties.containsSubSet(Utils::Property::BondOrderMatrix);
  bool atomicChargesRequired = requiredProperties.containsSubSet(Utils::Property::AtomicCharges);
  if (bondOrdersRequired && atomicChargesRequired) {
    out << "%output\nprint[P_Mayer] 1\nprint[P_Hirshfeld] 1\nend" << std::endl;
  }
  else if (bondOrdersRequired) {
    out << "%output\nprint[P_Mayer] 1\nend" << std::endl;
  }
  else if (atomicChargesRequired) {
    out << "%output\nprint[P_Hirshfeld] 1\nend" << std::endl;
  }
  // Set temperature if thermochemistry is to be calculated
  if (requiredProperties.containsSubSet(Utils::Property::Thermochemistry)) {
    out << "%freq\nTemp " << settings.getDouble(Utils::SettingsNames::temperature) << "\nend" << std::endl;
  }
  // Scf convergence and Scf max iterations
  out << "%SCF\nTolE " << settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) << std::endl
      << "MaxIter " << settings.getInt(Scine::Utils::SettingsNames::maxScfIterations) << "\nend" << std::endl;
  // Write name of point charges file if it was set
  auto pointChargesFile = settings.getString(SettingsNames::pointChargesFile);
  if (!pointChargesFile.empty()) {
    out << "%pointcharges \"" << pointChargesFile << "\"" << std::endl;
  }
}

void OrcaInputFileCreator::printTitle(std::ostream& out) {
  out << "# Orca calculation created by SCINE" << std::endl;
}

void OrcaInputFileCreator::printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings) {
  out << "*xyz " << settings.getInt(Scine::Utils::SettingsNames::molecularCharge) << " "
      << settings.getInt(Scine::Utils::SettingsNames::spinMultiplicity) << std::endl;
  for (const auto& a : atoms) {
    MolecularTrajectoryIO::writeXYZLine(out, a.getElementType(), a.getPosition());
  }
  out << "*" << std::endl;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
