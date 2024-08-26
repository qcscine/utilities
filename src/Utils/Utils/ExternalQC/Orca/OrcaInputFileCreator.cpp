/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaInputFileCreator.h"
#include "OrcaCalculatorSettings.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>

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
  printCalculationType(out, atoms, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void OrcaInputFileCreator::printCalculationType(std::ostream& out, const AtomCollection& atoms,
                                                const Settings& settings, const PropertyList& requiredProperties) {
  const std::string basisSet = settings.getString(Utils::SettingsNames::basisSet);
  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(
      settings.getString(Scine::Utils::SettingsNames::method));
  out << "! " << methodInput.first << " " << methodInput.second << " " << basisSet << std::endl;

  if ((boost::to_upper_copy<std::string>(methodInput.first).find("DLPNO") != std::string::npos) ||
      (boost::to_upper_copy<std::string>(methodInput.first).find("CC") != std::string::npos)) {
    auto auxCBasisSet = settings.getString(Scine::Utils::ExternalQC::SettingsNames::orcaAuxCBasisSet);
    if (!auxCBasisSet.empty()) {
      out << "! " << auxCBasisSet << "/C" << std::endl;
    }
    else {
      out << "! " << basisSet << "/C" << std::endl;
    }
  }

  if (boost::to_upper_copy<std::string>(methodInput.first).find("F12") != std::string::npos) {
    auto cabsBasisSet = settings.getString(Scine::Utils::ExternalQC::SettingsNames::orcaCabsBasisSet);
    if (!cabsBasisSet.empty()) {
      out << "! " << cabsBasisSet << std::endl;
    }
    else {
      out << "! " << basisSet << "-CABS" << std::endl;
    }
  }

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
  double epsilon = -1.0;
  double probeRadius = -1.0;
  if (!solvent.empty() && solvent != "none") {
    if (solvent.find("user_defined") != std::string::npos) {
      auto epsRadPair = interpretAsUserDefinedImplicitSolvation(solvent);
      epsilon = std::get<0>(epsRadPair);
      probeRadius = std::get<1>(epsRadPair);
    }
    else {
      out << "! CPCM(" << solvent << ")" << std::endl;
    }
  }
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    std::string gradType =
        (settings.getString(SettingsNames::gradientCalculationType) == "analytical") ? "EnGrad TightSCF" : "NumGrad";
    out << "! " << gradType << std::endl;
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
    if (solvent.find("user_defined") != std::string::npos) {
      out << "epsilon " << epsilon << std::endl;
      out << "rsolv " << probeRadius << std::endl;
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
      << "MaxIter " << settings.getInt(Scine::Utils::SettingsNames::maxScfIterations);
  if (settings.getBool(Utils::ExternalQC::SettingsNames::performBrokenSymmetryCalculation)) {
    // First check if the
    int numUnpairedElec = settings.getInt(Utils::SettingsNames::spinMultiplicity) - 1;
    int numUnpairedElecBeforeFlip = settings.getInt(Utils::ExternalQC::SettingsNames::initialSpinMultiplicity) - 1;
    bool numUnpairedElecBeforeFlipIsEven = numUnpairedElecBeforeFlip % 2 == 0;
    bool numUnpairedElecAfterFlipIsEven = numUnpairedElec % 2 == 0;

    // Check if combination of spin multiplicies is logic
    if ((numUnpairedElecBeforeFlipIsEven && !numUnpairedElecAfterFlipIsEven) ||
        (!numUnpairedElecBeforeFlipIsEven && numUnpairedElecAfterFlipIsEven)) {
      throw std::logic_error(
          "The final spin multiplicity cannot be generated from the initial spin multiplicity via spin flip!");
    }

    if (settings.getInt(Utils::ExternalQC::SettingsNames::initialSpinMultiplicity) == -1) {
      std::string errorString =
          "Please set both the initial (setting name: " +
          std::string(Scine::Utils::ExternalQC::SettingsNames::initialSpinMultiplicity) +
          ") and the final spin multiplicity (setting name: " + std::string(Scine::Utils::SettingsNames::spinMultiplicity) +
          ") if you want to perform a broken-symmetry calculation.";
      throw std::runtime_error(errorString);
    }
    if (settings.getIntList(Utils::ExternalQC::SettingsNames::spinFlipSites).empty()) {
      throw std::runtime_error("Please set the atom indices of all sites at which spin density should be flipped after "
                               "converging to the high-spin solution!");
    }
    out << "Flipspin ";
    auto spinFlipSites = settings.getIntList(Utils::ExternalQC::SettingsNames::spinFlipSites);
    for (unsigned long i = 0; i < spinFlipSites.size(); i++) {
      int site = spinFlipSites.at(i);
      if (i != spinFlipSites.size() - 1)
        out << site << ", ";
      else
        out << site;
    }
    out << std::endl;
    double ms = (settings.getInt(Scine::Utils::SettingsNames::spinMultiplicity) - 1.0) / 2.0;
    out << "FinalMs " << std::fixed << std::setprecision(1) << ms;
  }
  out << "\nend" << std::endl;
  if (settings.getBool(Utils::ExternalQC::SettingsNames::calculateMoessbauerParameter)) {
    if (Moessbauer::moessbauerNeededAndPossible(atoms, settings)) {
      // These settings correspond to the defaults recommended by the ORCA manual
      out << "%basis NewGTO 26 \"CP(PPP)\" end\nend";
      out << std::endl;
    }
    else {
      throw std::logic_error("You requested the calculation of 75-Fe Moessbauer parameters, but your structure does "
                             "not contain any Fe atoms!");
    }
  }

  // Write name of point charges file if it was set
  auto pointChargesFile = settings.getString(SettingsNames::pointChargesFile);
  if (!pointChargesFile.empty()) {
    boost::filesystem::path filePath(pointChargesFile);
    if (filePath.is_absolute()) {
      out << "%pointcharges \"" << pointChargesFile << "\"" << std::endl;
    }
    else {
      // The point charge file will not be in the orca calculation directory but one above.
      out << "%pointcharges \""
          << "../" << pointChargesFile << "\"" << std::endl;
    }
  }
}

void OrcaInputFileCreator::printTitle(std::ostream& out) {
  out << "# Orca calculation created by SCINE" << std::endl;
}

void OrcaInputFileCreator::printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings) {
  out << "*xyz " << settings.getInt(Scine::Utils::SettingsNames::molecularCharge) << " ";
  // If broken-symmetry calculation is enabled, the initial spin multiplicity must be written to that line
  if (settings.getBool(Utils::ExternalQC::SettingsNames::performBrokenSymmetryCalculation)) {
    out << settings.getInt(Scine::Utils::ExternalQC::SettingsNames::initialSpinMultiplicity) << std::endl;
  }
  else {
    out << settings.getInt(Scine::Utils::SettingsNames::spinMultiplicity) << std::endl;
  }
  for (const auto& a : atoms) {
    MolecularTrajectoryIO::writeXYZLine(out, a.getElementType(), a.getPosition());
  }
  out << "*" << std::endl;
  // The EPR NMR block needs to be below the coordinate block for the "all Fe" command to be recognized.
  if (Moessbauer::moessbauerNeededAndPossible(atoms, settings)) {
    // These settings correspond to the defaults recommended by the ORCA manual
    out << "%eprnmr nuclei = all Fe {rho, fgrad}";
    out << std::endl;
    out << "end";
  }
}

std::tuple<double, double> OrcaInputFileCreator::interpretAsUserDefinedImplicitSolvation(std::string solvent) {
  const std::string userDefined = "user_defined";
  auto startPosition = solvent.find("user_defined");
  // This should now look like (78.39 1.93)
  std::string epsilonAndProbeRadius = solvent.erase(startPosition, userDefined.length());
  // Remove brackets
  const bool firstIsABracket = epsilonAndProbeRadius[0] == '(';
  const bool lastIsABracket = epsilonAndProbeRadius[epsilonAndProbeRadius.length() - 1] == ')';
  if (!firstIsABracket || !lastIsABracket) {
    throw std::logic_error(
        "The solvent '" + solvent + "' is labeled as user defined but has a wrong format.\n" +
        "The format must be user_defined(<epsilon>, <probe_radius>). The brackets are misplaced or missing.");
  }
  epsilonAndProbeRadius.erase(0, 1);
  epsilonAndProbeRadius.erase(epsilonAndProbeRadius.length() - 1, 1);
  std::stringstream sstream(epsilonAndProbeRadius);
  double epsilon = -1.0;
  double probeRadius = -1.0;
  try {
    std::string epsilonString, probeRadiusString;
    std::getline(sstream, epsilonString, ',');
    std::getline(sstream, probeRadiusString, ',');
    epsilon = std::stod(epsilonString);
    probeRadius = std::stod(probeRadiusString);
  }
  catch (...) {
    throw std::logic_error("The solvent '" + solvent + "' is labeled as user defined but has a wrong format.\n" +
                           "The format must be user_defined(<epsilon>, <probe_radius>). Unable to convert the given\n" +
                           "dielectric constant or probe radius to a floating point number.");
  }
  if (sstream.rdbuf()->in_avail()) {
    throw std::logic_error("The solvent '" + solvent + "' is labeled as user defined but has a wrong format.\n" +
                           "The format must be user_defined(<epsilon>, <probe_radius>). The number of arguments in\n" +
                           "the brackets is larger than two.");
  }
  return std::make_tuple(epsilon, probeRadius);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
