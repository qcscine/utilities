/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaInputFileCreator.h"
#include "OrcaCalculatorSettings.h"
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
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
}

void OrcaInputFileCreator::createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                           const PropertyList& requiredProperties) {
  printCalculationType(out, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void OrcaInputFileCreator::printCalculationType(std::ostream& out, const Settings& settings,
                                                const PropertyList& requiredProperties) {
  out << "! " << settings.getString(Utils::SettingsNames::method) << " "
      << settings.getString(Utils::SettingsNames::basisSet) << std::endl;

  auto solvent = settings.getString(Utils::SettingsNames::solvent);
  if (!solvent.empty()) {
    out << "! CPCM(" << solvent << ")" << std::endl;
  }
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    out << "! EnGrad TightSCF" << std::endl;
  }
  if (requiredProperties.containsSubSet(Property::Hessian)) {
    out << "! AnFreq" << std::endl;
  }
  // Number of processes and memory per core
  int numProcs = settings.getInt(Utils::SettingsNames::externalProgramNProcs);
  out << "%maxcore " << settings.getInt(Utils::SettingsNames::externalProgramMemory) / numProcs << std::endl;
  if (numProcs != 1)
    out << "%pal\nnprocs " << numProcs << "\nend" << std::endl;

  // Print Mayer bond orders and/or Hirshfeld charges if required
  bool bondOrdersRequired = requiredProperties.containsSubSet(Utils::Property::BondOrderMatrix);
  bool atomicChargesRequired = requiredProperties.containsSubSet(Utils::Property::AtomicCharges);
  if (bondOrdersRequired && atomicChargesRequired)
    out << "%output\nprint[P_Mayer] 1\nprint[P_Hirshfeld] 1\nend" << std::endl;
  else if (bondOrdersRequired)
    out << "%output\nprint[P_Mayer] 1\nend" << std::endl;
  else if (atomicChargesRequired)
    out << "%output\nprint[P_Hirshfeld] 1\nend" << std::endl;
  // Set temperature if thermochemistry is to be calculated
  if (requiredProperties.containsSubSet(Utils::Property::Thermochemistry))
    out << "%freq\nTemp " << settings.getDouble(Utils::SettingsNames::temperature) << "\nend" << std::endl;
  // Scf convergence and Scf max iterations
  double tolE = settings.getDouble(Scine::Utils::SettingsNames::selfConsistanceCriterion);
  if (requiredProperties.containsSubSet(Property::Gradients) && tolE > 1e-8)
    tolE = 1e-8;
  out << "%SCF\nTolE " << tolE << std::endl
      << "MaxIter " << settings.getInt(Scine::Utils::SettingsNames::maxIterations) << "\nend" << std::endl;
  // Write name of point charges file if it was set
  auto pointChargesFile = settings.getString(SettingsNames::pointChargesFile);
  if (!pointChargesFile.empty())
    out << "%pointcharges \"" << pointChargesFile << "\"" << std::endl;
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
