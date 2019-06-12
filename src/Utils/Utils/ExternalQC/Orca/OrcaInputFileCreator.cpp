/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

void InputFileCreator::createInputFile(const std::string& filename, const AtomCollection& atoms,
                                       const Settings& settings, const PropertyList& requiredProperties) {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, atoms, settings, requiredProperties);
}

void InputFileCreator::createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                       const PropertyList& requiredProperties) {
  printCalculationType(out, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void InputFileCreator::printCalculationType(std::ostream& out, const Settings& settings, const PropertyList& requiredProperties) {
  out << "! " << settings.getString(SettingsNames::orcaMethod) << std::endl;
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    out << "! EnGrad" << std::endl;
  }
  if (requiredProperties.containsSubSet(Property::Hessian)) {
    out << "! AnFreq" << std::endl;
  }
  // Number of processes
  int numProcs = settings.getInt(SettingsNames::orcaNumProcs);
  if (numProcs != 1)
    out << "%pal\nnprocs " << numProcs << "\nend" << std::endl;
  // SCF convergence and SCF max iterations
  out << "%SCF\nsthresh " << settings.getDouble(SettingsNames::selfConsistanceCriterion) << std::endl
      << "MaxIter " << settings.getInt(SettingsNames::maxIterations) << "\nend" << std::endl;
}

void InputFileCreator::printTitle(std::ostream& out) {
  out << "# Orca calculation created by SCINE" << std::endl;
}

void InputFileCreator::printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings) {
  out << "*xyz " << settings.getInt(SettingsNames::molecularCharge) << " "
      << settings.getInt(SettingsNames::spinMultiplicity) << std::endl;
  for (const auto& a : atoms) {
    MolecularTrajectoryIO::writeXYZLine(out, a.getElementType(), a.getPosition());
  }
  out << "*" << std::endl;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
