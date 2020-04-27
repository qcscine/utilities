/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GaussianInputFileCreator.h"
#include "GaussianCalculatorSettings.h"
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <fstream>

namespace Scine {
namespace Utils {
namespace ExternalQC {

void GaussianInputFileCreator::createInputFile(const std::string& filename, const AtomCollection& atoms,
                                               const Settings& settings, const PropertyList& requiredProperties) {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, atoms, settings, requiredProperties);
}

void GaussianInputFileCreator::createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                               const PropertyList& requiredProperties) {
  printCalculationType(out, settings, requiredProperties);
  printTitle(out);
  printStructure(out, atoms, settings);
}

void GaussianInputFileCreator::printCalculationType(std::ostream& out, const Settings& settings,
                                                    const PropertyList& requiredProperties) {
  out << "%NProcShared=" << settings.getInt(SettingsNames::gaussianNumProcs) << std::endl;
  out << "%Mem=" << settings.getInt(SettingsNames::externalQCMemory) << "MB" << std::endl;

  out << "# " << settings.getString(Utils::SettingsNames::method) << "/" << settings.getString(Utils::SettingsNames::basisSet);
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
