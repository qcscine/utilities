/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H
#define UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H

#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include <ostream>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {

class AtomCollection;
class PropertyList;
class Settings;

namespace ExternalQC {

/**
 * @class Cp2kInputFileCreator Cp2kInputFileCreator.h
 * @brief This class creates CP2K input files.
 */
class Cp2kInputFileCreator {
 public:
  Cp2kInputFileCreator() = default;
  /**
   * @brief Create the CP2K input file with the filename 'filename'
   * @param filename Name of the input file.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  void createInputFile(const std::string& filename, const AtomCollection& atoms, const Settings& settings,
                       const PropertyList& requiredProperties, const std::string& basename) const;
  /**
   * @brief Create the CP2K input file with the filename 'filename'
   * @param out std::ostream.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  void createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                       const PropertyList& requiredProperties, const std::string& basename) const;

 private:
  // main sections
  void printGlobal(std::ostream& out, const PropertyList& requiredProperties, const std::string& basename) const;
  void printForceEval(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                      const PropertyList& requiredProperties) const;
  // structure related
  void printSubsys(std::ostream& out, const AtomCollection& atoms, const Settings& settings) const;
  void printCell(std::ostream& out, const Settings& settings) const;
  void printCoords(std::ostream& out, const AtomCollection& atoms) const;
  void printBasis(std::ostream& out, const AtomCollection& atoms, const Settings& settings) const;
  // calculation related
  void printDftInput(std::ostream& out, const Settings& settings, const PropertyList& requiredProperties) const;
  void printElectronicStructureBasics(std::ostream& out, const Settings& settings) const;
  void printDispersionCorrection(std::ostream& out, const Settings& settings) const;
  std::string determineCp2kSpinMode(SpinMode spinMode, int multiplicity) const;
  void printScfInput(std::ostream& out, const Settings& settings) const;
  void printPoissonSolver(std::ostream& out, const Settings& settings) const;
  void printGridInput(std::ostream& out, const Settings& settings) const;
  void printMatrixPrint(std::ostream& out, const Settings& settings, const PropertyList& requiredProperties) const;

  const std::vector<std::string> _availableZetaStrings = {"SZV", "DZVP", "TZVP", "TZV2P", "TZV2PX"};
  const std::vector<ElementType> _highQualityBasisElements = {ElementType::H, ElementType::C, ElementType::N,
                                                              ElementType::O, ElementType::F, ElementType::Si,
                                                              ElementType::P, ElementType::S, ElementType::Cl};
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H
