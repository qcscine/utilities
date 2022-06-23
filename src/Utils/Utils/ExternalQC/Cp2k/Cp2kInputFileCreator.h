/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H
#define UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include "Utils/Settings.h"
#include <map>
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
  /**
   * @brief Constructor of the writer.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   * @param isDft Whether a DFT input file wants to be written opposed to one for semiempirics
   */
  explicit Cp2kInputFileCreator(const AtomCollection& atoms, const Settings& settings,
                                const PropertyList& requiredProperties, bool isDft);
  /**
   * @brief Create the CP2K input file with the filename 'filename'
   * @param filename Name of the input file.
   * @param projectName Name of the project and basename of most output files
   */
  void createInputFile(const std::string& filename, const std::string& projectName) const;
  /**
   * @brief Write the CP2K input into the given outstream
   * @param out std::ostream.
   * @param projectName Name of the project and basename of most output files
   */
  void createInputFile(std::ostream& out, const std::string& projectName) const;

 private:
  // main sections
  void printGlobal(std::ostream& out, const std::string& projectName) const;
  void printForceEval(std::ostream& out) const;
  // structure related
  void printSubsys(std::ostream& out) const;
  void printCell(std::ostream& out) const;
  void printCoords(std::ostream& out) const;
  void printBasis(std::ostream& out) const;
  // calculation related
  void printDftInput(std::ostream& out) const;
  void printElectronicStructureBasics(std::ostream& out) const;
  void printFunctional(std::ostream& out) const;
  void printDispersionCorrection(std::ostream& out) const;
  void printDispersionCorrection(std::ostream& out, std::pair<std::string, std::string> methodInputs) const;
  void printSemiempiricalMethod(std::ostream& out) const;
  std::string determineCp2kSpinMode(SpinMode spinMode, int multiplicity) const;
  void printScfInput(std::ostream& out) const;
  void printPoissonSolver(std::ostream& out) const;
  void printGridInput(std::ostream& out) const;
  void printMatrixPrint(std::ostream& out) const;

  const std::vector<std::string> _availableZetaStrings = {"SZV", "DZVP", "TZVP", "TZV2P", "TZV2PX"};
  const std::vector<ElementType> _highQualityBasisElements = {ElementType::H, ElementType::C, ElementType::N,
                                                              ElementType::O, ElementType::F, ElementType::Si,
                                                              ElementType::P, ElementType::S, ElementType::Cl};
  const std::map<std::string, std::string> _dispersionSettingMap = {{"D3BJ", "DFTD3(BJ)"}, {"D3", "DFTD3"},
                                                                    {"D2", "DFTD2"},       {"DRSLL", "DRSLL"},
                                                                    {"LMKLL", "LMKLL"},    {"RVV10", "RVV10"}};
  AtomCollection _atoms;
  Settings _settings;
  PropertyList _requiredProperties;
  bool _isDft;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_CP2KINPUTFILECREATOR_H
