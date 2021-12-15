/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_TURBOMOLEINPUTFILECREATOR_H
#define UTILS_EXTERNALQC_TURBOMOLEINPUTFILECREATOR_H

#include <Utils/ExternalQC/Turbomole/TurbomoleFiles.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace Scine {
namespace Utils {

class AtomCollection;
class PropertyList;
class Settings;

namespace ExternalQC {

/**
 * @class TurbomoleInputFileCreator TurbomoleInputFileCreator.h
 * @brief This class creates Turbomole input files.
 */
class TurbomoleInputFileCreator {
 public:
  TurbomoleInputFileCreator(std::string& calculationDirectory, std::string& turbomoleExecutableBase, TurbomoleFiles& files);
  /**
   * @brief Create the Turbomole input file with the filename 'filename'
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  void createInputFiles(const AtomCollection& atoms, const Settings& settings);

 private:
  // Writes the coord file
  void writeCoordFile(const AtomCollection& atoms);
  // prepares the stdin for define
  void prepareDefineSession(const Settings& settings, const AtomCollection& atoms);
  // execution of preprocessing tool "define"
  void runDefine();
  // add different settings to control file if required
  void checkAndUpdateControlFile(const Settings& settings);
  // execution of "cosmoprep"
  void addSolvation(const Settings& settings);
  // helper function to check if charge and multiplicity are valid
  void checkValidityOfChargeAndMultiplicity(const Settings& settings, const AtomCollection& atoms);
  std::string& calculationDirectory_;
  std::string& turbomoleExecutableBase_;
  std::string defineExecutableBase_ = "define";
  // maps the available solvents to the corresponding dielectric constant
  // Constants have been taken from CRC Handbook of Chemistry and Physics, 101st edition,
  // section "Permittivity (Dielectric Constant) of Liquids"
  std::unordered_map<std::string, double> availableSolventModels_ = {
      {"acetone", 21.01}, {"ammonia", 16.61}, {"benzene", 2.2825}, {"chloroform", 4.8069}, {"dmso", 47.24},
      {"ethanol", 25.3},  {"hexane", 1.8865}, {"h2o", 80.1},       {"methanol", 33.0},     {"nitrobenzene", 35.6},
      {"thf", 7.52},      {"toluene", 2.379}, {"water", 80.1}};
  const std::vector<std::string> availableD3Params_ = {"D3", "D3BJ"};
  TurbomoleFiles files_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_TURBOMOLEINPUTFILECREATOR_H
