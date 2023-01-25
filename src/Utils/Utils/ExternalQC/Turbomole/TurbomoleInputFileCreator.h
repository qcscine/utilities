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
  // interpret user defined solvation settings
  void interpretAsUserDefinedImplicitSolvation(std::string solvent, double& epsilon, double& probeRadius);
  // correction for the define input if elements are present with multiple EHT parameter sets.
  std::string getMultipleEHTParameterCorrection(const AtomCollection& atoms);

  std::string& calculationDirectory_;
  std::string& turbomoleExecutableBase_;
  std::string defineExecutableBase_ = "define";
  // maps the available solvents to the corresponding dielectric constants and radii of the solvent molecules in
  // Angstrom. Constants have been taken from the ADF documentation (https://www.scm.com/doc/ADF/Input/COSMO.html, last
  // visited March 23, 2022).
  std::unordered_map<std::string, std::pair<double, double>> availableSolventModels_ = {
      {"aceticacid", std::make_pair(6.19, 2.83)},
      {"acetonitrile", std::make_pair(37.5, 2.76)},
      {"aniline", std::make_pair(6.8, 3.31)},
      {"benzylalcohol", std::make_pair(13.1, 3.45)},
      {"bromoform", std::make_pair(4.3, 3.26)},
      {"butanol", std::make_pair(17.5, 3.31)},
      {"isobutanol", std::make_pair(17.9, 3.33)},
      {"tertbutanol", std::make_pair(12.4, 3.35)},
      {"carbondisulfide", std::make_pair(2.6, 2.88)},
      {"carbontetrachloride", std::make_pair(2.2, 3.37)},
      {"cyclohexane", std::make_pair(2.0, 3.5)},
      {"cyclohexanone", std::make_pair(15.0, 3.46)},
      {"dichlorobenzene", std::make_pair(9.8, 3.54)},
      {"diethylether", std::make_pair(4.34, 3.46)},
      {"dioxane", std::make_pair(2.2, 3.24)},
      {"dmfa", std::make_pair(37.0, 3.13)},
      {"ethylacetate", std::make_pair(6.02, 3.39)},
      {"dichloroethane", std::make_pair(10.66, 3.15)},
      {"ethyleneglycol", std::make_pair(37.7, 2.81)},
      {"formicacid", std::make_pair(58.5, 2.47)},
      {"acetone", std::make_pair(20.7, 3.08)},
      {"ammonia", std::make_pair(16.9, 2.24)},
      {"benzene", std::make_pair(2.3, 3.28)},
      {"chloroform", std::make_pair(4.8, 3.17)},
      {"dmso", std::make_pair(46.7, 3.04)},
      {"ethanol", std::make_pair(24.55, 2.85)},
      {"hexane", std::make_pair(1.88, 3.74)},
      {"h2o", std::make_pair(78.39, 1.93)},
      {"methanol", std::make_pair(32.6, 2.53)},
      {"nitrobenzene", std::make_pair(34.8, 3.44)},
      {"thf", std::make_pair(7.58, 3.18)},
      {"toluene", std::make_pair(2.38, 3.48)},
      {"water", std::make_pair(78.39, 1.93)},
      {"isopropanol", std::make_pair(19.9, 3.12)}};

  const std::vector<std::string> availableDispersionParams_ = {"D3", "D3BJ", "D4"};

  TurbomoleFiles files_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_TURBOMOLEINPUTFILECREATOR_H
