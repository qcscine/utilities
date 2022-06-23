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
  std::string& calculationDirectory_;
  std::string& turbomoleExecutableBase_;
  std::string defineExecutableBase_ = "define";
  // maps the available solvents to the corresponding dielectric constants and radii of the solvent molecules in
  // Angstrom. Constants have been taken from the ADF documentation (https://www.scm.com/doc/ADF/Input/COSMO.html, last
  // visited Sept 06, 2021).
  std::unordered_map<std::string, std::pair<double, double>> availableSolventModels_ = {
      {"acetone", std::make_pair(21.7, 3.08)},  {"ammonia", std::make_pair(16.9, 2.24)},
      {"benzene", std::make_pair(2.3, 3.28)},   {"chloroform", std::make_pair(4.8, 3.17)},
      {"dmso", std::make_pair(46.7, 3.04)},     {"ethanol", std::make_pair(24.5, 2.85)},
      {"hexane", std::make_pair(1.88, 3.74)},   {"h2o", std::make_pair(78.4, 1.93)},
      {"methanol", std::make_pair(32.6, 2.53)}, {"nitrobenzene", std::make_pair(34.8, 3.44)},
      {"thf", std::make_pair(7.58, 3.18)},      {"toluene", std::make_pair(2.38, 3.48)},
      {"water", std::make_pair(78.4, 1.93)},    {"isopropanol", std::make_pair(19.9, 3.12)}};

  const std::vector<std::string> availableD3Params_ = {"D3", "D3BJ"};
  // Sets of damping parameters to aid SCF convergence. The old
  // Fock operator is added to the current one with a specific weight (first parameter). This weight is reduced in
  // specific steps (second parameter) until a minimum is reached (third parameter).
  std::unordered_map<std::string, std::vector<double>> availableDampingParameter_{
      {"low", {1.5, 0.05, 0.1}}, {"medium", {5.0, 0.1, 0.5}}, {"high", {8.5, 0.1, 0.5}}};
  TurbomoleFiles files_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_TURBOMOLEINPUTFILECREATOR_H
