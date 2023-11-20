/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TURBOMOLEHELPER_H
#define TURBOMOLEHELPER_H

#include <string>
#include <unordered_map>

namespace Scine {
namespace Utils {
namespace ExternalQC {
struct Filenames;

class TurbomoleHelper {
 public:
  TurbomoleHelper(std::string& calculationDirectory, std::string& turbomoleExecutableBase);
  /**
   * @brief Checks if a job was successful.
   *
   * @param out The name of the stream.
   * @param phrase The phrase to look for in the stream.
   */
  bool jobWasSuccessful(std::istream& out, std::string phrase = "(ended normally)");
  /**
   * @brief Helper function to execute any Turbomole process.
   *
   * @param binaryName The name of the binary.
   * @param outputToFile Whether the output should be printed to a file.
   */
  void execute(std::string binaryName, bool outputToFile, std::string outputFileName = "");
  /**
   * @brief Helper function to execute a Turbomole process that requires a stdin.
   */
  void execute(std::string binaryName, std::string stdInFile);
  /**
   * @brief Emptys a file if its present before the calculation starts.
   */
  void emptyFile(std::string file);
  /**
   * @brief Converts lower-case-only basis set strings to the case-sensitive file format required by Turbomole.
   * @note  Supports currently only Ahlrich's basis sets, Dunning basis set and Pople-style basis sets.
   */
  void mapBasisSetToTurbomoleStringRepresentation(std::string& basisSetString);
  /**
   * @brief Converts the DFT functional string to the turbomole input format if necessary (some DFT functional require a
   * "-" delimiter, e.g. b3-lyp).
   */
  void mapDftFunctionalToTurbomoleStringRepresentation(std::string& functionalString);

 private:
  std::string& calculationDirectory_;
  std::string& turbomoleExecutableBase_;
  std::unordered_map<std::string, std::string> correctedDftFunctionals_ = {{"b3lyp", "b3-lyp"},
                                                                           {"svwn", "s-vwn"},
                                                                           {"svwn_gaussian", "s-vwn_Gaussian"},
                                                                           {"b3lyp_gaussian", "b3-lyp_Gaussian"},
                                                                           {"blyp", "b-lyp"},
                                                                           {"m062x", "m06-2x"},
                                                                           {"b2plyp", "b2-plyp"},
                                                                           {"camb3lyp", "cam-b3lyp"},
                                                                           {"bvwn", "b-vwn"},
                                                                           {"bp", "b-p"},
                                                                           {"bp86", "b-p"},
                                                                           {"bhlyp", "bh-lyp"},
                                                                           {"m06l", "m06-l"},
                                                                           {"b97d", "b97-d"},
                                                                           {"pbeh3c", "pbeh-3c"},
                                                                           {"b973c", "b97-3c"}};
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // TURBOMOLEHELPER_H
