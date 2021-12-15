/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TURBOMOLEHELPER_H
#define TURBOMOLEHELPER_H

#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {
struct Filenames;

class TurbomoleHelper {
 public:
  TurbomoleHelper(std::string& calculationDirectory, std::string& turbomoleExecutableBase);
  /**
   * @brief Checks if a job was successful.
   */
  bool jobWasSuccessful(std::istream& out);
  /**
   * @brief Helper function to execute any Turbomole process.
   *
   * @param binaryName The name of the binary.
   * @param outputToFile Whether the output should be printed to a file.
   */
  void execute(std::string binaryName, bool outputToFile);
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

 private:
  std::string& calculationDirectory_;
  std::string& turbomoleExecutableBase_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // TURBOMOLEHELPER_H