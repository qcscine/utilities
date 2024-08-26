/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TURBOMOLEHELPER_H
#define TURBOMOLEHELPER_H

#include <algorithm>
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
  ///@brief List of all basis set notations supported by Turbomole. This is an optimistic list, because these
  /// are all basis sets for carbon in version 7.4.1 of Turbomole.
  /// in case the basis set is not supported for the actual element, we will detect the incorrect basis set written
  /// to the control file and throw an exception.
  const std::array<std::string, 130> basisSetNotations_ = {"CBSB5",
                                                           "SV",
                                                           "SVP",
                                                           "SV(P)",
                                                           "def-SVP",
                                                           "def-SV(P)",
                                                           "def2-SVP",
                                                           "dhf-SVP",
                                                           "dhf-SVP-2c",
                                                           "def2-SV(P)",
                                                           "dhf-SV(P)",
                                                           "dhf-SV(P)-2c",
                                                           "DZ",
                                                           "DZP",
                                                           "TZ",
                                                           "TZP",
                                                           "TZV",
                                                           "TZVP",
                                                           "def-TZVP",
                                                           "TZVPP",
                                                           "def-TZVPP",
                                                           "def2-TZVP",
                                                           "dhf-TZVP",
                                                           "dhf-TZVP-2c",
                                                           "def2-TZVPP",
                                                           "dhf-TZVPP",
                                                           "dhf-TZVPP-2c",
                                                           "TZVPPP",
                                                           "QZV",
                                                           "def-QZV",
                                                           "def2-QZV",
                                                           "QZVP",
                                                           "def-QZVP",
                                                           "def2-QZVP",
                                                           "dhf-QZVP",
                                                           "dhf-QZVP-2c",
                                                           "QZVPP",
                                                           "def-QZVPP",
                                                           "def2-QZVPP",
                                                           "dhf-QZVPP",
                                                           "dhf-QZVPP-2c",
                                                           "(7s3p)[3s2p]",
                                                           "(7s3p)[4s2p]",
                                                           "(9s5p)[5s3p]",
                                                           "(11s7p)[6s4p]",
                                                           "(13s8p)[8s5p]",
                                                           "minix",
                                                           "sto-3g hondo",
                                                           "sto-3g",
                                                           "4-31g hondo",
                                                           "4-31g",
                                                           "6-31g hondo",
                                                           "6-31g",
                                                           "3-21g hondo",
                                                           "3-21g",
                                                           "dz hondo",
                                                           "dz",
                                                           "dzp hondo",
                                                           "tz hondo",
                                                           "tz",
                                                           "10s6p-dun",
                                                           "tzp hondo",
                                                           "tzp",
                                                           "10s6p1d-dun",
                                                           "tz2p hondo",
                                                           "tz2p",
                                                           "10s6p2d-dun",
                                                           "ecp-2-sdf",
                                                           "ecp-02-sdf",
                                                           "2-sdf 2s2p",
                                                           "2-sdf 2s2p1d",
                                                           "6-31G",
                                                           "6-31G*",
                                                           "6-31G**",
                                                           "6-311G",
                                                           "6-311G*",
                                                           "6-311G**",
                                                           "6-311++G**",
                                                           "6-311G(2df,2pd)",
                                                           "cc-pVDZ",
                                                           "aug-cc-pVDZ",
                                                           "YP-aug-cc-pVDZ",
                                                           "d-aug-cc-pVDZ",
                                                           "cc-pwCVDZ",
                                                           "aug-cc-pwCVDZ",
                                                           "cc-pVTZ",
                                                           "cc-pVTZ-kernel",
                                                           "aug-cc-pVTZ",
                                                           "YP-aug-cc-pVTZ",
                                                           "d-aug-cc-pVTZ",
                                                           "d-aug-cc-pVTZ-oep",
                                                           "cc-pwCVTZ",
                                                           "aug-cc-pwCVTZ",
                                                           "cc-pVQZ",
                                                           "aug-cc-pVQZ",
                                                           "YP-aug-cc-pVQZ",
                                                           "d-aug-cc-pVQZ",
                                                           "cc-pwCVQZ",
                                                           "aug-cc-pwCVQZ",
                                                           "cc-pV5Z",
                                                           "aug-cc-pV5Z",
                                                           "YP-aug-cc-pV5Z",
                                                           "d-aug-cc-pV5Z",
                                                           "cc-pwCV5Z",
                                                           "aug-cc-pwCV5Z",
                                                           "cc-pV6Z",
                                                           "aug-cc-pV6Z",
                                                           "r12",
                                                           "cc-pVDZ-F12",
                                                           "cc-pVTZ-F12",
                                                           "cc-pVQZ-F12",
                                                           "def2-SVPD",
                                                           "def2-TZVPPD",
                                                           "def2-TZVPD",
                                                           "def2-QZVPPD",
                                                           "def2-QZVPD",
                                                           "pob-TZVP",
                                                           "x2c-SV(P)all",
                                                           "x2c-SV(P)all-s",
                                                           "x2c-SV(P)all-2c",
                                                           "x2c-SVPall",
                                                           "x2c-SVPall-s",
                                                           "x2c-SVPall-2c",
                                                           "x2c-TZVPall",
                                                           "x2c-TZVPall-s",
                                                           "x2c-TZVPall-2c",
                                                           "x2c-TZVPPall",
                                                           "x2c-TZVPPall-s",
                                                           "x2c-TZVPPall-2c",
                                                           "ANO-RCC-unc"};
  // generate another array with the same content but all lower-case with code
  const std::array<std::string, 130> basisSetNotationsLower_ = createLowerCaseNotations<130>(basisSetNotations_);

  template<unsigned int N>
  inline std::array<std::string, N> createLowerCaseNotations(const std::array<std::string, N>& contentList) {
    std::array<std::string, N> lowercaseArray;
    for (unsigned int i = 0; i < N; ++i) {
      lowercaseArray[i] = contentList[i];
      std::transform(lowercaseArray[i].begin(), lowercaseArray[i].end(), lowercaseArray[i].begin(), ::tolower);
    }
    return lowercaseArray;
  }
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // TURBOMOLEHELPER_H
