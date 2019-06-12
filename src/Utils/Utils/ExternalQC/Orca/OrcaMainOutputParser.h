/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H
#define UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class OrcaMainOutputParser OrcaMainOutputParser.h
 * @brief This class parses information out of the main ORCA output file.
 */
class OrcaMainOutputParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the output file.
   */
  explicit OrcaMainOutputParser(const std::string& outputFileName);

  /**
   * @brief Parse the energy from the ORCA output.
   * @return The energy as a double.
   */
  double getEnergy() const;
  /**
   * @brief Parse the number of atoms from the ORCA output.
   * @return The number of atoms.
   */
  int getNumberAtoms() const;
  /**
   * @brief Parse the gradients from the ORCA output.
   * @return GradientCollection
   */
  GradientCollection getGradients() const;

 private:
  void extractContent(const std::string& filename);

  std::string content_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H