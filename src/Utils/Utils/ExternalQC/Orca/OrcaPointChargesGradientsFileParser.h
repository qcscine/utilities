/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EXTERNALQC_ORCAPOINTCHARGESGRADIENTSFILEPARSER_H
#define UTILS_EXTERNALQC_ORCAPOINTCHARGESGRADIENTSFILEPARSER_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class OrcaPointChargesGradientsFileParser OrcaPointChargesGradientsFileParser.h
 * @brief Parses the point charges gradients from the ORCA .pcgrad file.
 */
class OrcaPointChargesGradientsFileParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the point charges gradients output file (*.pcgrad)
   */
  explicit OrcaPointChargesGradientsFileParser(const std::string& outputFileName);
  /**
   * @brief Parse the point charges gradients from the output file.
   * @return The point charges gradients.
   */
  GradientCollection getPointChargesGradients() const;

 private:
  // Extract the content of the file
  void extractContent(const std::string& filename);
  // The content of the file
  std::string content_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_ORCAPOINTCHARGESGRADIENTSFILEPARSER_H
