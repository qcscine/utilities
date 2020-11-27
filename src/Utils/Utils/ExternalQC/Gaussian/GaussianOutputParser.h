/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_GAUSSIANOUTPUTPARSER_H
#define UTILS_EXTERNALQC_GAUSSIANOUTPUTPARSER_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class GaussianOutputParser GaussianOutputParser.h
 * @brief This class parses information out of the Gaussian output file.
 */
class GaussianOutputParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the output file.
   */
  explicit GaussianOutputParser(const std::string& outputFileName);

  /**
   * @brief Parse the energy from the Gaussian output.
   * @return The energy as a double.
   */
  double getEnergy() const;
  /**
   * @brief Parse the number of atoms from the Gaussian output.
   * @return The number of atoms.
   */
  int getNumberAtoms() const;
  /**
   * @brief Parse the gradients from the Gaussian output.
   * @return GradientCollection
   */
  GradientCollection getGradients() const;
  /**
   * @brief Parse CM5 charges from the Gaussian output.
   * @return The CM5 charges.
   */
  std::vector<double> getCM5Charges() const;

 private:
  void extractContent(const std::string& filename);

  std::string content_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_GAUSSIANOUTPUTPARSER_H