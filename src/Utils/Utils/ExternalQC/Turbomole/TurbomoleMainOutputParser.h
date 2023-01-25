/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_TURBOMOLEMAINOUTPUTPARSER_H
#define UTILS_EXTERNALQC_TURBOMOLEMAINOUTPUTPARSER_H

#include <Utils/ExternalQC/Turbomole/TurbomoleFiles.h>
#include <Utils/Typenames.h>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
class BondOrderCollection;
namespace ExternalQC {

/**
 * @class TurbomoleMainOutputParser TurbomoleMainOutputParser.h
 * @brief This class parses information out of the main Turbomole output file.
 */
class TurbomoleMainOutputParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the output file.
   */
  explicit TurbomoleMainOutputParser(TurbomoleFiles& files);
  /**
   * @brief Parse the Turbomole output for errors, it is expected to do nothing if no error is present
   * @throws OutputFileParsingError if error is present
   */
  void checkForErrors() const;
  /**
   * @brief Parse the energy from the Turbomole output.
   * @return The energy as a double.
   */
  double getEnergy() const;
  /**
   * @brief Parse the energy of the excited state specified from the Turbomole output.
   * @param state The index of the excited state queried (the first excited state has index 1).
   * @return The energy of the excited state as a double.
   */
  double getExcitedStateEnergy(int state) const;
  /**
   * @brief Parse the Mayer bond orders from the Turbomole output.
   * @return The Utils::BondOrderCollection.
   */
  Utils::BondOrderCollection getBondOrders() const;
  /**
   * @brief Parse the Loewdin charges from the Turbomole output.
   * @return The Loewdin charges.
   */
  std::vector<double> getLoewdinCharges() const;
  /**
   * @brief Parse the number of atoms from the Turbomole output.
   * @return The number of atoms.
   */
  int getNumberAtoms() const;
  /**
   * @brief Parse the number of environment point charges from the Turbomole input with a nonzero charge.
   * @return The number of point charges.
   */
  int getNumberOfNonZeroPointCharges() const;
  /**
   * @brief Parse the gradients from the Turbomole output.
   * @return GradientCollection
   */
  GradientCollection getGradients() const;
  /**
   * @brief Parse the gradients with respect to coordinates of point charges from the Turbomole output.
   * @return GradientCollection
   */
  GradientCollection getPointChargesGradients() const;
  /**
   * @brief Parses the hessian matrix from the Turbomole output.
   * @return HessianMatrix
   */
  HessianMatrix getHessian() const;
  int getSymmetryNumber() const;

 private:
  void extractContent(const std::string& filename);

  std::string content_;
  TurbomoleFiles files_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_TURBOMOLEMAINOUTPUTPARSER_H